# Per-contingency reduced networks for a preventive SC-DCOPF, designed by
# adversarial dispatch generation. The formulation is in
# Adversarial_Injection_Generation_for_Safe_Network_Reduction.pdf.
#
# Network 0 is the full network, network c the full network without line c.
# Each gets its own merge mask. A mask is accepted when
#   - the reduced network accepts every hour's SC-DCOPF optimum (the witness), and
#   - no near-economic dispatch it accepts overloads a line of the network it
#     replaces (for an outage: among dispatches the intact network allows).
# Standalone: it does not use the greedy, KKT or proxy code.

module Adversarial

using LinearAlgebra, JuMP, Gurobi, Printf

export Net, net_from_base, with_ratings, add_line, clusters, ptdf, outages_of, lodf, post_ptdf,
       injection, sc_dcopf, adversary_set, limit!, unlimit!, reach, screen, check, design,
       exact_limits, bridges, unmerge, worst_loading

const GRB = Ref{Any}(nothing)
grb() = (isnothing(GRB[]) && (GRB[] = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))); GRB[])

function lp(; threads=1)
    m = Model(() -> Gurobi.Optimizer(grb()))
    set_silent(m)
    set_optimizer_attribute(m, "Threads", threads)
    set_optimizer_attribute(m, "FeasibilityTol", 1e-9)
    set_optimizer_attribute(m, "OptimalityTol", 1e-9)
    return m
end

# ---- network --------------------------------------------------------------------

"A DC network: lines (from, to, susceptance b, rating F) and generators."
struct Net
    N::Int
    L::Int
    from::Vector{Int}
    to::Vector{Int}
    b::Vector{Float64}
    F::Vector{Float64}
    ref::Int
    gen_bus::Vector{Int}
    pmin::Vector{Float64}
    pmax::Vector{Float64}
    c1::Vector{Float64}
    c0::Float64
end

"From the repo's case struct. The generator floor is relaxed to min(0, pmin) as in its OPFs."
net_from_base(base; relax_pmin=true) =
    Net(base.N, base.Ln, copy(base.Efrom), copy(base.Eto), copy(base.Dx), copy(base.frate),
        base.j0, copy(base.gen_bus), relax_pmin ? min.(0.0, base.pmin) : copy(base.pmin),
        copy(base.pmax), copy(base.c1), sum(base.c0))

"The same network with other ratings."
with_ratings(n::Net, F) = Net(n.N, n.L, n.from, n.to, n.b, F, n.ref, n.gen_bus, n.pmin, n.pmax, n.c1, n.c0)

add_line(n::Net, u, v, b, F) = Net(n.N, n.L + 1, [n.from; u], [n.to; v], [n.b; b], [n.F; F],
                                   n.ref, n.gen_bus, n.pmin, n.pmax, n.c1, n.c0)

"Cluster index of every bus, from the merged lines."
function clusters(n::Net, internal)
    parent = collect(1:n.N)
    root(x) = parent[x] == x ? x : (parent[x] = root(parent[x]))
    for l in 1:n.L
        internal[l] || continue
        a, b = root(n.from[l]), root(n.to[l])
        a == b || (parent[a] = b)
    end
    rep = [root(i) for i in 1:n.N]
    ids = Dict(r => k for (k, r) in enumerate(sort(unique(rep))))
    return [ids[r] for r in rep]
end

"""
Flow sensitivities with line `out` removed (0 = none) and the lines in
`internal` merged: flow on l = H[l, :] * p for any balanced p. Rows of merged
and removed lines are zero; `external` marks the others.
"""
function ptdf(n::Net; internal=falses(n.L), out=0)
    cl = clusters(n, internal)
    K = maximum(cl)
    ext = BitVector([l != out && cl[n.from[l]] != cl[n.to[l]] for l in 1:n.L])
    pos = zeros(Int, K)
    free = setdiff(1:K, [cl[n.ref]])
    pos[free] = 1:length(free)
    B = zeros(length(free), length(free))
    W = zeros(n.L, length(free))
    for l in findall(ext)
        a, b, y = pos[cl[n.from[l]]], pos[cl[n.to[l]]], n.b[l]
        a > 0 && (B[a, a] += y; W[l, a] += y)
        b > 0 && (B[b, b] += y; W[l, b] -= y)
        a > 0 && b > 0 && (B[a, b] -= y; B[b, a] -= y)
    end
    Hc = zeros(n.L, K)
    isempty(free) || (Hc[:, free] = transpose(lu(B) \ transpose(W)))
    return (H=Hc[:, cl], external=ext, cluster=cl, K=K)
end

"Lines whose outage keeps the network connected."
outages_of(n::Net) = [c for c in 1:n.L if maximum(clusters(n, [l != c for l in 1:n.L])) == 1]

"Line outage distribution factors for the listed outages (other columns zero)."
function lodf(n::Net, H, outs)
    Q = zeros(n.L, n.L)
    for c in outs
        a = H[:, n.from[c]] .- H[:, n.to[c]]
        Q[:, c] = a ./ (1 - a[c])
        Q[c, c] = -1.0
    end
    return Q
end

post_ptdf(H, Q, c) = c == 0 ? H : H .+ Q[:, c] * transpose(H[c, :])

function injection(n::Net, g, d)
    p = -Vector{Float64}(d)
    for k in eachindex(g)
        p[n.gen_bus[k]] += g[k]
    end
    return p
end

"Worst loading of the dispatch over the intact network and every listed outage."
function worst_loading(n::Net, H, Q, outs, g, d)
    f = H * injection(n, g, d)
    w = maximum(abs.(f) ./ n.F)
    for c in outs, l in 1:n.L
        l == c && continue
        w = max(w, abs(f[l] + Q[l, c] * f[c]) / n.F[l])
    end
    return w
end

# ---- full SC-DCOPF ----------------------------------------------------------------

"""
Preventive SC-DCOPF at every column of D by constraint generation: the
post-outage rows are true at every hour, so one pool is shared and only the
right-hand sides follow the demand. Returns cost and dispatch per hour (NaN
cost where infeasible) and the pool size.
"""
function sc_dcopf(n::Net, H, Q, outs, D)
    K, L, T = length(n.gen_bus), n.L, size(D, 2)
    G = H[:, n.gen_bus]
    m = lp()
    @variable(m, n.pmin[k] <= g[k=1:K] <= n.pmax[k])
    bal = @constraint(m, sum(g) == 0.0)
    up = @constraint(m, [l=1:L],  sum(G[l, k] * g[k] for k in 1:K) <= 0.0)
    dn = @constraint(m, [l=1:L], -sum(G[l, k] * g[k] for k in 1:K) <= 0.0)
    @objective(m, Min, sum(n.c1[k] * g[k] for k in 1:K) + n.c0)
    pool = Dict{Tuple{Int,Int},Any}()
    function rhs!(d)
        set_normalized_rhs(bal, sum(d))
        hd = H * d
        for l in 1:L
            set_normalized_rhs(up[l], n.F[l] + hd[l])
            set_normalized_rhs(dn[l], n.F[l] - hd[l])
        end
        for ((l, _), r) in pool
            w = dot(r.row, d)
            set_normalized_rhs(r.up, n.F[l] + w)
            set_normalized_rhs(r.dn, n.F[l] - w)
        end
    end
    Z, Gs = fill(NaN, T), zeros(K, T)
    for s in 1:T
        d = D[:, s]
        rhs!(d)
        while true
            optimize!(m)
            termination_status(m) == MOI.OPTIMAL || break
            f = H * injection(n, value.(g), d)
            added = 0
            for c in outs, l in 1:L
                (l == c || haskey(pool, (l, c))) && continue
                abs(f[l] + Q[l, c] * f[c]) > n.F[l] * (1 + 1e-7) || continue
                row = H[l, :] .+ Q[l, c] .* H[c, :]
                a, w = row[n.gen_bus], dot(row, d)
                pool[(l, c)] = (; row,
                    up=@constraint(m,  sum(a[k] * g[k] for k in 1:K) <= n.F[l] + w),
                    dn=@constraint(m, -sum(a[k] * g[k] for k in 1:K) <= n.F[l] - w))
                added += 1
            end
            if added == 0
                Z[s], Gs[:, s] = objective_value(m), value.(g)
                break
            end
        end
    end
    return Z, Gs, length(pool)
end

# ---- adversary set ------------------------------------------------------------------

"""
Dispatches an adversary may use: generator box, balance, demand anywhere in
the hull of D's columns, cost at most (1+x) times the interpolated SC-DCOPF
optimum Z and, when `intact` (the full PTDF) is given, feasible on the intact
network. One LP, re-solved with new objectives and reduced-network rows.
"""
function adversary_set(n::Net, D, Z, x; intact=nothing)
    K, T = length(n.gen_bus), size(D, 2)
    m = lp()
    @variable(m, n.pmin[k] <= g[k=1:K] <= n.pmax[k])
    @variable(m, lam[1:T] >= 0)
    @constraint(m, sum(lam) == 1)
    tot = vec(sum(D; dims=1))
    @constraint(m, sum(g) == sum(tot[t] * lam[t] for t in 1:T))
    isfinite(x) && @constraint(m, sum(n.c1[k] * g[k] for k in 1:K) + n.c0 <=
                                  (1 + x) * sum(Z[t] * lam[t] for t in 1:T))
    A = (; m, g, lam, D, net=n, rows=ConstraintRef[])
    isnothing(intact) || limit!(A, intact, trues(n.L); keep=true)
    return A
end

function flowexpr(A, row)
    a, w = row[A.net.gen_bus], transpose(A.D) * row
    return @expression(A.m, sum(a[k] * A.g[k] for k in eachindex(a)) -
                            sum(w[t] * A.lam[t] for t in eachindex(w)))
end

"""
Add |H_l p| <= F_l (or the given ratings) for the marked lines; an infinite
rating adds nothing. Rows added without `keep` go with `unlimit!`.
"""
function limit!(A, H, lines; keep=false, F=A.net.F)
    for l in findall(lines)
        isfinite(F[l]) || continue
        e = flowexpr(A, view(H, l, :))
        u = @constraint(A.m, e <= F[l])
        d = @constraint(A.m, -e <= F[l])
        keep || push!(A.rows, u, d)
    end
end
unlimit!(A) = (foreach(c -> delete(A.m, c), A.rows); empty!(A.rows))

"Max of sign * row' p over the set, and the injection attaining it (-Inf, nothing if empty)."
function reach(A, row, sign)
    @objective(A.m, Max, sign * flowexpr(A, row))
    optimize!(A.m)
    st = termination_status(A.m)
    st in (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED) && return -Inf, nothing
    st == MOI.OPTIMAL || error("adversary LP ended $st")
    return objective_value(A.m), injection(A.net, value.(A.g), A.D * value.(A.lam))
end

"""
Pairs (line, sign) of a network that some adversary can overload, with the
worst loading. `bound(l, sign)`, when given, is a cheap upper bound that clears
a pair without its own LP.
"""
function screen(A, Hc, lines; tau=1e-6, bound=nothing)
    crit = Tuple{Int,Int,Float64}[]
    for l in lines, s in (1, -1)
        !isnothing(bound) && bound(l, s) <= A.net.F[l] * (1 + tau) && continue
        v, _ = reach(A, view(Hc, l, :), s)
        v > A.net.F[l] * (1 + tau) && push!(crit, (l, s, v / A.net.F[l]))
    end
    return crit
end

# ---- master problem ---------------------------------------------------------------

"""
Bridges of the network without line `out`: lines whose removal disconnects it.
Merging a bridge changes no other line's flow, since the power beyond it has only
one way out. Iterative DFS with low-links; parallel lines are never bridges.
"""
function bridges(n::Net, out)
    adj = [Tuple{Int,Int}[] for _ in 1:n.N]          # (neighbour, line)
    for l in 1:n.L
        l == out && continue
        push!(adj[n.from[l]], (n.to[l], l))
        push!(adj[n.to[l]], (n.from[l], l))
    end
    disc, low = zeros(Int, n.N), zeros(Int, n.N)
    isbridge = falses(n.L)
    clock = 0
    for s in 1:n.N
        disc[s] > 0 && continue
        clock += 1
        disc[s] = low[s] = clock
        stack = [(s, 0, 1)]                           # (bus, line in, next neighbour)
        while !isempty(stack)
            u, lin, i = stack[end]
            if i <= length(adj[u])
                stack[end] = (u, lin, i + 1)
                v, l = adj[u][i]
                l == lin && continue
                if disc[v] == 0
                    clock += 1
                    disc[v] = low[v] = clock
                    push!(stack, (v, l, 1))
                else
                    low[u] = min(low[u], disc[v])
                end
            else
                pop!(stack)
                if !isempty(stack)
                    p = stack[end][1]
                    low[p] = min(low[p], low[u])
                    low[u] > disc[p] && (isbridge[lin] = true)
                end
            end
        end
    end
    return isbridge
end

"Unmerge the merged lines within `hops` of line `out`'s ends (hops counted on the full network)."
function unmerge(n::Net, mask, out, hops)
    adj = [Int[] for _ in 1:n.N]
    for l in 1:n.L
        push!(adj[n.from[l]], n.to[l])
        push!(adj[n.to[l]], n.from[l])
    end
    dist = fill(typemax(Int), n.N)
    frontier = unique([n.from[out], n.to[out]])
    dist[frontier] .= 0
    while !isempty(frontier)
        nxt = Int[]
        for u in frontier, v in adj[u]
            dist[v] == typemax(Int) || continue
            dist[v] = dist[u] + 1
            push!(nxt, v)
        end
        frontier = nxt
    end
    m = copy(mask)
    for l in 1:n.L
        m[l] && min(dist[n.from[l]], dist[n.to[l]]) < hops && (m[l] = false)
    end
    return m
end

# The master picks the lines to merge (z = 1) so the reduced network keeps as few
# lines as possible: maximize the merged lines. A line whose ends already share a
# cluster carries no flow either way, so counting it changes nothing physical.
# `keep` lines are fixed kept (the critical lines, which do the rejecting) and
# `merge` lines fixed merged (bridges). `cap` carries earlier rounds' bound while
# no neighbourhood is set; `near` keeps designs within a radius of a centre design.
# MIPFocus 1 and a gap: good designs fast, since the LP check certifies them anyway.
function master(n::Net, lines, Fr; threads=4, gap=0.05, keep=Int[], merge=Int[])
    m = Model(() -> Gurobi.Optimizer(grb()))
    set_silent(m)
    set_optimizer_attribute(m, "Threads", threads)
    # An adversary sits where some kept line is exactly at its rating, and the cut
    # asks for a little more. With the default 1e-6 tolerances that line would pass
    # as rejecting it and the cut would do nothing.
    set_optimizer_attribute(m, "FeasibilityTol", 1e-9)
    set_optimizer_attribute(m, "IntFeasTol", 1e-9)
    set_optimizer_attribute(m, "MIPFocus", 1)
    set_optimizer_attribute(m, "MIPGap", gap)
    # long runs: spill the branch-and-bound tree to disk past this many GB
    haskey(ENV, "ADV_NODEFILE_START") &&
        set_optimizer_attribute(m, "NodefileStart", parse(Float64, ENV["ADV_NODEFILE_START"]))
    @variable(m, z[lines], Bin)
    foreach(l -> fix(z[l], 0.0; force=true), keep)
    foreach(l -> fix(z[l], 1.0; force=true), merge)
    cap = @constraint(m, sum(z) <= length(lines))
    @objective(m, Max, sum(z))
    return (; m, z, lines, cap, Fr, near=Ref{Any}(nothing))
end

"Keep designs within `radius` line changes of `centre`; an infinite radius removes the limit."
function near!(M, centre, radius)
    isnothing(M.near[]) || (delete(M.m, M.near[]); M.near[] = nothing)
    isfinite(radius) || return
    M.near[] = @constraint(M.m, sum((centre[l] ? 1 - M.z[l] : M.z[l]) for l in M.lines) <= radius)
end

"""
Flows of the design at a fixed injection p: f on kept lines, free transfers t
inside clusters. With `limit` the kept lines carry their reduced ratings (a
witness); without, they carry none (an adversary the design must reject).

Bounds, so Gurobi can tighten the indicators: with positive susceptances a DC
flow has no cycle, so no line carries more than the positive injection P, and
routing each cluster's remainder over a tree needs at most (L+1)P per transfer.
"""
function flow_copy!(M, n::Net, p; limit, tau)
    m, z, lines = M.m, M.z, M.lines
    P = sum(max.(p, 0.0))
    fb = all(>(0), n.b) ? P : Inf
    th = @variable(m, [1:n.N])
    f = @variable(m, [lines], lower_bound=-fb, upper_bound=fb)
    t = @variable(m, [lines], lower_bound=-(n.L + 1) * P, upper_bound=(n.L + 1) * P)
    @constraint(m, th[n.ref] == 0)
    for l in lines
        @constraint(m, f[l] == n.b[l] * (th[n.from[l]] - th[n.to[l]]))
        @constraint(m, !z[l] => {t[l] == 0})
        if limit
            @constraint(m,  f[l] <= (M.Fr[l] + tau * n.F[l]) * (1 - z[l]))
            @constraint(m, -f[l] <= (M.Fr[l] + tau * n.F[l]) * (1 - z[l]))
        else
            @constraint(m, z[l] => {f[l] == 0})
        end
    end
    for i in 1:n.N
        @constraint(m, sum(f[l] + t[l] for l in lines if n.from[l] == i; init=0.0) -
                       sum(f[l] + t[l] for l in lines if n.to[l] == i; init=0.0) == p[i])
    end
    return f
end

"""
The design must reject the adversary injection. It overloads line l on the full
network, with flow v in direction sg; `keep` = Fr_l / F_l is the reduced rating's
share of the true one.
  :dominance  line l carries at least keep * v, so with the true flow above F its
              reduced flow is above keep * F = Fr: rejected by l itself
  :line       line l carries more than Fr_l
  :any        some kept candidate line carries more than its Fr (binaries)
"""
function reject!(M, n::Net, f, cand; tau, cut=:dominance, l=0, sg=1, v=0.0)
    cut === :dominance && return @constraint(M.m, sg * f[l] >= M.Fr[l] / n.F[l] * v)
    cut === :line && return @constraint(M.m, sg * f[l] >= M.Fr[l] + tau * n.F[l])
    y = @variable(M.m, [cand, [1, -1]], Bin)
    @constraint(M.m, sum(y) >= 1)
    for j in cand, s in (1, -1)
        @constraint(M.m, y[j, s] <= 1 - M.z[j])
        @constraint(M.m, y[j, s] => {s * f[j] >= M.Fr[j] + tau * n.F[j]})
    end
end

"""
Check a fixed mask of the network without line `out`, its kept lines rated
`rating`: the hours whose witness (columns of Pw) it rejects, and the worst
adversary's injection for every screened pair it fails. `worst` is the highest
loading any adversary reaches on the full network.
"""
function check(n::Net, out, Hc, A, crit, Pw, mask; tau=1e-6, rating=n.F)
    mp = ptdf(n; internal=mask, out)
    ext = findall(mp.external)
    Fw = mp.H * Pw
    cap = rating .+ tau .* n.F
    rejected = count(s -> any(abs(Fw[l, s]) > cap[l] for l in ext), axes(Pw, 2))
    # the worst hour of every kept line and sign the witnesses overload: accepting
    # those tends to accept the rest, since the accepted set is convex
    bad = Int[]
    for l in ext, sg in (1, -1)
        v, s = findmax(sg .* Fw[l, :])
        v > cap[l] && !(s in bad) && push!(bad, s)
    end
    limit!(A, mp.H, mp.external; F=rating)
    advs, worst = [], 0.0
    for (l, s, _) in crit
        v, p = reach(A, view(Hc, l, :), s)
        worst = max(worst, v / n.F[l])
        v > n.F[l] * (1 + tau) && push!(advs, (v / n.F[l], p, l, s, v))
    end
    unlimit!(A)
    sort!(advs; by=a -> -a[1])
    return (; bad, rejected, advs, worst, cluster=mp.cluster)
end

"""
Give the kept lines back as much rating as the check allows: r + a(F - r) with
the largest a in [0, 1] that still passes. Raising ratings only enlarges what
the reduced network accepts, so passing is monotone in a and bisection is exact.
"""
function raise_ratings(n::Net, out, Hc, A, crit, Pw, mask, rating; tau=1e-6, steps=12)
    at(a) = rating .+ a .* (n.F .- rating)
    ok(a) = (ck = check(n, out, Hc, A, crit, Pw, mask; tau, rating=at(a));
             isempty(ck.bad) && isempty(ck.advs))
    ok(1.0) && return at(1.0)
    lo, hi = 0.0, 1.0
    for _ in 1:steps
        mid = (lo + hi) / 2
        ok(mid) ? (lo = mid) : (hi = mid)
    end
    return at(lo)
end

"""
Limits for a finished mask, set by the critical lines alone. Each critical line
gets the highest rating at which it still rejects, by itself, every adversary
dispatch that loads the true line to its rating or beyond:

    r_l = min over its critical signs s of  min{ s phi_l(p) : p in A, s Hc_l p >= F_l }

less tau F_l; it may land above F_l. Every other line gets no limit (Inf):
screening showed it cannot overload for these dispatches, and a limit would only
reject good ones. Returns nothing when a critical line is merged, since some other
line must then do the rejecting. The result still has to pass `check`.
"""
function exact_limits(n::Net, out, Hc, A, crit, mask; tau=1e-6)
    mp = ptdf(n; internal=mask, out)
    rating = fill(Inf, n.L)
    for (l, s, _) in crit
        mp.external[l] || return nothing
        con = @constraint(A.m, s * flowexpr(A, view(Hc, l, :)) >= n.F[l])
        v, _ = reach(A, view(mp.H, l, :), -s)
        delete(A.m, con)
        isfinite(v) && (rating[l] = min(rating[l], -v - tau * n.F[l]))
    end
    return rating
end

"""
    design(n, out, Hc, A, crit, Pw; ...)

Reduced network for the network with line `out` removed (0 = none), keeping as
few lines as possible. `Hc` is that network's PTDF, `A` its adversary set, `crit`
its screened pairs and the columns of `Pw` the witnesses' injections. Kept lines
are rated (1 - derate) F.

Before any master solve: the critical lines are fixed kept and (with
`fix_radial`) every other bridge fixed merged. With a `start` design (the base
network's, for an outage) it tries that design, then the same with the merged
lines within `unmerge_hops` of the outage line unmerged; the first that passes
both checks is returned. Otherwise what they failed on goes into the master,
which then searches within `radius` line changes of the unmerged design,
doubling the radius whenever nothing is found there.

Each round solves the master, then checks the design: hours whose witness it
rejects become witness copies, and the worst adversary of every pair it fails
becomes a flow copy with a rejection cut. Stops when both checks pass, i.e. when
no adversary is left, then (with `raise`) gives ratings back as far as the check
allows. On a time or round limit it returns the unreduced network at full
ratings, which always passes when the witnesses are feasible at (1 - derate) F.
"""
function design(n::Net, out, Hc, A, crit, Pw; tau=1e-6, candidates=:all, cut=:dominance,
                derate=0.0, raise=true, witnesses_per_round=10, adversaries_per_round=3,
                master_time=20.0, time_limit=600.0, max_rounds=200, threads=4, gap=0.05,
                fix_radial=true, unmerge_hops=2, radius=20, start=nothing, log=nothing)
    lines = [l for l in 1:n.L if l != out]
    cand = candidates === :all ? lines : sort(unique(first.(crit)))
    Fr = (1 - derate) .* n.F
    keep = sort(unique(first.(crit)))
    merge = fix_radial ? [l for l in findall(bridges(n, out)) if !(l in keep)] : Int[]
    M = master(n, lines, Fr; threads, gap, keep, merge)
    started = time()
    rec = (rounds=Ref(0), adversaries=Ref(0), witnesses=Ref(0),
           master_s=Ref(0.0), check_s=Ref(0.0))
    finish(mask, status, rating=copy(n.F)) = (; mask, status, rating, rounds=rec.rounds[],
        adversaries=rec.adversaries[], witnesses=rec.witnesses[], master_s=rec.master_s[],
        check_s=rec.check_s[], seconds=time() - started)
    kept(ck) = count(l -> l != out && ck.cluster[n.from[l]] != ck.cluster[n.to[l]], 1:n.L)
    function accept(ck)
        closed = BitVector([l != out && ck.cluster[n.from[l]] == ck.cluster[n.to[l]] for l in 1:n.L])
        rating = copy(Fr)
        rating[closed] .= n.F[closed]
        raise && derate > 0 && (rating = raise_ratings(n, out, Hc, A, crit, Pw, closed, rating; tau))
        return closed, rating
    end
    function learn!(ck)
        for s in ck.bad[1:min(witnesses_per_round, end)]
            flow_copy!(M, n, Pw[:, s]; limit=true, tau)
            rec.witnesses[] += 1
        end
        for (_, p, l, sg, v) in ck.advs[1:min(adversaries_per_round, end)]
            reject!(M, n, flow_copy!(M, n, p; limit=false, tau), cand; tau, cut, l, sg, v)
            rec.adversaries[] += 1
        end
    end
    passes(ck) = isempty(ck.bad) && isempty(ck.advs)
    isnothing(log) || log(@sprintf("  fixed   : %d critical lines kept, %d bridges merged, %d of %d lines left to decide",
                                   length(keep), length(merge), length(lines) - length(keep) - length(merge), length(lines)))

    centre = nothing
    if !isnothing(start)
        s0 = BitVector([l != out && start[l] for l in 1:n.L])
        s0[keep] .= false
        s0[merge] .= true
        tries = [("start", s0)]
        out > 0 && unmerge_hops > 0 && push!(tries, ("unmerged", unmerge(n, s0, out, unmerge_hops)))
        for (label, mask) in tries
            t0 = time()
            ck = check(n, out, Hc, A, crit, Pw, mask; tau, rating=Fr)
            rec.check_s[] += time() - t0
            isnothing(log) || log(@sprintf("  %-8s: %3d buses, %3d lines, rejects %d witnesses, %d adversaries, worst loading %.4f",
                label, maximum(ck.cluster), kept(ck), ck.rejected, length(ck.advs), ck.worst))
            if passes(ck)
                closed, rating = accept(ck)
                return finish(closed, :certified, rating)
            end
            learn!(ck)
        end
        centre = tries[end][2]
        near!(M, centre, radius)
    end

    r = radius
    while true
        rec.rounds[] += 1
        left = time_limit - (time() - started)
        (left <= 0 || rec.rounds[] > max_rounds) && return finish(falses(n.L), :limit)
        set_time_limit_sec(M.m, min(master_time, left))
        # the centre as the start while the search is local; otherwise merging only the
        # fixed bridges, which changes no other flow and so meets every cut and witness
        local_search = !isnothing(M.near[])
        for l in lines
            set_start_value(M.z[l], local_search ? Float64(centre[l]) : Float64(l in merge))
        end
        optimize!(M.m)
        rec.master_s[] += solve_time(M.m)
        if primal_status(M.m) != MOI.FEASIBLE_POINT
            isnothing(M.near[]) && return finish(falses(n.L), :limit)
            r *= 2
            near!(M, centre, r >= length(lines) ? Inf : r)
            isnothing(log) || log("  nothing within the radius; widened to " *
                                  (isnothing(M.near[]) ? "the whole network" : string(r)))
            continue
        end
        mask = falses(n.L)
        for l in lines
            mask[l] = value(M.z[l]) > 0.5
        end
        master_s, master_status = solve_time(M.m), termination_status(M.m)
        isnothing(M.near[]) &&
            set_normalized_rhs(M.cap, min(normalized_rhs(M.cap), floor(objective_bound(M.m) + 1e-6)))

        t0 = time()
        ck = check(n, out, Hc, A, crit, Pw, mask; tau, rating=Fr)
        rec.check_s[] += time() - t0
        isnothing(log) || log(@sprintf("  round %3d: %3d buses, %3d lines (master %.2fs, %s), rejects %d witnesses, %d adversaries, worst loading %.4f",
            rec.rounds[], maximum(ck.cluster), kept(ck), master_s, master_status, ck.rejected,
            length(ck.advs), ck.worst))
        if passes(ck)
            closed, rating = accept(ck)
            return finish(closed, :certified, rating)
        end
        learn!(ck)
    end
end

end # module

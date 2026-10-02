# Exact SC-DCOPF methods to compare the reduced networks against. None reduces
# the network; each solves the full preventive SC-DCOPF.
#   ptdf_lazy        PTDF rows over the generators, added only when violated
#   btheta_lazy      B-theta network copies, a contingency's copy added only when violated
#   ptdf_decomposed  the decomposed PTDF of Alkhraijah et al. (PowerUp 2026): Kron
#                    reduction per area, boundary equivalent injections, lazy rows
# The two PTDF methods come in two forms. Substituted (flows=false): a limit row
# holds the whole flow expression, f_l + LODF f_c written out. Flow variables
# (flows=true, the paper's (4) and (9)): each line's base-case flow becomes a
# variable with one defining row, added the first time the line is needed, and a
# limit row holds just f_l, or f_l + LODF f_c.
# Each returns (; m, g, demand!, solve!) for Evaluate.solve_hours. The pools of
# lazy rows or copies are kept from one hour to the next.
# Needs adversarial/adversarial.jl loaded first (module Adversarial).

module Baselines

using LinearAlgebra, SparseArrays, Statistics, JuMP, Gurobi
using Main.Adversarial: Net, lp, injection, ptdf, clusters

export ptdf_lazy, btheta_lazy, ptdf_decomposed, partition, decomposed_gap,
       reduced_ptdf_lazy, reduced_btheta_lazy, compact_flows

const TOL = 1e-7
# a lazy solve past this time stops and marks its model unfinished (set by the runner)
const DEADLINE = Ref(Inf)

"Violated (line, outage) pairs of dispatch g at demand d; outage 0 = intact."
function violations(n::Net, H, Q, outs, g, d)
    f = H * injection(n, g, d)
    bad = Tuple{Int,Int}[]
    for l in 1:n.L
        abs(f[l]) > n.F[l] * (1 + TOL) && push!(bad, (l, 0))
    end
    for c in outs, l in 1:n.L
        l == c && continue
        abs(f[l] + Q[l, c] * f[c]) > n.F[l] * (1 + TOL) && push!(bad, (l, c))
    end
    return bad
end

function base_model(n::Net)
    m = lp()
    K = length(n.gen_bus)
    @variable(m, n.pmin[k] <= g[k=1:K] <= n.pmax[k])
    bal = @constraint(m, sum(g) == 0.0)
    @objective(m, Min, sum(n.c1[k] * g[k] for k in 1:K) + n.c0)
    return m, g, bal
end

# solve, add what the dispatch violates, repeat; returns solver seconds
function lazy_loop!(m, g, check, add!, unfinished)
    secs = 0.0
    while true
        time() > DEADLINE[] && (unfinished[] = true; return secs)
        optimize!(m)
        secs += solve_time(m)
        termination_status(m) == MOI.OPTIMAL || return secs
        bad = check(value.(g))
        isempty(bad) && return secs
        foreach(add!, bad)
    end
end

"""
Lazy flow limits over flow expressions: `flowexpr(l)` gives (expression, demand
weights) with f_l = expression - weights . d. See the header for the two forms.
"""
function lazy_limits(n::Net, H, Q, outs, m, flowexpr; flows)
    pool = Dict{Tuple{Int,Int},Any}()
    fv = Dict{Int,Any}()
    d = Ref(zeros(n.N))
    function flowvar(l)
        haskey(fv, l) && return fv[l].f
        ex, wd = flowexpr(l)
        f = @variable(m)
        fv[l] = (; f, wd, def=@constraint(m, f - ex == -dot(wd, d[])))
        return f
    end
    function add!((l, c))
        haskey(pool, (l, c)) && return
        if flows
            ex = c == 0 ? 1.0 * flowvar(l) : flowvar(l) + Q[l, c] * flowvar(c)
            pool[(l, c)] = (; wd=nothing, up=@constraint(m, ex <= n.F[l]), dn=@constraint(m, -ex <= n.F[l]))
        else
            ex, wd = flowexpr(l)
            if c != 0
                exc, wdc = flowexpr(c)
                ex, wd = ex + Q[l, c] * exc, wd + Q[l, c] * wdc
            end
            w = dot(wd, d[])
            pool[(l, c)] = (; wd, up=@constraint(m, ex <= n.F[l] + w), dn=@constraint(m, -ex <= n.F[l] - w))
        end
    end
    function update!(dd)
        d[] = dd
        for ((l, _), r) in pool
            isnothing(r.wd) && continue
            w = dot(r.wd, dd)
            set_normalized_rhs(r.up, n.F[l] + w)
            set_normalized_rhs(r.dn, n.F[l] - w)
        end
        foreach(x -> set_normalized_rhs(x.def, -dot(x.wd, dd)), values(fv))
    end
    check(gv) = filter(p -> !haskey(pool, p), violations(n, H, Q, outs, gv, d[]))
    return (; add!, update!, check)
end

"""
Compact rows with flow variables: the screened (line, outage) pairs, all present
from the start, over base-case flow variables (one defining row per line used).
"""
function compact_flows(n::Net, H, Q, outs, pairs)
    m, g, bal = base_model(n)
    K = length(n.gen_bus)
    fe(l) = (row = H[l, :]; (@expression(m, sum(row[n.gen_bus[k]] * g[k] for k in 1:K)), row))
    lim = lazy_limits(n, H, Q, outs, m, fe; flows=true)
    foreach(lim.add!, pairs)
    demand!(dd) = (set_normalized_rhs(bal, sum(dd)); lim.update!(dd))
    return (; m, g, demand!)
end

"PTDF model with every limit (intact and post-outage) added lazily."
function ptdf_lazy(n::Net, H, Q, outs; flows=false)
    m, g, bal = base_model(n)
    K = length(n.gen_bus)
    fe(l) = (row = H[l, :]; (@expression(m, sum(row[n.gen_bus[k]] * g[k] for k in 1:K)), row))
    lim = lazy_limits(n, H, Q, outs, m, fe; flows)
    demand!(dd) = (set_normalized_rhs(bal, sum(dd)); lim.update!(dd))
    unfinished = Ref(false)
    solve!() = lazy_loop!(m, g, lim.check, lim.add!, unfinished)
    return (; m, g, demand!, solve!, unfinished)
end

"B-theta model: the intact copy always, a full copy of an outage network once violated."
function btheta_lazy(n::Net, H, Q, outs)
    m, g, bal = base_model(n)
    at_from = [Int[] for _ in 1:n.N]
    at_to = [Int[] for _ in 1:n.N]
    for l in 1:n.L
        push!(at_from[n.from[l]], l)
        push!(at_to[n.to[l]], l)
    end
    gens = [findall(==(i), n.gen_bus) for i in 1:n.N]
    blocks = Dict{Int,Any}()
    d = Ref(zeros(n.N))
    function add!(c)
        haskey(blocks, c) && return
        th = @variable(m, [1:n.N])
        lines = [l for l in 1:n.L if l != c]
        f = @variable(m, [lines])
        @constraint(m, th[n.ref] == 0)
        for l in lines
            @constraint(m, f[l] == n.b[l] * (th[n.from[l]] - th[n.to[l]]))
            set_lower_bound(f[l], -n.F[l])
            set_upper_bound(f[l], n.F[l])
        end
        blocks[c] = @constraint(m, [i=1:n.N],
            sum(f[l] for l in at_from[i] if l != c; init=0.0) -
            sum(f[l] for l in at_to[i] if l != c; init=0.0) -
            sum(g[k] for k in gens[i]; init=0.0) == -d[][i])
    end
    add!(0)
    function demand!(dd)
        d[] = dd
        set_normalized_rhs(bal, sum(dd))
        for kcl in values(blocks), i in 1:n.N
            set_normalized_rhs(kcl[i], -dd[i])
        end
    end
    unfinished = Ref(false)
    solve!() = lazy_loop!(m, g,
        x -> unique(c for (_, c) in violations(n, H, Q, outs, x, d[]) if !haskey(blocks, c)), add!, unfinished)
    return (; m, g, demand!, solve!, unfinished)
end

"""
Areas by recursive spectral bisection: split the largest area at the median of
its Fiedler vector until there are k.
"""
function partition(n::Net, k)
    parts = [collect(1:n.N)]
    while length(parts) < k
        i = argmax(length.(parts))
        P = parts[i]
        pos = Dict(b => j for (j, b) in enumerate(P))
        Lp = zeros(length(P), length(P))
        for l in 1:n.L
            (haskey(pos, n.from[l]) && haskey(pos, n.to[l])) || continue
            a, b = pos[n.from[l]], pos[n.to[l]]
            a == b && continue
            Lp[a, a] += 1; Lp[b, b] += 1; Lp[a, b] -= 1; Lp[b, a] -= 1
        end
        v = eigen(Symmetric(Lp)).vectors[:, 2]
        left = v .<= median(v)
        (all(left) || !any(left)) && break
        parts[i] = P[left]
        push!(parts, P[.!left])
    end
    area = zeros(Int, n.N)
    for (a, P) in enumerate(parts)
        area[P] .= a
    end
    return area
end

"""
Decomposed PTDF (Alkhraijah, Sigler, Knueven, Maack, PowerUp 2026). For each
area a, the buses outside it and its boundary are Kron-reduced away:
    B^a = [B_II  B_IB; B_BI  B_BB + A^a B_EB],   A^a = -B_BE B_EE^-1,
boundary buses carry equivalent injections s^a = p_B + A^a p_E (the consistency
constraints), and each line's flow comes from its area's small PTDF over the
area's internal injections and s^a. Post-outage flows use the system LODF.
Every flow limit is added lazily. Exact: the flows equal H p.
"""
function ptdf_decomposed(n::Net, H, Q, outs; areas=4, flows=false)
    area = partition(n, areas)
    nA = maximum(area)
    m, g, bal = base_model(n)
    gens = [findall(==(i), n.gen_bus) for i in 1:n.N]
    B = spzeros(n.N, n.N)
    for l in 1:n.L
        u, v, y = n.from[l], n.to[l], n.b[l]
        B[u, u] += y; B[v, v] += y; B[u, v] -= y; B[v, u] -= y
    end
    flow = Vector{Any}(undef, n.L)           # (expression over g and s, demand weights)
    links = []                               # (constraints, boundary, A, external) per area
    for a in 1:nA
        I = findall(==(a), area)
        inI = falses(n.N); inI[I] .= true
        Bd = sort(unique(vcat([inI[n.from[l]] && !inI[n.to[l]] ? [n.to[l]] :
                               inI[n.to[l]] && !inI[n.from[l]] ? [n.from[l]] : Int[] for l in 1:n.L]...)))
        E = setdiff(1:n.N, I, Bd)
        R = [I; Bd]
        nI, nB = length(I), length(Bd)
        Ba = Matrix(B[R, R])
        if isempty(E)
            A = zeros(nB, 0)
        else
            X = Matrix(B[E, E]) \ Matrix(B[E, Bd])          # |E| x |B|
            A = -transpose(X)                                # -B_BE B_EE^-1 (B symmetric)
            Ba[nI+1:end, nI+1:end] .+= A * Matrix(B[E, Bd])
        end
        A[abs.(A) .< 1e-12] .= 0.0
        s = @variable(m, [1:nB])
        gE = [(t, k) for (t, e) in enumerate(E) for k in gens[e]]
        con = @constraint(m, [j=1:nB],
            s[j] - sum(g[k] for k in gens[Bd[j]]; init=0.0) -
            sum(A[j, t] * g[k] for (t, k) in gE if A[j, t] != 0; init=0.0) == 0.0)
        push!(links, (; con, Bd, A, E, s))
        # the area's PTDF, reference at its first bus
        Binv = zeros(length(R), length(R))
        Binv[2:end, 2:end] = inv(Ba[2:end, 2:end])
        pos = Dict(b => j for (j, b) in enumerate(R))
        for l in 1:n.L
            inI[n.from[l]] || continue                       # a line belongs to its from-bus's area
            h = n.b[l] .* (Binv[pos[n.from[l]], :] .- Binv[pos[n.to[l]], :])
            ex = AffExpr(0.0)
            wd = zeros(n.N)
            for j in 1:nI
                h[j] == 0 && continue
                foreach(k -> add_to_expression!(ex, h[j], g[k]), gens[R[j]])
                wd[R[j]] += h[j]
            end
            for j in 1:nB
                h[nI+j] == 0 || add_to_expression!(ex, h[nI+j], s[j])
            end
            flow[l] = (ex, sparsevec(wd))
        end
    end
    lim = lazy_limits(n, H, Q, outs, m, l -> flow[l]; flows)
    function demand!(dd)
        set_normalized_rhs(bal, sum(dd))
        for x in links
            rhs = -dd[x.Bd] .- (isempty(x.E) ? zeros(length(x.Bd)) : x.A * dd[x.E])
            for j in eachindex(x.Bd)
                set_normalized_rhs(x.con[j], rhs[j])
            end
        end
        lim.update!(dd)
    end
    unfinished = Ref(false)
    solve!() = lazy_loop!(m, g, lim.check, lim.add!, unfinished)
    return (; m, g, demand!, solve!, unfinished, area, flow, links)
end

"Largest gap between the decomposed flows and H p at dispatch gv and demand d (a check of exactness)."
function decomposed_gap(model, n::Net, H, gv, d)
    p = injection(n, gv, d)
    val = Dict{VariableRef,Float64}(model.g[k] => gv[k] for k in eachindex(gv))
    for x in model.links
        sv = p[x.Bd] .+ (isempty(x.E) ? zeros(length(x.Bd)) : x.A * p[x.E])
        foreach(j -> val[x.s[j]] = sv[j], eachindex(x.Bd))
    end
    f = H * p
    return maximum(abs(value(v -> val[v], model.flow[l][1]) - dot(model.flow[l][2], d) - f[l]) for l in 1:n.L)
end

# ---- the reduced networks in lazy forms ---------------------------------------------
# `blocks` lists the reduced networks as (out, mask, ratings). Both forms solve the
# same problem as the all-blocks B-theta model (Evaluate.network_scdcopf); only
# what is in the model at a time differs.

"Each reduced network's PTDF rows for its limited kept lines."
function reduced_ptdfs(n::Net, blocks)
    map(blocks) do (c, mask, rating)
        mp = ptdf(n; internal=mask, out=c)
        lim = [l for l in findall(mp.external) if isfinite(rating[l])]
        (; c, mask, rating, lim, Hl=mp.H[lim, :])
    end
end

"Limits of the reduced networks that dispatch g violates: (network index, line)."
function reduced_violations(R, n::Net, g, d)
    p = injection(n, g, d)
    bad = Tuple{Int,Int}[]
    for (i, r) in enumerate(R)
        f = r.Hl * p
        for (j, l) in enumerate(r.lim)
            abs(f[j]) > r.rating[l] * (1 + TOL) && push!(bad, (i, l))
        end
    end
    return bad
end

"""
The reduced networks in PTDF form: a network's line limit becomes a row over the
generators once violated. With `flows`, the line's flow is a variable with one
defining row and the limit is its bound (the paper's form).
"""
function reduced_ptdf_lazy(n::Net, blocks; flows=false)
    R = reduced_ptdfs(n, blocks)
    m, g, bal = base_model(n)
    K = length(n.gen_bus)
    pool = Dict{Tuple{Int,Int},Any}()
    d = Ref(zeros(n.N))
    function add!((i, l))
        haskey(pool, (i, l)) && return
        r = R[i]
        row = r.Hl[findfirst(==(l), r.lim), :]
        a, w, cap = row[n.gen_bus], dot(row, d[]), r.rating[l]
        if flows
            f = @variable(m, lower_bound=-cap, upper_bound=cap)
            pool[(i, l)] = (; row, cap, def=@constraint(m, f - sum(a[k] * g[k] for k in 1:K) == -w))
        else
            pool[(i, l)] = (; row, cap,
                up=@constraint(m,  sum(a[k] * g[k] for k in 1:K) <= cap + w),
                dn=@constraint(m, -sum(a[k] * g[k] for k in 1:K) <= cap - w))
        end
    end
    function demand!(dd)
        d[] = dd
        set_normalized_rhs(bal, sum(dd))
        for r in values(pool)
            w = dot(r.row, dd)
            if flows
                set_normalized_rhs(r.def, -w)
            else
                set_normalized_rhs(r.up, r.cap + w)
                set_normalized_rhs(r.dn, r.cap - w)
            end
        end
    end
    unfinished = Ref(false)
    solve!() = lazy_loop!(m, g, x -> filter(p -> !haskey(pool, p), reduced_violations(R, n, x, d[])), add!, unfinished)
    return (; m, g, demand!, solve!, unfinished)
end

"""
The reduced networks as B-theta blocks (cluster angles, kept lines, cluster
balance), the intact one from the start and an outage's block once its limits
are violated.
"""
function reduced_btheta_lazy(n::Net, blocks)
    R = reduced_ptdfs(n, blocks)
    m, g, bal = base_model(n)
    K = length(n.gen_bus)
    added = Dict{Int,Any}()
    d = Ref(zeros(n.N))
    function rhs!(b, dd)
        dc = zeros(length(b.kcl))
        foreach(i -> dc[b.cl[i]] += dd[i], 1:n.N)
        foreach(j -> set_normalized_rhs(b.kcl[j], -dc[j]), eachindex(dc))
    end
    function add!(i)
        haskey(added, i) && return
        c, mask, rating = blocks[i]
        cl = clusters(n, mask)
        nc = maximum(cl)
        ext = [l for l in 1:n.L if l != c && cl[n.from[l]] != cl[n.to[l]]]
        th = @variable(m, [1:nc])
        @constraint(m, th[cl[n.ref]] == 0)
        f = @variable(m, [ext])
        for l in ext
            @constraint(m, f[l] == n.b[l] * (th[cl[n.from[l]]] - th[cl[n.to[l]]]))
            isfinite(rating[l]) || continue
            set_lower_bound(f[l], -rating[l])
            set_upper_bound(f[l], rating[l])
        end
        out_of = [Int[] for _ in 1:nc]
        in_to = [Int[] for _ in 1:nc]
        foreach(l -> (push!(out_of[cl[n.from[l]]], l); push!(in_to[cl[n.to[l]]], l)), ext)
        gens = [[k for k in 1:K if cl[n.gen_bus[k]] == j] for j in 1:nc]
        kcl = @constraint(m, [j=1:nc], sum(f[l] for l in out_of[j]; init=0.0) -
                                      sum(f[l] for l in in_to[j]; init=0.0) -
                                      sum(g[k] for k in gens[j]; init=0.0) == 0.0)
        added[i] = (; cl, kcl)
        rhs!(added[i], d[])
    end
    i0 = findfirst(b -> b[1] == 0, blocks)
    isnothing(i0) || add!(i0)
    function demand!(dd)
        d[] = dd
        set_normalized_rhs(bal, sum(dd))
        foreach(b -> rhs!(b, dd), values(added))
    end
    unfinished = Ref(false)
    solve!() = lazy_loop!(m, g,
        x -> unique(i for (i, _) in reduced_violations(R, n, x, d[]) if !haskey(added, i)), add!, unfinished)
    return (; m, g, demand!, solve!, unfinished)
end

end # module

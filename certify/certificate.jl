# Certificate for a merged network: every dispatch the reduced network allows is
# deliverable on the full network, for every demand in a given set.
#
#   true flow on an external line = reduced flow + error,  error = M p
#   M = H - Hm   (full PTDF minus the merged network's PTDF, mapped back to buses)
#
# Over the injection set P = {C g - d : d in hull(D), g in its box, sum g = sum d}:
#   merged line l      safe if max |H_l p| <= F_l            (its limit is gone)
#   external line l    no change needed if max |H_l p| <= F_l, else its reduced
#                      rating becomes F_l - max |M_l p|        (triangle inequality)
# Every max is a small LP, so no OPF and no scenario loop.

module Certificate

using LinearAlgebra, JuMP, Gurobi

export merged_ptdf, ptdf, injection_set, extremes, certify, max_loading, minimal_derating,
       design_ratings, reach

"Cluster index of every bus, from the merged lines."
function clusters(base, internal)
    parent = collect(1:base.N)
    find(x) = parent[x] == x ? x : (parent[x] = find(parent[x]))
    for l in 1:base.Ln
        internal[l] || continue
        a, b = find(base.Efrom[l]), find(base.Eto[l])
        a == b || (parent[a] = b)
    end
    rep = [find(b) for b in 1:base.N]
    ids = Dict(r => i for (i, r) in enumerate(sort(unique(rep))))
    return [ids[rep[b]] for b in 1:base.N]
end

"""
Flow sensitivities of the merged network, one row per line and one column per
bus: reduced flow on line l = Hm[l, :] * p for any balanced injection p. Rows of
merged lines are zero. With nothing merged this is the full PTDF.
"""
function merged_ptdf(base, internal)
    cl = clusters(base, internal)
    K = maximum(cl)
    u, v, Ln = base.Efrom, base.Eto, base.Ln
    ext = BitVector([cl[u[l]] != cl[v[l]] for l in 1:Ln])
    ref = cl[base.j0]
    pos = zeros(Int, K)
    free = setdiff(1:K, [ref])
    for (i, k) in enumerate(free)
        pos[k] = i
    end
    B = zeros(length(free), length(free))
    W = zeros(Ln, length(free))
    for l in 1:Ln
        ext[l] || continue
        a, b, y = pos[cl[u[l]]], pos[cl[v[l]]], base.Dx[l]
        a > 0 && (B[a, a] += y; W[l, a] += y)
        b > 0 && (B[b, b] += y; W[l, b] -= y)
        a > 0 && b > 0 && (B[a, b] -= y; B[b, a] -= y)
    end
    Pc = zeros(Ln, K)
    # LU, not Cholesky: case300 has a negative-reactance branch, so B need not be definite
    isempty(free) || (Pc[:, free] = transpose(lu(B) \ transpose(W)))
    return (Hm=Pc[:, cl], external=ext, cluster=cl, K=K)
end

ptdf(base) = merged_ptdf(base, falses(base.Ln)).Hm

"""
The plausible injections: demand anywhere in the convex hull of the columns of
`D` (N x T), every generator anywhere in [pmin, pmax] with power balance. With
`relax_pmin` the lower bound is min(0, pmin), as in the OPFs of this repo.
Queries are LPs on one model, re-solved with a new objective.
"""
function injection_set(base, D::AbstractMatrix; relax_pmin::Bool=true)
    K, T = length(base.gen_bus), size(D, 2)
    pmin = relax_pmin ? min.(0.0, base.pmin) : base.pmin
    env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
    m = Model(() -> Gurobi.Optimizer(env))
    set_silent(m)
    set_optimizer_attribute(m, "FeasibilityTol", 1e-9)
    set_optimizer_attribute(m, "OptimalityTol", 1e-9)
    set_optimizer_attribute(m, "Threads", 1)
    @variable(m, pmin[k] <= g[k=1:K] <= base.pmax[k])
    @variable(m, lam[1:T] >= 0)
    @constraint(m, sum(lam) == 1)
    total = vec(sum(D; dims=1))
    @constraint(m, sum(g) == sum(total[t] * lam[t] for t in 1:T))
    return (; m, g, lam, D, base, rows=ConstraintRef[])
end

"Worst loading max |H_l p| / F_l of every line over the set, on the full network."
reach(base, S; H=ptdf(base)) =
    [begin lo, hi = extremes(S, view(H, l, :)); max(hi, -lo) / base.frate[l] end
     for l in 1:base.Ln]

"min and max of row' * p over the set."
function extremes(S, row::AbstractVector)
    a = [row[b] for b in S.base.gen_bus]
    w = transpose(S.D) * row
    obj = @expression(S.m, sum(a[k] * S.g[k] for k in eachindex(a)) -
                           sum(w[t] * S.lam[t] for t in eachindex(w)))
    out = Float64[]
    for sense in (MOI.MIN_SENSE, MOI.MAX_SENSE)
        @objective(S.m, sense, obj)
        optimize!(S.m)
        st = termination_status(S.m)
        st in (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED) && return NaN, NaN   # empty set
        st == MOI.OPTIMAL || error("injection-set LP ended $st")
        push!(out, objective_value(S.m))
    end
    return out[1], out[2]
end

"""
Restrict the set to what a reduced network allows: |Hm_l p| <= ratings[l] on its
external lines. Returns (line, upper row, lower row) per external line.
"""
function restrict!(S, Hm, external, ratings)
    rowmap = Tuple{Int,ConstraintRef,ConstraintRef}[]
    for l in findall(external)
        row = Hm[l, :]
        a = [row[b] for b in S.base.gen_bus]
        w = transpose(S.D) * row
        ex = @expression(S.m, sum(a[k] * S.g[k] for k in eachindex(a)) -
                              sum(w[t] * S.lam[t] for t in eachindex(w)))
        up, dn = @constraint(S.m, ex <= ratings[l]), @constraint(S.m, -ex <= ratings[l])
        push!(S.rows, up, dn)
        push!(rowmap, (l, up, dn))
    end
    return rowmap
end

"""
    design_ratings(base, S, internal; H=ptdf(base), lines, tol=1e-6)

Reduced ratings that pass the exact check on `lines` (the critical ones). First
scale down along the error bounds e until the check passes (F - alpha e), then
raise each derated line back as far as the check allows, largest derating first
(lines derated by less than `keep` of rating stay put). Every step is the exact
LP check, so what comes back is certified; certified = false when even alpha = 1
fails (a merged line that can still overload) or no dispatch is left.

The reduced rows stay in the model and a new set of ratings only changes their
right-hand sides, so each check re-solves from the previous basis.
"""
function design_ratings(base, S, internal; H=ptdf(base), lines, tol::Real=1e-6,
                        iters::Int=8, keep::Real=0.005)
    F = base.frate
    mp = merged_ptdf(base, internal)
    checks = Ref(0)
    done(r, a, reason="") = (rating=r, certified=isempty(reason), alpha=a,
                             checks=checks[], reason=reason)
    isempty(lines) && return done(copy(F), 0.0)
    cert = certify(base, internal, S; H)            # error bounds, on the unrestricted set
    e = F .- cert.rating
    rowmap = restrict!(S, mp.Hm, mp.external, F)
    rows = [(l, sgn, sgn .* H[l, :]) for l in lines for sgn in (1.0, -1.0)]
    function ok(r)
        checks[] += 1
        for (k, up, dn) in rowmap
            set_normalized_rhs(up, r[k]); set_normalized_rhs(dn, r[k])
        end
        for (l, _, row) in rows                     # stop at the first violation
            a = [row[b] for b in S.base.gen_bus]
            w = transpose(S.D) * row
            @objective(S.m, Max, sum(a[k] * S.g[k] for k in eachindex(a)) -
                                 sum(w[t] * S.lam[t] for t in eachindex(w)))
            optimize!(S.m)
            termination_status(S.m) == MOI.OPTIMAL || return false
            objective_value(S.m) <= F[l] * (1 + tol) || return false
        end
        return true
    end
    result = try
        if ok(F)
            done(copy(F), 0.0)
        elseif !ok(F .- e)
            done(F .- e, NaN, any(v -> v <= 0, F .- e) ?
                 "error bound exceeds a rating" : "a merged line can still overload")
        else
            lo, hi = 0.0, 1.0
            for _ in 1:iters
                mid = (lo + hi) / 2
                ok(F .- mid .* e) ? (hi = mid) : (lo = mid)
            end
            r = F .- hi .* e
            for k in sortperm(e; rev=true)
                F[k] - r[k] > keep * F[k] || continue
                trial = copy(r); trial[k] = F[k]
                if ok(trial)
                    r = trial
                    continue
                end
                a, b = r[k], F[k]
                for _ in 1:iters
                    m = (a + b) / 2
                    trial[k] = m
                    ok(trial) ? (a = m) : (b = m)
                end
                r[k] = a
            end
            done(r, hi)
        end
    finally
        unrestrict!(S)
    end
    return result
end

function unrestrict!(S)
    foreach(c -> delete(S.m, c), S.rows)
    empty!(S.rows)
    return S
end

"""
    certify(base, internal, S; H=ptdf(base), tol=1e-7)

Per line: kind, reach = max |H_l p| / F_l over S, error bound e_l (critical
external lines only) and the reduced rating. `certified` is false when a merged
line can overload or a derated rating would be <= 0.
"""
function certify(base, internal, S; H=ptdf(base), tol::Real=1e-7)
    mp = merged_ptdf(base, internal)
    F = base.frate
    rating = copy(F)
    rows = NamedTuple[]
    for l in 1:base.Ln
        lo, hi = extremes(S, view(H, l, :))
        reach = max(hi, -lo) / F[l]
        kind, e = mp.external[l] ? "external" : "merged", 0.0
        if mp.external[l] && reach > 1 + tol
            elo, ehi = extremes(S, H[l, :] .- mp.Hm[l, :])
            e = max(ehi, -elo, 0.0)
            rating[l] = F[l] - e
        end
        push!(rows, (line=l, kind=kind, reach=reach, critical=reach > 1 + tol,
                     error_bound=e, rating=F[l], reduced_rating=rating[l],
                     derating_pct=100 * (F[l] - rating[l]) / F[l]))
    end
    unsafe = [r.line for r in rows if r.kind == "merged" && r.critical]
    nonpos = [r.line for r in rows if r.kind == "external" && r.reduced_rating <= 0]
    return (; rows, rating, unsafe_merged=unsafe, nonpositive=nonpos,
            certified=isempty(unsafe) && isempty(nonpos), mp)
end

"""
Exact check: the largest true loading max |H_l p| / F_l over the set cut down
by the reduced network's own rows with `ratings`, for each line in `lines`. A
certified network keeps every line at or below 1. NaN when no dispatch is left.
"""
function max_loading(base, S, mp, ratings; H=ptdf(base), lines=1:base.Ln)
    restrict!(S, mp.Hm, mp.external, ratings)
    load = try
        [begin lo, hi = extremes(S, view(H, l, :)); max(hi, -lo) / base.frate[l] end
         for l in lines]
    finally
        unrestrict!(S)
    end
    return load
end

"""
Smallest derating that passes the exact check, along F - alpha * e (alpha in
[0, 1], e the certificate's error bounds). Lines whose box reach is <= 1 are
safe whatever the ratings, so only the critical ones are checked. Returns
alpha = NaN when even alpha = 1 fails, i.e. a merged line can still overload.
"""
function minimal_derating(base, S, cert; H=ptdf(base), tol::Real=1e-6, iters::Int=20)
    crit = [r.line for r in cert.rows if r.critical]
    F, e = base.frate, base.frate .- cert.rating
    ok(a) = all(x -> isfinite(x) && x <= 1 + tol,
                max_loading(base, S, cert.mp, F .- a .* e; H, lines=crit))
    (isempty(crit) || ok(0.0)) && return (alpha=0.0, rating=copy(F), lines=crit)
    ok(1.0) || return (alpha=NaN, rating=copy(cert.rating), lines=crit)
    lo, hi = 0.0, 1.0
    for _ in 1:iters
        mid = (lo + hi) / 2
        ok(mid) ? (hi = mid) : (lo = mid)
    end
    return (alpha=hi, rating=F .- hi .* e, lines=crit)
end

end # module

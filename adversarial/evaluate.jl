# The three SC-DCOPF forms the reduced networks are judged against:
#   network  - every block as a network copy (cluster angles, kept lines, ratings);
#              with no merges and every outage listed this is the full SC-DCOPF
#   compact  - only the screened (block, line) pairs, as PTDF rows over the generators
# Each model is built once; an hour only moves right-hand sides.
# Needs adversarial/adversarial.jl loaded first (module Adversarial).

module Evaluate

using LinearAlgebra, JuMP, Gurobi
using Main.Adversarial: Net, clusters, lp

export network_scdcopf, compact_scdcopf, solve_hours, model_size

"""
SC-DCOPF from network blocks. `blocks` lists (out, mask) or (out, mask, ratings):
the network without line `out` (0 = intact), with the lines in `mask` merged
and its kept lines rated `ratings` (default the network's own).
"""
function network_scdcopf(n::Net, blocks)
    m = lp()
    K = length(n.gen_bus)
    @variable(m, n.pmin[k] <= g[k=1:K] <= n.pmax[k])
    bal = @constraint(m, sum(g) == 0.0)
    parts = []
    for blk in blocks
        out, mask = blk[1], blk[2]
        F = length(blk) > 2 ? blk[3] : n.F
        cl = clusters(n, mask)
        nc = maximum(cl)
        ext = [l for l in 1:n.L if l != out && cl[n.from[l]] != cl[n.to[l]]]
        th = @variable(m, [1:nc])
        @constraint(m, th[cl[n.ref]] == 0)
        f = @variable(m, [ext])
        for l in ext
            @constraint(m, f[l] == n.b[l] * (th[cl[n.from[l]]] - th[cl[n.to[l]]]))
            isfinite(F[l]) || continue      # no limit on this line
            set_lower_bound(f[l], -F[l])
            set_upper_bound(f[l], F[l])
        end
        gens = [[k for k in 1:K if cl[n.gen_bus[k]] == i] for i in 1:nc]
        kcl = @constraint(m, [i=1:nc],
            sum(f[l] for l in ext if cl[n.from[l]] == i; init=0.0) -
            sum(f[l] for l in ext if cl[n.to[l]] == i; init=0.0) -
            sum(g[k] for k in gens[i]; init=0.0) == 0.0)
        push!(parts, (; cl, kcl))
    end
    @objective(m, Min, sum(n.c1[k] * g[k] for k in 1:K) + n.c0)
    function demand!(d)
        set_normalized_rhs(bal, sum(d))
        for P in parts
            dc = zeros(length(P.kcl))
            for i in eachindex(d)
                dc[P.cl[i]] += d[i]
            end
            for i in eachindex(dc)
                set_normalized_rhs(P.kcl[i], -dc[i])
            end
        end
    end
    return (; m, g, demand!)
end

"SC-DCOPF with one PTDF row per listed (row over buses, line) pair."
function compact_scdcopf(n::Net, rows)
    m = lp()
    K = length(n.gen_bus)
    @variable(m, n.pmin[k] <= g[k=1:K] <= n.pmax[k])
    bal = @constraint(m, sum(g) == 0.0)
    cons = map(rows) do (row, l)
        a = row[n.gen_bus]
        (; row, l,
         up=@constraint(m,  sum(a[k] * g[k] for k in 1:K) <= 0.0),
         dn=@constraint(m, -sum(a[k] * g[k] for k in 1:K) <= 0.0))
    end
    @objective(m, Min, sum(n.c1[k] * g[k] for k in 1:K) + n.c0)
    function demand!(d)
        set_normalized_rhs(bal, sum(d))
        for c in cons
            w = dot(c.row, d)
            set_normalized_rhs(c.up, n.F[c.l] + w)
            set_normalized_rhs(c.dn, n.F[c.l] - w)
        end
    end
    return (; m, g, demand!)
end

"""
Solve at every column of D: cost (NaN if infeasible), dispatch, total solver
seconds, total wall seconds and the hours attempted (all of them unless the
`time_limit` in seconds ran out first). A model with a `solve!` (a lazy method)
solves through it, which returns its solver seconds.
"""
function solve_hours(model, D; time_limit=Inf)
    T = size(D, 2)
    cost, disp, secs = fill(NaN, T), zeros(length(model.g), T), 0.0
    t0 = time()
    for s in 1:T
        time() - t0 > time_limit && return cost, disp, secs, time() - t0, s - 1
        model.demand!(D[:, s])
        if hasproperty(model, :solve!)
            secs += model.solve!()
            # a lazy solve cut short by its deadline has no valid answer for this hour
            hasproperty(model, :unfinished) && model.unfinished[] && return cost, disp, secs, time() - t0, s - 1
        else
            optimize!(model.m)
            secs += solve_time(model.m)
        end
        termination_status(model.m) == MOI.OPTIMAL || continue
        cost[s], disp[:, s] = objective_value(model.m), value.(model.g)
    end
    return cost, disp, secs, time() - t0, T
end

"Variables, constraints and nonzeros as Gurobi sees them (after a solve)."
model_size(m) = (vars=get_attribute(m, Gurobi.ModelAttribute("NumVars")),
                 cons=get_attribute(m, Gurobi.ModelAttribute("NumConstrs")),
                 nnz=get_attribute(m, Gurobi.ModelAttribute("NumNZs")))

end # module

# Final check and output files for a reduction, shared by kkt/, greedy/ and proxy/.
#
# Design demands: every one must have a deliverable cheapest reduced dispatch
# within the cost cap, and we also ask whether EVERY cheapest dispatch is safe.
# Held-out demands: the same check on demands the method never saw.
#
# summary.csv has the same leading columns for every approach, so runs can be
# stacked into one table.

module Results

using DelimitedFiles, Printf
import Main.Caps

export write_rows, finish, run_info, start_objective

function write_rows(path, rows)
    isempty(rows) && return
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(names, ","))
        for row in rows
            println(io, join((getproperty(row, n) for n in names), ","))
        end
    end
end

"""
Objective of the MIP start Gurobi accepted, read from its log. NaN when it
rejected the start or none was given. The last solve in the file counts.
"""
function start_objective(log_file)
    (isnothing(log_file) || !isfile(log_file)) && return NaN
    val = NaN
    for line in eachline(log_file)
        occursin("Gurobi Optimizer version", line) && (val = NaN)
        m = match(r"(?:Loaded user MIP start|User MIP start produced solution) with objective (\S+)", line)
        if !isnothing(m)
            val = something(tryparse(Float64, m.captures[1]), NaN)
        elseif occursin("User MIP start did not produce a new incumbent", line)
            val = NaN
        end
    end
    return val
end

"The run columns every summary.csv carries, in a fixed order."
function run_info(; mode, setting, hop_caps="none", status="", merged=-1, bound=NaN,
                  seed_merged=-1, seed_buses=-1, start_obj=NaN,
                  solve_seconds=NaN, seconds=NaN)
    gap = isfinite(bound) && merged > 0 ? round(100 * (bound - merged) / merged; digits=3) : NaN
    return (mode=mode, setting=setting, hop_caps=hop_caps, status=string(status),
            merged_bound=bound, gap_pct=gap, seed_merged=seed_merged, seed_buses=seed_buses,
            start_objective=start_obj,
            seed_used=seed_merged < 0 ? "" : string(start_obj >= seed_merged - 0.5),
            solve_seconds=round(solve_seconds; digits=1), seconds=round(seconds; digits=1))
end

"""
    finish(Audit, BL, design, heldout, internal, out; approach, case, cost_cap,
           info=run_info(...), extra=NamedTuple(), audit_seconds=600.0) -> summary row

`cost_cap = nothing` judges deliverability only (the proxy has no cost cap).
"""
function finish(Audit, BL, design, heldout, internal, out; approach, case, cost_cap,
                info, extra=NamedTuple(), audit_seconds::Real=600.0)
    mkpath(out)
    base = design.base
    N, Ln = base.N, base.Ln
    clus = BL.clustering_from_internal(base, internal)
    internal = BitVector([clus.rep_of[base.Efrom[l]] == clus.rep_of[base.Eto[l]] for l in 1:Ln])
    nb = length(clus.retained)

    # The audit runs without the cap, so a demand over the cap still reports its
    # overload; the cap is applied to the cost afterwards.
    function audit(c, all_optima)
        idx = collect(axes(c.load, 2))
        costs = Dict(s => BL.true_optimal_cost(base, c.load[:, s]; relax_pmin=true) for s in idx)
        res = Audit.audit_topology(c, internal, idx; all_optima, deadline=time() + audit_seconds)
        return [begin
                    full = costs[x.scenario]
                    pct = 100 * (x.reduced_cost - full) / max(abs(full), 1e-12)
                    within = isnothing(cost_cap) || x.reduced_cost <=
                        full + cost_cap / 100 * abs(full) + 1e-8 * max(1.0, abs(x.reduced_cost))
                    merge(Base.structdiff(x, NamedTuple{(:cost_gap_pass,)}),
                          (full_cost=full, cost_pct=pct,
                           deliverable=x.classification == :pass,
                           pass=x.classification == :pass && within))
                end for x in res.rows]
    end
    worst(rows, f) = isempty(rows) ? NaN : maximum(f, rows)
    t0 = time()
    a = audit(design, true)
    write_rows(joinpath(out, "design_audit.csv"), a)
    h = isnothing(heldout) ? NamedTuple[] : audit(heldout, false)
    isempty(h) || write_rows(joinpath(out, "heldout_audit.csv"), h)

    writedlm(joinpath(out, "internal.csv"), Int.(internal), ',')
    writedlm(joinpath(out, "assignment.csv"), clus.A, ',')
    open(joinpath(out, "bus_mapping.csv"), "w") do io
        println(io, "bus_id,representative_bus_id,retained")
        for b in 1:N
            println(io, design.bus_ids[b], ",", design.bus_ids[clus.rep_of[b]], ",",
                    clus.rep_of[b] == b ? 1 : 0)
        end
    end

    sizes = Caps.cluster_sizes(base, internal)
    row = (approach=approach, case=case, info...,
           buses=N, remaining_buses=nb, reduction_pct=round(100 * (N - nb) / N; digits=2),
           merged_lines=count(internal), lines=Ln,
           largest_cluster=maximum(values(sizes)),
           hop_diameter=Caps.hop_diameter(base, internal),
           cost_cap=something(cost_cap, NaN),
           design_count=length(a), design_pass=count(r -> r.pass, a),
           design_deliverable=count(r -> r.deliverable, a),
           design_all_optima_safe=count(r -> r.all_optima_checked && r.all_optima_safe, a),
           worst_design_overload_pct=round(100 * worst(a, r -> r.excess); digits=5),
           worst_design_cost_pct=round(worst(a, r -> r.cost_pct); digits=5),
           heldout_count=length(h), heldout_pass=count(r -> r.pass, h),
           heldout_deliverable=count(r -> r.deliverable, h),
           worst_heldout_overload_pct=round(100 * worst(h, r -> r.excess); digits=5),
           worst_heldout_cost_pct=round(worst(h, r -> r.cost_pct); digits=5),
           audit_seconds=round(time() - t0; digits=1), extra...)
    write_rows(joinpath(out, "summary.csv"), [row])

    println()
    println("=" ^ 72)
    @printf("%s  %s:  %d -> %d buses (%.1f%%), %d merged lines\n",
            approach, case, N, nb, row.reduction_pct, row.merged_lines)
    @printf("largest cluster %d buses, longest merged chain (BFS) %d\n",
            row.largest_cluster, row.hop_diameter)
    @printf("design demands   %d/%d pass   deliverable %d   every optimum safe on %d   worst overload %.4f%%   worst cost %+.4f%%\n",
            row.design_pass, row.design_count, row.design_deliverable,
            row.design_all_optima_safe, row.worst_design_overload_pct, row.worst_design_cost_pct)
    isempty(h) || @printf("held-out demands %d/%d pass   deliverable %d   worst overload %.4f%%   worst cost %+.4f%%\n",
                          row.heldout_pass, row.heldout_count, row.heldout_deliverable,
                          row.worst_heldout_overload_pct, row.worst_heldout_cost_pct)
    println("results -> ", out)
    println("=" ^ 72)
    return row
end

end # module

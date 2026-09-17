# Final check and output files for a reduction, shared by kkt/ and greedy/.
#
# Design demands: every one must have a deliverable cheapest reduced dispatch
# within the cost cap, and we also ask whether EVERY cheapest dispatch is safe.
# Held-out demands: the same check on demands the method never saw.

module Results

using DelimitedFiles, Printf
import Main.Caps

export write_rows, finish

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
    finish(Audit, BL, design, heldout, internal, out; approach, case, cost_cap,
           extra=NamedTuple(), audit_seconds=600.0) -> summary row
"""
function finish(Audit, BL, design, heldout, internal, out; approach, case, cost_cap,
                extra=NamedTuple(), audit_seconds::Real=600.0)
    mkpath(out)
    base = design.base
    N, Ln = base.N, base.Ln
    clus = BL.clustering_from_internal(base, internal)
    internal = BitVector([clus.rep_of[base.Efrom[l]] == clus.rep_of[base.Eto[l]] for l in 1:Ln])
    nb = length(clus.retained)

    function audit(c, all_optima)
        local idx = collect(axes(c.load, 2))
        local costs = Dict(s => BL.true_optimal_cost(base, c.load[:, s]; relax_pmin=true) for s in idx)
        local res = Audit.audit_topology(c, internal, idx; all_optima, full_costs=costs,
                                         cost_gap_pct=cost_cap, deadline=time() + audit_seconds)
        local worst = maximum(100 * (x.reduced_cost - costs[x.scenario]) /
                              max(abs(costs[x.scenario]), 1e-12) for x in res.rows)
        return res, worst
    end
    t0 = time()
    a, worst_design = audit(design, true)
    write_rows(joinpath(out, "design_audit.csv"), a.rows)
    ha, worst_held = nothing, NaN
    if !isnothing(heldout)
        ha, worst_held = audit(heldout, false)
        write_rows(joinpath(out, "heldout_audit.csv"), ha.rows)
    end

    writedlm(joinpath(out, "internal.csv"), Int.(internal), ',')
    writedlm(joinpath(out, "assignment.csv"), clus.A, ',')
    open(joinpath(out, "bus_mapping.csv"), "w") do io
        println(io, "bus_id,representative_bus_id,retained")
        for b in 1:N
            println(io, design.bus_ids[b], ",", design.bus_ids[clus.rep_of[b]], ",",
                    clus.rep_of[b] == b ? 1 : 0)
        end
    end

    pass(x) = isnothing(x) ? 0 : count(r -> r.classification == :pass, x.rows)
    sizes = Caps.cluster_sizes(base, internal)
    row = (approach=approach, case=case, buses=N, remaining_buses=nb,
           reduction_pct=round(100 * (N - nb) / N; digits=2),
           merged_lines=count(internal), lines=Ln,
           largest_cluster=maximum(values(sizes)),
           hop_diameter=Caps.hop_diameter(base, internal),
           design_pass=pass(a), design_count=size(design.load, 2),
           design_all_optima_safe=count(r -> r.all_optima_checked && r.all_optima_safe, a.rows),
           heldout_pass=pass(ha), heldout_count=isnothing(heldout) ? 0 : size(heldout.load, 2),
           worst_design_cost_pct=round(worst_design; digits=5),
           worst_heldout_cost_pct=isnan(worst_held) ? NaN : round(worst_held; digits=5),
           audit_seconds=round(time() - t0; digits=1), extra...)
    write_rows(joinpath(out, "summary.csv"), [row])

    println()
    println("=" ^ 72)
    @printf("%s  %s:  %d -> %d buses (%.1f%%), %d merged lines\n",
            approach, case, N, nb, row.reduction_pct, row.merged_lines)
    @printf("largest cluster %d buses, longest merged chain (BFS) %d\n",
            row.largest_cluster, row.hop_diameter)
    @printf("design demands   %d/%d pass   every optimum safe on %d/%d   worst cost %+.4f%%\n",
            row.design_pass, row.design_count, row.design_all_optima_safe,
            row.design_count, worst_design)
    isnothing(heldout) || @printf("held-out demands %d/%d pass   worst cost %+.4f%%\n",
                                  row.heldout_pass, row.heldout_count, worst_held)
    println("results -> ", out)
    println("=" ^ 72)
    return row
end

end # module

# Greedy network reduction.
#
#   1. Radial buses first: merging a degree-1 bus moves no other flow, so these
#      are taken without an LP when the absorbed line can never bind.
#   2. Remaining lines ranked by predicted flow change (flow_sensitivity.jl),
#      re-ranked after every accepted merge.
#   3. A merge is kept only if, at every design demand, some cheapest reduced
#      dispatch is deliverable on the full network within the cost cap.
#   4. Optional final KKT check on all demands jointly, rolling back merges
#      until it passes.

module Greedy

using LinearAlgebra, SparseArrays, Printf
using Main.FlowRanking
import Main.Caps

const FR = Main.FlowRanking

export greedy_reduction

"""
    greedy_reduction(Audit, BL, c, scope; kwargs...)

`Audit` is common/audit.jl, `BL` is kkt/kkt_model.jl.

Keyword options:
  time_limit      seconds; merging stops at 80%, the rest is for the KKT check
  cost_gap_pct    cost cap in percent over the full-network optimum (nothing = off)
  flow_tolerance  allowed overload, fraction of rating, in the acceptance test
  cost_tolerance  cost slack, relative, when looking for a deliverable dispatch
  line_budget     max merged lines (nothing = no limit)
  hop_cap         max merged chain length inside a cluster (nothing = no limit)
  size_cap        max buses per cluster (nothing = no limit)
  ordering        :flow (ranked, re-ranked) | :loading (static, least loaded first)
  norm            :headroom | :rating, what flow changes are divided by
  alpha           weight of the 2-norm against the max-norm in the score
  radial_first    take radial merges without an LP
  kkt_check       final joint KKT check with rollback
"""
function greedy_reduction(Audit, BL, c, scope; time_limit::Real=600.0,
                          cost_gap_pct=0.1, flow_tolerance::Real=1e-9,
                          cost_tolerance::Real=0.0, line_budget=nothing,
                          hop_cap=nothing, size_cap=nothing,
                          ordering::Symbol=:flow, norm::Symbol=:headroom,
                          alpha::Float64=0.5, radial_first::Bool=true,
                          kkt_check::Bool=true, relax_pmin::Bool=true,
                          threads::Int=1, verbose::Bool=true)
    ordering in (:flow, :loading) || error("ordering must be :flow or :loading")
    norm in (:headroom, :rating) || error("norm must be :headroom or :rating")
    base = c.base
    L, N = base.Ln, base.N
    pool = Int.(collect(scope))
    started = time()
    deadline = started + time_limit
    merge_deadline = started + (kkt_check ? 0.8 : 1.0) * time_limit
    buses(mask) = length(BL.clustering_from_internal(base, mask).retained)

    internal = falses(L)
    initial = Audit.audit_topology(c, internal, pool; relax_pmin, threads,
                                   cost_tolerance, flow_tolerance, deadline)
    fail = (covered=false, internal=internal, trace=NamedTuple[],
            elapsed_seconds=time() - started, reason=:full_network_not_verified,
            kkt_verified=false, rollbacks=0, radial_taken=0, rejected=0,
            lp_checks=0, buses=N)
    initial.covered || return fail
    reference = Dict(row.scenario => row.reduced_cost for row in initial.rows)

    check_order = sort(pool; by = s -> -sum(c.load[:, s]))
    trace = NamedTuple[]
    history = [copy(internal)]
    lp_checks = Ref(0)

    function accepts(proposal)
        oracle = Audit.topology_oracle(base, proposal; relax_pmin, threads,
                                       cost_tolerance, flow_tolerance)
        for s in check_order
            a = Audit.check_scenario!(oracle, c.load[:, s]; deadline=merge_deadline,
                                      full_cost=reference[s], cost_gap_pct)
            lp_checks[] += 1
            a.classification == :pass || return (ok=false, why=a.classification)
        end
        return (ok=true, why=:none)
    end

    # Close a line and every line its clustering makes internal.
    function close_line(cur, l)
        prop = copy(cur)
        prop[l] = true
        clus = BL.clustering_from_internal(base, prop)
        return BitVector([clus.rep_of[base.Efrom[k]] == clus.rep_of[base.Eto[k]] for k in 1:L])
    end
    allowed(prop, l) = Caps.within_caps(base, prop; budget=line_budget, hop_cap,
                                        size_cap, bus=base.Efrom[l])

    radial_taken = 0
    if radial_first
        dmax = vec(maximum(c.load; dims=2))
        dmin = vec(minimum(c.load; dims=2))
        for r in FR.radial_candidates(base, dmax, dmin)
            time() < merge_deadline || break
            internal[r.line] && continue
            r.capacity_safe || continue
            prop = close_line(internal, r.line)
            allowed(prop, r.line) || continue
            internal = prop
            push!(history, copy(internal))
            radial_taken += 1
            push!(trace, (line=r.line, accepted=true, kind=:radial, rejection=:none,
                          internal_lines=count(internal), buses=buses(internal),
                          elapsed_seconds=time() - started))
        end
        verbose && @printf(" radial merges: %d (no LP), %d buses left\n",
                           radial_taken, buses(internal))
    end

    fs = FR.build_sensitivity(base)
    for l in 1:L
        internal[l] && FR.apply_merge!(fs, base.Efrom[l], base.Eto[l])
    end
    # Injection held fixed while merging. The audit rows carry no dispatch, so
    # this is the base operating point for every demand, as in the runs so far.
    injections = [base.p for _ in pool]
    thetas = [fs.X * p for p in injections]
    flows = [fs.H * p for p in injections]
    loading_order = sortperm([maximum(abs, c.fhat[l, pool]) / base.frate[l] for l in 1:L])
    rejected = 0

    while time() < merge_deadline
        cands = [l for l in 1:L if !internal[l] && base.Efrom[l] != base.Eto[l]]
        isempty(cands) && break
        order = if ordering === :flow
            sc = FR.merge_scores(fs, thetas; alpha, skip = l -> internal[l],
                                 norm_mode=norm, flows = norm === :headroom ? flows : nothing)
            sort(cands; by = l -> sc[l])
        else
            [l for l in loading_order if !internal[l] && base.Efrom[l] != base.Eto[l]]
        end

        progressed = false
        for l in order
            time() < merge_deadline || break
            internal[l] && continue
            prop = close_line(internal, l)
            allowed(prop, l) || continue
            res = accepts(prop)
            if res.ok
                i, j = base.Efrom[l], base.Eto[l]
                zXz = fs.X[i, i] + fs.X[j, j] - 2 * fs.X[i, j]
                if zXz > 1e-12
                    Xz = fs.X[:, i] .- fs.X[:, j]
                    for (si, th) in enumerate(thetas)
                        thetas[si] = th .+ (-(th[i] - th[j]) / zXz) .* Xz
                    end
                end
                FR.apply_merge!(fs, i, j)
                for si in eachindex(flows)
                    flows[si] = fs.H * injections[si]
                end
                internal = prop
                push!(history, copy(internal))
                push!(trace, (line=l, accepted=true, kind=:ranked, rejection=:none,
                              internal_lines=count(internal), buses=buses(internal),
                              elapsed_seconds=time() - started))
                verbose && @printf(" merge L%-4d -> %d buses  (%.0f s)\n", l, buses(internal),
                                   time() - started)
                progressed = true
                break
            else
                rejected += 1
                push!(trace, (line=l, accepted=false, kind=:ranked, rejection=res.why,
                              internal_lines=count(internal), buses=buses(internal),
                              elapsed_seconds=time() - started))
            end
        end
        progressed || break               # a full pass with nothing accepted
    end

    rollbacks = 0
    kkt_verified = false
    if kkt_check
        pmin = relax_pmin ? min.(0.0, base.pmin) : base.pmin
        while time() < deadline && !isempty(history)
            candidate = last(history)
            point = BL.fixed_topology_kkt_start(base, [c.load[:, s] for s in pool], pmin,
                        candidate; time_limit=max(0.0, deadline - time()), solver_threads=threads)
            cost_ok = !isnothing(point) && (isnothing(cost_gap_pct) || all(
                sum(base.c1 .* point.g[:, j]) + sum(base.c0) <=
                reference[s] + cost_gap_pct / 100 * abs(reference[s]) + 1e-6
                for (j, s) in enumerate(pool)))
            if cost_ok
                internal = candidate
                kkt_verified = true
                break
            end
            pop!(history)
            rollbacks += 1
        end
        kkt_verified || (internal = falses(L))
    end

    return (; covered=true, internal, trace, elapsed_seconds=time() - started,
            reason = !kkt_check ? :lp_checked : kkt_verified ? :kkt_verified : :rolled_back_to_full,
            kkt_verified, rollbacks, radial_taken, rejected,
            lp_checks=lp_checks[], buses=buses(internal))
end

end # module

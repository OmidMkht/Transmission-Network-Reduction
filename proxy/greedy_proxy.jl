# Greedy reduction under the PROXY's feasibility test.
#
# greedy/greedy.jl accepts a merge when SOME cheapest re-dispatch is deliverable
# -- common/audit.jl leaves generation free and minimises cost. The proxy never
# re-dispatches: it holds injections at the reference point and asks only that the
# reduced network can still route them inside each line's eps window. Neither test
# implies the other (greedy may be leaning on a dispatch the proxy will not use;
# the proxy is looser on ratings by eps), so greedy/greedy.jl's answer is not a
# valid proxy warm start. This is the matching greedy.
#
# The oracle is the proxy's own rows with cl FIXED to the candidate mask and a
# slack on every window and balance row. Minimum slack 0 means the clustering is
# proxy-feasible; the value is how far off it is, in p.u. Building the rows from
# the proxy's own definitions is the point -- a reimplemented feasibility test
# that disagreed even slightly would defeat the purpose, since the seed has to be
# feasible for the solve it warm-starts.
#
# The LP is built once. Flipping a line updates two row right-hand sides and two
# variable bounds per scenario, so a sweep costs one LP re-solve per candidate
# from a warm basis rather than a rebuild.

module GreedyProxy

using JuMP, Gurobi, Printf
import Main.Caps

export greedy_proxy_reduction

"""
    greedy_proxy_reduction(TR, c, epsL; kwargs...)

`TR` is the module holding `multiscenario_windows` and
`multiscenario_internal_bounds` (common/preprocessing.jl, included inside the
proxy's data module). `c` is a MultiScenarioTxReductionCase, `epsL` the per-line
flow tolerance -- scalar or one entry per line, exactly as the proxy takes it.

Every setting that shapes feasibility mirrors `solve_reduction_edge_multiscenario`
so the result is a valid warm start for it: same `epsL`, `near_limit_threshold`,
scenario sets, caps and `internal_rating_bound`. A seed built at a tighter eps is
still valid for a looser solve (windows only widen); the reverse is not.

Keyword options:
  scenario_indices        scenarios the windows are written for
  protection_indices      scenarios that decide which lines are protected
  near_limit_threshold    fraction of rating above which a line is protected
  internal_bound_scale    per-line internal-transfer bound scale
  internal_rating_bound   cap the internal transfer by the line rating too
  line_budget/hop_cap/size_cap   caps, checked incrementally
  force_internal          lines an earlier pass merged and this one keeps
  ordering                :loading (least loaded first) | :none (line order)
  feas_tol                worst single slack, p.u., still counted as feasible
  time_limit              seconds

`feas_tol` is the WORST row, not the sum. A sum grows with the model -- case118
carries 186 x 20 window rows, and LP noise of 1e-9 apiece adds up past any
threshold tight enough to mean something -- so the per-row figure is the one that
keeps its meaning as the case grows. The default is 1e-5 p.u., a kilowatt at
baseMVA = 100, which sits above Gurobi's own feasibility tolerance and far below
anything physical.
"""
function greedy_proxy_reduction(TR, c, epsL;
                               scenario_indices=axes(c.p, 2),
                               protection_indices=axes(c.p, 2),
                               near_limit_threshold=nothing,
                               internal_bound_scale::Real=3.0,
                               internal_rating_bound::Bool=true,
                               line_budget=nothing, hop_cap=nothing,
                               size_cap=nothing,
                               force_internal=Int[],
                               ordering::Symbol=:loading,
                               feas_tol::Real=1e-5,
                               time_limit::Real=600.0,
                               threads::Int=1,
                               verbose::Bool=true)
    ordering in (:loading, :none) || error("ordering must be :loading or :none")
    base = c.base
    N, Ln = base.N, base.Ln
    u, v = base.Efrom, base.Eto
    started = time()
    deadline = started + time_limit

    win = TR.multiscenario_windows(c, epsL; near_limit_threshold=near_limit_threshold,
                                   scenario_indices=scenario_indices,
                                   protection_indices=protection_indices)
    selected = Int.(collect(scenario_indices))
    philo, phiup, protected = win.philo, win.phiup, win.protected
    S = length(selected)
    S == 0 && error("no scenarios selected")
    G = TR.multiscenario_internal_bounds(c, selected, philo, phiup, win.epsv;
                                         internal_bound_scale=internal_bound_scale)

    # Same two numbers the proxy uses: how much transfer a collapsed line may
    # carry, and how far f is allowed to travel on an external one.
    gcap = [internal_rating_bound ? min(base.frate[l], G[l, s]) : G[l, s]
            for l in 1:Ln, s in 1:S]
    fbound = [max(base.frate[l], abs(philo[l, s]), abs(phiup[l, s]))
              for l in 1:Ln, s in 1:S]

    env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
    m = Model(() -> Gurobi.Optimizer(env))
    set_silent(m)
    set_optimizer_attribute(m, "Threads", threads)
    @variable(m, th[1:N, 1:S])
    @variable(m, f[1:Ln, 1:S])
    @variable(m, gint[1:Ln, 1:S])
    @variable(m, sup[1:Ln, 1:S] >= 0)      # window overshoot, both directions
    @variable(m, slo[1:Ln, 1:S] >= 0)
    @variable(m, bsp[1:N, 1:S] >= 0)       # unroutable injection at a bus
    @variable(m, bsm[1:N, 1:S] >= 0)

    @constraint(m, [s=1:S], th[base.j0, s] == 0)
    @constraint(m, [l=1:Ln, s=1:S], f[l, s] == base.Dx[l] * (th[u[l], s] - th[v[l], s]))
    outgoing = [findall(==(b), u) for b in 1:N]
    incoming = [findall(==(b), v) for b in 1:N]
    @constraint(m, [b=1:N, s=1:S],
        sum(f[l, s] + gint[l, s] for l in outgoing[b]; init=0.0) -
        sum(f[l, s] + gint[l, s] for l in incoming[b]; init=0.0) +
        bsp[b, s] - bsm[b, s] == c.p[b, selected[s]])

    # The merge switch rides on these right-hand sides. cl = 0 gives the window;
    # cl = 1 sets them to 0 and f is bounded to 0, which makes the angles equal.
    # Slack never sits on the merge itself -- only on the window it replaces --
    # so a merge cannot be bought instead of made.
    rup = @constraint(m, [l=1:Ln, s=1:S],  f[l, s] - sup[l, s] <=  phiup[l, s])
    rlo = @constraint(m, [l=1:Ln, s=1:S], -f[l, s] - slo[l, s] <= -philo[l, s])
    @objective(m, Min, sum(sup) + sum(slo) + sum(bsp) + sum(bsm))

    for l in 1:Ln, s in 1:S
        set_lower_bound(f[l, s], -fbound[l, s]); set_upper_bound(f[l, s], fbound[l, s])
        set_lower_bound(gint[l, s], 0.0);        set_upper_bound(gint[l, s], 0.0)
    end

    function set_merged!(l, on)
        for s in 1:S
            if on
                set_lower_bound(f[l, s], 0.0); set_upper_bound(f[l, s], 0.0)
                set_normalized_rhs(rup[l, s], 0.0)
                set_normalized_rhs(rlo[l, s], 0.0)
                set_lower_bound(gint[l, s], -gcap[l, s])
                set_upper_bound(gint[l, s], gcap[l, s])
            else
                set_lower_bound(f[l, s], -fbound[l, s])
                set_upper_bound(f[l, s], fbound[l, s])
                set_normalized_rhs(rup[l, s],  phiup[l, s])
                set_normalized_rhs(rlo[l, s], -philo[l, s])
                set_lower_bound(gint[l, s], 0.0); set_upper_bound(gint[l, s], 0.0)
            end
        end
    end

    lp_checks = Ref(0)
    "Total slack in p.u., split into the window part and the balance part."
    function slack()
        optimize!(m)
        lp_checks[] += 1
        termination_status(m) == MOI.OPTIMAL ||
            return (total=Inf, window=Inf, balance=Inf, worst=Inf)
        w = sum(value.(sup)) + sum(value.(slo))
        b = sum(value.(bsp)) + sum(value.(bsm))
        worst = max(maximum(value.(sup); init=0.0), maximum(value.(slo); init=0.0),
                    maximum(value.(bsp); init=0.0), maximum(value.(bsm); init=0.0))
        return (total=w + b, window=w, balance=b, worst=worst)
    end

    closure(mask) = (rep = Caps.reps(base, mask);
                     BitVector([rep[u[l]] == rep[v[l]] for l in 1:Ln]))
    nbuses(mask) = length(unique(Caps.reps(base, mask)))

    internal = falses(Ln)
    force = Int.(collect(force_internal))
    if !isempty(force)
        clash = [l for l in force if protected[l]]
        isempty(clash) || error(
            "force_internal asks to collapse protected line(s) $clash. Protection " *
            "comes from near_limit_threshold, not eps, so the threshold or the " *
            "scenario set moved between passes.")
        seed = falses(Ln); seed[force] .= true
        internal = closure(seed)
        bad = [l for l in 1:Ln if internal[l] && protected[l]]
        isempty(bad) || error("held merges close a loop over protected line(s) $bad")
        verbose && @printf(" held internal from earlier passes: %d line(s), %d buses\n",
                           count(internal), nbuses(internal))
    end

    # Keep the LP in step with `internal` so a candidate only needs its own flip.
    cur = falses(Ln)
    function sync!(mask)
        for l in 1:Ln
            mask[l] == cur[l] && continue
            set_merged!(l, mask[l]); cur[l] = mask[l]
        end
    end
    sync!(internal)
    held = slack()
    isfinite(held.total) && held.worst <= feas_tol || error(
        "the held merges are already proxy-infeasible (worst slack $(held.worst) p.u.). " *
        "A ladder rung cannot start from a set the windows do not admit.")

    loading = [maximum(abs(c.fhat[l, s]) for s in selected) / max(base.frate[l], eps())
               for l in 1:Ln]
    order = ordering === :loading ? sortperm(loading) : collect(1:Ln)

    trace = NamedTuple[]
    accepted, rejected = 0, 0
    progress = true
    while progress && time() < deadline
        progress = false
        for l in order
            time() < deadline || break
            internal[l] && continue
            protected[l] && continue
            u[l] == v[l] && continue
            prop = closure(_with(internal, l))
            # Closing a loop may sweep other lines in with it; a protected one
            # among them makes the merge illegal, exactly as the proxy's cl == 0.
            if any(prop[k] && protected[k] for k in 1:Ln)
                rejected += 1
                push!(trace, (line=l, accepted=false, rejection=:protected_closure,
                              slack_pu=NaN, internal_lines=count(internal),
                              buses=nbuses(internal), elapsed_seconds=time() - started))
                continue
            end
            if !Caps.within_caps(base, prop; budget=line_budget, hop_cap=hop_cap,
                                 size_cap=size_cap, bus=u[l])
                rejected += 1
                push!(trace, (line=l, accepted=false, rejection=:caps,
                              slack_pu=NaN, internal_lines=count(internal),
                              buses=nbuses(internal), elapsed_seconds=time() - started))
                continue
            end
            sync!(prop)
            r = slack()
            if isfinite(r.total) && r.worst <= feas_tol
                internal = prop
                accepted += 1
                progress = true
                push!(trace, (line=l, accepted=true, rejection=:none,
                              slack_pu=r.total, internal_lines=count(internal),
                              buses=nbuses(internal), elapsed_seconds=time() - started))
                verbose && @printf(" merge L%-4d -> %d buses  (%.0f s)\n", l,
                                   nbuses(internal), time() - started)
            else
                rejected += 1
                sync!(internal)
                push!(trace, (line=l, accepted=false,
                              rejection = isfinite(r.total) ? :window : :infeasible,
                              slack_pu=r.total, internal_lines=count(internal),
                              buses=nbuses(internal), elapsed_seconds=time() - started))
            end
        end
    end

    sync!(internal)
    final = slack()
    verbose && @printf(" proxy greedy: %d merged line(s), %d of %d buses, %d LP check(s), %.1f s\n",
                       count(internal), nbuses(internal), N, lp_checks[], time() - started)
    return (; internal, buses=nbuses(internal), trace,
            elapsed_seconds=time() - started, lp_checks=lp_checks[],
            accepted, rejected, protected_lines=count(protected),
            slack_pu=final.total, slack_mw=final.total * base.baseMVA,
            reason = time() < deadline ? :converged : :time_limit)
end

"Copy of `mask` with line `l` also merged."
_with(mask, l) = (p = copy(BitVector(mask)); p[l] = true; p)

end # module

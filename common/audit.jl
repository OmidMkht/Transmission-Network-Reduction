module Audit

using JuMP, Gurobi, LinearAlgebra

export topology_oracle, check_scenario!, audit_topology

# Fixed topology: reuse two sparse LPs across demands. No KKT or SOS1 blocks
# are needed to ask whether SOME lower-level optimum is original-grid feasible.
function topology_oracle(base, internal::AbstractVector{Bool};
                         relax_pmin::Bool=true, threads::Int=1,
                         cost_tolerance::Real=1e-8, flow_tolerance::Real=1e-6)
    length(internal) == base.Ln || error("internal must have one entry per line")
    all(iszero, base.c2) || error("The oracle requires the same linear costs as the master")
    cost_tolerance >= 0 && flow_tolerance >= 0 || error("tolerances must be nonnegative")
    all(isfinite, base.frate) && all(>(0), base.frate) || error("ratings must be finite and positive")
    N, L, K = base.N, base.Ln, length(base.gen_bus)
    parent = collect(1:N)
    root(x) = parent[x] == x ? x : (parent[x] = root(parent[x]))
    for l in findall(internal)
        parent[root(base.Efrom[l])] = root(base.Eto[l])
    end
    rep = [root(i) for i in 1:N]
    retained = sort(unique(rep))
    index = Dict(r => i for (i, r) in enumerate(retained))
    bus = [index[r] for r in rep]
    R = length(retained)
    ext = [l for l in 1:L if bus[base.Efrom[l]] != bus[base.Eto[l]]]
    rf = [[l for l in ext if bus[base.Efrom[l]] == i] for i in 1:R]
    rt = [[l for l in ext if bus[base.Eto[l]] == i] for i in 1:R]
    rg = [[k for k in 1:K if bus[base.gen_bus[k]] == i] for i in 1:R]
    tf = [findall(==(i), base.Efrom) for i in 1:N]
    tt = [findall(==(i), base.Eto) for i in 1:N]
    tg = [findall(==(i), base.gen_bus) for i in 1:N]
    pmin = relax_pmin ? min.(0.0, base.pmin) : base.pmin
    env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
    function make_model(safety)
        m = Model(() -> Gurobi.Optimizer(env))
        set_silent(m)
        set_optimizer_attribute(m, "Threads", threads)
        set_optimizer_attribute(m, "FeasibilityTol", 1e-9)
        set_optimizer_attribute(m, "OptimalityTol", 1e-9)
        @variable(m, pmin[k] <= g[k=1:K] <= base.pmax[k])
        @variable(m, th[1:R])
        @variable(m, f[ext])
        @constraint(m, th[bus[base.j0]] == 0)
        @constraint(m, [l=ext], f[l] == base.Dx[l] *
            (th[bus[base.Efrom[l]]] - th[bus[base.Eto[l]]]))
        @constraint(m, [l=ext], -base.frate[l] <= f[l] <= base.frate[l])
        rb = @constraint(m, [i=1:R],
            sum(f[l] for l in rf[i]; init=0.0) - sum(f[l] for l in rt[i]; init=0.0) -
            sum(g[k] for k in rg[i]; init=0.0) == 0)
        @variable(m, tth[1:N])
        @variable(m, trueflow[1:L])
        @constraint(m, tth[base.j0] == 0)
        @constraint(m, [l=1:L], trueflow[l] == base.Dx[l] *
            (tth[base.Efrom[l]] - tth[base.Eto[l]]))
        tb = @constraint(m, [i=1:N],
            sum(trueflow[l] for l in tf[i]; init=0.0) - sum(trueflow[l] for l in tt[i]; init=0.0) -
            sum(g[k] for k in tg[i]; init=0.0) == 0)
        cost = @expression(m, sum(base.c1[k] * g[k] + base.c0[k] for k in 1:K))
        @variable(m, cap)
        @constraint(m, cost <= cap)
        violation = nothing
        if safety
            violation = @variable(m, lower_bound=0)
            @constraint(m, [l=1:L], trueflow[l] <= base.frate[l] * (1 + violation))
            @constraint(m, [l=1:L], -trueflow[l] <= base.frate[l] * (1 + violation))
            @objective(m, Min, violation)
        else
            @objective(m, Min, cost)
        end
        return (; m, g, rb, tb, cost, cap, trueflow, violation)
    end
    return (; base, bus, R, lower=make_model(false), safety=make_model(true),
            cost_tolerance=Float64(cost_tolerance), flow_tolerance=Float64(flow_tolerance))
end

function _demand!(oracle, d)
    length(d) == oracle.base.N && all(isfinite, d) || error("invalid demand")
    aggregate = zeros(oracle.R)
    for i in eachindex(d)
        aggregate[oracle.bus[i]] += d[i]
    end
    for part in (oracle.lower, oracle.safety)
        for i in eachindex(aggregate)
            set_normalized_rhs(part.rb[i], -aggregate[i])
        end
        for i in eachindex(d)
            set_normalized_rhs(part.tb[i], -d[i])
        end
    end
end

function _solve!(m, deadline)
    remaining = deadline - time()
    remaining > 0 || return MOI.TIME_LIMIT
    set_time_limit_sec(m, min(remaining, 1e10))
    optimize!(m)
    return termination_status(m)
end

"""
Check existence of an original-feasible reduced optimum, to stated cost/flow
tolerances. Optional `all_optima` maximizes both directions of EVERY original
flow over the slightly widened optimal cost band, including collapsed lines.
Nonoptimal LP termination never supplies a safety certificate or a no-good cut.
"""
function check_scenario!(oracle, demand; deadline::Real=Inf, all_optima::Bool=false,
                         full_cost=nothing, cost_gap_pct=nothing)
    _demand!(oracle, demand)
    lo, safe = oracle.lower, oracle.safety
    has_upper_bound(lo.cap) && delete_upper_bound(lo.cap)
    @objective(lo.m, Min, lo.cost)
    status = _solve!(lo.m, deadline)
    empty_result = (classification=:unknown, status=status, reduced_cost=NaN,
        excess=NaN, returned_utilization=NaN, worst_optimal_utilization=NaN,
        all_optima_checked=false, all_optima_safe=false, cost_slack=NaN,
        cost_gap_pass=false)
    status == MOI.INFEASIBLE && return merge(empty_result,
        (classification=:reduced_infeasible, excess=Inf))
    status == MOI.OPTIMAL || return empty_result
    z = value(lo.cost)
    slack = oracle.cost_tolerance * max(1.0, abs(z))
    util = maximum(abs.(value.(lo.trueflow)) ./ oracle.base.frate)
    cost_ok = if isnothing(cost_gap_pct)
        true
    else
        isnothing(full_cost) && error("full_cost is required for an economic gap check")
        cost_gap_pct >= 0 || error("cost_gap_pct must be nonnegative")
        z <= full_cost + cost_gap_pct / 100 * abs(full_cost) + slack
    end
    economic_cap = isnothing(cost_gap_pct) ? Inf :
        full_cost + cost_gap_pct/100 * abs(full_cost)
    if !cost_ok
        classification = objective_bound(lo.m) > economic_cap + slack ? :fail : :unknown
        return merge(empty_result,(; classification,status=MOI.OPTIMAL,
            reduced_cost=z,excess=Inf,returned_utilization=util,cost_slack=slack,cost_gap_pass=false))
    end
    set_upper_bound(safe.cap, min(z + slack,economic_cap))
    status = _solve!(safe.m, deadline)
    result = merge(empty_result, (status=status, reduced_cost=z,
        returned_utilization=util, cost_slack=slack, cost_gap_pass=cost_ok))
    status == MOI.OPTIMAL || return result
    excess = max(0.0, value(safe.violation))
    lower_excess = objective_bound(safe.m)
    classification = if !cost_ok
        :fail
    elseif excess <= oracle.flow_tolerance
        :pass
    elseif isfinite(lower_excess) && lower_excess > oracle.flow_tolerance
        :fail
    else
        :unknown
    end
    result = merge(result, (; classification, excess))
    all_optima || return result
    set_upper_bound(lo.cap, z + slack)
    worst = 0.0
    for l in 1:oracle.base.Ln, sign in (-1, 1)
        @objective(lo.m, Max, sign * lo.trueflow[l])
        status = _solve!(lo.m, deadline)
        status == MOI.OPTIMAL || return merge(result, (status=status,))
        bound = objective_bound(lo.m)
        isfinite(bound) || return result
        worst = max(worst, bound / oracle.base.frate[l])
    end
    return merge(result, (worst_optimal_utilization=worst,
        all_optima_checked=true, all_optima_safe=worst <= 1 + oracle.flow_tolerance))
end

"Check a finite pool; this pool is design data, never a held-out test set."
function audit_topology(c, internal, scope; relax_pmin=true, threads=1,
                        cost_tolerance=1e-8, flow_tolerance=1e-6,
                        deadline=Inf, all_optima=false, full_costs=nothing,
                        cost_gap_pct=nothing)
    oracle = topology_oracle(c.base, internal; relax_pmin, threads,
                             cost_tolerance, flow_tolerance)
    rows = NamedTuple[]
    for s in scope
        a = check_scenario!(oracle, c.load[:, s]; deadline, all_optima,
            full_cost=isnothing(full_costs) ? nothing : full_costs[s], cost_gap_pct)
        push!(rows, merge((scenario=s, scenario_id=c.scenario_ids[s]), a))
    end
    covered = !isempty(rows) && all(r -> r.classification == :pass, rows)
    return (; rows, covered,
        all_optima_covered=covered && all(r -> r.all_optima_checked && r.all_optima_safe, rows),
        guarantee=:optimistic_finite_pool, cost_tolerance, flow_tolerance)
end

end # module

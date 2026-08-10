# --------------------------------------------------------------------------- #
# POSTPROCESSING AND VALIDATION.
# --------------------------------------------------------------------------- #

# ============== 1. RECOVERING THE CLUSTERING FROM LINE STATUS =============== #
# Recover the assignment matrix from a line status vector
function assignment_from_line_status(c::TxReductionCase, cval; tol::Float64=0.5)
    N = c.N
    g = SimpleGraph(N)
    for l in 1:c.Ln
        cval[l] > tol && add_edge!(g, c.Efrom[l], c.Eto[l])
    end
    A = zeros(Int, N, N)
    for comp in connected_components(g)
        rep = c.j0 in comp ? c.j0 : minimum(comp)
        for j in comp
            A[rep, j] = 1
        end
    end
    return A
end

# Extract the reduction from a solved A: retained buses and j -> representative i.
function extract_reduction(Aval::AbstractMatrix)
    N = size(Aval, 1)
    retained = [i for i in 1:N if Aval[i, i] > 0.5]
    rep_of = [findfirst(i -> Aval[i, j] > 0.5, 1:N) for j in 1:N]
    return (retained=retained, rep_of=rep_of, n_retained=length(retained))
end

# ================== 2. FEASIBILITY OF THE RETURNED POINT ==================== #
"""Check the solver's unrounded point against every active model row."""
function model_feasibility_check_multiscenario(c::MultiScenarioTxReductionCase, r;
                                               tol::Float64=1e-6)
    base = c.base
    E = incidence_matrix(base)
    cr, f, g, th = r.c_raw, r.f, r.gint, r.vartheta
    selected = r.scenario_indices
    external = findall(cr .< 0.5)
    internal = findall(cr .>= 0.5)

    integrality = maximum(abs.(cr .- round.(cr)))
    objective_raw = sum(cr)
    objective_residue = abs(objective_raw - round(objective_raw))
    reference = maximum(abs.(th[base.j0, :]))
    flow_definition = maximum(abs.(f .-
        base.Dx .* (th[base.Efrom, :] .- th[base.Eto, :])))
    balance = maximum(abs.(E * (f + g) - c.p[:, selected]))
    window = isempty(external) ? 0.0 : maximum(max.(
        r.philo[external, :] .- f[external, :],
        f[external, :] .- r.phiup[external, :], 0.0))
    internal_flow = isempty(internal) ? 0.0 : maximum(abs.(f[internal, :]))
    external_transfer = isempty(external) ? 0.0 : maximum(abs.(g[external, :]))
    protected = isempty(findall(r.protected)) ? 0.0 : maximum(abs.(cr[findall(r.protected)]))
    # r.G is per (line, scenario) -- same shape as g, so no reshape is needed.
    g_bound = maximum(max.(abs.(g) .- r.G, 0.0))
    worst = maximum((reference, flow_definition, balance, window,
                     internal_flow, external_transfer, protected, g_bound))
    objective_is_integer = objective_residue <= 1e-12
    return (; integrality, objective_raw, objective_residue,
            objective_is_integer, reference, flow_definition, balance, window,
            internal_flow, external_transfer, protected, g_bound, worst,
            genuine=objective_is_integer && worst <= tol, tol)
end

function print_model_feasibility_check_multiscenario(chk)
    println("   objective (raw sum of c_l) = ", chk.objective_raw)
    println("   objective is exact integer = ", chk.objective_is_integer)
    println("   integrality residual       = ", chk.integrality)
    println("   flow-definition residual   = ", chk.flow_definition)
    println("   nodal-balance residual      = ", chk.balance)
    println("   external-window violation  = ", chk.window)
    println("   internal-flow violation     = ", chk.internal_flow)
    println("   external-transfer violation = ", chk.external_transfer)
    println("   protected-line violation    = ", chk.protected)
    println("   worst model residual        = ", chk.worst)
    println("   GENUINE                     = ", chk.genuine)
end

# =============== 3. WINDOW SCREENING AND MONTHLY BENCHMARK ================== #
"""
Recompute a clustering's reduced flows for arbitrary scenarios and report its
window violations.
"""
function screen_reduction_scenarios(c::MultiScenarioTxReductionCase, Aval, epsL;
                                    scenario_indices=axes(c.p, 2),
                                    protection_indices=axes(c.p, 2),
                                    near_limit_threshold=nothing,
                                    congestion_relaxation=0.0,
                                    congestion_relaxation_mode::Symbol=:none,
                                    tolerance::Float64=1e-7)
    base = c.base
    selected = Int.(collect(scenario_indices))
    # Screen against the SAME windows the MILP was given, so a window violation
    # means the clustering failed its own model rather than merely failing the
    # unrelaxed pin. True deviation from the base point is reported separately by
    # benchmark_reduction_scenarios, which compares flow to fhat directly and is
    # unaffected by any relaxation.
    win = multiscenario_windows(c, epsL;
        scenario_indices=selected,
        protection_indices=protection_indices,
        near_limit_threshold=near_limit_threshold,
        congestion_relaxation=congestion_relaxation,
        congestion_relaxation_mode=congestion_relaxation_mode)
    A = round.(Int, Aval)
    red = extract_reduction(A)
    cl = [red.rep_of[base.Efrom[l]] == red.rep_of[base.Eto[l]] ? 1 : 0 for l in 1:base.Ln]

    E = incidence_matrix(base)
    Bred = A * E * Diagonal(base.Dx) * E' * A'
    retained = red.retained
    ref = red.rep_of[base.j0]
    free = setdiff(retained, [ref])
    theta = zeros(base.N, length(selected))
    if !isempty(free)
        theta[free, :] = factorize(Bred[free, free]) \ (A * c.p[:, selected])[free, :]
    end
    flow = Diagonal(base.Dx) * E' * A' * theta

    external = findall(==(0), cl)
    line_violation = zeros(base.Ln, length(selected))
    # SIGNED violation, kept separately. The magnitude alone cannot tell scenario
    # generation which side of the window was breached, and the window is two
    # half-spaces: a support point for "flow too high" constrains nothing about
    # "flow too low". Adding only one side would just make the next iteration
    # rediscover the other.
    violation_above = zeros(base.Ln, length(selected))
    violation_below = zeros(base.Ln, length(selected))
    for l in external
        violation_above[l, :] = max.(flow[l, :] .- win.phiup[l, :], 0.0)
        violation_below[l, :] = max.(win.philo[l, :] .- flow[l, :], 0.0)
        line_violation[l, :] = max.(violation_above[l, :], violation_below[l, :])
    end
    protected_internal = findall(l -> win.protected[l] && cl[l] == 1, 1:base.Ln)
    if !isempty(protected_internal)
        line_violation[protected_internal, :] .= Inf
        violation_above[protected_internal, :] .= Inf
        violation_below[protected_internal, :] .= Inf
    end
    scenario_violation = vec(maximum(line_violation, dims=1))
    worst_value, worst_local = findmax(scenario_violation)
    worst_line = argmax(line_violation[:, worst_local])
    violating_local = findall(>(tolerance), scenario_violation)

    return (
        feasible=isempty(violating_local), flow=flow,
        scenario_indices=selected, scenario_ids=c.scenario_ids[selected],
        # The windows actually used, so a caller can measure overload against the
        # ADJUSTED cap (rating + delta on a relaxed congested line) as well as
        # against the plain rating.
        philo=win.philo, phiup=win.phiup, internal=cl,
        scenario_violation=scenario_violation,
        line_violation=line_violation,
        violation_above=violation_above, violation_below=violation_below,
        violating_scenarios=selected[violating_local],
        worst_violation=worst_value,
        worst_scenario=selected[worst_local],
        worst_scenario_id=c.scenario_ids[selected[worst_local]],
        worst_line=worst_line,
        protected_internal=protected_internal,
        tolerance=tolerance,
    )
end

"""Benchmark one clustering against every supplied monthly scenario."""
function benchmark_reduction_scenarios(c::MultiScenarioTxReductionCase,
                                       Aval, epsL, scenario_indices;
                                       congestion_threshold::Float64=0.9999,
                                       near_limit_threshold=congestion_threshold,
                                       protection_indices=scenario_indices,
                                       congestion_relaxation=0.0,
                                       congestion_relaxation_mode::Symbol=:none,
                                       tolerance::Float64=1e-7)
    scenario_indices = Int.(collect(scenario_indices))
    screen = screen_reduction_scenarios(c, Aval, epsL;
        scenario_indices=scenario_indices,
        protection_indices=protection_indices,
        near_limit_threshold=near_limit_threshold,
        congestion_relaxation=congestion_relaxation,
        congestion_relaxation_mode=congestion_relaxation_mode,
        tolerance=tolerance)
    base = c.base
    A = round.(Int, Aval)
    rep = extract_reduction(A).rep_of
    internal = [rep[base.Efrom[l]] == rep[base.Eto[l]] for l in 1:base.Ln]
    external = .!internal

    normalized_error = abs.(screen.flow - c.fhat[:, scenario_indices]) ./ base.frate
    absolute_error = abs.(screen.flow - c.fhat[:, scenario_indices])
    external_error = copy(normalized_error)
    external_error[internal, :] .= 0.0
    max_external_error = vec(maximum(external_error, dims=1))
    max_all_line_error = vec(maximum(normalized_error, dims=1))
    rating_overload = max.(abs.(screen.flow) .- base.frate, 0.0)

    # Overload of the SURVIVING (external) lines in the reduced network, scored
    # two ways. Internal lines are excluded because they carry no flow in the
    # reduced network by construction, so counting them would be meaningless.
    #
    #   vs RATING    the plain thermal limit -- a real violation either way.
    #   vs ADJUSTED  the cap the relaxation actually granted (rating + delta on a
    #                relaxed congested line). A line past its rating but inside
    #                its adjusted cap is doing exactly what the model permitted;
    #                one past the ADJUSTED cap has escaped even that, which means
    #                the clustering, not the relaxation, is at fault.
    absolute_flow = abs.(screen.flow)
    cap_rating = repeat(base.frate, 1, length(scenario_indices))
    cap_adjusted = max.(cap_rating, abs.(screen.philo), abs.(screen.phiup))
    overload_vs_rating = max.(absolute_flow .- cap_rating, 0.0)
    overload_vs_adjusted = max.(absolute_flow .- cap_adjusted, 0.0)
    overload_vs_rating[internal, :] .= 0.0
    overload_vs_adjusted[internal, :] .= 0.0
    over_rating = overload_vs_rating .> tolerance
    over_adjusted = overload_vs_adjusted .> tolerance
    lines_over_rating = findall(vec(any(over_rating, dims=2)))
    lines_over_adjusted = findall(vec(any(over_adjusted, dims=2)))
    scenarios_over_rating = findall(vec(any(over_rating, dims=1)))
    scenarios_over_adjusted = findall(vec(any(over_adjusted, dims=1)))

    utilization = abs.(c.fhat[:, scenario_indices]) ./ base.frate
    congested_mask = utilization .>= congestion_threshold
    congested_lines = findall(vec(any(congested_mask, dims=2)))
    congested_internal = [l for l in congested_lines if internal[l]]
    congested_errors = normalized_error[congested_mask]
    window_feasible = screen.scenario_violation .<= tolerance
    sorted_external = sort(max_external_error)
    p95_index = max(1, ceil(Int, 0.95 * length(sorted_external)))
    sorted_violation = sort(screen.scenario_violation)
    violation_p95 = sorted_violation[p95_index]
    relative_window_violation = screen.line_violation ./ base.frate

    return (
        screen=screen,
        n_scenarios=length(scenario_indices),
        n_window_feasible=count(window_feasible),
        window_feasible_fraction=mean(window_feasible),
        worst_window_violation=screen.worst_violation,
        mean_window_violation=mean(screen.scenario_violation),
        p95_window_violation=violation_p95,
        max_relative_window_violation=maximum(relative_window_violation),
        violating_scenario_line_pairs=count(>(tolerance), screen.line_violation),
        max_external_normalized_error=maximum(max_external_error),
        mean_external_normalized_error=mean(max_external_error),
        p95_external_normalized_error=sorted_external[p95_index],
        max_all_line_normalized_error=maximum(max_all_line_error),
        max_absolute_flow_error=maximum(absolute_error),
        max_rating_overload=maximum(rating_overload),
        # External-line overload, vs the plain rating and vs the adjusted cap.
        overload_vs_rating=overload_vs_rating,
        overload_vs_adjusted=overload_vs_adjusted,
        n_lines_over_rating=length(lines_over_rating),
        n_lines_over_adjusted=length(lines_over_adjusted),
        lines_over_rating=lines_over_rating,
        lines_over_adjusted=lines_over_adjusted,
        n_scenarios_over_rating=length(scenarios_over_rating),
        n_scenarios_over_adjusted=length(scenarios_over_adjusted),
        scenarios_over_rating=scenario_indices[scenarios_over_rating],
        scenarios_over_adjusted=scenario_indices[scenarios_over_adjusted],
        n_pairs_over_rating=count(over_rating),
        n_pairs_over_adjusted=count(over_adjusted),
        max_overload_vs_rating=maximum(overload_vs_rating),
        max_overload_vs_adjusted=maximum(overload_vs_adjusted),
        overload_hours_by_line=vec(sum(over_rating, dims=2)),
        congested_lines=congested_lines,
        congested_internal=congested_internal,
        max_congested_normalized_error=isempty(congested_errors) ? 0.0 : maximum(congested_errors),
        normalized_error=normalized_error,
        absolute_error=absolute_error,
        relative_window_violation=relative_window_violation,
        max_external_error_by_scenario=max_external_error,
        max_all_line_error_by_scenario=max_all_line_error,
        window_feasible_by_scenario=window_feasible,
    )
end

# ====================== 4. DC-OPF VALIDATION LAYER ========================= #
"""Solve one scenario's DC-OPF on the original (`A = I`) or reduced network."""
function solve_assigned_dc_opf_scenario(c::MultiScenarioTxReductionCase,
                                        Aval, scenario_index::Int;
                                        relax_pmin::Bool=true,
                                        time_limit=nothing,
                                        optimizer=Gurobi.Optimizer)
    base = c.base
    N, Ln = base.N, base.Ln
    1 <= scenario_index <= length(c.scenario_ids) ||
        error("scenario_index is outside the loaded scenario range")
    A = round.(Int, Aval)
    red = extract_reduction(A)
    retained, rep = red.retained, red.rep_of
    u = [rep[base.Efrom[l]] for l in 1:Ln]
    v = [rep[base.Eto[l]] for l in 1:Ln]
    gen_bus = [rep[base.gen_bus[g]] for g in eachindex(base.gen_bus)]
    demand = A * c.load[:, scenario_index]
    pmin = relax_pmin ? min.(0.0, base.pmin) : base.pmin
    nG = length(base.gen_bus)

    model = Model(optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "MIPGap", 0.00)
    isnothing(time_limit) || set_optimizer_attribute(model, "TimeLimit", time_limit)
    @variable(model, theta[1:N])
    @variable(model, pG[1:nG])
    @variable(model, flow[1:Ln])
    for i in 1:N
        i in retained || fix(theta[i], 0.0; force=true)
    end
    @constraint(model, theta[rep[base.j0]] == 0)
    @constraint(model, [g=1:nG], pmin[g] <= pG[g] <= base.pmax[g])
    @constraint(model, [l=1:Ln],
        flow[l] == base.Dx[l] * (theta[u[l]] - theta[v[l]]))
    @constraint(model, [l=1:Ln], -base.frate[l] <= flow[l] <= base.frate[l])

    balance = Dict{Int,ConstraintRef}()
    for i in retained
        balance[i] = @constraint(model,
            sum(pG[g] for g in 1:nG if gen_bus[g] == i; init=0.0) - demand[i] ==
            sum(flow[l] for l in 1:Ln if u[l] == i; init=0.0) -
            sum(flow[l] for l in 1:Ln if v[l] == i; init=0.0))
    end
    @objective(model, Min,
        sum(base.c2[g] * pG[g]^2 + base.c1[g] * pG[g] + base.c0[g]
            for g in 1:nG))
    optimize!(model)
    status = termination_status(model)
    optimal = status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED,
                         MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED)
    if !has_values(model)
        return (feasible=false, optimal=false, status=status, lmp=nothing)
    end

    lmp = nothing
    if optimal && has_duals(model)
        lmp = zeros(N)
        for i in retained
            lmp[i] = dual(balance[i]) / base.baseMVA
        end
        for j in 1:N
            lmp[j] = lmp[rep[j]]
        end
    end
    return (
        feasible=true, optimal=optimal, status=status, objective=objective_value(model),
        pG=value.(pG), theta=value.(theta), flow=value.(flow), lmp=lmp,
        solve_time=solve_time(model),
    )
end

"""Check a reduced-network dispatch on the original network for one scenario."""
function check_dispatch_on_original_scenario(c::MultiScenarioTxReductionCase,
                                             pG, scenario_index::Int;
                                             relative_tolerance::Float64=1e-3)
    base = c.base
    E = incidence_matrix(base)
    B = E * Diagonal(base.Dx) * E'
    injection = -copy(c.load[:, scenario_index])
    for g in eachindex(base.gen_bus)
        injection[base.gen_bus[g]] += pG[g]
    end
    theta = zeros(base.N)
    free = setdiff(1:base.N, [base.j0])
    theta[free] = B[free, free] \ injection[free]
    flow = Diagonal(base.Dx) * E' * theta
    balance_error = injection - E * flow
    overload = max.(abs.(flow) - base.frate, 0.0)
    relative_overload = overload ./ base.frate
    power_scale = max(maximum(abs.(injection)), 1e-8)
    relative_balance_error = maximum(abs.(balance_error)) / power_scale
    violating_lines = findall(relative_overload .> relative_tolerance)
    # A graded severity measure to sit beside the binary pass/fail. Once a
    # relaxation is in play the strict test is expected to fail, and "failed" on
    # its own says nothing about whether the dispatch missed by 0.05% or by 20%.
    # max_utilization is the worst |flow| / rating over the network, so
    # max_utilization - 1 is exactly the uniform rating headroom this dispatch
    # would have needed to be feasible.
    max_utilization = maximum(abs.(flow) ./ base.frate)
    return (
        feasible=relative_balance_error <= relative_tolerance && isempty(violating_lines),
        theta=theta, flow=flow, balance_error=balance_error,
        overload=overload, relative_overload=relative_overload,
        max_balance_error=maximum(abs.(balance_error)),
        max_relative_balance_error=relative_balance_error,
        max_overload=maximum(overload),
        max_relative_overload=maximum(relative_overload),
        total_overload=sum(overload),
        max_utilization=max_utilization,
        n_violating_lines=length(violating_lines),
        violating_lines=violating_lines,
    )
end

"""
Cheapest way to make a reduced-network dispatch work on the FULL network:
minimise the total generation move away from `pG_reference` subject to the full
network's DC power flow, generator limits and thermal ratings.

This is the fair economic comparison when the strict overload test fails. A
reduced dispatch that is 2% overloaded is not simply "infeasible" -- it is a plan
that costs a certain amount to repair, and that repair cost (in MW moved, and in
dollars against the true optimum) is the number worth comparing across
relaxation modes. `feasible=false` here means the FULL network cannot serve this
scenario at all, which is a different and much more serious failure.
"""
function repair_dispatch_on_original_scenario(c::MultiScenarioTxReductionCase,
                                              pG_reference, scenario_index::Int;
                                              relax_pmin::Bool=true,
                                              time_limit=nothing,
                                              optimizer=Gurobi.Optimizer)
    base = c.base
    N, Ln = base.N, base.Ln
    nG = length(base.gen_bus)
    pmin = relax_pmin ? min.(0.0, base.pmin) : base.pmin

    model = Model(optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "MIPGap", 0.00)
    isnothing(time_limit) || set_optimizer_attribute(model, "TimeLimit", time_limit)
    @variable(model, theta[1:N])
    @variable(model, pG[1:nG])
    @variable(model, flow[1:Ln])
    @variable(model, move[1:nG] >= 0)           # |pG - pG_reference|, linearised
    @constraint(model, theta[base.j0] == 0)
    @constraint(model, [g=1:nG], pmin[g] <= pG[g] <= base.pmax[g])
    @constraint(model, [g=1:nG],  pG[g] - pG_reference[g] <= move[g])
    @constraint(model, [g=1:nG], -pG[g] + pG_reference[g] <= move[g])
    @constraint(model, [l=1:Ln],
        flow[l] == base.Dx[l] * (theta[base.Efrom[l]] - theta[base.Eto[l]]))
    @constraint(model, [l=1:Ln], -base.frate[l] <= flow[l] <= base.frate[l])
    demand = c.load[:, scenario_index]
    @constraint(model, [i=1:N],
        sum(pG[g] for g in 1:nG if base.gen_bus[g] == i; init=0.0) - demand[i] ==
        sum(flow[l] for l in 1:Ln if base.Efrom[l] == i; init=0.0) -
        sum(flow[l] for l in 1:Ln if base.Eto[l] == i; init=0.0))
    @objective(model, Min, sum(move))
    optimize!(model)

    status = termination_status(model)
    optimal = status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED,
                         MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED)
    if !(optimal && has_values(model))
        return (feasible=false, status=status, redispatch=NaN, repaired_cost=NaN)
    end
    pGv = value.(pG)
    repaired_cost = sum(base.c2[g] * pGv[g]^2 + base.c1[g] * pGv[g] + base.c0[g]
                        for g in 1:nG)
    return (feasible=true, status=status, redispatch=sum(value.(move)),
            repaired_cost=repaired_cost, pG=pGv)
end

"""
Run full and reduced DC-OPFs for each supplied scenario and validate every
reduced dispatch back on the original network.
"""
function validate_reduced_dcopf_scenarios(c::MultiScenarioTxReductionCase,
                                          Aval, scenario_indices;
                                          relax_pmin::Bool=true,
                                          time_limit=nothing,
                                          relative_tolerance::Float64=1e-3,
                                          objective_tolerance_pct::Float64=0.1,
                                          lmp_tolerance::Float64=1e-3,
                                          measure_repair::Bool=false,
                                          training_indices=Int[],
                                          binding_tolerance::Float64=1e-3,
                                          progress_every::Int=100)
    scenarios = Int.(collect(scenario_indices))
    H = length(scenarios)
    H > 0 || error("At least one DC-OPF validation scenario is required")
    A = round.(Int, Aval)
    base = c.base
    # Which lines the clustering collapsed. An overload on an INTERNAL line and
    # one on an EXTERNAL line are different diagnoses: an internal line was
    # merged away, so the reduced network never modelled it at all; an external
    # line survived, so the reduced network modelled it and got the flow wrong.
    rep_of = extract_reduction(A).rep_of
    internal_line = [rep_of[base.Efrom[l]] == rep_of[base.Eto[l]] for l in 1:base.Ln]
    training_set = Set(Int.(collect(training_indices)))
    in_training = [s in training_set for s in scenarios]
    identity_assignment = Matrix{Int}(I, c.base.N, c.base.N)
    env_parameters = Dict{String,Any}("OutputFlag" => 0)
    isnothing(time_limit) || (env_parameters["TimeLimit"] = time_limit)
    env = Gurobi.Env(env_parameters)
    optimizer = () -> Gurobi.Optimizer(env)

    original_feasible = falses(H)
    reduced_feasible = falses(H)
    original_optimal = falses(H)
    reduced_optimal = falses(H)
    dispatch_feasible = falses(H)
    objective_within_tolerance = falses(H)
    original_status = Vector{String}(undef, H)
    reduced_status = Vector{String}(undef, H)
    original_objective = fill(NaN, H)
    reduced_objective = fill(NaN, H)
    objective_change_pct = fill(NaN, H)
    max_overload = fill(NaN, H)
    max_relative_overload = fill(NaN, H)
    max_balance_error = fill(NaN, H)
    max_lmp_error = fill(NaN, H)
    load_weighted_lmp_error = fill(NaN, H)
    # Graded severity, for when the binary overload test fails (the expected
    # case once a congestion relaxation is switched on).
    max_utilization = fill(NaN, H)
    total_overload = fill(NaN, H)
    n_violating_lines = fill(0, H)
    repair_redispatch = fill(NaN, H)
    repair_cost_pct = fill(NaN, H)
    repair_feasible = falses(H)
    # Failure anatomy: which KIND of line failed, how often each line failed, and
    # what the failing hours look like compared with the passing ones.
    n_violating_internal = zeros(Int, H)
    n_violating_external = zeros(Int, H)
    total_load = fill(NaN, H)
    line_violation_hours = zeros(Int, base.Ln)
    line_worst_overload = zeros(base.Ln)
    # Spuriously binding lines: at their rating in the reduced dispatch while the
    # true full-network optimum leaves them slack -- the opposite failure from an
    # overload (reduction being needlessly pessimistic, forcing a costlier
    # dispatch with no feasibility violation to show for it).
    n_spurious_binding = zeros(Int, H)
    line_spurious_hours = zeros(Int, base.Ln)
    line_worst_reduced_utilization = zeros(base.Ln)

    for (h, s) in enumerate(scenarios)
        original = solve_assigned_dc_opf_scenario(c, identity_assignment, s;
            relax_pmin=relax_pmin, time_limit=nothing, optimizer=optimizer)
        reduced = solve_assigned_dc_opf_scenario(c, A, s;
            relax_pmin=relax_pmin, time_limit=nothing, optimizer=optimizer)
        original_status[h] = string(original.status)
        reduced_status[h] = string(reduced.status)
        original_feasible[h] = original.feasible
        reduced_feasible[h] = reduced.feasible
        original_optimal[h] = original.optimal
        reduced_optimal[h] = reduced.optimal

        if original.optimal && reduced.optimal
            original_objective[h] = original.objective
            reduced_objective[h] = reduced.objective
            objective_change_pct[h] = 100 * (reduced.objective - original.objective) /
                max(abs(original.objective), 1e-8)
            check = check_dispatch_on_original_scenario(c, reduced.pG, s;
                relative_tolerance=relative_tolerance)
            dispatch_feasible[h] = check.feasible
            max_overload[h] = check.max_overload
            max_relative_overload[h] = check.max_relative_overload
            max_balance_error[h] = check.max_balance_error
            max_utilization[h] = check.max_utilization
            total_overload[h] = check.total_overload
            n_violating_lines[h] = check.n_violating_lines
            total_load[h] = sum(c.load[:, s])
            for l in check.violating_lines
                internal_line[l] ? (n_violating_internal[h] += 1) :
                                   (n_violating_external[h] += 1)
                line_violation_hours[l] += 1
            end
            line_worst_overload .= max.(line_worst_overload, check.overload)

            # Spuriously binding lines. An INTERNAL line is skipped: the reduced
            # network gives it flow == 0 identically, so it can never appear
            # binding there and its true loading is a separate question already
            # covered by the overload split above.
            for l in 1:base.Ln
                internal_line[l] && continue
                util_red = abs(reduced.flow[l]) / base.frate[l]
                util_full = abs(original.flow[l]) / base.frate[l]
                line_worst_reduced_utilization[l] =
                    max(line_worst_reduced_utilization[l], util_red)
                if util_red >= 1 - binding_tolerance && util_full < 1 - binding_tolerance
                    n_spurious_binding[h] += 1
                    line_spurious_hours[l] += 1
                end
            end
            # Cost accuracy and dispatch feasibility are DIFFERENT questions and
            # are now scored separately. Gating this on check.feasible (as it was)
            # made every overloaded hour also count as a cost failure, so the two
            # could never be compared independently -- which is exactly what is
            # needed to judge a relaxation that is expected to overload slightly.
            objective_within_tolerance[h] =
                abs(objective_change_pct[h]) <= objective_tolerance_pct

            if measure_repair && !check.feasible
                repair = repair_dispatch_on_original_scenario(c, reduced.pG, s;
                    relax_pmin=relax_pmin, time_limit=nothing, optimizer=optimizer)
                repair_feasible[h] = repair.feasible
                if repair.feasible
                    repair_redispatch[h] = repair.redispatch
                    repair_cost_pct[h] = 100 * (repair.repaired_cost - original.objective) /
                        max(abs(original.objective), 1e-8)
                end
            elseif measure_repair
                # Already feasible: nothing to repair, so the honest entries are
                # zero movement and the reduced dispatch's own cost gap.
                repair_feasible[h] = true
                repair_redispatch[h] = 0.0
                repair_cost_pct[h] = objective_change_pct[h]
            end

            if original.lmp !== nothing && reduced.lmp !== nothing
                price_error = abs.(reduced.lmp - original.lmp)
                max_lmp_error[h] = maximum(price_error)
                demand = c.load[:, s]
                load_weighted_lmp_error[h] = sum(demand) > 1e-12 ?
                    sum(demand .* price_error) / sum(demand) : 0.0
            end
        end
        if h == 1 || h % progress_every == 0 || h == H
            println("DC-OPF validation: $h/$H scenarios")
        end
    end

    finite_max(v) = any(isfinite, v) ? maximum(filter(isfinite, v)) : NaN
    finite_maxabs(v) = any(isfinite, v) ? maximum(abs, filter(isfinite, v)) : NaN
    all_dcopf_feasible = original_optimal .& reduced_optimal
    lmp_within_tolerance = isfinite.(max_lmp_error) .& (max_lmp_error .<= lmp_tolerance)
    return (
        scenario_indices=scenarios, scenario_ids=c.scenario_ids[scenarios],
        original_feasible=original_feasible, reduced_feasible=reduced_feasible,
        original_optimal=original_optimal, reduced_optimal=reduced_optimal,
        all_dcopf_feasible=all_dcopf_feasible,
        dispatch_feasible=dispatch_feasible,
        objective_within_tolerance=objective_within_tolerance,
        original_status=original_status, reduced_status=reduced_status,
        original_objective=original_objective, reduced_objective=reduced_objective,
        objective_change_pct=objective_change_pct,
        max_overload=max_overload, max_relative_overload=max_relative_overload,
        max_balance_error=max_balance_error,
        max_lmp_error=max_lmp_error,
        load_weighted_lmp_error=load_weighted_lmp_error,
        lmp_within_tolerance=lmp_within_tolerance,
        n_dcopf_feasible=count(all_dcopf_feasible),
        n_dispatch_feasible=count(dispatch_feasible),
        n_objective_within_tolerance=count(objective_within_tolerance),
        n_lmp_within_tolerance=count(lmp_within_tolerance),
        worst_abs_objective_change_pct=finite_maxabs(objective_change_pct),
        worst_overload=finite_max(max_overload),
        worst_relative_overload=finite_max(max_relative_overload),
        worst_lmp_error=finite_max(max_lmp_error),
        # Graded severity, comparable across relaxation modes even when the
        # strict pass/fail is uniformly "fail".
        max_utilization=max_utilization,
        total_overload=total_overload,
        n_violating_lines=n_violating_lines,
        # Failure anatomy
        internal_line=internal_line,
        in_training=in_training,
        n_violating_internal=n_violating_internal,
        n_violating_external=n_violating_external,
        line_violation_hours=line_violation_hours,
        line_worst_overload=line_worst_overload,
        total_load=total_load,
        # Spuriously binding: at the rating in the reduced dispatch, slack in the
        # true full-network optimum. The pessimistic counterpart of an overload.
        n_spurious_binding=n_spurious_binding,
        line_spurious_hours=line_spurious_hours,
        line_worst_reduced_utilization=line_worst_reduced_utilization,
        n_hours_with_spurious_binding=count(>(0), n_spurious_binding),
        max_spurious_binding=isempty(n_spurious_binding) ? 0 : maximum(n_spurious_binding),
        n_spurious_lines=count(>(0), line_spurious_hours),
        binding_tolerance=binding_tolerance,
        worst_max_utilization=finite_max(max_utilization),
        mean_max_utilization=any(isfinite, max_utilization) ?
            mean(filter(isfinite, max_utilization)) : NaN,
        worst_total_overload=finite_max(total_overload),
        max_violating_lines=isempty(n_violating_lines) ? 0 : maximum(n_violating_lines),
        # Repair economics (only populated when measure_repair=true).
        measure_repair=measure_repair,
        repair_feasible=repair_feasible,
        repair_redispatch=repair_redispatch,
        repair_cost_pct=repair_cost_pct,
        n_repair_feasible=count(repair_feasible),
        worst_repair_redispatch=finite_max(repair_redispatch),
        mean_repair_redispatch=any(isfinite, repair_redispatch) ?
            mean(filter(isfinite, repair_redispatch)) : NaN,
        worst_repair_cost_pct=finite_maxabs(repair_cost_pct),
        mean_repair_cost_pct=any(isfinite, repair_cost_pct) ?
            mean(filter(isfinite, repair_cost_pct)) : NaN,
        objective_tolerance_pct=objective_tolerance_pct,
        relative_tolerance=relative_tolerance,
        lmp_tolerance=lmp_tolerance,
        relax_pmin=relax_pmin,
    )
end

# --------------------------------------------------------------------------- #
# Canonical internal flows
# --------------------------------------------------------------------------- #
function shorted_internal_flows(c::TxReductionCase, Aval; short_factor::Float64=1e8)
    A = round.(Int, Aval)
    rep = extract_reduction(A).rep_of
    E = Matrix(incidence_matrix(c))
    scale = maximum(c.Dx) * short_factor
    Dxs = [rep[c.Efrom[l]] == rep[c.Eto[l]] ? scale : c.Dx[l] for l in 1:c.Ln]
    theta = pinv(E * Diagonal(Dxs) * E') * c.p
    return Diagonal(Dxs) * E' * theta          # total flow per line
end

# ================== 6. SOLVE-TIME BENCHMARK (full vs reduced) ============== #
"""Build one DC-OPF over the bus set `retained`, mapping every bus through `rep`."""
function _dcopf_on_partition(c::MultiScenarioTxReductionCase, scenario_index::Int,
                             rep::AbstractVector{Int}, retained::AbstractVector{Int},
                             keep_lines::AbstractVector{Int};
                             relax_pmin::Bool=true, threads::Int=1,
                             time_limit=nothing, optimizer=Gurobi.Optimizer,
                             direct::Bool=true)
    base = c.base
    nb = length(retained)
    idx = Dict{Int,Int}(b => k for (k, b) in enumerate(retained))
    demand = zeros(nb)
    for i in 1:base.N
        demand[idx[rep[i]]] += c.load[i, scenario_index]
    end
    gbus = [idx[rep[b]] for b in base.gen_bus]
    nG = length(base.gen_bus)
    pmin = relax_pmin ? min.(0.0, base.pmin) : base.pmin
    uu = [idx[rep[base.Efrom[l]]] for l in keep_lines]
    vv = [idx[rep[base.Eto[l]]] for l in keep_lines]
    nL = length(keep_lines)


    model = direct ? direct_model(optimizer()) : Model(optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "Threads", threads)
    set_optimizer_attribute(model, "MIPGap", 0.00)
    isnothing(time_limit) || set_optimizer_attribute(model, "TimeLimit", time_limit)

    frate_k = [base.frate[l] for l in keep_lines]
    @variable(model, theta[1:nb])
    @variable(model, pmin[g] <= pG[g=1:nG] <= base.pmax[g])
    @variable(model, -frate_k[k] <= flow[k=1:nL] <= frate_k[k])
    @constraint(model, theta[idx[rep[base.j0]]] == 0)
    @constraint(model, [k=1:nL],
        flow[k] == base.Dx[keep_lines[k]] * (theta[uu[k]] - theta[vv[k]]))
    outg = [Int[] for _ in 1:nb]
    inc  = [Int[] for _ in 1:nb]
    for k in 1:nL
        push!(outg[uu[k]], k); push!(inc[vv[k]], k)
    end
    gens_at = [Int[] for _ in 1:nb]
    for g in 1:nG
        push!(gens_at[gbus[g]], g)
    end
    @constraint(model, [i=1:nb],
        sum(pG[g] for g in gens_at[i]; init=0.0) - demand[i] ==
        sum(flow[k] for k in outg[i]; init=0.0) -
        sum(flow[k] for k in inc[i];  init=0.0))
    @objective(model, Min,
        sum(base.c2[g] * pG[g]^2 + base.c1[g] * pG[g] + base.c0[g] for g in 1:nG))
    return (model=model, n_bus=nb, n_line=nL, n_gen=nG)
end

"""
    benchmark_dcopf_solve_times(c, Aval, scenario_indices; repeats=5, ...)

Solve time of the DC-OPF on the FULL network versus the natively-built REDUCED
network, under identical solver settings. Returns per-scenario minimum and median
`solve_time(model)` for each, plus the model dimensions being compared.
"""
function benchmark_dcopf_solve_times(c::MultiScenarioTxReductionCase, Aval,
                                     scenario_indices;
                                     repeats::Int=5, relax_pmin::Bool=true,
                                     threads::Int=1, time_limit=nothing,
                                     optimizer=nothing,
                                     progress_every::Int=25)
    base = c.base
    scenarios = Int.(collect(scenario_indices))
    H = length(scenarios)
    repeats >= 1 || error("repeats must be at least 1")

    if optimizer === nothing
        env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
        optimizer = () -> Gurobi.Optimizer(env)
    end
    A = round.(Int, Aval)
    red = extract_reduction(A)
    rep_full = collect(1:base.N)
    all_buses = collect(1:base.N)
    all_lines = collect(1:base.Ln)
    # A line whose ends land in one cluster carries no power in the reduced
    # network and would only appear as a zero-flow self-loop. Drop it.
    keep = [l for l in 1:base.Ln
            if red.rep_of[base.Efrom[l]] != red.rep_of[base.Eto[l]]]

    tsol_full = fill(NaN, H); tsol_red = fill(NaN, H)     # solve_time(model)
    wall_full = fill(NaN, H); wall_red = fill(NaN, H)     # batched cold solve
    work_full = fill(NaN, H); work_red = fill(NaN, H)     # Gurobi work units
    iter_full = fill(NaN, H); iter_red = fill(NaN, H)
    obj_full  = fill(NaN, H); obj_red  = fill(NaN, H)
    ok_full = falses(H); ok_red = falses(H)
    dims_full = nothing; dims_red = nothing

    # metrics:
    function timed(build)
        warm = build()
        optimize!(warm.model)                  # absorbs JIT; never timed
        models = [build() for _ in 1:repeats]  # construction excluded from timing
        t0 = time_ns()
        for mm in models
            optimize!(mm.model)
        end
        wall = (time_ns() - t0) / 1e9 / repeats
        st = [solve_time(mm.model) for mm in models]
        wk = [MOI.get(mm.model, Gurobi.ModelAttribute("Work")) for mm in models]
        it = [MOI.get(mm.model, Gurobi.ModelAttribute("IterCount")) for mm in models]
        good = termination_status(warm.model) in
               (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL)
        return (tsol=minimum(st), wall=wall, work=minimum(wk), iter=minimum(it),
                ok=good, obj=has_values(warm.model) ? objective_value(warm.model) : NaN,
                dims=(n_bus=warm.n_bus, n_line=warm.n_line, n_gen=warm.n_gen))
    end

    for (h, s) in enumerate(scenarios)
        f = timed(() -> _dcopf_on_partition(c, s, rep_full, all_buses, all_lines;
                relax_pmin=relax_pmin, threads=threads, time_limit=time_limit,
                optimizer=optimizer))
        r = timed(() -> _dcopf_on_partition(c, s, red.rep_of, red.retained, keep;
                relax_pmin=relax_pmin, threads=threads, time_limit=time_limit,
                optimizer=optimizer))
        tsol_full[h], wall_full[h], work_full[h], iter_full[h] = f.tsol, f.wall, f.work, f.iter
        tsol_red[h],  wall_red[h],  work_red[h],  iter_red[h]  = r.tsol, r.wall, r.work, r.iter
        ok_full[h], obj_full[h] = f.ok, f.obj
        ok_red[h],  obj_red[h]  = r.ok, r.obj
        dims_full = f.dims; dims_red = r.dims
        (progress_every > 0 && h % progress_every == 0) &&
            println("   benchmarked $h / $H scenarios")
    end

    both = ok_full .& ok_red
    safe(x) = max(x, eps())
    ratio(a, b) = sum(a[both]; init=0.0) / safe(sum(b[both]; init=0.0))
    return (
        scenario_indices=scenarios, scenario_ids=c.scenario_ids[scenarios],
        solve_time_full=tsol_full, solve_time_reduced=tsol_red,
        wall_full=wall_full, wall_reduced=wall_red,
        work_full=work_full, work_reduced=work_red,
        iters_full=iter_full, iters_reduced=iter_red,
        objective_full=obj_full, objective_reduced=obj_red,
        optimal_full=ok_full, optimal_reduced=ok_red, both_optimal=both,
        wall_speedup=wall_full ./ wall_red, work_speedup=work_full ./ work_red,
        total_solve_time_full=sum(tsol_full[both]; init=0.0),
        total_solve_time_reduced=sum(tsol_red[both]; init=0.0),
        total_wall_full=sum(wall_full[both]; init=0.0),
        total_wall_reduced=sum(wall_red[both]; init=0.0),
        total_work_full=sum(work_full[both]; init=0.0),
        total_work_reduced=sum(work_red[both]; init=0.0),
        # BOTH totals must be non-zero, or the ratio is a division by the timer's
        # resolution floor and prints a meaningless 1e12x speedup.
        solve_time_resolved=sum(tsol_full[both]; init=0.0) > 0 &&
                            sum(tsol_red[both]; init=0.0) > 0,
        speedup_solve_time=ratio(tsol_full, tsol_red),
        speedup_wall=ratio(wall_full, wall_red),
        speedup_work=ratio(work_full, work_red),
        speedup_iters=ratio(iter_full, iter_red),
        dims_full=dims_full, dims_reduced=dims_red,
        repeats=repeats, threads=threads, relax_pmin=relax_pmin,
    )
end

"""Print the full-vs-reduced solve-time benchmark."""
function report_dcopf_solve_times(bench)
    df, dr = bench.dims_full, bench.dims_reduced
    println("\n", repeat("=", 78))
    println("DC-OPF SOLVE-TIME BENCHMARK  (", bench.repeats,
            " repeats/model, Threads=", bench.threads, ", min of repeats)")
    println(repeat("=", 78))
    println("  model size    full    ", df.n_bus, " buses, ", df.n_line,
            " lines, ", df.n_gen, " gens")
    println("                reduced ", dr.n_bus, " buses, ", dr.n_line,
            " lines, ", dr.n_gen, " gens",
            "   (", round(100 * (1 - dr.n_bus / df.n_bus), digits=1), "% fewer buses, ",
            round(100 * (1 - dr.n_line / df.n_line), digits=1), "% fewer lines)")
    n = count(bench.both_optimal)
    println("  scenarios both solved to optimality = ", n, " / ",
            length(bench.scenario_indices))
    if n == 0
        println("  no comparable scenarios -- nothing to report")
        return
    end
    if bench.solve_time_resolved
        println("  solve_time(model)   full = ", round(bench.total_solve_time_full, digits=3),
                " s     reduced = ", round(bench.total_solve_time_reduced, digits=3),
                " s     -> ", round(bench.speedup_solve_time, digits=2), "x")
    else
        println("  solve_time(model)   UNRESOLVED -- Gurobi's Runtime timer returns")
        println("                      0.0 on models this small (a solve is well under")
        println("                      1 ms). Use the wall and work rows below.")
    end
    println("  wall clock / solve  full = ", round(1000 * bench.total_wall_full / n, digits=3),
            " ms    reduced = ", round(1000 * bench.total_wall_reduced / n, digits=3),
            " ms    -> ", round(bench.speedup_wall, digits=2), "x")
    println("  work units          full = ", round(bench.total_work_full, digits=4),
            "      reduced = ", round(bench.total_work_reduced, digits=4),
            "      -> ", round(bench.speedup_work, digits=2), "x   <-- deterministic")
    println("  simplex iterations  full = ", round(Int, sum(bench.iters_full[bench.both_optimal])),
            "       reduced = ", round(Int, sum(bench.iters_reduced[bench.both_optimal])),
            "       -> ", round(bench.speedup_iters, digits=2), "x")
    sp = filter(isfinite, bench.work_speedup[bench.both_optimal])
    if length(sp) > 1
        println("  per-scenario work speedup: min ", round(minimum(sp), digits=2),
                "x   q25 ", round(quantile(sp, 0.25), digits=2),
                "x   median ", round(median(sp), digits=2),
                "x   q75 ", round(quantile(sp, 0.75), digits=2),
                "x   max ", round(maximum(sp), digits=2), "x")
    end
    gap = abs.(bench.objective_reduced .- bench.objective_full) ./
          max.(abs.(bench.objective_full), 1e-9) .* 100
    gapv = filter(isfinite, gap[bench.both_optimal])
    println("  objective difference (context, not a timing result): max ",
            round(maximum(gapv; init=0.0), digits=4), "%")
end

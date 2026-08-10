# --------------------------------------------------------------------------- #
# THE OPTIMIZATION PROBLEM -- edge-based, assignment-matrix-free reduction
# (see transmission_network_reduction.pdf, the edge-based model).
#
# c_l (internal/external) is shared across scenarios; angles, external flows
# and internal transfers are per scenario. A single-scenario run is just
# S = 1, same model.
#
# No assignment matrix, no connectivity constraint: a solved internal-line set
# is connected by construction, and clusters are recovered afterwards as
# connected components of that subgraph (`assignment_from_line_status`).
#
# The objective maximizes INTERNAL LINES, not retained buses -- two
# clusterings with the same bus count can score differently, and the denser
# one wins even if less bus-efficient. Read n_retained as a byproduct, not
# the target.
#
# Only the runners `include`; assumes tnr_preprocessing.jl is already loaded.
# --------------------------------------------------------------------------- #

"""
Solve one common edge-based reduction over the selected scenarios.

`protection_indices` defaults to every loaded scenario, so a line that is
binding (or above `near_limit_threshold`) anywhere in the horizon remains
external even when only a subset supplies full physics constraints.
"""
function solve_reduction_edge_multiscenario(c::MultiScenarioTxReductionCase, epsL;
                                            scenario_indices=axes(c.p, 2),
                                            protection_indices=axes(c.p, 2),
                                            near_limit_threshold=nothing,
                                            time_limit=nothing,
                                            mipgap::Float64=1e-3,
                                            threads::Int=8,
                                            numeric_focus=nothing,
                                            cycle_cut_lens=(),
                                            warm_c=nothing,
                                            congestion_relaxation=0.0,
                                            congestion_relaxation_mode::Symbol=:none,
                                            internal_bound_scale::Real=3.0,
                                            merge_exact_blocks::Bool=false,
                                            merge_exact_mode::Symbol=:fix,
                                            merge_leaf_blocks::Bool=true,
                                            lmp_separation::Bool=false,
                                            lmp_threshold::Real=5.0,
                                            lmp_relax_pmin::Bool=true,
                                            lmp_opf_time_limit=nothing,
                                            int_feas_tol=nothing,
                                            feasibility_tol=nothing,
                                            optimality_tol=nothing,
                                            log_file=nothing)
    base = c.base
    win = multiscenario_windows(c, epsL;
        near_limit_threshold=near_limit_threshold,
        scenario_indices=scenario_indices,
        protection_indices=protection_indices,
        congestion_relaxation=congestion_relaxation,
        congestion_relaxation_mode=congestion_relaxation_mode)
    selected, protected = win.selected, win.protected
    philo, phiup = win.philo, win.phiup
    
    # Per-(line, scenario) internal-transfer bound: min(x*H, G), x =
    # internal_bound_scale. G alone is proven but one loose number for the whole
    # network; H is a per-line heuristic. See multiscenario_internal_bounds.
    G = multiscenario_internal_bounds(c, selected, philo, phiup, win.epsv;
                                      internal_bound_scale=internal_bound_scale)

    N, Ln, S = base.N, base.Ln, length(selected)
    u, v = base.Efrom, base.Eto

    println("\n Multi-scenario reduction: $S active / $(length(c.scenario_ids)) loaded scenarios")
    println(" Protected lines across the screening horizon = ", count(protected))
    # Which mechanism actually reaches how many lines: a merely near-limit line
    # already carries an eps-wide window, so only the exactly binding ones are
    # point-pinned and only they can be freed by the relaxation.
    n_exact_lines = count(vec(any(win.exact_binding, dims=2)))
    println(" Pinned (line, scenario) pairs = ", count(win.pinned),
            "   of which exactly binding = ", count(win.exact_binding),
            "   (", n_exact_lines, " distinct lines)")
    if congestion_relaxation_mode !== :none
        println(" Congestion relaxation = ", congestion_relaxation_mode,
                "   max delta = ", base.baseMVA * maximum(win.delta), " MW",
                congestion_relaxation_mode === :symmetric ?
                "   <-- may UNDER-state congestion; check the full-network overloads" : "")
    end

    m = Model(Gurobi.Optimizer)
    set_optimizer_attribute(m, "MIPGap", mipgap)
    set_optimizer_attribute(m, "Threads", threads)
    isnothing(int_feas_tol)    || set_optimizer_attribute(m, "IntFeasTol", int_feas_tol)
    isnothing(feasibility_tol) || set_optimizer_attribute(m, "FeasibilityTol", feasibility_tol)
    isnothing(optimality_tol)  || set_optimizer_attribute(m, "OptimalityTol", optimality_tol)
    isnothing(numeric_focus) || set_optimizer_attribute(m, "NumericFocus", numeric_focus)
    isnothing(time_limit) || set_optimizer_attribute(m, "TimeLimit", time_limit)
    isnothing(log_file) || set_optimizer_attribute(m, "LogFile", log_file)

    @variable(m, cl[1:Ln], Bin)
    @variable(m, el[1:Ln], Bin)
    @variable(m, vartheta[1:N, 1:S])
    # A :conservative or :symmetric relaxation deliberately pushes a congested
    # line's window past its rating. A hard +-frate bound on f would silently
    # clip that window back to the rating -- re-collapsing an exactly binding
    # line's window to the single point fhat and making the whole relaxation a
    # no-op. Bound f by what the window actually needs, never tighter than the
    # rating. (The reduced DC-OPF downstream still enforces the true rating.)
    fbound = [max(base.frate[l], abs(philo[l, s]), abs(phiup[l, s]))
              for l in 1:Ln, s in 1:S]
    @variable(m, -fbound[l, s] <= f[l=1:Ln, s=1:S] <= fbound[l, s])
    @variable(m, -G[l, s] <= gint[l=1:Ln, s=1:S] <= G[l, s])
    @constraint(m, [l=1:Ln], el[l] == 1 - cl[l])
    for l in findall(protected)
        @constraint(m, cl[l] == 0)
    end

    # Bridges and protected-free leaf blocks: provably mergeable at ZERO flow
    # error anywhere else (see exactly_mergeable_lines). 
    merge_lines = Int[]
    merge_info = nothing
    if merge_exact_blocks
        merge_exact_mode in (:fix, :warm) ||
            error("merge_exact_mode must be :fix or :warm, got $merge_exact_mode")
        merge_info = exactly_mergeable_lines(base, protected;
                                             include_leaf_blocks=merge_leaf_blocks)
        merge_lines = merge_info.lines
        # PRECONDITION. Shorting l makes g^int absorb the flow f carried, so the
        # transfer bound must actually admit it in every modeled scenario. A
        # tight internal_bound_scale could violate this; drop such a line rather
        # than hand the solver an infeasible fixing.
        keep = Int[]
        dropped = Int[]
        for l in merge_lines
            (all(abs(c.fhat[l, s]) <= G[l, ss] + 1e-9
                 for (ss, s) in enumerate(selected)) ? push!(keep, l) :
                                                      push!(dropped, l))
        end
        merge_lines = keep
        println(" Exactly-mergeable lines = ", length(merge_lines), " / ", Ln,
                "   (", length(merge_info.bridges), " bridges, ",
                merge_info.n_leaf_blocks, " leaf blocks)",
                "   mode = ", merge_exact_mode)
        isempty(dropped) || println("   dropped ", length(dropped),
                " line(s) whose |fhat| exceeds the internal-transfer bound G",
                " -- raise internal_bound_scale to use them: ", dropped)
        if merge_exact_mode === :fix
            for l in merge_lines
                @constraint(m, cl[l] == 1)
            end
        end
    end

    # LMP SEPARATION. Keep buses whose prices differ by more than lmp_threshold
    # out of a common cluster, by forbidding the shortest path between them from
    # going fully internal. 
    lmp_paths = Vector{Vector{Int}}()
    lmp_pairs = Tuple{Int,Int}[]
    lmp_mat = nothing
    if lmp_separation
        lmp_mat = full_network_lmps(c, selected;
                    relax_pmin=lmp_relax_pmin, time_limit=lmp_opf_time_limit)
        sep = lmp_separation_paths(base, lmp_mat;
                    lmp_threshold=lmp_threshold, protected=protected)
        lmp_paths, lmp_pairs = sep.paths, sep.pairs
        for p in lmp_paths
            @constraint(m, sum(cl[l] for l in p) <= length(p) - 1)
        end
        println(" LMP separation: threshold = ", lmp_threshold, " \$/MWh",
                "   max pair gap = ", round(sep.max_gap, digits=2),
                "   violating pairs = ", sep.n_violating_pairs,
                "   rows added = ", length(lmp_paths),
                isempty(lmp_paths) ? "" :
                "   (path length: min " * string(minimum(length, lmp_paths)) *
                ", median " * string(round(median(length.(lmp_paths)), digits=1)) *
                ", max " * string(maximum(length, lmp_paths)) * ")")
        n_forced = count(p -> length(p) == 1, lmp_paths)
        n_forced > 0 && println("   of which ", n_forced,
                " adjacent pair(s) -> that line is pinned external outright")
    end

    @constraint(m, [l=1:Ln, s=1:S],
        f[l, s] == base.Dx[l] * (vartheta[u[l], s] - vartheta[v[l], s]))
    @constraint(m, [s=1:S], vartheta[base.j0, s] == 0)

    outgoing = [findall(==(b), u) for b in 1:N]
    incoming = [findall(==(b), v) for b in 1:N]
    @constraint(m, [b=1:N, ss=1:S],
        sum(f[l, ss] + gint[l, ss] for l in outgoing[b]; init=0.0) -
        sum(f[l, ss] + gint[l, ss] for l in incoming[b]; init=0.0) ==
        c.p[b, selected[ss]])

    # SOS1 
    @variable(m, sup[1:Ln, 1:S] >= 0)
    @variable(m, slo[1:Ln, 1:S] >= 0)
    @constraint(m, [l=1:Ln, s=1:S], f[l, s] - phiup[l, s] <= sup[l, s])
    @constraint(m, [l=1:Ln, s=1:S], philo[l, s] - f[l, s] <= slo[l, s])
    for l in 1:Ln, s in 1:S
        @constraint(m, [cl[l], f[l, s]] in MOI.SOS1([1.0, 2.0]))
        @constraint(m, [el[l], gint[l, s]] in MOI.SOS1([1.0, 2.0]))
        @constraint(m, [el[l], sup[l, s]] in MOI.SOS1([1.0, 2.0]))
        @constraint(m, [el[l], slo[l, s]] in MOI.SOS1([1.0, 2.0]))
    end

    n_cycle_cuts = 0
    if !isempty(cycle_cut_lens)
        for cyc in short_cycles(base; lens=cycle_cut_lens)
            k = length(cyc)
            for l in cyc
                @constraint(m, cl[l] >= sum(cl[ll] for ll in cyc if ll != l) - (k - 2))
                n_cycle_cuts += 1
            end
        end
    end

    # All-external is feasible for every active scenario. Radial merges provide
    # a stronger common warm start whenever the line is unprotected everywhere.
    protected_lines = findall(protected)
    if !isnothing(warm_c)
        length(warm_c) == Ln || error("warm_c must have length $Ln")
        cl_warm = round.(Int, warm_c)
        for l in 1:Ln
            set_start_value(cl[l], cl_warm[l])
            set_start_value(el[l], 1 - cl_warm[l])
        end
    else
        cl_warm, rep_warm = _radial_warm_start(base, protected_lines, Int[])
        for l in merge_lines
            cl_warm[l] = 1
        end
        for l in 1:Ln
            set_start_value(cl[l], cl_warm[l])
            set_start_value(el[l], 1 - cl_warm[l])
        end
        # External lines keep f = fhat; gint is left for Gurobi to complete from nodal balance.
        A_warm = assignment_from_line_status(base, cl_warm)
        red_warm = extract_reduction(A_warm)
        Ew = incidence_matrix(base)
        Bred_w = A_warm * Ew * Diagonal(base.Dx) * Ew' * A_warm'
        ref_w = red_warm.rep_of[base.j0]
        free_w = setdiff(red_warm.retained, [ref_w])
        Fw = isempty(free_w) ? nothing : factorize(Bred_w[free_w, free_w])
        for ss in 1:S
            s = selected[ss]
            th = zeros(N)
            if !isnothing(Fw)
                th[free_w] = Fw \ (A_warm * c.p[:, s])[free_w]
            end
            for b in 1:N
                set_start_value(vartheta[b, ss], th[red_warm.rep_of[b]])
            end
            for l in 1:Ln
                fw = base.Dx[l] * (th[red_warm.rep_of[u[l]]] - th[red_warm.rep_of[v[l]]])
                set_start_value(f[l, ss], fw)
                set_start_value(sup[l, ss], max(0.0, fw - phiup[l, ss]))
                set_start_value(slo[l, ss], max(0.0, philo[l, ss] - fw))
            end
        end
    end

    @objective(m, Max, sum(cl))
    optimize!(m)
    st = termination_status(m)
    has_values(m) || error("Multi-scenario edge reduction failed: $st")

    c_raw = value.(cl)
    cval = round.(Int, c_raw)
    A = assignment_from_line_status(base, cval)
    red = extract_reduction(A)

    n_lmp_sep = 0
    n_lmp_merged = 0
    lmp_merged_pairs = Tuple{Int,Int}[]
    if lmp_separation
        for (i, j) in lmp_pairs
            if red.rep_of[i] == red.rep_of[j]
                n_lmp_merged += 1
                length(lmp_merged_pairs) < 20 && push!(lmp_merged_pairs, (i, j))
            else
                n_lmp_sep += 1
            end
        end
        println(" LMP separation achieved on ", n_lmp_sep, " / ",
                length(lmp_pairs), " constrained pairs",
                n_lmp_merged == 0 ? "  (all separated)" :
                "   -- " * string(n_lmp_merged) *
                " still merged through another path")
    end

    return (
        n_retained=red.n_retained, retained=red.retained, rep_of=red.rep_of,
        A=Float64.(A), c=cval, c_raw=c_raw,
        f=Array(value.(f)), gint=Array(value.(gint)),
        vartheta=Array(value.(vartheta)), G=G,
        philo=philo, phiup=phiup, protected=protected,
        pinned=win.pinned, exact_binding=win.exact_binding,
        delta=win.delta, congestion_relaxation_mode=congestion_relaxation_mode,
        scenario_indices=selected, scenario_ids=c.scenario_ids[selected],
        n_internal_lines=sum(cval), n_external_lines=Ln-sum(cval),
        n_line_binaries=Ln, n_cycle_cuts=n_cycle_cuts,
        merge_lines=merge_lines, n_merge_lines=length(merge_lines),
        merge_exact_mode=merge_exact_blocks ? merge_exact_mode : :off,
        n_merge_bridges=isnothing(merge_info) ? 0 : length(merge_info.bridges),
        n_merge_leaf_blocks=isnothing(merge_info) ? 0 : merge_info.n_leaf_blocks,
        status=st, solve_time=solve_time(m),
    )
end

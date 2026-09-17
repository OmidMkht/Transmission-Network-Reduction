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
# Only the runners `include`; assumes common/preprocessing.jl is already loaded.
# --------------------------------------------------------------------------- #

"""
Solve one common edge-based reduction over the selected scenarios.

`protection_indices` defaults to every loaded scenario, so a line that is
binding (or above `near_limit_threshold`) anywhere in the horizon remains
external even when only a subset supplies full physics constraints.
"""
# Which constraint family does a held set actually conflict with? Guessing that
# from a bare INFEASIBLE wastes runs, so ask the solver instead: the conflict
# (IIS) is a minimal set of rows that cannot all hold, and naming the families
# in it points straight at the option that differs between the two solves.
function _conflict_report(m, force_cons, families...)
    try
        compute_conflict!(m)
        get_attribute(m, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND ||
            return "\n  No conflict certificate available; rerun with gurobi_log=true."
        inconflict(c) =
            get_attribute(c, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
        nheld = count(inconflict, force_cons)
        hit = [(nm, count(inconflict, cons)) for (nm, cons) in families]
        named = [nm * " (" * string(k) * " rows)" for (nm, k) in hit if k > 0]
        msg = "\n  Conflict involves $nheld of the held lines"
        msg *= isempty(named) ?
            ", and no other tagged family -- the window or balance rows are the " *
            "binding side, which means the held set is not feasible at this eps." :
            " together with: " * join(named, ", ") * "."
        return msg * "\n  A held set only transfers between solves that see the SAME" *
               " constraints -- compare that option against the solve that produced it."
    catch err
        return "\n  (conflict computation unavailable: $err)"
    end
end

# --------------------------------------------------------------------------- #
# EPS LADDER support: keep one model open across rungs.
#
# An iteration changes exactly two things -- some binaries become fixed to 1,
# and eps grows. Everything else (protection, LMP separation paths, cycle cuts,
# the Kron chains, the graph, the injections) is identical at every rung, so
# rebuilding the model to express that is both wasteful and unsafe: two builds
# can disagree, and then a held set gets replayed into constraints it was never
# solved against.
#
# Both changes only ever RELAX or RESTRICT in a known direction:
#   windows  [fhat-eps, fhat+eps] widen, so every previous solution stays feasible
#   G        grows with eps (verified: every entry, every rung)
#   held     fixes binaries to 1, the only restriction
# so a rung's answer is always a feasible point of the next rung.
# --------------------------------------------------------------------------- #

"""
    retune_windows!(mdl, epsL_new; verbose=false)

Widen the flow windows to a new eps, in place. Only the window rows and the
internal-transfer bounds move; nothing else in the model is touched.

Errors if the protected set moves, which would mean eps had leaked into
something it must not -- protection comes from near_limit_threshold, and a
protected line appearing or vanishing mid-ladder would invalidate every held
decision taken under the old set.
"""
function retune_windows!(mdl, epsL_new; verbose::Bool=false)
    mdl.switch_form === :hull ||
        error("the eps ladder needs switch_form = :hull -- :sos1 carries the " *
              "window in SOS1 slacks, which have no row to retune")
    c, base, Ln, S = mdl.case, mdl.base, mdl.Ln, mdl.S
    win = multiscenario_windows(c, epsL_new;
        near_limit_threshold=mdl.near_limit_threshold,
        scenario_indices=mdl.selected,
        protection_indices=mdl.protection_indices,
        congestion_relaxation=mdl.congestion_relaxation,
        congestion_relaxation_mode=mdl.congestion_relaxation_mode)
    win.protected == mdl.protected || error(
        "the protected set changed when eps moved. Protection is set by " *
        "near_limit_threshold, not eps, so this is a bug -- every decision held " *
        "so far was taken under the old set and cannot be trusted.")

    philo, phiup = win.philo, win.phiup
    all(phiup .>= mdl.phiup[] .- 1e-12) && all(philo .<= mdl.philo[] .+ 1e-12) || error(
        "retune_windows! was asked to NARROW a window. The ladder only holds " *
        "merges while the windows widen -- a smaller eps invalidates every line " *
        "already held. Check that eps_ladder is non-decreasing.")
    G = multiscenario_internal_bounds(c, mdl.selected, philo, phiup, win.epsv;
            internal_bound_scale=mdl.internal_bound_scale)
    all(G .>= mdl.G[] .- 1e-9) || error(
        "the internal-transfer bound SHRANK when eps grew. The ladder relies on " *
        "it only relaxing; a held solution may no longer be feasible.")

    # Delete and re-add rather than patching coefficients: the window sits in
    # both the coefficient and the RHS of these rows, and re-adding is obviously
    # correct where a normalized-coefficient update is obviously fiddly. It is
    # 2*Ln*S rows, which is noise next to the solve.
    m, cl, f, gint = mdl.m, mdl.cl, mdl.f, mdl.gint
    for l in 1:Ln, s in 1:S
        delete(m, mdl.win_lo[l, s])
        delete(m, mdl.win_up[l, s])
        mdl.win_lo[l, s] = @constraint(m, philo[l, s] * (1 - cl[l]) <= f[l, s])
        mdl.win_up[l, s] = @constraint(m, f[l, s] <= phiup[l, s] * (1 - cl[l]))
        # gint keeps its own bound; only the ceiling moves, and only upward.
        set_lower_bound(gint[l, s], -G[l, s])
        set_upper_bound(gint[l, s], G[l, s])
    end
    mdl.philo[] = philo
    mdl.phiup[] = phiup
    mdl.G[] = G
    mdl.win[] = win
    verbose && println(" retuned windows: max |phi| = ",
                       round(maximum(abs, phiup), digits=5),
                       ", max G = ", round(maximum(G), digits=5))
    return mdl
end

"""
    hold_internal!(mdl, lines)

Fix `lines` internal for the rest of the ladder. Uses `fix` rather than a new
constraint so the same call can be made every rung without accumulating rows.
"""
function hold_internal!(mdl, lines)
    held = Int.(collect(lines))
    clash = [l for l in held if mdl.protected[l]]
    isempty(clash) || error(
        "asked to hold protected line(s) $clash internal. Protection is set by " *
        "near_limit_threshold, not eps, so this means the threshold or the " *
        "scenario set moved between rungs.")
    for l in held
        fix(mdl.cl[l], 1; force=true)
    end
    mdl.force_lines[] = held
    return mdl
end

function solve_reduction_edge_multiscenario(c::MultiScenarioTxReductionCase, epsL;
                                            scenario_indices=axes(c.p, 2),
                                            protection_indices=axes(c.p, 2),
                                            near_limit_threshold=nothing,
                                            time_limit=nothing,
                                            mipgap::Real=1e-3,
                                            threads::Int=8,
                                            numeric_focus=nothing,
                                            cycle_cut_lens=(),
                                            warm_c=nothing,
                                            congestion_relaxation=0.0,
                                            congestion_relaxation_mode::Symbol=:none,
                                            internal_bound_scale::Real=3.0,
                                            internal_rating_bound::Bool=false,
                                            switch_form::Symbol=:hull,
                                            force_internal=Int[],
                                            return_model::Bool=false,
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
                                            line_budget=nothing,
                                            hop_cap=nothing,
                                            size_cap=nothing,
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
    # (:sos1 only -- :hull bounds f by the window itself, which respects a
    # relaxed window exactly and makes the workaround unnecessary.)
    fbound = [max(base.frate[l], abs(philo[l, s]), abs(phiup[l, s]))
              for l in 1:Ln, s in 1:S]
    @variable(m, -fbound[l, s] <= f[l=1:Ln, s=1:S] <= fbound[l, s])
    @variable(m, -G[l, s] <= gint[l=1:Ln, s=1:S] <= G[l, s])
    @constraint(m, [l=1:Ln], el[l] == 1 - cl[l])
    # A collapsed line still EXISTS and still has a rating, so the transfer it
    # carries cannot exceed it. gint is a variable, so this asks only that SOME
    # rating-feasible routing exists inside the cluster -- which is exactly the
    # condition for the merge to be physically deliverable. It also subsumes the
    # [el, gint] SOS1: cl = 0 forces gint = 0.
    if internal_rating_bound
        @constraint(m, [l=1:Ln, s=1:S],
            gint[l, s] <= min(base.frate[l], G[l, s]) * cl[l])
        @constraint(m, [l=1:Ln, s=1:S],
            -min(base.frate[l], G[l, s]) * cl[l] <= gint[l, s])
    end
    # Tagged so an infeasible solve can say WHICH family it conflicts with
    # rather than leaving the caller to guess -- see _conflict_report.
    prot_cons = [@constraint(m, cl[l] == 0) for l in findall(protected)]
    caps = Main.Caps.add_caps!(m, cl, base; budget=line_budget, hop_cap, size_cap,
                               protected=protected)

    # Lines an earlier pass already collapsed and this one may not reopen. A
    # solution feasible at some eps stays feasible at any larger eps (windows
    # only widen), so a ladder that grows eps can hold its earlier merges and
    # still be feasible -- see iterative/.
    force_lines = Int.(collect(force_internal))
    force_cons = ConstraintRef[]
    if !isempty(force_lines)
        clash = [l for l in force_lines if protected[l]]
        isempty(clash) || error(
            "force_internal asks to collapse protected line(s) $clash. Protection " *
            "is set by near_limit_threshold, not eps, so this means the threshold " *
            "or the scenario set moved between passes.")
        append!(force_cons, [@constraint(m, cl[l] == 1) for l in force_lines])
        println(" Held internal from earlier passes = ", length(force_lines), " / ", Ln)
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
    lmp_cons = ConstraintRef[]
    lmp_pairs = Tuple{Int,Int}[]
    lmp_mat = nothing
    if lmp_separation
        lmp_mat = full_network_lmps(c, selected;
                    relax_pmin=lmp_relax_pmin, time_limit=lmp_opf_time_limit)
        sep = lmp_separation_paths(base, lmp_mat;
                    lmp_threshold=lmp_threshold, protected=protected)
        lmp_paths, lmp_pairs = sep.paths, sep.pairs
        for p in lmp_paths
            push!(lmp_cons, @constraint(m, sum(cl[l] for l in p) <= length(p) - 1))
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

    win_lo, win_up = nothing, nothing
    # THE SWITCH. Collapsing a line is zeroing its impedance, so with
    # x_l = 1/D_l the physics reads (1-c_l) x_l f_l = theta_u - theta_v: the
    # disjunction leaves the flow equation and lands on the window, where
    #
    #     (1-c_l) philo <= f <= (1-c_l) phiup
    #
    # does merge-switch and thermal window in ONE pair. c=0 gives the window;
    # c=1 gives f=0, hence theta_u = theta_v, which is the merge. It is also
    # the CONVEX HULL of those two cases, so no tighter per-line form exists.
    #
    # :sos1 is the older equivalent, kept to reproduce earlier runs. Gurobi's
    # presolve turns it into the rows below anyway, so the bound is identical;
    # what :hull saves is 4 SOS families and the sup/slo slacks -- measured on
    # case118 at 5677 -> 5172 nodes, and 2684 with internal_rating_bound.
    if switch_form === :hull
        win_lo = @constraint(m, [l=1:Ln, s=1:S], philo[l, s] * (1 - cl[l]) <= f[l, s])
        win_up = @constraint(m, [l=1:Ln, s=1:S], f[l, s] <= phiup[l, s] * (1 - cl[l]))
        # gint = 0 when the line stays external. internal_rating_bound adds the
        # tighter rating version of exactly this pair.
        internal_rating_bound || begin
            @constraint(m, [l=1:Ln, s=1:S], -G[l, s] * cl[l] <= gint[l, s])
            @constraint(m, [l=1:Ln, s=1:S], gint[l, s] <= G[l, s] * cl[l])
        end
    elseif switch_form === :sos1
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
    else
        error("switch_form must be :hull or :sos1, got $switch_form")
    end

    n_cycle_cuts = 0
    cycle_cons = ConstraintRef[]
    cyc_list = Vector{Vector{Int}}()
    if !isempty(cycle_cut_lens)
        for cyc in short_cycles(base; lens=cycle_cut_lens)
            k = length(cyc)
            push!(cyc_list, collect(cyc))
            for l in cyc
                push!(cycle_cons,
                    @constraint(m, cl[l] >= sum(cl[ll] for ll in cyc if ll != l) - (k - 2)))
                n_cycle_cuts += 1
            end
        end
    end

    # All-external is feasible for every active scenario. Radial merges provide
    # a stronger common warm start whenever the line is unprotected everywhere.
    protected_lines = findall(protected)
    if !isnothing(warm_c)
        length(warm_c) == Ln || error("warm_c must have length $Ln")
        cl_warm = Int.(Main.Caps.trim_to_caps(base, round.(Int, warm_c);
                      budget=line_budget, hop_cap, size_cap))
        for l in 1:Ln
            set_start_value(cl[l], cl_warm[l])
            set_start_value(el[l], 1 - cl_warm[l])
        end
    else
        cl_warm, rep_warm = _radial_warm_start(base, protected_lines, Int[])
        for l in merge_lines
            cl_warm[l] = 1
        end
        cl_warm = Int.(Main.Caps.trim_to_caps(base, cl_warm; budget=line_budget,
                                             hop_cap, size_cap))
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
                if switch_form === :sos1
                    set_start_value(sup[l, ss], max(0.0, fw - phiup[l, ss]))
                    set_start_value(slo[l, ss], max(0.0, philo[l, ss] - fw))
                end
            end
        end
    end

    @objective(m, Max, sum(cl))

    # Everything the solve-and-extract step needs, in one bundle. The eps ladder
    # holds this open across rungs and mutates only what eps actually changes
    # (retune_windows!) instead of building a second model that can drift away
    # from this one -- which is how a held set ends up being replayed into
    # constraints it was never solved against.
    mdl = (m=m, cl=cl, el=el, f=f, gint=gint, vartheta=vartheta,
           case=c, base=base, selected=selected, Ln=Ln, S=S,
           switch_form=switch_form, win_lo=win_lo, win_up=win_up,
           philo=Ref(philo), phiup=Ref(phiup), G=Ref(G), win=Ref(win),
           protected=protected,
           force_lines=Ref(force_lines), force_cons=Ref(force_cons),
           prot_cons=prot_cons, lmp_cons=lmp_cons, cycle_cons=cycle_cons,
           lmp_paths=lmp_paths, lmp_pairs=lmp_pairs, lmp_separation=lmp_separation,
           cyc_list=cyc_list, n_cycle_cuts=n_cycle_cuts,
           merge_lines=merge_lines, merge_info=merge_info,
           merge_exact_blocks=merge_exact_blocks, merge_exact_mode=merge_exact_mode,
           congestion_relaxation_mode=congestion_relaxation_mode,
           congestion_relaxation=congestion_relaxation,
           internal_bound_scale=internal_bound_scale,
           internal_rating_bound=internal_rating_bound,
           near_limit_threshold=near_limit_threshold,
           protection_indices=protection_indices,
           caps=(line_budget=line_budget, hop_cap=hop_cap, size_cap=size_cap,
                 hop_lazy=caps.hop_lazy, n_hop_rows=caps.n_hop_rows))
    return solve_reduction_model!(mdl; return_model=return_model)
end

"""
    solve_reduction_model!(mdl; return_model=false)

Optimize the model in `mdl` and extract the answer. Split out of
`solve_reduction_edge_multiscenario` so the eps ladder can re-solve the SAME
model after retuning it, instead of building a second one.
"""
function solve_reduction_model!(mdl; return_model::Bool=false)
    m, cl, base, c = mdl.m, mdl.cl, mdl.base, mdl.case
    f, gint, vartheta = mdl.f, mdl.gint, mdl.vartheta
    Ln, selected = mdl.Ln, mdl.selected
    philo, phiup, G, win = mdl.philo[], mdl.phiup[], mdl.G[], mdl.win[]
    protected, lmp_paths, cyc_list = mdl.protected, mdl.lmp_paths, mdl.cyc_list
    lmp_pairs, lmp_separation = mdl.lmp_pairs, mdl.lmp_separation
    merge_lines, merge_info = mdl.merge_lines, mdl.merge_info
    merge_exact_blocks, merge_exact_mode = mdl.merge_exact_blocks, mdl.merge_exact_mode
    congestion_relaxation_mode = mdl.congestion_relaxation_mode
    n_cycle_cuts = mdl.n_cycle_cuts
    force_lines, force_cons = mdl.force_lines[], mdl.force_cons[]
    prot_cons, lmp_cons, cycle_cons = mdl.prot_cons, mdl.lmp_cons, mdl.cycle_cons

    optimize!(m)
    st = termination_status(m)

    if !has_values(m) && !isempty(force_lines)
        # Nothing came back AND lines were held: is the held set really in
        # conflict, or did presolve just decline to decide? INFEASIBLE_OR_
        # UNBOUNDED is a presolve verdict and this model is bounded by
        # construction (max sum(cl) <= Ln), so ask again with DualReductions
        # off. This is a DIAGNOSIS, not a second attempt at the problem: cap it,
        # and put the solver back as it was, because the ladder reuses this
        # model for the next rung.
        tl0 = get_optimizer_attribute(m, "TimeLimit")
        set_optimizer_attribute(m, "DualReductions", 0)
        set_optimizer_attribute(m, "TimeLimit", 120.0)
        optimize!(m)
        st = termination_status(m)
        set_optimizer_attribute(m, "DualReductions", 1)
        set_optimizer_attribute(m, "TimeLimit", tl0)
    end

    # TIME_LIMIT with an incumbent is a perfectly usable answer -- the ladder is
    # built on incumbents, and a rung that ran out of clock still hands a valid
    # clustering upward. Only a solve with NOTHING to return has failed.
    if !has_values(m)
        diag = isempty(force_lines) ? "" :
            "\n  $(length(force_lines)) line(s) were held internal by an earlier pass." *
            (st == MOI.INFEASIBLE ?
                _conflict_report(m, force_cons, ("protection", prot_cons),
                                 ("LMP separation", lmp_cons), ("cycle cuts", cycle_cons)) :
                "\n  The solve ended $st with no incumbent at all, so this is a " *
                "budget problem rather than a conflict -- raise solve_time_limit, " *
                "or make the step to this rung smaller.")
        error("Multi-scenario edge reduction failed: $st" * diag)
    end

    c_raw = value.(cl)
    cval = round.(Int, c_raw)

    # cval is what callers replay -- the eps ladder re-imposes it as cl == 1
    # equalities on the next rung. Gurobi only promises each binary is within
    # IntFeasTol of an integer, so a row satisfied at the margin by c_raw can be
    # VIOLATED by cval. Downstream that shows up as an inexplicable INFEASIBLE
    # on a model whose cuts never changed, so catch it at the source instead.
    let
        bad_lmp = count(p -> sum(cval[l] for l in p) > length(p) - 1, lmp_paths)
        bad_cyc = 0
        for cyc in cyc_list, l in cyc
            cval[l] < sum(cval[ll] for ll in cyc if ll != l) - (length(cyc) - 2) &&
                (bad_cyc += 1)
        end
        worst = maximum(abs.(c_raw .- cval); init=0.0)
        if bad_lmp + bad_cyc > 0
            @warn "The ROUNDED solution violates rows the solved one satisfied" *
                  " -- callers that replay it will see INFEASIBLE" bad_lmp bad_cyc worst
        end
    end
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
        status=st, solve_time=solve_time(m), caps=mdl.caps,
        model=return_model ? mdl : nothing,
    )
end

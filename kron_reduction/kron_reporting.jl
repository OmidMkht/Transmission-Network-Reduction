# --------------------------------------------------------------------------- #
# KRON-AWARE REPORTING -- mirrors tnr_reporting.jl's sweep_multiscenario /
# write_multiscenario_outputs closely enough that TxReport.print_relaxation_table,
# TxReport.write_relaxation_table, TxReport.report_reduction,
# TxReport.report_dcopf_failures and TxReport.plot_network all work UNMODIFIED
# against this module's output. tnr_reporting.jl itself is never touched --
# sweep_multiscenario has no seam to redirect its one hardcoded solve call at a
# smaller case, so this is a small, deliberate duplication of orchestration
# only; every validation/benchmark CALL inside it is the existing function,
# unmodified, run against the TRUE full case.
# --------------------------------------------------------------------------- #
module KronReport

using Printf, DelimitedFiles, Dates, Statistics

const TxReport = Main.TxReport   # already defined by the time this file loads

"""
    sweep_kron_multiscenario(TR, c_full, c_boundary, epsL_full, epsL_boundary,
                             kron_map, selected, month_indices, cfg) -> (rows, artifacts)

`c_boundary`/`epsL_boundary`/`kron_map` are computed ONCE by the caller before
this loop -- Kron reduction depends only on `near_limit_threshold` and the
base fhat/frate data, never on the relaxation (mode, delta), so re-deriving it
per setting would be pure waste. Only the boundary solve re-runs per setting.
"""
function sweep_kron_multiscenario(TR, c_full, c_boundary, epsL_full, epsL_boundary,
                                  kron_map, selected, month_indices, cfg)
    rows = NamedTuple[]
    artifacts = NamedTuple[]
    mva = c_full.base.baseMVA
    dc = cfg.dcopf

    for (k, (mode, delta)) in enumerate(cfg.relaxation_sweep)
        println()
        println(repeat("=", 78))
        println("KRON SETTING $k/$(length(cfg.relaxation_sweep)):  ",
                TxReport.relaxation_pretty(mode, delta))
        println(repeat("=", 78))

        log_file = nothing
        if get(cfg, :gurobi_log, false)
            log_dir = joinpath(cfg.output_dir, TxReport.relaxation_label(mode, delta))
            mkpath(log_dir)
            log_file = joinpath(log_dir, "gurobi.log")
            isfile(log_file) && rm(log_file)
            println("Gurobi log -> ", log_file)
        end

        gen = nothing
        if cfg.scenario_generation
            gen = TR.solve_reduction_with_scenario_generation(
                c_boundary, epsL_boundary, selected, month_indices;
                near_limit_threshold=cfg.near_limit_threshold,
                max_iterations=cfg.max_sc_generation_iterations,
                tolerance=cfg.screening_tolerance,
                time_limit=cfg.solve_time_limit,
                numeric_focus=cfg.numeric_focus,
                cycle_cut_lens=cfg.cycle_cut_lens,
                congestion_relaxation=delta,
                congestion_relaxation_mode=mode,
                internal_bound_scale=get(cfg, :internal_bound_scale, 3.0),
                merge_exact_blocks=get(cfg, :merge_exact_blocks, false),
                merge_exact_mode=get(cfg, :merge_exact_mode, :fix),
                merge_leaf_blocks=get(cfg, :merge_leaf_blocks, true),
                lmp_separation=get(cfg, :lmp_separation, false),
                lmp_threshold=get(cfg, :lmp_threshold, 5.0),
                lmp_relax_pmin=get(cfg, :lmp_relax_pmin, true),
                lmp_opf_time_limit=get(cfg, :opf_time_limit, nothing),
                int_feas_tol=get(cfg, :int_feas_tol, 1e-7),
                feasibility_tol=get(cfg, :feasibility_tol, 1e-7),
                optimality_tol=get(cfg, :optimality_tol, 1e-7),
                log_file=log_file)
            r_boundary = gen.r
            active = gen.active
        else
            r_boundary = TR.solve_reduction_edge_multiscenario(
                c_boundary, epsL_boundary;
                scenario_indices=selected,
                protection_indices=month_indices,
                near_limit_threshold=cfg.near_limit_threshold,
                time_limit=cfg.solve_time_limit,
                numeric_focus=cfg.numeric_focus,
                cycle_cut_lens=cfg.cycle_cut_lens,
                congestion_relaxation=delta,
                congestion_relaxation_mode=mode,
                internal_bound_scale=get(cfg, :internal_bound_scale, 3.0),
                merge_exact_blocks=get(cfg, :merge_exact_blocks, false),
                merge_exact_mode=get(cfg, :merge_exact_mode, :fix),
                merge_leaf_blocks=get(cfg, :merge_leaf_blocks, true),
                lmp_separation=get(cfg, :lmp_separation, false),
                lmp_threshold=get(cfg, :lmp_threshold, 5.0),
                lmp_relax_pmin=get(cfg, :lmp_relax_pmin, true),
                lmp_opf_time_limit=get(cfg, :opf_time_limit, nothing),
                int_feas_tol=get(cfg, :int_feas_tol, 1e-7),
                feasibility_tol=get(cfg, :feasibility_tol, 1e-7),
                optimality_tol=get(cfg, :optimality_tol, 1e-7),
                log_file=log_file)
            active = selected
        end

        A_full = TR.unfold_kron_assignment(c_full.base, kron_map, r_boundary.A)
        # Per-chain internal/external at the BOUNDARY level -- independent of
        # collapse_external_chains, which only changes what A_display folds,
        # never what the boundary MILP itself decided.
        red_b = TR.extract_reduction(round.(Int, r_boundary.A))
        chain_internal = [red_b.rep_of[kron_map.boundary_of[ch.x]] ==
                           red_b.rep_of[kron_map.boundary_of[ch.y]]
                           for ch in kron_map.chains]
        # collapse_external_chains (default true): chain-interior buses never
        # resurface in the REPORTED network, even when their equivalent line
        # stays external -- see kron_unfold.jl's kron_display_assignment.
        # bench/val/timing below always use the fully-reinserted A_full
        # regardless of this flag; only r_display/the plot/the CSVs change.
        A_display = get(cfg.kron_reduction, :collapse_external_chains, true) ?
            TR.kron_display_assignment(c_full.base, kron_map, r_boundary.A) : A_full
        r_display = TR.make_r_display(c_full, r_boundary, kron_map, A_display;
            near_limit_threshold=cfg.near_limit_threshold, protection_indices=month_indices)

        dcopf_indices = dc.scope === :month ? month_indices :
                        dc.scope === :selected ? active : Int[]
        # Boundary-solve feasibility: is the SOLVER'S OWN boundary-network point
        # feasible? A different, narrower question than everything below, which
        # all runs against the TRUE full case -- see make_r_display's contract.
        chk = TR.model_feasibility_check_multiscenario(c_boundary, r_boundary; tol=1e-6)
        congestion_threshold = isnothing(cfg.near_limit_threshold) ? 0.9999 :
                               Float64(cfg.near_limit_threshold)
        bench = TR.benchmark_reduction_scenarios(
            c_full, A_full, epsL_full, month_indices;
            congestion_threshold=congestion_threshold,
            near_limit_threshold=cfg.near_limit_threshold,
            protection_indices=month_indices,
            congestion_relaxation=delta,
            congestion_relaxation_mode=mode,
            tolerance=cfg.screening_tolerance)
        val = isempty(dcopf_indices) ? nothing : TR.validate_reduced_dcopf_scenarios(
            c_full, A_full, dcopf_indices;
            relax_pmin=dc.relax_pmin,
            time_limit=dc.time_limit,
            relative_tolerance=dc.relative_tolerance,
            objective_tolerance_pct=dc.objective_tolerance_pct,
            lmp_tolerance=dc.lmp_tolerance,
            measure_repair=dc.measure_repair,
            training_indices=active)
        timing = get(cfg, :solve_time_benchmark, true) && !isempty(dcopf_indices) ?
            TR.benchmark_dcopf_solve_times(
                c_full, A_full, dcopf_indices;
                repeats=get(cfg, :benchmark_repeats, 5),
                relax_pmin=dc.relax_pmin,
                threads=get(cfg, :benchmark_threads, 1),
                time_limit=dc.time_limit,
                progress_every=get(cfg, :benchmark_progress_every, 200)) : nothing

        push!(rows, (
            mode=mode, delta=delta,
            label=TxReport.relaxation_label(mode, delta),
            n_retained=r_display.n_retained,
            reduction_pct=100 * (c_full.base.N - r_display.n_retained) / c_full.base.N,
            n_internal=r_display.n_internal_lines,
            solve_time=r_boundary.solve_time,
            status=string(r_boundary.status),
            genuine=chk.genuine,
            window_feasible_pct=100 * bench.window_feasible_fraction,
            worst_window_mw=mva * bench.worst_window_violation,
            worst_external_error=bench.max_external_normalized_error,
            n_dcopf_feasible=isnothing(val) ? -1 : val.n_dispatch_feasible,
            n_dcopf_total=isnothing(val) ? 0 : length(val.scenario_indices),
            worst_utilization=isnothing(val) ? NaN : val.worst_max_utilization,
            worst_objective_pct=isnothing(val) ? NaN : val.worst_abs_objective_change_pct,
            worst_repair_mw=isnothing(val) ? NaN : mva * val.worst_repair_redispatch,
            mean_repair_cost_pct=isnothing(val) ? NaN : val.mean_repair_cost_pct,
            n_ext_lines_over_rating=bench.n_lines_over_rating,
            n_scenarios_over_rating=bench.n_scenarios_over_rating,
            n_ext_lines_over_adjusted=bench.n_lines_over_adjusted,
            n_scenarios_over_adjusted=bench.n_scenarios_over_adjusted,
            n_fail_in_training=isnothing(val) ? -1 :
                count((val.original_optimal .& val.reduced_optimal .&
                       .!val.dispatch_feasible) .& val.in_training),
            overload_internal_pairs=isnothing(val) ? -1 :
                sum(val.n_violating_internal),
            overload_external_pairs=isnothing(val) ? -1 :
                sum(val.n_violating_external),
            n_hours_spurious_binding=isnothing(val) ? -1 :
                val.n_hours_with_spurious_binding,
            n_spurious_lines=isnothing(val) ? -1 : val.n_spurious_lines,
            n_active=length(active),
            cg_iterations=isnothing(gen) ? 0 : gen.iterations,
            cg_converged=isnothing(gen) ? false : gen.converged,
            bench_buses_full=isnothing(timing) ? -1 : timing.dims_full.n_bus,
            bench_buses_reduced=isnothing(timing) ? -1 : timing.dims_reduced.n_bus,
            bench_lines_full=isnothing(timing) ? -1 : timing.dims_full.n_line,
            bench_lines_reduced=isnothing(timing) ? -1 : timing.dims_reduced.n_line,
            speedup_work=isnothing(timing) ? NaN : timing.speedup_work,
            speedup_wall=isnothing(timing) ? NaN : timing.speedup_wall,
            speedup_iters=isnothing(timing) ? NaN : timing.speedup_iters,
            speedup_solve_time=isnothing(timing) || !timing.solve_time_resolved ?
                NaN : timing.speedup_solve_time,
            wall_ms_full=isnothing(timing) ? NaN :
                1000 * timing.total_wall_full / max(count(timing.both_optimal), 1),
            wall_ms_reduced=isnothing(timing) ? NaN :
                1000 * timing.total_wall_reduced / max(count(timing.both_optimal), 1),
            # --- Kron-specific ---
            n_chains=length(kron_map.chains),
            n_kron_eliminated=kron_map.full_N - kron_map.boundary_N,
            boundary_n_bus=kron_map.boundary_N,
            boundary_n_line=c_boundary.base.Ln,
        ))
        push!(artifacts, (r=r_display, r_boundary=r_boundary, A=A_display, A_full=A_full,
                          chain_internal=chain_internal,
                          chk=chk, bench=bench, val=val, gen=gen,
                          timing=timing, active=active, kron_map=kron_map,
                          mode=mode, delta=delta,
                          label=TxReport.relaxation_label(mode, delta)))
    end
    return rows, artifacts
end

# --------------------------------------------------------------------------- #
# Kron-aware plotting: three panels (full network with Kron-chains highlighted,
# the Kron-reduced "boundary" network, and the final MILP-reduced network).
# Loaded lazily, same reason as TxReport.ensure_plots -- CairoMakie must never
# load on a plots-off run.
# --------------------------------------------------------------------------- #
const _KRON_PLOTS_LOADED = Ref(false)

function ensure_kron_plots(plots_file::AbstractString, kron_plots_file::AbstractString)
    TxReport.ensure_plots(plots_file)
    _KRON_PLOTS_LOADED[] && return nothing
    isfile(kron_plots_file) || error("Kron plot helper not found: $kron_plots_file")
    Base.include(@__MODULE__, kron_plots_file)
    _KRON_PLOTS_LOADED[] = true
    return nothing
end

"""
Render the three-panel Kron figure. `plots_file`/`kron_plots_file` are
included on first use.
"""
function plot_kron_network(plots_file, kron_plots_file, c_full, kron_map, c_boundary, A_full;
                           path, open_after::Bool=false, kwargs...)
    ensure_kron_plots(plots_file, kron_plots_file)
    f = getfield(@__MODULE__, :plot_kron_reduction)
    Base.invokelatest(f, c_full, kron_map, c_boundary, A_full; path=path, kwargs...)
    open_after && TxReport.open_externally(path)
    return path
end

"Standalone printout of chain counts/sizes -- run once, before the sweep."
function report_kron_chains(c_full, kron_map)
    TxReport.section("Kron reduction")
    println("  full network        = ", kron_map.full_N, " buses, ",
            c_full.base.Ln, " lines")
    println("  boundary network    = ", kron_map.boundary_N, " buses  ",
            "(", kron_map.full_N - kron_map.boundary_N, " eliminated)")
    println("  chains kept         = ", length(kron_map.chains))
    println("  rings (untouched)   = ", length(kron_map.rings))
    println("  self-loop chains (untouched) = ", length(kron_map.self_loop_chains))
    isempty(kron_map.chains) && return nothing
    lens = [length(ch.interior) for ch in kron_map.chains]
    println("  chain length: min=", minimum(lens), "  median=", sort(lens)[cld(end,2)],
            "  max=", maximum(lens))
    return nothing
end

# --------------------------------------------------------------------------- #
# Per-setting file output. Mirrors write_multiscenario_outputs, reading
# internal/protected/rep_of from art.r (= r_display, full-network-correct).
# selected_internal_transfer_mw.csv is DROPPED (r_boundary.gint has no
# full-network analogue for former chain-interior lines -- see kron_unfold.jl)
# and replaced with two Kron-specific files.
# --------------------------------------------------------------------------- #
function write_kron_multiscenario_outputs(dir, c, art, kron_map,
                                          selection, selected, month_indices, cfg)
    mkpath(dir)
    mva = c.base.baseMVA
    r, A, bench, val = art.r, art.A, art.bench, art.val
    has_calendar = haskey(cfg, :month) && haskey(cfg, :year)
    first_hour = DateTime(get(cfg, :year, 2017), 1, 1)

    writedlm(joinpath(dir, "assignment_matrix.csv"), A, ',')
    writedlm(joinpath(dir, "monthly_reduced_flow_mw.csv"), mva .* bench.screen.flow, ',')

    timing = get(art, :timing, nothing)
    if !isnothing(timing)
        open(joinpath(dir, "solve_time_benchmark.csv"), "w") do io
            println(io, "scenario_id,timestamp,optimal_full,optimal_reduced,",
                        "solve_time_full_s,solve_time_reduced_s,",
                        "wall_full_ms,wall_reduced_ms,work_full,work_reduced,",
                        "iters_full,iters_reduced,work_speedup,wall_speedup,",
                        "objective_full,objective_reduced")
            for (h, s) in enumerate(timing.scenario_indices)
                ts = has_calendar ? first_hour + Hour(c.scenario_ids[s] - 1) : ""
                println(io, c.scenario_ids[s], ",", ts, ",",
                        Int(timing.optimal_full[h]), ",", Int(timing.optimal_reduced[h]), ",",
                        timing.solve_time_full[h], ",", timing.solve_time_reduced[h], ",",
                        1000 * timing.wall_full[h], ",", 1000 * timing.wall_reduced[h], ",",
                        timing.work_full[h], ",", timing.work_reduced[h], ",",
                        timing.iters_full[h], ",", timing.iters_reduced[h], ",",
                        timing.work_speedup[h], ",", timing.wall_speedup[h], ",",
                        timing.objective_full[h], ",", timing.objective_reduced[h])
            end
        end
        open(joinpath(dir, "solve_time_summary.csv"), "w") do io
            println(io, "metric,full,reduced,speedup")
            println(io, "buses,", timing.dims_full.n_bus, ",", timing.dims_reduced.n_bus,
                    ",", timing.dims_full.n_bus / max(timing.dims_reduced.n_bus, 1))
            println(io, "lines,", timing.dims_full.n_line, ",", timing.dims_reduced.n_line,
                    ",", timing.dims_full.n_line / max(timing.dims_reduced.n_line, 1))
            println(io, "solve_time_s,", timing.total_solve_time_full, ",",
                    timing.total_solve_time_reduced, ",",
                    timing.solve_time_resolved ? timing.speedup_solve_time : NaN)
            println(io, "wall_s,", timing.total_wall_full, ",", timing.total_wall_reduced,
                    ",", timing.speedup_wall)
            println(io, "work_units,", timing.total_work_full, ",", timing.total_work_reduced,
                    ",", timing.speedup_work)
            println(io, "simplex_iterations,", sum(timing.iters_full), ",",
                    sum(timing.iters_reduced), ",", timing.speedup_iters)
        end
    end

    open(joinpath(dir, "bus_mapping.csv"), "w") do io
        println(io, "bus_id,representative_bus_id,retained")
        for j in 1:c.base.N
            rep = r.rep_of[j]
            println(io, "$(c.bus_ids[j]),$(c.bus_ids[rep]),$(j == rep ? 1 : 0)")
        end
    end

    congested_set = Set(selection.congested_lines)
    protected_set = Set(findall(r.protected))
    open(joinpath(dir, "line_status.csv"), "w") do io
        println(io, "line,from_bus_id,to_bus_id,internal,protected,congested_during_month")
        for l in 1:c.base.Ln
            println(io, "$l,$(c.bus_ids[c.base.Efrom[l]]),$(c.bus_ids[c.base.Eto[l]]),",
                    "$(r.c[l]),$(l in protected_set ? 1 : 0),",
                    "$(l in congested_set ? 1 : 0)")
        end
    end

    open(joinpath(dir, "monthly_benchmark.csv"), "w") do io
        println(io, "scenario_id,timestamp,selected,congested_line_count,window_feasible,",
                "max_window_violation_mw,max_external_normalized_error,",
                "max_all_line_normalized_error")
        selected_set = Set(selected)
        for (h, s) in enumerate(month_indices)
            println(io, "$(c.scenario_ids[s]),$(first_hour + Hour(c.scenario_ids[s] - 1)),",
                    "$(s in selected_set ? 1 : 0),",
                    "$(selection.congested_count_by_hour[h]),",
                    "$(bench.window_feasible_by_scenario[h] ? 1 : 0),",
                    "$(mva * bench.screen.scenario_violation[h]),",
                    "$(bench.max_external_error_by_scenario[h]),",
                    "$(bench.max_all_line_error_by_scenario[h])")
        end
    end

    open(joinpath(dir, "monthly_line_benchmark.csv"), "w") do io
        println(io, "line,max_absolute_flow_error_mw,max_normalized_flow_error,",
                "max_relative_window_violation")
        for l in 1:c.base.Ln
            println(io, "$l,$(mva * maximum(bench.absolute_error[l, :])),",
                    "$(maximum(bench.normalized_error[l, :])),",
                    "$(maximum(bench.relative_window_violation[l, :]))")
        end
    end

    if !isnothing(val)
        open(joinpath(dir, "monthly_dcopf_validation.csv"), "w") do io
            println(io, "scenario_id,original_status,reduced_status,",
                    "original_optimal,reduced_optimal,dispatch_feasible_on_original,",
                    "original_objective,reduced_objective,objective_change_pct,",
                    "objective_within_tolerance,max_overload_mw,max_relative_overload,",
                    "max_full_network_utilization,total_overload_mw,n_overloaded_lines,",
                    "max_lmp_error_per_mwh,load_weighted_lmp_error_per_mwh,",
                    "repair_feasible,repair_redispatch_mw,repair_cost_gap_pct,",
                    "in_milp_training_set,n_overloaded_internal_lines,",
                    "n_overloaded_external_lines,total_load_mw,",
                    "n_spuriously_binding_lines")
            for h in eachindex(val.scenario_ids)
                println(io, join((val.scenario_ids[h], val.original_status[h],
                    val.reduced_status[h], val.original_optimal[h] ? 1 : 0,
                    val.reduced_optimal[h] ? 1 : 0, val.dispatch_feasible[h] ? 1 : 0,
                    val.original_objective[h], val.reduced_objective[h],
                    val.objective_change_pct[h], val.objective_within_tolerance[h] ? 1 : 0,
                    mva * val.max_overload[h], val.max_relative_overload[h],
                    val.max_utilization[h], mva * val.total_overload[h],
                    val.n_violating_lines[h], val.max_lmp_error[h],
                    val.load_weighted_lmp_error[h], val.repair_feasible[h] ? 1 : 0,
                    mva * val.repair_redispatch[h], val.repair_cost_pct[h],
                    val.in_training[h] ? 1 : 0, val.n_violating_internal[h],
                    val.n_violating_external[h], mva * val.total_load[h],
                    val.n_spurious_binding[h]), ','))
            end
        end

        open(joinpath(dir, "line_overload_summary.csv"), "w") do io
            println(io, "line,from_bus_id,to_bus_id,internal,protected,",
                    "dcopf_overload_hours,dcopf_worst_overload_mw,",
                    "benchmark_overload_hours_vs_rating,",
                    "spuriously_binding_hours,worst_reduced_utilization")
            protected_set2 = Set(findall(r.protected))
            for l in 1:c.base.Ln
                println(io, join((l, c.bus_ids[c.base.Efrom[l]], c.bus_ids[c.base.Eto[l]],
                    val.internal_line[l] ? 1 : 0, l in protected_set2 ? 1 : 0,
                    val.line_violation_hours[l], mva * val.line_worst_overload[l],
                    bench.overload_hours_by_line[l],
                    val.line_spurious_hours[l],
                    val.line_worst_reduced_utilization[l]), ','))
            end
        end
    end

    open(joinpath(dir, "active_scenarios.csv"), "w") do io
        println(io, "scenario_id,timestamp,from_seed,total_load_mw,total_abs_flow_mw")
        seed_set = Set(selected)
        for s in art.active
            println(io, join((c.scenario_ids[s],
                first_hour + Hour(c.scenario_ids[s] - 1),
                s in seed_set ? 1 : 0,
                mva * sum(c.load[:, s]),
                mva * sum(abs, c.fhat[:, s])), ','))
        end
    end
    if !isnothing(art.gen)
        open(joinpath(dir, "scenario_generation_history.csv"), "w") do io
            println(io, "iteration,n_active,n_retained,n_internal_lines,solve_time_s,",
                    "worst_violation_above_mw,worst_violation_below_mw,",
                    "worst_line_above,worst_line_below,n_violating_hours,",
                    "both_sides,structural_fix,added_scenario_ids")
            for h in art.gen.history
                println(io, join((h.iteration, h.n_active, h.n_retained, h.n_internal,
                    h.solve_time, mva * h.worst_above, mva * h.worst_below,
                    h.worst_line_above, h.worst_line_below, h.n_violating_hours,
                    h.both_sides ? 1 : 0, h.structural ? 1 : 0,
                    "\"" * join(h.added_ids, ";") * "\""), ','))
            end
        end
    end

    # Kron-specific: one row per chain, and the boundary-scope internal
    # transfer (the closest available substitute for the dropped
    # selected_internal_transfer_mw.csv -- chain-interior lines never had an
    # individual transfer variable once folded into the Schur complement).
    open(joinpath(dir, "kron_chains.csv"), "w") do io
        println(io, "chain,x_bus_id,y_bus_id,n_interior,interior_bus_ids,",
                "Dx_eq,frate_eq_mw,epsL_eq_mw,driving_line,ended_up_internal")
        for (ci, ch) in enumerate(kron_map.chains)
            # BOUNDARY-level internal/external, from art.chain_internal --
            # NOT derived from r.rep_of, which (when collapse_external_chains
            # is on) always shows interior buses folded into x regardless of
            # the equivalent line's own internal/external status.
            internal_now = art.chain_internal[ci]
            println(io, join((ci, c.bus_ids[ch.x], c.bus_ids[ch.y], length(ch.interior),
                "\"" * join(c.bus_ids[ch.interior], ";") * "\"",
                ch.Dx_eq, mva * ch.frate_eq, mva * ch.epsL_eq, ch.driving_line,
                internal_now ? 1 : 0), ','))
        end
    end
    open(joinpath(dir, "kron_boundary_internal_transfer_mw.csv"), "w") do io
        println(io, "# boundary-network line index -> gint (MW). Scope: boundary")
        println(io, "# network only -- chain-interior lines have no individual")
        println(io, "# transfer variable once folded into the Schur complement.")
        writedlm(io, mva .* art.r_boundary.gint, ',')
    end

    open(joinpath(dir, "README.txt"), "w") do io
        println(io, has_calendar ?
            "Month: $(monthname(cfg.month)) $(cfg.year)" :
            "Case: single operating point (S = 1), no calendar")
        println(io, "Relaxation: $(TxReport.relaxation_pretty(art.mode, art.delta))")
        println(io, "Congestion definition: abs(flow)/rating >= $(cfg.near_limit_threshold)")
        println(io, "Kron reduction: $(kron_map.full_N) -> $(kron_map.boundary_N) buses ",
                "($(length(kron_map.chains)) chains)")
        println(io, "MILP scenario IDs: $(join(selection.scenario_ids, ", "))")
        println(io, "Those cover all $(length(selection.congested_lines)) congested lines.")
        println(io, "Every loaded scenario was used for the post-solve benchmark.")
        println(io)
        println(io, "line_status.csv / bus_mapping.csv / *.csv are all at FULL-network")
        println(io, "granularity, unfolded from the boundary (Kron-reduced) solve -- see")
        println(io, "kron_chains.csv for exactly which buses/lines were collapsed and")
        println(io, "whether each chain ended up internal (merged) or reinserted as-is.")
    end
    return dir
end

end # module KronReport

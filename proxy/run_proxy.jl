# Proxy reduction: one MILP with the injections fixed at each design demand.
#
#   julia --project=. --startup-file=no proxy/run_proxy.jl [key=value ...]
#
# Edit SETTINGS, or override any of them on the command line, e.g.
#   julia --project=. --startup-file=no proxy/run_proxy.jl case=case300 eps=0.05
#
# The proxy does not certify feasibility, so judge it by how far off the reduced
# dispatch is on the full network (overload MW / % of rating, cost and LMP error),
# reported per setting and for held-out demands.

SETTINGS = (
    # --- case and demands (see common/cases.jl) ---
    case              = :case118,
    demands           = :single,     # :single | :scaled | :hourly
    n_demands         = 9,           # :scaled load levels, or :hourly design hours (max)
    scale_range       = (0.9, 1.1),  # :scaled
    month             = 3,           # :hourly
    horizon_days      = 7,           # :hourly
    horizon_start_day = 1,           # :hourly
    linear_costs      = false,       # true = same cost model as kkt/ and greedy/

    # --- reduction ---
    eps               = 0.1,         # flow-error window, fraction of each rating
    near_limit        = 0.8,         # lines loaded at or above this stay external
    relaxation        = [(:conservative, 0.001)],  # (mode, delta) per solve: :none | :conservative | :symmetric
    lmp_separation    = true,
    lmp_threshold     = 3.0,         # $/MWh

    # --- limits (nothing = no limit) ---
    budget            = nothing,     # max merged lines
    hop_cap           = nothing,     # max merged chain length inside a cluster
    size_cap          = nothing,     # max buses per cluster

    # --- extras ---
    kron              = false,       # eliminate interior chain buses first (no caps)
    derate            = false,       # derate ratings of the reduced network afterwards
    scenario_generation = false,     # :hourly only, add violating hours and re-solve

    # --- run ---
    time_limit        = 600.0,       # seconds per solve
    opf_time_limit    = 60.0,
    plots             = true,
    open_plots        = false,
    export_matpower   = false,       # also write the reduced network as a .m file
    output_dir        = nothing,     # nothing = outputs/proxy/<case>/<tag>
)

using LinearAlgebra, Statistics, Dates, DelimitedFiles, Printf
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "common", "settings.jl"))
S = Main.Settings.apply_overrides(SETTINGS)
S.demands === :hourly && !isfile(joinpath(ROOT, "common", "multiscenario.jl")) &&
    error("demands = :hourly needs common/multiscenario.jl and the ACTIVSg scenario data, "
          * "which are not in the public repo")
S.demands in (:single, :scaled, :hourly) || error("demands must be :single, :scaled or :hourly")
S.kron && any(!isnothing, (S.budget, S.hop_cap, S.size_cap)) &&
    error("caps are not supported with kron = true")
S.scenario_generation && S.demands !== :hourly &&
    error("scenario_generation needs demands = :hourly")

include(joinpath(ROOT, "common", "cases.jl"))
include(joinpath(ROOT, "common", "caps.jl"))
@eval module TNR
    include(joinpath(dirname(@__DIR__), "common", "preprocessing.jl"))
    include(joinpath(@__DIR__, "proxy_model.jl"))
    include(joinpath(dirname(@__DIR__), "common", "postprocessing.jl"))
    Main.S.demands === :hourly &&
        include(joinpath(dirname(@__DIR__), "common", "multiscenario.jl"))
    if Main.S.kron
        include(joinpath(@__DIR__, "kron", "kron_core.jl"))
        include(joinpath(@__DIR__, "kron", "kron_unfold.jl"))
    end
end
TR = Main.TNR
include(joinpath(@__DIR__, "reporting.jl"))
include(joinpath(ROOT, "common", "matpower_export.jl"))
S.kron && include(joinpath(@__DIR__, "kron", "kron_reporting.jl"))
if S.derate
    isfile(joinpath(@__DIR__, "derating", "derate_core.jl")) ||
        error("derate = true needs proxy/derating/, which is not in the public repo")
    include(joinpath(@__DIR__, "derating", "derate_core.jl"))
end
DR = S.derate ? Main.Derating : nothing

function run_tag(s)
    parts = [string(s.demands, s.demands === :single ? "" : s.n_demands), "eps$(s.eps)"]
    s.kron && push!(parts, "kron")
    isnothing(s.budget) || push!(parts, "budget$(s.budget)")
    isnothing(s.hop_cap) || push!(parts, "hop$(s.hop_cap)")
    isnothing(s.size_cap) || push!(parts, "size$(s.size_cap)")
    s.derate && push!(parts, "derate")
    return join(parts, "_")
end
CASEFILE = Main.Cases.case_file(S.case)
OUTPUT_DIR = isnothing(S.output_dir) ?
    joinpath(ROOT, "outputs", "proxy", string(S.case), run_tag(S)) : abspath(S.output_dir)

# The reporting code reads this structure; :single/:scaled carry no calendar.
cfg = (
    normalized_error_threshold = S.eps,
    near_limit_threshold = S.near_limit,
    relaxation_sweep = S.relaxation,
    line_budget = S.budget, hop_cap = S.hop_cap, size_cap = S.size_cap,
    kron_reduction = (enabled=S.kron, min_chain_length=1, collapse_external_chains=true),
    derating = (certify=true, cost_cap_pct=0.5, empirical_scenarios=:all,
                certify_scenarios=:selected, max_empirical_iters=20, max_certify_iters=15),
    opf_time_limit = S.opf_time_limit,
    solve_time_limit = S.time_limit,
    numeric_focus = 3,
    screening_tolerance = 1e-6,
    cycle_cut_lens = (2, 3, 4),
    internal_bound_scale = 3.0,
    internal_rating_bound = true,
    switch_form = :hull,
    gurobi_log = true,
    merge_exact_blocks = false, merge_exact_mode = :fix, merge_leaf_blocks = true,
    lmp_separation = S.lmp_separation, lmp_threshold = S.lmp_threshold, lmp_relax_pmin = true,
    solve_time_benchmark = true, benchmark_repeats = 5, benchmark_threads = 1,
    dcopf = (scope=:month, time_limit=300.0, relax_pmin=true, relative_tolerance=1e-3,
             objective_tolerance_pct=1, lmp_tolerance=1e-3, measure_repair=true),
    show = (case=true, kron_chains=true, generation=true, reduction=true,
            model_check=true, benchmark=true, solve_time=true, dcopf=true,
            dcopf_failures=true, original_model=true, sweep_table=true,
            plots=S.plots, open_plots=S.open_plots, write_files=true,
            export_matpower=S.export_matpower),
    casefile = CASEFILE,
    output_dir = OUTPUT_DIR,
)
cfg = S.demands === :hourly ?
    merge(cfg, (case_name=string(S.case), month=S.month, year=Main.Cases.CASE_YEAR[S.case],
                matrix_dir=joinpath(ROOT, "outputs", "$(S.case)_dcopf"),
                operating_point_scale=1.0, line_limit_scale=1.0,
                seed_ranking=:flow, min_scenarios=S.n_demands, max_scenarios=S.n_demands,
                scenario_generation=S.scenario_generation, max_sc_generation_iterations=10)) :
    merge(cfg, (scenario_generation=false, max_sc_generation_iterations=0))

println("proxy reduction -> ", OUTPUT_DIR)
Main.Settings.print_settings(S)
Main.Settings.save_settings(joinpath(OUTPUT_DIR, "settings.txt"), S)

preamble = IOBuffer()
heldout = nothing
if S.demands in (:single, :scaled)
    c, heldout = Main.Cases.load_demands(TR, S)
    println("case = ", basename(CASEFILE), "   ", c.base.N, " buses, ", c.base.Ln,
            " lines, ", size(c.p, 2), " design demand(s)",
            isnothing(heldout) ? "" : ", $(size(heldout.p, 2)) held out")
    scenario_indices = collect(axes(c.p, 2))
    selection = TR.all_scenarios_selection(c, scenario_indices;
                                           congestion_threshold=something(S.near_limit, 0.9999))
else
    c = TR.build_multiscenario_tx_case(CASEFILE, cfg.matrix_dir;
            time_limit=cfg.solve_time_limit, validate_saved_ratings=false)
    month_idx = TR.month_scenario_indices(c, S.month; year=cfg.year, require_complete=true)
    first_hour = (S.horizon_start_day - 1) * 24 + 1
    last_hour = min(first_hour + S.horizon_days * 24 - 1, length(month_idx))
    c = TR.subset_multiscenario_case(c, month_idx[first_hour:last_hour])
    if S.linear_costs
        c.base.c2 .= 0.0
        c = TR.redispatch_dc_opf_scenarios(c; relax_pmin=true, time_limit=S.opf_time_limit)
    end
    scenario_indices = collect(axes(c.p, 2))
    selection = TR.select_seed_scenarios(c, scenario_indices;
            ranking=cfg.seed_ranking, congestion_threshold=something(S.near_limit, 0.9999),
            min_scenarios=cfg.min_scenarios, max_scenarios=cfg.max_scenarios)
    TxReport.report_case(c, selection, cfg, scenario_indices; io=preamble)
end

selected = selection.scenario_indices
epsL = cfg.normalized_error_threshold .* c.base.frate

if S.kron
    c_boundary, epsL_boundary, kron_map = TR.kron_reduce_case(c, epsL;
        near_limit_threshold=cfg.near_limit_threshold,
        eligibility_scenario_indices=scenario_indices,
        min_chain_length=cfg.kron_reduction.min_chain_length)
    KronReport.report_kron_chains(c, kron_map; io=preamble)
    rows, artifacts = KronReport.sweep_kron_multiscenario(
        TR, c, c_boundary, epsL, epsL_boundary, kron_map,
        selected, scenario_indices, cfg)
else
    rows, artifacts = TxReport.sweep_multiscenario(
        TR, c, epsL, selected, scenario_indices, cfg)
end

# One summary row per setting: size of the reduction and how far off the reduced
# dispatch is on the full network, for design and held-out demands.
function dispatch_error(val, mva)
    isnothing(val) && return (feasible="n/a", worst_overload_mw=NaN,
                              worst_overload_pct=NaN, worst_cost_change_pct=NaN,
                              worst_lmp_error=NaN)
    return (feasible="$(val.n_dispatch_feasible)/$(length(val.scenario_indices))",
            worst_overload_mw=round(mva * val.worst_overload; digits=3),
            worst_overload_pct=round(100 * val.worst_relative_overload; digits=3),
            worst_cost_change_pct=round(val.worst_abs_objective_change_pct; digits=4),
            worst_lmp_error=round(val.worst_lmp_error; digits=3))
end
summary_rows = NamedTuple[]

report_preamble = String(take!(preamble))
reports = String[]

for art in artifacts
    # ---- derating: a POST-reduction step, run BEFORE the reports ------------
    # The clustering never changes -- only the ratings the reduced network is
    # shipped with. It runs first so every section below describes the network
    # as it would actually be delivered, not the un-derated one.
    #
    # ONE rating per line, not one per scenario: the binding scenario sets it,
    # so a derating that holds in the worst hour holds in all of them.
    Frat, dinfo, val_d, dmaps = nothing, nothing, nothing, nothing
    if S.derate
        pick(w) = w === :all ? scenario_indices : selected
        emp_idx, cert_idx = pick(cfg.derating.empirical_scenarios),
                            pick(cfg.derating.certify_scenarios)
        dem(idx) = [c.load[:, s] for s in idx]
        # Use the validated/exported clustering. Older Kron artifacts stored
        # the display clustering in art.A, so retain the A_full fallback.
        Aval = get(art, :A_full, art.A)
        maps = DR.reduction_maps(Aval)
        dmaps = maps
        F0 = c.base.frate

        emp = DR.empirical_derate(c.base, maps.rep_of, maps.retained, dem(emp_idx),
                copy(F0); relax_pmin=cfg.dcopf.relax_pmin,
                max_iters=cfg.derating.max_empirical_iters, verbose=false)
        Frat = emp.Frat
        cert = cfg.derating.certify ?
            DR.certify_derating(c.base, maps.rep_of, maps.retained, dem(cert_idx),
                copy(Frat); relax_pmin=cfg.dcopf.relax_pmin,
                max_iters=cfg.derating.max_certify_iters,
                cost_cap_pct=cfg.derating.cost_cap_pct, verbose=false) : nothing
        isnothing(cert) || (Frat = cert.Frat)

        derated = [l for l in 1:c.base.Ln if Frat[l] < F0[l] - 1e-9]
        dinfo = (n_emp=length(emp_idx), n_cert=length(cert_idx),
                 empirical_ids=c.scenario_ids[emp_idx],
                 certificate_ids=c.scenario_ids[cert_idx],
                 emp=emp, cert=cert, derated=derated,
                 worst_pct=isempty(derated) ? 0.0 :
                     maximum(100 * (F0[l] - Frat[l]) / F0[l] for l in derated))

        # Re-validate with the delivered ratings, using the parent validator
        # unmodified apart from being told which ratings to use.
        val_d = isnothing(art.val) ? nothing :
            TR.validate_reduced_dcopf_scenarios(c, Aval, art.val.scenario_indices;
                line_ratings=Frat,
                relax_pmin=cfg.dcopf.relax_pmin, time_limit=cfg.dcopf.time_limit,
                relative_tolerance=cfg.dcopf.relative_tolerance,
                objective_tolerance_pct=cfg.dcopf.objective_tolerance_pct,
                lmp_tolerance=cfg.dcopf.lmp_tolerance,
                measure_repair=cfg.dcopf.measure_repair)
    end

    TxReport.report_digest(c, art, scenario_indices; derate=dinfo, val_derated=val_d)

    buf = IOBuffer()
    print(buf, report_preamble)
    println(buf)
    println(buf, repeat("=", 78))
    println(buf, "RESULTS:  ", TxReport.relaxation_pretty(art.mode, art.delta))
    println(buf, repeat("=", 78))
    cfg.show.generation && !isnothing(art.gen) &&
        TxReport.report_scenario_generation(c, art.gen; io=buf)
    cfg.show.reduction && TxReport.report_reduction(c, art.r, art.A; io=buf)
    if cfg.show.model_check
        if S.kron
            # The boundary point is the only one the solver itself produced;
            # everything below is checked on the true full network instead.
            TxReport.section("Boundary-solve feasibility (is the solver's own boundary-network point feasible?)"; io=buf)
            TR.print_model_feasibility_check_multiscenario(art.chk; io=buf)
        else
            TxReport.report_model_check(TR, art.chk; io=buf)
        end
    end
    cfg.show.benchmark  && TxReport.report_benchmark(c, art.bench, scenario_indices; io=buf)
    cfg.show.solve_time && !isnothing(art.timing) &&
        TR.report_dcopf_solve_times(art.timing; io=buf)
    cfg.show.dcopf && !isnothing(art.val) && TxReport.report_dcopf(c, art.val; io=buf)
    TxReport.report_derating_scope(c, dinfo; io=buf)
    if cfg.show.dcopf && !isnothing(val_d)
        TxReport.section("4) Full-vs-reduced DC-OPF validation  --  WITH derating"; io=buf)
        TxReport.report_dcopf(c, val_d; io=buf)
    end
    cfg.show.dcopf_failures && !isnothing(art.val) &&
        TxReport.report_dcopf_failures(c, art.val, art; io=buf)
    # The Kron sweep has no `cons`: the original A-based model is a statement
    # about the full network, and the Kron solve never formulates one.
    cfg.show.original_model && !isnothing(get(art, :cons, nothing)) &&
        TxReport.report_original_model_constraints(art.cons; io=buf)

    dir = !cfg.show.write_files ? nothing :
        S.kron ?
            KronReport.write_kron_multiscenario_outputs(
                joinpath(cfg.output_dir, art.label), c, art, art.kron_map,
                selection, selected, scenario_indices, cfg) :
            TxReport.write_multiscenario_outputs(
                joinpath(cfg.output_dir, art.label), c, art,
                selection, selected, scenario_indices, cfg)

    if !isnothing(dir)
        if S.derate
            open(joinpath(dir, "derating_certificate.txt"), "w") do certio
                TxReport.report_derating_scope(c, dinfo; io=certio)
            end
            open(joinpath(dir, "derated_ratings.csv"), "w") do io_
                println(io_, "line,from_bus_id,to_bus_id,internal,rating_mw,derated_mw")
                for l in 1:c.base.Ln
                    println(io_, join((l, c.bus_ids[c.base.Efrom[l]],
                        c.bus_ids[c.base.Eto[l]],
                        dmaps.rep_of[c.base.Efrom[l]] ==
                            dmaps.rep_of[c.base.Eto[l]] ? 1 : 0,
                        c.base.baseMVA * c.base.frate[l],
                        c.base.baseMVA * Frat[l]), ','))
                end
            end
        end

        get(cfg.show, :export_matpower, false) && export_reduced_matpower(
            cfg.casefile, c, art.r,
            joinpath(dir, "reduced_" * basename(cfg.casefile)))

        if cfg.show.plots
            plot_path = joinpath(dir, "full_and_reduced_network.png")
            plot_kwargs = (
                open_after=cfg.show.open_plots,
                # The model's OWN protected set. Without it the plot falls back
                # to the |fhat| = frate rule on c.base.fhat -- a stale single-OPF
                # snapshot that knows nothing about the monthly scenarios.
                binding_lines=findall(art.r.protected),
                congested_lines=selection.congested_lines,
                congested_line_labels=Dict(
                    l => "L$l  $(c.bus_ids[c.base.Efrom[l]])-$(c.bus_ids[c.base.Eto[l]])"
                    for l in selection.congested_lines),
                label_congested=true, label_retained=false, io=buf)
            if S.kron
                KronReport.plot_kron_network(
                    joinpath(ROOT, "common", "plots.jl"),
                    joinpath(@__DIR__, "kron", "kron_plots.jl"),
                    c, kron_map, c_boundary, art.A;
                    path=plot_path, plot_kwargs...)
            else
                TxReport.plot_network(
                    joinpath(ROOT, "common", "plots.jl"), c.base, art.A;
                    path=plot_path, plot_kwargs...)
            end
        end
    end



    # Held-out demands and the per-setting summary row.
    val_h = nothing
    if !isnothing(heldout)
        val_h = TR.validate_reduced_dcopf_scenarios(heldout, round.(Int, get(art, :A_full, art.A)),
            collect(axes(heldout.p, 2)); relax_pmin=cfg.dcopf.relax_pmin,
            time_limit=cfg.dcopf.time_limit, relative_tolerance=cfg.dcopf.relative_tolerance,
            objective_tolerance_pct=cfg.dcopf.objective_tolerance_pct,
            lmp_tolerance=cfg.dcopf.lmp_tolerance, measure_repair=cfg.dcopf.measure_repair)
        TxReport.section("Held-out demands: full-vs-reduced DC-OPF validation"; io=buf)
        TxReport.report_dcopf(heldout, val_h; io=buf)
    end
    let A = round.(Int, get(art, :A_full, art.A)), base = c.base, mva = base.baseMVA
        rep_of = [findfirst(==(1), A[:, b]) for b in 1:base.N]
        internal = [rep_of[base.Efrom[l]] == rep_of[base.Eto[l]] for l in 1:base.Ln]
        nb = length(unique(rep_of))
        d, h = dispatch_error(art.val, mva), dispatch_error(val_h, mva)
        row = (approach="proxy", case=string(S.case), setting=art.label, buses=base.N,
               remaining_buses=nb, reduction_pct=round(100 * (base.N - nb) / base.N; digits=2),
               merged_lines=count(internal), lines=base.Ln,
               largest_cluster=maximum(values(Main.Caps.cluster_sizes(base, internal))),
               hop_diameter=Main.Caps.hop_diameter(base, internal),
               solve_status=string(art.r.status),
               design_feasible=d.feasible, design_worst_overload_mw=d.worst_overload_mw,
               design_worst_overload_pct=d.worst_overload_pct,
               design_worst_cost_change_pct=d.worst_cost_change_pct,
               design_worst_lmp_error=d.worst_lmp_error,
               heldout_feasible=h.feasible, heldout_worst_overload_mw=h.worst_overload_mw,
               heldout_worst_overload_pct=h.worst_overload_pct,
               heldout_worst_cost_change_pct=h.worst_cost_change_pct,
               heldout_worst_lmp_error=h.worst_lmp_error)
        push!(summary_rows, row)
        isnothing(dir) || open(joinpath(dir, "summary.csv"), "w") do io_
            println(io_, join(keys(row), ","))
            println(io_, join(values(row), ","))
        end
    end
    text = String(take!(buf))
    push!(reports, text)
    isnothing(dir) || write(joinpath(dir, "report.txt"), text)
end

# Cross-setting comparison. Pairs with relaxation_comparison.csv, so it goes to
# a .txt next to it rather than to the console.
sweep_table = ""
if cfg.show.sweep_table
    buf = IOBuffer()
    TxReport.print_relaxation_table(rows; io=buf)
    sweep_table = String(take!(buf))
end

if cfg.show.write_files
    mkpath(cfg.output_dir)
    TxReport.write_relaxation_table(
        joinpath(cfg.output_dir, "relaxation_comparison.csv"), rows)
    isempty(sweep_table) ||
        write(joinpath(cfg.output_dir, "relaxation_comparison.txt"), sweep_table)
    println()
    println("Wrote all outputs to ", cfg.output_dir)
    println("Full report   -> <setting>/report.txt   (or print(reports[1]) from the REPL)")
end

if !isempty(summary_rows)
    mkpath(cfg.output_dir)
    open(joinpath(cfg.output_dir, "summary.csv"), "w") do io_
        println(io_, join(keys(first(summary_rows)), ","))
        foreach(r -> println(io_, join(values(r), ",")), summary_rows)
    end
    println()
    println("=" ^ 72)
    for r in summary_rows
        @printf("proxy  %s  %s:  %d -> %d buses (%.1f%%), %d merged lines, largest cluster %d\n",
                r.case, r.setting, r.buses, r.remaining_buses, r.reduction_pct,
                r.merged_lines, r.largest_cluster)
        @printf("  design   dispatch feasible %s   worst overload %.3f MW (%.3f%% of rating)   worst cost change %.4f%%\n",
                r.design_feasible, r.design_worst_overload_mw, r.design_worst_overload_pct,
                r.design_worst_cost_change_pct)
        isnothing(heldout) || @printf("  held-out dispatch feasible %s   worst overload %.3f MW (%.3f%% of rating)   worst cost change %.4f%%\n",
                r.heldout_feasible, r.heldout_worst_overload_mw, r.heldout_worst_overload_pct,
                r.heldout_worst_cost_change_pct)
    end
    println("results -> ", cfg.output_dir)
    println("=" ^ 72)
end

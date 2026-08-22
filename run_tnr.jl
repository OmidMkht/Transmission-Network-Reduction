# --------------------------------------------------------------------------- #
# THE runner. One file, four combinations:
#
#   julia --project=. --startup-file=no run_tnr.jl
#
#                       scenarios = :single        scenarios = :multi
#   kron = false        one .m operating point     hourly matrices, one month
#   kron = true         + chain elimination        + chain elimination
#
# Set the two switches in the RUN block below; everything after it is shared by
# all four, which is the point of the file. The model, preprocessing, validation
# and reporting are identical in every combination, so a change is tested the
# same way whichever one you run.
#
# :multi needs tnr_multiscenario.jl and the ACTIVSg scenario matrices, neither
# of which is published. A clean clone runs :single; asking for :multi without
# them stops with a message naming what is missing.
#
# The objective maximizes the number of INTERNAL LINES -- see tnr_model.jl.
# --------------------------------------------------------------------------- #

using LinearAlgebra, Statistics, Dates, DelimitedFiles

# ================================================== #
RUN = (
    scenarios = :multi,   # :single = one operating point | :multi = hourly matrices
    kron      = true,     # true = collapse degree-2 chains before the MILP

    # --- :single ---------------------------------------------------------- #
    casefile = joinpath(@__DIR__, "case studies", "pglib_opf_case118_ieee.m"),
    case_tag = "case118",  # names the output directory

    # --- :multi ----------------------------------------------------------- #
    case_name = "ACTIVSg200",   # "ACTIVSg200" | "ACTIVSg2000"
    month     = 3,
    # How much of the month to actually run. nothing = all of it. Every hour in
    # the horizon costs one full-network DC-OPF plus one reduced-network DC-OPF
    # in validation, so 7 days is ~4x cheaper end to end than 31. It does NOT
    # make the MILP smaller -- that is min/max_scenarios -- but it shrinks the
    # redispatch, screening and validation passes, which is where the time goes.
    horizon_days      = 7,
    horizon_start_day = 1,
)
# =========================================================================== #

RUN.scenarios in (:single, :multi) ||
    error("RUN.scenarios must be :single or :multi, got $(RUN.scenarios)")

# ACTIVSg2000's scenarios run on a LEAP calendar (8784 hours), ACTIVSg200's do
# not (8760). month_scenario_indices maps scenario_id -> date through this year,
# so the wrong one silently shifts every month after February by a day. One
# table, read by both the Kron and the non-Kron path: this used to be typed once
# per runner, and the two copies disagreed.
const CASE_YEAR = Dict("ACTIVSg200" => 2017, "ACTIVSg2000" => 2016)

const MULTI_FILE = joinpath(@__DIR__, "tnr_multiscenario.jl")
if RUN.scenarios === :multi
    isfile(MULTI_FILE) || error("""
        scenarios = :multi needs $(basename(MULTI_FILE)), which is not part of the
        published repo (it only runs against the ACTIVSg scenario matrices, which
        are not published either). Use scenarios = :single, or restore both.""")
    haskey(CASE_YEAR, RUN.case_name) || error(
        "No calendar year known for $(RUN.case_name). Add it to CASE_YEAR: a case " *
        "with 8784 hourly scenarios needs a leap year, one with 8760 a non-leap one.")
end

# The ONLY place include() is used. Preprocessing must come first: it defines
# the case structs every other file's methods dispatch on. The multi-scenario
# and Kron files are pulled in only when the run asks for them.
@eval module TNR
    include(joinpath(@__DIR__, "tnr_preprocessing.jl"))
    include(joinpath(@__DIR__, "tnr_model.jl"))
    include(joinpath(@__DIR__, "tnr_postprocessing.jl"))
    Main.RUN.scenarios === :multi &&
        include(joinpath(@__DIR__, "tnr_multiscenario.jl"))
    if Main.RUN.kron
        include(joinpath(@__DIR__, "kron_reduction", "kron_core.jl"))
        include(joinpath(@__DIR__, "kron_reduction", "kron_unfold.jl"))
    end
end
TR = Main.TNR
include(joinpath(@__DIR__, "tnr_reporting.jl"))
include(joinpath(@__DIR__, "matpower_export.jl"))
RUN.kron && include(joinpath(@__DIR__, "kron_reduction", "kron_reporting.jl"))

# ========================== CONFIGURATION ================================== #
# Shared by all four combinations. The handful of mode-specific fields are
# merged in below -- :single deliberately carries no month/year, which is how
# the reporting code knows the run has no calendar.
common = (
    # --- reduction -------------------------------------------------------- #
    normalized_error_threshold = 0.10,   # eps, as a fraction of each rating
    # nothing = only lines EXACTLY at their rating are protected. A value in
    # (0,1) also protects everything loaded at or above that fraction -- and
    # doubles as the Kron eligibility screen (E5): no chain line may reach it.
    near_limit_threshold = 0.80,

    # (mode, delta) per solve. delta is a fraction of each line's rating.
    #   :none          exact pin on congested lines -- the unrelaxed baseline
    #   :conservative  window opens toward MORE congestion only, so the reduced
    #                  network may over-state a corridor's loading, never
    #                  under-state it
    #   :symmetric     opens both ways -- twice the freedom, but it can
    #                  under-state congestion, which is what lets a reduced
    #                  dispatch overload the line back on the full network
    # Each entry is solved, benchmarked and validated on its own and gets its
    # own output subdirectory; all of them land in one comparison table.
    relaxation_sweep = [
        (:none,         0.00),
        # (:conservative, 0.01),
        # (:conservative, 0.03),
        # (:symmetric,    0.01),
        # (:symmetric,    0.03),
    ],

    # --- Kron preprocessing (ignored when RUN.kron is false) --------------- #
    #   min_chain_length  skip chains with fewer interior buses than this. 1
    #                     accepts every valid chain; raising it keeps the MILP's
    #                     freedom on short chains at the cost of size.
    #   collapse_external_chains  REPORTING ONLY -- fold a chain's interior into
    #                     its anchor even when the equivalent line stayed
    #                     external. Validation always uses the fully reinserted
    #                     network either way, so the safety numbers do not move;
    #                     n_retained and the plot do.
    kron_reduction = (
        enabled                  = RUN.kron,
        min_chain_length         = 1,
        collapse_external_chains = true,
    ),

    # --- solver ------------------------------------------------------------ #
    opf_time_limit      = 2 * 60.0,
    solve_time_limit    = 5 * 60.0,
    # NOT optional in practice. G spans several orders of magnitude against
    # constraint sensitivities of ~1e-4, so with NumericFocus off Gurobi can
    # prune wrongly and report a FALSE optimum. Measured on case118: it returned
    # "OPTIMAL" at 134 internal lines when 136 was achievable and provable, and
    # the model_feasibility_check passed it as GENUINE -- the answer was simply
    # 2 lines short. numeric_focus=3 returned 136 and was also 6x FASTER
    # (69 s vs 432 s). Set to nothing only to reproduce that failure.
    numeric_focus       = 3,
    screening_tolerance = 1e-6,
    cycle_cut_lens      = (2, 3, 4),   # () = off; 2 = parallel, 3 = triangles, 4 = 4-cycles
    # Internal-transfer bound: min(x * H_l, G), x = internal_bound_scale. G is
    # the proven but loose global bound; H_l is a tighter per-line heuristic.
    internal_bound_scale = 3.0,
    gurobi_log          = false,   # <output_dir>/<setting>/gurobi.log

    # --- exact-merge preprocessing (bridges, protected-free leaf blocks) ---- #
    #   merge_exact_blocks  false = off, the model decides these itself.
    #   merge_exact_mode    :fix  = add c_l == 1. A DOMINANCE fixing, so it
    #                               cannot change the optimum -- it only deletes
    #                               search. Fastest.
    #                       :warm = the same set as a MIP start only, leaving the
    #                               solver free to disagree. Slower, but it lets
    #                               you VERIFY the argument on a new case.
    merge_exact_blocks  = false,
    merge_exact_mode    = :fix,    # :fix | :warm
    merge_leaf_blocks   = true,    # false = bridges only

    # --- LMP separation ---------------------------------------------------- #
    # Two buses whose prices differ by more than lmp_threshold should not share
    # a cluster: merging them replaces two prices with one and no window
    # accuracy recovers the difference. Enforced the cheap way -- the SHORTEST
    # path between them may not go fully internal: sum_{l in P} c_l <= |P| - 1.
    # NECESSARY, NOT SUFFICIENT: blocking one path does not block every path, so
    # a constrained pair can still merge through a detour. The solve reports how
    # many pairs actually came out separated; read that, do not assume.
    lmp_separation      = false,
    lmp_threshold       = 1,       # $/MWh, worst gap over the ACTIVE scenarios
    lmp_relax_pmin      = true,    # match the DC-OPF validation's pmin handling

    # --- computational payoff: full vs reduced DC-OPF solve cost ----------- #
    # The reason to reduce a network is that the reduced one is cheaper to
    # solve, so measure it. Both networks are built by ONE builder, identical
    # solver settings, benchmark_repeats solves each after a warm-up, min taken.
    #   work units  deterministic, machine-independent -- TRUST THIS ONE
    #   wall clock  real time, but sensitive to other load on the machine
    #   solve_time  Gurobi's own Runtime; reads 0.0 under a millisecond
    solve_time_benchmark = true,
    benchmark_repeats    = 5,
    benchmark_threads    = 1,      # 1 keeps timings reproducible

    # --- DC-OPF validation ------------------------------------------------- #
    dcopf = (
        scope                   = :month,   # :month | :selected | :none
        time_limit              = 2 * 60.0,
        relax_pmin              = true,
        relative_tolerance      = 1e-3,
        objective_tolerance_pct = 0.1,
        lmp_tolerance           = 1e-3,
        measure_repair          = true,     # cheapest fix for an infeasible dispatch
    ),

    # --- which sections are produced --------------------------------------- #
    # The console always gets the digest only (reduction, window benchmark,
    # DC-OPF pass/fail). These decide what else is computed and written to
    # <setting>/report.txt -- printing is nearly free now, so leaving them on
    # costs only the sections that do real work: plots and original_model.
    show = (
        case           = true,   # :multi only -- month, congestion, seed choice
        kron_chains    = true,   # RUN.kron only -- what the elimination did
        generation     = true,   # per-iteration trace of scenario generation
        reduction      = true,
        model_check    = true,   # is the SOLVER'S OWN point feasible?
        benchmark      = true,
        solve_time     = true,   # full-vs-reduced DC-OPF cost
        dcopf          = true,
        dcopf_failures = true,   # anatomy of the hours that failed the strict test
        original_model = true,   # would the original assignment-based model accept this?
        sweep_table    = true,
        plots          = true,
        open_plots     = false,
        write_files    = true,
        # write the reduced network as a MATPOWER .m file too -- see
        # matpower_export.jl. Off by default: most runs are exploratory, not
        # meant to be published as an example case.
        export_matpower = false,
    ),
)

kron_tag = RUN.kron ? "kron_edge" : "edge"
cfg = RUN.scenarios === :single ?
    merge(common, (
        casefile   = RUN.casefile,
        output_dir = joinpath(@__DIR__, "outputs", "$(RUN.case_tag)_$(kron_tag)"),
        # S = 1, so there is nothing for scenario generation to add.
        scenario_generation          = false,
        max_sc_generation_iterations = 0,
    )) :
    merge(common, (
        case_name  = RUN.case_name,
        month      = RUN.month,
        year       = CASE_YEAR[RUN.case_name],
        casefile   = joinpath(@__DIR__, "case studies", RUN.case_name,
                              "case_$(RUN.case_name).m"),
        matrix_dir = joinpath(@__DIR__, "outputs", "$(RUN.case_name)_dcopf"),
        output_dir = joinpath(@__DIR__, "outputs",
                              "$(RUN.case_name)_$(kron_tag)_month$(RUN.month)"),
        # --- operating point ---
        operating_point_scale = 1.0,   # scales load/generation before redispatch
        line_limit_scale      = 1.0,   # scales every rating, then redispatches
        # --- scenario selection ---
        # SEED: hours ranked biggest-first by :flow (sum |line flow|) or :load
        # (total demand), unioned with a minimum cover of the congested lines.
        # The cover is a hard requirement -- a line congested anywhere in the
        # month must be congested in some seed hour, or it will not be
        # protected. Ranking fills the rest up to min_scenarios, never past max.
        seed_ranking  = :flow,   # :flow | :load
        min_scenarios = 20,
        max_scenarios = 20,
        # GENERATION: solve on the active set, screen the clustering against all
        # monthly hours (a linear solve each -- no MILP), add the worst-violating
        # hour(s), re-solve warm-started. On convergence the clustering is
        # window-feasible for the WHOLE month, which no fixed set can promise.
        scenario_generation          = false,
        max_sc_generation_iterations = 10,
    ))
# =========================================================================== #

# Preamble to every setting's report.txt: what the case is and what the Kron
# step did. Both happen once, before the sweep, so they are the same for all
# settings.
preamble = IOBuffer()

# ---------------------------- build the case -------------------------------- #
if RUN.scenarios === :single
    c = TR.build_single_scenario_case(cfg.casefile; time_limit=cfg.opf_time_limit)
    println("case = ", basename(cfg.casefile), "   ",
            c.base.N, " buses, ", c.base.Ln, " lines, ", size(c.p, 2), " scenario")
    scenario_indices = collect(axes(c.p, 2))
    selection = TR.all_scenarios_selection(c, scenario_indices;
            congestion_threshold=isnothing(cfg.near_limit_threshold) ? 0.9999 :
                                 cfg.near_limit_threshold)
else
    # Hours that must have been solved for the requested month to be complete:
    # the day-of-year of that month's LAST day, times 24.
    hours_needed = Dates.dayofyear(Date(cfg.year, cfg.month,
                        daysinmonth(Date(cfg.year, cfg.month, 1)))) * 24
    isdir(cfg.matrix_dir) || error("""
        No scenario matrices for $(cfg.case_name) at $(cfg.matrix_dir).
        Build them with:
            julia --startup-file=no run_scenario_dcopf.jl $(cfg.case_name) $(hours_needed)
        ($(hours_needed) hours reaches the end of month $(cfg.month) in $(cfg.year).)""")
    println("Loading network and scenario matrices from ", cfg.matrix_dir)
    c = TR.build_multiscenario_tx_case(cfg.casefile, cfg.matrix_dir;
            time_limit=cfg.solve_time_limit, validate_saved_ratings=false)
    month_idx = TR.month_scenario_indices(c, cfg.month;
            year=cfg.year, require_complete=true)
    # Trim to horizon_days BEFORE subsetting, so every later pass -- redispatch,
    # screening, DC-OPF validation -- sees only the hours being studied. The
    # completeness check above still runs on the full month first, so a
    # truncated matrix build is caught rather than silently shortening this.
    if !isnothing(RUN.horizon_days)
        first_hour = (RUN.horizon_start_day - 1) * 24 + 1
        last_hour = first_hour + RUN.horizon_days * 24 - 1
        first_hour >= 1 && first_hour <= length(month_idx) ||
            error("horizon_start_day $(RUN.horizon_start_day) is outside month $(cfg.month)")
        last_hour > length(month_idx) && (last_hour = length(month_idx))
        println("Horizon: days $(RUN.horizon_start_day)-$(fld(last_hour, 24)) of ",
                monthname(cfg.month), " $(cfg.year)  ",
                "($(last_hour - first_hour + 1) of $(length(month_idx)) hours)")
        month_idx = month_idx[first_hour:last_hour]
    end
    c = TR.subset_multiscenario_case(c, month_idx)
    c = TR.scale_operating_points(c, cfg.operating_point_scale;
            relax_pmin=cfg.dcopf.relax_pmin, time_limit=cfg.dcopf.time_limit)
    c = TR.scale_line_limits_and_redispatch(c, cfg.line_limit_scale;
            relax_pmin=cfg.dcopf.relax_pmin, time_limit=cfg.dcopf.time_limit)
    scenario_indices = collect(axes(c.p, 2))
    selection = TR.select_seed_scenarios(c, scenario_indices;
            ranking=cfg.seed_ranking,
            congestion_threshold=cfg.near_limit_threshold,
            min_scenarios=cfg.min_scenarios, max_scenarios=cfg.max_scenarios)
    cfg.show.case && TxReport.report_case(c, selection, cfg, scenario_indices; io=preamble)
end

selected = selection.scenario_indices
epsL = cfg.normalized_error_threshold .* c.base.frate

# ------------------------------ solve --------------------------------------- #
if RUN.kron
    c_boundary, epsL_boundary, kron_map = TR.kron_reduce_case(c, epsL;
        near_limit_threshold=cfg.near_limit_threshold,
        eligibility_scenario_indices=scenario_indices,
        min_chain_length=cfg.kron_reduction.min_chain_length)
    cfg.show.kron_chains &&
        KronReport.report_kron_chains(c, kron_map; io=preamble)
    rows, artifacts = KronReport.sweep_kron_multiscenario(
        TR, c, c_boundary, epsL, epsL_boundary, kron_map,
        selected, scenario_indices, cfg)
else
    rows, artifacts = TxReport.sweep_multiscenario(
        TR, c, epsL, selected, scenario_indices, cfg)
end

# ----------------------------- report --------------------------------------- #
# The console gets the digest only -- reduction, window benchmark, DC-OPF
# pass/fail. Every other section is still produced in full; it just goes to a
# string, which lands in <setting>/report.txt and stays here in `reports`, so
# print(reports[1]) brings the whole thing back without re-running anything.
report_preamble = String(take!(preamble))
reports = String[]

for art in artifacts
    TxReport.report_digest(c, art, scenario_indices)

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
        if RUN.kron
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
    cfg.show.dcopf_failures && !isnothing(art.val) &&
        TxReport.report_dcopf_failures(c, art.val, art; io=buf)
    # The Kron sweep has no `cons`: the original A-based model is a statement
    # about the full network, and the Kron solve never formulates one.
    cfg.show.original_model && !isnothing(get(art, :cons, nothing)) &&
        TxReport.report_original_model_constraints(art.cons; io=buf)

    dir = !cfg.show.write_files ? nothing :
        RUN.kron ?
            KronReport.write_kron_multiscenario_outputs(
                joinpath(cfg.output_dir, art.label), c, art, art.kron_map,
                selection, selected, scenario_indices, cfg) :
            TxReport.write_multiscenario_outputs(
                joinpath(cfg.output_dir, art.label), c, art,
                selection, selected, scenario_indices, cfg)

    if !isnothing(dir)
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
            if RUN.kron
                KronReport.plot_kron_network(
                    joinpath(@__DIR__, "transmission_plots.jl"),
                    joinpath(@__DIR__, "kron_reduction", "kron_plots.jl"),
                    c, kron_map, c_boundary, art.A;
                    path=plot_path, plot_kwargs...)
            else
                TxReport.plot_network(
                    joinpath(@__DIR__, "transmission_plots.jl"), c.base, art.A;
                    path=plot_path, plot_kwargs...)
            end
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

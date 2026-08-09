# --------------------------------------------------------------------------- #
# Single-scenario runner: one pglib/MATPOWER test case (case118, case300, ...).
#
#   julia --project=. --startup-file=no run_tnr.jl
#
# There is no separate single-scenario model: a test case is built as a case
# with S = 1 (`build_single_scenario_case`) and goes through the same model,
# preprocessing, validation and reporting as a multi-hour horizon would.
#
# READ THIS BEFORE INTERPRETING n_retained. The objective maximizes the number
# of INTERNAL LINES, not the number of buses removed -- see tnr_model.jl.
# --------------------------------------------------------------------------- #

using LinearAlgebra, Statistics, Dates, DelimitedFiles

# The ONLY place include() is used. Preprocessing must come first: it defines
# the case structs every other file's methods dispatch on.
@eval module TNR
    include(joinpath(@__DIR__, "tnr_preprocessing.jl"))
    include(joinpath(@__DIR__, "tnr_model.jl"))
    include(joinpath(@__DIR__, "tnr_postprocessing.jl"))
end
TR = Main.TNR
include(joinpath(@__DIR__, "tnr_reporting.jl"))

# ========================== CONFIGURATION ================================== #
cfg = (
    # --- inputs ---
    casefile   = joinpath(@__DIR__, "case studies", "pglib_opf_case118_ieee.m"),
    output_dir = joinpath(@__DIR__, "outputs", "case118_edge"),

    # --- reduction ---
    normalized_error_threshold = 0.10,   # eps, as a fraction of each rating
    # nothing = only lines EXACTLY at their rating are protected. A value in
    # (0,1) also protects everything loaded at or above that fraction.
    near_limit_threshold = 0.8,

    relaxation_sweep = [(:none, 0.00)],  # a single, unrelaxed solve

    # --- solver ---
    opf_time_limit      = 2 * 60.0,
    solve_time_limit    = 5 * 60.0,
    cycle_cut_lens      = (2, 3, 4),   # () = off; 2 = parallel lines, 3 = triangles, 4 = chordless 4-cycles
    # NOT optional. Without it Gurobi's presolve can prune wrongly on a
    # numerically wide model and report a false optimum -- see tnr_model.jl.
    numeric_focus       = 3,
    screening_tolerance = 1e-6,
    gurobi_log          = false,   # <output_dir>/<setting>/gurobi.log

    # --- scenario settings (fixed: there is only one operating point) ---
    scenario_generation = false,

    # --- computational payoff: full vs reduced DC-OPF solve cost -------------
    # The reason to reduce a network is that the reduced one is cheaper to
    # solve, so measure it: both networks built by one builder, identical
    # solver settings, best of a few repeats after a warm-up.
    solve_time_benchmark = true,

    # --- DC-OPF validation ---
    dcopf = (
        scope                  = :month,   # :month = the one scenario | :none
        time_limit             = 2 * 60.0,
        relax_pmin              = true,
        relative_tolerance     = 1e-3,
        objective_tolerance_pct = 0.1,
        lmp_tolerance          = 1e-3,
        measure_repair         = true,
    ),

    # --- which validations to display / write ---
    show = (
        reduction      = true,
        model_check    = true,   # is the SOLVER'S OWN point feasible?
        benchmark      = true,
        solve_time     = true,   # full-vs-reduced DC-OPF cost
        dcopf          = true,
        dcopf_failures = true,
        original_model = true,   # would the original assignment-based model accept this?
        sweep_table    = true,
        plots          = true,
        open_plots     = false,
        write_files    = true,
    ),
)
# =========================================================================== #

c = TR.build_single_scenario_case(cfg.casefile; time_limit=cfg.opf_time_limit)
println("case = ", basename(cfg.casefile), "   ",
        c.base.N, " buses, ", c.base.Ln, " lines, ", size(c.p, 2), " scenario")

scenario_indices = collect(axes(c.p, 2))
selection = TR.all_scenarios_selection(c, scenario_indices;
        congestion_threshold=isnothing(cfg.near_limit_threshold) ? 0.9999 :
                             cfg.near_limit_threshold)
selected = selection.scenario_indices
epsL = cfg.normalized_error_threshold .* c.base.frate

rows, artifacts = TxReport.sweep_multiscenario(TR, c, epsL, selected, scenario_indices, cfg)

for art in artifacts
    println()
    println(repeat("=", 78))
    println("RESULTS:  ", TxReport.relaxation_pretty(art.mode, art.delta))
    println(repeat("=", 78))
    cfg.show.reduction   && TxReport.report_reduction(c, art.r, art.A)
    # This has to come BEFORE the checks below: they discard the solver's own
    # (f, vartheta, g) and recompute from the recovered clustering, so none of
    # them would notice that the returned point was itself infeasible.
    cfg.show.model_check && TxReport.report_model_check(TR, art.chk)
    cfg.show.benchmark   && TxReport.report_benchmark(c, art.bench, scenario_indices)
    cfg.show.solve_time && !isnothing(art.timing) && TR.report_dcopf_solve_times(art.timing)
    cfg.show.dcopf && !isnothing(art.val) && TxReport.report_dcopf(c, art.val)
    cfg.show.dcopf_failures && !isnothing(art.val) &&
        TxReport.report_dcopf_failures(c, art.val, art)
    cfg.show.original_model && !isnothing(art.cons) &&
        TxReport.report_original_model_constraints(art.cons)

    cfg.show.write_files || continue
    dir = TxReport.write_multiscenario_outputs(
        joinpath(cfg.output_dir, art.label), c, art,
        selection, selected, scenario_indices, cfg)

    cfg.show.plots || continue
    TxReport.plot_network(
        joinpath(@__DIR__, "transmission_plots.jl"), c.base, art.A;
        path=joinpath(dir, "full_and_reduced_network.png"),
        open_after=cfg.show.open_plots,
        title="$(basename(cfg.casefile)):  $(c.base.N) -> $(art.r.n_retained) buses  ·  " *
              "eps = $(cfg.normalized_error_threshold)",
        binding_lines=findall(art.r.protected),
        congested_lines=selection.congested_lines,
        label_congested=true, label_retained=false)
end

cfg.show.sweep_table && TxReport.print_relaxation_table(rows)
if cfg.show.write_files
    mkpath(cfg.output_dir)
    TxReport.write_relaxation_table(
        joinpath(cfg.output_dir, "relaxation_comparison.csv"), rows)
    println()
    println("Wrote all outputs to ", cfg.output_dir)
end

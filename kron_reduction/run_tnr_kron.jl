# --------------------------------------------------------------------------- #
# runner: Kron preprocessing + reduction, one pglib/MATPOWER test case.
#
# julia --project=. --startup-file=no kron_reduction/run_tnr_kron.jl
#
# Same case and same operating point as ../run_tnr.jl, with one step added in
# front of the MILP: every maximal chain of degree-2, generator-free,
# uncongested buses is collapsed by a Kron reduction (Schur complement) into a
# single equivalent line, the MILP solves on that smaller "boundary" network,
# and the clustering is unfolded back onto the full bus set before anything is
# validated. See reference/kron_preprocessing.pdf for the formulation.
#
# Set kron.enabled = false to get ../run_tnr.jl's answer through this same
# path, which is how the two are meant to be compared.
# --------------------------------------------------------------------------- #

using LinearAlgebra, Statistics, Dates, DelimitedFiles

@eval module TNR
    include(joinpath(@__DIR__, "..", "tnr_preprocessing.jl"))
    include(joinpath(@__DIR__, "..", "tnr_model.jl"))
    include(joinpath(@__DIR__, "..", "tnr_postprocessing.jl"))
    include(joinpath(@__DIR__, "kron_core.jl"))
    include(joinpath(@__DIR__, "kron_unfold.jl"))
end
TR = Main.TNR
include(joinpath(@__DIR__, "..", "tnr_reporting.jl"))
include(joinpath(@__DIR__, "kron_reporting.jl"))

# ========================== CONFIGURATION ================================== #
cfg = (
    # --- inputs
    casefile   = joinpath(@__DIR__, "..", "case studies", "pglib_opf_case118_ieee.m"),
    output_dir = joinpath(@__DIR__, "..", "outputs", "case118_kron"),

    # --- reduction
    normalized_error_threshold = 0.10,   # eps, as a fraction of each rating
    # nothing = only lines EXACTLY at their rating are protected. A value in
    # (0,1) also protects everything loaded at or above that fraction -- and
    # doubles as the Kron eligibility screen (E5): no chain line may reach it.
    near_limit_threshold = 0.8,

    relaxation_sweep = [(:none, 0.00)],  # a single, unrelaxed solve

    # --- Kron preprocessing
    #   enabled           false = no elimination, i.e. run_tnr.jl's problem.
    #   min_chain_length  skip chains with fewer interior buses than this.
    #                     1 accepts every valid chain; raising it keeps the
    #                     MILP's freedom on short chains at the cost of size.
    #   collapse_external_chains  reporting only -- fold a chain's interior
    #                     into its anchor even when the equivalent line stayed
    #                     external. Validation always uses the fully
    #                     reinserted network either way.
    kron_reduction = (
        enabled                  = true,
        min_chain_length         = 1,
        collapse_external_chains = true,
    ),

    # --- solver ---
    opf_time_limit      = 2 * 60.0,
    solve_time_limit    = 5 * 60.0,
    numeric_focus       = 3,
    screening_tolerance = 1e-6,

    # cuts
    cycle_cut_lens      = (2, 3, 4),   # () = off; 2 = parallel lines, 3 = triangles, 4 = chordless 4-cycles

    gurobi_log          = false,   # <output_dir>/<setting>/gurobi.log

    # --- scenario settings (fixed: there is only one operating point) ---
    scenario_generation = false,

    # --- computational payoff
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
        kron_chains    = true,
        reduction      = true,
        model_check    = true,   # is the SOLVER'S OWN point feasible?
        benchmark      = true,
        solve_time     = true,   # full-vs-reduced DC-OPF cost
        dcopf          = true,
        dcopf_failures = true,
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

if cfg.kron_reduction.enabled
    c_boundary, epsL_boundary, kron_map = TR.kron_reduce_case(c, epsL;
        near_limit_threshold=cfg.near_limit_threshold,
        eligibility_scenario_indices=scenario_indices,
        min_chain_length=cfg.kron_reduction.min_chain_length)
    cfg.show.kron_chains && KronReport.report_kron_chains(c, kron_map)
    rows, artifacts = KronReport.sweep_kron_multiscenario(
        TR, c, c_boundary, epsL, epsL_boundary, kron_map,
        selected, scenario_indices, cfg)
else
    rows, artifacts = TxReport.sweep_multiscenario(TR, c, epsL, selected, scenario_indices, cfg)
end

for art in artifacts
    println()
    println(repeat("=", 78))
    println("RESULTS:  ", TxReport.relaxation_pretty(art.mode, art.delta))
    println(repeat("=", 78))
    cfg.show.reduction && TxReport.report_reduction(c, art.r, art.A)
    if cfg.show.model_check
        if cfg.kron_reduction.enabled
            # The boundary point is the only one the solver itself produced;
            # everything below is checked on the true full network instead.
            TxReport.section("Boundary-solve feasibility (is the solver's own boundary-network point feasible?)")
            TR.print_model_feasibility_check_multiscenario(art.chk)
        else
            TxReport.report_model_check(TR, art.chk)
        end
    end
    cfg.show.benchmark  && TxReport.report_benchmark(c, art.bench, scenario_indices)
    cfg.show.solve_time && !isnothing(art.timing) && TR.report_dcopf_solve_times(art.timing)
    cfg.show.dcopf && !isnothing(art.val) && TxReport.report_dcopf(c, art.val)
    cfg.show.dcopf_failures && !isnothing(art.val) &&
        TxReport.report_dcopf_failures(c, art.val, art)

    cfg.show.write_files || continue
    dir = cfg.kron_reduction.enabled ?
        KronReport.write_kron_multiscenario_outputs(
            joinpath(cfg.output_dir, art.label), c, art, art.kron_map,
            selection, selected, scenario_indices, cfg) :
        TxReport.write_multiscenario_outputs(
            joinpath(cfg.output_dir, art.label), c, art,
            selection, selected, scenario_indices, cfg)

    cfg.show.plots || continue
    plot_path = joinpath(dir, "full_and_reduced_network.png")
    plot_kwargs = (
        open_after=cfg.show.open_plots,
        binding_lines=findall(art.r.protected),
        congested_lines=selection.congested_lines,
        label_congested=true, label_retained=false)
    if cfg.kron_reduction.enabled
        KronReport.plot_kron_network(
            joinpath(@__DIR__, "..", "transmission_plots.jl"),
            joinpath(@__DIR__, "kron_plots.jl"),
            c, kron_map, c_boundary, art.A;
            path=plot_path, plot_kwargs...)
    else
        TxReport.plot_network(
            joinpath(@__DIR__, "..", "transmission_plots.jl"), c.base, art.A;
            path=plot_path, plot_kwargs...)
    end
end

cfg.show.sweep_table && TxReport.print_relaxation_table(rows)
if cfg.show.write_files
    mkpath(cfg.output_dir)
    TxReport.write_relaxation_table(
        joinpath(cfg.output_dir, "relaxation_comparison.csv"), rows)
    println()
    println("Wrote all outputs to ", cfg.output_dir)
end

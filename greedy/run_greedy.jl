# Greedy network reduction.
#
#   julia --project=. --startup-file=no greedy/run_greedy.jl [key=value ...]
#
# Edit SETTINGS, or override any of them on the command line, e.g.
#   julia --project=. --startup-file=no greedy/run_greedy.jl case=case300 hop_cap=5

SETTINGS = (
    # --- case and demands (see common/cases.jl) ---
    case              = :case118,
    demands           = :scaled,     # :single | :scaled | :hourly
    n_demands         = 9,           # :scaled load levels, or :hourly design hours
    scale_range       = (0.9, 1.1),  # :scaled
    month             = 3,           # :hourly
    horizon_days      = 7,           # :hourly
    horizon_start_day = 1,           # :hourly

    # --- acceptance test ---
    cost_cap          = 0.1,         # % over the full-network optimum; nothing = off
    flow_tol          = 1e-9,        # allowed overload, fraction of rating
    cost_tol          = 0.0,         # relative cost slack for a deliverable dispatch
    kkt_check         = true,        # final joint KKT check, rolls merges back if it fails

    # --- limits (nothing = no limit) ---
    budget            = nothing,     # max merged lines
    hop_cap           = nothing,     # max merged chain length inside a cluster
    size_cap          = nothing,     # max buses per cluster

    # --- search ---
    ordering          = :flow,       # :flow (ranked, re-ranked) | :loading (static)
    norm              = :headroom,   # :headroom | :rating
    alpha             = 0.5,         # 2-norm weight against the max-norm
    radial_first      = true,

    # --- run ---
    time_limit        = 600.0,       # seconds
    threads           = 1,           # per LP; the checks are many small LPs
    opf_time_limit    = 60.0,
    output_dir        = nothing,     # nothing = outputs/greedy/<case>/<tag>
)

using Dates, Printf, LinearAlgebra, DelimitedFiles
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "common", "settings.jl"))
S = merge(Main.Settings.apply_overrides(SETTINGS), (linear_costs=true,))
S.demands === :hourly && !isfile(joinpath(ROOT, "common", "multiscenario.jl")) &&
    error("demands = :hourly needs common/multiscenario.jl and the ACTIVSg scenario data, "
          * "which are not in the public repo")

@eval module Data
    include(joinpath(dirname(@__DIR__), "common", "preprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "postprocessing.jl"))
    Main.S.demands === :hourly &&
        include(joinpath(dirname(@__DIR__), "common", "multiscenario.jl"))
end
include(joinpath(ROOT, "common", "cases.jl"))
include(joinpath(ROOT, "common", "caps.jl"))
include(joinpath(ROOT, "common", "audit.jl"))
include(joinpath(ROOT, "common", "results.jl"))
include(joinpath(ROOT, "kkt", "kkt_model.jl"))
include(joinpath(@__DIR__, "flow_sensitivity.jl"))
include(joinpath(@__DIR__, "greedy.jl"))

function run_tag(s)
    parts = [string(s.demands, s.demands === :single ? "" : s.n_demands),
             string(s.ordering, s.ordering === :flow ? "-$(s.norm)" : ""),
             "tol$(s.flow_tol)"]
    isnothing(s.budget) || push!(parts, "budget$(s.budget)")
    isnothing(s.hop_cap) || push!(parts, "hop$(s.hop_cap)")
    isnothing(s.size_cap) || push!(parts, "size$(s.size_cap)")
    s.kkt_check || push!(parts, "nokkt")
    return join(parts, "_")
end

out = isnothing(S.output_dir) ?
    joinpath(ROOT, "outputs", "greedy", string(S.case), run_tag(S)) : abspath(S.output_dir)
println("greedy reduction -> ", out)
Main.Settings.print_settings(S)
Main.Settings.save_settings(joinpath(out, "settings.txt"), S)

design, heldout = Main.Cases.load_demands(Main.Data, S)
scope = collect(axes(design.load, 2))
@printf("\n%s: %d buses, %d lines, %d design demand(s), %d held out\n\n",
        S.case, design.base.N, design.base.Ln, length(scope),
        isnothing(heldout) ? 0 : size(heldout.load, 2))

r = Main.Greedy.greedy_reduction(Main.Audit, Main.BilevelReduction, design, scope;
        time_limit=S.time_limit, cost_gap_pct=S.cost_cap, flow_tolerance=S.flow_tol,
        cost_tolerance=S.cost_tol, line_budget=S.budget, hop_cap=S.hop_cap,
        size_cap=S.size_cap, ordering=S.ordering, norm=S.norm, alpha=S.alpha,
        radial_first=S.radial_first, kkt_check=S.kkt_check, threads=S.threads)
r.covered || error("the full network itself fails the acceptance test ($(r.reason)); " *
                   "try a looser flow_tol or cost_tol")
Main.Results.write_rows(joinpath(out, "trace.csv"), r.trace)

Main.Results.finish(Main.Audit, Main.BilevelReduction, design, heldout, r.internal, out;
    approach="greedy", case=string(S.case), cost_cap=S.cost_cap,
    extra=(status=r.reason, seconds=round(r.elapsed_seconds; digits=1),
           lp_checks=r.lp_checks, radial_merges=r.radial_taken,
           rejected=r.rejected, rollbacks=r.rollbacks))

# Greedy network reduction with a hop limit: no merged chain inside a cluster is
# longer than hop_cap lines.
#
#   julia --project=. --startup-file=no greedy/run_greedy_hop.jl [key=value ...]
#
# hop_cap is one number, or a ladder: one greedy pass per entry, each starting
# from what the previous one merged, e.g.
#   julia --project=. --startup-file=no greedy/run_greedy_hop.jl case=case300 hop_cap=5
#   julia --project=. --startup-file=no greedy/run_greedy_hop.jl case=case300 'hop_cap=[5,10,nothing]'
# Each step writes a row to steps.csv; the audit runs on the last step.

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
    kkt_check         = true,        # final joint KKT check per step, rolls merges back if it fails

    # --- hop limit ---
    hop_cap           = 5,           # a number, or a ladder like [5, 10, nothing]

    # --- search ---
    ordering          = :flow,       # :flow (ranked, re-ranked) | :loading (static)
    norm              = :headroom,   # :headroom | :rating
    alpha             = 0.5,         # 2-norm weight against the max-norm
    radial_first      = true,

    # --- run ---
    time_limit        = 600.0,       # seconds per step
    threads           = 1,           # per LP; the checks are many small LPs
    opf_time_limit    = 60.0,
    output_dir        = nothing,     # nothing = outputs/greedy/<case>/<tag>
)

HOP_LADDER = true
include(joinpath(@__DIR__, "pipeline.jl"))

# Proxy reduction with a hop limit: no merged chain inside a cluster is longer
# than hop_cap lines.
#
#   julia --project=. --startup-file=no proxy/run_proxy_hop.jl [key=value ...]
#
# hop_cap is one number, or a ladder: one solve per entry, each keeping what the
# previous one merged, e.g.
#   julia --project=. --startup-file=no proxy/run_proxy_hop.jl case=case300 hop_cap=5
#   julia --project=. --startup-file=no proxy/run_proxy_hop.jl case=case300 'hop_cap=[5,10,nothing]'
# A ladder takes one relaxation setting; steps.csv has a row per step and the
# last step gets the full report.
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

    # --- hop limit ---
    hop_cap           = 5,           # a number, or a ladder like [5, 10, nothing]

    # --- extras ---
    derate            = false,       # derate ratings of the reduced network afterwards
    scenario_generation = false,     # :hourly only, add violating hours and re-solve

    # --- run ---
    time_limit        = 600.0,       # seconds per step
    opf_time_limit    = 60.0,
    plots             = true,
    open_plots        = false,
    export_matpower   = false,       # also write the reduced network as a .m file
    output_dir        = nothing,     # nothing = outputs/proxy/<case>/<tag>
)

HOP_LADDER = true
include(joinpath(@__DIR__, "pipeline.jl"))

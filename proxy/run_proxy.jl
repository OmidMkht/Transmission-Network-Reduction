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
#
# One solve per relaxation setting. For a hop ladder (hop 5, then 10, then free)
# use proxy/run_proxy_hop.jl.

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
    cycle_cuts        = (2, 3, 4),   # short-cycle closure lengths; () turns them off

    # --- greedy warm start ---
    greedy_seed       = false,       # seed the solve with proxy/greedy_proxy.jl
    greedy_seed_time_limit = 600.0,  # seconds per seed (once per ladder rung)
    greedy_seed_verbose = false,     # print every accepted merge

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

HOP_LADDER = false
include(joinpath(@__DIR__, "pipeline.jl"))

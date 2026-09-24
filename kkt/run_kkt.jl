# Bilevel reduction solved as one MILP through the lower level's KKT conditions.
#
#   julia --project=. --startup-file=no kkt/run_kkt.jl [key=value ...]
#
# Edit SETTINGS, or override any of them on the command line, e.g.
#   julia --project=. --startup-file=no kkt/run_kkt.jl case=case300 hop_cap=5
#
# One solve. For a hop ladder (hop 5, then 10, then free) use kkt/run_kkt_hop.jl.

SETTINGS = (
    # --- case and demands (see common/cases.jl) ---
    case              = :case118,
    demands           = :scaled,     # :single | :scaled | :hourly
    n_demands         = 9,           # :scaled load levels, or :hourly design hours
    scale_range       = (0.9, 1.1),  # :scaled
    month             = 3,           # :hourly
    horizon_days      = 7,           # :hourly
    horizon_start_day = 1,           # :hourly

    # --- target ---
    cost_cap          = 0.1,         # % over the full-network optimum; nothing = off
    objective         = :lines,      # :lines (max merged lines) | :clusters (min buses)

    # --- limits (nothing = no limit) ---
    budget            = nothing,     # max merged lines
    hop_cap           = nothing,     # max merged chain length inside a cluster
    size_cap          = nothing,     # max buses per cluster

    # --- start ---
    start_from        = nothing,     # internal.csv of another run (e.g. greedy), kept merged
    radial            = :warm,       # :warm | :enforce | :none, safe radial merges
    greedy_seed       = false,       # seed the solve with greedy/greedy.jl
    greedy_seed_time_limit = 600.0,  # seconds per seed (once per ladder rung)
    greedy_flow_tol   = 1e-9,        # greedy's overload tolerance while seeding
    greedy_seed_verbose = false,     # print every accepted merge

    # --- run ---
    time_limit        = 600.0,       # seconds
    threads           = parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", string(Sys.CPU_THREADS))),
    mipgap            = 1e-4,
    opf_time_limit    = 60.0,
    output_dir        = nothing,     # nothing = outputs/kkt/<case>/<tag>
)

HOP_LADDER = false
include(joinpath(@__DIR__, "pipeline.jl"))

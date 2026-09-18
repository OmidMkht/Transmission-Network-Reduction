# Bilevel KKT reduction with a hop limit: no merged chain inside a cluster is
# longer than hop_cap lines.
#
#   julia --project=. --startup-file=no kkt/run_kkt_hop.jl [key=value ...]
#
# hop_cap is one number, or a ladder: one solve per entry, each keeping what the
# previous one merged, e.g.
#   julia --project=. --startup-file=no kkt/run_kkt_hop.jl case=case300 hop_cap=5
#   julia --project=. --startup-file=no kkt/run_kkt_hop.jl case=ACTIVSg200 'hop_cap=[10,20,nothing]'
# Each step writes a row to steps.csv; the audit runs on the last feasible step.

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

    # --- hop limit ---
    hop_cap           = 5,           # a number, or a ladder like [5, 10, nothing]

    # --- start ---
    start_from        = nothing,     # internal.csv of another run (e.g. greedy), kept merged
    radial            = :warm,       # :warm | :enforce | :none, safe radial merges

    # --- run ---
    time_limit        = 600.0,       # seconds per step
    threads           = parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", string(Sys.CPU_THREADS))),
    mipgap            = 1e-4,
    opf_time_limit    = 60.0,
    output_dir        = nothing,     # nothing = outputs/kkt/<case>/<tag>
)

HOP_LADDER = true
include(joinpath(@__DIR__, "pipeline.jl"))

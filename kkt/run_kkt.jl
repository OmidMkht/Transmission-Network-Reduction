# Bilevel reduction solved as one MILP through the lower level's KKT conditions.
#
#   julia --project=. --startup-file=no kkt/run_kkt.jl [key=value ...]
#
# Edit SETTINGS, or override any of them on the command line, e.g.
#   julia --project=. --startup-file=no kkt/run_kkt.jl case=case300 hop_cap=5
#
# A cap given as a list is a schedule: one solve per entry, each keeping what
# the previous one merged, e.g. hop_cap=[5,10,nothing] runs hop 5, hop 10, free.

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

    # --- limits (nothing = no limit; a list = schedule) ---
    budget            = nothing,     # max merged lines
    hop_cap           = nothing,     # max merged chain length inside a cluster
    size_cap          = nothing,     # max buses per cluster

    # --- start ---
    start_from        = nothing,     # internal.csv of another run (e.g. greedy), kept merged
    radial            = :warm,       # :warm | :enforce | :none, safe radial merges

    # --- run ---
    time_limit        = 600.0,       # seconds per solve (per step in a schedule)
    threads           = parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", string(Sys.CPU_THREADS))),
    mipgap            = 1e-4,
    opf_time_limit    = 60.0,
    output_dir        = nothing,     # nothing = outputs/kkt/<case>/<tag>
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
include(joinpath(@__DIR__, "kkt_model.jl"))
const BL = Main.BilevelReduction

# Caps as one entry per step.
as_steps(x) = x isa AbstractVector ? collect(x) : [x]
hops, sizes, budgets = as_steps(S.hop_cap), as_steps(S.size_cap), as_steps(S.budget)
nsteps = maximum(length.((hops, sizes, budgets)))
for (name, v) in (("hop_cap", hops), ("size_cap", sizes), ("budget", budgets))
    length(v) in (1, nsteps) || error("$name has $(length(v)) entries, other schedules have $nsteps")
end
at_step(v, k) = length(v) == 1 ? v[1] : v[k]
label(v) = join((isnothing(x) ? "free" : string(x) for x in v), "-")

function run_tag(s)
    parts = [string(s.demands, s.demands === :single ? "" : s.n_demands), string(s.objective)]
    any(!isnothing, hops) && push!(parts, "hop" * label(hops))
    any(!isnothing, sizes) && push!(parts, "size" * label(sizes))
    any(!isnothing, budgets) && push!(parts, "budget" * label(budgets))
    isnothing(s.start_from) || push!(parts, "start-" * basename(dirname(abspath(s.start_from))))
    return join(parts, "_")
end

out = isnothing(S.output_dir) ?
    joinpath(ROOT, "outputs", "kkt", string(S.case), run_tag(S)) : abspath(S.output_dir)
println("KKT bilevel reduction -> ", out)
Main.Settings.print_settings(S)
Main.Settings.save_settings(joinpath(out, "settings.txt"), S)

design, heldout = Main.Cases.load_demands(Main.Data, S)
base = design.base
scope = collect(axes(design.load, 2))
@printf("\n%s: %d buses, %d lines, %d design demand(s), %d held out, %d step(s)\n\n",
        S.case, base.N, base.Ln, length(scope),
        isnothing(heldout) ? 0 : size(heldout.load, 2), nsteps)

held = falses(base.Ln)
if !isnothing(S.start_from)
    held = BitVector(vec(readdlm(S.start_from, ',', Int)) .== 1)
    length(held) == base.Ln || error("start_from has $(length(held)) lines, case has $(base.Ln)")
    println("starting from ", S.start_from, ": ", count(held), " merged lines kept")
end
radial = BL.radial_internal_mask(base).internal
warm = Main.Caps.trim_to_caps(base, radial .| held; budget=at_step(budgets, 1),
                              hop_cap=at_step(hops, 1), size_cap=at_step(sizes, 1)) .| held

started = time()
best = nothing
trace = NamedTuple[]
for k in 1:nsteps
    capkw = (line_budget=at_step(budgets, k), hop_cap=at_step(hops, k),
             cluster_size_cap=at_step(sizes, k))
    r = BL.solve_bilevel_reduction(design, scope; capkw...,
        cost_gap_pct=S.cost_cap, objective=S.objective, radial_mode=S.radial,
        held_internal=held, warm_internal=warm, time_limit=S.time_limit,
        solver_threads=S.threads, mipgap=S.mipgap, t1_screening=true,
        cycle_cut_lens=(), solver_seed=0, output_flag=1, start_time_limit=60.0,
        log_file=joinpath(out, "gurobi_step$(k).log"))
    push!(trace, (step=k, budget=something(capkw.line_budget, "free"),
                  hop_cap=something(capkw.hop_cap, "free"),
                  size_cap=something(capkw.cluster_size_cap, "free"),
                  status=string(r.status), feasible=r.feasible,
                  buses=r.feasible ? r.n_retained : -1,
                  merged_lines=r.feasible ? count(r.internal) : -1,
                  bound=r.bound, seconds=round(time() - started; digits=1)))
    Main.Results.write_rows(joinpath(out, "steps.csv"), trace)
    @printf("step %d: %s, %s buses\n", k, r.status, r.feasible ? r.n_retained : "no")
    if r.feasible
        global best = r
        global held = copy(r.internal)
        global warm = copy(r.internal)
    end
end
isnothing(best) && error("no feasible reduction found; raise time_limit")

Main.Results.finish(Main.Audit, BL, design, heldout, best.internal, out;
    approach="kkt", case=string(S.case), cost_cap=S.cost_cap,
    extra=(status=string(best.status), bound=best.bound, steps=nsteps,
           seconds=round(time() - started; digits=1)))

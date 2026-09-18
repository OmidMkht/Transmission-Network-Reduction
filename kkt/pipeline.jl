# Shared by kkt/run_kkt.jl and kkt/run_kkt_hop.jl: load the case, solve, audit.
# The runner defines SETTINGS and HOP_LADDER before including this file.
#
# hop_cap as a list is a ladder: one solve per entry, each keeping what the
# previous one merged, e.g. [5, 10, nothing] runs hop 5, hop 10, then free.

using Dates, Printf, LinearAlgebra, DelimitedFiles
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "common", "settings.jl"))
S = merge(Main.Settings.apply_overrides(SETTINGS), (linear_costs=true,))
S.demands === :hourly && !isfile(joinpath(ROOT, "common", "multiscenario.jl")) &&
    error("demands = :hourly needs common/multiscenario.jl and the ACTIVSg scenario data, "
          * "which are not in the public repo")
!HOP_LADDER && S.hop_cap isa AbstractVector &&
    error("hop_cap is a list; for a hop ladder use kkt/run_kkt_hop.jl")

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

budget, size_cap = get(S, :budget, nothing), get(S, :size_cap, nothing)
(budget isa AbstractVector || size_cap isa AbstractVector) &&
    error("budget and size_cap take one value; only hop_cap has a ladder")
hops = S.hop_cap isa AbstractVector ? collect(S.hop_cap) : [S.hop_cap]
isempty(hops) && error("hop_cap ladder is empty")
nsteps = length(hops)
label(v) = join((isnothing(x) ? "free" : string(x) for x in v), "-")

function run_tag(s)
    parts = [string(s.demands, s.demands === :single ? "" : s.n_demands), string(s.objective)]
    any(!isnothing, hops) && push!(parts, "hop" * label(hops))
    isnothing(size_cap) || push!(parts, "size$(size_cap)")
    isnothing(budget) || push!(parts, "budget$(budget)")
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
warm = Main.Caps.trim_to_caps(base, radial .| held; budget=budget,
                              hop_cap=hops[1], size_cap=size_cap) .| held

started = time()
best = nothing
trace = NamedTuple[]
for k in 1:nsteps
    r = BL.solve_bilevel_reduction(design, scope;
        line_budget=budget, hop_cap=hops[k], cluster_size_cap=size_cap,
        cost_gap_pct=S.cost_cap, objective=S.objective, radial_mode=S.radial,
        held_internal=held, warm_internal=warm, time_limit=S.time_limit,
        solver_threads=S.threads, mipgap=S.mipgap, t1_screening=true,
        cycle_cut_lens=(), solver_seed=0, output_flag=1, start_time_limit=60.0,
        log_file=joinpath(out, nsteps == 1 ? "gurobi.log" : "gurobi_step$(k).log"))
    push!(trace, (step=k, hop_cap=something(hops[k], "free"),
                  status=string(r.status), feasible=r.feasible,
                  buses=r.feasible ? r.n_retained : -1,
                  merged_lines=r.feasible ? count(r.internal) : -1,
                  bound=r.bound, seconds=round(time() - started; digits=1)))
    Main.Results.write_rows(joinpath(out, "steps.csv"), trace)
    @printf("step %d (hop %s): %s, %s buses\n", k, something(hops[k], "free"),
            r.status, r.feasible ? r.n_retained : "no")
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

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
include(joinpath(ROOT, "greedy", "flow_sensitivity.jl"))
include(joinpath(ROOT, "greedy", "greedy.jl"))
const BL = Main.BilevelReduction

budget, size_cap = get(S, :budget, nothing), get(S, :size_cap, nothing)
(budget isa AbstractVector || size_cap isa AbstractVector) &&
    error("budget and size_cap take one value; only hop_cap has a ladder")
hops = S.hop_cap isa AbstractVector ? collect(S.hop_cap) : [S.hop_cap]
isempty(hops) && error("hop_cap ladder is empty")
nsteps = length(hops)
label(v) = join((isnothing(x) ? "free" : string(x) for x in v), "-")
# time_limit is one number for every step, or one per ladder step.
limits = S.time_limit isa AbstractVector ? Float64.(S.time_limit) : fill(Float64(S.time_limit), nsteps)
length(limits) == nsteps || error("time_limit has $(length(limits)) entries for $nsteps step(s)")

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

# Greedy seed. Same cost cap and the same caps as the rung it seeds, so the
# topology is admissible to topology_start_allowed instead of being silently
# swapped for the radial mask. Greedy's own kkt_check calls
# fixed_topology_kkt_start -- the very routine that builds the KKT start -- so a
# seed that verifies here is one the solve can complete.
function kkt_seed(force, hop)
    get(S, :greedy_seed, false) || return nothing
    g = Main.Greedy.greedy_reduction(Main.Audit, BL, design, scope;
        time_limit=get(S, :greedy_seed_time_limit, 600.0),
        cost_gap_pct=S.cost_cap, flow_tolerance=get(S, :greedy_flow_tol, 1e-9),
        cost_tolerance=1e-8, line_budget=budget, hop_cap=hop, size_cap=size_cap,
        ordering=:flow, norm=:headroom, alpha=0.5, radial_first=true,
        kkt_check=true, threads=S.threads,
        verbose=get(S, :greedy_seed_verbose, false),
        held_internal=any(force) ? force : nothing)
    @printf("greedy seed (hop %s): %d merged line(s), %d of %d buses, %s, %.1f s\n",
            something(hop, "free"), count(g.internal), g.buses, base.N,
            g.reason, g.elapsed_seconds)
    g.covered || error("greedy seed failed: $(g.reason)")
    return (internal=BitVector(g.internal), merged=count(g.internal), buses=g.buses)
end
SEEDED = get(S, :greedy_seed, false)
nothing_seed = (internal=nothing, merged=-1, buses=-1)

seed = something(kkt_seed(held, hops[1]), nothing_seed)
warm = isnothing(seed.internal) ?
    Main.Caps.trim_to_caps(base, radial .| held; budget=budget,
                           hop_cap=hops[1], size_cap=size_cap) .| held :
    seed.internal .| held

# Hold-forward fixes every rung's merges into the next. The hop cap only
# loosens, so the held set stays feasible. Warm-forward seeds the next rung with
# it instead and lets the solver undo a merge.
HOLD_FORWARD = get(S, :hold_forward, true)
HOLD_FORWARD || println("ladder: warm-forward (rungs are seeded, not forced)")

started = time()
best = nothing
trace = NamedTuple[]
prev_point = nothing   # the last rung's full solution, a fallback start for the next
for k in 1:nsteps
    log_file = joinpath(out, nsteps == 1 ? "gurobi.log" : "gurobi_step$(k).log")
    isfile(log_file) && rm(log_file)   # Gurobi appends; a rerun starts a clean log
    r = BL.solve_bilevel_reduction(design, scope;
        line_budget=budget, hop_cap=hops[k], cluster_size_cap=size_cap,
        cost_gap_pct=S.cost_cap, objective=S.objective, radial_mode=S.radial,
        held_internal=held, warm_internal=warm, time_limit=limits[k],
        solver_threads=S.threads, mipgap=S.mipgap, t1_screening=true,
        cycle_cut_lens=(), solver_seed=0, output_flag=1, start_time_limit=3600.0,
        require_warm_start=SEEDED, start_point=prev_point, log_file=log_file)
    start_obj = Main.Results.start_objective(log_file)
    row = (step=k, hop_cap=something(hops[k], "free"), time_limit=limits[k],
           seed_merged=seed.merged, seed_buses=seed.buses,
           start_from=!SEEDED ? "" : get(r, :start_from_point, false) ? "previous rung" : "greedy seed",
           start_objective=start_obj,
           status=string(r.status), feasible=r.feasible,
           buses=r.feasible ? r.n_retained : -1,
           merged_lines=r.feasible ? count(r.internal) : -1,
           bound=r.bound, solve_seconds=round(r.solve_time; digits=1),
           seconds=round(time() - started; digits=1))
    push!(trace, row)
    Main.Results.write_rows(joinpath(out, "steps.csv"), trace)
    @printf("step %d (hop %s): %s, %s buses\n", k, something(hops[k], "free"),
            r.status, r.feasible ? r.n_retained : "no")
    if r.feasible
        # Warm-forward lets a rung come back worse than an earlier one, so take
        # the best over rungs rather than the last. Under hold-forward the bus
        # count only falls, so this picks the last one anyway.
        (isnothing(best) || r.n_retained < best.r.n_retained) && (global best = (r=r, row=row))
        nsteps > 1 && writedlm(joinpath(out, "internal_step$(k).csv"), Int.(r.internal), ',')
        global prev_point = r.point
        HOLD_FORWARD && (global held = copy(r.internal))
        if k < nsteps
            nxt = kkt_seed(copy(r.internal), hops[k+1])
            global seed = something(nxt, nothing_seed)
            global warm = isnothing(nxt) ? copy(r.internal) : nxt.internal
        end
    end
end
isnothing(best) && error("no feasible reduction found; raise time_limit")

Main.Results.finish(Main.Audit, BL, design, heldout, best.r.internal, out;
    approach="kkt", case=string(S.case), cost_cap=S.cost_cap,
    info=Main.Results.run_info(mode=nsteps > 1 ? "ladder" : "noladder",
        setting="cap=$(S.cost_cap)", hop_caps=label(hops), status=best.r.status,
        merged=best.row.merged_lines, bound=best.r.bound,
        seed_merged=best.row.seed_merged, seed_buses=best.row.seed_buses,
        start_obj=best.row.start_objective,
        solve_seconds=sum(t.solve_seconds for t in trace),
        seconds=time() - started),
    extra=(steps=nsteps,))

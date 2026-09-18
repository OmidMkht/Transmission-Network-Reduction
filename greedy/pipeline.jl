# Shared by greedy/run_greedy.jl and greedy/run_greedy_hop.jl: load the case,
# merge, audit. The runner defines SETTINGS and HOP_LADDER before including this.
#
# hop_cap as a list is a ladder: one greedy pass per entry, each starting from
# what the previous one merged, e.g. [5, 10, nothing] runs hop 5, hop 10, then free.

using Dates, Printf, LinearAlgebra, DelimitedFiles
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "common", "settings.jl"))
S = merge(Main.Settings.apply_overrides(SETTINGS), (linear_costs=true,))
S.demands === :hourly && !isfile(joinpath(ROOT, "common", "multiscenario.jl")) &&
    error("demands = :hourly needs common/multiscenario.jl and the ACTIVSg scenario data, "
          * "which are not in the public repo")
!HOP_LADDER && S.hop_cap isa AbstractVector &&
    error("hop_cap is a list; for a hop ladder use greedy/run_greedy_hop.jl")

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

budget, size_cap = get(S, :budget, nothing), get(S, :size_cap, nothing)
(budget isa AbstractVector || size_cap isa AbstractVector) &&
    error("budget and size_cap take one value; only hop_cap has a ladder")
hops = S.hop_cap isa AbstractVector ? collect(S.hop_cap) : [S.hop_cap]
isempty(hops) && error("hop_cap ladder is empty")
nsteps = length(hops)
label(v) = join((isnothing(x) ? "free" : string(x) for x in v), "-")

function run_tag(s)
    parts = [string(s.demands, s.demands === :single ? "" : s.n_demands),
             string(s.ordering, s.ordering === :flow ? "-$(s.norm)" : ""),
             "tol$(s.flow_tol)"]
    isnothing(budget) || push!(parts, "budget$(budget)")
    any(!isnothing, hops) && push!(parts, "hop" * label(hops))
    isnothing(size_cap) || push!(parts, "size$(size_cap)")
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
@printf("\n%s: %d buses, %d lines, %d design demand(s), %d held out, %d step(s)\n\n",
        S.case, design.base.N, design.base.Ln, length(scope),
        isnothing(heldout) ? 0 : size(heldout.load, 2), nsteps)

started = time()
held = nothing
r = nothing
steps, trace = NamedTuple[], NamedTuple[]
totals = (lp_checks=0, radial_merges=0, rejected=0, rollbacks=0)
for k in 1:nsteps
    nsteps > 1 && @printf("\n--- step %d: hop %s ---\n", k, something(hops[k], "free"))
    global r = Main.Greedy.greedy_reduction(Main.Audit, Main.BilevelReduction, design, scope;
            time_limit=S.time_limit, cost_gap_pct=S.cost_cap, flow_tolerance=S.flow_tol,
            cost_tolerance=S.cost_tol, line_budget=budget, hop_cap=hops[k],
            size_cap=size_cap, ordering=S.ordering, norm=S.norm, alpha=S.alpha,
            radial_first=S.radial_first, kkt_check=S.kkt_check, threads=S.threads,
            held_internal=held)
    r.covered || error("the full network itself fails the acceptance test ($(r.reason)); " *
                       "try a looser flow_tol or cost_tol")
    append!(trace, [merge((step=k,), row) for row in r.trace])
    push!(steps, (step=k, hop_cap=something(hops[k], "free"), status=r.reason,
                  buses=r.buses, merged_lines=count(r.internal),
                  seconds=round(time() - started; digits=1)))
    Main.Results.write_rows(joinpath(out, "trace.csv"), trace)
    nsteps > 1 && Main.Results.write_rows(joinpath(out, "steps.csv"), steps)
    global totals = (lp_checks=totals.lp_checks + r.lp_checks,
                     radial_merges=totals.radial_merges + r.radial_taken,
                     rejected=totals.rejected + r.rejected,
                     rollbacks=totals.rollbacks + r.rollbacks)
    global held = copy(r.internal)
    @printf("step %d (hop %s): %s, %d buses\n", k, something(hops[k], "free"), r.reason, r.buses)
end

Main.Results.finish(Main.Audit, Main.BilevelReduction, design, heldout, r.internal, out;
    approach="greedy", case=string(S.case), cost_cap=S.cost_cap,
    extra=(status=r.reason, steps=nsteps, seconds=round(time() - started; digits=1),
           totals...))

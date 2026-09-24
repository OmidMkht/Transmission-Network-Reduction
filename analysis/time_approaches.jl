# Run AND time the three reduction approaches over a settings sweep.
#
#   julia --project=. --startup-file=no analysis/time_approaches.jl case=14bus
#   julia --project=. --startup-file=no analysis/time_approaches.jl case=case118 demands=scaled
#
# One pass does both jobs: it writes each setting's clustering to
# outputs/<approach>/<caseid>/<setting>/internal.csv, which analysis/review_toy.jl
# then scores, and outputs/<case>_review/compute_time.csv with the timings.
#
# Why timing needs its own harness: a runner is a fresh Julia process that calls
# its algorithm exactly once, so its reported seconds are mostly JIT -- greedy's
# "5.8 s" is 5.15 s of compilation and 0.09 s of work -- while the proxy reports
# JuMP's solver time alone and never includes JIT at all. Here the case is loaded
# once, every path is warmed up, those timings are thrown away, and only then is
# each setting timed. All three run on the same THREADS so it is like for like.

using JuMP, Gurobi, Printf, DelimitedFiles, Statistics

const ROOT = dirname(@__DIR__)
const THREADS = 8
const CASEID = Dict("6busww" => :case6ww, "14bus" => :case14toy, "case118" => :case118)

casename = "14bus"; demands = :hourly; n_train = 20; time_limit = 300
for a in ARGS
    k, v = split(a, Char(61); limit=2)
    if k == "case"; global casename = String(v)
    elseif k == "demands"; global demands = Symbol(v)
    elseif k == "n_train"; global n_train = parse(Int, v)
    elseif k == "time_limit"; global time_limit = parse(Int, v)
    else; error("accepts case=, demands=, n_train=, time_limit=")
    end
end
const TIME_LIMIT = time_limit
caseid = CASEID[casename]

SETTINGS = (case=caseid, demands=demands, n_demands=n_train, scale_range=(0.9, 1.1),
            month=[9, 10, 11], horizon_days=92, horizon_start_day=1,
            opf_time_limit=60.0, linear_costs=true)
include(joinpath(ROOT, "common", "settings.jl"))
S = merge(Main.Settings.apply_overrides(SETTINGS, String[]), (linear_costs=true,))

# proxy_model.jl is written to live INSIDE the data module (it refers to
# MultiScenarioTxReductionCase unqualified), exactly as proxy/pipeline.jl
# arranges it, so it is included there rather than at top level.
@eval module Data
    include(joinpath(dirname(@__DIR__), "common", "preprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "proxy", "proxy_model.jl"))
    include(joinpath(dirname(@__DIR__), "common", "postprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "multiscenario.jl"))
end
include(joinpath(ROOT, "common", "cases.jl"))
include(joinpath(ROOT, "common", "caps.jl"))
include(joinpath(ROOT, "common", "audit.jl"))
include(joinpath(ROOT, "kkt", "kkt_model.jl"))
include(joinpath(ROOT, "greedy", "flow_sensitivity.jl"))
include(joinpath(ROOT, "greedy", "greedy.jl"))

const BL = Main.BilevelReduction
const TR = Main.Data

design, heldout = Main.Cases.load_demands(Main.Data, S)
scope = collect(axes(design.load, 2))
base = design.base
@printf("[%s] %s: %d buses, %d lines, %d train / %d test scenarios   %d threads\n",
        casename, demands, base.N, base.Ln, length(scope),
        isnothing(heldout) ? 0 : size(heldout.load, 2), THREADS)

kkt_call(cap; tl=TIME_LIMIT) = BL.solve_bilevel_reduction(design, scope;
    cost_gap_pct=cap, objective=:lines, radial_mode=:warm, time_limit=tl,
    solver_threads=THREADS, mipgap=1e-4, t1_screening=true, cycle_cut_lens=(),
    solver_seed=0, output_flag=0, start_time_limit=60.0)

greedy_call(cap; tol=1e-9) = Main.Greedy.greedy_reduction(Main.Audit, BL, design, scope;
    time_limit=TIME_LIMIT, cost_gap_pct=cap, flow_tolerance=tol, cost_tolerance=1e-8,
    ordering=:flow, norm=:headroom, alpha=0.5, radial_first=true,
    kkt_check=true, threads=THREADS)

proxy_call(eps; near=0.8) = TR.solve_reduction_edge_multiscenario(design, eps .* base.frate;
    scenario_indices=scope, protection_indices=scope,
    near_limit_threshold=near, time_limit=TIME_LIMIT, mipgap=1e-3, threads=THREADS,
    numeric_focus=3, cycle_cut_lens=(), internal_rating_bound=true,
    switch_form=:hull, lmp_separation=false, merge_leaf_blocks=true)

println("\nwarming up (compiling every path; these timings are discarded)")
try; kkt_call(3.0; tl=20); catch e; println("  kkt warmup: ", e); end
try; greedy_call(3.0); catch e; println("  greedy warmup: ", e); end
try; proxy_call(0.25); catch e; println("  proxy warmup: ", e); end

"0/1 merge mask, whatever the approach calls it."
mask_of(r) = hasproperty(r, :internal) ? Int.(r.internal) :
             hasproperty(r, :c) ? round.(Int, r.c) : Int[]

rows = NamedTuple[]
function bench(approach, dir, setting, f)
    t = time(); r = f(); el = time() - t
    buses = hasproperty(r, :n_retained) ? r.n_retained :
            hasproperty(r, :buses) ? r.buses : -1
    st = hasproperty(r, :status) ? string(r.status) :
         hasproperty(r, :reason) ? string(r.reason) : ""
    # a run that stopped on the clock is not a solve time and must not read as one
    hit = el >= TIME_LIMIT - 1 || occursin("TIME_LIMIT", st)
    out = joinpath(ROOT, "outputs", dir, string(caseid), setting)
    mkpath(out)
    writedlm(joinpath(out, "internal.csv"), mask_of(r), ',')
    push!(rows, (approach=approach, setting=setting, buses=buses,
                 compute_s=round(el; digits=3), hit_time_limit=hit))
    @printf("%-8s %-18s %6d buses %10.3f s%s\n", approach, setting, buses, el,
            hit ? "   TIME LIMIT (not proved optimal)" : "")
    return nothing
end

println("\ntimed (case already loaded, code already compiled)")
@printf("%-8s %-18s %13s %10s\n", "approach", "setting", "buses", "compute s")
println("-"^52)
for cap in (1.0, 2.0, 3.0)
    bench("KKT", "kkt", "cap$(Int(cap))pct", () -> kkt_call(cap))
end
tag(x) = replace(@sprintf("%.2f", x), "." => "p")
for eps in (0.10, 0.25, 0.50)
    bench("Proxy", "proxy", "eps$(tag(eps))", () -> proxy_call(eps))
end
for nl in (0.95, 0.99)
    bench("Proxy", "proxy", "eps0p25_nl$(tag(nl))", () -> proxy_call(0.25; near=nl))
end
for cap in (1.0, 2.0, 3.0)
    bench("Greedy", "greedy", "cap$(Int(cap))pct", () -> greedy_call(cap))
end
bench("Greedy", "greedy", "cap3pct_tol1e-4", () -> greedy_call(3.0; tol=1e-4))

OUT = joinpath(ROOT, "outputs", casename * "_review")
mkpath(OUT)
hdr = collect(string.(keys(rows[1])))
writedlm(joinpath(OUT, "compute_time.csv"),
         vcat(permutedims(hdr), [getfield(r, Symbol(c)) for r in rows, c in hdr]), ',')
println("\nwrote compute_time.csv and per-setting internal.csv -> ", OUT)

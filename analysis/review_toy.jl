# Train/test review of the three reduction approaches on a toy case.
#
#   julia --project=. --startup-file=no analysis/review_toy.jl case=6busww month=9,10,11 n_train=20
#
# `month` takes one or several months. The pipeline's seed selection ranks every
# hour of the WHOLE window by total |line flow| and takes the n_train highest as
# the TRAINING set; every other hour is the TEST set, which no approach ever
# saw. The flow ranks of the chosen hours are printed so that is checkable.
#
# For each clustering it solves the reduced DC-OPF at each hour, puts that
# dispatch back on the full network, and measures the overload. This is the
# dispatch a solver actually RETURNS -- not the optimistic "some optimum is
# safe" reading the certificates use -- because that is what you get if you
# hand the reduced network to a downstream OPF.
#
# Everything is reported as a percentage: overload as % of the line's rating,
# cost as % over the full network's own optimum, violations as % of hours.

using JuMP, Gurobi, Printf, DelimitedFiles, Statistics, Dates

const ROOT = dirname(@__DIR__)
# Overloads are in % of rating. The repo's audit calls a line safe within
# flow_tolerance = 1e-6 of its rating, which is 1e-4 % here. Anything tighter
# sits below the LP solver's own feasibility tolerance and counts numerical
# noise as a violation -- at 1e-7 the proxy's case118 clustering 'violated'
# 60% of scenarios by margins that round to 0.000%.
const TOL = 1e-4

casename = "6busww"; months = [11]; n_train = 20; demands = :hourly
for a in ARGS
    k, v = split(a, Char(61); limit=2)
    if k == "case"; global casename = String(v)
    elseif k == "month"; global months = [parse(Int, x) for x in split(v, Char(44))]
    elseif k == "n_train"; global n_train = parse(Int, v)
    elseif k == "demands"; global demands = Symbol(v)
    elseif k == "all_settings"; nothing          # handled below, after ROOT setup
    else; error("accepts case=, month=, n_train=, demands=, all_settings=")
    end
end
const CASEID = Dict("6busww" => :case6ww, "14bus" => :case14toy,
                    "case118" => :case118)
caseid = CASEID[casename]

SETTINGS = (case=caseid, demands=demands, n_demands=n_train, scale_range=(0.9, 1.1),
            month=months, horizon_days=32*length(months), horizon_start_day=1,
            opf_time_limit=60.0, linear_costs=true)
include(joinpath(ROOT, "common", "settings.jl"))
S = merge(Main.Settings.apply_overrides(SETTINGS, String[]), (linear_costs=true,))

@eval module Data
    include(joinpath(dirname(@__DIR__), "common", "preprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "postprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "multiscenario.jl"))
end
include(joinpath(ROOT, "common", "cases.jl"))

design, heldout = Main.Cases.load_demands(Main.Data, S)
base = design.base
N, L, K = base.N, base.Ln, length(base.gen_bus)
ntr, nte = size(design.load, 2), isnothing(heldout) ? 0 : size(heldout.load, 2)
@printf("[%s] %d buses, %d lines, %d gens   months %s: %d train / %d test hours\n",
        casename, N, L, K, join(months, ","), ntr, nte)

# Check that the training set really is the high-flow end of the whole window.
# select_seed_scenarios ranks by total |line flow| and then unions in a minimum
# cover of the congested lines, so a couple of picks can sit outside the top n.
if nte > 0 && demands === :hourly
    flowscore = vcat([sum(abs, view(design.fhat, :, s)) for s in 1:ntr],
                     [sum(abs, view(heldout.fhat, :, s)) for s in 1:nte])
    rank = sortperm(sortperm(flowscore; rev=true))        # 0-based rank
    trrank = sort(rank[1:ntr] .+ 1)
    permonth = Dict{Int,Int}()
    for s in 1:ntr
        m = Dates.month(DateTime(2015,1,1) + Dates.Hour(design.scenario_ids[s] - 1))
        permonth[m] = get(permonth, m, 0) + 1
    end
    @printf("training hours: flow ranks %s ... %s of %d   (%d in the top %d)\n",
            join(trrank[1:min(5,end)], ","), string(trrank[end]), ntr + nte,
            count(<=(ntr), trrank), ntr)
    @printf("                by month: %s\n",
            join(["$m:$(permonth[m])" for m in sort(collect(keys(permonth)))], "  "))
end

env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
pmin = min.(0.0, base.pmin)

function reduced_model(internal)
    m = Model(() -> Gurobi.Optimizer(env)); set_silent(m)
    @variable(m, g[k=1:K]); @variable(m, th[1:N])
    @variable(m, f[1:L]); @variable(m, t[1:L])
    @variable(m, shed[1:N] >= 0); @variable(m, curt[1:N] >= 0)
    @constraint(m, [k=1:K], pmin[k] <= g[k] <= base.pmax[k])
    @constraint(m, th[base.j0] == 0)
    @constraint(m, [l=1:L], f[l] == base.Dx[l] * (th[base.Efrom[l]] - th[base.Eto[l]]))
    for l in 1:L
        if internal[l]
            @constraint(m, f[l] == 0)
        else
            @constraint(m, -base.frate[l] <= f[l] <= base.frate[l])
            @constraint(m, t[l] == 0)
        end
    end
    @constraint(m, bal[b=1:N],
        sum(g[k] for k in 1:K if base.gen_bus[k] == b; init=0.0) + shed[b] - curt[b]
        - sum(f[l] + t[l] for l in 1:L if base.Efrom[l] == b; init=0.0)
        + sum(f[l] + t[l] for l in 1:L if base.Eto[l]   == b; init=0.0) == 0.0)
    @objective(m, Min, sum(base.c1[k] * g[k] for k in 1:K) + 1e5*sum(shed) + 1e4*sum(curt))
    return (m=m, g=g, bal=bal, shed=shed, curt=curt)
end

# the true network carrying a fixed dispatch
fm = Model(() -> Gurobi.Optimizer(env)); set_silent(fm)
@variable(fm, tth[1:N]); @variable(fm, tf[1:L])
@constraint(fm, tth[base.j0] == 0)
@constraint(fm, [l=1:L], tf[l] == base.Dx[l] * (tth[base.Efrom[l]] - tth[base.Eto[l]]))
@constraint(fm, tbal[b=1:N],
    - sum(tf[l] for l in 1:L if base.Efrom[l] == b; init=0.0)
    + sum(tf[l] for l in 1:L if base.Eto[l]   == b; init=0.0) == 0.0)
@objective(fm, Min, 0.0)

"-> (max overload % of rating, cost, solved) per column of `loads`"
function run_scope(R, loads)
    S = size(loads, 2)
    ovl = zeros(S); cost = zeros(S); ok = trues(S)
    for s in 1:S
        for b in 1:N
            set_normalized_rhs(R.bal[b], loads[b, s])
            set_upper_bound(R.shed[b], max(loads[b, s], 0.0))
            set_upper_bound(R.curt[b], max(-loads[b, s], 0.0) + sum(base.pmax))
        end
        optimize!(R.m)
        if termination_status(R.m) != MOI.OPTIMAL
            ok[s] = false; ovl[s] = Inf; continue
        end
        gv, sv, cv = value.(R.g), value.(R.shed), value.(R.curt)
        cost[s] = sum(base.c1[k] * gv[k] for k in 1:K)
        for b in 1:N
            inj = sum(gv[k] for k in 1:K if base.gen_bus[k] == b; init=0.0) +
                  sv[b] - cv[b] - loads[b, s]
            set_normalized_rhs(tbal[b], -inj)
        end
        optimize!(fm)
        if termination_status(fm) != MOI.OPTIMAL
            ok[s] = false; ovl[s] = Inf; continue
        end
        tfv = value.(tf)
        ovl[s] = maximum(100 * max(abs(tfv[l]) - base.frate[l], 0.0) / base.frate[l]
                         for l in 1:L)
    end
    return ovl, cost, ok
end

# baseline: the full network's own optimum, for the cost gap
B = reduced_model(falses(L))
ovl0_tr, cost0_tr, _ = run_scope(B, design.load)
ovl0_te, cost0_te, _ = nte > 0 ? run_scope(B, heldout.load) : (Float64[], Float64[], Bool[])

clusters(internal) = begin
    p = collect(1:N); rt(x) = p[x] == x ? x : (p[x] = rt(p[x]))
    for l in findall(internal); p[rt(base.Efrom[l])] = rt(base.Eto[l]); end
    length(unique(rt(b) for b in 1:N))
end

"""Newest clustering written by `approach` for this case, as a Bool mask.

kkt/ and greedy/ go through Results.finish and write internal.csv. The proxy has
its own reporting and writes line_status.csv one level deeper, with the mask in
an `internal` column -- same information, different file.
"""
EVERY = any(a -> a == "all_settings=true", ARGS)

"""Every clustering `approach` has written for this case: [(label, path)].

The label is the run directory, which is what carries the setting -- so the
runners must be given an explicit output_dir when sweeping a knob that their
own run_tag ignores (cost_cap for kkt and greedy; the proxy's tag has eps).
Sorted by directory name so the table order is stable across runs.
"""
function all_clusterings(approach)
    root = joinpath(ROOT, "outputs", approach, string(caseid))
    isdir(root) || return Tuple{String,String}[]
    out = Tuple{String,String}[]
    for d in sort(readdir(root; join=true))
        isdir(d) || continue
        f = joinpath(d, "internal.csv")
        isfile(f) && push!(out, (basename(d), f))
        for sub in sort(readdir(d; join=true))     # proxy: <tag>/<setting>/
            isdir(sub) || continue
            g = joinpath(sub, "line_status.csv")
            isfile(g) && push!(out, (basename(d), g))
        end
    end
    # default to the newest run only; all_settings=true reports the whole sweep
    (EVERY || isempty(out)) && return out
    newest = argmax(p -> mtime(p[2]), out)
    return [newest]
end

function read_mask(path)
    if endswith(path, "line_status.csv")
        tab, hd = readdlm(path, ',', Any, '\n'; header=true)
        col = findfirst(==("internal"), vec(hd))
        isnothing(col) && error("no `internal` column in $path")
        return Int.(tab[:, col]) .== 1
    end
    return vec(readdlm(path, ',', Int)) .== 1
end

"""Measured compute time, keyed by (approach, setting).

Written by analysis/time_approaches.jl, which loads the case once, calls every
approach once to compile, throws those timings away, and only then times each
setting. That is the only fair number here: a runner is a fresh process calling
its algorithm exactly once, so its own reported seconds are mostly Julia JIT --
greedy's 5.8 s is 5.15 s of compilation and 0.09 s of work -- while the proxy
reports JuMP's solver time alone and never includes JIT at all.
"""
const COMPUTE = let f = joinpath(ROOT, "outputs", casename * "_review", "compute_time.csv")
    d = Dict{Tuple{String,String},Tuple{Float64,Bool}}()
    if isfile(f)
        tab, hd = readdlm(f, ',', Any, '\n'; header=true)
        ia = findfirst(==("approach"), vec(hd))
        is = findfirst(==("setting"), vec(hd))
        ic = findfirst(==("compute_s"), vec(hd))
        il = findfirst(==("hit_time_limit"), vec(hd))
        if !any(isnothing, (ia, is, ic))
            for r in axes(tab, 1)
                hit = isnothing(il) ? false : (string(tab[r, il]) in ("true", "1"))
                d[(strip(String(tab[r, ia])), strip(String(tab[r, is])))] =
                    (Float64(tab[r, ic]), hit)
            end
        end
    else
        @warn "no compute_time.csv; run analysis/time_approaches.jl case=$casename"
    end
    d
end

compute_of(approach, setting) = get(COMPUTE, (approach, setting), (NaN, false))

rows = NamedTuple[]
function score(label, internal, srcnote, t=(NaN, false))
    R = reduced_model(internal)
    otr, ctr, oktr = run_scope(R, design.load)
    ote, cte, okte = nte > 0 ? run_scope(R, heldout.load) : (Float64[], Float64[], Bool[])
    gap(c, c0) = isempty(c) ? 0.0 :
        maximum(c0[i] > 1e-9 ? 100*(c[i]-c0[i])/c0[i] : 0.0 for i in eachindex(c))
    vtr = count(i -> !oktr[i] || otr[i] > TOL, 1:ntr)
    vte = nte == 0 ? 0 : count(i -> !okte[i] || ote[i] > TOL, 1:nte)
    fin(v) = isempty(v) ? 0.0 : maximum(x -> isfinite(x) ? x : 0.0, v)
    push!(rows, (approach=label, setting=srcnote, buses=clusters(internal),
        reduction_pct = round(100*(N - clusters(internal))/N; digits=2),
        train_violating_pct = round(100*vtr/max(ntr,1); digits=2),
        test_violating_pct  = round(100*vte/max(nte,1); digits=2),
        train_max_overload_pct = round(fin(otr); digits=5),
        test_max_overload_pct  = round(fin(ote); digits=5),
        train_max_cost_gap_pct = round(gap(ctr, cost0_tr); digits=4),
        test_max_cost_gap_pct  = round(gap(cte, cost0_te); digits=4),
        compute_s = isnan(t[1]) ? "" : round(t[1]; digits=3),
        time_limit_hit = t[2] ? "yes" : ""))
end

score("Full network (no reduction)", falses(L), "-")
for (lab, dir) in (("KKT", "kkt"), ("Proxy", "proxy"), ("Greedy", "greedy"))
    found = all_clusterings(dir)
    if isempty(found)
        @printf("  %-6s nothing under outputs/%s/%s -- run it first\n", lab, dir, caseid)
        continue
    end
    for (setting, f) in found
        v = read_mask(f)
        length(v) == L || error("$f has $(length(v)) lines, case has $L")
        score(lab, v, setting, compute_of(lab, setting))
    end
end

OUT = joinpath(ROOT, "outputs", casename * "_review")
mkpath(OUT)
hdr = collect(string.(keys(rows[1])))
writedlm(joinpath(OUT, "summary.csv"),
         vcat(permutedims(hdr), [getfield(r, Symbol(c)) for r in rows, c in hdr]), ',')

println()
@printf("%-8s %-22s %6s %8s | %8s %8s | %8s %8s | %8s %8s | %10s\n",
        "approach", "setting", "buses", "reduc %", "train v%", "test v%",
        "tr ovl%", "te ovl%", "tr cost%", "te cost%", "compute s")
println("-"^140)
for r in rows
    cv = r.compute_s == "" ? "-" :
         (r.compute_s < 1 ? @sprintf("%.3f", r.compute_s) : @sprintf("%.1f", r.compute_s)) *
         (r.time_limit_hit == "yes" ? "*" : "")
    @printf("%-8s %-22s %6d %8.2f | %8.2f %8.2f | %8.3f %8.3f | %8.4f %8.4f | %10s\n",
            r.approach, first(r.setting, 22), r.buses, r.reduction_pct,
            r.train_violating_pct, r.test_violating_pct,
            r.train_max_overload_pct, r.test_max_overload_pct,
            r.train_max_cost_gap_pct, r.test_max_cost_gap_pct, cv)
end
@printf("\nwrote summary.csv -> %s\n", OUT)

# Do the networks we have survive losing one generator? Every hour of the March
# week, with each generator out in turn (re-dispatched after the outage, on the
# full network and on the reduced one alike).
#
#   julia --project=. --startup-file=no analysis/contingency_profile.jl
#
# Per network: deliverable = some cheapest reduced dispatch works on the full
# network; within cap = deliverable and cost within the cap it was built with.
# Writes outputs/contingency/ACTIVSg200/profile.csv and failures.csv.

using Printf, DelimitedFiles, JuMP
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "common", "settings.jl"))
S = (case=:ACTIVSg200, demands=:hourly, n_demands=168, month=3, horizon_days=7,
     horizon_start_day=1, opf_time_limit=60.0, linear_costs=true)
@eval module Data
    include(joinpath(dirname(@__DIR__), "common", "preprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "postprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "multiscenario.jl"))
end
include(joinpath(ROOT, "common", "cases.jl"))
include(joinpath(ROOT, "common", "audit.jl"))
include(joinpath(ROOT, "common", "results.jl"))
const AU = Main.Audit

OUT = joinpath(ROOT, "outputs", "contingency", "ACTIVSg200")
mkpath(OUT)
c, _ = Main.Cases.load_demands(Main.Data, S)
base, D = c.base, c.load
K, T = length(base.gen_bus), size(D, 2)
pmin = min.(0.0, base.pmin)
gens = [k for k in 1:K if base.pmax[k] > 1e-9]
@printf("ACTIVSg200: %d buses, %d generators (%d with capacity), %d hours\n", base.N, K, length(gens), T)

# scenario (t, k): hour t with generator k out; k = 0 is no outage
function outage!(oracle, k, on)
    k == 0 && return
    for part in (oracle.lower, oracle.safety)
        set_upper_bound(part.g[k], on ? 0.0 : base.pmax[k])
        set_lower_bound(part.g[k], on ? 0.0 : pmin[k])
    end
end
function sweep(oracle)
    res = Matrix{Any}(undef, T, length(gens) + 1)
    for (j, k) in enumerate(vcat(0, gens))
        outage!(oracle, k, true)
        for t in 1:T
            res[t, j] = AU.check_scenario!(oracle, D[:, t])
        end
        outage!(oracle, k, false)
    end
    return res
end

t0 = time()
full = sweep(AU.topology_oracle(base, falses(base.Ln)))
fullcost = map(a -> a.status == MOI.OPTIMAL ? a.reduced_cost : NaN, full)
infeasible = count(isnan, fullcost)
@printf("full network: %d of %d outage scenarios infeasible (not enough capacity) (%.0fs)\n",
        infeasible, T * length(gens), time() - t0)

readmask(dir) = vec(readdlm(joinpath(dir, "internal.csv"), ',', Int)) .== 1
camp = joinpath(ROOT, "outputs", "campaign", "ACTIVSg200")
step2 = joinpath(ROOT, "outputs", "step2", "ACTIVSg200")
nets = [("greedy week", 0.1, joinpath(step2, "greedy_week_cap0.1")),
        ("greedy week", 0.5, joinpath(step2, "greedy_week_cap0.5")),
        ("greedy week", 1.0, joinpath(step2, "greedy_week_cap1")),
        ("greedy 10h", 0.1, joinpath(camp, "greedy_noladder_cap0.1")),
        ("greedy 10h", 0.5, joinpath(camp, "greedy_noladder_cap0.5")),
        ("kkt 10h", 0.1, joinpath(camp, "kkt_noladder_cap0.1")),
        ("proxy 10h eps1%", NaN, joinpath(camp, "proxy_noladder_eps1"))]

rows, fails = NamedTuple[], NamedTuple[]
for (name, cap, dir) in nets
    isfile(joinpath(dir, "internal.csv")) || (println("missing ", dir); continue)
    internal = readmask(dir)
    oracle = AU.topology_oracle(base, internal)
    t1 = time()
    res = sweep(oracle)
    for (label, cols) in (("no outage", 1:1), ("one generator out", 2:size(res, 2)))
        idx = [(t, j) for t in 1:T for j in cols if isfinite(fullcost[t, j])]
        ok = [res[t, j].classification == :pass for (t, j) in idx]
        gap = [100 * (res[t, j].reduced_cost - fullcost[t, j]) / abs(fullcost[t, j]) for (t, j) in idx]
        incap = [ok[i] && (isnan(cap) || gap[i] <= cap + 1e-6) for i in eachindex(idx)]
        exc = [res[t, j].excess for (t, j) in idx]
        r = (network="$name cap $(cap)%", buses=oracle.R, scenarios=label, total=length(idx),
             deliverable=count(ok), within_cap=count(incap),
             reduced_infeasible=count(i -> res[idx[i]...].classification == :reduced_infeasible, eachindex(idx)),
             worst_overload_pct=100 * maximum(e -> isnan(e) ? 0.0 : e, exc; init=0.0),
             worst_cost_pct=maximum(filter(isfinite, gap); init=-Inf))
        push!(rows, r)
        @printf("%-22s %4d buses  %-18s deliverable %5d  within cap %5d / %5d  worst overload %6.2f%%  worst cost %7.3f%%\n",
                r.network, r.buses, label, r.deliverable, r.within_cap, r.total, r.worst_overload_pct, r.worst_cost_pct)
        label == "no outage" && continue
        for (j, k) in enumerate(gens)
            bad = count(t -> isfinite(fullcost[t, j+1]) && res[t, j+1].classification != :pass, 1:T)
            bad > 0 && push!(fails, (network=r.network, generator=k, bus=base.gen_bus[k],
                                     pmax=base.pmax[k], failed_hours=bad))
        end
    end
    worst = sort([f for f in fails if f.network == "$name cap $(cap)%"]; by=f -> -f.failed_hours)
    isempty(worst) || println("    outages that break it most: ",
        join((@sprintf("gen %d (%.0f MW) %d h", f.generator, 100 * f.pmax, f.failed_hours) for f in first(worst, 5)), ", "))
    @printf("    (%.0fs)\n", time() - t1)
    flush(stdout)
    Main.Results.write_rows(joinpath(OUT, "profile.csv"), rows)
    isempty(fails) || Main.Results.write_rows(joinpath(OUT, "failures.csv"), fails)
end
println("-> ", OUT)

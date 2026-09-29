# Step 2: greedy built on all 168 hours of a week against the certified reduction
# built on the same week, both tested on that week and on three unseen weeks.
#
#   julia --project=. --startup-file=no analysis/step2_week_compare.jl
#
# Networks
#   greedy week    greedy/run_greedy.jl with every hour of Mar 1-7 as a design demand
#                  (outputs/step2/ACTIVSg200/greedy_week_cap*/internal.csv)
#   greedy week +  the same networks, certified afterwards (ratings from the
#   certificate    exact check over the week's demand hull)
#   certified      sensitivity clustering + ratings, over the same hull, beta sweep
#   reference      the campaign's 10-hour greedy and KKT networks
# Every hour of every week: deliverable (some cheapest reduced dispatch works on
# the full network), safe (the dispatch the solver returned works), cost vs full.
# Writes outputs/step2/ACTIVSg200/comparison.csv.

using Printf, DelimitedFiles, LinearAlgebra
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "common", "settings.jl"))
week(month, day) = (case=:ACTIVSg200, demands=:hourly, n_demands=168, month=month,
                    horizon_days=7, horizon_start_day=day, opf_time_limit=60.0, linear_costs=true)
S = week(3, 1)
@eval module Data
    include(joinpath(dirname(@__DIR__), "common", "preprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "postprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "multiscenario.jl"))
end
include(joinpath(ROOT, "common", "cases.jl"))
include(joinpath(ROOT, "common", "caps.jl"))
include(joinpath(ROOT, "common", "audit.jl"))
include(joinpath(ROOT, "common", "results.jl"))
include(joinpath(ROOT, "kkt", "kkt_model.jl"))
include(joinpath(ROOT, "certify", "certificate.jl"))
include(joinpath(ROOT, "certify", "coarsen.jl"))
const BL, CT, CO, AU = Main.BilevelReduction, Main.Certificate, Main.Coarsen, Main.Audit

OUT = joinpath(ROOT, "outputs", "step2", "ACTIVSg200")
mkpath(OUT)
weeks = [("Mar 1-7 (built on)", week(3, 1)), ("Mar 8-14", week(3, 8)),
         ("Jul 17-23", week(7, 17)), ("Jan 16-22", week(1, 16))]
data = map(weeks) do (label, s)
    c, _ = Main.Cases.load_demands(Main.Data, s)
    D = c.load
    (label=label, D=D, full=[BL.true_optimal_cost(c.base, D[:, t]; relax_pmin=true) for t in axes(D, 2)],
     base=c.base)
end
base = data[1].base
for d in data
    d.base.frate == base.frate && d.base.Dx == base.Dx || error("weeks load different networks")
end
@printf("ACTIVSg200: %d buses; weeks: %s\n", base.N, join(("$(d.label) $(size(d.D, 2)) h" for d in data), ", "))
flush(stdout)

H = CT.ptdf(base)
set = CT.injection_set(base, data[1].D)
crit = findall(>(1 + 1e-7), CT.reach(base, set; H))
@printf("lines that can overload over the Mar 1-7 hull: %s\n", crit)

function evaluate(internal, ratings, d)
    oracle = AU.topology_oracle(base, internal; reduced_ratings=ratings)
    ev = [let a = AU.check_scenario!(oracle, d.D[:, t])
              (ok=a.classification == :pass, safe=a.returned_utilization <= 1 + 1e-6,
               excess=a.excess, infeasible=a.classification == :reduced_infeasible,
               cost=100 * (a.reduced_cost - d.full[t]) / abs(d.full[t]))
          end for t in axes(d.D, 2)]
    fin(x) = filter(isfinite, x)
    (ok=count(e -> e.ok, ev), safe=count(e -> e.safe, ev), hours=length(ev),
     infeasible=count(e -> e.infeasible, ev),
     overload=100 * maximum(e -> isnan(e.excess) ? Inf : e.excess, ev),
     mean_cost=(x = fin([e.cost for e in ev]); isempty(x) ? NaN : sum(x) / length(x)),
     worst_cost=maximum(fin([e.cost for e in ev]); init=-Inf))
end

readmask(dir) = vec(readdlm(joinpath(dir, "internal.csv"), ',', Int)) .== 1
build_seconds(dir) = (hdr = split(readlines(joinpath(dir, "summary.csv"))[1], ',');
                      val = split(readlines(joinpath(dir, "summary.csv"))[2], ',');
                      parse(Float64, val[findfirst(==("seconds"), hdr)]))

nets = NamedTuple[]
for cap in ("0.1", "0.5", "1")
    dir = joinpath(OUT, "greedy_week_cap$cap")
    isfile(joinpath(dir, "internal.csv")) || (println("missing ", dir); continue)
    m = readmask(dir)
    push!(nets, (method="greedy week", name="cap $cap%", internal=m, ratings=nothing,
                 seconds=build_seconds(dir)))
    t0 = time()
    dr = CT.design_ratings(base, set, m; H, lines=crit)
    push!(nets, (method="greedy week + certificate", name="cap $cap%", internal=m,
                 ratings=dr.certified ? dr.rating : missing, seconds=build_seconds(dir) + time() - t0))
end
for beta in (0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0)
    t0 = time()
    m = CO.coarsen(base, crit, set, beta)
    dr = CT.design_ratings(base, set, m; H, lines=crit)
    push!(nets, (method="certified", name="beta $beta", internal=m,
                 ratings=dr.certified ? dr.rating : missing, seconds=time() - t0))
end
camp = joinpath(ROOT, "outputs", "campaign", "ACTIVSg200")
for (run, label) in (("greedy_noladder_cap0.1", "greedy 10h cap 0.1%"),
                     ("greedy_noladder_cap0.5", "greedy 10h cap 0.5%"),
                     ("kkt_noladder_cap0.1", "kkt 10h cap 0.1%"))
    isdir(joinpath(camp, run)) || continue
    push!(nets, (method="reference", name=label, internal=readmask(joinpath(camp, run)),
                 ratings=nothing, seconds=NaN))
end

rows = NamedTuple[]
for n in nets
    buses = CT.merged_ptdf(base, n.internal).K
    if ismissing(n.ratings)
        @printf("%-26s %-20s %4d buses  not certifiable\n", n.method, n.name, buses)
        push!(rows, (method=n.method, network=n.name, buses=buses, certified=false,
                     build_seconds=round(n.seconds; digits=1), week="", ok=-1, safe=-1, hours=0,
                     infeasible=-1, worst_overload_pct=NaN, mean_cost_pct=NaN, worst_cost_pct=NaN))
        continue
    end
    for d in data
        r = evaluate(n.internal, n.ratings, d)
        push!(rows, (method=n.method, network=n.name, buses=buses,
                     certified=n.method in ("certified", "greedy week + certificate"),
                     build_seconds=round(n.seconds; digits=1), week=d.label, ok=r.ok, safe=r.safe,
                     hours=r.hours, infeasible=r.infeasible, worst_overload_pct=r.overload,
                     mean_cost_pct=r.mean_cost, worst_cost_pct=r.worst_cost))
        @printf("%-26s %-20s %4d buses  %-19s deliverable %3d  safe %3d /%d  worst overload %8.4f%%  cost %7.3f%% (worst %7.3f%%)\n",
                n.method, n.name, buses, d.label, r.ok, r.safe, r.hours, r.overload,
                r.mean_cost, r.worst_cost)
    end
    flush(stdout)
    Main.Results.write_rows(joinpath(OUT, "comparison.csv"), rows)
end
println("-> ", joinpath(OUT, "comparison.csv"))

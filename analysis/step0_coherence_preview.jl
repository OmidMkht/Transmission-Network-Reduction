# Preview for step 1: cluster ACTIVSg200 by sensitivity to the lines that can
# overload, then certify with the same machinery as step0_certificate.jl.
#
#   julia --project=. --startup-file=no analysis/step0_coherence_preview.jl
#
# A line is merged when its two ends have PTDF entries within delta on every
# critical line; critical lines themselves are never merged. Buses with equal
# sensitivity to line l can merge without changing l's flow for any injection.

using Printf, DelimitedFiles, LinearAlgebra
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "common", "settings.jl"))
S = (case=:ACTIVSg200, demands=:hourly, n_demands=10, month=3, horizon_days=7,
     horizon_start_day=1, opf_time_limit=60.0, linear_costs=true)
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
const BL, CT, AU = Main.BilevelReduction, Main.Certificate, Main.Audit

OUT = joinpath(ROOT, "outputs", "step0_ACTIVSg200")
design, heldout = Main.Cases.load_demands(Main.Data, S)
base = design.base
Dall = hcat(design.load, heldout.load)
nd = size(design.load, 2)
H = CT.ptdf(base)
set = CT.injection_set(base, Dall)

crit = Int[]
for l in 1:base.Ln
    lo, hi = CT.extremes(set, H[l, :])
    max(hi, -lo) > base.frate[l] * (1 + 1e-7) && push!(crit, l)
end
println("lines that can overload over the week: ", crit)

full_cost = [BL.true_optimal_cost(base, Dall[:, s]; relax_pmin=true) for s in 1:size(Dall, 2)]
function evaluate(internal, ratings)
    oracle = AU.topology_oracle(base, internal; reduced_ratings=ratings)
    [let a = AU.check_scenario!(oracle, Dall[:, s])
         (ok=a.classification == :pass, util=a.returned_utilization,
          cost=100 * (a.reduced_cost - full_cost[s]) / abs(full_cost[s]))
     end for s in 1:size(Dall, 2)]
end
closure(mask) = (cl = CT.clusters(base, mask);
                 BitVector([cl[base.Efrom[l]] == cl[base.Eto[l]] for l in 1:base.Ln]))

rows = NamedTuple[]
for delta in (1e-6, 1e-3, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
    t0 = time()
    seed = BitVector([!(l in crit) &&
                      maximum(abs(H[k, base.Efrom[l]] - H[k, base.Eto[l]]) for k in crit) <= delta
                      for l in 1:base.Ln])
    internal = closure(seed)
    cert = CT.certify(base, internal, set; H)
    mind = CT.minimal_derating(base, set, cert; H)
    secs = time() - t0
    orig, cer = evaluate(internal, nothing), evaluate(internal, mind.rating)
    extra = [cer[s].cost - orig[s].cost for s in eachindex(cer)]
    der = isempty(crit) ? 0.0 :
        maximum(100 * (base.frate[l] - mind.rating[l]) / base.frate[l] for l in crit)
    r = (delta=delta, buses=cert.mp.K, merged_lines=count(internal),
         critical_merged=count(l -> internal[l], crit), certified=isfinite(mind.alpha),
         max_derating_pct=der, orig_ok=count(e -> e.ok, orig), cert_ok=count(e -> e.ok, cer),
         cert_max_util=maximum(e -> e.util, cer), orig_mean_cost_pct=sum(e -> e.cost, orig) / length(orig),
         cert_mean_extra_cost_pct=sum(extra) / length(extra), cert_worst_cost_pct=maximum(e -> e.cost, cer),
         seconds=round(secs; digits=2))
    push!(rows, r)
    @printf("delta %-6g %3d buses  crit merged %d  certified %-5s derate %5.1f%%  deliverable %3d -> %3d of 168  util %.4f  extra cost %.3f%%  (%.1fs)\n",
            delta, r.buses, r.critical_merged, r.certified, r.max_derating_pct, r.orig_ok,
            r.cert_ok, r.cert_max_util, r.cert_mean_extra_cost_pct, secs)
end
Main.Results.write_rows(joinpath(OUT, "coherence_preview.csv"), rows)

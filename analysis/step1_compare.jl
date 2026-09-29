# Step 1: certified reduction (sensitivity clustering + dual-based ratings) against
# the campaign's greedy / KKT / proxy networks, all certified the same way.
#
#   julia --project=. --startup-file=no analysis/step1_compare.jl ACTIVSg200
#   julia --project=. --startup-file=no analysis/step1_compare.jl case300
#
# Demands: ACTIVSg200 uses the campaign week (10 design + 158 held-out hours);
# the pglib cases use loads scaled 0.9 / 1.0 / 1.1 (0.95 and 1.05 held out).
# Every network is certified over a demand hull with any dispatch in the
# generator box, and then checked hour by hour with its certified ratings.
# Writes outputs/step1/<case>.csv.

using Printf, DelimitedFiles, LinearAlgebra
const ROOT = dirname(@__DIR__)
CASE = Symbol(get(ARGS, 1, "ACTIVSg200"))
include(joinpath(ROOT, "common", "settings.jl"))
S = CASE === :ACTIVSg200 ?
    (case=CASE, demands=:hourly, n_demands=10, month=3, horizon_days=7, horizon_start_day=1,
     opf_time_limit=60.0, linear_costs=true) :
    (case=CASE, demands=:scaled, n_demands=3, scale_range=(0.9, 1.1),
     opf_time_limit=60.0, linear_costs=true)
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
include(joinpath(ROOT, "certify", "certificate.jl"))
include(joinpath(ROOT, "certify", "coarsen.jl"))
const BL, CT, CO, AU = Main.BilevelReduction, Main.Certificate, Main.Coarsen, Main.Audit

OUT = joinpath(ROOT, "outputs", "step1")
mkpath(OUT)
design, heldout = Main.Cases.load_demands(Main.Data, S)
base = design.base
Dd, Dh = design.load, heldout.load
Dall = hcat(Dd, Dh)
nd, nt = size(Dd, 2), size(Dall, 2)
H = CT.ptdf(base)
full_cost = [BL.true_optimal_cost(base, Dall[:, s]; relax_pmin=true) for s in 1:nt]
@printf("%s: %d buses, %d lines, %d generators, %d design + %d held-out demands\n",
        CASE, base.N, base.Ln, length(base.gen_bus), nd, nt - nd)
flush(stdout)

function evaluate(internal, ratings)
    oracle = AU.topology_oracle(base, internal; reduced_ratings=ratings)
    [let a = AU.check_scenario!(oracle, Dall[:, s])
         (ok=a.classification == :pass, util=a.returned_utilization,
          infeasible=a.classification == :reduced_infeasible,
          cost=100 * (a.reduced_cost - full_cost[s]) / abs(full_cost[s]))
     end for s in 1:nt]
end
fin(x) = filter(isfinite, x)
mean_or_nan(x) = isempty(x) ? NaN : sum(x) / length(x)

campaign_dir = joinpath(ROOT, "outputs", "campaign", string(CASE))
campaign = isdir(campaign_dir) ?
    sort([d for d in readdir(campaign_dir) if isfile(joinpath(campaign_dir, d, "internal.csv"))]) : String[]
betas = (0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0)
# ACTIVSg200 is certified over the whole week and over the design hours alone;
# for the pglib cases the design hull already contains the held-out scalings.
variants = CASE === :ACTIVSg200 ? (("week", Dall), ("design", Dd)) : (("design", Dd),)

rows = NamedTuple[]
for (variant, D) in variants
    set = CT.injection_set(base, D)
    t0 = time()
    rch = CT.reach(base, set; H)
    crit = findall(>(1 + 1e-7), rch)
    t_reach = time() - t0
    @printf("\n[%s] %d of %d lines can overload over the set (%.1fs)\n", variant, length(crit), base.Ln, t_reach)
    nets = vcat([("coherent", @sprintf("beta=%g", b), b) for b in betas],
                [("campaign", r, NaN) for r in campaign])
    for (method, name, beta) in nets
        t1 = time()
        internal = method == "coherent" ? CO.coarsen(base, crit, set, beta) :
            vec(readdlm(joinpath(campaign_dir, name, "internal.csv"), ',', Int)) .== 1
        t_cluster = time() - t1
        mp = CT.merged_ptdf(base, internal)
        t2 = time()
        dr = CT.design_ratings(base, set, internal; H, lines=crit)
        t_cert = time() - t2
        if dr.certified && !isempty(crit)   # independent re-check of the claim
            chk = maximum(CT.max_loading(base, set, mp, dr.rating; H, lines=crit))
            chk <= 1 + 1e-5 || error("$name: certified but exact check gives $chk")
        end
        orig = evaluate(internal, nothing)
        cer = evaluate(internal, dr.rating)
        derated = [l for l in 1:base.Ln if dr.rating[l] < base.frate[l] * (1 - 1e-6)]
        r = (case=string(CASE), variant=variant, method=method, network=name, beta=beta,
             buses=mp.K, merged_lines=count(internal), critical_lines=length(crit),
             critical_merged=count(l -> internal[l], crit),
             certified=dr.certified, reason=dr.reason, checks=dr.checks,
             derated_lines=length(derated),
             max_derating_pct=isempty(derated) ? 0.0 :
                 maximum(100 * (base.frate[l] - dr.rating[l]) / base.frate[l] for l in derated),
             orig_ok=count(e -> e.ok, orig), cert_ok=count(e -> e.ok, cer), hours=nt,
             cert_infeasible=count(e -> e.infeasible, cer),
             cert_max_util=maximum(fin([e.util for e in cer]); init=-Inf),
             orig_mean_cost_pct=mean_or_nan(fin([e.cost for e in orig])),
             cert_mean_cost_pct=mean_or_nan(fin([e.cost for e in cer])),
             cert_worst_cost_pct=maximum(fin([e.cost for e in cer]); init=-Inf),
             cluster_seconds=round(t_cluster; digits=2), certify_seconds=round(t_cert; digits=2))
        push!(rows, r)
        @printf("  %-8s %-24s %4d buses  cert %-5s checks %3d  derated %2d (max %5.1f%%)  ok %3d -> %3d/%d  util %.4f  cost %7.3f%% -> %7.3f%% (worst %7.3f%%)  %.1fs\n",
                method, name, r.buses, r.certified, r.checks, r.derated_lines, r.max_derating_pct,
                r.orig_ok, r.cert_ok, nt, r.cert_max_util, r.orig_mean_cost_pct,
                r.cert_mean_cost_pct, r.cert_worst_cost_pct, t_cluster + t_cert)
        flush(stdout)
        Main.Results.write_rows(joinpath(OUT, "$(CASE).csv"), rows)
    end
end

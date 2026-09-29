# Step 0 of the certified-reduction idea, on the campaign's ACTIVSg200 networks.
#
#   julia --project=. --startup-file=no analysis/step0_certificate.jl
#
# For each network (greedy, kkt, proxy; outputs/campaign/ACTIVSg200/*/internal.csv)
# and two demand sets -- the hull of all 168 hours of the week, and the hull of
# the 10 design hours only:
#   1. certificate: which merged lines can overload, derated ratings F - e
#   2. exact check: worst true loading over every dispatch the reduced network
#      allows, before and after derating (an LP per line)
#   3. every hour: reduced OPF with original and with derated ratings --
#      deliverable?, overload, cost against the full OPF
# Writes outputs/step0_ACTIVSg200/{checks.txt, networks.csv, lines.csv, hours.csv}.

using Printf, DelimitedFiles, LinearAlgebra, Random
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

CAMPAIGN = joinpath(ROOT, "outputs", "campaign", "ACTIVSg200")
OUT = joinpath(ROOT, "outputs", "step0_ACTIVSg200")
mkpath(OUT)
log_io = open(joinpath(OUT, "checks.txt"), "w")
say(args...) = (println(args...); println(log_io, args...); flush(log_io))

design, heldout = Main.Cases.load_demands(Main.Data, S)
base = design.base
Dd, Dh = design.load, heldout.load
Dall = hcat(Dd, Dh)
nd, nh = size(Dd, 2), size(Dh, 2)
hour_ids = vcat(design.scenario_ids, heldout.scenario_ids)
pmin = min.(0.0, base.pmin)
say(@sprintf("ACTIVSg200: %d buses, %d lines, %d generators, %d design + %d held-out hours",
             base.N, base.Ln, length(base.gen_bus), nd, nh))
runs = sort([d for d in readdir(CAMPAIGN) if isfile(joinpath(CAMPAIGN, d, "internal.csv"))])
internal_of(run) = vec(readdlm(joinpath(CAMPAIGN, run, "internal.csv"), ',', Int)) .== 1
H = CT.ptdf(base)

# ---- check 1: merged-network sensitivities against the physics ---------------
# A merged line is a line of zero reactance. Give merged lines a huge susceptance
# in the FULL network and its flows on the other lines must match the merged
# network's, for the injections the full OPF actually produced.
let run = "kkt_noladder_cap1", internal = internal_of(run)
    mp = CT.merged_ptdf(base, internal)
    Dx2 = [internal[l] ? 1e7 * base.Dx[l] : base.Dx[l] for l in 1:base.Ln]
    base2 = typeof(base)((f === :Dx ? Dx2 : getfield(base, f) for f in fieldnames(typeof(base)))...)
    H2 = CT.ptdf(base2)
    worst = 0.0
    for s in 1:nd
        p = design.p[:, s]
        d = abs.((H2 * p .- mp.Hm * p)[mp.external]) ./ base.frate[mp.external]
        worst = max(worst, maximum(d))
    end
    say(@sprintf("check 1  merged PTDF vs 1e7-susceptance network (%s, %d clusters): worst gap %.2e of rating",
                 run, mp.K, worst))
    worst < 1e-4 || error("check 1 failed")
end

# ---- check 2: worst-case LP against T1's closed form -------------------------
let S1 = CT.injection_set(base, Dd[:, 1:1]), worst = 0.0
    total = sum(Dd[:, 1])
    for l in 1:base.Ln
        a = [H[l, b] for b in base.gen_bus]
        off = dot(H[l, :], Dd[:, 1])
        hi = BL._extreme_over_dispatch(a, pmin, base.pmax, total; maximize=true) - off
        lo = BL._extreme_over_dispatch(a, pmin, base.pmax, total; maximize=false) - off
        lo2, hi2 = CT.extremes(S1, H[l, :])
        worst = max(worst, abs(hi - hi2), abs(lo - lo2))
    end
    t1 = BL.t1_redundant_limits(base, [Dd[:, s] for s in 1:nd], pmin)
    reach = 0
    for s in 1:nd
        Ss = CT.injection_set(base, Dd[:, s:s])
        for l in 1:base.Ln
            lo, hi = CT.extremes(Ss, H[l, :])
            reach += max(hi, -lo) > base.frate[l] ? 1 : 0
        end
    end
    say(@sprintf("check 2  LP vs knapsack extremes on design hour 1: worst gap %.2e p.u.;  reachable line-hours: LP %d, T1 %d",
                 worst, reach, count(.!t1)))
    (worst < 1e-6 && reach == count(.!t1)) || error("check 2 failed")
end

# ---- the per-hour evaluation ---------------------------------------------------
full_cost = [BL.true_optimal_cost(base, Dall[:, s]; relax_pmin=true) for s in 1:size(Dall, 2)]
function evaluate(internal, ratings)
    oracle = AU.topology_oracle(base, internal; reduced_ratings=ratings)
    map(1:size(Dall, 2)) do s
        a = AU.check_scenario!(oracle, Dall[:, s])
        (ok=a.classification == :pass, excess=a.excess, util=a.returned_utilization,
         cost_pct=100 * (a.reduced_cost - full_cost[s]) / abs(full_cost[s]),
         infeasible=a.classification == :reduced_infeasible)
    end
end
tally(ev, idx) = (ok=count(e -> e.ok, ev[idx]),
                  overload=100 * maximum(e -> isfinite(e.excess) ? e.excess : Inf, ev[idx]),
                  util=maximum(e -> isnan(e.util) ? Inf : e.util, ev[idx]),
                  cost=maximum(e -> isnan(e.cost_pct) ? Inf : e.cost_pct, ev[idx]),
                  infeasible=count(e -> e.infeasible, ev[idx]))
D_idx, H_idx = 1:nd, nd+1:nd+nh

# ---- check 3: original ratings reproduce the campaign's audit ------------------
orig = Dict(run => evaluate(internal_of(run), nothing) for run in runs)
for run in runs
    camp = split(readlines(joinpath(CAMPAIGN, run, "summary.csv"))[2], ',')
    hdr = split(readlines(joinpath(CAMPAIGN, run, "summary.csv"))[1], ',')
    col(n) = parse(Int, camp[findfirst(==(n), hdr)])
    t_d, t_h = tally(orig[run], D_idx), tally(orig[run], H_idx)
    (t_d.ok == col("design_deliverable") && t_h.ok == col("heldout_deliverable")) ||
        error("check 3 failed on $run: $(t_d.ok)/$(t_h.ok) vs campaign " *
              "$(col("design_deliverable"))/$(col("heldout_deliverable"))")
end
say("check 3  original ratings reproduce the campaign's design and held-out counts on all $(length(runs)) networks")

# ---- main ----------------------------------------------------------------------
net_rows, line_rows, hour_rows = NamedTuple[], NamedTuple[], NamedTuple[]
for (variant, D) in (("hull_all168", Dall), ("hull_design10", Dd))
    set = CT.injection_set(base, D)
    for run in runs
        internal = internal_of(run)
        t0 = time()
        cert = CT.certify(base, internal, set; H)
        mind = CT.minimal_derating(base, set, cert; H)
        t_cert = time() - t0
        crit_lines = mind.lines   # box reach > 1; every other line is safe whatever the ratings
        wmax(v) = isempty(v) ? 1.0 : maximum(x -> isnan(x) ? Inf : x, v)
        load0 = CT.max_loading(base, set, cert.mp, base.frate; H, lines=crit_lines)
        load1 = CT.max_loading(base, set, cert.mp, cert.rating; H, lines=crit_lines)
        load2 = CT.max_loading(base, set, cert.mp, mind.rating; H, lines=crit_lines)
        der = evaluate(internal, cert.rating)
        dmin = evaluate(internal, mind.rating)
        o_d, o_h = tally(orig[run], D_idx), tally(orig[run], H_idx)
        d_d, d_h = tally(der, D_idx), tally(der, H_idx)
        m_d, m_h = tally(dmin, D_idx), tally(dmin, H_idx)
        crit = [r for r in cert.rows if r.kind == "external" && r.critical]
        extra(ev) = (x = [ev[s].cost_pct - orig[run][s].cost_pct for s in 1:size(Dall, 2)
                          if isfinite(ev[s].cost_pct) && isfinite(orig[run][s].cost_pct)];
                     isempty(x) ? NaN : sum(x) / length(x))
        mind_pct = isempty(crit) ? 0.0 :
            maximum(100 * (base.frate[r.line] - mind.rating[r.line]) / base.frate[r.line] for r in crit)
        push!(net_rows, (variant=variant, run=run, buses=cert.mp.K,
            merged_lines=count(internal), critical_lines=length(crit_lines),
            unsafe_merged=length(cert.unsafe_merged), critical_external=length(crit),
            worst_loading_orig=wmax(load0),
            bound_derating_pct=isempty(crit) ? 0.0 : maximum(r.derating_pct for r in crit),
            bound_worst_loading=wmax(load1),
            alpha=mind.alpha, min_derating_pct=mind_pct, min_worst_loading=wmax(load2),
            certified=isfinite(mind.alpha),
            orig_design_ok=o_d.ok, orig_heldout_ok=o_h.ok, orig_heldout_overload_pct=o_h.overload,
            orig_worst_cost_pct=max(o_d.cost, o_h.cost), orig_mean_cost_pct=
                sum(e.cost_pct for e in orig[run]) / length(orig[run]),
            bound_design_ok=d_d.ok, bound_heldout_ok=d_h.ok,
            bound_max_returned_util=max(d_d.util, d_h.util),
            bound_worst_cost_pct=max(d_d.cost, d_h.cost), bound_mean_extra_cost_pct=extra(der),
            min_design_ok=m_d.ok, min_heldout_ok=m_h.ok, min_heldout_overload_pct=m_h.overload,
            min_max_returned_util=max(m_d.util, m_h.util),
            min_infeasible_hours=m_d.infeasible + m_h.infeasible,
            min_worst_cost_pct=max(m_d.cost, m_h.cost), min_mean_extra_cost_pct=extra(dmin),
            certificate_seconds=round(t_cert; digits=2)))
        for (i, l) in enumerate(crit_lines)
            r = cert.rows[l]
            push!(line_rows, merge((variant=variant, run=run), r,
                                   (min_rating=mind.rating[l], loading_orig=load0[i],
                                    loading_bound=load1[i], loading_min=load2[i])))
        end
        for s in 1:size(Dall, 2)
            push!(hour_rows, (variant=variant, run=run, hour=hour_ids[s],
                              set=s <= nd ? "design" : "heldout",
                              orig_ok=orig[run][s].ok, orig_overload_pct=100 * orig[run][s].excess,
                              orig_cost_pct=orig[run][s].cost_pct,
                              bound_ok=der[s].ok, bound_util=der[s].util, bound_cost_pct=der[s].cost_pct,
                              min_ok=dmin[s].ok, min_overload_pct=100 * dmin[s].excess,
                              min_util=dmin[s].util, min_cost_pct=dmin[s].cost_pct))
        end
        r = last(net_rows)
        say(@sprintf("%-13s %-24s %3d buses  crit %d (merged %d)  orig load %.3f  held %3d | bound: derate %5.1f%% held %3d cost+%6.3f%% | min: alpha %.3f derate %5.1f%% load %.3f held %3d util %.4f cost+%6.3f%%  (%.1fs)",
                     variant, run, r.buses, r.critical_lines, r.unsafe_merged, r.worst_loading_orig,
                     r.orig_heldout_ok, r.bound_derating_pct, r.bound_heldout_ok, r.bound_mean_extra_cost_pct,
                     r.alpha, r.min_derating_pct, r.min_worst_loading, r.min_heldout_ok,
                     r.min_max_returned_util, r.min_mean_extra_cost_pct, t_cert))
    end
end
Main.Results.write_rows(joinpath(OUT, "networks.csv"), net_rows)
Main.Results.write_rows(joinpath(OUT, "lines.csv"), line_rows)
Main.Results.write_rows(joinpath(OUT, "hours.csv"), hour_rows)
say("wrote ", OUT)
close(log_io)

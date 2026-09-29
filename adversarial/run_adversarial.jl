# Per-contingency reduced networks for a preventive SC-DCOPF.
#
#   julia --project=. --startup-file=no adversarial/run_adversarial.jl [key=value ...]
#
# Designs one reduced network for the intact case and one per line outage, then
# solves the SC-DCOPF forms on the hull hours and on a window outside the hull:
#   full     every block a full network copy
#   lazy     PTDF rows added only when violated (constraint generation)
#   compact  the screened (block, line) pairs as PTDF rows
#   reduced  the designed networks; redundant blocks dropped
# With designs=false only the first three are compared.
# With `planning`, each of the case's candidate lines is added in turn: the
# reduced networks take the line as it is, the compact rows go stale.
#
# On a cluster the networks can be designed in parallel (hpc/submit_adversarial.sh):
#   network_index=k   design only network k of [base; outages] into output_dir, then stop
#   start_from=path   try the base network's design in that masks.csv first
#   load_designs=dir  read every part's designs from dir/*/ and run the evaluation

SETTINGS = (
    case              = :case14toy,
    month             = 11,          # the demand hull: every hour of this window
    horizon_start_day = 1,
    horizon_days      = 30,
    test_month        = 12,          # a window outside the hull
    test_start_day    = 1,
    test_days         = 31,
    designs           = true,        # false: compare full, lazy and compact only
    band              = 0.05,        # adversaries cost at most (1+band) x the loading's SC-DCOPF optimum
    tau               = 1e-6,        # overload the checks let pass, fraction of rating
    cut               = :dominance,  # rejection: :dominance | :line (the overloaded line) | :any (some kept line)
    derate            = 0.0,         # kept lines of a reduced network are rated (1 - derate) F
    margin            = nothing,     # witness: SC-DCOPF optimum at ratings (1 - margin) F; nothing = derate
    raise             = true,        # after certification, raise ratings back as far as the check allows
    candidates        = :all,        # lines that may reject an adversary under :any: :all | :screened
    adversaries_per_round = 3,       # worst adversaries added to the master per round
    master_threads    = 4,
    master_time       = 20.0,        # seconds per master solve
    master_gap        = 0.05,        # relative gap at which a master solve may stop
    fix_radial        = true,        # merge every non-critical bridge before the master
    unmerge_hops      = 2,           # outage start: unmerge merged lines this close to the outage line
    radius            = 20,          # the master may change at most this many lines of that start
    network_time      = 600.0,       # seconds per network
    max_rounds        = 200,
    networks          = nothing,     # nothing = all; e.g. [0, 3] designs only those
    network_index     = nothing,     # design only network k of [base; outages] (1-based)
    start_from        = nothing,     # masks.csv whose base design is tried first
    load_designs      = nothing,     # folder of per-network parts to evaluate instead of designing
    exact_limits      = false,       # also evaluate each design with exact limits on its critical lines only
    full_base         = false,       # evaluation: keep the intact network whole, reduce only the outages
    planning          = true,        # add each candidate line of case studies/14bus in turn
    output_dir        = nothing,     # nothing = outputs/adversarial/<case>
)

using Printf, LinearAlgebra, DelimitedFiles, JuMP
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "common", "settings.jl"))
S = Main.Settings.apply_overrides(SETTINGS)
@eval module Data
    include(joinpath(dirname(@__DIR__), "common", "preprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "postprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "multiscenario.jl"))
end
include(joinpath(ROOT, "common", "cases.jl"))
include(joinpath(@__DIR__, "adversarial.jl"))
include(joinpath(@__DIR__, "evaluate.jl"))
using .Adversarial, .Evaluate

OUT = isnothing(S.output_dir) ? joinpath(ROOT, "outputs", "adversarial", string(S.case)) : S.output_dir
mkpath(OUT)
Main.Settings.save_settings(joinpath(OUT, "settings.txt"), S)
log_io = open(joinpath(OUT, "summary.txt"), "w")
say(args...) = (println(args...); println(log_io, args...); flush(log_io))
function csv(name, header, rows)
    open(joinpath(OUT, name), "w") do io
        println(io, header)
        foreach(r -> println(io, join(r, ",")), rows)
    end
end
Main.Settings.print_settings(S)

# ---- data ---------------------------------------------------------------------------
function hours(month, day, days)
    s = (case=S.case, demands=:hourly, n_demands=1, month, horizon_days=days,
         horizon_start_day=day, opf_time_limit=60.0, linear_costs=true)
    design, held = Main.Cases.load_demands(Main.Data, s)
    isnothing(held) && return design.base, design.load, design.scenario_ids
    return design.base, hcat(design.load, held.load), vcat(design.scenario_ids, held.scenario_ids)
end
base, Dhull, ids_hull = hours(S.month, S.horizon_start_day, S.horizon_days)
_, Dtest, ids_test = hours(S.test_month, S.test_start_day, S.test_days)
net = net_from_base(base)
N, L = net.N, net.L
H = ptdf(net).H
outs = outages_of(net)
Q = lodf(net, H, outs)
say(@sprintf("%s: %d buses, %d lines, %d generators; hull %d hours, test %d hours; %d of %d outages keep it connected",
             S.case, N, L, length(net.gen_bus), size(Dhull, 2), size(Dtest, 2), length(outs), L))
all(bridges(net, 0) .== BitVector([!(l in outs) for l in 1:L])) || error("bridge search disagrees with the outage list")
# at most 20 outages, spread over the list: each check is a full PTDF
let gap = maximum(maximum(abs.(ptdf(net; out=c).H .- post_ptdf(H, Q, c)))
                  for c in outs[unique(round.(Int, range(1, length(outs); length=min(20, length(outs)))))])
    say(@sprintf("check    LODF against the PTDF with the line removed: worst gap %.1e", gap))
    gap < 1e-8 || error("LODF check failed")
end

# ---- full SC-DCOPF ------------------------------------------------------------------
t0 = time()
Zh, Gh, pool = sc_dcopf(net, H, Q, outs, Dhull)
lazy_hull = time() - t0
t0 = time()
Zt, _, pool_t = sc_dcopf(net, H, Q, outs, Dtest)
lazy_test = time() - t0
hull, test = findall(isfinite, Zh), findall(isfinite, Zt)
say(@sprintf("SC-DCOPF: hull %d of %d hours feasible, test %d of %d; %d post-outage rows (%.1fs + %.1fs)",
             length(hull), length(Zh), length(test), length(Zt), pool, lazy_hull, lazy_test))
isempty(hull) && error("no feasible hull hour")
D, Z = Dhull[:, hull], Zh[hull]
# The witness every reduced network must accept, one per hour. With a margin it is
# the SC-DCOPF optimum under ratings lowered by the margin: it costs a little more,
# and its slack lets a reduced network lower ratings the optimum would pin.
Gw = Gh[:, hull]
margin = isnothing(S.margin) ? S.derate : S.margin
margin >= S.derate || error("margin must be at least derate, or the unreduced network rejects the witness")
if margin > 0
    Zm, Gm, _ = sc_dcopf(with_ratings(net, net.F .* (1 - margin)), H, Q, outs, D)
    okm = findall(isfinite, Zm)
    wgap = 100 .* (Zm[okm] ./ Z[okm] .- 1)
    say(@sprintf("witness (margin %g%%): %d of %d hours, cost over the SC-DCOPF optimum mean %.3f%%, max %.3f%%",
                 100 * margin, length(okm), length(hull), sum(wgap) / length(wgap), maximum(wgap)))
    maximum(wgap) <= 100 * S.band || say("warning: the witness costs more than the band covers")
    # an hour with no secure dispatch at the margin has no witness a derated
    # network can accept, so it leaves the design hull
    if length(okm) < length(hull)
        say("         ", length(hull) - length(okm), " hours have none and leave the design hull")
        hull, D, Z = hull[okm], D[:, okm], Z[okm]
    end
    Gw = Gm[:, okm]
end
Pw = reduce(hcat, [injection(net, Gw[:, s], D[:, s]) for s in axes(D, 2)])
csv("hours.csv", "set,hour,hour_id,sc_dcopf_cost",
    vcat([("hull", s, ids_hull[s], Zh[s]) for s in eachindex(Zh)],
         [("test", s, ids_test[s], Zt[s]) for s in eachindex(Zt)]))

# ---- screening ----------------------------------------------------------------------
# Outages share one adversary set (the intact network's), so the extremes of the
# intact flows bound every post-outage flow and clear most pairs without an LP.
function screen_all(n, H, Q, outs, D, Z; only=[0; outs])
    A0 = adversary_set(n, D, Z, S.band)
    Ac = adversary_set(n, D, Z, S.band; intact=H)
    # the bounds need every intact line's extremes: skipped when only the base is screened
    some = any(!=(0), only)
    hi = some ? [reach(Ac, H[l, :], 1)[1] for l in 1:n.L] : Float64[]
    lo = some ? [-reach(Ac, H[l, :], -1)[1] for l in 1:n.L] : Float64[]
    nets = Dict{Int,Any}()
    for c in only
        Hc = post_ptdf(H, Q, c)
        bound = c == 0 ? nothing : (l, s) -> (q = Q[l, c];
            s == 1 ? hi[l] + (q >= 0 ? q * hi[c] : q * lo[c]) :
                     -(lo[l] + (q >= 0 ? q * lo[c] : q * hi[c])))
        A = c == 0 ? A0 : Ac
        nets[c] = (; Hc, A, crit=screen(A, Hc, [l for l in 1:n.L if l != c]; tau=S.tau, bound))
    end
    return nets
end
chosen = S.networks
if !isnothing(S.network_index)
    all_nets = [0; outs]
    if S.network_index > length(all_nets)
        say("network_index ", S.network_index, " is past the last of ", length(all_nets), " networks")
        close(log_io)
        exit()
    end
    chosen = [all_nets[S.network_index]]
end
# a run that designs only some networks screens only those (the evaluation needs all)
only = S.designs && isnothing(S.load_designs) && !isnothing(chosen) ? chosen : [0; outs]
t0 = time()
scr = screen_all(net, H, Q, outs, D, Z; only)
redundant = [c for c in only if isempty(scr[c].crit)]
say(@sprintf("screening (band %g%%): %d of %d networks redundant (%.1fs)",
             100 * S.band, length(redundant), length(only), time() - t0))

# ---- designs ------------------------------------------------------------------------
name(c) = c == 0 ? "base" : string(c)
designs = Dict{Int,Any}()
net_rows, mask_rows = [], []
start = nothing
if !isnothing(S.start_from)
    rows, _ = readdlm(S.start_from, ',', header=true)
    start = falses(L)
    for r in eachrow(rows)
        string(r[1]) == "base" && (start[Int(r[2])] = Int(r[3]) == 1)
    end
    say("start design from ", S.start_from, ": ", maximum(clusters(net, start)), " buses")
end
t0 = time()
for c in (S.designs && isnothing(S.load_designs) ? [0; outs] : Int[])
    isnothing(chosen) || c in chosen || continue
    s = scr[c]
    if isempty(s.crit)
        push!(net_rows, (name(c), 0, true, 1, 0, 0.0, "redundant", 0, 0, 0, 0.0, 0.0, 0.0))
        continue
    end
    say("network ", name(c), ": ", length(s.crit), " critical pairs")
    r = design(net, c, s.Hc, s.A, s.crit, Pw; tau=S.tau, candidates=S.candidates, cut=S.cut,
               derate=S.derate, raise=S.raise,
               adversaries_per_round=S.adversaries_per_round, threads=S.master_threads,
               master_time=S.master_time, time_limit=S.network_time, max_rounds=S.max_rounds,
               gap=S.master_gap, fix_radial=S.fix_radial, unmerge_hops=S.unmerge_hops,
               radius=S.radius, log=say,
               # outages start from the base design: a given one, else this run's
               start=(c == 0 ? nothing : !isnothing(start) ? start :
                      haskey(designs, 0) ? designs[0].mask : nothing))
    designs[c] = r
    nb = maximum(clusters(net, r.mask))
    kept = count(l -> l != c && !r.mask[l], 1:L)
    derate = 100 * (1 - minimum(r.rating ./ net.F))
    push!(net_rows, (name(c), length(s.crit), false, nb, kept, round(derate; digits=3), r.status, r.rounds, r.adversaries,
                     r.witnesses, round(r.master_s; digits=3), round(r.check_s; digits=3),
                     round(r.seconds; digits=3)))
    append!(mask_rows, [(name(c), l, Int(r.mask[l]), r.rating[l] / net.F[l]) for l in 1:L if l != c])
    say(@sprintf("network %4s: %2d critical pairs -> %2d buses, %2d lines, max derating %.2f%%, %s after %d rounds (%d adversaries, %d witness hours), %.1fs",
                 name(c), length(s.crit), nb, kept, derate, r.status, r.rounds, r.adversaries, r.witnesses, r.seconds))
end
S.designs && say(@sprintf("designs: %.1fs", time() - t0))
csv("networks.csv", "network,critical_pairs,redundant,buses,kept_lines,max_derating_pct,status,rounds,adversaries,witness_hours,master_s,check_s,seconds", net_rows)
csv("masks.csv", "network,line,internal,rating_over_F", mask_rows)
if !isnothing(chosen)
    say("only some networks designed: evaluation skipped")
    close(log_io)
    exit()
end

# designs written by separate runs, one folder each (see network_index)
if !isnothing(S.load_designs)
    net_rows = []
    for part in sort(readdir(S.load_designs; join=true))
        isfile(joinpath(part, "networks.csv")) && isfile(joinpath(part, "masks.csv")) || continue
        countlines(joinpath(part, "networks.csv")) > 1 || continue
        nrows, _ = readdlm(joinpath(part, "networks.csv"), ',', Any, header=true)
        mrows = countlines(joinpath(part, "masks.csv")) > 1 ?
                readdlm(joinpath(part, "masks.csv"), ',', Any, header=true)[1] : Matrix{Any}(undef, 0, 4)
        for r in eachrow(nrows)
            push!(net_rows, Tuple(r))
            (r[3] == true || string(r[3]) == "true") && continue
            c = string(r[1]) == "base" ? 0 : Int(r[1])
            mask, rating = falses(L), copy(net.F)
            for m in eachrow(mrows)
                string(m[1]) == string(r[1]) || continue
                mask[Int(m[2])] = Int(m[3]) == 1
                rating[Int(m[2])] = Float64(m[4]) * net.F[Int(m[2])]
            end
            designs[c] = (; mask, rating, status=Symbol(r[7]))
        end
    end
    missing_nets = [c for c in [0; outs] if !isempty(scr[c].crit) && !haskey(designs, c)]
    for c in missing_nets
        designs[c] = (; mask=falses(L), rating=copy(net.F), status=:missing)
    end
    say(@sprintf("loaded %d designs from %s; %d non-redundant networks missing (kept unreduced)",
                 length(designs) - length(missing_nets), S.load_designs, length(missing_nets)))
    csv("networks.csv", "network,critical_pairs,redundant,buses,kept_lines,max_derating_pct,status,rounds,adversaries,witness_hours,master_s,check_s,seconds", net_rows)
end

# ---- check: the adversary again, with the reduced network in angle form ---------------
function angle_worst(A, n, out, mask, rating, Hc, crit)
    m, lines, K = A.m, [l for l in 1:n.L if l != out], length(n.gen_bus)
    th = @variable(m, [1:n.N])
    f = @variable(m, [lines])
    t = @variable(m, [lines])
    cons = ConstraintRef[@constraint(m, th[n.ref] == 0)]
    for l in lines
        push!(cons, @constraint(m, f[l] == n.b[l] * (th[n.from[l]] - th[n.to[l]])))
        if mask[l]
            push!(cons, @constraint(m, f[l] == 0))
        else
            push!(cons, @constraint(m, t[l] == 0), @constraint(m, f[l] <= rating[l]),
                  @constraint(m, -f[l] <= rating[l]))
        end
    end
    for i in 1:n.N
        push!(cons, @constraint(m,
            sum(f[l] + t[l] for l in lines if n.from[l] == i; init=0.0) -
            sum(f[l] + t[l] for l in lines if n.to[l] == i; init=0.0) ==
            sum(A.g[k] for k in 1:K if n.gen_bus[k] == i; init=0.0) -
            sum(A.D[i, s] * A.lam[s] for s in axes(A.D, 2))))
    end
    w = maximum(reach(A, view(Hc, l, :), s)[1] / n.F[l] for (l, s, _) in crit)
    foreach(c -> delete(m, c), cons)
    delete(m, th); delete(m, f.data); delete(m, t.data)
    return w
end
S.designs && let gap = 0.0, worst = 0.0
    for (c, r) in designs
        r.status === :certified || continue
        s = scr[c]
        w = check(net, c, s.Hc, s.A, s.crit, Pw, r.mask; tau=S.tau, rating=r.rating).worst
        gap = max(gap, abs(w - angle_worst(s.A, net, c, r.mask, r.rating, s.Hc, s.crit)))
        worst = max(worst, w)
    end
    say(@sprintf("check    certified designs: worst adversary loading %.6f; angle form vs PTDF form gap %.1e", worst, gap))
end

# ---- exact limits: critical lines only -------------------------------------------------
# Each certified design gets the limits its critical lines need (exact_limits) and
# none elsewhere, and keeps them only if the result still passes the check;
# otherwise it keeps its own ratings. Also counts the designs that accept the exact
# SC-DCOPF optimum at every hour, which makes the reduced cost exact.
Pexact = reduce(hcat, [injection(net, Gh[:, s], Dhull[:, s]) for s in hull])
designs_x = Dict{Int,Any}()
S.designs && S.exact_limits && let changed = 0, merged_crit = 0, failed = 0, exact_ok = 0, lim_rows = []
    for c in sort(collect(keys(designs)))
        r, s = designs[c], scr[c]
        rt = r.status === :certified ? exact_limits(net, c, s.Hc, s.A, s.crit, r.mask; tau=S.tau) : nothing
        designs_x[c] = r
        if isnothing(rt)
            merged_crit += r.status === :certified
        else
            ck = check(net, c, s.Hc, s.A, s.crit, Pw, r.mask; tau=S.tau, rating=rt)
            if isempty(ck.bad) && isempty(ck.advs)
                changed += 1
                designs_x[c] = (; mask=r.mask, rating=rt, status=r.status)
            else
                failed += 1
            end
        end
        x = designs_x[c]
        mp = ptdf(net; internal=x.mask, out=c)
        ext = findall(mp.external)
        Fx = mp.H[ext, :] * Pexact
        acc = all(abs(Fx[i, t]) <= x.rating[ext[i]] + S.tau * net.F[ext[i]] for i in eachindex(ext), t in axes(Fx, 2))
        exact_ok += acc
        push!(lim_rows, (name(c), x === r ? "design" : "exact", count(l -> isfinite(x.rating[l]), ext), length(ext),
                         join([@sprintf("%d:%.4f", l, x.rating[l] / net.F[l]) for l in unique(first.(s.crit))], " "), acc))
    end
    say(@sprintf("exact limits: %d designs changed, %d kept their own (%d with a merged critical line, %d failing the check); %d of %d accept the exact SC-DCOPF optimum at every hour",
                 changed, merged_crit + failed, merged_crit, failed, exact_ok, length(designs)))
    csv("exact_limits.csv", "network,limits,limited_lines,kept_lines,critical_limit_over_F,accepts_exact_optimum", lim_rows)
end

# ---- the SC-DCOPF forms ---------------------------------------------------------------
blocks_full = [(c, falses(L)) for c in [0; outs]]
# with full_base the reduced forms keep the intact network whole and reduce only the outages
with_base(blocks) = S.full_base ? [(0, falses(L), copy(net.F)); [b for b in blocks if b[1] != 0]] : blocks
blocks_red = with_base([(c, designs[c].mask, designs[c].rating) for c in sort(collect(keys(designs)))])
rows_compact = [(scr[c].Hc[l, :], l) for c in [0; outs] for l in unique(first.(scr[c].crit))]
forms = [("full", () -> network_scdcopf(net, blocks_full)),
         ("compact", () -> compact_scdcopf(net, rows_compact))]
S.designs && push!(forms, ("reduced", () -> network_scdcopf(net, blocks_red)))
isempty(designs_x) || push!(forms, ("reduced_exact", () -> network_scdcopf(net,
    with_base([(c, designs_x[c].mask, designs_x[c].rating) for c in sort(collect(keys(designs_x)))]))))
"Worst loading of a dispatch and where it happens: (loading, network (0 = intact), line)."
function worst_where(g, d)
    f = H * injection(net, g, d)
    v, l = findmax(abs.(f) ./ net.F)
    best = (v, 0, l)
    for c in outs, k in 1:L
        k == c && continue
        x = abs(f[k] + Q[k, c] * f[c]) / net.F[k]
        x > best[1] && (best = (x, c, k))
    end
    return best
end
# lazy: one generator vector, the base rows, the balance and the pooled rows it needed
eval_rows = Any[("lazy", "hull", length(hull), length(net.gen_bus), 2L + 2pool + 1, "", "", round(lazy_hull; digits=4), length(Zh) - length(hull), 0.0, 0.0),
                ("lazy", "test", length(test), length(net.gen_bus), 2L + 2pool_t + 1, "", "", round(lazy_test; digits=4), length(Zt) - length(test), 0.0, 0.0)]
say(@sprintf("compact rows: %d (block, line) pairs from screening, of %d possible", length(rows_compact), (length(outs) + 1) * L))
# how far each test hour's demand is from the design hull: L1 distance over total
# demand, 0 = inside (where the reduced networks are guaranteed)
hull_dist = let m = Adversarial.lp(), T = size(D, 2)
    @variable(m, lam[1:T] >= 0)
    @variable(m, up[1:N] >= 0)
    @variable(m, dn[1:N] >= 0)
    @constraint(m, sum(lam) == 1)
    fit = @constraint(m, D * lam .+ up .- dn .== 0.0)
    @objective(m, Min, sum(up) + sum(dn))
    map(test) do s
        set_normalized_rhs.(fit, Dtest[:, s])
        optimize!(m)
        objective_value(m) / sum(Dtest[:, s])
    end
end
inside = findall(<=(1e-6), hull_dist)
say(@sprintf("test hours inside the design hull: %d of %d (the rest up to %.1f%% of demand away)",
             length(inside), length(test), 100 * maximum(hull_dist; init=0.0)))
hour_rows = []
for (form, build) in forms, (set, Ds, Zs, ids, dist) in (("hull", D, Z, ids_hull[hull], zeros(length(hull))),
                                                          ("test", Dtest[:, test], Zt[test], ids_test[test], hull_dist))
    t1 = time()
    model = build()
    tb = time() - t1
    cost, disp, secs = solve_hours(model, Ds)
    sz = model_size(model.m)
    ok = findall(isfinite, cost)
    gap = isempty(ok) ? NaN : maximum(abs(100 * (cost[s] - Zs[s]) / Zs[s]) for s in ok)
    over = isempty(ok) ? NaN : maximum(worst_loading(net, H, Q, outs, disp[:, s], Ds[:, s]) for s in ok)
    push!(eval_rows, (form, set, size(Ds, 2), sz.vars, sz.cons, sz.nnz, round(tb; digits=4),
                      round(secs; digits=4), size(Ds, 2) - length(ok), gap, 100 * max(0.0, over - 1)))
    say(@sprintf("%-8s %-4s: %5d vars %5d rows %6d nnz | build %.3fs, solve %.3fs over %d hours | %d infeasible, worst cost gap %.2e%%, worst overload %.2e%%",
                 form, set, sz.vars, sz.cons, sz.nnz, tb, secs, size(Ds, 2), size(Ds, 2) - length(ok),
                 gap, 100 * max(0.0, over - 1)))
    # per hour: signed cost gap (below 0 = cheaper than the secure optimum) and worst loading
    places = [isfinite(cost[s]) ? worst_where(disp[:, s], Ds[:, s]) : (NaN, "", "") for s in axes(Ds, 2)]
    append!(hour_rows, [(form, set, ids[s], 100 * dist[s], Zs[s], cost[s], 100 * (cost[s] - Zs[s]) / Zs[s], places[s]...)
                        for s in axes(Ds, 2)])
    if startswith(form, "reduced") && !isempty(ok)
        bad = [s for s in ok if places[s][1] > 1 + 1e-6]
        w = isempty(bad) ? nothing : places[argmax(s -> places[s][1], bad)]
        isempty(bad) || say(@sprintf("         %d hours overload (%d inside the hull); worst %.3f%% on line %d after %s",
            length(bad), count(s -> dist[s] <= 1e-6, bad), 100 * (w[1] - 1), w[3],
            w[2] == 0 ? "no outage (intact network)" : "outage $(w[2])"))
    end
end
csv("forms.csv", "form,hours,n_hours,vars,rows,nnz,build_s,solve_s,infeasible_hours,worst_cost_gap_pct,worst_overload_pct", eval_rows)
csv("hours_eval.csv", "form,hours,hour_id,hull_distance_pct,sc_dcopf_cost,cost,cost_gap_pct,worst_loading,worst_outage,worst_line", hour_rows)

# ---- network form: add each candidate line ---------------------------------------------
if S.planning && S.designs
    file = joinpath(ROOT, "case studies", "14bus", "network_planning.csv")
    cand, _ = readdlm(file, ',', header=true)
    # the compact rows were computed without the new line, so its dispatch does not change
    cost_c, disp_c, _ = solve_hours(compact_scdcopf(net, rows_compact), D)
    plan_rows = []
    for r in eachrow(cand)
        id, u, v, x, fmax = Int(r[1]), Int(r[2]), Int(r[3]), Float64(r[4]), Float64(r[5])
        n2 = add_line(net, u, v, 1 / x, fmax / base.baseMVA)
        H2 = ptdf(n2).H
        Q2 = lodf(n2, H2, outs)
        Z2, G2, _ = sc_dcopf(n2, H2, Q2, outs, D)
        ok2 = findall(isfinite, Z2)
        cr, dr, _ = solve_hours(network_scdcopf(n2, [(c, [m; false], [rt; n2.F[end]]) for (c, m, rt) in blocks_red]), D)
        Z2m, G2m, _ = margin > 0 ? sc_dcopf(with_ratings(n2, n2.F .* (1 - margin)), H2, Q2, outs, D) : (Z2, G2, 0)
        Gw2 = [isfinite(Z2m[s]) ? G2m[:, s] : G2[:, s] for s in ok2]
        Pw2 = reduce(hcat, [injection(n2, Gw2[i], D[:, s]) for (i, s) in enumerate(ok2)])
        scr2 = screen_all(n2, H2, Q2, outs, D[:, ok2], Z2[ok2])
        held = count(keys(designs)) do c
            s2 = scr2[c]
            ck = check(n2, c, s2.Hc, s2.A, s2.crit, Pw2, [designs[c].mask; false]; tau=S.tau,
                       rating=[designs[c].rating; n2.F[end]])
            isempty(ck.bad) && isempty(ck.advs)
        end
        woke = count(c -> !haskey(designs, c) && !isempty(scr2[c].crit), [0; outs])
        need = sum(length(unique(first.(scr2[c].crit))) for c in [0; outs])
        for (form, cost, disp) in (("reduced", cr, dr), ("compact_stale", cost_c, disp_c))
            ok = [s for s in ok2 if isfinite(cost[s])]
            gap = isempty(ok) ? NaN : maximum(abs(100 * (cost[s] - Z2[s]) / Z2[s]) for s in ok)
            over = isempty(ok) ? NaN : maximum(worst_loading(n2, H2, Q2, outs, disp[:, s], D[:, s]) for s in ok)
            push!(plan_rows, (id, u, v, form, length(ok2) - length(ok), gap, 100 * max(0.0, over - 1),
                              held, length(designs), woke, need, length(rows_compact)))
            say(@sprintf("line %d (%d-%d) %-13s: %d infeasible, worst cost gap %.3f%%, worst overload %.3f%% | designs still certified %d of %d, dropped blocks now critical %d, compact rows %d needed vs %d stale",
                         id, u, v, form, length(ok2) - length(ok), gap, 100 * max(0.0, over - 1),
                         held, length(designs), woke, need, length(rows_compact)))
        end
    end
    csv("planning.csv", "candidate,from,to,form,infeasible_hours,worst_cost_gap_pct,worst_overload_pct,designs_still_certified,designs,dropped_now_critical,compact_rows_needed,compact_rows_stale", plan_rows)
end
say("wrote ", OUT)
close(log_io)

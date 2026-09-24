# Comprehensive review of network reduction on the 6-bus toy case.
#
#   julia --project=. --startup-file=no analysis/review_6bus.jl [scale=0.018]
#
# For each candidate clustering it solves the REDUCED DC-OPF at every one of the
# 8760 hours, puts the resulting dispatch back on the FULL network, and measures
# how far the true flows exceed their ratings. Design hours (the window the
# reduction was built from) and the remaining out-of-sample hours are scored
# separately, which is the only way to see whether a reduction generalises.
#
# Violations are reported three ways because they answer different questions:
#   MW   how much copper is missing
#   %    overload relative to that line's own rating  (severity)
#   p.u. overload on the 100 MVA system base          (system-wide scale)
#
# Writes outputs/6busww/review/*.csv; analysis/xlsx_write.py turns them into a
# workbook and analysis/plot_6bus_review.jl draws the figures.

using JuMP, Gurobi, Printf, DelimitedFiles, Statistics

const ROOT = dirname(@__DIR__)
const DATA = joinpath(ROOT, "case studies", "6busww")
const OUT  = joinpath(ROOT, "outputs", "6busww", "review")
const BASEMVA = 100.0
const TOL = 1e-6

scale = 0.018
for a in ARGS
    k, v = split(a, Char(61); limit=2)
    k == "scale" || error("only scale=<float> is accepted")
    global scale = parse(Float64, v)
end
const DATA_SF = 0.08

rd(f) = readdlm(joinpath(DATA, f), ',', Any, '\n'; header=true)[1]
G, Lr = rd("network_gens.csv"), rd("network_lines.csv")
gbus = Int.(G[:, 2]); gpmax = Float64.(G[:, 3]); gcost = Float64.(G[:, 5])
lfrom = Int.(Lr[:, 2]); lto = Int.(Lr[:, 3])
lfmax = Float64.(Lr[:, 5]); sus = BASEMVA ./ Float64.(Lr[:, 4])

raw = rd("hourly_demand.csv")
stamp = Int.(raw[:, 2:5])
dem = Float64.(raw[:, 6:end]) .* (scale / DATA_SF)
ren = Float64.(rd("hourly_wind.csv")[:, 6:end]) .+ Float64.(rd("hourly_solar.csv")[:, 6:end])
D = dem .- ren
H, N = size(D); L, K = length(lfrom), length(gbus)
lbl(l) = "$(lfrom[l])-$(lto[l])"

# The design window: month 11, days 8..14, matching the runners.
inwin = [stamp[h, 2] == 11 && 8 <= stamp[h, 3] <= 14 for h in 1:H]
window = findall(inwin)
@printf("scale %.3f   %d hours   design window = Nov 8-14 (%d h)\n", scale, H, length(window))

# --------------------------------------------------------------------------- #
# Two models, built once. `red` is the reduced OPF for a given clustering;
# `full` puts a fixed dispatch back on the original network.
# --------------------------------------------------------------------------- #
env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))

function reduced_model(internal::AbstractVector{Bool})
    m = Model(() -> Gurobi.Optimizer(env)); set_silent(m)
    @variable(m, 0 <= g[k=1:K] <= gpmax[k])
    @variable(m, th[1:N]); @variable(m, f[1:L]); @variable(m, t[1:L])
    @variable(m, shed[1:N] >= 0)
    # Net load goes negative at the renewable buses, so the reduced OPF needs a
    # spill variable or it is simply infeasible in every high-renewable hour --
    # which would show up as a "violation" that has nothing to do with the
    # reduction. Bounded per hour by the local renewable output.
    @variable(m, curt[1:N] >= 0)
    @constraint(m, th[1] == 0)
    @constraint(m, [l=1:L], f[l] == sus[l] * (th[lfrom[l]] - th[lto[l]]))
    for l in 1:L
        if internal[l]
            @constraint(m, f[l] == 0)          # merged: angles tied, transfer free
        else
            @constraint(m, -lfmax[l] <= f[l] <= lfmax[l])
            @constraint(m, t[l] == 0)
        end
    end
    @constraint(m, bal[i=1:N],
        sum(g[k] for k in 1:K if gbus[k] == i; init=0.0) + shed[i] - curt[i]
        - sum(f[l] + t[l] for l in 1:L if lfrom[l] == i; init=0.0)
        + sum(f[l] + t[l] for l in 1:L if lto[l]   == i; init=0.0) == 0.0)
    @objective(m, Min, sum(gcost[k]*g[k] for k in 1:K) + 10_000*sum(shed) + 1_000*sum(curt))
    return (m=m, g=g, bal=bal, shed=shed, curt=curt)
end

# true network under a FIXED dispatch: is it deliverable?
fullm = Model(() -> Gurobi.Optimizer(env)); set_silent(fullm)
@variable(fullm, tth[1:N]); @variable(fullm, tf[1:L])
@constraint(fullm, tth[1] == 0)
@constraint(fullm, [l=1:L], tf[l] == sus[l] * (tth[lfrom[l]] - tth[lto[l]]))
@constraint(fullm, tbal[i=1:N],
    - sum(tf[l] for l in 1:L if lfrom[l] == i; init=0.0)
    + sum(tf[l] for l in 1:L if lto[l]   == i; init=0.0) == 0.0)
@objective(fullm, Min, 0.0)

"Solve the reduced OPF at every hour, then measure the true-network overload."
function evaluate(internal::AbstractVector{Bool})
    R = reduced_model(internal)
    cost = zeros(H); shedmw = zeros(H)
    ovl_mw = zeros(H); ovl_pct = zeros(H); worst_line = zeros(Int, H)
    infeas = falses(H)
    for h in 1:H
        for i in 1:N
            set_normalized_rhs(R.bal[i], D[h, i])
            set_upper_bound(R.shed[i], max(D[h, i], 0.0))
            set_upper_bound(R.curt[i], ren[h, i])
        end
        optimize!(R.m)
        if termination_status(R.m) != MOI.OPTIMAL
            infeas[h] = true; ovl_mw[h] = Inf; continue
        end
        gv = value.(R.g); sv = value.(R.shed); cv = value.(R.curt)
        cost[h] = sum(gcost[k]*gv[k] for k in 1:K); shedmw[h] = sum(sv)
        # put that dispatch on the full network
        for i in 1:N
            inj = sum(gv[k] for k in 1:K if gbus[k] == i; init=0.0) + sv[i] - cv[i] - D[h, i]
            set_normalized_rhs(tbal[i], -inj)
        end
        optimize!(fullm)
        if termination_status(fullm) != MOI.OPTIMAL
            infeas[h] = true; ovl_mw[h] = Inf; continue
        end
        tfv = value.(tf)
        for l in 1:L
            over = abs(tfv[l]) - lfmax[l]
            if over > ovl_mw[h]
                ovl_mw[h] = over; worst_line[h] = l; ovl_pct[h] = 100 * over / lfmax[l]
            end
        end
        ovl_mw[h] = max(ovl_mw[h], 0.0); ovl_pct[h] = max(ovl_pct[h], 0.0)
    end
    return (cost=cost, shed=shedmw, ovl_mw=ovl_mw, ovl_pct=ovl_pct,
            worst_line=worst_line, infeas=infeas)
end

clusters(internal) = begin
    p = collect(1:N)
    rt(x) = p[x] == x ? x : (p[x] = rt(p[x]))
    for l in findall(internal); p[rt(lfrom[l])] = rt(lto[l]); end
    length(unique(rt(b) for b in 1:N))
end

# --------------------------------------------------------------------------- #
mask(ls) = (v = falses(L); for l in ls; v[l] = true; end; v)
cases = [
    ("base (no reduction)",        Int[]),
    ("KKT optimum {2-4, 3-6}",     [5, 9]),
    ("single {2-4}",               [5]),
    ("single {3-6}",               [9]),
    ("single {2-3}",               [4]),
    ("single {2-5}",               [6]),
    ("single {2-6}",               [7]),
    ("single {1-5}",               [3]),
    # Each of 2-3, 2-5, 2-6 is feasible on its own. Their unions test the other
    # direction of monotonicity: does a set of individually safe merges stay safe?
    ("union {2-3, 2-5}",           [4, 6]),
    ("union {2-3, 2-5, 2-6}",      [4, 6, 7]),
    ("copper plate (all 11)",      collect(1:L)),
]

base_eval = evaluate(falses(L))       # the full network's own optimum, per hour

rows = NamedTuple[]
hourly = Dict{String,Any}()
println()
@printf("%-24s %6s %8s  | %-28s | %-28s\n", "", "buses", "merged",
        "DESIGN WINDOW (168 h)", "OUT OF SAMPLE (8592 h)")
@printf("%-24s %6s %8s  | %6s %8s %8s | %6s %8s %8s %8s\n", "case", "", "",
        "viol h", "max MW", "max %", "viol h", "max MW", "max %", "max p.u.")
println("-"^122)
for (name, ls) in cases
    im = mask(ls)
    e = evaluate(im)
    hourly[name] = e
    out = setdiff(1:H, window)
    vw = [h for h in window if e.ovl_mw[h] > TOL || e.infeas[h]]
    vo = [h for h in out    if e.ovl_mw[h] > TOL || e.infeas[h]]
    fin(x) = isempty(x) ? 0.0 : maximum(y -> isfinite(y) ? y : 0.0, x)
    mw_w  = fin([e.ovl_mw[h]  for h in vw]); pc_w = fin([e.ovl_pct[h] for h in vw])
    mw_o  = fin([e.ovl_mw[h]  for h in vo]); pc_o = fin([e.ovl_pct[h] for h in vo])
    dcost = sum(e.cost) - sum(base_eval.cost)
    # Per-HOUR cost gap against the full network's own optimum. The KKT cap is a
    # per-scenario condition, so an annual average can sit inside 0.1% while a
    # single hour blows straight through it.
    gap(h) = base_eval.cost[h] > 1e-9 ?
             100 * (e.cost[h] - base_eval.cost[h]) / base_eval.cost[h] : 0.0
    gw = isempty(window) ? 0.0 : maximum(gap, window)
    go = isempty(out) ? 0.0 : maximum(gap, out)
    push!(rows, (case=name, buses=clusters(im), merged=length(ls),
        design_violating_hours=length(vw), design_max_ovl_mw=round(mw_w; digits=4),
        design_max_ovl_pct=round(pc_w; digits=3),
        design_max_ovl_pu=round(mw_w/BASEMVA; digits=6),
        oos_violating_hours=length(vo), oos_max_ovl_mw=round(mw_o; digits=4),
        oos_max_ovl_pct=round(pc_o; digits=3),
        oos_max_ovl_pu=round(mw_o/BASEMVA; digits=6),
        oos_violating_pct_of_hours=round(100*length(vo)/length(out); digits=3),
        annual_cost=round(sum(e.cost); digits=2),
        annual_cost_change_pct=round(100*dcost/sum(base_eval.cost); digits=4),
        design_max_cost_gap_pct=round(gw; digits=4),
        oos_max_cost_gap_pct=round(go; digits=4),
        unserved_mwh=round(sum(e.shed); digits=2)))
    @printf("%-24s %6d %8d  | %6d %8.3f %8.2f | %6d %8.3f %8.2f %8.4f\n",
            name, clusters(im), length(ls), length(vw), mw_w, pc_w,
            length(vo), mw_o, pc_o, mw_o/BASEMVA)
end

println("\nWORST PER-HOUR COST GAP vs the full network's own optimum")
@printf("%-24s %14s %14s   %s\n", "case", "design max %", "out-of-sample %",
        "(KKT ran with a 0.1% cap)")
println("-"^78)
for r in rows
    flag = r.design_max_cost_gap_pct > 0.1 ? "  <- over the 0.1% cap" : ""
    @printf("%-24s %14.4f %14.4f%s\n", r.case, r.design_max_cost_gap_pct,
            r.oos_max_cost_gap_pct, flag)
end

mkpath(OUT)
hdr = collect(string.(keys(rows[1])))
writedlm(joinpath(OUT, "summary.csv"),
         vcat(permutedims(hdr), [getfield(r, Symbol(c)) for r in rows, c in hdr]), ',')

# per-line overload frequency, out of sample, for the infeasible cases
lrows = NamedTuple[]
for (name, ls) in cases
    e = hourly[name]
    out = setdiff(1:H, window)
    for l in 1:L
        hits = count(h -> e.worst_line[h] == l && e.ovl_mw[h] > TOL, out)
        hits == 0 && continue
        mws = [e.ovl_mw[h] for h in out if e.worst_line[h] == l && e.ovl_mw[h] > TOL]
        push!(lrows, (case=name, line=l, path=lbl(l), rating_mw=lfmax[l],
            hours_worst=hits, max_ovl_mw=round(maximum(mws); digits=4),
            max_ovl_pct=round(100*maximum(mws)/lfmax[l]; digits=3),
            max_ovl_pu=round(maximum(mws)/BASEMVA; digits=6),
            mean_ovl_mw=round(mean(mws); digits=4)))
    end
end
if !isempty(lrows)
    h2 = collect(string.(keys(lrows[1])))
    writedlm(joinpath(OUT, "per_line_violations.csv"),
             vcat(permutedims(h2), [getfield(r, Symbol(c)) for r in lrows, c in h2]), ',')
end

# hourly trace for the two headline cases
for (name, tag) in (("KKT optimum {2-4, 3-6}", "kkt_opt"), ("single {3-6}", "single_3_6"))
    e = hourly[name]
    writedlm(joinpath(OUT, "hourly_$(tag).csv"),
        vcat(permutedims(["hour","month","day","hour_of_day","in_window","cost",
                          "unserved_mw","ovl_mw","ovl_pct","ovl_pu","worst_line"]),
             hcat(1:H, stamp[:,2], stamp[:,3], stamp[:,4], Int.(inwin), e.cost,
                  e.shed, e.ovl_mw, e.ovl_pct, e.ovl_mw ./ BASEMVA, e.worst_line)), ',')
end

@printf("\nwrote summary.csv, per_line_violations.csv, 2 hourly traces -> %s\n", OUT)

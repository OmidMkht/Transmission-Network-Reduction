# Profile for the per-contingency reduction (Adversarial_Injection_Generation_for_Safe_Network_Reduction.pdf).
#
#   julia --project=. --startup-file=no analysis/sc_profile.jl
#
# ACTIVSg200, Mar 1-7 (168 hours), normal ratings, preventive N-1 over every line
# whose outage keeps the network connected.
#   1. full SC-DCOPF per hour (constraint generation): Z_SC and the base-OPF cost Z0
#   2. for each band x: which lines some adversary can overload
#        base network:     box, balance, demand hull, cost <= (1+x) * interpolated Z_SC
#        outage c:         the same, plus feasible on the full intact network
#      A network with no such line is redundant (its reduced network is one bus).
# Writes outputs/sc_profile/ACTIVSg200/{hours.csv, networks.csv, pairs.csv, summary.txt}.

using Printf, LinearAlgebra, Statistics, JuMP, Gurobi
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
include(joinpath(ROOT, "certify", "certificate.jl"))
const CT = Main.Certificate

const BANDS = [0.0, 0.01, 0.05, Inf]
const TAU = 1e-6

OUT = joinpath(ROOT, "outputs", "sc_profile", "ACTIVSg200")
mkpath(OUT)
log_io = open(joinpath(OUT, "summary.txt"), "w")
say(args...) = (println(args...); println(log_io, args...); flush(log_io))

design, heldout = Main.Cases.load_demands(Main.Data, S)
base = design.base
D = hcat(design.load, heldout.load)
hour_ids = vcat(design.scenario_ids, heldout.scenario_ids)
N, L, K, T = base.N, base.Ln, length(base.gen_bus), size(D, 2)
u, v, F, gb = base.Efrom, base.Eto, base.frate, base.gen_bus
pmin, pmax = min.(0.0, base.pmin), base.pmax
c1, c0 = base.c1, sum(base.c0)
H = CT.ptdf(base)
say(@sprintf("ACTIVSg200: %d buses, %d lines, %d generators, %d hours", N, L, K, T))

# ---- outages that keep the network connected, and their LODFs ------------------
function connected_without(c)
    parent = collect(1:N)
    root(x) = parent[x] == x ? x : (parent[x] = root(parent[x]))
    for l in 1:L
        l == c && continue
        a, b = root(u[l]), root(v[l])
        a == b || (parent[a] = b)
    end
    r = root(1)
    return all(root(i) == r for i in 1:N)
end
outages = [c for c in 1:L if connected_without(c)]
LODF = zeros(L, L)
for c in outages
    a = H[:, u[c]] .- H[:, v[c]]
    LODF[:, c] = a ./ (1 - a[c])
    LODF[c, c] = -1.0
end
postrow(l, c) = H[l, :] .+ LODF[l, c] .* H[c, :]
say(@sprintf("outages: %d of %d lines keep the network connected", length(outages), L))

# LODF check against the PTDF of the network with the line removed
let worst = 0.0
    for c in outages[1:min(10, end)]
        Dx2 = copy(base.Dx); Dx2[c] = 0.0
        base2 = typeof(base)((f === :Dx ? Dx2 : getfield(base, f) for f in fieldnames(typeof(base)))...)
        Hc = CT.ptdf(base2)
        worst = max(worst, maximum(abs.(Hc .- (H .+ LODF[:, c] * transpose(H[c, :])))))
    end
    say(@sprintf("check    LODF vs direct PTDF on 10 outages: worst gap %.2e", worst))
    worst < 1e-8 || error("LODF check failed")
end

env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
function lp()
    m = Model(() -> Gurobi.Optimizer(env))
    set_silent(m)
    set_optimizer_attribute(m, "Threads", 1)
    set_optimizer_attribute(m, "FeasibilityTol", 1e-9)
    set_optimizer_attribute(m, "OptimalityTol", 1e-9)
    return m
end

# ---- 1. full SC-DCOPF per hour ----------------------------------------------------
# Rows are true constraints at every hour, so the pool is shared and only the
# right-hand sides move with the demand.
function injection(g, d)
    p = -copy(d)
    for k in 1:K
        p[gb[k]] += g[k]
    end
    return p
end
G = H[:, gb]
function opf_model()
    m = lp()
    @variable(m, pmin[k] <= g[k=1:K] <= pmax[k])
    bal = @constraint(m, sum(g) == 0.0)
    up = @constraint(m, [l=1:L],  sum(G[l, k] * g[k] for k in 1:K) <= 0.0)
    dn = @constraint(m, [l=1:L], -sum(G[l, k] * g[k] for k in 1:K) <= 0.0)
    @objective(m, Min, sum(c1[k] * g[k] for k in 1:K) + c0)
    return (; m, g, bal, up, dn, rows=Dict{Tuple{Int,Int},Any}())
end
function set_demand!(M, d)
    set_normalized_rhs(M.bal, sum(d))
    hd = H * d
    for l in 1:L
        set_normalized_rhs(M.up[l], F[l] + hd[l])
        set_normalized_rhs(M.dn[l], F[l] - hd[l])
    end
    for ((l, c), r) in M.rows
        w = dot(r.row, d)
        set_normalized_rhs(r.up, F[l] + w)
        set_normalized_rhs(r.dn, F[l] - w)
    end
end
function add_row!(M, l, c, d)
    row = postrow(l, c)
    a = row[gb]
    w = dot(row, d)
    up = @constraint(M.m,  sum(a[k] * M.g[k] for k in 1:K) <= F[l] + w)
    dn = @constraint(M.m, -sum(a[k] * M.g[k] for k in 1:K) <= F[l] - w)
    M.rows[(l, c)] = (; row, up, dn)
end

base_opf, scopf = opf_model(), opf_model()
Z0, Zsc, rounds = fill(NaN, T), fill(NaN, T), zeros(Int, T)
t0 = time()
for s in 1:T
    d = D[:, s]
    set_demand!(base_opf, d)
    optimize!(base_opf.m)
    termination_status(base_opf.m) == MOI.OPTIMAL || (say("hour $s: base OPF ", termination_status(base_opf.m)); continue)
    Z0[s] = objective_value(base_opf.m)
    set_demand!(scopf, d)
    while true
        optimize!(scopf.m)
        rounds[s] += 1
        termination_status(scopf.m) == MOI.OPTIMAL || (say("hour $s: SC-DCOPF ", termination_status(scopf.m)); break)
        f = H * injection(value.(scopf.g), d)
        added = 0
        for c in outages, l in 1:L
            (l == c || haskey(scopf.rows, (l, c))) && continue
            abs(f[l] + LODF[l, c] * f[c]) > F[l] * (1 + 1e-7) || continue
            add_row!(scopf, l, c, d)
            added += 1
        end
        if added == 0
            Zsc[s] = objective_value(scopf.m)
            break
        end
    end
end
ok = findall(isfinite, Zsc)
secpct = 100 .* (Zsc .- Z0) ./ abs.(Z0)
say(@sprintf("SC-DCOPF: %d of %d hours feasible, %d contingency rows in the pool, %.0fs",
             length(ok), T, length(scopf.rows), time() - t0))
isempty(ok) && error("no feasible hour")
say(@sprintf("          security cost over the base OPF: mean %.3f%%, max %.3f%%",
             mean(secpct[ok]), maximum(secpct[ok])))
open(joinpath(OUT, "hours.csv"), "w") do io
    println(io, "hour,hour_id,base_opf_cost,sc_dcopf_cost,security_cost_pct,rounds")
    for s in 1:T
        @printf(io, "%d,%d,%.6f,%.6f,%.6f,%d\n", s, hour_ids[s], Z0[s], Zsc[s], secpct[s], rounds[s])
    end
end

# ---- 2. screening per band ------------------------------------------------------
Dk, Zk = D[:, ok], Zsc[ok]
Tk = length(ok)
tot = vec(sum(Dk; dims=1))
HD = H * Dk
function adversary_set(x; intact)
    m = lp()
    @variable(m, pmin[k] <= g[k=1:K] <= pmax[k])
    @variable(m, lam[1:Tk] >= 0)
    @constraint(m, sum(lam) == 1)
    @constraint(m, sum(g) == sum(tot[t] * lam[t] for t in 1:Tk))
    isfinite(x) && @constraint(m, sum(c1[k] * g[k] for k in 1:K) + c0 <=
                                  (1 + x) * sum(Zk[t] * lam[t] for t in 1:Tk))
    if intact
        @constraint(m, [l=1:L],  sum(G[l, k] * g[k] for k in 1:K) - sum(HD[l, t] * lam[t] for t in 1:Tk) <= F[l])
        @constraint(m, [l=1:L], -sum(G[l, k] * g[k] for k in 1:K) + sum(HD[l, t] * lam[t] for t in 1:Tk) <= F[l])
    end
    return (; m, g, lam, lps=Ref(0))
end
"max of sense * row' p over the set"
function reach(A, row, sense)
    a, w = row[gb], transpose(Dk) * row
    @objective(A.m, Max, sense * (sum(a[k] * A.g[k] for k in 1:K) - sum(w[t] * A.lam[t] for t in 1:Tk)))
    optimize!(A.m)
    A.lps[] += 1
    termination_status(A.m) == MOI.OPTIMAL || error("adversary LP ended $(termination_status(A.m))")
    return objective_value(A.m)
end

pairs_io = open(joinpath(OUT, "pairs.csv"), "w")
println(pairs_io, "band_pct,network,line,sign,loading")
nets_io = open(joinpath(OUT, "networks.csv"), "w")
println(nets_io, "band_pct,network,critical_lines,worst_checked_loading,redundant")
bandname(x) = isfinite(x) ? @sprintf("%g", 100x) : "inf"

for x in BANDS
    t1 = time()
    # base network
    A0 = adversary_set(x; intact=false)
    crit0 = Int[]; worst0 = 0.0
    for l in 1:L, sg in (1, -1)
        r = reach(A0, H[l, :], sg) / F[l]
        worst0 = max(worst0, r)
        if r > 1 + TAU
            println(pairs_io, "$(bandname(x)),base,$l,$sg,$r")
            l in crit0 || push!(crit0, l)
        end
    end
    println(nets_io, "$(bandname(x)),base,$(length(crit0)),$worst0,$(isempty(crit0))")
    # no band and no intact network is the certificate's injection set (2 lines in step 0)
    isinf(x) && say("check    base network with no band: $(length(crit0)) critical lines (step 0 found 2)")

    # outages: bound each post-outage flow by the pre-outage extremes first
    Ac = adversary_set(x; intact=true)
    hi = [reach(Ac, H[l, :], 1) for l in 1:L]
    lo = [-reach(Ac, H[l, :], -1) for l in 1:L]
    ncrit = Int[]; exact = 0
    for c in outages
        crit = Int[]; worst = 0.0
        for l in 1:L
            l == c && continue
            q = LODF[l, c]
            bhi = hi[l] + (q >= 0 ? q * hi[c] : q * lo[c])
            blo = lo[l] + (q >= 0 ? q * lo[c] : q * hi[c])
            # worst is over the pairs the bound could not clear
            for (sg, b) in ((1, bhi), (-1, -blo))
                b / F[l] > 1 + TAU || continue
                exact += 1
                r = reach(Ac, postrow(l, c), sg) / F[l]
                worst = max(worst, r)
                if r > 1 + TAU
                    println(pairs_io, "$(bandname(x)),$c,$l,$sg,$r")
                    l in crit || push!(crit, l)
                end
            end
        end
        push!(ncrit, length(crit))
        println(nets_io, "$(bandname(x)),$c,$(length(crit)),$worst,$(isempty(crit))")
    end
    flush(pairs_io); flush(nets_io)
    nonred = filter(>(0), ncrit)
    say(@sprintf("band %-4s base: %d critical lines (worst %.3f) | outages: %d of %d redundant; the rest %s critical lines (mean %.1f, max %d) | %d exact LPs, %.0fs",
                 bandname(x) * (isfinite(x) ? "%" : ""), length(crit0), worst0,
                 count(==(0), ncrit), length(outages),
                 isempty(nonred) ? "have no" : "have",
                 isempty(nonred) ? 0.0 : mean(nonred), isempty(nonred) ? 0 : maximum(nonred),
                 exact, time() - t1))
end
close(pairs_io); close(nets_io)
say("wrote ", OUT)
close(log_io)

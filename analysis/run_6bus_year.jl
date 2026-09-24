# Full-year hourly DC-OPF on the 6-bus Wood & Wollenberg case study.
#
#   julia --project=. --startup-file=no analysis/run_6bus_year.jl [pmin=relaxed|enforced]
#
# Reads the CSVs written by analysis/extract_6bus_hourly.py. Wind and solar are
# already folded into the net load there, so this file sees load only.
#
# Every hour is solved with unserved energy and curtailment variables priced
# far above any generator, so the LP never goes infeasible and the hours that
# WOULD have been infeasible are reported by their shortfall instead of by a
# status code. That is what makes a year's worth of infeasibility measurable.

using JuMP, Gurobi, Printf, DelimitedFiles, Statistics

const ROOT = dirname(@__DIR__)
const DATA = joinpath(ROOT, "case studies", "6busww")
const OUT  = joinpath(ROOT, "outputs", "6busww")
const BASEMVA = 100.0
const VOLL = 10_000.0          # $/MWh on unserved energy
const DUMP = 1_000.0           # $/MWh on curtailment (Pmin floor relief)
const TOL  = 1e-6

const DATA_SF = 0.08           # the Scaling Factor the Nodes sheet actually carries

global pmin_mode = :relaxed
global scale = DATA_SF         # demand-side scaling factor to solve at
global renew = true            # false = drop wind/solar too, so scale=0 is a dead system
for a in ARGS
    k, v = split(a, Char(61); limit=2)
    if k == "pmin"
        global pmin_mode = Symbol(v)
    elseif k == "scale"
        global scale = parse(Float64, v)
    elseif k == "renewables"
        global renew = parse(Bool, v)
    else
        error("accepts pmin=relaxed|enforced, scale=<float>, renewables=true|false")
    end
end
pmin_mode in (:relaxed, :enforced) || error("pmin must be relaxed or enforced")

readcsv(f) = readdlm(joinpath(DATA, f), ',', Any, '\n'; header=true)

gens_raw, _ = readcsv("network_gens.csv")
lines_raw, _ = readcsv("network_lines.csv")
load_raw, load_hdr = readcsv("hourly_net_load.csv")

gbus  = Int.(gens_raw[:, 2]);  gpmax = Float64.(gens_raw[:, 3])
gpmin_data = Float64.(gens_raw[:, 4]);  gcost = Float64.(gens_raw[:, 5])
gpmin = pmin_mode === :relaxed ? zeros(length(gpmin_data)) : gpmin_data

lfrom = Int.(lines_raw[:, 2]);  lto = Int.(lines_raw[:, 3])
lx    = Float64.(lines_raw[:, 4]);  lfmax = Float64.(lines_raw[:, 5])

buses = sort(unique(vcat(lfrom, lto)))
N, L, K = length(buses), length(lfrom), length(gbus)
idx = Dict(b => i for (i, b) in enumerate(buses))
slack = 1                                   # bus 1 is the Type 3 bus
suscept = BASEMVA ./ lx                     # MW per radian

# hourly net load, MW, in bus order. The stored series already has the data's
# own 0.08 applied; solving at another scaling rebuilds it from the components,
# because only the DEMAND side scales -- wind and solar are fixed nameplate MW.
stamp = Int.(load_raw[:, 2:5])
Dmat = Float64.(load_raw[:, 6:end])
Renew = Float64.(readcsv("hourly_wind.csv")[1][:, 6:end]) .+
        Float64.(readcsv("hourly_solar.csv")[1][:, 6:end])
if scale != DATA_SF || !renew
    dem = Float64.(readcsv("hourly_demand.csv")[1][:, 6:end]) .* (scale / DATA_SF)
    Dmat = renew ? dem .- Renew : dem
    renew || (Renew = zero(Renew))
end
H = size(Dmat, 1)
size(Dmat, 2) == N || error("net-load CSV has $(size(Dmat,2)) bus columns, network has $N")

@printf("6-bus W&W: %d buses, %d lines, %d generators, %d hours   pmin=%s\n",
        N, L, K, H, pmin_mode)
@printf("  gen capacity %.1f MW   net load %.1f .. %.1f MW (mean %.1f)\n",
        sum(gpmax), minimum(sum(Dmat, dims=2)), maximum(sum(Dmat, dims=2)),
        mean(sum(Dmat, dims=2)))

# --------------------------------------------------------------------------- #
# One model, re-solved 8760 times with only the load RHS moving. Rebuilding it
# per hour costs more in JuMP overhead than the LP itself takes to solve.
# --------------------------------------------------------------------------- #
env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
m = Model(() -> Gurobi.Optimizer(env))
set_optimizer_attribute(m, "OutputFlag", 0)
set_silent(m)

@variable(m, g[k=1:K])
@variable(m, th[1:N])
@variable(m, f[l=1:L])
@variable(m, shed[1:N] >= 0)
@variable(m, curt[1:N] >= 0)

@constraint(m, [k=1:K], gpmin[k] <= g[k] <= gpmax[k])
@constraint(m, th[slack] == 0)
@constraint(m, flowdef[l=1:L], f[l] == suscept[l] * (th[idx[lfrom[l]]] - th[idx[lto[l]]]))
@constraint(m, fup[l=1:L],  f[l] <=  lfmax[l])
@constraint(m, flo[l=1:L], -f[l] <=  lfmax[l])

gens_at = [findall(==(buses[i]), gbus) for i in 1:N]
gpmin_bus = [sum(gpmin[k] for k in findall(==(buses[i]), gbus); init=0.0) for i in 1:N]
out_of  = [findall(==(buses[i]), lfrom) for i in 1:N]
in_to   = [findall(==(buses[i]), lto)   for i in 1:N]

# LMP = dual of this row, so the load sits alone on the RHS.
@constraint(m, balance[i=1:N],
    sum(g[k] for k in gens_at[i]; init=0.0) + shed[i] - curt[i]
    - sum(f[l] for l in out_of[i]; init=0.0)
    + sum(f[l] for l in in_to[i];  init=0.0) == 0.0)

@objective(m, Min, sum(gcost[k] * g[k] for k in 1:K) +
                   VOLL * sum(shed) + DUMP * sum(curt))

# --------------------------------------------------------------------------- #
hours_cost   = zeros(H)
hours_shed   = zeros(H)
hours_curt   = zeros(H)
hours_spread = zeros(H)
lmp          = zeros(H, N)
flow         = zeros(H, L)
disp         = zeros(H, K)
bind_hi      = zeros(Int, L)
bind_lo      = zeros(Int, L)
bind_ok      = zeros(Int, L)   # binding in hours that were actually servable
cong_rent    = zeros(L)        # rent over servable hours only: during a
                               # shortfall the price is VOLL, not a real
                               # congestion price, and it would swamp the sum
nfail        = 0

t0 = time()
for h in 1:H
    for i in 1:N
        d = Dmat[h, i]
        set_normalized_rhs(balance[i], d)
        set_upper_bound(shed[i], max(d, 0.0))
        # Curtailment spills what is produced AT this bus: local renewable
        # output, plus the must-run floor when Pmin is enforced. Unbounded, the
        # LP dumps wherever the network reaches and the location is meaningless.
        set_upper_bound(curt[i], Renew[h, i] + gpmin_bus[i])
    end
    optimize!(m)
    if termination_status(m) != MOI.OPTIMAL
        global nfail += 1
        continue
    end
    hours_cost[h] = objective_value(m)
    sv, cv = value.(shed), value.(curt)
    hours_shed[h] = sum(sv)
    hours_curt[h] = sum(cv)
    fv = value.(f)
    flow[h, :] = fv
    disp[h, :] = value.(g)
    for i in 1:N
        lmp[h, i] = dual(balance[i])
    end
    hours_spread[h] = maximum(view(lmp, h, :)) - minimum(view(lmp, h, :))
    servable = hours_shed[h] <= TOL
    for l in 1:L
        atlimit = false
        if fv[l] >= lfmax[l] - 1e-5
            bind_hi[l] += 1; atlimit = true
        elseif fv[l] <= -lfmax[l] + 1e-5
            bind_lo[l] += 1; atlimit = true
        end
        if servable
            atlimit && (bind_ok[l] += 1)
            cong_rent[l] += (abs(dual(fup[l])) + abs(dual(flo[l]))) * lfmax[l]
        end
    end
    h % 2000 == 0 && @printf("  ... %d / %d hours (%.1fs)\n", h, H, time() - t0)
end
@printf("solved %d hours in %.1f s   (%d solver failures)\n\n", H, time() - t0, nfail)

# --------------------------------------------------------------------------- #
mkpath(OUT)
tag = (scale == DATA_SF && renew) ? string(pmin_mode) :
      "$(pmin_mode)_scale$(scale)" * (renew ? "" : "_norenew")

shed_hours = findall(>(TOL), hours_shed)
curt_hours = findall(>(TOL), hours_curt)
util = [maximum(abs.(flow[h, :]) ./ lfmax) for h in 1:H]

open(joinpath(OUT, "report_$(tag).txt"), "w") do io
    for sink in (stdout, io)
        println(sink, "="^74)
        println(sink, "6-bus W&W -- full-year hourly DC-OPF   (pmin $pmin_mode)")
        println(sink, "="^74)
        @printf(sink, "\nhours %d   buses %d   lines %d   generators %d\n", H, N, L, K)
        tl = vec(sum(Dmat, dims=2))
        @printf(sink, "net load  min %.1f  mean %.1f  max %.1f MW\n",
                minimum(tl), mean(tl), maximum(tl))
        @printf(sink, "gen Pmax %.1f MW   Pmin(data) %.1f MW\n\n",
                sum(gpmax), sum(gpmin_data))

        println(sink, "-- FEASIBILITY " * "-"^58)
        @printf(sink, "hours with unserved energy : %4d  (%.2f%% of the year)\n",
                length(shed_hours), 100 * length(shed_hours) / H)
        if !isempty(shed_hours)
            @printf(sink, "  total unserved          : %.1f MWh\n", sum(hours_shed))
            @printf(sink, "  worst hour              : %.1f MW (hour %d, %d-%02d-%02d h%02d)\n",
                    maximum(hours_shed), argmax(hours_shed),
                    stamp[argmax(hours_shed), 1], stamp[argmax(hours_shed), 2],
                    stamp[argmax(hours_shed), 3], stamp[argmax(hours_shed), 4])
            bym = zeros(Int, 12)
            for h in shed_hours; bym[stamp[h, 2]] += 1; end
            println(sink, "  by month                : ", join(bym, " "))
        end
        @printf(sink, "hours with curtailment     : %4d  (%.2f%%)\n",
                length(curt_hours), 100 * length(curt_hours) / H)
        isempty(curt_hours) || @printf(sink, "  total curtailed         : %.1f MWh\n",
                                       sum(hours_curt))

        nok = H - length(shed_hours)
        println(sink, "\n-- CONGESTION " * "-"^59)
        println(sink, "hrs+/hrs- count every hour; the last two columns cover only")
        @printf(sink, "the %d servable hours, where the price is a real congestion price.\n\n", nok)
        @printf(sink, "%-5s %-7s %6s %7s %7s %8s %9s %12s\n",
                "line", "path", "Fmax", "hrs +", "hrs -", "%% of yr",
                "bind(ok)", "rent \$")
        order = sortperm(bind_hi .+ bind_lo; rev=true)
        for l in order
            nb = bind_hi[l] + bind_lo[l]
            @printf(sink, "%-5d %-7s %6.0f %7d %7d %7.1f%% %9d %12.0f\n",
                    l, "$(lfrom[l])-$(lto[l])", lfmax[l], bind_hi[l], bind_lo[l],
                    100 * nb / H, bind_ok[l], cong_rent[l])
        end
        ncong = count(h -> util[h] >= 1 - 1e-5, 1:H)
        @printf(sink, "\nhours with >=1 line at its limit : %d  (%.1f%%)\n",
                ncong, 100 * ncong / H)
        @printf(sink, "mean / max worst-line utilisation : %.3f / %.3f\n",
                mean(util), maximum(util))
        never = [l for l in 1:L if bind_hi[l] + bind_lo[l] == 0]
        isempty(never) || @printf(sink, "never binding                     : %s\n",
            join(["$(lfrom[l])-$(lto[l])" for l in never], ", "))

        println(sink, "\n-- PRICES " * "-"^63)
        ok = setdiff(1:H, shed_hours)          # VOLL hours would swamp the stats
        @printf(sink, "LMP stats over the %d hours without unserved energy\n", length(ok))
        for i in 1:N
            col = lmp[ok, i]
            @printf(sink, "  bus %d   min %8.2f   mean %8.2f   max %8.2f  \$/MWh\n",
                    buses[i], minimum(col), mean(col), maximum(col))
        end
        sep = count(h -> hours_spread[h] > 1e-4, ok)
        @printf(sink, "hours with price separation (>1e-4) : %d  (%.1f%% of those hours)\n",
                sep, 100 * sep / length(ok))
        @printf(sink, "mean / max LMP spread               : %.2f / %.2f \$/MWh\n",
                mean(hours_spread[ok]), maximum(hours_spread[ok]))

        println(sink, "\n-- GENERATION " * "-"^59)
        @printf(sink, "%-5s %-5s %8s %10s %10s %8s\n",
                "gen", "bus", "Pmax", "MWh", "CF", "\$/MWh")
        for k in 1:K
            e = sum(view(disp, :, k))
            @printf(sink, "%-5d %-5d %8.0f %10.0f %9.1f%% %8.2f\n",
                    k, gbus[k], gpmax[k], e, 100 * e / (gpmax[k] * H), gcost[k])
        end
        @printf(sink, "\nannual cost (incl. penalties) : \$%.3g\n", sum(hours_cost))
        energy_cost = sum(gcost[k] * sum(view(disp, :, k)) for k in 1:K)
        @printf(sink, "annual generation cost         : \$%.3g\n", energy_cost)
        println(sink)
    end
end

writedlm(joinpath(OUT, "hourly_$(tag).csv"),
    vcat(permutedims(vcat(["hour", "month", "day", "hour_of_day", "cost",
                           "shed_mw", "curt_mw", "lmp_spread", "worst_util"],
                          ["lmp_bus$(b)" for b in buses],
                          ["f_$(lfrom[l])_$(lto[l])" for l in 1:L])),
         hcat(1:H, stamp[:, 2], stamp[:, 3], stamp[:, 4], hours_cost,
              hours_shed, hours_curt, hours_spread, util, lmp, flow)), ',')

writedlm(joinpath(OUT, "lines_$(tag).csv"),
    vcat(permutedims(["line", "from", "to", "fmax", "hours_at_pos_limit",
                      "hours_at_neg_limit", "pct_of_year", "congestion_rent"]),
         hcat(1:L, lfrom, lto, lfmax, bind_hi, bind_lo,
              100 .* (bind_hi .+ bind_lo) ./ H, cong_rent)), ',')

println("wrote report_$(tag).txt, hourly_$(tag).csv, lines_$(tag).csv -> $OUT")

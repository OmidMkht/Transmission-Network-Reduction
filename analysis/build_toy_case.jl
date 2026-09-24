# Build a reduction toy case from an extracted planning workbook: a MATPOWER
# file plus the hourly scenario matrices the :hourly demand path expects.
#
#   julia --project=. --startup-file=no analysis/build_toy_case.jl case=6busww scale=0.018
#   julia --project=. --startup-file=no analysis/build_toy_case.jl case=14bus  scale=0.05
#
# Writes
#   case studies/<case>/<case>_toy.m           network, linear costs, base load
#   outputs/<caseid>_dcopf/{bus_ids,load_mw,generation_mw,scenario_summary}.csv
#
# Only the EXISTING network is used: network_lines.csv already excludes the
# IntI = 0 planning candidates, and network_gens.csv the candidate units. Wind
# and solar stay folded into the load as negative injection, because the
# reduction models dispatch the thermal units only.
#
# load_mw is the SERVED net load (net - shed + curtailment) so generation
# balances it exactly: build_multiscenario_tx_case rejects any imbalance, and
# the reduction MILP is infeasible for every clustering when sum(p) != 0.

using JuMP, Gurobi, Printf, DelimitedFiles, Statistics

const ROOT = dirname(@__DIR__)
const BASEMVA = 100.0
const VOLL, DUMP, TOL = 10_000.0, 1_000.0, 1e-6
const CASEID = Dict("6busww" => "case6ww", "14bus" => "case14toy")

casename = "6busww"
scale = 0.018
for a in ARGS
    k, v = split(a, Char(61); limit=2)
    if k == "case"; global casename = String(v)
    elseif k == "scale"; global scale = parse(Float64, v)
    else; error("accepts case=6busww|14bus and scale=<float>")
    end
end
haskey(CASEID, casename) || error("unknown case $casename")

DATA = joinpath(ROOT, "case studies", casename)
MATDIR = joinpath(ROOT, "outputs", "$(CASEID[casename])_dcopf")

rd(f) = readdlm(joinpath(DATA, f), ',', Any, '\n'; header=true)[1]
G, Lr, Bs = rd("network_gens.csv"), rd("network_lines.csv"), rd("network_buses.csv")
gbus = Int.(G[:, 2]); gpmax = Float64.(G[:, 3]); gcost = Float64.(G[:, 5])
lfrom = Int.(Lr[:, 2]); lto = Int.(Lr[:, 3])
lx = Float64.(Lr[:, 4]); lfmax = Float64.(Lr[:, 5])
allbus = Int.(Bs[:, 1]); btype = Int.(Bs[:, 2])

raw = rd("hourly_demand.csv")
stamp = Int.(raw[:, 2:5])
# the stored demand already carries the workbook's own scaling factor
data_sf = maximum(Float64.(Bs[:, 6]))
dem = Float64.(raw[:, 6:end]) .* (scale / data_sf)
wnd = Float64.(rd("hourly_wind.csv")[:, 6:end])
sol = Float64.(rd("hourly_solar.csv")[:, 6:end])
ren = wnd .+ sol
D = dem .- ren
H, N = size(D); L, K = length(lfrom), length(gbus)
sus = BASEMVA ./ lx
slack = something(findfirst(==(3), btype), 1)

@printf("[%s] scale %.4g (workbook %.4g)   %d buses, %d lines, %d gens, %d hours\n",
        casename, scale, data_sf, N, L, K, H)
@printf("   net load %.1f .. %.1f MW (mean %.1f)   gen Pmax %.1f MW\n",
        minimum(sum(D, dims=2)), maximum(sum(D, dims=2)), mean(sum(D, dims=2)), sum(gpmax))

env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
m = Model(() -> Gurobi.Optimizer(env)); set_silent(m)
@variable(m, 0 <= g[k=1:K] <= gpmax[k])          # Pmin relaxed, as the repo does
@variable(m, th[1:N]); @variable(m, f[1:L])
@variable(m, shed[1:N] >= 0); @variable(m, curt[1:N] >= 0)
@constraint(m, th[slack] == 0)
@constraint(m, [l=1:L], f[l] == sus[l] * (th[lfrom[l]] - th[lto[l]]))
@constraint(m, [l=1:L], -lfmax[l] <= f[l] <= lfmax[l])
@constraint(m, bal[i=1:N],
    sum(g[k] for k in 1:K if gbus[k] == allbus[i]; init=0.0) + shed[i] - curt[i]
    - sum(f[l] for l in 1:L if lfrom[l] == allbus[i]; init=0.0)
    + sum(f[l] for l in 1:L if lto[l]   == allbus[i]; init=0.0) == 0.0)
@objective(m, Min, sum(gcost[k]*g[k] for k in 1:K) + VOLL*sum(shed) + DUMP*sum(curt))

genmat = zeros(N, H); loadmat = zeros(N, H)
nshed = 0; ncurt = 0
for h in 1:H
    for i in 1:N
        set_normalized_rhs(bal[i], D[h, i])
        set_upper_bound(shed[i], max(D[h, i], 0.0))
        set_upper_bound(curt[i], ren[h, i])   # spill only what is produced here
    end
    optimize!(m)
    termination_status(m) == MOI.OPTIMAL || error("hour $h did not solve")
    gv, sv, cv = value.(g), value.(shed), value.(curt)
    sum(sv) > TOL && (global nshed += 1)
    sum(cv) > TOL && (global ncurt += 1)
    for k in 1:K
        genmat[findfirst(==(gbus[k]), allbus), h] += gv[k]
    end
    loadmat[:, h] = D[h, :] .- sv .+ cv
end
imb = maximum(abs.(vec(sum(genmat, dims=1) .- sum(loadmat, dims=1))))
@printf("   solved %d hours   shed in %d, curtailment in %d   imbalance %.2e MW\n",
        H, nshed, ncurt, imb)
imb < 1e-6 || error("generation does not balance served load")

mkpath(MATDIR)
writedlm(joinpath(MATDIR, "bus_ids.csv"), allbus, ',')
writedlm(joinpath(MATDIR, "load_mw.csv"), loadmat, ',')
writedlm(joinpath(MATDIR, "generation_mw.csv"), genmat, ',')
open(joinpath(MATDIR, "scenario_summary.csv"), "w") do io
    println(io, "scenario_id,month,day,hour")
    for h in 1:H
        @printf(io, "%d,%d,%d,%d\n", h, stamp[h, 2], stamp[h, 3], stamp[h, 4])
    end
end

basePd = vec(mean(loadmat, dims=2))
mfile = joinpath(DATA, "$(casename)_toy.m")
open(mfile, "w") do io
    println(io, "function mpc = $(casename)_toy")
    println(io, "%$(uppercase(casename))_TOY  reduction toy case, existing network only.")
    println(io, "%   Generated by analysis/build_toy_case.jl from the planning workbook.")
    @printf(io, "%%   Demand scaling %.4g; wind and solar are negative load.\n", scale)
    println(io, "%   Planning-candidate lines and generators are NOT included.")
    println(io, "%   Pmin is 0: it is a commitment floor, every solver here relaxes it,")
    println(io, "%   and build_tx_case's base DC-OPF does not. Real floor in network_gens.csv.")
    println(io, "mpc.version = '2';")
    println(io, "mpc.baseMVA = 100;\n")
    println(io, "%% bus data")
    println(io, "%\tbus_i\ttype\tPd\tQd\tGs\tBs\tarea\tVm\tVa\tbaseKV\tzone\tVmax\tVmin")
    println(io, "mpc.bus = [")
    for i in 1:N
        @printf(io, "\t%d\t%d\t%.4f\t0\t0\t0\t1\t1\t0\t230\t1\t1.1\t0.9;\n",
                allbus[i], btype[i], basePd[i])
    end
    println(io, "];\n")
    println(io, "%% generator data")
    println(io, "%\tbus\tPg\tQg\tQmax\tQmin\tVg\tmBase\tstatus\tPmax\tPmin\tPc1\tPc2\tQc1min\tQc1max\tQc2min\tQc2max\tramp_agc\tramp_10\tramp_30\tramp_q\tapf")
    println(io, "mpc.gen = [")
    for k in 1:K
        @printf(io, "\t%d\t0\t0\t0\t0\t1\t100\t1\t%.4f\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;\n",
                gbus[k], gpmax[k])
    end
    println(io, "];\n")
    println(io, "%% branch data")
    println(io, "%\tfbus\ttbus\tr\tx\tb\trateA\trateB\trateC\tratio\tangle\tstatus\tangmin\tangmax")
    println(io, "mpc.branch = [")
    for l in 1:L
        @printf(io, "\t%d\t%d\t0\t%.6f\t0\t%.4f\t%.4f\t%.4f\t0\t0\t1\t-360\t360;\n",
                lfrom[l], lto[l], lx[l], lfmax[l], lfmax[l], lfmax[l])
    end
    println(io, "];\n")
    println(io, "%% generator cost data")
    println(io, "mpc.gencost = [")
    for k in 1:K
        @printf(io, "\t2\t0\t0\t3\t0\t%.6f\t0;\n", gcost[k])
    end
    println(io, "];")
end
@printf("   wrote %s and 4 matrices -> %s\n", basename(mfile), MATDIR)

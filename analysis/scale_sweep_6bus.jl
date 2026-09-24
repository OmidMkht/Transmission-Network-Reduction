# How the 6-bus system behaves as the Nodes-sheet Scaling Factor is varied.
#
#   julia --project=. --startup-file=no analysis/scale_sweep_6bus.jl
#
# The data ships Scaling Factor = 0.08. Only the DEMAND side scales with it;
# wind and solar are fixed nameplate MW, so net load = demand*s/0.08 - wind - solar.
# At s = 0 the net load is pure renewable injection and every price collapses to
# the curtailment penalty, which is why that point reports nothing useful.

using JuMP, Gurobi, Printf, DelimitedFiles, Statistics

const ROOT = dirname(@__DIR__)
const DATA = joinpath(ROOT, "case studies", "6busww")
const OUT  = joinpath(ROOT, "outputs", "6busww")
const BASEMVA = 100.0
const VOLL, DUMP, TOL = 10_000.0, 1_000.0, 1e-6
const DATA_SF = 0.08

rd(f) = readdlm(joinpath(DATA, f), ',', Any, '\n'; header=true)[1]
G, Lr = rd("network_gens.csv"), rd("network_lines.csv")
gbus = Int.(G[:, 2]); gpmax = Float64.(G[:, 3]); gcost = Float64.(G[:, 5])
lfrom = Int.(Lr[:, 2]); lto = Int.(Lr[:, 3])
lfmax = Float64.(Lr[:, 5]); sus = BASEMVA ./ Float64.(Lr[:, 4])

dem = Float64.(rd("hourly_demand.csv")[:, 6:end])
wnd = Float64.(rd("hourly_wind.csv")[:, 6:end])
sol = Float64.(rd("hourly_solar.csv")[:, 6:end])
H, N = size(dem); L, K = length(lfrom), length(gbus)

env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
m = Model(() -> Gurobi.Optimizer(env)); set_silent(m)
@variable(m, 0 <= g[k=1:K] <= gpmax[k])
@variable(m, th[1:N]); @variable(m, f[1:L])
@variable(m, shed[1:N] >= 0); @variable(m, curt[1:N] >= 0)
@constraint(m, th[1] == 0)
@constraint(m, [l=1:L], f[l] == sus[l] * (th[lfrom[l]] - th[lto[l]]))
@constraint(m, fup[l=1:L],  f[l] <=  lfmax[l])
@constraint(m, flo[l=1:L], -f[l] <=  lfmax[l])
@constraint(m, bal[i=1:N],
    sum(g[k] for k in 1:K if gbus[k] == i; init=0.0) + shed[i] - curt[i]
    - sum(f[l] for l in 1:L if lfrom[l] == i; init=0.0)
    + sum(f[l] for l in 1:L if lto[l]   == i; init=0.0) == 0.0)
@objective(m, Min, sum(gcost[k] * g[k] for k in 1:K) + VOLL*sum(shed) + DUMP*sum(curt))

function year(s)
    D = dem .* (s / DATA_SF) .- wnd .- sol
    shed_h = 0; shed_e = 0.0; curt_h = 0; cong_h = 0
    bind = zeros(Int, L); lmpsum = zeros(N); spread = 0.0; nok = 0; gen = zeros(K)
    for h in 1:H
        for i in 1:N
            set_normalized_rhs(bal[i], D[h, i])
            set_upper_bound(shed[i], max(D[h, i], 0.0))
            set_upper_bound(curt[i], sum(gpmax) + 500.0)
        end
        optimize!(m)
        termination_status(m) == MOI.OPTIMAL || continue
        sv, cv, fv = sum(value.(shed)), sum(value.(curt)), value.(f)
        sv > TOL && (shed_h += 1; shed_e += sv)
        cv > TOL && (curt_h += 1)
        gen .+= value.(g)
        any(l -> abs(fv[l]) >= lfmax[l] - 1e-5, 1:L) && (cong_h += 1)
        for l in 1:L
            abs(fv[l]) >= lfmax[l] - 1e-5 && (bind[l] += 1)
        end
        if sv <= TOL && cv <= TOL        # a price-meaningful hour
            lv = [dual(bal[i]) for i in 1:N]
            lmpsum .+= lv
            spread += maximum(lv) - minimum(lv)
            nok += 1
        end
    end
    tot = vec(sum(D, dims=2))
    return (peak=maximum(tot), mean=mean(tot), shed_h=shed_h, shed_e=shed_e,
            curt_h=curt_h, cong_h=cong_h, bind=bind, nok=nok,
            lmp=nok > 0 ? lmpsum ./ nok : fill(NaN, N),
            spread=nok > 0 ? spread / nok : NaN, gen=gen)
end

scales = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.10]
res = Dict{Float64,Any}()

open(joinpath(OUT, "scale_sweep.txt"), "w") do io
  for sink in (stdout, io)
    println(sink, "="^100)
    println(sink, "6-bus W&W -- full-year DC-OPF across the Nodes-sheet Scaling Factor")
    println(sink, "(data ships 0.08; only demand scales, wind/solar are fixed nameplate MW)")
    println(sink, "="^100)
    @printf(sink, "\n%-6s %8s %8s %8s %10s %8s %8s %10s\n",
            "scale", "peak MW", "mean MW", "shed hrs", "shed GWh",
            "curt hrs", "cong hrs", "spread \$")
    println(sink, "-"^76)
    for s in scales
      r = get!(res, s) do; year(s); end
      @printf(sink, "%-6.2f %8.1f %8.1f %8d %10.1f %8d %8d %10.2f\n",
              s, r.peak, r.mean, r.shed_h, r.shed_e/1000, r.curt_h, r.cong_h, r.spread)
    end

    println(sink, "\n\nMEAN LMP BY BUS  (\$/MWh, over hours with no shortfall and no curtailment)")
    @printf(sink, "%-6s %6s %9s %9s %9s %9s %9s %9s\n",
            "scale", "hrs", "bus1", "bus2", "bus3", "bus4", "bus5", "bus6")
    println(sink, "-"^76)
    for s in scales
      r = res[s]
      @printf(sink, "%-6.2f %6d", s, r.nok)
      for i in 1:N
        isnan(r.lmp[i]) ? @printf(sink, " %9s", "-") : @printf(sink, " %9.2f", r.lmp[i])
      end
      println(sink)
    end

    println(sink, "\n\nHOURS AT LIMIT, BY LINE  (out of 8760)")
    @printf(sink, "%-6s", "scale")
    for l in 1:L; @printf(sink, " %7s", "$(lfrom[l])-$(lto[l])"); end
    println(sink); println(sink, "-"^(6 + 8L))
    for s in scales
      @printf(sink, "%-6.2f", s)
      for l in 1:L; @printf(sink, " %7d", res[s].bind[l]); end
      println(sink)
    end

    println(sink, "\n\nGENERATION, GWh")
    @printf(sink, "%-6s %12s %12s %12s   %s\n", "scale",
            "g1 \$30.45", "g2 \$47.16", "g3 \$96.12", "(bus 1 / 2 / 3)")
    println(sink, "-"^60)
    for s in scales
      r = res[s]
      @printf(sink, "%-6.2f %12.0f %12.0f %12.0f\n", s,
              r.gen[1]/1000, r.gen[2]/1000, r.gen[3]/1000)
    end
    println(sink)
  end
end
println("wrote scale_sweep.txt -> $OUT")

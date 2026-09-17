# --------------------------------------------------------------------------- #
# How infeasible is a proxy reduction? Replays one proxy/run_proxy.jl result and
# measures it line by line.
#
#   julia --project=. --startup-file=no analysis/proxy_infeasibility.jl \
#         [casefile] [setting_dir]
#
#   casefile     name under "case studies/", default pglib_opf_case118_ieee.m
#   setting_dir  one <output>/<setting> folder written by proxy/run_proxy.jl, default
#                outputs/proxy/case118/single_eps0.1/conservative_0p1pct
#
# Uses the run's own assignment_matrix.csv, so the clustering is exactly the one
# that was reported, and checks both DC-OPF objectives against the recorded
# ones before trusting anything.
#
# Three dispatches are compared:
#   full    DC-OPF on the original network (the truth)
#   reduced DC-OPF on the reduced network
#   replay  the reduced dispatch pushed through the ORIGINAL network's physics
#
# Writes <setting_dir>/infeasibility/{summary.txt, lines.csv,
# congested_lines.csv, generators.csv, lmp.csv}.
# --------------------------------------------------------------------------- #

using LinearAlgebra, Printf, DelimitedFiles, Statistics

const ROOT = dirname(@__DIR__)
const CASE = length(ARGS) >= 1 ? ARGS[1] : "pglib_opf_case118_ieee.m"
const SETTING = length(ARGS) >= 2 ? abspath(ARGS[2]) :
    joinpath(ROOT, "outputs", "proxy", "case118", "single_eps0.1", "conservative_0p1pct")

# Same tolerances as the pipeline's validator, so the counts agree with it.
const BIND_TOL = 1e-3     # at its limit:  |f| >= (1 - tol) * rating
const VIOL_TOL = 1e-3     # violated:      overload > tol * rating

@eval module TR
    include(joinpath(dirname(@__DIR__), "common", "preprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "postprocessing.jl"))
end

c = TR.build_single_scenario_case(joinpath(ROOT, "case studies", CASE); time_limit=60.0)
base = c.base
N, L, G = base.N, base.Ln, length(base.gen_bus)
mva = base.baseMVA
bus = c.bus_ids

A = round.(Int, readdlm(joinpath(SETTING, "assignment_matrix.csv"), ','))
size(A) == (N, N) || error("assignment matrix is $(size(A)), expected ($N, $N)")
rep = TR.extract_reduction(A).rep_of
nret = length(TR.extract_reduction(A).retained)
internal = [rep[base.Efrom[l]] == rep[base.Eto[l]] for l in 1:L]

# Protected flags come from the run; the internal flags are recomputed and must agree.
protected = falses(L)
let path = joinpath(SETTING, "line_status.csv")
    if isfile(path)
        M, hdr = readdlm(path, ','; header=true)
        col(name) = findfirst(==(name), vec(hdr))
        for k in axes(M, 1)
            l = Int(M[k, col("line")])
            protected[l] = Int(M[k, col("protected")]) == 1
            (Int(M[k, col("internal")]) == 1) == internal[l] ||
                error("line $l: internal flag disagrees with line_status.csv")
        end
    end
end

# ---- solve -----------------------------------------------------------------
full = TR.solve_assigned_dc_opf_scenario(c, Matrix{Int}(I, N, N), 1; relax_pmin=true)
red = TR.solve_assigned_dc_opf_scenario(c, A, 1; relax_pmin=true)
full.optimal && red.optimal ||
    error("DC-OPF not optimal: full $(full.status), reduced $(red.status)")
replay = TR.check_dispatch_on_original_scenario(c, red.pG, 1; relative_tolerance=VIOL_TOL)
repair = TR.repair_dispatch_on_original_scenario(c, red.pG, 1; relax_pmin=true)

# Did we reproduce the reported run?
reproduced = let path = joinpath(SETTING, "monthly_dcopf_validation.csv")
    if isfile(path)
        M, hdr = readdlm(path, ','; header=true)
        col(name) = findfirst(==(name), vec(hdr))
        f0, r0 = Float64(M[1, col("original_objective")]), Float64(M[1, col("reduced_objective")])
        (abs(full.objective - f0) / abs(f0) < 1e-6 && abs(red.objective - r0) / abs(r0) < 1e-6,
         f0, r0)
    else
        (missing, NaN, NaN)
    end
end

# ---- line quantities -------------------------------------------------------
rate = base.frate
u_full = abs.(full.flow) ./ rate
u_red = abs.(red.flow) ./ rate
u_rep = abs.(replay.flow) ./ rate
over = replay.overload                         # p.u.
over_rel = replay.relative_overload            # fraction of rating

external = findall(.!internal)
bind_full = findall(u_full .>= 1 - BIND_TOL)
bind_red = [l for l in external if u_red[l] >= 1 - BIND_TOL]
bind_rep = findall(u_rep .>= 1 - BIND_TOL)
viol = findall(over_rel .> VIOL_TOL)
both = intersect(bind_full, bind_red)
missed = setdiff(bind_full, bind_red)          # at limit in truth, not in reduced model
spurious = setdiff(bind_red, bind_full)        # at limit in reduced model, not in truth

band(u, idx, t) = count(l -> u[l] >= t, idx)

# ---- generators and prices -------------------------------------------------
dpg = (red.pG .- full.pG) .* mva
moved = findall(abs.(dpg) .> 0.1)
load = c.load[:, 1]
lmp_err = abs.(red.lmp .- full.lmp)
lmp_lw = sum(load .* lmp_err) / sum(load)
levels(v) = length(unique(round.(v; digits=2)))

# ---- report ----------------------------------------------------------------
io = IOBuffer()
p(args...) = println(io, args...)
f(fmt, args...) = println(io, Printf.format(Printf.Format(fmt), args...))

p("=" ^ 78)
p("PROXY REDUCTION INFEASIBILITY  --  ", basename(CASE))
p("=" ^ 78)
p("setting    : ", SETTING)
f("network    : %d buses, %d lines, %d generators, baseMVA %.0f", N, L, G, mva)
f("reduction  : %d buses retained (%.1f%%), %d internal / %d external lines",
  nret, 100 * (N - nret) / N, count(internal), length(external))
f("protected  : %d lines", count(protected))
f("tolerances : at-limit |f| >= (1 - %.0e) rating;  violated overload > %.0e rating",
  BIND_TOL, VIOL_TOL)
if reproduced[1] === missing
    p("reproduced : no recorded objectives to compare against")
else
    f("reproduced : %s  (full %.4f vs recorded %.4f;  reduced %.4f vs recorded %.4f)",
      reproduced[1] ? "YES" : "NO -- numbers below may not match the run",
      full.objective, reproduced[2], red.objective, reproduced[3])
end

p()
p("1) CONGESTED LINES ON THE FULL NETWORK  (full DC-OPF)")
p("-" ^ 78)
f("   at limit                 = %d of %d lines", length(bind_full), L)
f("   loading >= 99%%           = %d", band(u_full, 1:L, 0.99))
f("   loading >= 90%%           = %d", band(u_full, 1:L, 0.90))
f("   loading >= 80%%           = %d", band(u_full, 1:L, 0.80))
f("   of the at-limit lines: %d internal (merged away), %d external, %d protected",
  count(l -> internal[l], bind_full), count(l -> !internal[l], bind_full),
  count(l -> protected[l], bind_full))

p()
p("2) CONGESTED LINES ON THE REDUCED NETWORK  (reduced DC-OPF, external lines)")
p("-" ^ 78)
f("   at limit                 = %d of %d external lines", length(bind_red), length(external))
f("   loading >= 99%%           = %d", band(u_red, external, 0.99))
f("   loading >= 90%%           = %d", band(u_red, external, 0.90))
f("   loading >= 80%%           = %d", band(u_red, external, 0.80))
f("   at limit in BOTH         = %d", length(both))
f("   MISSED  (full only)      = %d   truly congested, reduced model thinks slack", length(missed))
f("   SPURIOUS (reduced only)  = %d   reduced model thinks full, truly slack", length(spurious))

p()
p("3) THERMAL VIOLATIONS  (reduced dispatch replayed on the original network)")
p("-" ^ 78)
f("   feasible on original     = %s", replay.feasible ? "yes" : "NO")
f("   violated lines           = %d   (%d internal, %d external, %d protected)",
  length(viol), count(l -> internal[l], viol), count(l -> !internal[l], viol),
  count(l -> protected[l], viol))
if !isempty(viol)
    w = viol[argmax(over[viol])]
    f("   worst violation          = %.3f MW  =  %.5f p.u.  =  %.3f%% of rating",
      over[w] * mva, over[w], 100 * over_rel[w])
    f("                              on line %d (%d -> %d), rating %.1f MW, true flow %.1f MW",
      w, bus[base.Efrom[w]], bus[base.Eto[w]], rate[w] * mva, abs(replay.flow[w]) * mva)
    wr = viol[argmax(over_rel[viol])]
    wr == w || f("   worst in %% of rating     = %.3f%% on line %d (%.3f MW)",
                 100 * over_rel[wr], wr, over[wr] * mva)
    f("   total overload           = %.3f MW  =  %.5f p.u.", sum(over) * mva, sum(over))
    f("   worst loading            = %.2f%% of rating", 100 * maximum(u_rep))
    f("   uniform rating uplift needed for feasibility = %.3f%%", 100 * (maximum(u_rep) - 1))
    p()
    p("   line   from->to   kind      prot  rating MW  reduced-model MW  true MW  over MW   over %  over p.u.")
    for l in sort(viol; by = l -> -over[l])
        f("   %-5d %4d->%-4d  %-8s  %-4s  %9.1f  %16.2f  %7.2f  %7.3f  %6.3f  %9.5f",
          l, bus[base.Efrom[l]], bus[base.Eto[l]], internal[l] ? "internal" : "external",
          protected[l] ? "yes" : "no", rate[l] * mva, red.flow[l] * mva,
          replay.flow[l] * mva, over[l] * mva, 100 * over_rel[l], over[l])
    end
end
f("   sub-tolerance residuals  = %d lines over rating by less than %.0e (solver dust)",
  count(l -> over[l] > 1e-9 && over_rel[l] <= VIOL_TOL, 1:L), VIOL_TOL)

p()
p("4) FLOW ACCURACY ON SURVIVING LINES  (reduced model vs true replay flow)")
p("-" ^ 78)
ferr = [abs(replay.flow[l] - red.flow[l]) for l in external]
ferr_n = [abs(replay.flow[l] - red.flow[l]) / rate[l] for l in external]
f("   max flow error           = %.2f MW  (%.2f%% of that line's rating)",
  maximum(ferr) * mva, 100 * ferr_n[argmax(ferr)])
f("   max normalised error     = %.2f%% of rating", 100 * maximum(ferr_n))
f("   mean normalised error    = %.3f%% of rating", 100 * mean(ferr_n))
f("   lines with error > 1%% rating = %d,  > 5%% = %d", count(>(0.01), ferr_n), count(>(0.05), ferr_n))
f("   internal lines, true loading: max %.1f%%, %d above 80%%, %d at limit",
  100 * maximum(u_rep[internal]; init=0.0), count(l -> internal[l] && u_rep[l] >= 0.8, 1:L),
  count(l -> internal[l] && u_rep[l] >= 1 - BIND_TOL, 1:L))

p()
p("5) COST")
p("-" ^ 78)
f("   full-network optimum     = %.2f \$/h", full.objective)
f("   reduced-network optimum  = %.2f \$/h   (%+.4f%%)", red.objective,
  100 * (red.objective - full.objective) / full.objective)
p("   ", red.objective < full.objective ?
      "reduced is CHEAPER than the truth: optimistic -- it has dropped a real constraint" :
      "reduced is at least as expensive as the truth: pessimistic or exact")
if repair.feasible
    f("   repair: minimum redispatch to make the reduced dispatch feasible = %.2f MW", repair.redispatch * mva)
    f("           (%.3f%% of total load, %.5f p.u.)", 100 * repair.redispatch / sum(load), repair.redispatch)
    f("           repaired dispatch cost = %.2f \$/h  (%+.4f%% vs full optimum)",
      repair.repaired_cost, 100 * (repair.repaired_cost - full.objective) / full.objective)
else
    f("   repair: INFEASIBLE (%s) -- the full network cannot absorb this dispatch", repair.status)
end

p()
p("6) GENERATOR DISPATCH  (reduced vs full optimum)")
p("-" ^ 78)
f("   generators moved > 0.1 MW = %d of %d", length(moved), G)
f("   total |dispatch change|   = %.2f MW  (%.3f%% of load)", sum(abs, dpg), 100 * sum(abs, dpg) / (sum(load) * mva))
if !isempty(moved)
    g = moved[argmax(abs.(dpg[moved]))]
    f("   largest single move       = %+.2f MW at generator %d (bus %d)", dpg[g], g, bus[base.gen_bus[g]])
end

p()
p("7) PRICES  (LMP, \$/MWh)")
p("-" ^ 78)
f("   max LMP error            = %.3f   at bus %d", maximum(lmp_err), bus[argmax(lmp_err)])
f("   load-weighted LMP error  = %.3f", lmp_lw)
f("   buses with error > 0.1   = %d,  > 1 = %d  (of %d)", count(>(0.1), lmp_err), count(>(1), lmp_err), N)
f("   distinct price levels    = %d full  ->  %d reduced", levels(full.lmp), levels(red.lmp))
f("   LMP range full           = %.2f .. %.2f", minimum(full.lmp), maximum(full.lmp))
f("   LMP range reduced        = %.2f .. %.2f", minimum(red.lmp), maximum(red.lmp))

report = String(take!(io))
print(report)

# ---- files -----------------------------------------------------------------
out = joinpath(SETTING, "infeasibility")
mkpath(out)
write(joinpath(out, "summary.txt"), report)

line_header = "line,from_bus,to_bus,internal,protected,rating_mw,full_flow_mw,full_loading," *
    "reduced_flow_mw,reduced_loading,true_flow_mw,true_loading,overload_mw,overload_pct," *
    "overload_pu,flow_error_mw,at_limit_full,at_limit_reduced,violated"
line_row(l) = join((l, bus[base.Efrom[l]], bus[base.Eto[l]], Int(internal[l]), Int(protected[l]),
    rate[l] * mva, full.flow[l] * mva, u_full[l], red.flow[l] * mva, u_red[l],
    replay.flow[l] * mva, u_rep[l], over[l] * mva, 100 * over_rel[l], over[l],
    internal[l] ? "" : abs(replay.flow[l] - red.flow[l]) * mva,
    Int(l in bind_full), Int(l in bind_red), Int(l in viol)), ",")
open(joinpath(out, "lines.csv"), "w") do fh
    println(fh, line_header)
    foreach(l -> println(fh, line_row(l)), 1:L)
end
open(joinpath(out, "congested_lines.csv"), "w") do fh
    println(fh, line_header)
    foreach(l -> println(fh, line_row(l)), sort(union(bind_full, bind_red, viol)))
end
open(joinpath(out, "generators.csv"), "w") do fh
    println(fh, "generator,bus,pmin_mw,pmax_mw,full_mw,reduced_mw,change_mw,repair_mw")
    for g in 1:G
        println(fh, join((g, bus[base.gen_bus[g]], base.pmin[g] * mva, base.pmax[g] * mva,
            full.pG[g] * mva, red.pG[g] * mva, dpg[g],
            repair.feasible ? repair.pG[g] * mva : ""), ","))
    end
end
open(joinpath(out, "lmp.csv"), "w") do fh
    println(fh, "bus,representative_bus,load_mw,lmp_full,lmp_reduced,abs_error")
    for i in 1:N
        println(fh, join((bus[i], bus[rep[i]], load[i] * mva, full.lmp[i], red.lmp[i], lmp_err[i]), ","))
    end
end
println("\nwrote ", out)

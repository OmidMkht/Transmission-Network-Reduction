# --------------------------------------------------------------------------- #
# Reporting, file output, and the relaxation sweep driver.
#
# Everything here is PRESENTATION and ORCHESTRATION: each function takes results
# that have already been computed and prints them, writes them, or loops a solve
# over a list of settings. No model is formulated in this file.
#
# It exists so the runners can stay what a runner should be -- a configuration
# block, one call to run, and a set of toggles for which validations to show.
#
#   include("transmission_reduction_reporting.jl")
#   TxReport.section("1) Reduction summary")
#
# The module is deliberately agnostic about WHICH reduction module it is driving:
# the solving module is handed in as an argument (`TR`), so the same sweep driver
# serves both the single-scenario and the multi-scenario runner.
# --------------------------------------------------------------------------- #

module TxReport

using Printf
using DelimitedFiles
using Dates
using Statistics
using LinearAlgebra
using Graphs

# --------------------------------------------------------------------------- #
# Small formatting helpers
# --------------------------------------------------------------------------- #

section(title; io::IO=stdout) =
    (println(io); println(io, title); println(io, repeat("-", length(title))))
section(io::IO, title) = section(title; io=io)

"Round for display without turning a meaningful tiny number into 0.0."
function sig(x, digits::Int=4)
    (x isa Real && isfinite(x)) || return x
    x == 0 && return 0.0
    return round(x, sigdigits=digits)
end

pct(x, digits::Int=3) = (x isa Real && isfinite(x)) ? string(round(x, digits=digits), "%") : "n/a"

"Label a (mode, delta) relaxation setting for tables, filenames and titles."
function relaxation_label(mode::Symbol, delta::Real)
    mode === :none && return "none"
    return string(mode, "_", replace(string(round(100 * delta, digits=3)), "." => "p"), "pct")
end

relaxation_pretty(mode::Symbol, delta::Real) =
    mode === :none ? "none (exact pin)" :
    string(mode, "  δ=", round(100 * delta, digits=3), "% of rating")

"Fixed-width form for table rows -- `relaxation_pretty` overflows the column."
relaxation_short(mode::Symbol, delta::Real) =
    mode === :none ? "none (exact pin)" :
    string(mode, " ", round(100 * delta, digits=3), "%")

"""
Open a file for viewing. Prefers the VS Code CLI (`code -r`, opens as a tab
in the running IDE window) when it's on PATH; falls back to the OS default
application otherwise. Never throws.
"""
function open_externally(path::AbstractString)
    try
        code_cli = Sys.which("code")
        if !isnothing(code_cli)
            run(ignorestatus(`$code_cli -r $path`))
        elseif Sys.iswindows()
            # PowerShell's -Command joins every remaining argument and re-parses
            # the result, so a path containing spaces must arrive as ONE
            # already-quoted command string.
            # Invoke-Item, NOT Start-Process: Windows PowerShell 5.1's
            # Start-Process has no -LiteralPath parameter and fails with
            # "A parameter cannot be found that matches parameter name
            # 'LiteralPath'". Invoke-Item does have it, and -LiteralPath is what
            # stops [ ] in a path being treated as a wildcard.
            run(ignorestatus(`powershell -NoProfile -Command "Invoke-Item -LiteralPath '$path'"`))
        elseif Sys.isapple()
            run(ignorestatus(`open $path`); wait=false)
        elseif Sys.islinux()
            run(ignorestatus(`xdg-open $path`); wait=false)
        end
    catch err
        @warn "Could not open the file automatically; use the saved copy" path exception=err
    end
    return nothing
end

# --------------------------------------------------------------------------- #
# Plotting is loaded lazily: CairoMakie is by far the slowest include in the
# project, so a run with plots switched off must never pay for it.
# --------------------------------------------------------------------------- #
const _PLOTS_LOADED = Ref(false)

function ensure_plots(plots_file::AbstractString)
    _PLOTS_LOADED[] && return nothing
    isfile(plots_file) || error("Plot helper not found: $plots_file")
    Base.include(@__MODULE__, plots_file)
    _PLOTS_LOADED[] = true
    return nothing
end

"""
Render the before/after network figure. `plots_file` is included on first use.
Returns the path, or `nothing` when plotting is disabled.
"""
function plot_network(plots_file, case, A; path, open_after::Bool=false, kwargs...)
    ensure_plots(plots_file)
    f = getfield(@__MODULE__, :plot_reduction)
    Base.invokelatest(f, case, A; path=path, kwargs...)
    open_after && open_externally(path)
    return path
end

# --------------------------------------------------------------------------- #
# Multi-scenario report sections
# --------------------------------------------------------------------------- #

function report_case(c, selection, cfg, month_indices; io::IO=stdout)
    section(io, "Case")
    println(io, "  month                     = ", monthname(cfg.month), " ", cfg.year,
            "   (", length(month_indices), " hourly scenarios)")
    println(io, "  operating-point/load scale= ", cfg.operating_point_scale,
            "   line-limit scale = ", cfg.line_limit_scale)
    println(io, "  congestion definition     = |flow| / rating >= ", cfg.near_limit_threshold)
    worst_util, worst_at = findmax(selection.utilization)
    println(io, "  max monthly utilization   = ", sig(worst_util),
            "  (scenario ID ", c.scenario_ids[month_indices[worst_at[2]]],
            ", line ", worst_at[1], ")")
    println(io, "  congested-line union      = ", length(selection.congested_lines), " lines")
    println(io, "  seed ranking              = ", get(selection, :ranking, :none),
            "  (biggest first; :flow = sum |line flow|, :load = total demand)")
    println(io, "  min set-cover / seed size = ", selection.minimum_cover_size,
            " / ", selection.selected_count,
            "   (cover is a hard requirement, rank fills the rest)")
    println(io, "  seed scenario IDs         = ", join(selection.scenario_ids, ", "))
    isempty(selection.congested_lines) && println(io, 
        "  WARNING: no line meets the congestion threshold this month; " *
        "congestion coverage is vacuous and the selector simply took the highest-stress hours.")
    return nothing
end

function report_reduction(c, r, A; io::IO=stdout)
    section(io, "1) Reduction summary")
    sizes = vec(sum(A, dims=2))
    sizes = sizes[sizes .> 0]
    println(io, "  status / solve time = ", r.status, "  /  ", sig(r.solve_time), " s")
    println(io, "  reduction           = ",
            pct(100 * (c.base.N - r.n_retained) / c.base.N, 1),
            "  (", r.n_retained, " / ", c.base.N, " buses retained)")
    println(io, "  internal / external = ", r.n_internal_lines, " / ", r.n_external_lines, " lines")
    println(io, "  largest cluster     = ", maximum(sizes), " buses")
    if get(r, :merge_exact_mode, :off) !== :off
        println(io, "  exactly-mergeable   = ", r.n_merge_lines, " lines ",
                r.merge_exact_mode === :fix ? "FIXED internal" : "warm-started",
                "   (", r.n_merge_bridges, " bridges, ",
                r.n_merge_leaf_blocks, " leaf blocks)")
        println(io, "                        fixed-injection flow preservation; ",
                r.merge_exact_mode === :fix ?
                "fixed contractions restrict the search" :
                "warm starts leave the search unrestricted")
    end
    if get(r, :n_lmp_rows, 0) > 0
        println(io, "  LMP separation      = ", r.n_lmp_rows, " shortest-path rows at ",
                r.lmp_threshold, " \$/MWh")
        println(io, "                        ", r.n_lmp_separated, " / ",
                length(r.lmp_pairs), " constrained pairs came out separated",
                r.n_lmp_still_merged == 0 ? "" :
                ";  " * string(r.n_lmp_still_merged) *
                " merged anyway via another path, e.g. " *
                string(r.lmp_merged_examples[1:min(5, end)]))
    end
    return nothing
end

function report_model_check(TR, chk; io::IO=stdout)
    section(io, "2) Feasibility of the returned MILP point")
    TR.print_model_feasibility_check_multiscenario(chk; io=io)
    return chk
end

function report_benchmark(c, bench, month_indices; io::IO=stdout)
    section(io, "3) Benchmark over all $(length(month_indices)) monthly scenarios")
    mva = c.base.baseMVA
    println(io, "  scenarios inside reduction windows = ", bench.n_window_feasible,
            " / ", bench.n_scenarios,
            "  (", pct(100 * bench.window_feasible_fraction, 2), ")")
    println(io, "  worst window violation             = ", sig(mva * bench.worst_window_violation), " MW")
    println(io, "  p95 / mean window violation        = ", sig(mva * bench.p95_window_violation),
            " / ", sig(mva * bench.mean_window_violation), " MW")
    println(io, "  worst external flow error          = ", sig(bench.max_external_normalized_error),
            "  (normalized by rating)")
    println(io, "  worst all-line flow difference     = ", sig(bench.max_all_line_normalized_error),
            "  (includes lines made internal)")
    println(io, "  max absolute flow difference       = ", sig(mva * bench.max_absolute_flow_error), " MW")
    println(io, "  max reduced-network overload       = ", sig(mva * bench.max_rating_overload), " MW")
    println(io)
    println(io, "  external (surviving) lines overloaded in the reduced network:")
    plural(n, word) = string(n, " ", word, n == 1 ? "" : "s")
    println(io, "    vs plain RATING    = ", plural(bench.n_lines_over_rating, "line"), " in ",
            bench.n_scenarios_over_rating, " / ", bench.n_scenarios, " scenarios",
            "   (", plural(bench.n_pairs_over_rating, "line-hour pair"), ", worst ",
            sig(mva * bench.max_overload_vs_rating), " MW)")
    println(io, "    vs ADJUSTED cap    = ", plural(bench.n_lines_over_adjusted, "line"), " in ",
            bench.n_scenarios_over_adjusted, " / ", bench.n_scenarios, " scenarios",
            "   (", plural(bench.n_pairs_over_adjusted, "line-hour pair"), ", worst ",
            sig(mva * bench.max_overload_vs_adjusted), " MW)")
    println(io, "    ADJUSTED cap = rating + delta on a relaxed congested line. A line past")
    println(io, "    the rating but inside its adjusted cap is doing what the model allowed;")
    println(io, "    one past the ADJUSTED cap escaped even that -- the clustering is at fault.")
    if !isempty(bench.lines_over_rating)
        worst_lines = sort(bench.lines_over_rating;
                           by=l -> -bench.overload_hours_by_line[l])
        println(io, "    most frequently overloaded: ",
                join(["L$l ($(bench.overload_hours_by_line[l]) h)"
                      for l in first(worst_lines, 8)], ", "))
    end
    println(io, "  congested lines made internal      = ", length(bench.congested_internal),
            isempty(bench.congested_internal) ? "" : "   <-- coverage/protection failure")
    println(io, "  worst scenario ID / line           = ", bench.screen.worst_scenario_id,
            " / ", bench.screen.worst_line)
    return nothing
end

function report_dcopf(c, val; io::IO=stdout)
    section(io, "4) Full-vs-reduced DC-OPF validation")
    println(io, "  scope: returned dispatches at listed demands; alternate optima and unseen demands are not certified")
    println(io, "  scenario IDs                       = ", join(val.scenario_ids, ", "))
    mva = c.base.baseMVA
    H = length(val.scenario_indices)
    println(io, "  Pmin relaxed to zero               = ", val.relax_pmin)
    println(io, "  both DC-OPFs solved                = ", val.n_dcopf_feasible, " / ", H)
    println(io, "  line limits used here              = TRUE ratings (frate), NOT the relaxed cap")
    println(io, "    the congestion relaxation delta enters only the MILP window and the")
    println(io, "    window benchmark. It reaches this section solely through the clustering")
    println(io, "    it produced, so these numbers are directly comparable across settings.")
    println(io)
    println(io, "  4a) strict pass/fail")
    println(io, "      reduced dispatch feasible on original network = ",
            val.n_dispatch_feasible, " / ", H,
            "   (tolerance ", val.relative_tolerance, " of rating)")
    println(io, "      objective gap <= ", val.objective_tolerance_pct, "%                      = ",
            val.n_objective_within_tolerance, " / ", H)
    println(io, "      max LMP difference <= ", val.lmp_tolerance, " \$/MWh          = ",
            val.n_lmp_within_tolerance, " / ", H)
    println(io, "      worst absolute objective change              = ",
            pct(val.worst_abs_objective_change_pct))
    println(io, "      worst LMP difference                         = ",
            sig(val.worst_lmp_error), " \$/MWh")
    println(io)
    println(io, "  4b) graded severity  (stays comparable when the strict test is uniformly \"fail\")")
    println(io, "      worst / mean full-network utilization        = ",
            sig(val.worst_max_utilization), " / ", sig(val.mean_max_utilization),
            "   (1.0 = exactly at rating)")
    println(io, "      worst original-network overload              = ",
            sig(mva * val.worst_overload), " MW")
    println(io, "      worst total overload over all lines          = ",
            sig(mva * val.worst_total_overload), " MW")
    println(io, "      most overloaded lines in a single hour       = ", val.max_violating_lines)
    println(io)
    println(io, "  4c) spuriously binding lines  (the PESSIMISTIC error, opposite of an overload)")
    println(io, "      a line at its rating in the REDUCED dispatch but slack in the true")
    println(io, "      full-network optimum. Nothing is violated -- the reduction simply")
    println(io, "      believes a corridor is full when it is not, pushing the OPF off the")
    println(io, "      cheap dispatch. This is where a cost gap with NO feasibility failure")
    println(io, "      comes from.")
    println(io, "      hours with >= 1 spurious binding line   = ",
            val.n_hours_with_spurious_binding, " / ", H)
    println(io, "      distinct lines ever spuriously binding  = ", val.n_spurious_lines)
    println(io, "      most in a single hour                   = ", val.max_spurious_binding,
            "   (binding tolerance ", val.binding_tolerance, " of rating)")
    if val.n_spurious_lines > 0
        worst = sort(findall(>(0), val.line_spurious_hours);
                     by=l -> -val.line_spurious_hours[l])
        println(io, "      most frequent: ",
                join(["L$l ($(val.line_spurious_hours[l]) h)" for l in first(worst, 8)], ", "))
        println(io, "      => these corridors are over-constrained by the clustering; a better")
        println(io, "         clustering around them would recover cost, not feasibility")
    end
    if val.measure_repair
        println(io)
        println(io, "  4d) repair economics  (cheapest generation move that makes the reduced")
        println(io, "      dispatch work on the full network)")
        println(io, "      hours repairable                             = ",
                val.n_repair_feasible, " / ", H,
                val.n_repair_feasible == H ? "" :
                "   <-- a failed repair solve does not by itself prove full-network infeasibility")
        println(io, "      worst / mean redispatch                      = ",
                sig(mva * val.worst_repair_redispatch), " / ",
                sig(mva * val.mean_repair_redispatch), " MW")
        println(io, "      worst / mean repaired-cost gap vs optimum    = ",
                pct(val.worst_repair_cost_pct), " / ", pct(val.mean_repair_cost_pct))
    end
    return nothing
end

"Per-iteration trace of the scenario-generation loop."
# --------------------------------------------------------------------------- #
# Console digest -- the only thing a run prints by default.
#
# Three blocks, because a run is really asked three questions: how much did it
# reduce, did the reduced flows stay inside their own windows, and does the
# reduced dispatch survive on the true full network. Everything else -- the
# residual check, the solve-time benchmark, the graded severity, the overload
# anatomy, the sweep table -- is still produced in full and goes to report.txt.
# --------------------------------------------------------------------------- #
"Persist the same certificate outcome and scope in console and saved reports."
function report_derating_scope(c, derate; io::IO=stdout)
    isnothing(derate) && return nothing
    section(io, "Derating outcome and certification scope")
    println(io, "  empirical outcome = ", derate.emp.reason)
    println(io, "  empirical scenario IDs = ", join(get(derate, :empirical_ids, []), ";"))
    println(io, "  requested certificate scenario IDs = ", join(get(derate, :certificate_ids, []), ";"))
    cert = derate.cert
    println(io, "  certification = ", isnothing(cert) ? "not run" :
        cert.certified ? "passed within stated scope and tolerances" : "not established: $(cert.reason)")
    if !isnothing(cert)
        println(io, "  dispatch scope = ", get(cert, :dispatch_scope, "unrecorded"))
        println(io, "  cost band percent above reduced optimum = ", get(cert, :cost_cap_pct, "unrecorded"))
        println(io, "  effective cost caps in supplied-scenario order = ", get(cert, :cost_caps, "unrecorded"))
        println(io, "  monitored original lines = ", get(cert, :n_monitored_lines, "unrecorded"))
        println(io, "  absolute flow tolerance (pu) = ", get(cert, :tolerance, "unrecorded"))
        println(io, "  generation minimum relaxed = ", get(cert, :relax_pmin, "unrecorded"))
        println(io, "  This does not cover other demands or prove that unsuccessful derating is impossible.")
    end
    println(io, "  derated original line IDs = ", join(derate.derated, ";"))
    return nothing
end

"Console summary, with optional empirical validation under derated ratings."
function report_digest(c, art, month_indices; io::IO=stdout,
                       derate=nothing, val_derated=nothing)
    println(io)
    println(io, repeat("=", 78))
    println(io, "RESULTS:  ", relaxation_pretty(art.mode, art.delta))
    println(io, repeat("=", 78))

    report_reduction(c, art.r, art.A; io=io)

    bench = art.bench
    section(io, "3) Benchmark over all $(length(month_indices)) monthly scenarios")
    println(io, "  scenarios inside reduction windows = ", bench.n_window_feasible,
            " / ", bench.n_scenarios,
            "  (", pct(100 * bench.window_feasible_fraction, 2), ")")
    println(io, "  worst external flow error          = ",
            sig(bench.max_external_normalized_error), "  (normalized by rating)")

    if !isnothing(derate)
        section(io, "Derating (post-reduction; the clustering is unchanged)")
        println(io, "  hours checked       = ", derate.n_emp, " empirical / ",
                derate.n_cert, " requested for certification",
                "   (one rating per line, set by the binding hour)")
        println(io, "  empirical descent   = ", derate.emp.reason)
        println(io, "  certified           = ", isnothing(derate.cert) ? "not run" :
                derate.cert.certified ?
                    "YES -- within the reported cost band and checked scenarios" :
                    "NO (" * string(derate.cert.reason) * ")")
        println(io, "  lines derated       = ", length(derate.derated), " / ", c.base.Ln,
                isempty(derate.derated) ? "" :
                "   worst " * string(round(derate.worst_pct, digits=1)) * "% off its rating")
        if !isnothing(derate.cert)
            println(io, "  certificate scope   = ", get(derate.cert, :n_scenarios, derate.n_cert),
                    " checked demands; cost cap (%): ", get(derate.cert, :cost_cap_pct, "unrecorded"),
                    "; Pmin relaxed: ", get(derate.cert, :relax_pmin, "unrecorded"))
            println(io, "  certificate tolerance (pu) = ", get(derate.cert, :tolerance, "unrecorded"))
            derate.cert.certified || println(io,
                "  failure means this search did not certify the ratings; other deratings may still work")
        end
    end

    val = art.val
    isnothing(val) && return nothing
    H = length(val.scenario_indices)

    "One 4a block for a given validation result."
    function pass_fail(v, title)
        section(io, title)
        println(io, "  scope: one returned dispatch per checked demand; this is empirical validation")
        println(io, "  reduced dispatch feasible on original network = ",
                v.n_dispatch_feasible, " / ", H,
                "   (tolerance ", v.relative_tolerance, " of rating)")
        println(io, "  objective gap <= ", v.objective_tolerance_pct, "%                      = ",
                v.n_objective_within_tolerance, " / ", H)
        println(io, "  max LMP difference <= ", v.lmp_tolerance, " \$/MWh          = ",
                v.n_lmp_within_tolerance, " / ", H)
        println(io, "  worst absolute objective change              = ",
                pct(v.worst_abs_objective_change_pct))
        println(io, "  worst LMP difference                         = ",
                sig(v.worst_lmp_error), " \$/MWh")
    end

    if isnothing(val_derated)
        pass_fail(val, "4a) Full-vs-reduced DC-OPF: strict pass/fail")
    else
        pass_fail(val, "4a) Full-vs-reduced DC-OPF: strict pass/fail  --  WITHOUT derating")
        pass_fail(val_derated, "4a) Full-vs-reduced DC-OPF: strict pass/fail  --  WITH derating")
    end
    return nothing
end

function report_scenario_generation(c, gen; io::IO=stdout)
    section(io, "0) Scenario generation")
    mva = c.base.baseMVA
    println(io, "  seed scenarios   = ", length(gen.seed),
            "   final active = ", length(gen.active),
            "   iterations = ", gen.iterations)
    println(io, "  converged        = ", gen.converged,
            gen.converged ?
            "   (no monthly hour violates its window)" :
            "   <-- some monthly hours still violate")
    n_structural_fixes = count(h -> h.structural, gen.history)
    n_structural_fixes > 0 && println(io, "  shorted-protected-line fixes applied = ", n_structural_fixes,
            "  (a protected line's endpoints got merged via other internal lines;",
            " see the trace below)")
    if !isempty(gen.structural_failure)
        println(io)
        println(io, "  UNRESOLVED -- protected line(s) ", join(gen.structural_failure, ", "),
                " are still shorted by a path of internal lines, and every hour that")
        println(io, "  stresses them is already active. That combination should be provably")
        println(io, "  infeasible under an exact solve, so this points to solver tolerance")
        println(io, "  (IntFeasTol/FeasibilityTol) rather than a genuine model limit.")
    end
    header = @sprintf("%5s %8s %7s %10s %11s %11s %9s %s",
                      "iter", "active", "buses", "int.lines",
                      "worst above", "worst below", "bad hrs", "added IDs")
    println(io, "  ", header)
    println(io, "  ", repeat("-", length(header)))
    for h in gen.history
        finite(x) = isfinite(x) ? @sprintf("%.4f", mva * x) : "Inf"
        tag = h.structural ? "  (fixing a shorted protected line)" :
              h.both_sides  ? "  (both sides)" : ""
        println(io, "  ", @sprintf("%5d %8d %7d %10d %11s %11s %9d %s",
            h.iteration, h.n_active, h.n_retained, h.n_internal,
            finite(h.worst_above), finite(h.worst_below), h.n_violating_hours,
            isempty(h.added_ids) ?
                (gen.converged ? "-- converged --" :
                 isempty(gen.structural_failure) ? "-- stopped: max iterations --" :
                                                   "-- stopped: unresolved short --") :
                join(h.added_ids, ", ") * tag))
    end
    println(io)
    println(io, "  worst above/below are in MW, on opposite sides of the flow window.")
    println(io, "  Two hours are added when BOTH sides are breached: the window is two")
    println(io, "  half-spaces, and a support point for one side constrains nothing")
    println(io, "  about the other. A row tagged \"fixing a shorted protected line\" instead")
    println(io, "  targets the hour that forces a merged-but-protected line apart.")
    return nothing
end

"""
Anatomy of the hours where the reduced dispatch was NOT feasible on the full
network. "How many failed" says nothing about how to fix it; these four cuts do.

  INTERNAL vs EXTERNAL overloads. An internal line was merged away, so the
  reduced network never modelled it and could not have respected it -- the fix is
  to stop collapsing it (protect it, or tighten eps). An external line survived,
  so the reduced network did model it and got the flow wrong -- the fix is a
  better clustering around it, not more protection.

  IN-TRAINING vs OUT-OF-TRAINING. If failures land only on hours the MILP never
  saw, the scenario selection is too narrow -- widen the cover. If hours the MILP
  was trained on also fail, the model itself is too loose and no amount of extra
  scenarios will help.

  PER-LINE FREQUENCY. A handful of lines failing in most hours is a targeted fix
  (protect those corridors). Many lines each failing rarely is a systemic
  accuracy problem.

  LOAD LEVEL. Failures concentrated in high-load hours mean the reduction is fine
  off-peak and the binding cases are simply the stressed ones.
"""
function report_dcopf_failures(c, val, art; top_n::Int=10, io::IO=stdout)
    section(io, "5) Anatomy of the DC-OPF failures")
    mva = c.base.baseMVA
    evaluated = val.original_optimal .& val.reduced_optimal
    failed = evaluated .& .!val.dispatch_feasible
    nfail = count(failed)
    if nfail == 0
        println(io, "  every evaluated hour's reduced dispatch was feasible on the full network")
        return nothing
    end
    H = length(val.scenario_indices)
    println(io, "  infeasible hours = ", nfail, " / ", count(evaluated), " evaluated (of $H)")

    # (b) training-set membership
    fail_train = count(failed .& val.in_training)
    n_train = count(val.in_training)
    println(io)
    println(io, "  a) MILP training coverage")
    println(io, "     failing hours that were IN the MILP training set  = ", fail_train,
            " / ", n_train, " training hours present")
    println(io, "     failing hours the MILP never saw                  = ", nfail - fail_train)
    println(io, "     => ", fail_train == 0 ?
            "failures are ALL out-of-sample: widen the scenario cover" :
            "hours the MILP was trained on also fail: the model itself is too loose, " *
            "more scenarios will not fix it")

    # (a) internal vs external
    tot_int = sum(val.n_violating_internal[failed])
    tot_ext = sum(val.n_violating_external[failed])
    println(io)
    println(io, "  b) which KIND of line is overloaded  (line-hour pairs over all failing hours)")
    println(io, "     INTERNAL (merged away)  = ", tot_int)
    println(io, "     EXTERNAL (survived)     = ", tot_ext)
    println(io, "     => ", tot_int > tot_ext ?
            "mostly collapsed lines: the clustering is merging corridors that matter -- " *
            "protect them or lower the congestion threshold" :
            "mostly surviving lines: the reduced network models them but gets the flow " *
            "wrong -- an accuracy problem, not a protection problem")
    worst_h = argmax(ifelse.(failed, val.n_violating_lines, -1))
    println(io, "     worst single hour: ID ", val.scenario_ids[worst_h], " with ",
            val.n_violating_lines[worst_h], " overloaded lines (",
            val.n_violating_internal[worst_h], " internal / ",
            val.n_violating_external[worst_h], " external), utilization ",
            sig(val.max_utilization[worst_h]))

    # (c) per-line frequency
    println(io)
    println(io, "  c) most frequently overloaded lines")
    offenders = sort(findall(>(0), val.line_violation_hours);
                     by=l -> -val.line_violation_hours[l])
    protected_set = Set(findall(art.r.protected))
    println(io, "     ", @sprintf("%-6s %-9s %-10s %8s %12s", "line", "kind", "protected",
                              "hours", "worst MW"))
    for l in first(offenders, top_n)
        println(io, "     ", @sprintf("%-6s %-9s %-10s %8d %12s",
            "L$l",
            val.internal_line[l] ? "internal" : "external",
            l in protected_set ? "yes" : "no",
            val.line_violation_hours[l],
            @sprintf("%.3f", mva * val.line_worst_overload[l])))
    end

    # (d) load level
    println(io)
    println(io, "  d) are failures concentrated at peak load?")
    lf = val.total_load[failed]
    lp = val.total_load[evaluated .& val.dispatch_feasible]
    println(io, "     mean total load, failing hours = ",
            sig(mva * (isempty(lf) ? NaN : mean(lf))), " MW")
    println(io, "     mean total load, passing hours = ",
            sig(mva * (isempty(lp) ? NaN : mean(lp))), " MW")
    if !isempty(lf) && !isempty(lp)
        println(io, "     => ", mean(lf) > mean(lp) ?
                "failures skew toward HIGH load: the reduction holds off-peak and breaks under stress" :
                "failures are NOT load-driven: the cause is structural, not stress")
    end
    return nothing
end

# --------------------------------------------------------------------------- #
# Relaxation sweep: one solve per (mode, delta) setting, then a table.
#
# Every setting gets its own output subdirectory, so the per-setting CSVs and
# figures never overwrite one another and can be diffed directly.
# --------------------------------------------------------------------------- #

"""
    sweep_multiscenario(TR, c, epsL, selected, month_indices, cfg) -> (rows, artifacts)

Solve the reduction once per entry of `cfg.relaxation_sweep`, benchmarking and
(optionally) DC-OPF-validating each. Returns the comparison rows and the full
per-setting artifacts, so the caller can write detailed outputs for any of them.
"""
function sweep_multiscenario(TR, c, epsL, selected, month_indices, cfg)
    rows = NamedTuple[]
    artifacts = NamedTuple[]
    mva = c.base.baseMVA
    dc = cfg.dcopf

    for (k, (mode, delta)) in enumerate(cfg.relaxation_sweep)
        println()
        println(repeat("=", 78))
        println("SETTING $k/$(length(cfg.relaxation_sweep)):  ", relaxation_pretty(mode, delta))
        println(repeat("=", 78))

        # Optional Gurobi log, one file per relaxation setting so the settings do
        # not interleave. Gurobi APPENDS to an existing LogFile, which is what we
        # want WITHIN a run (scenario generation solves repeatedly and every
        # iteration lands in the same file, giving one continuous trace) but not
        # ACROSS runs, so a stale file from a previous run is cleared first.
        # `get` keeps this optional: a cfg without the field simply gets nothing.
        log_file = nothing
        if get(cfg, :gurobi_log, false)
            log_dir = joinpath(cfg.output_dir, relaxation_label(mode, delta))
            mkpath(log_dir)
            log_file = joinpath(log_dir, "gurobi.log")
            isfile(log_file) && rm(log_file)
            println("Gurobi log -> ", log_file)
        end

        # Scenario generation adapts the active set to THIS relaxation setting --
        # a looser window binds at different hours -- so `active` is per setting,
        # not shared. It is also what counts as the MILP training set downstream.
        gen = nothing
        if cfg.scenario_generation
            gen = TR.solve_reduction_with_scenario_generation(
                c, epsL, selected, month_indices;
                near_limit_threshold=cfg.near_limit_threshold,
                max_iterations=cfg.max_sc_generation_iterations,
                tolerance=cfg.screening_tolerance,
                time_limit=cfg.solve_time_limit,
                numeric_focus=cfg.numeric_focus,
                cycle_cut_lens=cfg.cycle_cut_lens,
                congestion_relaxation=delta,
                congestion_relaxation_mode=mode,
                internal_bound_scale=get(cfg, :internal_bound_scale, 3.0),
                internal_rating_bound=get(cfg, :internal_rating_bound, false),
                switch_form=get(cfg, :switch_form, :hull),
                force_internal=get(cfg, :force_internal, Int[]),
                merge_exact_blocks=get(cfg, :merge_exact_blocks, false),
                merge_exact_mode=get(cfg, :merge_exact_mode, :fix),
                merge_leaf_blocks=get(cfg, :merge_leaf_blocks, true),
                lmp_separation=get(cfg, :lmp_separation, false),
                lmp_threshold=get(cfg, :lmp_threshold, 5.0),
                lmp_relax_pmin=get(cfg, :lmp_relax_pmin, true),
                lmp_opf_time_limit=get(cfg, :opf_time_limit, nothing),
                int_feas_tol=get(cfg, :int_feas_tol, nothing),
                feasibility_tol=get(cfg, :feasibility_tol, nothing),
                optimality_tol=get(cfg, :optimality_tol, nothing),
                line_budget=get(cfg, :line_budget, nothing),
                hop_cap=get(cfg, :hop_cap, nothing),
                size_cap=get(cfg, :size_cap, nothing),
                log_file=log_file)
            r = gen.r
            active = gen.active
        else
            r = TR.solve_reduction_edge_multiscenario(
                c, epsL;
                scenario_indices=selected,
                protection_indices=month_indices,
                near_limit_threshold=cfg.near_limit_threshold,
                time_limit=cfg.solve_time_limit,
                numeric_focus=cfg.numeric_focus,
                cycle_cut_lens=cfg.cycle_cut_lens,
                congestion_relaxation=delta,
                congestion_relaxation_mode=mode,
                internal_bound_scale=get(cfg, :internal_bound_scale, 3.0),
                internal_rating_bound=get(cfg, :internal_rating_bound, false),
                switch_form=get(cfg, :switch_form, :hull),
                force_internal=get(cfg, :force_internal, Int[]),
                merge_exact_blocks=get(cfg, :merge_exact_blocks, false),
                merge_exact_mode=get(cfg, :merge_exact_mode, :fix),
                merge_leaf_blocks=get(cfg, :merge_leaf_blocks, true),
                lmp_separation=get(cfg, :lmp_separation, false),
                lmp_threshold=get(cfg, :lmp_threshold, 5.0),
                lmp_relax_pmin=get(cfg, :lmp_relax_pmin, true),
                lmp_opf_time_limit=get(cfg, :opf_time_limit, nothing),
                int_feas_tol=get(cfg, :int_feas_tol, nothing),
                feasibility_tol=get(cfg, :feasibility_tol, nothing),
                optimality_tol=get(cfg, :optimality_tol, nothing),
                line_budget=get(cfg, :line_budget, nothing),
                hop_cap=get(cfg, :hop_cap, nothing),
                size_cap=get(cfg, :size_cap, nothing),
                log_file=log_file)
            active = selected
        end
        dcopf_indices = dc.scope === :month ? month_indices :
                        dc.scope === :selected ? active : Int[]
        A = round.(Int, r.A)
        chk = TR.model_feasibility_check_multiscenario(c, r; tol=1e-6)
        # near_limit_threshold = nothing means "protect only EXACTLY binding
        # lines". That is a valid model setting, but the benchmark still needs a
        # concrete number to label a line congested, so fall back to the same
        # ~1.0 definition `binding_lines` uses for exact binding.
        congestion_threshold = isnothing(cfg.near_limit_threshold) ? 0.9999 :
                               Float64(cfg.near_limit_threshold)
        bench = TR.benchmark_reduction_scenarios(
            c, A, epsL, month_indices;
            congestion_threshold=congestion_threshold,
            near_limit_threshold=cfg.near_limit_threshold,
            protection_indices=month_indices,
            congestion_relaxation=delta,
            congestion_relaxation_mode=mode,
            tolerance=cfg.screening_tolerance)
        val = isempty(dcopf_indices) ? nothing : TR.validate_reduced_dcopf_scenarios(
            c, A, dcopf_indices;
            relax_pmin=dc.relax_pmin,
            time_limit=dc.time_limit,
            relative_tolerance=dc.relative_tolerance,
            objective_tolerance_pct=dc.objective_tolerance_pct,
            lmp_tolerance=dc.lmp_tolerance,
            measure_repair=dc.measure_repair,
            training_indices=active)
        # Legacy A-based cross-check. Only meaningful for ONE operating point
        # (it re-derives a single angle vector from A), so the single-scenario
        # runner enables it and the multi-scenario runner leaves it off.
        # THE COMPUTATIONAL PAYOFF. Everything above measures whether the reduced
        # network is ACCURATE; this measures whether it is CHEAPER, which is the
        # entire point of reducing it. Full vs natively-rebuilt reduced DC-OPF,
        # one builder for both, Threads=1, min of `repeats` after a warm-up.
        timing = get(cfg, :solve_time_benchmark, true) && !isempty(dcopf_indices) ?
            TR.benchmark_dcopf_solve_times(
                c, A, dcopf_indices;
                repeats=get(cfg, :benchmark_repeats, 5),
                relax_pmin=dc.relax_pmin,
                threads=get(cfg, :benchmark_threads, 1),
                time_limit=dc.time_limit,
                progress_every=get(cfg, :benchmark_progress_every, 200)) : nothing
        cons = get(cfg.show, :original_model, false) ?
            check_original_model_constraints(
                TR, c.base, A, r, epsL;
                near_limit_threshold=cfg.near_limit_threshold,
                congestion_relaxation=delta,
                congestion_relaxation_mode=mode) : nothing

        push!(rows, (
            mode=mode, delta=delta,
            label=relaxation_label(mode, delta),
            n_retained=r.n_retained,
            reduction_pct=100 * (c.base.N - r.n_retained) / c.base.N,
            n_internal=r.n_internal_lines,
            solve_time=r.solve_time,
            status=string(r.status),
            genuine=chk.genuine,
            window_feasible_pct=100 * bench.window_feasible_fraction,
            worst_window_mw=mva * bench.worst_window_violation,
            worst_external_error=bench.max_external_normalized_error,
            n_dcopf_feasible=isnothing(val) ? -1 : val.n_dispatch_feasible,
            n_dcopf_total=isnothing(val) ? 0 : length(val.scenario_indices),
            worst_utilization=isnothing(val) ? NaN : val.worst_max_utilization,
            worst_objective_pct=isnothing(val) ? NaN : val.worst_abs_objective_change_pct,
            worst_repair_mw=isnothing(val) ? NaN : mva * val.worst_repair_redispatch,
            mean_repair_cost_pct=isnothing(val) ? NaN : val.mean_repair_cost_pct,
            n_ext_lines_over_rating=bench.n_lines_over_rating,
            n_scenarios_over_rating=bench.n_scenarios_over_rating,
            n_ext_lines_over_adjusted=bench.n_lines_over_adjusted,
            n_scenarios_over_adjusted=bench.n_scenarios_over_adjusted,
            n_fail_in_training=isnothing(val) ? -1 :
                count((val.original_optimal .& val.reduced_optimal .&
                       .!val.dispatch_feasible) .& val.in_training),
            overload_internal_pairs=isnothing(val) ? -1 :
                sum(val.n_violating_internal),
            overload_external_pairs=isnothing(val) ? -1 :
                sum(val.n_violating_external),
            n_hours_spurious_binding=isnothing(val) ? -1 :
                val.n_hours_with_spurious_binding,
            n_spurious_lines=isnothing(val) ? -1 : val.n_spurious_lines,
            n_active=length(active),
            cg_iterations=isnothing(gen) ? 0 : gen.iterations,
            cg_converged=isnothing(gen) ? false : gen.converged,
            # --- computational payoff of the reduction ---
            bench_buses_full=isnothing(timing) ? -1 : timing.dims_full.n_bus,
            bench_buses_reduced=isnothing(timing) ? -1 : timing.dims_reduced.n_bus,
            bench_lines_full=isnothing(timing) ? -1 : timing.dims_full.n_line,
            bench_lines_reduced=isnothing(timing) ? -1 : timing.dims_reduced.n_line,
            speedup_work=isnothing(timing) ? NaN : timing.speedup_work,
            speedup_wall=isnothing(timing) ? NaN : timing.speedup_wall,
            speedup_iters=isnothing(timing) ? NaN : timing.speedup_iters,
            speedup_solve_time=isnothing(timing) || !timing.solve_time_resolved ?
                NaN : timing.speedup_solve_time,
            wall_ms_full=isnothing(timing) ? NaN :
                1000 * timing.total_wall_full / max(count(timing.both_optimal), 1),
            wall_ms_reduced=isnothing(timing) ? NaN :
                1000 * timing.total_wall_reduced / max(count(timing.both_optimal), 1),
        ))
        push!(artifacts, (r=r, A=A, chk=chk, bench=bench, val=val, gen=gen,
                          cons=cons, timing=timing, active=active,
                          mode=mode, delta=delta,
                          label=relaxation_label(mode, delta)))
    end
    return rows, artifacts
end

"Side-by-side comparison of every relaxation setting in the sweep."
function print_relaxation_table(rows; io::IO=stdout)
    section(io, "RELAXATION COMPARISON")
    header = @sprintf("%-22s %6s %5s %7s %8s %9s %10s %11s %11s %10s",
                      "setting", "scen", "iter", "buses", "reduct.", "int.lines",
                      "windowOK%", "worstUtil", "repair MW", "cost gap%")
    println(io, header)
    println(io, repeat("-", length(header)))
    for row in rows
        println(io, @sprintf("%-22s %6d %5s %7d %7.1f%% %9d %10.1f %11s %11s %10s",
            relaxation_short(row.mode, row.delta),
            row.n_active,
            row.cg_iterations == 0 ? "-" :
                string(row.cg_iterations, row.cg_converged ? "" : "*"),
            row.n_retained, row.reduction_pct, row.n_internal,
            row.window_feasible_pct,
            isfinite(row.worst_utilization) ? @sprintf("%.4f", row.worst_utilization) : "-",
            isfinite(row.worst_repair_mw) ? @sprintf("%.2f", row.worst_repair_mw) : "-",
            isfinite(row.mean_repair_cost_pct) ? @sprintf("%.4f", row.mean_repair_cost_pct) : "-"))
    end
    println(io)
    println(io, "  scen       scenarios the MILP was finally solved on")
    println(io, "  iter       scenario-generation iterations; * = hit max without converging")
    println(io, "  buses      retained after reduction (lower = more reduction)")
    println(io, "  windowOK%  monthly hours whose reduced flows stay inside their own windows")
    println(io, "  worstUtil  worst |flow| / rating on the FULL network under the reduced")
    println(io, "             dispatch; 1.0 = exactly at rating, excess is the headroom needed")
    println(io, "  repair MW  worst-hour generation move needed to make that dispatch feasible")
    println(io, "  cost gap%  mean cost of the repaired dispatch vs the true full-network optimum")
    base = first(rows)
    println(io)
    for row in rows[2:end]
        Δbuses = base.n_retained - row.n_retained
        println(io, "  ", relaxation_short(row.mode, row.delta), ":  ",
                Δbuses > 0 ? "$Δbuses fewer buses" :
                Δbuses == 0 ? "NO extra reduction" : "$(-Δbuses) MORE buses",
                " than the unrelaxed baseline",
                isfinite(row.worst_utilization) ?
                    "   (worst utilization $(sig(row.worst_utilization)))" : "")
    end
    return nothing
end

function write_relaxation_table(path, rows)
    open(path, "w") do io
        println(io, "mode,delta_fraction_of_rating,n_retained,reduction_pct,",
                "n_internal_lines,solve_time_s,status,genuine,",
                "window_feasible_pct,worst_window_violation_mw,worst_external_error,",
                "n_dcopf_dispatch_feasible,n_dcopf_scenarios,worst_full_network_utilization,",
                "worst_objective_change_pct,worst_repair_redispatch_mw,mean_repair_cost_pct,",
                "n_external_lines_over_rating,n_scenarios_over_rating,",
                "n_external_lines_over_adjusted,n_scenarios_over_adjusted,",
                "n_failing_hours_in_training,overload_internal_line_hours,",
                "overload_external_line_hours,",
                "n_hours_with_spuriously_binding_lines,n_spuriously_binding_lines,",
                "n_active_scenarios,generation_iterations,generation_converged")
        for r in rows
            println(io, join((r.mode, r.delta, r.n_retained, r.reduction_pct,
                              r.n_internal, r.solve_time, r.status, r.genuine,
                              r.window_feasible_pct, r.worst_window_mw,
                              r.worst_external_error, r.n_dcopf_feasible,
                              r.n_dcopf_total, r.worst_utilization,
                              r.worst_objective_pct, r.worst_repair_mw,
                              r.mean_repair_cost_pct,
                              r.n_ext_lines_over_rating, r.n_scenarios_over_rating,
                              r.n_ext_lines_over_adjusted, r.n_scenarios_over_adjusted,
                              r.n_fail_in_training, r.overload_internal_pairs,
                              r.overload_external_pairs,
                              r.n_hours_spurious_binding, r.n_spurious_lines,
                              r.n_active,
                              r.cg_iterations, r.cg_converged), ','))
        end
    end
    return path
end

# --------------------------------------------------------------------------- #
# Per-setting file output
# --------------------------------------------------------------------------- #

function write_multiscenario_outputs(dir, c, art, selection, selected, month_indices, cfg)
    mkpath(dir)
    mva = c.base.baseMVA
    r, A, bench, val = art.r, art.A, art.bench, art.val
    # A single-scenario case has no calendar: its scenario IDs are plain indices,
    # not hour-of-year offsets. Timestamps are still emitted (so the column set
    # never changes) but they are only meaningful when a calendar was supplied.
    has_calendar = haskey(cfg, :month) && haskey(cfg, :year)
    first_hour = DateTime(get(cfg, :year, 2017), 1, 1)

    writedlm(joinpath(dir, "assignment_matrix.csv"), A, ',')
    writedlm(joinpath(dir, "monthly_reduced_flow_mw.csv"), mva .* bench.screen.flow, ',')
    writedlm(joinpath(dir, "selected_internal_transfer_mw.csv"), mva .* r.gint, ',')

    # Per-scenario solve cost, full vs reduced -- the computational payoff the
    # whole reduction exists to deliver. Work units are the deterministic column;
    # wall clock is real time but load-sensitive, and solve_time is Gurobi's own
    # Runtime, which reads 0.0 on networks small enough to solve in under a ms.
    timing = get(art, :timing, nothing)
    if !isnothing(timing)
        open(joinpath(dir, "solve_time_benchmark.csv"), "w") do io
            println(io, "scenario_id,timestamp,optimal_full,optimal_reduced,",
                        "solve_time_full_s,solve_time_reduced_s,",
                        "wall_full_ms,wall_reduced_ms,work_full,work_reduced,",
                        "iters_full,iters_reduced,work_speedup,wall_speedup,",
                        "objective_full,objective_reduced")
            for (h, s) in enumerate(timing.scenario_indices)
                ts = has_calendar ? first_hour + Hour(c.scenario_ids[s] - 1) : ""
                println(io, c.scenario_ids[s], ",", ts, ",",
                        Int(timing.optimal_full[h]), ",", Int(timing.optimal_reduced[h]), ",",
                        timing.solve_time_full[h], ",", timing.solve_time_reduced[h], ",",
                        1000 * timing.wall_full[h], ",", 1000 * timing.wall_reduced[h], ",",
                        timing.work_full[h], ",", timing.work_reduced[h], ",",
                        timing.iters_full[h], ",", timing.iters_reduced[h], ",",
                        timing.work_speedup[h], ",", timing.wall_speedup[h], ",",
                        timing.objective_full[h], ",", timing.objective_reduced[h])
            end
        end
        open(joinpath(dir, "solve_time_summary.csv"), "w") do io
            println(io, "metric,full,reduced,speedup")
            println(io, "buses,", timing.dims_full.n_bus, ",", timing.dims_reduced.n_bus,
                    ",", timing.dims_full.n_bus / max(timing.dims_reduced.n_bus, 1))
            println(io, "lines,", timing.dims_full.n_line, ",", timing.dims_reduced.n_line,
                    ",", timing.dims_full.n_line / max(timing.dims_reduced.n_line, 1))
            println(io, "solve_time_s,", timing.total_solve_time_full, ",",
                    timing.total_solve_time_reduced, ",",
                    timing.solve_time_resolved ? timing.speedup_solve_time : NaN)
            println(io, "wall_s,", timing.total_wall_full, ",", timing.total_wall_reduced,
                    ",", timing.speedup_wall)
            println(io, "work_units,", timing.total_work_full, ",", timing.total_work_reduced,
                    ",", timing.speedup_work)
            println(io, "simplex_iterations,", sum(timing.iters_full), ",",
                    sum(timing.iters_reduced), ",", timing.speedup_iters)
        end
    end

    open(joinpath(dir, "bus_mapping.csv"), "w") do io
        println(io, "bus_id,representative_bus_id,retained")
        for j in 1:c.base.N
            rep = r.rep_of[j]
            println(io, "$(c.bus_ids[j]),$(c.bus_ids[rep]),$(j == rep ? 1 : 0)")
        end
    end

    congested_set = Set(selection.congested_lines)
    protected_set = Set(findall(r.protected))
    open(joinpath(dir, "line_status.csv"), "w") do io
        println(io, "line,from_bus_id,to_bus_id,internal,protected,congested_during_month")
        for l in 1:c.base.Ln
            println(io, "$l,$(c.bus_ids[c.base.Efrom[l]]),$(c.bus_ids[c.base.Eto[l]]),",
                    "$(r.c[l]),$(l in protected_set ? 1 : 0),",
                    "$(l in congested_set ? 1 : 0)")
        end
    end

    open(joinpath(dir, "monthly_benchmark.csv"), "w") do io
        println(io, "scenario_id,timestamp,selected,congested_line_count,window_feasible,",
                "max_window_violation_mw,max_external_normalized_error,",
                "max_all_line_normalized_error")
        selected_set = Set(selected)
        for (h, s) in enumerate(month_indices)
            println(io, "$(c.scenario_ids[s]),$(first_hour + Hour(c.scenario_ids[s] - 1)),",
                    "$(s in selected_set ? 1 : 0),",
                    "$(selection.congested_count_by_hour[h]),",
                    "$(bench.window_feasible_by_scenario[h] ? 1 : 0),",
                    "$(mva * bench.screen.scenario_violation[h]),",
                    "$(bench.max_external_error_by_scenario[h]),",
                    "$(bench.max_all_line_error_by_scenario[h])")
        end
    end

    open(joinpath(dir, "monthly_line_benchmark.csv"), "w") do io
        println(io, "line,max_absolute_flow_error_mw,max_normalized_flow_error,",
                "max_relative_window_violation")
        for l in 1:c.base.Ln
            println(io, "$l,$(mva * maximum(bench.absolute_error[l, :])),",
                    "$(maximum(bench.normalized_error[l, :])),",
                    "$(maximum(bench.relative_window_violation[l, :]))")
        end
    end

    if !isnothing(val)
        open(joinpath(dir, "monthly_dcopf_validation.csv"), "w") do io
            println(io, "scenario_id,original_status,reduced_status,",
                    "original_optimal,reduced_optimal,dispatch_feasible_on_original,",
                    "original_objective,reduced_objective,objective_change_pct,",
                    "objective_within_tolerance,max_overload_mw,max_relative_overload,",
                    "max_full_network_utilization,total_overload_mw,n_overloaded_lines,",
                    "max_lmp_error_per_mwh,load_weighted_lmp_error_per_mwh,",
                    "repair_feasible,repair_redispatch_mw,repair_cost_gap_pct,",
                    "in_milp_training_set,n_overloaded_internal_lines,",
                    "n_overloaded_external_lines,total_load_mw,",
                    "n_spuriously_binding_lines")
            for h in eachindex(val.scenario_ids)
                println(io, join((val.scenario_ids[h], val.original_status[h],
                    val.reduced_status[h], val.original_optimal[h] ? 1 : 0,
                    val.reduced_optimal[h] ? 1 : 0, val.dispatch_feasible[h] ? 1 : 0,
                    val.original_objective[h], val.reduced_objective[h],
                    val.objective_change_pct[h], val.objective_within_tolerance[h] ? 1 : 0,
                    mva * val.max_overload[h], val.max_relative_overload[h],
                    val.max_utilization[h], mva * val.total_overload[h],
                    val.n_violating_lines[h], val.max_lmp_error[h],
                    val.load_weighted_lmp_error[h], val.repair_feasible[h] ? 1 : 0,
                    mva * val.repair_redispatch[h], val.repair_cost_pct[h],
                    val.in_training[h] ? 1 : 0, val.n_violating_internal[h],
                    val.n_violating_external[h], mva * val.total_load[h],
                    val.n_spurious_binding[h]), ','))
            end
        end

        # Per-line failure frequency: which corridors to protect next.
        open(joinpath(dir, "line_overload_summary.csv"), "w") do io
            println(io, "line,from_bus_id,to_bus_id,internal,protected,",
                    "dcopf_overload_hours,dcopf_worst_overload_mw,",
                    "benchmark_overload_hours_vs_rating,",
                    "spuriously_binding_hours,worst_reduced_utilization")
            protected_set = Set(findall(r.protected))
            for l in 1:c.base.Ln
                println(io, join((l, c.bus_ids[c.base.Efrom[l]], c.bus_ids[c.base.Eto[l]],
                    val.internal_line[l] ? 1 : 0, l in protected_set ? 1 : 0,
                    val.line_violation_hours[l], mva * val.line_worst_overload[l],
                    bench.overload_hours_by_line[l],
                    val.line_spurious_hours[l],
                    val.line_worst_reduced_utilization[l]), ','))
            end
        end
    end

    # The scenario set the MILP was finally solved on, and how it got there.
    open(joinpath(dir, "active_scenarios.csv"), "w") do io
        println(io, "scenario_id,timestamp,from_seed,total_load_mw,total_abs_flow_mw")
        seed_set = Set(selected)
        for s in art.active
            println(io, join((c.scenario_ids[s],
                first_hour + Hour(c.scenario_ids[s] - 1),
                s in seed_set ? 1 : 0,
                mva * sum(c.load[:, s]),
                mva * sum(abs, c.fhat[:, s])), ','))
        end
    end
    if !isnothing(art.gen)
        open(joinpath(dir, "scenario_generation_history.csv"), "w") do io
            println(io, "iteration,n_active,n_retained,n_internal_lines,solve_time_s,",
                    "worst_violation_above_mw,worst_violation_below_mw,",
                    "worst_line_above,worst_line_below,n_violating_hours,",
                    "both_sides,structural_fix,added_scenario_ids")
            for h in art.gen.history
                println(io, join((h.iteration, h.n_active, h.n_retained, h.n_internal,
                    h.solve_time, mva * h.worst_above, mva * h.worst_below,
                    h.worst_line_above, h.worst_line_below, h.n_violating_hours,
                    h.both_sides ? 1 : 0, h.structural ? 1 : 0,
                    "\"" * join(h.added_ids, ";") * "\""), ','))
            end
        end
    end

    open(joinpath(dir, "README.txt"), "w") do io
        println(io, has_calendar ?
            "Month: $(monthname(cfg.month)) $(cfg.year)" :
            "Case: single operating point (S = 1), no calendar")
        println(io, "Relaxation: $(relaxation_pretty(art.mode, art.delta))")
        println(io, "Congestion definition: abs(flow)/rating >= $(cfg.near_limit_threshold)")
        println(io, "Operating-point/load scale: $(get(cfg, :operating_point_scale, 1.0))")
        println(io, "Line-limit scale: $(get(cfg, :line_limit_scale, 1.0))")
        println(io, "MILP scenario IDs: $(join(selection.scenario_ids, ", "))")
        println(io, "Those cover all $(length(selection.congested_lines)) congested lines.")
        println(io, "Fixed-injection benchmark scenario IDs: ", join(c.scenario_ids[month_indices], ", "))
        println(io, "DC-OPF replay scenario IDs: ",
                isnothing(val) ? "not run" : join(val.scenario_ids, ", "))
        println(io, "Replay checks returned dispatches only; alternate optima and other demands are not certified.")
        println(io)
        println(io, "A congestion relaxation loosens the flow pin on congested lines to buy")
        println(io, "reduction. :conservative opens the window only toward MORE congestion, so")
        println(io, "the reduced network can over-state but never under-state a corridor's")
        println(io, "loading. :symmetric opens it both ways -- more reduction, but it may")
        println(io, "under-state congestion and so overload the line on the full network.")
        println(io, "delta is far larger than the DC-OPF overload tolerance, so the strict")
        println(io, "pass/fail is EXPECTED to fail; compare max_full_network_utilization and")
        println(io, "repair_redispatch_mw across settings instead.")
    end
    return dir
end

# --------------------------------------------------------------------------- #
# Legacy A-based cross-check.
#
# Kept because it is an INDEPENDENT test: the edge model has no assignment
# matrix and no connectivity constraints, so nothing structurally guarantees its
# answer would satisfy the older A-based model (36)/(42). This re-derives each
# of those requirements from the solved (A, c) instead of assuming them.
#
# Single-operating-point only -- it recovers ONE angle vector from A, so it is
# wired up by the single-scenario runner and left off by the multi-scenario one.
# --------------------------------------------------------------------------- #
"""
Would this (A, c) be accepted by the ORIGINAL A-based model (36)/(42)?

The edge model has no assignment matrix and no connectivity constraints, so
nothing structurally guarantees that its answer satisfies what an A-based model
would demand. This checks each requirement directly instead of assuming it.
"""
function check_original_model_constraints(TR, c, A, r, epsL;
                                          near_limit_threshold=nothing,
                                          congestion_relaxation=0.0,
                                          congestion_relaxation_mode::Symbol=:none,
                                          tol::Real=1e-6)
    rep = TR.extract_reduction(A).rep_of
    cl_from_A = [rep[c.Efrom[l]] == rep[c.Eto[l]] ? 1 : 0 for l in 1:c.Ln]

    assignment_ok = all(vec(sum(A, dims=1)) .== 1) &&
                    all(A[i, i] in (0, 1) for i in 1:c.N) &&
                    all(A[i, j] <= A[i, i] for i in 1:c.N, j in 1:c.N)

    members = Dict{Int,Vector{Int}}()
    for j in 1:c.N
        push!(get!(members, rep[j], Int[]), j)
    end
    function connected(mem)
        length(mem) <= 1 && return true
        mset = Set(mem)
        idx = Dict(b => i for (i, b) in enumerate(mem))
        g = SimpleGraph(length(mem))
        for l in 1:c.Ln
            u, v = c.Efrom[l], c.Eto[l]
            (u in mset && v in mset) && add_edge!(g, idx[u], idx[v])
        end
        return length(connected_components(g)) == 1
    end
    connectivity_ok = all(connected(mem) for mem in values(members))

    linking_mismatches = findall(l -> cl_from_A[l] != r.c[l], 1:c.Ln)
    Lplus, Lminus = TR.binding_lines(c; near_limit_threshold=near_limit_threshold)
    binding = sort(unique(vcat(Lplus, Lminus)))
    binding_ok = all(cl_from_A[l] == 0 for l in binding)

    philo, phiup = TR.phi_window(c, epsL, Lplus, Lminus;
        congestion_relaxation=congestion_relaxation,
        congestion_relaxation_mode=congestion_relaxation_mode)
    E = TR.incidence_matrix(c)
    Dx = Diagonal(c.Dx)
    theta = pinv(Matrix(A * E * Dx * E' * A')) * (A * c.p)
    flow = Dx * E' * A' * theta
    ext = findall(cl_from_A .== 0)
    int = findall(cl_from_A .== 1)
    window_ok = all(philo[l] - tol <= flow[l] <= phiup[l] + tol for l in ext) &&
                all(abs(flow[l]) <= tol for l in int)
    g_canon = TR.shorted_internal_flows(c, A)
    # r.G is per (line, scenario). Take each line's loosest bound over the
    # scenarios actually solved -- that is the bound the model applied to it
    # somewhere, so a canonical flow inside it is not evidence of a violation.
    Gline = r.G isa AbstractMatrix ? vec(maximum(r.G, dims=2)) : r.G
    g_ok = all(abs(g_canon[l]) <= Gline[l] + tol for l in int)

    return (assignment_ok=assignment_ok, connectivity_ok=connectivity_ok,
            linking_ok=isempty(linking_mismatches),
            n_linking_mismatches=length(linking_mismatches),
            binding_ok=binding_ok, window_ok=window_ok, g_ok=g_ok,
            all_ok=assignment_ok && connectivity_ok && isempty(linking_mismatches) &&
                   binding_ok && window_ok && g_ok)
end

function report_original_model_constraints(chk; io::IO=stdout)
    section(io, "4) Feasible under the ORIGINAL A-based model's constraints?")
    println(io, "  assignment (5)-(6)              = ", chk.assignment_ok)
    println(io, "  cluster connectivity            = ", chk.connectivity_ok)
    println(io, "  c-linking (14)-(16)             = ", chk.linking_ok,
            chk.linking_ok ? "" : "  ($(chk.n_linking_mismatches) lines mismatched)")
    println(io, "  binding lines stay external (17)= ", chk.binding_ok)
    println(io, "  physical-limit window (22)      = ", chk.window_ok)
    println(io, "  internal-transfer bound (25)    = ", chk.g_ok)
    println(io, "  ALL CONSTRAINTS SATISFIED       = ", chk.all_ok)
    return nothing
end

end # module

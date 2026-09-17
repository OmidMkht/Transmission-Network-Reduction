# --------------------------------------------------------------------------- #
# UNFOLD -- compose a Kron elimination (kron_core.jl) with a solved boundary
# clustering into a full-size assignment matrix, so every EXISTING validation/
# reporting function (model_feasibility_check_multiscenario,
# benchmark_reduction_scenarios, validate_reduced_dcopf_scenarios, ...) can run
# completely unmodified against the true original case.
# --------------------------------------------------------------------------- #

"""
    unfold_kron_assignment(base_full, kron_map, A_boundary) -> A_full::Matrix{Float64}

Composes the elimination map with a solved boundary clustering into a full
N_full x N_full assignment matrix, in the same convention
assignment_from_line_status/extract_reduction use elsewhere (A[i,i]=1 for
every final representative, exactly one A[rep,j]=1 per bus j).

Rule: if a chain's equivalent line ended up INTERNAL (its two anchors landed
in the same boundary cluster), every interior bus joins that cluster too
(lossless -- internal structure is already arbitrary, absorbed by g^int). If
it stayed EXTERNAL, every interior bus becomes its own singleton cluster
(lossless -- nothing was actually discarded, only algebraically folded for
the solve; both original lines per bus are simply never merged, exactly as
if Kron reduction had never touched that chain).

Buses untouched by Kron reduction -- including ring and self-loop-chain
buses, which were never removed from the boundary set in the first place --
need no special case: they are ordinary boundary buses under `full_of`.
"""
function unfold_kron_assignment(base_full::TxReductionCase, kron_map, A_boundary)
    base_full.N == kron_map.full_N ||
        error("kron_map was built from a different case ($(kron_map.full_N) buses) " *
              "than base_full ($(base_full.N) buses)")
    red_b = extract_reduction(round.(Int, A_boundary))
    N = kron_map.full_N
    rep_of_full = zeros(Int, N)
    for bb in 1:kron_map.boundary_N
        rep_of_full[kron_map.full_of[bb]] = kron_map.full_of[red_b.rep_of[bb]]
    end
    for ch in kron_map.chains
        x_b, y_b = kron_map.boundary_of[ch.x], kron_map.boundary_of[ch.y]
        if red_b.rep_of[x_b] == red_b.rep_of[y_b]
            rep = rep_of_full[ch.x]
            for b in ch.interior
                rep_of_full[b] = rep
            end
        else
            for b in ch.interior
                rep_of_full[b] = b
            end
        end
    end
    A_full = zeros(Float64, N, N)
    for j in 1:N
        A_full[rep_of_full[j], j] = 1.0
    end
    return A_full
end

"""
    kron_display_assignment(base_full, kron_map, A_boundary) -> A_display::Matrix{Float64}

Same shape and convention as `unfold_kron_assignment`, but every chain's
interior buses ALWAYS fold into `ch.x`'s representative, regardless of
whether the equivalent line ended up internal or external at the boundary
level. A chain's last original leg (its final interior bus -> `y`) is then
automatically the one surviving line connecting `x`'s cluster to `y`'s --
no separate bookkeeping needed, since `full_line_view` derives internal/
external purely from `rep_of`.

CONTRACT: for REPORTING/DISPLAY only (n_retained, the plot, the CSVs) --
never pass this to `benchmark_reduction_scenarios` or
`validate_reduced_dcopf_scenarios`. Those hard-require `unfold_kron_assignment`'s
full reinsertion, because a chain's individual original lines only have a
well-defined flow/window once actually reinstated; this collapsed version
relies entirely on chain eligibility already having proven those lines are
never near their rating (see kron_core.jl's kron_eligible_buses), which is
true by construction but is not something this function re-derives.
"""
function kron_display_assignment(base_full::TxReductionCase, kron_map, A_boundary)
    base_full.N == kron_map.full_N ||
        error("kron_map was built from a different case ($(kron_map.full_N) buses) " *
              "than base_full ($(base_full.N) buses)")
    red_b = extract_reduction(round.(Int, A_boundary))
    N = kron_map.full_N
    rep_of_full = zeros(Int, N)
    for bb in 1:kron_map.boundary_N
        rep_of_full[kron_map.full_of[bb]] = kron_map.full_of[red_b.rep_of[bb]]
    end
    for ch in kron_map.chains
        rep = rep_of_full[ch.x]
        for b in ch.interior
            rep_of_full[b] = rep
        end
    end
    A_display = zeros(Float64, N, N)
    for j in 1:N
        A_display[rep_of_full[j], j] = 1.0
    end
    return A_display
end

"""
    full_line_view(c_full, rep_of_full; near_limit_threshold, protection_indices)
      -> (internal::BitVector, protected::BitVector)

Per-FULL-line status, at full-Ln granularity.

`internal` is a pure function of the final cluster assignment (exactly the
same way screen_reduction_scenarios/benchmark_reduction_scenarios already
independently recompute line status from A rather than trusting a solver's
raw output) -- correct for every full line, including former chain-interior
lines that never existed as boundary MILP lines at all.

`protected` is recomputed fresh on the FULL network via multiscenario_windows'
own criterion, not looked up through the boundary solve -- so no full-line/
boundary-line bookkeeping is needed, and it is automatically correct
everywhere. A chain-touched line is provably always false here: chain
eligibility already required every constituent line's utilization under
`near_limit_threshold`.
"""
function full_line_view(c_full::MultiScenarioTxReductionCase, rep_of_full::AbstractVector;
                        near_limit_threshold, protection_indices=axes(c_full.p, 2))
    base = c_full.base
    internal = falses(base.Ln)
    for l in 1:base.Ln
        internal[l] = rep_of_full[base.Efrom[l]] == rep_of_full[base.Eto[l]]
    end
    # epsL is irrelevant to `protected` (computed from fhat/frate/threshold
    # alone) -- a placeholder vector avoids duplicating that criterion here.
    win = multiscenario_windows(c_full, ones(base.Ln);
        near_limit_threshold=near_limit_threshold, protection_indices=protection_indices)
    return internal, win.protected
end

"""
    make_r_display(c_full, r_boundary, kron_map, A_full; near_limit_threshold,
                   protection_indices=axes(c_full.p,2))

r_display = merge(r_boundary, full-granularity overrides for
n_retained/retained/rep_of/A/c/n_internal_lines/n_external_lines/protected).

CONTRACT: for full-network DISPLAY only -- `report_reduction`,
`report_dcopf_failures`, `plot_network`'s
`binding_lines=findall(r_display.protected)`. Never pass this to
`model_feasibility_check_multiscenario` or `check_original_model_constraints`:
the remaining fields (f, gint, vartheta, philo, phiup, pinned, exact_binding,
G, scenario_indices) stay BOUNDARY-shaped after the merge and are only
self-consistent paired with `c_boundary` -- use `(c_boundary, r_boundary)`
for that check instead, reported separately as "boundary-solve feasibility".
"""
function make_r_display(c_full::MultiScenarioTxReductionCase, r_boundary, kron_map, A_full;
                        near_limit_threshold, protection_indices=axes(c_full.p, 2))
    red = extract_reduction(round.(Int, A_full))
    internal, protected = full_line_view(c_full, red.rep_of;
        near_limit_threshold=near_limit_threshold, protection_indices=protection_indices)
    return merge(r_boundary, (
        n_retained=red.n_retained, retained=red.retained, rep_of=red.rep_of,
        A=A_full, c=Int.(internal),
        n_internal_lines=count(internal), n_external_lines=length(internal) - count(internal),
        protected=protected,
    ))
end

"""
    solve_and_unfold_kron(c_full, epsL_full, c_boundary, epsL_boundary, kron_map;
                          solve_kwargs...) -> (r_boundary, A_full, r_display)

Thin driver: solves the boundary MILP with every keyword forwarded verbatim
(no remapping needed -- Kron reduction only removes buses/lines, never
scenario columns, so `c_boundary.p` has exactly `c_full.p`'s S columns in the
same order), then unfolds. `near_limit_threshold`/`protection_indices` are
read out of `solve_kwargs` (without removing them, so they still reach the
solve call) to keep the display step's protection criterion consistent with
what the boundary solve itself used.
"""
function solve_and_unfold_kron(c_full::MultiScenarioTxReductionCase, epsL_full,
                               c_boundary::MultiScenarioTxReductionCase, epsL_boundary,
                               kron_map; solve_kwargs...)
    r_boundary = solve_reduction_edge_multiscenario(c_boundary, epsL_boundary; solve_kwargs...)
    A_full = unfold_kron_assignment(c_full.base, kron_map, r_boundary.A)
    near_limit_threshold = get(solve_kwargs, :near_limit_threshold, nothing)
    protection_indices = get(solve_kwargs, :protection_indices, axes(c_full.p, 2))
    r_display = make_r_display(c_full, r_boundary, kron_map, A_full;
        near_limit_threshold=near_limit_threshold, protection_indices=protection_indices)
    return r_boundary, A_full, r_display
end

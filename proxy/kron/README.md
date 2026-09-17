# Kron preprocessing

The reduction MILP carries one binary per line. Eliminating eligible degree-two
chains reduces its size, although solve difficulty also depends on the
formulation, scenarios and relaxation strength. This preprocessing groups the
chain's contraction decisions into one decision.

This directory collapses every such chain into one equivalent series line
before the MILP runs, solves the MILP on the smaller "boundary" network, and
unfolds the result back onto the full bus set. The elimination is a Schur
complement of the susceptance matrix, so boundary angles and boundary line
flows are unchanged in every scenario -- and because the clustering is
unfolded to full size, the feasibility check, the flow-window benchmark and
the reduced DC-OPF validation all run against the true original network,
unmodified. Nothing is taken on trust from the elimination.

`reference/kron_preprocessing.pdf` is the formulation: the Schur complement,
the closed-form injection split, the eligibility rules, the unfolding, and
what the whole composition costs.

## Running it

Set `kron = true` in `proxy/run_proxy.jl`, or pass it on the command line:

```
julia --project=. --startup-file=no proxy/run_proxy.jl kron=true
```

Writing to `outputs/case118_kron_edge/`. Setting it back to `false` gives the
plain run through the same code path and the same settings, which is how the
two are meant to be compared -- one flag is the only difference between them.
It composes with `scenarios = :multi` the same way.

## Same scope limit as the parent

The Schur complement preserves the DC equations for the specified injections;
this does not establish thermal-limit preservation or OPF dispatch safety for
the subsequent clustering. Full-network flow-window validation and independent
OPF redispatch checks must pass after unfolding. Each check covers only its
specified scenarios and returned dispatches. The eligibility screen uses the
supplied fixed injections, so eliminated lines can become congested after
redispatch, even at the same demand.

## Which buses are eliminated

A bus qualifies only if all of the following hold:

| Rule | Why |
|---|---|
| degree exactly 2 | the Schur complement is then a single series line, with no fill-in |
| no generator attached | nothing with a dispatch variable is ever eliminated |
| not the angle reference | avoids remapping the reference into the boundary network |
| its two lines reach two distinct neighbours | a parallel circuit ending at one neighbour has no two-anchor series equivalent |
| both lines under `near_limit_threshold` in every supplied scenario | preserves the proxy's protected-line set at the screened injections; it is not a redispatch-safety certificate |

Cycles lying entirely inside the eligible set, and paths whose two ends are
the same bus, have no series equivalent and are left untouched.

## The trade-off

A chain becomes one binary, so it is merged or not as a unit -- the MILP can
no longer absorb a chain's first bus and leave the second. The reachable set
of clusterings is restricted. Equivalent-line ratings and windows do not in
general preserve all original chain limits, so only results that pass the same
full-network checks should be compared. What is bought is a smaller MILP.

Which side wins depends on the case and the time budget, and it does change
sign between cases, so it is measured rather than assumed. Run the same case
both ways -- `kron = true` and `kron = false` change nothing else -- and
compare `relaxation_comparison.csv`.

A smaller model can reduce construction time or reach a better incumbent under
the same time budget. At a proved optimum under the same bus-count objective
and original-network constraints, restricting the available clusterings cannot
improve the optimum. The current line-count objective and equivalent-line
windows require care when interpreting this comparison.

The primary `n_retained`, `art.A`, and `assignment_matrix.csv` now refer to the
unfolded network used for validation and timing. `A_display` and
`display_assignment_matrix.csv` preserve the optional compact display separately.
Older outputs with `collapse_external_chains=true` could report fewer buses than
the validated network actually had; those counts need to be recomputed before
comparison.

## Files

| File | Role |
|---|---|
| `kron_core.jl` | Eligibility, chain detection, Schur complement, boundary case |
| `kron_unfold.jl` | Boundary clustering -> full-size assignment matrix |
| `kron_reporting.jl` | Sweep driver, console reports, CSV output |
| `kron_plots.jl` | Three-panel figure: chains, boundary network, reduced network |

Every validation and benchmark call inside `kron_reporting.jl` is the existing
function from the parent directory, unmodified, run against the true full
case. Nothing in the parent directory is modified by anything here.

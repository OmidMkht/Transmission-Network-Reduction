# Kron preprocessing

The reduction MILP carries one binary per line, so its difficulty is set by
the size of the network it is handed -- and much of that network carries no
decision. A bus of degree two with no generator and no nearby congestion
contributes a binary whose only effect is to lengthen a chain.

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

Set `kron = true` in the `RUN` block of `../run_tnr.jl`:

```
julia --project=. --startup-file=no run_tnr.jl
```

Writing to `outputs/case118_kron_edge/`. Setting it back to `false` gives the
plain run through the same code path and the same settings, which is how the
two are meant to be compared -- one flag is the only difference between them.
It composes with `scenarios = :multi` the same way.

## Which buses are eliminated

A bus qualifies only if all of the following hold:

| Rule | Why |
|---|---|
| degree exactly 2 | the Schur complement is then a single series line, with no fill-in |
| no generator attached | nothing with a dispatch variable is ever eliminated |
| not the angle reference | avoids remapping the reference into the boundary network |
| its two lines reach two distinct neighbours | a dead-end into the same neighbour twice is a parallel circuit, and eliminating it would need a shunt |
| both lines under `near_limit_threshold` in every scenario | this is the safety rule: no eliminated line is ever one the reduction had to keep external, so a chain cannot hide a congested corridor |

Cycles lying entirely inside the eligible set, and paths whose two ends are
the same bus, have no series equivalent and are left untouched.

## The trade-off

A chain becomes one binary, so it is merged or not as a unit -- the MILP can
no longer absorb a chain's first bus and leave the second. The reachable set
of clusterings is a strict subset of the full model's, so the reduction at a
given tolerance is weakly worse. What is bought is a smaller MILP.

Whether that pays depends entirely on the case and the time budget, so it is
measured rather than assumed. On the bundled case118 it does not pay:

| case118, eps = 10%, no relaxation | plain | Kron |
|---|---|---|
| buses retained | 33 / 118 (72.0%) | 40 / 118 (66.1%) |
| MILP solve time | 12.3 s | 6.1 s |
| reduced DC-OPF speedup (work units) | 3.03x | 1.98x |

case118 halves a MILP time that was never the bottleneck, and pays with a
coarser reduction and a weaker downstream speedup. The elimination earns its
place on networks where the full MILP does not finish -- on ACTIVSg2000 the
boundary MILP reaches a materially better incumbent inside the same one-hour
budget. Run both and read `relaxation_comparison.csv` before assuming either
way.

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

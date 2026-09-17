# Reduced network examples

Three networks reduced by the edge-based MILP in this repo (`reference/transmission_network_reduction.pdf`
-- plain aggregation, no Kron preprocessing), exported as standalone MATPOWER
`.m` files. Each loads directly in MATPOWER or PowerModels; none needs
anything else from this repo.

| File | Source case | Buses retained | MILP | DC-OPF check |
|---|---|---|---|---|
| `reduced_pglib_opf_case118_ieee.m` | pglib_opf_case118_ieee | 33 / 118 (72.0%) | `OPTIMAL`, 12 s | 1/1 feasible |
| `reduced_pglib_opf_case500_goc.m` | pglib_opf_case500_goc | 16 / 500 (96.8%) | `OPTIMAL`, 919 s | 0/1 feasible -- see caveat below |
| `reduced_case_ACTIVSg200.m` | ACTIVSg200, March 2017, 100 seed hours | 128 / 200 (36.0%) | `TIME_LIMIT` (20 min), incumbent | 168/168 returned dispatches feasible (one week) |

Same settings for all three: `normalized_error_threshold = 0.10` (each external
line's flow window is +-10% of its rating), `near_limit_threshold = 0.80`, and
`relaxation_sweep = [(:none, 0.00)]` -- the strict, unrelaxed pin: a line at or
above 80% of its rating anywhere in the scenario set is forced external and
reproduced exactly. Other external lines retain their configured flow windows.

## Caveat on case500_goc

The MILP solved to a *proven* optimum -- not a time-out -- and the reduced
network is 100% window-feasible. But at 96.8% reduction it is aggressive
enough that the strict DC-OPF re-dispatch overloads one line by about 2.1%
above rating. It therefore fails the stated dispatch-feasibility check, despite
satisfying the proxy model. The repair-economics check recorded 24 MW of
redispatch and a 0.043% increase in generation cost. That separate repair does
not make the exported network's unmodified OPF dispatch safe.
`pglib_opf_case300_ieee` was attempted at the same
settings and dropped from this batch: its MILP did not converge in 20 minutes
and its DC-OPF gap was neither small nor proven-optimal, unlike case500_goc's.

## Scenario and validation scope

The MILP constrains fixed-injection flow errors. The DC-OPF column reports
independently optimized dispatches checked on the original network; it does not
certify all optimal or all feasible dispatches. The case500 failure demonstrates
that even a modeled scenario can fail this independent check.

`reduced_case_ACTIVSg200.m` was reduced against 100 seed hours and validated
across 168 hours (one week), with every returned dispatch feasible on the true
full network. These are sampled checks, not continuous-envelope certification.

`reduced_pglib_opf_case118_ieee.m` and `reduced_pglib_opf_case500_goc.m` were
reduced against a **single** operating point, because a pglib `.m` file carries
exactly one. Case118 passed its one scenario; case500 failed it. Measured
on case118: hold the clustering fixed and scale the load, and the strict DC-OPF
check fails at +-2.5%, worst utilisation 1.056. Reducing the same case against
two operating points spanning +-10% gives 108 / 118 buses instead of 33 / 118,
and that clustering passed the tested load levels within the range. No
continuous-interval guarantee was established by those samples.

Use the two single-point files as examples of what the exporter produces and of
what a reduced network looks like -- not as networks to run studies on.

## What's preserved, what isn't

Every line NOT collapsed by the reduction is copied verbatim from the source
`.m` file -- resistance, reactance, susceptance, tap ratio, phase shift,
rating, in-service status. Nothing about a surviving line's data is touched,
because reduction never modifies a surviving line, only which buses its ends
attach to. Bus shunts and loads belonging to an eliminated bus are moved onto
its cluster's representative bus (summed, not dropped); generators keep their
own row, limits and cost, only their bus is redirected. Buses keep their
original external MATPOWER numbers (not renumbered), so `bus_i` in these
files is directly comparable to the source case.

One pre-existing limitation of the whole pipeline, not something this export
introduces: this project's own DC model (`common/preprocessing.jl`) ignores
transformer tap ratio and phase shift when computing flows (`Dx = 1/br_x`
only). A tap/shift-aware DC or AC solve run on these files elsewhere may
therefore see slightly different flows than what this repo validated
internally.

## How these were produced

`common/matpower_export.jl`'s `export_reduced_matpower`, called from `proxy/run_proxy.jl`
with `show.export_matpower = true` -- `scenarios = :single` for case118 and
case500_goc, `scenarios = :multi` for ACTIVSg200. ACTIVSg200's own scenario-generation
pipeline is not published in this repo (see the root README's Kron section
for the same pattern: derived artifact published, generating pipeline kept
local) -- only the reduced network itself.

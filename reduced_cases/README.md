# Reduced network examples

Three networks reduced by the edge-based MILP in this repo (`reference/transmission_network_reduction.pdf`
-- plain aggregation, no Kron preprocessing), exported as standalone MATPOWER
`.m` files. Each loads directly in MATPOWER or PowerModels; none needs
anything else from this repo.

| File | Source case | Buses retained | MILP | DC-OPF check |
|---|---|---|---|---|
| `reduced_pglib_opf_case118_ieee.m` | pglib_opf_case118_ieee | 33 / 118 (72.0%) | `OPTIMAL`, 12 s | 1/1 feasible |
| `reduced_pglib_opf_case500_goc.m` | pglib_opf_case500_goc | 16 / 500 (96.8%) | `OPTIMAL`, 919 s | 0/1 feasible -- see caveat below |
| `reduced_case_ACTIVSg200.m` | ACTIVSg200, March 2017, 100 seed hours | 128 / 200 (36.0%) | `TIME_LIMIT` (20 min), genuine | 168/168 feasible (whole month) |

Same settings for all three: `normalized_error_threshold = 0.10` (each external
line's flow window is +-10% of its rating), `near_limit_threshold = 0.80`, and
`relaxation_sweep = [(:none, 0.00)]` -- the strict, unrelaxed pin: a line at or
above 80% of its rating anywhere in the scenario set is forced external and
reproduced exactly, no window slack is granted anywhere.

## Caveat on case500_goc

The MILP solved to a *proven* optimum -- not a time-out -- and the reduced
network is 100% window-feasible. But at 96.8% reduction it is aggressive
enough that the strict DC-OPF re-dispatch overloads one line by about 2.1%
above rating. This is not evidence the reduction is broken: the codebase's own
repair-economics check quantifies it as a cheap fix -- 24 MW of redispatch, a
0.043% increase in generation cost -- and every other line and every window
constraint holds. Read it as "reduced networks built with an exact congestion
pin can leave a small, repairable gap at very high reduction ratios," not as a
validated dispatch. `pglib_opf_case300_ieee` was attempted at the same
settings and dropped from this batch: its MILP did not converge in 20 minutes
and its DC-OPF gap was neither small nor proven-optimal, unlike case500_goc's.

case118 and the ACTIVSg200 month are both fully validated: every scenario's
reduced-network optimal dispatch is feasible on the true full network, not
merely inside the MILP's own flow windows.

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
introduces: this project's own DC model (`tnr_preprocessing.jl`) ignores
transformer tap ratio and phase shift when computing flows (`Dx = 1/br_x`
only). A tap/shift-aware DC or AC solve run on these files elsewhere may
therefore see slightly different flows than what this repo validated
internally.

## How these were produced

`matpower_export.jl`'s `export_reduced_matpower`, called from `run_tnr.jl`
(case118, case500_goc) and the local multi-scenario runner (ACTIVSg200) with
`cfg.show.export_matpower = true`. ACTIVSg200's own scenario-generation
pipeline is not published in this repo (see the root README's Kron section
for the same pattern: derived artifact published, generating pipeline kept
local) -- only the reduced network itself.

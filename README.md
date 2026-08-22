# Transmission Network Reduction

A Julia/JuMP implementation of an edge-based, assignment-matrix-free MILP for
reducing a transmission network: it chooses which lines to collapse (merging
the buses at their ends) so that the reduced network's DC power flow stays
within an error window of the full network's, on every scenario supplied,
while maximizing the number of internal (collapsed) lines.

The reason to reduce a network is computational: a reduced DC-OPF is cheaper
to solve than the full one. `tnr_postprocessing.jl` benchmarks both under
identical solver settings (Gurobi work units, not just wall clock) so that
claim is measured, not assumed, on every run.

## Requirements

- [Julia](https://julialang.org/) 1.11
- [Gurobi](https://www.gurobi.com/) with a valid license (academic licenses
  are free) and the `GUROBI_HOME` / `GRB_LICENSE_FILE` environment variables
  set up per Gurobi's own instructions
- The Julia packages pinned in `Project.toml` / `Manifest.toml`

## Setup

```
git clone <this-repo>
cd "Transmission Network Reduction"
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

`Pkg.instantiate()` installs exactly the package versions recorded in
`Manifest.toml`, so a fresh clone reproduces the environment this code was
developed and tested against.

## Quick start

Two small PGLib-OPF cases (`case14`, `case118`) are bundled so the pipeline
runs with no data download:

```
julia --project=. --startup-file=no run_tnr.jl
```

This builds `case studies/pglib_opf_case118_ieee.m` as a one-scenario case,
solves the reduction MILP, validates the result against a DC-OPF, benchmarks
full-vs-reduced solve time, and writes everything to `outputs/case118_edge/`.

`run_tnr.jl` is the only runner. Two switches at the top of it pick what runs:

```julia
RUN = (
    scenarios = :single,   # :single = one operating point | :multi = hourly matrices
    kron      = false,     # true = collapse degree-2 chains before the MILP
    ...
)
```

Everything below that block -- the model, the tolerances, the sweep, the
validation -- is shared by all four combinations, so a change is tested the
same way whichever one you run. The console prints a short digest; the full
report for each setting lands in `<output_dir>/<setting>/report.txt`.

`scenarios = :multi` runs a month of hourly scenarios instead of one operating
point. It needs `tnr_multiscenario.jl` and the ACTIVSg scenario matrices,
neither of which is published, so a fresh clone runs `:single`; asking for
`:multi` stops with a message naming what is missing.

## Kron preprocessing

`kron_reduction/` collapses every chain of degree-2, generator-free,
uncongested buses into one equivalent series line before the MILP runs, then
unfolds the result back onto the full bus set so every validation still runs
against the true original network. Set `kron = true` in `run_tnr.jl` to use
it; outputs go to `outputs/case118_kron_edge/`.

See [`kron_reduction/README.md`](kron_reduction/README.md) for the eligibility
rules and the measured trade-off, and `reference/kron_preprocessing.pdf` for
the formulation.

## Reduced network examples

[`reduced_cases/`](reduced_cases/) has three networks already reduced by this
repo, exported as standalone MATPOWER `.m` files that load directly in
MATPOWER or PowerModels -- no dependency on this repo or on the (much larger)
original case files. See its README for exactly what settings produced each
one and a caveat on one of them.

## Case data

`case14` and `case118` ship in `case studies/` for the quick start above.
Larger PGLib-OPF cases (`case300`, `case500`, `case2000`, `case6515`, ...) can
be added the same way -- download the `.m` file from
[power-grid-lib/pglib-opf](https://github.com/power-grid-lib/pglib-opf) and
place it under `case studies/<filename>.m`.

## File overview

| File | Role |
|---|---|
| `tnr_preprocessing.jl` | Case structs, case building, DC-OPF, redispatch, windows |
| `tnr_model.jl` | The reduction MILP itself (edge-based model) |
| `tnr_postprocessing.jl` | Feasibility checks, DC-OPF validation, solve-time benchmark |
| `tnr_reporting.jl` | Console reports, plots, CSV output, the sweep driver |
| `transmission_plots.jl` | Before/after network figure |
| `run_tnr.jl` | The runner: single/multi scenario, with/without Kron |
| `matpower_export.jl` | Write a reduced network back out as a MATPOWER `.m` file |
| `kron_reduction/` | Chain elimination before the MILP |
| `reduced_cases/` | Reduced networks published as standalone MATPOWER `.m` files |
| `reference/` | Compiled papers describing the formulations this code implements |

`outputs/` is where every runner writes results -- it's gitignored and
reproducible from the case data and `cfg` blocks above.

## License

MIT -- see `LICENSE`. This covers the code only; third-party case data under
`case studies/` and the paper under `reference/` carry their own terms.

# Transmission Network Reduction

Julia/JuMP code for reducing a transmission network by merging buses, so that a
DC-OPF on the reduced network gives dispatches that still work on the full one.

Three approaches, each with its own runner:

| Folder | Approach | Guarantee |
|---|---|---|
| `kkt/` | Bilevel model solved as one MILP through the lower level's KKT conditions | Some reduced optimum is feasible on the full network at every design demand |
| `proxy/` | Proxy MILP: flows of the reduced network stay within a window of the full network's at fixed injections | None on the re-optimized dispatch; judged by how far off it is |
| `greedy/` | Greedy merges ranked by predicted flow change, each checked by LPs | Same as `kkt/`, per accepted merge |

## Requirements

- [Julia](https://julialang.org/) 1.11+
- [Gurobi](https://www.gurobi.com/) with a license (free for academics)

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## Running

Each runner starts with a `SETTINGS` block. Edit it, or override any setting on
the command line:

```
julia --project=. --startup-file=no kkt/run_kkt.jl
julia --project=. --startup-file=no greedy/run_greedy.jl case=case300 hop_cap=5
julia --project=. --startup-file=no proxy/run_proxy.jl case=case14 demands=scaled
```

Values can be numbers, `true`/`false`, `nothing` or words (`case300`).

### Hop limit

Each approach has a second runner that limits how long a merged chain inside a
cluster may be. `hop_cap` is one number, or a ladder: one solve per entry, each
keeping what the previous one merged. Quote lists in the shell.

```
julia --project=. --startup-file=no kkt/run_kkt_hop.jl case=ACTIVSg200 'hop_cap=[10,20,nothing]'
julia --project=. --startup-file=no greedy/run_greedy_hop.jl case=case300 hop_cap=5
julia --project=. --startup-file=no proxy/run_proxy_hop.jl case=case300 'hop_cap=[5,10,nothing]'
```

A ladder writes `steps.csv`, one row per step, and the usual outputs for the
last step. `time_limit` is per step.

### Settings shared by all three

| Setting | Values |
|---|---|
| `case` | `case14`, `case118`, `case300`, `case500`, `case2000`, `case6515`, `ACTIVSg200`, `ACTIVSg2000` |
| `demands` | `:single` base demand; `:scaled` `n_demands` load levels over `scale_range`, midpoints held out; `:hourly` ACTIVSg hours |
| `budget` | max merged lines, `nothing` = no limit |
| `hop_cap` | max chain of merged lines inside a cluster |
| `size_cap` | max buses per cluster |
| `time_limit` | seconds |
| `output_dir` | `nothing` = `outputs/<approach>/<case>/<tag>/` |

The `_hop` runners drop `budget` and `size_cap` and take `hop_cap` as a number
or a ladder.

### Approach-specific

- **kkt:** `cost_cap` (% over the full-network optimum), `objective`
  (`:lines` or `:clusters`), `start_from` (an `internal.csv` from another run,
  e.g. the greedy, kept merged).
- **greedy:** `cost_cap`, `flow_tol` (allowed overload, fraction of rating),
  `cost_tol`, `kkt_check` (final joint check with rollback), `ordering`,
  `norm`, `alpha`, `radial_first`.
- **proxy:** `eps` (flow window, fraction of rating), `near_limit` (lines
  loaded above this stay unmerged), `relaxation`, `lmp_separation`, `kron`,
  `plots`, `export_matpower`.

### Outputs

Every run writes `settings.txt`, `summary.csv`, the reduction (`internal.csv`,
`assignment.csv`, `bus_mapping.csv`) and its checks. `kkt/` and `greedy/` check
the design demands (and whether every reduced optimum is safe) and the held-out
demands. `proxy/` replays the reduced dispatch on the full network and reports
overloads in MW and % of rating, cost and LMP errors, per setting in
`report.txt`.

## Scope

All checks use a lossless DC model (taps and phase shifts ignored) and hold
only for the demands checked. A reduced OPF with several optimal dispatches may
return one that was not checked; `kkt/` and `greedy/` report how many design
demands have every optimum safe.

## Other folders

| Folder | Contents |
|---|---|
| `common/` | Case table and demands, caps, audit, pre/postprocessing, plots, MATPOWER export |
| `analysis/` | Line-by-line infeasibility of a proxy result, network plot |
| `reduced_cases/` | Reduced networks as standalone MATPOWER files |
| `reference/` | Papers describing the formulations |

## Case data

`case14` and `case118` are included. Other PGLib-OPF cases go in
`case studies/` (from [pglib-opf](https://github.com/power-grid-lib/pglib-opf)).
`demands=hourly` needs the ACTIVSg hourly data and the scripts that build its
scenario matrices, which are not in this repo; `single` and `scaled` work on any
case above.

## License

MIT, see `LICENSE`. Covers the code only; case data and papers carry their own
terms.

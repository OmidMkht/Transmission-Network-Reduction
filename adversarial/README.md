# Per-contingency network reduction for SC-DCOPF

A preventive SC-DCOPF has one network copy for the intact case and one per line
outage. This folder replaces each copy with its own reduced network (buses
merged), designed so that **any dispatch the reduced networks accept is secure on
the full network**, for every demand in the design range and every dispatch
within 5% of the secure optimal cost.

Each reduced network is found by a loop:
1. A MIP proposes the most-merged network the evidence so far allows.
2. An LP searches for an **adversary**: a near-economic dispatch the proposal
   accepts but the full network cannot carry.
3. If one exists, it is added to the MIP with a cut that makes the next proposal
   reject it, and the loop repeats.

**Status (29 September 2026).**
- **ACTIVSg200:** every network is certified; median 4 of 200 buses. Designed on
  July, the reduced SC-DCOPF is exact on July and has no overload on August and
  March.
- **ACTIVSg2000:** runs, but the intact network does not reduce yet (1,550 of
  2,000 buses) with either the dominance cut or the weaker cut
  (`docs/results.md`, Sections 8–9).

## Documents

| File | Contents |
|---|---|
| [docs/method.md](docs/method.md) | the method: every step and constraint, the guarantee and its scope |
| [docs/results.md](docs/results.md) | every run with its numbers |
| [docs/decisions.md](docs/decisions.md) | decisions and reasons, rejected ideas, proposals, open questions |
| `hpc/VACC.md`, section 12 | running on the VACC (local only; `hpc/` is not in git) |

The PDFs `Adversarial_Injection_Generation_for_Safe_Network_Reduction.pdf` and
`Adversarial_Reduction_Report.pdf` are older snapshots (before master v2), kept
out of git. These markdown files replace them.

## Files

| File | Contents |
|---|---|
| `adversarial.jl` | module `Adversarial`: network data, PTDF and LODF, lazy SC-DCOPF, adversary set, screening, bridges, master, checks, design loop, exact limits |
| `evaluate.jl` | module `Evaluate`: SC-DCOPF from network blocks (full or reduced) and from compact PTDF rows; per-hour solves |
| `run_adversarial.jl` | the runner: data, screening, designs, checks, exact limits, evaluation, output files |

## Running

```bash
julia --project=. --startup-file=no adversarial/run_adversarial.jl [key=value ...]

# ACTIVSg200, design on July, test on August
julia --project=. --startup-file=no adversarial/run_adversarial.jl case=ACTIVSg200 \
    month=7 horizon_days=31 test_month=8 test_days=31 planning=false exact_limits=true
```

It needs the hourly scenarios: `common/multiscenario.jl` and
`outputs/<case>_dcopf/`, both gitignored and built by
`common/make_hourly_scenarios.jl` from the ACTIVSg data. ACTIVSg2000's hourly
data covers January–March only.

## Settings

| Setting | Default | Meaning |
|---|---|---|
| `case` | `:case14toy` | test case |
| `month`, `horizon_start_day`, `horizon_days` | 11, 1, 30 | design window (the demand hull); `month` may list several months |
| `test_month`, `test_start_day`, `test_days` | 12, 1, 31 | test window |
| `designs` | `true` | `false`: compare the full, lazy and compact forms only |
| `contingencies` | `true` | `false`: the intact network alone; the cost cap and the witness come from the plain DCOPF, and the evaluation runs in the same job |
| `search` | `:master` | `:master` (master + polishing) or `:greedy` (merge in batches while the exact check passes; no master, every step certified) |
| `band` | 0.05 | cost cap: at most (1 + band) × the interpolated SC-DCOPF optimum |
| `tau` | 1e-6 | overload the checks let pass, fraction of the rating |
| `cut` | `:dominance` | rejection cut: `:dominance` (reduced flow beyond full flow), `:line` (weaker: above the rating), `:any` (any kept line; binaries) |
| `tolerance` | 0.0 | accepted dispatches may overload the full network by at most this fraction; critical pairs and adversaries count only above it, and the dominance cut gets that much room |
| `limits` | `:all` | ratings during the design: `:all` kept lines, or `:critical` lines only (the others unlimited, as after exact limits) |
| `hop_limit` | `nothing` | outage networks: fix every line beyond this many hops from the outage and critical lines to the start design (doubles when nothing fits); `nothing` uses `radius` |
| `cluster_hops` | `[3, 6, nothing]` | searches without a seed (the intact network; an outage network once its neighbourhood or hop region covers everything): one rung per entry, each capping every cluster's chains of merged lines at that many hops, a certified rung's merges held into the next; `nothing` turns it off |
| `ladder_base` | `false` | also use the ladder for the intact network |
| `derate` | 0.0 | kept lines rated (1 − derate) F during the design |
| `margin` | `nothing` | witness at ratings (1 − margin) F; `nothing` = `derate` |
| `raise` | `true` | after certification, raise lowered ratings back as far as the check allows |
| `candidates` | `:all` | lines that may reject under `:any`: `:all` or `:screened` |
| `adversaries_per_round` | 3 | worst adversaries added to the master per round |
| `master_threads` | 4 | Gurobi threads for the master |
| `master_time` | 20 | seconds per master solve |
| `master_gap` | 0.05 | relative gap at which a master solve may stop |
| `polish_time` | 300 | seconds per round to add merges the master missed, checked exactly without the MIP; 0 turns it off |
| `fix_radial` | `true` | fix every non-critical bridge merged |
| `unmerge_hops` | 2 | outage start: unmerge merged lines this close to the outage line |
| `radius` | 20 | the master may change at most this many lines of the outage start (doubles when nothing is found) |
| `network_time` | 600 | seconds per network |
| `max_rounds` | 200 | rounds per network |
| `networks` | `nothing` | design only these networks, e.g. `[0, 3]` (0 = intact) |
| `network_index` | `nothing` | design only network k of [intact; outages], or a list such as `[2,3,4]`, then stop (cluster mode) |
| `start_from` | `nothing` | `masks.csv` whose intact design outage networks start from |
| `load_designs` | `nothing` | folder of per-network parts to evaluate instead of designing |
| `exact_limits` | `false` | also evaluate each design with exact limits on its critical lines only |
| `full_base` | `false` | evaluation: keep the intact network whole, reduce only the outages |
| `planning` | `true` | add each candidate line of `case studies/14bus` in turn (14-bus only) |
| `output_dir` | `nothing` | default `outputs/adversarial/<case>` |

## Output files

| File | Columns / contents |
|---|---|
| `summary.txt` | the log: screening, every design round, every form |
| `settings.txt` | the settings used |
| `hours.csv` | `set, hour, hour_id, sc_dcopf_cost` for the design (`hull`) and test hours |
| `networks.csv` | `network, critical_pairs, redundant, buses, kept_lines, max_derating_pct, status, rounds, adversaries, witness_hours, master_s, check_s, seconds` |
| `masks.csv` | `network, line, internal, rating_over_F`: the designs (internal = 1 is merged) |
| `exact_limits.csv` | `network, limits, limited_lines, kept_lines, critical_limit_over_F, accepts_exact_optimum` |
| `forms.csv` | `form, hours, n_hours, vars, rows, nnz, build_s, solve_s, infeasible_hours, worst_cost_gap_pct, worst_overload_pct` |
| `hours_eval.csv` | `form, hours, hour_id, hull_distance_pct, sc_dcopf_cost, cost, cost_gap_pct, worst_loading, worst_outage, worst_line`: one row per form and hour; the gap is signed, the hull distance is 0 inside the design hull |
| `planning.csv` | the planning-line test (`planning=true`) |
| `critical_<network>.csv`, `master_copies_<network>.csv` | the critical pairs, and every injection the master held (kind, round, line, sign, full loading, hour, then the injection per bus), so its answers can be checked without the MIP |
| `gurobi_master_<network>.log` | Gurobi's log of every master solve for that network (incumbent, bound, gap); the round lines in `summary.txt` also show merges and bound |

## Cluster mode

Every network is designed on its own, so on a cluster each one is a separate job:

1. `network_index=1` designs the intact network into its own folder.
2. `network_index=k start_from=<intact masks.csv>` designs outage network k,
   starting from the intact design. A job that designs only some networks
   screens only those.
3. `load_designs=<parts folder>` reads every part, then checks and evaluates.

`hpc/submit_adversarial.sh` chains the three as SLURM jobs. Environment knobs:
- `CASE`, `WINDOW` (design and test windows as settings);
- `TESTS` (test months as `month:days`, one evaluation each in `eval_m<month>/`), `CAMPAIGN` (runs go to `outputs/adversarial/<case>_<campaign>/<tag>/`);
- `CHUNK` (networks per outage job; SLURM counts every job against a 1,000-job submit limit);
- `NETWORK_TIME`, `MASTER_TIME`;
- `ARRAY` (default `2-174`; SLURM caps indices at 1000) and `ARRAY_LIMIT`;
- `RES` (extra sbatch resources), `NO_EVAL`;
- `EVAL_NAME` for `hpc/adversarial_job.sh` (the evaluation folder name).

Details in `hpc/VACC.md`, section 12.

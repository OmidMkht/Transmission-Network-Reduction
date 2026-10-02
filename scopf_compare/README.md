# Reduced networks against exact SC-DCOPF methods

A standalone comparison: the per-contingency reduced networks of `adversarial/`
against methods that solve the full preventive SC-DCOPF (line outages) without
reducing anything. The methods are separate approaches; none is combined with the
reduction.

| Method | What it is | Exact? |
|---|---|---|
| `btheta_full` | every network (intact and each outage) as a B-θ copy sharing one dispatch | yes |
| `btheta_lazy` | the intact copy; an outage's copy is added once the dispatch violates it | yes |
| `ptdf_lazy` | PTDF rows over the generator outputs, each limit added once violated ("typical PTDF") | yes |
| `ptdf_decomposed` | decomposed PTDF of Alkhraijah, Sigler, Knueven and Maack, "A Sparsification Method for Security-Constrained Optimal Power Flow" (PowerUp 2026): Kron reduction per area, boundary equivalent injections, lazy rows | yes |
| `ptdf_lazy_flows`, `ptdf_decomposed_flows` | the same two in the paper's form: a flow variable per needed line with one defining row, limit rows of one or two flows (the substituted forms write the whole flow expression into every limit row) | yes |
| `compact` | only the (network, line) pairs screened on the design range, as PTDF rows | inside the design range |
| `red_btheta:<run>` | the reduced networks of a design run as B-θ blocks, all present | inside the design range |
| `red_btheta_lazy:<run>` | the same, an outage's block added once its limits are violated | inside the design range |
| `red_ptdf_lazy:<run>` | each reduced network's PTDF rows, added once violated | inside the design range |
| `...:<run>+x` | any of the three with exact limits on the critical lines (from the run's `masks_exact.csv`) | inside the design range |

## The decomposed PTDF (`baselines.jl`)
- Buses are split into areas by recursive spectral bisection (`partition`).
- For each area, the buses beyond its boundary are Kron-reduced away:
  $B^a = \begin{bmatrix} B_{II} & B_{IB} \\ B_{BI} & B_{BB} + A^a B_{EB}\end{bmatrix}$ with $A^a = -B_{BE}B_{EE}^{-1}$.
- The boundary buses carry equivalent injections $s^a = p_B + A^a p_E$; these are the consistency constraints.
- Each line belongs to the area of its from-bus. Its flow comes from that area's small PTDF over the internal injections and $s^a$.
- Post-outage flows use the system LODF.
- Every limit is added lazily.

The runner checks the construction: the decomposed flows equal $Hp$ to about $10^{-14}$.

## Running

```bash
julia --project=. --startup-file=no scopf_compare/run_compare.jl \
    'months=[7,6,8,9]' 'designs=["outputs/adversarial/ACTIVSg200_julbench/dom"]'
```

| Setting | Default | Meaning |
|---|---|---|
| `case` | `:ACTIVSg200` | test case |
| `months` | `[7, 6, 8, 9]` | whole months, solved in this order (the first is the design month) |
| `max_hours` | `nothing` | hours per month, for quick tests |
| `designs` | `[]` | design run folders (`parts/`, and `eval*/masks_exact.csv` for `exact:`); the first also gives the compact rows |
| `methods` | the exact methods, compact and the three reduced forms | subset to run (`:ptdf_lazy_flows` and `:ptdf_decomposed_flows` on request) |
| `areas` | 4 | areas of the decomposed PTDF, or a list such as `[4, 6]` |
| `time_limit` | 3600 | seconds per method and month; unsolved hours are reported |

On VACC: `sbatch hpc/compare_job.sh <out_dir> [settings]`.

Each method is built once and solves the months in order; the lazy methods keep
their rows or copies from one hour to the next. Every dispatch is judged against
the exact secure optimum, and by its worst loading on the full network over the
intact case and every outage.

## Output
- `forms.csv`: per method and month
  - `hours`, `unsolved_hours`
  - `vars`, `rows`, `nnz` (lazy methods: at the end)
  - `build_s`, `solver_s`, `wall_s`
  - `infeasible_hours`, `gap_mean_pct`, `gap_max_pct`
  - `overload_hours`, `worst_overload_pct`
- `hours.csv`: per method and hour — cost, gap and worst loading.
- `summary.txt`: the log.

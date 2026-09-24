# 6-bus Wood & Wollenberg case study

The 6-bus system from Wood & Wollenberg, with a year of hourly data attached.
`6busww.xlsx` is a capacity-expansion dataset: alongside the existing network it
carries candidate generators, candidate lines, candidate renewables and storage,
each with an investment cost.

## Files

| file | what it holds |
|---|---|
| `case6ww.m` | the original MATPOWER case |
| `6busww.xlsx` | source workbook (network, costs, 8760 h series) |
| `network_buses.csv` | bus, type, base kV, load factor, zone, scaling factor |
| `network_gens.csv` | the 3 **existing** generators (`IntI = 1`) |
| `network_lines.csv` | the 11 **existing** lines (`IntI = 1`) |
| `hourly_demand.csv` | gross bus demand, MW, 8760 h |
| `hourly_wind.csv` | wind output by bus, MW |
| `hourly_solar.csv` | solar output by bus, MW |
| `hourly_net_load.csv` | demand − wind − solar, the series the OPF uses |

Regenerate everything with `python analysis/extract_case_hourly.py 6busww` (standard
library only; xlsx is a zip of XML).

## How the workbook maps onto the network

**Load.** The `Demand_z1` sheet is a zonal series of 3228–8061 MW, far larger
than this system. The `Nodes` sheet carries the mapping:

    bus_demand[b,h] = Demand_z1[h] * ScalingFactor[b] * LoadFactor[b]

`LoadFactor` is 0.3 / 0.3 / 0.4 at buses 4, 5, 6 and sums to 1, so it is the bus
split; `ScalingFactor` is a uniform 0.08 and sets the system size. Buses 1–3
carry no load. Result: **258–645 MW** of gross demand.

**Renewables as load.** Wind and solar enter as *negative load* at their own
bus, not as generators:

    net_load[b,h] = demand[b,h] - wind[b,h] - solar[b,h]

Capacity comes from `Renewables.MaxCap` and the shape from the capacity-factor
sheet named by that row's **`T_index`** — 4 → `Wind(z1)`, 5 → `Solar(z3)`. Not by
`Zone`: on the 14-bus, renewable 6 sits at bus 3 with `T_index = 4` (wind) but
`Zone = 3`, so zone would mislabel it as solar. The two agree on this case.
Placement: 70 MW wind at bus 4; 50 MW wind + 80 MW solar at bus 5; 35 MW solar
at bus 6. Net load runs **75.7–592.7 MW**, mean 308.7 MW.

**Existing vs candidate.** Only `IntI = 1` rows are kept for generators and
lines, which reproduces `case6ww.m`'s topology exactly. The renewables are all
`IntI = 0` but are still used, because here they are load rather than
generation. Storage is ignored.

**Costs.** `Gen_Cost` (linear, $/MWh) is used, not the quadratic `gencost` in
the `.m`. The workbook's `info` sheet states the existing units' costs and
technical parameters were replaced with actual values.

## Planning data

The workbook is a generation-and-transmission expansion dataset. Every
candidate carries `Inv_Cost`, which is the **annualised** cost per year in
dollars (capital recovery plus fixed O&M), derived on the `Cost_Calculations`
sheet from overnight cost, life span and interest rate.

| candidate | where | size | `Inv_Cost` $/yr |
|---|---|---|---|
| thermal 4 (coal) | bus 1 | 100 MW | 20,480,046 |
| thermal 5 (CCGT) | bus 2 | 100 MW | 16,106,144 |
| thermal 6 (CCGT) | bus 6 | 50 MW | 8,053,072 |
| thermal 7 (SCGT) | bus 5 | 50 MW | 4,931,551 |
| line 12 | 1–5 | 80 MW, 50 km | 318,823 |
| line 13 | 4–5 | 70 MW, 30 km | 167,382 |
| line 14 | 5–6 | 100 MW, 70 km | 557,941 |
| wind 1 | bus 4 | 70 MW | 15,767,648 |
| wind 2 | bus 5 | 50 MW | 11,262,605 |
| solar 3 | bus 6 | 35 MW | 3,725,382 |
| solar 4 | bus 5 | 80 MW | 8,515,160 |
| battery 1 | bus 4 | 50 MWh / 25 MW | 3,586,844 |
| battery 2 | bus 5 | 34 MWh / 20 MW | 2,869,475 |

Technology index (`T_index`): 1 coal, 2 CCGT, 3 SCGT, 4 wind, 5 solar,
6 Li-ion storage, 7 HV AC line. Life spans are 30 y thermal, 25 y renewable,
15 y storage, 50 y line; interest 10% on everything except lines at 7%. Lines
are priced at \$1100/MW/km.

Operating cost is heat rate x fuel price + variable O&M, which reproduces
`Gen_Cost` exactly — e.g. unit 1: 8.0 x 2.89 + 7.33 = \$30.45/MWh.

Unit-commitment parameters are on the `Generators` sheet and repeated at the
foot of `Cost_Calculations`: `IntS` initial status (h), `IntP` initial output,
`MinU`/`MinD` minimum up/down time, `RU`/`RD` ramp rates, `SU`/`SD` start-up and
shut-down ramps, plus `Startup_Cost` in \$/start.

Note that the three candidate lines are **1–5, 4–5 and 5–6** — the first is the
line the full-year DC-OPF finds binding 92.6% of the time, and the third feeds
the bus whose import limit causes every shortfall.

## Running the year

    julia --project=. --startup-file=no analysis/run_6bus_year.jl pmin=relaxed
    julia --project=. --startup-file=no analysis/run_6bus_year.jl pmin=enforced

8760 DC-OPFs, about 3 s. Each hour carries unserved-energy and curtailment
variables priced above any generator, so the LP never reports INFEASIBLE and
the hours that would have been infeasible are measured by their shortfall
instead. Reports and per-hour CSVs land in `outputs/6busww/`.

## The reduction toy case

`analysis/build_toy_case.jl` turns this data into a case the three reduction
approaches can run on directly:

    julia --project=. --startup-file=no analysis/build_toy_case.jl case=6busww scale=0.018

It writes `6busww_toy.m` here and the hourly scenario matrices to
`outputs/case6ww_dcopf/` (`bus_ids`, `load_mw`, `generation_mw`,
`scenario_summary`). The case is registered as `:case6ww` in `common/cases.jl`
with `CASE_YEAR = 2015`, so `demands = :hourly` works.

Two things the generated `.m` does deliberately:

* **`Pmin` is written as 0.** It is a unit-commitment floor; every solver here
  runs `relax_pmin = true`, but `build_tx_case`'s base DC-OPF does *not* relax
  it, and `sum(Pmin) = 190 MW` against a 163 MW mean load makes the case fail to
  load at all. The real floor stays in `network_gens.csv`.
* **`load_mw` is the SERVED net load** (net load − shed + curtailment), so
  generation balances it exactly. `build_multiscenario_tx_case` rejects any
  imbalance, and the reduction MILP is infeasible for *every* clustering when
  `sum(p) != 0`.

### Running the approaches

The window is November (720 h). The pipeline's seed selection takes the 20
highest-flow hours as training; the other 700 are the test set.

    W="case=case6ww demands=hourly month=11 horizon_start_day=1 horizon_days=31 n_demands=20"
    julia --project=. --startup-file=no kkt/run_kkt.jl       $W cost_cap=0.1 time_limit=120
    julia --project=. --startup-file=no greedy/run_greedy.jl $W cost_cap=0.1 time_limit=120
    julia --project=. --startup-file=no proxy/run_proxy.jl   $W eps=0.05 time_limit=120 plots=false

### Reviewing a reduction out of sample

    julia --project=. --startup-file=no analysis/review_toy.jl case=6busww month=11 n_train=20
    julia --project=. --startup-file=no analysis/plot_toy_review.jl case=6busww
    python analysis/xlsx_write.py outputs/6busww_review/6busww_review.xlsx \
        Summary=outputs/6busww_review/summary.csv ...

`review_toy.jl` reads each approach's own clustering (`internal.csv` for kkt and
greedy, `line_status.csv` for the proxy), solves the reduced DC-OPF at every hour
of the month, puts that dispatch back on the full network, and reports
everything as a percentage — overload as % of the line's rating, cost as % over
the full network's optimum — split into the 20 training hours and the 700 test
hours.

`analysis/review_6bus.jl` is the separate year-long study that sweeps *candidate*
clusterings (rather than the approaches' answers) across all 8760 hours; it is
what produced the non-monotonicity results below.

### What the toy shows

At scaling 0.018, with a 0.1% cost cap and the November 8–14 window:

| clustering | buses | design viol. | out-of-sample viol. | worst overload | worst hourly cost gap |
|---|---|---|---|---|---|
| KKT optimum {2-4, 3-6} | 4 | 0 / 168 | **15 / 8592** | 0.99% (0.0040 p.u.) | 0.091% |
| {2-3, 2-5, 2-6} | **3** | 0 / 168 | **0 / 8592** | none | 4.12% |
| greedy | 6 | — | — | — | — |
| proxy | 6 | — | — | — | — |

Three properties, all reproducible in under a second:

1. **Merging more can *fix* feasibility.** `{3-6}` alone fails a design
   scenario; `{2-4, 3-6}` passes. The feasible family is not downward closed, so
   an infeasible set yields no valid cut.
2. **Merging more can *lower* the cost gap.** `{2-4}` alone is 0.158% over the
   full optimum — outside the cap — while `{2-4, 3-6}` is 0.091%, inside it.
3. **The cost cap costs reduction and robustness.** `{2-3, 2-5, 2-6}` reaches 3
   buses and never violates in the whole year, but its worst hourly cost gap is
   4.12%, so the 0.1% cap excludes it. The KKT's admissible answer is both
   *larger* (4 buses) and *less robust* (15 out-of-sample violations).

## Headline result

The existing network cannot serve this demand series. **3,633 h (41.5% of the
year) end with unserved energy, 161,384 MWh in total**, identical whether or not
`Pmin` is enforced — so it is a network limit, not a commitment artefact. Line
**3–6** is at its 80 MW limit in every one of those hours and in none of the
servable ones. Line **1–5** (40 MW) binds 92.6% of the year and is what bottles
up the cheapest unit at bus 1.

Delivering to one load bus at a time, with the others at zero, the network tops
out at roughly 139 MW (bus 4), 183 MW (bus 5) and 209 MW (bus 6); loaded
together it saturates near 350 MW against a 592.7 MW peak. The candidate lines,
generators and renewables in the workbook are what the dataset expects you to
build to close that gap.

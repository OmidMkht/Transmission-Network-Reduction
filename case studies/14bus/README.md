# 14-bus case study

The IEEE 14-bus system with a year of hourly data, packaged the same way as
[`../6busww`](../6busww/README.md). `14bus.xlsx` is a generation-and-transmission
expansion dataset: the existing network plus candidate generators, lines,
renewables and storage, each with an annualised investment cost.

Regenerate the CSVs with `python analysis/extract_case_hourly.py 14bus`.

## Load factor

The `Nodes` sheet's `Load Factor` column is the **bus share of total demand**,
and it reproduces the IEEE case14 load split exactly:

| bus | case14 Pd (MW) | Pd / 259 | workbook `Load Factor` |
|---|---|---|---|
| 2 | 21.7 | 0.083784 | 0.0837838 |
| 3 | 94.2 | 0.363707 | 0.363707 |
| 4 | 47.8 | 0.184556 | 0.184556 |
| 5 | 7.6 | 0.029344 | 0.0293436 |
| 9 | 29.5 | 0.113900 | 0.1139 |

The eleven load factors sum to **1.000000**, so:

    bus_demand[b,h] = Demand_z1[h] * ScalingFactor * LoadFactor[b]

with `ScalingFactor = 0.1563` uniform on the load buses. Buses 1, 7 and 8 carry
no load. That gives 774–1806 MW of gross demand against 770 MW of existing
generation, so the workbook's own scaling is far past what this network can
serve — the toy case below uses a lower scaling.

## Existing network vs planning candidates

Only `IntI = 1` rows are part of the network:

| | existing (`IntI = 1`) | planning-only (`IntI = 0`) |
|---|---|---|
| lines | **20** (matching IEEE case14) | **7** — lines 21–27 |
| generators | **5** (buses 1, 2, 3, 6, 8) | 6 |
| renewables | 0 | 6 (used anyway, as negative load) |
| storage | 0 | 4 (ignored) |

**The 7 planning-only lines are excluded from both the full and the reduced
network.** They are written to `network_planning.csv` for reference only:

| line | path | Fmax | `Inv_Cost` $/yr |
|---|---|---|---|
| 21 | 1–2 | 100 | 398,529 |
| 22 | 2–3 | 100 | 239,118 |
| 23 | 1–3 | 150 | 1,076,029 |
| 24 | 10–9 | 70 | 223,176 |
| 25 | 10–14 | 150 | 836,911 |
| 26 | 6–8 | 90 | 143,471 |
| 27 | 8–9 | 100 | 278,970 |

Several duplicate an existing corridor (21 duplicates line 1, 22 duplicates line
3), so treating them as part of the network would silently double those
capacities.

## Renewables

Wind and solar enter as **negative load** at their own bus. The capacity-factor
series is chosen by `T_index` (4 → `Wind(z1)`, 5 → `Solar(z3)`), **not** by
`Zone`: renewable 6 sits at bus 3 with `T_index = 4` but `Zone = 3`, so zone
would mislabel it as solar.

270 MW of wind (buses 6, 12, 7, 3) and 200 MW of solar (buses 9, 14).

## Toy reduction case

    julia --project=. --startup-file=no analysis/build_toy_case.jl case=14bus scale=0.05

Registered as `:case14toy`. At scaling 0.05 the net load runs −95 … 535 MW
against 770 MW of capacity, no hour is short, and the reduction problem has a
genuine interior — unlike the 6-bus, which jumps straight from "merge
everything" to "merge nothing".

| scaling | KKT result |
|---|---|
| 0.02 | 14 → 1 bus (everything merges) |
| 0.03 | 14 → 1 bus |
| 0.04 | 14 → 7 buses |
| **0.05** | **14 → 8 buses** |

### Running the approaches

    W="case=case14toy demands=hourly month=11 horizon_start_day=1 horizon_days=31 n_demands=20"
    julia --project=. --startup-file=no kkt/run_kkt.jl       $W cost_cap=0.1 time_limit=240
    julia --project=. --startup-file=no greedy/run_greedy.jl $W cost_cap=0.1 time_limit=240
    julia --project=. --startup-file=no proxy/run_proxy.jl   $W eps=0.05 time_limit=240 plots=false
    julia --project=. --startup-file=no analysis/review_toy.jl   case=14bus month=11 n_train=20
    julia --project=. --startup-file=no analysis/plot_toy_review.jl case=14bus

### Result (November, 20 training hours by flow, 700 test hours)

| approach | buses | reduction | test hours violating | worst test overload | worst test cost gap |
|---|---|---|---|---|---|
| Full network | 14 | — | 0% | — | — |
| **KKT** | 8 | **42.9%** | **5.43%** | 0.396% of rating | 0.076% |
| **Proxy** | 13 | 7.1% | **0%** | none | 0% |
| **Greedy** | 8 | **42.9%** | **3.00%** | 0.382% of rating | 0.085% |

All three pass every training hour. The differences appear only out of sample,
which is the point of the split: KKT merges one more line than greedy for the
same bus count and fails nearly twice as many test hours, while the proxy buys
its clean test record by barely reducing at all.

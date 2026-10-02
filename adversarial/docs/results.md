# Results

Every run of the per-contingency reduction so far, in order (28–29 September
2026). The method is in `method.md`; the reasons behind each choice are in
`decisions.md`.

## How to read the numbers

| Term | Meaning |
|---|---|
| network | the intact network (base) or the network after one line outage |
| critical pair | a (line, direction) that some dispatch in the design range can overload |
| redundant | a network with no critical pair; dropped from the reduced SC-DCOPF |
| buses / kept lines | size of a reduced network |
| cost gap | reduced SC-DCOPF cost minus the full SC-DCOPF (secure) optimum, in % of the optimum; above 0 means dearer |
| overload | worst loading of the returned dispatch on the full network, over the intact network and every outage, above 100% |
| inside the hull | the hour's demand is a convex combination of the design hours (where the guarantee holds) |
| design limits / exact limits | reduced ratings as designed (all kept lines at $F$) / exact limits on critical lines and none elsewhere |

Unless stated otherwise: 5% cost cap, no derating, exact witness, dominance cut.

## 1. Cases and data

| Case | Buses / lines / generators | Outages that keep it connected | Hourly data |
|---|---|---|---|
| 14-bus toy (`case14toy`) | 14 / 20 / 5 | 19 of 20 | November design (683 of 720 hours feasible under N-1), December test (705 of 744) |
| ACTIVSg200 | 200 / 245 / 38 | 173 of 245 | full year (2017, 8,760 h) |
| ACTIVSg2000 | 2,000 / 3,206 / 432 | 2,756 of 3,206 | January–March only (2016, 2,184 h) |

ACTIVSg200 demand is allocated to buses by zone, so the demand vectors span only
6 dimensions. Which months fall inside the July demand hull:

| Month | Jan–Apr | May | Jun | Jul | Aug | Sep | Oct | Nov–Dec |
|---|---|---|---|---|---|---|---|---|
| Hours inside the July hull | 0 | 146 | 220 | 744 | 320 | 194 | 72 | 0 |

No month has an hour above July's peak total demand (2,178 MW); months fall
outside the hull because their zonal mix differs, not only their total.

## 2. Profiles

**ACTIVSg200, March 1–7** (`analysis/sc_profile.jl`, `outputs/sc_profile/ACTIVSg200`):
- The full SC-DCOPF is feasible at all 168 hours with 79–81 post-outage rows.
- Security costs 9.3% over the plain DC-OPF on average (10.8% at most).
- At a 5% cap, 32 of 173 outage networks are redundant; the others have 1–2
  critical lines, almost always line 184.

**ACTIVSg2000, March 1–7:**
- The full SC-DCOPF is feasible at all 168 hours with 3,943 post-outage rows
  (59 s per 168 hours on VACC).
- The intact network has 398 critical pairs on 289 lines. 391 of them can be
  overloaded by more than 1%, 324 by more than 10%.

## 3. The 14-bus toy case

| Variant | Mean buses (of 14) | Reduction | Smallest | Cost gap | December: infeasible hours |
|---|---|---|---|---|---|
| Dominance, exact witness | 12.1 | 13.6% | 9 | 0 | 0 of 705 |
| Dominance, 2% or 5% witness margin, no derating | 12.1 | 13.6% | 9 | 0 | 0 |
| Dominance, 2% derating | 10.75 | 23.2% | 5 | ≤1.61% | 5 |
| … with exact limits | 10.75 | 23.2% | 5 | ≤0.98% | 5 |

Other cuts (`outputs/adversarial/case14toy_line`, `case14toy_any_base`):
- **Weaker cut** (`:line`): same reduction (12.1 buses), 248 s of design against
  13 s, 9–27 adversaries per network.
- **Any-line cut** (`:any`): did not converge in 10 minutes on the base network.

**Network form (planning lines).** Each candidate line is added to the full
network and to every reduced block:
- With the exact-witness designs the reduced SC-DCOPF stays exact for lines 21–25.
- Lines 26 and 27 attach to bus 8, which every design merged because it was
  radial, and break the designs (up to 9.6% overload).
- Compact rows computed beforehand are wrong for every candidate (cost errors up
  to 16.6%, overloads up to 8%).

## 4. ACTIVSg200, March 1–7 design, first master

The first master minimised buses (root and commodity-flow variables) and started
from "merge nothing". Outage networks started from the base design, one VACC job
each (`outputs/adversarial/vacc/ACTIVSg200_{d0,d2,d5,mix}`).

| | d0 (your rule) | d2 | d5 | Combined (d2, d0 for 3) |
|---|---|---|---|---|
| Rule on the critical line | reduced ≥ full | reduced ≥ 98% of full | reduced ≥ 95% of full | – |
| Witness | exact optimum | optimum at 98% ratings | optimum at 95% ratings | exact optimum |
| Certified | 142 of 142 | 142 of 142 | 142 of 142 | 142 of 142 |
| Buses: median / mean (of 200) | 3 / 6.6 | 3 / 4.2 | 6 / 7.6 | 3 / 5.0 |
| Networks with ≤4 buses | 134 | 139 | 16 | 139 |
| Larger networks (buses) | 128, 128, 74, 56, 47, 45, 18, 6 | 91, 51, 22 | 95, 56, 49, 43, 16, 16, 12 | 128, 128, 22 |
| Ratings lowered at the end | none | 7, ≤1.93% | 100, ≤4.83% | – |
| Certified straight from the base design | 121 | 129 | 117 | – |
| Design time: total / slowest | 3.5 h / 20 min | 2.3 h / 18 min | 4.3 h / 18 min | – |

The combination takes d2's design everywhere except outages 27, 73 and 242,
whose d2 designs reject the exact optimum.

**SC-DCOPF forms** (reduced forms with exact limits; test week July 17–23):

| Form | Variables | Rows | Solve, 168 h | Design week: gap / overload | July: infeasible h / worst gap / overload |
|---|---|---|---|---|---|
| Full network copies | 77,295 | 77,432 | 9–10 s | exact / 0 | 0 / exact / 0 |
| Lazy PTDF rows | 38 | 653 | 0.8 s (wall, with the checks) | exact / 0 | 0 / exact / 0 |
| Compact PTDF rows | 38 | 299 | 0.02–0.03 s | exact / 0 | 0 / ≈0 / 0.004% |
| Reduced, d0 | 3,551 | 3,656 | 0.39 s | exact / 0 | 32 / 4.6% / 0 |
| Reduced, d2 | 2,286 | 2,391 | 0.24 s | ≤0.24% / 0 | 14 / 3.2% / 0 |
| Reduced, d5 (own ratings) | 3,665 | 3,770 | 0.56 s | ≤0.55% / 0 | 91 / 0.9% / 0 |
| Reduced, combined | 2,540 | 2,645 | 0.31 s | exact / 0 | 14 / 3.2% / 0 |

**Exact limits.** In d2 only 1.3 of the 11.6 kept lines per network keep a
limit. The critical lines' exact limits have a median of 1.053 F (range
0.967–1.69 F), above F for 130 of 139: the reduced networks overstate line 184's
flow, so they can safely accept more.

## 5. ACTIVSg200, March 1–7 design, master v2

Master v2 (`method.md`, Section 7): maximise merged lines, critical lines fixed
kept, bridges fixed merged, outage start with local unmerge and a neighbourhood,
MIPFocus 1 and a 5% gap. VACC run `v2` (`outputs/adversarial/vacc/ACTIVSg200_v2`).

**Designs:**
- All 142 non-redundant networks certified (32 redundant); 112 straight from the base design.
- Buses: median 3, mean 5.3, 127 networks with ≤4 buses. Larger ones: 128
  (outages 73 and 242), 21 (195), 13, 12, 8.
- 13.1 kept lines per network on average.
- Design time 59 minutes summed over networks (first master, d0: 3.5 hours);
  slowest network 9 minutes.

**Against the first master** (local runs on the same settings): 3–7 times faster
on the hard networks, with leaner designs; 31 minutes in total locally
(`outputs/adversarial/ACTIVSg200_v2_b5`) against 3.5 hours.

**Without a cost cap** (local run, stopped part-way, `outputs/adversarial/ACTIVSg200_v2`):
about 15 times slower, with more looping networks, more rounds, 47 s master
solves of which 76% hit the time limit. The base network came out at 33 buses
instead of 3.

**SC-DCOPF forms**, July 17–23 as the test week (`eval_fullbase_false`, `eval_fullbase_true`):

| Reduced SC-DCOPF | Variables | Design week | July: infeasible h | July: worst overload (hours) |
|---|---|---|---|---|
| Reduced base, design limits | 2,644 | exact, no overload | 0 | 0.23% (2) |
| Reduced base, exact limits | 2,644 | exact, no overload | 0 | 1.66% (4) |
| Full base, design limits | 3,076 | exact, no overload | 0 | 1.32% (4) |
| Full base, exact limits | 3,076 | exact, no overload | 0 | 1.66% (1) |

Every overload is on line 185 after outage 195. Keeping the intact network
whole changes the cost at no hour. The different overload figures come from the
solver picking a different dispatch at the same cost.

### Why July overloads: diagnosis

The overloads come from July lying outside the design range, not from the
reduction method:

- **All 168 July hours lie outside the March-week hull**, by 3.6–17.9% of total
  demand (L1). Total demand is 14.0–17.9 p.u. in the design week and 12.7–21.8
  p.u. in the July week.
- **Every overloaded hour is a July peak hour**, with a secure cost 27% above
  March's highest.
- **Line 185 after outage 195 was not critical in March.** No dispatch in the
  March range can push it past 0.95 F, so no cut protects it. Outage 195's design
  (21 buses, 38 kept lines) keeps it with limit F, but the merges elsewhere shift
  flow. At hour 4767 the reduced flow is −0.975 F against −1.002 F on the full
  network: the reduced network accepts, the full one overloads. With exact
  limits the line has no limit at all, hence the larger overload (hour 4744:
  −1.0049 F reduced, −1.0166 F full).
- **The problem is wider than line 185.** Screening March and July together gives:
  - 168 critical pairs instead of 150;
  - 18 new pairs, including line 185 after outage 195, which now reaches 1.25 F;
  - 4 outage networks that March dropped as redundant (76, 175, 176, 182) now critical;
  - 106 of 146 designs fail the check over that range (132 with exact limits), by up to 25%.

  The July overloads were small only because an optimal dispatch is not a worst case.
- **Compact PTDF rows fail the same way:** one July hour overloads by 0.004%,
  because those rows were also chosen by March screening.

### Optimality, March-week design

Signed cost gap per hour, 168 hours each:

| Form | Design week | July: mean | July: median | July: max | July hours above the optimum | July hours overloaded |
|---|---|---|---|---|---|---|
| Compact PTDF rows | 0 | 0 | 0 | 0 | 0 | 1 (0.004%) |
| Reduced, design limits | 0 | +0.30% | +0.13% | +1.40% | 104 | 2 (4 with full base) |
| Reduced, exact limits | 0 | +0.09% | 0 | +0.87% | 49 | 4 (1 with full base) |

No hour is cheaper than the secure optimum. By July load level:

| July hours | Count | Design limits: mean / max | Exact limits: mean / max |
|---|---|---|---|
| Secure cost above March's peak | 88 | +0.56% / +1.40% | +0.17% / +0.87% |
| Within March's cost range | 64 | +0.015% / +0.26% | 0 / 0 |
| Below March's lowest | 16 | 0 | 0 |

In July's peak hours the reduced networks are wrong both ways at once: too strict
elsewhere (cost above the optimum) and too loose on line 185 (overload).

## 6. ACTIVSg200, July design (all of July)

VACC run `jul` (`outputs/adversarial/vacc/ACTIVSg200_jul`): design on July 1–31
(744 hours), tested on August 1–31 (`eval/`) and March 1–31 (`eval_mar/`), both
full months.

**Designs:**
- All 145 non-redundant networks certified (29 redundant); 102 straight from the base design.
- Buses: median 4, mean 9.3, 114 networks with ≤4 buses. Larger ones: 128 (outages
  242, 195, 193), 127 (73), 123 (7), 46 (53).
- 18.1 kept lines per network on average. The base network: 4 buses, 11 lines, 39 s.
- Design time 2.7 hours summed; slowest network 21 minutes; 19 networks over a minute.
- Exact limits: 136 designs changed, 9 kept their own (6 with a merged critical
  line, 3 failing the check). All 145 accept the exact optimum at every hour.

**SC-DCOPF forms** (744 hours per month):

| Form | Variables | Rows | Solve, 744 h |
|---|---|---|---|
| Full network copies | 77,295 | 77,432 | 55–63 s |
| Compact PTDF rows | 38 | 337 | 0.2 s |
| Reduced | 4,011 | 4,119 | 2.6–3.5 s |

| Month | Hours inside July's hull | Form | Infeasible h | Worst overload | Cost gap: mean / max | Hours above 0.001% |
|---|---|---|---|---|---|---|
| July (design) | 744 | reduced, both limits | 0 | 0 | 0 / 0.00001% | 0 |
| August | 321 (rest up to 4.4% away) | reduced, design limits | 0 | 0 | 0.0002% / 0.036% | 16 |
| August | | reduced, exact limits | 0 | 0 | 0.0001% / 0.007% | 10 |
| March | 0 (up to 9.4% away) | reduced, both limits | 0 | 0 | 0.00001% / 0.00001% | 0 |
| August, March | | compact | 0 | 0 | 0 | 0 |

- Inside July's hull the cost gap is at most 0.00001%, i.e. exact. Every August
  hour with a gap lies outside the hull.
- The March-week design had given 0.2–1.7% overloads and up to 1.4% extra cost in
  July. Designing on the peak month removed both, even in March, which lies
  entirely outside July's hull. That is observed, not guaranteed: July's range
  covers the flows that matter in March.

## 7. Screening statistics

Measured with a scratch script (`screen_stats.jl`, not in the repo) that repeats the runner's screening with counters.

**ACTIVSg200, July design, all 173 outage networks:**

| Step | (line, direction) pairs |
|---|---|
| All outages × lines × 2 | 84,424 |
| Cleared by the LODF bound | 84,224 (99.76%) |
| LPs solved | 200 |
| … already critical at a stored dispatch (the intact-extreme LP optima and the hourly optima) | 167 |
| … truly undecided | 33 |
| Critical | 168, in 144 networks |

- Critical lines per network: median 1, maximum 5.
- Worst loading of the critical pairs: 144 above 0.1%, 126 above 0.5%, 112 above
  1%, 93 above 2%, 60 above 5%, 40 above 10%.
- Exact duplicates (pairs that load identically at every injection): none.

**ACTIVSg2000, intact network:** 6,412 LPs in 11 s give 398 critical pairs on 289
lines. 5 of them duplicate another. Dropping pairs whose worst overload is at most
ε would remove only 7 at ε = 1% and 74 at ε = 10%.

## 8. ACTIVSg2000

**Local smoke test** (1 design day, 30 s master solves): 370 critical pairs; the
second master solve hit its limit and certified at 1,547 of 2,000 buses.

**VACC run `s1`** (`outputs/adversarial/vacc/ACTIVSg2000_s1`): design March 1–7,
test February 1–7, 600 s master solves, 8 CPUs, 32 GB. The intact network plus a
sample of 20 outage networks (array indices 2, 52, …, 952), no evaluation.

The intact network:
- 289 critical lines fixed kept, 450 bridges fixed merged, 2,467 lines left to decide.

| Round | Design | Master | Adversaries left | Worst loading |
|---|---|---|---|---|
| 1 | 103 buses, 201 lines | 0 s | 76 | 2.04 |
| 2 | 1,505 buses, 2,703 lines | 600 s (time limit) | 79 | 1.12 |
| 3 | 1,550 buses, 2,756 lines | 516 s (within 5% gap) | 0: certified | 1.00 |

- 22 minutes in total, but only the 450 bridges were merged (22.5% of buses
  removed).
- The last master reports that, with its 6 adversaries, at most about 22 more
  merges could exist.
- Cause: with the dominance cut and no derating, each cut is met exactly by the
  unreduced network, and almost every merge lowers one of the 6 cut flows
  slightly (`decisions.md`, "The weaker cut").

The outage sample:
- 19 networks designed, all certified straight from the base design at 1,543–1,550
  buses, with 17–89 critical pairs each. One (network 792) is redundant.
- Each design took 22–75 s, but each job took 14–25 minutes, about 1,000 s of it
  screening: every job re-solves the 6,412 intact-extreme LPs.

**Hop locality of the intact network** (how far each undecided line is from the
nearest critical line, and how much merging that line alone can shift a critical
line's flow, at the merged line's rating):

| Hops | Undecided lines | Cumulative share | Shift of a critical flow (% of its rating): median / 90th pct / max | Lines above 1% |
|---|---|---|---|---|
| 1 | 538 | 21.8% | 26.2 / 76.1 / 206 | 529 |
| 2 | 503 | 42.2% | 10.3 / 40.1 / 176 | 484 |
| 3 | 391 | 58.0% | 3.7 / 16.3 / 80 | 321 |
| 4 | 375 | 73.2% | 1.1 / 7.0 / 134 | 201 |
| 5 | 277 | 84.5% | 0.39 / 2.2 / 29 | 67 |
| 6 | 156 | 90.8% | 0.15 / 1.5 / 7.4 | 25 |
| 7 | 125 | 95.9% | 0.12 / 0.72 / 3.6 | 7 |
| 8–11 | 102 | 100% | ≤0.09 / ≤0.39 / 4.4 | 3 |

Merging line $j = (a, b)$ alone shifts line $\ell$'s flow by
$-f_j\, \mathrm{PTDF}_\ell(a \to b) / \mathrm{PTDF}_j(a \to b)$ (rank-one update
with infinite susceptance on $j$).

## 9. Runs in progress (29 September 2026)

| Run | What | Status when written |
|---|---|---|
| `ACTIVSg200_jul_line` (jobs 5508740–43) | ACTIVSg200 July design with the weaker cut, evaluated on August and March | finished, below |

**`ACTIVSg200_jul_line`, finished:** the July design with the weaker cut, against
the dominance cut (`jul`).

| | Dominance (`jul`) | Weaker cut (`jul_line`) |
|---|---|---|
| Buses: median / mean / max | 4 / 9.3 / 128 | 3 / 9.7 / 128 |
| Networks with ≤4 buses | 114 | 116 |
| Kept lines per network | 18.1 | 19.3 |
| Certified straight from the base design | 102 | 103 |
| Adversaries added in total | 333 | 755 |
| Design time: total / slowest | 2.7 h / 21 min | 6.8 h / 43 min |
| Reduced SC-DCOPF variables | 4,011 | 4,247 |
| August: worst cost gap, design / exact limits | 0.036% / 0.007% | 0.046% / 0.0000% |
| March: worst cost gap | 0.00001% | 0.0000% |
| Overload, July / August / March | 0 | 0 |

Same reduction and safety; the weaker cut is 2.5 times slower.

**`ACTIVSg2000_line`, finished** (job 5508700): the intact network with the weaker
cut, otherwise as `s1`.

| Round | Design | Master | Witness hours rejected | Adversaries left | Worst loading |
|---|---|---|---|---|---|
| 1 | 103 buses, 201 lines | 0 s | 0 | 76 | 2.04 |
| 2 | 1,046 buses, 2,062 lines | 600 s (time limit) | 121 | 60 | 1.48 |
| 3 | 1,550 buses, 2,756 lines | 51 s (within 5% gap) | 0 | 0: certified | 1.00 |

It ends with the same design as the dominance cut (only the bridges merged),
after 6 adversaries and 1 witness hour, in 14 minutes. The weaker cut did not
give the master room on this network.

**Master bound, rerun with Gurobi logs saved** (`ACTIVSg2000_line_bound`,
`ACTIVSg2000_s1_bound`; the first runs logged no bound). The objective counts
every merged line, including the 450 fixed bridges.

| Run | Round | Incumbent (merged lines) | Best bound | Gap | Status |
|---|---|---|---|---|---|
| weaker cut | 1 | 2,917 | 2,917 | 0% | optimal |
| weaker cut | 2 | 1,319 (925 buses) | 2,915 | 121% | time limit (600 s) |
| weaker cut | 3 | 450 (1,550 buses) | 450 | 0% | optimal (195 s) |
| dominance | 1 | 2,917 | 2,917 | 0% | optimal |
| dominance | 2 | 503 (1,505 buses) | 2,915 | 480% | time limit (600 s) |

In round 3 of the weaker-cut run the root bound fell from 2,915 to 450 in one
step through cuts, and Gurobi warned "max constraint violation (1.78e-08) exceeds
tolerance" (tolerance 1e-9). The dominance run's root LP passed through an
objective of 1.9e31 and a primal infeasibility of 9.5e36. Both are signs of
numerical trouble, so the "no free line can be merged" proof is not yet trusted.

**The proof is wrong.** A rerun that saves every injection the master held
(`ACTIVSg2000_line_copies`, which repeated the run above exactly: 1,046 buses in
round 2, then 450 merges with bound 450 in 51 s) was checked without the MIP. Merging
one line $j=(a,b)$ changes line $\ell$'s flow by exactly
$-f_j\,T_{\ell j}/T_{jj}$ with $T_{:,j}=H_{:,a}-H_{:,b}$ (confirmed against the
contracted network's PTDF to 10 digits).

| Master (round) | Copies held | Slack of "bridges only" | Single extra merges that satisfy every copy |
|---|---|---|---|
| weaker cut, round 3, "optimal, gap 0%" | 6 adversaries, 1 witness hour | 23.8% of rating | 2,460 of 2,467 (many with about 25% slack) |
| dominance, `ACTIVSg2000_s1_copies`, round 2, time limit with no improvement | 3 adversaries | 0 (exactly tight) | 294 of 2,467 |

Any of those 2,460 designs has 451 merges, above Gurobi's "best bound" of 450, so
the bound is invalid: a numerical error in the master. The dominance master
did not claim optimality; it simply found nothing better than its start in
600 s although single merges were available. Safety is not affected (the check
is exact); only the reduction is.

**Far-merge check** (scratch script `far_merge_check.jl`): merge the bridges plus
every undecided line at least 7 hops from any critical line (227 lines, 1,366
buses left) and run the exact check.

| Limits on the reduced network | Witness hours rejected | Adversaries | Worst loading |
|---|---|---|---|
| F on every kept line (as in the design) | 168 of 168 | 52 | 1.26 |
| F on critical lines only | 0 of 168 | 63 | 1.37 |

- **The limits on non-critical kept lines alone make the reduced network reject
  the secure optimum at every hour.** Without them, all witnesses pass.
- **Far merges together do hide real overloads on critical lines.** 52–63
  adversaries remain even though no merged line is within 6 hops of a critical
  line. Individually small shifts (Section 8) add up when whole regions are
  contracted.
- At the secure optimum no intact line sits at its rating, on ACTIVSg200 (744
  July hours) or ACTIVSg2000 (168 March hours).

**Reruns with the master fixes and polishing** (intact network only, otherwise as
`s1`; 600 s master, 300 s polishing per round, the 3 worst adversaries added per
round; VACC `ACTIVSg2000_s1_polish`, `ACTIVSg2000_line_polish`):

| Run | Cut | Rounds | Result |
|---|---|---|---|
| `s1_polish` (job 5510081) | dominance | 8 | certified at 1,543 buses, 2,741 lines, in 1.3 h (21 adversaries, 2 witness hours) |
| `line_polish` (job 5510080) | weaker | 107 after 19 h, still going | not certified; designs of 103–291 buses, each failed by the next check (54–103 violations, worst loading 1.55–2.04) |

- **Dominance:** only 7 merges beyond the bridges. In round 7, polishing
  from 1,542 buses found no single merge that satisfies all 21 adversaries.
- **Weaker cut:**
  - From round 4 the master never improves on its start. By round 107 it has
    1.76 M rows and 2.77 M columns, and its bound is still 2,917 (merge
    everything). Every design comes from polishing alone.
  - Polishing merges as much as the ~320 adversaries held allow, and the next
    check always finds new overloads. The worst loading has stayed at 1.55–1.9
    since round 30, with no sign of converging.
- 103 buses is the floor: all 289 critical lines kept, every other line merged.
- `line_polish` ended at its 20 h budget after 113 rounds (336 adversaries),
  uncertified, so it returned the unreduced network.

**Greedy on the intact network** (`search = :greedy`, 1 October 2026; VACC
`ACTIVSg2000_greedy_sc`, `ACTIVSg2000_greedy_dc`; design March 1–7, test February
1–7):

| Run | Range | Critical pairs | Result | Single merges that failed (adversary / witness) | Time |
|---|---|---|---|---|---|
| `greedy_sc` (job 5539125) | 5% over the SC-DCOPF optimum, as `s1` | 398 on 289 lines | 1,540 buses, 2,738 lines | 4,896 / 12 over 2 sweeps | 17 min |
| `greedy_dc` (job 5539126) | 5% over the plain DCOPF optimum (`contingencies=false`) | 372 on 276 lines | 1,540 buses, 2,739 lines | 1,585 / 3,348 over 2 sweeps | 27 min |

- Both certified 10 merges beyond the bridges, all within the first 25 lines
  tried. After that, no single remaining merge passes the exact check.
- This is the exact condition, so the stall is the condition, not the search. It
  matches the dominance master (1,543 buses).
- Any merge shifts critical flows by an amount whose sign depends on the
  dispatch. A dispatch just above a critical line's rating, shifted just below
  it, is then accepted by the reduced network: an adversary, however small the
  shift.
- On the plain DCOPF range the optimum sits at line ratings, so most merges also
  reject the witness.
- `greedy_dc` evaluation (B-θ DCOPF):

  | Form | Variables | Solve, 168 h | Cost gap | Overload |
  |---|---|---|---|---|
  | full | 5,638 | 3.0 s | — | 0 |
  | reduced | 4,710 | 2.2 s | ≤ 7e-8 % | 0 |

  - The cost gap and overload hold on both weeks.
  - Every February hour lies outside the March hull (up to 18.9% of demand
    away).
  - The full B-θ form left 2 test hours unsolved.
- The faster check (lazy rows, stop at the first failure) takes about 0.2–0.3 s
  per failed single merge at 2,000 buses, against about 30 s per check in the
  polish runs.

On ACTIVSg200 (intact network, March week) the same greedy reaches 5 buses in 3 s
on the SC range. On the plain DCOPF range it stops at the bridges (128 buses).

## 10. July benchmark of the master approaches (ACTIVSg200)

Fresh designs on all of July (744 h), 4 h per network, master 60 s, exact
witness, 5% cap (VACC campaign `ACTIVSg200_julbench`, one folder per run;
`analysis/adversarial_benchmark.py`). Design statistics:

| Run | Master | Certified (of 145) | Buses: median / mean / max | Adversaries | Design time: total / slowest |
|---|---|---|---|---|---|
| `dom` | dominance cut | 145 | 4 / 8.9 / 128 | 329 | 2.6 h / 15 min |
| `dom_polish` | + polishing | 145 | 4 / 7.9 / 128 | 347 | 2.8 h / 22 min |
| `line` | weaker cut | 145 | 3 / 10.4 / 128 | 796 | 7.4 h / 50 min |
| `line_polish` | weaker cut + polishing | 144 (1 hit 4 h) | 3 / 7.0 / 200 | 1,611 | 19.2 h / 4 h |
| `any_crit` | any critical line rejects | 145 | 4 / 10.2 / 128 | 778 | 7.1 h / 56 min |
| `any_all` | any kept line rejects | 143 (2 hit 4 h) | 5 / 19.9 / 200 | 916 | 8.9 h / 45 min |
| `tol1` | dominance, 1% tolerance | 93 (fewer critical networks) | 3 / 8.7 / 112 | 302 | 2.2 h / 16 min |
| `critlim` | critical-only ratings | 145 | 5 / 9.3 / 128 | 399 | 3.0 h / 22 min |
| `hops` | hop-limited outage search | 145 | 3 / 8.0 / 128 | 306 | 1.9 h / 13 min |
| `combo` | polishing + critical-only + hops | 145 | 5 / 8.2 / 128 | 406 | 2.7 h / 17 min |

Kept from this benchmark: dominance cut + hop-limited outage search (`hops`),
with polishing as an optional add-on (`decisions.md`).

### Cluster hop ladder (1–2 October)

Same setup as `hops` (polishing off), with the ladder `[3, 6, free]`, hold-forward.
`hops_v2` reruns `hops` on the current code. `hops_ladder_out` copies the intact
network's design from `hops_v2`, so only the outage networks differ.

| Run | Ladder on | Buses: median / mean | ≤4-bus networks | Kept lines | Variables | Design time: total / slowest | Worst test cost gap: design / exact limits |
|---|---|---|---|---|---|---|---|
| `hops` | – | 3 / 8.0 | 121 | 17.1 | 3,685 | 1.9 h / 12.5 min | 0.0012% / 0.00001% |
| `hops_v2` | – | 3 / 7.7 | 122 | 17.5 | 3,701 | 2.5 h / 16 min | 0.12% / 0 |
| `hops_ladder_out` | seedless outage networks | 3 / **7.3** | **124** | **17.0** | **3,562** | 2.7 h / 21 min | 0.12% / 0 |
| `hops_ladder` | intact network too | 4 / 8.4 | 112 | 22.3 | 4,485 | 2.9 h / 15 min | 0.37% / 0.35% |

The variables are the reduced SC-DCOPF's. All 145 networks were certified in
every run. The design month (July) is exact in every run, with no overload or
infeasible hour in any month.

- **Intact network:** the ladder gave 4 buses and 17 lines in 145 s; the free
  master gave 3 buses and 11 lines in 12 s.
  - The rungs went 31 → 23 → 4 buses, and the merges held from the capped rungs
    cost one bus and six lines.
  - Every outage network starts from that design, hence `hops_ladder`'s larger
    networks (outage 6: 28 buses against 5) and cost gaps outside the hull.
- **Outage networks:** the 10 that lost their seed (whose hop region grew to the
  whole network) went from 666 to 620 buses in total, and from 100 to 116
  minutes.

  | Outage | Without the ladder | With it |
  |---|---|---|
  | 7 | 13 | 6 |
  | 53 | 49 | 6 |
  | 91 | 10 | 5 |
  | 122 | 10 | 9 |
  | 86 | 63 | 72 |
  | 140 | 9 | 10 |
  | 73, 193, 195, 242 | 128 | 128 |

  128 buses means only the bridges are merged; those four stay there in every run.
- **The differences are about as large as the noise between runs.** The two runs
  without the ladder (`hops` and `hops_v2`) differ by 0.3 buses on the mean, by
  0.6 h in design time, and by 0.0012% against 0.12% in test cost gap with design
  limits, from slightly different intact designs.
- **Default since:** the ladder only for outage networks (`ladder_base = false`).

## 11. Reduced networks against exact SC-DCOPF methods

`scopf_compare/` (see its README), run from scratch
(`outputs/scopf_compare/ACTIVSg200_fair_r1`, `_r2`: the same job twice, times
averaged; `analysis/scopf_compare_table.py`). The July designs of Section 10
against the full network, over July (design month), June, August and September:
2,928 hours. The reduced networks are solved in every formulation the full
network is, so each pair is like for like:

- **B-θ, every block:** a B-θ block per network.
- **B-θ, lazy blocks:** an outage's block is added once violated.
- **PTDF, lazy rows, substituted:** each limit row writes the whole flow out over the generators.
- **PTDF, lazy rows, flow variables:** the paper's form. A flow variable per line with one defining row. The reduced version bounds each reduced network's flow variable.

The forms of one design return the same costs (within $10^{-9}$). The two
replicates differ by up to about 15-20% in wall time on these short totals, so
smaller differences are noise.

**Full vs reduced, formulation by formulation** (per hour; reduced: `hops`,
exact limits):

| Formulation | Full network: size, solver / wall | Reduced networks: size, solver / wall | Speed-up (wall) |
|---|---|---|---|
| B-θ, every block | 77,295 vars; 70.3 / 98.9 ms | 3,685 vars; 3.3 / 6.7 ms | 15× |
| B-θ, lazy blocks | 36,447 vars; 31.7 / 40.9 ms | 1,675 vars; 1.4 / 1.8 ms | 23× |
| PTDF, lazy rows, substituted | 38 vars, 351 rows; 0.22 / 0.68 ms | 38 vars, 293 rows; 0.20 / 0.53 ms | 1.3× (solver 1.1×) |
| PTDF, lazy rows, flow variables | 211 vars, 524 rows; 0.21 / 0.50 ms | 184 vars, 147 rows; 0.12 / 0.35 ms | 1.4× (solver 1.75×) |

Other exact methods on the full network, and compact rows (screened on the
design range, not a network), per hour wall: decomposed PTDF substituted 1.00 ms
(4 areas) / 0.84 ms (6); decomposed PTDF with flow variables 0.57 / 0.60 ms;
compact rows 0.40 ms, with flow variables 0.33 ms.

**Every design run** (per hour wall, exact limits):

| Run | B-θ, every block | B-θ, lazy blocks | PTDF, substituted | PTDF, flow variables | Worst cost gap, June/Aug/Sep (design / exact limits) | Infeasible / overloaded hours |
|---|---|---|---|---|---|---|
| `hops` | 6.7 ms | 1.8 ms | 0.53 ms | 0.35 ms | 0.0012% / 0.00001% | 0 / 0 |
| `dom_polish` | 6.3 ms | 2.1 ms | 0.54 ms | 0.35 ms | 0.016% / 0.00001% | 0 / 0 |
| `dom` | 6.8 ms | 2.0 ms | 0.54 ms | 0.35 ms | 0.035% / 0.00001% | 0 / 0 |
| `critlim` | 6.7 ms | 2.1 ms | 0.52 ms | 0.33 ms | 0.011% / 0.0072% | 0 / 0 |
| `line` | 7.2 ms | 2.0 ms | 0.52 ms | 0.35 ms | 0.025% / 0.025% | 0 / 0 |
| `line_polish` | 6.0 ms | 2.0 ms | 0.52 ms | 0.36 ms | 0.076% / 0.035% | 0 / 0 |
| `any_crit` | 7.0 ms | 2.1 ms | 0.51 ms | 0.34 ms | 0.11% / 0.0028% | 0 / 0 |
| `any_all` | 11.7 ms | 6.1 ms | 0.53 ms | 0.36 ms | 0.0082% / exact | 0 / 0 |
| `combo` | 6.2 ms | 2.1 ms | 0.51 ms | 0.33 ms | 0.55% / 0.15% | 0 / 0 |
| `tol1` | 4.3 ms | 1.3 ms | 0.38 ms | 0.27 ms | 0.069% / 0.025% | 0 / every hour, ≤0.73% / ≤1.00% (its tolerance) |
| *full network* | *98.9 ms* | *40.9 ms* | *0.68 ms* | *0.50 ms* | *exact* | *0 / 0* |

In July every run is exact to 0.00001% (`tol1`: cheaper by design). With design
ratings the lazy forms are a little faster still (`hops`: 1.4 ms lazy B-θ,
0.44 ms PTDF with flow variables).

- Like for like, the reduced networks are faster in every formulation: 15-23
  times in the B-θ forms, 1.3-1.4 times in the PTDF forms.
- The fastest exact method on the full network is PTDF with flow variables
  (0.50 ms). The reduced networks in that form (0.35 ms) tie with compact rows
  in the same form (0.33 ms), and beat the decomposed PTDF, which does not pay
  off at 200 buses.
- An earlier, less careful comparison (`ACTIVSg200_julbench9b`) used the
  substituted form only, with one job; its PTDF figures are superseded by these.

ACTIVSg2000, March 1-7 (168 h): lazy B-θ copies solved 1 hour in 4.5 h (6.95 M
variables, about 1,330 outage copies).

**ACTIVSg2000, January (the highest-load month: mean 34.2 GW, peak 45.8 GW),
744 hours, exact methods only** (`outputs/scopf_compare/ACTIVSg2000_jan`). Both
PTDF methods in both forms: *substituted* (each limit row holds the whole flow
expression, over all generators) and *flow variables* (the paper's (4) and (9):
a flow variable per needed line with one defining row, and limit rows of one or
two flows).

| Method | Variables | Rows | Nonzeros | Solver | Wall |
|---|---|---|---|---|---|
| PTDF, lazy rows, substituted | 432 | 49,819 | 21.4 M | 543 s | 646 s |
| Decomposed PTDF, 4 areas, substituted | 557 | 49,944 | 11.9 M | 373 s | 472 s |
| Decomposed PTDF, 6 areas, substituted | 603 | 49,990 | 10.1 M | 330 s | 422 s |
| PTDF, lazy rows, flow variables | 3,188 | 52,575 | 1.28 M | 42 s | 66 s |
| Decomposed PTDF, 4 areas, flow variables (the paper) | 3,313 | 52,700 | 0.51 M | 25 s | 48 s |
| Decomposed PTDF, 6 areas, flow variables | 3,359 | 52,746 | 0.46 M | 24 s | 46 s |

- All six return the secure optimum at every hour (gaps below $10^{-12}$, no
  overload).
- Flow variables are what matter most: 17 times fewer nonzeros and 13 times less
  solver time for the plain PTDF, 23 times fewer nonzeros for the decomposed one.
- On top of that the decomposition saves 1.7 times the solver time (42 s to
  24-25 s) and about 1.4 times the wall time. 6 areas are slightly better than 4.

## 12. Where the results are

| Folder | Contents |
|---|---|
| `outputs/adversarial/case14toy*` | 14-bus runs (`_line`, `_any_base`, `_d2`, `_rate2`, `_v2`, …) |
| `outputs/sc_profile/ACTIVSg200` | ACTIVSg200 profile |
| `outputs/adversarial/vacc/ACTIVSg200_{d0,d2,d5,mix}` | first master, March week |
| `outputs/adversarial/vacc/ACTIVSg200_v2` | master v2, March week; `eval_fullbase_{false,true}` with `hours_eval.csv` and `diagnosis.txt` |
| `outputs/adversarial/ACTIVSg200_v2_b5`, `ACTIVSg200_v2` | local v2 run; the stopped no-cap run |
| `outputs/adversarial/vacc/ACTIVSg200_jul` | July design; `eval/` (August), `eval_mar/` (March) |
| `outputs/adversarial/vacc/ACTIVSg2000_s1` | ACTIVSg2000 intact network and outage sample |
| VACC `~/tnr_hpc/outputs/adversarial/ACTIVSg2000_line`, `ACTIVSg200_jul_line` | weaker-cut runs (not yet copied back) |
| VACC `~/tnr_hpc/outputs/adversarial/ACTIVSg2000_s1_polish`, `ACTIVSg2000_line_polish` | ACTIVSg2000 intact network, master fixes and polishing (not yet copied back) |

Each run folder has `parts/net_<k>/` (one per network: `summary.txt`,
`networks.csv`, `masks.csv`) and evaluation folders with `summary.txt`,
`forms.csv`, `hours_eval.csv`, `exact_limits.csv` and `networks.csv`
(file formats in `../README.md`).

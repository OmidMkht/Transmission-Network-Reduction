# Decisions, rejected ideas and open questions

What was decided for the per-contingency reduction, and why (28–29 September
2026). Method in `method.md`, numbers in `results.md`.

## Settled

| Decision | Reason |
|---|---|
| **Goal: containment.** Any dispatch feasible on the reduced network must be feasible on the full network. | The earlier feasibility test only checked around the optimum. |
| **Separate from proxy/greedy, no KKT.** Own folder `adversarial/`. | A different approach; KKT is not needed. |
| **Application: preventive SC-DCOPF with N-1 line outages.** Each contingency and the intact case get their own reduced network. | Each contingency adds a network copy $\lvert A x\rvert \le b$; each is replaced by $\lvert A_{\mathrm{red}} x_{\mathrm{red}}\rvert \le b_{\mathrm{red}}$. |
| Preventive: the same dispatch before and after the outage. | |
| Post-contingency ratings equal normal ratings. | |
| Outages that island the network are skipped. | |
| The intact network is reduced too (option b). | |
| **Reference:** the full SC-DCOPF optimum. | |
| **Demand:** the convex hull of the design hours. | "Any demand in the window" as an LP-friendly set. |
| **Cost cap:** 5% above each loading's own SC-DCOPF optimum, interpolated over the hull. | Without a cap the design is about 15 times slower and bigger. A cap on the maximum plain DC-OPF cost was considered and rejected: it misses the security premium and is loose at low load. |
| **Outage adversaries respect the intact limits.** | The intact network's reduced block already guarantees that; it makes every network independent. |
| **The rule: "reduced flow beyond full flow", sign-aware** (dominance cut, no derating, exact witness). | Your rule. It keeps every rating at F and gives exact cost on the design range. |
| **No ratings or derating to choose.** | Asked twice ("Why do we need ratings and deratings at all?"; "we discussed not to consider the flow limits"). The d2/d5 runs remain as experiments only. |
| "Useful" = reduction above 50%. | |
| Convergence: stop when no adversary can violate. | Each round adds a constraint the current design breaks and none is removed (`method.md`, Section 8). |
| One worst adversary per round is acceptable. The code takes the 3 worst overall per round (not per line). | Per-line adversaries would add too many copies when many lines are critical. |
| No early stopping of the adversary checks. | Parallel LP checks were asked for instead (not implemented yet). Only the greedy's pass/fail checks stop at the first failure; that cannot change its answer. |
| **Master v2:** maximise merged lines; critical lines fixed kept; bridges fixed merged; outage start = base design → local unmerge → neighbourhood; MIPFocus 1, 5% gap. | "Make the master more scalable", items 1–4. 3–7 times faster on ACTIVSg200's hard networks. |
| **Design on the peak month.** | The March-week designs overloaded in July because July lies outside the March hull. Designing on all of July removed the overloads on August and March. |
| **Master approaches kept** (1 October, from the July benchmark): dominance cut + hop-limited outage search (`hops`); polishing as an optional add-on; greedy with the exact check as the alternative to the master. | `hops` was smallest, fastest and most accurate. The weaker and any-line cuts, critical-only limits and `combo` cost more time for no smaller networks. |
| **Cluster hop ladder for seedless outage networks:** rungs `[3, 6, free]` by default (adjustable), hold-forward. Not for the intact network (`ladder_base = false`). | Outage networks whose hop region grows to the whole network have no seed to search near. On the July benchmark, the 10 such networks shrank from 666 to 620 buses at 16% more time, with the same accuracy. On the intact network the ladder lost: 4 buses against 3, 12 times slower (`results.md`, Section 10). |

## Removed from the code (2 October)

The state before is commit `155ab51`.

| Removed | Why |
|---|---|
| Weaker cut (`cut=:line`; runs `line`, `line_polish`) | Same reduction as the dominance cut at 2.9 times the design time; with polishing 7.4 times, and one network uncertified at 4 h. |
| Any-line cut (`cut=:any`, `candidates`; runs `any_crit`, `any_all`) | 2.7–3.4 times slower, and no smaller networks. |

Kept for sure: the dominance cut with the hop-limited outage search (`hops`),
polishing (`dom_polish`), the cluster hop ladder, and any combination of them.
Everything else stays in the code for now; to be decided later.

## Rejected

| Idea | Why not |
|---|---|
| Adversary pool (keep 1–3 adversaries per line active, the rest in a pool, re-activate on violation) | Rejected: "dont use the pool idea!". The master keeps every adversary it adds. |
| Derating to give the master slack on ACTIVSg2000 | Brings back the ratings you asked to drop. The weaker cut was chosen instead. |
| Keeping the intact network whole ("full base") | Changes the cost at no hour; overloads outside the range stay the same. |
| A pure hop limit for the intact network | On ACTIVSg2000 the 289 critical lines are spread everywhere: 3–4 hops still cover 58–73% of the undecided lines. Every copy stays full size, and far lines still shift critical flows (201 lines at 4 hops by more than 1%). |
| The first master (bus count with root and commodity flows) | Replaced by v2: slower and larger. |
| Free per-line rating variables | They crept down one hair per round and the witness set grew. |
| The any-line cut (`:any`) | Exact, but did not converge in 10 minutes on the 14-bus base network. |

## The weaker cut

**Why the intact network of ACTIVSg2000 did not reduce.** The dominance cut asks
the reduced flow on the overloaded line to be at least the full flow. The
unreduced network meets that with equality, so there is no room: any merge that
lowers one of the cut flows even slightly is forbidden. On ACTIVSg200 every
intact-network adversary was on line 184, and the master could find merges that
push more power through that one line. On ACTIVSg2000 the 6 cuts sit on
different lines, and almost every merge lowers at least one of them a little. The
design kept everything except the bridges (1,550 of 2,000 buses). The number of
critical lines only sets a floor (289 kept lines would still allow a few hundred
buses).

**Update: this explanation was incomplete.** The weaker cut, which leaves plenty
of room at the unreduced network, ended at the same 1,550 buses
(`results.md`, Section 9). Two other things block merging on ACTIVSg2000:
- **Unneeded limits.** Kept lines that are not critical carry limit F during the
  design. Merging far regions pushes some of them over F at the secure optimum,
  so the witness is rejected. With limits on critical lines only, every witness
  passes.
- **Real hidden overloads.** Merging only lines 7 or more hops from any critical
  line still lets 52–63 overloading dispatches through. That is the safety
  condition itself; no choice of cut relaxes it.

**Ways to give the master room** (all keep the guarantee on the design range):

| Option | Room for merges | Guarantee | Cost | Downside |
|---|---|---|---|---|
| **Weaker cut** (`cut=:line`): the reduced flow only has to exceed the rating | how far each adversary overloads (12–104% on ACTIVSg2000) | exact | exact | each cut rejects only that dispatch, so more rounds (20 times slower on 14 buses, same reduction) |
| Tolerance ε: the reduced flow may fall short of the full flow by up to ε of the rating | ε | any accepted dispatch overloads the full network by at most ε | exact | the guarantee is "within ε" |
| Derating δ: kept lines rated $(1-\delta)F$ | δ | exact | slightly higher | changes the ratings |

**Chosen: the weaker cut** (29 September), running on ACTIVSg2000 and on
ACTIVSg200 (July design) for comparison.

## Finding critical lines

The approach in the code matches the one you described:

- **Per contingency, 2L LPs** maximising and minimising each post-outage flow over
  the adversary set, with the cost cap.
- **Each contingency keeps only its own critical lines.**
- **Most LPs are skipped** by the LODF bound. On ACTIVSg200 it clears 99.76% of
  84,424 pairs, leaving 200 LPs.

Improvements not implemented yet:
- **Stored dispatches.** Before solving an LP, test the dispatches already at hand
  (the intact-extreme LP optima and the hourly optima). If one overloads the
  pair, it is critical with no LP. On ACTIVSg200 that leaves 33 LPs instead of 200.
- **Solve the 2L intact LPs once and share them.** At 2,000 buses they are the
  heaviest part (6,412 LPs with 6,412 dense rows, about 1,000 s per job), and
  every outage job currently re-solves them.

## Protecting fewer critical pairs

For an exact guarantee every critical pair must be protected: the reduced
network must reject every dispatch that overloads it. It does not need its own
cut. Today only the pairs whose adversaries break a design get one (6 cuts for
398 pairs on the ACTIVSg2000 intact network); the check tests all of them.

| Option | Guarantee | ACTIVSg200 (168 pairs) | ACTIVSg2000 intact (398) |
|---|---|---|---|
| Drop exact duplicates | exact | removes 0 | removes 5 |
| Drop a pair whenever its overload implies another protected pair's (one LP per candidate) | exact | not tested | not tested |
| Drop pairs whose worst possible overload is ≤ ε | overload ≤ ε | ε=1%: drops 56; ε=2%: 75 | ε=1%: drops 7; ε=10%: 74 |
| Protect only pairs that bind in some hour's secure OPF | none | – | – |

## Proposed, not implemented

**Hop-limited outage networks.** Start from the base design and make only the
lines within h hops of the outage line and of the network's own critical lines
free decisions. Fix everything else to the base design, and grow h if no design
is found. This replaces the `radius` rule, which keeps every binary in the model
and only caps how many change. With the far lines fixed, presolve shrinks every
copy. The guarantee is unchanged, since the checks run on the full network.

**Region plus exterior equivalent (master v3).**
1. **Region:** the critical lines and every line within h hops of them and of the
   outage line (or every line whose merge can shift a critical flow by more than
   a threshold, measured from the PTDF). Only the region's lines are decisions;
   the exterior keeps the start design's choice.
2. **Exact equivalent of the exterior:** with its merges fixed, the exterior is a
   fixed DC network. Kron reduction onto the boundary buses gives
   - equivalent lines among the boundary buses,
     $Y_{eq} = B^{ext}_{BB} - B_{BE} B_{EE}^{-1} B_{EB}$;
   - a fixed map moving exterior injections onto the boundary,
     $W = -B_{BE} B_{EE}^{-1}$.

   Both are computed once per network with one sparse factorisation.
3. **Region-sized copies:** each adversary or witness copy has only the region's
   buses and lines plus $Y_{eq}$ on the boundary, with injection
   $p_{\mathrm{region}} + W p_{\mathrm{exterior}}$. For a fixed exterior this is
   exact, not an approximation: DC flow is linear.
4. **Checks stay on the full network**, so the certificate is unchanged.
5. **Growth:** if no design exists inside the region, grow it and recompute the
   equivalent. In the worst case it becomes today's master.
6. **Requirement:** limits only on critical lines during the design (the witness
   copies then need only in-region flows).

Costs:
- merges outside the region are only reconsidered when the region grows;
- $Y_{eq}$ is dense among the boundary buses;
- many scattered critical lines make a large region.

## Open questions and next steps

As of 2 October 2026:

1. **ACTIVSg2000: the condition, not the search, blocks the reduction.**
   - With original ratings and zero tolerance, the exact-check greedy merges only
     10 lines beyond the bridges on the intact network (1,540 of 2,000 buses), in
     both the SC-DCOPF and the plain DCOPF range (`results.md`, Section 9).
   - Options for room: an overload tolerance ε (about 1%, with the greedy), or a
     narrower cost cap. Not decided.
2. **Screening at scale:** share the intact extremes across jobs. Each 2,000-bus
   outage job spends about 1,000 s re-solving the same 6,412 LPs. Add the
   stored-dispatch test too.
3. **Check speed:**
   - Parallel adversary LPs (one model per thread).
   - Rank-one updates of the reduced PTDF instead of rebuilding it (1–2 s per
     check at 2,000 buses).
   - Already done: lazy reduced rows, and stopping at the first failure for
     pass/fail.
4. **Ladder at scale:** lazy hop-cap rows. ACTIVSg2000 has 437,000 paths at k = 4
   and over 2 million at k = 6; today rungs past 200,000 paths are skipped.
5. **Master scalability beyond hops:** region plus exterior equivalent
   ("Proposed, not implemented").
6. **Design range:** all of July covered August and March on ACTIVSg200. For
   ACTIVSg2000 only January–March hourly data exists; July would need new
   scenarios (`common/make_hourly_scenarios.jl`).
7. **ACTIVSg200 outages 73, 193, 195 and 242** stay at 128 buses (bridges only) in
   every run.
8. **Still to decide** whether to keep in the code: critical-only limits, the
   radius search, the full base, derating, the tolerance ε, greedy, the
   main-network mode, compare-only mode and the planning test.

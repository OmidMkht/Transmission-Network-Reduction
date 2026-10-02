# Method: per-contingency reduced networks for a preventive SC-DCOPF

This is the method as implemented in `adversarial/` (29 September 2026). Code
references point to `adversarial/adversarial.jl` unless another file is named.

## 1. The problem

### Network and notation

| Symbol | Meaning |
|---|---|
| $N, L, K$ | buses, lines, generators |
| $u_\ell, v_\ell$ | from and to bus of line $\ell$ |
| $b_\ell = 1/x_\ell$, $F_\ell$ | susceptance and rating of line $\ell$ |
| $C$ | generator-to-bus map |
| $H$ | PTDF, $H = B E^\top (E B E^\top)^\dagger$ with a reference bus (`ptdf`) |
| $\mathcal{C}$ | outages: lines whose removal keeps the network connected (`outages_of`); islanding outages are skipped |
| $H^c$ | PTDF after outage $c$: $H^c = H + \mathrm{LODF}_{:,c} H_{c,:}$ (`post_ptdf`) |
| $\mathrm{LODF}_{\ell c}$ | $a_\ell / (1 - a_c)$ with $a = H(e_{u_c} - e_{v_c})$, and $\mathrm{LODF}_{cc} = -1$ (`lodf`) |
| $\mathcal{G}$ | generator box $[\min(0, g^{\min}),\, g^{\max}]$ (minimum outputs relaxed to 0) |
| $c^\top g + c_0$ | linear generation cost |
| $d_1, \dots, d_T$ | the design hours (bus demand vectors), columns of $D$ |
| $p = Cg - d$ | net injection |

The LODF is checked against a PTDF rebuilt with the line removed (worst gap
about $10^{-13}$; the runner checks up to 20 outages spread over the list).
Post-outage ratings equal normal ratings.

### Preventive SC-DCOPF

The same dispatch must hold before and after every outage:

$$
Z(d) = \min_g\ c^\top g + c_0 \quad\text{s.t.}\quad g \in \mathcal{G},\ \ \mathbf{1}^\top g = \mathbf{1}^\top d,\ \
|H_\ell p| \le F_\ell\ \forall \ell,\ \ |H^c_\ell p| \le F_\ell\ \ \forall c \in \mathcal{C},\ \ell \ne c.
$$

$Z_s$ and $g^\star_s$ are its optimal cost and dispatch at design hour $s$. They are
computed once by constraint generation (`sc_dcopf`): start with the intact rows,
solve, compute every post-outage flow $f_\ell + \mathrm{LODF}_{\ell c} f_c$, add the
violated rows, repeat. One pool of rows serves all hours.

### What we want

The SC-DCOPF has one block per network $c \in \{0\} \cup \mathcal{C}$ ($c = 0$ is
the intact network). Each block is replaced by its own **reduced network**. The
requirement is containment:

> every dispatch that the reduced networks accept is secure on the full network

for every demand in the design range and every dispatch within the cost cap
(Section 2). "Secure" means all intact and post-outage limits hold on the full
network.

## 2. The design range (adversary set)

The guarantee is made for a set of dispatches, the **adversary set**. Its
variables are the generator outputs $g$ and demand weights $\lambda$
(`adversary_set`):

- generators: $g \in \mathcal{G}$
- demand anywhere in the convex hull of the design hours: $d = D\lambda$, $\lambda \ge 0$, $\mathbf{1}^\top\lambda = 1$
- balance: $\mathbf{1}^\top g = \sum_t \lambda_t\, \mathbf{1}^\top d_t$
- **cost cap**: $c^\top g + c_0 \le (1 + x) \sum_t \lambda_t Z_t$, with $x = 5\%$ (`band`)

This is $\mathcal{A}^0$, used for the intact network. For an outage network the
dispatch must also respect every intact limit:
$\mathcal{A}^c = \mathcal{A}^0 \cap \{|H_\ell p| \le F_\ell\ \forall \ell\}$.

Why each piece is there:

- **Cost cap.** The SC-DCOPF only ever returns near-economic dispatches. Without
  the cap, the adversary can pick absurd dispatches that force the reduced
  network to stay large; a run without it was about 15 times slower and bigger
  (see `results.md`). The cap removes only dispatches dearer than the secure
  optimum plus 5%; cheaper, insecure ones stay in.
- **Interpolated cap.** $Z$ is convex in demand, so $\sum_t \lambda_t Z_t \ge Z(D\lambda)$.
  At a mixed demand the cap is therefore never tighter than 5% above that
  demand's own secure optimum; it errs on the side of more adversaries.
- **Intact limits for outages.** The reduced SC-DCOPF also contains the intact
  network's block, which is certified on $\mathcal{A}^0$. Any dispatch the reduced
  SC-DCOPF accepts is therefore intact-feasible, so outage networks may assume
  it. This lets every network be designed on its own.

## 3. Reduced networks

A **design** for network $c$ marks lines as merged. Clusters are the connected
components of the merged lines; a line whose two ends end up in the same cluster
is closed as well (it carries no flow). Kept lines keep their reactance and get a
rating $r_\ell$: $F_\ell$ by default, $(1-\delta)F_\ell$ with derating, or none
($\infty$) after exact limits.

For a fixed design the kept-line flows are $\varphi(p) = H^{\mathrm{red}} p$: the
PTDF of the cluster network (cluster injection = sum of its buses' injections),
mapped back to buses (`ptdf(n; internal=mask, out=c)`). The reduced network
**accepts** $p$ when $|\varphi_\ell(p)| \le r_\ell$ on every kept line.

In the reduced SC-DCOPF (`evaluate.jl`, `network_scdcopf`) a block has one angle
per cluster, a flow per kept line with its rating, and one balance row per
cluster. The generator outputs are shared by all blocks.

## 4. Conditions and the guarantee

Each network's reduced design must meet two conditions:

1. **No hidden overload.** Every $p \in \mathcal{A}^c$ that the reduced network
   accepts satisfies $|H^c_\ell p| \le F_\ell$ for all $\ell \ne c$ (for $c = 0$:
   $|H_\ell p| \le F_\ell$).
2. **Witness.** The reduced network accepts the witness injection
   $\hat p_s = C \hat g_s - d_s$ at every design hour $s$. The witness is the exact
   secure optimum, $\hat g_s = g^\star_s$ (with derating: the SC-DCOPF optimum at
   ratings $(1-\delta)F$, the "margin witness").

**Proposition.** If every network meets both conditions, then:

- the reduced SC-DCOPF is feasible at every demand in the hull;
- every dispatch it accepts inside the cost cap is secure on the full network;
- at each design hour its optimum costs exactly $Z_s$ (with the exact witness).

*Why.* The accepted sets are convex, so mixing the hour witnesses
($\sum_s \lambda_s \hat g_s$) gives an accepted dispatch at demand $D\lambda$ with
cost $\sum_s \lambda_s Z_s$, inside the cap. An accepted dispatch in the cap lies
in $\mathcal{A}^0$; the intact network's condition 1 makes it intact-feasible, so
it lies in every $\mathcal{A}^c$, and each outage network's condition 1 makes it
secure. At a design hour $g^\star_s$ is accepted, so the reduced optimum costs at
most $Z_s$; it is secure, so it costs at least $Z_s$.

**Scope.** The guarantee holds only inside the design range. Outside it the
reduced networks can be wrong in both directions (too loose: overloads; too
strict: cost above the optimum). The March-week designs tested on July showed
both (`results.md`, Section 5).

## 5. Finding the critical lines (screening)

A line and direction $(\ell, \sigma)$ of network $c$ is **critical** if some
dispatch in the range can overload it:

$$
u^c_{\ell\sigma} = \max\{\sigma H^c_\ell p : p \in \mathcal{A}^c\} > (1 + \tau) F_\ell, \qquad \tau = 10^{-6}.
$$

Only critical pairs need protecting; the others cannot overload anywhere in the
range. A network with no critical pair is **redundant** and dropped from the
reduced SC-DCOPF. Each network keeps only its own critical lines.

That is up to $2L$ LPs per network, $2L^2$ in all. Most are skipped
(`screen_all` in the runner, `screen`):

- **Intact extremes.** $2L$ LPs give $\overline h_\ell = \max H_\ell p$ and
  $\underline h_\ell = \min H_\ell p$ over $\mathcal{A}^c$. They are shared by every
  outage.
- **LODF bound.** After outage $c$ the flow is $f_\ell + q f_c$ with
  $q = \mathrm{LODF}_{\ell c}$, so
  $u^c_{\ell,+} \le \overline h_\ell + \max(q\overline h_c,\, q\underline h_c)$ (and
  the mirror for $\sigma = -1$). A pair whose bound is within the rating is
  cleared without an LP.
- An LP is solved only for pairs the bound cannot clear.

A run that designs only some networks (`network_index`, `networks`) screens only
those networks, and skips the intact extremes when only the intact network is
screened.

## 6. Checks for a fixed design

`check(n, c, Hc, A, crit, Pw, mask; rating)`:

- **Witness check.** For every design hour and kept line:
  $|H^{\mathrm{red}}_\ell \hat p_s| \le r_\ell + \tau F_\ell$. For every kept line and
  sign the worst failing hour is collected (distinct hours).
- **Adversary check.** For every critical pair $(\ell, \sigma)$, one LP:
  $$v = \max\{\sigma H^c_\ell p : p \in \mathcal{A}^c,\ |H^{\mathrm{red}}_m p| \le r_m\ \forall \text{ kept } m\}.$$
  If $v > (1 + \tau) F_\ell$, the maximiser is an **adversary**: a dispatch the
  reduced network accepts that overloads the full network. Adversaries are
  sorted by loading $v / F_\ell$. Taking the max over all pairs gives the single
  worst adversary.
- A design is **certified** when no witness is rejected and no adversary exists.
- The reduced network's rows go into each LP lazily: a kept line's row is added
  once the maximiser breaks it, and stays for the remaining pairs. Exact, and each
  LP holds only the rows that bind.
- With `stop=true` the check returns at the first failure (rejected witnesses, or
  one adversary), for a pass/fail answer.

The runner repeats the adversary check with the reduced network in angle form
(one angle per bus, merged lines forced to zero angle difference) for every
certified design; the two forms agree to about $10^{-14}$.

## 7. The master problem (v2)

`master`, `flow_copy!`, `reject!`, `near!`. One MIP per network, kept across
rounds.

**Decision.** $z_\ell \in \{0,1\}$ for every line $\ell \ne c$; $z_\ell = 1$ merges it.

**Fixed before any solve.**
- Critical lines are fixed kept ($z = 0$): they are the ones that reject adversaries.
- With `fix_radial`, every non-critical bridge is fixed merged ($z = 1$): merging a
  bridge changes no other line's flow, since the power beyond it has one way out
  (`bridges`, iterative Tarjan).

**Objective.** Maximise $\sum_\ell z_\ell$, i.e. keep as few lines as possible. A
line closed by other merges counts as merged too; physically it makes no
difference.

**No bound carried between rounds.** An earlier version capped $\sum_\ell z_\ell$ at
the previous round's bound. On ACTIVSg2000 Gurobi reported bounds that designs it
had not found beat (`results.md`, Section 9), so the cap was removed.

**Neighbourhood** (outage networks, see Section 8):
$\sum_{\ell:\ \text{centre merges } \ell} (1 - z_\ell) + \sum_{\ell:\ \text{centre keeps } \ell} z_\ell \le \text{radius}$.

**Flow copy** at a fixed injection $p$ (one per adversary and per witness hour).
Variables: angle $\theta_i$ per bus, flow $f_\ell$ and free transfer $t_\ell$ per line.

| Constraint | Meaning |
|---|---|
| $\theta_{\mathrm{ref}} = 0$ | reference angle |
| $f_\ell = b_\ell(\theta_{u_\ell} - \theta_{v_\ell})$ | DC flow |
| $\sum_{u_\ell = i}(f_\ell + t_\ell) - \sum_{v_\ell = i}(f_\ell + t_\ell) = p_i$ | balance at every bus |
| $z_\ell = 0 \Rightarrow t_\ell = 0$ (indicator) | only merged lines carry free transfers |
| $\lvert f_\ell \rvert \le P$, $\lvert t_\ell \rvert \le P$, $P = \sum_i \max(p_i, 0)$ | bounds that let Gurobi tighten the indicators: with positive susceptances the contracted network's DC flow has no cycle, and routing inside each cluster over a tree adds none, so some solution carries at most $P$ on every line (was $(L+1)P$ for $t$) |
| adversary copy: $z_\ell = 1 \Rightarrow f_\ell = 0$ (indicator), no limits | merged lines carry no DC flow |
| witness copy: $\lvert f_\ell \rvert \le (r_\ell + \tau F_\ell)(1 - z_\ell)$ | merged lines carry no DC flow; kept lines respect their rating |

With $z$ fixed, the copy's $f$ on kept lines equals $\varphi(p)$: merged lines
force equal angles inside a cluster, and the transfers move power freely inside
it.

**Dominance cut** for adversary $k$ that overloads line $\ell$ in direction
$\sigma$ with full flow $v_k$: $\sigma f^k_\ell \ge (r_\ell / F_\ell)\, v_k$, i.e.
"reduced flow beyond full flow", sign-aware; with $r = F$ it is
$\sigma f^k_\ell \ge v_k$. No binaries.

The reduced network then rejects $p^k$ by line $\ell$ itself. The cut asks for more
than rejection, so it also rejects nearby dispatches, which means fewer rounds. Its
price is that the unreduced network meets it with equality, leaving merges no room
(`decisions.md`, "The weaker cut"). The weaker cut (`:line`) and the any-line cut
(`:any`) were tried and removed (`results.md`, Section 10).

**Master variants** (settings; all keep the exact check):

| Setting | What changes |
|---|---|
| `tolerance` ε | The guarantee becomes "accepted dispatches overload the full network by at most ε". A pair is critical only if it can exceed $(1+\varepsilon)F$, adversaries count only above $(1+\varepsilon)F$, the dominance cut becomes $\sigma f^k_\ell \ge (r_\ell/F_\ell)(v_k - \varepsilon F_\ell)$, and exact limits use $(1+\varepsilon)F$. |
| `limits = :critical` | During the design only the critical lines carry a rating. The others are unlimited, so witness copies limit only critical lines. |
| `hop_limit` h | For outage networks, instead of the `radius` rule: every line more than h hops from the outage line and from the critical lines is fixed to the start design; h doubles whenever the master finds nothing. |
| `polish_time` | Seconds per round for polishing (below); 0 turns it off. |

**Solver settings.**

| Setting | Value | Why |
|---|---|---|
| FeasibilityTol, IntFeasTol | $10^{-8}$ | cuts and witness limits sit at ratings and the check allows $\tau F$ more (at least $2\times10^{-7}$ on ACTIVSg2000), so the tolerances must stay well below that. $10^{-9}$, Gurobi's minimum, gave invalid bounds at 2,000 buses |
| NumericFocus | 2 | same reason |
| MIPFocus | 1 | good designs fast; the LP check certifies them anyway |
| MIPGap | 5% (`master_gap`) | same reason |
| time per solve | `master_time` (20 s local, 60 s on VACC, 600 s for ACTIVSg2000) | |
| Threads | `master_threads` | |
| NodefileStart | env `ADV_NODEFILE_START` (GB) | long solves spill the tree to disk |
| start values | the centre design during a local search; otherwise only the bridges merged | "bridges merged" meets every cut and witness |

## 8. The design loop, per network

`design(n, c, Hc, A, crit, Pw; ...)`:

1. Screen; a redundant network is dropped.
2. Fix the critical lines kept and the non-critical bridges merged.
3. **Outage start** (outage networks, when the base network's design is given):
   - try the base design (with the fixings applied); if it passes both checks, done;
   - otherwise try it with the merged lines within `unmerge_hops` (2) of the
     outage line's ends unmerged; if that passes, done;
   - otherwise whatever they failed on (witness hours, adversaries) goes into the
     master, and the unmerged design becomes the **centre**: the master may change
     at most `radius` (20) lines of it.
4. **Round.**
   - Solve the master within its time limit.
   - No feasible point: during a local search double the radius (once it reaches
     the number of lines the neighbourhood is dropped); during a global search
     return the unreduced network.
   - Otherwise read the design and **polish** it (below, up to `polish_time`
     seconds), then check it.
   - Certified: stop.
   - Otherwise add up to `witnesses_per_round` (10) witness copies and the
     `adversaries_per_round` (3) worst adversaries with their cuts, and repeat.
5. **Limits.** Past `network_time` or `max_rounds`, return the unreduced network
   (status `limit`). It always passes, since it is the full network.
6. **Derating only.** Raise the kept ratings back towards $F$ by bisection as far
   as the check allows (`raise_ratings`).

**Polishing** (`polish`, `merge_slack`, `copies_slack`). The master can miss
merges, and with bad numerics claim none exist. Polishing adds them without the
MIP:
- Every single merge from the current design is scored exactly against every
  copy the master holds. Merging line $j=(a,b)$ moves line $\ell$'s flow by
  $-f_j\,T_{\ell j}/T_{jj}$, with $T_{:,j} = H^{\mathrm{red}}_{:,a} - H^{\mathrm{red}}_{:,b}$
  (infinite susceptance on $j$).
- The merges that keep every copy satisfied are taken best first in one batch.
  The batch is re-checked exactly on the contracted network and halved until it
  passes.
- Scoring repeats on the new design until no single merge fits.

Critical lines stay kept. The polished design satisfies everything the master
knows, so it is at least as good as the master's; the check then certifies it or
yields new copies. On the ACTIVSg2000 round the master called optimal at 1,550
buses, polishing reached 110 buses in 25 s with every copy still satisfied.

**Why the loop ends.** Each round adds a constraint that the current design
breaks: the design accepted the adversary, so its flow on the line was within the
rating, below what the cut asks; likewise for a rejected witness. Constraints are
never removed, so no design repeats, and there are finitely many designs. The
neighbourhood only delays this: it doubles until it covers the whole network.
The loop stops when no adversary is left.

### Cluster hop ladder (`cluster_hops`)

For searches with no seed to stay near: the intact network, or an outage network
whose neighbourhood or hop region has grown to the whole network.

- **Rung with hop cap k:** no chain of more than k merged lines inside a cluster.
  - One row per simple path P of k+1 lines: $\sum_{l \in P} z_l \le k$
    (`common/caps.jl`, shared with greedy, KKT and proxy).
  - Bridges and critical lines are left out of the paths: bridges are fixed
    merged and change no flow, critical lines are fixed kept.
  - Rows on ACTIVSg200: about 1,400 at k = 3 and 9,800 at k = 6.
  - A rung with more than 200,000 paths is skipped (ACTIVSg2000 from k = 4); lazy
    rows are not implemented.
- **Each rung runs the usual loop** (master, polishing within the cap, check)
  until certified. Its merges are then fixed merged in the next rung
  (hold-forward). The adversaries and witnesses carry over.
- **Time:** each rung gets an equal share of what is left; the last rung gets the
  rest.
  - A rung that runs out of time, or whose master finds nothing, passes the last
    certified design on.
  - On the time or round limit, the last certified rung's design is returned
    instead of the unreduced network.
- Default rungs `[3, 6, nothing]`, the last one free.
- By default only outage networks use it; `ladder_base = true` adds the intact
  network. On ACTIVSg200 the free master did better there.

The master's witness copies hold the rating itself, while the check allows
$\tau F$ more. Otherwise a design that met a witness copy only within solver
error was rejected round after round.

### Greedy instead of the master (`search = :greedy`)

No master, no copies. Every merge is kept only if the exact check passes, so each
design on the way is certified.
- Start: every non-critical bridge merged.
- Order: by how far merging line $j$ alone can shift a kept critical line's flow
  at $j$'s rating, $\max_\ell |T_{\ell j}| F_j / (T_{jj} F_\ell)$ (Section 7's
  rank-one formula), smallest first.
- Batches of lines: the size doubles after a pass and halves after a failure. A
  line that fails alone waits for the next sweep.
- Each check stops at the first failure, and the critical pair that failed last
  goes first.
- Sweeps repeat, re-ranked, until one merges nothing.

## 9. After the design

**Exact limits** (`exact_limits=true`, `exact_limits`). Each critical line gets the
highest rating at which it still rejects, by itself, every dispatch that loads
the true line to its rating or beyond:

$$
r_\ell = \min_\sigma\ \min\{\sigma \varphi_\ell(p) : p \in \mathcal{A}^c,\ \sigma H^c_\ell p \ge F_\ell\} - \tau F_\ell .
$$

It can land above or below $F_\ell$. Every other kept line gets no limit, since
screening showed it cannot overload in the range. The result must pass the check
again; if it does not, or a critical line is merged, the design keeps its own
ratings. The runner also reports whether each design accepts the exact optimum
at every hour.

**Combining designs.** Every network is certified on its own against the same
adversary set, so designs from different runs can be mixed network by network,
provided all of them accept one common witness.

**Full base** (`full_base=true`, evaluation only). The intact block is kept as
the full network and only the outages are reduced.

## 10. Evaluation

The runner solves every SC-DCOPF form at every design hour and every test hour
(`evaluate.jl`: each model is built once and only the right-hand sides change per
hour):

| Form | What it is |
|---|---|
| `full` | every block (intact and every outage) as a full network copy: an angle per bus, a flow per line, limits |
| `lazy` | PTDF rows added only when violated (the reference; also gives $Z$) |
| `compact` | only the screened (block, line) pairs, as PTDF rows over the generator outputs |
| `reduced` | the designed networks, redundant ones dropped, design ratings |
| `reduced_exact` | the same with exact limits |

Per hour it records (`hours_eval.csv`):
- the cost and the signed gap to the secure optimum (negative would mean cheaper
  than possible, i.e. insecure);
- the worst loading of the returned dispatch on the full network, over the intact
  network and every outage (via the LODF), and where it occurs;
- the demand's distance from the design hull: the L1 distance to the nearest
  $D\lambda$, as a percentage of total demand (0 means inside, where the guarantee holds).

With `planning=true` (14-bus only) each candidate line of
`case studies/14bus/network_planning.csv` is added to the full network and to
every reduced block, to test the network form against compact rows computed
beforehand.

# transmission_plots.jl

Before/after network figure shared by both active runners
([transmission_reduction_edge.md](transmission_reduction_edge.md),
[transmission_reduction_edge_multiscenario.md](transmission_reduction_edge_multiscenario.md))
via `TxReport.plot_network` in
[transmission_reduction_reporting.md](transmission_reduction_reporting.md).

This is a **rewrite** of the version documented in `others/transmission_plots.md`
— that file describes an earlier, differently-colored, differently-labeled
version and applied to the older assignment-based drivers in `others/`; it
was not updated when this file changed and should be treated as historical
only. This doc describes the file as it exists at the top level now.

Deliberately a **separate file, not part of either model module**: CairoMakie
is the slowest thing in the directory to load and no solving needs it, so a
plain solve run never pays for it — `TxReport.ensure_plots` loads it lazily,
on first actual use. It also reaches into nothing but its arguments — the
retained set and the `bus → super-node` map are re-derived from `A` — so it
works with the case and matrix from either runner.

```julia
include(joinpath(@__DIR__, "transmission_plots.jl"))
plot_reduction(c, A; path="outputs/case118.png", title="case118")
```

## Design

**Visual priority.** The reduced network *is* the external lines — they are
what survives. So external lines are drawn dark and thick, and the (usually
far more numerous) internal lines are demoted to a soft cluster-colored web:
a wide low-alpha halo plus a thin core line, tracing that cluster's *actual*
footprint. A convex hull was considered and rejected: on a 200-bus case one
cluster can hold 185 buses spread across the whole layout, and its hull
would visually claim every other cluster's buses too. Stroking the cluster's
own internal edges can never over-claim.

**Color.** Cluster hues come from a validated categorical palette (7 hues;
the reference palette's 8th slot, red, is held out and reserved for
congestion so a cluster color can never be mistaken for a congestion
marker). Hues are assigned by **greedy graph coloring of the
cluster-adjacency graph** — two clusters are adjacent only if a surviving
external line joins them — largest cluster first. This means only clusters
that actually touch need distinct hues, so the palette never has to grow
with the cluster count; non-adjacent clusters legitimately reuse a hue, the
ordinary map-coloring convention.

## `plot_reduction(c, A; path=nothing, kwargs...)`

| Panel | Shows |
|---|---|
| **Left (original)** | every bus, every line. Each merged cluster's own internal lines get the soft colored web described above. External (surviving) lines are drawn dark and thick on top. Unmerged single buses stay neutral grey. |
| **Right (reduced)** | retained buses only, marker area ∝ `√(cluster size / largest cluster size)` (bounded, so one dominant cluster can't swallow the panel), external lines only. |

### Protected lines <a name="binding_lines"></a>

Lines the model forbids from ever becoming internal are stroked near-black on
**both** panels — they must be visible surviving on the right. Pass the
model's own protected set explicitly via **`binding_lines`**:

```julia
binding_lines = findall(r.protected)                    # multi-scenario
binding_lines = sort(unique(vcat(r.Lplus, r.Lminus)))    # single-scenario
```

Without it, the fallback is `|fhat| == frate` on `c.fhat`/`c.frate` directly
— correct **only** for a genuine single-hour `TxReductionCase`. Passing
`c.base` from a `MultiScenarioTxReductionCase` without also passing
`binding_lines` would use `c.base.fhat`, the stale single-OPF snapshot
`build_multiscenario_tx_case` builds once for topology/generator data only —
unrelated to any modeled hour and never refreshed by the operating-point
scale or redispatch. Always pass `binding_lines` explicitly for a
multi-scenario plot.

### Parallel lines

Contracting endpoints routinely leaves several external lines joining the
same pair of super-nodes. They remain distinct lines with their own ratings
and flows — the model never merges them — but drawn straight they'd lie
exactly on top of one another, silently under-reporting how many lines
survive. The reduced panel instead **fans** them into separate arcs
(quadratic Bézier, symmetric curvature spread) and labels the bundle
`×n`. The panel title's "`N lines over M connections`" is the same count
made explicit.

### Congestion

Lines passed via `congested_lines` are stroked in the reserved red on both
panels (with a white halo so the stroke reads over a colored cluster web or
another line), and labeled by line number (`congested_line_labels` for a
custom label per line, else `L<n>`) — capped at `max_congested_labels`
before the midpoint annotations are dropped as unreadable clutter. A
congested line that ended up **internal** (a coverage/protection failure —
see the multiscenario doc's scenario-selection cover requirement) triggers
an explicit warning in the figure's footer, not just a silent miscoloring.

### Text readability

All labels (`L`-numbers, `×n` bundle counts, retained-bus numbers) use a
`_halo_text!` helper rather than Makie's built-in `text!(strokewidth=...)`.
Makie centers a text stroke *on* the glyph path, so half of it lands inside
the letter — a halo wide enough to read over a busy, colored background
(3–4px) is then wider than the stems of a 14–17px bold glyph, and the label
renders as an unreadable white smear. `_halo_text!` instead draws eight
copies of the same text offset in a ring **behind** the real text, leaving
the glyph itself untouched no matter how wide the halo needs to be.

### One layout, one frame

Coordinates come from a stress layout of the **original** graph
(`network_layout`, `algorithm=:stress` default — reads better on meshed
transmission graphs than `:spring`), computed once, and a super-node keeps
its own bus's coordinate — so a node can be traced from left to right by
eye. That alone isn't enough: each `Axis` otherwise autoscales to its own
data, and the right panel holds a subset of the nodes, so it would render at
a different zoom and defeat the comparison. Both axes are pinned to the full
layout extent. pglib cases carry no geographic data, so positions are
topological only. Pass `pos=` to reuse one layout across several figures —
the only way to compare two *different* reductions of the same case side by
side.

### Keywords

`path` (save location; `nothing` returns the `Figure` without writing),
`title`, `algorithm`/`seed` (layout), `pos` (reuse a layout), `resolution`
(default `(1700, 950)`), `node_size`, `label_retained` (annotate super-node
bus numbers — useful up to ~30 buses), `congested_lines`,
`congested_line_labels`, `label_congested`, `binding_lines` (see above),
`max_congested_labels`.

## In the runners

Both runners route through `TxReport.plot_network(plots_file, case, A;
path, open_after, kwargs...)` rather than calling `plot_reduction` directly
— that wrapper is what defers the `include` (see
[transmission_reduction_reporting.md](transmission_reduction_reporting.md#lazy-plotting)).
`cfg.show.plots = false` skips the figure (and the CairoMakie load)
entirely; `cfg.show.open_plots = true` also opens each saved figure in the OS
default viewer via `TxReport.open_externally`.

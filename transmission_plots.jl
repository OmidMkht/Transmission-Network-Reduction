# --------------------------------------------------------------------------- #
# Before/after graph plots for the reduction. Loaded lazily by tnr_reporting.jl
# (CairoMakie is the slowest include in the project, so a run with plots off
# never pays for it) -- or standalone:
#
#   include("transmission_plots.jl")
#   plot_reduction(c, A; path="outputs/case118.png")
#
# `A` is the assignment matrix (A[i,j]=1 => bus j is in cluster i); the maps
# are re-derived from it, so this works with any runner's output.
#
# The two panels share one layout (coordinates computed once on the original
# graph, a super-node keeps its bus's coordinate) so they're comparable by eye.
#
# External lines are drawn dark and thick -- they're what the reduced network
# actually is. Internal lines become a soft cluster-coloured halo+core instead
# of a convex hull, because a hull lies on a sprawling cluster (one ACTIVSg200
# cluster holds 185/200 buses and its hull would cover everyone else's too);
# stroking the cluster's own edges traces its real footprint.
#
# Cluster hues come from Colors.jl's distinguishable_colors, one per cluster --
# every cluster gets its own colour, not just adjacent ones. White, black and
# the reserved congestion red are seeded in so no cluster is confusable with
# the background, text, or a congestion marker.
# --------------------------------------------------------------------------- #

using CairoMakie
using Colors
using Graphs
using NetworkLayout
using GeometryBasics: Point2f

const LEGEND_SWATCH = colorant"#2a78d6"   # representative example only, not an actual cluster's color
const SURFACE      = colorant"#fcfcfb"
const INK_PRIMARY  = colorant"#0b0b0b"
const INK_SECOND   = colorant"#52514e"
const INK_MUTED    = colorant"#898781"
const CONGESTED_COLOR = colorant"#d03b3b"   # status:critical, reserved
const SINGLETON_COLOR = colorant"#898781"

# retained buses + bus -> super-node map, straight from A (see extract_reduction)
function _reduction_maps(A::AbstractMatrix)
    N = size(A, 1)
    retained = [i for i in 1:N if A[i, i] > 0.5]
    rep = zeros(Int, N)
    for j in 1:N
        i = findfirst(i -> A[i, j] > 0.5, 1:N)
        i === nothing && error("bus $j is unassigned in A")
        rep[j] = i
    end
    return retained, rep
end

function _simple_graph(c)
    g = SimpleGraph(c.N)
    for l in 1:c.Ln
        add_edge!(g, c.Efrom[l], c.Eto[l])
    end
    return g
end

"""
    network_layout(c; algorithm=:stress, seed=1)

2-D coordinates for every bus, from the network topology alone (pglib cases carry
no geographic data). `:stress` reads better on meshed transmission graphs;
`:spring` is the force-directed alternative.
"""
function network_layout(c; algorithm::Symbol=:stress, seed::Int=1)
    adj = adjacency_matrix(_simple_graph(c))
    pos = algorithm === :spring ? Spring(; seed=seed)(adj) : Stress()(adj)
    return [Point2f(p[1], p[2]) for p in pos]
end

# Line endpoints as a flat point vector for linesegments!.
_segments(pos, from, to, lines) =
    [p for l in lines for p in (pos[from[l]], pos[to[l]])]

_midpoints(pos, from, to, lines) =
    [Point2f((pos[from[l]][1] + pos[to[l]][1]) / 2,
             (pos[from[l]][2] + pos[to[l]][2]) / 2) for l in lines]

# Quadratic-Bezier arc from p1 to p2, bulging perpendicular to the chord by
# `curvature` (a fraction of the chord length). curvature = 0 gives the straight
# segment. Used to FAN OUT parallel connections on the reduced panel: contracting
# endpoints routinely leaves several external lines joining the same super-node
# pair, and drawn straight they lie exactly on top of one another, so the panel
# silently under-reports how many lines survive.
function _arc(p1, p2, curvature::Real; n::Int=48)
    abs(curvature) < 1e-9 && return [p1, p2]
    dx = p2[1] - p1[1]
    dy = p2[2] - p1[2]
    cx = (p1[1] + p2[1]) / 2 - curvature * dy
    cy = (p1[2] + p2[2]) / 2 + curvature * dx
    return [Point2f((1 - t)^2 * p1[1] + 2 * (1 - t) * t * cx + t^2 * p2[1],
                    (1 - t)^2 * p1[2] + 2 * (1 - t) * t * cy + t^2 * p2[2])
            for t in range(0, 1; length=n)]
end

_fan(k::Int) = k == 1 ? [0.0] : collect(range(-0.20, 0.20; length=k))

# Readable text over a busy network.
#
# Makie's text `strokewidth` centres the outline ON the glyph path, so half of it
# lands INSIDE the letter. A halo wide enough to be visible over coloured cluster
# webs (3-4 px) is then wider than the stems of a 14-17 px bold glyph, and the
# label renders as an unreadable white smear -- which is exactly what the
# L-numbers and the parallel-line counts were doing. Drawing the halo as copies
# of the text offset in a ring BEHIND it leaves the glyph itself untouched, so
# the halo can be as wide as it needs to be without eating the letter.
function _halo_text!(ax, positions; text, fontsize, color, font=:bold,
                     align=(:center, :bottom), offset=(0, 0),
                     halo::Real=1.7, halo_color=:white)
    isempty(positions) && return nothing
    for dx in (-halo, 0.0, halo), dy in (-halo, 0.0, halo)
        (dx == 0.0 && dy == 0.0) && continue
        text!(ax, positions; text=text, fontsize=fontsize, font=font,
              color=halo_color, align=align,
              offset=(offset[1] + dx, offset[2] + dy))
    end
    text!(ax, positions; text=text, fontsize=fontsize, font=font,
          color=color, align=align, offset=offset)
    return nothing
end

# One distinguishable colour per cluster (not just per adjacent pair), largest
# cluster first so the most visually prominent regions get the most-separated
# hues. Seeded with white/black/congestion-red so a cluster colour is never
# confusable with the background, text, or a congestion marker; those three
# come back as the first entries of distinguishable_colors and are dropped.
function _cluster_colors(merged_reps, csize)
    isempty(merged_reps) && return Dict{Int,RGB{Float64}}()
    seed = [RGB(1.0, 1.0, 1.0), RGB(0.0, 0.0, 0.0), convert(RGB, CONGESTED_COLOR)]
    palette = distinguishable_colors(length(merged_reps) + length(seed), seed)[length(seed)+1:end]
    ordered = sort(merged_reps; by = r -> (-csize[r], r))
    return Dict(i => palette[k] for (k, i) in enumerate(ordered))
end

"""
    plot_reduction(c, A; path=nothing, kwargs...)

Two-panel before/after figure on a shared layout.

**Left (original).** Every bus and every line. Each merged cluster is drawn as a
soft coloured web tracing its own internal lines -- those are the edges being
collapsed, so the coloured regions are exactly the clusters. External lines,
which are what actually survives the reduction, are drawn dark and thick on top.
Unmerged single buses stay neutral grey.

**Right (reduced).** Retained buses only, marker area growing with cluster size
(bounded, so one huge cluster cannot swallow the panel), and external lines only.
Parallel external lines joining the same pair of super-nodes are fanned into
separate arcs and the bundle is labelled with its line count, so the surviving
line count is visible rather than merely asserted in the panel title. Same
colours, same coordinates as the left panel.

Protected lines -- those the model forbids from ever becoming internal -- are
stroked in near-black on both panels, so they must be visible on the right. Pass
the model's own protected set via `binding_lines`; with no such argument the set
falls back to lines at their rating in `c.fhat`, which is only meaningful for a
single-scenario case.

Congested lines are stroked in reserved red on top of everything and labelled by
line number.

Keywords: `path` (save location; `nothing` shows/returns only), `title`,
`algorithm`/`seed` (layout), `pos` (reuse a layout), `resolution`,
`node_size`, `label_retained` (annotate super-node bus numbers),
`congested_lines` (line indices highlighted in red on both panels),
`congested_line_labels` (optional line-index-to-label dictionary),
`binding_lines` (protected/forced-external line indices; defaults to the
`|f̂| = f̄` rule), and `max_congested_labels` (label cap before the midpoint
annotations are dropped as unreadable clutter).
"""
function plot_reduction(c, A;
                        path=nothing,
                        title::AbstractString="",
                        algorithm::Symbol=:stress,
                        seed::Int=1,
                        pos=nothing,
                        resolution=(1700, 950),
                        node_size::Real=9,
                        label_retained::Bool=false,
                        congested_lines=Int[],
                        congested_line_labels=nothing,
                        label_congested::Bool=true,
                        binding_lines=nothing,
                        max_congested_labels::Int=30)
    Aint = round.(Int, A)
    retained, rep = _reduction_maps(Aint)
    pos === nothing && (pos = network_layout(c; algorithm=algorithm, seed=seed))

    csize = zeros(Int, c.N)                        # cluster size seen by each bus
    for j in 1:c.N
        csize[rep[j]] += 1
    end
    bus_csize = [csize[rep[j]] for j in 1:c.N]

    internal = [l for l in 1:c.Ln if rep[c.Efrom[l]] == rep[c.Eto[l]]]
    external = setdiff(1:c.Ln, internal)

    congested = sort!(unique!(Int.(collect(congested_lines))))
    all((1 .<= congested) .& (congested .<= c.Ln)) ||
        error("congested_lines must lie in 1:$(c.Ln)")
    congested_set = Set(congested)
    congested_ext = intersect(congested, external)
    congested_internal = intersect(congested, internal)
    congested_label(l) = isnothing(congested_line_labels) ? "L$l" :
        get(congested_line_labels, l, "L$l")

    # The model's own protected set when supplied. The |fhat| = frate fallback is
    # only correct for a single-scenario case: a multi-scenario run protects a
    # line that binds in ANY studied hour, which no single fhat snapshot knows.
    binding = if isnothing(binding_lines)
        [l for l in 1:c.Ln
         if abs(abs(c.fhat[l]) - c.frate[l]) <= 1e-6 * max(c.frate[l], 1.0)]
    else
        b = sort!(unique!(Int.(collect(binding_lines))))
        all((1 .<= b) .& (b .<= c.Ln)) ||
            error("binding_lines must lie in 1:$(c.Ln)")
        b
    end
    binding_set = Set(binding)
    binding_ext = intersect(binding, external)

    merged_reps = [i for i in retained if csize[i] > 1]
    cluster_color = _cluster_colors(merged_reps, csize)
    color_of(i) = get(cluster_color, i, SINGLETON_COLOR)

    internal_by_cluster = Dict{Int,Vector{Int}}()
    for l in internal
        push!(get!(internal_by_cluster, rep[c.Efrom[l]], Int[]), l)
    end

    merged = [j for j in 1:c.N if bus_csize[j] > 1]
    single = [j for j in 1:c.N if bus_csize[j] == 1]

    fig = Figure(size=resolution, backgroundcolor=SURFACE)
    isempty(title) || Label(fig[0, 1:2], title;
                            fontsize=21, font=:bold, color=INK_PRIMARY)

    axL = Axis(fig[1, 1], aspect=DataAspect(), backgroundcolor=SURFACE,
               titlesize=15, titlecolor=INK_SECOND,
               title="ORIGINAL  ·  $(c.N) buses, $(c.Ln) lines\n" *
                     "$(length(internal)) lines collapse into $(length(merged_reps)) clusters  ·  " *
                     "$(length(external)) lines survive")
    # Contracting endpoints can leave several external lines joining the SAME pair
    # of super-nodes. They are all still distinct lines with their own ratings and
    # flows -- the model never merges them -- so they are fanned into separate
    # arcs on the right panel and the bundle carries its own count.
    pair_lines = Dict{Tuple{Int,Int},Vector{Int}}()
    for l in external
        push!(get!(pair_lines, minmax(rep[c.Efrom[l]], rep[c.Eto[l]]), Int[]), l)
    end
    distinct_pairs = length(pair_lines)

    axR = Axis(fig[1, 2], aspect=DataAspect(), backgroundcolor=SURFACE,
               titlesize=15, titlecolor=INK_SECOND,
               title="REDUCED  ·  $(length(retained)) buses, $(length(external)) lines " *
                     "over $(distinct_pairs) connections\n" *
                     "largest cluster $(maximum(csize[retained])) buses  ·  " *
                     "$(length(congested_ext)) congested, $(length(binding_ext)) protected")
    # Shared coordinates are only half the job: each Axis otherwise autoscales to
    # its OWN data, and the right panel holds a subset of the nodes, so it would
    # render at a different zoom and defeat the comparison. Pin both to the full
    # layout extent.
    xs = [p[1] for p in pos]; ys = [p[2] for p in pos]
    padx = 0.06 * max(maximum(xs) - minimum(xs), eps())
    pady = 0.06 * max(maximum(ys) - minimum(ys), eps())
    for ax in (axL, axR)
        hidedecorations!(ax); hidespines!(ax)
        limits!(ax, minimum(xs) - padx, maximum(xs) + padx,
                    minimum(ys) - pady, maximum(ys) + pady)
    end

    # ---- left panel ----
    # Soft cluster webs first (largest cluster underneath, so a small cluster
    # sitting inside a sprawling one still reads on top).
    for i in sort(merged_reps; by = r -> -csize[r])
        ls = get(internal_by_cluster, i, Int[])
        isempty(ls) && continue
        segs = _segments(pos, c.Efrom, c.Eto, ls)
        linesegments!(axL, segs; color=(color_of(i), 0.16), linewidth=11)
        linesegments!(axL, segs; color=(color_of(i), 0.90), linewidth=1.4)
    end
    # Surviving lines are the point of the figure: dark and thick, over the webs.
    isempty(external) ||
        linesegments!(axL, _segments(pos, c.Efrom, c.Eto, external);
                      color=INK_SECOND, linewidth=2.4)
    isempty(binding) ||
        linesegments!(axL, _segments(pos, c.Efrom, c.Eto, binding);
                      color=INK_PRIMARY, linewidth=3.4)
    # A white halo keeps congested corridors legible wherever they cross a
    # cluster web; the L-number label means the state never rests on colour alone.
    isempty(congested) || begin
        segs = _segments(pos, c.Efrom, c.Eto, congested)
        linesegments!(axL, segs; color=(:white, 0.95), linewidth=7.5)
        linesegments!(axL, segs; color=CONGESTED_COLOR, linewidth=4.0)
    end

    isempty(single) || scatter!(axL, pos[single];
                                color=SINGLETON_COLOR, markersize=node_size * 0.6,
                                strokewidth=0)
    isempty(merged) || scatter!(axL, pos[merged];
                                color=[color_of(rep[j]) for j in merged],
                                markersize=node_size * 0.85,
                                strokewidth=0.4, strokecolor=SURFACE)
    # Super-nodes must be findable on the LEFT too: the bus each cluster collapses
    # onto is drawn as a larger dark-ringed diamond, so the right panel's markers
    # can be traced back to their originals.
    scatter!(axL, pos[retained];
             color=[color_of(i) for i in retained], marker=:diamond,
             markersize=[node_size * (csize[i] > 1 ? 1.7 : 1.1) for i in retained],
             strokewidth=1.1, strokecolor=INK_PRIMARY)

    # ---- right panel: super-nodes + surviving lines, same coordinates ----
    # Build every arc first, then draw by priority, so a later bundle's white
    # halo can never paint over an already-drawn congested line.
    plain_arcs = Vector{Vector{Point2f}}()
    protected_arcs = Vector{Vector{Point2f}}()
    congested_arcs = Vector{Vector{Point2f}}()
    bundle_labels = Point2f[]
    bundle_text = String[]
    for (a, b) in sort(collect(keys(pair_lines)))
        ls = sort(pair_lines[(a, b)])
        curvatures = _fan(length(ls))
        for (t, l) in enumerate(ls)
            arc = _arc(pos[a], pos[b], curvatures[t])
            if l in congested_set
                push!(congested_arcs, arc)
            elseif l in binding_set
                push!(protected_arcs, arc)
            else
                push!(plain_arcs, arc)
            end
        end
        if length(ls) > 1
            # Anchor the count just OUTSIDE the widest arc of its own fan, or it
            # lands in the middle of the bundle it is describing. Put it at the
            # QUARTER point rather than the midpoint: a congested line in the
            # same bundle already claims the midpoint for its L-number, and the
            # two collide there.
            outer = _arc(pos[a], pos[b], maximum(curvatures) + 0.12)
            push!(bundle_labels, outer[cld(length(outer), 4)])
            push!(bundle_text, "×$(length(ls))")
        end
    end
    for arc in plain_arcs
        lines!(axR, arc; color=INK_SECOND, linewidth=1.8)
    end
    for arc in protected_arcs
        lines!(axR, arc; color=INK_PRIMARY, linewidth=3.0)
    end
    for arc in congested_arcs
        lines!(axR, arc; color=(:white, 0.95), linewidth=7.5)
        lines!(axR, arc; color=CONGESTED_COLOR, linewidth=4.0)
    end
    _halo_text!(axR, bundle_labels; text=bundle_text, fontsize=15,
                color=INK_PRIMARY, align=(:center, :center))

    # Marker area grows with cluster size but is bounded: scaling by sqrt of the
    # RELATIVE size keeps a 185-of-200 cluster from covering the panel.
    maxc = maximum(csize[retained])
    scatter!(axR, pos[retained];
             color=[color_of(i) for i in retained], marker=:diamond,
             markersize=[node_size * (1.0 + 2.2 * sqrt(csize[i] / maxc))
                         for i in retained],
             strokewidth=1.2, strokecolor=INK_PRIMARY)

    if label_retained
        _halo_text!(axR, pos[retained]; text=string.(retained), fontsize=12,
                    color=INK_PRIMARY, align=(:center, :bottom), offset=(0, 8))
    end

    if label_congested && !isempty(congested) && length(congested) <= max_congested_labels
        _halo_text!(axL, _midpoints(pos, c.Efrom, c.Eto, congested);
                    text=["L$l" for l in congested], fontsize=17,
                    color=CONGESTED_COLOR, align=(:center, :bottom), offset=(0, 11))
        if !isempty(congested_ext)
            reduced_from = [rep[u] for u in c.Efrom]
            reduced_to = [rep[v] for v in c.Eto]
            _halo_text!(axR, _midpoints(pos, reduced_from, reduced_to, congested_ext);
                        text=["L$l" for l in congested_ext], fontsize=17,
                        color=CONGESTED_COLOR, align=(:center, :bottom), offset=(0, 11))
        end
    end

    legend_elements = Any[
        PolyElement(color=(LEGEND_SWATCH, 0.16),
                    strokecolor=LEGEND_SWATCH, strokewidth=1.2),
        MarkerElement(marker=:diamond, color=LEGEND_SWATCH,
                      strokecolor=INK_PRIMARY, strokewidth=1.2, markersize=15),
        MarkerElement(marker=:circle, color=SINGLETON_COLOR,
                      strokecolor=SURFACE, strokewidth=0.5, markersize=9),
        LineElement(color=INK_SECOND, linewidth=2.4),
    ]
    legend_labels = String[
        "merged cluster (its own internal lines)",
        "super-node (both panels; area ∝ cluster size)",
        "unmerged single bus",
        "external line — survives the reduction",
    ]
    if !isempty(binding)
        push!(legend_elements, LineElement(color=INK_PRIMARY, linewidth=3.4))
        push!(legend_labels, "protected — forced to stay external")
    end
    if !isempty(congested)
        push!(legend_elements, LineElement(color=CONGESTED_COLOR, linewidth=4.0))
        push!(legend_labels, "congested line (labelled by number)")
    end
    if any(length(v) > 1 for v in values(pair_lines))
        push!(legend_elements,
              MarkerElement(marker=:hline, color=INK_SECOND, markersize=15))
        push!(legend_labels, "×n — parallel lines fanned on the reduced panel")
    end
    Legend(fig[2, 1:2], legend_elements, legend_labels;
           orientation=:horizontal, nbanks=2, framevisible=true,
           framecolor=colorant"#e1e0d9", labelsize=12.5, labelcolor=INK_SECOND,
           patchsize=(26, 12), colgap=18, padding=(12, 12, 8, 8))

    if !isempty(congested)
        details = [congested_label(l) for l in congested]
        detail_rows = [join(details[i:min(i + 3, end)], "   ·   ")
                       for i in 1:4:length(details)]
        footer = "CONGESTED   " * join(detail_rows, "\nCONGESTED   ")
        color = CONGESTED_COLOR
        if !isempty(congested_internal)
            footer *= "\nWARNING: $(length(congested_internal)) congested line(s) became internal"
        end
        Label(fig[3, 1:2], footer; fontsize=13, font=:bold, color=color)
    end

    rowgap!(fig.layout, 6)

    if path !== nothing
        mkpath(dirname(abspath(path)))
        save(path, fig)
        println("wrote ", path)
    end
    return fig
end

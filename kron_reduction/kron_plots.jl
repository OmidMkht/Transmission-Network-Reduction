# --------------------------------------------------------------------------- #
# Three-panel Kron figure. Loaded lazily by kron_reporting.jl's
# ensure_kron_plots (mirrors ../transmission_plots.jl's own laziness -- never
# pay for CairoMakie on a plots-off run). Reuses transmission_plots.jl's
# low-level helpers (network_layout, _reduction_maps, _cluster_colors,
# _segments, _arc, _fan, _halo_text!, _midpoints, colour constants) via
# TxReport.* qualification rather than redefining them -- nothing in the
# parent directory is modified by this file.
#
# Left, middle and right panels share ONE layout, computed once on the full
# network, so all three are comparable by eye:
#
#   LEFT   full network, every Kron-eligible chain drawn as a coloured web
#          (one colour per chain) over its own interior buses/lines -- what
#          gets collapsed before the reduction MILP ever runs.
#   MIDDLE the "boundary" network the MILP actually solves on, in the SAME
#          visual language as plot_reduction's own left panel: a soft
#          coloured web per MILP cluster, dark external survivors, grey
#          singles. Coloured by final MILP cluster (NOT by Kron chain --
#          that's the left panel's job), using the SAME colours as the right
#          panel's super-nodes. Equivalent lines that survive as external are
#          dashed, so you can still see which boundary lines stand in for a
#          collapsed chain.
#   RIGHT  the final reduced network -- retained super-nodes and external
#          lines only, from the unfolded full-size assignment matrix. Same
#          rendering as transmission_plots.jl's plot_reduction right panel.
# --------------------------------------------------------------------------- #

using CairoMakie
using Colors
using GeometryBasics: Point2f

"""
    plot_kron_reduction(c_full, kron_map, c_boundary, A_full; path=nothing, kwargs...)

`c_full`/`c_boundary` are `MultiScenarioTxReductionCase`s, `A_full` the
unfolded full-size assignment matrix (`kron_reduction/kron_unfold.jl`'s
`unfold_kron_assignment`). Keyword arguments mirror `plot_reduction`'s:
`path`, `title`, `algorithm`/`seed`/`pos` (layout), `resolution`, `node_size`,
`label_retained`, `congested_lines`, `congested_line_labels`,
`label_congested`, `binding_lines`, `max_congested_labels`.
"""
function plot_kron_reduction(c_full, kron_map, c_boundary, A_full;
                             path=nothing,
                             title::AbstractString="",
                             algorithm::Symbol=:stress,
                             seed::Int=1,
                             pos=nothing,
                             resolution=(2500, 950),
                             node_size::Real=9,
                             label_retained::Bool=false,
                             congested_lines=Int[],
                             congested_line_labels=nothing,
                             label_congested::Bool=true,
                             binding_lines=nothing,
                             max_congested_labels::Int=30)
    c = c_full.base
    Aint = round.(Int, A_full)
    retained, rep = TxReport._reduction_maps(Aint)
    pos === nothing && (pos = TxReport.network_layout(c; algorithm=algorithm, seed=seed))

    csize = zeros(Int, c.N)
    for j in 1:c.N
        csize[rep[j]] += 1
    end

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

    binding = if isnothing(binding_lines)
        Int[]
    else
        b = sort!(unique!(Int.(collect(binding_lines))))
        all((1 .<= b) .& (b .<= c.Ln)) ||
            error("binding_lines must lie in 1:$(c.Ln)")
        b
    end
    binding_set = Set(binding)
    binding_ext = intersect(binding, external)

    # ---- one colour per Kron chain (independent of the final-cluster palette
    # used on the right -- chains and final clusters are different concepts) ----
    chains = kron_map.chains
    nchains = length(chains)
    chain_color = Dict{Int,RGB{Float64}}()
    if nchains > 0
        seed_colors = [RGB(1.0, 1.0, 1.0), RGB(0.0, 0.0, 0.0), convert(RGB, TxReport.CONGESTED_COLOR)]
        palette = distinguishable_colors(nchains + length(seed_colors), seed_colors)[length(seed_colors)+1:end]
        order = sortperm(chains; by = ch -> -length(ch.interior))
        for (k, ci) in enumerate(order)
            chain_color[ci] = palette[k]
        end
    end
    color_of_chain(ci) = get(chain_color, ci, TxReport.SINGLETON_COLOR)

    merged_reps = [i for i in retained if csize[i] > 1]
    cluster_color = TxReport._cluster_colors(merged_reps, csize)
    color_of_cluster(i) = get(cluster_color, i, TxReport.SINGLETON_COLOR)

    fig = Figure(size=resolution, backgroundcolor=TxReport.SURFACE)
    isempty(title) || Label(fig[0, 1:3], title;
                            fontsize=21, font=:bold, color=TxReport.INK_PRIMARY)

    xs = [p[1] for p in pos]; ys = [p[2] for p in pos]
    padx = 0.06 * max(maximum(xs) - minimum(xs), eps())
    pady = 0.06 * max(maximum(ys) - minimum(ys), eps())

    # ============================= LEFT: full network, chains highlighted ============================= #
    axL = Axis(fig[1, 1], aspect=DataAspect(), backgroundcolor=TxReport.SURFACE,
               titlesize=15, titlecolor=TxReport.INK_SECOND,
               title="$(nchains) Kron-eligible chain$(nchains == 1 ? "" : "s")  ·  " *
                     "$(kron_map.full_N - kron_map.boundary_N) buses eliminable")

    # Neutral backdrop first, so every original line is visible even where no
    # chain touches it.
    c.Ln > 0 && linesegments!(
        axL, TxReport._segments(pos, c.Efrom, c.Eto, 1:c.Ln);
        color=(TxReport.INK_MUTED, 0.45), linewidth=1.2)
    for (ci, ch) in enumerate(chains)
        segs = TxReport._segments(pos, c.Efrom, c.Eto, ch.chain_lines)
        col = color_of_chain(ci)
        linesegments!(axL, segs; color=(col, 0.20), linewidth=11)
        linesegments!(axL, segs; color=(col, 0.95), linewidth=1.8)
    end
    isempty(binding) ||
        linesegments!(axL, TxReport._segments(pos, c.Efrom, c.Eto, binding);
                      color=TxReport.INK_PRIMARY, linewidth=3.4)
    isempty(congested) || begin
        segs = TxReport._segments(pos, c.Efrom, c.Eto, congested)
        linesegments!(axL, segs; color=(:white, 0.95), linewidth=7.5)
        linesegments!(axL, segs; color=TxReport.CONGESTED_COLOR, linewidth=4.0)
    end

    chain_interior_set = Set(b for ch in chains for b in ch.interior)
    other_buses = [j for j in 1:c.N if !(j in chain_interior_set)]
    isempty(other_buses) || scatter!(axL, pos[other_buses];
                                     color=TxReport.SINGLETON_COLOR, markersize=node_size * 0.6,
                                     strokewidth=0)
    for (ci, ch) in enumerate(chains)
        isempty(ch.interior) && continue
        scatter!(axL, pos[ch.interior]; color=color_of_chain(ci),
                 markersize=node_size * 0.85, strokewidth=0.4, strokecolor=TxReport.SURFACE)
    end
    anchor_buses = unique(vcat([ch.x for ch in chains], [ch.y for ch in chains]))
    isempty(anchor_buses) || scatter!(axL, pos[anchor_buses];
                                      color=TxReport.SURFACE, markersize=node_size * 0.85,
                                      strokewidth=1.6, strokecolor=TxReport.INK_PRIMARY)

    if label_congested && !isempty(congested) && length(congested) <= max_congested_labels
        TxReport._halo_text!(axL, TxReport._midpoints(pos, c.Efrom, c.Eto, congested);
                    text=["L$l" for l in congested], fontsize=17,
                    color=TxReport.CONGESTED_COLOR, align=(:center, :bottom), offset=(0, 11))
    end

    # ============================ MIDDLE: Kron-reduced (boundary) network, MILP clusters shown ============================ #
    # Same visual language as plot_reduction's own left panel (soft coloured
    # webs for merged clusters, dark external survivors, grey singles) --
    # restricted to the boundary network, and using the SAME cluster colours
    # as the right panel's super-nodes (color_of_cluster, keyed by the
    # full-index representative -- identical for a boundary bus whether you
    # look it up here or on the right, since both come from the one A_full).
    # This is deliberately NOT coloured by Kron chain -- that's the left
    # panel's job; this panel's job is showing what the MILP itself grouped.
    cb = c_boundary.base
    boundary_pos = [pos[kron_map.full_of[bb]] for bb in 1:kron_map.boundary_N]
    boundary_rep = [rep[kron_map.full_of[bb]] for bb in 1:kron_map.boundary_N]
    boundary_csize = [csize[r] for r in boundary_rep]
    chain_of_boundary_line = Dict(ch.boundary_line => ci for (ci, ch) in enumerate(chains))
    equiv_b = collect(keys(chain_of_boundary_line))

    b_internal = [l for l in 1:cb.Ln if boundary_rep[cb.Efrom[l]] == boundary_rep[cb.Eto[l]]]
    b_external = setdiff(1:cb.Ln, b_internal)
    b_equiv_external = intersect(equiv_b, b_external)
    b_plain_external = setdiff(b_external, b_equiv_external)
    b_merged_reps = unique(boundary_rep[bb] for bb in 1:kron_map.boundary_N if boundary_csize[bb] > 1)

    axM = Axis(fig[1, 2], aspect=DataAspect(), backgroundcolor=TxReport.SURFACE,
               titlesize=15, titlecolor=TxReport.INK_SECOND,
               title="$(length(b_internal)) lines collapse into $(length(b_merged_reps)) clusters  ·  " *
                     "$(length(b_external)) lines survive  ($(length(equiv_b)) equivalent)")

    internal_by_cluster_b = Dict{Int,Vector{Int}}()
    for l in b_internal
        push!(get!(internal_by_cluster_b, boundary_rep[cb.Efrom[l]], Int[]), l)
    end
    for i in sort(b_merged_reps; by = r -> -csize[r])
        ls = get(internal_by_cluster_b, i, Int[])
        isempty(ls) && continue
        segs = TxReport._segments(boundary_pos, cb.Efrom, cb.Eto, ls)
        linesegments!(axM, segs; color=(color_of_cluster(i), 0.18), linewidth=11)
        linesegments!(axM, segs; color=(color_of_cluster(i), 0.90), linewidth=1.4)
    end
    isempty(b_plain_external) ||
        linesegments!(axM, TxReport._segments(boundary_pos, cb.Efrom, cb.Eto, b_plain_external);
                      color=TxReport.INK_SECOND, linewidth=2.4)
    isempty(b_equiv_external) ||
        linesegments!(axM, TxReport._segments(boundary_pos, cb.Efrom, cb.Eto, b_equiv_external);
                      color=TxReport.INK_SECOND, linewidth=2.4, linestyle=:dash)

    # Protected/congested overlay, same as the left panel -- both are provably
    # never chain lines (chain eligibility already excludes near-congested
    # lines), so boundary_line_of_full is always nonzero here; the filter is
    # just defensive.
    binding_b = [kron_map.boundary_line_of_full[l] for l in binding
                 if kron_map.boundary_line_of_full[l] != 0]
    isempty(binding_b) ||
        linesegments!(axM, TxReport._segments(boundary_pos, cb.Efrom, cb.Eto, binding_b);
                      color=TxReport.INK_PRIMARY, linewidth=3.4)
    congested_b = Int[]
    congested_b_labels = String[]
    for l in congested
        bl = kron_map.boundary_line_of_full[l]
        bl == 0 && continue
        push!(congested_b, bl)
        push!(congested_b_labels, "L$l")
    end
    isempty(congested_b) || begin
        segs = TxReport._segments(boundary_pos, cb.Efrom, cb.Eto, congested_b)
        linesegments!(axM, segs; color=(:white, 0.95), linewidth=7.5)
        linesegments!(axM, segs; color=TxReport.CONGESTED_COLOR, linewidth=4.0)
    end
    if label_congested && !isempty(congested_b) && length(congested) <= max_congested_labels
        TxReport._halo_text!(axM, TxReport._midpoints(boundary_pos, cb.Efrom, cb.Eto, congested_b);
                    text=congested_b_labels, fontsize=17,
                    color=TxReport.CONGESTED_COLOR, align=(:center, :bottom), offset=(0, 11))
    end

    # Interior-bus count for each equivalent line that survived as its own
    # external line -- chains folded into a cluster's web have no individual
    # label, exactly like any other internal line in this scheme.
    eq_labels = Point2f[]
    eq_text = String[]
    for l in b_equiv_external
        ci = chain_of_boundary_line[l]
        length(chains[ci].interior) <= 1 && continue
        push!(eq_labels, TxReport._midpoints(boundary_pos, cb.Efrom, cb.Eto, [l])[1])
        push!(eq_text, "×$(length(chains[ci].interior))")
    end
    TxReport._halo_text!(axM, eq_labels; text=eq_text, fontsize=12,
                         color=TxReport.INK_SECOND, align=(:center, :center))

    boundary_single = [bb for bb in 1:kron_map.boundary_N if boundary_csize[bb] == 1]
    boundary_merged = [bb for bb in 1:kron_map.boundary_N if boundary_csize[bb] > 1]
    isempty(boundary_single) || scatter!(axM, boundary_pos[boundary_single];
                                         color=TxReport.SINGLETON_COLOR, markersize=node_size * 0.6,
                                         strokewidth=0)
    isempty(boundary_merged) || scatter!(axM, boundary_pos[boundary_merged];
                                         color=[color_of_cluster(boundary_rep[bb]) for bb in boundary_merged],
                                         markersize=node_size * 0.85,
                                         strokewidth=0.4, strokecolor=TxReport.SURFACE)

    # ================================ RIGHT: final reduced network ================================ #
    axR = Axis(fig[1, 3], aspect=DataAspect(), backgroundcolor=TxReport.SURFACE,
               titlesize=15, titlecolor=TxReport.INK_SECOND,
               title="largest cluster $(maximum(csize[retained])) buses  ·  " *
                     "$(length(congested_ext)) congested, $(length(binding_ext)) protected")

    pair_lines = Dict{Tuple{Int,Int},Vector{Int}}()
    for l in external
        push!(get!(pair_lines, minmax(rep[c.Efrom[l]], rep[c.Eto[l]]), Int[]), l)
    end

    plain_arcs = Vector{Vector{Point2f}}()
    protected_arcs = Vector{Vector{Point2f}}()
    congested_arcs = Vector{Vector{Point2f}}()
    bundle_labels = Point2f[]
    bundle_text = String[]
    for (a, b) in sort(collect(keys(pair_lines)))
        ls = sort(pair_lines[(a, b)])
        curvatures = TxReport._fan(length(ls))
        for (t, l) in enumerate(ls)
            arc = TxReport._arc(pos[a], pos[b], curvatures[t])
            if l in congested_set
                push!(congested_arcs, arc)
            elseif l in binding_set
                push!(protected_arcs, arc)
            else
                push!(plain_arcs, arc)
            end
        end
        if length(ls) > 1
            outer = TxReport._arc(pos[a], pos[b], maximum(curvatures) + 0.12)
            push!(bundle_labels, outer[cld(length(outer), 4)])
            push!(bundle_text, "×$(length(ls))")
        end
    end
    for arc in plain_arcs
        lines!(axR, arc; color=TxReport.INK_SECOND, linewidth=1.8)
    end
    for arc in protected_arcs
        lines!(axR, arc; color=TxReport.INK_PRIMARY, linewidth=3.0)
    end
    for arc in congested_arcs
        lines!(axR, arc; color=(:white, 0.95), linewidth=7.5)
        lines!(axR, arc; color=TxReport.CONGESTED_COLOR, linewidth=4.0)
    end
    TxReport._halo_text!(axR, bundle_labels; text=bundle_text, fontsize=15,
                color=TxReport.INK_PRIMARY, align=(:center, :center))

    maxc = maximum(csize[retained])
    scatter!(axR, pos[retained];
             color=[color_of_cluster(i) for i in retained], marker=:diamond,
             markersize=[node_size * (1.0 + 2.2 * sqrt(csize[i] / maxc))
                         for i in retained],
             strokewidth=1.2, strokecolor=TxReport.INK_PRIMARY)

    if label_retained
        TxReport._halo_text!(axR, pos[retained]; text=string.(retained), fontsize=12,
                    color=TxReport.INK_PRIMARY, align=(:center, :bottom), offset=(0, 8))
    end

    if label_congested && !isempty(congested_ext) && length(congested) <= max_congested_labels
        reduced_from = [rep[u] for u in c.Efrom]
        reduced_to = [rep[v] for v in c.Eto]
        TxReport._halo_text!(axR, TxReport._midpoints(pos, reduced_from, reduced_to, congested_ext);
                    text=["L$l" for l in congested_ext], fontsize=17,
                    color=TxReport.CONGESTED_COLOR, align=(:center, :bottom), offset=(0, 11))
    end

    for ax in (axL, axM, axR)
        hidedecorations!(ax); hidespines!(ax)
        limits!(ax, minimum(xs) - padx, maximum(xs) + padx,
                    minimum(ys) - pady, maximum(ys) + pady)
    end

    # ---------------------------------- legend + footer ---------------------------------- #
    legend_elements = Any[
        PolyElement(color=(TxReport.LEGEND_SWATCH, 0.20),
                    strokecolor=TxReport.LEGEND_SWATCH, strokewidth=1.2),
        MarkerElement(marker=:circle, color=TxReport.SURFACE,
                      strokecolor=TxReport.INK_PRIMARY, strokewidth=1.6, markersize=9),
        PolyElement(color=(TxReport.LEGEND_SWATCH, 0.18),
                    strokecolor=TxReport.LEGEND_SWATCH, strokewidth=1.2),
        LineElement(color=TxReport.INK_SECOND, linewidth=2.4, linestyle=:dash),
        LineElement(color=TxReport.INK_SECOND, linewidth=1.8),
        MarkerElement(marker=:diamond, color=TxReport.LEGEND_SWATCH,
                      strokecolor=TxReport.INK_PRIMARY, strokewidth=1.2, markersize=15),
        MarkerElement(marker=:circle, color=TxReport.SINGLETON_COLOR,
                      strokecolor=TxReport.SURFACE, strokewidth=0.5, markersize=9),
        LineElement(color=TxReport.INK_SECOND, linewidth=2.4),
    ]
    legend_labels = String[
        "Kron chain (left) — collapsible bus web, coloured per chain",
        "chain anchor bus (left) — kept in the boundary network",
        "MILP cluster (middle web / right super-node) — same colours both places",
        "equivalent line (middle, dashed) — stands in for one collapsed chain",
        "untouched line (middle)",
        "super-node (right; area ∝ cluster size)",
        "unmerged single bus",
        "external line — survives the reduction (right)",
    ]
    if !isempty(binding)
        push!(legend_elements, LineElement(color=TxReport.INK_PRIMARY, linewidth=3.4))
        push!(legend_labels, "protected — forced to stay external (right)")
    end
    if !isempty(congested)
        push!(legend_elements, LineElement(color=TxReport.CONGESTED_COLOR, linewidth=4.0))
        push!(legend_labels, "congested line (labelled by number)")
    end
    if any(length(v) > 1 for v in values(pair_lines))
        push!(legend_elements,
              MarkerElement(marker=:hline, color=TxReport.INK_SECOND, markersize=15))
        push!(legend_labels, "×n — parallel lines fanned out (right)")
    end
    Legend(fig[2, 1:3], legend_elements, legend_labels;
           orientation=:horizontal, nbanks=3, framevisible=true,
           framecolor=colorant"#e1e0d9", labelsize=12, labelcolor=TxReport.INK_SECOND,
           patchsize=(24, 12), colgap=16, padding=(12, 12, 8, 8))

    if !isempty(congested)
        details = [congested_label(l) for l in congested]
        detail_rows = [join(details[i:min(i + 3, end)], "   ·   ")
                       for i in 1:4:length(details)]
        footer = "CONGESTED   " * join(detail_rows, "\nCONGESTED   ")
        if !isempty(congested_internal)
            footer *= "\nWARNING: $(length(congested_internal)) congested line(s) became internal"
        end
        Label(fig[3, 1:3], footer; fontsize=13, font=:bold, color=TxReport.CONGESTED_COLOR)
    end

    rowgap!(fig.layout, 6)

    if path !== nothing
        mkpath(dirname(abspath(path)))
        save(path, fig)
        println("wrote ", path)
    end
    return fig
end

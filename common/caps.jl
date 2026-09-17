# Line budget, hop cap and cluster-size cap, shared by all three approaches.
#
#   hop cap   no chain of more than k merged lines inside a cluster
#   size cap  no cluster with more than k buses
#   budget    no more than k merged lines

module Caps

using JuMP
import MathOptInterface as MOI

export capped_paths, long_chain, size_witness, merged_path, cluster_sizes, hop_diameter,
       within_caps, trim_to_caps, add_caps!

# Every simple path whose length crosses `limit` (edges, or path weight when
# `weight` is given). Errors past max_paths rather than returning a partial list.
function capped_paths(base, limit::Float64; weight=nothing,
                      forbidden=falses(base.Ln), max_paths::Int=200_000)
    N, Ln = base.N, base.Ln
    u, v = base.Efrom, base.Eto
    inc = [Int[] for _ in 1:N]
    for l in 1:Ln
        forbidden[l] && continue
        u[l] == v[l] && continue
        push!(inc[u[l]], l)
        push!(inc[v[l]], l)
    end
    w = isnothing(weight) ? ones(Ln) : collect(weight)
    out = Vector{Vector{Int}}()
    seen = Set{Vector{Int}}()
    onpath = falses(N)
    path = Int[]
    function walk(b, acc)
        for l in inc[b]
            l in path && continue
            o = u[l] == b ? v[l] : u[l]
            onpath[o] && continue
            push!(path, l)
            acc2 = acc + w[l]
            if acc2 > limit
                key = sort(path)
                if !(key in seen)
                    push!(seen, key)
                    push!(out, copy(path))
                    length(out) <= max_paths || error(
                        "path cap exceeds max_paths=$max_paths; no exact cap " *
                        "was constructed. Use a tighter cap or a smaller graph.")
                end
            else
                onpath[o] = true
                walk(o, acc2)
                onpath[o] = false
            end
            pop!(path)
        end
    end
    for b in 1:N
        onpath[b] = true
        walk(b, 0.0)
        onpath[b] = false
    end
    return out
end

# A connected Kmax+1-bus subset of a cluster and the Kmax merged lines spanning
# it (BFS), or empty when the cluster is small enough.
function size_witness(u, v, clv, ls, Kmax::Int)
    adj = Dict{Int,Vector{Tuple{Int,Int}}}()
    for l in ls
        clv[l] > 0.5 || continue
        push!(get!(adj, u[l], Tuple{Int,Int}[]), (l, v[l]))
        push!(get!(adj, v[l], Tuple{Int,Int}[]), (l, u[l]))
    end
    isempty(adj) && return (Int[], Int[])
    src = first(keys(adj))
    seen = Set{Int}([src])
    W, T, queue = [src], Int[], [src]
    while !isempty(queue) && length(W) < Kmax + 1
        b = popfirst!(queue)
        for (l, o) in get(adj, b, Tuple{Int,Int}[])
            o in seen && continue
            push!(seen, o); push!(W, o); push!(T, l); push!(queue, o)
            length(W) >= Kmax + 1 && break
        end
    end
    return length(W) == Kmax + 1 ? (W, T) : (Int[], Int[])
end

# A merged simple path of cap+1 lines, or Int[] if none is found. BFS first
# (exact on trees), then a bounded DFS for clusters with cycles; `exhausted` is
# set if the DFS budget runs out.
function long_chain(u, v, clv, ls, cap::Int; max_steps::Int=1_000_000,
                    exhausted::Ref{Bool}=Ref(false))
    adj = Dict{Int,Vector{Tuple{Int,Int}}}()
    for l in ls
        clv[l] > 0.5 || continue
        u[l] == v[l] && continue
        push!(get!(adj, u[l], Tuple{Int,Int}[]), (l, v[l]))
        push!(get!(adj, v[l], Tuple{Int,Int}[]), (l, u[l]))
    end
    length(adj) >= cap + 2 || return Int[]
    for src in keys(adj)
        dist = Dict(src => 0)
        pred = Dict{Int,Tuple{Int,Int}}()
        queue = [src]; head = 1
        while head <= length(queue)
            b = queue[head]; head += 1
            for (l, o) in adj[b]
                haskey(dist, o) && continue
                dist[o] = dist[b] + 1
                pred[o] = (b, l)
                push!(queue, o)
                if dist[o] == cap + 1
                    path = Int[]; x = o
                    while x != src
                        p, pl = pred[x]; push!(path, pl); x = p
                    end
                    return reverse!(path)
                end
            end
        end
    end
    steps = 0
    path = Int[]
    onpath = Set{Int}()
    function dfs(b)
        length(path) > cap && return true
        for (l, o) in adj[b]
            o in onpath && continue
            (steps += 1) > max_steps && return false
            push!(path, l); push!(onpath, o)
            dfs(o) && return true
            pop!(path); delete!(onpath, o)
        end
        return false
    end
    for src in keys(adj)
        push!(onpath, src)
        dfs(src) && return copy(path)
        delete!(onpath, src)
        steps > max_steps && (exhausted[] = true; break)
    end
    return Int[]
end

# Merged lines joining buses a and b inside a cluster (BFS), or Int[] if none.
function merged_path(u, v, clv, ls, a::Int, b::Int)
    adj = Dict{Int,Vector{Tuple{Int,Int}}}()
    for l in ls
        clv[l] > 0.5 || continue
        push!(get!(adj, u[l], Tuple{Int,Int}[]), (l, v[l]))
        push!(get!(adj, v[l], Tuple{Int,Int}[]), (l, u[l]))
    end
    pred = Dict{Int,Tuple{Int,Int}}(a => (0, 0))
    queue = [a]; head = 1
    while head <= length(queue)
        x = queue[head]; head += 1
        x == b && break
        for (l, o) in get(adj, x, Tuple{Int,Int}[])
            haskey(pred, o) && continue
            pred[o] = (x, l); push!(queue, o)
        end
    end
    haskey(pred, b) && a != b || return Int[]
    path = Int[]; x = b
    while x != a
        p, l = pred[x]; push!(path, l); x = p
    end
    return path
end

# Cluster representative of every bus, from the merged-line mask.
function reps(base, internal)
    parent = collect(1:base.N)
    find(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    for l in 1:base.Ln
        internal[l] > 0.5 || continue
        a, b = find(base.Efrom[l]), find(base.Eto[l])
        a != b && (parent[a] = b)
    end
    return [find(b) for b in 1:base.N]
end

"Bus count per cluster, keyed by representative."
function cluster_sizes(base, internal)
    size = Dict{Int,Int}()
    for r in reps(base, internal)
        size[r] = get(size, r, 0) + 1
    end
    return size
end

"Longest BFS distance through merged lines inside any cluster."
function hop_diameter(base, internal)
    adj = [Int[] for _ in 1:base.N]
    for l in 1:base.Ln
        internal[l] > 0.5 || continue
        a, b = base.Efrom[l], base.Eto[l]
        a == b && continue
        push!(adj[a], b); push!(adj[b], a)
    end
    worst = 0
    for s in 1:base.N
        isempty(adj[s]) && continue
        dist = Dict(s => 0); q = [s]; h = 1
        while h <= length(q)
            x = q[h]; h += 1
            for y in adj[x]
                haskey(dist, y) && continue
                dist[y] = dist[x] + 1; push!(q, y)
            end
        end
        worst = max(worst, maximum(values(dist)))
    end
    return worst
end

"""
    within_caps(base, internal; budget, hop_cap, size_cap, bus=nothing)

True when the merged-line mask respects every cap. With `bus`, only the cluster
holding that bus is checked (enough after a single merge).
"""
function within_caps(base, internal; budget=nothing, hop_cap=nothing,
                     size_cap=nothing, bus=nothing)
    isnothing(budget) || count(>(0.5), internal) <= budget || return false
    isnothing(hop_cap) && isnothing(size_cap) && return true
    rep = reps(base, internal)
    target = isnothing(bus) ? nothing : rep[bus]
    if !isnothing(size_cap)
        n = isnothing(target) ? maximum(count(==(r), rep) for r in unique(rep)) :
                                count(==(target), rep)
        n <= size_cap || return false
    end
    if !isnothing(hop_cap)
        clv = Float64.(internal)
        groups = Dict{Int,Vector{Int}}()
        for l in 1:base.Ln
            internal[l] > 0.5 || continue
            r = rep[base.Efrom[l]]
            (isnothing(target) || r == target) || continue
            push!(get!(groups, r, Int[]), l)
        end
        for ls in values(groups)
            length(ls) > hop_cap || continue
            isempty(long_chain(base.Efrom, base.Eto, clv, ls, hop_cap)) || return false
        end
    end
    return true
end

"Drop lines from a start mask, in order, until it respects the caps."
function trim_to_caps(base, internal; budget=nothing, hop_cap=nothing, size_cap=nothing)
    kept = falses(base.Ln)
    for l in findall(>(0.5), internal)
        kept[l] = true
        within_caps(base, kept; budget, hop_cap, size_cap, bus=base.Efrom[l]) ||
            (kept[l] = false)
    end
    return kept
end

"""
    add_caps!(m, cl, base; budget, hop_cap, size_cap, protected, hop_mode=:auto)

Caps on a JuMP model with one merge binary per line. Hop rows are listed when
they fit (they tighten the LP) and cut lazily otherwise; the size cap is lazy.
With a budget, lines closing a loop inside a cluster are forced merged (lazily),
so the budget counts every line inside a cluster. Registers a lazy-constraint
callback, so the model must not already have one.
"""
function add_caps!(m, cl, base; budget=nothing, hop_cap=nothing, size_cap=nothing,
                   protected=falses(base.Ln), hop_mode::Symbol=:auto)
    hop_mode in (:auto, :static, :lazy) || error("hop_mode must be :auto, :static or :lazy")
    Ln, u, v = base.Ln, base.Efrom, base.Eto
    isnothing(budget) || @constraint(m, sum(cl[l] for l in 1:Ln) <= budget)
    hop_lazy, n_rows = false, 0
    if !isnothing(hop_cap)
        ps = hop_mode === :lazy ? nothing : try
            capped_paths(base, float(hop_cap); forbidden=protected)
        catch err
            hop_mode === :auto && err isa ErrorException &&
                occursin("max_paths", err.msg) || rethrow()
            nothing
        end
        if isnothing(ps)
            hop_lazy = true
        else
            for P in ps
                @constraint(m, sum(cl[l] for l in P) <= length(P) - 1)
            end
            n_rows = length(ps)
        end
        println(" Hop cap = ", hop_cap, hop_lazy ? "   lazy" : "   rows = $n_rows")
    end
    isnothing(size_cap) || println(" Cluster size cap = ", size_cap, "   lazy")
    isnothing(budget) || println(" Line budget = ", budget)
    if hop_lazy || !isnothing(size_cap) || !isnothing(budget)
        set_optimizer_attribute(m, "LazyConstraints", 1)
        function caps_callback(cb_data)
            callback_node_status(cb_data, m) == MOI.CALLBACK_NODE_STATUS_INTEGER || return
            clv = [callback_value(cb_data, cl[l]) for l in 1:Ln]
            rep = reps(base, clv)
            groups = Dict{Int,Vector{Int}}()
            for l in 1:Ln
                rep[u[l]] == rep[v[l]] || continue
                push!(get!(groups, rep[u[l]], Int[]), l)
            end
            for (r, ls) in groups
                if !isnothing(budget)
                    for l in ls
                        clv[l] > 0.5 && continue
                        P = merged_path(u, v, clv, ls, u[l], v[l])
                        isempty(P) || MOI.submit(m, MOI.LazyConstraint(cb_data),
                            @build_constraint(cl[l] >= sum(cl[k] for k in P) - (length(P) - 1)))
                    end
                end
                if hop_lazy && length(ls) > hop_cap
                    P = long_chain(u, v, clv, ls, hop_cap)
                    isempty(P) || MOI.submit(m, MOI.LazyConstraint(cb_data),
                        @build_constraint(sum(cl[l] for l in P) <= length(P) - 1))
                end
                if !isnothing(size_cap) && count(==(r), rep) > size_cap
                    _, T = size_witness(u, v, clv, ls, size_cap)
                    isempty(T) || MOI.submit(m, MOI.LazyConstraint(cb_data),
                        @build_constraint(sum(cl[l] for l in T) <= size_cap - 1))
                end
            end
        end
        MOI.set(m, MOI.LazyConstraintCallback(), caps_callback)
    end
    return (hop_lazy=hop_lazy, n_hop_rows=n_rows)
end

end # module

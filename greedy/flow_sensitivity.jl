# Flow-sensitivity ranking for greedy merges.
#
# Merging the ends of line (i,j) at a fixed injection moves every line flow by
#     df = t * (H[:,i] - H[:,j]),   t = -(theta_i - theta_j) / (X_ii + X_jj - 2X_ij)
# so scoring a candidate is O(L), no LP. After a merge X and H get exact rank-one
# updates (X' = X - (Xz)(z'X)/(z'Xz), H' = H - (Hz)(z'X)/(z'Xz)).
#
#     score = alpha * ||df ./ w||_2 + (1 - alpha) * ||df ./ w||_inf,  worst demand
# with w = rating or remaining headroom. It only orders candidates; acceptance
# is still the exact LP check.

module FlowRanking

using LinearAlgebra, SparseArrays
using Base.Threads: @threads, nthreads

export FlowSensitivity, build_sensitivity, merge_scores, apply_merge!,
       predicted_delta, radial_candidates

"""
Impedance and PTDF matrices for the network as currently merged, maintained by
rank-one updates so they stay exact without re-inversion.
"""
mutable struct FlowSensitivity
    N::Int
    L::Int
    Efrom::Vector{Int}
    Eto::Vector{Int}
    frate::Vector{Float64}
    X::Matrix{Float64}      # reference-removed inverse of B
    H::Matrix{Float64}      # L x N PTDF
    ref::Int
end

function build_sensitivity(base; ref::Int=base.j0)
    N, L = base.N, base.Ln
    I = vcat(base.Efrom, base.Eto)
    J = vcat(1:L, 1:L)
    V = vcat(fill(1.0, L), fill(-1.0, L))
    E = sparse(I, J, V, N, L)
    B = Matrix(E * Diagonal(base.Dx) * E')
    keep = [b for b in 1:N if b != ref]
    X = zeros(N, N)
    X[keep, keep] = inv(B[keep, keep])
    H = Diagonal(base.Dx) * Matrix(E') * X
    return FlowSensitivity(N, L, copy(base.Efrom), copy(base.Eto), copy(base.frate),
                           X, H, ref)
end

"z = e_i - e_j applied to a matrix column-wise: returns M[:,i] - M[:,j]."
@inline coldiff(M, i, j) = @view(M[:, i]) .- @view(M[:, j])

"""
Flow change on every line if buses i and j are shorted, at the given angles.
Returns the vector df and the transfer t; t is NaN when i and j are already the
same node (nothing to merge).
"""
function predicted_delta(fs::FlowSensitivity, theta::AbstractVector, i::Int, j::Int)
    i == j && return (df=zeros(fs.L), t=NaN)
    zXz = fs.X[i, i] + fs.X[j, j] - 2 * fs.X[i, j]
    # Already electrically identical: shorting them moves nothing.
    zXz <= 1e-12 && return (df=zeros(fs.L), t=0.0)
    d = theta[i] - theta[j]
    t = -d / zXz
    return (df=t .* coldiff(fs.H, i, j), t=t)
end

"""
Score every candidate line, lowest = least disturbance. `thetas` holds one angle
vector per scenario; the worst scenario decides, as in the acceptance test.

`skip(l)` lets the caller exclude lines already absorbed.
"""
function merge_scores(fs::FlowSensitivity, thetas::Vector{<:AbstractVector};
                      alpha::Float64=0.5, skip = l -> false, threaded::Bool=true,
                      norm_mode::Symbol=:rating,
                      flows::Union{Nothing,Vector{<:AbstractVector}}=nothing)
    scores = fill(Inf, fs.L)
    # Two ways to normalise the flow change before taking the norm.
    #
    #   :rating    df_k / fbar_k -- how big the disturbance is relative to the
    #              line. Treats a lightly loaded line and a nearly full one the
    #              same, which is wrong when the point is to avoid violations.
    #
    #   :headroom  df_k / (fbar_k - |f_k|) -- how big it is relative to what the
    #              line can still absorb. A 1 MW shift onto a line with 0.5 MW
    #              spare scores worse than 50 MW onto an empty one, which is the
    #              ordering that actually predicts rejection.
    #
    # Headroom needs the current flows, one vector per scenario, matching
    # `thetas`. Without them it falls back to :rating rather than guessing.
    use_headroom = norm_mode === :headroom && flows !== nothing &&
                   length(flows) == length(thetas)
    inv_rate = [1.0 / max(1e-9, fs.frate[l]) for l in 1:fs.L]
    # Floor the headroom at a small fraction of the rating so a line sitting
    # exactly on its limit gives a large-but-finite score instead of Inf.
    inv_head = if use_headroom
        [[1.0 / max(fs.frate[l] - abs(fv[l]), 1e-3 * max(fs.frate[l], 1e-9))
          for l in 1:fs.L] for fv in flows]
    else
        Vector{Vector{Float64}}()
    end

    # Each line's score depends only on read-only state (X, H, thetas) and is
    # written to its own slot, so this parallelises with no locking. Allocating
    # the delta vector per candidate would dominate, so the norms are computed
    # from the scalar and the PTDF column difference directly.
    function score_line(l)
        skip(l) && return Inf
        i, j = fs.Efrom[l], fs.Eto[l]
        i == j && return Inf
        zXz = fs.X[i, i] + fs.X[j, j] - 2 * fs.X[i, j]
        zXz <= 1e-12 && return Inf
        worst = 0.0
        for (si, theta) in enumerate(thetas)
            t = -(theta[i] - theta[j]) / zXz
            w = use_headroom ? inv_head[si] : inv_rate
            n2 = 0.0; ninf = 0.0
            @inbounds for k in 1:fs.L
                v = abs(t * (fs.H[k, i] - fs.H[k, j])) * w[k]
                n2 += v * v
                v > ninf && (ninf = v)
            end
            worst = max(worst, alpha * sqrt(n2) + (1 - alpha) * ninf)
        end
        return worst
    end

    if threaded && nthreads() > 1
        @threads for l in 1:fs.L
            scores[l] = score_line(l)
        end
    else
        for l in 1:fs.L
            scores[l] = score_line(l)
        end
    end
    return scores
end

"""
Lines whose collapse cannot move any OTHER line's flow, so they need no ranking
and no LP check to order.

A bus of degree 1 sits on exactly one line, so every unit injected there must
traverse that line whatever the rest of the network does. Formally
H[:,i] - H[:,j] is zero on every line except l itself, hence df = 0 on all
RETAINED lines: merging a radial bus is invisible to the rest of the system.

The catch, which is why this returns a flag rather than a blanket permission:
absorbing line l also DELETES its capacity constraint. A radial line's flow is
exactly the net injection at its leaf bus, so whether that matters is an O(1)
question rather than an LP -- if the leaf can never inject more than the rating,
the constraint was slack anyway and dropping it is free.

`dmax`/`dmin` are the worst-case demands at each bus over the scenario pool.
Returns one entry per radial line: its index, the leaf bus, and `capacity_safe`.
"""
function radial_candidates(base, dmax::AbstractVector, dmin::AbstractVector)
    N, L = base.N, base.Ln
    deg = zeros(Int, N)
    incident = [Int[] for _ in 1:N]
    for l in 1:L
        base.Efrom[l] == base.Eto[l] && continue
        deg[base.Efrom[l]] += 1; deg[base.Eto[l]] += 1
        push!(incident[base.Efrom[l]], l); push!(incident[base.Eto[l]], l)
    end
    gmax = zeros(N); gmin = zeros(N)
    for (k, b) in enumerate(base.gen_bus)
        gmax[b] += max(0.0, base.pmax[k])
        gmin[b] += min(0.0, base.pmin[k])
    end
    out = NamedTuple[]
    for i in 1:N
        deg[i] == 1 || continue
        l = incident[i][1]
        # Net injection at the leaf spans generation range minus demand range.
        hi = gmax[i] - dmin[i]
        lo = gmin[i] - dmax[i]
        safe = max(abs(hi), abs(lo)) <= base.frate[l] + 1e-9
        push!(out, (line=l, leaf=i, inj_hi=hi, inj_lo=lo,
                    rating=base.frate[l], capacity_safe=safe))
    end
    return out
end

"""
Fold a merge of buses i and j into X and H exactly, by rank-one update. Call
after a merge is ACCEPTED so later scores reflect the new network.
"""
function apply_merge!(fs::FlowSensitivity, i::Int, j::Int)
    i == j && return fs
    zXz = fs.X[i, i] + fs.X[j, j] - 2 * fs.X[i, j]
    zXz <= 1e-12 && return fs          # already shorted
    Xz = fs.X[:, i] .- fs.X[:, j]      # X z
    zX = fs.X[i, :] .- fs.X[j, :]      # z' X  (X symmetric, but keep it explicit)
    Hz = fs.H[:, i] .- fs.H[:, j]      # H z
    fs.X .-= (Xz * zX') ./ zXz
    fs.H .-= (Hz * zX') ./ zXz
    return fs
end

end # module

# Clustering for a certified reduction.
#
# Merging clusters a and b shorts them. By the rank-one (Sherman-Morrison) update,
# the flow on line l then changes by
#
#   -(H_l(a) - H_l(b)) * (current through the short),   current = (theta_a - theta_b) / R_ab
#
# with H the current reduced network's sensitivities and R_ab its effective
# resistance. So a merge costs, on every line that can overload,
#
#   |H_l(a) - H_l(b)| * max over the set |current|
#
# and is free when a and b are equally sensitive to l. Merge cheapest first,
# keeping each critical line's summed cost within beta * F_l; recompute on the
# merged network after every pass. Critical lines are never merged, not even by
# closing a loop. The cost only steers the clustering; certify/certificate.jl
# computes the exact error and the ratings.

module Coarsen

using LinearAlgebra
import Main.Certificate as CT

export coarsen

"Grounded inverse of the reduced network's susceptance matrix, K x K."
function reduced_inverse(base, cl, K)
    ref = cl[base.j0]
    free = setdiff(1:K, [ref])
    pos = zeros(Int, K)
    for (i, k) in enumerate(free)
        pos[k] = i
    end
    B = zeros(length(free), length(free))
    for l in 1:base.Ln
        a, b = pos[cl[base.Efrom[l]]], pos[cl[base.Eto[l]]]
        cl[base.Efrom[l]] == cl[base.Eto[l]] && continue
        y = base.Dx[l]
        a > 0 && (B[a, a] += y)
        b > 0 && (B[b, b] += y)
        a > 0 && b > 0 && (B[a, b] -= y; B[b, a] -= y)
    end
    X = zeros(K, K)
    isempty(free) || (X[free, free] = inv(Symmetric(B)))
    return X
end

"""
    coarsen(base, crit, S, beta; maxpass=200) -> internal line mask

`crit` the lines that can overload over the injection set `S`
(certificate.jl's injection_set). beta = 0 allows only free merges.
"""
function coarsen(base, crit, S, beta::Real; maxpass::Int=200)
    Ln, u, v = base.Ln, base.Efrom, base.Eto
    budget = beta .* base.frate[crit] .+ 1e-9
    total = zeros(length(crit))
    is_crit = falses(Ln)
    is_crit[crit] .= true
    internal = falses(Ln)
    for _ in 1:maxpass
        mp = CT.merged_ptdf(base, internal)
        cl, K = mp.cluster, mp.K
        rep = [findfirst(==(k), cl) for k in 1:K]
        Hc = mp.Hm[crit, rep]                   # sensitivity of each critical line to each cluster
        X = reduced_inverse(base, cl, K)
        allowed = Dict{Tuple{Int,Int},Bool}()
        for l in 1:Ln
            a, b = cl[u[l]], cl[v[l]]
            a == b && continue
            key = minmax(a, b)
            allowed[key] = get(allowed, key, true) && !is_crit[l]
        end
        cands = Tuple{Float64,Int,Int,Vector{Float64}}[]
        for ((a, b), ok) in allowed
            ok || continue
            R = X[a, a] + X[b, b] - 2 * X[a, b]
            row = (X[a, :] .- X[b, :]) ./ R      # current through a short a-b, per cluster injection
            lo, hi = CT.extremes(S, row[cl])
            I = max(hi, -lo)
            c = abs.(Hc[:, a] .- Hc[:, b]) .* I
            push!(cands, (isempty(c) ? 0.0 : maximum(c ./ budget), a, b, c))
        end
        sort!(cands; by=first)
        touched = Set{Int}()
        merged = 0
        for (_, a, b, c) in cands
            (a in touched || b in touched) && continue
            all(total .+ c .<= budget) || continue
            total .+= c
            push!(touched, a, b)
            merged += 1
            for l in 1:Ln                         # merge a and b: every line between them
                (cl[u[l]] == a && cl[v[l]] == b || cl[u[l]] == b && cl[v[l]] == a) &&
                    (internal[l] = true)
            end
        end
        merged == 0 && break
    end
    cl = CT.clusters(base, internal)
    return BitVector([cl[u[l]] == cl[v[l]] for l in 1:Ln])
end

end # module

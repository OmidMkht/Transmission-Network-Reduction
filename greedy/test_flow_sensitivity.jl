# --------------------------------------------------------------------------- #
# Validate the flow-sensitivity ordering against brute force.
#
#   julia --project=. --startup-file=no greedy/test_flow_sensitivity.jl [case]
#
# Two things must hold, and neither is obvious enough to take on trust:
#
#   1. PREDICTION. The closed form
#          df = -[(theta_i - theta_j)/(X_ii + X_jj - 2X_ij)] * (H[:,i] - H[:,j])
#      must equal the flow change actually obtained by rebuilding the network
#      with i and j shorted and recomputing flows at the SAME injection.
#
#   2. RANK-ONE UPDATE. After apply_merge!, X and H must equal the matrices of
#      the genuinely merged network -- otherwise predictions drift as merges
#      accumulate, which is precisely the failure the new ordering exists to fix.
#
# Brute force here means: build the shorted network explicitly by adding a very
# stiff line between i and j, invert, and compare.
# --------------------------------------------------------------------------- #

using LinearAlgebra, SparseArrays, Printf

const ROOT = dirname(@__DIR__)
const CASE = length(ARGS) >= 1 ? ARGS[1] : "pglib_opf_case118_ieee.m"

RUN = (scenarios = :single,)
@eval module TData
    include(joinpath(dirname(@__DIR__), "common", "preprocessing.jl"))
end
const TR = Main.TData
include(joinpath(@__DIR__, "flow_sensitivity.jl"))
const SO = Main.FlowRanking

c0 = TR.build_single_scenario_case(joinpath(ROOT, "case studies", CASE))
base = c0.base
base.c2 .= 0.0
N, L = base.N, base.Ln

"Flows and angles of a network whose susceptances are `dx`, at injection p."
function solve_flows(dx, p, ref)
    I = vcat(base.Efrom, base.Eto)
    J = vcat(1:L, 1:L)
    V = vcat(fill(1.0, L), fill(-1.0, L))
    E = sparse(I, J, V, N, L)
    B = Matrix(E * Diagonal(dx) * E')
    keep = [b for b in 1:N if b != ref]
    X = zeros(N, N)
    X[keep, keep] = inv(B[keep, keep])
    th = X * p
    f = Diagonal(dx) * Matrix(E') * th
    return (f=f, theta=th, X=X, H=Diagonal(dx) * Matrix(E') * X)
end

# A physically consistent injection: the base case's own.
p = base.p
fs = SO.build_sensitivity(base)
b0 = solve_flows(base.Dx, p, base.j0)

@printf("case %s: %d buses, %d lines\n", basename(CASE), N, L)
@printf("reference bus %d\n\n", base.j0)

# --------------------------------------------------------------------------- #
# 1. Prediction accuracy on single merges.
# --------------------------------------------------------------------------- #
# Shorting is approximated by a very stiff line; STIFF must be large enough that
# the residual angle difference is numerically zero but not so large that the
# matrix becomes ill-conditioned.
const STIFF = 1e7

println("1) closed-form prediction vs rebuilding the shorted network")
@printf("%6s %6s %6s %14s %14s %12s\n", "line", "i", "j", "max|df_pred|", "max|df_true|", "rel err")

function test_prediction(fs, b0, p, ntest)
    worst = 0.0
    tested = 0
    for l in 1:L
        tested >= ntest && break
        i, j = base.Efrom[l], base.Eto[l]
        i == j && continue
        r = SO.predicted_delta(fs, b0.theta, i, j)
        (isnan(r.t) || all(iszero, r.df)) && continue
        dx2 = copy(base.Dx)
        dx2[l] = STIFF                      # short this line
        b1 = solve_flows(dx2, p, base.j0)
        # The shorted line's own flow is meaningless (it absorbs the transfer),
        # so compare the OTHER lines, which is what the ordering scores.
        others = [k for k in 1:L if k != l]
        dtrue = b1.f[others] .- b0.f[others]
        dpred = r.df[others]
        # A merge that moves nothing carries no information about accuracy: the
        # relative error would just divide float noise by float noise. Skip it.
        if maximum(abs, dtrue) < 1e-4
            continue
        end
        rel = maximum(abs, dpred .- dtrue) / maximum(abs, dtrue)
        worst = max(worst, rel)
        tested += 1
        tested <= 8 && @printf("%6d %6d %6d %14.4f %14.4f %12.2e\n",
                               l, i, j, maximum(abs, dpred), maximum(abs, dtrue), rel)
    end
    return worst, tested
end

worst_pred, ntested = test_prediction(fs, b0, p, 40)
@printf("\nworst relative error over %d merges: %.3e  -> %s\n\n", ntested, worst_pred,
        worst_pred < 1e-4 ? "PASS" : "FAIL")

# --------------------------------------------------------------------------- #
# 2. Rank-one update stays exact over a sequence of merges.
# --------------------------------------------------------------------------- #
println("2) rank-one update of X and H vs full re-inversion, over repeated merges")
@printf("%6s %8s %16s %16s\n", "merge", "line", "max|dX|", "max|dH|")

function test_updates(nmerge)
    fs2 = SO.build_sensitivity(base)
    dx = copy(base.Dx)
    shorted = Int[]                 # lines whose own PTDF row is meaningless
    worstX = 0.0; worstH = 0.0
    done = 0
    for l in 1:L
        done >= nmerge && break
        i, j = base.Efrom[l], base.Eto[l]
        i == j && continue
        # Skip if already electrically shorted by an earlier merge.
        (fs2.X[i, i] + fs2.X[j, j] - 2 * fs2.X[i, j]) <= 1e-9 && continue
        SO.apply_merge!(fs2, i, j)
        dx[l] = STIFF
        push!(shorted, l)
        ref = solve_flows(dx, p, base.j0)
        # Compare only the rows that still mean something. A shorted line takes
        # whatever transfer the merge needs, so the reference gives it a PTDF row
        # near 1 while the update keeps its original susceptance -- a difference
        # of definition, not an error. The ordering never scores those rows.
        keep = [k for k in 1:L if !(k in shorted)]
        dX = maximum(abs, fs2.X .- ref.X)
        dH = isempty(keep) ? 0.0 : maximum(abs, fs2.H[keep, :] .- ref.H[keep, :])
        worstX = max(worstX, dX); worstH = max(worstH, dH)
        done += 1
        done <= 8 && @printf("%6d %8d %16.3e %16.3e\n", done, l, dX, dH)
    end
    return worstX, worstH, done
end

wX, wH, nmerged = test_updates(25)
# Successive rank-one updates accumulate float error: |dH| grows steadily from
# ~3e-7 after one merge to ~1e-5 after 25. That is ordinary drift, not a wrong
# formula, and PTDF entries are O(1) so this is a ~1e-5 RELATIVE error. It is
# immaterial here because the score only ORDERS candidates and acceptance is
# decided by an exact LP. Rebuild the sensitivity from scratch if a run ever
# needs many hundreds of merges.
@printf("\nworst |dX| %.3e, worst |dH| %.3e over %d merges -> %s\n\n",
        wX, wH, nmerged, (wX < 1e-4 && wH < 1e-4) ? "PASS" : "FAIL")

# --------------------------------------------------------------------------- #
# 3. Does the new score order candidates differently from utilisation?
# --------------------------------------------------------------------------- #
println("3) new ordering vs the existing utilisation ordering")
util = [abs(base.fhat[l]) / base.frate[l] for l in 1:L]
old_order = sortperm(util)
scores = SO.merge_scores(fs, [b0.theta]; alpha=0.5)
new_order = sortperm(scores)
overlap10 = length(intersect(Set(old_order[1:min(10, L)]), Set(new_order[1:min(10, L)])))
overlap25 = length(intersect(Set(old_order[1:min(25, L)]), Set(new_order[1:min(25, L)])))
@printf("  first 10 candidates shared: %d of 10\n", overlap10)
@printf("  first 25 candidates shared: %d of 25\n", overlap25)
@printf("  top 5 by utilisation : %s\n", join(old_order[1:5], ", "))
@printf("  top 5 by flow delta  : %s\n", join(new_order[1:5], ", "))
println()
println(overlap10 <= 7 ? ">>> the orderings genuinely differ, so this can change what the seed finds." :
                         ">>> the orderings mostly agree; expect little behavioural change.")

# --------------------------------------------------------------------------- #
# 4. The point of the proposal: does RE-RANKING after merges diverge from the
#    static order? Comparing only at iteration 0 understates the idea, because
#    the existing heuristic cannot re-rank at all.
# --------------------------------------------------------------------------- #
println("4) does re-ranking drift away from the static order as merges accumulate?")
@printf("%8s %10s %14s %14s\n", "merges", "cand left", "top-10 shared", "top-1 same?")

function test_drift(nsteps)
    fs3 = SO.build_sensitivity(base)
    merged = falses(L)
    theta = copy(b0.theta)
    dx = copy(base.Dx)
    static_rank = sortperm([abs(base.fhat[l]) / base.frate[l] for l in 1:L])
    rows = NamedTuple[]
    for step in 0:nsteps
        remaining = [l for l in 1:L if !merged[l] && base.Efrom[l] != base.Eto[l]]
        length(remaining) < 12 && break
        sc = SO.merge_scores(fs3, [theta]; alpha=0.5, skip = l -> merged[l])
        dyn = sort(remaining; by = l -> sc[l])
        stat = [l for l in static_rank if !merged[l]]
        shared = length(intersect(Set(dyn[1:10]), Set(stat[1:10])))
        push!(rows, (step=step, left=length(remaining), shared=shared,
                     same_top = dyn[1] == stat[1]))
        @printf("%8d %10d %11d/10 %14s\n", step, length(remaining), shared,
                dyn[1] == stat[1] ? "yes" : "NO")
        # Take the dynamically-best merge and move the operating point with it.
        l = dyn[1]
        i, j = base.Efrom[l], base.Eto[l]
        r = SO.predicted_delta(fs3, theta, i, j)
        SO.apply_merge!(fs3, i, j)
        merged[l] = true
        dx[l] = STIFF
        theta = solve_flows(dx, p, base.j0).theta
    end
    return rows
end

drift = test_drift(10)
println()
if !isempty(drift)
    later = [r for r in drift if r.step >= 3]
    avg = isempty(later) ? NaN : sum(r.shared for r in later) / length(later)
    ndiff = count(r -> !r.same_top, drift)
    @printf("after 3+ merges the two orders share %.1f of their top 10 on average\n", avg)
    @printf("the top candidate differed at %d of %d steps\n", ndiff, length(drift))
    println(avg < 8 || ndiff > 0 ?
        ">>> re-ranking DOES diverge once merges accumulate -- the dynamic score is doing work." :
        ">>> re-ranking tracks the static order closely -- little to gain here.")
end

# --------------------------------------------------------------------------- #
# 5. Radial buses: collapsing one must leave every OTHER line's flow untouched,
#    so they can be taken without ranking. Verify against brute force, and
#    report how many are also free of a capacity concern.
# --------------------------------------------------------------------------- #
println("5) radial (degree-1) buses: collapse must not move any other flow")
dmax = vec(maximum(c0.load; dims=2))
dmin = vec(minimum(c0.load; dims=2))
rad = SO.radial_candidates(base, dmax, dmin)
@printf("  %d radial buses of %d; %d of those are also capacity-safe\n",
        length(rad), N, count(r -> r.capacity_safe, rad))

function test_radial(rad)
    worst = 0.0
    shown = 0
    for r in rad
        l = r.line
        dx2 = copy(base.Dx); dx2[l] = STIFF
        b1 = solve_flows(dx2, p, base.j0)
        others = [k for k in 1:L if k != l]
        moved = maximum(abs, b1.f[others] .- b0.f[others])
        worst = max(worst, moved)
        shown += 1
        shown <= 6 && @printf("  line %4d leaf %4d  inj [%.1f, %.1f] vs rating %.1f  %-12s max|df_other| = %.2e\n",
                              l, r.leaf, r.inj_lo, r.inj_hi, r.rating,
                              r.capacity_safe ? "cap-safe" : "CAP RISK", moved)
    end
    return worst
end

worst_radial = isempty(rad) ? 0.0 : test_radial(rad)
@printf("\n  worst flow moved on any other line by a radial collapse: %.3e -> %s\n\n",
        worst_radial, worst_radial < 1e-6 ? "PASS (free merges confirmed)" : "FAIL")

# --------------------------------------------------------------------------- #
# 6. Threading.
# --------------------------------------------------------------------------- #
println("6) threaded scoring")
@printf("  julia threads available: %d\n", Threads.nthreads())
let th = [b0.theta]
    s1 = SO.merge_scores(fs, th; alpha=0.5, threaded=false)
    t1 = @elapsed SO.merge_scores(fs, th; alpha=0.5, threaded=false)
    s2 = SO.merge_scores(fs, th; alpha=0.5, threaded=true)
    t2 = @elapsed SO.merge_scores(fs, th; alpha=0.5, threaded=true)
    agree = maximum(abs, replace(s1, Inf => 0.0) .- replace(s2, Inf => 0.0))
    @printf("  serial %.4f s, threaded %.4f s (%.2fx)\n", t1, t2, t1 / max(t2, 1e-9))
    @printf("  identical scores: %s (max diff %.2e)\n", agree < 1e-12 ? "yes" : "NO", agree)
end
println()

ok = worst_pred < 1e-4 && wX < 1e-4 && wH < 1e-4 && worst_radial < 1e-6
println()
println(ok ? ">>> ALL CHECKS PASS" : ">>> CHECKS FAILED -- do not use this ordering")
exit(ok ? 0 : 1)

# --------------------------------------------------------------------------- #
# KRON REDUCTION -- collapse chains of "boring" buses (degree 2, no generator,
# not near any congested line) into one equivalent line each, before the
# reduction MILP ever sees them. Included inside the same TNR module as the
# rest of the pipeline, so everything here dispatches on TxReductionCase /
# MultiScenarioTxReductionCase and reuses incidence_matrix, extract_reduction,
# etc. directly. Nothing in the parent directory is modified by this file.
#
# See kron_reduction/../.claude-adjacent design discussion for the two rules
# this all serves: an eliminated bus rejoins whichever cluster its equivalent
# line's two endpoints end up in (if they merge), or is reinserted with its
# true original topology and injection (if the equivalent line stays
# external) -- see kron_unfold.jl for that half.
# --------------------------------------------------------------------------- #

# Per-bus incident line-index lists. NOT Graphs.SimpleGraph, which collapses
# parallel circuits into one edge -- exactly the pitfall exactly_mergeable_lines
# already documents. A degree-2 bus with two lines to the SAME neighbour is a
# parallel-circuit dead-end, not a pass-through, and needs both lines visible
# to be excluded correctly (see kron_eligible_buses).
function kron_incident_lines(base::TxReductionCase)
    incident = [Int[] for _ in 1:base.N]
    for l in 1:base.Ln
        push!(incident[base.Efrom[l]], l)
        push!(incident[base.Eto[l]], l)
    end
    return incident
end

other_endpoint(base::TxReductionCase, l::Int, b::Int) =
    base.Efrom[l] == b ? base.Eto[l] : base.Efrom[l]

"""
    kron_worst_utilization(c::MultiScenarioTxReductionCase, scenario_indices) -> Vector{Float64}

`|fhat_l| / frate_l`, maximized over `scenario_indices` -- the same worst-case
aggregation `multiscenario_windows`'s `near` test uses per scenario, reduced
to one number per line so chain eligibility is a single static screen.
"""
function kron_worst_utilization(c::MultiScenarioTxReductionCase, scenario_indices)
    base = c.base
    s = Int.(collect(scenario_indices))
    isempty(s) && error("At least one scenario is required")
    return vec(maximum(abs.(c.fhat[:, s]) ./ base.frate; dims=2))
end

"""
    kron_eligible_buses(base, worst_utilization; near_limit_threshold) -> BitVector

A bus is eligible for Kron elimination iff: degree exactly 2, no generator
attached, is not the reference bus (base.j0 is the angle reference every
solve/validator pins to 0 -- excluding it outright avoids inventing
reference-bus remapping), its two incident lines lead to two DISTINCT
neighbours (a degree-2 dead-end into the same neighbour twice is a parallel
circuit, not a pass-through), and BOTH incident lines have
worst_utilization < near_limit_threshold. This last condition is also what
guarantees, structurally, that a Kron-eliminated chain can never mask a truly
near-congested original line: nothing in it was ever close to its limit.
"""
function kron_eligible_buses(base::TxReductionCase, worst_utilization::AbstractVector;
                             near_limit_threshold::Real)
    incident = kron_incident_lines(base)
    gen_buses = Set(base.gen_bus)
    eligible = falses(base.N)
    for b in 1:base.N
        b == base.j0 && continue
        b in gen_buses && continue
        length(incident[b]) == 2 || continue
        l1, l2 = incident[b]
        n1, n2 = other_endpoint(base, l1, b), other_endpoint(base, l2, b)
        n1 == n2 && continue
        (worst_utilization[l1] < near_limit_threshold &&
         worst_utilization[l2] < near_limit_threshold) || continue
        eligible[b] = true
    end
    return eligible
end

"""
    kron_chains(base, eligible, incident) -> (chains, rings, self_loop_chains)

The eligible-bus subgraph is degree-2 everywhere in it, so it decomposes into
disjoint simple paths and simple cycles.

  * Each maximal path is a `chain`: `(x, y, interior, chain_lines)`, x/y the
    anchor (non-eligible) buses at each end, `interior` the eliminated buses
    in walk order, `chain_lines` the `length(interior)+1` original line
    indices in walk order (leg x-b1, interior legs, leg bk-y).
  * A maximal cycle entirely inside the eligible set (never touches a
    non-eligible bus) has no anchor pair and cannot be represented as one
    equivalent line -- returned in `rings`, left untouched.
  * A path whose two ends are the SAME anchor bus would require a shunt at
    that bus, which TxReductionCase cannot express (series lines only) --
    returned in `self_loop_chains`, left untouched.

Two chains sharing one anchor (a hub with several pendant chains) or sharing
both anchors (parallel chains, or a chain whose anchor pair already has a
direct line) both fall out for free from this construction: `chains` is just
a flat list, more than one entry may share an (x,y) pair, exactly like
ACTIVSg200's existing double circuits.
"""
function kron_chains(base::TxReductionCase, eligible::BitVector, incident)
    N = base.N
    visited = falses(N)
    chains = NamedTuple[]
    rings = Vector{NamedTuple}()
    self_loop_chains = NamedTuple[]

    # Walk from the first interior bus `first` away from anchor `anchor`,
    # along `first_line`, until hitting a non-eligible bus or (for a ring)
    # back to a visited eligible bus.
    function walk(anchor::Int, first::Int, first_line::Int)
        interior = Int[first]
        chain_lines = Int[first_line]
        visited[first] = true
        prev, cur = anchor, first
        while eligible[cur]
            l1, l2 = incident[cur]
            nb1, nb2 = other_endpoint(base, l1, cur), other_endpoint(base, l2, cur)
            nxt, lnxt = nb1 == prev ? (nb2, l2) : (nb1, l1)
            push!(chain_lines, lnxt)
            if eligible[nxt]
                visited[nxt] = true
                push!(interior, nxt)
            end
            prev, cur = cur, nxt
        end
        return interior, chain_lines, cur
    end

    # Pass 1: start from every anchor -> eligible-bus edge.
    for l in 1:base.Ln
        a, b = base.Efrom[l], base.Eto[l]
        a_elig, b_elig = eligible[a], eligible[b]
        if !a_elig && b_elig && !visited[b]
            interior, chain_lines, y = walk(a, b, l)
            push!(y == a ? self_loop_chains : chains,
                  (x=a, y=y, interior=interior, chain_lines=chain_lines))
        elseif a_elig && !b_elig && !visited[a]
            interior, chain_lines, y = walk(b, a, l)
            push!(y == b ? self_loop_chains : chains,
                  (x=b, y=y, interior=interior, chain_lines=chain_lines))
        end
    end

    # Pass 2: any eligible bus not yet visited belongs to a pure ring (never
    # adjacent to any anchor), walked in an arbitrary direction back to itself.
    for start in 1:N
        (eligible[start] && !visited[start]) || continue
        l1, _ = incident[start]
        interior, chain_lines, closer = walk(start, other_endpoint(base, l1, start), l1)
        @assert closer == start "ring walk did not close"
        push!(rings, (start=start, interior=vcat(start, interior), chain_lines=chain_lines))
    end

    return chains, rings, self_loop_chains
end

"""
    kron_chain_equivalent(base, chain) -> (Dx_eq::Float64, transfer::Matrix{Float64})

Builds the local (k+2)x(k+2) susceptance-weighted Laplacian for
{x, interior..., y} from ONLY the chain's own `chain_lines` (exact: chain
eligibility already guarantees an interior bus has no other incident line),
partitions it into boundary rows/cols r=[x,y] and interior rows/cols
kk=interior, and returns:

  Dx_eq     the Schur complement collapsed to the scalar series susceptance
            between x and y (equals the classic 1/sum(1/Dx_i) for a plain
            chain -- kept as a Schur complement rather than that closed form
            directly so the same computation also produces `transfer`).
  transfer  2 x k matrix T = B_rk * inv(B_kk). For any scenario's interior
            injection vector p_k (net, generation - load; generation is
            always 0 here by eligibility), the boundary injection add-on is
            `p_r_eq = p_r - T * p_k` (standard Schur-complement elimination
            of `[[B_rr,B_rk];[B_kr,B_kk]] * [th_r;th_k] = [p_r;p_k]`).
"""
function kron_chain_equivalent(base::TxReductionCase, chain)
    interior = chain.interior
    k = length(interior)
    local_bus = Dict{Int,Int}(chain.x => 1, chain.y => k + 2)
    for (i, b) in enumerate(interior)
        local_bus[b] = i + 1
    end
    B = zeros(k + 2, k + 2)
    for l in chain.chain_lines
        p = local_bus[base.Efrom[l]]
        q = local_bus[base.Eto[l]]
        w = base.Dx[l]
        B[p, p] += w; B[q, q] += w
        B[p, q] -= w; B[q, p] -= w
    end
    r = [1, k + 2]
    kk = 2:(k + 1)
    B_rr = B[r, r]
    B_rk = B[r, kk]
    B_kk = B[kk, kk]
    B_kr = B[kk, r]
    F = factorize(B_kk)
    T = B_rk * (F \ Matrix(I, k, k))
    B_sc = B_rr - T * B_kr
    Dx_eq = B_sc[1, 1]
    return Dx_eq, T
end

"""
    kron_chain_ratings(base, chain, epsL_full, worst_utilization) -> (frate_eq, epsL_eq, driving_line)

`frate_eq = max(frate)` over the chain -- safe once `epsL_eq` is tight,
since the window clipping `min(frate_eq, fhat_eq+epsL_eq)` is then dominated
by the tight epsL term regardless.

`epsL_eq = min(epsL)` over the chain. Not recomputed as a fraction of
`frate_eq` (that would silently loosen it back up). Originally this used the
epsL of the highest-utilization ("driving") line instead of the minimum;
measured on case118 and ACTIVSg200, that let a reinserted (non-merged)
chain's OTHER lines end up checked -- by benchmark_reduction_scenarios,
against their OWN epsL after unfolding -- against a window effectively wider
than what the equivalent line itself was allowed, so a reinserted line could
fail its own window even though the boundary solve was fully valid on its
own terms. Taking the minimum means no reinserted line is ever asked to fit
inside less room than the equivalent line had. `driving_line` (highest
utilization, kept for diagnostics/reporting) is tracked separately from
whichever line actually sets `epsL_eq`.
"""
function kron_chain_ratings(base::TxReductionCase, chain, epsL_full::AbstractVector,
                            worst_utilization::AbstractVector)
    frate_eq = maximum(base.frate[l] for l in chain.chain_lines)
    driving_line = chain.chain_lines[1]
    best_u = worst_utilization[driving_line]
    for l in chain.chain_lines
        if worst_utilization[l] > best_u
            best_u = worst_utilization[l]
            driving_line = l
        end
    end
    epsL_eq = minimum(epsL_full[l] for l in chain.chain_lines)
    return frate_eq, epsL_eq, driving_line
end

"""
    kron_reduce_case(c, epsL_full; near_limit_threshold,
                     eligibility_scenario_indices=axes(c.p,2), min_chain_length=1)
      -> (c_boundary::MultiScenarioTxReductionCase, epsL_boundary::Vector{Float64}, kron_map)

Top-level entry point. `eligibility_scenario_indices` defaults to every
loaded scenario (mirroring how `protection_indices` defaults elsewhere), so
eligibility reflects the whole horizon, not whatever subset later becomes the
MILP's active set.

`c_boundary`'s singular (non-scenario) TxReductionCase fields
(p/thetahat/fhat/Pd) are populated from scenario column 1 as vestigial
placeholders -- confirmed by inspection that no function this module calls
with `c_boundary` reads them; never pass `c_boundary` to
`report_case`/`plot_network`/`write_multiscenario_outputs`.
"""
function kron_reduce_case(c::MultiScenarioTxReductionCase, epsL_full::AbstractVector;
                          near_limit_threshold::Real,
                          eligibility_scenario_indices=axes(c.p, 2),
                          min_chain_length::Int=1)
    base = c.base
    S = size(c.p, 2)
    length(epsL_full) == base.Ln ||
        error("epsL_full must have $(base.Ln) entries")
    min_chain_length >= 1 || error("min_chain_length must be at least 1")

    worst_u = kron_worst_utilization(c, eligibility_scenario_indices)
    eligible = kron_eligible_buses(base, worst_u; near_limit_threshold=near_limit_threshold)
    incident = kron_incident_lines(base)
    chains_all, rings, self_loop_chains = kron_chains(base, eligible, incident)
    chains = filter(ch -> length(ch.interior) >= min_chain_length, chains_all)

    # Buses removed from the boundary set: interior buses of every kept chain.
    removed = falses(base.N)
    for ch in chains, b in ch.interior
        removed[b] = true
    end

    boundary_of = zeros(Int, base.N)   # full_idx -> boundary_idx, 0 if eliminated
    full_of = Int[]                    # boundary_idx -> full_idx
    for b in 1:base.N
        removed[b] && continue
        push!(full_of, b)
        boundary_of[b] = length(full_of)
    end
    boundary_N = length(full_of)

    # Lines: every line with neither endpoint removed is copied through
    # unchanged; every kept chain contributes exactly one new equivalent line.
    chain_of_line = Dict{Int,Int}()   # original chain_lines[l] -> chain index, for bookkeeping
    for (ci, ch) in enumerate(chains), l in ch.chain_lines
        chain_of_line[l] = ci
    end
    removed_lines = Set(keys(chain_of_line))

    Efrom_b = Int[]; Eto_b = Int[]; Dx_b = Float64[]; frate_b = Float64[]; epsL_b = Float64[]
    boundary_line_of_full = zeros(Int, base.Ln)   # full line idx -> boundary line idx, 0 if collapsed into a chain
    for l in 1:base.Ln
        l in removed_lines && continue
        push!(Efrom_b, boundary_of[base.Efrom[l]])
        push!(Eto_b, boundary_of[base.Eto[l]])
        push!(Dx_b, base.Dx[l])
        push!(frate_b, base.frate[l])
        push!(epsL_b, epsL_full[l])
        boundary_line_of_full[l] = length(Efrom_b)
    end
    chain_boundary_line = zeros(Int, length(chains))
    chain_Dx_eq = zeros(length(chains)); chain_frate_eq = zeros(length(chains))
    chain_epsL_eq = zeros(length(chains)); chain_driving_line = zeros(Int, length(chains))
    chain_transfer = Vector{Matrix{Float64}}(undef, length(chains))
    for (ci, ch) in enumerate(chains)
        Dx_eq, T = kron_chain_equivalent(base, ch)
        frate_eq, epsL_eq, driving_line = kron_chain_ratings(base, ch, epsL_full, worst_u)
        push!(Efrom_b, boundary_of[ch.x]); push!(Eto_b, boundary_of[ch.y])
        push!(Dx_b, Dx_eq); push!(frate_b, frate_eq); push!(epsL_b, epsL_eq)
        chain_boundary_line[ci] = length(Efrom_b)
        chain_Dx_eq[ci] = Dx_eq; chain_frate_eq[ci] = frate_eq
        chain_epsL_eq[ci] = epsL_eq; chain_driving_line[ci] = driving_line
        chain_transfer[ci] = T
    end
    Ln_b = length(Efrom_b)

    gen_bus_b = [boundary_of[g] for g in base.gen_bus]
    any(==(0), gen_bus_b) && error("a generator bus was unexpectedly eliminated")

    # Injection redistribution. generation is unchanged (eligibility forbids a
    # generator on any interior bus, so nothing to move); load picks up each
    # chain's transferred contribution: load_r_eq = load_r - T*load_interior
    # (equivalently p_r_eq = p_r - T*p_interior, since p = generation - load
    # and generation_interior == 0 identically -- see kron_chain_equivalent).
    load_b = c.load[full_of, :]
    generation_b = c.generation[full_of, :]
    for (ci, ch) in enumerate(chains)
        x_b, y_b = boundary_of[ch.x], boundary_of[ch.y]
        load_interior = c.load[ch.interior, :]
        correction = chain_transfer[ci] * load_interior   # 2 x S
        load_b[x_b, :] .-= correction[1, :]
        load_b[y_b, :] .-= correction[2, :]
    end
    p_b = generation_b .- load_b

    # thetahat/fhat: one fresh linear DC solve on the boundary topology for all
    # S columns, factorized once -- exact, by the Schur-complement guarantee
    # that boundary-bus angles are unchanged by eliminating zero-residual
    # interior structure. Mirrors redispatch_dc_opf_scenarios's own
    # angle re-derivation pattern.
    boundary_base_stub = TxReductionCase(
        boundary_N, Ln_b, Efrom_b, Eto_b, Dx_b, frate_b,
        zeros(boundary_N), zeros(boundary_N), zeros(Ln_b), zeros(boundary_N),
        gen_bus_b, copy(base.pmin), copy(base.pmax),
        copy(base.c2), copy(base.c1), copy(base.c0),
        boundary_of[base.j0], base.baseMVA,
    )
    Eb = incidence_matrix(boundary_base_stub)
    Bb = Eb * Diagonal(Dx_b) * Eb'
    freebus = setdiff(1:boundary_N, [boundary_of[base.j0]])
    Fb = factorize(Bb[freebus, freebus])
    thetahat_b = zeros(boundary_N, S)
    fhat_b = zeros(Ln_b, S)
    for s in 1:S
        th = zeros(boundary_N)
        th[freebus] = Fb \ p_b[freebus, s]
        thetahat_b[:, s] = th
        fhat_b[:, s] = Dx_b .* (th[Efrom_b] .- th[Eto_b])
    end

    base_boundary = TxReductionCase(
        boundary_N, Ln_b, Efrom_b, Eto_b, Dx_b, frate_b,
        copy(p_b[:, 1]), copy(thetahat_b[:, 1]), copy(fhat_b[:, 1]), copy(load_b[:, 1]),
        gen_bus_b, copy(base.pmin), copy(base.pmax),
        copy(base.c2), copy(base.c1), copy(base.c0),
        boundary_of[base.j0], base.baseMVA,
    )
    c_boundary = MultiScenarioTxReductionCase(
        base_boundary, copy(c.scenario_ids), c.bus_ids[full_of],
        load_b, generation_b, p_b, thetahat_b, fhat_b,
    )

    chains_full = [merge(ch, (
        boundary_line=chain_boundary_line[ci], Dx_eq=chain_Dx_eq[ci],
        frate_eq=chain_frate_eq[ci], epsL_eq=chain_epsL_eq[ci],
        driving_line=chain_driving_line[ci], transfer=chain_transfer[ci],
    )) for (ci, ch) in enumerate(chains)]

    kron_map = (
        full_N=base.N, boundary_N=boundary_N, full_of=full_of, boundary_of=boundary_of,
        chains=chains_full, rings=rings, self_loop_chains=self_loop_chains,
        eligible=eligible, near_limit_threshold=near_limit_threshold,
        eligibility_scenario_indices=Int.(collect(eligibility_scenario_indices)),
        boundary_line_of_full=boundary_line_of_full,
    )
    return c_boundary, epsL_b, kron_map
end

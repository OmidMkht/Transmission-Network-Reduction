# --------------------------------------------------------------------------- #
# Target-directed reduction: the bilevel program of README.md, reformulated to
# a single MILP through the lower level's KKT conditions.
#
# Placement I -- full-network feasibility (T1) is an UPPER-level constraint on
# the lower level's solution, so the model establishes existence of a reduced
# optimum that is feasible on the true network on the modelled scenarios. A
# separately deployed solver may choose another optimum. The
# older cost-gap condition (T2) remains available as an explicit option.
#
# Needs common/caps.jl loaded first (module Caps).
# --------------------------------------------------------------------------- #

module BilevelReduction

using JuMP, Gurobi, LinearAlgebra, Printf
using Main.Caps: capped_paths, size_witness, long_chain, merged_path

export solve_bilevel_reduction

"Canonical preprocessing helper defined beside the case type."
function case_helper(c, name::Symbol)
    mod = parentmodule(typeof(c))
    isdefined(mod, name) || error(
        "$(nameof(mod)).$name is required. Include common/preprocessing.jl and " *
        "common/postprocessing.jl before kkt_model.jl.")
    return getfield(mod, name)
end

# --------------------------------------------------------------------------- #
# The true full-network optimum, per scenario. This is LL(c = 0): setting every
# line external makes the lower level the full DC-OPF exactly, so Z* is computed
# by the same equations the bilevel uses rather than by a second, parallel
# implementation that could drift from it.
# --------------------------------------------------------------------------- #
function true_optimal_cost(base, demand::AbstractVector; relax_pmin::Bool=true,
                           env=nothing)
    N, Ln = base.N, base.Ln
    K = length(base.gen_bus)
    pmin = relax_pmin ? min.(0.0, base.pmin) : base.pmin

    m = isnothing(env) ? Model(Gurobi.Optimizer) : Model(() -> Gurobi.Optimizer(env))
    set_optimizer_attribute(m, "OutputFlag", 0)
    @variable(m, th[1:N])
    @variable(m, g[1:K])
    @variable(m, f[1:Ln])
    @constraint(m, th[base.j0] == 0)
    @constraint(m, [k=1:K], pmin[k] <= g[k] <= base.pmax[k])
    @constraint(m, [l=1:Ln], f[l] == base.Dx[l] * (th[base.Efrom[l]] - th[base.Eto[l]]))
    @constraint(m, [l=1:Ln], -base.frate[l] <= f[l] <= base.frate[l])
    @constraint(m, [b=1:N],
        sum(f[l] for l in 1:Ln if base.Efrom[l] == b; init=0.0) -
        sum(f[l] for l in 1:Ln if base.Eto[l]   == b; init=0.0) ==
        sum(g[k] for k in 1:K if base.gen_bus[k] == b; init=0.0) - demand[b])
    @objective(m, Min, sum(base.c1[k] * g[k] + base.c0[k] for k in 1:K))
    optimize!(m)
    termination_status(m) in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) ||
        error("Full-network DC-OPF failed: $(termination_status(m))")
    return objective_value(m)
end

"""Generator-free radial lines that can be contracted without changing DC-OPF."""
function radial_internal_mask(base; forbidden=falses(base.Ln))
    length(forbidden) == base.Ln || error("forbidden must have length $(base.Ln)")
    rep = collect(1:base.N)
    find(x) = (rep[x] == x ? x : (rep[x] = find(rep[x])))
    has_gen = falses(base.N)
    for b in base.gen_bus
        has_gen[b] = true
    end
    alive = trues(base.Ln)
    while true
        deg = zeros(Int, base.N)
        for l in 1:base.Ln
            alive[l] || continue
            a, b = find(base.Efrom[l]), find(base.Eto[l])
            a == b && continue
            deg[a] += 1
            deg[b] += 1
        end
        merged = false
        for l in 1:base.Ln
            (alive[l] && !forbidden[l]) || continue
            a, b = find(base.Efrom[l]), find(base.Eto[l])
            a == b && continue
            leaf = deg[a] == 1 && !has_gen[a] ? a :
                   deg[b] == 1 && !has_gen[b] ? b : 0
            leaf == 0 && continue
            other = leaf == a ? b : a
            rep[leaf] = other
            has_gen[other] |= has_gen[leaf]
            alive[l] = false
            merged = true
        end
        merged || break
    end
    rep_of = [find(b) for b in 1:base.N]
    internal = [!forbidden[l] && rep_of[base.Efrom[l]] == rep_of[base.Eto[l]]
                for l in 1:base.Ln]
    return (internal=internal, rep_of=rep_of, lines=findall(internal))
end

"""A complete primal-dual KKT point for one fixed reduction topology."""
function fixed_topology_kkt_start(base, demand, pmin, internal::AbstractVector{Bool};
                                  status_out=nothing, time_limit=nothing,
                                  solver_threads=nothing)
    N, Ln = base.N, base.Ln
    K, S = length(base.gen_bus), length(demand)
    length(internal) == Ln || error("internal must have length $Ln")
    u, v, F = base.Efrom, base.Eto, base.frate
    incident_from = [findall(==(b), u) for b in 1:N]
    incident_to   = [findall(==(b), v) for b in 1:N]
    gens_at       = [findall(==(b), base.gen_bus) for b in 1:N]

    env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))

    m = Model(() -> Gurobi.Optimizer(env))
    set_optimizer_attribute(m, "OutputFlag", 0)
    isnothing(time_limit) || set_optimizer_attribute(m, "TimeLimit", time_limit)
    isnothing(solver_threads) || set_optimizer_attribute(m, "Threads", solver_threads)

    @variable(m, g[1:K, 1:S])
    @variable(m, th[1:N, 1:S])
    @variable(m, f[1:Ln, 1:S])
    @variable(m, t[1:Ln, 1:S])
    @variable(m, lam[1:N, 1:S])
    @variable(m, sig[1:Ln, 1:S])
    @variable(m, pi0[1:S])
    @variable(m, mup[1:Ln, 1:S] >= 0)
    @variable(m, mum[1:Ln, 1:S] >= 0)
    @variable(m, nu[1:Ln, 1:S])
    @variable(m, tau[1:Ln, 1:S])
    @variable(m, rhop[1:K, 1:S] >= 0)
    @variable(m, rhom[1:K, 1:S] >= 0)
    @variable(m, tth[1:N, 1:S])
    @variable(m, tf[1:Ln, 1:S])

    for l in 1:Ln, s in 1:S
        if internal[l]
            @constraint(m, f[l, s] == 0)
            fix(mup[l, s], 0.0; force=true)
            fix(mum[l, s], 0.0; force=true)
            fix(tau[l, s], 0.0; force=true)
        else
            @constraint(m, f[l, s] <= F[l])
            @constraint(m, -f[l, s] <= F[l])
            @constraint(m, t[l, s] == 0)
            fix(nu[l, s], 0.0; force=true)
        end
    end
    @constraint(m, [k=1:K, s=1:S], g[k, s] <= base.pmax[k])
    @constraint(m, [k=1:K, s=1:S], pmin[k] <= g[k, s])
    @constraint(m, [l=1:Ln, s=1:S],
        f[l, s] == base.Dx[l] * (th[u[l], s] - th[v[l], s]))
    @constraint(m, [s=1:S], th[base.j0, s] == 0)
    @constraint(m, [b=1:N, s=1:S],
        sum(f[l, s] + t[l, s] for l in incident_from[b]; init=0.0) -
        sum(f[l, s] + t[l, s] for l in incident_to[b]; init=0.0) ==
        sum(g[k, s] for k in gens_at[b]; init=0.0) - demand[s][b])

    @constraint(m, [k=1:K, s=1:S],
        base.c1[k] - lam[base.gen_bus[k], s] + rhop[k, s] - rhom[k, s] == 0)
    @constraint(m, [l=1:Ln, s=1:S],
        lam[u[l], s] - lam[v[l], s] + sig[l, s] +
        mup[l, s] - mum[l, s] + nu[l, s] == 0)
    @constraint(m, [l=1:Ln, s=1:S],
        lam[u[l], s] - lam[v[l], s] + tau[l, s] == 0)
    @constraint(m, [b=1:N, s=1:S],
        -sum(base.Dx[l] * sig[l, s] for l in incident_from[b]; init=0.0) +
         sum(base.Dx[l] * sig[l, s] for l in incident_to[b]; init=0.0) +
        (b == base.j0 ? pi0[s] : 0.0) == 0)

    # The same dispatch on the true network. This makes a failed radial-safety
    # assumption visible while constructing the start instead of inside the MIP.
    @constraint(m, [l=1:Ln, s=1:S],
        tf[l, s] == base.Dx[l] * (tth[u[l], s] - tth[v[l], s]))
    @constraint(m, [l=1:Ln, s=1:S], -F[l] <= tf[l, s] <= F[l])
    @constraint(m, [s=1:S], tth[base.j0, s] == 0)
    @constraint(m, [b=1:N, s=1:S],
        sum(tf[l, s] for l in incident_from[b]; init=0.0) -
        sum(tf[l, s] for l in incident_to[b]; init=0.0) ==
        sum(g[k, s] for k in gens_at[b]; init=0.0) - demand[s][b])

    # For a fixed topology this is an ordinary LP primal-dual equality. It lets a
    # continuous LP construct every KKT multiplier without branching on SOS1.
    @constraint(m, [s=1:S],
        sum(base.c1[k] * g[k, s] for k in 1:K) ==
        sum(demand[s][b] * lam[b, s] for b in 1:N) -
        sum(F[l] * (mup[l, s] + mum[l, s]) for l in 1:Ln) -
        sum(base.pmax[k] * rhop[k, s] for k in 1:K) +
        sum(pmin[k] * rhom[k, s] for k in 1:K))
    @objective(m, Min,
        sum(mup) + sum(mum) + sum(rhop) + sum(rhom))
    optimize!(m)
    start_status = termination_status(m)
    isnothing(status_out) || (status_out[] = start_status)
    start_status == MOI.OPTIMAL || return nothing

    mupv, mumv = value.(mup), value.(mum)
    nuv = value.(nu)
    # Map the internal f=0 equality multiplier into the scaled-bound model's
    # nonnegative pair. Externally, select the ray-free representative.
    for l in 1:Ln, s in 1:S
        if internal[l]
            mupv[l, s] = max(nuv[l, s], 0.0)
            mumv[l, s] = max(-nuv[l, s], 0.0)
        else
            shift = min(mupv[l, s], mumv[l, s])
            mupv[l, s] -= shift
            mumv[l, s] -= shift
        end
    end
    return (g=value.(g), th=value.(th), f=value.(f), t=value.(t),
            lam=value.(lam), sig=value.(sig), pi0=value.(pi0),
            mup=mupv, mum=mumv, tau=value.(tau),
            rhop=value.(rhop), rhom=value.(rhom),
            tth=value.(tth), tf=value.(tf))
end

all_external_kkt_start(base, demand, pmin) =
    fixed_topology_kkt_start(base, demand, pmin, falses(base.Ln))

"""Extreme value of `a'g` over the T1 dispatch set, by continuous knapsack.

The set is `{g : pmin <= g <= pmax, sum(g) = total}` -- a box plus one equality,
so the extreme is reached by loading the most (or least) favourable generators
first. Returns `nothing` when the set is empty.
"""
function _extreme_over_dispatch(a::AbstractVector, pmin, pmax, total::Real;
                                maximize::Bool)
    K = length(a)
    headroom = pmax .- pmin
    remaining = total - sum(pmin)
    (remaining < -1e-9 || remaining > sum(headroom) + 1e-9) && return nothing
    remaining = max(remaining, 0.0)
    value = sum(a[k] * pmin[k] for k in 1:K; init=0.0)
    for k in sortperm(a; rev=maximize)
        remaining <= 0 && break
        take = min(headroom[k], remaining)
        value += a[k] * take
        remaining -= take
    end
    return value
end

"""T1 thermal limits that no admissible dispatch can reach.

T1 is written on the FULL network, which never changes with `c`, so a bound
computed once there is valid for every topology the upper level might choose.
The lower level is NOT screenable this way: collapsing a line redistributes
flow, so a full-network bound proves nothing about the reduced network.

The bound relaxes T1's dispatch to its box and the power balance, dropping the
requirement that `g` is what the reduced OPF returns. A relaxation can only
overstate the reachable flow, so a line it clears is genuinely unreachable and
its limit rows are redundant.

Returns an `Ln x S` mask of the (line, scenario) pairs whose rows can go, or
`nothing` if no screening is possible.
"""
function t1_redundant_limits(base, demand, pmin; margin::Float64=0.0)
    N, Ln = base.N, base.Ln
    K, S = length(base.gen_bus), length(demand)
    free = setdiff(1:N, [base.j0])
    pos = zeros(Int, N)
    for (i, b) in enumerate(free)
        pos[b] = i
    end
    nf = length(free)
    # Reduced Laplacian and the per-line susceptance rows, built directly: this
    # module takes `base` as data and does not include the preprocessing.
    B = zeros(nf, nf)
    W = zeros(Ln, nf)
    for l in 1:Ln
        ia, ib, y = pos[base.Efrom[l]], pos[base.Eto[l]], base.Dx[l]
        ia > 0 && (B[ia, ia] += y; W[l, ia] += y)
        ib > 0 && (B[ib, ib] += y; W[l, ib] -= y)
        if ia > 0 && ib > 0
            B[ia, ib] -= y
            B[ib, ia] -= y
        end
    end
    factor = try
        factorize(B)
    catch
        return nothing
    end
    P = zeros(Ln, N)
    P[:, free] = transpose(factor \ transpose(W))

    redundant = falses(Ln, S)
    a = zeros(K)
    for s in 1:S
        total = sum(demand[s])
        for l in 1:Ln
            for k in 1:K
                a[k] = P[l, base.gen_bus[k]]
            end
            offset = dot(view(P, l, :), demand[s])
            hi = _extreme_over_dispatch(a, pmin, base.pmax, total; maximize=true)
            lo = _extreme_over_dispatch(a, pmin, base.pmax, total; maximize=false)
            (isnothing(hi) || isnothing(lo)) && return nothing
            limit = base.frate[l] * (1 - margin)
            redundant[l, s] = (hi - offset <= limit) && (lo - offset >= -limit)
        end
    end
    return redundant
end

"Clusters induced by the internal lines: union-find over {l : cl_l = 1}."
function clustering_from_internal(base, internal::AbstractVector{Bool})
    parent = collect(1:base.N)
    find(x) = (parent[x] == x ? x : (parent[x] = find(parent[x])))
    for l in 1:base.Ln
        internal[l] || continue
        a, b = find(base.Efrom[l]), find(base.Eto[l])
        a == b || (parent[a] = b)
    end
    rep = [find(b) for b in 1:base.N]
    retained = sort(unique(rep))
    # Canonical representative = the smallest bus index in the cluster, matching
    # the parent code's extract_reduction so the two are directly comparable.
    canon = Dict(r => minimum(b for b in 1:base.N if rep[b] == r) for r in retained)
    rep_of = [canon[rep[b]] for b in 1:base.N]
    A = zeros(Int, base.N, base.N)
    for b in 1:base.N
        A[rep_of[b], b] = 1
    end
    return (rep_of=rep_of, retained=sort(unique(rep_of)), A=A)
end

"""
    solve_bilevel_reduction(c, scenario_indices; ...)

Maximize the number of collapsed lines subject to: for every modelled scenario,
some reduced-network optimal dispatch is feasible on the true full network
(optimistic T1). A finite `cost_gap_pct` additionally enforces the economic-gap
condition (T2); pass `nothing` for the feasibility-only target.

`kkt_form=:standard` reproduces the original SOS1 complementarity system.
`kkt_form=:normalized` adds the exact normalization that at most one of the
two thermal multipliers is nonzero, removing their collapsed-line common ray.

`radial_mode=:enforce` contracts every generator-free radial subtree to a
fixpoint and fixes those safe lines internal. Use `:warm` to use the same
reduction only as a start, or `:none` to disable radial preprocessing.

`warm_internal` supplies a topology seed without fixing its lines; in contrast,
`held_internal` fixes selected lines internal. An infeasible warm seed falls
back to the topology implied by the required fixed lines and radial mode.

`cycle_cut_lens=(2,3,4)` adds topology-closure cuts. They preserve the follower
OPF of a partition, but can restrict attainable partitions in combination with
line budgets, protected lines, or caps on selected paths. Setting
`lmp_separation=true` adds the normal shortest-path LMP rows computed from the
full-network OPFs. LMP separation is a selectable heuristic restriction and
can reduce the best attainable reduction.

`c` is a MultiScenarioTxReductionCase. Returns the clustering plus the pieces
needed to hand it to the parent directory's validators unchanged.
"""
function solve_bilevel_reduction(c, scenario_indices;
                                 cost_gap_pct::Union{Nothing,Real}=nothing,
                                 relax_pmin::Bool=true,
                                 protected=nothing,
                                 held_internal=nothing,
                                 warm_internal=nothing,
                                 neighborhood_center=nothing,
                                 neighborhood_radius=nothing,
                                 line_budget=nothing,
                                 time_limit=nothing,
                                 mipgap::Real=1e-4,
                                 numeric_focus=3,
                                 solver_threads=nothing,
                                 solver_seed=nothing,
                                 output_flag=nothing,
                                 start_time_limit=nothing,
                                 log_file=nothing,
                                 kkt_form::Symbol=:standard,
                                 complete_warm_start::Bool=true,
                                 radial_mode::Symbol=:enforce,
                                 objective::Symbol=:lines,
                                 t1_screening::Bool=true,
                                 t1_screen_margin::Float64=0.0,
                                 cycle_cut_lens=(),
                                 hop_cap=nothing,
                                 hop_cap_mode::Symbol=:auto,
                                 elec_cap=nothing,
                                 cluster_size_cap=nothing,
                                 lmp_separation::Bool=false,
                                 lmp_threshold::Real=3.0,
                                 lmp_relax_pmin::Bool=true,
                                 lmp_opf_time_limit=nothing)
    base = c.base
    N, Ln = base.N, base.Ln
    K = length(base.gen_bus)
    selected = Int.(collect(scenario_indices))
    S = length(selected)
    S > 0 || error("scenario_indices must not be empty")
    length(unique(selected)) == S || error("scenario_indices must be unique")
    all((1 .<= selected) .& (selected .<= size(c.p, 2))) ||
        error("scenario_indices must lie in 1:$(size(c.p, 2))")
    u, v = base.Efrom, base.Eto
    F = base.frate
    pmin = relax_pmin ? min.(0.0, base.pmin) : base.pmin

    isnothing(neighborhood_center) == isnothing(neighborhood_radius) ||
        error("neighborhood_center and neighborhood_radius must be supplied together")
    center = if isnothing(neighborhood_center)
        nothing
    else
        neighborhood_center isa AbstractVector{Bool} && length(neighborhood_center)==Ln ||
            error("neighborhood_center must be a Bool mask with $Ln entries")
        neighborhood_radius isa Integer && 0<=neighborhood_radius<=Ln ||
            error("neighborhood_radius must be an integer between 0 and $Ln")
        BitVector(neighborhood_center)
    end

    kkt_form in (:standard, :normalized) || error(
        "kkt_form must be :standard or :normalized, got $kkt_form")
    radial_mode in (:none, :warm, :enforce) || error(
        "radial_mode must be :none, :warm, or :enforce, got $radial_mode")
    objective in (:lines, :clusters) || error(
        "objective must be :lines or :clusters, got $objective")
    lmp_threshold >= 0 || error("lmp_threshold must be nonnegative")
    isfinite(t1_screen_margin) && t1_screen_margin >= 0 ||
        error("t1_screen_margin must be finite and nonnegative")
    isnothing(hop_cap) || (hop_cap isa Integer && hop_cap >= 1) ||
        error("hop_cap must be a positive integer edge count or nothing")
    hop_cap_mode in (:auto, :static, :lazy) ||
        error("hop_cap_mode must be :auto, :static or :lazy, got $hop_cap_mode")
    isnothing(elec_cap) || elec_cap > 0 ||
        error("elec_cap must be positive or nothing")
    isnothing(cluster_size_cap) ||
        (cluster_size_cap isa Integer && cluster_size_cap >= 2) ||
        error("cluster_size_cap must be an integer >= 2 or nothing")
    isnothing(cost_gap_pct) || cost_gap_pct >= 0 ||
        error("cost_gap_pct must be nonnegative or nothing")

    all(iszero, base.c2) || error(
        "Lower level assumes a LINEAR generation cost (c2 = 0); the runner " *
        "drops the quadratic term before the case reaches this point.")

    demand = [c.load[:, s] for s in selected]
    Zstar = if isnothing(cost_gap_pct)
        println(" Cost-gap constraint disabled (feasibility-only target)")
        nothing
    else
        env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
        z = [true_optimal_cost(base, d; relax_pmin=relax_pmin, env=env) for d in demand]
        println(" True full-network optimum per scenario = ", round.(z, digits=2))
        z
    end

    incident_from = [findall(==(b), u) for b in 1:N]
    incident_to   = [findall(==(b), v) for b in 1:N]
    gens_at       = [findall(==(b), base.gen_bus) for b in 1:N]

    m = Model(Gurobi.Optimizer)
    set_optimizer_attribute(m, "MIPGap", mipgap)
    isnothing(numeric_focus) || set_optimizer_attribute(m, "NumericFocus", numeric_focus)
    isnothing(time_limit) || set_optimizer_attribute(m, "TimeLimit", time_limit)
    isnothing(log_file) || set_optimizer_attribute(m, "LogFile", log_file)
    isnothing(solver_threads) || set_optimizer_attribute(m, "Threads", solver_threads)
    isnothing(solver_seed) || set_optimizer_attribute(m, "Seed", solver_seed)
    isnothing(output_flag) || set_optimizer_attribute(m, "OutputFlag", Int(output_flag))

    # ---------------------- upper level: c alone --------------------------- #
    @variable(m, cl[1:Ln], Bin)
    # el is pinned to 1 - cl by the row below, so cl's integrality already
    # makes it integral: declaring it Bin adds a binary per line that carries
    # no information. Presolve cannot substitute it out either, because it
    # sits in the SOS1 pairs, so the declaration is the only place to fix it.
    @variable(m, 0 <= el[1:Ln] <= 1)
    @constraint(m, [l=1:Ln], el[l] == 1 - cl[l])
    if !isnothing(center)
        @constraint(m,sum(center[l] ? 1-cl[l] : cl[l] for l in 1:Ln)<=neighborhood_radius)
        println(" Reversible neighborhood: at most ",neighborhood_radius," changed edge indicators")
    end
    protected_mask = if isnothing(protected)
        falses(Ln)
    elseif protected isa AbstractVector{Bool}
        length(protected) == Ln ||
            error("protected Bool mask must have $Ln entries")
        BitVector(protected)
    else
        ids = sort!(unique!(Int.(collect(protected))))
        all((1 .<= ids) .& (ids .<= Ln)) ||
            error("protected line IDs must lie in 1:$Ln")
        mask = falses(Ln)
        mask[ids] .= true
        mask
    end
    forbidden = copy(protected_mask)
    if any(protected_mask)
        for l in findall(protected_mask)
            forbidden[l] = true
            @constraint(m, cl[l] == 0)
        end
        println(" Protected lines pinned external = ", count(protected_mask),
                "   (optional here -- T1 checks safety directly)")
    end
    radial = radial_internal_mask(base; forbidden=forbidden)
    if radial_mode === :enforce
        @constraint(m, [l=radial.lines], cl[l] == 1)
    end
    println(" Safe generator-free radial lines = ", length(radial.lines),
            " / ", Ln, "   mode = ", radial_mode)

    # ------------------- ladder rung: held set and budget ------------------ #
    # Two changes and no others, exactly as iterative/ladder_core.jl does for
    # the proxy: some binaries are fixed to 1 because an earlier rung collapsed
    # them, and the cardinality cap moves up. Both leave every earlier rung's
    # solution feasible, so a rung's answer is always a feasible point of the
    # next one.
    held_mask = if isnothing(held_internal)
        falses(Ln)
    elseif held_internal isa AbstractVector{Bool}
        length(held_internal) == Ln ||
            error("held_internal Bool mask must have $Ln entries")
        BitVector(held_internal)
    else
        ids = sort!(unique!(Int.(collect(held_internal))))
        all((1 .<= ids) .& (ids .<= Ln)) ||
            error("held_internal line IDs must lie in 1:$Ln")
        mask = falses(Ln)
        mask[ids] .= true
        mask
    end
    clash = findall(held_mask .& protected_mask)
    isempty(clash) || error(
        "lines $clash are both protected external and held internal")
    if any(held_mask)
        @constraint(m, [l=findall(held_mask)], cl[l] == 1)
        println(" Held internal from earlier rungs = ", count(held_mask), " / ", Ln)
    end
    if !isnothing(line_budget)
        budget = floor(Int, line_budget)
        required_internal = copy(held_mask)
        radial_mode === :enforce && (required_internal .|= radial.internal)
        budget >= count(required_internal) || error(
            "line_budget $budget is below the $(count(required_internal)) " *
            "lines fixed by held_internal and enforced radial preprocessing")
        @constraint(m, sum(cl) <= budget)
        println(" Line budget = ", budget, " / ", Ln, "   (RR = ",
                round(budget / Ln, digits=3), ")")
    end

    # Topology closure. If every edge but one on a cycle is internal, the last
    # edge has endpoints in the same cluster and is redundant as an external
    # self-loop. These rows preserve the follower OPF of the closed partition,
    # but closing it also changes the selected-line count and can conflict with
    # a line budget, protected lines, or caps on selected paths.
    cycle_list = Vector{Vector{Int}}()
    n_cycle_cuts = 0
    if !isempty(cycle_cut_lens)
        cycles = case_helper(c, :short_cycles)(base; lens=cycle_cut_lens)
        for cyc in cycles
            push!(cycle_list, collect(cyc))
            k = length(cyc)
            for l in cyc
                @constraint(m,
                    cl[l] >= sum(cl[ll] for ll in cyc if ll != l) - (k - 2))
                n_cycle_cuts += 1
            end
        end
        println(" Short-cycle closure: lengths = ", cycle_cut_lens,
                "   cycles = ", length(cycle_list),
                "   rows added = ", n_cycle_cuts)
    end

    # ------------------ cluster spread: path caps -------------------------- #
    # Collapsing every line of a simple path puts its two ends in one cluster,
    # so capping the path caps how far apart two buses in a cluster can be.
    # These rows are written in cl alone -- no rank variable -- and are exact
    # under both objectives. The size cap below likewise uses cl witnesses. They
    # bite at the root, where the line budget provably does not: the relaxation
    # cannot set every cl to 1 along any capped path.
    #
    # What they enforce is one notch stronger than a diameter: no collapsed
    # SIMPLE PATH longer than the cap, which also rules out a long cycle inside
    # a cluster even though its diameter is half its length.
    path_rows = Vector{Vector{Int}}()
    n_hop_rows, n_elec_rows = 0, 0
    # :auto lists the paths up front when they fit (the rows tighten the root
    # LP, which is what made hop caps useful) and falls back to separating them
    # lazily when they do not -- a hop cap of 10 on a meshed network has far too
    # many paths to list. Same forbidden set either way.
    hop_lazy = false
    if !isnothing(hop_cap)
        ps = if hop_cap_mode === :lazy
            nothing
        else
            try
                capped_paths(base, float(hop_cap); forbidden=protected_mask)
            catch err
                hop_cap_mode === :auto && err isa ErrorException &&
                    occursin("max_paths", err.msg) || rethrow()
                nothing
            end
        end
        if isnothing(ps)
            hop_lazy = true
            println(" Hop cap = ", hop_cap, " edge(s)   lazy (paths separated on incumbents)")
        else
            append!(path_rows, ps)
            n_hop_rows = length(ps)
            println(" Hop cap = ", hop_cap, " edge(s)   rows = ", n_hop_rows)
        end
    end
    if !isnothing(elec_cap)
        # Merging pins theta_u = theta_v, and the angle it destroys across a
        # path is what sum(x_l) measures -- so this is the metric the cap is
        # really about, with hops as its unweighted stand-in.
        ps = capped_paths(base, float(elec_cap); weight=1 ./ base.Dx,
                           forbidden=protected_mask)
        append!(path_rows, ps)
        n_elec_rows = length(ps)
        println(" Electrical cap = ", elec_cap, " p.u.   rows = ", n_elec_rows)
    end
    for P in path_rows
        @constraint(m, sum(cl[l] for l in P) <= length(P) - 1)
    end

    # A pinned line cannot be released, so a cap its own set already breaks is
    # infeasible by construction rather than by physics -- and reads in the log
    # as the rung failing. Name the cause instead.
    if !isempty(path_rows) || !isnothing(cluster_size_cap) || hop_lazy
        pinned = copy(held_mask)
        radial_mode === :enforce && (pinned[radial.lines] .= true)
        bad = findfirst(P -> all(l -> pinned[l], P), path_rows)
        isnothing(bad) || error(
            "the pinned lines already collapse a path of " *
            "$(length(path_rows[bad])) edges $(path_rows[bad]), which the cap " *
            "forbids. Raise the cap, or set radial_mode = :warm.")
        if hop_lazy
            chain = long_chain(u, v, Float64.(pinned), 1:Ln, hop_cap)
            isempty(chain) || error(
                "the pinned lines already collapse a path of $(length(chain)) " *
                "edges $chain, which the hop cap forbids. Raise the cap, or set " *
                "radial_mode = :warm.")
        end
        if !isnothing(cluster_size_cap)
            pc = clustering_from_internal(base, pinned)
            big = maximum(count(==(r), pc.rep_of) for r in pc.retained)
            big <= cluster_size_cap || error(
                "the pinned lines already form a cluster of $big buses, above " *
                "the cap of $cluster_size_cap. Raise the cap, or set " *
                "radial_mode = :warm.")
        end
    end

    # The same baseline-LMP heuristic used by the regular TNR model. It is not
    # needed for T1 and is not an ideal-problem valid inequality: it deliberately
    # forbids some reductions to preserve economically distinct buses. Keep it
    # explicit and report it separately from the exact cycle cuts.
    lmp_paths = Vector{Vector{Int}}()
    lmp_pairs = Tuple{Int,Int}[]
    lmp_max_gap = 0.0
    n_lmp_skipped_radial = 0
    if lmp_separation
        lmp = case_helper(c, :full_network_lmps)(c, selected;
            relax_pmin=lmp_relax_pmin, time_limit=lmp_opf_time_limit)
        sep = case_helper(c, :lmp_separation_paths)(base, lmp;
            lmp_threshold=lmp_threshold, protected=forbidden)
        lmp_max_gap = sep.max_gap
        for (p, pair) in zip(sep.paths, sep.pairs)
            # Guaranteed radial elimination has priority over a baseline-price
            # heuristic. This only matters for a price-separated pair whose
            # complete chosen path lies inside the enforced radial forest.
            if radial_mode === :enforce && all(radial.internal[l] for l in p)
                n_lmp_skipped_radial += 1
                continue
            end
            push!(lmp_paths, p)
            push!(lmp_pairs, pair)
            @constraint(m, sum(cl[l] for l in p) <= length(p) - 1)
        end
        println(" LMP separation: threshold = ", lmp_threshold, " \$/MWh",
                "   max pair gap = ", round(lmp_max_gap, digits=2),
                "   violating pairs = ", sep.n_violating_pairs,
                "   rows added = ", length(lmp_paths),
                "   skipped for enforced radial = ", n_lmp_skipped_radial)
    end

    # ------------------- lower level LL(c, s): primal ---------------------- #
    # Written on the FULL bus set so that c enters only as a variable bound and
    # never as a coefficient. t is the intra-cluster transfer (gint upstream):
    # free on an internal line, zero on an external one, which is what makes a
    # cluster behave as a single node without changing the constraint matrix.
    @variable(m, g[1:K, 1:S])
    @variable(m, th[1:N, 1:S])
    @variable(m, f[1:Ln, 1:S])
    @variable(m, t[1:Ln, 1:S])

    # cl = 1 forces f = 0 (hence th_u = th_v: the buses are merged); cl = 0
    # leaves the plain thermal limit. One constraint pair does both, with no
    # big-M: F is the natural bound either way.
    @variable(m, slp[1:Ln, 1:S] >= 0)      # slack of  f <=  F(1-cl)
    @variable(m, slm[1:Ln, 1:S] >= 0)      # slack of -f <=  F(1-cl)
    @constraint(m, [l=1:Ln, s=1:S], slp[l, s] == F[l] * (1 - cl[l]) - f[l, s])
    @constraint(m, [l=1:Ln, s=1:S], slm[l, s] == F[l] * (1 - cl[l]) + f[l, s])
    @constraint(m, [l=1:Ln, s=1:S], [el[l], t[l, s]] in MOI.SOS1([1.0, 2.0]))

    @variable(m, sgp[1:K, 1:S] >= 0)       # slack of  g <= pmax
    @variable(m, sgm[1:K, 1:S] >= 0)       # slack of -g <= -pmin
    @constraint(m, [k=1:K, s=1:S], sgp[k, s] == base.pmax[k] - g[k, s])
    @constraint(m, [k=1:K, s=1:S], sgm[k, s] == g[k, s] - pmin[k])

    @constraint(m, flowdef[l=1:Ln, s=1:S],
        f[l, s] == base.Dx[l] * (th[u[l], s] - th[v[l], s]))
    @constraint(m, [s=1:S], th[base.j0, s] == 0)
    @constraint(m, balance[b=1:N, s=1:S],
        sum(f[l, s] + t[l, s] for l in incident_from[b]; init=0.0) -
        sum(f[l, s] + t[l, s] for l in incident_to[b];   init=0.0) ==
        sum(g[k, s] for k in gens_at[b]; init=0.0) - demand[s][b])

    # --------------------- lower level LL(c, s): dual ---------------------- #
    @variable(m, lam[1:N, 1:S])            # nodal price
    @variable(m, sig[1:Ln, 1:S])           # flow definition
    @variable(m, pi0[1:S])                 # angle reference
    @variable(m, mup[1:Ln, 1:S] >= 0)      # f <=  F(1-cl)
    @variable(m, mum[1:Ln, 1:S] >= 0)      # -f <= F(1-cl)
    @variable(m, tau[1:Ln, 1:S])           # t = 0 on an external line
    @variable(m, rhop[1:K, 1:S] >= 0)      # g <= pmax
    @variable(m, rhom[1:K, 1:S] >= 0)      # g >= pmin

    # Stationarity.
    @constraint(m, [k=1:K, s=1:S],
        base.c1[k] - lam[base.gen_bus[k], s] + rhop[k, s] - rhom[k, s] == 0)
    @constraint(m, [l=1:Ln, s=1:S],
        lam[u[l], s] - lam[v[l], s] + sig[l, s] + mup[l, s] - mum[l, s] == 0)
    # d/dt: on an internal line tau is forced to zero below, so this collapses to
    # lam_u == lam_v -- merging two buses equalizes their prices. The parent
    # model's lmp_separation assumes the contrapositive; here it is derived.
    @constraint(m, [l=1:Ln, s=1:S],
        lam[u[l], s] - lam[v[l], s] + tau[l, s] == 0)
    @constraint(m, [b=1:N, s=1:S],
        -sum(base.Dx[l] * sig[l, s] for l in incident_from[b]; init=0.0) +
         sum(base.Dx[l] * sig[l, s] for l in incident_to[b];   init=0.0) +
        (b == base.j0 ? pi0[s] : 0.0) == 0)

    # Complementary slackness. Every pair is SOS1, the same constraint type the
    # parent model is built from -- no big-M anywhere in this file.
    for l in 1:Ln, s in 1:S
        @constraint(m, [mup[l, s], slp[l, s]] in MOI.SOS1([1.0, 2.0]))
        @constraint(m, [mum[l, s], slm[l, s]] in MOI.SOS1([1.0, 2.0]))
        if kkt_form === :normalized
            # On a collapsed line both flow-bound slacks are zero, so the two
            # multipliers otherwise contain the cost-free ray
            # (mup,mum) -> (mup+k,mum+k).  At most one nonzero multiplier is
            # WLOG: every signed difference mup-mum has such a representation.
            # On an external line this is already implied by complementary
            # slackness whenever F > 0.  The extra SOS1 is therefore an exact
            # dual normalization, not an assumption about congestion.
            @constraint(m, [mup[l, s], mum[l, s]] in MOI.SOS1([1.0, 2.0]))
        end
        @constraint(m, [cl[l], tau[l, s]] in MOI.SOS1([1.0, 2.0]))
    end
    for k in 1:K, s in 1:S
        @constraint(m, [rhop[k, s], sgp[k, s]] in MOI.SOS1([1.0, 2.0]))
        @constraint(m, [rhom[k, s], sgm[k, s]] in MOI.SOS1([1.0, 2.0]))
    end

    # ------------- upper-level coupling on the lower solution -------------- #
    # T1: the dispatch g that the REDUCED network returns must be deliverable on
    # the TRUE network. Written in angle space (N variables, Ln rows, two
    # nonzeros each) rather than as dense PTDF rows -- algebraically identical,
    # and it keeps the model sparse.
    @variable(m, tth[1:N, 1:S])
    @variable(m, tf[1:Ln, 1:S])
    @constraint(m, [l=1:Ln, s=1:S],
        tf[l, s] == base.Dx[l] * (tth[u[l], s] - tth[v[l], s]))
    n_t1_screened = 0
    if t1_screening
        screened = t1_redundant_limits(base, demand, pmin; margin=t1_screen_margin)
        if isnothing(screened)
            println(" T1 screening: not applicable on this case; all rows kept")
            @constraint(m, [l=1:Ln, s=1:S], -F[l] <= tf[l, s] <= F[l])
        else
            for l in 1:Ln, s in 1:S
                screened[l, s] && continue
                @constraint(m, -F[l] <= tf[l, s] <= F[l])
            end
            n_t1_screened = count(screened)
            @printf(" T1 screening: %d of %d limit rows unreachable, dropped (%.1f%%)\n",
                    n_t1_screened, Ln * S, 100 * n_t1_screened / (Ln * S))
        end
    else
        @constraint(m, [l=1:Ln, s=1:S], -F[l] <= tf[l, s] <= F[l])
    end
    @constraint(m, [s=1:S], tth[base.j0, s] == 0)
    @constraint(m, [b=1:N, s=1:S],
        sum(tf[l, s] for l in incident_from[b]; init=0.0) -
        sum(tf[l, s] for l in incident_to[b];   init=0.0) ==
        sum(g[k, s] for k in gens_at[b]; init=0.0) - demand[s][b])

    # T2 is optional.  The feasibility-only target agreed in the main problem
    # uses T1 alone; a finite cost_gap_pct reproduces the earlier experiment.
    if !isnothing(cost_gap_pct)
        for s in 1:S
            @constraint(m,
                sum(base.c1[k] * g[k, s] + base.c0[k] for k in 1:K) <=
                Zstar[s] + cost_gap_pct / 100 * abs(Zstar[s]))
        end
    end

    # Objective.
    #
    # :lines    max sum(cl) -- the number of lines removed, the parent model's
    #           own objective, so the two differ only in the criterion.
    #
    # :clusters min the number of clusters. That count is a RANK, not a sum:
    #           #clusters = N - rank{l : cl_l = 1}, and an internal line joining
    #           two clusters removes a bus while one closing a cycle removes
    #           none. Continuous y auxiliaries bounded by cl and component ranks
    #           give max sum(y) = rank for each integral cl: a spanning forest
    #           attains the bound. The auxiliaries need not themselves describe
    #           a forest at every feasible or time-limited incumbent.
    #
    # Measured on case118 S=1: :clusters proves the bound (2 clusters, exactly
    # optimal) but its incumbent stalls around 95, where :lines walks to 2. Rank
    # is a plateau -- most single-line flips leave it unchanged -- so branching
    # gets no gradient. Keep both and compare; do not assume :clusters wins
    # because its objective is the one you care about.
    if objective === :clusters
        @variable(m, 0 <= y[1:Ln] <= 1)
        @constraint(m, [l=1:Ln], y[l] <= cl[l])
        # The whole-network rank row has to be present at the ROOT. Left to the
        # lazy callback, which only fires at integer nodes, the LP claims rank Ln
        # and the cluster bound pins at 1 for the entire solve.
        @constraint(m, sum(y) <= N - 1)
        if !isnothing(center)
            # Deleting center edges cannot increase graph rank. Each added edge
            # can increase it by at most one. This valid inequality also gives
            # the local relaxation useful information before integer callbacks.
            center_rank = N - length(clustering_from_internal(base, center).retained)
            @constraint(m, sum(y) <= center_rank + sum(cl[l] for l in 1:Ln if !center[l]))
        end
        # The short cycles are already enumerated for topology closure. Their
        # forest inequalities are free and, unlike the callback, bite at
        # fractional nodes. Valid because C is a subset of E(V(C)) and V(C)
        # induces a connected subgraph, so the row is a weakened rank inequality.
        for cyc in cycle_list
            @constraint(m, sum(y[l] for l in cyc) <= length(cyc) - 1)
        end
        if !isnothing(cluster_size_cap)
            # Clusters of at most Kmax buses need at least ceil(N/Kmax) of them,
            # so the rank cannot reach N-1. Static, exact, and it moves the root
            # bound -- the one thing the line budget never did.
            @constraint(m, sum(y) <= N - cld(N, cluster_size_cap))
            println(" Cluster size cap = ", cluster_size_cap,
                    " buses   rank bound = ", N - cld(N, cluster_size_cap))
        end

        # A connected internal topology is a copper plate whichever connected
        # edge set represents it: every angle is equal, every line carries zero
        # flow, and the follower's dispatch set is the same. So one KKT/T1 test
        # settles whether rank N-1 is reachable at all.
        #
        # The cap is sound only on a PROVEN infeasibility. A start that merely
        # failed to solve says nothing, and capping on it would cut off the true
        # optimum and then report the wrong answer as proven. INFEASIBLE_OR_
        # UNBOUNDED counts as proof here: that model minimises a sum of
        # nonnegative variables, so it is bounded below by 0 and cannot be
        # unbounded -- only the infeasible half of the status is reachable.
        one_cluster_status = Ref{Any}(nothing)
        one_cluster_start = fixed_topology_kkt_start(
            base, demand, pmin, trues(Ln); status_out=one_cluster_status,
            time_limit=start_time_limit, solver_threads=solver_threads)
        one_cluster_feasible = !isnothing(one_cluster_start)
        if one_cluster_feasible
            println(" One-cluster precheck: feasible")
        elseif one_cluster_status[] in (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED) &&
               N >= 2
            @constraint(m, sum(y) <= N - 2)
            println(" One-cluster precheck: proven infeasible; rank capped at ", N - 2)
        else
            println(" One-cluster precheck: INCONCLUSIVE (",
                    one_cluster_status[], "); no rank cap applied")
        end

        @objective(m, Max, sum(y))
    else
        @objective(m, Max, sum(cl))
    end

    # A binary-only start is not enough: Gurobi may fail to complete hundreds of
    # continuous SOS1 pairs for several minutes. Construct the full primal-dual
    # KKT point for the safe radial topology with a small LP instead.
    fallback_internal = radial_mode === :none ? falses(Ln) : copy(radial.internal)
    fallback_internal .|= held_mask
    start_internal = if !isnothing(warm_internal)
        seed = if warm_internal isa AbstractVector{Bool}
            length(warm_internal) == Ln || error("warm_internal Bool mask must have $Ln entries")
            BitVector(warm_internal)
        else
            ids = sort!(unique!(Int.(collect(warm_internal))))
            all((1 .<= ids) .& (ids .<= Ln)) || error("warm_internal line IDs must lie in 1:$Ln")
            mask = falses(Ln)
            mask[ids] .= true
            mask
        end
        seed[protected_mask] .= false
        seed .|= held_mask
        radial_mode === :enforce && (seed .|= radial.internal)
        seed
    elseif !isnothing(center)
        copy(center)
    elseif objective === :clusters && one_cluster_feasible && !any(protected_mask) &&
                        isempty(lmp_paths) && isnothing(line_budget) &&
                        isempty(path_rows) && isnothing(cluster_size_cap)
        # The copper plate IS the optimum when it is feasible, so hand it over.
        # Not under a budget: the copper plate blows it and the start is dropped.
        trues(Ln)
    else
        # A rung starts from what it is holding, which is the previous rung's
        # topology -- so the KKT start is built for a point known to be feasible.
        copy(fallback_internal)
    end
    # Repair a warm-only radial seed if a heuristic LMP row cuts it off. The
    # enforced case already omitted precisely those conflicting heuristic rows.
    for p in lmp_paths
        all(start_internal[l] for l in p) || continue
        j = findfirst(l -> !held_mask[l], p)
        isnothing(j) || (start_internal[p[j]] = false)
    end
    function topology_start_allowed(seed)
        isnothing(center) || count(seed .!= center)<=neighborhood_radius || return false
        isnothing(line_budget) || count(seed) <= floor(Int, line_budget) || return false
        any(P -> all(l -> seed[l], P), path_rows) && return false
        hop_lazy && !isempty(long_chain(u, v, Float64.(seed), 1:Ln, hop_cap)) && return false
        any(P -> all(l -> seed[l], P), lmp_paths) && return false
        any(cyc -> count(l -> seed[l], cyc) == length(cyc) - 1, cycle_list) && return false
        if !isnothing(cluster_size_cap)
            sc = clustering_from_internal(base, seed)
            maximum(count(==(r), sc.rep_of) for r in sc.retained) <= cluster_size_cap || return false
        end
        return true
    end
    topology_start_allowed(start_internal) || (start_internal = copy(fallback_internal))
    make_start(seed) = topology_start_allowed(seed) ?
        fixed_topology_kkt_start(base, demand, pmin, seed;
                                time_limit=start_time_limit, solver_threads=solver_threads) : nothing
    warm = complete_warm_start ? make_start(start_internal) : nothing
    if complete_warm_start && isnothing(warm) && start_internal != fallback_internal
        println(" Warm topology did not yield a KKT start; retrying required/radial topology")
        start_internal = copy(fallback_internal)
        warm = make_start(start_internal)
    end
    if complete_warm_start && isnothing(warm) && radial_mode === :warm &&
       start_internal != held_mask
        println(" Radial warm topology did not yield a KKT start; retrying held lines only")
        start_internal = copy(held_mask)
        warm = make_start(start_internal)
    end
    for l in 1:Ln
        cv = start_internal[l] ? 1.0 : 0.0
        set_start_value(cl[l], cv)
        set_start_value(el[l], 1.0 - cv)
    end
    if objective === :clusters
        # A spanning forest of the start topology, so the start satisfies every
        # rank row it will ever be asked about and carries the rank as its value.
        parent0 = collect(1:N)
        find0(x) = (parent0[x] == x ? x : (parent0[x] = find0(parent0[x])))
        for l in 1:Ln
            keep = false
            if start_internal[l]
                a, b = find0(u[l]), find0(v[l])
                if a != b
                    parent0[a] = b
                    keep = true
                end
            end
            set_start_value(y[l], keep ? 1.0 : 0.0)
        end
    end
    if !isnothing(warm)
        for s in 1:S
            set_start_value(pi0[s], warm.pi0[s])
            for b in 1:N
                set_start_value(th[b, s], warm.th[b, s])
                set_start_value(lam[b, s], warm.lam[b, s])
                set_start_value(tth[b, s], warm.tth[b, s])
            end
            for l in 1:Ln
                fv = warm.f[l, s]
                ev = start_internal[l] ? 0.0 : 1.0
                set_start_value(f[l, s], fv)
                set_start_value(t[l, s], warm.t[l, s])
                set_start_value(slp[l, s], F[l] * ev - fv)
                set_start_value(slm[l, s], F[l] * ev + fv)
                set_start_value(sig[l, s], warm.sig[l, s])
                set_start_value(mup[l, s], warm.mup[l, s])
                set_start_value(mum[l, s], warm.mum[l, s])
                set_start_value(tau[l, s], warm.tau[l, s])
                set_start_value(tf[l, s], warm.tf[l, s])
            end
            for k in 1:K
                gv = warm.g[k, s]
                set_start_value(g[k, s], gv)
                set_start_value(sgp[k, s], base.pmax[k] - gv)
                set_start_value(sgm[k, s], gv - pmin[k])
                set_start_value(rhop[k, s], warm.rhop[k, s])
                set_start_value(rhom[k, s], warm.rhom[k, s])
            end
        end
        println(" Complete primal-dual MIP start constructed (",
                count(start_internal), " internal lines)")
    else
        println(" No complete primal-dual start available; using binary start only")
    end

    # Component rank bounds for y, separated lazily. These rows do not describe
    # the entire forest polytope: a subcomponent cycle may still carry too much y.
    # They nevertheless bound sum(y) by the rank of integral cl, and a spanning
    # forest attains that value. Thus they give the correct optimization objective
    # without requiring full matroid separation. E(S) spans every line inside the
    # component, not only the selected ones: y <= cl zeroes the rest, and the
    # wider row is the stronger cut. The static rows above carry the fractional
    # nodes, which this callback cannot reach.
    hop_unchecked = Ref(0)
    if objective === :clusters || !isnothing(cluster_size_cap) || hop_lazy ||
            !isnothing(line_budget)
        set_optimizer_attribute(m, "LazyConstraints", 1)
        function forest_callback(cb_data)
            callback_node_status(cb_data, m) == MOI.CALLBACK_NODE_STATUS_INTEGER ||
                return
            clv = [callback_value(cb_data, cl[l]) for l in 1:Ln]
            yv  = objective === :clusters ?
                  [callback_value(cb_data, y[l]) for l in 1:Ln] : Float64[]
            parent = collect(1:N)
            find(x) = (parent[x] == x ? x : (parent[x] = find(parent[x])))
            for l in 1:Ln
                clv[l] > 0.5 || continue
                a, b = find(u[l]), find(v[l])
                a == b || (parent[a] = b)
            end
            rep = [find(b) for b in 1:N]
            size_of = zeros(Int, N)
            for b in 1:N
                size_of[rep[b]] += 1
            end
            lines_of = Dict{Int,Vector{Int}}()
            load_of = Dict{Int,Float64}()
            for l in 1:Ln
                rep[u[l]] == rep[v[l]] || continue
                r = rep[u[l]]
                push!(get!(lines_of, r, Int[]), l)
                objective === :clusters &&
                    (load_of[r] = get(load_of, r, 0.0) + yv[l])
            end
            # Under a budget, a line closing a loop inside a cluster must count.
            if !isnothing(line_budget)
                for (r, ls) in lines_of, l in ls
                    clv[l] > 0.5 && continue
                    P = merged_path(u, v, clv, ls, u[l], v[l])
                    isempty(P) && continue
                    MOI.submit(m, MOI.LazyConstraint(cb_data),
                        @build_constraint(cl[l] >= sum(cl[k] for k in P) - (length(P) - 1)))
                end
            end
            if objective === :clusters
                for (r, ls) in lines_of
                    n_i = size_of[r]
                    n_i >= 2 || continue
                    load_of[r] <= n_i - 1 + 1e-6 && continue
                    MOI.submit(m, MOI.LazyConstraint(cb_data),
                               @build_constraint(sum(y[l] for l in ls) <= n_i - 1))
                end
            end
            # Hop cap, lazily: one forbidden path per offending cluster. The
            # row is the static one, sum over a path of cap+1 lines <= cap.
            if hop_lazy
                for (r, ls) in lines_of
                    length(ls) > hop_cap || continue
                    ran_out = Ref(false)
                    P = long_chain(u, v, clv, ls, hop_cap; exhausted=ran_out)
                    ran_out[] && (hop_unchecked[] += 1)
                    isempty(P) && continue
                    MOI.submit(m, MOI.LazyConstraint(cb_data),
                        @build_constraint(sum(cl[l] for l in P) <= length(P) - 1))
                end
            end
            # Oversized cluster: cut the actual selected topology under BOTH
            # objectives. Restricting only y would let the same oversized cl
            # component survive with a smaller auxiliary forest. Any connected
            # Kmax+1-vertex witness must lose a selected tree edge.
            isnothing(cluster_size_cap) && return
            for (r, ls) in lines_of
                size_of[r] > cluster_size_cap || continue
                W, T = size_witness(u, v, clv, ls, cluster_size_cap)
                isempty(T) && continue
                MOI.submit(m, MOI.LazyConstraint(cb_data),
                    @build_constraint(sum(cl[l] for l in T) <=
                                      cluster_size_cap - 1))
            end
            return
        end
        MOI.set(m, MOI.LazyConstraintCallback(), forest_callback)
        objective === :clusters || isnothing(cluster_size_cap) || println(
            " Cluster size cap = ", cluster_size_cap,
            " buses   (:lines has no rank variable, so lazy cl cuts only --",
            " no static bound)")
    end

    println(" Bilevel MILP [", kkt_form, "]: ", num_variables(m), " variables, ",
            S, " scenario(s), ", N, " buses, ", Ln, " lines, ", K, " generators")
    t0 = time()
    optimize!(m)
    solve_time = time() - t0

    status = termination_status(m)
    nodes = try
        MOI.get(m, Gurobi.ModelAttribute("NodeCount"))
    catch
        NaN
    end
    work = try
        MOI.get(m, Gurobi.ModelAttribute("Work"))
    catch
        NaN
    end
    bound_value = try
        objective_bound(m)
    catch
        Inf # A solve stopped before its root LP may have no bound attribute.
    end
    has_values(m) || return (status=status, solve_time=solve_time, feasible=false,
                             neighborhood_center=center,neighborhood_radius=neighborhood_radius,
                             bound_scope=isnothing(center) ? :full_problem : :neighborhood,
                             kkt_form=kkt_form, node_count=nodes, work=work,
                             objective_kind=objective,
                             n_path_rows=length(path_rows),
                             n_hop_rows=n_hop_rows, n_elec_rows=n_elec_rows,
                             n_variables=num_variables(m),
                             bound=bound_value,
                             scenario_indices=selected,
                             scenario_ids=c.scenario_ids[selected],
                             relax_pmin=relax_pmin,
                             cost_gap_pct=cost_gap_pct,
                             protected=protected_mask)

    cv, evv = value.(cl), value.(el)
    gv, thv, fv, tv = value.(g), value.(th), value.(f), value.(t)
    slpv, slmv = value.(slp), value.(slm)
    sgpv, sgmv = value.(sgp), value.(sgm)
    lamv, sigv, pi0v = value.(lam), value.(sig), value.(pi0)
    mupv, mumv, tauv = value.(mup), value.(mum), value.(tau)
    rhopv, rhomv = value.(rhop), value.(rhom)
    tthv, tfv = value.(tth), value.(tf)

    internal = cv .> 0.5
    clus = clustering_from_internal(base, internal)
    lmp_merged_pairs = [(i, j) for (i, j) in lmp_pairs
                        if clus.rep_of[i] == clus.rep_of[j]]
    return (status            = status,
            neighborhood_center=center, neighborhood_radius=neighborhood_radius,
            bound_scope=isnothing(center) ? :full_problem : :neighborhood,
            feasible          = true,
            solve_time        = solve_time,
            objective_kind    = objective,
            objective         = objective_value(m),
            bound             = bound_value,
            # The cluster count the bound stands for, whichever objective ran.
            # For :clusters the MILP maximises a rank, so N - bound is the best
            # cluster count still reachable.
            n_clusters_bound  = objective === :clusters ?
                                max(N - floor(bound_value + 1e-7), 1.0) : NaN,
            node_count        = nodes,
            work              = work,
            n_variables       = num_variables(m),
            n_sos1            = S * (4 * Ln + 2 * K +
                                      (kkt_form === :normalized ? Ln : 0)),
            internal          = internal,
            n_internal_lines  = count(internal),
            n_external_lines  = count(!, internal),
            rep_of            = clus.rep_of,
            retained          = clus.retained,
            n_retained        = length(clus.retained),
            A                 = clus.A,
            dispatch          = gv,
            scenario_indices  = selected,
            scenario_ids      = c.scenario_ids[selected],
            Zstar             = Zstar,
            kkt_form          = kkt_form,
            relax_pmin        = relax_pmin,
            cost_gap_pct      = cost_gap_pct,
            n_t1_screened     = n_t1_screened,
            cost_gap_enforced = !isnothing(cost_gap_pct),
            protected         = protected_mask,
            held_internal     = held_mask,
            line_budget       = line_budget,
            path_rows         = path_rows,
            radial_mode       = radial_mode,
            radial_lines      = radial.lines,
            n_radial_lines    = length(radial.lines),
            hop_cap           = hop_cap,
            hop_cap_lazy      = hop_lazy,
            n_hop_unchecked   = hop_unchecked[],
            elec_cap          = elec_cap,
            cluster_size_cap  = cluster_size_cap,
            n_path_rows       = length(path_rows),
            n_hop_rows        = n_hop_rows,
            n_elec_rows       = n_elec_rows,
            max_cluster_size  = maximum(count(==(r), clus.rep_of)
                                        for r in clus.retained),
            cycle_cut_lens    = cycle_cut_lens,
            cycles            = cycle_list,
            n_cycles          = length(cycle_list),
            n_cycle_cuts      = n_cycle_cuts,
            lmp_separation    = lmp_separation,
            lmp_threshold     = lmp_threshold,
            lmp_paths         = lmp_paths,
            lmp_pairs         = lmp_pairs,
            n_lmp_rows        = length(lmp_paths),
            n_lmp_separated   = length(lmp_pairs) - length(lmp_merged_pairs),
            n_lmp_merged      = length(lmp_merged_pairs),
            lmp_merged_pairs  = lmp_merged_pairs,
            n_lmp_skipped_radial = n_lmp_skipped_radial,
            lmp_max_gap       = lmp_max_gap,
            # Raw point for the bilevel-specific post-solve residual audit.
            c_raw             = cv,
            e_raw             = evv,
            theta             = thv,
            reduced_flow      = fv,
            internal_transfer = tv,
            slack_flow_upper  = slpv,
            slack_flow_lower  = slmv,
            slack_gen_upper   = sgpv,
            slack_gen_lower   = sgmv,
            lambda            = lamv,
            sigma             = sigv,
            pi0               = pi0v,
            mu_upper          = mupv,
            mu_lower          = mumv,
            tau               = tauv,
            rho_upper         = rhopv,
            rho_lower         = rhomv,
            true_theta        = tthv,
            true_flow         = tfv,
            pmin_used         = pmin,
            cost              = [sum(base.c1[k] * gv[k, s] + base.c0[k]
                                     for k in 1:K) for s in 1:S])
end

end # module

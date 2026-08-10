# --------------------------------------------------------------------------- #
# PREPROCESSING -- everything that turns raw inputs into a solvable case.
#
# Every function here works for BOTH the single-scenario and the multi-scenario
# path: a single-scenario case is simply a MultiScenarioTxReductionCase with
# S = 1 (see `build_single_scenario_case`), so there is one code path, one set
# of windows, and one internal-transfer bound for both runners.
#
# Load order: this file FIRST -- it defines the two case structs that every
# other file's methods dispatch on. Only the runners call `include`.
# --------------------------------------------------------------------------- #

using PowerModels
PowerModels.silence()
using JuMP
using Gurobi
import MathOptInterface as MOI
using LinearAlgebra
using SparseArrays
using Graphs
using DelimitedFiles
using Dates
using Statistics

# =========================== 1. CASE DATA ================================== #
# --------------------------------------------------------------------------- #
# 1. Case data. Pure DC (b_l = 1/x_l, taps/shifts ignored -- the paper's Dx is a
# plain diagonal susceptance). Only in-service branches become lines l.
# --------------------------------------------------------------------------- #
struct TxReductionCase
    N::Int                      # buses
    Ln::Int                     # in-service lines
    Efrom::Vector{Int}          # u_l  (from-bus of line l), 1-based bus index
    Eto::Vector{Int}            # v_l  (to-bus of line l)
    Dx::Vector{Float64}         # diagonal susceptance b_l = 1/x_l  (Dx_ll)
    frate::Vector{Float64}      # thermal rating (>0, finite); f_lo=-frate, f_up=+frate
    # base operating point (consistent: p = E Dx E' thetahat, fhat = Dx E' thetahat)
    p::Vector{Float64}          # nodal injection gen-load (zero-sum)
    thetahat::Vector{Float64}
    fhat::Vector{Float64}
    Pd::Vector{Float64}
    gen_bus::Vector{Int}
    pmin::Vector{Float64}
    pmax::Vector{Float64}
    c2::Vector{Float64}
    c1::Vector{Float64}
    c0::Vector{Float64}
    j0::Int                     # reference bus
    baseMVA::Float64
end

# Signed node-line incidence E (N x Ln): +1 at from, -1 at to.
function incidence_matrix(c::TxReductionCase)
    I = vcat(c.Efrom, c.Eto)
    J = vcat(1:c.Ln, 1:c.Ln)
    V = vcat(fill(1.0, c.Ln), fill(-1.0, c.Ln))
    return sparse(I, J, V, c.N, c.Ln)
end


# Parse a pglib .m case and build a consistent DC base operating point via a
# DC-OPF (min gen cost). p is taken as the dispatch injection so p = E fhat holds
# exactly; thetahat/fhat come from the same solve.
function build_tx_case(casefile::AbstractString;
                       big_rate::Float64=1e3, time_limit=nothing)
    raw = PowerModels.parse_file(casefile)
    data = make_basic_network(raw)
    baseMVA = data["baseMVA"]

    buses = data["bus"]
    bkeys = sort(collect(keys(buses)), by = x -> parse(Int, x))
    bus_i = [buses[k]["bus_i"] for k in bkeys]
    bus_type = [buses[k]["bus_type"] for k in bkeys]
    N = length(bkeys)
    idx_of(b) = findfirst(==(b), bus_i)      # external bus number -> 1..N

    # generators
    gens = data["gen"]
    gk = sort(collect(keys(gens)), by = x -> parse(Int, x))
    nG = length(gk)
    gen_bus = [idx_of(gens[k]["gen_bus"]) for k in gk]
    pmin = [Float64(gens[k]["pmin"]) for k in gk]
    pmax = [Float64(gens[k]["pmax"]) for k in gk]
    cost = [gens[k]["cost"] for k in gk]     # [c2,c1,c0] or shorter
    c2 = [length(cst) >= 3 ? Float64(cst[end-2]) : 0.0 for cst in cost]
    c1 = [length(cst) >= 2 ? Float64(cst[end-1]) : 0.0 for cst in cost]
    c0 = [length(cst) >= 1 ? Float64(cst[end]) : 0.0 for cst in cost]

    # loads -> per-bus demand
    Pd = zeros(N)
    for (_, ld) in data["load"]
        Pd[idx_of(ld["load_bus"])] += Float64(ld["pd"])
    end

    # in-service branches -> lines
    brs = data["branch"]
    brk = sort(collect(keys(brs)), by = x -> parse(Int, x))
    Efrom = Int[]; Eto = Int[]; Dx = Float64[]; frate = Float64[]
    for k in brk
        br = brs[k]
        br["br_status"] == 0 && continue
        push!(Efrom, idx_of(br["f_bus"]))
        push!(Eto,   idx_of(br["t_bus"]))
        push!(Dx, 1.0 / br["br_x"])
        r = get(br, "rate_a", 0.0)
        push!(frate, (isfinite(r) && r > 0) ? Float64(r) : big_rate)
    end
    Ln = length(Efrom)
    j0 = something(findfirst(==(3), bus_type), 1)

    # base DC-OPF for a consistent (p, thetahat, fhat). OutputFlag is set on the
    # Env before it starts, not on the model after -- set_optimizer_attribute
    # here would be too late to suppress the license banner.
    env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
    m = Model(() -> Gurobi.Optimizer(env))
    set_optimizer_attribute(m, "MIPGap", 0.00)
    isnothing(time_limit) || set_optimizer_attribute(m, "TimeLimit", time_limit)
    @variable(m, th[1:N]); @variable(m, pG[1:nG]); @variable(m, pF[1:Ln])
    @constraint(m, th[j0] == 0)
    @constraint(m, [g=1:nG], pmin[g] <= pG[g] <= pmax[g])
    for l in 1:Ln
        @constraint(m, pF[l] == Dx[l] * (th[Efrom[l]] - th[Eto[l]]))
        @constraint(m, -frate[l] <= pF[l] <= frate[l])
    end
    for b in 1:N
        @constraint(m,
            sum(pG[g] for g in 1:nG if gen_bus[g] == b; init=0.0) - Pd[b] ==
            sum(pF[l] for l in 1:Ln if Efrom[l] == b; init=0.0) -
            sum(pF[l] for l in 1:Ln if Eto[l]   == b; init=0.0))
    end
    @objective(m, Min,
        sum(c2[g] * pG[g]^2 + c1[g] * pG[g] + c0[g] for g in 1:nG))
    optimize!(m)
    st = termination_status(m)
    st in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) ||
        error("Base DC-OPF failed on $casefile: $st (infeasible ratings?).")

    thetahat = value.(th)
    fhat = value.(pF)
    p = zeros(N)
    for g in 1:nG; p[gen_bus[g]] += value(pG[g]); end
    p .-= Pd

    return TxReductionCase(
        N, Ln, Efrom, Eto, Dx, frate,
        p, thetahat, fhat,
        Pd, gen_bus, pmin, pmax, c2, c1, c0,
        j0, baseMVA,
    )
end

struct MultiScenarioTxReductionCase
    base::TxReductionCase
    scenario_ids::Vector{Int}
    bus_ids::Vector{Int}
    load::Matrix{Float64}        # N x S, per unit
    generation::Matrix{Float64}  # N x S, per unit
    p::Matrix{Float64}           # generation - load, N x S, per unit
    thetahat::Matrix{Float64}    # N x S
    fhat::Matrix{Float64}        # Ln x S, per unit
end

"""
Build a ONE-SCENARIO case straight from a MATPOWER/pglib `.m` file.

This is what lets a plain test case (case118, case300, ...) run through exactly
the same model, benchmark and validation path as a 744-hour month: there is no
separate single-scenario formulation, only S = 1. No scenario matrices, no
`matrix_dir`, and no calendar are involved -- the single operating point is the
base DC-OPF that `build_tx_case` already solves.

Everything is read back out of that base point so the case is self-consistent by
construction (p = generation - load, and p = E*fhat exactly), which is the same
invariant `build_multiscenario_tx_case` checks after reconstructing flows from
saved matrices.
"""
function build_single_scenario_case(casefile::AbstractString;
                                    scenario_id::Int=1, time_limit=nothing)
    base = build_tx_case(casefile; time_limit=time_limit)

    # bus_ids is only carried for labelling (plots, CSV output). build_tx_case
    # indexes buses 1..N internally, so recover the external MATPOWER numbers in
    # the same sorted order it used.
    raw = PowerModels.parse_file(casefile)
    basic = make_basic_network(raw)
    buses = basic["bus"]
    bus_keys = sort(collect(keys(buses)), by=x -> parse(Int, x))
    bus_ids = Int[buses[k]["bus_i"] for k in bus_keys]
    length(bus_ids) == base.N ||
        error("Bus count mismatch building a single-scenario case from $casefile")

    return MultiScenarioTxReductionCase(
        base, [scenario_id], bus_ids,
        reshape(copy(base.Pd), base.N, 1),              # load
        reshape(base.p .+ base.Pd, base.N, 1),          # generation = p + load
        reshape(copy(base.p), base.N, 1),
        reshape(copy(base.thetahat), base.N, 1),
        reshape(copy(base.fhat), base.Ln, 1),
    )
end

"""
Scale every operating point consistently. A factor of 0.8 reduces load,
generation, nodal injections, DC angles, and reference line flows by 20%.
Topology, generator limits/costs, and line ratings are unchanged.
"""
function scale_operating_points(c::MultiScenarioTxReductionCase, factor::Real)
    isfinite(factor) && factor > 0 ||
        error("operating-point scale must be positive and finite")
    alpha = Float64(factor)
    return MultiScenarioTxReductionCase(
        c.base, copy(c.scenario_ids), copy(c.bus_ids),
        alpha .* c.load, alpha .* c.generation, alpha .* c.p,
        alpha .* c.thetahat, alpha .* c.fhat,
    )
end

"""Return a copy containing only the requested scenario columns."""
function subset_multiscenario_case(c::MultiScenarioTxReductionCase, scenario_indices)
    scenarios = sort!(unique!(Int.(collect(scenario_indices))))
    isempty(scenarios) && error("At least one scenario is required")
    all((1 .<= scenarios) .& (scenarios .<= length(c.scenario_ids))) ||
        error("scenario_indices must lie in 1:$(length(c.scenario_ids))")
    return MultiScenarioTxReductionCase(
        c.base, copy(c.scenario_ids[scenarios]), copy(c.bus_ids),
        copy(c.load[:, scenarios]), copy(c.generation[:, scenarios]),
        copy(c.p[:, scenarios]), copy(c.thetahat[:, scenarios]),
        copy(c.fhat[:, scenarios]),
    )
end

function _case_with_scaled_line_limits(c::MultiScenarioTxReductionCase, factor::Real)
    isfinite(factor) && factor > 0 ||
        error("line-limit scale must be positive and finite")
    alpha = Float64(factor)
    b = c.base
    scaled_base = TxReductionCase(
        b.N, b.Ln, copy(b.Efrom), copy(b.Eto), copy(b.Dx), alpha .* b.frate,
        copy(b.p), copy(b.thetahat), copy(b.fhat), copy(b.Pd),
        copy(b.gen_bus), copy(b.pmin), copy(b.pmax),
        copy(b.c2), copy(b.c1), copy(b.c0), b.j0, b.baseMVA,
    )
    return MultiScenarioTxReductionCase(
        scaled_base, copy(c.scenario_ids), copy(c.bus_ids),
        copy(c.load), copy(c.generation), copy(c.p),
        copy(c.thetahat), copy(c.fhat),
    )
end

"""
Select EVERY supplied scenario, in the same shape `select_seed_scenarios`
returns.

With one operating point (or any horizon small enough to solve on in full)
there is nothing to choose: the active set is all of it, so the set-cover MILP
that `select_seed_scenarios` runs would be a no-op. This produces the same
descriptor -- utilization, congested-line union, per-hour congestion counts --
so the reporting and file-output code cannot tell the two apart.

`minimum_cover_size` is reported as the number of scenarios taken, because
that IS the cover here: everything is included.
"""
function all_scenarios_selection(c::MultiScenarioTxReductionCase, scenario_indices;
                                 congestion_threshold::Float64=0.9999)
    idx = Int.(collect(scenario_indices))
    isempty(idx) && error("At least one scenario is required")
    base = c.base
    utilization = abs.(c.fhat[:, idx]) ./ base.frate
    congested = utilization .>= congestion_threshold
    return (
        scenario_indices=idx,
        scenario_ids=c.scenario_ids[idx],
        month_indices=idx,
        month_scenario_ids=c.scenario_ids[idx],
        congestion_threshold=congestion_threshold,
        congested_lines=findall(vec(any(congested, dims=2))),
        congested_mask=congested,
        utilization=utilization,
        congested_count_by_hour=vec(sum(congested, dims=1)),
        chosen_local_indices=collect(eachindex(idx)),
        minimum_cover_size=length(idx),
        selected_count=length(idx),
        ranking=:all,
        score=vec(sum(abs.(c.fhat[:, idx]), dims=1)),
    )
end

# ==================== 2. WINDOWS AND THE G BOUND =========================== #

# --------------------------------------------------------------------------- #
# 2. Windows and the data-only internal-transfer bound G (Lemma 1, eq 29).
#
# Fixed-tolerance window (eq 20-21): binding lines collapse to a point (Remark 2).
#   L+ (upper-binding): fhat == +frate ;  L- (lower-binding): fhat == -frate.
#
# `near_limit_threshold` (optional) extends the SAME protection to lines that
# are heavily loaded but not exactly at their rating: any line with
# |fhat| >= near_limit_threshold * frate is folded into L+/L- too, so it also
# gets forced external (17)/(43c) and, via phi_window, has its window pinned so
# it cannot look LESS loaded than it already is (only headroom AWAY from the
# limit is closed off; it can still move further toward the limit). This is a
# threshold on ORIGINAL loading, independent of the tol used to detect EXACT
# binding -- e.g. near_limit_threshold=0.9 additionally protects every line
# already carrying >= 90% of its rating.
function binding_lines(c::TxReductionCase; tol::Float64=1e-6, near_limit_threshold=nothing)
    Lplus  = [l for l in 1:c.Ln if abs(c.fhat[l] - c.frate[l])  <= tol * max(1.0, c.frate[l])]
    Lminus = [l for l in 1:c.Ln if abs(c.fhat[l] + c.frate[l])  <= tol * max(1.0, c.frate[l])]
    if !isnothing(near_limit_threshold)
        thr = near_limit_threshold
        Lplus  = sort(union(Lplus,  [l for l in 1:c.Ln if c.fhat[l]  >=  thr * c.frate[l]]))
        Lminus = sort(union(Lminus, [l for l in 1:c.Ln if c.fhat[l]  <= -thr * c.frate[l]]))
    end
    return Lplus, Lminus
end

# phi window (eq 20-21) for a given per-line tolerance eps_L. Single-scenario
# version of multiscenario_windows below -- see its docstring for what
# congestion_relaxation_mode does.
function phi_window(c::TxReductionCase, epsL::AbstractVector, Lplus, Lminus;
                    congestion_relaxation=0.0,
                    congestion_relaxation_mode::Symbol=:none)
    congestion_relaxation_mode in (:none, :conservative, :symmetric) ||
        error("congestion_relaxation_mode must be :none, :conservative or " *
              ":symmetric, got $congestion_relaxation_mode")
    relax = congestion_relaxation isa AbstractVector ?
        Vector{Float64}(congestion_relaxation) :
        fill(Float64(congestion_relaxation), c.Ln)
    length(relax) == c.Ln ||
        error("congestion_relaxation must be scalar or have $(c.Ln) entries")
    all(>=(0.0), relax) || error("congestion_relaxation must be nonnegative")
    delta = congestion_relaxation_mode === :none ? zeros(c.Ln) : relax .* c.frate

    philo = similar(c.frate); phiup = similar(c.frate)
    inp = Set(Lplus); inm = Set(Lminus)
    for l in 1:c.Ln
        in_p = l in inp
        in_m = l in inm
        # Only the congested SIDE may run past the rating, and only for a
        # protected line; everything else keeps its ordinary cap.
        lo_cap = -c.frate[l]
        hi_cap =  c.frate[l]
        in_p && (hi_cap += delta[l])
        in_m && (lo_cap -= delta[l])
        lo = max(lo_cap, c.fhat[l] - epsL[l])
        hi = min(hi_cap, c.fhat[l] + epsL[l])
        # keep congestion visible
        if in_p
            lo = max(lo, congestion_relaxation_mode === :symmetric ?
                         c.fhat[l] - delta[l] : c.fhat[l])
        end
        if in_m
            hi = min(hi, congestion_relaxation_mode === :symmetric ?
                         c.fhat[l] + delta[l] : c.fhat[l])
        end
        philo[l] = lo; phiup[l] = hi
    end
    return philo, phiup
end


"""
Scenario-specific flow windows and the union of protected lines.

`congestion_relaxation` trades base-point accuracy on CONGESTED lines for extra
reduction. It is a fraction of each line's rating (scalar, or a per-line vector)
and applies only to (line, scenario) pairs that are actually pinned -- exactly
binding, or above `near_limit_threshold`. `congestion_relaxation_mode` picks the
direction, and the direction is the whole argument:

  :none          the pin is exact (philo = fhat for a positively loaded binding
                 line). For an EXACTLY binding line the base window collapses to
                 the single point fhat, which pins an angle difference outright
                 and is the strongest single obstacle to merging.

  :conservative  the window opens only in the MORE-congested direction, to
                 [fhat, fhat + delta]. The reduced network may then over-state a
                 corridor's loading but never under-state it, so a reduced
                 dispatch can never be built on headroom that does not exist.
                 This necessarily lets the window exceed the thermal rating --
                 legitimate here because the reduction model is a FIT to the base
                 point, not a dispatch; the reduced DC-OPF still enforces the
                 true rating downstream.

  :symmetric     the window opens both ways, to [fhat - delta, fhat + delta].
                 Twice the freedom, but the model may now UNDER-state a congested
                 corridor by up to delta -- and an under-statement is precisely
                 what lets a reduced dispatch overload that line back on the full
                 network. Provided for comparison against :conservative, not as a
                 safe default.

A merely near-limit line already carries a window `eps` wide, so delta changes
almost nothing there. The relaxation is aimed at the exactly binding lines whose
window has collapsed to a point. `pinned` and `exact_binding` are returned so a
caller can report how many lines each mechanism actually reaches.
"""
function multiscenario_windows(c::MultiScenarioTxReductionCase, epsL;
                               near_limit_threshold=nothing,
                               scenario_indices=axes(c.p, 2),
                               protection_indices=axes(c.p, 2),
                               binding_tolerance::Float64=1e-6,
                               congestion_relaxation=0.0,
                               congestion_relaxation_mode::Symbol=:none)
    base, Ln = c.base, c.base.Ln
    if !isnothing(near_limit_threshold)
        0.0 < near_limit_threshold <= 1.0 ||
            error("near_limit_threshold must lie in (0, 1]")
    end
    congestion_relaxation_mode in (:none, :conservative, :symmetric) ||
        error("congestion_relaxation_mode must be :none, :conservative or " *
              ":symmetric, got $congestion_relaxation_mode")
    selected = Int.(collect(scenario_indices))
    protected_from = Int.(collect(protection_indices))
    epsv = epsL isa AbstractVector ? Vector{Float64}(epsL) : fill(Float64(epsL), Ln)
    length(epsv) == Ln || error("epsL must be scalar or have $Ln entries")

    relaxation = congestion_relaxation isa AbstractVector ?
        Vector{Float64}(congestion_relaxation) :
        fill(Float64(congestion_relaxation), Ln)
    length(relaxation) == Ln ||
        error("congestion_relaxation must be scalar or have $Ln entries")
    all(>=(0.0), relaxation) ||
        error("congestion_relaxation must be nonnegative")
    # delta is absolute (p.u.); the input is a fraction of each line's rating.
    delta = congestion_relaxation_mode === :none ? zeros(Ln) :
            relaxation .* base.frate

    protected = falses(Ln)
    for s in protected_from, l in 1:Ln
        flow = c.fhat[l, s]
        exact = abs(abs(flow) - base.frate[l]) <=
                binding_tolerance * max(1.0, base.frate[l])
        near = !isnothing(near_limit_threshold) &&
               abs(flow) >= near_limit_threshold * base.frate[l]
        protected[l] |= exact || near
    end

    philo = Matrix{Float64}(undef, Ln, length(selected))
    phiup = similar(philo)
    pinned = falses(Ln, length(selected))
    exact_binding = falses(Ln, length(selected))
    for (ss, s) in enumerate(selected), l in 1:Ln
        center = c.fhat[l, s]
        exact = abs(abs(center) - base.frate[l]) <=
                binding_tolerance * max(1.0, base.frate[l])
        near = !isnothing(near_limit_threshold) &&
               abs(center) >= near_limit_threshold * base.frate[l]
        # Pin only scenarios in which this line is actually protected. The
        # shared line is still forced external if protection occurs anywhere.
        is_pinned = exact || near
        pinned[l, ss] = is_pinned
        exact_binding[l, ss] = exact

        # Only the congested SIDE may run past the rating, and only for a pinned
        # pair -- an unpinned line keeps its ordinary [-frate, frate] cap.
        lo_cap = -base.frate[l]
        hi_cap =  base.frate[l]
        if is_pinned
            center >= 0 ? (hi_cap += delta[l]) : (lo_cap -= delta[l])
        end
        philo[l, ss] = max(lo_cap, center - epsv[l])
        phiup[l, ss] = min(hi_cap, center + epsv[l])
        if is_pinned
            if center >= 0
                floor_value = congestion_relaxation_mode === :symmetric ?
                    center - delta[l] : center
                philo[l, ss] = max(philo[l, ss], floor_value)
            else
                ceil_value = congestion_relaxation_mode === :symmetric ?
                    center + delta[l] : center
                phiup[l, ss] = min(phiup[l, ss], ceil_value)
            end
        end

        # The window must contain its own centre, or c = 0 (reduce nothing) isn't
        # feasible -- and it always must be. A pinned line's philo/phiup can land
        # exactly at center with zero slack, but f is only reproduced to solver
        # precision (~1e-5 p.u. after an OPF redispatch), so round-off alone can
        # make it infeasible (this is what broke ACTIVSg2000: 363 pinned pairs,
        # ~half rounding the wrong way). Floor the slack at 1e-5 relative --
        # three orders above the measured ~1.6e-6 p.u. reconstruction gap.
        pin_slack = max(binding_tolerance, 1e-5) * max(1.0, base.frate[l])
        philo[l, ss] = min(philo[l, ss], center - pin_slack)
        phiup[l, ss] = max(phiup[l, ss], center + pin_slack)
    end
    return (; selected, protected, philo, phiup, epsv, pinned, exact_binding,
            delta, congestion_relaxation_mode)
end

# Internal-transfer bound, per (line, scenario): min(x*H, G).
#   G  proven (0.5||p||_1 + sum F_m), but one number for the whole network --
#      loose, and numerically nasty since it spans orders of magnitude against
#      ~1e-4 constraint sensitivities.
#   H  per-line heuristic: |fhat_l| plus the `deg` largest tolerances eps_m in
#      the network (deg = max degree of l's endpoints) -- a line can only
#      absorb transfers through edges incident to it. Not proven; used alone
#      (x=1) it can silently under-bound a feasible clustering.
# x (internal_bound_scale) trades between them: large x -> falls back to G
# (safe, loose), x=3 default keeps most of H's benefit with headroom. The min
# with G means no x makes the bound looser than the proven one.
function multiscenario_internal_bounds(c::MultiScenarioTxReductionCase,
                                       selected, philo, phiup, epsv;
                                       internal_bound_scale::Real=3.0)
    internal_bound_scale > 0 ||
        error("internal_bound_scale must be positive, got $internal_bound_scale")
    base = c.base
    Ln = base.Ln
    deg = zeros(Int, base.N)
    for l in 1:Ln
        deg[base.Efrom[l]] += 1
        deg[base.Eto[l]] += 1
    end
    eps_desc = sort(collect(Float64, epsv); rev=true)
    S = length(selected)
    Gl = zeros(Ln, S)
    for (ss, s) in enumerate(selected)
        F = max.(abs.(philo[:, ss]), abs.(phiup[:, ss]))
        Gglobal = 0.5 * sum(abs, c.p[:, s]) + sum(F)
        for l in 1:Ln
            nbr = min(max(deg[base.Efrom[l]], deg[base.Eto[l]]), Ln)
            H = abs(c.fhat[l, s]) + sum(@view eps_desc[1:nbr])
            Gl[l, ss] = min(internal_bound_scale * H, Gglobal)
        end
    end
    return Gl
end

# DC power-flow consistency -- check this FIRST when a solve is infeasible. If
# (p, Dx, E, thetahat, fhat) satisfy the DC equations, c = 0 (reduce nothing)
# is feasible by construction, so an INFEASIBLE result means the operating
# point itself is broken, not that the search failed. Three residuals, each
# failing for a different reason: sum(p) != 0 is fatal (no c can fix it, since
# summing the nodal rows cancels every flow term); E*fhat != p means fhat and p
# came from different injections; fhat != Dx*E'*th means fhat and th came from
# different solves.
"""
    dc_consistency(c; scenario_indices=axes(c.p, 2))

Residuals of the three DC power-flow identities over the requested scenarios, in
p.u. Everything downstream assumes these are ~0; when a reduction is infeasible,
check this before anything else.
"""
function dc_consistency(c::MultiScenarioTxReductionCase;
                        scenario_indices=axes(c.p, 2))
    base = c.base
    s = Int.(collect(scenario_indices))
    E = incidence_matrix(base)
    P = c.p[:, s]
    Fh = c.fhat[:, s]
    Th = c.thetahat[:, s]
    global_imbalance = vec(sum(P, dims=1))
    nodal = P - E * Fh
    flowdef = Fh - Diagonal(base.Dx) * E' * Th
    return (
        scenarios=s,
        max_global_imbalance=maximum(abs, global_imbalance; init=0.0),
        max_nodal_residual=maximum(abs, nodal; init=0.0),
        max_flow_residual=maximum(abs, flowdef; init=0.0),
        worst_scenario=isempty(s) ? 0 : s[argmax(abs.(global_imbalance))],
        global_imbalance=global_imbalance,
    )
end

"""Print `dc_consistency`, in MW, with a verdict on whether c = 0 is feasible."""
function report_dc_consistency(c::MultiScenarioTxReductionCase, chk; tol::Float64=1e-9)
    mva = c.base.baseMVA
    println("\nDC power-flow consistency over ", length(chk.scenarios), " scenario(s)")
    println("   sum(p) = 0        max ", chk.max_global_imbalance * mva, " MW")
    println("   E fhat = p        max ", chk.max_nodal_residual * mva, " MW")
    println("   fhat = Dx E' th   max ", chk.max_flow_residual * mva, " MW")
    worst = max(chk.max_global_imbalance, chk.max_nodal_residual, chk.max_flow_residual)
    if worst <= tol
        println("   OK -- the data describes a DC power flow, so c = 0 is feasible ",
                "and the reduction model cannot be infeasible.")
    else
        println("   BROKEN -- worst residual ", worst * mva, " MW exceeds ", tol * mva,
                " MW.")
        chk.max_global_imbalance > tol && println(
            "   sum(p) != 0 is the fatal one: it makes the reduction model ",
            "infeasible for EVERY clustering, c = 0 included.")
    end
    return worst <= tol
end

# ============ 3. OPERATING-POINT REDISPATCH AND LIMIT SCALING ============== #
"""
Re-solve the full-network DC-OPF for the requested scenario columns and
atomically replace their generation, injection, angle, and flow data. Loads
and line ratings are taken from `c`, so this can be used after applying an
operating-point/load scale and a line-limit scale.
"""
function redispatch_dc_opf_scenarios(c::MultiScenarioTxReductionCase,
                                     scenario_indices=axes(c.p, 2);
                                     relax_pmin::Bool=true,
                                     time_limit=nothing,
                                     progress_every::Int=100,
                                     balance_tolerance::Float64=1e-6,
                                     rating_relative_tolerance::Float64=1e-3)
    scenarios = sort!(unique!(Int.(collect(scenario_indices))))
    isempty(scenarios) && error("At least one redispatch scenario is required")
    all((1 .<= scenarios) .& (scenarios .<= length(c.scenario_ids))) ||
        error("scenario_indices must lie in 1:$(length(c.scenario_ids))")
    progress_every > 0 || error("progress_every must be positive")
    balance_tolerance >= 0 || error("balance_tolerance must be nonnegative")
    rating_relative_tolerance >= 0 ||
        error("rating_relative_tolerance must be nonnegative")

    base = c.base
    identity_assignment = Matrix{Int}(I, base.N, base.N)
    env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
    optimizer = () -> Gurobi.Optimizer(env)

    # Work on copies and construct the new immutable case only after every OPF
    # succeeds, so a failed hour cannot leave a partially refreshed data set.
    generation = copy(c.generation)
    p = copy(c.p)
    thetahat = copy(c.thetahat)
    fhat = copy(c.fhat)
    started = time()
    println("Automatically redispatching $(length(scenarios)) scenario(s) under the adjusted line limits")

    for (k, s) in enumerate(scenarios)
        result = solve_assigned_dc_opf_scenario(
            c, identity_assignment, s;
            relax_pmin=relax_pmin,
            time_limit=time_limit,
            optimizer=optimizer,
        )
        result.optimal || error(
            "Automatic DC-OPF redispatch failed for scenario ID $(c.scenario_ids[s]) " *
            "with status $(result.status). The requested load/operating-point and " *
            "line-limit scales may be physically infeasible.")

        generation[:, s] .= 0.0
        for g in eachindex(base.gen_bus)
            generation[base.gen_bus[g], s] += result.pG[g]
        end
        p[:, s] .= generation[:, s] .- c.load[:, s]

        # sum(p) must be EXACTLY 0: summing every nodal balance row cancels the
        # flow terms and leaves sum_b p_b = 0, a hard requirement independent of
        # clustering -- if it's violated, even c = 0 is infeasible. A DC-OPF only
        # satisfies its own balance to solver tolerance (~1.2e-5 p.u. measured),
        # so push the residual onto the reference bus, whose balance row is
        # redundant anyway (theta_ref is pinned).
        imbalance = sum(view(p, :, s))
        p[base.j0, s] -= imbalance
        generation[base.j0, s] -= imbalance

        thetahat[:, s] .= result.theta
        fhat[:, s] .= result.flow

        if k == 1 || k % progress_every == 0 || k == length(scenarios)
            println("Automatic DC-OPF redispatch: $k/$(length(scenarios)) scenarios (",
                    round(time() - started, digits=1), " s)")
        end
    end

    # Re-derive angles/flows from the corrected p rather than trusting the OPF's
    # own theta/flow (only solver-tolerance accurate, and p just moved at the
    # reference bus). This makes (p, fhat, thetahat) self-consistent to ~1e-12,
    # which matters because a pinned line's window is one-sided at fhat -- any
    # gap between fhat and what p actually implies makes the model infeasible.
    E = incidence_matrix(base)
    let B = E * Diagonal(base.Dx) * E',
        freebus = setdiff(1:base.N, [base.j0])
        F = factorize(B[freebus, freebus])
        for s in scenarios
            th = zeros(base.N)
            th[freebus] = F \ p[freebus, s]
            thetahat[:, s] .= th
            fhat[:, s] .= Diagonal(base.Dx) * E' * th
        end
    end
    selected_p = p[:, scenarios]
    selected_theta = thetahat[:, scenarios]
    selected_flow = fhat[:, scenarios]
    max_system_imbalance = maximum(abs, vec(sum(selected_p, dims=1)))
    max_nodal_residual = maximum(abs, selected_p - E * selected_flow)
    max_flow_residual = maximum(abs,
        selected_flow - Diagonal(base.Dx) * E' * selected_theta)
    overload = max.(abs.(selected_flow) .- base.frate, 0.0)
    allowed_overload = balance_tolerance .+
        rating_relative_tolerance .* base.frate
    max_rating_excess = maximum(overload .- allowed_overload)

    # System imbalance sums over all N buses, so its round-off grows with system
    # size, unlike the per-bus/per-line residuals below (bounded by the solver's
    # FeasibilityTol, not accumulating) -- give only this one a size-relative
    # floor (a fixed 1e-6 is fine for ACTIVSg200 but 12x too tight for ACTIVSg2000).
    peak_demand = maximum(sum(view(c.load, :, scenarios), dims=1); init=0.0)
    system_tolerance = max(balance_tolerance, 1e-7 * peak_demand)
    max_system_imbalance <= system_tolerance || error(
        "Automatic redispatch power imbalance is $max_system_imbalance p.u., " *
        "over the $system_tolerance p.u. limit " *
        "($(round(100 * max_system_imbalance / max(peak_demand, eps()), digits=8))% of peak demand)")
    # nodal/flow residuals are differences of large flow terms (nodal_i =
    # p_i - sum(+/-f_l), flow_l = f_l - Dx_l(th_u-th_v)), so round-off scales
    # with flow magnitude, not bus count -- scale the tolerance by the largest
    # flow present. A genuine defect (mis-mapped bus, transposed incidence)
    # still stands out orders of magnitude above this.
    flow_scale = max(1.0, maximum(abs, selected_flow; init=0.0))
    residual_tolerance = max(balance_tolerance, 1e-6 * flow_scale)
    max_nodal_residual <= residual_tolerance || error(
        "Automatic redispatch nodal residual is $max_nodal_residual p.u., " *
        "over the $residual_tolerance p.u. limit (peak flow $flow_scale p.u.)")
    max_flow_residual <= residual_tolerance || error(
        "Automatic redispatch DC-flow residual is $max_flow_residual p.u., " *
        "over the $residual_tolerance p.u. limit (peak flow $flow_scale p.u.)")
    max_rating_excess <= 0.0 || error(
        "Automatic redispatch exceeds an adjusted line limit beyond tolerance by " *
        "$(base.baseMVA * max_rating_excess) MW")

    return MultiScenarioTxReductionCase(
        base, copy(c.scenario_ids), copy(c.bus_ids), copy(c.load),
        generation, p, thetahat, fhat,
    )
end

"""
Scale all line limits and automatically recompute economically dispatched,
limit-feasible operating points for the requested scenarios.
"""
function scale_line_limits_and_redispatch(
        c::MultiScenarioTxReductionCase, factor::Real;
        scenario_indices=axes(c.p, 2), kwargs...)
    scaled = _case_with_scaled_line_limits(c, factor)
    return redispatch_dc_opf_scenarios(
        scaled, scenario_indices; kwargs...)
end

# ================= 4. CYCLE ENUMERATION AND WARM START ===================== #
# Warm start: merge every radial line (a degree-1 endpoint in the original
# topology) that isn't congested. A degree-1 bus's whole injection already
# flows through its one line, so shorting it is free (f_l=0, g_l=fhat_l) and
# disturbs nothing else -- every other line keeps f = fhat.
#
# The representative of a merged (leaf, hub) pair must be the hub, not the
# smaller index: vartheta_warm is read off as thetahat[rep_warm[bus]], and the
# hub's OTHER lines need vartheta there to stay thetahat[hub]. Picking the leaf
# breaks those rows (surfaced as a real Gurobi constraint-violation report on
# the delivered MIP start). j0 is never treated as a leaf -- its angle is
# hard-fixed to 0.
function _radial_warm_start(c::TxReductionCase, Lplus, Lminus)
    deg = zeros(Int, c.N)
    for l in 1:c.Ln
        deg[c.Efrom[l]] += 1
        deg[c.Eto[l]]   += 1
    end
    congested = falses(c.Ln)
    for l in Lplus;  congested[l] = true; end
    for l in Lminus; congested[l] = true; end

    cl_warm = zeros(Int, c.Ln)
    rep_warm = collect(1:c.N)
    for l in 1:c.Ln
        congested[l] && continue
        u, v = c.Efrom[l], c.Eto[l]
        leaf_is_u = deg[u] == 1 && u != c.j0
        leaf_is_v = deg[v] == 1 && v != c.j0
        (leaf_is_u || leaf_is_v) || continue
        cl_warm[l] = 1
        if leaf_is_u && !leaf_is_v
            rep_warm[u] = v
        elseif leaf_is_v && !leaf_is_u
            rep_warm[v] = u
        else
            rep_warm[u] = v   # isolated (leaf,leaf) pair -- either choice is fine
        end
    end
    return cl_warm, rep_warm
end

# --------------------------------------------------------------------------- #
# EXACTLY-MERGEABLE LINES: bridges and protected-free leaf blocks.
#
# A bridge (u,v) is the only path between the two sides of the network, so its
# flow is fixed by one side's net injection; shorting it leaves every other
# line's flow unchanged, in every scenario -- zero window error. Same argument
# for a leaf block reached through a single articulation point. So for any
# unprotected line like this, c_l = 1 is a dominance fixing (free objective
# gain, zero accuracy cost) -- fixing it deletes search, it can't change the
# answer. (Measured: 9/186 lines on case118, 94/411 on case300, 71/245 on
# ACTIVSg200 March.) Presolve can't find this itself: it's a global
# connectivity + DC-physics property, not a local bound implication.
#
# SimpleGraph collapses parallel circuits, so a double-circuit corridor looks
# like one edge and Tarjan would call it a bridge -- wrong, since the pair
# stays connected if either circuit opens. Candidates are filtered by bus-pair
# multiplicity to avoid this (ACTIVSg200 has double circuits).
#
# One precondition is checked, not assumed: shorting l needs g^int to absorb
# its flow, i.e. |fhat_l| <= G_l.
# --------------------------------------------------------------------------- #
"""
    exactly_mergeable_lines(base, protected; include_leaf_blocks=true)

Line indices that can be made internal with provably zero flow error anywhere
else in the network. `protected` is a Boolean mask or index list of lines that
must stay external; those are never returned, nor is any leaf block containing
one.
"""
function exactly_mergeable_lines(base::TxReductionCase, protected;
                                 include_leaf_blocks::Bool=true)
    prot = Set(protected isa AbstractVector{Bool} ? findall(protected) :
               Int.(collect(protected)))
    N, Ln = base.N, base.Ln

    pair_lines = Dict{Tuple{Int,Int},Vector{Int}}()
    for l in 1:Ln
        a, b = minmax(base.Efrom[l], base.Eto[l])
        a == b && continue                       # ignore any self-loop
        push!(get!(pair_lines, (a, b), Int[]), l)
    end
    g = SimpleGraph(N)
    for ((a, b), _) in pair_lines
        add_edge!(g, a, b)
    end

    bridge_lines = Int[]
    for e in Graphs.bridges(g)
        a, b = minmax(src(e), dst(e))
        ls = get(pair_lines, (a, b), Int[])
        length(ls) == 1 || continue              # parallel pair -> NOT a bridge
        l = ls[1]
        l in prot || push!(bridge_lines, l)
    end

    leaf_lines = Int[]
    n_leaf_blocks = 0
    if include_leaf_blocks
        arts = Set(Graphs.articulation(g))
        for blk in Graphs.biconnected_components(g)
            verts = Set{Int}()
            for e in blk
                push!(verts, src(e)); push!(verts, dst(e))
            end
            length(intersect(verts, arts)) == 1 || continue   # leaf of block-cut tree
            lines = Int[]
            clean = true
            for e in blk
                a, b = minmax(src(e), dst(e))
                for l in get(pair_lines, (a, b), Int[])
                    l in prot && (clean = false)
                    push!(lines, l)
                end
            end
            (clean && length(lines) > 1) || continue          # 1 line = already a bridge
            n_leaf_blocks += 1
            append!(leaf_lines, lines)
        end
    end

    lines = sort!(unique!(vcat(bridge_lines, leaf_lines)))
    return (lines=lines, bridges=sort(bridge_lines),
            leaf_block_lines=sort(unique(leaf_lines)),
            n_leaf_blocks=n_leaf_blocks)
end

# --------------------------------------------------------------------------- #
# LMP SEPARATION. Two buses with materially different prices shouldn't end up
# in one cluster -- merging replaces two prices with one, and the reduction
# objective (count internal lines) has no reason to care about that itself.
#
# "i and j in different clusters" means every i-j path has an external line --
# exponentially many paths. Instead this imposes it on just the shortest one:
#     sum_{l in P(i,j)} c_l <= |P(i,j)| - 1
# Necessary, not sufficient: never excludes a genuinely separated clustering,
# but admits some that merge i,j via a different path (reported as the
# achieved-separation count). |P|=1 degenerates to c_l = 0.
#
# Minimality: on one BFS tree rooted at i, the path to any k on j's parent
# chain is a prefix of the path to j, so if (i,k) also violates, its row
# already implies j's and the (i,j) row is dropped. Only prefixes on the same
# tree qualify -- the mirrored "suffix" rule isn't sound, since the shortest
# k-to-j path need not be the tail of i's path to j.
# --------------------------------------------------------------------------- #
"""
    lmp_separation_paths(base, lmp; lmp_threshold, protected=Int[])

Shortest-path separation rows for every bus pair whose price gap exceeds
`lmp_threshold` in at least one scenario. `lmp` is bus x scenario, in \$/MWh.

Returns `paths` (each a vector of line indices, to be constrained to at most
`length - 1` internal lines), plus the pair and gap diagnostics.
"""
function lmp_separation_paths(base::TxReductionCase, lmp::AbstractMatrix;
                              lmp_threshold::Real, protected=Int[])
    N, Ln = base.N, base.Ln
    # Worst price gap each pair reaches ANYWHERE in the horizon. A pair that
    # separates in one hour must stay separated in the shared clustering, so the
    # max over scenarios -- not the mean -- is the protective aggregation.
    gap = zeros(N, N)
    for s in axes(lmp, 2), i in 1:N, j in i+1:N
        d = abs(lmp[i, s] - lmp[j, s])
        if d > gap[i, j]
            gap[i, j] = d
            gap[j, i] = d
        end
    end
    viol = gap .> lmp_threshold
    for i in 1:N
        viol[i, i] = false
    end
    # A line that is already pinned external satisfies any row it appears in, so
    # every path through one is redundant and is dropped rather than added.
    prot = Set(protected isa AbstractVector{Bool} ? findall(protected) :
               Int.(collect(protected)))

    adj = [Int[] for _ in 1:N]          # neighbour -> connecting line
    adjl = [Int[] for _ in 1:N]
    for l in 1:Ln
        a, b = base.Efrom[l], base.Eto[l]
        a == b && continue
        push!(adj[a], b); push!(adjl[a], l)
        push!(adj[b], a); push!(adjl[b], l)
    end

    seen = Set{Vector{Int}}()
    paths = Vector{Vector{Int}}()
    pairs = Tuple{Int,Int}[]
    parent = zeros(Int, N)
    pline = zeros(Int, N)
    depth = fill(-1, N)
    queue = Int[]
    for i in 1:N
        any(@view viol[i, :]) || continue
        fill!(parent, 0); fill!(pline, 0); fill!(depth, -1)
        empty!(queue); push!(queue, i); depth[i] = 0
        head = 1
        while head <= length(queue)
            b = queue[head]; head += 1
            for (t, nb) in enumerate(adj[b])
                depth[nb] == -1 || continue
                depth[nb] = depth[b] + 1
                parent[nb] = b
                pline[nb] = adjl[b][t]
                push!(queue, nb)
            end
        end
        for j in 1:N
            (viol[i, j] && depth[j] > 0) || continue
            # Prefix-minimality: a violating pair strictly inside the chain to j
            # already implies j's row on this same tree.
            k = parent[j]
            minimal = true
            while k != i && k != 0
                if viol[i, k]
                    minimal = false
                    break
                end
                k = parent[k]
            end
            minimal || continue
            path = Int[]
            node = j
            while node != i
                push!(path, pline[node])
                node = parent[node]
            end
            any(in(prot), path) && continue      # already satisfied by protection
            sort!(path)
            path in seen && continue
            push!(seen, path)
            push!(paths, path)
            push!(pairs, (i, j))
        end
    end
    return (paths=paths, pairs=pairs, gap=gap, n_violating_pairs=count(viol) ÷ 2,
            max_gap=maximum(gap; init=0.0))
end

"""Full-network LMPs (bus x selected scenario, \$/MWh) for the separation cuts."""
function full_network_lmps(c::MultiScenarioTxReductionCase, selected;
                           relax_pmin::Bool=true, time_limit=nothing)
    N = c.base.N
    lmp = zeros(N, length(selected))
    A0 = Matrix{Float64}(I, N, N)
    # One shared, silenced environment for every scenario -- each fresh Gurobi.Env
    # reprints the license banner (OutputFlag set after the fact doesn't suppress
    # it, only a shared env avoids repeating it).
    env = Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
    optimizer = () -> Gurobi.Optimizer(env)
    for (ss, s) in enumerate(selected)
        res = solve_assigned_dc_opf_scenario(c, A0, s;
                relax_pmin=relax_pmin, time_limit=time_limit, optimizer=optimizer)
        isnothing(res.lmp) &&
            error("full-network DC-OPF gave no duals for scenario $s ($(res.status))")
        lmp[:, ss] = res.lmp
    end
    return lmp
end

# Short cycles as line-index sets, for the closure cuts. Graphs.jl's
# cycle_basis only returns a basis of the cycle space, not every cycle, so
# filtering it to length 3/4 would miss most short cycles -- and it runs on a
# SimpleGraph, which collapses parallel lines (the strongest case: two lines on
# one bus pair force c_a = c_b). Enumerate directly instead. Multi-line bus
# pairs expand combinatorially: a triangle with 1/2/1 parallel edges yields 2
# distinct cycles.
# `lens` selects which cycle lengths to emit, independently: () = none,
# (3,) = triangles only, (3,4) = triangles + chordless 4-cycles, (2,3,4) = also
# parallel-line pairs.
function short_cycles(c::TxReductionCase; lens=(2, 3, 4))
    want = Set(Int.(lens))
    N, Ln = c.N, c.Ln
    cycles = Vector{Vector{Int}}()
    isempty(want) && return cycles

    pair_lines = Dict{Tuple{Int,Int},Vector{Int}}()
    for l in 1:Ln
        a, b = minmax(c.Efrom[l], c.Eto[l])
        a == b && continue                       # ignore any self-loop
        push!(get!(pair_lines, (a, b), Int[]), l)
    end
    nbr = [Int[] for _ in 1:N]
    for ((a, b), _) in pair_lines
        push!(nbr[a], b); push!(nbr[b], a)
    end
    for i in 1:N; sort!(unique!(nbr[i])); end

    # length 2: parallel lines on the same bus pair. c_a = c_b (both directions
    # of the generic cut collapse to equality).
    if 2 in want
        for (_, ls) in pair_lines
            length(ls) < 2 && continue
            for i in 1:length(ls)-1, j in i+1:length(ls)
                push!(cycles, [ls[i], ls[j]])
            end
        end
    end

    # length 3: triangles. Emit once per {a,b,w} by requiring a < b < w.
    if 3 in want
        for ((a, b), lab) in pair_lines
            for w in nbr[a]
                w <= b && continue
                haskey(pair_lines, minmax(b, w)) || continue
                for l1 in lab, l2 in pair_lines[minmax(b, w)], l3 in pair_lines[minmax(a, w)]
                    push!(cycles, [l1, l2, l3])
                end
            end
        end
    end
    4 in want || return cycles

    # length 4: u - x - w - y - u. Collect the midpoints of every 2-path, then
    # pair them up. CHORDLESS only: if u-w is itself a line the 4-cycle carries a
    # chord and its cut is dominated by the two triangle cuts, so it is skipped
    # to keep the row count down.
    mids = Dict{Tuple{Int,Int},Vector{Int}}()
    for x in 1:N
        nb = nbr[x]
        for i in 1:length(nb)-1, j in i+1:length(nb)
            push!(get!(mids, (nb[i], nb[j]), Int[]), x)
        end
    end
    for ((u, w), xs) in mids
        length(xs) < 2 && continue
        haskey(pair_lines, (u, w)) && continue       # chorded -> skip
        sort!(xs)
        for i in 1:length(xs)-1, j in i+1:length(xs)
            x, y = xs[i], xs[j]
            for l1 in pair_lines[minmax(u, x)], l2 in pair_lines[minmax(x, w)],
                l3 in pair_lines[minmax(w, y)], l4 in pair_lines[minmax(y, u)]
                push!(cycles, [l1, l2, l3, l4])
            end
        end
    end
    return cycles
end

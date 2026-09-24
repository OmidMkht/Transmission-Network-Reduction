# Test cases and demand sets, shared by the runners.
#
#   demands = :single   the case's base demand
#             :scaled   n_demands uniform load scales over scale_range; the
#                       midpoints between them are held out
#             :hourly   ACTIVSg hourly scenarios (month / horizon); n_demands
#                       hours are chosen as design demands, the rest held out.
#                       Needs common/multiscenario.jl and the scenario matrices
#                       from common/make_hourly_scenarios.jl.

module Cases

export CASES, case_file, load_demands

const CASES = Dict(
    :case6ww     => "6busww/6busww_toy.m",
    :case14toy   => "14bus/14bus_toy.m",
    :case14      => "pglib_opf_case14_ieee.m",
    :case118     => "pglib_opf_case118_ieee.m",
    :case300     => "pglib_opf_case300_ieee.m",
    :case500     => "pglib_opf_case500_goc.m",
    :case2000    => "pglib_opf_case2000_goc.m",
    :case6515    => "pglib_opf_case6515_rte.m",
    :ACTIVSg200  => "ACTIVSg200/case_ACTIVSg200.m",
    :ACTIVSg2000 => "ACTIVSg2000/case_ACTIVSg2000.m",
)

# ACTIVSg2000's hourly data use a leap year (8784 h), ACTIVSg200's do not.
const CASE_YEAR = Dict(:ACTIVSg200 => 2017, :ACTIVSg2000 => 2016,
                       :case6ww => 2015, :case14toy => 2015)

const ROOT = dirname(@__DIR__)

function case_file(case::Symbol)
    haskey(CASES, case) || error("unknown case :$case. Known: " *
                                 join(sort(string.(keys(CASES))), ", "))
    path = joinpath(ROOT, "case studies", CASES[case])
    isfile(path) || error("case file not found: $path (see README, Case data)")
    return path
end

function stack(TR, c0, cases)
    TR.MultiScenarioTxReductionCase(c0.base, collect(eachindex(cases)), c0.bus_ids,
        hcat((p.load for p in cases)...), hcat((p.generation for p in cases)...),
        hcat((p.p for p in cases)...), hcat((p.thetahat for p in cases)...),
        hcat((p.fhat for p in cases)...))
end

"""
    load_demands(TR, s) -> (design, heldout)

`TR` is the module holding common/preprocessing.jl (and multiscenario.jl for
:hourly). `s` is a runner's settings. `heldout` is nothing when there is none.
"""
function load_demands(TR, s)
    file = case_file(s.case)
    tl = s.opf_time_limit
    if s.demands in (:single, :scaled)
        c0 = TR.build_single_scenario_case(file; time_limit=tl)
        if s.linear_costs
            c0.base.c2 .= 0.0
            c0 = TR.redispatch_dc_opf_scenarios(c0, 1:1; relax_pmin=true, time_limit=tl)
        end
        s.demands === :single && return c0, nothing
        n = s.n_demands
        n >= 2 || error("demands = :scaled needs n_demands >= 2")
        lo, hi = s.scale_range
        scales = collect(range(lo, hi; length=n))
        held = [(scales[i] + scales[i+1]) / 2 for i in 1:n-1]
        build(ss) = stack(TR, c0,
            [TR.scale_operating_points(c0, x; relax_pmin=true, time_limit=tl) for x in ss])
        return build(scales), build(held)
    elseif s.demands === :hourly
        isdefined(TR, :build_multiscenario_tx_case) ||
            error("demands = :hourly needs common/multiscenario.jl")
        name = string(s.case)
        matrix_dir = joinpath(ROOT, "outputs", "$(name)_dcopf")
        isdir(matrix_dir) || error("no hourly scenarios at $matrix_dir; build them with " *
                                   "common/make_hourly_scenarios.jl")
        c = TR.build_multiscenario_tx_case(file, matrix_dir; time_limit=tl,
                                           validate_saved_ratings=false)
        # `month` takes one number or several, e.g. month=[9,10,11] for a
        # quarter. The months are concatenated in the order given and
        # horizon_start_day / horizon_days then slice that whole stretch, so a
        # 3-month window needs horizon_days >= 92.
        months = s.month isa Integer ? [Int(s.month)] : Int.(collect(s.month))
        isempty(months) && error("month must not be empty")
        idx = reduce(vcat, (TR.month_scenario_indices(c, m; year=CASE_YEAR[s.case],
                                                      require_complete=true)
                            for m in months))
        first_h = (s.horizon_start_day - 1) * 24 + 1
        last_h = min(first_h + s.horizon_days * 24 - 1, length(idx))
        first_h <= last_h || error("horizon_start_day is past the end of the window")
        c = TR.subset_multiscenario_case(c, idx[first_h:last_h])
        if s.linear_costs
            c.base.c2 .= 0.0
            c = TR.redispatch_dc_opf_scenarios(c; relax_pmin=true, time_limit=tl)
        end
        all_h = collect(axes(c.p, 2))
        n = min(s.n_demands, length(all_h))
        sel = TR.select_seed_scenarios(c, all_h; ranking=:flow, min_scenarios=n,
                                       max_scenarios=n, congestion_threshold=0.8)
        design = sel.scenario_indices
        rest = setdiff(all_h, design)
        return TR.subset_multiscenario_case(c, design),
               isempty(rest) ? nothing : TR.subset_multiscenario_case(c, rest)
    end
    error("demands must be :single, :scaled or :hourly, got :$(s.demands)")
end

end # module

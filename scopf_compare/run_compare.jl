# The reduced networks against exact SC-DCOPF methods, on the same hours.
#
#   julia --project=. --startup-file=no scopf_compare/run_compare.jl [key=value ...]
#
# Methods (every one solves the preventive SC-DCOPF with line outages):
#   btheta_full      every network copy (intact and each outage), B-theta     exact
#   btheta_lazy      intact copy, an outage copy added once violated          exact
#   ptdf_lazy        PTDF rows added once violated ("typical PTDF")           exact
#   ptdf_decomposed  decomposed PTDF of Alkhraijah et al., lazy rows          exact
#   ptdf_lazy_flows, ptdf_decomposed_flows   the same two in the paper's form: a flow
#                    variable per needed line, limit rows of one or two flows     exact
#   compact          the screened (network, line) rows only                   exact inside the design range
# and the reduced networks of each design run, in three forms, with design ratings
# (<form>:<run>) and with exact limits (<form>:<run>+x, from the run's masks_exact.csv):
#   red_btheta       every reduced network as a B-theta block (cluster angles, kept lines)
#   red_btheta_lazy  the intact reduced network, an outage's block added once violated
#   red_ptdf_lazy    each reduced network's PTDF rows, added once violated
#   red_ptdf_lazy_flows   the same, each row's flow a bounded variable (the paper's form)
# compact_flows: the compact rows over base-case flow variables.
# Each method is built once and solves the months in order. Accuracy is judged
# against the exact secure optimum and by the worst loading of the returned
# dispatch on the full network (intact and every outage).

SETTINGS = (
    case         = :ACTIVSg200,
    months       = [7, 6, 8, 9],     # whole months, solved in this order (7 = the design month)
    max_hours    = nothing,          # hours per month (nothing = all), for quick tests
    designs      = String[],         # design run folders (with parts/, and eval_m*/masks_exact.csv for exact:)
    methods      = [:btheta_full, :btheta_lazy, :ptdf_lazy, :ptdf_decomposed, :compact,
                    :reduced, :reduced_lazy, :reduced_ptdf],   # the last three: red_btheta, red_btheta_lazy, red_ptdf_lazy
    areas        = 4,                # areas of the decomposed PTDF: a number or a list, e.g. [4, 6]
    time_limit   = 3600.0,           # seconds per method and month; hours left unsolved count as such
    output_dir   = nothing,          # nothing = outputs/scopf_compare/<case>
)

using Printf, LinearAlgebra, DelimitedFiles, JuMP, Statistics
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "common", "settings.jl"))
S = Main.Settings.apply_overrides(SETTINGS)
@eval module Data
    include(joinpath(dirname(@__DIR__), "common", "preprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "postprocessing.jl"))
    include(joinpath(dirname(@__DIR__), "common", "multiscenario.jl"))
end
include(joinpath(ROOT, "common", "cases.jl"))
include(joinpath(ROOT, "adversarial", "adversarial.jl"))
include(joinpath(ROOT, "adversarial", "evaluate.jl"))
include(joinpath(@__DIR__, "baselines.jl"))
using .Adversarial, .Evaluate, .Baselines

OUT = isnothing(S.output_dir) ? joinpath(ROOT, "outputs", "scopf_compare", string(S.case)) : S.output_dir
mkpath(OUT)
Main.Settings.save_settings(joinpath(OUT, "settings.txt"), S)
log_io = open(joinpath(OUT, "summary.txt"), "w")
say(args...) = (println(args...); println(log_io, args...); flush(log_io))
function csv(name, header, rows)
    open(joinpath(OUT, name), "w") do io
        println(io, header)
        foreach(r -> println(io, join(r, ",")), rows)
    end
end

# ---- data ---------------------------------------------------------------------------
function month_hours(month)
    s = (case=S.case, demands=:hourly, n_demands=1, month, horizon_days=31,
         horizon_start_day=1, opf_time_limit=60.0, linear_costs=true)
    design, held = Main.Cases.load_demands(Main.Data, s)
    D = isnothing(held) ? design.load : hcat(design.load, held.load)
    ids = isnothing(held) ? design.scenario_ids : vcat(design.scenario_ids, held.scenario_ids)
    order = sortperm(ids)                         # chronological
    k = isnothing(S.max_hours) ? length(order) : min(S.max_hours, length(order))
    return design.base, D[:, order[1:k]], ids[order[1:k]]
end
data = Dict(m => month_hours(m) for m in S.months)
net = net_from_base(data[first(S.months)][1])
N, L = net.N, net.L
H = ptdf(net).H
outs = outages_of(net)
Q = lodf(net, H, outs)
say(@sprintf("%s: %d buses, %d lines, %d generators, %d outages; months %s",
             S.case, N, L, length(net.gen_bus), length(outs), join(S.months, ", ")))
# the exact secure optimum every method is judged against
ref = Dict(m => sc_dcopf(net, H, Q, outs, data[m][2])[1] for m in S.months)

# ---- designs --------------------------------------------------------------------------
name_of(s) = s == "base" ? 0 : parse(Int, s)
function read_masks(file)
    rows = readdlm(file, ',', Any, header=true)[1]
    designs = Dict{Int,Any}()
    for r in eachrow(rows)
        c = name_of(string(r[1]))
        haskey(designs, c) || (designs[c] = (mask=falses(L), rating=copy(net.F)))
        designs[c].mask[Int(r[2])] = Int(r[3]) == 1
        designs[c].rating[Int(r[2])] = Float64(r[4]) * net.F[Int(r[2])]
    end
    return designs
end
function read_designs(run)
    designs = Dict{Int,Any}()
    for part in readdir(joinpath(run, "parts"); join=true)
        f = joinpath(part, "masks.csv")
        isfile(f) && countlines(f) > 1 && merge!(designs, read_masks(f))
    end
    return designs
end
blocks(designs) = [(c, designs[c].mask, designs[c].rating) for c in sort(collect(keys(designs)))]
"Critical (line, network) pairs the design run screened."
function compact_pairs(run)
    pairs = Tuple{Int,Int}[]
    for f in vcat([filter(x -> startswith(basename(x), "critical_"), readdir(p; join=true))
                   for p in readdir(joinpath(run, "parts"); join=true)]...)
        countlines(f) > 1 || continue
        c = name_of(replace(basename(f), "critical_" => "", ".csv" => ""))
        foreach(l -> push!(pairs, (l, c)), unique(Int.(readdlm(f, ',', Any, header=true)[1][:, 1])))
    end
    return pairs
end
"Critical (network, line) pairs the design run screened, as compact PTDF rows."
function compact_rows(run)
    rows = []
    for f in vcat([filter(x -> startswith(basename(x), "critical_"), readdir(p; join=true))
                   for p in readdir(joinpath(run, "parts"); join=true)]...)
        c = name_of(replace(basename(f), "critical_" => "", ".csv" => ""))
        countlines(f) > 1 || continue
        Hc = post_ptdf(H, Q, c)
        for l in unique(Int.(readdlm(f, ',', Any, header=true)[1][:, 1]))
            push!(rows, (Hc[l, :], l))
        end
    end
    return rows
end

methods = Tuple{String,Function}[]
for m in S.methods
    m === :btheta_full && push!(methods, ("btheta_full", () -> network_scdcopf(net, [(c, falses(L)) for c in [0; outs]])))
    m === :btheta_lazy && push!(methods, ("btheta_lazy", () -> btheta_lazy(net, H, Q, outs)))
    m === :ptdf_lazy && push!(methods, ("ptdf_lazy", () -> ptdf_lazy(net, H, Q, outs)))
    m === :ptdf_lazy_flows && push!(methods, ("ptdf_lazy_flows", () -> ptdf_lazy(net, H, Q, outs; flows=true)))
    for a in vcat(S.areas)
        m === :ptdf_decomposed && push!(methods, ("ptdf_decomposed/$a", () -> ptdf_decomposed(net, H, Q, outs; areas=a)))
        m === :ptdf_decomposed_flows &&
            push!(methods, ("ptdf_decomposed_flows/$a", () -> ptdf_decomposed(net, H, Q, outs; areas=a, flows=true)))
    end
    if m === :compact && !isempty(S.designs)
        rows = compact_rows(first(S.designs))
        isempty(rows) || push!(methods, ("compact", () -> compact_scdcopf(net, rows)))
    end
    if m === :compact_flows && !isempty(S.designs)
        pairs = compact_pairs(first(S.designs))
        isempty(pairs) || push!(methods, ("compact_flows", () -> compact_flows(net, H, Q, outs, pairs)))
    end
    if m in (:reduced, :reduced_lazy, :reduced_ptdf, :reduced_ptdf_flows)
        form, build = m === :reduced ? ("red_btheta", b -> network_scdcopf(net, b)) :
                       m === :reduced_lazy ? ("red_btheta_lazy", b -> reduced_btheta_lazy(net, b)) :
                       m === :reduced_ptdf ? ("red_ptdf_lazy", b -> reduced_ptdf_lazy(net, b)) :
                                             ("red_ptdf_lazy_flows", b -> reduced_ptdf_lazy(net, b; flows=true))
        for run in S.designs
            label = basename(rstrip(run, ['/', '\\']))
            ds = read_designs(run)
            m === :reduced && say(@sprintf("design run %s: %d reduced networks, median %d buses", label, length(ds),
                                           round(Int, median([maximum(clusters(net, x.mask)) for x in values(ds)]))))
            push!(methods, ("$form:$label", () -> build(blocks(ds))))
            ex = filter(isfile, [joinpath(e, "masks_exact.csv") for e in readdir(run; join=true)
                                 if startswith(basename(e), "eval")])
            if !isempty(ex)
                dx = read_masks(first(ex))
                push!(methods, ("$form:$label+x", () -> build(blocks(dx))))
            end
        end
    end
end

# ---- solve ------------------------------------------------------------------------------
form_rows, hour_rows = [], []
for (label, build) in methods
    t0 = time()
    model = build()
    tb = time() - t0
    if startswith(label, "ptdf_decomposed")
        # a balanced dispatch: every generator at the same share of its range
        d1 = data[first(S.months)][2][:, 1]
        a = (sum(d1) - sum(net.pmin)) / sum(net.pmax .- net.pmin)
        gv = net.pmin .+ a .* (net.pmax .- net.pmin)
        say(@sprintf("ptdf_decomposed: %d areas; decomposed flows vs H p: worst gap %.1e",
                     maximum(model.area), decomposed_gap(model, net, H, gv, d1)))
    end
    for m in S.months
        _, D, ids = data[m]
        Z = ref[m]
        Baselines.DEADLINE[] = time() + S.time_limit
        hasproperty(model, :unfinished) && (model.unfinished[] = false)
        cost, disp, secs, wall, done = solve_hours(model, D; time_limit=S.time_limit)
        sz = model_size(model.m)
        ok = findall(s -> s <= done && isfinite(cost[s]) && isfinite(Z[s]), 1:done)
        gap = [100 * (cost[s] - Z[s]) / Z[s] for s in ok]
        load = [worst_loading(net, H, Q, outs, disp[:, s], D[:, s]) for s in ok]
        for (i, s) in enumerate(ok)
            push!(hour_rows, (label, m, ids[s], Z[s], cost[s], gap[i], load[i]))
        end
        push!(form_rows, (label, m, size(D, 2), size(D, 2) - done, sz.vars, sz.cons, sz.nnz, round(tb; digits=3),
                          round(secs; digits=3), round(wall; digits=3), done - length(ok),
                          isempty(gap) ? NaN : mean(gap), isempty(gap) ? NaN : maximum(gap),
                          count(>(1 + 1e-6), load), isempty(load) ? NaN : 100 * max(0.0, maximum(load) - 1)))
        say(@sprintf("%-26s month %2d: %7d vars %7d rows %9d nnz | build %6.2fs, solver %7.2fs, wall %7.2fs over %d of %d h | %d infeasible, cost gap mean %.2e%% max %.2e%%, %d h overloaded, worst %.3f%%",
                     label, m, sz.vars, sz.cons, sz.nnz, tb, secs, wall, done, size(D, 2), done - length(ok),
                     form_rows[end][12], form_rows[end][13], form_rows[end][14], form_rows[end][15]))
        tb = 0.0
    end
end
csv("forms.csv", "method,month,hours,unsolved_hours,vars,rows,nnz,build_s,solver_s,wall_s,infeasible_hours,gap_mean_pct,gap_max_pct,overload_hours,worst_overload_pct", form_rows)
csv("hours.csv", "method,month,hour_id,sc_dcopf_cost,cost,cost_gap_pct,worst_loading", hour_rows)
say("wrote ", OUT)
close(log_io)

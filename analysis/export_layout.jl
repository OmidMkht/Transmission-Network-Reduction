# --------------------------------------------------------------------------- #
# Write the bus coordinates the pipeline's plots use, so other plotting tools
# can draw on the same layout.
#
#   julia --project=. --startup-file=no analysis/export_layout.jl [casefile] [out.csv]
#
# Same graph and algorithm as network_layout() in common/plots.jl (stress
# layout on the simple bus graph), without loading CairoMakie.
# --------------------------------------------------------------------------- #

using Graphs, NetworkLayout

const ROOT = dirname(@__DIR__)
const CASE = length(ARGS) >= 1 ? ARGS[1] : "pglib_opf_case118_ieee.m"
const OUT = length(ARGS) >= 2 ? abspath(ARGS[2]) :
    joinpath(ROOT, "outputs", "proxy", "case118", "single_eps0.1", "conservative_0p1pct",
             "infeasibility", "bus_positions.csv")

@eval module TR
    include(joinpath(dirname(@__DIR__), "common", "preprocessing.jl"))
end

c = TR.build_single_scenario_case(joinpath(ROOT, "case studies", CASE); time_limit=60.0)
base = c.base

g = SimpleGraph(base.N)
for l in 1:base.Ln
    add_edge!(g, base.Efrom[l], base.Eto[l])
end
pos = Stress()(adjacency_matrix(g))

mkpath(dirname(OUT))
open(OUT, "w") do io
    println(io, "bus_index,bus_id,x,y")
    for i in 1:base.N
        println(io, join((i, c.bus_ids[i], pos[i][1], pos[i][2]), ","))
    end
end
println("wrote ", OUT)

# Figures for the 6-bus reduction review.
#
#   julia --project=. --startup-file=no analysis/plot_6bus_review.jl
#
# Reads outputs/6busww/review/*.csv written by analysis/review_6bus.jl.

using CairoMakie, DelimitedFiles, Printf, Statistics

const ROOT = dirname(@__DIR__)
const DATA = joinpath(ROOT, "case studies", "6busww")
const REV  = joinpath(ROOT, "outputs", "6busww", "review")
const OUT  = joinpath(ROOT, "outputs", "6busww")

rd(f) = readdlm(f, ',', Any, '\n'; header=true)
Lr = rd(joinpath(DATA, "network_lines.csv"))[1]
lfrom = Int.(Lr[:, 2]); lto = Int.(Lr[:, 3]); lfmax = Float64.(Lr[:, 5])
L = length(lfrom)

S, Sh = rd(joinpath(REV, "summary.csv"))
col(n) = S[:, findfirst(==(n), vec(Sh))]
cases = String.(col("case"))
buses = Int.(col("buses"))
oos_h = Int.(col("oos_violating_hours"))
oos_mw = Float64.(col("oos_max_ovl_mw"))
oos_pct = Float64.(col("oos_max_ovl_pct"))
dsg_h = Int.(col("design_violating_hours"))
gapd = Float64.(col("design_max_cost_gap_pct"))

pos = Dict(1 => (0.0, 1.0), 2 => (1.0, 1.0), 3 => (2.0, 1.0),
           4 => (0.0, 0.0), 5 => (1.0, 0.0), 6 => (2.0, 0.0))

function clusters_of(ls)
    p = collect(1:6)
    rt(x) = p[x] == x ? x : (p[x] = rt(p[x]))
    for l in ls; p[rt(lfrom[l])] = rt(lto[l]); end
    return [rt(b) for b in 1:6]
end

PAL = [RGBf(0.20,0.45,0.75), RGBf(0.93,0.55,0.18), RGBf(0.25,0.62,0.35),
       RGBf(0.80,0.30,0.30), RGBf(0.55,0.40,0.70), RGBf(0.45,0.45,0.48)]

function draw_net!(ax, merged, title)
    ax.title = title; ax.titlesize = 13
    hidedecorations!(ax); hidespines!(ax)
    xlims!(ax, -0.35, 2.35); ylims!(ax, -0.3, 1.3)
    rep = clusters_of(merged)
    ids = sort(unique(rep)); cidx = Dict(r => i for (i, r) in enumerate(ids))
    for l in 1:L
        (x1,y1) = pos[lfrom[l]]; (x2,y2) = pos[lto[l]]
        ismerged = l in merged
        lines!(ax, [x1,x2], [y1,y2];
               color = ismerged ? RGBf(0.85,0.25,0.20) : RGBf(0.78,0.80,0.83),
               linewidth = ismerged ? 5 : 1.5,
               linestyle = ismerged ? :solid : :solid)
    end
    for b in 1:6
        (x,y) = pos[b]
        scatter!(ax, [x],[y]; markersize=34, color=PAL[cidx[rep[b]]],
                 strokecolor=:black, strokewidth=1.2)
        text!(ax, x, y; text=string(b), align=(:center,:center),
              fontsize=14, color=:white, font=:bold)
    end
end

fig = Figure(size = (1180, 830), backgroundcolor = :white)

draw_net!(Axis(fig[1,1]), [5,9],
          "KKT optimum {2-4, 3-6}\n4 buses — 15 out-of-sample violations")
draw_net!(Axis(fig[1,2]), [4,6,7],
          "union {2-3, 2-5, 2-6}\n3 buses — ZERO violations all year")

# (c) reduction vs out-of-sample violations
ax3 = Axis(fig[2,1], xlabel="buses in the reduced network",
           ylabel="out-of-sample hours with a violation",
           title="Smaller is not automatically worse", titlesize=13)
keep = [i for i in eachindex(cases) if !occursin("copper", cases[i])]
scatter!(ax3, Float64.(buses[keep]), Float64.(oos_h[keep]);
         markersize=15, color=[oos_h[i] == 0 ? RGBf(0.25,0.62,0.35) : RGBf(0.85,0.30,0.25) for i in keep])
# Several cases land on the same (buses, violations) point; one label each
# would print them on top of one another, so group them first.
short(i) = replace(replace(replace(replace(cases[i], "single " => ""),
           "union " => ""), "KKT optimum " => "KKT "), " (no reduction)" => "")
grouped = Dict{Tuple{Int,Int},Vector{Int}}()
for i in keep
    push!(get!(grouped, (buses[i], oos_h[i]), Int[]), i)
end
for ((b, v), is) in grouped
    text!(ax3, b, v; text = " " * join(short.(is), ", "),
          align = (:left, :center), fontsize = 9.5)
end
xlims!(ax3, 2.5, 7.4)

# (d) violation magnitude, worst line overload
ax4 = Axis(fig[2,2], ylabel="worst out-of-sample overload (% of rating)",
           title="How infeasible, when it fails", titlesize=13,
           xticks=(1:length(keep), [replace(replace(cases[i], "single " => ""),
                    "union " => "") for i in keep]),
           xticklabelrotation=pi/4, xticklabelsize=9)
barplot!(ax4, 1:length(keep), [oos_pct[i] for i in keep];
         color=[oos_pct[i] == 0 ? RGBf(0.25,0.62,0.35) : RGBf(0.85,0.30,0.25) for i in keep])
ylims!(ax4, 0, maximum(oos_pct[keep]) * 1.28)
for (j,i) in enumerate(keep)
    oos_pct[i] == 0 && continue
    text!(ax4, j, oos_pct[i]; text=@sprintf(" %.2f%%\n %.4f p.u.", oos_pct[i], oos_mw[i]/100),
          align=(:center,:bottom), fontsize=8.5)
end

Label(fig[0, :], "6-bus toy — reduction quality in and out of sample (8760 h, scaling 0.018)",
      fontsize = 16, font = :bold)

save(joinpath(OUT, "review_6bus.png"), fig; px_per_unit = 2)
println("wrote ", joinpath(OUT, "review_6bus.png"))

# Network diagram of the 6-bus Wood & Wollenberg case, annotated with how often
# each line is at its limit over the full year.
#
#   julia --project=. --startup-file=no analysis/plot_6bus_network.jl [scale=0.05]
#
# Congestion frequency comes from outputs/6busww/lines_relaxed*.csv, written by
# analysis/run_6bus_year.jl at the same scale.

using CairoMakie, DelimitedFiles, Printf

const ROOT = dirname(@__DIR__)
const DATA = joinpath(ROOT, "case studies", "6busww")
const OUT  = joinpath(ROOT, "outputs", "6busww")

scale = 0.05
for a in ARGS
    k, v = split(a, Char(61); limit=2)
    k == "scale" || error("only scale=<float> is accepted")
    global scale = parse(Float64, v)
end

rd(f) = readdlm(joinpath(DATA, f), ',', Any, '\n'; header=true)[1]
Lr, G = rd("network_lines.csv"), rd("network_gens.csv")
lfrom = Int.(Lr[:, 2]); lto = Int.(Lr[:, 3])
lx = Float64.(Lr[:, 4]); lfmax = Float64.(Lr[:, 5])
gbus = Int.(G[:, 2]); gpmax = Float64.(G[:, 3]); gcost = Float64.(G[:, 5])
L = length(lfrom)

# congestion, if the matching year run exists
tag = scale == 0.08 ? "relaxed" : "relaxed_scale$(scale)"
cfile = joinpath(OUT, "lines_$(tag).csv")
pct = zeros(L)
if isfile(cfile)
    cl = readdlm(cfile, ',', Any, '\n'; header=true)[1]
    pct = Float64.(cl[:, 7])
    println("congestion from $(basename(cfile))")
else
    @warn "no $(basename(cfile)); drawing without congestion shading"
end

# Gens on the top row, loads on the bottom. 7 of the 11 lines become grid
# edges this way and only 4 are diagonals, which keeps the picture readable.
pos = Dict(1 => (0.0, 1.0), 2 => (1.0, 1.0), 3 => (2.0, 1.0),
           4 => (0.0, 0.0), 5 => (1.0, 0.0), 6 => (2.0, 0.0))

renew = Dict(4 => "70 MW wind", 5 => "50 wind + 80 solar", 6 => "35 MW solar")
gen_at = Dict(b => (gpmax[k], gcost[k]) for (k, b) in enumerate(gbus))

fig = Figure(size = (1020, 620), backgroundcolor = :white)
ax = Axis(fig[1, 1],
          title = "6-bus Wood & Wollenberg — line loading over 8760 h " *
                  "(demand scaling $(scale))",
          titlesize = 17)
hidedecorations!(ax); hidespines!(ax)
xlims!(ax, -0.5, 2.55); ylims!(ax, -0.38, 1.32)

maxpct = maximum(pct; init = 0.0)
shade(p) = maxpct <= 0 ? RGBf(0.72, 0.75, 0.78) :
           (p <= 0.05 ? RGBf(0.72, 0.75, 0.78) :
            RGBf(0.85, 0.30 - 0.22 * (p / maxpct), 0.20 - 0.14 * (p / maxpct)))

for l in 1:L
    (x1, y1) = pos[lfrom[l]]
    (x2, y2) = pos[lto[l]]
    p = pct[l]
    lw = 1.8 + 6.0 * (maxpct > 0 ? p / maxpct : 0.0)
    lines!(ax, [x1, x2], [y1, y2]; color = shade(p), linewidth = lw)
    # The four diagonals cross at two points; staggering the label along each
    # line keeps the two texts at a crossing apart.
    t = get(Dict((1, 5) => 0.27, (2, 6) => 0.27, (2, 4) => 0.73, (3, 5) => 0.73),
            (lfrom[l], lto[l]), 0.5)
    mx, my = x1 + t * (x2 - x1), y1 + t * (y2 - y1)
    dx, dy = x2 - x1, y2 - y1
    n = sqrt(dx^2 + dy^2)
    ox, oy = -dy / n * 0.08, dx / n * 0.08
    hot = p >= 0.5
    txt = hot ? @sprintf("%.0f MW\n%.0f%% at limit", lfmax[l], p) :
                @sprintf("%.0f MW", lfmax[l])
    text!(ax, mx + ox, my + oy; text = txt, align = (:center, :center),
          fontsize = 10, color = hot ? RGBf(0.6, 0.1, 0.05) : RGBf(0.35, 0.35, 0.38))
end

for (b, (x, y)) in sort(collect(pos))
    isgen = haskey(gen_at, b)
    scatter!(ax, [x], [y]; markersize = 44,
             color = isgen ? RGBf(0.20, 0.45, 0.75) : RGBf(0.96, 0.72, 0.25),
             strokecolor = :black, strokewidth = 1.4)
    text!(ax, x, y; text = string(b), align = (:center, :center),
          fontsize = 17, color = isgen ? :white : :black, font = :bold)
    if isgen
        pmax, c = gen_at[b]
        text!(ax, x, y + 0.155; text = @sprintf("G%d  %.0f MW\n\$%.2f/MWh", b, pmax, c),
              align = (:center, :bottom), fontsize = 11, color = RGBf(0.12, 0.3, 0.55))
    end
    if haskey(renew, b)
        text!(ax, x, y - 0.155; text = "load\n" * renew[b],
              align = (:center, :top), fontsize = 10.5, color = RGBf(0.45, 0.35, 0.05))
    end
end

Legend(fig[2, 1],
    [LineElement(color = shade(maxpct), linewidth = 7),
     LineElement(color = RGBf(0.72, 0.75, 0.78), linewidth = 2),
     MarkerElement(color = RGBf(0.20, 0.45, 0.75), marker = :circle, markersize = 15),
     MarkerElement(color = RGBf(0.96, 0.72, 0.25), marker = :circle, markersize = 15)],
    ["congested (thicker = more hours at limit)", "never at its limit",
     "generator bus", "load bus"],
    orientation = :horizontal, framevisible = false, labelsize = 11.5)

mkpath(OUT)
out = joinpath(OUT, "network_scale$(scale).png")
save(out, fig; px_per_unit = 2)
println("wrote ", out)

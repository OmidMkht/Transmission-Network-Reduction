# Figures for a toy-case train/test review.
#
#   julia --project=. --startup-file=no analysis/plot_toy_review.jl case=14bus
#
# One network panel per approach, buses coloured by cluster and merged lines
# highlighted, plus a bar chart of reduction against out-of-sample violations.

using CairoMakie, DelimitedFiles, Printf, Random

const ROOT = dirname(@__DIR__)
const CASEID = Dict("6busww" => "case6ww", "14bus" => "case14toy")

casename = "14bus"
for a in ARGS
    k, v = split(a, Char(61); limit=2)
    k == "case" || error("only case= is accepted")
    global casename = String(v)
end
DATA = joinpath(ROOT, "case studies", casename)
REV  = joinpath(ROOT, "outputs", casename * "_review")
caseid = CASEID[casename]

rd(f) = readdlm(f, ',', Any, '\n'; header=true)
Lr = rd(joinpath(DATA, "network_lines.csv"))[1]
Gg = rd(joinpath(DATA, "network_gens.csv"))[1]
lfrom = Int.(Lr[:, 2]); lto = Int.(Lr[:, 3]); L = length(lfrom)
gbus = Set(Int.(Gg[:, 2]))
buses = sort(unique(vcat(lfrom, lto))); N = length(buses)
bidx = Dict(b => i for (i, b) in enumerate(buses))

S, Sh = rd(joinpath(REV, "summary.csv"))
col(n) = S[:, findfirst(==(n), vec(Sh))]
approaches = String.(col("approach"))
settings = String.(col("setting"))
red_pct = Float64.(col("reduction_pct"))
tev_pct = Float64.(col("test_violating_pct"))
teo_pct = Float64.(col("test_max_overload_pct"))

# --- layout: fixed hand layout for the 6-bus, spring otherwise --------------- #
function layout()
    if casename == "6busww"
        p = Dict(1 => (0.0,1.0), 2 => (1.0,1.0), 3 => (2.0,1.0),
                 4 => (0.0,0.0), 5 => (1.0,0.0), 6 => (2.0,0.0))
        return [p[b] for b in buses]
    end
    Random.seed!(7)                       # deterministic across runs
    pos = [(cos(2pi*i/N), sin(2pi*i/N)) for i in 1:N]
    adj = [Int[] for _ in 1:N]
    for l in 1:L
        i, j = bidx[lfrom[l]], bidx[lto[l]]
        i == j && continue
        push!(adj[i], j); push!(adj[j], i)
    end
    k = 1.2 / sqrt(N)
    for it in 1:600
        temp = 0.10 * (1 - it/600) + 0.002
        fx = zeros(N); fy = zeros(N)
        for i in 1:N, j in 1:N
            i == j && continue
            dx, dy = pos[i][1]-pos[j][1], pos[i][2]-pos[j][2]
            d2 = max(dx*dx + dy*dy, 1e-4)
            fx[i] += k*k*dx/d2; fy[i] += k*k*dy/d2          # repulsion
        end
        for i in 1:N, j in adj[i]
            dx, dy = pos[i][1]-pos[j][1], pos[i][2]-pos[j][2]
            d = sqrt(max(dx*dx + dy*dy, 1e-6))
            fx[i] -= d*dx/k; fy[i] -= d*dy/k                # attraction
        end
        for i in 1:N
            m = sqrt(fx[i]^2 + fy[i]^2)
            m < 1e-9 && continue
            s = min(m, temp) / m
            pos[i] = (pos[i][1] + fx[i]*s, pos[i][2] + fy[i]*s)
        end
    end
    xs = [p[1] for p in pos]; ys = [p[2] for p in pos]
    sx, sy = maximum(xs)-minimum(xs), maximum(ys)-minimum(ys)
    return [((p[1]-minimum(xs))/max(sx,1e-9)*2, (p[2]-minimum(ys))/max(sy,1e-9)*2)
            for p in pos]
end
P = layout()

"Mask for one (approach, setting); `setting` is the run directory name."
function find_mask(approach, setting)
    d = joinpath(ROOT, "outputs", approach, caseid, setting)
    isdir(d) || return nothing
    f = joinpath(d, "internal.csv")
    if !isfile(f)                                  # proxy: <setting>/<mode>/
        for sub in sort(readdir(d; join=true))
            isdir(sub) || continue
            g = joinpath(sub, "line_status.csv")
            isfile(g) && (f = g; break)
        end
    end
    isfile(f) || return nothing
    if endswith(f, "line_status.csv")
        tab, hd = readdlm(f, ',', Any, '\n'; header=true)
        return Int.(tab[:, findfirst(==("internal"), vec(hd))]) .== 1
    end
    return vec(readdlm(f, ',', Int)) .== 1
end

PAL = [RGBf(0.20,0.45,0.75), RGBf(0.93,0.55,0.18), RGBf(0.25,0.62,0.35),
       RGBf(0.80,0.30,0.30), RGBf(0.55,0.40,0.70), RGBf(0.35,0.62,0.68),
       RGBf(0.85,0.45,0.60), RGBf(0.50,0.50,0.25), RGBf(0.45,0.45,0.48),
       RGBf(0.65,0.35,0.15), RGBf(0.30,0.35,0.60), RGBf(0.60,0.60,0.30),
       RGBf(0.20,0.60,0.55), RGBf(0.70,0.25,0.45)]

function draw!(ax, mask, title)
    ax.title = title; ax.titlesize = 12.5
    hidedecorations!(ax); hidespines!(ax)
    xlims!(ax, -0.18, 2.18); ylims!(ax, -0.18, 2.18)
    p = collect(1:N); rt(x) = p[x] == x ? x : (p[x] = rt(p[x]))
    for l in findall(mask); p[rt(bidx[lfrom[l]])] = rt(bidx[lto[l]]); end
    rep = [rt(i) for i in 1:N]
    ids = sort(unique(rep)); ci = Dict(r => i for (i, r) in enumerate(ids))
    for l in 1:L
        i, j = bidx[lfrom[l]], bidx[lto[l]]
        lines!(ax, [P[i][1], P[j][1]], [P[i][2], P[j][2]];
               color = mask[l] ? RGBf(0.85,0.25,0.20) : RGBf(0.80,0.82,0.85),
               linewidth = mask[l] ? 4.5 : 1.2)
    end
    for i in 1:N
        scatter!(ax, [P[i][1]], [P[i][2]];
                 markersize = buses[i] in gbus ? 26 : 20,
                 color = PAL[mod1(ci[rep[i]], length(PAL))],
                 marker = buses[i] in gbus ? :rect : :circle,
                 strokecolor = :black, strokewidth = 1.0)
        text!(ax, P[i][1], P[i][2]; text = string(buses[i]),
              align = (:center, :center), fontsize = 9, color = :white, font = :bold)
    end
end

# one network panel per approach, at whichever setting reduced the most
best = Tuple{String,Int}[]
for a in ("KKT", "Proxy", "Greedy")
    idx = [i for i in eachindex(approaches) if approaches[i] == a &&
           !isnothing(find_mask(lowercase(a), settings[i]))]
    isempty(idx) || push!(best, (a, idx[argmax(red_pct[idx])]))
end

fig = Figure(size = (max(420 * length(best), 900), 800), backgroundcolor = :white)
for (i, (a, r)) in enumerate(best)
    draw!(Axis(fig[1, i]), find_mask(lowercase(a), settings[r]),
          @sprintf("%s  [%s]\n%d buses (%.0f%% reduction) — %.1f%% of test hours violate",
                   a, settings[r], Int(round(N*(1-red_pct[r]/100))), red_pct[r], tev_pct[r]))
end

# the bar chart covers EVERY setting, which is where the sensitivity shows
lab = [approaches[i] == "Full network (no reduction)" ? "full" :
       approaches[i] * "\n" * settings[i] for i in eachindex(approaches)]
ax = Axis(fig[2, 1:max(length(best), 1)],
          ylabel = "percent",
          title = "$(casename): every setting — reduction vs out-of-sample failure " *
                  "(one month, 20 train / rest test)",
          titlesize = 13,
          xticks = (1:length(approaches), lab), xticklabelsize = 8.5)
w = 0.24
barplot!(ax, (1:length(approaches)) .- w, red_pct; width=w*2,
         color=RGBf(0.20,0.45,0.75), label="bus reduction %")
barplot!(ax, (1:length(approaches)) .+ w, tev_pct; width=w*2,
         color=RGBf(0.85,0.30,0.25), label="test hours violating %")
for i in eachindex(approaches)
    red_pct[i] > 0 && text!(ax, i - w, red_pct[i]; text=@sprintf("%.0f%%", red_pct[i]),
        align=(:center,:bottom), fontsize=8.5)
    tev_pct[i] > 0 && text!(ax, i + w, tev_pct[i];
        text=@sprintf("%.1f%%\nmax %.1f%%", tev_pct[i], teo_pct[i]),
        align=(:center,:bottom), fontsize=8)
end
axislegend(ax; position=:lt, framevisible=false, labelsize=10)
ylims!(ax, 0, max(maximum(red_pct), maximum(tev_pct)) * 1.35 + 1)

out = joinpath(ROOT, "outputs", "$(casename)_review", "review_$(casename).png")
save(out, fig; px_per_unit = 2)
println("wrote ", out)

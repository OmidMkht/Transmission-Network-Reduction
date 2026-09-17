"""Full vs reduced network for a proxy reduction, coloured by what the DC-OPF
replay found.

    python analysis/plot_proxy_reduction.py [setting_dir] [title]

setting_dir is a proxy/run_proxy.jl <output>/<setting> folder that already holds
infeasibility/lines.csv (analysis/proxy_infeasibility.jl) and
infeasibility/bus_positions.csv (analysis/export_layout.jl). Both panels use the
same bus coordinates; a super-node sits at its representative bus.
"""
import csv
import math
import os
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SETTING = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.join(
    ROOT, "outputs", "proxy", "case118", "single_eps0.1", "conservative_0p1pct")
TITLE = sys.argv[2] if len(sys.argv) > 2 else "case118, proxy"
INF = os.path.join(SETTING, "infeasibility")

# Palette: reference chart palette (light surface). Blue = merged; the two
# status colours carry line state, always with a label and a width step.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
BLUE = "#2a78d6"
BLUE_LINE = "#86b6ef"
BLUE_WASH = "#dbe9fb"
AT_LIMIT = "#ec835a"     # status: serious
OVERLOAD = "#d03b3b"     # status: critical

plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans"],
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


def rows(name):
    with open(os.path.join(name), newline="") as fh:
        return list(csv.DictReader(fh))


pos = {int(r["bus_id"]): (float(r["x"]), float(r["y"]))
       for r in rows(os.path.join(INF, "bus_positions.csv"))}
rep = {int(r["bus_id"]): int(r["representative_bus_id"])
       for r in rows(os.path.join(SETTING, "bus_mapping.csv"))}
lines = rows(os.path.join(INF, "lines.csv"))
for l in lines:
    l["id"] = int(l["line"])
    l["u"], l["v"] = int(l["from_bus"]), int(l["to_bus"])
    for k in ("internal", "protected", "at_limit_full", "at_limit_reduced", "violated"):
        l[k] = l[k] == "1"
    for k in ("full_loading", "reduced_loading", "true_loading", "overload_pct",
              "overload_mw", "overload_pu", "rating_mw"):
        l[k] = float(l[k])

size = defaultdict(int)
for b, r in rep.items():
    size[r] += 1
n_full, n_red = len(pos), len(size)
external = [l for l in lines if not l["internal"]]
viol = [l for l in lines if l["violated"]]
lim_full = [l for l in lines if l["at_limit_full"]]
lim_red = [l for l in lines if l["at_limit_reduced"]]
missed = [l for l in lim_full if not l["at_limit_reduced"]]


def seg(ax, p, q, **kw):
    ax.plot([p[0], q[0]], [p[1], q[1]], solid_capstyle="round", **kw)


def arc(p, q, bend, n=40):
    """Quadratic Bezier from p to q bulging `bend` x chord length sideways."""
    mx, my = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2
    dx, dy = q[0] - p[0], q[1] - p[1]
    cx, cy = mx - dy * bend, my + dx * bend
    ts = [i / n for i in range(n + 1)]
    xs = [(1 - t) ** 2 * p[0] + 2 * (1 - t) * t * cx + t ** 2 * q[0] for t in ts]
    ys = [(1 - t) ** 2 * p[1] + 2 * (1 - t) * t * cy + t ** 2 * q[1] for t in ts]
    return xs, ys, (xs[n // 2], ys[n // 2])


HALO = [pe.withStroke(linewidth=3.2, foreground=SURFACE)]


def place_labels(ax, fig, items, fontsize=8.5, obstacles=()):
    """Greedy label placement in pixel space: try 12 spots around each anchor,
    keep the one that overlaps least with labels already placed and with the
    other anchors; draw a leader when the label sits away from its line."""
    fig.canvas.draw()
    to_px = ax.transData.transform
    to_data = ax.transData.inverted().transform
    dpi_scale = fig.dpi / 72.0
    placed = []
    anchors_px = [to_px(it["anchor"]) for it in items]
    nodes_px = [to_px(o) for o in obstacles]
    for idx, it in enumerate(sorted(items, key=lambda t: -len(t["text"]))):
        ax_px = to_px(it["anchor"])
        nl = it["text"].count("\n") + 1
        w = max(len(s) for s in it["text"].split("\n")) * fontsize * 0.56 * dpi_scale + 6
        h = nl * fontsize * 1.35 * dpi_scale + 4
        best = None
        for ring in (1.0, 1.9, 2.8):
            for k in range(12):
                ang = 2 * math.pi * k / 12
                r = ring * (h * 0.9 + 10)
                cx = ax_px[0] + math.cos(ang) * (r + w * 0.35 * abs(math.cos(ang)))
                cy = ax_px[1] + math.sin(ang) * r
                box = (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
                cost = ring * 40.0
                for b in placed:
                    ox = max(0, min(box[2], b[2]) - max(box[0], b[0]))
                    oy = max(0, min(box[3], b[3]) - max(box[1], b[1]))
                    cost += ox * oy
                for a in anchors_px:
                    if box[0] - 6 < a[0] < box[2] + 6 and box[1] - 6 < a[1] < box[3] + 6:
                        cost += 400.0
                for a in nodes_px:
                    if box[0] - 4 < a[0] < box[2] + 4 and box[1] - 4 < a[1] < box[3] + 4:
                        cost += 120.0
                if best is None or cost < best[0]:
                    best = (cost, cx, cy, box)
        _, cx, cy, box = best
        placed.append(box)
        tx, ty = to_data((cx, cy))
        ax.plot([it["anchor"][0], tx], [it["anchor"][1], ty], color=MUTED,
                lw=0.7, zorder=6)
        ax.text(tx, ty, it["text"], fontsize=fontsize, color=INK, ha="center",
                va="center", zorder=9, path_effects=HALO,
                fontweight="bold" if it.get("bold") else "normal")


fig, (axL, axR) = plt.subplots(1, 2, figsize=(17, 9.2), dpi=200)
for ax in (axL, axR):
    ax.set_aspect("equal")
    ax.axis("off")
xs = [p[0] for p in pos.values()]
ys = [p[1] for p in pos.values()]
pad = 0.6
for ax in (axL, axR):
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)

# ---------------- left: full network ----------------
for l in lines:
    if l["internal"]:
        seg(axL, pos[l["u"]], pos[l["v"]], color=BLUE_WASH, lw=7, zorder=1)
for l in lines:
    p, q = pos[l["u"]], pos[l["v"]]
    if l["internal"]:
        seg(axL, p, q, color=BLUE_LINE, lw=1.0, zorder=2)
    else:
        seg(axL, p, q, color=INK2, lw=1.0, zorder=3)
for l in lim_full:
    if not l["violated"]:
        seg(axL, pos[l["u"]], pos[l["v"]], color=AT_LIMIT, lw=3.0, zorder=5)
for l in viol:
    seg(axL, pos[l["u"]], pos[l["v"]], color=OVERLOAD, lw=4.4, zorder=5)
merged_bus = [b for b in pos if size[rep[b]] > 1]
single_bus = [b for b in pos if size[rep[b]] == 1]
axL.scatter([pos[b][0] for b in merged_bus], [pos[b][1] for b in merged_bus], s=10,
            color=BLUE, edgecolors=SURFACE, linewidths=0.5, zorder=4)
axL.scatter([pos[b][0] for b in single_bus], [pos[b][1] for b in single_bus], s=13,
            color=INK2, edgecolors=SURFACE, linewidths=0.5, zorder=4)
mid = lambda l, P: ((P[l["u"]][0] + P[l["v"]][0]) / 2, (P[l["u"]][1] + P[l["v"]][1]) / 2)
place_labels(axL, fig, [
    dict(anchor=mid(l, pos),
         text=("L{} +{:.2f}%".format(l["id"], l["overload_pct"]) if l["violated"]
               else "L{}".format(l["id"])),
         bold=l["violated"]) for l in lim_full], obstacles=list(pos.values()))
axL.set_title("Full network  ·  {} buses, {} lines".format(n_full, len(lines)),
              fontsize=14, color=INK, loc="left", pad=26)
axL.text(0.0, 1.01, "{} at limit at the full DC-OPF optimum; {} overloaded when the reduced "
         "dispatch runs here".format(len(lim_full), len(viol)),
         transform=axL.transAxes, fontsize=9.5, color=INK2, va="bottom")

# ---------------- right: reduced network ----------------
rpos = {r: pos[r] for r in size}
groups = defaultdict(list)
for l in external:
    a, b = rep[l["u"]], rep[l["v"]]
    groups[(min(a, b), max(a, b))].append(l)
mids = {}
for key, ls in groups.items():
    n = len(ls)
    for k, l in enumerate(sorted(ls, key=lambda t: t["id"])):
        bend = 0.0 if n == 1 else (k - (n - 1) / 2) * 0.14
        p, q = rpos[rep[l["u"]]], rpos[rep[l["v"]]]
        if rep[l["u"]] > rep[l["v"]]:
            p, q = q, p
        cx, cy, m = arc(p, q, bend)
        mids[l["id"]] = m
        if l["violated"]:
            axR.plot(cx, cy, color=OVERLOAD, lw=4.4, zorder=5, solid_capstyle="round")
        elif l["at_limit_reduced"]:
            axR.plot(cx, cy, color=AT_LIMIT, lw=3.0, zorder=5, solid_capstyle="round")
        else:
            axR.plot(cx, cy, color=INK2, lw=1.0, zorder=3, solid_capstyle="round")
big = [r for r in size if size[r] > 1]
solo = [r for r in size if size[r] == 1]
axR.scatter([rpos[r][0] for r in big], [rpos[r][1] for r in big],
            s=[22 + 11 * size[r] for r in big], color=BLUE, alpha=0.9,
            edgecolors=SURFACE, linewidths=1.6, zorder=7)
axR.scatter([rpos[r][0] for r in solo], [rpos[r][1] for r in solo], s=16,
            color=INK2, edgecolors=SURFACE, linewidths=0.6, zorder=7)
largest = max(size, key=lambda r: size[r])
labels_r = []
for l in lim_red:
    if l["violated"]:
        txt = "L{}  model {:.0f}%, true {:.2f}%".format(
            l["id"], 100 * l["reduced_loading"], 100 * l["true_loading"])
    else:
        txt = "L{}".format(l["id"])
    labels_r.append(dict(anchor=mids[l["id"]], text=txt, bold=l["violated"]))
for l in missed:
    labels_r.append(dict(anchor=mids[l["id"]],
                         text="L{}  {:.0f}% here, at limit in full".format(
                             l["id"], 100 * l["reduced_loading"])))
place_labels(axR, fig, labels_r, obstacles=list(rpos.values()) * 3)
axR.set_title("Reduced network  ·  {} buses, {} lines".format(n_red, len(external)),
              fontsize=14, color=INK, loc="left", pad=26)
axR.text(0.0, 1.01, "{} at limit in the reduced DC-OPF; {} lines at limit in full are not "
         "at limit here".format(len(lim_red), len(missed)),
         transform=axR.transAxes, fontsize=9.5, color=INK2, va="bottom")

# ---------------- legend + summary ----------------
handles = [
    Line2D([], [], color=OVERLOAD, lw=4.4, label="Overloaded on the full network by the reduced dispatch"),
    Line2D([], [], color=AT_LIMIT, lw=3.0, label="At its thermal limit"),
    Line2D([], [], color=INK2, lw=1.2, label="Surviving (external) line"),
    Line2D([], [], color=BLUE_LINE, lw=6, alpha=0.6, label="Merged (internal) line, cluster footprint"),
    Line2D([], [], marker="o", ls="", color=BLUE, markeredgecolor=SURFACE, markersize=10,
           label="Cluster / super-node (area ~ buses merged; largest {})".format(size[largest])),
    Line2D([], [], marker="o", ls="", color=INK2, markeredgecolor=SURFACE, markersize=5,
           label="Bus kept on its own"),
]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9.5,
           labelcolor=INK2, bbox_to_anchor=(0.5, 0.045))
w = max(viol, key=lambda l: l["overload_mw"]) if viol else None
summary = ("Reduced dispatch run on the full network:  {} lines overloaded".format(len(viol)))
if w:
    summary += ("  ·  worst L{} +{:.2f}% of rating = {:.2f} MW = {:.3f} p.u.  ·  total {:.2f} MW"
                .format(w["id"], w["overload_pct"], w["overload_mw"], w["overload_pu"],
                        sum(l["overload_mw"] for l in viol)))
fig.text(0.5, 0.018, summary, ha="center", fontsize=10, color=INK)
fig.suptitle("{}:  {} → {} buses ({:.1f}% reduction)".format(
    TITLE, n_full, n_red, 100 * (n_full - n_red) / n_full),
    x=0.02, ha="left", y=0.975, fontsize=17, color=INK, fontweight="bold")
fig.subplots_adjust(left=0.02, right=0.98, top=0.855, bottom=0.12, wspace=0.04)

out = os.path.join(INF, "full_vs_reduced_dcopf.png")
fig.savefig(out, dpi=200)
print("wrote", out)

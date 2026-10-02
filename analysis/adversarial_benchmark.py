# One comparison of the reduction approaches in a campaign folder
# (hpc/submit_adversarial.sh with CAMPAIGN set).
#
#   python analysis/adversarial_benchmark.py outputs/adversarial/vacc/ACTIVSg200_julbench
#
# Reads every run's parts/*/networks.csv (the designs) and eval_m*/ (forms.csv,
# hours_eval.csv, one folder per test month). Writes designs.csv and tests.csv
# next to the runs and prints both as tables. Runs that are not finished show
# what they have so far.
import csv
import glob
import os
import statistics as st
import sys


def settings(run):
    out = {}
    # an outage part: the base part may be copied from another run
    for f in sorted(glob.glob(os.path.join(run, "parts", "*", "settings.txt")))[-1:]:
        for line in open(f, encoding="utf-8"):
            if "=" in line:
                k, v = line.split("=", 1)
                out[k.strip()] = v.strip().strip('"')
    return out


def approach(s):
    cut = s.get("cut", ":dominance").lstrip(":")
    parts = [cut + (" (" + s.get("candidates", ":all").lstrip(":") + ")" if cut == "any" else "")]
    if float(s.get("tolerance", "0") or 0) > 0:
        parts.append("tolerance %g%%" % (100 * float(s["tolerance"])))
    if s.get("limits", ":all").lstrip(":") == "critical":
        parts.append("critical-only limits")
    if s.get("hop_limit", "nothing") not in ("nothing", ""):
        parts.append("hops %s" % s["hop_limit"])
    if float(s.get("polish_time", "0") or 0) > 0:
        parts.append("polish")
    ladder = s.get("cluster_hops", "nothing")
    if ladder not in ("nothing", ""):
        parts.append("ladder " + ladder[ladder.find("["):].replace("nothing", "free"))
    return ", ".join(parts)


def designs(run):
    rows = []
    for f in glob.glob(os.path.join(run, "parts", "*", "networks.csv")):
        rows += list(csv.DictReader(open(f, encoding="utf-8")))
    live = [r for r in rows if r["redundant"] not in ("true", "True")]
    if not live:
        return None
    b = [int(r["buses"]) for r in live]
    sec = [float(r["seconds"]) for r in live]
    return dict(
        networks=len(rows), designed=len(live),
        certified=sum(r["status"] == "certified" for r in live),
        limit=sum(r["status"] == "limit" for r in live),
        from_start=sum(r["status"] == "certified" and int(r["rounds"]) == 0 for r in live),
        buses_median=st.median(b), buses_mean=round(st.mean(b), 1), buses_max=max(b),
        small=sum(x <= 4 for x in b),
        kept_mean=round(st.mean(int(r["kept_lines"]) for r in live), 1),
        rounds_mean=round(st.mean(int(r["rounds"]) for r in live), 1),
        adversaries=sum(int(r["adversaries"]) for r in live),
        witness_hours=sum(int(r["witness_hours"]) for r in live),
        design_h=round(sum(sec) / 3600, 2), slowest_min=round(max(sec) / 60, 1),
        master_h=round(sum(float(r["master_s"]) for r in live) / 3600, 2),
        check_h=round(sum(float(r["check_s"]) for r in live) / 3600, 2))


def tests(run):
    out = []
    for ev in sorted(glob.glob(os.path.join(run, "eval_m*"))):
        month = os.path.basename(ev)[len("eval_m"):]
        hf, ff = os.path.join(ev, "hours_eval.csv"), os.path.join(ev, "forms.csv")
        if not (os.path.isfile(hf) and os.path.isfile(ff)):
            continue
        forms = {(r["form"], r["hours"]): r for r in csv.DictReader(open(ff, encoding="utf-8"))}
        hours = list(csv.DictReader(open(hf, encoding="utf-8")))
        for form in ("reduced", "reduced_exact", "compact", "full"):
            for part in ("hull", "test"):
                r = [x for x in hours if x["form"] == form and x["hours"] == part]
                if not r:
                    continue
                ok = [x for x in r if x["cost"] not in ("NaN", "")]
                gap = [float(x["cost_gap_pct"]) for x in ok]
                load = [float(x["worst_loading"]) for x in ok]
                inside = [x for x in ok if float(x["hull_distance_pct"]) <= 1e-4]
                fm = forms.get((form, part), {})
                out.append(dict(
                    month="7 (design)" if part == "hull" else month, form=form, hours=len(r),
                    inside_hull=len(inside), infeasible=len(r) - len(ok),
                    gap_mean_pct=round(st.mean(gap), 5) if gap else "",
                    gap_max_pct=round(max(gap), 5) if gap else "",
                    overload_hours=sum(v > 1 + 1e-6 for v in load),
                    overload_inside=sum(float(x["worst_loading"]) > 1 + 1e-6 for x in inside),
                    worst_overload_pct=round(100 * max(0.0, max(load) - 1), 4) if load else "",
                    vars=fm.get("vars", ""), solve_s=fm.get("solve_s", "")))
    return out


def table(rows, cols):
    print("| " + " | ".join(cols) + " |")
    print("|" + "---|" * len(cols))
    for r in rows:
        print("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")


def main(root):
    runs = sorted(d for d in glob.glob(os.path.join(root, "*")) if os.path.isdir(os.path.join(d, "parts")))
    drows, trows = [], []
    for run in runs:
        name, s = os.path.basename(run), settings(run)
        d = designs(run)
        if d:
            drows.append(dict(run=name, approach=approach(s), **d))
        for t in tests(run):
            trows.append(dict(run=name, **t))
    for fname, rows in (("designs.csv", drows), ("tests.csv", trows)):
        if rows:
            with open(os.path.join(root, fname), "w", newline="", encoding="utf-8") as io:
                w = csv.DictWriter(io, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
    print("\nDesigns (non-redundant networks)\n")
    table(drows, ["run", "approach", "designed", "certified", "limit", "from_start", "buses_median",
                  "buses_mean", "buses_max", "small", "kept_mean", "rounds_mean", "adversaries",
                  "design_h", "slowest_min"])
    for form in ("reduced", "reduced_exact"):
        print("\nTests, form %s (design limits)\n" % form if form == "reduced" else "\nTests, form %s\n" % form)
        table([t for t in trows if t["form"] == form],
              ["run", "month", "hours", "inside_hull", "infeasible", "gap_mean_pct", "gap_max_pct",
               "overload_hours", "overload_inside", "worst_overload_pct", "vars", "solve_s"])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "outputs/adversarial/vacc/ACTIVSg200_julbench")

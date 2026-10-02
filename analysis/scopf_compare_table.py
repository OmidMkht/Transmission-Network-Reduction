# One table from scopf_compare runs (scopf_compare/run_compare.jl), averaging the
# times over replicate runs of the same comparison.
#
#   python analysis/scopf_compare_table.py outputs/scopf_compare/ACTIVSg200_fair_r1 outputs/scopf_compare/ACTIVSg200_fair_r2
#
# Per method, over all months: model size (largest), build, solver and wall
# seconds (mean over replicates, with the spread), per-hour wall, infeasible
# hours, worst cost gap in the first month (the design month) and in the others,
# overloaded hours and the worst overload. Writes table.csv next to the first run.
import csv
import os
import statistics as st
import sys


def read(run):
    return list(csv.DictReader(open(os.path.join(run, "forms.csv"), encoding="utf-8")))


def main(runs):
    reps = [read(r) for r in runs]
    first_month = reps[0][0]["month"]
    out = []
    for method in dict.fromkeys(r["method"] for r in reps[0]):
        per = []
        for rows in reps:
            mine = [r for r in rows if r["method"] == method]
            if mine:
                per.append(mine)
        rows = per[0]
        hours = sum(int(r["hours"]) for r in rows)
        wall = [sum(float(r["wall_s"]) for r in m) for m in per]
        solver = [sum(float(r["solver_s"]) for r in m) for m in per]
        design = [float(r["gap_max_pct"]) for r in rows if r["month"] == first_month]
        test = [float(r["gap_max_pct"]) for r in rows if r["month"] != first_month]
        out.append(dict(
            method=method, hours=hours,
            unsolved=sum(int(r["unsolved_hours"]) for r in rows),
            vars=max(int(r["vars"]) for r in rows), rows=max(int(r["rows"]) for r in rows),
            nnz=max(int(r["nnz"]) for r in rows),
            build_s=round(st.mean(sum(float(r["build_s"]) for r in m) for m in per), 2),
            solver_s=round(st.mean(solver), 2), wall_s=round(st.mean(wall), 2),
            wall_spread_s=round(max(wall) - min(wall), 2),
            ms_per_hour=round(1000 * st.mean(wall) / hours, 3),
            solver_ms_per_hour=round(1000 * st.mean(solver) / hours, 3),
            infeasible=sum(int(r["infeasible_hours"]) for r in rows),
            gap_design_pct=max(design) if design else "",
            gap_test_pct=max(test) if test else "",
            overload_hours=sum(int(r["overload_hours"]) for r in rows),
            worst_overload_pct=max(float(r["worst_overload_pct"]) for r in rows)))
    with open(os.path.join(runs[0], "table.csv"), "w", newline="", encoding="utf-8") as io:
        w = csv.DictWriter(io, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    cols = ["method", "vars", "rows", "nnz", "build_s", "solver_ms_per_hour", "ms_per_hour", "wall_spread_s",
            "infeasible", "gap_design_pct", "gap_test_pct", "overload_hours", "worst_overload_pct"]
    print("| " + " | ".join(cols) + " |")
    print("|" + "---|" * len(cols))
    for r in out:
        print("| " + " | ".join(("%.2g" % r[c]) if isinstance(r[c], float) and c.startswith("gap") else str(r[c])
                                for c in cols) + " |")


if __name__ == "__main__":
    main(sys.argv[1:])

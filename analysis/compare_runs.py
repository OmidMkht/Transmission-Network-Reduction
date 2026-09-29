# Collect finished runs and compare them the way a reduction actually has to be
# compared.
#
#   python analysis/compare_runs.py                     # everything under outputs/
#   python analysis/compare_runs.py outputs/bench       # one tree
#   python analysis/compare_runs.py --csv runs.csv      # also write a csv
#
# HOW RUNS ARE COMPARED
#
# Gap and DCOPF response are ALWAYS reported. What they decide depends on
# whether every run in the pair had hours to work with.
#
# Both runs had hours (>= LONG_SECONDS) -- convergence decides:
#
#   both converged (gap 0)   -> compare on reduction and DCOPF response.
#                               The search is settled; only the answer matters.
#   exactly one converged    -> the converged one wins outright. Still report
#                               reduction and DCOPF response for both.
#   neither converged        -> no winner. At this budget both are lower bounds.
#
# Either run is shorter than that -- compare on reduction, as before. A short
# run is not expected to converge, so its gap is not evidence against it; it is
# reported for context, not used as a veto.

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LONG_SECONDS = 3600.0        # 'enough time' = hours; below this, compare on reduction
GAP_ZERO = 1e-6


def _f(row, *keys, default=None):
    for k in keys:
        v = row.get(k)
        if v not in (None, "", "NA"):
            try:
                return float(v)
            except ValueError:
                pass
    return default


def _s(row, *keys, default=""):
    for k in keys:
        v = row.get(k)
        if v not in (None, ""):
            return v
    return default


def gap_from_log(run_dir):
    """Last MIP gap Gurobi printed, as a fraction. The proxy summary has no bound."""
    logs = sorted(run_dir.glob("gurobi*.log"))
    if not logs:
        return None
    text = logs[-1].read_text(errors="ignore")
    m = re.findall(r"gap ([0-9.]+)%", text)
    if m:
        return float(m[-1]) / 100.0
    if "Best objective" in text and "best bound" in text:
        mm = re.findall(r"Best objective ([0-9.e+-]+), best bound ([0-9.e+-]+)", text)
        if mm:
            obj, bnd = float(mm[-1][0]), float(mm[-1][1])
            if obj not in (0.0,):
                return abs(bnd - obj) / abs(obj)
    return None


def load(run_dir):
    f = run_dir / "summary.csv"
    if not f.is_file():
        return None
    rows = list(csv.DictReader(f.open()))
    if not rows:
        return None
    r = rows[-1]

    merged = _f(r, "merged_lines")
    bound = _f(r, "bound")
    gap = None
    if bound is not None and merged:
        gap = abs(bound - merged) / abs(merged)      # objective is :lines
    if gap is None:
        gap = gap_from_log(run_dir)

    status = _s(r, "status", "solve_status", default="?")
    converged = status.upper().startswith("OPTIMAL") or (gap is not None and gap <= GAP_ZERO)
    seconds = _f(r, "seconds", "audit_seconds", default=0.0) or 0.0

    # DCOPF response: the KKT audit reports pass counts, the proxy reports
    # feasibility plus how far off it is. Normalise to one shape.
    dp, dc = _f(r, "design_pass"), _f(r, "design_count")
    hp, hc = _f(r, "heldout_pass"), _f(r, "heldout_count")
    # The proxy writes design_feasible as "20/20", not a boolean. Accept both.
    def pass_count(key):
        v = _s(r, key).strip().lower()
        if not v:
            return (None, None)
        if "/" in v:
            a, _, b = v.partition("/")
            try:
                return (float(a), float(b))
            except ValueError:
                return (None, None)
        return (1.0, 1.0) if v in ("true", "yes") else (0.0, 1.0)

    if dp is None:
        dp, dc = pass_count("design_feasible")
    if hp is None:
        hp, hc = pass_count("heldout_feasible")

    return dict(
        run=run_dir.name,
        path=str(run_dir.relative_to(ROOT)),
        approach=_s(r, "approach"),
        case=_s(r, "case"),
        status=status,
        converged=converged,
        gap=gap,
        short=seconds < LONG_SECONDS,
        seconds=seconds,
        buses=_f(r, "remaining_buses", "buses"),
        reduction=_f(r, "reduction_pct"),
        merged=merged,
        bound=bound,
        design=(dp, dc),
        heldout=(hp, hc),
        worst_cost=_f(r, "worst_design_cost_pct", "design_worst_cost_change_pct"),
        worst_overload=_f(r, "design_worst_overload_pct"),
        ho_cost=_f(r, "worst_heldout_cost_pct", "heldout_worst_cost_change_pct"),
        ho_overload=_f(r, "heldout_worst_overload_pct"),
    )


def verdict(a, b):
    """Apply the rule to a pair. Returns (winner_or_None, reason)."""
    def by_reduction(note):
        if a["reduction"] is None or b["reduction"] is None:
            return None, note + "; reduction missing"
        if abs(a["reduction"] - b["reduction"]) < 1e-9:
            return None, note + "; equal reduction -> separate on DCOPF response"
        w = a if a["reduction"] > b["reduction"] else b
        return w["run"], note

    # Short runs are not expected to converge, so the gap is context, not a veto.
    if a["short"] or b["short"]:
        return by_reduction("short run(s) -> compared on reduction (gap/DCOPF reported)")

    if a["converged"] and b["converged"]:
        return by_reduction("both converged -> more reduction wins (check DCOPF response)")
    if a["converged"] != b["converged"]:
        w = a if a["converged"] else b
        return w["run"], "only this run converged -> it wins outright"
    return None, "neither converged after hours -> no winner; both are lower bounds"


def fmt(x, spec=".1f", dash="-"):
    return dash if x is None else format(x, spec)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    trees = [Path(a) if Path(a).is_absolute() else ROOT / a for a in args] or [ROOT / "outputs"]
    csv_out = None
    if "--csv" in sys.argv:
        csv_out = sys.argv[sys.argv.index("--csv") + 1]

    runs = []
    for t in trees:
        for f in sorted(t.rglob("summary.csv")):
            r = load(f.parent)
            if r:
                runs.append(r)
    if not runs:
        print("no runs with a summary.csv found")
        return

    hdr = f"{'run':<34}{'status':<13}{'gap':>8}{'buses':>7}{'red%':>7}{'design':>9}{'held':>8}{'cost%':>8}{'ovl%':>8}{'sec':>8}"
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(runs, key=lambda r: (r["approach"], r["case"], r["run"])):
        flag = "" if not r["short"] else "  SHORT"
        d = "%d/%d" % (r["design"][0], r["design"][1]) if r["design"][0] is not None else "-"
        h = "%d/%d" % (r["heldout"][0], r["heldout"][1]) if r["heldout"][0] is not None else "-"
        print(f"{r['run'][:33]:<34}{r['status'][:12]:<13}"
              f"{(fmt(None) if r['gap'] is None else format(100*r['gap'], '.1f')+'%'):>8}"
              f"{fmt(r['buses'], '.0f'):>7}{fmt(r['reduction']):>7}"
              f"{d:>9}{h:>8}{fmt(r['worst_cost'], '.3f'):>8}"
              f"{fmt(r['worst_overload'], '.3f'):>8}{fmt(r['seconds'], '.0f'):>8}{flag}")

    conv = [r for r in runs if r["converged"]]
    print(f"\n{len(conv)} of {len(runs)} run(s) converged to a 0% gap.")
    if len(runs) > 1:
        print("\npairwise verdicts (same case and approach only):")
        for i, a in enumerate(runs):
            for b in runs[i + 1:]:
                if (a["case"], a["approach"]) != (b["case"], b["approach"]):
                    continue
                w, why = verdict(a, b)
                print(f"  {a['run'][:26]:<27} vs {b['run'][:26]:<27} "
                      f"-> {w or 'no winner':<27} {why}")

    if csv_out:
        keys = ["run", "path", "approach", "case", "status", "converged", "gap",
                "short", "seconds", "buses", "reduction", "merged", "bound",
                "worst_cost", "worst_overload", "ho_cost", "ho_overload"]
        with open(csv_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            for r in runs:
                w.writerow(r)
        print("\nwrote", csv_out)


if __name__ == "__main__":
    main()

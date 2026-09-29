# One table from a campaign folder (hpc/submit_campaign.sh).
#
#   python analysis/collect_campaign.py outputs/campaign
#
# Reads <case>/<run>/summary.csv (same columns for greedy, kkt and proxy), writes
# results.csv next to the case folders and prints one block per case. Runs with
# no summary.csv are listed with the last line of their SLURM log.

import csv
import sys
from pathlib import Path

CASE_ORDER = ["case118", "ACTIVSg200", "case300", "case500"]
APPROACH_ORDER = {"greedy": 0, "kkt": 1, "proxy": 2}
MODE_ORDER = {"noladder": 0, "ladder": 1}

COLUMNS = [
    "case", "approach", "mode", "setting", "status", "remaining_buses", "buses",
    "reduction_pct", "merged_lines", "merged_bound", "gap_pct",
    "seed_merged", "seed_buses", "seed_used",
    "design_deliverable", "design_pass", "design_count",
    "worst_design_overload_pct", "worst_design_cost_pct",
    "heldout_deliverable", "heldout_pass", "heldout_count",
    "worst_heldout_overload_pct", "worst_heldout_cost_pct",
    "largest_cluster", "hop_diameter", "solve_seconds", "seconds", "run",
]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def setting_value(s):
    return num(s.split("=", 1)[1]) if "=" in s else float("nan")


def setting_label(s):
    """cap is already in %; the proxy's eps is a fraction of rating, shown in %."""
    name, _, value = s.partition("=")
    x = num(value) * (100 if name == "eps" else 1)
    return f"{name} {x:g}%" if x == x else s


def last_log_line(logs, run_name):
    hits = sorted(logs.glob(run_name + "_*.log"))
    if not hits:
        return "no SLURM log"
    lines = [l.strip() for l in hits[-1].read_text(errors="replace").splitlines() if l.strip()]
    return lines[-1][:120] if lines else "empty log"


def fmt(x, digits=2):
    v = num(x)
    if v != v:
        return "-"
    return f"{v:.{digits}f}" if abs(v) < 1e6 else "inf"


def main(root):
    root = Path(root)
    logs = root / "logs"
    rows, missing = [], []
    for case_dir in sorted(p for p in root.iterdir() if p.is_dir() and p.name != "logs"):
        for run_dir in sorted(p for p in case_dir.iterdir() if p.is_dir()):
            f = run_dir / "summary.csv"
            if not f.exists():
                missing.append((case_dir.name, run_dir.name,
                                last_log_line(logs, f"{case_dir.name}_{run_dir.name}")))
                continue
            with open(f, newline="") as fh:
                r = next(csv.DictReader(fh))
            r["run"] = f"{case_dir.name}/{run_dir.name}"
            rows.append(r)

    def key(r):
        c = r["case"]
        return (CASE_ORDER.index(c) if c in CASE_ORDER else 99, c,
                APPROACH_ORDER.get(r["approach"], 9), MODE_ORDER.get(r["mode"], 9),
                setting_value(r["setting"]))
    rows.sort(key=key)

    with open(root / "results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in COLUMNS})

    head = (f"{'approach':8} {'mode':9} {'setting':10} {'status':12} {'buses':>9} "
            f"{'red%':>6} {'merged':>6} {'gap%':>7} {'seed':>6} {'used':>5} "
            f"{'design ok':>9} {'ovl%':>7} {'cost%':>7} {'held ok':>9} {'ovl%':>7} "
            f"{'cost%':>7} {'solve h':>7}")
    case = None
    for r in rows:
        if r["case"] != case:
            case = r["case"]
            print(f"\n{case}  ({r['buses']} buses)")
            print(head)
            print("-" * len(head))
        held = (f"{r['heldout_deliverable']}/{r['heldout_count']}"
                if r["heldout_count"] not in ("", "0") else "-")
        print(f"{r['approach']:8} {r['mode']:9} {setting_label(r['setting']):10} {r['status'][:12]:12} "
              f"{r['remaining_buses'] + '/' + r['buses']:>9} {fmt(r['reduction_pct'], 1):>6} "
              f"{r['merged_lines']:>6} {fmt(r['gap_pct'], 1):>7} "
              f"{(r['seed_merged'] if r['seed_merged'] != '-1' else '-'):>6} "
              f"{(r['seed_used'] or '-')[:5]:>5} "
              f"{r['design_deliverable'] + '/' + r['design_count']:>9} "
              f"{fmt(r['worst_design_overload_pct'], 4):>7} {fmt(r['worst_design_cost_pct'], 3):>7} "
              f"{held:>9} {fmt(r['worst_heldout_overload_pct'], 4):>7} "
              f"{fmt(r['worst_heldout_cost_pct'], 3):>7} "
              f"{fmt(num(r['solve_seconds']) / 3600, 2):>7}")

    if missing:
        print(f"\nno result yet ({len(missing)}):")
        for c, run, why in missing:
            print(f"  {c}/{run}: {why}")
    print(f"\n{len(rows)} run(s) -> {root / 'results.csv'}")
    print("design ok / held ok = demands whose reduced OPF dispatch is deliverable on the "
          "full network (overload <= 1e-6 of rating); ovl% = worst overload in % of rating "
          "(best reduced optimum); "
          "cost% = worst reduced-vs-full cost; seed = merged lines in the greedy seed, "
          "used = Gurobi accepted it.")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "outputs/campaign")

# One comparison table per test case, all in one csv (outputs/campaign/tables.csv).
#
#   python3 analysis/campaign_tables.py outputs/campaign
#
# Run it where the campaign ran (VACC): wall times come from sacct, using the
# job id in each run's SLURM log name (the latest job per run). Without sacct the
# wall-time column is left empty. Column meanings are written under the tables.

import csv
import re
import subprocess
import sys
from pathlib import Path

CASES = ["case118", "ACTIVSg200", "case300", "case500"]
ORDER = [("greedy", "noladder"), ("kkt", "noladder"), ("kkt", "ladder"),
         ("proxy", "noladder"), ("proxy", "ladder")]
NAME = {("greedy", "noladder"): "Greedy",
        ("kkt", "noladder"): "KKT",
        ("kkt", "ladder"): "KKT + hop ladder (3-6-free)",
        ("proxy", "noladder"): "Proxy",
        ("proxy", "ladder"): "Proxy + hop ladder (3-6-free)"}
STATUS = {"kkt_verified": "heuristic (KKT-checked)", "TIME_LIMIT": "time limit",
          "OPTIMAL": "optimal"}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def ok(x):
    return x == x and abs(x) < 1e12


def f(x, d):
    x = num(x)
    if not ok(x):
        return ""
    s = f"{x:.{d}f}"
    return s[1:] if s.startswith("-") and float(s) == 0 else s   # no "-0.000"


def hms(seconds):
    s = num(seconds)
    if not ok(s):
        return ""
    s = int(round(s))
    return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"


def slurm_seconds(elapsed):
    """sacct Elapsed: [D-]HH:MM:SS or MM:SS."""
    days, _, rest = elapsed.partition("-") if "-" in elapsed else ("0", "", elapsed)
    parts = [int(p) for p in rest.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    return int(days) * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def wall_times(logs):
    """run name -> (job id, wall seconds, state) of its latest SLURM job."""
    latest = {}
    for p in logs.glob("*.log"):
        m = re.match(r"(.+)_(\d+)\.log$", p.name)
        if m and int(m.group(2)) > latest.get(m.group(1), 0):
            latest[m.group(1)] = int(m.group(2))
    if not latest:
        return {}
    try:
        out = subprocess.run(["sacct", "-X", "-P", "-n", "--format=JobID,Elapsed,State",
                              "-j", ",".join(str(j) for j in latest.values())],
                             capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return {}
    info = {}
    for line in out.splitlines():
        jid, elapsed, state = line.split("|")[:3]
        if jid.isdigit():
            info[int(jid)] = (slurm_seconds(elapsed), state.split()[0])
    return {run: (jid,) + info[jid] for run, jid in latest.items() if jid in info}


def setting(approach, s):
    name, _, value = s.partition("=")
    v = num(value)
    return f"cost cap {v:g}%" if name == "cap" else f"flow window {100 * v:g}% of rating"


def main(root):
    root = Path(root)
    walls = wall_times(root / "logs")
    runs = {}
    for case in CASES:
        for d in sorted((root / case).glob("*/summary.csv")):
            with open(d, newline="") as fh:
                r = next(csv.DictReader(fh))
            r["wall"] = walls.get(f"{case}_{d.parent.name}")
            steps = d.parent / "steps.csv"
            r["steps"] = list(csv.DictReader(open(steps, newline=""))) if steps.exists() else []
            runs.setdefault(case, []).append(r)

    out = []
    for i, case in enumerate(CASES):
        rows = runs.get(case, [])
        if not rows:
            continue
        n = rows[0]["buses"]
        nd, nh = int(num(rows[0]["design_count"])), int(num(rows[0]["heldout_count"]) or 0)
        held = nh > 0
        # best proven bound on merged lines at each cap: a network within cap c is
        # also within any larger cap, so every KKT (no ladder) bound at cap >= c holds
        kb = {num(r["setting"].split("=")[1]): num(r["merged_bound"])
              for r in rows if r["approach"] == "kkt" and r["mode"] == "noladder"}

        def best_bound(cap):
            b = [v for c, v in kb.items() if c >= cap - 1e-12 and ok(v)]
            return min(b) if b else float("nan")

        demand = (f"{nd} design hours + {nh} held-out hours of Mar 1-7" if held
                  else "one design loading (the case's base load); no held-out set")
        out.append([f"Table {i + 1}: {case} ({n} buses, {demand})"])
        head = ["Approach", "Setting", "Buses kept", "Reduction %", "Lines merged",
                f"Design: DC-OPF passed (of {nd})", "Design: worst overload %",
                "Design: worst cost increase %"]
        if held:
            head += [f"Held-out: DC-OPF passed (of {nh})", "Held-out: worst overload %",
                     "Held-out: worst cost increase %"]
        head += ["Solver status", "MIP gap % (own model)", "Gap to best bound %",
                 "Greedy warm start (buses)", "Warm start accepted",
                 "Wall time (h:mm:ss)", "Solve time (h:mm:ss)", "Remarks"]
        out.append(head)
        key = lambda r: (ORDER.index((r["approach"], r["mode"])), num(r["setting"].split("=")[1]))
        for r in sorted(rows, key=key):
            kind = (r["approach"], r["mode"])
            is_cap = r["setting"].startswith("cap")
            cap = num(r["setting"].split("=")[1])
            merged = num(r["merged_lines"])
            bb = best_bound(cap) if is_cap else float("nan")
            to_bound = 100 * (bb - merged) / merged if ok(bb) and merged > 0 else float("nan")
            greedy = r["approach"] == "greedy"
            # warm start: the run's greedy start, or for a ladder the last rung's
            seed_buses, seed_merged = num(r["seed_buses"]), num(r["seed_merged"])
            accepted = r["seed_used"] == "true"
            gap = f(r["gap_pct"], 1)
            remarks = []
            if r["steps"]:
                last = r["steps"][-1]
                seed_buses, seed_merged = num(last["seed_buses"]), num(last["seed_merged"])
                accepted = ok(num(last["start_objective"]))
                if last.get("feasible", "true") == "false":
                    good = [s for s in r["steps"] if s.get("feasible", "true") == "true"]
                    gap = "no solution in last rung"
                    remarks.append("free rung: Gurobi rejected the greedy start and found no solution "
                                   f"in its time limit; the network is from the hop-{good[-1]['hop_cap']} rung")
            if not greedy and accepted and merged <= seed_merged and not (r["approach"] == "proxy" and r["steps"]):
                remarks.append("no improvement over the warm start" if not r["steps"]
                               else "last rung did not improve on its warm start")
            if num(r["design_pass"]) < num(r["design_count"]):
                remarks.append("fails its own design check")
            line = [NAME[kind], setting(r["approach"], r["setting"]),
                    r["remaining_buses"], f(r["reduction_pct"], 1), r["merged_lines"],
                    r["design_pass"], f(r["worst_design_overload_pct"], 4),
                    f(r["worst_design_cost_pct"], 3)]
            if held:
                line += [r["heldout_pass"], f(r["worst_heldout_overload_pct"], 4),
                         f(r["worst_heldout_cost_pct"], 3)]
            line += [STATUS.get(r["status"], r["status"]),
                     "n/a (heuristic)" if greedy else gap,
                     f"{max(to_bound, 0.0):.1f}" if ok(to_bound) else "n/a",
                     "" if greedy or not seed_buses >= 0 else f"{seed_buses:g}",
                     "" if greedy else ("yes" if accepted else "no"),
                     hms(r["wall"][1]) if r["wall"] else "",
                     hms(r["solve_seconds"]), "; ".join(remarks)]
            out.append(line)
        out.append([])

    out += [["Notes"],
            ["DC-OPF passed", "demands where the reduced network's OPF has a dispatch that is feasible on the full network (overload <= 1e-6 of rating) and, for greedy/KKT, costs no more than the cap above the full OPF"],
            ["Worst overload %", "largest line overload on the full network of that dispatch, in % of rating (0 = feasible)"],
            ["Worst cost increase %", "largest reduced-OPF cost above the full-network OPF cost, in %"],
            ["MIP gap % (own model)", "Gurobi's final gap on the model it solved, in merged lines; for a ladder this is the last rung, whose merges from earlier rungs are fixed, so it is not a gap for the unrestricted problem"],
            ["Gap to best bound %", "(best proven upper bound on merged lines - lines merged) / lines merged; the bound is the KKT (no ladder) bound at this cap or any larger cap. n/a for proxy, which has no cost cap"],
            ["Greedy warm start (buses)", "buses in the greedy network given to Gurobi as a start (for a ladder, the last rung's start)"],
            ["Warm start accepted", "whether Gurobi accepted that start"],
            ["Wall time", "SLURM elapsed time of the whole job on VACC: Julia start-up, data, greedy warm starts, solve, validation"],
            ["Solve time", "time inside the reduction itself (greedy search, or Gurobi summed over rungs)"],
            ["Time limits", "MIP: 19 h without ladder, rungs 3/6/10 h with it (case500: 38 h; 6/12/20 h). Greedy warm start: up to 3 h, 1 h per rung (case500: 4 h; 2 h). Jobs: 16 CPUs, 24 h wall (case500 48 h)"]]
    dest = root / "tables.csv"
    with open(dest, "w", newline="") as fh:
        csv.writer(fh).writerows(out)
    print(f"{sum(len(v) for v in runs.values())} runs, "
          f"{sum(1 for v in runs.values() for r in v if r['wall'])} with wall times -> {dest}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "outputs/campaign")

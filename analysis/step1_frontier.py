# One table per case from outputs/step1/<case>.csv (analysis/step1_compare.jl).
#
#   python analysis/step1_frontier.py
#
# "usable" = certified and the reduced OPF serves every evaluated demand with its
# certified ratings. Costs are the certified network's OPF cost against the full
# OPF, mean and worst over the evaluated demands. Networks that cannot be
# certified are listed with their uncertified record instead.

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "outputs" / "step1"
rows_out = []
for case in ("ACTIVSg200", "case118", "case300", "case500"):
    f = ROOT / f"{case}.csv"
    if not f.exists():
        continue
    rows = list(csv.DictReader(open(f)))
    for variant in dict.fromkeys(r["variant"] for r in rows):
        rs = [r for r in rows if r["variant"] == variant]
        n = rs[0]["hours"]
        print(f"\n{case} ({variant} set, {rs[0]['critical_lines']} lines can overload, {n} demands evaluated)")
        print(f"  {'buses':>5}  {'network':26} {'status':12} {'mean cost%':>10} {'worst%':>8} {'max derate%':>11}  {'uncertified: ok':>15}  {'secs':>6}")
        for r in sorted(rs, key=lambda r: -int(r["buses"])):
            usable = r["certified"] == "true" and r["cert_infeasible"] == "0"
            status = "certified" if usable else ("no dispatch" if r["certified"] == "true" else "not certif.")
            name = ("new " if r["method"] == "coherent" else "") + r["network"]
            secs = float(r["cluster_seconds"]) + float(r["certify_seconds"])
            cost = f"{float(r['cert_mean_cost_pct']):10.3f} {float(r['cert_worst_cost_pct']):8.3f} {float(r['max_derating_pct']):11.1f}" \
                if usable else f"{'-':>10} {'-':>8} {'-':>11}"
            print(f"  {int(r['buses']):5d}  {name:26} {status:12} {cost}  {r['orig_ok'] + '/' + n:>15}  {secs:6.1f}")
            rows_out.append(dict(case=case, variant=variant, network=name, buses=r["buses"],
                                 status=status, mean_cost_pct=r["cert_mean_cost_pct"] if usable else "",
                                 worst_cost_pct=r["cert_worst_cost_pct"] if usable else "",
                                 max_derating_pct=r["max_derating_pct"] if usable else "",
                                 uncertified_ok=r["orig_ok"], demands=n, seconds=round(secs, 1)))
with open(ROOT / "frontier.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
    w.writeheader()
    w.writerows(rows_out)
print(f"\n-> {ROOT / 'frontier.csv'}")

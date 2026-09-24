"""Extract a planning workbook (6busww.xlsx / 14bus.xlsx) into plain CSVs.

    python analysis/extract_case_hourly.py 6busww
    python analysis/extract_case_hourly.py 14bus

Writes into "case studies/<case>/":
    network_buses.csv  network_gens.csv  network_lines.csv  network_planning.csv
    hourly_demand.csv  hourly_wind.csv   hourly_solar.csv   hourly_net_load.csv

Modelling rules, all visible here rather than buried in a solver:

  load       bus_demand[b,h] = Demand_z1[h] * ScalingFactor[b] * LoadFactor[b].
             LoadFactor sums to 1 over the load buses, so it is the bus share;
             ScalingFactor sets the system size. Verified on 14bus, whose load
             factors reproduce the IEEE case14 split exactly (21.7/259 = .0838,
             94.2/259 = .3637, ...).

  renewables NEGATIVE LOAD at their own bus: MaxCap * capacity factor. The
             series is chosen by T_index (4 = wind -> Wind(z1), 5 = solar ->
             Solar(z3)), NOT by Zone. 14bus renewable 6 sits at bus 3 with
             T_index 4 but Zone 3, so zone would mislabel it as solar.

  existing   only IntI = 1 rows are part of the network, for BOTH generators and
             lines. The IntI = 0 lines are planning candidates and must not
             appear in the full or the reduced network; they are written to
             network_planning.csv for reference only. Renewables are all
             IntI = 0 but are still used, because they enter as load.

  storage    ignored.

xlsx is a zip of XML, so this needs no third-party package.
"""
import csv
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RELNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
PKGREL = "{http://schemas.openxmlformats.org/package/2006/relationships}"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = {"6busww": "6busww.xlsx", "14bus": "14bus.xlsx"}
HOURS_PER_DAY = 24
WIND_T, SOLAR_T = 4, 5


def _col(ref):
    n = 0
    for ch in re.match(r"([A-Z]+)", ref).group(1):
        n = n * 26 + (ord(ch) - 64)
    return n


def _strings(z):
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    return ["".join(t.text or "" for t in si.iter(f"{NS}t"))
            for si in root.findall(f"{NS}si")]


def _sheets(z):
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    tgt = {r.get("Id"): r.get("Target") for r in rels.findall(f"{PKGREL}Relationship")}
    out = {}
    for sh in wb.find(f"{NS}sheets").findall(f"{NS}sheet"):
        t = tgt[sh.get(f"{RELNS}id")].lstrip("/")
        out[sh.get("name")] = t if t.startswith("xl/") else "xl/" + t
    return out


def _rows(z, path, strings):
    root = ET.fromstring(z.read(path))
    out = []
    for r in root.find(f"{NS}sheetData").findall(f"{NS}row"):
        cells = {}
        for c in r.findall(f"{NS}c"):
            t, v = c.get("t"), c.find(f"{NS}v")
            if t == "inlineStr":
                el = c.find(f"{NS}is")
                val = "".join(x.text or "" for x in el.iter(f"{NS}t")) if el is not None else None
            elif v is None or v.text is None:
                val = None
            elif t == "s":
                val = strings[int(v.text)]
            elif t in ("str", "e"):
                val = v.text
            else:
                try:
                    val = float(v.text)
                except ValueError:
                    val = v.text
            if val is not None:
                cells[_col(c.get("r"))] = val
        w = max(cells) if cells else 0
        out.append([cells.get(i) for i in range(1, w + 1)])
    return out


def _series(rows):
    vals, stamps = [], []
    for r in rows[1:]:
        if not r or r[0] is None:
            continue
        y, m, d = int(r[0]), int(r[1]), int(r[2])
        for h in range(HOURS_PER_DAY):
            vals.append(float(r[3 + h]))
            stamps.append((y, m, d, h + 1))
    return vals, stamps


def _write(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    return path


def extract(case):
    data = os.path.join(ROOT, "case studies", case)
    src = os.path.join(data, CASES[case])
    with zipfile.ZipFile(src) as z:
        strings = _strings(z)
        sheets = _sheets(z)
        grab = lambda n: [r for r in _rows(z, sheets[n], strings)[1:]
                          if r and r[0] is not None]
        nodes, gens = grab("Nodes"), grab("Generators")
        lines, rens = grab("Lines"), grab("Renewables")
        demand, stamps = _series(_rows(z, sheets["Demand_z1"], strings))
        wind, _ = _series(_rows(z, sheets["Wind(z1)"], strings))
        solar, _ = _series(_rows(z, sheets["Solar(z3)"], strings))

    H = len(demand)
    buses = [int(r[0]) for r in nodes]
    lf = {int(r[0]): float(r[3]) for r in nodes}
    sf = {int(r[0]): float(r[5]) for r in nodes}
    lfsum = sum(lf.values())

    _write(os.path.join(data, "network_buses.csv"),
           ["bus", "type", "base_kv", "load_factor", "zone", "scaling_factor"],
           [[int(r[0]), int(r[1]), r[2], r[3], int(r[4]), r[5]] for r in nodes])
    _write(os.path.join(data, "network_gens.csv"),
           ["gen", "bus", "pmax", "pmin", "cost"],
           [[int(r[0]), int(r[1]), r[3], r[4], r[16]] for r in gens if int(r[13]) == 1])
    _write(os.path.join(data, "network_lines.csv"),
           ["line", "from", "to", "x", "fmax"],
           [[int(r[0]), int(r[1]), int(r[2]), r[3], r[4]] for r in lines if int(r[6]) == 1])
    # planning candidates: reference only, never part of any network here
    _write(os.path.join(data, "network_planning.csv"),
           ["line", "from", "to", "x", "fmax", "inv_cost"],
           [[int(r[0]), int(r[1]), int(r[2]), r[3], r[4], r[7]]
            for r in lines if int(r[6]) == 0])

    cf = {WIND_T: wind, SOLAR_T: solar}
    gross = {b: [demand[h] * sf[b] * lf[b] for h in range(H)] for b in buses}
    windb = {b: [0.0] * H for b in buses}
    solarb = {b: [0.0] * H for b in buses}
    for r in rens:
        b, tech, cap = int(r[1]), int(r[2]), float(r[3])
        tech in cf or _fail(f"renewable {int(r[0])} has T_index {tech}; expected 4 or 5")
        tgt = windb if tech == WIND_T else solarb
        s = cf[tech]
        for h in range(H):
            tgt[b][h] += cap * s[h]
    net = {b: [gross[b][h] - windb[b][h] - solarb[b][h] for h in range(H)] for b in buses}

    head = ["hour", "year", "month", "day", "hour_of_day"] + [f"bus{b}" for b in buses]
    for name, tab in (("hourly_demand.csv", gross), ("hourly_wind.csv", windb),
                      ("hourly_solar.csv", solarb), ("hourly_net_load.csv", net)):
        _write(os.path.join(data, name), head,
               [[h + 1, *stamps[h], *[f"{tab[b][h]:.6f}" for b in buses]] for h in range(H)])

    tot = lambda t: [sum(t[b][h] for b in buses) for h in range(H)]
    tn = tot(net)
    pmax = sum(float(r[3]) for r in gens if int(r[13]) == 1)
    nex = sum(1 for r in lines if int(r[6]) == 1)
    npl = sum(1 for r in lines if int(r[6]) == 0)
    print(f"[{case}] {len(buses)} buses   lines: {nex} existing + {npl} planning-only (excluded)")
    print(f"        generators: {sum(1 for r in gens if int(r[13]) == 1)} existing "
          f"(Pmax {pmax:.0f} MW), {sum(1 for r in gens if int(r[13]) == 0)} candidates")
    print(f"        load factors sum to {lfsum:.6f}   scaling factor "
          f"{max(sf.values()):.4f}")
    wcap = sum(float(r[3]) for r in rens if int(r[2]) == WIND_T)
    scap = sum(float(r[3]) for r in rens if int(r[2]) == SOLAR_T)
    print(f"        renewables: {wcap:.0f} MW wind + {scap:.0f} MW solar (negative load)")
    for lab, v in (("gross demand", tot(gross)), ("NET load", tn)):
        print(f"        {lab:<13} min {min(v):8.1f}  mean {sum(v)/len(v):8.1f}  max {max(v):8.1f} MW")
    print(f"        hours net load > Pmax: {sum(1 for x in tn if x > pmax)} / {H}")


def _fail(msg):
    raise SystemExit("ERROR: " + msg)


if __name__ == "__main__":
    which = sys.argv[1:] or list(CASES)
    for c in which:
        c in CASES or _fail(f"unknown case {c}; known: {', '.join(CASES)}")
        extract(c)

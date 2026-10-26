#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Validation of the numerical parameters (issue #127): how the observables
of a Pb+Pb event depend on the lattice spacing a, the number of layers Ny
and the IR regulator m, at a fixed geometry.

With a fixed `seed` and b = 0 (or an impact parameter that is a lattice
vector for all sizes, e.g. a multiple of L/128), the geometry of an event
(nucleon positions, hot spots, T_p) does not depend on `size`, `Ny` or
`m`; only the color charges change. The color charges of each setting are
varied with `colorChargeSeed`, which draws them from their own stream, so
that the spread at a fixed geometry gives the statistical error of each
setting.

The script runs the example input without JIMWLK, with b fixed and
L = 30 fm, for three scans:

  * size (a = L/size) at the Ny and m of the input,
  * Ny at a fixed size,
  * m at a fixed size,

each with several color-charge seeds, and computes dE/deta = tau int T^tautau
d^2x at tau = 0.2 fm and at the final time (0.4 fm), dN/dy, eps_2, eps_3
and their angles. It checks that the geometry is the same in all runs
(N_part and N_coll lists, b; <Q_s> of the IP-Sat input is reported) and that different
color-charge seeds give different events, and reports the mean and its
statistical error of each observable against a, Ny and m. The physical
Q_s of the Wilson lines at fixed QsMuRatio, which changes with Ny and m,
is measured by the dipole setup (validations/dipole_correlator, #123).

Usage (from validations/numerical_parameters/):
    # quick test: sizes 64 and 128, Ny and m scans at 64, 2 seeds
    python3 numerical_parameters.py --ipglasma-path ../.. --quick
    # the setup of the paper
    python3 numerical_parameters.py --ipglasma-path ../.. --max-workers 4

--ipglasma-path is a directory with the ipglasma executable,
qs2Adj_vs_Tp_vs_Y_240.in, nucleusConfigurations/ and tables/. The results
go to <datadir>/numerical_parameters_summary.json, the figures to
<datadir>/numerical_parameters_*.pdf.
"""

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "utilities"))
from read_tmunu import COMPONENTS, read_tmunu  # noqa: E402

OBSERVABLES = ["dE/deta(0.2)", "dE/deta(final)", "dN/dy", "eps_2", "eps_3",
               "Psi_2", "Psi_3"]


def write_input(template, overrides, path):
    remaining = dict(overrides)
    lines = []
    with open(template) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) >= 2 and tokens[0] in remaining:
                lines.append(f"{tokens[0]} {remaining.pop(tokens[0])}\n")
            else:
                lines.append(line)
    if remaining:
        raise KeyError(f"keys not in the template: {sorted(remaining)}")
    with open(path, "w") as f:
        f.writelines(lines)


def run_tag(cfg):
    return (f"size{cfg['size']}_ny{cfg['ny']}_m{cfg['m']:g}"
            f"_cc{cfg['color_seed']}")


def read_run(workdir):
    """The observables and the geometry of one run."""
    out = {}
    for path in glob.glob(os.path.join(workdir, "Tmunu-t*-0.ipgt")):
        tmunu, meta = read_tmunu(path)
        tau = meta["tau_fm"]
        de = tau * tmunu[..., COMPONENTS.index("T00")].sum() \
            * meta["dx_fm"] * meta["dy_fm"]
        key = "dE/deta(0.2)" if abs(tau - 0.2) < 0.05 else "dE/deta(final)"
        out[key] = float(de)
    with open(glob.glob(os.path.join(workdir, "NpartdNdy-t*-0.dat"))[0]) as f:
        values = f.read().split()
    out["dN/dy"] = float(values[1])
    ecc = np.loadtxt(os.path.join(workdir, "eccentricities0.dat"), ndmin=2)[-1]
    out.update({"eps_2": ecc[3], "Psi_2": ecc[4], "eps_3": ecc[5],
                "Psi_3": ecc[6]})
    geometry = {"b": ecc[17]}
    for name in ("NpartList0.dat", "NcollList0.dat"):
        with open(os.path.join(workdir, name)) as f:
            geometry[name] = f.read()
    with open(os.path.join(workdir, "run.log")) as f:
        log = f.read()
    for which in ("min", "avg", "max"):
        match = re.search(rf"Average Qs\({which}\) = ([0-9.eE+-]+)", log)
        geometry[f"Qs_{which}"] = float(match.group(1)) if match else None
    return {k: float(v) for k, v in out.items()}, geometry


def run_one(cfg, settings):
    workdir = os.path.join(settings["datadir"], run_tag(cfg))
    cache = os.path.join(workdir, "result.json")
    if os.path.exists(cache):
        with open(cache) as f:
            data = json.load(f)
        return cfg, data["observables"], data["geometry"], None
    os.makedirs(workdir, exist_ok=True)
    tables = os.path.join(workdir, "tables")
    if os.path.isdir(settings["tables"]) and not os.path.exists(tables):
        os.symlink(settings["tables"], tables)
    overrides = dict(settings["overrides"], size=cfg["size"],
                     sizeOutput=cfg["size"], Ny=cfg["ny"], m=cfg["m"],
                     colorChargeSeed=cfg["color_seed"])
    write_input(settings["template"], overrides,
                os.path.join(workdir, "input"))
    env = dict(os.environ, OMP_NUM_THREADS=str(settings["threads"]))
    with open(os.path.join(workdir, "run.log"), "w") as log:
        result = subprocess.run([settings["ipglasma"], "input"], cwd=workdir,
                                stdout=log, stderr=subprocess.STDOUT, env=env)
    if result.returncode != 0:
        return cfg, None, None, f"ipglasma failed (exit {result.returncode})"
    observables, geometry = read_run(workdir)
    with open(cache, "w") as f:
        json.dump({"observables": observables, "geometry": geometry}, f)
    for path in glob.glob(os.path.join(workdir, "*.ipgt")):
        os.remove(path)
    return cfg, observables, geometry, None


def statistics(values):
    values = np.array(values, dtype=float)
    err = values.std(ddof=1) / np.sqrt(len(values)) if len(values) > 1 \
        else float("nan")
    return {"mean": float(values.mean()), "error": float(err),
            "values": values.tolist()}


def plot(summary, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    scans = [("size", "a [fm]", lambda s: summary["L"] / s),
             ("ny", r"$N_y$", float), ("m", r"$m$ [GeV]", float)]
    names = ["dE/deta(0.2)", "dE/deta(final)", "dN/dy", "eps_2", "eps_3"]
    fig, axes = plt.subplots(len(names), 3, figsize=(10, 2.2 * len(names)),
                             sharex="col")
    for col, (scan, label, x_of) in enumerate(scans):
        rows = summary["scans"][scan]
        x = np.array([x_of(r["value"]) for r in rows])
        for row, name in enumerate(names):
            ax = axes[row, col]
            y = np.array([r["observables"][name]["mean"] for r in rows])
            e = np.array([r["observables"][name]["error"] for r in rows])
            ax.errorbar(x, y, yerr=e, fmt="o-", ms=4)
            if col == 0:
                ax.set_ylabel(name)
            if row == len(names) - 1:
                ax.set_xlabel(label)
            if scan == "ny":
                ax.set_xscale("log")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Dependence on a, Ny and m at a fixed Pb+Pb geometry "
                    "(issue #127).")
    parser.add_argument("--ipglasma-path", default="../..")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--template", default=os.path.join(HERE, "..", "..",
                                                           "input"))
    parser.add_argument("--datadir", default="./numerical_parameters_output")
    parser.add_argument("--L", type=float, default=30., help="[fm]")
    parser.add_argument("--b", type=float, default=0.,
                        help="impact parameter [fm]; a multiple of L/128 "
                             "is a lattice vector for all sizes")
    parser.add_argument("--seed", type=int, default=3,
                        help="the seed of the geometry")
    parser.add_argument("--sizes", type=int, nargs="+",
                        default=[128, 256, 512, 1024])
    parser.add_argument("--ny", type=int, nargs="+",
                        default=[1, 5, 10, 25, 50, 100])
    parser.add_argument("--m", type=float, nargs="+",
                        default=[0.1, 0.2, 0.4, 0.8], help="[GeV]")
    parser.add_argument("--scan-size", type=int, default=256,
                        help="the size of the Ny and m scans")
    parser.add_argument("--color-seeds", type=int, default=4)
    parser.add_argument("--maxtime", type=float, default=0.4, help="[fm/c]")
    parser.add_argument("--threads", type=int, default=1,
                        help="OpenMP threads per run")
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--quick", action="store_true",
                        help="sizes 64 and 128, Ny and m scans at 64 with "
                             "three values each, 2 color-charge seeds")
    args = parser.parse_args()
    if args.quick:
        args.sizes, args.scan_size, args.color_seeds = [64, 128], 64, 2
        args.ny, args.m = [5, 25, 50], [0.2, 0.4, 0.8]

    with open(args.template) as f:
        template_values = dict(line.split()[:2] for line in f
                               if len(line.split()) >= 2
                               and not line.startswith("#"))
    ny0, m0 = int(template_values["Ny"]), float(template_values["m"])
    path = os.path.abspath(args.ipglasma_path)
    settings = {
        "ipglasma": os.path.join(path, args.ipglasma_cmd),
        "template": os.path.abspath(args.template),
        "datadir": os.path.abspath(args.datadir),
        "threads": args.threads,
        "overrides": {
            "L": args.L, "LOutput": args.L, "seed": args.seed,
            "maxTime": args.maxtime, "bMin": args.b, "bMax": args.b,
            "rotateReactionPlane": 0, "useJIMWLK": 0, "useFluctuatingX": 0,
            "writeWilsonLines": 0, "writeTmunu": 1, "writeTmunuBinary": 1,
            "writeHydro": 0, "outputTimes": "0.2",
            "computeEccentricities": 1, "computeGluonMultiplicity": 1,
            "nucleusQsTableFileName": os.path.join(
                path, "qs2Adj_vs_Tp_vs_Y_240.in"),
            "nuclearConfigurationsPath": os.path.join(
                path, "nucleusConfigurations"),
        },
    }
    settings["tables"] = os.path.join(path, "tables")
    os.makedirs(settings["datadir"], exist_ok=True)

    color_seeds = range(1, args.color_seeds + 1)
    points = {"size": [dict(size=s, ny=ny0, m=m0) for s in args.sizes],
              "ny": [dict(size=args.scan_size, ny=n, m=m0) for n in args.ny],
              "m": [dict(size=args.scan_size, ny=ny0, m=m) for m in args.m]}
    configs = {}
    for scan in points.values():
        for point in scan:
            for cs in color_seeds:
                cfg = dict(point, color_seed=cs)
                configs[run_tag(cfg)] = cfg
    print(f"{len(configs)} runs (b = {args.b} fm, L = {args.L} fm, seed "
          f"{args.seed}) with up to {args.max_workers} workers ...",
          flush=True)
    start, results, failures = time.time(), {}, []
    with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
        futures = [pool.submit(run_one, cfg, settings)
                   for cfg in configs.values()]
        for i, future in enumerate(as_completed(futures), 1):
            cfg, obs, geometry, error = future.result()
            if error:
                failures.append(f"{run_tag(cfg)}: {error}")
            else:
                results[run_tag(cfg)] = (cfg, obs, geometry)
            if i % max(1, len(configs) // 10) == 0 or i == len(configs):
                print(f"  {i}/{len(configs)} runs done "
                      f"({time.time() - start:.0f} s)", flush=True)
    for failure in failures:
        print("FAILED:", failure)

    # the geometry must be the same in all runs: the N_part and N_coll lists
    # and b (<Q_s> is averaged on the lattice, so it depends a little on a)
    geometries = [g for _, _, g in results.values()]
    keys = ("b", "NpartList0.dat", "NcollList0.dat")
    same_geometry = all(all(g[k] == geometries[0][k] for k in keys)
                        for g in geometries)
    # different color-charge seeds must give different events
    differ = all(
        len({results[run_tag(dict(p, color_seed=cs))][1]["dN/dy"]
             for cs in color_seeds if run_tag(dict(p, color_seed=cs))
             in results}) == len(color_seeds)
        for scan in points.values() for p in scan) if len(color_seeds) > 1 \
        else True

    summary = {"settings": vars(args), "L": args.L,
               "geometry": geometries[0] if geometries else None,
               "checks": {"same_geometry": bool(same_geometry),
                          "color_seeds_differ": bool(differ),
                          "all_runs_succeeded": not failures},
               "scans": {}}
    for scan, scan_points in points.items():
        rows = []
        for p in scan_points:
            runs = [results[run_tag(dict(p, color_seed=cs))][1]
                    for cs in color_seeds
                    if run_tag(dict(p, color_seed=cs)) in results]
            if not runs:
                continue
            rows.append({"value": p[scan],
                         "observables": {name: statistics(
                             [r[name] for r in runs]) for name in OBSERVABLES}})
        summary["scans"][scan] = rows

    print(f"\nGeometry (seed {args.seed}, b = {args.b} fm): <Q_s>(avg) = "
          f"{summary['geometry']['Qs_avg']} GeV")
    for scan, rows in summary["scans"].items():
        print(f"\n{scan} scan:")
        print(f"  {scan:>6s}" + "".join(f"{n:>22s}" for n in OBSERVABLES[:5]))
        for row in rows:
            cells = "".join(
                f"{row['observables'][n]['mean']:13.5g} "
                f"({row['observables'][n]['error']:.1g})"
                for n in OBSERVABLES[:5])
            print(f"  {row['value']:>6g}{cells}")
    passed = all(summary["checks"].values())
    for name, ok in summary["checks"].items():
        print(f"{'OK  ' if ok else 'FAIL'} {name}")
    summary["passed"] = bool(passed)
    with open(os.path.join(settings["datadir"],
                           "numerical_parameters_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    plot(summary, os.path.join(settings["datadir"],
                               "numerical_parameters.pdf"))
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{settings['datadir']}/numerical_parameters_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

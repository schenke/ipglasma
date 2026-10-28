#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
The two definitions of the gluon number (gluonSpectrumEstimator, issue
#140) in the example input: Pb+Pb with hot spots and JIMWLK from x = 0.01
to 1e-4, at b = 0 (--b), 768^2 at L = 22.5 fm, the geometry of --seed and
--configs color-charge configurations (colorChargeSeed 1, 2, ...).

The gluon multiplicity is measured at the final time, so each time of
--times and configuration is a run of its own. The runs use
gluonSpectrumEstimator 1, the projection of the Coulomb-gauge potentials
and momenta onto the free boost-invariant modes, and gluonMultiplicity0.json
then also gives dN/dy and dE/dy of 2|E|^2/omega (gluonSpectrumEstimator 0),
from the same gauge-fixed fields. For both, the script reports dN/dy and
<k_T> = (dE/dy)/(dN/dy) and their ratio, and from the log the longitudinal
part of the gauge-fixed potential, sqrt(|A_L|^2/|A|^2), a check of the gauge
fixing. 2|E|^2/omega counts the gluons of a mode only on average over the
oscillation of its fields, i.e. for k tau >> 1, so the two must agree once
k tau >> 1 at the peak of the spectrum, near Q_s. The script checks that

  1. at the times from --agree-from (0.4 fm) on, the mean over the
     configurations of dN/dy and of <k_T> of 2|E|^2/omega over the
     projection is 1 within --tolerance (3%).

Check 1 is skipped with --quick. The comparison in the dilute limit, where
the lowest modes dominate and the two differ at all times, is part of
validations/dilute_limit (--gluon-spectrum-estimator 1).

Usage (from validations/gluon_spectrum_estimator/):
    # quick test: 128^2, one configuration at 0.2 and 0.4 fm
    python3 gluon_spectrum_estimator.py --ipglasma-path ../.. --quick
    # 768^2, 2 configurations at 0.2, 0.4 and 1 fm (6 runs, about 5 min
    # each with 10 threads)
    python3 gluon_spectrum_estimator.py --ipglasma-path ../.. --threads 10
"""

import argparse
import itertools
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "numerical_parameters"))
from numerical_parameters import write_input  # noqa: E402
sys.path.insert(0, os.path.join(HERE, "..", "dilute_limit"))
from leading_order import pi_line_fraction  # noqa: E402

LONGITUDINAL_RE = re.compile(
    r"longitudinal part of the gauge-fixed potential: ([0-9.eE+-]+)")
QUANTITIES = ["dN/dy", "<k_T>"]


def run_tag(cfg):
    return f"tau{cfg['tau']:g}_cc{cfg['config']}"


def run_one(cfg, settings):
    """Run the example input to the time tau for one color-charge
    configuration; both gluon numbers (cached in <workdir>/result.json)."""
    workdir = os.path.join(settings["datadir"], run_tag(cfg))
    cache = os.path.join(workdir, "result.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return cfg, json.load(f), None
    os.makedirs(workdir, exist_ok=True)
    tables = os.path.join(workdir, "tables")
    if os.path.isdir(settings["tables"]) and not os.path.exists(tables):
        os.symlink(settings["tables"], tables)
    overrides = dict(settings["overrides"], maxTime=f"{cfg['tau']:g}",
                     colorChargeSeed=cfg["config"])
    write_input(settings["template"], overrides,
                os.path.join(workdir, "input"))
    env = dict(os.environ, OMP_NUM_THREADS=str(settings["threads"]))
    with open(os.path.join(workdir, "run.log"), "w") as log:
        returncode = subprocess.run(
            [settings["ipglasma"], "input"], cwd=workdir, stdout=log,
            stderr=subprocess.STDOUT, env=env).returncode
    path = os.path.join(workdir, "gluonMultiplicity0.json")
    if returncode != 0 or not os.path.exists(path):
        return cfg, None, (f"ipglasma failed (exit {returncode}, no "
                           f"gluonMultiplicity0.json), see {workdir}")
    with open(path) as f:
        gluons = json.load(f)
    if "dN_electric_estimator" not in gluons:
        return cfg, None, (f"{path} has no dN_electric_estimator: this "
                           "ipglasma has no gluonSpectrumEstimator 1")
    with open(os.path.join(workdir, "run.log")) as f:
        match = LONGITUDINAL_RE.search(f.read())
    out = {"tau_fm": gluons["tau_fm"], "estimator": gluons["estimator"],
           "projection": {"dN/dy": gluons["dN"],
                          "dE/dy": gluons["dE_GeV"]},
           "electric": {"dN/dy": gluons["dN_electric_estimator"],
                        "dE/dy": gluons["dE_electric_estimator_GeV"]},
           "longitudinal_part": float(match[1]) if match else None,
           # the modes on the lines k = -pi/a, which only the projection
           # counts (estimated from the spectrum)
           "pi_line_fraction": pi_line_fraction(
               gluons, settings["size"], settings["L"])}
    for which in ("projection", "electric"):
        out[which]["<k_T>"] = out[which]["dE/dy"] / out[which]["dN/dy"]
    with open(cache, "w") as f:
        json.dump(out, f)
    return cfg, out, None


def mean_and_error(values):
    values = np.asarray(values, dtype=float)
    error = values.std(ddof=1) / np.sqrt(len(values)) if len(values) > 1 \
        else np.nan
    return float(values.mean()), float(error)


def analyze(runs, settings):
    by_key = {(cfg["tau"], cfg["config"]): out for cfg, out in runs}
    result = {"times": [], "configs": settings["configs"]}
    for tau in sorted(settings["times"]):
        outs = [by_key[(tau, c)] for c in range(1, settings["configs"] + 1)
                if (tau, c) in by_key]
        if not outs:
            continue
        entry = {"tau_fm": tau, "configurations": len(outs),
                 "longitudinal_part": [o["longitudinal_part"] for o in outs],
                 "pi_line_fraction": [o.get("pi_line_fraction")
                                      for o in outs]}
        for name in QUANTITIES:
            projection = [o["projection"][name] for o in outs]
            electric = [o["electric"][name] for o in outs]
            ratio = [e / p for e, p in zip(electric, projection)]
            entry[name] = {
                "projection": projection, "electric": electric,
                "ratio": ratio,
                "projection_mean": mean_and_error(projection),
                "electric_mean": mean_and_error(electric),
                "ratio_mean": mean_and_error(ratio)}
        result["times"].append(entry)
    return result


def checks_for(result, settings):
    checks = {}
    if settings["quick"]:
        return checks
    tolerance = settings["tolerance"]
    for entry in result["times"]:
        if entry["tau_fm"] < settings["agree_from"] - 1e-9:
            continue
        for name in QUANTITIES:
            ratio = entry[name]["ratio_mean"][0]
            checks[f"{name} of 2|E|^2/omega over the projection is 1 within "
                   f"{tolerance:g} at tau = {entry['tau_fm']:g} fm"] = bool(
                       abs(ratio - 1.) < tolerance)
    return checks


def print_result(result):
    print(f"\n{result['configs']} color-charge configurations; 2|E|^2/omega "
          "against the projection onto the free Bjorken modes")
    for entry in result["times"]:
        print(f"tau = {entry['tau_fm']:g} fm:")
        for name in QUANTITIES:
            q = entry[name]
            print(f"  {name:6s} projection "
                  + " ".join(f"{v:9.4g}" for v in q["projection"])
                  + "   2|E|^2/omega "
                  + " ".join(f"{v:9.4g}" for v in q["electric"])
                  + "   ratio " + " ".join(f"{v:.4f}" for v in q["ratio"]))
        longitudinal = [v for v in entry["longitudinal_part"]
                        if v is not None]
        if longitudinal:
            print("  longitudinal part of the gauge-fixed potential: "
                  + " ".join(f"{v:.4f}" for v in longitudinal))
        lines = [v for v in entry["pi_line_fraction"] if v is not None]
        if lines:
            print("  part of dN/dy on the lines k = -pi/a (estimate): "
                  + " ".join(f"{v:.1e}" for v in lines))


def plot(result, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"projection": "#2a78d6", "electric": "#eb6834"}
    labels = {"projection": "projection onto the Bjorken modes (1)",
              "electric": r"$2|E|^2/\omega$ (0)"}
    tau = np.array([e["tau_fm"] for e in result["times"]])
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, name in zip(axes, QUANTITIES):
        for which, color in colors.items():
            mean = [e[name][f"{which}_mean"][0] for e in result["times"]]
            ax.plot(tau, mean, "-", color=color, lw=1., label=labels[which])
            for e in result["times"]:
                ax.plot([e["tau_fm"]] * len(e[name][which]), e[name][which],
                        "o", color=color, ms=3)
        ax.set_xlabel(r"$\tau$ [fm]")
        ax.set_ylabel("dN/dy" if name == "dN/dy"
                      else r"$\langle k_T \rangle$ [GeV]")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="The two definitions of the gluon number in the example "
                    "input (issue #140).")
    parser.add_argument("--ipglasma-path", default="../..")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--template", default=os.path.join(HERE, "..", "..",
                                                           "input"))
    parser.add_argument("--datadir",
                        default="./gluon_spectrum_estimator_output")
    parser.add_argument("--size", type=int, default=768)
    parser.add_argument("--L", type=float, default=22.5,
                        help="[fm]; 22.5 fm suffices at b = 0")
    parser.add_argument("--b", type=float, default=0., help="[fm]")
    parser.add_argument("--seed", type=int, default=3,
                        help="the seed of the geometry")
    parser.add_argument("--configs", type=int, default=2,
                        help="color-charge configurations (colorChargeSeed "
                             "1, 2, ...)")
    parser.add_argument("--times", type=float, nargs="+",
                        default=[0.2, 0.4, 1.0], help="[fm]")
    parser.add_argument("--agree-from", type=float, default=0.4,
                        help="the earliest time [fm] of check 1")
    parser.add_argument("--tolerance", type=float, default=0.03)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--quick", action="store_true",
                        help="128^2, one configuration at 0.2 and 0.4 fm")
    args = parser.parse_args()
    if args.quick:
        args.size, args.configs, args.times = 128, 1, [0.2, 0.4]

    path = os.path.abspath(args.ipglasma_path)
    settings = dict(vars(args))
    settings.update({
        "ipglasma": os.path.join(path, args.ipglasma_cmd),
        "template": os.path.abspath(args.template),
        "datadir": os.path.abspath(args.datadir),
        "tables": os.path.join(path, "tables"),
        "overrides": {
            "size": args.size, "L": args.L, "LOutput": args.L,
            "sizeOutput": args.size, "bMin": args.b, "bMax": args.b,
            "seed": args.seed, "rotateReactionPlane": 0,
            "jimwlkSaveSnapshots": 0, "writeWilsonLines": 0,
            "computeGluonMultiplicity": 1, "gluonSpectrumEstimator": 1,
            "nucleusQsTableFileName": os.path.join(
                path, "qs2Adj_vs_Tp_vs_Y_240.in"),
            "nuclearConfigurationsPath": os.path.join(
                path, "nucleusConfigurations")}})
    os.makedirs(settings["datadir"], exist_ok=True)

    # the long runs first
    configs = [{"tau": tau, "config": c} for tau, c in itertools.product(
        sorted(args.times, reverse=True), range(1, args.configs + 1))]
    print(f"{len(configs)} runs of size {args.size}, L = {args.L:g} fm, with "
          f"up to {args.max_workers} workers ...", flush=True)
    start, runs, failures = time.time(), [], []
    with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
        futures = [pool.submit(run_one, cfg, settings) for cfg in configs]
        for future in as_completed(futures):
            cfg, out, error = future.result()
            if error:
                failures.append(f"{run_tag(cfg)}: {error}")
            else:
                runs.append((cfg, out))
    print(f"done in {time.time() - start:.0f} s")
    for failure in failures:
        print("FAILED:", failure)

    summary = {"settings": vars(args), "checks": {}}
    if runs:
        result = analyze(runs, settings)
        summary["result"] = result
        summary["checks"] = checks_for(result, settings)
        print_result(result)
    summary["checks"]["all runs succeeded"] = not failures
    print()
    for name, ok in summary["checks"].items():
        print(f"{'OK  ' if ok else 'FAIL'} {name}")
    passed = all(summary["checks"].values())
    summary["passed"] = passed
    with open(os.path.join(settings["datadir"],
                           "gluon_spectrum_estimator_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    if runs:
        plot(summary["result"], os.path.join(
            settings["datadir"], "gluon_spectrum_estimator.pdf"))
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{settings['datadir']}/gluon_spectrum_estimator_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Step-size convergence of the JIMWLK evolution (issue #128).

The MV model with a constant color-charge density (useNucleus 0) is evolved
with JIMWLK from x0 = jimwlkInitialX, with running coupling (jimwlkAlphaS 0)
and with fixed coupling, for four step sizes jimwlkDs each. All step sizes
use the same seeds, so they start from the same initial Wilson lines, and
both nuclei of a run are evolved independently (two configurations per
run). The evolution time s (ln(x0/x)/pi^2 with running and
alpha_s ln(x0/x)/pi^2 with fixed coupling) is a multiple of the largest
step, and so are the snapshots, so all step sizes end and take their
snapshots at exactly the same x.

From the Wilson lines at x0, at the snapshots and at the final x, the
script computes the dipole correlator D(r) (with the functions of
validations/dipole_correlator) and Q_s from D(r_s) = exp(-1/2),
Q_s^2 = 2/r_s^2. Per configuration, it takes the differences to the
smallest step size (the initial Wilson lines and with them most of the
fluctuations are common to all step sizes), and averages them over the
configurations. Each step size draws its own Langevin noise, so these
differences fluctuate by 1-2% per configuration at 512^2. A first-order
scheme gives Q_s(ds) - Q_s(ds_min) = c (ds - ds_min); the script fits c
with the covariance of the differences over the configurations, tests the
fit with its chi^2, extrapolates Q_s(0) = Q_s(ds_min) - c ds_min, and
measures the order p from a fit c (ds^p - ds_min^p). The same is done for
D at the r where D = exp(-1/2) at the smallest step. It checks, at the
final x, that

  1. the initial Wilson lines are the same for all step sizes,
  2. Q_s and D converge monotonically: every successive difference has the
     sign of c or is zero within two standard deviations,
  3. Q_s and D are consistent with first order: the chi^2 of the linear
     fit has a p-value above --min-p-value,
  4. the deviation c ds / Q_s(0) of Q_s (and of D) at the recommended step
     sizes of the README (0.005 with running and 0.0005 with fixed
     coupling) is below --max-deviation; otherwise the recommendation
     should be revised.

Usage (from validations/jimwlk_step_size/):
    python3 jimwlk_step_size.py --ipglasma-path ../.. --quick
    python3 jimwlk_step_size.py --ipglasma-path ../.. --max-workers 4
"""

import argparse
import glob
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
sys.path.insert(0, os.path.join(HERE, "..", "dipole_correlator"))
from dipole_correlator import (HBARC, crossing, dipole_correlator,  # noqa
                               profiles, write_input)
from read_wilson_lines import read_wilson_lines  # noqa: E402

FILE_RE = re.compile(r"WilsonLine_x_([0-9.eE+-]+)_(\d+)$")


def x_values(mode, settings):
    """The x of the initial Wilson lines, the snapshots and the final
    Wilson lines of a coupling mode."""
    s = settings[mode]["s"]
    alphas = settings[mode]["alphas"]
    # s = alpha_s ln(x0/x) / pi^2, with alpha_s -> 1 for running coupling
    y = s * np.pi ** 2 / (alphas if alphas > 0. else 1.)
    n = settings["snapshots"] + 1
    return [settings["x0"] * np.exp(-y * k / n) for k in range(n + 1)]


def run_tag(cfg):
    return f"{cfg['mode']}_ds{cfg['ds']:g}_seed{cfg['seed']}"


def run_one(cfg, settings):
    """Run IP-Glasma with JIMWLK for one step size and seed and return the
    D(r) profiles (axis) of both nuclei at all x (cached in
    <workdir>/profiles.npz)."""
    workdir = os.path.join(settings["datadir"], run_tag(cfg))
    cache = os.path.join(workdir, "profiles.npz")
    if os.path.exists(cache):
        with np.load(cache) as data:
            return cfg, {k: data[k] for k in data.files}, None
    if settings["analyze_only"]:
        return cfg, None, "no cached result (--analyze-only)"
    os.makedirs(workdir, exist_ok=True)
    for old in glob.glob(os.path.join(workdir, "WilsonLine_*")):
        os.remove(old)
    xs = x_values(cfg["mode"], settings)
    mode = settings[cfg["mode"]]
    overrides = {
        "size": settings["size"], "L": settings["L"],
        "g2muGeV": settings["g2mu"], "Ny": settings["ny"],
        "m": settings["m"], "UVDamp": 0, "seed": cfg["seed"],
        "useNucleus": 0, "runEvolution": 0, "writeWilsonLines": 2,
        "writeWilsonLineSnapshot": 0, "useFluctuatingX": 0,
        "nucleusQsTableFileName": settings["qs_table"],
        "useJIMWLK": 1, "jimwlkInitialX": f"{xs[0]:.12g}",
        "projectileX": f"{xs[-1]:.12g}", "targetX": f"{xs[-1]:.12g}",
        "jimwlkDs": f"{cfg['ds']:.12g}", "jimwlkAlphaS": mode["alphas"],
        "jimwlkMass": settings["jimwlk_mass"], "jimwlkSaveSnapshots": 1,
        "jimwlkXSnapshotList": ",".join(f"{x:.12g}" for x in xs[1:-1]),
    }
    write_input(settings["template"], overrides,
                os.path.join(workdir, "input"))
    env = dict(os.environ, OMP_NUM_THREADS=str(settings["threads"]))
    with open(os.path.join(workdir, "run.log"), "w") as log:
        result = subprocess.run([settings["ipglasma"], "input"], cwd=workdir,
                                stdout=log, stderr=subprocess.STDOUT, env=env)
    files = glob.glob(os.path.join(workdir, "WilsonLine_x_*"))
    if result.returncode != 0 or len(files) != 2 * len(xs):
        return cfg, None, (f"ipglasma failed (exit {result.returncode}, "
                           f"{len(files)} of {2 * len(xs)} Wilson-line "
                           f"files), see {workdir}")
    n = settings["size"] // 2 + 1
    out = {"axis": np.zeros((len(xs), 2, n))}
    for path in files:
        match = FILE_RE.search(os.path.basename(path))
        ix = int(np.argmin(np.abs(np.log(np.array(xs) / float(match[1])))))
        nucleus = (int(match[2]) - 1) % 2
        V, meta = read_wilson_lines(path)
        out["axis"][ix, nucleus] = profiles(dipole_correlator(V))["axis"]
        out["a_fm"] = np.array(meta["a_fm"])
        if not settings["keep_files"]:
            os.remove(path)
    np.savez(cache, **out)
    return cfg, out, None


def fit_step_size(values, ds):
    """Fit the step-size dependence of one observable.

    values has one row per configuration and one column per step size, ds
    decreasing. The fit uses the differences to the smallest step,
    d_i = O(ds_i) - O(ds_min), which share the initial Wilson lines, with
    their covariance over the configurations. A first-order scheme gives
    d_i = c (ds_i - ds_min) (generalized least squares; chi^2 tests it,
    with the F distribution for a covariance estimated from the
    configurations); a scan of d_i = c (ds_i^p - ds_min^p) over p gives
    the order p and its interval at Delta chi^2 = 1."""
    from scipy.stats import f as f_distribution

    ds = np.asarray(ds, dtype=float)
    v = np.asarray(values, dtype=float)
    n_conf = len(v)
    d = v[:, :-1] - v[:, -1:]
    mean_d = d.mean(axis=0)
    inverse = np.linalg.inv(np.atleast_2d(np.cov(d, rowvar=False)) / n_conf)

    def gls(h):
        weight = h @ inverse @ h
        c = (h @ inverse @ mean_d) / weight
        residual = mean_d - c * h
        return c, 1. / np.sqrt(weight), float(residual @ inverse @ residual)

    c, c_error, chi2 = gls(ds[:-1] - ds[-1])
    dof = len(ds) - 2
    p_value = float("nan")
    if dof > 0 and n_conf > dof:
        p_value = float(f_distribution.sf(
            chi2 * (n_conf - dof) / (dof * (n_conf - 1)), dof,
            n_conf - dof))
    powers = np.linspace(0.05, 4., 791)
    chi2_p = np.array([gls(ds[:-1] ** p - ds[-1] ** p)[2] for p in powers])
    inside = powers[chi2_p <= chi2_p.min() + 1.]
    smallest = v[:, -1]
    zero = smallest.mean() - c * ds[-1]
    zero_error = np.hypot(smallest.std(ddof=1) / np.sqrt(n_conf),
                          ds[-1] * c_error)
    # the measured deviation from the extrapolation, O(ds_i) - O(0) =
    # d_i + c ds_min
    d_error = np.sqrt(np.diag(np.atleast_2d(np.cov(d, rowvar=False)))
                      / n_conf)
    measured = np.append(mean_d, 0.) + c * ds[-1]
    measured_error = np.hypot(np.append(d_error, 0.), ds[-1] * c_error)
    steps = v[:, :-1] - v[:, 1:]
    return {
        "zero": float(zero), "zero_error": float(zero_error),
        "slope": float(c), "slope_error": float(c_error),
        "chi2": chi2, "dof": dof, "p_value": p_value,
        "order": float(powers[np.argmin(chi2_p)]),
        "order_interval": [float(inside.min()), float(inside.max())],
        "deviation_model": (c * ds / zero).tolist(),
        "deviation_model_error": (c_error * ds / zero).tolist(),
        "deviation": (measured / zero).tolist(),
        "deviation_error": (measured_error / zero).tolist(),
        "steps": steps.mean(axis=0).tolist(),
        "steps_error": (steps.std(axis=0, ddof=1)
                        / np.sqrt(n_conf)).tolist(),
    }


def analyze_mode(mode, runs, settings):
    ds_list = sorted(settings[mode]["ds"], reverse=True)
    xs = x_values(mode, settings)
    by_key = {(cfg["ds"], cfg["seed"]): out for cfg, out in runs}
    seeds = sorted({cfg["seed"] for cfg, _ in runs})
    seeds = [s for s in seeds if all((ds, s) in by_key for ds in ds_list)]
    a = float(runs[0][1]["a_fm"])
    r = np.arange(settings["size"] // 2 + 1)
    # [configuration, ds, x, r]
    profile = np.array([[by_key[(ds, s)]["axis"][:, nucleus]
                         for ds in ds_list]
                        for s in seeds for nucleus in (0, 1)])
    result = {"ds": ds_list, "x": xs, "configurations": len(profile),
              "recommended_ds": settings[mode]["recommended"], "points": []}
    # the initial Wilson lines are the same for all step sizes
    result["initial_spread"] = float(np.abs(
        profile[:, :, 0] - profile[:, -1:, 0]).max())
    for ix, x in enumerate(xs[1:], start=1):
        qs = np.array([[np.sqrt(2.) / (crossing(r, p) * a) * HBARC
                        for p in conf[:, ix]] for conf in profile])
        # D at the r where the mean D of the smallest step is exp(-1/2)
        r_star = crossing(r, profile[:, -1, ix].mean(axis=0))
        lower = int(np.floor(r_star))
        weight = r_star - lower
        d_star = ((1. - weight) * profile[:, :, ix, lower]
                  + weight * profile[:, :, ix, lower + 1])
        result["points"].append({
            "x": x, "y": float(np.log(settings["x0"] / x)),
            "r_star": r_star, "Qs_GeV": float(qs.mean()),
            "Qs": fit_step_size(qs, ds_list),
            "D": fit_step_size(d_star, ds_list)})
    return result


def checks_for(mode, result, settings):
    checks = {}
    tag = f"{mode} coupling"
    checks[f"{tag}: same initial Wilson lines for all step sizes"] = bool(
        result["initial_spread"] < 1e-10)
    final = result["points"][-1]
    i_rec = result["ds"].index(result["recommended_ds"])
    for name in ("Qs", "D"):
        fit = final[name]
        # each successive difference has the sign of the slope, or is zero
        # within two standard deviations
        sign = np.sign(fit["slope"])
        checks[f"{tag}: {name} converges monotonically"] = bool(all(
            sign * step > -2. * error for step, error
            in zip(fit["steps"], fit["steps_error"])))
        checks[f"{tag}: {name} is consistent with first order"] = bool(
            fit["p_value"] > settings["min_p_value"])
        checks[f"{tag}: {name} deviates by less than "
               f"{settings['max_deviation']:g} at jimwlkDs "
               f"{result['recommended_ds']:g}"] = bool(
            abs(fit["deviation_model"][i_rec]) < settings["max_deviation"])
    return checks


def print_mode(mode, result):
    print(f"\n{mode} coupling ({result['configurations']} configurations, "
          f"initial Wilson lines equal to {result['initial_spread']:.1e}):")
    i_rec = result["ds"].index(result["recommended_ds"])
    for point in result["points"]:
        print(f"  y = {point['y']:.3f} (x = {point['x']:.4e}), "
              f"Q_s = {point['Qs_GeV']:.4f} GeV, r* = {point['r_star']:.2f} a")
        for name in ("Qs", "D"):
            fit = point[name]
            print(f"    {name:2s}: ds -> 0: {fit['zero']:.5f} +- "
                  f"{fit['zero_error']:.5f}; deviation at ds = "
                  + ", ".join(f"{ds:g}: {v:+.4f} +- {e:.4f}" for ds, v, e
                              in zip(result["ds"], fit["deviation"],
                                     fit["deviation_error"]))
                  + f"; first order: chi2 {fit['chi2']:.2f}/{fit['dof']} "
                  f"(p = {fit['p_value']:.2f}), c ds/O(0) = "
                  f"{fit['deviation_model'][i_rec]:+.4f} +- "
                  f"{fit['deviation_model_error'][i_rec]:.4f} at ds = "
                  f"{result['recommended_ds']:g}; order "
                  f"{fit['order']:.2f} [{fit['order_interval'][0]:.2f}, "
                  f"{fit['order_interval'][1]:.2f}]")


def plot(results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(results), 2,
                             figsize=(10, 3.8 * len(results)), squeeze=False)
    for (mode, result), (left, right) in zip(results.items(), axes):
        ds = np.array(result["ds"])
        line = np.linspace(0., ds.max(), 50)
        for i, point in enumerate(result["points"]):
            for ax, name in ((left, "Qs"), (right, "D")):
                fit = point[name]
                ax.errorbar(ds, 100 * np.array(fit["deviation"]),
                            yerr=100 * np.array(fit["deviation_error"]),
                            fmt="o", ms=3, color=f"C{i}",
                            label=f"y = {point['y']:.2f}")
                ax.plot(line, 100 * fit["slope"] * line / fit["zero"], "--",
                        lw=0.8, color=f"C{i}")
        for ax, name in ((left, r"$Q_s$"), (right, r"$D(r^*)$")):
            ax.axvline(result["recommended_ds"], color="k", lw=0.6, ls=":")
            ax.axhline(0., color="k", lw=0.6)
            ax.set_xlim(left=0.)
            ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
            ax.set_xlabel("jimwlkDs")
            ax.set_ylabel(f"{name}: deviation from ds = 0 [%]")
            ax.set_title(f"{mode} coupling", fontsize=9)
        left.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path)


def main():
    parser = argparse.ArgumentParser(
        description="Step-size convergence of the JIMWLK evolution (issue "
                    "#128).")
    parser.add_argument("--ipglasma-path", default="../..")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--template", default=os.path.join(HERE, "..", "..",
                                                           "input"))
    parser.add_argument("--datadir", default="./jimwlk_step_size_output")
    parser.add_argument("--size", type=int, default=512)
    parser.add_argument("--L", type=float, default=10., help="[fm]")
    parser.add_argument("--g2mu", type=float, default=2., help="[GeV]")
    parser.add_argument("--ny", type=int, default=50)
    parser.add_argument("--m", type=float, default=0.4,
                        help="IR regulator of the initial condition [GeV]")
    parser.add_argument("--jimwlk-mass", type=float, default=0.4,
                        help="IR regulator of the JIMWLK kernel [GeV]")
    parser.add_argument("--x0", type=float, default=0.01,
                        help="jimwlkInitialX")
    parser.add_argument("--running-ds", type=float, nargs="+",
                        default=[0.01, 0.005, 0.0025, 0.00125])
    parser.add_argument("--running-s", type=float, default=0.48,
                        help="evolution time ln(x0/x)/pi^2, a multiple of "
                             "the largest step times (snapshots + 1)")
    parser.add_argument("--fixed-ds", type=float, nargs="+",
                        default=[0.002, 0.001, 0.0005, 0.00025])
    parser.add_argument("--fixed-s", type=float, default=0.04,
                        help="evolution time alpha_s ln(x0/x)/pi^2")
    parser.add_argument("--fixed-alphas", type=float, default=0.2)
    parser.add_argument("--modes", nargs="+", default=["running", "fixed"],
                        choices=["running", "fixed"])
    parser.add_argument("--snapshots", type=int, default=3,
                        help="equally spaced snapshots before the final x")
    parser.add_argument("--seeds", type=int, default=8,
                        help="runs per step size (two configurations each)")
    parser.add_argument("--max-deviation", type=float, default=0.03,
                        help="largest relative deviation at the recommended "
                             "step sizes")
    parser.add_argument("--min-p-value", type=float, default=0.01,
                        help="smallest p-value of the first-order fit")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--quick", action="store_true",
                        help="128^2 with the same lattice spacing, four "
                             "seeds")
    parser.add_argument("--analyze-only", action="store_true",
                        help="only analyze the results in --datadir")
    parser.add_argument("--keep-files", action="store_true",
                        help="keep the Wilson-line files")
    args = parser.parse_args()
    if args.quick:
        args.L *= 128 / args.size
        args.size, args.seeds = 128, 4

    path = os.path.abspath(args.ipglasma_path)
    settings = {
        "ipglasma": os.path.join(path, args.ipglasma_cmd),
        "qs_table": os.path.join(path, "qs2Adj_vs_Tp_vs_Y_240.in"),
        "template": os.path.abspath(args.template),
        "datadir": os.path.abspath(args.datadir),
        "size": args.size, "L": args.L, "g2mu": args.g2mu, "ny": args.ny,
        "m": args.m, "jimwlk_mass": args.jimwlk_mass, "x0": args.x0,
        "snapshots": args.snapshots, "threads": args.threads,
        "analyze_only": args.analyze_only, "keep_files": args.keep_files,
        "max_deviation": args.max_deviation,
        "min_p_value": args.min_p_value,
        "running": {"ds": args.running_ds, "s": args.running_s,
                    "alphas": 0., "recommended": 0.005},
        "fixed": {"ds": args.fixed_ds, "s": args.fixed_s,
                  "alphas": args.fixed_alphas, "recommended": 0.0005},
    }
    for mode in args.modes:
        steps = settings[mode]["s"] / max(settings[mode]["ds"]) \
            / (args.snapshots + 1)
        if abs(steps - round(steps)) > 1e-9:
            parser.error(f"--{mode}-s is not a multiple of the largest step "
                         f"times {args.snapshots + 1}")
        if settings[mode]["recommended"] not in settings[mode]["ds"]:
            parser.error(f"--{mode}-ds must contain "
                         f"{settings[mode]['recommended']:g}")
    os.makedirs(settings["datadir"], exist_ok=True)

    configs = [{"mode": mode, "ds": ds, "seed": seed}
               for mode in args.modes
               for ds, seed in itertools.product(settings[mode]["ds"],
                                                 range(1, args.seeds + 1))]
    # the long runs first
    configs.sort(key=lambda c: c["ds"] * (settings["fixed"]["alphas"]
                                          if c["mode"] == "fixed" else 1.))
    print(f"{len(configs)} runs of size {args.size}, L = {args.L:g} fm, "
          f"with up to {args.max_workers} workers ...", flush=True)
    start, results, failures = time.time(), [], []
    with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
        futures = [pool.submit(run_one, cfg, settings) for cfg in configs]
        for future in as_completed(futures):
            cfg, out, error = future.result()
            if error:
                failures.append(f"{run_tag(cfg)}: {error}")
            else:
                results.append((cfg, out))
    print(f"done in {time.time() - start:.0f} s")
    for failure in failures:
        print("FAILED:", failure)

    summary = {"settings": {k: v for k, v in vars(args).items()},
               "modes": {}, "checks": {}}
    for mode in args.modes:
        runs = [(cfg, out) for cfg, out in results if cfg["mode"] == mode]
        if not runs:
            continue
        result = analyze_mode(mode, runs, settings)
        summary["modes"][mode] = result
        summary["checks"].update(checks_for(mode, result, settings))
        print_mode(mode, result)
    summary["checks"]["all runs succeeded"] = not failures
    print()
    for name, ok in summary["checks"].items():
        print(f"{'OK  ' if ok else 'FAIL'} {name}")
    passed = all(summary["checks"].values())
    summary["passed"] = passed
    with open(os.path.join(settings["datadir"],
                           "jimwlk_step_size_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    if summary["modes"]:
        plot(summary["modes"], os.path.join(settings["datadir"],
                                            "jimwlk_step_size.pdf"))
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{settings['datadir']}/jimwlk_step_size_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

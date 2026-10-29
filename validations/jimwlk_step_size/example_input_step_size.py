#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
The JIMWLK step size in the example input (issue #128): Pb+Pb with hot
spots, the running coupling and IR regulator of the input, evolved from
x0 ~ 0.01 to projectileX = targetX = 1e-4, at b = 0 (--b), for several step
sizes jimwlkDs from the same initial Wilson lines.

The evolution is a multiple of the largest step, with --snapshots
equally spaced snapshots in between (by default 93 steps of 0.005 from
x0 = 0.0098, snapshots after 31 and 62), so all step sizes end and take
their snapshots at exactly the same x. The seed fixes the geometry and the
initial Wilson lines (both nuclei, two configurations per run).

A nucleus is not translation invariant, so Q_s is measured in regions:
the sites with an initial g^2 mu above 5% of its maximum (from the
geometry file), and the thirds of these by g^2 mu (dilute, medium, dense;
the dilute third includes the edges and the gaps between hot spots, where
the coupling is strongest). In each region, D(r) is the average of
Re tr V(x) V^+(x + r)/Nc over the sites x of the region, along the lattice
axes, and Q_s = sqrt(2)/r_s with D(r_s) = e^{-1/2}. After the Yang-Mills
evolution, the script reads dE/deta at 0.2 fm and at the final time, dN/dy,
eps_2 and eps_3 (with the functions of validations/numerical_parameters).

ln Q_s and ln O (eps_n itself) are fitted linearly in the step, which gives
the deviation at each step from zero step. The fit is linear in the data,
so the fit of the averages over the configurations is the average of the
fits of the single configurations, and the spread of the latter gives the
error. All step sizes start from the same initial Wilson lines, so the
fluctuations of the initial condition cancel in this error, while those of
the Langevin noise, which differs between step sizes, remain part of it
and are not taken for a step-size effect. (Three runs give only three
configurations of the observables, too few for the covariance of
jimwlk_step_size.py.) A quadratic fit tests the first order. The script
checks that

  1. the initial Wilson lines are the same for all step sizes,
  2. Q_s of every region and dN/dy deviate by less than --max-deviation at
     the recommended step (0.005),
  3. dE/deta deviates by less than --max-deviation at --energy-ds (0.0025):
     its step error grows on finer lattices, since it sits in the hard tail
     of the evolved spectrum (0.2% at 256^2, 5.6% at 768^2 with 0.005),
  4. ln Q_s of the whole nucleus is consistent with first order (the
     quadratic term within three standard deviations of zero).

Checks 2 and 3 are skipped with --quick.

Usage (from validations/jimwlk_step_size/):
    # quick test: 128^2, L = 30 fm, three step sizes, 2 seeds
    python3 example_input_step_size.py --ipglasma-path ../.. --quick
    # 768^2 at L = 22.5 fm, 3 seeds (52 minutes with 2 workers of 8
    # threads)
    python3 example_input_step_size.py --ipglasma-path ../.. --threads 8 \\
        --max-workers 2
    # 256^2 at L = 30 fm, the lattice of the example input: 16 seeds (with
    # four, the error of dN/dy is 6%, twice the tolerance)
    python3 example_input_step_size.py --ipglasma-path ../.. --size 256 \\
        --L 30 --seeds 16 --threads 6 --max-workers 4
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
sys.path.insert(0, os.path.join(HERE, "..", "numerical_parameters"))
from numerical_parameters import (HBARC, read_g2mu, read_run,  # noqa: E402
                                  write_input)
sys.path.insert(0, os.path.join(HERE, "..", "..", "utilities"))
from read_wilson_lines import read_wilson_lines  # noqa: E402

FILE_RE = re.compile(r"WilsonLine_x_([0-9.eE+-]+)_(\d+)$")
REGIONS = ["all", "dilute", "medium", "dense"]
OBSERVABLES = ["dN/dy", "dE/deta(0.2)", "dE/deta(final)", "eps_2", "eps_3"]
QS_LEVEL = np.exp(-0.5)


def evolution(settings):
    """x0, the number of largest steps, and the x of the initial Wilson
    lines, the snapshots and the final ones."""
    dy = np.pi ** 2 * max(settings["ds"])
    parts = settings["snapshots"] + 1
    n = parts * max(1, round(np.log(settings["x0"] / settings["x"])
                             / dy / parts))
    # a little more than n steps, so that rounding cannot drop one
    x0 = settings["x"] * np.exp((n + 1e-6) * dy)
    xs = [x0 * np.exp(-k * (n // parts) * dy) for k in range(parts)]
    return x0, n, xs + [settings["x"]]


def regional_dipole(V, masks, rmax):
    """D(r) along the axes (mean of +-x, +-y), r = 0..rmax, averaged over
    the sites x of each mask."""
    out = np.ones((len(masks), rmax + 1))
    for r in range(1, rmax + 1):
        d = 0.
        for axis in (0, 1):
            # Re tr V(x) V^+(x + r)/Nc, and at x - r the same field shifted
            plus = np.einsum("xyij,xyij->xy", V,
                             np.roll(V, -r, axis=axis).conj()).real / 3.
            d = d + plus + np.roll(plus, r, axis=axis)
        out[:, r] = [d[m].mean() / 4. for m in masks]
    return out


def crossing(d):
    below = np.nonzero(d < QS_LEVEL)[0]
    if below.size == 0 or below[0] == 0:
        return np.nan
    i = below[0]
    return i - 1 + (d[i - 1] - QS_LEVEL) / (d[i - 1] - d[i])


def run_tag(cfg):
    return f"ds{cfg['ds']:g}_seed{cfg['seed']}"


def run_one(cfg, settings):
    """Run the example input for one step size and seed; Q_s of the regions
    at every snapshot and the observables (cached in <workdir>/result.json)."""
    workdir = os.path.join(settings["datadir"], run_tag(cfg))
    cache = os.path.join(workdir, "result.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return cfg, json.load(f), None
    if settings["analyze_only"]:
        return cfg, None, "no cached result (--analyze-only)"
    os.makedirs(workdir, exist_ok=True)
    for old in glob.glob(os.path.join(workdir, "WilsonLine*")):
        os.remove(old)
    tables = os.path.join(workdir, "tables")
    if os.path.isdir(settings["tables"]) and not os.path.exists(tables):
        os.symlink(settings["tables"], tables)
    x0, _, xs = evolution(settings)
    overrides = dict(settings["overrides"], seed=cfg["seed"],
                     jimwlkInitialX=f"{x0:.15g}",
                     jimwlkDs=f"{cfg['ds']:.12g}",
                     jimwlkXSnapshotList=",".join(f"{x:.15g}"
                                                  for x in xs[1:-1]))
    write_input(settings["template"], overrides,
                os.path.join(workdir, "input"))
    env = dict(os.environ, OMP_NUM_THREADS=str(settings["threads"]))
    with open(os.path.join(workdir, "run.log"), "w") as log:
        returncode = subprocess.run(
            [settings["ipglasma"], "input"], cwd=workdir, stdout=log,
            stderr=subprocess.STDOUT, env=env).returncode
    files = glob.glob(os.path.join(workdir, "WilsonLine_x_*"))
    if returncode != 0 or len(files) != 2 * len(xs):
        return cfg, None, (f"ipglasma failed (exit {returncode}, "
                           f"{len(files)} of {2 * len(xs)} Wilson-line "
                           f"files), see {workdir}")
    size = settings["size"]
    rmax = min(size // 4, int(np.ceil(1.5 / (settings["L"] / size))))
    qs = np.full((2, len(REGIONS), len(xs)), np.nan)
    D = np.zeros((2, len(REGIONS), len(xs), rmax + 1))
    sites = np.zeros((2, len(REGIONS)), dtype=int)
    for path in sorted(files):
        match = FILE_RE.search(os.path.basename(path))
        ix = int(np.argmin(np.abs(np.log(np.array(xs) / float(match[1])))))
        nucleus = (int(match[2]) - 1) % 2
        g2mu = read_g2mu(os.path.join(
            workdir, f"WilsonLineGeometry_{match[2]}")).reshape(size, size)
        inside = g2mu > 0.05 * g2mu.max()
        low, high = np.quantile(g2mu[inside], [1. / 3., 2. / 3.])
        masks = [inside, inside & (g2mu < low),
                 inside & (g2mu >= low) & (g2mu < high), g2mu >= high]
        sites[nucleus] = [m.sum() for m in masks]
        V, meta = read_wilson_lines(path)
        D[nucleus, :, ix] = regional_dipole(V, masks, rmax)
        for k in range(len(REGIONS)):
            qs[nucleus, k, ix] = np.sqrt(2.) / (
                crossing(D[nucleus, k, ix]) * meta["a_fm"]) * HBARC
    observables, geometry = read_run(workdir)
    out = {"qs_GeV": qs.tolist(), "D": D.tolist(), "sites": sites.tolist(),
           "observables": observables, "b": geometry["b"]}
    with open(cache, "w") as f:
        json.dump(out, f)
    if not settings["keep_files"]:
        for path in glob.glob(os.path.join(workdir, "WilsonLine*")) \
                + glob.glob(os.path.join(workdir, "*.ipgt")):
            os.remove(path)
    return cfg, out, None


def fit_linear(values, ds):
    """values[configuration, ds]: the fit linear in the step of the
    averages over the configurations, with the error from the spread of the
    fits of the single configurations (which share the initial Wilson
    lines): the deviation from zero step at each step, zero step, and the
    quadratic term at the largest step from a quadratic fit."""
    ds = np.asarray(ds, dtype=float)
    v = np.asarray(values, dtype=float)
    n = len(v)
    linear = np.linalg.lstsq(np.vander(ds, 2), v.T, rcond=None)[0]
    deviation = linear[0][:, None] * ds[None, :]
    quadratic = np.linalg.lstsq(np.vander(ds, 3), v.T, rcond=None)[0]
    q = quadratic[0] * ds.max() ** 2
    return {"deviation": deviation.mean(axis=0).tolist(),
            "deviation_error": (deviation.std(axis=0, ddof=1)
                                / np.sqrt(n)).tolist(),
            "zero": float(linear[1].mean()),
            "zero_error": float(linear[1].std(ddof=1) / np.sqrt(n)),
            "quadratic": float(q.mean()),
            "quadratic_error": float(q.std(ddof=1) / np.sqrt(n))}


def analyze(runs, settings):
    ds_list = sorted(settings["ds"], reverse=True)
    by_key = {(cfg["ds"], cfg["seed"]): out for cfg, out in runs}
    seeds = sorted({cfg["seed"] for cfg, _ in runs})
    seeds = [s for s in seeds if all((ds, s) in by_key for ds in ds_list)]
    _, n, xs = evolution(settings)
    dy = np.pi ** 2 * max(settings["ds"])
    parts = settings["snapshots"] + 1
    y = [k * (n // parts) * dy for k in range(parts)] + [n * dy]
    # [configuration (seed, nucleus), ds, region, y]
    qs = np.array([[np.array(by_key[(ds, s)]["qs_GeV"])[nucleus]
                    for ds in ds_list] for s in seeds for nucleus in (0, 1)])
    result = {"ds": ds_list, "y": y, "x": xs, "seeds": seeds,
              "configurations": len(qs),
              "initial_spread": float(np.nanmax(np.abs(
                  qs[:, :, :, 0] - qs[:, -1:, :, 0]))),
              "sites": by_key[(ds_list[0], seeds[0])]["sites"],
              "qs": {}, "observables": {}}
    for k, region in enumerate(REGIONS):
        lnq = np.log(qs[:, :, k])
        fits = [fit_linear(lnq[:, :, iy], ds_list)
                for iy in range(1, len(y))]
        for fit in fits:
            # Q_s(ds)/Q_s(0) - 1
            fit["deviation"] = np.expm1(fit["deviation"]).tolist()
            fit["zero"] = float(np.exp(fit["zero"]))
        result["qs"][region] = {
            "GeV": np.exp(lnq).mean(axis=0).tolist(), "fits": fits}
    for name in OBSERVABLES:
        v = np.array([[by_key[(ds, s)]["observables"][name]
                       for ds in ds_list] for s in seeds])
        relative = not name.startswith("eps")
        fit = fit_linear(np.log(v) if relative else v, ds_list)
        if relative:
            fit["deviation"] = np.expm1(fit["deviation"]).tolist()
            fit["zero"] = float(np.exp(fit["zero"]))
        else:
            # relative to the value at zero step
            for key in ("deviation", "deviation_error"):
                fit[key] = (np.array(fit[key]) / fit["zero"]).tolist()
        fit["values"] = v.mean(axis=0).tolist()
        result["observables"][name] = fit
    return result


def checks_for(result, settings):
    checks = {"same initial Wilson lines for all step sizes":
              bool(result["initial_spread"] < 1e-10)}
    fit = result["qs"]["all"]["fits"][-1]
    checks["ln Q_s (all) is consistent with first order"] = bool(
        abs(fit["quadratic"]) < 3. * fit["quadratic_error"])
    if settings["quick"]:
        return checks
    ds = result["ds"]
    tolerance = settings["max_deviation"]
    i_rec = ds.index(settings["recommended"])
    for region in REGIONS:
        dev = result["qs"][region]["fits"][-1]["deviation"][i_rec]
        checks[f"Q_s ({region}) deviates by less than {tolerance:g} at "
               f"jimwlkDs {settings['recommended']:g}"] = bool(
            abs(dev) < tolerance)
    checks[f"dN/dy deviates by less than {tolerance:g} at jimwlkDs "
           f"{settings['recommended']:g}"] = bool(abs(
               result["observables"]["dN/dy"]["deviation"][i_rec])
               < tolerance)
    if settings["energy_ds"] in ds:
        i_energy = ds.index(settings["energy_ds"])
        checks[f"dE/deta deviates by less than {tolerance:g} at jimwlkDs "
               f"{settings['energy_ds']:g}"] = bool(abs(
                   result["observables"]["dE/deta(final)"]["deviation"]
                   [i_energy]) < tolerance)
    return checks


def print_result(result):
    ds = result["ds"]
    print(f"\n{result['configurations']} configurations of Q_s (seeds "
          f"{', '.join(map(str, result['seeds']))}), initial Wilson lines "
          f"equal to {result['initial_spread']:.1e}")
    print("Q_s [GeV] at y = " + ", ".join(f"{v:.2f}" for v in result["y"])
          + f" (jimwlkDs {ds[-1]:g}); deviation from zero step at the final "
          "y at jimwlkDs " + ", ".join(f"{d:g}" for d in ds))
    for region in REGIONS:
        q = np.array(result["qs"][region]["GeV"])[-1]
        fit = result["qs"][region]["fits"][-1]
        print(f"  {region:7s}" + " ".join(f"{v:6.3f}" for v in q) + "   "
              + "  ".join(f"{100 * d:+5.2f}% +- {100 * e:.2f}%" for d, e
                          in zip(fit["deviation"], fit["deviation_error"])))
    print("deviation of Q_s (all) at jimwlkDs " f"{ds[0]:g} against y: "
          + ", ".join(f"{100 * f['deviation'][0]:+.2f}% +- "
                      f"{100 * f['deviation_error'][0]:.2f}% (y = {yy:.2f})"
                      for f, yy in zip(result["qs"]["all"]["fits"],
                                       result["y"][1:])))
    print("observables (per seed), deviation from zero step:")
    for name, fit in result["observables"].items():
        print(f"  {name:15s} {fit['zero']:10.4g}   "
              + "  ".join(f"{100 * d:+6.2f}% +- {100 * e:.2f}%" for d, e
                          in zip(fit["deviation"], fit["deviation_error"])))


def plot(result, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a8a85"]
    ds = np.array(result["ds"])
    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 3.8))
    for (name, fits), color in zip(
            ((r, result["qs"][r]["fits"][-1]) for r in REGIONS), colors):
        left.errorbar(ds, 100 * np.array(fits["deviation"]),
                      100 * np.array(fits["deviation_error"]), fmt="o-",
                      ms=3, color=color, label=f"Q_s, {name}")
    for (name, fit), color in zip(result["observables"].items(), colors):
        right.errorbar(ds, 100 * np.array(fit["deviation"]),
                       100 * np.array(fit["deviation_error"]), fmt="o-",
                       ms=3, color=color, label=name)
    for ax in (left, right):
        ax.axhline(0., color="k", lw=0.6)
        ax.set_xlim(left=0.)
        ax.set_xlabel("jimwlkDs")
        ax.set_ylabel("deviation from zero step [%]")
        ax.legend(fontsize=7)
    left.set_title(f"y = {result['y'][-1]:.2f}", fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="The JIMWLK step size in the example input (issue "
                    "#128).")
    parser.add_argument("--ipglasma-path", default="../..")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--template", default=os.path.join(HERE, "..", "..",
                                                           "input"))
    parser.add_argument("--datadir", default="./example_input_step_size_output")
    parser.add_argument("--size", type=int, default=768)
    parser.add_argument("--L", type=float, default=22.5,
                        help="[fm]; 22.5 fm suffices at b = 0")
    parser.add_argument("--b", type=float, default=0., help="[fm]")
    parser.add_argument("--ds", type=float, nargs="+",
                        default=[0.005, 0.0025, 0.00125, 0.000625],
                        help="step sizes; the largest must be a multiple of "
                             "each")
    parser.add_argument("--x0", type=float, default=0.01,
                        help="the nominal jimwlkInitialX; the evolution is "
                             "rounded to multiples of the largest step")
    parser.add_argument("--x", type=float, default=1e-4,
                        help="projectileX and targetX")
    parser.add_argument("--snapshots", type=int, default=2)
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--recommended", type=float, default=0.005)
    parser.add_argument("--energy-ds", type=float, default=0.0025,
                        help="the step at which dE/deta must be within "
                             "--max-deviation")
    parser.add_argument("--max-deviation", type=float, default=0.03)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--quick", action="store_true",
                        help="128^2 at L = 30 fm, step sizes 0.005, 0.0025 "
                             "and 0.00125, 2 seeds")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--keep-files", action="store_true")
    args = parser.parse_args()
    if args.quick:
        args.size, args.L, args.seeds = 128, 30., 2
        args.ds = [0.005, 0.0025, 0.00125]
    ratios = [max(args.ds) / ds for ds in args.ds]
    if any(abs(q - round(q)) > 1e-9 for q in ratios):
        parser.error("the largest step must be a multiple of each step")

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
            "rotateReactionPlane": 0, "useFluctuatingX": 0,
            "useJIMWLK": 1, "projectileX": args.x, "targetX": args.x,
            "jimwlkSaveSnapshots": 1, "writeWilsonLines": 2,
            "writeWilsonLineGeometry": 1, "writeTmunu": 1,
            "writeTmunuBinary": 1, "writeHydro": 0, "outputTimes": "0.2",
            "computeEccentricities": 1, "computeGluonMultiplicity": 1,
            "nucleusQsTableFileName": os.path.join(
                path, "qs2Adj_vs_Tp_vs_Y_240.in"),
            "nuclearConfigurationsPath": os.path.join(
                path, "nucleusConfigurations")}})
    os.makedirs(settings["datadir"], exist_ok=True)
    x0, n, _ = evolution(settings)
    print(f"evolution: {n} steps of {max(args.ds):g} from x0 = {x0:.5g} to "
          f"{args.x:g}")

    # the long runs first
    configs = [{"ds": ds, "seed": seed} for ds, seed in
               itertools.product(sorted(args.ds), range(1, args.seeds + 1))]
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
                           "example_input_step_size_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    if runs:
        plot(summary["result"], os.path.join(
            settings["datadir"], "example_input_step_size.pdf"))
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{settings['datadir']}/example_input_step_size_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

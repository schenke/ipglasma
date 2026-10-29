#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
The JIMWLK evolution in the setup of Lappi and Mantysaari, arXiv:1212.4825
(issue #128): the step-size dependence, the evolution speed against the
curves of that paper, and BK with the same kernel.

The setup of the paper (its curves with the lower initial Q_s): the MV
model (useNucleus 0) with g^2 mu L = 31, Ny = 100, m = 0 and no UV damping
on a 1024^2 lattice, evolved with the square-root running coupling, the
prescription IP-Glasma implements (jimwlkAlphaS 0), with Lambda L = 6,
mu0 L = 15, c = 0.2, N_f = 3 and no IR regulator in the kernel
(jimwlkMass 0). The coupling of the paper has the scale 4 e^{-2 gamma_E}/r^2
instead of 4/r^2, so jimwlkLambdaQCD and jimwlkMu0 are e^{gamma_E} times its
Lambda and mu0. Everything is in units of L; L only sets the units of the
input.

All step sizes evolve the same seeds, i.e. the same initial Wilson lines
(both nuclei of a run, two configurations), to y = ln(x0/x) = --ymax, with
snapshots at every unit of y. These are rounded to multiples of the largest
step, so all step sizes take them at the same x. From the dipole
correlator, Q_s = sqrt(2)/r_s with D(r_s) = e^{-1/2}, in units of Lambda.

At each y, ln Q_s is fitted against the step with fit_step_size of
jimwlk_step_size.py (first order, with the covariance over the
configurations); the configurations at zero step are those of the smallest
step minus the fitted slope times it. The evolution speed lambda =
d ln Q_s^2/dy comes from the differences between successive snapshots. As a
reference, bk_sqrt.py solves the BK equation with the same kernel from the
mean dipole correlator of the initial Wilson lines; lambda of BK is computed
with the same Delta y. The curves of the paper (lappi_mantysaari_2012/) are
Q_s(y) from Fig. 6 and lambda of JIMWLK (Fig. 8) and BK (Fig. 11), all with
the square-root prescription and c = 0.2.

The script checks that

  1. the initial Wilson lines are the same for all step sizes,
  2. ln Q_s at the final y is consistent with a first-order step (p-value
     above --min-p-value),
  3. the BK solution reproduces lambda of Fig. 11 within --bk-tolerance
     (this tests the coupling and its scale, independent of the lattice),
  4. at zero step, JIMWLK is slower than BK: Q_s at the final y is below
     that of BK, as in the paper,
  5. Q_s at zero step agrees with Fig. 6 within --fig6-tolerance (not with
     --quick, whose lattice is coarser than that of the paper).

With the defaults (step sizes 0.001, 0.0005 and 0.00025 to y = 8, 2 seeds),
the step error is about 30 times larger per unit of step and of rapidity
than in jimwlk_step_size.py (no IR regulator, a coupling that freezes at
0.76 instead of 0.36); the paper's curves match the smallest step.

Usage (from validations/jimwlk_step_size/):
    # quick test: 256^2, to y = 4, step sizes 0.002, 0.001, 0.0005, 2 seeds
    python3 lappi_mantysaari.py --ipglasma-path ../.. --quick
    # the setup of the paper (6 runs, 3.6 hours with 4 workers of 6
    # threads)
    python3 lappi_mantysaari.py --ipglasma-path ../.. --threads 6 \\
        --max-workers 4
"""

import argparse
import csv
import glob
import itertools
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# (jimwlk_step_size puts validations/dipole_correlator on the path)
from jimwlk_step_size import FILE_RE, fit_step_size  # noqa: E402
from dipole_correlator import (HBARC, crossing, dipole_correlator,  # noqa
                               profiles, write_input)
from read_wilson_lines import read_wilson_lines  # noqa: E402
import bk_sqrt  # noqa: E402

REFERENCE = os.path.join(HERE, "lappi_mantysaari_2012")


def snapshot_y(settings):
    """y of the initial Wilson lines, the snapshots and the final ones: the
    integers up to ymax, rounded to multiples of the largest step."""
    dy = np.pi ** 2 * max(settings["ds"])
    return np.array([round(k / dy) * dy for k in range(settings["ymax"] + 1)])


def run_tag(cfg):
    return f"ds{cfg['ds']:g}_seed{cfg['seed']}"


def run_one(cfg, settings):
    """Run IP-Glasma for one step size and seed and return the D(r)
    profiles (axis) of both nuclei at every snapshot (cached in
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
    xs = settings["x0"] * np.exp(-snapshot_y(settings))
    L = settings["L"]
    overrides = {
        "size": settings["size"], "L": L,
        "g2muGeV": f"{settings['g2mu_L'] * HBARC / L:.12g}",
        "Ny": settings["ny"], "m": 0, "UVDamp": 0, "seed": cfg["seed"],
        "useNucleus": 0, "runEvolution": 0, "writeWilsonLines": 2,
        "writeWilsonLineSnapshot": 0, "useFluctuatingX": 0,
        "nucleusQsTableFileName": settings["qs_table"],
        "useJIMWLK": 1, "jimwlkInitialX": f"{xs[0]:.12g}",
        "projectileX": f"{xs[-1]:.12g}", "targetX": f"{xs[-1]:.12g}",
        "jimwlkDs": f"{cfg['ds']:.12g}", "jimwlkAlphaS": 0,
        "jimwlkMass": 0, "jimwlkC": settings["c"], "nFlavors": settings["nf"],
        "jimwlkLambdaQCD": f"{settings['lambda_L'] * HBARC / L * np.exp(np.euler_gamma):.12g}",
        "jimwlkMu0": f"{settings['mu0_L'] * HBARC / L * np.exp(np.euler_gamma):.12g}",
        "jimwlkSaveSnapshots": 1,
        "jimwlkXSnapshotList": ",".join(f"{x:.12g}" for x in xs[1:-1]),
    }
    write_input(settings["template"], overrides,
                os.path.join(workdir, "input"))
    env = dict(os.environ, OMP_NUM_THREADS=str(settings["threads"]))
    with open(os.path.join(workdir, "run.log"), "w") as log:
        result = subprocess.run([settings["ipglasma"], "input"], cwd=workdir,
                                stdout=log, stderr=subprocess.STDOUT,
                                env=env).returncode
    files = glob.glob(os.path.join(workdir, "WilsonLine_x_*"))
    if result != 0 or len(files) != 2 * len(xs):
        return cfg, None, (f"ipglasma failed (exit {result}, {len(files)} of "
                           f"{2 * len(xs)} Wilson-line files), see {workdir}")
    out = {"axis": np.zeros((len(xs), 2, settings["size"] // 2 + 1))}
    for path in files:
        match = FILE_RE.search(os.path.basename(path))
        ix = int(np.argmin(np.abs(np.log(xs / float(match[1])))))
        nucleus = (int(match[2]) - 1) % 2
        V, meta = read_wilson_lines(path)
        out["axis"][ix, nucleus] = profiles(dipole_correlator(V))["axis"]
        out["a_fm"] = np.array(meta["a_fm"])
        if not settings["keep_files"]:
            os.remove(path)
    np.savez(cache, **out)
    return cfg, out, None


def read_curves(name):
    curves = {}
    with open(os.path.join(REFERENCE, name)) as f:
        for row in csv.DictReader(f):
            curves.setdefault(row["curve"], []).append(
                (float(row["x"]), float(row["y"])))
    return {k: np.array(sorted(v)) for k, v in curves.items()}


def reference_curves(c):
    """The square-root curves of arXiv:1212.4825 for c = 0.2: Q_s/Lambda at
    y = 0, 4, ..., 20 from N(r_s) = 1 - e^{-1/2} of Fig. 6, and lambda
    against Q_s/Lambda of JIMWLK (Fig. 8) and BK (Fig. 11)."""
    if abs(c - 0.2) > 1e-12:
        return None
    level = 1. - np.exp(-0.5)
    fig6 = {}
    for name, curve in read_curves("fig6_dipole.csv").items():
        if not name.startswith("sqrtas_JIMWLK_y"):
            continue
        i = int(np.argmax(curve[:, 1] > level))
        lnr = np.interp(level, curve[i - 1:i + 1, 1],
                        np.log(curve[i - 1:i + 1, 0]))
        fig6[int(name.split("_y")[1])] = float(np.sqrt(2.) / np.exp(lnr))
    return {"fig6_qs": dict(sorted(fig6.items())),
            "fig8_jimwlk": read_curves("fig8_lambda.csv")[
                "sqrtas_JIMWLK_c0.2_lowQs0"].tolist(),
            "fig11_bk": read_curves("fig11_lambda.csv")[
                "BK_sqrtas_c0.2"].tolist()}


def speed(y, lnq):
    """lambda = d ln Q_s^2/dy between successive snapshots, per
    configuration (lnq[configuration, y]), at the geometric mean of Q_s;
    mean and error over the configurations."""
    lam = 2. * np.diff(lnq, axis=1) / np.diff(y)
    q = np.exp(0.5 * (lnq[:, 1:] + lnq[:, :-1])).mean(axis=0)
    n = len(lnq)
    return {"qs": q.tolist(), "lambda": lam.mean(axis=0).tolist(),
            "lambda_error": (lam.std(axis=0, ddof=1) / np.sqrt(n)).tolist()}


def analyze(runs, settings):
    ds_list = sorted(settings["ds"], reverse=True)
    by_key = {(cfg["ds"], cfg["seed"]): out for cfg, out in runs}
    seeds = sorted({cfg["seed"] for cfg, _ in runs})
    seeds = [s for s in seeds if all((ds, s) in by_key for ds in ds_list)]
    y = snapshot_y(settings)
    a = float(runs[0][1]["a_fm"])
    lam_unit = settings["lambda_L"] * HBARC / settings["L"]
    r = np.arange(settings["size"] // 2 + 1)
    # [configuration, ds, y, r]
    profile = np.array([[by_key[(ds, s)]["axis"][:, nucleus]
                         for ds in ds_list]
                        for s in seeds for nucleus in (0, 1)])
    lnq = np.log(np.array([[[np.sqrt(2.) / (crossing(r, p) * a) * HBARC
                             / lam_unit for p in conf_ds]
                            for conf_ds in conf] for conf in profile]))
    result = {"y": y.tolist(), "ds": ds_list, "configurations": len(lnq),
              "initial_spread": float(np.abs(
                  profile[:, :, 0] - profile[:, -1:, 0]).max()),
              "ln_qs": {f"{ds:g}": lnq[:, i].tolist()
                        for i, ds in enumerate(ds_list)},
              "fits": []}
    zero = lnq[:, -1].copy()
    for iy in range(1, len(y)):
        fit = fit_step_size(lnq[:, :, iy], ds_list)
        zero[:, iy] -= fit["slope"] * ds_list[-1]
        result["fits"].append({
            "y": float(y[iy]), "ln_qs_zero": fit["zero"],
            "ln_qs_zero_error": fit["zero_error"], "slope": fit["slope"],
            "slope_error": fit["slope_error"], "chi2": fit["chi2"],
            "dof": fit["dof"], "p_value": fit["p_value"],
            "order": fit["order"], "order_interval": fit["order_interval"],
            # Q_s(ds)/Q_s(0) - 1 from the fit
            "deviation": [float(np.expm1(fit["slope"] * ds))
                          for ds in ds_list]})
    result["ln_qs_zero"] = zero.tolist()
    result["speed"] = {f"{ds:g}": speed(y, lnq[:, i])
                       for i, ds in enumerate(ds_list)}
    result["speed"]["0"] = speed(y, zero)

    # BK with the same kernel from the mean initial dipole correlator
    D0 = profile[:, -1, 0].mean(axis=0)
    yb, qb = bk_sqrt.solve(r[1:] * a * lam_unit / HBARC, 1. - D0[1:],
                           y[-1] + 0.5, settings["c"],
                           settings["mu0_L"] / settings["lambda_L"],
                           settings["nf"], **settings["bk_grid"])
    lnqb = np.log(qb)
    # lambda with the Delta y of the snapshots, centered at each y of BK
    half = 0.5 * np.mean(np.diff(y))
    inner = (yb >= half) & (yb <= yb[-1] - half)
    lam_bk = (2. * (np.interp(yb[inner] + half, yb, lnqb)
                    - np.interp(yb[inner] - half, yb, lnqb)) / (2. * half))
    result["bk"] = {"y": yb[::5].tolist(), "qs": qb[::5].tolist(),
                    "qs_at_snapshots": np.exp(np.interp(y, yb, lnqb)).tolist(),
                    "speed": {"qs": qb[inner][::5].tolist(),
                              "lambda": lam_bk[::5].tolist()}}
    result["reference"] = reference_curves(settings["c"])
    return result


def checks_for(result, settings):
    checks = {"same initial Wilson lines for all step sizes":
              bool(result["initial_spread"] < 1e-10)}
    final = result["fits"][-1]
    checks["ln Q_s at the final y is consistent with first order"] = bool(
        final["p_value"] > settings["min_p_value"])
    zero_final = np.exp(np.array(result["ln_qs_zero"])[:, -1])
    checks["at zero step, JIMWLK is slower than BK"] = bool(
        zero_final.mean() < result["bk"]["qs_at_snapshots"][-1])
    ref = result["reference"]
    if ref is None:
        return checks
    # BK against Fig. 11 where both have points
    fig11 = np.array(ref["fig11_bk"])
    q = np.array(result["bk"]["speed"]["qs"])
    lam = np.array(result["bk"]["speed"]["lambda"])
    inside = (q > fig11[:, 0].min()) & (q < fig11[:, 0].max())
    deviation = lam[inside] / np.interp(np.log(q[inside]),
                                        np.log(fig11[:, 0]), fig11[:, 1]) - 1.
    result["bk_fig11_max_deviation"] = float(np.abs(deviation).max())
    checks[f"BK reproduces Fig. 11 within {settings['bk_tolerance']:g}"] = \
        bool(result["bk_fig11_max_deviation"] < settings["bk_tolerance"])
    if not settings["quick"]:
        y = np.array(result["y"])
        zero = np.exp(np.array(result["ln_qs_zero"])).mean(axis=0)
        devs = {k: float(zero[np.argmin(np.abs(y - k))] / v - 1.)
                for k, v in ref["fig6_qs"].items()
                if 0 < k <= y[-1] + 0.5}
        result["fig6_deviation"] = devs
        checks[f"Q_s at zero step agrees with Fig. 6 within "
               f"{settings['fig6_tolerance']:g}"] = bool(all(
                   abs(d) < settings["fig6_tolerance"] for d in devs.values()))
    return checks


def print_result(result):
    print(f"\n{result['configurations']} configurations, initial Wilson "
          f"lines equal to {result['initial_spread']:.1e}")
    ds = result["ds"]
    print("Q_s/Lambda at y = " + ", ".join(f"{v:.2f}" for v in result["y"]))
    for k in ds:
        q = np.exp(np.array(result["ln_qs"][f"{k:g}"])).mean(axis=0)
        print(f"  ds {k:<9g}" + " ".join(f"{v:7.2f}" for v in q))
    zero = np.exp(np.array(result["ln_qs_zero"])).mean(axis=0)
    print(f"  {'ds -> 0':12s}" + " ".join(f"{v:7.2f}" for v in zero))
    print(f"  {'BK':12s}" + " ".join(f"{v:7.2f}" for v in
                                     result["bk"]["qs_at_snapshots"]))
    final = result["fits"][-1]
    print(f"at y = {final['y']:.3f}: Q_s(ds)/Q_s(0) - 1 = "
          + ", ".join(f"{d:+.4f} (ds {k:g})"
                      for d, k in zip(final["deviation"], ds))
          + f"; first order: chi2 {final['chi2']:.2f}/{final['dof']} "
          f"(p = {final['p_value']:.2f}), order {final['order']:.2f} "
          f"[{final['order_interval'][0]:.2f}, "
          f"{final['order_interval'][1]:.2f}]")
    ref = result["reference"]
    fig8 = np.array(ref["fig8_jimwlk"]) if ref else None
    fig11 = np.array(ref["fig11_bk"]) if ref else None
    bk = result["bk"]["speed"]
    print("lambda = d ln Q_s^2/dy:  Q_s/Lambda   ds -> 0          "
          + "".join(f"ds {k:<8g}" for k in ds)
          + "BK     " + ("Fig. 8  Fig. 11" if ref else ""))
    for i, q in enumerate(result["speed"]["0"]["qs"]):
        def at(curve):
            return np.interp(np.log(q), np.log(curve[:, 0]), curve[:, 1],
                             left=np.nan, right=np.nan)
        z = result["speed"]["0"]
        print(f"  {q:24.2f}   {z['lambda'][i]:.3f} +- "
              f"{z['lambda_error'][i]:.3f}  "
              + "".join(f"{result['speed'][f'{k:g}']['lambda'][i]:<11.3f}"
                        for k in ds)
              + f"{np.interp(np.log(q), np.log(bk['qs']), bk['lambda'], left=np.nan, right=np.nan):.3f}  "
              + (f"{at(fig8):.3f}   {at(fig11):.3f}" if ref else ""))
    if ref:
        print("Fig. 6: Q_s/Lambda = " + ", ".join(
            f"{v:.2f} (y = {k})" for k, v in ref["fig6_qs"].items()))
    if "bk_fig11_max_deviation" in result:
        print(f"BK against Fig. 11: largest deviation in lambda "
              f"{100 * result['bk_fig11_max_deviation']:.2f}%")
    if "fig6_deviation" in result:
        print("Q_s at zero step against Fig. 6: " + ", ".join(
            f"{100 * v:+.1f}% (y = {k})"
            for k, v in result["fig6_deviation"].items()))


def plot(result, path):
    """(a) lambda against Q_s/Lambda, (b) Q_s/Lambda at the final y
    against the step."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ramp = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
    ref = result["reference"]
    if ref:
        fig11 = np.array(ref["fig11_bk"])
        fig8 = np.array(ref["fig8_jimwlk"])
        left.plot(fig11[:, 0], fig11[:, 1], color="#1baf7a", lw=3,
                  label="arXiv:1212.4825 BK (Fig. 11)")
        left.plot(fig8[:, 0], fig8[:, 1], color="#eb6834", lw=1.5,
                  label="arXiv:1212.4825 JIMWLK (Fig. 8)")
    bk = result["bk"]["speed"]
    left.plot(bk["qs"], bk["lambda"], color="#eda100", lw=1.2, label="BK")
    for k, color in zip(result["ds"], ramp):
        s = result["speed"][f"{k:g}"]
        left.plot(s["qs"], s["lambda"], "o-", ms=3, lw=1, color=color,
                  label=f"jimwlkDs {k:g}")
    z = result["speed"]["0"]
    left.errorbar(z["qs"], z["lambda"], z["lambda_error"], fmt="D",
                  color="#2a78d6", ms=4, label="ds -> 0")
    left.set_xscale("log")
    left.set_xlabel(r"$Q_s/\Lambda$")
    left.set_ylabel(r"$\lambda = d\ln Q_s^2/dy$")
    left.legend(fontsize=7)
    ds = np.array(result["ds"])
    lnq = [np.array(result["ln_qs"][f"{k:g}"])[:, -1] for k in ds]
    q = np.array([np.exp(v).mean() for v in lnq])
    e = np.array([np.exp(v).std(ddof=1) / np.sqrt(len(v)) for v in lnq])
    final = result["fits"][-1]
    line = np.linspace(0., ds.max(), 50)
    right.plot(line, np.exp(final["ln_qs_zero"] + final["slope"] * line),
               color="#8a8a85", lw=0.8)
    right.errorbar(ds, q, e, fmt="o", color="#1c5cab")
    zero = np.exp(np.array(result["ln_qs_zero"])[:, -1])
    right.errorbar([0.], [zero.mean()], [zero.std(ddof=1) / np.sqrt(len(zero))],
                   fmt="D", color="#2a78d6", label="ds -> 0")
    right.axhline(result["bk"]["qs_at_snapshots"][-1], color="#eda100",
                  label="BK")
    y_final = result["y"][-1]
    if ref and any(abs(k - y_final) < 0.5 for k in ref["fig6_qs"]):
        k = min(ref["fig6_qs"], key=lambda k: abs(k - y_final))
        right.axhline(ref["fig6_qs"][k], color="#eb6834", ls="--",
                      label="arXiv:1212.4825 (Fig. 6)")
    right.set_xlabel("jimwlkDs")
    right.set_ylabel(rf"$Q_s/\Lambda$ at $y = {y_final:.2f}$")
    right.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="JIMWLK in the setup of arXiv:1212.4825: step size, "
                    "evolution speed and BK (issue #128).")
    parser.add_argument("--ipglasma-path", default="../..")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--template", default=os.path.join(HERE, "..", "..",
                                                           "input"))
    parser.add_argument("--datadir", default="./lappi_mantysaari_output")
    parser.add_argument("--size", type=int, default=1024)
    parser.add_argument("--L", type=float, default=10.,
                        help="[fm]; only sets the units of the input")
    parser.add_argument("--g2mu-L", type=float, default=31.)
    parser.add_argument("--ny", type=int, default=100)
    parser.add_argument("--lambda-L", type=float, default=6.,
                        help="Lambda_QCD L in the convention of the paper")
    parser.add_argument("--mu0-L", type=float, default=15.)
    parser.add_argument("--c", type=float, default=0.2)
    parser.add_argument("--nf", type=int, default=3)
    parser.add_argument("--x0", type=float, default=0.01,
                        help="jimwlkInitialX; only y = ln(x0/x) matters")
    parser.add_argument("--ds", type=float, nargs="+",
                        default=[0.001, 0.0005, 0.00025],
                        help="step sizes; the largest must be a multiple of "
                             "each")
    parser.add_argument("--ymax", type=int, default=8)
    parser.add_argument("--seeds", type=int, default=2,
                        help="runs per step size (two configurations each)")
    parser.add_argument("--bk-tolerance", type=float, default=0.02,
                        help="largest relative deviation of lambda of BK "
                             "from Fig. 11")
    parser.add_argument("--fig6-tolerance", type=float, default=0.1,
                        help="largest relative deviation of Q_s at zero "
                             "step from Fig. 6")
    parser.add_argument("--min-p-value", type=float, default=0.01)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--quick", action="store_true",
                        help="256^2 with the same L, to y = 4, step sizes "
                             "0.002, 0.001 and 0.0005")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--keep-files", action="store_true")
    args = parser.parse_args()
    if args.quick:
        args.size, args.ymax, args.ds = 256, 4, [0.002, 0.001, 0.0005]
    ratios = [max(args.ds) / ds for ds in args.ds]
    if any(abs(q - round(q)) > 1e-9 for q in ratios):
        parser.error("the largest step must be a multiple of each step")

    path = os.path.abspath(args.ipglasma_path)
    settings = dict(vars(args))
    settings.update({
        "ipglasma": os.path.join(path, args.ipglasma_cmd),
        "qs_table": os.path.join(path, "qs2Adj_vs_Tp_vs_Y_240.in"),
        "template": os.path.abspath(args.template),
        "datadir": os.path.abspath(args.datadir),
        # a coarser grid for --quick (9 s instead of 3 min to y = 4.5; Q_s
        # within 0.2%)
        "bk_grid": dict(nr=200, ntheta=16, dy=0.04) if args.quick else {}})
    os.makedirs(settings["datadir"], exist_ok=True)

    configs = [{"ds": ds, "seed": seed} for ds, seed in
               itertools.product(sorted(args.ds), range(1, args.seeds + 1))]
    print(f"{len(configs)} runs of size {args.size} to y = {args.ymax} with "
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
                           "lappi_mantysaari_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    if runs:
        plot(summary["result"], os.path.join(settings["datadir"],
                                             "lappi_mantysaari.pdf"))
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{settings['datadir']}/lappi_mantysaari_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

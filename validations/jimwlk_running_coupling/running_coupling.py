#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Validation of the JIMWLK running coupling: the "noise" coupling
(jimwlkRunningCoupling 2) against the "square root" coupling
(jimwlkRunningCoupling 1), reproducing Fig. 6 of T. Lappi and
H. Mantysaari, Eur. Phys. J. C 73 (2013) 2307, arXiv:1212.4825.

The square-root coupling multiplies the kernel by sqrt(alpha_s(r)); the
noise coupling instead gives the noise the correlator
<eta eta> = int_k e^{ik.(x-y)} alpha_s(k) (Eqs. (24)-(25) there), which
the paper finds to evolve significantly slower. The setup is that of the
paper, in units of the lattice size L:

  * 1024^2 lattice, MV model with g^2 mu L = 31 in Ny = 100 layers and no
    infrared mass (Q_s L ~ 22 initially),
  * Lambda_QCD L = 6, mu0 L = 15 (alpha_s frozen at 0.76), c = 0.2, Nf = 3,
  * alpha_s(k) of Eq. (39) for the noise coupling (the code uses the
    lattice momentum), alpha_s(r) of Eq. (40), with the scale
    4 e^{-2 gamma_E} / r^2, for the square-root coupling.

The code's alpha_s(r) has the scale 4 / r^2; the scale of Eq. (40) is
obtained by running the square-root coupling with Lambda_QCD and mu0 both
multiplied by e^{gamma_E}, which leaves alpha_s(r -> infinity) unchanged.
The prescription "sqrt4" runs the code's scale instead, for comparison.

Each run evolves the two nuclei, two independent configurations, over
`--ymax` units of rapidity and writes the Wilson lines every
`--dy-snapshot` units. The script computes the dipole correlator
D(r) = Re < tr V^+(x) V(x + r) > / Nc along the lattice axes, averaged over
all sites, deletes the Wilson-line files, and checks that

  1. both prescriptions start from the same Q_s,
  2. the noise coupling gives a smaller Q_s than the square-root coupling
     at every later rapidity, and
  3. a smaller evolution speed lambda = d ln Q_s^2 / dy, averaged over the
     evolution,

with D(r_s) = exp(-1/2) and Q_s^2 = 2/r_s^2 as in the paper. The results
go to <datadir>/running_coupling_summary.json, the figures to
<datadir>/running_coupling_*.pdf.

Usage (from validations/jimwlk_running_coupling/):
    # a quick test on a 256^2 lattice over 8 units, 4 minutes on 8 cores
    python3 running_coupling.py --ipglasma-path ../.. --quick
    # the setup of the paper, about 2 hours per run on 8 cores
    python3 running_coupling.py --ipglasma-path ../.. --seeds 2

The default step jimwlkDs = 0.0005 is converged on 256^2: with the
0.005 recommended in the README, Q_s at y = 8 comes out about 25% too
large there (the growth of ln Q_s^2 10-15% too fast) for both
prescriptions, here with alpha_s frozen at 0.76.

--ipglasma-path is a directory with the ipglasma executable and
qs2Adj_vs_Tp_vs_Y_240.in (e.g. the repository root after
./compile_IPGlasma.sh).
"""

import argparse
import glob
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "utilities"))
from read_wilson_lines import read_wilson_lines  # noqa: E402

HBARC = 0.1973269804  # GeV fm
NC = 3
X0 = 0.01  # jimwlkInitialX, y = ln(X0 / x)
EULER_GAMMA = 0.5772156649015329

# jimwlkRunningCoupling and the factor multiplying Lambda_QCD and mu0
PRESCRIPTIONS = {
    "noise": {"jimwlkRunningCoupling": 2, "scale": 1.},
    "sqrt": {"jimwlkRunningCoupling": 1, "scale": math.exp(EULER_GAMMA)},
    "sqrt4": {"jimwlkRunningCoupling": 1, "scale": 1.},
}
LABELS = {
    "noise": r"Noise $\alpha_s$",
    "sqrt": r"$\sqrt{\alpha_s}$",
    "sqrt4": r"$\sqrt{\alpha_s}$, scale $4/r^2$",
}


# --------------------------------------------------------------------------
# physics
# --------------------------------------------------------------------------

def dipole_axis_profile(V):
    """D(r) = Re sum_x tr V^+(x) V(x + r) / (Nc N^2) along the lattice axes
    (mean of +-x, +-y), for r = 0 .. N/2 in lattice units; V has shape
    (N, N, 3, 3)."""
    n = V.shape[0]
    total = np.zeros((n, n))
    for a in range(NC):
        for b in range(NC):
            f = np.fft.fft2(V[:, :, a, b])
            # ifft2(|f_k|^2)[r] = sum_x conj(V(x)) V(x + r)
            total += np.fft.ifft2(np.conj(f) * f).real
    d2 = total / (NC * n * n)
    r = np.arange(n // 2 + 1)
    return 0.25 * (d2[r, 0] + d2[-r, 0] + d2[0, r] + d2[0, -r])


def crossing(r, d, level=math.exp(-0.5)):
    """The r at which the decreasing profile d first falls below level,
    linearly interpolated; NaN if it does not."""
    below = np.nonzero(d < level)[0]
    if below.size == 0 or below[0] == 0:
        return float("nan")
    i = below[0]
    return float(r[i - 1] + (d[i - 1] - level) / (d[i - 1] - d[i]))


# --------------------------------------------------------------------------
# running IP-Glasma
# --------------------------------------------------------------------------

def write_input(template, overrides, path):
    """Copy the template, replacing the values of the keys in overrides."""
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


def snapshot_rapidities(settings):
    n = int(round(settings["ymax"] / settings["dy_snapshot"]))
    return [settings["dy_snapshot"] * k for k in range(1, n)]


def run_one(prescription, seed, settings):
    """Run IP-Glasma for one prescription and seed and return the D(r)
    profiles of its two nuclei at all rapidities (cached in
    <workdir>/profiles.npz)."""
    workdir = os.path.join(settings["datadir"], f"{prescription}_seed{seed}")
    cache = os.path.join(workdir, "profiles.npz")
    if os.path.exists(cache):
        with np.load(cache) as data:
            return {k: data[k] for k in data.files}
    if settings["analyze_only"]:
        raise RuntimeError(f"no cached result in {workdir} (--analyze-only)")
    os.makedirs(workdir, exist_ok=True)
    for old in glob.glob(os.path.join(workdir, "WilsonLine*")):
        os.remove(old)

    L = settings["lambda_L"] * HBARC / settings["lambda_GeV"]  # fm
    scale = PRESCRIPTIONS[prescription]["scale"]
    x_final = X0 * math.exp(-settings["ymax"])
    snapshots = ",".join(f"{X0 * math.exp(-y):.10e}"
                         for y in snapshot_rapidities(settings))
    overrides = {
        "size": settings["size"], "L": f"{L:.10g}", "Ny": settings["ny"],
        "g2muGeV": f"{settings['g2mu_L'] * HBARC / L:.10g}",
        "seed": seed,
        "nucleusQsTableFileName": settings["qs_table"],
        "projectileX": f"{x_final:.10e}", "targetX": f"{x_final:.10e}",
        "jimwlkLambdaQCD": f"{settings['lambda_GeV'] * scale:.10g}",
        "jimwlkMu0": f"{settings['mu0_L'] * HBARC / L * scale:.10g}",
        "jimwlkC": settings["c"], "jimwlkDs": settings["ds"],
        "jimwlkInitialX": X0,
        "jimwlkRunningCoupling":
            PRESCRIPTIONS[prescription]["jimwlkRunningCoupling"],
        "jimwlkXSnapshotList": snapshots,
    }
    write_input(settings["template"], overrides,
                os.path.join(workdir, "input"))
    env = dict(os.environ)
    if settings["threads"]:
        env["OMP_NUM_THREADS"] = str(settings["threads"])
    start = time.time()
    with open(os.path.join(workdir, "run.log"), "w") as log:
        result = subprocess.run([settings["ipglasma"], "input"], cwd=workdir,
                                stdout=log, stderr=subprocess.STDOUT, env=env)
    files = sorted(glob.glob(os.path.join(workdir, "WilsonLine_x_*")))
    expected = 2 * (len(snapshot_rapidities(settings)) + 2)
    if result.returncode != 0 or len(files) != expected:
        raise RuntimeError(f"ipglasma failed (exit {result.returncode}, "
                           f"{len(files)} of {expected} Wilson-line files), "
                           f"see {workdir}")
    print(f"  {prescription} seed {seed}: ran in {time.time() - start:.0f} s",
          flush=True)

    # group the files by nucleus (the trailing file number) and rapidity
    by_nucleus = {}
    for path in files:
        V, meta = read_wilson_lines(path)
        y = math.log(X0 / meta["x"])
        number = path.rsplit("_", 1)[1]
        by_nucleus.setdefault(number, []).append((y, dipole_axis_profile(V)))
        a_fm = meta["a_fm"]
        if not settings["keep_files"]:
            os.remove(path)
    out = {"a_fm": np.array(a_fm)}
    for i, number in enumerate(sorted(by_nucleus)):
        rows = sorted(by_nucleus[number], key=lambda t: t[0])
        out["y"] = np.array([t[0] for t in rows])
        out[f"D{i}"] = np.array([t[1] for t in rows])
    np.savez(cache, **out)
    return out


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------

def evolved_rapidity(y, ds):
    """The rapidity actually reached at the snapshot requested at y: the
    closest whole number of steps of pi^2 ds (see JIMWLK::snapshotSteps)."""
    dy = math.pi ** 2 * ds
    return np.floor(np.asarray(y) / dy + 0.5) * dy


def analyze(results, settings):
    """Mean N(r), Q_s(y) and lambda(y) of each prescription."""
    summary = {}
    lam_qcd_lat = settings["lambda_L"] / settings["size"]  # Lambda a
    for prescription, runs in results.items():
        y = evolved_rapidity(runs[0]["y"], settings["ds"])
        D = np.array([run[f"D{i}"] for run in runs for i in (0, 1)])
        mean = D.mean(axis=0)
        err = D.std(axis=0, ddof=1) / np.sqrt(len(D))
        r = np.arange(mean.shape[1])
        # Q_s of each configuration, for the error of Q_s(y)
        qs = np.array([[math.sqrt(2.) / crossing(r, d) / lam_qcd_lat
                        for d in conf] for conf in D])  # Q_s / Lambda
        log_qs2 = np.log(qs ** 2)
        lam = np.gradient(log_qs2, y, axis=1)
        summary[prescription] = {
            "configurations": len(D),
            "y": y.tolist(),
            "r_Lambda": (r * lam_qcd_lat).tolist(),
            "N": (1. - mean).tolist(),
            "N_error": err.tolist(),
            "Qs_over_Lambda": qs.mean(axis=0).tolist(),
            "Qs_over_Lambda_error":
                (qs.std(axis=0, ddof=1) / np.sqrt(len(qs))).tolist(),
            "lambda": lam.mean(axis=0).tolist(),
            "lambda_error":
                (lam.std(axis=0, ddof=1) / np.sqrt(len(lam))).tolist(),
            # the mean evolution speed, ln(Q_s^2(ymax)/Q_s^2(0)) / ymax
            "lambda_mean": float(np.mean(log_qs2[:, -1] - log_qs2[:, 0])
                                 / (y[-1] - y[0])),
        }
    checks = {}
    if "noise" in summary and "sqrt" in summary:
        noise, sqrt = summary["noise"], summary["sqrt"]
        qn = np.array(noise["Qs_over_Lambda"])
        qr = np.array(sqrt["Qs_over_Lambda"])
        en = np.array(noise["Qs_over_Lambda_error"])
        er = np.array(sqrt["Qs_over_Lambda_error"])
        checks["same_initial_Qs"] = bool(
            abs(qn[0] - qr[0]) < 3. * math.hypot(en[0], er[0]) + 1e-6 * qr[0])
        checks["noise_Qs_smaller"] = bool(np.all(qn[1:] < qr[1:]))
        checks["noise_slower"] = bool(noise["lambda_mean"] < sqrt["lambda_mean"])
        # the rapidity the noise coupling needs to reach the Q_s that the
        # square-root coupling reaches at ymax (interpolating ln Q_s^2)
        target = math.log(qr[-1] ** 2)
        ln_qn2 = np.log(qn ** 2)
        y = np.array(noise["y"])
        if ln_qn2[-1] >= target:
            y_needed = float(np.interp(target, ln_qn2, y))
        else:  # extrapolate with the last evolution speed
            y_needed = float(y[-1] + (target - ln_qn2[-1])
                             / max(noise["lambda"][-1], 1e-6))
        checks["noise_rapidity_for_sqrt_final_Qs"] = y_needed
    return {"settings": {k: v for k, v in settings.items()
                         if isinstance(v, (int, float, str, bool))},
            "prescriptions": summary, "checks": checks}


STYLES = {"noise": dict(color="k", ls="-"),
          "sqrt": dict(color="tab:blue", ls=":", lw=2.5),
          "sqrt4": dict(color="tab:orange", ls="--")}


def plot_fig6(summary, path, dy_plot):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    for prescription, s in summary["prescriptions"].items():
        y = np.array(s["y"])
        r = np.array(s["r_Lambda"])
        N = np.array(s["N"])
        shown = [i for i, yy in enumerate(y)
                 if abs(yy / dy_plot - round(yy / dy_plot)) < 0.05]
        for j, i in enumerate(shown):
            ax.plot(r[1:], N[i, 1:], **STYLES[prescription],
                    label=LABELS[prescription] if j == 0 else None)
    ax.set_xscale("log")
    ax.set_xlim(0.003, 1.5)
    ax.set_ylim(0., 1.)
    ax.set_xlabel(r"$r\Lambda_\mathrm{QCD}$")
    ax.set_ylabel(r"$N = 1 - \langle D(r) \rangle$")
    ax.set_title(fr"$y = 0, {dy_plot:g}, {2 * dy_plot:g}, \ldots$ "
                 "(right to left)", fontsize=10)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_speed(summary, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for prescription, s in summary["prescriptions"].items():
        style = dict(STYLES[prescription])
        style.pop("lw", None)
        left.errorbar(s["y"], s["Qs_over_Lambda"],
                      yerr=s["Qs_over_Lambda_error"], marker="o", ms=3,
                      label=LABELS[prescription], **style)
        right.errorbar(s["Qs_over_Lambda"], s["lambda"],
                       yerr=s["lambda_error"], marker="o", ms=3,
                       label=LABELS[prescription], **style)
    left.set_yscale("log")
    left.set_xlabel(r"$y$")
    left.set_ylabel(r"$Q_s/\Lambda_\mathrm{QCD}$")
    left.legend(frameon=False)
    right.set_xscale("log")
    right.set_xlabel(r"$Q_s/\Lambda_\mathrm{QCD}$")
    right.set_ylabel(r"$\lambda = d\ln Q_s^2/dy$")
    right.set_ylim(bottom=0.)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="JIMWLK noise vs square-root running coupling, Fig. 6 "
                    "of arXiv:1212.4825.")
    parser.add_argument("--ipglasma-path", default="../..",
                        help="directory with the ipglasma executable and "
                             "qs2Adj_vs_Tp_vs_Y_240.in")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--datadir", default="./running_coupling_output")
    parser.add_argument("--template",
                        default=os.path.join(HERE, "input_running_coupling"))
    parser.add_argument("--prescriptions", default="noise,sqrt",
                        help="comma-separated, of " + ", ".join(PRESCRIPTIONS))
    parser.add_argument("--seeds", type=int, default=1,
                        help="runs per prescription (two configurations "
                             "each)")
    parser.add_argument("--size", type=int, default=1024)
    parser.add_argument("--ny", type=int, default=100)
    parser.add_argument("--g2mu-L", type=float, default=31.)
    parser.add_argument("--lambda-L", type=float, default=6.)
    parser.add_argument("--mu0-L", type=float, default=15.)
    parser.add_argument("--lambda-GeV", type=float, default=0.2,
                        help="sets the units; the results depend only on "
                             "the ratios")
    parser.add_argument("--c", type=float, default=0.2)
    parser.add_argument("--ds", type=float, default=0.0005,
                        help="jimwlkDs; the recommended 0.005 evolves "
                             "gives a 25%% too large Q_s at y = 8 on 256^2")
    parser.add_argument("--ymax", type=float, default=20.)
    parser.add_argument("--dy-snapshot", type=float, default=1.)
    parser.add_argument("--dy-plot", type=float, default=4.)
    parser.add_argument("--threads", type=int, default=0,
                        help="OMP_NUM_THREADS of ipglasma (0: unset)")
    parser.add_argument("--quick", action="store_true",
                        help="256^2 lattice with the same physical size, "
                             "8 units of rapidity")
    parser.add_argument("--keep-files", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    args = parser.parse_args()

    if args.quick:
        args.size = 256
        args.ymax = 8.
        args.dy_snapshot = 2.
        args.dy_plot = 2.
    settings = {
        "ipglasma": os.path.abspath(os.path.join(args.ipglasma_path,
                                                 args.ipglasma_cmd)),
        "qs_table": os.path.abspath(os.path.join(
            args.ipglasma_path, "qs2Adj_vs_Tp_vs_Y_240.in")),
        "datadir": os.path.abspath(args.datadir),
        "template": os.path.abspath(args.template),
        "size": args.size, "ny": args.ny, "g2mu_L": args.g2mu_L,
        "lambda_L": args.lambda_L, "mu0_L": args.mu0_L,
        "lambda_GeV": args.lambda_GeV, "c": args.c, "ds": args.ds,
        "ymax": args.ymax, "dy_snapshot": args.dy_snapshot,
        "threads": args.threads, "keep_files": args.keep_files,
        "analyze_only": args.analyze_only,
    }
    os.makedirs(settings["datadir"], exist_ok=True)
    prescriptions = args.prescriptions.split(",")
    for p in prescriptions:
        if p not in PRESCRIPTIONS:
            parser.error(f"unknown prescription {p}")

    results = {}
    for p in prescriptions:
        results[p] = [run_one(p, seed, settings)
                      for seed in range(1, args.seeds + 1)]
    summary = analyze(results, settings)
    with open(os.path.join(settings["datadir"],
                           "running_coupling_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    plot_fig6(summary, os.path.join(settings["datadir"],
                                    "running_coupling_fig6.pdf"), args.dy_plot)
    plot_speed(summary, os.path.join(settings["datadir"],
                                     "running_coupling_speed.pdf"))

    for p, s in summary["prescriptions"].items():
        print(f"{p}: {s['configurations']} configurations, mean lambda "
              f"{s['lambda_mean']:.3f}")
        for y, q, e in zip(s["y"], s["Qs_over_Lambda"],
                           s["Qs_over_Lambda_error"]):
            print(f"  y = {y:5.2f}  Q_s/Lambda = {q:7.3f} +- {e:.3f}")
    failed = [k for k, v in summary["checks"].items()
              if isinstance(v, bool) and not v]
    for k, v in summary["checks"].items():
        print(f"{k}: {v}")
    if failed:
        print("FAILED: " + ", ".join(failed))
        sys.exit(1)
    print("All checks passed.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Validation of the forward light cone, the evolution, T^munu and the gluon
spectrum in their perturbative limit (issue #125): the dilute glasma.

With a constant color-charge density (useNucleus 0) and an IR regulator m
large compared with g^2 mu, both nuclei enter to second order, so the
energy density at tau -> 0+, tau eps at late times and dN/dy scale as
(g^2 mu)^4 / g^2, and at m << k << 1/a the gluon spectrum falls like
(g^2 mu)^4 ln(k/m) / k^4, the leading-order result (Kovner, McLerran,
Weigert; Gyulassy, McLerran). The deviation from the (g^2 mu)^4 scaling as
g^2 mu approaches m shows the onset of saturation.

The dilute range is bounded below by numerics: the field produced in the
collision is of second order, and for g^2 mu a and m a small enough (e.g.
g^2 mu = 0.01 GeV at a = 0.0125 fm) it falls below what the light-cone
solver resolves with its tolerance |F| < 1e-6 (eps(0+) comes out 13% low),
and the late energy density below the floor of the T^munu output
(1e-16 GeV/fm^3). The default scan therefore starts at 0.05 GeV.

The script runs IP-Glasma (the example input with useNucleus 0, no
JIMWLK) for a logarithmic scan of g2muGeV, all with the same seed (so the
random numbers are the same and only scaled), plus a run with g = 2, and
computes eps(tau -> 0+) (at the first time step) and tau eps at the final
time, averaged over the interior of the lattice (away from the open
boundaries of the evolution), dN/dy per area and the spectrum. It checks
that

  1. the exponents of eps(0+) and tau eps in g^2 mu are 4 within
     --exponent-tolerance (1%) over the dilute range,
  2. eps(0+) g^2 does not depend on g.

dN/dy and the spectrum (a fit A + B ln k of k^4 dN/d^2k in the window
[--k-min-over-m m, --k-max-times-a / a]) are only reported: they are
computed from the fields of the whole lattice, and the open boundaries of
the evolution turn the pure-gauge fields of the nuclei, which are of first
order in g^2 mu, into physical fields in a band of width ~tau near the
edges. In the dilute limit this band dominates the whole-lattice
quantities, which then scale like (g^2 mu)^2; the script reports its share
of the energy. (In collisions of nuclei the fields vanish near the edges.)

Usage (from validations/dilute_limit/):
    # quick test, 64^2 (about 10 s)
    python3 dilute_limit.py --ipglasma-path ../.. --quick
    # the setup of the paper: 256^2, L = 3.2 fm, m = 1 GeV
    python3 dilute_limit.py --ipglasma-path ../.. --max-workers 4

--ipglasma-path is a directory with the ipglasma executable,
qs2Adj_vs_Tp_vs_Y_240.in and tables/. The results go to
<datadir>/dilute_summary.json and <datadir>/dilute_limit.pdf.
"""

import argparse
import glob
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "utilities"))
from read_tmunu import COMPONENTS, read_tmunu  # noqa: E402

HBARC = 0.1973269804


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


def interior_mean(tmunu, meta, margin_fm):
    """The mean of T^tautau over the sites farther than margin_fm from the
    edges of the lattice."""
    t00 = tmunu[..., COMPONENTS.index("T00")]
    n = t00.shape[0]
    x = -0.5 * n * meta["dx_fm"] + meta["dx_fm"] * np.arange(n)
    inside = np.abs(x) < 0.5 * n * meta["dx_fm"] - margin_fm
    return float(t00[np.ix_(inside, inside)].mean())


def run_one(cfg, settings):
    tag = f"g2mu{cfg['g2mu']:g}_g{cfg['g']:g}"
    workdir = os.path.join(settings["datadir"], tag)
    cache = os.path.join(workdir, "result.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return cfg, json.load(f), None
    os.makedirs(workdir, exist_ok=True)
    tables = os.path.join(workdir, "tables")
    if os.path.isdir(settings["tables"]) and not os.path.exists(tables):
        os.symlink(settings["tables"], tables)
    write_input(settings["template"],
                dict(settings["overrides"], g2muGeV=cfg["g2mu"], g=cfg["g"]),
                os.path.join(workdir, "input"))
    env = dict(os.environ, OMP_NUM_THREADS=str(settings["threads"]))
    with open(os.path.join(workdir, "run.log"), "w") as log:
        result = subprocess.run([settings["ipglasma"], "input"], cwd=workdir,
                                stdout=log, stderr=subprocess.STDOUT, env=env)
    if result.returncode != 0:
        return cfg, None, f"ipglasma failed (exit {result.returncode})"
    out = {}
    files = sorted(glob.glob(os.path.join(workdir, "Tmunu-t*-0.ipgt")),
                   key=lambda p: float(p.split("Tmunu-t")[1].split("-")[0]))
    early, meta_early = read_tmunu(files[0])
    late, meta_late = read_tmunu(files[-1])
    margin = settings["margin_fm"]
    out["tau_early"] = meta_early["tau_fm"]
    out["tau_late"] = meta_late["tau_fm"]
    out["eps0"] = interior_mean(early, meta_early, margin)  # GeV/fm^3
    out["tau_eps"] = meta_late["tau_fm"] * interior_mean(late, meta_late,
                                                         margin)
    # the share of the energy outside the interior, where the open
    # boundaries turn the pure-gauge fields of the nuclei into physical ones
    t00 = late[..., COMPONENTS.index("T00")]
    n = t00.shape[0]
    inside_sites = (n - 2 * int(np.ceil(margin / meta_late["dx_fm"]))) ** 2
    out["edge_share"] = float(1. - interior_mean(late, meta_late, margin)
                              * inside_sites / t00.sum())
    with open(os.path.join(workdir, "gluonMultiplicity0.json")) as f:
        gluons = json.load(f)
    out["dN_dy_per_area"] = gluons["dN"] / settings["L"] ** 2  # 1/fm^2
    out["kt_GeV"] = gluons["kt_GeV"]
    out["dN_d2k"] = gluons["dN_d2k_GeV_minus2"]
    with open(cache, "w") as f:
        json.dump(out, f)
    for path in files:
        os.remove(path)
    return cfg, out, None


def power_fit(x, y):
    """The exponent p and its error of y = C x^p, fitted in log-log."""
    lx, ly = np.log(x), np.log(y)
    (p, c), cov = np.polyfit(lx, ly, 1, cov=True) if len(x) > 2 else \
        (np.polyfit(lx, ly, 1), np.zeros((2, 2)))
    return float(p), float(np.sqrt(cov[0, 0])), float(np.exp(c))


def main():
    parser = argparse.ArgumentParser(
        description="The dilute limit of the glasma (issue #125).")
    parser.add_argument("--ipglasma-path", default="../..")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--template", default=os.path.join(HERE, "..", "..",
                                                           "input"))
    parser.add_argument("--datadir", default="./dilute_output")
    parser.add_argument("--size", type=int, default=256)
    parser.add_argument("--L", type=float, default=3.2, help="[fm]")
    parser.add_argument("--m", type=float, default=1., help="[GeV]")
    parser.add_argument("--g2mu", type=float, nargs="+",
                        default=[0.05, 0.1, 0.2, 0.5, 1., 2.],
                        help="g^2 mu [GeV]")
    parser.add_argument("--dilute-max", type=float, default=0.2,
                        help="largest g^2 mu [GeV] of the dilute range")
    parser.add_argument("--maxtime", type=float, default=0.5,
                        help="the late time [fm/c]")
    parser.add_argument("--seed", type=int, default=5)
    parser.add_argument("--exponent-tolerance", type=float, default=0.04,
                        help="|p - 4| accepted (1%% of 4)")
    parser.add_argument("--k-min-over-m", type=float, default=2.)
    parser.add_argument("--k-max-times-a", type=float, default=0.3,
                        help="k_max a / hbar c of the spectrum window")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--quick", action="store_true",
                        help="64^2, L = 3.2 fm, maxtime 0.4, five g^2 mu "
                             "(the spectrum window is then too narrow)")
    args = parser.parse_args()
    if args.quick:
        args.size, args.maxtime = 64, 0.4
        args.g2mu = [0.01, 0.03, 0.1, 0.3, 1.]

    a = args.L / args.size
    path = os.path.abspath(args.ipglasma_path)
    settings = {
        "ipglasma": os.path.join(path, args.ipglasma_cmd),
        "template": os.path.abspath(args.template),
        "datadir": os.path.abspath(args.datadir),
        "tables": os.path.join(path, "tables"),
        "threads": args.threads, "L": args.L,
        # the open boundaries reach tau into the lattice
        "margin_fm": args.maxtime + 2 * a,
        "overrides": {
            "size": args.size, "L": args.L, "sizeOutput": args.size,
            "LOutput": args.L, "useNucleus": 0, "m": args.m,
            "seed": args.seed, "maxTime": args.maxtime, "useJIMWLK": 0,
            "useFluctuatingX": 0, "writeWilsonLines": 0, "writeTmunu": 1,
            "writeTmunuBinary": 1, "writeHydro": 0,
            # the first time step (an output time is rounded down to it)
            "outputTimes": f"{1.5 * a / 10:.6g}",
            "timeStepsPerLatticeSpacing": 10,
            "computeGluonMultiplicity": 1, "computeEccentricities": 0,
            "nucleusQsTableFileName": os.path.join(
                path, "qs2Adj_vs_Tp_vs_Y_240.in"),
        },
    }
    os.makedirs(settings["datadir"], exist_ok=True)
    dilute = [g for g in args.g2mu if g <= args.dilute_max]
    configs = [{"g2mu": g, "g": 1.} for g in args.g2mu]
    configs.append({"g2mu": dilute[len(dilute) // 2], "g": 2.})
    print(f"{len(configs)} runs of size {args.size}, L = {args.L} fm "
          f"(a = {a:.4f} fm, hbar c/a = {HBARC / a:.1f} GeV), m = {args.m} "
          f"GeV ...", flush=True)
    start, results, failures = time.time(), {}, []
    with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
        futures = [pool.submit(run_one, cfg, settings) for cfg in configs]
        for future in as_completed(futures):
            cfg, out, error = future.result()
            if error:
                failures.append(f"g2mu {cfg['g2mu']}, g {cfg['g']}: {error}")
            else:
                results[(cfg["g2mu"], cfg["g"])] = out
    print(f"done in {time.time() - start:.0f} s")
    for failure in failures:
        print("FAILED:", failure)

    summary = {"settings": vars(args), "a_fm": a, "runs": {}, "checks": {}}
    g2mu = np.array([g for g in args.g2mu if (g, 1.) in results])
    for g in g2mu:
        summary["runs"][f"{g:g}"] = results[(g, 1.)]
    names = ["eps0", "tau_eps", "dN_dy_per_area"]
    # dN/dy and the spectrum come from the fields of the whole lattice,
    # including the band near the open boundaries, where the pure-gauge
    # fields of the nuclei (first order in g^2 mu) become physical ones;
    # in the dilute limit this band dominates, so they are only reported
    checked = ["eps0", "tau_eps"]
    use = g2mu <= args.dilute_max
    print(f"\nScaling with g^2 mu (dilute range: g^2 mu <= "
          f"{args.dilute_max} GeV, {use.sum()} points):")
    fits = {}
    for name in names:
        y = np.array([results[(g, 1.)][name] for g in g2mu])
        p, dp, c = power_fit(g2mu[use], y[use])
        fits[name] = {"exponent": p, "error": dp, "coefficient": c,
                      "ratio_to_dilute": (y / (c * g2mu ** p)).tolist()}
        ok = abs(p - 4.) <= args.exponent_tolerance
        if name in checked:
            summary["checks"][f"exponent_{name}"] = bool(ok)
        label = ("OK  " if ok else "FAIL") if name in checked else "info"
        print(f"  {label} {name:15s} exponent "
              f"{p:.4f} +- {dp:.4f}; value / dilute fit: "
              + " ".join(f"{r:.3f}" for r in fits[name]["ratio_to_dilute"]))
    summary["fits"] = fits

    # k^4 dN/d^2k ~ A + B ln k in the window
    k_lo, k_hi = args.k_min_over_m * args.m, args.k_max_times_a * HBARC / a
    spectra = {}
    print(f"\nSpectrum, k^4 dN/d^2k = A + B ln k for {k_lo:.2f} < k < "
          f"{k_hi:.2f} GeV:")
    for g in g2mu[use]:
        out = results[(g, 1.)]
        k = np.array(out["kt_GeV"])
        dn = np.array(out["dN_d2k"])
        window = (k > k_lo) & (k < k_hi) & (dn > 0)
        if window.sum() < 4:
            print(f"  info g^2 mu {g:g}: fewer than 4 bins in the window")
            continue
        y = k[window] ** 4 * dn[window] / g ** 4
        B, A = np.polyfit(np.log(k[window]), y, 1)
        residual = float(np.sqrt(np.mean((y - A - B * np.log(k[window])) ** 2))
                         / np.mean(y))
        spectra[f"{g:g}"] = {"A": float(A), "B": float(B),
                             "relative_residual": residual,
                             "bins": int(window.sum())}
        print(f"  info g^2 mu {g:g} GeV: B/A = "
              f"{B / A:+.3f}, relative residual {residual:.3f} "
              f"({window.sum()} bins)")
    summary["spectra"] = spectra
    print("\nShare of the energy at the final time in the band near the open "
          "boundaries:\n  " + ", ".join(
              f"g^2 mu {g:g}: {results[(g, 1.)]['edge_share']:.3f}"
              for g in g2mu))
    print("  (dN/dy and the spectrum include this band, whose fields are of "
          "first order in g^2 mu, so they are not checked)")

    gc = configs[-1]["g2mu"]
    if (gc, 2.) in results and (gc, 1.) in results:
        ratio = 4. * results[(gc, 2.)]["eps0"] / results[(gc, 1.)]["eps0"]
        ok = abs(ratio - 1.) < 1e-3
        summary["checks"]["g_scaling"] = bool(ok)
        print(f"\n{'OK  ' if ok else 'FAIL'} eps(0+) g^2 at g = 2 over g = 1: "
              f"{ratio:.6f}")
    summary["checks"]["all_runs_succeeded"] = not failures
    passed = all(summary["checks"].values())
    summary["passed"] = bool(passed)
    with open(os.path.join(settings["datadir"], "dilute_summary.json"),
              "w") as f:
        json.dump(summary, f, indent=1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
    for name in names:
        left.plot(g2mu, fits[name]["ratio_to_dilute"], "o-", label=name)
    left.axhline(1., color="k", lw=0.5)
    left.set_xscale("log")
    left.set_xlabel(r"$g^2\mu$ [GeV]")
    left.set_ylabel(r"value / $C (g^2\mu)^4$ (fit to the dilute range)")
    left.legend(fontsize=8)
    for g in g2mu[use]:
        out = results[(g, 1.)]
        k = np.array(out["kt_GeV"])
        dn = np.array(out["dN_d2k"])
        good = dn > 0
        right.plot(k[good], k[good] ** 4 * dn[good] / g ** 4, "-",
                   label=fr"$g^2\mu = {g:g}$ GeV")
    right.axvspan(k_lo, k_hi, color="0.9")
    right.set_xscale("log")
    right.set_xlabel(r"$k$ [GeV]")
    right.set_ylabel(r"$k^4\,dN/dy\,d^2k\,/\,(g^2\mu)^4$")
    right.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(settings["datadir"], "dilute_limit.pdf"))
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{settings['datadir']}/dilute_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

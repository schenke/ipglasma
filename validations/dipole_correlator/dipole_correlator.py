#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Validation of the Wilson lines (issue #123): the dipole correlator of the
MV model against its exact lattice prediction.

IP-Glasma builds the Wilson lines of a constant color-charge density
(useNucleus 0) as V = prod_k exp(-i Lambda_k^a t^a) over Ny layers. Each
Lambda_k is Gaussian, with the two-point function (summed over the layers)

    G(r) = (g^2 mu a)^2 / N^2 sum_k K(k)^2 cos(k.r),
    K(k) = exp(-|k~| u) / [(k~ a)^2 + (m a)^2],  K(0) = 0 for m = 0,

on the periodic lattice (k~ the lattice momentum, u the UV damping), so
that the dipole correlator

    D(r) = Re < tr V^+(x) V(x + r) > / Nc

tends to D_G(r) = exp(-C_F [G(0) - G(r)]) for Ny -> infinity, with
corrections of order 1/Ny. No continuum or infinite-volume approximation
enters. This script runs IP-Glasma for a scan of g^2 mu, m, Ny (and
optionally UVDamp), computes D(r) from the binary Wilson lines along the
axes and the diagonal, averaged over all sites and over the two nuclei of
each run (two independent configurations), and checks that

  1. D approaches D_G for large Ny,
  2. Ny (D - D_G) tends to a constant,
  3. the extrapolation of D - D_G to Ny -> infinity vanishes within its
     statistical error,
  4. g drops out: a run with g = 2 gives the same D as with g = 1.

It also reports Q_s / g^2 mu, with D(r_s) = exp(-1/2) and Q_s^2 = 2/r_s^2
(T. Lappi, Eur. Phys. J. C 55 (2008) 285, arXiv:0711.3039). The results go
to <datadir>/dipole_summary.json, the figures to <datadir>/dipole_*.pdf.

Usage (from validations/dipole_correlator/):
    # a quick test on small lattices, about a minute
    python3 dipole_correlator.py --ipglasma-path ../.. --quick
    # the scan for the paper
    python3 dipole_correlator.py --ipglasma-path ../.. --max-workers 8

--ipglasma-path is a directory with the ipglasma executable and
qs2Adj_vs_Tp_vs_Y_240.in (e.g. the repository root after
./compile_IPGlasma.sh).
"""

import argparse
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
sys.path.insert(0, os.path.join(HERE, "..", "..", "utilities"))
from read_wilson_lines import read_wilson_lines  # noqa: E402

CF = 4. / 3.
HBARC = 0.1973269804  # GeV fm
NC = 3


# --------------------------------------------------------------------------
# physics
# --------------------------------------------------------------------------

def dipole_correlator(V):
    """D(r) = Re sum_x tr V^+(x) V(x + r) / (Nc N^2) for all lattice
    vectors r = (rx, ry), periodic; V has shape (N, N, 3, 3)."""
    n = V.shape[0]
    total = np.zeros((n, n))
    for a in range(NC):
        for b in range(NC):
            f = np.fft.fft2(V[:, :, a, b])
            # ifft2(|f_k|^2)[r] = sum_x conj(V(x)) V(x + r)
            total += np.fft.ifft2(np.conj(f) * f).real
    return total / (NC * n * n)


def lattice_two_point(n, g2mu_lat, m_lat, uv_lat):
    """G(r) of the summed layers, for all r, in the code's conventions
    (Init::computeWilsonLineMomentumKernel); lattice units."""
    k = 2. * np.pi * np.fft.fftfreq(n)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    kt2 = 4. * (np.sin(kx / 2.) ** 2 + np.sin(ky / 2.) ** 2)
    damping = np.exp(-np.sqrt(kt2) * uv_lat)
    with np.errstate(divide="ignore"):
        kernel = damping / (kt2 + m_lat * m_lat)
    if m_lat == 0.:
        kernel[0, 0] = 0.
    # ifft2 includes the 1/N^2
    return g2mu_lat ** 2 * np.fft.ifft2(kernel ** 2).real


def profiles(d2):
    """Axis (mean of +-x, +-y) and diagonal profiles of a 2D field, for
    r = 0 .. N/2 in lattice units."""
    n = d2.shape[0]
    r = np.arange(n // 2 + 1)
    axis = 0.25 * (d2[r, 0] + d2[-r, 0] + d2[0, r] + d2[0, -r])
    diag = 0.25 * (d2[r, r] + d2[-r, -r] + d2[r, -r] + d2[-r, r])
    return {"axis": axis, "diag": diag}


def crossing(r, d, level=np.exp(-0.5)):
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


def run_tag(cfg):
    return (f"g2mu{cfg['g2mu']:g}_m{cfg['m']:g}_uv{cfg['uvdamp']:g}"
            f"_ny{cfg['ny']}_g{cfg['g']:g}_seed{cfg['seed']}")


def run_one(cfg, settings):
    """Run IP-Glasma for one configuration and return the D(r) profiles of
    its two nuclei (cached in <workdir>/profiles.npz)."""
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
    overrides = {
        "size": settings["size"], "L": settings["L"], "Ny": cfg["ny"],
        "m": cfg["m"], "UVDamp": cfg["uvdamp"], "g2muGeV": cfg["g2mu"],
        "g": cfg["g"], "seed": cfg["seed"],
        "nucleusQsTableFileName": settings["qs_table"],
        "writeWilsonLines": 2, "runEvolution": 0, "useNucleus": 0,
    }
    write_input(settings["template"], overrides,
                os.path.join(workdir, "input"))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with open(os.path.join(workdir, "run.log"), "w") as log:
        result = subprocess.run([settings["ipglasma"], "input"], cwd=workdir,
                                stdout=log, stderr=subprocess.STDOUT, env=env)
    files = sorted(glob.glob(os.path.join(workdir, "WilsonLine_*")))
    if result.returncode != 0 or len(files) != 2:
        return cfg, None, (f"ipglasma failed (exit {result.returncode}, "
                           f"{len(files)} Wilson-line files), see {workdir}")
    out = {}
    for i, path in enumerate(files):
        V, meta = read_wilson_lines(path)
        prof = profiles(dipole_correlator(V))
        out[f"axis{i}"] = prof["axis"]
        out[f"diag{i}"] = prof["diag"]
        out["a_fm"] = np.array(meta["a_fm"])
        if not settings["keep_files"]:
            os.remove(path)
    np.savez(cache, **out)
    return cfg, out, None


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------

def weighted_fit(ny, dev, sigma):
    """Fit dev = d_inf + c / ny; returns d_inf, its error, c."""
    w = 1. / np.maximum(sigma, 1e-15) ** 2
    A = np.vstack([np.ones_like(ny), 1. / ny]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    d_inf, c = cov @ (A.T @ (w * dev))
    chi2 = float(np.sum(w * (dev - d_inf - c / ny) ** 2))
    dof = max(len(ny) - 2, 1)
    # inflate the error if the fit is worse than the statistical errors
    scale = max(1., np.sqrt(chi2 / dof))
    return float(d_inf), float(np.sqrt(cov[0, 0]) * scale), float(c)


def analyze_group(key, runs, settings):
    g2mu, m, uv = key
    n = settings["size"]
    a = settings["L"] / n
    g2mu_lat = g2mu * a / HBARC
    G = lattice_two_point(n, g2mu_lat, m * a / HBARC, uv / a * HBARC)
    exact = profiles(np.exp(-CF * (G[0, 0] - G)))
    r = np.arange(n // 2 + 1)
    # compare where D varies: up to the r where D_G has fallen to 0.05, at
    # most N/4 (away from the periodic images)
    small = np.nonzero(exact["axis"] < 0.05)[0]
    r_max = int(min(n // 4, small[0] if small.size else n // 4))
    r_star = int(np.argmin(np.abs(exact["axis"] - np.exp(-0.5))))
    rows = []
    for ny in sorted({cfg["ny"] for cfg, _ in runs}):
        samples = {"axis": [], "diag": []}
        for cfg, out in runs:
            if cfg["ny"] != ny:
                continue
            for i in (0, 1):
                samples["axis"].append(out[f"axis{i}"])
                samples["diag"].append(out[f"diag{i}"])
        row = {"ny": ny, "configurations": len(samples["axis"])}
        for direction in ("axis", "diag"):
            s = np.array(samples[direction])
            mean = s.mean(axis=0)
            err = s.std(axis=0, ddof=1) / np.sqrt(len(s)) if len(s) > 1 \
                else np.full_like(mean, np.nan)
            dev = mean - exact[direction]
            row[direction] = {"D": mean.tolist(), "error": err.tolist(),
                              "max_abs_dev": float(np.abs(dev[:r_max + 1]).max())}
        dev_star = row["axis"]["D"][r_star] - exact["axis"][r_star]
        row["dev_rstar"] = float(dev_star)
        row["error_rstar"] = float(row["axis"]["error"][r_star])
        row["ny_dev_rstar"] = float(ny * dev_star)
        r_s = crossing(r, np.array(row["axis"]["D"]))
        row["Qs_GeV"] = float(np.sqrt(2.) / (r_s * a) * HBARC)
        row["Qs_over_g2mu"] = row["Qs_GeV"] / g2mu
        rows.append(row)

    group = {"g2mu_GeV": g2mu, "m_GeV": m, "UVDamp_invGeV": uv,
             "a_fm": a, "r_max": r_max, "r_star": r_star,
             "exact": {d: exact[d].tolist() for d in exact},
             "Qs_exact_over_g2mu": float(
                 np.sqrt(2.) / (crossing(r, exact["axis"]) * a) * HBARC / g2mu),
             "rows": rows}

    # acceptance
    ny = np.array([row["ny"] for row in rows], dtype=float)
    max_dev = np.array([row["axis"]["max_abs_dev"] for row in rows])
    dev = np.array([row["dev_rstar"] for row in rows])
    err = np.array([row["error_rstar"] for row in rows])
    checks = {}
    if len(rows) >= 2:
        checks["approaches_exact"] = bool(max_dev[-1] < 0.2 * max_dev[0])
        # Ny (D - D_G) at the two largest Ny, within three combined errors
        a1, a2 = ny[-2] * dev[-2], ny[-1] * dev[-1]
        e = np.hypot(ny[-2] * err[-2], ny[-1] * err[-1])
        checks["ny_dev_constant"] = bool(abs(a1 - a2) < 3. * e + 0.05 * abs(a2))
    use = ny >= settings["ny_fit_min"]
    if use.sum() >= 3:
        d_inf, d_err, c = weighted_fit(ny[use], dev[use], err[use])
        group["fit"] = {"d_inf": d_inf, "error": d_err, "c": c,
                        "ny_min": settings["ny_fit_min"]}
        checks["extrapolation_vanishes"] = bool(abs(d_inf) < 3. * d_err)
    group["checks"] = checks
    return group


def plot_group(group, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    a = group["a_fm"]
    r_max = group["r_max"]
    rows = group["rows"]
    r = np.arange(len(group["exact"]["axis"]))
    show = rows if len(rows) <= 4 else [rows[0], rows[len(rows) // 3],
                                        rows[2 * len(rows) // 3], rows[-1]]
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(5.5, 7), sharex=True)
    top.plot(r[:r_max + 1] * a, group["exact"]["axis"][:r_max + 1], "k-",
             label=r"$D_\mathrm{G}$ (exact, $N_y\to\infty$)")
    for row in show:
        top.errorbar(r[:r_max + 1] * a, row["axis"]["D"][:r_max + 1],
                     yerr=row["axis"]["error"][:r_max + 1], fmt="o", ms=3,
                     label=fr"$N_y = {row['ny']}$")
        # the statistical error is scaled by Ny as well
        bottom.errorbar(
            r[:r_max + 1] * a,
            row["ny"] * (np.array(row["axis"]["D"][:r_max + 1])
                         - np.array(group["exact"]["axis"][:r_max + 1])),
            yerr=row["ny"] * np.array(row["axis"]["error"][:r_max + 1]),
            fmt="o-", ms=3, lw=0.8, label=fr"$N_y = {row['ny']}$")
    top.set_ylabel(r"$D(r)$")
    top.set_title(fr"$g^2\mu = {group['g2mu_GeV']:g}$ GeV, "
                  fr"$m = {group['m_GeV']:g}$ GeV, "
                  fr"UVDamp $= {group['UVDamp_invGeV']:g}$ GeV$^{{-1}}$",
                  fontsize=10)
    top.legend(fontsize=8)
    bottom.axhline(0., color="k", lw=0.5)
    bottom.set_xlabel(r"$r$ [fm]")
    bottom.set_ylabel(r"$N_y\,(D - D_\mathrm{G})$")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Dipole correlator of the MV model against its exact "
                    "lattice prediction (issue #123).")
    parser.add_argument("--ipglasma-path", default="../..",
                        help="directory with the ipglasma executable and "
                             "qs2Adj_vs_Tp_vs_Y_240.in")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--datadir", default="./dipole_output")
    parser.add_argument("--template", default=os.path.join(HERE, "input_dipole"))
    parser.add_argument("--size", type=int, default=512)
    parser.add_argument("--L", type=float, default=25.6, help="[fm]")
    parser.add_argument("--g2mu", type=float, nargs="+", default=[0.5, 1., 2.],
                        help="g^2 mu [GeV]")
    parser.add_argument("--ny", type=int, nargs="+",
                        default=[1, 2, 5, 10, 25, 50, 100])
    parser.add_argument("--m", type=float, nargs="+", default=[0., 0.2, 0.4],
                        help="IR regulator [GeV]")
    parser.add_argument("--uvdamp", type=float, nargs="*", default=[0.1],
                        help="UVDamp values [1/GeV] for an extra Ny scan at "
                             "the middle g^2 mu and m")
    parser.add_argument("--seeds", type=int, default=4,
                        help="runs per setting (two configurations each)")
    parser.add_argument("--first-seed", type=int, default=1)
    parser.add_argument("--ny-fit-min", type=int, default=5,
                        help="smallest Ny in the 1/Ny extrapolation")
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--quick", action="store_true",
                        help="small test: size 64, L 8, g2mu 1, m 0.2, "
                             "Ny 1..50, 2 seeds, no UVDamp scan")
    parser.add_argument("--analyze-only", action="store_true",
                        help="only analyze the results in --datadir")
    parser.add_argument("--keep-files", action="store_true",
                        help="keep the Wilson-line files")
    args = parser.parse_args()
    if args.quick:
        args.size, args.L, args.g2mu, args.m = 64, 8., [1.], [0.2]
        args.ny, args.seeds, args.uvdamp = [1, 2, 5, 10, 25, 50], 2, []

    ipglasma_path = os.path.abspath(args.ipglasma_path)
    settings = {
        "ipglasma": os.path.join(ipglasma_path, args.ipglasma_cmd),
        "qs_table": os.path.join(ipglasma_path, "qs2Adj_vs_Tp_vs_Y_240.in"),
        "template": os.path.abspath(args.template),
        "datadir": os.path.abspath(args.datadir),
        "size": args.size, "L": args.L, "ny_fit_min": args.ny_fit_min,
        "analyze_only": args.analyze_only, "keep_files": args.keep_files,
    }
    if not args.analyze_only:
        for path in (settings["ipglasma"], settings["qs_table"]):
            if not os.path.exists(path):
                sys.exit(f"not found: {path} (see --ipglasma-path)")
    os.makedirs(settings["datadir"], exist_ok=True)

    seeds = range(args.first_seed, args.first_seed + args.seeds)
    m_ref = 0.2 if 0.2 in args.m else args.m[0]
    g2mu_ref = sorted(args.g2mu)[len(args.g2mu) // 2]
    groups = [(g, m, 0.) for g, m in itertools.product(args.g2mu, args.m)]
    groups += [(g2mu_ref, m_ref, uv) for uv in args.uvdamp if uv > 0.]
    configs = [{"g2mu": g, "m": m, "uvdamp": uv, "ny": ny, "g": 1.,
                "seed": seed}
               for (g, m, uv), ny, seed in itertools.product(
                   groups, args.ny, seeds)]
    # g must drop out: the same run with g = 2
    g_check = {"g2mu": g2mu_ref, "m": m_ref, "uvdamp": 0., "ny": args.ny[-1],
               "g": 2., "seed": seeds[0]}
    configs.append(g_check)

    print(f"{len(configs)} runs of size {args.size}, L = {args.L} fm, "
          f"with up to {args.max_workers} workers ...", flush=True)
    start = time.time()
    results, failures = {}, []
    with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
        futures = [pool.submit(run_one, cfg, settings) for cfg in configs]
        for i, future in enumerate(as_completed(futures), 1):
            cfg, out, error = future.result()
            if error:
                failures.append(f"{run_tag(cfg)}: {error}")
            else:
                results[run_tag(cfg)] = (cfg, out)
            if i % max(1, len(configs) // 10) == 0 or i == len(configs):
                print(f"  {i}/{len(configs)} runs done "
                      f"({time.time() - start:.0f} s)", flush=True)
    for failure in failures:
        print("FAILED:", failure)

    summary = {"settings": {k: v for k, v in vars(args).items()},
               "groups": [], "failures": failures}
    for key in groups:
        runs = [(cfg, out) for cfg, out in results.values()
                if (cfg["g2mu"], cfg["m"], cfg["uvdamp"]) == key
                and cfg["g"] == 1.]
        if not runs:
            continue
        group = analyze_group(key, runs, settings)
        summary["groups"].append(group)
        tag = f"g2mu{key[0]:g}_m{key[1]:g}_uv{key[2]:g}"
        plot_group(group, os.path.join(settings["datadir"],
                                       f"dipole_{tag}.pdf"))

    reference = results.get(run_tag(dict(g_check, g=1.)))
    other = results.get(run_tag(g_check))
    if reference and other:
        diff = max(float(np.abs(reference[1][k] - other[1][k]).max())
                   for k in ("axis0", "axis1", "diag0", "diag1"))
        summary["g_check"] = {"max_abs_difference": diff,
                              "passed": bool(diff < 1e-10)}

    # report
    print("\nDipole correlator against its exact lattice prediction:")
    passed = not failures
    for group in summary["groups"]:
        print(f"\n g2mu = {group['g2mu_GeV']:g} GeV, m = {group['m_GeV']:g} "
              f"GeV, UVDamp = {group['UVDamp_invGeV']:g}/GeV "
              f"(a = {group['a_fm']:.4f} fm, r* = {group['r_star']} a, "
              f"Qs/g2mu exact = {group['Qs_exact_over_g2mu']:.4f})")
        print("    Ny  max|D-D_G|  D-D_G at r*    Ny(D-D_G)   Qs/g2mu")
        for row in group["rows"]:
            print(f"  {row['ny']:4d}  {row['axis']['max_abs_dev']:10.2e}  "
                  f"{row['dev_rstar']:+.2e}({row['error_rstar']:.0e})  "
                  f"{row['ny_dev_rstar']:+9.4f}  {row['Qs_over_g2mu']:.4f}")
        if "fit" in group:
            fit = group["fit"]
            print(f"  extrapolation (Ny >= {fit['ny_min']}): D - D_G -> "
                  f"{fit['d_inf']:+.2e} +- {fit['error']:.1e}, "
                  f"c = {fit['c']:+.3f}")
        for name, ok in group["checks"].items():
            print(f"  {'OK  ' if ok else 'FAIL'} {name}")
            passed &= ok
    if "g_check" in summary:
        g = summary["g_check"]
        print(f"\n {'OK  ' if g['passed'] else 'FAIL'} g drops out: "
              f"max |D(g=2) - D(g=1)| = {g['max_abs_difference']:.1e}")
        passed &= g["passed"]
    summary["passed"] = bool(passed)
    with open(os.path.join(settings["datadir"], "dipole_summary.json"),
              "w") as f:
        json.dump(summary, f, indent=1)
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{settings['datadir']}/dipole_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

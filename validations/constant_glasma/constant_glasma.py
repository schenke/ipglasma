#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Validation of the classical Yang-Mills part against the published results
of the same model (issue #126): the glasma of a constant color-charge
density, as in T. Lappi, Phys. Rev. C 67 (2003) 054903 [hep-ph/0303076].

SU(3), useNucleus 0, one layer (Ny 1) and no IR regulator (m 0), on the
lattices of Lappi's Fig. 3: N = 128, 256 and 512 sites at fixed g^2 mu L,
so g^2 mu a = g^2 mu L / N. The script runs IP-Glasma, writes T^munu at
dense output times and averages it over the whole lattice, which for
useNucleus 0 is periodic, as the lattices of the references. It computes

  * f_E = g^2 tau eps / (g^2 mu)^3 against g^2 mu tau, and its plateau,
    the value at the last time (g^2 mu tau = 20 by default). After its
    maximum at g^2 mu tau ~ 1, tau eps decreases slowly (d(tau eps)/dtau
    = -P_L, and P_L is slightly positive). Lappi's Table I gives the
    asymptote A of a fit A + B exp(-tau/tau0); his Fig. 4 reaches it at
    g^2 mu tau ~ 20 (tau = 2 fm there).
  * P_L/eps = tau^2 T^etaeta / T^tautau and P_T/eps = (T^xx + T^yy) /
    (2 T^tautau) against tau,
  * eps(tau -> 0+), at the first time step, against the lattice spacing,
  * f_N = g^2 (dN/dy) / (g^2 mu L)^2, the gluon number per unit rapidity
    and transverse area, from the multiplicity of the whole lattice at the
    last time. The runs compute it in transverse Coulomb gauge with the
    free lattice dispersion and the electric fields
    (gluonSpectrumEstimator 0), as Eqs. (27) and (28) of Lappi 2003, the
    definition of his Table I,

and checks that

  1. P_L/eps = -1 and P_T/eps = 1 within 0.02 at the first time step
     (purely longitudinal fields; the first step is at g^2 mu tau =
     g^2 mu a / 10, up to 0.094 here),
  2. P_L/eps is about 0 and P_T/eps about 1/2 at the last time,
  3. f_E reaches its plateau: it changes by less than --plateau-tolerance
     between half the last time and the last time,
  4. eps(0+) grows as a -> 0, while the plateau of f_E changes by less
     than --spacing-tolerance between the two finest lattices, beyond two
     combined standard deviations (T. Lappi, Phys. Lett. B 643 (2006) 11
     [hep-ph/0606207]),
  5. on the 512^2 lattices of Lappi's Table I (g^2 mu L = 72, 120, 192),
     the plateau of f_E and f_N agree with his values within three
     combined standard deviations (his quoted error and the error of the
     mean over the seeds).

The f_E of Krasnitz, Nara and Venugopalan, Phys. Rev. Lett. 87 (2001)
192302 [hep-ph/0108092], as quoted by Lappi and halved as in their
erratum, Nucl. Phys. A 727 (2003) 427 [hep-ph/0305112], are printed for
information; their lattice spacing is not quoted, so they are not checked.

Usage (from validations/constant_glasma/):
    python3 constant_glasma.py --ipglasma-path ../.. --quick
    python3 constant_glasma.py --ipglasma-path ../.. --max-workers 2
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
from read_tmunu import COMPONENTS, read_tmunu  # noqa: E402

HBARC = 0.1973269804
STEPS_PER_SPACING = 10

# Lappi 2003, Table I: f_E and f_N (SU(3), f_N with the free dispersion)
# on the 512^2 lattices nearest to the continuum limit in his Fig. 3, by
# g^2 mu L
LAPPI_SIZE = 512
LAPPI_FE = {72.: (0.265, 0.005), 120.: (0.227, 0.003), 192.: (0.238, 0.005)}
LAPPI_FN = {72.: (0.297, 0.006), 120.: (0.289, 0.003), 192.: (0.329, 0.006)}
# Krasnitz, Nara, Venugopalan 2001, as quoted in Lappi 2003 (Sec. V) for
# Lambda_s R_A = 25 and 83.7 with pi R_A^2 = L^2, halved per their erratum
KNV_FE = {np.sqrt(np.pi) * 25.: 0.537 / 2, np.sqrt(np.pi) * 83.7: 0.497 / 2}


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


def lattice_mean(tmunu):
    """Means of T^tautau, T^xx + T^yy and tau^2 T^etaeta over the whole
    (periodic) lattice."""
    t = tmunu.astype(np.float64)
    c = {name: COMPONENTS.index(name)
         for name in ("T00", "Txx", "Tyy", "tau2_Tetaeta")}
    return (float(t[..., c["T00"]].mean()),
            float((t[..., c["Txx"]] + t[..., c["Tyy"]]).mean()),
            float(t[..., c["tau2_Tetaeta"]].mean()))


def run_one(cfg, settings):
    tag = f"g2muL{cfg['g2muL']:g}_N{cfg['size']}_seed{cfg['seed']}"
    workdir = os.path.join(settings["datadir"], tag)
    cache = os.path.join(workdir, "result.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return cfg, json.load(f), None
    os.makedirs(workdir, exist_ok=True)
    g2mu = settings["g2mu"]
    L = cfg["g2muL"] * HBARC / g2mu
    size = cfg["size"]
    a = L / size
    dt = a / STEPS_PER_SPACING
    taus = [x * HBARC / g2mu for x in settings["g2mu_tau"]]
    taus = [t for t in taus if 2 * dt <= t]
    max_time = taus[-1]
    # the first time step (an output time is rounded down to it) for
    # eps(0+)
    outputs = [1.5 * dt] + taus[:-1]
    overrides = dict(settings["overrides"], size=size, L=f"{L:.10g}",
                     sizeOutput=size, LOutput=f"{L:.10g}", seed=cfg["seed"],
                     maxTime=f"{max_time:.10g}",
                     outputTimes=",".join(f"{t:.10g}" for t in outputs))
    write_input(settings["template"], overrides,
                os.path.join(workdir, "input"))
    env = dict(os.environ, OMP_NUM_THREADS=str(settings["threads"]))
    with open(os.path.join(workdir, "run.log"), "w") as log:
        result = subprocess.run([settings["ipglasma"], "input"], cwd=workdir,
                                stdout=log, stderr=subprocess.STDOUT, env=env)
    if result.returncode != 0:
        return cfg, None, f"ipglasma failed (exit {result.returncode})"
    rows = []
    for path in glob.glob(os.path.join(workdir, "Tmunu-t*-0.ipgt")):
        tmunu, meta = read_tmunu(path)
        tau = meta["tau_fm"]
        eps, pt2, pl = lattice_mean(tmunu)
        os.remove(path)
        rows.append({"tau_fm": tau, "g2mu_tau": tau * g2mu / HBARC,
                     # f_E = g^2 tau eps / (g^2 mu)^3, tau eps in GeV^3
                     "f_E": settings["g"] ** 2 * tau * eps * HBARC ** 2
                     / g2mu ** 3,
                     # g^2 eps / (g^2 mu)^4, eps in GeV^4
                     "eps_scaled": settings["g"] ** 2 * eps * HBARC ** 3
                     / g2mu ** 4,
                     "PL_over_eps": pl / eps, "PT_over_eps": 0.5 * pt2 / eps})
    rows.sort(key=lambda r: r["tau_fm"])
    with open(os.path.join(workdir, "gluonMultiplicity0.json")) as f:
        dN = json.load(f)["dN"]
    out = {"L_fm": L, "a_fm": a, "rows": rows,
           # f_N = g^2 (dN/dy) / (g^2 mu)^2 / L^2
           "f_N": settings["g"] ** 2 * dN / cfg["g2muL"] ** 2}
    with open(cache, "w") as f:
        json.dump(out, f)
    return cfg, out, None


def mean_error(values):
    v = np.array(values, dtype=float)
    error = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.
    return float(v.mean()), float(error)


def main():
    parser = argparse.ArgumentParser(
        description="Constant-g^2 mu glasma against published classical "
                    "Yang-Mills results (issue #126).")
    parser.add_argument("--ipglasma-path",
                        default=os.path.join(HERE, "..", ".."),
                        help="directory with the ipglasma executable and "
                             "qs2Adj_vs_Tp_vs_Y_240.in (default: the "
                             "repository; build/src is tried as well)")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--template", default=os.path.join(HERE, "..", "..",
                                                           "input"))
    parser.add_argument("--datadir", default="./constant_glasma_output")
    parser.add_argument("--g2mu", type=float, default=1., help="[GeV]")
    parser.add_argument("--g", type=float, default=1.)
    parser.add_argument("--g2muL", type=float, nargs="+", default=[72., 120.])
    parser.add_argument("--sizes", type=int, nargs="+",
                        default=[128, 256, 512],
                        help="lattice sizes N, g^2 mu a = g^2 mu L / N")
    parser.add_argument("--g2mu-tau", type=float, nargs="+",
                        default=[0.1, 0.2, 0.3, 0.5, 0.75, 1., 1.5, 2., 2.5,
                                 3., 4., 5., 6., 8., 10., 12., 14., 16., 18.,
                                 20.],
                        help="output times g^2 mu tau; the last is the "
                             "final time, where f_N is measured")
    parser.add_argument("--seeds", type=int, default=4)
    parser.add_argument("--plateau-tolerance", type=float, default=0.1)
    parser.add_argument("--spacing-tolerance", type=float, default=0.1)
    parser.add_argument("--pressure-tolerance", type=float, default=0.1,
                        help="for P_L/eps -> 0 and P_T/eps -> 1/2 at the "
                             "last time")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--quick", action="store_true",
                        help="g^2 mu L = 36 on 128^2 and 256^2, four "
                             "seeds")
    parser.add_argument("--analyze-only", action="store_true",
                        help="only analyze the cached results")
    args = parser.parse_args()
    if args.quick:
        args.g2muL, args.sizes = [36.], [128, 256]

    path = os.path.abspath(args.ipglasma_path)
    executable = os.path.join(path, args.ipglasma_cmd)
    if not args.analyze_only and not os.path.exists(executable):
        # a CMake build without `make install` leaves it in build/src
        built = os.path.join(path, "build", "src", args.ipglasma_cmd)
        if not os.path.exists(built):
            sys.exit(f"no {args.ipglasma_cmd} in {path} or {path}/build/src;"
                     " build IP-Glasma or set --ipglasma-path")
        executable = built
    settings = {
        "ipglasma": executable,
        "template": os.path.abspath(args.template),
        "datadir": os.path.abspath(args.datadir),
        "threads": args.threads, "g2mu": args.g2mu, "g": args.g,
        "g2mu_tau": sorted(args.g2mu_tau),
        "overrides": {
            "useNucleus": 0, "g2muGeV": args.g2mu, "g": args.g, "Ny": 1,
            "m": 0, "UVDamp": 0, "useJIMWLK": 0, "useFluctuatingX": 0,
            "runEvolution": 1, "inverseQsForMaxTime": 0,
            "timeStepsPerLatticeSpacing": STEPS_PER_SPACING,
            "writeWilsonLines": 0, "writeWilsonLineGeometry": 0,
            "writeTmunu": 1, "writeTmunuBinary": 1, "writeHydro": 0,
            "computeGluonMultiplicity": 1, "computeEccentricities": 0,
            # f_N of Lappi counts gluons with the electric fields
            "gluonSpectrumEstimator": 0,
            "nucleusQsTableFileName": os.path.join(
                path, "qs2Adj_vs_Tp_vs_Y_240.in"),
        },
    }
    os.makedirs(settings["datadir"], exist_ok=True)
    configs = [{"g2muL": L, "size": n, "seed": s} for L, n, s in
               itertools.product(args.g2muL, args.sizes,
                                 range(1, args.seeds + 1))]
    results, failures = [], []
    if args.analyze_only:
        for cfg in configs:
            cache = os.path.join(
                settings["datadir"],
                f"g2muL{cfg['g2muL']:g}_N{cfg['size']}_seed{cfg['seed']}",
                "result.json")
            if os.path.exists(cache):
                with open(cache) as f:
                    results.append((cfg, json.load(f)))
            else:
                failures.append(f"{cfg}: no cached result")
    else:
        print(f"{len(configs)} runs, g^2 mu = {args.g2mu} GeV ...",
              flush=True)
        start = time.time()
        with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
            futures = [pool.submit(run_one, cfg, settings)
                       for cfg in configs]
            for future in as_completed(futures):
                cfg, out, error = future.result()
                if error:
                    failures.append(f"{cfg}: {error}")
                else:
                    results.append((cfg, out))
        print(f"done in {time.time() - start:.0f} s")
    for failure in failures:
        print("FAILED:", failure)

    summary = {"settings": vars(args), "groups": [], "checks": {}}
    checks = summary["checks"]
    for L_, n in itertools.product(args.g2muL, args.sizes):
        runs = [out for cfg, out in results
                if cfg["g2muL"] == L_ and cfg["size"] == n]
        if not runs:
            continue
        # the seeds share the output steps; average them time by time
        rows = []
        for i, first in enumerate(runs[0]["rows"]):
            row = {"g2mu_tau": first["g2mu_tau"]}
            for key in ("f_E", "eps_scaled", "PL_over_eps", "PT_over_eps"):
                row[key], row[key + "_error"] = mean_error(
                    [out["rows"][i][key] for out in runs])
            rows.append(row)
        tag = f"g2muL {L_:g}, N {n}"
        first, last = rows[0], rows[-1]
        middle = min(rows, key=lambda r: abs(r["g2mu_tau"]
                                             - 0.5 * last["g2mu_tau"]))
        group = {"g2muL": L_, "size": n, "g2mua": L_ / n, "seeds": len(runs),
                 "L_fm": runs[0]["L_fm"], "rows": rows,
                 "eps0": first["eps_scaled"],
                 "eps0_error": first["eps_scaled_error"],
                 "fE": last["f_E"], "fE_error": last["f_E_error"]}
        group["fN"], group["fN_error"] = mean_error(
            [out["f_N"] for out in runs])
        checks[f"{tag}: P_L/eps = -1, P_T/eps = 1 at the first step"] = bool(
            abs(first["PL_over_eps"] + 1.) < 0.02
            and abs(first["PT_over_eps"] - 1.) < 0.02)
        checks[f"{tag}: P_L/eps -> 0, P_T/eps -> 1/2"] = bool(
            abs(last["PL_over_eps"]) < args.pressure_tolerance
            and abs(last["PT_over_eps"] - 0.5) < args.pressure_tolerance)
        checks[f"{tag}: f_E reaches its plateau"] = bool(
            abs(last["f_E"] / middle["f_E"] - 1.) < args.plateau_tolerance)
        reference = LAPPI_FE.get(L_) if n == LAPPI_SIZE else None
        if reference is not None:
            group["lappi_fE"] = reference
            sigma = np.hypot(group["fE_error"], reference[1])
            checks[f"{tag}: f_E = {reference[0]} +- {reference[1]} "
                   "(Lappi 2003, Table I)"] = bool(
                abs(group["fE"] - reference[0]) <= 3 * sigma)
            reference_fN = LAPPI_FN[L_]
            group["lappi_fN"] = reference_fN
            sigma = np.hypot(group["fN_error"], reference_fN[1])
            checks[f"{tag}: f_N = {reference_fN[0]} +- {reference_fN[1]} "
                   "(Lappi 2003, Table I)"] = bool(
                abs(group["fN"] - reference_fN[0]) <= 3 * sigma)
        summary["groups"].append(group)
        print(f"\n{tag} (g2mu a {L_ / n:.3g}, L = {runs[0]['L_fm']:.3f} fm, "
              f"{len(runs)} seeds):")
        print("   g2mu tau   f_E               P_L/eps   P_T/eps")
        for row in rows:
            print(f"   {row['g2mu_tau']:8.4g}   {row['f_E']:.4f} +- "
                  f"{row['f_E_error']:.4f}   {row['PL_over_eps']:+.3f}"
                  f"    {row['PT_over_eps']:.3f}")
        print(f"   f_N at g2mu tau {last['g2mu_tau']:.4g}: "
              f"{group['fN']:.4f} +- {group['fN_error']:.4f}")
        if reference is not None:
            print(f"   Lappi 2003: f_E = {reference[0]} +- {reference[1]}, "
                  f"f_N = {LAPPI_FN[L_][0]} +- {LAPPI_FN[L_][1]}")

    # the lattice-spacing dependence at fixed g^2 mu L
    for L_ in args.g2muL:
        groups = sorted([g for g in summary["groups"] if g["g2muL"] == L_],
                        key=lambda g: g["size"])
        if len(groups) < 2:
            continue
        eps0 = [g["eps0"] for g in groups]
        checks[f"g2muL {L_:g}: eps(0+) grows as a -> 0"] = bool(
            all(b > a for a, b in zip(eps0, eps0[1:])))
        fine, coarse = groups[-1], groups[-2]
        sigma = np.hypot(fine["fE_error"], coarse["fE_error"])
        checks[f"g2muL {L_:g}: f_E plateau changes by less than "
               f"{args.spacing_tolerance:g} between the finest lattices"] = \
            bool(abs(fine["fE"] - coarse["fE"]) - 2 * sigma
                 < args.spacing_tolerance * fine["fE"])
        print(f"\ng2muL {L_:g}:  g2mu a    g^2 eps(0+)/(g^2 mu)^4   f_E "
              "plateau         f_N")
        for g in groups:
            print(f"             {g['g2mua']:.4f}   {g['eps0']:8.4f}"
                  f"                 {g['fE']:.4f} +- {g['fE_error']:.4f}"
                  f"   {g['fN']:.4f} +- {g['fN_error']:.4f}")
    print("\nKrasnitz, Nara, Venugopalan (2001), halved per their erratum "
          "(lattice spacing not quoted, not checked): "
          + ", ".join(f"f_E = {v:.3f} at g2muL {k:.1f}"
                      for k, v in KNV_FE.items()))
    checks["all runs succeeded"] = not failures
    print()
    for name, ok in checks.items():
        print(f"{'OK  ' if ok else 'FAIL'} {name}")
    passed = all(checks.values())
    summary["passed"] = bool(passed)
    with open(os.path.join(settings["datadir"],
                           "constant_glasma_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (left, middle, right) = plt.subplots(1, 3, figsize=(14, 4))
    for i, g in enumerate(summary["groups"]):
        color = f"C{i}"
        x = np.array([r["g2mu_tau"] for r in g["rows"]])
        label = fr"$g^2\mu L = {g['g2muL']:g}$, $N = {g['size']}$"
        left.errorbar(x, [r["f_E"] for r in g["rows"]],
                      yerr=[r["f_E_error"] for r in g["rows"]], fmt="o-",
                      ms=3, color=color, label=label)
        if "lappi_fE" in g:
            value, error = g["lappi_fE"]
            left.errorbar([x[-1] * 1.15], [value], yerr=[error], fmt="D",
                          ms=5, color=color, mfc="none", capsize=3)
        middle.plot(x, [r["PL_over_eps"] for r in g["rows"]], "o-", ms=3,
                    color=color)
        middle.plot(x, [r["PT_over_eps"] for r in g["rows"]], "s--", ms=3,
                    color=color)
    left.set_ylabel(r"$f_E$ (open diamonds: Lappi 2003)")
    middle.set_ylabel(r"$P_L/\varepsilon$ (o), $P_T/\varepsilon$ (s)")
    for ax in (left, middle):
        ax.set_xscale("log")
        ax.set_xlabel(r"$g^2\mu\tau$")
    left.legend(fontsize=7)
    twin = right.twinx()
    for i, L_ in enumerate(args.g2muL):
        groups = [g for g in summary["groups"] if g["g2muL"] == L_]
        x = [g["g2mua"] for g in groups]
        right.errorbar(x, [g["fE"] for g in groups],
                       yerr=[g["fE_error"] for g in groups], fmt="o-",
                       ms=4, color=f"C{i}", label=fr"$g^2\mu L = {L_:g}$")
        right.errorbar(x, [g["fN"] for g in groups],
                       yerr=[g["fN_error"] for g in groups], fmt="D:",
                       ms=4, color=f"C{i}")
        twin.plot(x, [g["eps0"] for g in groups], "s--", ms=4,
                  color=f"C{i}")
    right.set_xlabel(r"$g^2\mu a$")
    right.set_ylabel(r"plateau of $f_E$ (o), $f_N$ (D)")
    twin.set_ylabel(r"$g^2\varepsilon(0^+)/(g^2\mu)^4$ (s)")
    right.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(settings["datadir"], "constant_glasma.pdf"))
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{settings['datadir']}/constant_glasma_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

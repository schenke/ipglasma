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

The script runs the example input without JIMWLK, with b fixed (0 by
default) and L = 22.5 fm (which suffices at b = 0; use 30 fm for other
impact parameters), for three scans:

  * size (a = L/size, by default 0.234 to 0.018 fm) at the Ny and m of the
    input,
  * Ny at a fixed size (768, a = 0.029 fm),
  * m at a fixed size,

each with several color-charge seeds, and computes dE/deta = tau int T^tautau
d^2x at tau = 0.2 fm and at the final time (0.4 fm), dN/dy, eps_2, eps_3
and their angles. It checks that the geometry is the same in all runs
(N_part and N_coll lists, b; <Q_s> of the IP-Sat input is reported) and that different
color-charge seeds give different events, and reports the mean and its
statistical error of each observable against a, Ny and m. The physical
Q_s of the Wilson lines at fixed QsMuRatio, which changes with Ny and m,
is measured by the dipole setup (validations/dipole_correlator, #123).

The lattice converges once it resolves the color-charge scale, g^2 mu a
well below 1. A run that only builds the Wilson lines gives the g^2 mu =
Q_s/QsMuRatio of each site of the geometry, and the script reports its
median and maximum (the median weighted with (g^2 mu)^2, i.e. with the
multiplicity) and g^2 mu a for each size. For seed 3 the hot spots of the
example input (BGq = 0.3 GeV^-2, a width of 0.11 fm) overlap to g^2 mu of
up to 13 GeV, with a median of 5.2 GeV (1/g^2 mu = 0.038 fm), so that only
the smallest lattice spacings resolve most of the event. As a reference,
a fourth scan runs a constant g^2 mu (useNucleus 0, L = 3 fm, on its
periodic lattice) of about this median at the Ny and m of the input,
against g^2 mu a: dN/dy is within about 4% of its limit for g^2 mu a below
about 0.6, while dE/deta still grows by 5% from g^2 mu a = 0.3 to 0.15.
The energy sits in modes up to several g^2 mu, where the gluon spectrum
falls like 1/k^4, and converges only like a, with a deficit of about
0.3 g^2 mu a (the energy density at tau = 0 even diverges, as ln^2 a;
T. Lappi, Phys. Lett. B 643 (2006) 11). An independent classical
Yang-Mills implementation gives the same dependence on a for the same
Wilson lines. With JIMWLK the tail of the spectrum is less steep (the
anomalous dimension of the evolved Wilson lines), so the energy converges
even more slowly.

The script checks the convergence in the lattice spacing for each
observable: the relative change between the two finest lattices must be
below --convergence-tolerance (1%) plus two standard deviations. For an
observable that does not converge, it prints a prominent warning, with an
estimate of the remaining deviation from an extrapolation linear in a,
marks it in the figure and in the summary. With the defaults, dN/dy and
the eccentricities converge, while dE/deta still changes by about 4% between
the two finest lattices; the linear extrapolation puts its continuum value
about 16% above the finest one, which would need lattice spacings far below
0.018 fm to reach.

Where the energy sits in k_T: for each lattice of the size scans, the
script reads the gluon spectrum (gluonMultiplicity0.json) and reports the
energy below and above 2 and 4 times the g^2 mu scale (the binned energy
below, the total dE minus it above), an estimate of the continuum energy
(the binned energy below k_c plus a tail k^4 dN/d^2k = C beyond it,
2 pi C/k_c, for k_c = 0.3 and 0.4 times the largest k_T; meaningful once k_c
lies in the 1/k^4 tail) and the exponent p of k^4 dN/d^2k ~ k^p between 6
and 12 GeV (--tail-fit); the mean spectra go to the summary and to
numerical_parameters_spectra.pdf.

With --jimwlk, all runs evolve the Wilson lines with the JIMWLK parameters
of the template (the example input: from x = 0.01 to 1e-4), and the g^2 mu
scale is that of the initial Wilson lines; --x sets projectileX and targetX
(the x of the Wilson lines, or with --jimwlk the x the evolution ends at).
The runs, the summary and the figures of these variants are tagged
(jimwlk_, x0.01_, ...), so they never mix with the default ones. --ny, --m
and --uniform-sizes without values skip those scans.

Usage (from validations/numerical_parameters/):
    # quick test: sizes 64 and 128, Ny and m scans at 64, the constant
    # g^2 mu at sizes 32 and 64, 2 seeds
    python3 numerical_parameters.py --ipglasma-path ../.. --quick
    # the setup of the paper (about 25 min with 2 workers of 8 threads)
    python3 numerical_parameters.py --ipglasma-path ../.. --threads 8 \\
        --max-workers 2
    # the size scan with JIMWLK (the example input from x = 0.01 to 1e-4)
    python3 numerical_parameters.py --ipglasma-path ../.. --jimwlk \\
        --sizes 768 1024 1280 1536 --ny --m --uniform-sizes --color-seeds 2

--ipglasma-path is a directory with the ipglasma executable,
qs2Adj_vs_Tp_vs_Y_240.in, nucleusConfigurations/ and tables/. The results
go to <datadir>/numerical_parameters[_<variant>]_summary.json, the figures
to <datadir>/numerical_parameters[_<variant>]*.pdf.
"""

import argparse
import glob
import json
import os
import re
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
# the constant g^2 mu has no geometry, so no eccentricities
UNIFORM_OBSERVABLES = OBSERVABLES[:3]
HBARC = 0.1973269804  # GeV fm
# the energy is split at these k_T, in units of the g^2 mu scale; the tail
# estimate assumes k^4 dN/d^2k constant above these fractions of the largest
# k_T of the spectrum
ENERGY_CUTS = (2., 4.)
TAIL_FRACTIONS = (0.3, 0.4)


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
    # e.g. jimwlk_ for the runs with --jimwlk, so that they never reuse the
    # cached results of runs without it
    variant = f"{cfg['variant']}_" if cfg.get("variant") else ""
    uniform = (f"uniform_g2mu{cfg['g2mu']:g}_L{cfg['L']:g}_"
               if "g2mu" in cfg else "")
    return (f"{variant}{uniform}size{cfg['size']}_ny{cfg['ny']}"
            f"_m{cfg['m']:g}_cc{cfg['color_seed']}")


def read_run(workdir, with_geometry=True):
    """The observables and (with_geometry) the geometry of one run."""
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
    if not with_geometry:
        return {k: float(v) for k, v in out.items()}, None
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


def read_spectrum(workdir):
    """The gluon spectrum of gluonMultiplicity0.json (OUTPUT.md): the k_T
    bins, dN/d^2k and dE/d^2k, and the totals dN and dE, which include the
    modes above the largest bin; None without the file."""
    path = os.path.join(workdir, "gluonMultiplicity0.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        data = json.load(f)
    return {"kt": data["kt_GeV"], "dN_d2k": data["dN_d2k_GeV_minus2"],
            "dE_d2k": data["dE_d2k_GeV_minus1"], "dN": data["dN"],
            "dE": data["dE_GeV"]}


def run_ipglasma(settings, overrides, workdir):
    """Runs ipglasma in workdir with the template and the overrides; the
    error message if it failed, else None."""
    os.makedirs(workdir, exist_ok=True)
    tables = os.path.join(workdir, "tables")
    if os.path.isdir(settings["tables"]) and not os.path.exists(tables):
        os.symlink(settings["tables"], tables)
    write_input(settings["template"], overrides,
                os.path.join(workdir, "input"))
    env = dict(os.environ, OMP_NUM_THREADS=str(settings["threads"]))
    with open(os.path.join(workdir, "run.log"), "w") as log:
        result = subprocess.run([settings["ipglasma"], "input"], cwd=workdir,
                                stdout=log, stderr=subprocess.STDOUT, env=env)
    if result.returncode != 0:
        return f"ipglasma failed (exit {result.returncode})"
    return None


def read_g2mu(path):
    """g^2 mu [GeV] of each site of a WilsonLineGeometry file (OUTPUT.md)."""
    with open(path, "rb") as f:
        data = f.read()
    if data[:8] != b"IPGGEO1\0":
        raise ValueError(f"{path} is not a Wilson-line geometry file")
    length = int.from_bytes(data[8:16], "little")
    meta = json.loads(data[16:16 + length])
    sites = meta["N"] * meta["N"]
    # the nucleons (x, y, z, proton) come first, then g^2 mu^2 of each site,
    # (g^2 mu a)^2 / g^2 in lattice units
    g2mu2 = np.frombuffer(data, "<f8", sites,
                          16 + length + meta["nucleons"] * 4 * 8)
    return meta["g"] * np.sqrt(g2mu2) * HBARC / (meta["L_fm"] / meta["N"])


def g2mu_scale(settings, size):
    """The scale g^2 mu = Q_s/QsMuRatio of the geometry, from a run that
    only builds the Wilson lines (one layer) and writes the geometry files:
    the median weighted with (g^2 mu)^2 and the maximum, over the sites of
    both nuclei."""
    variant = f"{settings['variant']}_" if settings["variant"] else ""
    workdir = os.path.join(settings["datadir"], f"{variant}g2mu_size{size}")
    cache = os.path.join(workdir, "g2mu.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return json.load(f)
    # the initial Wilson lines, without JIMWLK (with --jimwlk at its
    # starting x)
    overrides = dict(settings["overrides"], size=size, sizeOutput=size,
                     Ny=1, runEvolution=0, writeWilsonLines=2,
                     writeWilsonLineGeometry=1, colorChargeSeed=1,
                     useJIMWLK=0, **settings["initial_x"])
    error = run_ipglasma(settings, overrides, workdir)
    if error:
        print(f"FAILED: the g^2 mu run in {workdir}: {error}")
        return None
    g2mu = np.concatenate([read_g2mu(path) for path in glob.glob(
        os.path.join(workdir, "WilsonLineGeometry_*"))])
    order = np.argsort(g2mu)
    weight = np.cumsum(g2mu[order] ** 2)
    scale = {"median": float(g2mu[order][np.searchsorted(
                 weight, 0.5 * weight[-1])]),
             "max": float(g2mu.max())}
    with open(cache, "w") as f:
        json.dump(scale, f)
    for path in glob.glob(os.path.join(workdir, "WilsonLine*")):
        os.remove(path)
    return scale


def run_one(cfg, settings):
    workdir = os.path.join(settings["datadir"], run_tag(cfg))
    cache = os.path.join(workdir, "result.json")
    if os.path.exists(cache):
        with open(cache) as f:
            data = json.load(f)
        return (cfg, data["observables"], data["geometry"],
                data.get("spectrum"), None)
    overrides = dict(settings["overrides"], size=cfg["size"],
                     sizeOutput=cfg["size"], Ny=cfg["ny"], m=cfg["m"],
                     colorChargeSeed=cfg["color_seed"])
    uniform = "g2mu" in cfg
    if uniform:
        overrides.update(useNucleus=0, g2muGeV=cfg["g2mu"], L=cfg["L"],
                         LOutput=cfg["L"], computeEccentricities=0)
    error = run_ipglasma(settings, overrides, workdir)
    if error:
        return cfg, None, None, None, error
    observables, geometry = read_run(workdir, with_geometry=not uniform)
    spectrum = read_spectrum(workdir)
    with open(cache, "w") as f:
        json.dump({"observables": observables, "geometry": geometry,
                   "spectrum": spectrum}, f)
    for path in glob.glob(os.path.join(workdir, "*.ipgt")):
        os.remove(path)
    return cfg, observables, geometry, spectrum, None


def statistics(values):
    values = np.array(values, dtype=float)
    err = values.std(ddof=1) / np.sqrt(len(values)) if len(values) > 1 \
        else float("nan")
    return {"mean": float(values.mean()), "error": float(err),
            "values": values.tolist()}


# observables whose convergence in the lattice spacing is checked (not the
# angles)
CONVERGENCE_OBSERVABLES = ["dE/deta(0.2)", "dE/deta(final)", "dN/dy",
                           "eps_2", "eps_3"]


def convergence(rows, L, tolerance):
    """Per observable of a size scan: the relative change between the two
    finest lattices with its error, whether it is below tolerance plus two
    standard deviations, and, from an extrapolation linear in a, the
    relative deviation of the finest lattice from a = 0."""
    if len(rows) < 2:
        return {}
    second, finest = sorted(rows, key=lambda r: r["value"])[-2:]
    a1, a2 = L / finest["value"], L / second["value"]
    result = {}
    for name in CONVERGENCE_OBSERVABLES:
        if name not in finest["observables"]:
            continue
        x1 = finest["observables"][name]["mean"]
        x2 = second["observables"][name]["mean"]
        e1 = finest["observables"][name]["error"]
        e2 = second["observables"][name]["error"]
        change = x1 / x2 - 1.
        error = abs(x1 / x2) * float(np.hypot(e1 / x1, e2 / x2))
        result[name] = {
            "a_fm": [a2, a1], "change": change, "error": error,
            "converged": bool(abs(change) <= tolerance + 2. * error),
            # x(a) = x(0) + c a through both points
            "extrapolated_deviation": (x1 - x2) / (a2 - a1) * a1 / x1}
    return result


def energy_by_kt(spectra, scale):
    """Where dE/deta sits in k_T, from the spectra of the color-charge seeds
    of one setting and the g^2 mu scale [GeV]: the energy below and above
    ENERGY_CUTS times the scale (the binned energy below the cut, and the
    total dE, which includes the modes above the largest bin, minus it), and
    an estimate of the continuum energy, the binned energy below a cutoff
    k_c plus a tail k^4 dN/d^2k = C beyond it, 2 pi C/k_c, with C the mean of
    k^4 dN/d^2k between 0.85 k_c and k_c, for k_c = TAIL_FRACTIONS times the
    largest binned k_T. The bins are equally spaced."""
    kt = np.array(spectra[0]["kt"])
    weight = 2. * np.pi * kt * (kt[1] - kt[0])
    rows = {"total": [], "binned": []}
    for cut in ENERGY_CUTS:
        rows[f"below_{cut:g}"], rows[f"above_{cut:g}"] = [], []
    for fraction in TAIL_FRACTIONS:
        rows[f"tail_estimate_{fraction:g}"] = []
    for spectrum in spectra:
        de = np.array(spectrum["dE_d2k"]) * weight
        dn = np.array(spectrum["dN_d2k"])
        rows["total"].append(spectrum["dE"])
        rows["binned"].append(de.sum())
        for cut in ENERGY_CUTS:
            below = de[kt < cut * scale].sum()
            rows[f"below_{cut:g}"].append(below)
            rows[f"above_{cut:g}"].append(spectrum["dE"] - below)
        for fraction in TAIL_FRACTIONS:
            kc = fraction * kt[-1]
            band = (kt > 0.85 * kc) & (kt < kc)
            c = (kt ** 4 * dn)[band].mean()
            rows[f"tail_estimate_{fraction:g}"].append(
                de[kt < kc].sum() + 2. * np.pi * c / kc)
    return {name: statistics(values) for name, values in rows.items()}


def tail_exponent(spectrum, kmin, kmax):
    """The exponent p of k^4 dN/d^2k ~ k^p between kmin and kmax [GeV] (a
    straight line in log-log); p = 0 is the 1/k^4 tail of the dilute
    limit."""
    kt = np.array(spectrum["kt"])
    y = kt ** 4 * np.array(spectrum["dN_d2k"])
    inside = (kt >= kmin) & (kt <= kmax) & (y > 0.)
    if inside.sum() < 3:
        return None
    return float(np.polyfit(np.log(kt[inside]), np.log(y[inside]), 1)[0])


def mean_spectrum(spectra):
    return {"kt": spectra[0]["kt"],
            "dN_d2k": np.mean([s["dN_d2k"] for s in spectra],
                              axis=0).tolist(),
            "dE_d2k": np.mean([s["dE_d2k"] for s in spectra],
                              axis=0).tolist()}


def plot_spectra(summary, path):
    """k^4 dN/d^2k of the size scan and of the constant g^2 mu, one line
    per lattice (dashed above 0.4 times the largest k_T, where the lattice
    cutoff suppresses the spectrum)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    scans = [s for s in ("size", "uniform")
             if any("spectrum" in r for r in summary["scans"].get(s, []))]
    if not scans:
        return
    fig, axes = plt.subplots(1, len(scans), figsize=(5 * len(scans), 3.8),
                             squeeze=False)
    for ax, scan in zip(axes[0], scans):
        rows = [r for r in summary["scans"][scan] if "spectrum" in r]
        for i, row in enumerate(rows):
            kt = np.array(row["spectrum"]["kt"])
            y = kt ** 4 * np.array(row["spectrum"]["dN_d2k"])
            color = plt.cm.Blues(0.35 + 0.65 * i / max(1, len(rows) - 1))
            cut = 0.4 * kt[-1]
            ax.plot(kt[kt <= cut], y[kt <= cut], color=color,
                    label=f"{row['value']:g}$^2$")
            ax.plot(kt[kt >= cut], y[kt >= cut], color=color, ls="--")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$k_T$ [GeV]")
        ax.set_ylabel(r"$k_T^4\, dN/dy\, d^2k_T$ [GeV$^2$]")
        ax.set_title("constant $g^2\\mu$" if scan == "uniform"
                     else summary.get("variant") or "size scan", fontsize=9)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot(summary, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    scans = [("size", "a [fm]", lambda r: summary["L"] / r["value"]),
             ("ny", r"$N_y$", lambda r: r["value"]),
             ("m", r"$m$ [GeV]", lambda r: r["value"]),
             ("uniform", r"$g^2\mu\, a$ (constant $g^2\mu$)",
              lambda r: r["g2mu_a"])]
    scans = [s for s in scans if summary["scans"].get(s[0])]
    names = ["dE/deta(0.2)", "dE/deta(final)", "dN/dy", "eps_2", "eps_3"]
    fig, axes = plt.subplots(len(names), len(scans),
                             figsize=(3.4 * len(scans), 2.2 * len(names)),
                             sharex="col", squeeze=False)
    for col, (scan, label, x_of) in enumerate(scans):
        rows = summary["scans"][scan]
        x = np.array([x_of(r) for r in rows])
        for row, name in enumerate(names):
            ax = axes[row, col]
            if name not in rows[0]["observables"]:
                ax.axis("off")
                continue
            y = np.array([r["observables"][name]["mean"] for r in rows])
            e = np.array([r["observables"][name]["error"] for r in rows])
            ax.errorbar(x, y, yerr=e, fmt="o-", ms=4)
            if col == 0:
                ax.set_ylabel(name)
            # the constant g^2 mu has no eccentricity rows below dN/dy
            if row == len(names) - 1 or (scan == "uniform" and row == 2):
                ax.set_xlabel(label)
                ax.xaxis.set_tick_params(labelbottom=True)
            if scan in ("ny", "uniform"):
                ax.set_xscale("log")
            check = summary.get("convergence", {}).get(scan, {}).get(name)
            if check and not check["converged"]:
                # the values grow towards small a: the lower left is free
                ax.text(0.03, 0.06, f"NOT CONVERGED: "
                        f"{100 * check['change']:+.1f}% in the last step",
                        transform=ax.transAxes, color="#c0392b",
                        fontsize=7, fontweight="bold", va="bottom",
                        bbox={"facecolor": "white", "edgecolor": "none",
                              "alpha": 0.8})
    scale = (summary["geometry"] or {}).get("g2mu_GeV")
    if scale and scans[0][0] == "size":
        # the median g^2 mu a of the geometry on top of the size scan
        top = axes[0, 0].secondary_xaxis(
            "top", functions=(lambda a: a * scale["median"] / HBARC,
                              lambda x: x * HBARC / scale["median"]))
        top.set_xlabel(r"median $g^2\mu\, a$")
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
    parser.add_argument("--L", type=float, default=22.5,
                        help="[fm]; 22.5 fm suffices for b = 0, use 30 fm "
                             "for other impact parameters")
    parser.add_argument("--b", type=float, default=0.,
                        help="impact parameter [fm]; a multiple of L/128 "
                             "is a lattice vector for all sizes")
    parser.add_argument("--seed", type=int, default=3,
                        help="the seed of the geometry")
    parser.add_argument("--sizes", type=int, nargs="+",
                        default=[96, 192, 384, 768, 1024, 1280])
    parser.add_argument("--ny", type=int, nargs="*",
                        default=[1, 5, 10, 25, 50, 100],
                        help="the Ny scan; none skips it")
    parser.add_argument("--m", type=float, nargs="*",
                        default=[0.1, 0.2, 0.4, 0.8],
                        help="the m scan [GeV]; none skips it")
    parser.add_argument("--scan-size", type=int, default=768,
                        help="the size of the Ny and m scans, and of the "
                             "run that gives g^2 mu of the geometry")
    parser.add_argument("--uniform-g2mu", type=float, default=5.,
                        help="the constant g^2 mu of the reference scan "
                             "[GeV]")
    parser.add_argument("--uniform-L", type=float, default=3.,
                        help="L of the reference scan [fm]")
    parser.add_argument("--uniform-sizes", type=int, nargs="*",
                        default=[32, 64, 128, 256, 512],
                        help="the sizes of the reference scan; none skips "
                             "it")
    parser.add_argument("--jimwlk", action="store_true",
                        help="evolve the Wilson lines with the JIMWLK "
                             "parameters of the template (from "
                             "jimwlkInitialX to projectileX and targetX)")
    parser.add_argument("--x", type=float, default=None,
                        help="projectileX and targetX: the x of the Wilson "
                             "lines, with --jimwlk the x the evolution ends "
                             "at (default: those of the template)")
    parser.add_argument("--tail-fit", type=float, nargs=2, default=[6., 12.],
                        metavar=("KMIN", "KMAX"),
                        help="k_T range [GeV] of the fit k^4 dN/d^2k ~ k^p")
    parser.add_argument("--color-seeds", type=int, default=4)
    parser.add_argument("--maxtime", type=float, default=0.4, help="[fm/c]")
    parser.add_argument("--convergence-tolerance", type=float,
                        default=0.01,
                        help="largest relative change between the two "
                             "finest lattices (plus two standard deviations) "
                             "of a converged observable")
    parser.add_argument("--threads", type=int, default=4,
                        help="OpenMP threads per run")
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--quick", action="store_true",
                        help="sizes 64 and 128, Ny and m scans at 64 with "
                             "three values each, the constant g^2 mu at "
                             "sizes 32 and 64, 2 color-charge seeds")
    args = parser.parse_args()
    if args.quick:
        args.sizes, args.scan_size, args.color_seeds = [64, 128], 64, 2
        args.ny, args.m = [5, 25, 50], [0.2, 0.4, 0.8]
        args.uniform_sizes = [32, 64]

    with open(args.template) as f:
        template_values = dict(line.split()[:2] for line in f
                               if len(line.split()) >= 2
                               and not line.startswith("#"))
    ny0, m0 = int(template_values["Ny"]), float(template_values["m"])
    # the tag of the runs with --jimwlk or --x, e.g. jimwlk or x0.01
    variant = "_".join((["jimwlk"] if args.jimwlk else [])
                       + ([f"x{args.x:g}"] if args.x is not None else []))
    # numerical_parameters[_<variant>]_summary.json and the figures
    stem = os.path.join(os.path.abspath(args.datadir), "numerical_parameters"
                        + (f"_{variant}" if variant else ""))
    x_final = ({"projectileX": args.x, "targetX": args.x}
               if args.x is not None else {})
    path = os.path.abspath(args.ipglasma_path)
    settings = {
        "ipglasma": os.path.join(path, args.ipglasma_cmd),
        "template": os.path.abspath(args.template),
        "datadir": os.path.abspath(args.datadir),
        "threads": args.threads,
        "variant": variant,
        # the x of the initial Wilson lines, for the g^2 mu of the geometry
        "initial_x": ({"projectileX": template_values["jimwlkInitialX"],
                       "targetX": template_values["jimwlkInitialX"]}
                      if args.jimwlk else x_final),
        "overrides": dict({
            "L": args.L, "LOutput": args.L, "seed": args.seed,
            "maxTime": args.maxtime, "bMin": args.b, "bMax": args.b,
            "rotateReactionPlane": 0, "useJIMWLK": int(args.jimwlk),
            "jimwlkSaveSnapshots": 0, "useFluctuatingX": 0,
            "writeWilsonLines": 0, "writeTmunu": 1, "writeTmunuBinary": 1,
            "writeHydro": 0, "outputTimes": "0.2",
            "computeEccentricities": 1, "computeGluonMultiplicity": 1,
            "nucleusQsTableFileName": os.path.join(
                path, "qs2Adj_vs_Tp_vs_Y_240.in"),
            "nuclearConfigurationsPath": os.path.join(
                path, "nucleusConfigurations"),
        }, **x_final),
    }
    settings["tables"] = os.path.join(path, "tables")
    os.makedirs(settings["datadir"], exist_ok=True)

    color_seeds = range(1, args.color_seeds + 1)
    points = {"size": [dict(size=s, ny=ny0, m=m0) for s in args.sizes],
              "ny": [dict(size=args.scan_size, ny=n, m=m0) for n in args.ny],
              "m": [dict(size=args.scan_size, ny=ny0, m=m) for m in args.m],
              "uniform": [dict(size=s, ny=ny0, m=m0, g2mu=args.uniform_g2mu,
                               L=args.uniform_L)
                          for s in args.uniform_sizes]}
    for scan in points.values():
        for point in scan:
            point["variant"] = variant
    configs = {}
    for scan in points.values():
        for point in scan:
            for cs in color_seeds:
                cfg = dict(point, color_seed=cs)
                configs[run_tag(cfg)] = cfg
    scale = g2mu_scale(settings, args.scan_size)
    print(f"{len(configs)} runs (b = {args.b} fm, L = {args.L} fm, seed "
          f"{args.seed}) with up to {args.max_workers} workers ...",
          flush=True)
    start, results, failures = time.time(), {}, []
    with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
        futures = [pool.submit(run_one, cfg, settings)
                   for cfg in configs.values()]
        for i, future in enumerate(as_completed(futures), 1):
            cfg, obs, geometry, spectrum, error = future.result()
            if error:
                failures.append(f"{run_tag(cfg)}: {error}")
            else:
                results[run_tag(cfg)] = (cfg, obs, geometry, spectrum)
            if i % max(1, len(configs) // 10) == 0 or i == len(configs):
                print(f"  {i}/{len(configs)} runs done "
                      f"({time.time() - start:.0f} s)", flush=True)
    for failure in failures:
        print("FAILED:", failure)

    # the geometry must be the same in all runs: the N_part and N_coll lists
    # and b (<Q_s> is averaged on the lattice, so it depends a little on a);
    # the constant g^2 mu has none
    geometries = [g for _, _, g, _ in results.values() if g is not None]
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

    summary = {"settings": vars(args), "L": args.L, "variant": variant,
               "geometry": dict(geometries[0], g2mu_GeV=scale)
               if geometries else None,
               "checks": {"same_geometry": bool(same_geometry),
                          "color_seeds_differ": bool(differ),
                          "all_runs_succeeded": not failures
                          and scale is not None},
               "scans": {}}
    for scan, scan_points in points.items():
        rows = []
        names = UNIFORM_OBSERVABLES if scan == "uniform" else OBSERVABLES
        for p in scan_points:
            runs = [results[run_tag(dict(p, color_seed=cs))][1]
                    for cs in color_seeds
                    if run_tag(dict(p, color_seed=cs)) in results]
            if not runs:
                continue
            row = {"value": p[scan] if scan != "uniform" else p["size"],
                   "observables": {name: statistics([r[name] for r in runs])
                                   for name in names}}
            # g^2 mu a: of the constant g^2 mu, or the median of the geometry
            if scan == "uniform":
                row["g2mu_a"] = p["g2mu"] * p["L"] / p["size"] / HBARC
            elif scale:
                row["g2mu_a"] = scale["median"] * args.L / p["size"] / HBARC
            # where the energy sits in k_T, for the scans in the lattice
            # spacing
            spectra = [results[run_tag(dict(p, color_seed=cs))][3]
                       for cs in color_seeds
                       if run_tag(dict(p, color_seed=cs)) in results]
            g2mu = p["g2mu"] if scan == "uniform" else \
                (scale or {}).get("median")
            if scan in ("size", "uniform") and g2mu and \
                    all(s is not None for s in spectra):
                row["spectrum"] = mean_spectrum(spectra)
                row["energy_by_kt"] = energy_by_kt(spectra, g2mu)
                row["tail_exponent"] = tail_exponent(row["spectrum"],
                                                     *args.tail_fit)
            rows.append(row)
        summary["scans"][scan] = rows

    print(f"\nGeometry (seed {args.seed}, b = {args.b} fm): <Q_s>(avg) = "
          f"{summary['geometry']['Qs_avg']} GeV")
    if scale:
        print(f"g^2 mu of the geometry: median {scale['median']:.3g} GeV "
              f"(1/g^2 mu = {HBARC / scale['median']:.3g} fm, weighted with "
              f"(g^2 mu)^2), max {scale['max']:.3g} GeV")
    for scan, rows in summary["scans"].items():
        if not rows:
            continue
        names = list(rows[0]["observables"])[:5]
        title = {"uniform": f"constant g^2 mu = {args.uniform_g2mu:g} GeV, "
                            f"L = {args.uniform_L:g} fm, size"}.get(scan, scan)
        print(f"\n{title} scan:")
        print(f"  {scan:>7s}{'g2mu a':>8s}"
              + "".join(f"{n:>22s}" for n in names))
        for row in rows:
            cells = "".join(
                f"{row['observables'][n]['mean']:13.5g} "
                f"({row['observables'][n]['error']:.1g})"
                for n in names)
            g2mu_a = f"{row['g2mu_a']:8.3g}" if "g2mu_a" in row else " " * 8
            print(f"  {row['value']:>7g}{g2mu_a}{cells}")
        if any("energy_by_kt" in row for row in rows):
            print(f"  dE/deta by k_T at the final time [GeV] (cuts in units "
                  f"of the g^2 mu scale; tail estimate: the binned energy "
                  f"below k_c plus a 1/k^4 tail, k_c = "
                  f"{' and '.join(f'{f:g}' for f in TAIL_FRACTIONS)} "
                  f"k_max; p: k^4 dN/d^2k ~ k^p for k_T = "
                  f"{args.tail_fit[0]:g}-{args.tail_fit[1]:g} GeV):")
            columns = (["total"] + [f"below_{c:g}" for c in ENERGY_CUTS]
                       + [f"above_{c:g}" for c in ENERGY_CUTS]
                       + [f"tail_estimate_{f:g}" for f in TAIL_FRACTIONS])
            print(f"  {scan:>7s}" + "".join(f"{c:>18s}" for c in columns)
                  + f"{'p':>7s}")
            for row in rows:
                if "energy_by_kt" not in row:
                    continue
                e = row["energy_by_kt"]
                p_tail = row["tail_exponent"]
                print(f"  {row['value']:>7g}" + "".join(
                    f"{e[c]['mean']:11.4g} ({e[c]['error']:.1g})"
                    for c in columns)
                    + (f"{p_tail:7.2f}" if p_tail is not None else ""))
    summary["convergence"] = {
        scan: convergence(summary["scans"][scan],
                          args.uniform_L if scan == "uniform" else args.L,
                          args.convergence_tolerance)
        for scan in ("size", "uniform") if summary["scans"].get(scan)}
    print(f"\nconvergence in the lattice spacing (change between the two "
          f"finest lattices; converged if below "
          f"{100 * args.convergence_tolerance:g}% plus two standard "
          f"deviations):")
    not_converged = []
    for scan, checks in summary["convergence"].items():
        for name, check in checks.items():
            state = "ok" if check["converged"] else "NOT CONVERGED"
            print(f"  {scan:8s}{name:16s}{100 * check['change']:+7.2f}% "
                  f"(+- {100 * check['error']:.2f}%)  {state}")
            if not check["converged"]:
                not_converged.append((scan, name, check))
    if not_converged:
        line = "!" * 76
        print(f"\n{line}\nWARNING: NOT CONVERGED IN THE LATTICE SPACING")
        for scan, name, check in not_converged:
            print(f"  {name} ({scan} scan) changes by "
                  f"{100 * check['change']:+.1f}% "
                  f"(+- {100 * check['error']:.1f}%) from a = "
                  f"{check['a_fm'][0]:.3g} to {check['a_fm'][1]:.3g} fm; "
                  f"an extrapolation linear in a puts a = 0 "
                  f"{100 * check['extrapolated_deviation']:+.0f}% away.")
        print("The gluon number converges like a^2 once g^2 mu a < 1, the "
              "energy only like a\n(a deficit of about 0.3 g^2 mu a from the "
              "1/k^4 tail of the spectrum; with\nJIMWLK, whose tail is less "
              "steep, even more slowly), so these values still\ndepend on "
              "the lattice spacing; see validations/README.md.\n" + line)
    summary["not_converged"] = [f"{scan}: {name}"
                                for scan, name, _ in not_converged]
    passed = all(summary["checks"].values())
    for name, ok in summary["checks"].items():
        print(f"{'OK  ' if ok else 'FAIL'} {name}")
    summary["passed"] = bool(passed)
    with open(f"{stem}_summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    plot(summary, f"{stem}.pdf")
    plot_spectra(summary, f"{stem}_spectra.pdf")
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{stem}_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

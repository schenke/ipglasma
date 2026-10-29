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
enters. The expansion parameter of these corrections is G(0)/Ny: each
layer also carries the local fluctuations of the field, so D approaches
D_G as 1/Ny only for Ny >~ G(0). This script runs IP-Glasma for a scan of
g^2 mu, m, Ny (and optionally UVDamp), computes D(r) from the binary
Wilson lines along the axes and the diagonal, averaged over all sites and
over the two nuclei of each run (two independent configurations), and
checks, at the r_s where D_G = exp(-1/2), that

  1. D approaches D_G for large Ny: |D - D_G| at the largest Ny is below a
     fifth of that at the smallest,
  2. Ny (D - D_G) tends to a constant: the two largest Ny agree within
     three combined errors (plus 5%),
  3. the extrapolation of D - D_G = d_inf + c/Ny to Ny -> infinity, fitted
     for Ny >= 2 G(0) (rounded to an integer), vanishes within three
     statistical errors,
  4. g drops out: a run with g = 2 gives the same D as with g = 1.

A check that cannot be made (e.g. with fewer than three Ny in the fit
range) counts as failed. For G(0) <~ 1, D - D_G changes sign near r_s and is
not resolved there, while it is at smaller r; the checks are therefore made
at r_s where D - D_G is resolved there at the smallest Ny (more than five
standard errors), and otherwise at the r <= r_max where it is most
significant. Where fewer than three Ny of the scan satisfy Ny >= 2 G(0)
(m = 0 at large g^2 mu L), the scan of that setting is extended by
doubling the largest Ny (up to --ny-max), and settings in which D_G does
not fall to exp(-1/2) within N/4 are skipped. Each Ny has its own seeds,
so the configurations of different Ny are independent.

It also reports Q_s / g^2 mu, with D(r_s) = exp(-1/2) and Q_s^2 = 2/r_s^2,
the fundamental-representation saturation scale in coordinate space.

With --lappi, the script instead runs the setup of T. Lappi, Eur. Phys. J.
C 55 (2008) 285, arXiv:0711.3039 (g^2 mu L = 100, g^2 mu a = 0.5, m = 0
and 0.1 g^2 mu) and compares Q_s / g^2 mu with his Fig. 4 and Fig. 7,
with the definitions he used (issue #148):
  * the adjoint correlator in the factorized form, |S(r)|^2 - 1 with
    S(r) = sum_x tr V^+(x) V(x + r) / N^2 of each configuration (his code
    takes 2|S|^2 - 1, twice that plus a constant, with the same maximum of
    k~^2 C(k~)), and Q_s from that maximum, averaged over the polar angle;
  * Q_s as his scripts obtain it (lappi_procedure()): the spectra in the
    bins of his code, the maximum of k~^2 C of each configuration from a
    parabola around its largest bin, and their mean and its error over the
    configurations (his error bars), for the factorized adjoint and the
    fundamental spectrum;
  * the k = 0 mode of the kernel removed also for m > 0
    (removeZeroMode 1); the runs with m > 0 are made with and without it;
  * in coordinate space, C(r_s) = C(0) exp(-1/2), Q_s^2 = 2/r_s^2, with
    C(r) the 2D momentum integral of the binned spectrum of 2|S|^2 - 1,
    with k~ as the momentum and the k = 0 mode by its magnitude, as in his
    spectrum files (lappi_coordinate_profile()); that mode is negative, so
    its magnitude adds about 1.9 to C(r) at every r.
It checks that, for m = 0 and for m > 0 with the zero mode removed, the
adjoint and the fundamental sqrt(C_A/C_F) Q~_s of his procedure agree with
Fig. 4 within --lappi-n-sigma (3) combined errors at every Ny, and the
coordinate-space Q_s with Fig. 7 within --lappi-coord-tolerance (3%). It
also reports the adjoint correlator as written in his paper,
< |tr V^+(x) V(x + r)|^2 > - 1, from the adjoint Wilson lines
U_ab = 2 tr(t^a V^+ t^b V) with Fourier transforms like D, the coordinate
space of the lattice correlators, the maxima of the averaged spectra in
bins of width 2 pi/L (peak_momentum()) and in his bins (kfit.pl), and the
Gaussian fits of his scripts (runn.pl), which are not those of Fig. 4.

The results go to <datadir>/dipole_summary.json (with --lappi,
<datadir>/dipole_lappi_summary.json), the figures to <datadir>/dipole_*.pdf.

Usage (from validations/dipole_correlator/):
    # a quick test on small lattices, a few seconds
    python3 dipole_correlator.py --ipglasma-path ../.. --quick
    # the scan for the paper
    python3 dipole_correlator.py --ipglasma-path ../.. --max-workers 8
    # the comparison with Lappi (2008)
    python3 dipole_correlator.py --ipglasma-path ../.. --max-workers 8 --lappi

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
CA = 3.
HBARC = 0.1973269804  # GeV fm
NC = 3
DA = NC * NC - 1
# D - D_G at r_s counts as resolved at the smallest Ny above this many
# standard errors; otherwise the checks are made where it is most significant
SIGNIFICANCE = 5.

# the generators t^a = lambda^a / 2
GELL_MANN = np.zeros((8, 3, 3), dtype=complex)
GELL_MANN[0][0, 1] = GELL_MANN[0][1, 0] = 1.
GELL_MANN[1][0, 1], GELL_MANN[1][1, 0] = -1j, 1j
GELL_MANN[2][0, 0], GELL_MANN[2][1, 1] = 1., -1.
GELL_MANN[3][0, 2] = GELL_MANN[3][2, 0] = 1.
GELL_MANN[4][0, 2], GELL_MANN[4][2, 0] = -1j, 1j
GELL_MANN[5][1, 2] = GELL_MANN[5][2, 1] = 1.
GELL_MANN[6][1, 2], GELL_MANN[6][2, 1] = -1j, 1j
GELL_MANN[7] = np.diag([1., 1., -2.]) / np.sqrt(3.)
TA = GELL_MANN / 2.

# Q_s / g^2 mu of arXiv:0711.3039 (g^2 mu L = 100, g^2 mu a = 0.5), read
# off the vector graphics of its figures: the marker centres and, for
# Fig. 4, the error bars (half their length). Fig. 4: momentum space,
# adjoint Q_s and fundamental sqrt(C_A/C_F) Q~_s, {Ny: (value, error)};
# Fig. 7: coordinate space, adjoint, {Ny: value} (no error bars drawn).
# Keys: m / g^2 mu.
LAPPI_FIG4_ADJ = {
    0.: {1: (0.5697, 0.0019), 2: (0.6200, 0.0020), 4: (0.6949, 0.0017),
         10: (0.8422, 0.0011), 20: (0.9605, 0.0011), 50: (1.0765, 0.0011),
         100: (1.1264, 0.0016), 150: (1.1434, 0.0015)},
    0.1: {1: (0.5190, 0.0021), 2: (0.6167, 0.0008), 4: (0.7034, 0.0016),
          10: (0.7939, 0.0014), 20: (0.8308, 0.0017), 50: (0.8563, 0.0017),
          150: (0.8678, 0.0022)},
}
LAPPI_FIG4_FUND = {
    0.: {1: (0.5833, 0.0055), 2: (0.6228, 0.0060), 4: (0.6901, 0.0101),
         10: (0.8449, 0.0057), 20: (0.9724, 0.0066), 50: (1.0748, 0.0081),
         100: (1.1341, 0.0134), 150: (1.1475, 0.0141)},
    0.1: {1: (0.5024, 0.0088), 2: (0.5922, 0.0041), 4: (0.6738, 0.0076),
          10: (0.7548, 0.0123), 20: (0.7936, 0.0098), 50: (0.8152, 0.0103),
          150: (0.8142, 0.0199)},
}
LAPPI_FIG7_COORD = {
    0.: {1: 0.6367, 10: 0.9209, 50: 1.1659},
    0.1: {1: 0.6227, 10: 0.9340, 50: 1.0031},
}


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


def fundamental_spectrum(V):
    """C~(k) = sum_ij |sum_x e^{-ik.x} V_ij(x)|^2 / N^2, the Fourier
    transform of tr V^+(x) V(x + r) averaged over x (lattice units)."""
    n = V.shape[0]
    f = np.fft.fft2(V, axes=(0, 1))
    return (np.abs(f) ** 2).sum(axis=(2, 3)) / (n * n)


def adjoint_lines(V):
    """The adjoint Wilson lines U_ab(x) = 2 tr[t^a V^+(x) t^b V(x)], real
    orthogonal 8x8 matrices; shape (N, N, 8, 8)."""
    Vd = np.conj(np.swapaxes(V, -1, -2))
    # V^+ t^b V for all b, shape (N, N, 8, 3, 3)
    M = np.matmul(Vd[:, :, None], np.matmul(TA, V[:, :, None]))
    return 2. * np.einsum("aji,xybij->xyab", TA, M).real


def adjoint_spectrum(V):
    """C(k) = sum_ab |sum_x e^{-ik.x} U_ab(x)|^2 / N^2, the Fourier
    transform of sum_ab U_ab(x) U_ab(x + r) = |tr V^+(x) V(x + r)|^2 - 1
    averaged over x (lattice units)."""
    U = adjoint_lines(V)
    n = U.shape[0]
    f = np.fft.fft2(U, axes=(0, 1))
    return (np.abs(f) ** 2).sum(axis=(2, 3)) / (n * n)


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


def momentum_bins(n):
    """The lattice momentum k~ a = 2 [sin^2(kx a/2) + sin^2(ky a/2)]^(1/2)
    of all modes, their bin (width 2 pi/N, centred on multiples of it)
    and the mean k~ a of each bin; bin 0 is k = 0."""
    k = 2. * np.pi * np.fft.fftfreq(n)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    kt = 2. * np.sqrt(np.sin(kx / 2.) ** 2 + np.sin(ky / 2.) ** 2)
    index = np.floor(kt * n / (2. * np.pi) + 0.5).astype(int)
    count = np.bincount(index.ravel())
    with np.errstate(invalid="ignore"):
        mean = np.bincount(index.ravel(), kt.ravel()) / count
    return kt, index, mean


def radial_spectrum(ck, kt, index):
    """k~^2 C(k), averaged over the modes of each bin (polar angle)."""
    count = np.bincount(index.ravel())
    with np.errstate(invalid="ignore"):
        return np.bincount(index.ravel(), (kt ** 2 * ck).ravel()) / count


def peak_momentum(k, f):
    """The k~ of the maximum of the binned k~^2 C: the vertex of the
    parabola through the largest bin and two bins on either side."""
    ok = np.isfinite(f)
    k, f = k[ok][1:], f[ok][1:]
    i = int(np.argmax(f))
    s = slice(max(i - 2, 0), i + 3)
    c = np.polyfit(k[s], f[s], 2)
    if c[0] >= 0.:
        return float("nan")
    return float(-c[1] / (2. * c[0]))


def lappi_bins(kt, bins):
    """The bins of the spectrum files of the code of arXiv:0711.3039
    (#148): bins of width 2 sqrt(2)/bins in k~ a (N bins for the adjoint,
    2N for the fundamental spectrum); the bin of each mode and the mean
    k~ a of each bin, NaN for an empty bin."""
    index = np.floor(kt / (2. * np.sqrt(2.) + 1e-4) * bins).astype(int)
    return index, binned_mean(kt, index, bins)


def binned_mean(values, index, bins):
    """The mean of values over the modes of each bin; NaN for an empty
    bin."""
    count = np.bincount(index.ravel(), minlength=bins)
    with np.errstate(invalid="ignore"):
        return np.bincount(index.ravel(), values.ravel(),
                           minlength=bins) / count


def lappi_gaussian_fit(k, c, err, q, kmin, kmax, iterations=5):
    """Q_s from the fit of ln C = const - k~^2 / Q_s^2 to the spectrum C
    averaged over the configurations, weighted with its relative errors,
    as runn.pl of arXiv:0711.3039 does (#148): over kmin < k~ < kmax, then
    over 0.1 Q_s < k~ < 1.3 Q_s, iterations times in all; q is the start
    value, kept while the fit has no maximum. Returns Q_s, its error and
    the last upper limit of the fit."""
    ok = np.isfinite(c) & (c > 0.) & (err > 0.)
    k, c, err = k[ok], c[ok], err[ok]
    for _ in range(iterations):
        s = (k > kmin) & (k < kmax)
        x, y, w = k[s] ** 2, np.log(c[s]), (c[s] / err[s]) ** 2
        S, Sx, Sxx = w.sum(), (w * x).sum(), (w * x * x).sum()
        delta = S * Sxx - Sx * Sx
        b = (S * (w * x * y).sum() - Sx * (w * y).sum()) / delta
        if b < 0.:
            q = np.sqrt(-1. / b)
        kmin, kmax = 0.1 * q, 1.3 * q
    return float(q), float(-0.5 * np.sqrt(S / delta) * q / b), float(kmax)


def lappi_peak_fit(k, f, err=None):
    """The k~ of the maximum of f = k~^2 C from a parabola around its
    largest bin, as arXiv:0711.3039 does (#148). With err (kfit.pl, for the
    averaged spectrum), the fit is weighted with 1/err^2 and its window,
    from 0.95 to 1.05 times the k~ of the largest bin, is widened by 5% on
    either side until it holds at least 3 bins and either 10 bins or an
    upper limit of 1.5 times that k~. Without err (kfitone.pl, for one
    configuration), the weight is k~ and the window grows by 10% from that
    k~ until it holds at least 4 bins and the same; the result is NaN if
    the parabola has no maximum or it lies more than 50% away from the
    largest bin."""
    ok = np.isfinite(f) & (k > 0.)
    k, f = k[ok], f[ok]
    kpeak = k[int(np.argmax(f))]
    if err is not None:
        w, step, lower, upper, least = 1. / err[ok] ** 2, 0.05, 0.95, 1.05, 3
    else:
        w, step, lower, upper, least = k, 0.1, 1., 1., 4
    lower, upper = lower * kpeak, upper * kpeak
    while True:
        lower, upper = lower * (1. - step), upper * (1. + step)
        s = (k > lower) & (k < upper)
        if s.sum() >= least and (s.sum() >= 10 or upper >= 1.5 * kpeak):
            break
    c = np.polyfit(k[s], f[s], 2, w=np.sqrt(w[s]))
    peak = -c[1] / (2. * c[0])
    if err is None and not (c[0] < 0. and abs(peak / kpeak - 1.) < 0.5):
        return float("nan")
    return float(peak)


def lappi_coordinate_profile(k, count, spectrum, c0, n, r):
    """The coordinate-space adjoint correlator as arXiv:0711.3039 computes
    it (#148): the 2D momentum integral of the binned spectrum of
    2|S(r)|^2 - 1, with the lattice momentum k~ as the momentum,
    sum_bins count C(k~) J_0(k~ r) / n^2, and the k = 0 mode by its
    magnitude, as in the spectrum files of that code. spectrum is the
    binned k~^2 C of |S|^2 - 1 (lappi_observables), c0 its k = 0 mode."""
    from scipy.special import j0
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(k > 0., spectrum / k ** 2, 0.)
    c = np.nan_to_num(c)
    integral = (count[1:, None] * c[1:, None]
                * j0(k[1:, None] * r[None, :])).sum(axis=0)
    return (abs(2. * c0 + n * n) + 2. * integral) / (n * n)


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


def seed_of(ny, i, first_seed):
    """The seed of the i-th run with Ny layers: each Ny has its own seeds,
    so the configurations of different Ny are independent."""
    return 1000 * ny + first_seed + i


def run_tag(cfg, settings):
    return (f"N{settings['size']}_L{settings['L']:g}"
            f"_g2mu{cfg['g2mu']:g}_m{cfg['m']:g}_uv{cfg['uvdamp']:g}"
            f"_ny{cfg['ny']}_g{cfg['g']:g}_seed{cfg['seed']}"
            + ("_nozero" if cfg.get("remove_zero_mode") else ""))


def lappi_observables(V, kt, index):
    """The adjoint correlator of one configuration (profiles in r) and the
    binned k~^2 C(k~) of the adjoint, factorized and fundamental ones;
    also the fundamental C(k~) and the magnitude of the factorized one in
    the bins of arXiv:0711.3039 (lappi_bins)."""
    out = {}
    cka = adjoint_spectrum(V)
    adj = profiles(np.fft.ifft2(cka).real)
    ckf = fundamental_spectrum(V)
    # S(r) = sum_x tr V^+(x) V(x + r) / N^2, complex
    s = np.fft.ifft2(ckf)
    fact = np.abs(s) ** 2 - 1.
    ck_fact = np.fft.fft2(fact).real
    out["adj_axis"], out["adj_diag"] = adj["axis"], adj["diag"]
    out["fact_axis"] = profiles(fact)["axis"]
    # the k = 0 modes, which k~^2 C(k~) does not show
    out["c0_adj"] = np.array(cka[0, 0])
    out["c0_fact"] = np.array(ck_fact[0, 0])
    out["spec_adj"] = radial_spectrum(cka, kt, index)
    out["spec_fund"] = radial_spectrum(ckf, kt, index)
    out["spec_fact"] = radial_spectrum(ck_fact, kt, index)
    # the spectrum files of that code: the mean of tr U(k) U^+(k) and of
    # the magnitude of the transform of 2|S|^2 - 1 (twice that of
    # |S|^2 - 1 for k != 0) over each bin
    n = V.shape[0]
    out["lappi_fund"] = binned_mean(ckf, lappi_bins(kt, 2 * n)[0], 2 * n)
    out["lappi_fact"] = binned_mean(np.abs(ck_fact), lappi_bins(kt, n)[0], n)
    return out


def run_one(cfg, settings):
    """Run IP-Glasma for one configuration and return the D(r) profiles of
    its two nuclei (cached in <workdir>/profiles.npz); with --lappi also
    the correlators of lappi_observables."""
    workdir = os.path.join(settings["datadir"], run_tag(cfg, settings))
    cache = os.path.join(workdir, "profiles.npz")
    if os.path.exists(cache):
        with np.load(cache) as data:
            cached = {k: data[k] for k in data.files}
        if not settings["lappi"] or "lappi_fact0" in cached:
            return cfg, cached, None
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
    if cfg.get("remove_zero_mode"):
        overrides["removeZeroMode"] = 1
    write_input(settings["template"], overrides,
                os.path.join(workdir, "input"))
    env = dict(os.environ, OMP_NUM_THREADS=str(settings["threads"]))
    with open(os.path.join(workdir, "run.log"), "w") as log:
        result = subprocess.run([settings["ipglasma"], "input"], cwd=workdir,
                                stdout=log, stderr=subprocess.STDOUT, env=env)
    files = sorted(glob.glob(os.path.join(workdir, "WilsonLine_*")))
    if result.returncode != 0 or len(files) != 2:
        return cfg, None, (f"ipglasma failed (exit {result.returncode}, "
                           f"{len(files)} Wilson-line files), see {workdir}")
    out = {}
    if settings["lappi"]:
        kt, index, kmean = momentum_bins(settings["size"])
        out["k_bins"] = kmean
        out["k_lappi_fund"] = lappi_bins(kt, 2 * settings["size"])[1]
        out["k_lappi_fact"] = lappi_bins(kt, settings["size"])[1]
    for i, path in enumerate(files):
        V, meta = read_wilson_lines(path)
        prof = profiles(dipole_correlator(V))
        out[f"axis{i}"] = prof["axis"]
        out[f"diag{i}"] = prof["diag"]
        out["a_fm"] = np.array(meta["a_fm"])
        if settings["lappi"]:
            for name, value in lappi_observables(V, kt, index).items():
                out[f"{name}{i}"] = value
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


def setting_two_point(key, settings):
    """G(r) of a setting (g2mu, m, UVDamp) and its lattice spacing."""
    g2mu, m, uv = key
    a = settings["L"] / settings["size"]
    return lattice_two_point(settings["size"], g2mu * a / HBARC,
                             m * a / HBARC, uv / a * HBARC), a


def setting_ny(base, g0, fit_from, ny_max):
    """The Ny of a setting: those of --ny, extended by doubling the
    largest until three of them lie in the fit range Ny >= fit_from G(0)
    (at most up to ny_max)."""
    nys = sorted(base)
    while (sum(ny >= fit_from * g0 for ny in nys) < 3
           and 2 * nys[-1] <= ny_max):
        nys.append(2 * nys[-1])
    return nys


def analyze_group(key, runs, settings):
    g2mu, m, uv = key
    n = settings["size"]
    G, a = setting_two_point(key, settings)
    g0 = float(G[0, 0])
    exact = profiles(np.exp(-CF * (g0 - G)))
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

    # where the checks are made: at r_s if D - D_G is resolved there at the
    # smallest Ny, else (for G(0) <~ 1, where it changes sign near r_s) at
    # the r <= r_max where it is most significant
    first = rows[0]["axis"]
    with np.errstate(invalid="ignore", divide="ignore"):
        significance = np.nan_to_num(np.abs(
            np.array(first["D"]) - exact["axis"]) / np.array(first["error"]))
    resolved = bool(significance[r_star] > SIGNIFICANCE)
    r_check = r_star if resolved else \
        int(1 + np.argmax(significance[1:r_max + 1]))
    for row in rows:
        row["dev_rcheck"] = float(row["axis"]["D"][r_check]
                                  - exact["axis"][r_check])
        row["error_rcheck"] = float(row["axis"]["error"][r_check])
    # Ny >= fit_from G(0), rounded to an integer
    ny_from = max(1., float(round(settings["fit_from"] * g0)))
    group = {"g2mu_GeV": g2mu, "m_GeV": m, "UVDamp_invGeV": uv,
             "a_fm": a, "G0": g0, "resolved": resolved,
             "significance_rstar": float(significance[r_star]),
             "r_max": r_max, "r_star": r_star, "r_check": r_check,
             "exact": {d: exact[d].tolist() for d in exact},
             "Qs_exact_over_g2mu": float(
                 np.sqrt(2.) / (crossing(r, exact["axis"]) * a) * HBARC / g2mu),
             "rows": rows}

    # acceptance
    # (a check that cannot be made is None, which counts as a failure)
    ny = np.array([row["ny"] for row in rows], dtype=float)
    dev = np.array([row["dev_rcheck"] for row in rows])
    err = np.array([row["error_rcheck"] for row in rows])
    checks = {"approaches_exact": None, "ny_dev_constant": None,
              "extrapolation_vanishes": None}
    if len(rows) >= 2:
        checks["approaches_exact"] = bool(abs(dev[-1]) < 0.2 * abs(dev[0]))
        # Ny (D - D_G) at the two largest Ny, within three combined errors
        a1, a2 = ny[-2] * dev[-2], ny[-1] * dev[-1]
        e = np.hypot(ny[-2] * err[-2], ny[-1] * err[-1])
        checks["ny_dev_constant"] = bool(abs(a1 - a2) < 3. * e + 0.05 * abs(a2))
    use = ny >= ny_from
    group["fit"] = {"ny_from": ny_from, "ny": ny[use].astype(int).tolist()}
    if use.sum() >= 3:
        d_inf, d_err, c = weighted_fit(ny[use], dev[use], err[use])
        group["fit"].update({"d_inf": d_inf, "error": d_err, "c": c})
        checks["extrapolation_vanishes"] = bool(abs(d_inf) < 3. * d_err)
    group["checks"] = checks
    return group


def jackknife(samples, estimate):
    """estimate(mean of samples) and its delete-one jackknife error;
    samples has the configurations along axis 0."""
    samples = np.asarray(samples)
    n = len(samples)
    total = samples.sum(axis=0)
    value = estimate(total / n)
    partial = np.array([estimate((total - s) / (n - 1)) for s in samples])
    # an estimate that fails (no maximum, no crossing) for one of the
    # leave-one-out samples makes the error undefined
    if np.isnan(partial).any():
        return float(value), float("nan")
    return float(value), float(np.sqrt((n - 1) * np.var(partial)))


def lappi_procedure(outs, g2mu_lat, m_lat):
    """Q_s / g^2 mu as the scripts of arXiv:0711.3039 obtain it (#148),
    from the spectra in its bins (lappi_observables), averaged over the
    configurations, with errors sqrt(var / n) in each bin:
      * runn.pl: Gaussian fits (lappi_gaussian_fit), first of the
        fundamental spectrum, from 0.1 m to 0.65 g^2 mu with the start value
        0.5 g^2 mu, then of the factorized adjoint one, from 0.05 / a to
        the last upper limit of the fundamental fit; the fundamental Q~_s
        multiplied by sqrt(C_A / C_F);
      * kfit.pl: the maximum of k~^2 C of the averaged spectra
        (lappi_peak_fit with errors);
      * kfitone.pl: the mean of the maxima of k~^2 C of the single
        configurations (lappi_peak_fit without errors) that do not fail,
        with the error sqrt(var / n)."""
    row = {}
    spectra = {}
    for name in ("fund", "fact"):
        k = outs[0][f"k_lappi_{name}"]
        c = np.array([out[f"lappi_{name}{i}"] for out in outs
                      for i in (0, 1)])
        spectra[name] = (k, c.mean(axis=0), c.std(axis=0) / np.sqrt(len(c)),
                         c)
    k, c, err, _ = spectra["fund"]
    q, dq, kmax = lappi_gaussian_fit(k, c, err, 0.5 * g2mu_lat,
                                     0.1 * m_lat, 0.65 * g2mu_lat)
    scale = np.sqrt(CA / CF) / g2mu_lat
    row["Qs_runn_fund_scaled"] = q * scale
    row["Qs_runn_fund_scaled_error"] = dq * scale
    k, c, err, _ = spectra["fact"]
    q, dq, _ = lappi_gaussian_fit(k, c, err, q, 0.05, kmax)
    row["Qs_runn"], row["Qs_runn_error"] = q / g2mu_lat, dq / g2mu_lat
    for name, suffix, factor in (("fact", "", 1.),
                                 ("fund", "_fund_scaled", np.sqrt(CA / CF))):
        k, c, err, single = spectra[name]
        scale = factor / g2mu_lat
        row[f"Qs_kfit{suffix}"] = scale * lappi_peak_fit(k, k ** 2 * c,
                                                         k ** 2 * err)
        peaks = np.array([lappi_peak_fit(k, k ** 2 * f) for f in single])
        peaks = scale * peaks[np.isfinite(peaks)]
        row[f"Qs_kfitone{suffix}"] = float(peaks.mean())
        row[f"Qs_kfitone{suffix}_error"] = float(peaks.std()
                                                 / np.sqrt(len(peaks)))
        row[f"kfitone{suffix}_rejected"] = int(len(single) - len(peaks))
    return row


def analyze_lappi(key, runs, settings):
    """Q_s / g^2 mu of the Lappi setup for one m (and removeZeroMode), for
    each Ny."""
    g2mu, m, uv, removed = key
    g2mu_lat = g2mu * settings["L"] / settings["size"] / HBARC
    m_lat = m * settings["L"] / settings["size"] / HBARC
    m_rel = round(m / g2mu, 6)
    r = np.arange(settings["size"] // 2 + 1)
    count = np.bincount(settings["bin_index"].ravel()).astype(float)
    n = settings["size"]
    rows = []
    for ny in sorted({cfg["ny"] for cfg, _ in runs}):
        outs = [out for cfg, out in runs if cfg["ny"] == ny]
        k = outs[0]["k_bins"]

        def collect(name):
            return [out[f"{name}{i}"] for out in outs for i in (0, 1)]

        def momentum(f):
            return peak_momentum(k, f) / g2mu_lat

        def coordinate(c):
            return np.sqrt(2.) / crossing(r, c / DA) / g2mu_lat

        row = {"ny": ny, "configurations": 2 * len(outs)}
        row["Qs_adj"], row["Qs_adj_error"] = jackknife(
            collect("spec_adj"), momentum)
        row["Qs_fact"], row["Qs_fact_error"] = jackknife(
            collect("spec_fact"), momentum)
        q, e = jackknife(collect("spec_fund"), momentum)
        row["Qs_fund_scaled"] = float(np.sqrt(CA / CF) * q)
        row["Qs_fund_scaled_error"] = float(np.sqrt(CA / CF) * e)
        row["Qs_coord_adj"], row["Qs_coord_adj_error"] = jackknife(
            collect("adj_axis"), coordinate)
        row["Qs_coord_fact"], row["Qs_coord_fact_error"] = jackknife(
            collect("fact_axis"), coordinate)
        # the coordinate space of Lappi: the 2D momentum integral of the
        # binned spectrum, normalized at r = 0

        def coordinate_lappi(sample):
            c = lappi_coordinate_profile(k, count, sample[:-1], sample[-1],
                                         n, r)
            return np.sqrt(2.) / crossing(r, c / c[0]) / g2mu_lat

        samples = [np.append(out[f"spec_fact{i}"], out[f"c0_fact{i}"])
                   for out in outs for i in (0, 1)]
        row["Qs_coord_lappi"], row["Qs_coord_lappi_error"] = jackknife(
            samples, coordinate_lappi)
        row.update(lappi_procedure(outs, g2mu_lat, m_lat))
        for name, table in (("adj", LAPPI_FIG4_ADJ),
                            ("fund", LAPPI_FIG4_FUND)):
            if ny in table.get(m_rel, {}):
                row[f"lappi_{name}"] = list(table[m_rel][ny])
        if ny in LAPPI_FIG7_COORD.get(m_rel, {}):
            row["lappi_coord"] = LAPPI_FIG7_COORD[m_rel][ny]
        rows.append(row)
    return {"g2mu_GeV": g2mu, "m_GeV": m, "m_over_g2mu": m_rel,
            "remove_zero_mode": bool(removed),
            "g2mu_L": g2mu * settings["L"] / HBARC, "g2mu_a": g2mu_lat,
            "rows": rows}


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
                  fr"UVDamp $= {group['UVDamp_invGeV']:g}$ GeV$^{{-1}}$, "
                  fr"$G(0) = {group['G0']:.3g}$", fontsize=9)
    top.legend(fontsize=8)
    bottom.axhline(0., color="k", lw=0.5)
    bottom.set_xlabel(r"$r$ [fm]")
    bottom.set_ylabel(r"$N_y\,(D - D_\mathrm{G})$")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_lappi(groups, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(groups), figsize=(5. * len(groups), 4.),
                             squeeze=False)
    for ax, group in zip(axes[0], groups):
        rows = group["rows"]
        ny = np.array([row["ny"] for row in rows])
        for name, fmt, label in (
                ("Qs_adj", "o-", r"adjoint, $\langle|\mathrm{tr}|^2\rangle-1$"),
                ("Qs_kfitone", "s-", r"factorized, $|S|^2-1$"),
                ("Qs_kfitone_fund_scaled", "^-",
                 r"fundamental $\times\sqrt{C_A/C_F}$"),
                ("Qs_coord_lappi", "v-",
                 "coordinate space, momentum integral of $|S|^2-1$")):
            ax.errorbar(ny, [row[name] for row in rows],
                        yerr=[row[name + "_error"] for row in rows],
                        fmt=fmt, ms=4, lw=0.8, label=label)
        for name, fmt, label in (("lappi_adj", "k*", "Lappi Fig. 4, adjoint"),
                                 ("lappi_fund", "kx",
                                  "Lappi Fig. 4, fundamental")):
            sel = [row for row in rows if name in row]
            ax.errorbar([row["ny"] for row in sel],
                        [row[name][0] for row in sel],
                        yerr=[row[name][1] for row in sel], fmt=fmt, ms=6,
                        label=label)
        sel = [row for row in rows if "lappi_coord" in row]
        ax.plot([row["ny"] for row in sel],
                [row["lappi_coord"] for row in sel], "k+", ms=8,
                label="Lappi Fig. 7, coordinate space")
        ax.set_xscale("log")
        ax.set_xlabel(r"$N_y$")
        ax.set_ylabel(r"$Q_s/g^2\mu$")
        zero = "" if group["m_over_g2mu"] == 0. else (
            ", zero mode removed" if group["remove_zero_mode"]
            else ", zero mode kept")
        ax.set_title(fr"$g^2\mu L = {group['g2mu_L']:.0f}$, "
                     fr"$g^2\mu a = {group['g2mu_a']:.2f}$, "
                     fr"$m = {group['m_over_g2mu']:g}\,g^2\mu${zero}",
                     fontsize=9)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def run_all(configs, settings, max_workers):
    print(f"{len(configs)} runs of size {settings['size']}, "
          f"L = {settings['L']:g} fm, with up to {max_workers} workers ...",
          flush=True)
    start = time.time()
    results, failures = {}, []
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        # the longest runs first
        futures = [pool.submit(run_one, cfg, settings)
                   for cfg in sorted(configs, key=lambda c: -c["ny"])]
        for i, future in enumerate(as_completed(futures), 1):
            cfg, out, error = future.result()
            if error:
                failures.append(f"{run_tag(cfg, settings)}: {error}")
            else:
                results[run_tag(cfg, settings)] = (cfg, out)
            if i % max(1, len(configs) // 10) == 0 or i == len(configs):
                print(f"  {i}/{len(configs)} runs done "
                      f"({time.time() - start:.0f} s)", flush=True)
    for failure in failures:
        print("FAILED:", failure)
    return results, failures


def main_lappi(args, settings):
    """The setup of arXiv:0711.3039: g^2 mu L = 100, g^2 mu a = 0.5,
    m = 0 and 0.1 g^2 mu."""
    # for m > 0 with and without the zero mode of the kernel, which
    # arXiv:0711.3039 leaves out also there (#148)
    groups = [(args.g2mu[0], m, 0., removed) for m in args.m
              for removed in ((False, True) if m > 0. else (False,))]
    configs = [{"g2mu": g, "m": m, "uvdamp": uv, "ny": ny, "g": 1.,
                "remove_zero_mode": removed,
                "seed": seed_of(ny, i, args.first_seed)}
               for (g, m, uv, removed), ny, i in itertools.product(
                   groups, args.ny, range(args.seeds))]
    results, failures = run_all(configs, settings, args.max_workers)
    settings["bin_index"] = momentum_bins(args.size)[1]
    summary = {"settings": {k: v for k, v in vars(args).items()},
               "groups": [], "failures": failures}
    for key in groups:
        runs = [(cfg, out) for cfg, out in results.values()
                if (cfg["g2mu"], cfg["m"], cfg["uvdamp"],
                    cfg.get("remove_zero_mode", False)) == key]
        if runs:
            summary["groups"].append(analyze_lappi(key, runs, settings))

    print("\nQ_s / g^2 mu in the setup of Lappi, arXiv:0711.3039; momentum "
          "space: maximum of k~^2 C(k~),\ncoordinate space: C(r_s) = 8 "
          "exp(-1/2), Q_s^2 = 2/r_s^2, along the axes:")
    for group in summary["groups"]:
        print(f"\n g2mu L = {group['g2mu_L']:.1f}, g2mu a = "
              f"{group['g2mu_a']:.3f}, m = {group['m_over_g2mu']:g} g2mu"
              + (", zero mode removed" if group["remove_zero_mode"]
                 and group["m_over_g2mu"] > 0. else ""))
        print("    Ny conf  adjoint         factorized      Lappi Fig. 4 adj"
              "  sqrt(CA/CF) fund  Lappi Fig. 4 fund  adj, coordinate"
              "  Lappi Fig. 7")
        for row in group["rows"]:
            def entry(name):
                return f"{row[name]:.4f}({row[name + '_error']:.4f})"

            def lappi(name):
                if name not in row:
                    return f"{'-':>15}"
                return f"{row[name][0]:.4f}({row[name][1]:.4f})"
            coord = (f"{row['lappi_coord']:.4f}" if "lappi_coord" in row
                     else "-")
            print(f"  {row['ny']:4d} {row['configurations']:4d}  "
                  f"{entry('Qs_adj')}  {entry('Qs_fact')}  "
                  f"{lappi('lappi_adj')}   {entry('Qs_fund_scaled')}   "
                  f"{lappi('lappi_fund')}    {entry('Qs_coord_adj')}  "
                  f"{coord:>12}")
    print("\nWith the procedure of arXiv:0711.3039 (#148): the mean of the "
          "maxima of the single\nconfigurations in his bins, adjoint "
          "(factorized) and sqrt(CA/CF) fundamental, with\ntheir deviation "
          "from Fig. 4 in combined errors; the Gaussian fits of his "
          "runn.pl;\nthe coordinate space from the 2D momentum integral of "
          "the binned spectrum:")
    checks, deviations = {}, {}
    for group in summary["groups"]:
        lappi_like = (group["m_over_g2mu"] == 0.
                      or group["remove_zero_mode"])
        label = f"m = {group['m_over_g2mu']:g} g2mu" + (
            "" if group["m_over_g2mu"] == 0. else
            ", zero mode removed" if group["remove_zero_mode"]
            else ", zero mode kept")
        print(f" {label}")
        pulls = {"adj": [], "fund": []}
        for row in group["rows"]:
            line = f"  {row['ny']:4d}"
            for name, ours in (("adj", "Qs_kfitone"),
                               ("fund", "Qs_kfitone_fund_scaled")):
                line += f"  {name} {row[ours]:.4f}({row[ours + '_error']:.4f})"
                if f"lappi_{name}" in row:
                    value, error = row[f"lappi_{name}"]
                    pull = (row[ours] - value) / np.hypot(
                        row[ours + "_error"], error)
                    row[f"pull_{name}"] = float(pull)
                    pulls[name].append(abs(pull))
                    line += f" {pull:+5.1f}"
                else:
                    line += "      "
            ref = (f"{row['lappi_coord']:.4f}" if "lappi_coord" in row
                   else "-")
            print(f"{line}  Gaussian fit {row['Qs_runn']:.4f};  coordinate "
                  f"{row['Qs_coord_lappi']:.4f}"
                  f"({row['Qs_coord_lappi_error']:.4f}), Fig. 7 {ref}")
        if not lappi_like:
            continue
        # the checks, where the definitions are those of arXiv:0711.3039
        coordinate = [abs(row["Qs_coord_lappi"] / row["lappi_coord"] - 1.)
                      for row in group["rows"] if "lappi_coord" in row]
        for name, text in (("adj", "adjoint"), ("fund", "fundamental")):
            if pulls[name]:
                check = (f"{label}: the {text} Q_s is within "
                         f"{args.lappi_n_sigma:g} errors of Fig. 4 at every "
                         f"Ny")
                checks[check] = bool(max(pulls[name]) < args.lappi_n_sigma)
                deviations[check] = f"{max(pulls[name]):.2f} errors"
        if coordinate:
            check = (f"{label}: the coordinate-space Q_s is within "
                     f"{args.lappi_coord_tolerance:g} of Fig. 7")
            checks[check] = bool(max(coordinate)
                                 < args.lappi_coord_tolerance)
            deviations[check] = f"{max(coordinate):.4f}"
    summary["checks"] = checks
    summary["largest_deviations"] = deviations
    summary["passed"] = bool(all(checks.values()) and not failures)
    print()
    for name, ok in checks.items():
        print(f"{'OK  ' if ok else 'FAIL'} {name} (largest deviation "
              f"{deviations[name]})")
    plot_lappi(summary["groups"],
               os.path.join(settings["datadir"], "dipole_lappi.pdf"))
    path = os.path.join(settings["datadir"], "dipole_lappi_summary.json")
    with open(path, "w") as f:
        json.dump(summary, f, indent=1)
    print(f"\n{'PASSED' if summary['passed'] else 'FAILED'}; summary in "
          f"{path}")
    return 0 if summary["passed"] else 1


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
    parser.add_argument("--uvdamp", type=float, nargs="*", default=[],
                        help="UVDamp values [1/GeV] for an extra Ny scan at "
                             "the middle g^2 mu and m (none by default: "
                             "UVDamp 0 is the MV model used in practice)")
    parser.add_argument("--seeds", type=int, default=4,
                        help="runs per setting and Ny (two configurations "
                             "each)")
    parser.add_argument("--first-seed", type=int, default=1)
    parser.add_argument("--fit-from", type=float, default=2.,
                        help="the 1/Ny extrapolation uses Ny >= FIT_FROM "
                             "G(0)")
    parser.add_argument("--ny-max", type=int, default=3200,
                        help="largest Ny to which the scan of a setting is "
                             "extended")
    parser.add_argument("--threads", type=int, default=1,
                        help="OpenMP threads per run")
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--quick", action="store_true",
                        help="small test: size 64, L 8, g2mu 1, m 0.2, "
                             "Ny 1..50, 2 seeds")
    parser.add_argument("--lappi", action="store_true",
                        help="the setup of arXiv:0711.3039 instead of the "
                             "scan: size 200, g2mu 1, g2mu L = 100, m = 0 "
                             "and 0.1 g2mu, Ny 1..150, 32 seeds")
    parser.add_argument("--lappi-n-sigma", type=float, default=3.,
                        help="with --lappi: the largest deviation of Q_s "
                             "from Fig. 4 of arXiv:0711.3039, in combined "
                             "errors")
    parser.add_argument("--lappi-coord-tolerance", type=float, default=0.03,
                        help="with --lappi: the same for the "
                             "coordinate-space Q_s and its Fig. 7")
    parser.add_argument("--analyze-only", action="store_true",
                        help="only analyze the results in --datadir")
    parser.add_argument("--keep-files", action="store_true",
                        help="keep the Wilson-line files")
    args = parser.parse_args()
    if args.quick:
        args.size, args.L, args.g2mu, args.m = 64, 8., [1.], [0.2]
        args.ny, args.seeds = [1, 2, 5, 10, 25, 50], 2
    if args.lappi:
        # g^2 mu a = 0.5 and g^2 mu L = 100; only these combinations enter
        args.size, args.g2mu = 200, [1.]
        args.L = 100. * HBARC / args.g2mu[0]
        args.m = [f * args.g2mu[0] for f in sorted(LAPPI_FIG4_ADJ)]
        args.ny = [1, 2, 4, 10, 20, 50, 100, 150]
        if args.seeds == parser.get_default("seeds"):
            args.seeds = 32
    if args.seeds > 1000:
        sys.exit("--seeds must be at most 1000 (see seed_of)")

    ipglasma_path = os.path.abspath(args.ipglasma_path)
    settings = {
        "ipglasma": os.path.join(ipglasma_path, args.ipglasma_cmd),
        "qs_table": os.path.join(ipglasma_path, "qs2Adj_vs_Tp_vs_Y_240.in"),
        "template": os.path.abspath(args.template),
        "datadir": os.path.abspath(args.datadir),
        "size": args.size, "L": args.L, "fit_from": args.fit_from,
        "threads": args.threads, "lappi": args.lappi,
        "analyze_only": args.analyze_only, "keep_files": args.keep_files,
    }
    if not args.analyze_only:
        for path in (settings["ipglasma"], settings["qs_table"]):
            if not os.path.exists(path):
                sys.exit(f"not found: {path} (see --ipglasma-path)")
    os.makedirs(settings["datadir"], exist_ok=True)
    if args.lappi:
        return main_lappi(args, settings)

    m_ref = 0.2 if 0.2 in args.m else args.m[0]
    g2mu_ref = sorted(args.g2mu)[len(args.g2mu) // 2]
    keys = [(g, m, 0.) for g, m in itertools.product(args.g2mu, args.m)]
    keys += [(g2mu_ref, m_ref, uv) for uv in args.uvdamp if uv > 0.]
    # the Ny of each setting; settings without r_s are skipped
    groups, skipped = {}, []
    for key in keys:
        G, _ = setting_two_point(key, settings)
        g0 = float(G[0, 0])
        exact = profiles(np.exp(-CF * (g0 - G)))["axis"]
        if exact[:args.size // 4 + 1].min() > np.exp(-0.5):
            skipped.append({"g2mu_GeV": key[0], "m_GeV": key[1],
                            "UVDamp_invGeV": key[2], "G0": g0,
                            "reason": "D_G > exp(-1/2) for r <= N/4"})
            continue
        groups[key] = setting_ny(args.ny, g0, args.fit_from, args.ny_max)
    configs = [{"g2mu": g, "m": m, "uvdamp": uv, "ny": ny, "g": 1.,
                "seed": seed_of(ny, i, args.first_seed)}
               for (g, m, uv), nys in groups.items()
               for ny, i in itertools.product(nys, range(args.seeds))]
    # g must drop out: the same run with g = 2
    ny_ref = sorted(args.ny)[-1]
    g_check = {"g2mu": g2mu_ref, "m": m_ref, "uvdamp": 0., "ny": ny_ref,
               "g": 2., "seed": seed_of(ny_ref, 0, args.first_seed)}
    configs.append(g_check)
    if (g2mu_ref, m_ref, 0.) not in groups:
        # its setting is skipped: run the g = 1 reference on its own
        configs.append(dict(g_check, g=1.))
    for s in skipped:
        print(f"skipped: g2mu = {s['g2mu_GeV']:g} GeV, m = {s['m_GeV']:g} "
              f"GeV, UVDamp = {s['UVDamp_invGeV']:g}/GeV (G(0) = "
              f"{s['G0']:.3f}): {s['reason']}")

    results, failures = run_all(configs, settings, args.max_workers)

    summary = {"settings": {k: v for k, v in vars(args).items()},
               "groups": [], "skipped": skipped, "failures": failures}
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

    reference = results.get(run_tag(dict(g_check, g=1.), settings))
    other = results.get(run_tag(g_check, settings))
    if reference and other:
        diff = max(float(np.abs(reference[1][k] - other[1][k]).max())
                   for k in ("axis0", "axis1", "diag0", "diag1"))
        summary["g_check"] = {"max_abs_difference": diff,
                              "passed": bool(diff < 1e-10)}
    else:
        summary["g_check"] = {"max_abs_difference": None, "passed": None,
                              "reason": "the runs with g = 1 and 2 failed"}

    # report
    print("\nDipole correlator against its exact lattice prediction:")
    passed = not failures
    for group in summary["groups"]:
        fit = group["fit"]
        print(f"\n g2mu = {group['g2mu_GeV']:g} GeV, m = {group['m_GeV']:g} "
              f"GeV, UVDamp = {group['UVDamp_invGeV']:g}/GeV "
              f"(a = {group['a_fm']:.4f} fm, G(0) = {group['G0']:.3g}, "
              f"r* = {group['r_star']} a, "
              f"Qs/g2mu exact = {group['Qs_exact_over_g2mu']:.4f})")
        if not group["resolved"]:
            print(f"  D - D_G at r* is not resolved at the smallest Ny "
                  f"({group['significance_rstar']:.1f} standard errors); "
                  f"the checks are made at r = {group['r_check']} a, where "
                  "it is most significant")
        rc = group["r_check"]
        print(f"    Ny  max|D-D_G|  D-D_G at r={rc:<3d}  Ny(D-D_G)   Qs/g2mu")
        for row in group["rows"]:
            mark = "*" if row["ny"] in fit["ny"] else " "
            print(f" {mark}{row['ny']:4d}  {row['axis']['max_abs_dev']:10.2e}  "
                  f"{row['dev_rcheck']:+.2e}({row['error_rcheck']:.0e})  "
                  f"{row['ny'] * row['dev_rcheck']:+9.4f}  "
                  f"{row['Qs_over_g2mu']:.4f}")
        if "d_inf" in fit:
            print(f"  extrapolation (* : Ny >= {fit['ny_from']:.3g}): "
                  f"D - D_G -> {fit['d_inf']:+.2e} +- {fit['error']:.1e}, "
                  f"c = {fit['c']:+.3f}")
        else:
            print(f"  no extrapolation: fewer than three Ny >= "
                  f"{fit['ny_from']:.3g} (see --ny-max)")
        for name, ok in group["checks"].items():
            status = "INCOMPLETE" if ok is None else ("OK  " if ok
                                                      else "FAIL")
            print(f"  {status} {name}")
            passed &= ok is True
    g = summary["g_check"]
    if g["passed"] is None:
        print(f"\n INCOMPLETE g drops out: {g['reason']}")
    else:
        print(f"\n {'OK  ' if g['passed'] else 'FAIL'} g drops out: "
              f"max |D(g=2) - D(g=1)| = {g['max_abs_difference']:.1e}")
    passed &= g["passed"] is True
    summary["passed"] = bool(passed)
    with open(os.path.join(settings["datadir"], "dipole_summary.json"),
              "w") as f:
        json.dump(summary, f, indent=1)
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{settings['datadir']}/dipole_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

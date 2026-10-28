# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
The BK equation with the square-root running coupling of Lappi and
Mantysaari (arXiv:1212.4825), the mean-field limit of the JIMWLK evolution
with the kernel sqrt(alpha_s(r)) K(r) (Section 10 of the paper on
IP-Glasma 2.0), for a translation-invariant dipole amplitude N(r):

  dS(r)/dy = Nc/(2 pi^2) int d^2z Kt(X, Y) [S(X) S(Y) - S(r)],
  Kt = a(X)/X^2 + a(Y)/Y^2 - 2 sqrt(a(X) a(Y)) X.Y/(X^2 Y^2),

with S = 1 - N, X = x - z, Y = y - z and r = x - y. All lengths are in units
of 1/Lambda_QCD. The coupling is Eq. (40) of arXiv:1212.4825, with the scale
4 e^{-2 gamma_E}/r^2:

  alpha_s(r) = 4 pi / (beta c ln[(mu0^2/Lambda^2)^(1/c)
                                  + (4 e^{-2 gamma_E}/(r^2 Lambda^2))^(1/c)]),

beta = 11 - 2 N_f/3. The integrand is symmetric under X <-> Y, so z is
restricted to X < Y (twice that region); in polar coordinates around x the
angle runs over [theta_min, pi] and its mirror image, so the point Y = 0
never enters. ln N is evolved on a logarithmic grid in r with RK4, the
integral is done with the trapezoidal rule in ln X and Gauss-Legendre in the
angle, and N(Y) is interpolated linearly in ln Y (as N ~ r^2 below the
grid and constant above it).
"""

import numpy as np

NC = 3
QS_LEVEL = np.exp(-0.5)


def alpha_s(r, c, mu0, nf):
    """The coupling of Eq. (40) of arXiv:1212.4825; r in 1/Lambda, mu0 in
    Lambda."""
    beta = 11. - 2. * nf / 3.
    return 4. * np.pi / (beta * c * np.log(
        mu0 ** (2. / c)
        + (4. * np.exp(-2. * np.euler_gamma) / r ** 2) ** (1. / c)))


class BK:
    """The right-hand side of the BK equation for ln N on the grid lnr."""

    def __init__(self, lnr, c, mu0, nf, ntheta=32):
        self.lnr = lnr
        r = np.exp(lnr)
        dln = lnr[1] - lnr[0]
        wx = np.full(len(r), dln)
        wx[[0, -1]] *= 0.5
        t, wt = np.polynomial.legendre.leggauss(ntheta)
        R = r[:, None, None]
        X = r[None, :, None]
        # X < Y  <=>  cos(theta) > r/(2 X)
        tmin = np.arccos(np.clip(R / (2. * X), -1., 1.))
        theta = tmin + (np.pi - tmin) * (t[None, None, :] + 1.) / 2.
        wtheta = (np.pi - tmin) / 2. * wt[None, None, :]
        cos = np.cos(theta)
        Y2 = X ** 2 + R ** 2 - 2. * X * R * cos
        aX, aY = alpha_s(X, c, mu0, nf), alpha_s(np.sqrt(Y2), c, mu0, nf)
        kernel = (aX / X ** 2 + aY / Y2
                  - 2. * np.sqrt(aX * aY) * (X ** 2 - X * R * cos)
                  / (X ** 2 * Y2))
        # factor 2 for X <-> Y and 2 for the mirror image of the angle
        self.w = NC / (2. * np.pi ** 2) * 4. * kernel * X ** 2 \
            * wx[None, :, None] * wtheta
        u = (0.5 * np.log(Y2) - lnr[0]) / dln
        self.i = np.clip(np.floor(u).astype(np.int64), 0, len(lnr) - 2)
        self.f = u - self.i
        self.below = u < 0.
        self.above = u > len(lnr) - 1
        self.ratio2 = np.exp(2. * u[self.below] * dln)      # (Y/r_min)^2

    def rhs(self, lnN):
        """d ln N/dy."""
        N = np.exp(lnN)
        S = 1. - N
        NY = np.exp(lnN[self.i] * (1. - self.f) + lnN[self.i + 1] * self.f)
        NY[self.below] = N[0] * self.ratio2
        NY[self.above] = N[-1]
        integrand = S[None, :, None] * (1. - NY) - S[:, None, None]
        dS = np.einsum("ijk,ijk->i", self.w, integrand)
        return -dS / N


def saturation_scale(lnr, lnN):
    """Q_s = sqrt(2)/r_s with N(r_s) = 1 - e^{-1/2}, from a cubic in
    ln r(S) around the crossing."""
    S = 1. - np.exp(lnN)
    k = int(np.argmax(S < QS_LEVEL))
    j = slice(max(k - 2, 0), k + 2)
    return np.sqrt(2.) / np.exp(np.polyval(
        np.polyfit(S[j] - QS_LEVEL, lnr[j], 3), 0.))


def initial_condition(r, N, lnr):
    """N on the grid lnr from N at the points r (e.g. the lattice dipole
    correlator, N = 1 - D), extended as N ~ r^2 to smaller and as a
    constant to larger r."""
    keep = (r > 0.) & (N > 0.)
    r, N = r[keep], N[keep]
    lnN = np.interp(lnr, np.log(r), np.log(N))
    small = lnr < np.log(r[0])
    lnN[small] = np.log(N[0]) + 2. * (lnr[small] - np.log(r[0]))
    lnN[lnr > np.log(r[-1])] = np.log(N[-1])
    return lnN


def solve(r, N, ymax, c, mu0, nf, dy=0.02, nr=400, rmin=1e-7, rmax=1e2,
          ntheta=32):
    """Evolve N(r) (r in 1/Lambda) to ymax; returns y and Q_s/Lambda at
    every step."""
    lnr = np.linspace(np.log(rmin), np.log(rmax), nr)
    bk = BK(lnr, c, mu0, nf, ntheta)
    lnN = initial_condition(np.asarray(r), np.asarray(N), lnr)
    steps = int(round(ymax / dy))
    ys, qs = [0.], [saturation_scale(lnr, lnN)]
    for n in range(1, steps + 1):
        k1 = bk.rhs(lnN)
        k2 = bk.rhs(lnN + dy / 2. * k1)
        k3 = bk.rhs(lnN + dy / 2. * k2)
        k4 = bk.rhs(lnN + dy * k3)
        lnN = np.minimum(lnN + dy / 6. * (k1 + 2. * k2 + 2. * k3 + k4), 0.)
        ys.append(n * dy)
        qs.append(saturation_scale(lnr, lnN))
    return np.array(ys), np.array(qs)

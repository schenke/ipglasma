# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
The gluon spectrum of two MV nuclei to leading order in g^2 mu, in the
conventions of IP-Glasma, for the comparison of the dilute limit.

Color charges <g rho^a(x) g rho^b(y)> = (g^2 mu)^2 delta^ab delta^2(x - y)
(summed over the Ny layers), (-nabla^2 + m^2) Lambda = g rho and
V = exp(-i Lambda^a t^a), so to first order A^i = -(1/g) d_i Lambda for
each nucleus, and <Lambda^a(p) Lambda^b(-p)> = (g^2 mu)^2/(p^2 + m^2)^2 per
area. The fields at tau = 0+ (A^eta = (ig/2)[A^i_A, A^i_B]) are, to order
rho_A rho_B, the longitudinal fields

  E^eta_c = (1/g) f^abc d_i Lambda_A^a d_i Lambda_B^b,
  B^eta_c = (1/g) f^abc eps_ij d_i Lambda_A^a d_j Lambda_B^b,

and, since the nonlinear terms are of higher order, the field strengths
then evolve as free boost-invariant fields: the E^eta mode as
E_z = c J0(k tau) with a transverse magnetic field c J1(k tau), the B^eta
mode as B_z = d J0(k tau) with a transverse electric field d J1(k tau).
Averaged over the color charges (sum_abc f^abc f^abc = Nc (Nc^2 - 1)),
per transverse area and with q = k - p,

  <|E^eta(k)|^2> = Nc (Nc^2 - 1)/g^2 (g^2 mu_A)^2 (g^2 mu_B)^2 I_E(k),
  I_E(k) = int d^2p/(2 pi)^2 (p.q)^2 / ((p^2 + m^2)^2 (q^2 + m^2)^2),

and likewise <|B^eta(k)|^2> with (p x q)^2, I_B. IP-Glasma counts gluons
with the electric fields in Coulomb gauge, n = 2 |E|^2/omega, which at the
time tau of the measurement gives

  dN/dy d^2k = Nc (Nc^2 - 1)/g^2 (g^2 mu_A)^2 (g^2 mu_B)^2 / (2 pi)^2
               (tau/k) [I_E(k) J0(k tau)^2 + I_B(k) J1(k tau)^2]

per area, and at late times, with J0^2, J1^2 -> 1/(pi k tau) on average,

  dN/dy d^2k = Nc (Nc^2 - 1)/g^2 (g^2 mu_A)^2 (g^2 mu_B)^2
               (I_E + I_B)/(4 pi^3 k^2),

the classical result of Kovner, McLerran and Weigert and of Gyulassy and
McLerran, which falls like ln(k/m)/k^4 for k >> m. Lattice artifacts (the
lattice momenta of the kernel, the lattice E^eta and the dispersion
omega = k~ of the estimator) are not included, so the comparison holds
for k a << 1.
"""

import numpy as np
from scipy.special import j0, j1

HBARC = 0.1973269804  # GeV fm
# the k_T bins of GluonMultiplicity: width 2.83/100 in lattice units
BIN_WIDTH = 2.83 / 100
BINS = 100


def integrals(k, m, n_p=1500, n_theta=300):
    """I_E(k) and I_B(k) [GeV^-4 per (2 pi)^2 ... see above], for k and m in
    GeV: the integral over d^2p in polar coordinates around p = 0, on a
    logarithmic grid in p and Gauss-Legendre in the angle."""
    lnp = np.linspace(np.log(1e-5 * m), np.log(1e4 * m), n_p)
    p = np.exp(lnp)
    t, wt = np.polynomial.legendre.leggauss(n_theta)
    theta = 0.5 * np.pi * (t + 1.)          # [0, pi], twice for [0, 2 pi]
    wtheta = np.pi * wt
    P, TH = np.meshgrid(p, theta, indexing="ij")
    q2 = k * k + P * P - 2. * k * P * np.cos(TH)
    denominator = (P * P + m * m) ** 2 * (q2 + m * m) ** 2
    p_dot_q = P * k * np.cos(TH) - P * P
    p_cross_q = P * k * np.sin(TH)
    out = []
    for numerator in (p_dot_q ** 2, p_cross_q ** 2):
        inner = (numerator / denominator * wtheta[None, :]).sum(axis=1) \
            * p * p                          # d^2p = p^2 dln p dtheta
        out.append(float(np.trapz(inner, lnp)) / (2. * np.pi) ** 2)
    return out


def spectrum(k, ie, ib, tau=None, nc=3, g=1.):
    """dN/dy d^2k per area and per (g^2 mu_A)^2 (g^2 mu_B)^2, for k [GeV]
    and the integrals at k: at the time tau [GeV^-1] for the electric
    estimator of IP-Glasma, or at late times for tau=None."""
    k = np.asarray(k, dtype=float)
    color = nc * (nc * nc - 1) / g ** 2
    if tau is None:
        return color * (ie + ib) / (4. * np.pi ** 3 * k * k)
    return color / (2. * np.pi) ** 2 * (tau / k) \
        * (ie * j0(k * tau) ** 2 + ib * j1(k * tau) ** 2)


def lattice_prediction(size, L, m, tau, nc=3, g=1., grid_points=120):
    """The prediction for one lattice, per (g^2 mu_A)^2 (g^2 mu_B)^2: the
    k_T bins of GluonMultiplicity [GeV], dN/dy d^2k of the whole lattice
    [GeV^-2] at the time tau [fm] and at late times, binned like the code
    (the modes of the periodic lattice binned in the lattice momentum k~,
    the modes on the lines k_x = -pi/a or k_y = -pi/a left out, and the
    average over the modes of a bin), and dN/dy per area [fm^-2] at tau,
    the sum over the same modes without the zero mode."""
    a = L / size
    area = (L / HBARC) ** 2                  # GeV^-2
    index = np.arange(size)
    kx = 2. * np.pi * (-0.5 + index / size)
    KX, KY = np.meshgrid(kx, kx, indexing="ij")
    used = (index[:, None] != 0) & (index[None, :] != 0)
    ktilde = 2. * np.sqrt(np.sin(KX / 2.) ** 2 + np.sin(KY / 2.) ** 2)[used]
    k = np.sqrt(KX ** 2 + KY ** 2)[used] * HBARC / a
    nonzero = ktilde > 0.
    ktilde, k = ktilde[nonzero], k[nonzero]
    grid = np.geomspace(0.5 * k.min(), 1.01 * k.max(), grid_points)
    ie, ib = np.array([integrals(kk, m) for kk in grid]).T
    at = {}
    for name, t in (("tau", tau / HBARC), ("late", None)):
        values = spectrum(grid, ie, ib, t, nc, g) * area
        if t is None:
            at[name] = np.exp(np.interp(np.log(k), np.log(grid),
                                        np.log(values)))
        else:
            # the Bessel functions oscillate: interpolate the integrals
            iek = np.exp(np.interp(np.log(k), np.log(grid), np.log(ie)))
            ibk = np.exp(np.interp(np.log(k), np.log(grid), np.log(ib)))
            at[name] = spectrum(k, iek, ibk, t, nc, g) * area
    kt = (np.arange(BINS) + 0.5) * BIN_WIDTH / a * HBARC
    binned = {name: np.full(BINS, np.nan) for name in at}
    for ik in range(BINS):
        inside = (ktilde > ik * BIN_WIDTH) & (ktilde <= (ik + 1) * BIN_WIDTH)
        if inside.any():
            for name in at:
                binned[name][ik] = at[name][inside].mean()
    # each mode is a cell (2 pi/L)^2 of k space
    dn_dy = float(at["tau"].sum() * (2. * np.pi / (L / HBARC)) ** 2
                  / area / HBARC ** 2)   # fm^-2
    return kt, binned["tau"], binned["late"], dn_dy

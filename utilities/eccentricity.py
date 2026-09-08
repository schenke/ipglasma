#!/usr/bin/env python3
"""
Energy-density-weighted spatial eccentricities epsilon_n, computed
directly from a 2D energy density grid (e.g. the T00 component of an
IP-Glasma Tmunu snapshot, see read_tmunu.py).

    eps_n * exp(i n Psi_n) = -<r^n exp(i n phi)> / <r^n>

with <O> = sum_cells e(x,y) O(x,y) / sum_cells e(x,y), and (x, y) measured
relative to the energy-density-weighted centroid of the event. This is the
standard participant-plane eccentricity definition used throughout the
heavy-ion literature (e.g. Alver & Roland, PRC 81, 054905 (2010)); the
overall sign convention does not matter for the magnitude eps_n.
"""

import numpy as np


def grid_coordinates(nx, ny, dx, dy):
    """Return (X, Y) meshgrids of shape (ny, nx) with physical spacing
    dx, dy. The absolute origin is arbitrary -- only used differences
    between grid points matter for eps_n, since the centroid is
    subtracted before computing eps_n."""
    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    X, Y = np.meshgrid(x, y)  # shapes (ny, nx), matches energy_density
    return X, Y


def eccentricity(energy_density, dx, dy, n=2, clip_negative=True):
    """Compute the energy-density-weighted eccentricity eps_n and
    participant-plane angle Psi_n for a 2D energy density grid.

    Parameters
    ----------
    energy_density : 2D array, shape (ny, nx)
    dx, dy : float
        Lattice spacing in the x and y directions (physical units, e.g.
        fm). Only the ratio/product of dx, dy to r^n matters; an overall
        unit choice cancels in eps_n since it appears in both numerator
        and denominator.
    n : int
        Harmonic order (2 for the standard ellipticity eps_2).
    clip_negative : bool
        Clip tiny negative energy-density values (numerical noise) to
        zero before using them as weights.

    Returns
    -------
    eps_n : float
    psi_n : float
        Participant-plane angle in radians.
    """
    e = np.asarray(energy_density, dtype=np.float64)
    if clip_negative:
        e = np.clip(e, 0.0, None)

    ny, nx = e.shape
    X, Y = grid_coordinates(nx, ny, dx, dy)

    total = e.sum()
    if total <= 0:
        return float("nan"), float("nan")

    xbar = (e * X).sum() / total
    ybar = (e * Y).sum() / total
    xr = X - xbar
    yr = Y - ybar

    r_n = (xr**2 + yr**2) ** (n / 2.0)
    phi = np.arctan2(yr, xr)

    num_c = (e * r_n * np.cos(n * phi)).sum()
    num_s = (e * r_n * np.sin(n * phi)).sum()
    denom = (e * r_n).sum()

    eps_n = np.sqrt(num_c**2 + num_s**2) / denom
    psi_n = np.arctan2(num_s, num_c) / n + np.pi / n
    return float(eps_n), float(psi_n)


def eccentricity_from_tmunu_file(path, n=2, clip_negative=True):
    """Load a Tmunu snapshot file (.ipgt or .dat, see read_tmunu.py) and
    return (eps_n, psi_n, tau_fm)."""
    from read_tmunu import get_energy_density  # local import, avoids a
    # hard dependency on read_tmunu.py for callers that only need the
    # eccentricity() function above.

    energy_density, dx, dy, tau_fm = get_energy_density(path)
    eps_n, psi_n = eccentricity(
        energy_density, dx, dy, n=n, clip_negative=clip_negative)
    return eps_n, psi_n, tau_fm

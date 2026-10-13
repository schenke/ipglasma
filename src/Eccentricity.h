// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#ifndef SRC_ECCENTRICITY_H_
#define SRC_ECCENTRICITY_H_

#include "Lattice.h"
#include "Parameters.h"

/**
 * The spatial eccentricities and the energy-momentum-tensor anisotropy of
 * the evolved fields.
 */
namespace Eccentricity {

/**
 * The participant-plane angle \f$\Psi_n\f$ of the moments
 * \f$\langle r^m \sin n\phi\rangle\f$ and
 * \f$\langle r^m \cos n\phi\rangle\f$, defined by
 * \f$\varepsilon_n e^{in\Psi_n} = -\langle r^m e^{in\phi}\rangle /
 * \langle r^m\rangle\f$, i.e. \f$\Psi_n = [\operatorname{atan2}(
 * \langle r^m \sin n\phi\rangle, \langle r^m \cos n\phi\rangle)
 * + \pi]/n\f$, which holds for either sign of the cosine moment.
 * \param[in] n Harmonic, at least 1.
 * \param[in] sinMoment \f$\langle r^m \sin n\phi\rangle\f$ (any
 * normalization).
 * \param[in] cosMoment \f$\langle r^m \cos n\phi\rangle\f$ (the same
 * normalization).
 * \return \f$\Psi_n \in [0, 2\pi/n)\f$ [rad].
 */
double participantPlaneAngle(int n, double sinMoment, double cosMoment);

/**
 * Computes the spatial eccentricities \f$\varepsilon_1,\ldots,
 * \varepsilon_6\f$ and their participant-plane angles
 * \f$\Psi_n \in [0, 2\pi/n)\f$, defined by
 * \f$\varepsilon_n e^{in\Psi_n} = -\langle r^m e^{in\phi}\rangle /
 * \langle r^m\rangle\f$ (\f$m = 3\f$ for \f$n = 1\f$, \f$m = n\f$
 * otherwise), from the energy-density-weighted
 * (running-coupling-corrected, `getEpsilon()*getutau()`-weighted,
 * cells below \p cutoff excluded) spatial moments, first recentering
 * on the energy-weighted centroid. If \p doAniso is `0`, appends one
 * row to `eccentricities<id>.dat`. If \p doAniso is `1`, instead
 * appends to `anisotropy<id>.dat` the \f$T^{xx}-T^{yy}\f$ spatial
 * anisotropy (via the anonymous-namespace `computeRotatedAnisotropy`
 * helper) resampled at ten angles around the flow-velocity event
 * plane \f$\Psi_U\f$ (computed from `getux()`/`getuy()`); the program only
 * calls it with \p doAniso `0`, the anisotropy is a debugging helper.
 * The file is described in \ref md_OUTPUT "OUTPUT.md".
 * \param[in] lat Lattice to read the energy density (and, for
 * \p doAniso==1, `Txx`/`Txy`/`Tyy`/flow velocity) from.
 * \param[in] param Simulation parameters.
 * \param[in] it Current time step index, used for the output file's
 * time column and (on `it==1`) to store \f$\Psi_2\f$ in
 * `param->event.psi`.
 * \param[in] cutoff Energy-density cutoff [GeV/fm\f$^3\f$] below which a
 * cell is excluded from the weighted averages (`eccentricityCutoff`).
 * \param[in] doAniso Selects which output file/quantity is written
 * (`0`: eccentricities, `1`: rotated-tensor anisotropy).
 */
void compute(
    Lattice *lat, Parameters *param, int it, double cutoff, int doAniso);

}  // namespace Eccentricity

#endif  // SRC_ECCENTRICITY_H_

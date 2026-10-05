// Eccentricity.h is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.

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
 * Computes the spatial eccentricities \f$\varepsilon_1,\ldots,
 * \varepsilon_6\f$ and their event-plane angles
 * \f$\Psi_1,\ldots,\Psi_6\f$ from the energy-density-weighted
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
 * time column and (on `it==1`) to store \f$\Psi\f$ in
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

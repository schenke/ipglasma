// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#ifndef SRC_ENERGYMOMENTUMTENSOR_H_
#define SRC_ENERGYMOMENTUMTENSOR_H_

#include "Lattice.h"
#include "Parameters.h"

/**
 * The energy-momentum tensor \f$T^{\mu\nu}\f$ of the classical
 * Yang-Mills fields on the lattice.
 */
namespace EnergyMomentumTensor {

/**
 * Computes every component of the energy-momentum tensor
 * \f$T^{\mu\nu}\f$ at every cell (the spatial plaquette, the diagonal
 * electric and magnetic/gradient contributions, and the six
 * off-diagonal components), storing the results into `lat->cells`. The
 * value of cell \f$(i, j)\f$ belongs to the center of its plaquette,
 * \f$x + (a/2, a/2)\f$, where the fields are averaged to; the output
 * writers (MyEigen) interpolate it from there. The
 * outermost ring of cells is a guard region and gets \f$T^{\mu\nu}=0\f$.
 * \f$T^{\tau x}\f$ and \f$T^{x\eta}\f$ need the fields two sites ahead in
 * \f$x\f$, and \f$T^{\tau y}\f$ and \f$T^{y\eta}\f$ two sites ahead in
 * \f$y\f$, so they are also 0 on the last interior column
 * (\f$i = N-2\f$) and row (\f$j = N-2\f$), respectively.
 * Runs its team kernels in sequence inside one shared
 * `#pragma omp parallel` region.
 * \param[in,out] lat Lattice to read the fields from and write
 * \f$T^{\mu\nu}\f$ into (via `lat->cells`); \c lat->Uy1 is overwritten
 * with the spatial plaquettes.
 * \param[in] param Simulation parameters.
 * \param[in] it Current time step index (\f$\tau=it\cdot d\tau\f$),
 * used to convert the momentum fields' lattice normalization to
 * physical units.
 */
void compute(Lattice *lat, Parameters *param, int it);

}  // namespace EnergyMomentumTensor

#endif  // SRC_ENERGYMOMENTUMTENSOR_H_

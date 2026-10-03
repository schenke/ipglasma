// EnergyMomentumTensor.h is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.

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
 * outermost ring of cells is a guard region and gets \f$T^{\mu\nu}=0\f$.
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

// Evolution.h is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_EVOLUTION_H_
#define SRC_EVOLUTION_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include "GluonMultiplicity.h"
#include "Group.h"
#include "Lattice.h"
#include "Parameters.h"
#include "PrettyOstream.h"

/**
 * Advances the classical Yang-Mills fields built by Init via a leapfrog
 * time evolution: the coordinate fields (\f$\phi\f$, the transverse
 * gauge links `Ux`/`Uy`) live at integer-\f$\tau\f$ steps and are
 * updated by evolvePhi()/evolveU(), while the momentum fields (the
 * electric fields `U`/`U2`, and \f$\pi\f$) live at half-integer steps
 * and are updated by evolvePi()/evolveE(). Also computes the
 * energy-momentum tensor (EnergyMomentumTensor::compute()), derived observables
 * (Eccentricity::compute(), u()), and an optional final-time gluon
 * spectrum/multiplicity estimate (GluonMultiplicity::compute()) via
 * Coulomb-gauge fixing and FFT.
 */
class Evolution {
  private:
    /// Gluon spectrum/multiplicity measurement at the final time (owns
    /// the FFT used for it).
    GluonMultiplicity gluonMultiplicity_;
    /// Log sink for progress/warning/error messages.
    PrettyOstream messager_;

  public:
    /**
     * Constructs an Evolution for a lattice of the given dimensions.
     * \param[in] nn Two-element array `{N, N}`, the transverse lattice
     * dimensions, forwarded to the gluon-multiplicity FFT.
     */
    explicit Evolution(const int nn[]) : gluonMultiplicity_(nn) {}

    /**
     * Top-level evolution driver: after an initial coordinate half-step
     * (evolvePhi()/evolveU() from \f$\tau=0\f$ to \f$\tau=d\tau\f$,
     * since the momenta at \f$\tau=d\tau/2\f$ are taken equal to their
     * \f$\tau=0\f$ values), leapfrog-advances the fields
     * (evolveStepPersistent(), an anonymous-namespace helper in
     * Evolution.cpp that runs one shared `#pragma omp parallel` team
     * through evolvePi()/evolveE()/evolvePhi()/evolveU() in order) up to
     * `param->evolution.maxTime` (or, if `evolution.inverseQsForMaxTime`, up to
     * \f$1/Q_s\f$). At the final time step (and, if `output.writeOutputs
     * == 5`, at four additional fixed intermediate times), temporarily
     * recenters the momenta from \f$\tau_{n-1/2}\f$ to \f$\tau_n\f$ to
     * measure EnergyMomentumTensor::compute() and either u() or
     * `MyEigen::writeTmunu4D()` (depending on `output.writeEpsilonUHydro`),
     * then restores the unmodified momenta so the measurement cannot perturb
     * the trajectory. At the very end, runs checkGaussLaw(), and -- if
     * `output.computeGluonMultiplicity` -- Eccentricity::compute() and
     * GluonMultiplicity::compute(), stopping early if it reports no
     * collision.
     * \param[in,out] lat Lattice to evolve in place.
     * \param[in] group Group instance, forwarded to
     * GluonMultiplicity::compute().
     * \param[in,out] param Simulation parameters.
     */
    void run(Lattice *lat, Group *group, Parameters *param);
    /**
     * Leapfrog coordinate update for the transverse gauge links: rotates
     * `Ux`/`Uy` by \f$\exp(i g^2 d\tau/(\tau+d\tau/2)\,U)\f$ (a
     * second-order Padé approximant of the exponential), using the
     * electric fields `U`/`U2` as the generator. Runs
     * `evolveUTeam` (an anonymous-namespace helper in Evolution.cpp)
     * inside its own `#pragma omp parallel` region.
     * \param[in,out] lat Lattice whose `Ux`/`Uy` are updated in place.
     * \param[in] param Simulation parameters.
     * \param[in] dtau Time step [lattice units].
     * \param[in] tau Current proper time [lattice units], evaluated at
     * the half-step \f$\tau+d\tau/2\f$ inside the update.
     */
    void evolveU(Lattice *lat, Parameters *param, double dtau, double tau);
    /**
     * Leapfrog coordinate update for the scalar field: \f$\phi \mathrel{+}=
     * (\tau+d\tau/2)\,d\tau\,\pi\f$. Runs `evolvePhiTeam` (an
     * anonymous-namespace helper in Evolution.cpp) inside its own
     * `#pragma omp parallel` region.
     * \param[in,out] lat Lattice whose \c Uy2 (\f$\phi\f$) is updated in
     * place.
     * \param[in] param Simulation parameters.
     * \param[in] dtau Time step [lattice units].
     * \param[in] tau Current proper time [lattice units].
     */
    void evolvePhi(Lattice *lat, Parameters *param, double dtau, double tau);
    /**
     * Leapfrog momentum update for \f$\pi\f$: adds
     * \f$(d\tau/\tau)\f$ times the covariant discrete Laplacian of
     * \f$\phi\f$ (parallel-transported via `Ux`/`Uy` to its four
     * neighbors). Runs `evolvePiTeam` (an anonymous-namespace helper in
     * Evolution.cpp) inside its own `#pragma omp parallel` region.
     * \param[in,out] lat Lattice whose \c Ux2 (\f$\pi\f$) is updated in
     * place.
     * \param[in] param Simulation parameters.
     * \param[in] dtau Time step [lattice units].
     * \param[in] tau Current proper time [lattice units].
     */
    void evolvePi(Lattice *lat, Parameters *param, double dtau, double tau);
    /**
     * Leapfrog momentum update for the electric fields `U`/`U2`: adds
     * the traceless anti-Hermitian plaquette force (from the four
     * spatial plaquettes touching each link) and the \f$\phi\f$-\f$\pi\f$
     * commutator force. Runs `evolveETeam` (an anonymous-namespace
     * helper in Evolution.cpp, using the shared `addEForceSU3` kernel)
     * inside its own `#pragma omp parallel` region.
     * \param[in,out] lat Lattice whose `U`/`U2` are updated in place.
     * \param[in] param Simulation parameters.
     * \param[in] dtau Time step [lattice units].
     * \param[in] tau Current proper time [lattice units].
     */
    void evolveE(Lattice *lat, Parameters *param, double dtau, double tau);
    /**
     * Diagnostic: computes, at every cell, the discrete Gauss-law
     * constraint violation
     * \f$U_x(x-\hat x)^\dagger E_1(x-\hat x) U_x(x-\hat x) - E_1(x) +
     * U_y(x-\hat y)^\dagger E_2(x-\hat y) U_y(x-\hat y) - E_2(x) -
     * i[\phi,\pi]\f$ (evaluated using the momenta at their current
     * \f$\tau-d\tau/2\f$), and logs the largest value found.
     * \param[in] lat Lattice to read `U`/`U2`/`Ux`/`Uy`/`Uy2`/\c
     * Ux2 from.
     * \param[in] param Simulation parameters.
     */
    void checkGaussLaw(Lattice *lat, Parameters *param);
    /**
     * Writes a binary snapshot (`evolvedFields<id>_it<it>.ipgf`) of the
     * six evolved fields (\f$\phi\f$, \f$\pi\f$, \c U (\f$E_1\f$), \c U2
     * (\f$E_2\f$), \c Ux, \c Uy) as full \f$3\times3\f$ complex
     * matrices per cell, little-endian single-precision, preceded by an
     * 8-byte magic string, an 8-byte metadata length, and a JSON
     * metadata header describing the layout/units/staggering.
     * The file is described in \ref md_OUTPUT "OUTPUT.md".
     * \param[in] lat Lattice to read the fields from.
     * \param[in] param Simulation parameters.
     * \param[in] it Current time step index, used for the file name and
     * the metadata's time/tau fields.
     */
    void writeEvolvedFields(Lattice *lat, Parameters *param, int it);
    /**
     * Computes the hydrodynamic flow-velocity field \f$u^\mu\f$ by
     * delegating to `MyEigen::flowVelocity4D`.
     * \param[in,out] lat Lattice to read \f$T^{\mu\nu}\f$ from and write
     * the flow velocity into.
     * \param[in] param Simulation parameters.
     * \param[in] it Current time step index.
     * \param[in] finalFlag Forwarded to `MyEigen::flowVelocity4D`
     * (selects final- vs. intermediate-time output naming/behavior).
     */
    void u(Lattice *lat, Parameters *param, int it, bool finalFlag);
    /**
     * The final-time flow measurement of run(), after
     * EnergyMomentumTensor::compute(): with `output.writeEpsilonUHydro` the
     * full u() solve and hydro output; otherwise only the raw \f$T^{\mu\nu}\f$
     * output, preceded by the flow-velocity solve if
     * `output.computeGluonMultiplicity` (so that Eccentricity::compute(), which
     * weights by \f$\epsilon u^\tau\f$, sees the solved fields).
     * \param[in,out] lat Lattice holding \f$T^{\mu\nu}\f$; receives
     * \f$\epsilon\f$ and \f$u^\mu\f$ when the solve runs.
     * \param[in] param Simulation parameters.
     * \param[in] it Final time step index.
     */
    void finalFlowMeasurement(Lattice *lat, Parameters *param, int it);
};

#endif  // SRC_EVOLUTION_H_

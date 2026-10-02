// Evolution.cpp is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Evolution.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "EnergyMomentumTensor.h"
#include "Fragmentation.h"
#include "GaugeFix.h"
#include "Instrumentation.h"
#include "MyEigen.h"
#include "PhysConst.h"
#include "RunningCoupling.h"
#include "SU3.h"

using Fragmentation::kkp;
using PhysConst::hbarc;
using PhysConst::m_kaon;
using PhysConst::m_pion;
using PhysConst::m_proton;
using PhysConst::Nc;

using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;

//**************************************************************************
// Evolution class.

namespace {

/// Scratch matrices for evolveUTeam(), reused across cells to avoid
/// reallocating.
struct EvolveUScratch {
    /// Initializes \c one to the identity matrix.
    EvolveUScratch() : one(1.) {}

    /// Padé-approximant intermediate for the \c Ux-side rotation.
    Matrix E1;
    /// Padé-approximant intermediate for the \c Uy-side rotation.
    Matrix E2;
    /// General scratch used while building `E1`/`E2`.
    Matrix temp1;
    /// General scratch used while building `E1`/`E2`.
    Matrix temp2;
    /// Reusable identity matrix.
    Matrix one;
};

/// Scratch matrices for evolvePhiTeam(), reused across cells to avoid
/// reallocating.
struct EvolvePhiScratch {
    /// Current cell's \f$\phi\f$ (\c Uy2), updated in place.
    Matrix phi;
    /// Current cell's \f$\pi\f$ (\c Ux2).
    Matrix pi;
};

/// Scratch matrices for evolvePiTeam(), reused across cells to avoid
/// reallocating.
struct EvolvePiScratch {
    /// Current cell's \f$U_x\f$.
    Matrix Ux;
    /// Current cell's \f$U_y\f$.
    Matrix Uy;
    /// \f$U_x\f$ at the \f$-\hat x\f$ neighbor.
    Matrix UxXm1;
    /// \f$U_y\f$ at the \f$-\hat y\f$ neighbor.
    Matrix UyYm1;
    /// Current cell's \f$\phi\f$.
    Matrix phi;
    /// \f$\phi\f$ parallel-transported from the \f$+\hat x\f$ neighbor.
    Matrix phiX;
    /// \f$\phi\f$ parallel-transported from the \f$+\hat y\f$ neighbor.
    Matrix phiY;
    /// \f$\phi\f$ parallel-transported from the \f$-\hat x\f$ neighbor.
    Matrix phimX;
    /// \f$\phi\f$ parallel-transported from the \f$-\hat y\f$ neighbor.
    Matrix phimY;
    /// The covariant discrete Laplacian of \f$\phi\f$,
    /// \f$\phi_X+\phi_{-X}+\phi_Y+\phi_{-Y}-4\phi\f$.
    Matrix bracket;
    /// Current cell's \f$\pi\f$, updated in place.
    Matrix pi;
};

/// Scratch matrices for evolveETeam(), reused across cells to avoid
/// reallocating.
struct EvolveEScratch {
    /// Current cell's \f$U_x\f$.
    Matrix Ux;
    /// Current cell's \f$U_y\f$.
    Matrix Uy;
    /// General scratch used while building `U12`/`U1m2`/`U2m1`.
    Matrix temp1;
    /// General scratch used while building `U12`/`U1m2`/`U2m1`.
    Matrix temp2;
    /// Electric field (\c U or \c U2) being updated, passed to
    /// addEForceSU3().
    Matrix En;
    /// \f$\phi\f$ parallel-transported to the plaquette's far corner.
    Matrix phiN;
    /// Current cell's \f$\phi\f$.
    Matrix phi;
    /// The \f$(x,y)\f$-oriented spatial plaquette touching this link.
    Matrix U12;
    /// The \f$(x,-y)\f$-oriented spatial plaquette touching this link.
    Matrix U1m2;
    /// The \f$(-x,y)\f$-oriented spatial plaquette touching this link.
    Matrix U2m1;
};

/**
 * Adds one electric field's plaquette and \f$[\phi_N,\phi]\f$
 * commutator force to \p En in place: the traceless anti-Hermitian
 * part of \f$a+\text{bSign}\cdot b\f$ (scaled by \p coeffPlaq) plus
 * \f$[\phi_N,\phi]\f$ (scaled by \p coeffComm), then re-projects \p En
 * traceless. Shared by both calls in evolveETeam() (once for \c U with
 * \f$a=U_{12}\f$/\f$b=U_{1m2}\f$/`bSign=+1`, once for \c U2 with
 * \f$a=U_{2m1}\f$/\f$b=U_{12}\f$/`bSign=-1`), which previously
 * duplicated this derivation.
 * \param[in,out] En Electric field to add the force to.
 * \param[in] a First plaquette combination.
 * \param[in] b Second plaquette combination.
 * \param[in] bSign Sign \p b enters the plaquette combination with
 * (`+1` or `-1`).
 * \param[in] phiN \f$\phi\f$ parallel-transported to this plaquette's
 * far corner.
 * \param[in] phi This cell's \f$\phi\f$.
 * \param[in] coeffPlaq Prefactor multiplying the plaquette force
 * (\f$i\tau d\tau/(2g^2)\f$).
 * \param[in] coeffComm Prefactor multiplying the commutator force
 * (\f$i d\tau/\tau\f$).
 */
inline void addEForceSU3(
    Matrix &En, const Matrix &a, const Matrix &b, double bSign,
    const Matrix &phiN, const Matrix &phi, const complex<double> coeffPlaq,
    const complex<double> coeffComm) {
    complex<double> *E = En.data();
    const complex<double> *A = a.data();
    const complex<double> *B = b.data();
    const su3::Matrix3 comm = su3::commutator(phiN, phi);

    // The plaquette force is the traceless anti-Hermitian part of
    // M = a + bSign*b.  Form it directly instead of materializing M, M^dagger,
    // an identity-matrix scale, and the associated Matrix temporaries.
    const complex<double> m00 = A[0] + bSign * B[0];
    const complex<double> m11 = A[4] + bSign * B[4];
    const complex<double> m22 = A[8] + bSign * B[8];
    const complex<double> traceThird =
        ((m00 - std::conj(m00)) + (m11 - std::conj(m11))
         + (m22 - std::conj(m22)))
        / 3.0;

    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            const int idx = 3 * row + col;
            const int tidx = 3 * col + row;
            const complex<double> mij = A[idx] + bSign * B[idx];
            const complex<double> mji = A[tidx] + bSign * B[tidx];
            complex<double> plaq = mij - std::conj(mji);
            if (row == col) plaq -= traceThird;
            E[idx] += coeffPlaq * plaq + coeffComm * comm.e[idx];
        }
    }

    // E is constrained to be traceless.  The old code subtracts
    // trace(E)/3 times the identity Matrix; only the diagonal entries change.
    const complex<double> eTraceThird = (E[0] + E[4] + E[8]) / 3.0;
    E[0] -= eTraceThird;
    E[4] -= eTraceThird;
    E[8] -= eTraceThird;
}

/**
 * `Evolution::evolveU()`'s per-cell kernel, run inside an existing
 * `#pragma omp parallel` region (its own `#pragma omp for` provides
 * the worksharing): rotates `Ux`/`Uy` by a second-order Padé
 * approximant of \f$\exp(i g^2 d\tau/(\tau+d\tau/2)\,U)\f$, using the
 * electric fields `U`/`U2` as the generator.
 * \param[in,out] lat Lattice whose `Ux`/`Uy` are updated in place.
 * \param[in] N Lattice side length.
 * \param[in] g Coupling \f$g\f$.
 * \param[in] dtau Time step [lattice units].
 * \param[in] tau Current proper time [lattice units].
 * \param[in,out] scratch Thread-local scratch storage.
 */
void evolveUTeam(
    Lattice *lat, int N, double g, double dtau, double tau,
    EvolveUScratch &scratch) {
    const int n = 2;
    const complex<double> iOmega(0., g * g * dtau / (tau + dtau / 2.));

#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
        scratch.E1 = iOmega * lat->U[pos];

        scratch.temp2 = scratch.one + 1. / (double)n * scratch.E1;
        for (int in = 0; in < n - 1; in++) {
            scratch.temp1 = scratch.E1 * scratch.temp2;
            scratch.temp2 =
                scratch.one + 1. / (double)(n - 1 - in) * scratch.temp1;
        }

        scratch.E1 = scratch.temp2;

        scratch.E2 = iOmega * lat->U2[pos];

        scratch.temp2 = scratch.one + 1. / (double)n * scratch.E2;
        for (int in = 0; in < n - 1; in++) {
            scratch.temp1 = scratch.E2 * scratch.temp2;
            scratch.temp2 =
                scratch.one + 1. / (double)(n - 1 - in) * scratch.temp1;
        }

        scratch.E2 = scratch.temp2;

        lat->Ux[pos] = (scratch.E1 * lat->Ux[pos]);
        lat->Uy[pos] = (scratch.E2 * lat->Uy[pos]);
    }
}

/**
 * `Evolution::evolvePhi()`'s per-cell kernel, run inside an existing
 * `#pragma omp parallel` region: \f$\phi \mathrel{+}=
 * (\tau+d\tau/2)\,d\tau\,\pi\f$.
 * \param[in,out] lat Lattice whose \c Uy2 (\f$\phi\f$) is updated in
 * place.
 * \param[in] N Lattice side length.
 * \param[in] dtau Time step [lattice units].
 * \param[in] tau Current proper time [lattice units].
 * \param[in,out] scratch Thread-local scratch storage.
 */
void evolvePhiTeam(
    Lattice *lat, int N, double dtau, double tau, EvolvePhiScratch &scratch) {
#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
        scratch.phi = lat->Uy2[pos];
        scratch.pi = lat->Ux2[pos];

        scratch.phi = scratch.phi + (tau + dtau / 2.) * dtau * scratch.pi;

        lat->Uy2[pos] = (scratch.phi);
    }
}

/**
 * `Evolution::evolvePi()`'s per-cell kernel, run inside an existing
 * `#pragma omp parallel` region: adds \f$(d\tau/\tau)\f$ times the
 * covariant discrete Laplacian of \f$\phi\f$ (parallel-transported via
 * `Ux`/`Uy` to its four neighbors) to \f$\pi\f$.
 * \param[in,out] lat Lattice whose \c Ux2 (\f$\pi\f$) is updated in
 * place.
 * \param[in] N Lattice side length.
 * \param[in] dtau Time step [lattice units].
 * \param[in] tau Current proper time [lattice units].
 * \param[in,out] scratch Thread-local scratch storage.
 */
void evolvePiTeam(
    Lattice *lat, int N, double dtau, double tau, EvolvePiScratch &scratch) {
    const double dtauOverTau = dtau / tau;

#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
        scratch.Ux = lat->Ux[pos];
        scratch.Uy = lat->Uy[pos];
        scratch.pi = lat->Ux2[pos];
        scratch.phi = lat->Uy2[pos];

        scratch.phiX =
            scratch.Ux
            * scratch.Ux.prodABconj(lat->Uy2[lat->pospX[pos]], scratch.Ux);
        scratch.phiY =
            scratch.Uy
            * scratch.Uy.prodABconj(lat->Uy2[lat->pospY[pos]], scratch.Uy);

        scratch.UxXm1 = lat->Ux[lat->posmX[pos]];
        scratch.UyYm1 = lat->Uy[lat->posmY[pos]];

        scratch.phimX =
            scratch.Ux.prodAconjB(scratch.UxXm1, lat->Uy2[lat->posmX[pos]])
            * scratch.UxXm1;
        scratch.phimY =
            scratch.Ux.prodAconjB(scratch.UyYm1, lat->Uy2[lat->posmY[pos]])
            * scratch.UyYm1;

        scratch.bracket = scratch.phiX + scratch.phimX + scratch.phiY
                          + scratch.phimY - 4. * scratch.phi;

        scratch.pi += dtauOverTau * scratch.bracket;

        lat->Ux2[pos] = (scratch.pi);
    }
}

/**
 * `Evolution::evolveE()`'s per-cell kernel, run inside an existing
 * `#pragma omp parallel` region: adds the traceless anti-Hermitian
 * plaquette force (from the four spatial plaquettes touching each
 * link) and the \f$[\phi_N,\phi]\f$ commutator force to `U`/`U2`,
 * via the shared addEForceSU3() kernel.
 * \param[in,out] lat Lattice whose `U`/`U2` are updated in place.
 * \param[in] N Lattice side length.
 * \param[in] g Coupling \f$g\f$.
 * \param[in] dtau Time step [lattice units].
 * \param[in] tau Current proper time [lattice units].
 * \param[in,out] scratch Thread-local scratch storage.
 */
void evolveETeam(
    Lattice *lat, int N, double g, double dtau, double tau,
    EvolveEScratch &scratch) {
    const complex<double> coeffPlaq =
        complex<double>(0., 1.) * tau * dtau / (2. * g * g);
    const complex<double> coeffComm = complex<double>(0., 1.) * dtau / tau;

#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
        const int posmXpY = lat->posmXpY[pos];
        const int pospXmY = lat->pospXmY[pos];

        scratch.En = lat->U[pos];
        scratch.phi = lat->Uy2[pos];
        scratch.phiN = lat->Uy2[lat->pospX[pos]];
        scratch.Ux = lat->Ux[pos];
        scratch.phiN =
            scratch.Ux * scratch.Ux.prodABconj(scratch.phiN, scratch.Ux);

        scratch.Uy = lat->Uy[pos];
        scratch.temp1 = lat->Ux[lat->pospY[pos]];
        scratch.temp1.conjg();
        scratch.U12 = (scratch.Ux * lat->Uy[lat->pospX[pos]])
                      * (scratch.Ux.prodABconj(scratch.temp1, scratch.Uy));

        scratch.temp1 = lat->Ux[lat->posmY[pos]];
        scratch.temp2 = lat->Uy[pospXmY];
        scratch.U1m2 =
            (scratch.Ux.prodABconj(scratch.Ux, scratch.temp2))
            * (scratch.Ux.prodAconjB(scratch.temp1, lat->Uy[lat->posmY[pos]]));

        scratch.temp1 = lat->Uy[lat->posmX[pos]];
        scratch.temp2 = lat->Ux[posmXpY];
        scratch.U2m1 =
            (scratch.Ux.prodABconj(scratch.Uy, scratch.temp2))
            * (scratch.Ux.prodAconjB(scratch.temp1, lat->Ux[lat->posmX[pos]]));

        // U12 + U1m2 - U12^dagger - U1m2^dagger is the
        // anti-Hermitian part of U12 + U1m2.
        addEForceSU3(
            scratch.En, scratch.U12, scratch.U1m2, 1.0, scratch.phiN,
            scratch.phi, coeffPlaq, coeffComm);
        lat->U[pos] = scratch.En;

        scratch.phiN = lat->Uy2[lat->pospY[pos]];
        scratch.phiN =
            scratch.Uy * scratch.Uy.prodABconj(scratch.phiN, scratch.Uy);

        scratch.En = lat->U2[pos];
        // U12^dagger + U2m1 - U12 - U2m1^dagger is the
        // anti-Hermitian part of U2m1 - U12.
        addEForceSU3(
            scratch.En, scratch.U2m1, scratch.U12, -1.0, scratch.phiN,
            scratch.phi, coeffPlaq, coeffComm);
        lat->U2[pos] = scratch.En;
    }
}

/**
 * Records the elapsed wall time since \p started under \p phase in
 * the global profiler, without restarting the clock. Called from
 * inside `#pragma omp single` in evolveStepPersistent(), for phases
 * whose elapsed time is the final one measured in that scope.
 * \param[in] phase Profiler phase name to add the elapsed time to.
 * \param[in] started Wall-clock time the phase began (from
 * `ipg::wallSeconds()`).
 */
void addTeamPhase(const char *phase, double started) {
    ipg::Profiler::instance().add(phase, ipg::wallSeconds() - started);
}

/**
 * Same as addTeamPhase(), but also resets \p started to the current
 * time, so the next phase's elapsed time is measured from here.
 * \param[in] phase Profiler phase name to add the elapsed time to.
 * \param[in,out] started Wall-clock time the phase began; updated to
 * now on return.
 */
void addPhaseAndRestart(const char *phase, double &started) {
    const double now = ipg::wallSeconds();
    ipg::Profiler::instance().add(phase, now - started);
    started = now;
}

/**
 * Runs one leapfrog step (\f$\pi\f$, then `U`/`U2`, then --
 * if \p updateCoordinates -- \f$\phi\f$ and `Ux`/`Uy`) inside a
 * single shared `#pragma omp parallel` region: each team function's
 * own `#pragma omp for` provides the worksharing, and the implicit
 * barrier at the end of each keeps the update order correct without
 * tearing down and rebuilding the thread team between phases.
 * \param[in,out] lat Lattice to evolve in place.
 * \param[in] param Simulation parameters.
 * \param[in] dtau Time step [lattice units].
 * \param[in] tau Current proper time [lattice units].
 * \param[in] updateCoordinates Whether to also update the coordinate
 * fields (\f$\phi\f$, `Ux`/`Uy`) this step, or only the momenta
 * (\f$\pi\f$, `U`/`U2`) -- `false` for a measurement-only half-step
 * that must not advance the trajectory.
 */
void evolveStepPersistent(
    Lattice *lat, Parameters *param, double dtau, double tau,
    bool updateCoordinates) {
    IPG_PROFILE_SCOPE("evolution.parallel_step");
    const int N = param->lattice.size;
    const double g = param->coupling.g;
    double phaseStart = 0.0;

#pragma omp parallel shared(phaseStart)
    {
        EvolveUScratch uScratch;
        EvolvePhiScratch phiScratch;
        EvolvePiScratch piScratch;
        EvolveEScratch eScratch;

#pragma omp single
        {
            phaseStart = ipg::wallSeconds();
        }
        evolvePiTeam(lat, N, dtau, tau, piScratch);
#pragma omp single
        {
            addTeamPhase("evolution.evolvePi", phaseStart);
            phaseStart = ipg::wallSeconds();
        }

        evolveETeam(lat, N, g, dtau, tau, eScratch);
#pragma omp single
        {
            addTeamPhase("evolution.evolveE", phaseStart);
            phaseStart = ipg::wallSeconds();
        }

        if (updateCoordinates) {
            evolvePhiTeam(lat, N, dtau, tau, phiScratch);
#pragma omp single
            {
                addTeamPhase("evolution.evolvePhi", phaseStart);
                phaseStart = ipg::wallSeconds();
            }

            evolveUTeam(lat, N, g, dtau, tau, uScratch);
#pragma omp single
            {
                addTeamPhase("evolution.evolveU", phaseStart);
            }
        }
    }
}

/**
 * Computes the local-coupling factor \f$g^2/(4\pi\alpha_s)\f$ at one
 * cell, used to rescale \f$T^{\mu\nu}\f$/\f$\epsilon\f$-derived
 * quantities when running coupling is enabled (\f$\alpha_s\f$ runs
 * with either the local \f$Q_s\f$ at this cell or one of the
 * event-averaged \f$Q_s\f$ choices, per
 * `param->coupling.runWithLocalQs`/`coupling.runWithQs`).
 * \param[in] lat Lattice to read \f$g^2\mu_A^2\f$/\f$g^2\mu_B^2\f$
 * from (only used when `coupling.runWithLocalQs`).
 * \param[in] param Simulation parameters.
 * \param[in] pos Flat cell index.
 * \param[in] N Lattice side length.
 * \param[in] a Lattice spacing [fm].
 * \param[in] g Coupling \f$g\f$.
 * \param[in] c Running-coupling shape parameter.
 * \param[in] muZero \f$\mu_0\f$ in the running-coupling formula.
 * \return The local-coupling factor; `1` if running coupling is
 * disabled.
 * \see computeRunningCouplingGfactorFromScale(), which this calls
 * with either the local or an event-averaged \f$Q_s\f$ as the scale.
 */
double computeRunningCouplingGfactor(
    Lattice *lat, Parameters *param, int pos, int N, double a, double g,
    double c, double muZero) {
    if (!param->coupling.runningCoupling) return 1.;

    const double lambdaQCD = param->coupling.LambdaQCD;
    const int nFlavors = param->coupling.nFlavors;

    if (param->coupling.runWithLocalQs) {
        // run with the local (in transverse plane) coupling
        const bool inBounds =
            lat->xFromPosition(pos) > 0 && lat->xFromPosition(pos) < N - 1
            && lat->yFromPosition(pos) > 0 && lat->yFromPosition(pos) < N - 1;
        const double g2mu2A = inBounds ? lat->cells[pos]->getg2mu2A() : 0;
        const double g2mu2B = inBounds ? lat->cells[pos]->getg2mu2B() : 0;

        double Qs = 0.;
        if (param->coupling.runWithQs == 2) {
            Qs = sqrt(
                std::max(g2mu2A, g2mu2B) * param->colorCharge.QsMuRatio
                * param->colorCharge.QsMuRatio / a / a * hbarc * hbarc
                * param->coupling.g * param->coupling.g);
        } else if (param->coupling.runWithQs == 0) {
            Qs = sqrt(
                std::min(g2mu2A, g2mu2B) * param->colorCharge.QsMuRatio
                * param->colorCharge.QsMuRatio / a / a * hbarc * hbarc
                * param->coupling.g * param->coupling.g);
        } else if (param->coupling.runWithQs == 1) {
            Qs = sqrt(
                (g2mu2A + g2mu2B) / 2. * param->colorCharge.QsMuRatio
                * param->colorCharge.QsMuRatio / a / a * hbarc * hbarc
                * param->coupling.g * param->coupling.g);
        }

        return computeRunningCouplingGfactorFromScale(
            g, muZero, c, lambdaQCD, nFlavors,
            param->coupling.runningCouplingQsFactor * Qs);
    } else {
        double averageQs = 0.;
        if (param->coupling.runWithQs == 0)
            averageQs = param->event.averageQsmin;
        else if (param->coupling.runWithQs == 1)
            averageQs = param->event.averageQsAvg;
        else if (param->coupling.runWithQs == 2)
            averageQs = param->event.averageQs;

        return computeRunningCouplingGfactorFromScale(
            g, muZero, c, lambdaQCD, nFlavors,
            param->coupling.runningCouplingQsFactor * averageQs);
    }
}

/**
 * Fills \p E1[pos] from \p sourceField[pos] (one of `lat->U`/`U2`/
 * `Ux2`), scaling by the square root of the local running-coupling
 * g-factor (computeRunningCouplingGfactor()) unless \f$\alpha_s\f$
 * runs with \f$k_T\f$ (in which case the \f$k_T\f$-dependent factor
 * is applied later, per-mode, in accumulateGluonSpectrum() instead).
 * \param[in] lat Lattice to read \f$g^2\mu_A^2\f$/\f$g^2\mu_B^2\f$
 * from, forwarded to computeRunningCouplingGfactor().
 * \param[in] param Simulation parameters.
 * \param[in] N Lattice side length.
 * \param[in] a Lattice spacing [fm].
 * \param[in] g Coupling \f$g\f$.
 * \param[in] c Running-coupling shape parameter.
 * \param[in] muZero \f$\mu_0\f$ in the running-coupling formula.
 * \param[in] sourceField Field to copy from (`lat->U`, `lat->U2`, or
 * `lat->Ux2`).
 * \param[out] E1 Filled with the (optionally rescaled) field, as a
 * pointer-per-cell view ready for `FFT::fftn`.
 */
void prepareSpectrumField(
    Lattice *lat, Parameters *param, int N, double a, double g, double c,
    double muZero, const std::vector<Matrix> &sourceField,
    std::vector<Matrix *> &E1) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int pos = lat->positionFromXY(i, j);
            double gfactor = computeRunningCouplingGfactor(
                lat, param, pos, N, a, g, c, muZero);
            if (!param->coupling.runWithKt) {
                *E1[pos] = sourceField[pos] * sqrt(gfactor);
            } else {
                *E1[pos] = sourceField[pos];
            }
        }
    }
}

/**
 * Accumulates \p E1's (already FFT'd) momentum-space spectrum into
 * `dNdeta`/`dEdeta` and the `n`/`E`/`n2` \f$k_T\f$ bins. Called
 * once each for the \c U, \c U2, and \c Ux2 (\f$\pi\f$) fields by
 * `Evolution::multiplicity()`.
 * \param[in] param Simulation parameters.
 * \param[in] N Lattice side length.
 * \param[in] it Current time step index.
 * \param[in] dtau Time step [lattice units].
 * \param[in] g Coupling \f$g\f$.
 * \param[in] a Lattice spacing [fm].
 * \param[in] c Running-coupling shape parameter.
 * \param[in] muZero \f$\mu_0\f$ in the running-coupling formula.
 * \param[in] dkt Momentum-bin width [lattice units].
 * \param[in] bins Number of \f$k_T\f$ bins in `n`/`E`/`n2`/
 * \p counter.
 * \param[in] E1 The field's FFT'd momentum-space values, as a
 * pointer-per-cell view.
 * \param[in] useElectricNormalization Selects \c nkt's electric-field
 * normalization (`true`, \f$g^2/((it-0.5)d\tau)\f$) vs. the
 * \f$\pi\f$-field normalization (`false`, \f$(it-0.5)d\tau\f$, no
 * \f$g^2\f$).
 * \param[in] accumulateCounter Whether to also increment \p counter;
 * only one of the three per-field passes needs to, since all three
 * share the same \f$k_T\f$ grid.
 * \param[in,out] dNdeta Running sum of \f$dN/dy\f$ (or \f$d\eta\f$),
 * incremented by this field's contribution.
 * \param[in,out] dEdeta Running sum of \f$dE/dy\f$ (or \f$d\eta\f$),
 * incremented by this field's contribution.
 * \param[in,out] n Binned \f$dN/d^2k_T\f$, incremented by this
 * field's contribution, length \p bins.
 * \param[in,out] E Binned \f$dE/d^2k_T\f$, incremented by this
 * field's contribution, length \p bins.
 * \param[in,out] n2 Alternate binned \f$dN/d^2k_T\f$ normalization
 * (used for a cross-check), incremented by this field's contribution,
 * length \p bins.
 * \param[in,out] counter Number of lattice momentum modes falling
 * into each bin, incremented only if \p accumulateCounter, length
 * \p bins.
 */
void accumulateGluonSpectrum(
    Parameters *param, int N, int it, double dtau, double g, double a, double c,
    double muZero, double dkt, int bins, const std::vector<Matrix *> &E1,
    bool useElectricNormalization, bool accumulateCounter, double &dNdeta,
    double &dEdeta, double *n, double *E, double *n2, int *counter) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double nkt = 0.;
            int pos = latticeIndex(i, j, N);
            int npos = latticeIndex(N - i, N - j, N);

            double kx =
                2. * M_PI
                * (-0.5 + static_cast<double>(i) / static_cast<double>(N));
            double ky =
                2. * M_PI
                * (-0.5 + static_cast<double>(j) / static_cast<double>(N));
            double kt2 =
                4.
                * (sin(kx / 2.) * sin(kx / 2.) + sin(ky / 2.) * sin(ky / 2.));
            double omega2 =
                4.
                * (sin(kx / 2.) * sin(kx / 2.)
                   + sin(ky / 2.) * sin(ky / 2.));  // lattice dispersion
                                                    // relation (this is
                                                    // omega squared)

            // i=0 or j=0 have no negative k_T value available
            if (i != 0 && j != 0) {
                if (omega2 != 0) {
                    if (useElectricNormalization) {
                        nkt =
                            2. / sqrt(omega2) / static_cast<double>(N * N)
                            * (g * g / ((it - 0.5) * dtau)
                               * ((((*E1[pos]) * (*E1[npos])).trace()).real()));
                    } else {
                        nkt =
                            2. / sqrt(omega2) / static_cast<double>(N * N)
                            * (((it - 0.5) * dtau)
                               * ((((*E1[pos]) * (*E1[npos])).trace()).real()));
                    }
                    if (param->coupling.runWithKt) {
                        nkt *= computeRunningCouplingGfactorFromScale(
                            g, muZero, c, param->coupling.LambdaQCD,
                            param->coupling.nFlavors,
                            param->coupling.runningCouplingQsFactor * sqrt(kt2)
                                * hbarc / a);
                    }
                }

                dNdeta += nkt;
                dEdeta += nkt * sqrt(omega2) * hbarc / a;

                for (int ik = 0; ik < bins; ik++) {
                    if (abs(sqrt(kt2)) > ik * dkt
                        && abs(sqrt(kt2)) <= (ik + 1) * dkt) {
                        n[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2) * 2 * M_PI
                                 * sqrt(kt2) * dkt * N * N / M_PI / M_PI / 2.
                                 / 2.;
                        E[ik] += sqrt(omega2) * hbarc / a * nkt / dkt / 2 / M_PI
                                 / sqrt(kt2) * 2 * M_PI * sqrt(kt2) * dkt * N
                                 * N / M_PI / M_PI / 2. / 2.;
                        n2[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2);
                        // dividing by bin size; bin is dkt times Jacobian
                        // k(=ik*dkt) times 2Pi in phi times the correct
                        // number of counts for an infinite lattice: area in
                        // bin divided by total area
                        if (accumulateCounter) {
                            counter[ik] += 1;  // number of entries in n[ik]
                        }
                    }
                }
            }
        }
    }
}

/**
 * `Evolution::multiplicity()`'s per-bin \f$dN/dy\f$, \f$dE/dy\f$
 * weight: the phase-space factor \f$(ik+0.5)\,dk_T^2\,2\pi\f$, times
 * a Jacobian ratio when the rapidity input is actually a
 * pseudorapidity (the same factor previously computed identically
 * three times -- unconditionally, and again inside the \f$k_T>3\f$
 * and \f$k_T>6\f$ GeV cuts -- for both `usePseudoRapidity` branches).
 * \param[in] param Simulation parameters.
 * \param[in] m Jacobian mass term [GeV] (`param->colorCharge.jacobianMass`).
 * \param[in] ik Bin index.
 * \param[in] dkt Momentum-bin width [lattice units].
 * \param[in] a Lattice spacing [fm].
 * \return The per-bin weight.
 */
double computeMultiplicityBinWeight(
    Parameters *param, double m, int ik, double dkt, double a) {
    const double base = (ik + 0.5) * dkt * dkt * 2. * M_PI;
    if (!param->colorCharge.usePseudoRapidity) {
        return base;
    }
    return base * cosh(param->colorCharge.rapidity())
           / (sqrt(
               pow(cosh(param->colorCharge.rapidity()), 2.)
               + m * m
                     / (((ik + 0.5) * dkt / a * hbarc)
                        * ((ik + 0.5) * dkt / a * hbarc))));
}

}  // namespace

void Evolution::evolveU(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    IPG_PROFILE_SCOPE("evolution.evolveU");
    const int N = param->lattice.size;
    const double g = param->coupling.g;

#pragma omp parallel
    {
        EvolveUScratch scratch;
        evolveUTeam(lat, N, g, dtau, tau, scratch);
    }
}

void Evolution::evolvePhi(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    IPG_PROFILE_SCOPE("evolution.evolvePhi");
    const int N = param->lattice.size;

#pragma omp parallel
    {
        EvolvePhiScratch scratch;
        evolvePhiTeam(lat, N, dtau, tau, scratch);
    }
}

void Evolution::evolvePi(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    IPG_PROFILE_SCOPE("evolution.evolvePi");
    const int N = param->lattice.size;

#pragma omp parallel
    {
        EvolvePiScratch scratch;
        evolvePiTeam(lat, N, dtau, tau, scratch);
    }
}

void Evolution::evolveE(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    IPG_PROFILE_SCOPE("evolution.evolveE");
    const int N = param->lattice.size;
    const double g = param->coupling.g;

#pragma omp parallel
    {
        EvolveEScratch scratch;
        evolveETeam(lat, N, g, dtau, tau, scratch);
    }
}

void Evolution::checkGaussLaw(Lattice *lat, Parameters *param) {
    IPG_PROFILE_SCOPE("diagnostics.gauss_law");
    const int N = param->lattice.size;

    Matrix Ux;
    Matrix UxXm1;
    Matrix UxYm1;

    Matrix Uy;
    Matrix UyXm1;
    Matrix UyYm1;

    Matrix UxDag;
    Matrix UxXm1Dag;
    Matrix UxYm1Dag;

    Matrix UyDag;
    Matrix UyXm1Dag;
    Matrix UyYm1Dag;

    Matrix E1;
    Matrix E2;
    Matrix E1mX;
    Matrix E2mY;
    Matrix phi;
    Matrix pi;

    Matrix Gauss;
    double largest = 0;

    for (int pos = 0; pos < N * N; pos++) {
        // retrieve current Ux and Uy
        Ux = lat->Ux[pos];
        Uy = lat->Uy[pos];
        UxDag = Ux;
        UxDag.conjg();
        UyDag = Uy;
        UyDag.conjg();

        UxXm1 = lat->Ux[lat->posmX[pos]];
        UxYm1 = lat->Ux[lat->posmY[pos]];
        UxXm1Dag = UxXm1;
        UxXm1Dag.conjg();
        UxYm1Dag = UxYm1;
        UxYm1Dag.conjg();

        UyXm1 = lat->Uy[lat->posmX[pos]];
        UyYm1 = lat->Uy[lat->posmY[pos]];
        UyXm1Dag = UyXm1;
        UyXm1Dag.conjg();
        UyYm1Dag = UyYm1;
        UyYm1Dag.conjg();

        // retrieve current E1 and E2 (that's the one defined at tau-dtau/2)
        E1 = lat->U[pos];
        E2 = lat->U2[pos];
        E1mX = lat->U[lat->posmX[pos]];
        E2mY = lat->U2[lat->posmY[pos]];
        // retrieve current phi (at time tau) at this x_T
        phi = lat->Uy2[pos];
        // retrieve current pi
        pi = lat->Ux2[pos];

        Gauss = UxXm1Dag * E1mX * UxXm1 - E1 + UyYm1Dag * E2mY * UyYm1 - E2
                - complex<double>(0., 1.) * (phi * pi - pi * phi);

        if (Gauss.square() > largest) largest = Gauss.square();
    }
    messager_ << "[Evolution::checkGaussLaw]: Gauss violation=" << largest;
    messager_.flush("info");
}

void Evolution::writeEvolvedFields(Lattice *lat, Parameters *param, int it) {
    IPG_PROFILE_SCOPE("output.evolved_fields");
    const int N = param->lattice.size;
    const double a = param->lattice.L / static_cast<double>(N);
    const double dtau = param->run.dtau;
    const double tauLattice = static_cast<double>(it) * dtau;
    const double tauFm = a * tauLattice;
    const double momentumTauFm =
        (it == 0) ? 0.0 : a * (static_cast<double>(it) - 0.5) * dtau;

    // The payload layout is [field, real_or_imag, x, y, row, col], C-order.
    // Full matrices are stored so that no color information is discarded.
    constexpr int nFields = 6;
    const std::size_t matrixElements =
        static_cast<std::size_t>(N) * N * Nc * Nc;
    const std::size_t payloadElements =
        static_cast<std::size_t>(nFields) * 2 * matrixElements;
    std::vector<float> payload(payloadElements);

    auto matrixAt = [lat](const int field, const int pos) -> const Matrix & {
        switch (field) {
            case 0:
                return lat->Uy2[pos];
            case 1:
                return lat->Ux2[pos];
            case 2:
                return lat->U[pos];
            case 3:
                return lat->U2[pos];
            case 4:
                return lat->Ux[pos];
            case 5:
                return lat->Uy[pos];
            default:
                throw std::runtime_error("invalid evolved-field index");
        }
    };

    for (int field = 0; field < nFields; ++field) {
        const std::size_t realOffset =
            static_cast<std::size_t>(2 * field) * matrixElements;
        const std::size_t imagOffset = realOffset + matrixElements;
        for (int x = 0; x < N; ++x) {
            for (int y = 0; y < N; ++y) {
                const int pos = lat->positionFromXY(x, y);
                const Matrix &matrix = matrixAt(field, pos);
                const std::complex<double> *elements = matrix.data();
                const std::size_t siteOffset =
                    static_cast<std::size_t>(pos) * Nc * Nc;
                for (int row = 0; row < Nc; ++row) {
                    for (int col = 0; col < Nc; ++col) {
                        const std::size_t element =
                            static_cast<std::size_t>(row) * Nc + col;
                        payload[realOffset + siteOffset + element] =
                            static_cast<float>(elements[element].real());
                        payload[imagOffset + siteOffset + element] =
                            static_cast<float>(elements[element].imag());
                    }
                }
            }
        }
    }

    // This binary format is explicitly little-endian. IP-Glasma production
    // platforms are normally little-endian; fail loudly rather than emit an
    // ambiguous file on another architecture.
    const std::uint16_t endianProbe = 1;
    if (*reinterpret_cast<const unsigned char *>(&endianProbe) != 1) {
        throw std::runtime_error(
            "writeEvolvedFields currently requires a little-endian host");
    }

    std::stringstream metadata;
    metadata << std::setprecision(17)
             << "{\"format\":\"ipglasma-evolved-fields\","
             << "\"version\":1,"
             << "\"dtype\":\"<f4\","
             << "\"shape\":[6,2," << N << "," << N << "," << Nc << "," << Nc
             << "],"
             << "\"axis_order\":[\"field\",\"complex_part\",\"x\",\"y\","
                "\"row\",\"col\"],"
             << "\"fields\":[\"phi\",\"pi\",\"E1\",\"E2\",\"Ux\","
                "\"Uy\"],"
             << "\"complex_part\":[\"real\",\"imag\"],"
             << "\"native_site_index\":\"pos=x*N+y\","
             << "\"event_id\":" << param->event.eventId << ","
             << "\"step\":" << it << ","
             << "\"tau_lattice\":" << tauLattice << ","
             << "\"tau_fm\":" << tauFm << ","
             << "\"momentum_tau_fm\":" << momentumTauFm << ","
             << "\"a_fm\":" << a << ","
             << "\"dtau_lattice\":" << dtau << ","
             << "\"staggering\":\"Ux,Uy,phi at tau; E1,E2,pi at tau-dtau/2 "
                "for step>0; all variables are the initialized tau=0+ values "
                "for step=0\"}";
    const std::string metadataString = metadata.str();

    stringstream filename;
    filename << "evolvedFields" << param->event.eventId << "_it" << std::setw(8)
             << std::setfill('0') << it << ".ipgf";

    ofstream output(
        filename.str().c_str(),
        std::ios::out | std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "could not open evolved-field snapshot " + filename.str());
    }

    const char magic[8] = {'I', 'P', 'G', 'F', 'L', 'D', '1', '\0'};
    const std::uint64_t metadataBytes =
        static_cast<std::uint64_t>(metadataString.size());
    output.write(magic, sizeof(magic));
    output.write(
        reinterpret_cast<const char *>(&metadataBytes), sizeof(metadataBytes));
    output.write(metadataString.data(), metadataString.size());
    output.write(
        reinterpret_cast<const char *>(payload.data()),
        static_cast<std::streamsize>(payload.size() * sizeof(float)));
    output.close();

    if (!output) {
        throw std::runtime_error(
            "failed while writing evolved-field snapshot " + filename.str());
    }
    messager_ << "[Evolution::writeEvolvedFields]: Wrote evolved fields "
                 "at tau="
              << tauFm << " fm/c to " << filename.str();
    messager_.flush("info");
}

void Evolution::writeGluonMultiplicityTarget(
    Parameters *param, int it, double a, double dtau, double dNPrimary,
    double dNBinned, double dEPrimary, double dEBinned, double dNCut3,
    double dECut3, double dNCut6, double dECut6, const double *spectrumN,
    const double *spectrumE, const int *spectrumCounts, int bins, double dkt) {
    IPG_PROFILE_SCOPE("output.gluon_target");
    stringstream filename;
    filename << "gluonMultiplicity" << param->event.eventId << ".json";
    ofstream output(filename.str().c_str(), std::ios::out | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "could not open gluon-multiplicity target " + filename.str());
    }

    const char *rapidityVariable =
        (!param->colorCharge.usePseudoRapidity) ? "y" : "eta";
    const double meanKt = (dNPrimary != 0.0) ? dEPrimary / dNPrimary : 0.0;
    const double spectrumUnitFactor = (a / hbarc) * (a / hbarc);

    output << std::setprecision(17) << "{\n"
           << "  \"format\": \"ipglasma-gluon-target\",\n"
           << "  \"version\": 1,\n"
           << "  \"event_id\": " << param->event.eventId << ",\n"
           << "  \"step\": " << it << ",\n"
           << "  \"tau_fm\": " << static_cast<double>(it) * dtau * a << ",\n"
           << "  \"rapidity_variable\": \"" << rapidityVariable << "\",\n"
           << "  \"dN\": " << dNPrimary << ",\n"
           << "  \"dN_binned_check\": " << dNBinned << ",\n"
           << "  \"dE_GeV\": " << dEPrimary << ",\n"
           << "  \"dE_binned_check_GeV\": " << dEBinned << ",\n"
           << "  \"mean_kT_GeV\": " << meanKt << ",\n"
           << "  \"dN_kT_gt_3_GeV\": " << dNCut3 << ",\n"
           << "  \"dE_kT_gt_3_GeV\": " << dECut3 << ",\n"
           << "  \"dN_kT_gt_6_GeV\": " << dNCut6 << ",\n"
           << "  \"dE_kT_gt_6_GeV\": " << dECut6 << ",\n"
           << "  \"Npart\": " << param->event.Npart << ",\n"
           << "  \"Tpp\": " << param->event.Tpp << ",\n"
           << "  \"impact_parameter_fm\": " << param->event.b << ",\n"
           << "  \"random_seed\": " << param->run.randomSeed << ",\n"
           << "  \"spectrum_definition\": \"azimuthally averaged Coulomb-gauge "
              "gluon spectrum used by Evolution::multiplicity\",\n"
           << "  \"kt_GeV\": [";

    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << (static_cast<double>(ik) + 0.5) * dkt / a * hbarc;
    }
    output << "],\n  \"dN_d2k_GeV_minus2\": [";
    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << spectrumN[ik] * spectrumUnitFactor;
    }
    output << "],\n  \"dE_d2k_GeV_minus1\": [";
    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << spectrumE[ik] * spectrumUnitFactor;
    }
    output << "],\n  \"lattice_bin_counts\": [";
    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << spectrumCounts[ik];
    }
    output << "]\n}\n";
    output.close();

    if (!output) {
        throw std::runtime_error(
            "failed while writing gluon-multiplicity target " + filename.str());
    }
    messager_ << "[Evolution::writeGluonMultiplicityTarget]: Wrote gluon "
                 "target dN/d"
              << rapidityVariable << "=" << dNPrimary << " to "
              << filename.str();
    messager_.flush("info");
}

void Evolution::run(Lattice *lat, Group *group, Parameters *param) {
    IPG_PROFILE_SCOPE("evolution.total");
    int N = param->lattice.size;
    double a = param->lattice.L / N;  // lattice spacing in fm

    // do the first half step of the momenta (E1,E2,pi)
    // for now I use the \tau=0 value at \tau=d\tau/2.
    double dtau = param->run.dtau;  // dtau is in lattice units

    double maxtime = param->evolution.maxTime;  // maxtime is in fm
    if (param->evolution.inverseQsForMaxTime) {
        maxtime = 1. / param->event.averageQs * hbarc;
        messager_ << "[Evolution::run]: maximal evolution time = " << maxtime
                  << " fm";
        messager_.flush("info");
    }

    // E and Pi at tau=dtau/2 are equal to the initial ones (at tau=0)
    // now evolve phi and U to time tau=dtau.
    {
        IPG_PROFILE_SCOPE("evolution.initial_coordinate_half_step");
        evolvePhi(lat, param, dtau, 0.);
        evolveU(lat, param, dtau, 0.);
    }

    int itmax = static_cast<int>(maxtime / (a * dtau) + 0.00000000001);
    int it0 = static_cast<int>(0.1 / (a * dtau) + 0.0000000001);
    int it1 = static_cast<int>(0.2 / (a * dtau) + 0.0000000001);
    int it2 = static_cast<int>(0.3 / (a * dtau) + 0.0000000001);
    int it3 = static_cast<int>(0.4 / (a * dtau) + 0.0000000001);

    // Tmunu is defined at the integer coordinate time tau_n, while the
    // leapfrog momenta E1, E2, and pi live at tau_{n-1/2}.  Keep reusable
    // backups so measurements can temporarily center the momenta without
    // changing the actual evolution trajectory or making it output-dependent.
    const std::size_t latticeSites = static_cast<std::size_t>(N) * N;
    std::vector<Matrix> tmunuE1Backup(latticeSites);
    std::vector<Matrix> tmunuE2Backup(latticeSites);
    std::vector<Matrix> tmunuPiBackup(latticeSites);

    messager_ << "[Evolution::run]: Starting evolution: num of time steps="
              << itmax;
    messager_.flush("info");
    if ((param->output.writeOutputs == 5)) {
        messager_ << "[Evolution::run]: Measuring at times " << it0 * a * dtau
                  << ", " << it1 * a * dtau << ", " << it2 * a * dtau << ", "
                  << it3 * a * dtau << ", " << itmax * a * dtau << ". ";
        messager_.flush("info");
    }
    messager_ << "[Evolution::run]:  a = " << a;
    messager_.flush("info");
    messager_ << "[Evolution::run]:  dtau = " << dtau;
    messager_.flush("info");
    messager_ << "[Evolution::run]:  it0 = " << it0;
    messager_.flush("info");

    // do evolution
    for (int it = 1; it <= itmax; it++) {
        const bool finalTmunuMeasurement = (it == itmax);
        const bool intermediateTmunuMeasurement =
            (param->output.writeOutputs == 5)
            && (it == it0 || it == it1 || it == it2 || it == it3);
        const bool measureTmunu =
            finalTmunuMeasurement || intermediateTmunuMeasurement;

        if (measureTmunu) {
            std::copy(lat->U.begin(), lat->U.end(), tmunuE1Backup.begin());
            std::copy(lat->U2.begin(), lat->U2.end(), tmunuE2Backup.begin());
            std::copy(lat->Ux2.begin(), lat->Ux2.end(), tmunuPiBackup.begin());

            // Temporarily move E1, E2, and pi from tau_{n-1/2} to tau_n.
            // Coordinates are not advanced.  The original momenta are restored
            // after output, so enabling Tmunu measurements cannot change the
            // subsequent leapfrog trajectory.
            evolveStepPersistent(lat, param, dtau / 2., it * dtau, false);
        }

        if (finalTmunuMeasurement) {
            EnergyMomentumTensor::compute(lat, param, it);
            finalFlowMeasurement(lat, param, it);
        }

        if (intermediateTmunuMeasurement) {
            EnergyMomentumTensor::compute(lat, param, it);
            //  Preserve the historical intermediate-time finalFlag=false path
            //  when hydro output is enabled.
            if (param->output.writeEpsilonUHydro) {
                u(lat, param, it, false);
            } else {
                MyEigen myeigen;
                myeigen.writeTmunu4D(lat, param, it);
            }
        }

        if (measureTmunu) {
            std::copy(
                tmunuE1Backup.begin(), tmunuE1Backup.end(), lat->U.begin());
            std::copy(
                tmunuE2Backup.begin(), tmunuE2Backup.end(), lat->U2.begin());
            std::copy(
                tmunuPiBackup.begin(), tmunuPiBackup.end(), lat->Ux2.begin());
        }

        if (it % 10 == 1) {
            messager_ << "[Evolution::run]: Evolving to time " << it * a * dtau
                      << " fm/c";
            messager_.flush("info");
        }

        // Keep one OpenMP team alive across all leapfrog kernels in this
        // time step. The worksharing loops retain their implicit barriers,
        // preserving the Pi -> E -> phi -> U update order.
        if (it < itmax) {
            evolveStepPersistent(lat, param, dtau, (it)*dtau, true);
        } else {
            evolveStepPersistent(lat, param, dtau / 2., (it)*dtau, false);
        }

        if (it == itmax) {
            checkGaussLaw(lat, param);
        }

        int success = 1;
        if (param->output.computeGluonMultiplicity) {
            if (it == itmax) {
                eccentricity(lat, param, it, 0.0, 0);
                // eccentricity(lat, param, it, 0.1, 0);
                // eccentricity(lat, param, it, 1., 0);
                // eccentricity(lat, param, it, 10., 0);

                success = multiplicity(lat, group, param, it);
            }
        }

        if (success == 0) break;
    }
}

void Evolution::finalFlowMeasurement(Lattice *lat, Parameters *param, int it) {
    // Hydro flow fields are optional. Tmunu output remains available
    // through the lightweight writer when the expensive eigen solve is
    // disabled.
    if (param->output.writeEpsilonUHydro) {
        u(lat, param, it, true);
    } else {
        MyEigen myeigen;
        // eccentricity() weights by epsilon * u^tau, which only the
        // flow-velocity solve sets.
        if (param->output.computeGluonMultiplicity) {
            myeigen.solveFlowVelocity(lat, param, it);
        }
        myeigen.writeTmunu4D(lat, param, it);
    }
}

void Evolution::u(Lattice *lat, Parameters *param, int it, bool finalFlag) {
    IPG_PROFILE_SCOPE("observables.flow_velocity");
    MyEigen myeigen;
    myeigen.flowVelocity4D(lat, param, it, finalFlag);
}

namespace {
/// Result of computeRotatedAnisotropy(): the rotated- and
/// unrotated-frame \f$T^{xx}-T^{yy}\f$ spatial anisotropy sums.
struct AnisotropyResult {
    /// \f$\sum (T^{xx}_{\text{rot}}-T^{yy}_{\text{rot}})\f$ at the
    /// sampled angle \f$\Psi\f$.
    double num;
    /// \f$\sum (T^{xx}_{\text{rot}}+T^{yy}_{\text{rot}})\f$ at the
    /// sampled angle \f$\Psi\f$.
    double den;
    /// \f$\sum (T^{xx}-T^{yy})\f$ in the unrotated frame.
    double num2;
    /// \f$\sum (T^{xx}+T^{yy})\f$ in the unrotated frame.
    double den2;
};

/**
 * `Evolution::eccentricity()`'s `doAniso==1` branch samples this at
 * ten values of \p Psi (previously ten copy-pasted ~25-line blocks,
 * differing only in \p Psi): sums \f$T^{xx}-T^{yy}\f$ and
 * \f$T^{xx}+T^{yy}\f$ over the whole lattice, both in the unrotated
 * frame and after rotating \f$T^{xx}\f$/\f$T^{xy}\f$/\f$T^{yy}\f$ by
 * \p Psi.
 * \param[in] lat Lattice to read `Txx`/`Txy`/`Tyy` from.
 * \param[in] N Lattice side length.
 * \param[in] Psi Rotation angle [rad].
 * \return The summed rotated- and unrotated-frame anisotropy
 * numerators/denominators.
 */
AnisotropyResult computeRotatedAnisotropy(Lattice *lat, int N, double Psi) {
    double num = 0., den = 0., num2 = 0., den2 = 0.;
    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            int pos = lat->positionFromXY(ix, iy);

            double TxxRot = cos(Psi)
                                * (cos(Psi) * lat->cells[pos]->getTxx()
                                   - sin(Psi) * lat->cells[pos]->getTxy())
                            - sin(Psi)
                                  * (cos(Psi) * lat->cells[pos]->getTxy()
                                     - sin(Psi) * lat->cells[pos]->getTyy());
            double TyyRot = sin(Psi)
                                * (sin(Psi) * lat->cells[pos]->getTxx()
                                   + cos(Psi) * lat->cells[pos]->getTxy())
                            + cos(Psi)
                                  * (sin(Psi) * lat->cells[pos]->getTxy()
                                     + cos(Psi) * lat->cells[pos]->getTyy());

            num2 += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
            den2 += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();

            num += TxxRot - TyyRot;
            den += TxxRot + TyyRot;
        }
    }
    return {num, den, num2, den2};
}
}  // namespace

void Evolution::eccentricity(
    Lattice *lat, Parameters *param, int it, double cutoff, int doAniso) {
    IPG_PROFILE_SCOPE("observables.eccentricity");
    stringstream strecc_name;
    strecc_name << "eccentricities" << param->event.eventId << ".dat";
    string ecc_name;
    ecc_name = strecc_name.str();

    // cutoff on energy density is 'cutoff' times Lambda_QCD^4
    int N = param->lattice.size;
    int pos;
    double rA, phiA, x, y;
    double L = param->lattice.L;
    double a = L / N;  // lattice spacing in fm
    double eccentricity1, eccentricity2, eccentricity3, eccentricity4,
        eccentricity5, eccentricity6;
    double avcos, avsin, avcos1, avsin1, avcos3, avsin3, avrSq, avxSq, avySq,
        avr1, avr3, avcos4, avsin4, avr4, avcos5, avsin5, avr5, avcos6, avsin6,
        avr6;
    double Rbar;
    double Psi1, Psi2, Psi3, Psi4, Psi5, Psi6;
    double maxEps = 0;
    double g = param->coupling.g;

    double g2mu2A, g2mu2B, gfactor, alphas = 0., Qs = 0.;
    double c = param->coupling.c;
    double muZero = param->coupling.mu0;

    double weight;

    double area = 0.;
    double avgeden = 0.;
    int sum = 0;

    avrSq = 0.;
    avr3 = 0.;

    double avx = 0.;
    double avy = 0.;
    double toteps = 0.;
    int xshift;
    int yshift;
    double maxX = 0.;
    double maxY = 0.;

    double smallestX = 0.;
    double smallestY = 0.;
    double avgQs2AQs2B = 0.;

    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            pos = lat->positionFromXY(ix, iy);
            maxEps = std::max(lat->cells[pos]->getEpsilon(), maxEps);
        }
    }

    // first shift to the center
    for (int ix = 0; ix < N; ix++) {
        x = -L / 2. + a * ix;
        for (int iy = 0; iy < N; iy++) {
            y = -L / 2. + a * iy;
            pos = lat->positionFromXY(ix, iy);

            gfactor = computeRunningCouplingGfactor(
                lat, param, pos, N, a, g, c, muZero);

            if (lat->cells[pos]->getEpsilon() * gfactor
                < cutoff)  // this is 1/fm^4, so Lambda_QCD^{-4} (because
                           // \Lambda_QCD is roughly 1/fm)
            {
                weight = 0.;
            } else {
                weight =
                    (lat->cells[pos]->getEpsilon() * lat->cells[pos]->getutau()
                     * gfactor);
                area += a * a;
                sum += 1;
                avgeden += lat->cells[pos]->getEpsilon() * hbarc
                           * gfactor;  // GeV/fm^3
                g2mu2A = lat->cells[pos]->getg2mu2A();
                g2mu2B = lat->cells[pos]->getg2mu2B();
                avgQs2AQs2B += g2mu2A * param->colorCharge.QsMuRatio
                               * param->colorCharge.QsMuRatio * g2mu2B
                               * param->colorCharge.QsMuRatio
                               * param->colorCharge.QsMuRatio / a / a / a / a;
            }
            avx += x * weight;
            avy += y * weight;
            toteps += weight;
        }
    }

    avx /= toteps;
    avy /= toteps;
    avgeden /= double(sum);
    avgQs2AQs2B /= double(sum);
    param->event.area = area;

    xshift = static_cast<int>(floor(avx / a + 0.00000000001));
    yshift = static_cast<int>(floor(avy / a + 0.00000000001));

    avcos1 = 0.;
    avsin1 = 0.;
    avcos = 0.;
    avsin = 0.;
    avcos3 = 0.;
    avsin3 = 0.;
    avcos4 = 0.;
    avsin4 = 0.;
    avcos5 = 0.;
    avsin5 = 0.;
    avcos6 = 0.;
    avsin6 = 0.;
    avr1 = 0.;
    avrSq = 0.;
    avxSq = 0.;
    avySq = 0.;
    avr3 = 0.;
    avr4 = 0.;
    avr5 = 0.;
    avr6 = 0.;

    for (int ix = 2; ix < N - 2; ix++) {
        x = -L / 2. + a * ix - avx;
        for (int iy = 2; iy < N - 2; iy++) {
            pos = lat->positionFromXY(ix, iy);
            y = -L / 2. + a * iy - avy;
            if (x >= 0) {
                phiA = atan(y / x);
                if (x == 0) {
                    if (y >= 0)
                        phiA = M_PI / 2.;
                    else if (y < 0)
                        phiA = 3. * M_PI / 2.;
                }
            } else {
                phiA = atan(y / x) + M_PI;
            }

            gfactor = computeRunningCouplingGfactor(
                lat, param, pos, N, a, g, c, muZero);

            if (lat->cells[pos]->getEpsilon() * gfactor
                < cutoff)  // this is 1/fm^4, so Lambda_QCD^{-4}
            {
                weight = 0.;
            } else {
                weight =
                    (lat->cells[pos]->getEpsilon() * lat->cells[pos]->getutau()
                     * gfactor);
            }

            rA = sqrt(x * x + y * y);
            avr1 += rA * rA * rA * (weight);
            avrSq += rA * rA * (weight);  // compute average r^2
            avr3 += rA * rA * rA * (weight);
            avr4 += rA * rA * rA * rA * (weight);
            avr5 += rA * rA * rA * rA * rA * (weight);
            avr6 += rA * rA * rA * rA * rA * rA * (weight);

            avcos1 += rA * rA * rA * cos(phiA) * (weight);
            avsin1 += rA * rA * rA * sin(phiA) * (weight);
            avcos += rA * rA * cos(2. * phiA) * (weight);
            avsin += rA * rA * sin(2. * phiA) * (weight);
            avcos3 += rA * rA * rA * cos(3. * phiA) * (weight);
            avsin3 += rA * rA * rA * sin(3. * phiA) * (weight);
            avcos4 += rA * rA * rA * rA * cos(4. * phiA) * (weight);
            avsin4 += rA * rA * rA * rA * sin(4. * phiA) * (weight);
            avcos5 += rA * rA * rA * rA * rA * cos(5. * phiA) * (weight);
            avsin5 += rA * rA * rA * rA * rA * sin(5. * phiA) * (weight);
            avcos6 += rA * rA * rA * rA * rA * rA * cos(6. * phiA) * (weight);
            avsin6 += rA * rA * rA * rA * rA * rA * sin(6. * phiA) * (weight);

            if (weight > cutoff && iy == N / 2 + yshift) {
                maxX = x;
            }
            if (weight > cutoff && ix == N / 2 + xshift) {
                maxY = y;
            }

            if (weight < cutoff && iy == N / 2 + yshift && ix > N / 2 + xshift
                && smallestX == 0) {
                smallestX = x;
            }
            if (weight < cutoff && ix == N / 2 + xshift && iy > N / 2 + yshift
                && smallestY == 0) {
                smallestY = y;
            }
        }
    }

    // compute and print eccentricity and angles:
    Psi1 = (atan(avsin1 / avcos1) + M_PI) / 1.;
    Psi2 = (atan(avsin / avcos) + M_PI) / 2.;
    Psi3 = (atan(avsin3 / avcos3) + M_PI) / 3.;
    Psi4 = (atan(avsin4 / avcos4) + M_PI) / 4.;
    Psi5 = (atan(avsin5 / avcos5) + M_PI) / 5.;
    Psi6 = (atan(avsin6 / avcos6) + M_PI) / 6.;
    eccentricity1 = sqrt(avcos1 * avcos1 + avsin1 * avsin1) / avr1;
    eccentricity2 = sqrt(avcos * avcos + avsin * avsin) / avrSq;
    eccentricity3 = sqrt(avcos3 * avcos3 + avsin3 * avsin3) / avr3;
    eccentricity4 = sqrt(avcos4 * avcos4 + avsin4 * avsin4) / avr4;
    eccentricity5 = sqrt(avcos5 * avcos5 + avsin5 * avsin5) / avr5;
    eccentricity6 = sqrt(avcos6 * avcos6 + avsin6 * avsin6) / avr6;

    double avx2 = avx;
    double avy2 = avy;
    avx = 0.;
    avy = 0.;
    toteps = 0.;

    for (int ix = 0; ix < N; ix++) {
        x = -L / 2. + a * ix - avx2;
        for (int iy = 0; iy < N; iy++) {
            pos = lat->positionFromXY(ix, iy);

            gfactor = computeRunningCouplingGfactor(
                lat, param, pos, N, a, g, c, muZero);

            if (lat->cells[pos]->getEpsilon() * gfactor
                < cutoff)  // this is 1/fm^4, so Lambda_QCD^{-4}
            {
                weight = 0.;
            } else {
                weight =
                    (lat->cells[pos]->getEpsilon() * lat->cells[pos]->getutau()
                     * gfactor);
            }

            y = -L / 2. + a * iy - avy2;
            avx += x * weight;
            avy += y * weight;
            avxSq += x * x * weight;
            avySq += y * y * weight;
            toteps += weight;
        }
    }
    avx /= toteps;
    avy /= toteps;
    avxSq /= toteps;
    avySq /= toteps;
    avrSq /= toteps;
    Rbar = 1. / sqrt(1. / avxSq + 1. / avySq);
    if (it == 1) param->event.psi = Psi2;

    if (doAniso == 0) {
        ofstream foutEcc(ecc_name.c_str(), std::ios::app);
        foutEcc << it * a * param->run.dtau << " " << eccentricity1 << " "
                << Psi1 << " " << eccentricity2 << " " << Psi2 << " "
                << eccentricity3 << " " << Psi3 << " " << eccentricity4 << " "
                << Psi4 << " " << eccentricity5 << " " << Psi5 << " "
                << eccentricity6 << " " << Psi6 << " " << cutoff << " "
                << sqrt(avrSq) << " " << maxX << " " << maxY << " "
                << param->event.b << " " << param->event.Tpp << " "
                << param->event.area << " " << Rbar << " " << avgeden << " "
                << avgQs2AQs2B * hbarc << endl;
        foutEcc.close();
    }

    if (doAniso == 1) {
        stringstream straniso_name;
        straniso_name << "anisotropy" << param->event.eventId << ".dat";
        string aniso_name;
        aniso_name = straniso_name.str();

        ofstream foutAniso(aniso_name.c_str(), std::ios::app);

        double ux, uy, PsiU;
        double unum = 0., uden = 0.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = lat->positionFromXY(ix, iy);
                ux = lat->cells[pos]->getux();
                uy = lat->cells[pos]->getuy();
                unum += sqrt(ux * ux + uy * uy) * sin(2. * atan2(uy, ux));
                uden += sqrt(ux * ux + uy * uy) * cos(2. * atan2(uy, ux));
            }
        }

        PsiU = atan2(unum, uden) / 2.;

        foutAniso << "Psi2=" << Psi2 << ", cos(Psi2)=" << cos(Psi2)
                  << ", sin(Psi2)=" << sin(Psi2) << endl;
        foutAniso << "PsiU=" << PsiU << ", cos(PsiU)=" << cos(PsiU)
                  << ", sin(PsiU)=" << sin(PsiU) << endl;

        // Sample the rotated-tensor anisotropy at Psi = PsiU + k*Pi/8 for
        // k=0..9 (k=0: Psi = PsiU;  // param->event.psi;//-Pi/2.;).
        for (int k = 0; k < 10; ++k) {
            const double Psi = PsiU + static_cast<double>(k) * M_PI / 8.;
            const AnisotropyResult result =
                computeRotatedAnisotropy(lat, N, Psi);
            foutAniso << it * a * param->run.dtau << " "
                      << result.num / result.den << " "
                      << result.num2 / result.den2 << " angle=" << Psi << endl;
        }

        foutAniso.close();
    }
}

void Evolution::readNkt(Parameters *param) {
    messager_ << "[Evolution::readNkt]: Reading n(k_T) from file ";
    messager_.flush("info");
    string Npart, dummy;
    string kt, nkt, Tpp, b;
    double dkt = 0.;
    double dNdeta = 0.;

    // open file

    ifstream fin;
    stringstream strmult_name;
    strmult_name << "multiplicity" << param->event.eventId << ".dat";
    string mult_name;
    mult_name = strmult_name.str();
    fin.open(mult_name.c_str());
    messager_ << "[Evolution::readNkt]: File " << mult_name.c_str() << " ... ";
    messager_.flush("info");

    // open file

    ifstream fin2;
    stringstream strmult_name2;
    strmult_name2 << "NpartdNdy" << param->event.eventId << ".dat";
    string mult_name2;
    mult_name2 = strmult_name2.str();
    fin2.open(mult_name2.c_str());
    messager_ << "[Evolution::readNkt]: File " << mult_name2.c_str() << " ... ";
    messager_.flush("info");

    // read file

    if (fin) {
        for (int ikt = 0; ikt < 100; ikt++) {
            if (!fin.eof()) {
                fin >> dummy;
                fin >> kt;
                fin >> nkt;
                nIn_[ikt] = atof(nkt.c_str());
                fin >> dummy >> Tpp >> b >> Npart;
                if (ikt == 0) dkt = atof(kt.c_str());
                if (ikt == 1) dkt = dkt - atof(kt.c_str());
            }
            messager_ << "[Evolution::readNkt]: " << nIn_[ikt];
            messager_.flush("info");
        }
        fin.close();
        messager_ << "[Evolution::readNkt]:  done.";
        messager_.flush("info");
    } else {
        messager_ << "[Evolution::readNkt]: File " << mult_name.c_str()
                  << " does not exist. Exiting.";
        messager_.flush("error");
        exit(1);
    }

    if (fin2) {
        if (!fin2.eof()) {
            fin2 >> Npart;
            fin2 >> nkt;
            fin2 >> Tpp;
            fin2 >> b;

            dNdeta = atof(nkt.c_str());
        }
        fin2.close();
        messager_ << "[Evolution::readNkt]:  done.";
        messager_.flush("info");
    } else {
        messager_ << "[Evolution::readNkt]: File " << mult_name2.c_str()
                  << " does not exist. Exiting.";
        messager_.flush("error");
        exit(1);
    }

    double m, P;
    m = param->colorCharge.jacobianMass;                           // in GeV
    P = 0.13 + 0.32 * pow(param->collision.sqrtS / 1000., 0.115);  // in GeV
    double dNdeta2;
    dNdeta2 = 0.;

    for (int ik = 0; ik < 100; ik++) {
        if (!param->colorCharge.usePseudoRapidity) {
            dNdeta2 += nIn_[ik] * (ik + 0.5) * dkt * dkt * 2.
                       * M_PI;  // integrate, gives a ik*dkt*2pi*dkt
        } else {
            dNdeta2 +=
                nIn_[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI
                * cosh(param->colorCharge.rapidity())
                / (sqrt(
                    pow(cosh(param->colorCharge.rapidity()), 2.)
                    + m * m / (((ik + 0.5) * dkt) * ((ik + 0.5) * dkt))));
        }
    }

    dNdeta *=
        cosh(param->colorCharge.rapidity())
        / (sqrt(pow(cosh(param->colorCharge.rapidity()), 2.) + m * m / P / P));

    ofstream foutNN("NpartdNdy-mod.dat", std::ios::out);
    foutNN << Npart << " " << dNdeta << " " << dNdeta2 << " "
           << atof(Tpp.c_str()) << " " << atof(b.c_str()) << endl;
    foutNN.close();

    exit(1);
}

void Evolution::hadronizeAndWriteMultiplicity(
    Parameters *param, double a, double dkt, int bins, const double *n,
    double *Nhgsl, int hbins) {
    const double hadronizationStart = ipg::wallSeconds();
    messager_ << "[Evolution::multiplicity]:  Hadronizing ... ";
    messager_.flush("info");
    double z, frac;
    double mypt, kt, Ng;
    int ik;
    const int steps = 6000;
    double dz = 0.95 / static_cast<double>(steps);
    double zValues[steps + 1];
    double zintegrand[steps + 1];
    gsl_interp_accel *zacc = gsl_interp_accel_alloc();
    gsl_spline *zspline = gsl_spline_alloc(gsl_interp_cspline, steps + 1);

    for (int ih = 0; ih <= hbins; ih++) {
        mypt = ih * (20. / static_cast<double>(hbins));

        for (int iz = 0; iz <= steps; iz++) {
            z = 0.05 + iz * dz;
            zValues[iz] = z;

            kt = mypt / z;

            ik = static_cast<int>(
                floor(kt * a / hbarc / dkt - 0.5 + 0.00000001));

            frac = (kt - (ik + 0.5) * dkt / a * hbarc) / (dkt / a * hbarc);

            if (ik + 1 < bins && ik >= 0)
                Ng = ((1. - frac) * n[ik] + frac * n[ik + 1]) * a / hbarc * a
                     / hbarc;  // to make dN/d^2k_T fo k_T in GeV
            else
                Ng = 0.;

            if (!param->colorCharge.usePseudoRapidity) {
                zintegrand[iz] = 1. / (z * z) * Ng * kkp(7, 1, z, kt);
            } else {
                zintegrand[iz] =
                    1. / (z * z) * Ng * 2.
                    * (kkp(1, 1, z, kt) * cosh(param->colorCharge.rapidity())
                           / (sqrt(
                               pow(cosh(param->colorCharge.rapidity()), 2.)
                               + m_pion * m_pion / (mypt * mypt)))
                       + kkp(2, 1, z, kt) * cosh(param->colorCharge.rapidity())
                             / (sqrt(
                                 pow(cosh(param->colorCharge.rapidity()), 2.)
                                 + m_kaon * m_kaon / (mypt * mypt)))
                       + kkp(4, 1, z, kt) * cosh(param->colorCharge.rapidity())
                             / (sqrt(
                                 pow(cosh(param->colorCharge.rapidity()), 2.)
                                 + m_proton * m_proton / (mypt * mypt))));
            }
        }

        zValues[steps] = 1.;  // set exactly 1

        gsl_spline_init(zspline, zValues, zintegrand, steps + 1);
        Nhgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);
    }

    gsl_spline_free(zspline);
    gsl_interp_accel_free(zacc);

    stringstream strmultHad_name;
    strmultHad_name << "multiplicityHadrons" << param->event.eventId << ".dat";
    string multHad_name;
    multHad_name = strmultHad_name.str();

    ofstream foutdNdpt(multHad_name.c_str(), std::ios::out);
    for (int ih = 0; ih <= hbins; ih++) {
        if (ih % 10 == 0)
            foutdNdpt << ih * 20. / static_cast<double>(hbins) << " "
                      << Nhgsl[ih] << " " << 0. << " " << 0. << " "
                      << param->event.Tpp << " " << param->event.b
                      << endl;  // leaving out the L and H ones for now
    }
    foutdNdpt.close();

    messager_ << "[Evolution::multiplicity]:  done.";
    messager_.flush("info");

    ipg::Profiler::instance().add(
        "observables.gluon_multiplicity.hadronization",
        ipg::wallSeconds() - hadronizationStart);
}

int Evolution::multiplicity(
    Lattice *lat, Group *group, Parameters *param, int it) {
    IPG_PROFILE_SCOPE("observables.gluon_multiplicity");
    int N = param->lattice.size;
    int npos, pos;
    double L = param->lattice.L;
    double a = L / N;  // lattice spacing in fm
    double kx, ky, kt2, omega2;
    double g = param->coupling.g;
    int nn[2];
    nn[0] = N;
    nn[1] = N;
    double dtau = param->run.dtau;
    double nkt;
    const int bins = 100;
    double n[bins];   // k_T array
    double E[bins];   // k_T array
    double n2[bins];  // k_T array
    int counter[bins];
    double dkt = 2.83 / static_cast<double>(bins);
    double dNdeta = 0.;
    double dNdeta2 = 0.;
    double dNdetaCut = 0.;
    double dNdetaCut2 = 0.;
    double dEdetaCut = 0.;
    double dEdetaCut2 = 0.;
    double dEdeta = 0.;
    double dEdeta2 = 0.;

    stringstream strNpartdNdy_name;
    strNpartdNdy_name << "NpartdNdy-t" << it * dtau * a << "-"
                      << param->event.eventId << ".dat";
    string NpartdNdy_name;
    NpartdNdy_name = strNpartdNdy_name.str();
    messager_ << "[Evolution::multiplicity]: Measuring multiplicity ... ";
    messager_.flush("info");

    // fix transverse Coulomb gauge
    GaugeFix gaugefix;

    double maxtime;
    if (param->evolution.inverseQsForMaxTime) {
        maxtime = 1. / param->event.averageQs * hbarc;
        messager_ << "[Evolution::multiplicity]: maximal evolution time = "
                  << maxtime << " fm";
        messager_.flush("info");
    } else {
        maxtime = param->evolution.maxTime;  // maxtime is in fm
    }

    int itmax = static_cast<int>(floor(maxtime / (a * dtau) + 1e-10));

    double multiplicityPhaseStart = ipg::wallSeconds();
    gaugefix.fftChi(fft_, lat, group, param, 4000);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.gauge_fix", multiplicityPhaseStart);
    // gauge is fixed

    // E1Storage owns the N*N scratch matrices; E1 is a pointer-per-cell
    // view over it for FFT::fftn's T** interface.
    std::vector<Matrix> E1Storage(N * N, Matrix(0.));
    std::vector<Matrix *> E1(N * N);
    for (int i = 0; i < N * N; i++) {
        E1[i] = &E1Storage[i];
    }
    addPhaseAndRestart(
        "observables.gluon_multiplicity.allocate", multiplicityPhaseStart);

    double c = param->coupling.c;
    double muZero = param->coupling.mu0;

    prepareSpectrumField(lat, param, N, a, g, c, muZero, lat->U, E1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.prepare_E1", multiplicityPhaseStart);

    // do Fourier transforms
    fft_->fftn(E1.data(), E1.data(), nn, 1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.fft_E1", multiplicityPhaseStart);

    for (int ik = 0; ik < bins; ik++) {
        n[ik] = 0.;
        E[ik] = 0.;
        n2[ik] = 0.;
        counter[ik] = 0;
    }

    const int hbins = 2000;
    double Nhgsl[hbins + 1];

    addPhaseAndRestart(
        "observables.gluon_multiplicity.setup_bins", multiplicityPhaseStart);

    accumulateGluonSpectrum(
        param, N, it, dtau, g, a, c, muZero, dkt, bins, E1, true, true, dNdeta,
        dEdeta, n, E, n2, counter);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.spectrum_E1", multiplicityPhaseStart);

    /// -------- 2 ---------

    prepareSpectrumField(lat, param, N, a, g, c, muZero, lat->U2, E1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.prepare_E2", multiplicityPhaseStart);

    fft_->fftn(E1.data(), E1.data(), nn, 1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.fft_E2", multiplicityPhaseStart);

    accumulateGluonSpectrum(
        param, N, it, dtau, g, a, c, muZero, dkt, bins, E1, true, false, dNdeta,
        dEdeta, n, E, n2, counter);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.spectrum_E2", multiplicityPhaseStart);

    /// ------3 --------

    prepareSpectrumField(lat, param, N, a, g, c, muZero, lat->Ux2, E1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.prepare_pi", multiplicityPhaseStart);

    // do Fourier transforms
    fft_->fftn(E1.data(), E1.data(), nn, 1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.fft_pi", multiplicityPhaseStart);

    accumulateGluonSpectrum(
        param, N, it, dtau, g, a, c, muZero, dkt, bins, E1, false, false,
        dNdeta, dEdeta, n, E, n2, counter);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.spectrum_pi", multiplicityPhaseStart);

    double m, P;
    m = param->colorCharge.jacobianMass;                           // in GeV
    P = 0.13 + 0.32 * pow(param->collision.sqrtS / 1000., 0.115);  // in GeV

    for (int ik = 0; ik < bins; ik++) {
        if (counter[ik] > 0) {
            n[ik] = n[ik] / static_cast<double>(counter[ik]);
            E[ik] = E[ik] / static_cast<double>(counter[ik]);
            // integrate, gives a ik*dkt*2pi*dkt
            const double weight =
                computeMultiplicityBinWeight(param, m, ik, dkt, a);
            dNdeta2 += n[ik] * weight;
            dEdeta2 += E[ik] * weight;
            if (ik * dkt / a * hbarc > 3.)  //
            {
                dNdetaCut += n[ik] * weight;
                dEdetaCut += E[ik] * weight;
            }
            if (ik * dkt / a * hbarc > 6.)  // large cut
            {
                dNdetaCut2 += n[ik] * weight;
                dEdetaCut2 += E[ik] * weight;
            }
        }
    }

    addPhaseAndRestart(
        "observables.gluon_multiplicity.bin_postprocess",
        multiplicityPhaseStart);

    // compute hadrons using fragmentation function
    if (it == itmax && param->output.writeOutputs == 3) {
        hadronizeAndWriteMultiplicity(param, a, dkt, bins, n, Nhgsl, hbins);
        multiplicityPhaseStart = ipg::wallSeconds();
    }

    if (!param->colorCharge.usePseudoRapidity && param->run.MPIRank == 0) {
        messager_ << "[Evolution::multiplicity]: dN/dy 1 = " << dNdeta
                  << ", dE/dy 1 = " << dEdeta;
        messager_.flush("info");
        messager_ << "[Evolution::multiplicity]: dN/dy 2 = " << dNdeta2
                  << ", dE/dy 2 = " << dEdeta2;
        messager_.flush("info");
        messager_ << "[Evolution::multiplicity]: gluon <p_T> = "
                  << dEdeta / dNdeta;
        messager_.flush("info");
    } else if (param->colorCharge.usePseudoRapidity) {
        m = param->colorCharge.jacobianMass;                           // in GeV
        P = 0.13 + 0.32 * pow(param->collision.sqrtS / 1000., 0.115);  // in GeV
        dNdeta *= cosh(param->colorCharge.rapidity())
                  / (sqrt(
                      pow(cosh(param->colorCharge.rapidity()), 2.)
                      + m * m / (P * P)));
        dEdeta *= cosh(param->colorCharge.rapidity())
                  / (sqrt(
                      pow(cosh(param->colorCharge.rapidity()), 2.)
                      + m * m / (P * P)));

        if (param->run.MPIRank == 0) {
            messager_ << "[Evolution::multiplicity]: dN/deta 1 = " << dNdeta
                      << ", dE/deta 1 = " << dEdeta;
            messager_.flush("info");
            messager_ << "[Evolution::multiplicity]: dN/deta 2 = " << dNdeta2
                      << ", dE/deta 2 = " << dEdeta2;
            messager_.flush("info");
            messager_ << "[Evolution::multiplicity]: dN/deta_cut 1 = "
                      << dNdetaCut;
            messager_.flush("info");
            messager_ << "[Evolution::multiplicity]: dN/deta_cut 2 = "
                      << dNdetaCut2;
            messager_.flush("info");
            messager_ << "[Evolution::multiplicity]: gluon <p_T> = "
                      << dEdeta / dNdeta;
            messager_.flush("info");
        }
    }

    addPhaseAndRestart(
        "observables.gluon_multiplicity.report", multiplicityPhaseStart);

    if (dNdeta == 0.) {
        messager_ << "[Evolution::multiplicity]: No collision happened on "
                     "rank "
                  << param->run.MPIRank
                  << ". Restarting with new random number...";
        messager_.flush("warning");
        addPhaseAndRestart(
            "observables.gluon_multiplicity.cleanup", multiplicityPhaseStart);
        return 0;
    }

    if (it == itmax) {
        ofstream foutNN(NpartdNdy_name.c_str(), std::ios::out);
        foutNN << param->event.Npart << " " << dNdeta << " " << param->event.Tpp
               << " " << param->event.b << " " << dEdeta << " "
               << param->run.randomSeed << " "
               << "N/A"
               << " "
               << "N/A"
               << " "
               << "N/A"
               << " " << dNdetaCut << " " << dEdetaCut << " " << dNdetaCut2
               << " " << dEdetaCut2 << " "
               << computeRunningCouplingGfactorFromScale(
                      g, muZero, c, param->coupling.LambdaQCD,
                      param->coupling.nFlavors,
                      param->coupling.runningCouplingQsFactor
                          * param->event.averageQs)
               << endl;
        foutNN.close();
        writeGluonMultiplicityTarget(
            param, it, a, dtau, dNdeta, dNdeta2, dEdeta, dEdeta2, dNdetaCut,
            dEdetaCut, dNdetaCut2, dEdetaCut2, n, E, counter, bins, dkt);
    }
    addPhaseAndRestart("output.gluon_multiplicity", multiplicityPhaseStart);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.cleanup", multiplicityPhaseStart);

    messager_ << "[Evolution::multiplicity]:  done.";
    messager_.flush("info");
    param->event.success = 1;
    return 1;
}

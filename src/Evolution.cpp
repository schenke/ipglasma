// Evolution.cpp is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Evolution.h"

#include <algorithm>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "Eccentricity.h"
#include "EnergyMomentumTensor.h"
#include "Instrumentation.h"
#include "MyEigen.h"
#include "PhysConst.h"
#include "SU3.h"

using PhysConst::hbarc;
using PhysConst::Nc;

using std::endl;
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
                Eccentricity::compute(
                    lat, param, it, param->output.eccentricityCutoff, 0);

                success = gluonMultiplicity_.compute(lat, group, param, it);
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
        // Eccentricity::compute() weights by epsilon * u^tau, which only the
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

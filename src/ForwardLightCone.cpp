// ForwardLightCone.cpp is part of the IP-Glasma solver.

#include "ForwardLightCone.h"

#include <cmath>
#include <cstdint>
#include <vector>

#include "PhysConst.h"
#include "gsl/gsl_linalg.h"

using PhysConst::Nc2m1;
using std::complex;

namespace {

/**
 * SplitMix64 finalizer/counter step, used only to construct a
 * deterministic, stateless retry stream for the forward-lightcone
 * solver (see forwardLightconeRetrySeed()/deterministicRetryGaussian()).
 * \param[in] x Input state.
 * \return The mixed 64-bit output.
 */
inline std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

/**
 * Derives a deterministic seed for
 * ForwardLightCone::findU()'s restart stream, unique per
 * cell/direction/event/run without touching the shared Random state
 * (so retries are reproducible even though they're computed in
 * parallel across cells). Chains splitmix64() over \p runSeed and each
 * of `eventId`/`pos`/`direction` in turn.
 * \param[in] runSeed Run-wide seed (typically `param->run.randomSeed`).
 * \param[in] eventId Current event id.
 * \param[in] pos Flat cell index.
 * \param[in] direction Distinguishes the \f$x\f$ and \f$y\f$ link
 * solves at the same cell (e.g. `0`/`1`).
 * \return A 64-bit seed to drive deterministicRetryGaussian() for this
 * cell/direction/event/run.
 */
inline std::uint64_t forwardLightconeRetrySeed(
    std::uint64_t runSeed, int eventId, int pos, int direction) {
    std::uint64_t key = splitmix64(runSeed);
    key = splitmix64(
        key
        ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(eventId))
           + 0xD1B54A32D192ED03ULL));
    key = splitmix64(
        key
        ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(pos))
           + 0x94D049BB133111EBULL));
    key = splitmix64(
        key
        ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(direction))
           + 0xBF58476D1CE4E5B9ULL));
    return key;
}

/**
 * Draws one standard-normal sample from a deterministic, stateless
 * stream keyed by \p seed and \p drawIndex (see
 * forwardLightconeRetrySeed()), via the Box-Muller transform applied
 * to two independent SplitMix64 counters. No shared Random state is
 * touched, so this is safe to call from multiple threads/cells at
 * once and reproducible given the same `seed`/`drawIndex`.
 * \param[in] seed Stream seed (see forwardLightconeRetrySeed()).
 * \param[in] drawIndex Draw index within the stream; increment for
 * each successive sample requested with the same \p seed.
 * \return One standard-normal (\f$\mu=0\f$, \f$\sigma=1\f$) sample.
 */
inline double deterministicRetryGaussian(
    std::uint64_t seed, std::uint64_t drawIndex) {
    // Build two open-interval uniform doubles from independent SplitMix64
    // counters, then use Box-Muller.  No shared RNG state is touched.
    const std::uint64_t counter = 2ULL * drawIndex;
    const std::uint64_t r1 =
        splitmix64(seed + 0x9E3779B97F4A7C15ULL * (counter + 1ULL));
    const std::uint64_t r2 =
        splitmix64(seed + 0x9E3779B97F4A7C15ULL * (counter + 2ULL));

    constexpr double invTwo53 = 1.0 / 9007199254740992.0;
    const double u1 = (static_cast<double>(r1 >> 11) + 0.5) * invTwo53;
    const double u2 = (static_cast<double>(r2 >> 11) + 0.5) * invTwo53;
    constexpr double twoPi = 6.283185307179586476925286766559;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(twoPi * u2);
}

}  // namespace

void ForwardLightCone::solveAxb(
    double *Jab, double *Fa, std::vector<double> &xvec) {
    gsl_matrix_view m = gsl_matrix_view_array(Jab, Nc2m1, Nc2m1);
    gsl_vector_view c = gsl_vector_view_array(Fa, Nc2m1);
    gsl_vector *x = gsl_vector_alloc(Nc2m1);

    int s;
    gsl_permutation *p = gsl_permutation_alloc(Nc2m1);
    gsl_linalg_LU_decomp(&m.matrix, p, &s);
    gsl_linalg_LU_solve(&m.matrix, p, &c.vector, x);
    gsl_permutation_free(p);

    for (int i = 0; i < Nc2m1; i++) {
        xvec[i] = gsl_vector_get(x, i);
    }

    gsl_vector_free(x);
}

void ForwardLightCone::initialize(Lattice *lat, Parameters *param) {
    messager_.info(
        "[ForwardLightCone::initialize]: Finding fields in forward "
        "lightcone...");
    const int N = param->lattice.size;
    const int N2 = N * N;
#pragma omp parallel
    {
        sanitizeU(lat, N2);

        LinkScratch linkScratch;
        computeLinksTeam(lat, N2, linkScratch);

        UScratch uScratch;
        computeUxUyTeam(lat, param, N2, uScratch);

        // compute initial electric field
        ElectricFieldScratch electricScratch;
        // with minus ax, ay
        computeElectricFieldTeam(
            lat, N2, lat->posmX, lat->posmY, lat->U, electricScratch);
        // with plus ax, ay
        computeElectricFieldTeam(
            lat, N2, lat->pospX, lat->pospY, lat->U2, electricScratch);

        // compute the plaquette
        PlaquetteScratch plaquetteScratch;
        computePlaquetteTeam(lat, N2, plaquetteScratch);

        computePiTeam(lat, param, N2);

        resetFieldsTeam(lat, N2);
    }  // omp block
}

void ForwardLightCone::sanitizeU(Lattice *lat, int N2) {
// compute Ux(3) Uy(3) after the collision
#pragma omp for
    for (int pos = 0; pos < N2; pos++) {
        // loops over all cells
        auto checkU = lat->U[pos].trace();
        if (checkU != checkU) {
            lat->U[pos] = (one_);
        }

        checkU = lat->U2[pos].trace();
        if (checkU != checkU) {
            lat->U2[pos] = (one_);
        }
    }
}

void ForwardLightCone::computeLinksTeam(
    Lattice *lat, int N2, LinkScratch &scratch) {
#pragma omp for
    for (int pos = 0; pos < N2; pos++) {
        // loops over all cells
        scratch.UDx = lat->U[lat->pospX[pos]];
        scratch.UDx.conjg();
        lat->Ux1[pos] = (lat->U[pos] * scratch.UDx);

        scratch.UDy = lat->U[lat->pospY[pos]];
        scratch.UDy.conjg();
        lat->Uy1[pos] = (lat->U[pos] * scratch.UDy);

        scratch.UDx = lat->U2[lat->pospX[pos]];
        scratch.UDx.conjg();
        lat->Ux2[pos] = (lat->U2[pos] * scratch.UDx);

        scratch.UDy = lat->U2[lat->pospY[pos]];
        scratch.UDy.conjg();
        lat->Uy2[pos] = (lat->U2[pos] * scratch.UDy);
    }
}

// from Ux(1,2) and Uy(1,2) compute Ux(3) and Uy(3):
void ForwardLightCone::computeUxUyTeam(
    Lattice *lat, Parameters *param, int N2, UScratch &scratch) {
#pragma omp for
    for (int pos = 0; pos < N2; pos++) {
        // loops over all cells
        scratch.UDx1 = lat->Ux1[pos];
        scratch.UDx2 = lat->Ux2[pos];
        const std::uint64_t retrySeedX = forwardLightconeRetrySeed(
            param->run.randomSeed, param->event.eventId, pos, 0);
        bool status =
            findU(scratch.UDx1, scratch.UDx2, scratch.temp2, retrySeedX);
        lat->Ux[pos] = (scratch.temp2);
        if (!status) {
            // Inside an omp parallel/for region: use a fresh,
            // stack-local instance rather than sharing messager_.
            PrettyOstream localMessager;
            localMessager << "[ForwardLightCone::initialize]: Failed "
                             "to converge finding Ux in "
                             "the forward lightcone at pos x = "
                          << pos / param->lattice.size
                          << " y = " << pos % param->lattice.size;
            localMessager.flush("warning");
        }

        scratch.UDy1 = lat->Uy1[pos];
        scratch.UDy2 = lat->Uy2[pos];
        const std::uint64_t retrySeedY = forwardLightconeRetrySeed(
            param->run.randomSeed, param->event.eventId, pos, 1);
        status = findU(scratch.UDy1, scratch.UDy2, scratch.temp2, retrySeedY);
        lat->Uy[pos] = (scratch.temp2);
        if (!status) {
            // Inside an omp parallel/for region: use a fresh,
            // stack-local instance rather than sharing messager_.
            PrettyOstream localMessager;
            localMessager << "[ForwardLightCone::initialize]: Failed "
                             "to converge finding Uy in "
                             "the forward lightcone at pos x = "
                          << pos / param->lattice.size
                          << " y = " << pos % param->lattice.size;
            localMessager.flush("warning");
        }
    }
}

void ForwardLightCone::computeElectricFieldTeam(
    Lattice *lat, int N2, const std::vector<int> &neighborX,
    const std::vector<int> &neighborY, std::vector<Matrix> &outputField,
    ElectricFieldScratch &scratch) {
#pragma omp for
    for (int pos = 0; pos < N2; pos++) {
        // x part in sum:
        scratch.Ux1mUx2 = lat->Ux1[pos] - lat->Ux2[pos];
        scratch.UDx1 = lat->Ux1[pos];
        scratch.UDx1.conjg();
        scratch.UDx2 = lat->Ux2[pos];
        scratch.UDx2.conjg();
        scratch.UDx1mUDx2 = scratch.UDx1 - scratch.UDx2;

        scratch.Ux = lat->Ux[pos];
        scratch.UDx = scratch.Ux;
        scratch.UDx.conjg();

        scratch.temp2 = scratch.Ux1mUx2 * scratch.UDx - scratch.Ux1mUx2
                        - scratch.Ux * scratch.UDx1mUDx2 + scratch.UDx1mUDx2;

        scratch.Ux1mUx2 = lat->Ux1[neighborX[pos]] - lat->Ux2[neighborX[pos]];
        scratch.UDx1 = lat->Ux1[neighborX[pos]];
        scratch.UDx1.conjg();
        scratch.UDx2 = lat->Ux2[neighborX[pos]];
        scratch.UDx2.conjg();
        scratch.UDx1mUDx2 = scratch.UDx1 - scratch.UDx2;

        scratch.Ux = lat->Ux[neighborX[pos]];
        scratch.UDx = scratch.Ux;
        scratch.UDx.conjg();

        scratch.temp2 = scratch.temp2 - scratch.UDx * scratch.Ux1mUx2
                        + scratch.Ux1mUx2 + scratch.UDx1mUDx2 * scratch.Ux
                        - scratch.UDx1mUDx2;

        // y part in sum
        scratch.Uy1mUy2 = lat->Uy1[pos] - lat->Uy2[pos];
        scratch.UDy1 = lat->Uy1[pos];
        scratch.UDy1.conjg();
        scratch.UDy2 = lat->Uy2[pos];
        scratch.UDy2.conjg();
        scratch.UDy1mUDy2 = scratch.UDy1 - scratch.UDy2;

        scratch.Uy = lat->Uy[pos];
        scratch.UDy = scratch.Uy;
        scratch.UDy.conjg();

        // y part of the sum:
        scratch.temp2 = scratch.temp2 + scratch.Uy1mUy2 * scratch.UDy
                        - scratch.Uy1mUy2 - scratch.Uy * scratch.UDy1mUDy2
                        + scratch.UDy1mUDy2;

        scratch.Uy1mUy2 = lat->Uy1[neighborY[pos]] - lat->Uy2[neighborY[pos]];
        scratch.UDy1 = lat->Uy1[neighborY[pos]];
        scratch.UDy1.conjg();
        scratch.UDy2 = lat->Uy2[neighborY[pos]];
        scratch.UDy2.conjg();
        scratch.UDy1mUDy2 = scratch.UDy1 - scratch.UDy2;

        scratch.Uy = lat->Uy[neighborY[pos]];
        scratch.UDy = scratch.Uy;
        scratch.UDy.conjg();

        scratch.temp2 = scratch.temp2 - scratch.UDy * scratch.Uy1mUy2
                        + scratch.Uy1mUy2 + scratch.UDy1mUDy2 * scratch.Uy
                        - scratch.UDy1mUDy2;

        outputField[pos] = ((1. / 8.) * scratch.temp2);
    }
}

void ForwardLightCone::computePlaquetteTeam(
    Lattice *lat, int N2, PlaquetteScratch &scratch) {
#pragma omp for
    for (int pos = 0; pos < N2; pos++) {
        scratch.UDx = lat->Ux[lat->pospY[pos]];
        scratch.UDy = lat->Uy[pos];

        scratch.UDx.conjg();
        scratch.UDy.conjg();

        scratch.Uplaq =
            lat->Ux[pos]
            * (lat->Uy[lat->pospX[pos]] * (scratch.UDx * scratch.UDy));
        lat->Uy1[pos] = (scratch.Uplaq);
    }
}

void ForwardLightCone::computePiTeam(Lattice *lat, Parameters *param, int N2) {
#pragma omp for
    for (int pos = 0; pos < N2; pos++) {
        // this is pi in lattice units as needed for the evolution. (later,
        // the a^4 gives the right units for the energy density
        lat->Ux2[pos] =
            (complex<double>(0., -2. / param->coupling.g) * (lat->U[pos]));
        // factor -2 because I have A^eta (note the 1/8 before)
        // but want \pi (E^z).
    }
}

void ForwardLightCone::resetFieldsTeam(Lattice *lat, int N2) {
    const Matrix zero(0.);
#pragma omp for
    for (int pos = 0; pos < N2; pos++) {
        lat->U[pos] = (zero);
        lat->U2[pos] = (zero);
        lat->Uy2[pos] = (zero);

        // reset the Ux1 to be used for other purposes later
        lat->Ux1[pos] = (one_);
    }
}

double ForwardLightCone::computeResidual(
    const Matrix &U1pU2, const Matrix &U1pU2dagger, const Matrix &Usol,
    const Matrix &Usoldagger, const std::vector<complex<double>> &traceCache,
    double *Fa) {
    Matrix Mtemp = U1pU2 * Usoldagger - Usol * U1pU2dagger;
    double Fzero = 0.;
    for (int ai = 0; ai < Nc2m1; ai++) {
        complex<double> traceLoc =
            Mtemp.traceOfProductOfMatrix(group_->getT(ai), Mtemp);
        // minus trace if temp gives -F_ai
        auto traceRes = (-1.) * (traceCache[ai] + traceLoc);
        Fa[ai] = imag(traceRes);
        Fzero += std::abs(Fa[ai]);
    }
    return Fzero;
}

void ForwardLightCone::computeJacobian(
    const Matrix &U0, const Matrix &U1pU2, const Matrix &Usoldagger,
    std::vector<Matrix> &MtempArr, std::vector<double> &alpha, double *Jab) {
    Matrix Mtemp(0.);

    // numerical formula
    bool JabGood = true;
    for (int bi = 0; bi < Nc2m1; bi++) {
        double dalpha_bi =
            (std::max(0.001, std::min(10., 0.01 * std::abs(alpha[bi]))));
        alpha[bi] = alpha[bi] + dalpha_bi;
        Mtemp = Matrix::fromAlgebraExponent(alpha) * U0;
        Mtemp.conjg();
        Mtemp = U1pU2 * (Mtemp - Usoldagger);
        double Mcheck = 0.;
        for (int ai = 0; ai < Nc2m1; ai++) {
            int countMe = ai * Nc2m1 + bi;
            complex<double> traceLoc =
                Mtemp.traceOfProductOfMatrix(group_->getT(ai), Mtemp);
            Jab[countMe] = 2. * imag(traceLoc) / dalpha_bi;
            Mcheck += std::abs(Jab[countMe]);
            if (Mcheck < 1e-15) {
                // avoid matrix to be singular
                JabGood = false;
            }
        }
        alpha[bi] = alpha[bi] - dalpha_bi;
    }
    if (!JabGood) {
        // analytical approximated formula
        for (int bi = 0; bi < Nc2m1; bi++) {
            Mtemp = group_->getT(bi) * Usoldagger;
            for (int ai = 0; ai < Nc2m1; ai++) {
                int countMe = ai * Nc2m1 + bi;
                complex<double> traceLoc =
                    Mtemp.traceOfProductOfMatrix(MtempArr[ai], Mtemp);
                auto traceRes = -2. * real(traceLoc);
                Jab[countMe] = traceRes;
            }
        }
    }
}

bool ForwardLightCone::findU(
    Matrix &U1, Matrix &U2, Matrix &Usol, std::uint64_t retrySeed) {
    const int maxIterations = 2000;
    const int maxRetrys = 200;

    Matrix U0 = U1 * U2;
    Matrix U1pU2 = U1 + U2;
    Matrix U1pU2dagger = U1pU2;
    U1pU2dagger.conjg();

    Matrix Mtemp(0.);
    std::vector<Matrix> MtempArr;
    MtempArr.resize(Nc2m1);
    std::vector<complex<double>> traceCache(Nc2m1, 0.);
    for (int ai = 0; ai < Nc2m1; ai++) {
        Mtemp = group_->getT(ai) * (U1pU2 - U1pU2dagger);
        traceCache[ai] = Mtemp.trace();
        MtempArr[ai] = group_->getT(ai) * U1pU2;
    }

    // JabData/FaData own the storage; Jab/Fa are raw-pointer views into it
    // for solveAxb's GSL interface (gsl_matrix_view_array/
    // gsl_vector_view_array need a raw contiguous buffer, not a
    // std::vector).
    std::vector<double> JabData(Nc2m1 * Nc2m1);
    std::vector<double> FaData(Nc2m1);
    double *Jab = JabData.data();
    double *Fa = FaData.data();

    double Fzero = 10.;
    double FzeroMin = 1e6;
    Matrix UsolBestEst(1.);

    // set up initial guess
    std::vector<double> alpha(Nc2m1, 0.);  // solution
    std::vector<double> Dalpha(Nc2m1, 0.);
    Usol = Matrix::fromAlgebraExponent(alpha) * U0;
    Matrix Usoldagger = Usol;
    Usoldagger.conjg();

    int iter = 0;
    int nRestart = 0;
    std::uint64_t retryDraw = 0;
    auto nextRetryGaussian = [&]() {
        return deterministicRetryGaussian(retrySeed, retryDraw++);
    };
    while (Fzero > 1e-6 && iter < maxIterations && nRestart < maxRetrys) {
        iter++;

        Fzero = computeResidual(
            U1pU2, U1pU2dagger, Usol, Usoldagger, traceCache, Fa);

        computeJacobian(U0, U1pU2, Usoldagger, MtempArr, alpha, Jab);

        solveAxb(Jab, Fa, Dalpha);

        bool DalphaCheck = true;
        for (int ai = 0; ai < Nc2m1; ai++) {
            if (std::abs(Dalpha[ai]) > 1e5) {
                DalphaCheck = false;
                break;
            }
        }
        if (DalphaCheck) {
            for (int ai = 0; ai < Nc2m1; ai++) {
                alpha[ai] = alpha[ai] + Dalpha[ai];
            }
        } else {
            for (int ai = 0; ai < Nc2m1; ai++) {
                alpha[ai] = nextRetryGaussian();
            }
        }
        Usol = Matrix::fromAlgebraExponent(alpha) * U0;
        Usoldagger = Usol;
        Usoldagger.conjg();

        if (Fzero < FzeroMin) {
            FzeroMin = Fzero;
            UsolBestEst = Usol;
        }
        if (iter == maxIterations) {
            for (int ai = 0; ai < Nc2m1; ai++) {
                alpha[ai] = nextRetryGaussian();
            }
            Usol = Matrix::fromAlgebraExponent(alpha) * U0;
            Usoldagger = Usol;
            Usoldagger.conjg();

            nRestart++;
            iter = 0;
        }
    }
    bool success = true;
    if (nRestart == maxRetrys) {
        // findU() is called from inside an omp
        // parallel/for region (initialize), so use a
        // fresh, stack-local instance rather than sharing messager_.
        PrettyOstream localMessager;
        localMessager << "[ForwardLightCone::findU]: Did not converge, "
                      << "Fzero: " << FzeroMin;
        localMessager.flush("warning");
        Usol = UsolBestEst;  // return the best estimate
        success = false;
    }
    return (success);
}

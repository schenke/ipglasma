// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include "ForwardLightCone.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "Instrumentation.h"
#include "PhysConst.h"
#include "gsl/gsl_linalg.h"

using PhysConst::Nc2m1;
using std::complex;

namespace {

/**
 * Whether all nine elements of \p M are finite.
 * \param[in] M The matrix.
 * \return `false` if an element has a NaN or infinite part.
 */
bool isFinite(const Matrix &M) {
    for (int i = 0; i < 9; i++) {
        if (!std::isfinite(M(i).real()) || !std::isfinite(M(i).imag())) {
            return false;
        }
    }
    return true;
}

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

bool ForwardLightCone::solveAxb(
    double *Jab, double *Fa, std::vector<double> &xvec) {
    gsl_matrix_view m = gsl_matrix_view_array(Jab, Nc2m1, Nc2m1);
    gsl_vector_view c = gsl_vector_view_array(Fa, Nc2m1);
    gsl_vector_view x = gsl_vector_view_array(xvec.data(), Nc2m1);

    int s;
    gsl_permutation *p = gsl_permutation_alloc(Nc2m1);
    gsl_linalg_LU_decomp(&m.matrix, p, &s);
    // with a zero pivot gsl_linalg_LU_solve calls the GSL error handler,
    // which aborts the program by default (and JIMWLK turns it off, so
    // the solution would be left undefined)
    bool solved = true;
    for (int i = 0; i < Nc2m1; i++) {
        if (gsl_matrix_get(&m.matrix, i, i) == 0.) solved = false;
    }
    if (solved) {
        solved = gsl_linalg_LU_solve(&m.matrix, p, &c.vector, &x.vector)
                 == GSL_SUCCESS;
    }
    gsl_permutation_free(p);
    if (!solved) std::fill(xvec.begin(), xvec.end(), 0.);
    return solved;
}

void ForwardLightCone::initialize(Lattice *lat, Parameters *param) {
    IPG_PROFILE_SCOPE("initialization.forward_light_cone");
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

        // compute initial electric field (with minus ax, ay)
        ElectricFieldScratch electricScratch;
        computeElectricFieldTeam(
            lat, N2, lat->posmX, lat->posmY, lat->U, electricScratch);

        computePiTeam(lat, param, N2);

        resetFieldsTeam(lat, N2);
    }  // omp block
}

void ForwardLightCone::sanitizeU(Lattice *lat, int N2) {
#pragma omp for
    for (int pos = 0; pos < N2; pos++) {
        if (!isFinite(lat->U[pos])) lat->U[pos] = one_;
        if (!isFinite(lat->U2[pos])) lat->U2[pos] = one_;
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

void ForwardLightCone::computePiTeam(Lattice *lat, Parameters *param, int N2) {
#pragma omp for
    for (int pos = 0; pos < N2; pos++) {
        // this is pi in lattice units as needed for the evolution. (later,
        // the a^4 gives the right units for the energy density
        // factor -2 because I have A^eta (note the 1/8 before)
        // but want \pi (E^z).
        Matrix pi = complex<double>(0., -2. / param->coupling.g) * lat->U[pos];
        // The lattice expression has a trace part, a U(1) component that
        // vanishes in the continuum limit (a commutator there): keep pi
        // traceless, like the transverse electric fields. The evolution
        // keeps it traceless then.
        const complex<double> traceThird = pi.trace() / 3.;
        for (int i = 0; i < 3; i++) pi.set(i, i, pi(i, i) - traceThird);
        lat->Ux2[pos] = pi;
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

std::vector<complex<double>> ForwardLightCone::residualTraceCache(
    const Matrix &U1pU2) const {
    Matrix U1pU2dagger = U1pU2;
    U1pU2dagger.conjg();
    std::vector<complex<double>> traceCache(Nc2m1, 0.);
    for (int ai = 0; ai < Nc2m1; ai++) {
        traceCache[ai] = (group_->getT(ai) * (U1pU2 - U1pU2dagger)).trace();
    }
    return traceCache;
}

double ForwardLightCone::computeResidual(
    const Matrix &M, const std::vector<complex<double>> &traceCache,
    double *Fa) {
    // M - M^dagger = (U1 + U2) U^dagger - U (U1 + U2)^dagger
    Matrix Mdagger = M;
    Mdagger.conjg();
    Matrix Mtemp = M - Mdagger;
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

void ForwardLightCone::computeJacobian(const Matrix &M, double *Jab) {
    // dF_a/d delta_b for U -> e^{i delta_b t_b} U is
    // 2 Re tr[t_a (U1 + U2) U^dagger t_b]; solveAxb() solves J dx = F, so
    // J = -dF/d delta. t_a M is computed once per a, and the trace summed
    // in the order of su3::traceABC(t_a, M, t_b), which gives the same
    // numbers with a third of the operations.
    const complex<double> *B = M.data();
    complex<double> tM[9];
    for (int ai = 0; ai < Nc2m1; ai++) {
        const complex<double> *A = group_->getT(ai).data();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                tM[3 * i + j] = A[3 * i] * B[j] + A[3 * i + 1] * B[3 + j]
                                + A[3 * i + 2] * B[6 + j];
            }
        }
        for (int bi = 0; bi < Nc2m1; bi++) {
            const complex<double> *C = group_->getT(bi).data();
            complex<double> tr = 0.;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) tr += tM[3 * i + j] * C[3 * j + i];
            }
            Jab[ai * Nc2m1 + bi] = -2. * tr.real();
        }
    }
}

bool ForwardLightCone::findU(
    Matrix &U1, Matrix &U2, Matrix &Usol, std::uint64_t retrySeed) {
    const int maxIterations = 2000;
    const int maxRetrys = 200;

    // the iteration needs finite links; sanitizeU() has replaced non-finite
    // Wilson lines already, so this only guards other callers
    if (!isFinite(U1) || !isFinite(U2)) {
        Usol = one_;
        return false;
    }

    const Matrix U0 = U1 * U2;
    const Matrix U1pU2 = U1 + U2;
    const std::vector<complex<double>> traceCache = residualTraceCache(U1pU2);

    // JabData/FaData own the storage; Jab/Fa are raw-pointer views into it
    // for solveAxb's GSL interface (gsl_matrix_view_array/
    // gsl_vector_view_array need a raw contiguous buffer, not a
    // std::vector).
    std::vector<double> JabData(Nc2m1 * Nc2m1);
    std::vector<double> FaData(Nc2m1);
    double *Jab = JabData.data();
    double *Fa = FaData.data();

    // the lowest residual so far and the estimate that has it
    double FzeroMin = 1e6;
    Matrix UsolBestEst(1.);

    // Newton's method on the group: start from U1 U2 and rotate the
    // estimate from the left, U -> e^{i Dalpha_a t_a} U, with the exact
    // Jacobian of that rotation
    std::vector<double> Dalpha(Nc2m1, 0.);
    Usol = U0;
    Matrix Usoldagger = Usol;
    Usoldagger.conjg();

    std::uint64_t retryDraw = 0;
    // a fresh start e^{i alpha_a t_a} U1 U2 with random alpha
    auto restart = [&]() {
        std::vector<double> alpha(Nc2m1);
        for (double &a : alpha) {
            a = deterministicRetryGaussian(retrySeed, retryDraw++);
        }
        Usol = Matrix::fromAlgebraExponent(alpha) * U0;
        Usoldagger = Usol;
        Usoldagger.conjg();
    };

    int iter = 0;
    int nRestart = 0;
    bool converged = false;
    while (nRestart < maxRetrys) {
        iter++;

        // shared by the residual and the Jacobian
        const Matrix M = U1pU2 * Usoldagger;
        const double Fzero = computeResidual(M, traceCache, Fa);
        if (Fzero < FzeroMin) {
            FzeroMin = Fzero;
            UsolBestEst = Usol;
        }
        // false for a NaN residual
        converged = Fzero <= 1e-6;

        computeJacobian(M, Jab);
        bool stepIsGood = solveAxb(Jab, Fa, Dalpha);
        for (int ai = 0; ai < Nc2m1; ai++) {
            // also catches a NaN step
            if (!(std::abs(Dalpha[ai]) <= 1e5)) stepIsGood = false;
        }
        if (stepIsGood) {
            // also once more after convergence, which refines the solution
            Usol = Matrix::fromAlgebraExponent(Dalpha) * Usol;
            Usoldagger = Usol;
            Usoldagger.conjg();
        } else if (!converged) {
            restart();
        }
        // a converged estimate is kept, not replaced by a restart
        if (converged) break;

        if (iter == maxIterations) {
            restart();
            nRestart++;
            iter = 0;
        }
    }
    if (!converged) {
        // findU() is called from inside an omp
        // parallel/for region (initialize), so use a
        // fresh, stack-local instance rather than sharing messager_.
        PrettyOstream localMessager;
        localMessager << "[ForwardLightCone::findU]: Did not converge, "
                      << "Fzero: " << FzeroMin;
        localMessager.flush("warning");
        Usol = UsolBestEst;  // return the best estimate
    }
    return converged;
}

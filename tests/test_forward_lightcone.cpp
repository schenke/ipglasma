// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <vector>

#include "ForwardLightCone.h"
#include "Group.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "doctest.h"
#include "test_helpers.h"

TEST_CASE(
    "ForwardLightCone::sanitizeU replaces a U/U2 with a NaN or infinite "
    "element with identity") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    Lattice lat(&param, N);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    lat.U[5].set(0, 0, std::complex<double>(nan, 0.0));
    lat.U2[7].set(4, std::complex<double>(0.0, nan));
    // off the diagonal, where the trace does not see it
    lat.U[9].set(0, 1, std::complex<double>(nan, 0.0));
    lat.U2[11].set(2, 1, std::complex<double>(inf, 0.0));
    lat.U[13] = makeTestMatrix(13);

    Group group;
    ForwardLightCone flc(&group);
    flc.sanitizeU(&lat, N * N);

    const Matrix identity(1.0);
    CHECK(matricesClose(lat.U[5], identity, 1e-14));
    CHECK(matricesClose(lat.U2[7], identity, 1e-14));
    CHECK(matricesClose(lat.U[9], identity, 1e-14));
    CHECK(matricesClose(lat.U2[11], identity, 1e-14));
    // finite cells are not touched
    CHECK(matricesClose(lat.U[0], identity, 1e-14));
    CHECK(matricesClose(lat.U[13], makeTestMatrix(13), 1e-15));
}

TEST_CASE(
    "ForwardLightCone::computeLinksTeam computes Ux1/Uy1/Ux2/Uy2 as "
    "U * conjg(U at the +x/+y neighbor)") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    Lattice lat(&param, N);
    for (int pos = 0; pos < N * N; ++pos) {
        lat.U[pos] = makeTestMatrix(pos + 1);
        lat.U2[pos] = makeTestMatrix(pos + 101);
    }

    Group group;
    ForwardLightCone flc(&group);
    ForwardLightCone::LinkScratch scratch;
    flc.computeLinksTeam(&lat, N * N, scratch);

    for (int pos = 0; pos < N * N; ++pos) {
        Matrix expectedUx1Dagger = lat.U[lat.pospX[pos]];
        expectedUx1Dagger.conjg();
        Matrix expectedUx1 = lat.U[pos] * expectedUx1Dagger;
        CHECK(matricesClose(lat.Ux1[pos], expectedUx1, 1e-10));

        Matrix expectedUy1Dagger = lat.U[lat.pospY[pos]];
        expectedUy1Dagger.conjg();
        Matrix expectedUy1 = lat.U[pos] * expectedUy1Dagger;
        CHECK(matricesClose(lat.Uy1[pos], expectedUy1, 1e-10));

        Matrix expectedUx2Dagger = lat.U2[lat.pospX[pos]];
        expectedUx2Dagger.conjg();
        Matrix expectedUx2 = lat.U2[pos] * expectedUx2Dagger;
        CHECK(matricesClose(lat.Ux2[pos], expectedUx2, 1e-10));
    }
}

namespace {
// M = (U1 + U2) U^dagger, the argument of computeResidual()/
// computeJacobian()
Matrix residualProduct(const Matrix &U1, const Matrix &U2, const Matrix &U) {
    Matrix Udagger = U;
    Udagger.conjg();
    return (U1 + U2) * Udagger;
}

// The residual F_a of ForwardLightCone::computeResidual() for U1, U2 and U.
std::vector<double> residual(
    ForwardLightCone &flc, const Matrix &U1, const Matrix &U2,
    const Matrix &U) {
    std::vector<double> F(8);
    flc.computeResidual(
        residualProduct(U1, U2, U), flc.residualTraceCache(U1 + U2), F.data());
    return F;
}
}  // namespace

TEST_CASE(
    "ForwardLightCone::computeJacobian is the derivative of the residual for "
    "a rotation from the left (issue #107)") {
    Group group;
    ForwardLightCone flc(&group);
    std::mt19937 engine(9);
    for (int sample = 0; sample < 5; sample++) {
        const Matrix U1 = randomSU3(engine, 0.8);
        const Matrix U2 = randomSU3(engine, 0.8);
        const Matrix U = randomSU3(engine, 0.8);
        std::vector<double> J(64);
        flc.computeJacobian(residualProduct(U1, U2, U), J.data());
        // central difference of F under U -> e^{+-i eps t_b} U; findU()'s
        // convention is J = -dF/d delta
        const double eps = 1e-5;
        for (int b = 0; b < 8; b++) {
            std::vector<double> step(8, 0.);
            step[b] = eps;
            const Matrix Uplus = Matrix::fromAlgebraExponent(step) * U;
            step[b] = -eps;
            const Matrix Uminus = Matrix::fromAlgebraExponent(step) * U;
            const std::vector<double> Fplus = residual(flc, U1, U2, Uplus);
            const std::vector<double> Fminus = residual(flc, U1, U2, Uminus);
            for (int a = 0; a < 8; a++) {
                CAPTURE(a);
                CAPTURE(b);
                CHECK(
                    J[a * 8 + b]
                    == doctest::Approx(-(Fplus[a] - Fminus[a]) / (2. * eps))
                           .epsilon(1e-7)
                           .scale(1.));
            }
        }
    }
}

TEST_CASE(
    "ForwardLightCone::findU solves the matching condition for random "
    "links") {
    Group group;
    ForwardLightCone flc(&group);
    std::mt19937 engine(4);
    for (int sample = 0; sample < 20; sample++) {
        CAPTURE(sample);
        Matrix U1 = randomSU3(engine, 1.);
        Matrix U2 = randomSU3(engine, 1.);
        Matrix U;
        REQUIRE(flc.findU(U1, U2, U, /*retrySeed=*/sample));
        double Fzero = 0.;
        for (double F : residual(flc, U1, U2, U)) Fzero += std::abs(F);
        CHECK(Fzero < 1e-6);
        Matrix Udagger = U;
        Udagger.conjg();
        CHECK(matricesClose(U * Udagger, Matrix(1.), 1e-13));
        CHECK(std::abs(U.det() - 1.) < 1e-13);
    }
}

TEST_CASE(
    "ForwardLightCone::findU: U2 = 1 gives U = U1, also with a singular "
    "Jacobian, and the solution is kept") {
    // with U2 = 1 the matching condition is solved by U = U1, which is the
    // starting point U1 U2
    Group group;
    ForwardLightCone flc(&group);
    std::mt19937 engine(6);
    std::vector<Matrix> links;
    for (int sample = 0; sample < 5; sample++) {
        links.push_back(randomSU3(engine, 1.));
    }
    // (U1 + 1)(1 + U1^dagger) is singular here, and so is the Jacobian at
    // U = U1: the step after convergence fails and must not replace the
    // solution (it used to abort in GSL's error handler)
    Matrix singular(1.);
    singular.set(0, 0, -1.);
    singular.set(1, 1, -1.);
    links.push_back(singular);
    for (std::size_t i = 0; i < links.size(); i++) {
        CAPTURE(i);
        Matrix U1 = links[i];
        Matrix U2(1.);
        Matrix U;
        CHECK(flc.findU(U1, U2, U, /*retrySeed=*/i));
        CHECK(matricesClose(U, links[i], 1e-13));
    }
    // U1 = U2 = diag(-1, -1, 1): U1 + U2 is Hermitian, so U = U1 U2 = 1
    // solves the condition, and the Jacobian there is singular as well
    Matrix U1 = singular;
    Matrix U2 = singular;
    Matrix U;
    CHECK(flc.findU(U1, U2, U, /*retrySeed=*/7));
    double Fzero = 0.;
    for (double F : residual(flc, U1, U2, U)) Fzero += std::abs(F);
    CHECK(Fzero < 1e-6);
}

TEST_CASE(
    "ForwardLightCone::findU returns the identity and false for non-finite "
    "links") {
    Group group;
    ForwardLightCone flc(&group);
    Matrix U1(1.);
    // off the diagonal, where sanitizeU() used to miss it
    U1.set(
        0, 1,
        std::complex<double>(std::numeric_limits<double>::quiet_NaN(), 0.));
    Matrix U2(1.);
    Matrix U;
    CHECK_FALSE(flc.findU(U1, U2, U, /*retrySeed=*/0));
    CHECK(matricesClose(U, Matrix(1.), 1e-15));
}

TEST_CASE(
    "ForwardLightCone::solveAxb solves J x = F and reports a singular J") {
    Group group;
    ForwardLightCone flc(&group);
    // a regular J: diagonal 1..8 with one off-diagonal element
    std::vector<double> J(64, 0.), F(8), x(8, -1.);
    for (int a = 0; a < 8; a++) {
        J[a * 8 + a] = a + 1.;
        F[a] = 1.;
    }
    J[0 * 8 + 7] = 2.;
    CHECK(flc.solveAxb(J.data(), F.data(), x));
    CHECK(x[7] == doctest::Approx(1. / 8.));
    CHECK(x[0] == doctest::Approx(1. - 2. / 8.));
    for (int a = 1; a < 7; a++) CHECK(x[a] == doctest::Approx(1. / (a + 1.)));
    // a singular J (a zero row and column): false and x = 0 instead of the
    // GSL error handler
    std::fill(J.begin(), J.end(), 0.);
    for (int a = 0; a < 7; a++) J[a * 8 + a] = 1.;
    std::fill(x.begin(), x.end(), -1.);
    CHECK_FALSE(flc.solveAxb(J.data(), F.data(), x));
    for (double value : x) CHECK(value == 0.);
}

TEST_CASE(
    "ForwardLightCone::resetFieldsTeam zeroes U/U2/Uy2 and resets Ux1 to "
    "the identity") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    Lattice lat(&param, N);
    for (int pos = 0; pos < N * N; ++pos) {
        lat.U[pos] = makeTestMatrix(pos + 1);
        lat.U2[pos] = makeTestMatrix(pos + 2);
        lat.Uy2[pos] = makeTestMatrix(pos + 3);
        lat.Ux1[pos] = makeTestMatrix(pos + 4);
    }

    Group group;
    ForwardLightCone flc(&group);
    flc.resetFieldsTeam(&lat, N * N);

    const Matrix zero(0.0);
    const Matrix identity(1.0);
    for (int pos = 0; pos < N * N; ++pos) {
        CHECK(matricesClose(lat.U[pos], zero, 1e-14));
        CHECK(matricesClose(lat.U2[pos], zero, 1e-14));
        CHECK(matricesClose(lat.Uy2[pos], zero, 1e-14));
        CHECK(matricesClose(lat.Ux1[pos], identity, 1e-14));
    }
}

TEST_CASE(
    "ForwardLightCone::computePiTeam gives the traceless part of "
    "-(2i/g) A^eta (issue #102)") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    param.coupling.g = 1.3;
    Lattice lat(&param, N);
    // A^eta with a trace, as the lattice expression of
    // computeElectricFieldTeam() gives it
    for (int pos = 0; pos < N * N; ++pos) {
        lat.U[pos] = makeTestMatrix(pos + 1);
        REQUIRE(std::abs(lat.U[pos].trace()) > 1e-3);
    }

    Group group;
    ForwardLightCone flc(&group);
    flc.computePiTeam(&lat, &param, N * N);

    const std::complex<double> factor(0., -2. / param.coupling.g);
    for (int pos = 0; pos < N * N; ++pos) {
        CHECK(std::abs(lat.Ux2[pos].trace()) < 1e-14);
        Matrix expected = factor * lat.U[pos];
        const std::complex<double> traceThird = expected.trace() / 3.;
        for (int i = 0; i < 3; i++)
            expected.set(i, i, expected(i, i) - traceThird);
        CHECK(matricesClose(lat.Ux2[pos], expected, 1e-14));
    }
}

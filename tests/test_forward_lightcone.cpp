// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

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

TEST_CASE("ForwardLightCone::sanitizeU replaces a NaN U/U2 with identity") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    Lattice lat(&param, N);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    lat.U[5].set(0, 0, std::complex<double>(nan, 0.0));
    lat.U2[7].set(4, std::complex<double>(0.0, nan));

    Group group;
    ForwardLightCone flc(&group);
    flc.sanitizeU(&lat, N * N);

    const Matrix identity(1.0);
    CHECK(matricesClose(lat.U[5], identity, 1e-14));
    CHECK(matricesClose(lat.U2[7], identity, 1e-14));
    // An untouched cell must remain the identity it started as.
    CHECK(matricesClose(lat.U[0], identity, 1e-14));
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
// A random SU(3) matrix e^{i q_a t_a}, q_a Gaussian with width `width`.
Matrix randomSU3(std::mt19937 &engine, double width) {
    std::normal_distribution<double> gauss(0., width);
    std::vector<double> q(8);
    for (double &value : q) value = gauss(engine);
    return Matrix::fromAlgebraExponent(q);
}

// The residual F_a of ForwardLightCone::computeResidual() for U1, U2 and U.
std::vector<double> residual(
    ForwardLightCone &flc, Group &group, const Matrix &U1, const Matrix &U2,
    const Matrix &U) {
    const Matrix W = U1 + U2;
    Matrix Wdagger = W;
    Wdagger.conjg();
    std::vector<std::complex<double>> traceCache(8);
    for (int a = 0; a < 8; a++) {
        traceCache[a] = (group.getT(a) * (W - Wdagger)).trace();
    }
    Matrix Udagger = U;
    Udagger.conjg();
    std::vector<double> F(8);
    flc.computeResidual(W, Wdagger, U, Udagger, traceCache, F.data());
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
        Matrix Udagger = U;
        Udagger.conjg();
        std::vector<double> J(64);
        flc.computeJacobian(U1 + U2, Udagger, J.data());
        // central difference of F under U -> e^{+-i eps t_b} U; findU()'s
        // convention is J = -dF/d delta
        const double eps = 1e-5;
        for (int b = 0; b < 8; b++) {
            std::vector<double> step(8, 0.);
            step[b] = eps;
            const Matrix Uplus = Matrix::fromAlgebraExponent(step) * U;
            step[b] = -eps;
            const Matrix Uminus = Matrix::fromAlgebraExponent(step) * U;
            const std::vector<double> Fplus =
                residual(flc, group, U1, U2, Uplus);
            const std::vector<double> Fminus =
                residual(flc, group, U1, U2, Uminus);
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
        for (double F : residual(flc, group, U1, U2, U)) Fzero += std::abs(F);
        CHECK(Fzero < 1e-6);
        Matrix Udagger = U;
        Udagger.conjg();
        CHECK(matricesClose(U * Udagger, Matrix(1.), 1e-13));
        CHECK(std::abs(U.det() - 1.) < 1e-13);
    }
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

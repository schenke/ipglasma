#include <complex>
#include <limits>

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

TEST_CASE(
    "ForwardLightCone::computePlaquetteTeam wires the four links in the "
    "documented order") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    Lattice lat(&param, N);
    for (int pos = 0; pos < N * N; ++pos) {
        lat.Ux[pos] = makeTestMatrix(pos + 1);
        lat.Uy[pos] = makeTestMatrix(pos + 51);
    }

    Group group;
    ForwardLightCone flc(&group);
    ForwardLightCone::PlaquetteScratch scratch;
    flc.computePlaquetteTeam(&lat, N * N, scratch);

    for (int pos = 0; pos < N * N; ++pos) {
        Matrix UDx = lat.Ux[lat.pospY[pos]];
        UDx.conjg();
        Matrix UDy = lat.Uy[pos];
        UDy.conjg();
        Matrix expected = lat.Ux[pos] * (lat.Uy[lat.pospX[pos]] * (UDx * UDy));
        CHECK(matricesClose(lat.Uy1[pos], expected, 1e-10));
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

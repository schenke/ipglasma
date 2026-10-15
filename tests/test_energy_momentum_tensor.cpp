// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <cmath>
#include <random>
#include <vector>

#include "Cell.h"
#include "EnergyMomentumTensor.h"
#include "Group.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "SU3.h"
#include "doctest.h"
#include "test_helpers.h"

namespace {
const int N = 6;
const int it = 4;

void makeTmunuTestParam(Parameters &param) {
    makeLatticeParam(param, N);
    param.lattice.L = 3.;  // a = 0.5 fm
    param.coupling.g = 1.3;
    param.run.dtau = 0.1;
}

// Vacuum: unit links, vanishing electric fields, pi and phi.
void setVacuum(Lattice &lat) {
    const Matrix zero(0.);
    for (int pos = 0; pos < N * N; pos++) {
        lat.Ux[pos] = Matrix(1.);
        lat.Uy[pos] = Matrix(1.);
        lat.U[pos] = zero;
        lat.U2[pos] = zero;
        lat.Ux2[pos] = zero;
        lat.Uy2[pos] = zero;
    }
}

Matrix algebraElement(const Group &group, std::mt19937 &rng) {
    std::normal_distribution<double> gauss(0., 0.5);
    Matrix m(0.);
    for (int a = 0; a < 8; a++) m += gauss(rng) * group.getT(a);
    return m;
}

Matrix randomLink(std::mt19937 &rng) {
    std::normal_distribution<double> gauss(0., 0.4);
    std::vector<double> Q(8);
    for (double &q : Q) q = gauss(rng);
    return Matrix::fromAlgebraExponent(Q);
}

bool isGuard(int ix, int iy) {
    return ix == 0 || iy == 0 || ix == N - 1 || iy == N - 1;
}

void checkAllComponentsZero(const Cell &cell) {
    CHECK(cell.getEpsilon() == doctest::Approx(0.));
    CHECK(cell.getTtautau() == doctest::Approx(0.));
    CHECK(cell.getTxx() == doctest::Approx(0.));
    CHECK(cell.getTyy() == doctest::Approx(0.));
    CHECK(cell.getTetaeta() == doctest::Approx(0.));
    CHECK(cell.getTxy() == doctest::Approx(0.));
    CHECK(cell.getTtaux() == doctest::Approx(0.));
    CHECK(cell.getTtauy() == doctest::Approx(0.));
    CHECK(cell.getTtaueta() == doctest::Approx(0.));
    CHECK(cell.getTxeta() == doctest::Approx(0.));
    CHECK(cell.getTyeta() == doctest::Approx(0.));
}
}  // namespace

TEST_CASE("EnergyMomentumTensor::compute vanishes in the vacuum") {
    Parameters param;
    makeTmunuTestParam(param);
    Lattice lat(&param, N);
    setVacuum(lat);
    EnergyMomentumTensor::compute(&lat, &param, it);
    for (int pos = 0; pos < N * N; pos++) {
        CAPTURE(pos);
        checkAllComponentsZero(*lat.cells[pos]);
    }
}

TEST_CASE(
    "EnergyMomentumTensor::compute: a uniform E_1 field gives "
    "epsilon = -T^xx = T^yy = tau^2 T^etaeta and no off-diagonal "
    "components") {
    Parameters param;
    makeTmunuTestParam(param);
    Lattice lat(&param, N);
    setVacuum(lat);
    Group group;
    std::mt19937 rng(3);
    const Matrix E1 = algebraElement(group, rng);
    for (int pos = 0; pos < N * N; pos++) lat.U[pos] = E1;
    EnergyMomentumTensor::compute(&lat, &param, it);

    const double a = param.lattice.L / N;
    const double tauLattice = it * param.run.dtau;
    const double tau = tauLattice * a;  // fm
    // g^2/tau^2 (Tr E_1^2 at the cell and at its +y neighbor) / 2, in fm^-4
    const double expected = param.coupling.g * param.coupling.g
                            / (tauLattice * tauLattice)
                            * su3::traceSquare(E1).real() / std::pow(a, 4.);
    REQUIRE(expected > 0.);
    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            CAPTURE(ix);
            CAPTURE(iy);
            const Cell &cell = *lat.cells[lat.positionFromXY(ix, iy)];
            if (isGuard(ix, iy)) {
                checkAllComponentsZero(cell);
                continue;
            }
            CHECK(cell.getEpsilon() == doctest::Approx(expected));
            CHECK(cell.getTtautau() == doctest::Approx(expected));
            CHECK(cell.getTxx() == doctest::Approx(-expected));
            CHECK(cell.getTyy() == doctest::Approx(expected));
            CHECK(tau * tau * cell.getTetaeta() == doctest::Approx(expected));
            CHECK(cell.getTxy() == doctest::Approx(0.));
            CHECK(cell.getTtaux() == doctest::Approx(0.));
            CHECK(cell.getTtauy() == doctest::Approx(0.));
            CHECK(cell.getTtaueta() == doctest::Approx(0.));
            CHECK(cell.getTxeta() == doctest::Approx(0.));
            CHECK(cell.getTyeta() == doctest::Approx(0.));
        }
    }
}

TEST_CASE(
    "EnergyMomentumTensor::compute is traceless with a positive energy "
    "density for random fields") {
    Parameters param;
    makeTmunuTestParam(param);
    Lattice lat(&param, N);
    Group group;
    std::mt19937 rng(17);
    for (int pos = 0; pos < N * N; pos++) {
        lat.Ux[pos] = randomLink(rng);
        lat.Uy[pos] = randomLink(rng);
        lat.U[pos] = algebraElement(group, rng);
        lat.U2[pos] = algebraElement(group, rng);
        lat.Ux2[pos] = algebraElement(group, rng);
        lat.Uy2[pos] = algebraElement(group, rng);
    }
    EnergyMomentumTensor::compute(&lat, &param, it);

    const double tau = it * param.run.dtau * param.lattice.L / N;
    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            CAPTURE(ix);
            CAPTURE(iy);
            const Cell &cell = *lat.cells[lat.positionFromXY(ix, iy)];
            if (isGuard(ix, iy)) {
                checkAllComponentsZero(cell);
                continue;
            }
            CHECK(cell.getEpsilon() > 0.);
            CHECK(cell.getEpsilon() == cell.getTtautau());
            // conformal: T^mu_mu = T^tautau - T^xx - T^yy - tau^2 T^etaeta = 0
            CHECK(
                cell.getTxx() + cell.getTyy() + tau * tau * cell.getTetaeta()
                == doctest::Approx(cell.getTtautau()));
        }
    }
}

TEST_CASE(
    "EnergyMomentumTensor::compute does not use the links that leave the "
    "lattice, and T^{tau x}, T^{x eta} (T^{tau y}, T^{y eta}) vanish on the "
    "last interior column (row)") {
    // U_x at i = N - 1 and U_y at j = N - 1 point to sites outside the
    // lattice. T^{tau x} and T^{x eta} reach two sites ahead in x, so on
    // the last interior column i = N - 2 they would need them (issue #112);
    // likewise T^{tau y} and T^{y eta} on the row j = N - 2.
    Parameters param;
    makeTmunuTestParam(param);
    Lattice lat(&param, N);
    Group group;
    std::mt19937 rng(23);
    for (int pos = 0; pos < N * N; pos++) {
        lat.Ux[pos] = randomLink(rng);
        lat.Uy[pos] = randomLink(rng);
        lat.U[pos] = algebraElement(group, rng);
        lat.U2[pos] = algebraElement(group, rng);
        lat.Ux2[pos] = algebraElement(group, rng);
        lat.Uy2[pos] = algebraElement(group, rng);
    }
    auto components = [&](int pos) {
        const Cell &cell = *lat.cells[pos];
        return std::vector<double> {
            cell.getEpsilon(), cell.getTtautau(), cell.getTxx(),
            cell.getTyy(),     cell.getTetaeta(), cell.getTxy(),
            cell.getTtaux(),   cell.getTtauy(),   cell.getTtaueta(),
            cell.getTxeta(),   cell.getTyeta()};
    };
    EnergyMomentumTensor::compute(&lat, &param, it);
    std::vector<std::vector<double>> before(N * N);
    for (int pos = 0; pos < N * N; pos++) before[pos] = components(pos);

    for (int k = 0; k < N; k++) {
        lat.Ux[lat.positionFromXY(N - 1, k)] = randomLink(rng);
        lat.Uy[lat.positionFromXY(k, N - 1)] = randomLink(rng);
    }
    EnergyMomentumTensor::compute(&lat, &param, it);
    for (int ix = 1; ix < N - 1; ix++) {
        for (int iy = 1; iy < N - 1; iy++) {
            CAPTURE(ix);
            CAPTURE(iy);
            const int pos = lat.positionFromXY(ix, iy);
            CHECK(components(pos) == before[pos]);
            const Cell &cell = *lat.cells[pos];
            if (ix == N - 2) {
                CHECK(cell.getTtaux() == 0.);
                CHECK(cell.getTxeta() == 0.);
            } else {
                CHECK(cell.getTtaux() != 0.);
                CHECK(cell.getTxeta() != 0.);
            }
            if (iy == N - 2) {
                CHECK(cell.getTtauy() == 0.);
                CHECK(cell.getTyeta() == 0.);
            } else {
                CHECK(cell.getTtauy() != 0.);
                CHECK(cell.getTyeta() != 0.);
            }
            // their stencil reaches only one site ahead
            CHECK(cell.getTtaueta() != 0.);
            CHECK(cell.getTxy() != 0.);
        }
    }
}

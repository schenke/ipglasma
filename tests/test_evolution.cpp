// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "Cell.h"
#include "EnergyMomentumTensor.h"
#include "Evolution.h"
#include "ForwardLightCone.h"
#include "Group.h"
#include "Lattice.h"
#include "Parameters.h"
#include "doctest.h"
#include "test_helpers.h"

namespace {
// Matches test_lattice.cpp's makeLatticeParam: Parameters holds a
// PrettyOstream member, which holds a non-copyable/non-movable
// std::ostringstream, so this takes an out-parameter rather than
// returning by value.
void makeEvolutionTestParam(Parameters &param, int size) {
    param.lattice.size = size;
    param.lattice.L = static_cast<double>(size);  // a = 1 fm
    param.run.MPIRank = 0;
    param.event.eventId = 0;
    param.random.seed = 0;
    param.run.MPISize = 1;
    param.colorCharge.projectileX = 0.01;
    param.colorCharge.targetX = 0.01;
    param.coupling.g = 1.0;
    param.run.dtau = 0.1;
    param.coupling.runningCoupling =
        false;  // gfactor == 1 everywhere; see below
}

}  // namespace

TEST_CASE(
    "Evolution::finalFlowMeasurement fills epsilon and u^tau for "
    "Eccentricity::compute() also without hydro output") {
    // Every cell holds the T^{mu nu} of an ideal fluid with energy density
    // e, pressure p = e/3 and velocity v along x, so the flow solve must
    // return epsilon = e and u^tau = gamma.
    const int N = 4;
    const int it = 5;
    const double e = 1.0;
    const double p = e / 3.;
    const double v = 0.2;
    const double gamma = 1. / std::sqrt(1. - v * v);

    struct Case {
        bool writeHydro;
        bool computeEccentricities;
        bool expectSolve;
    };
    for (const Case &c : std::vector<Case> {
             {false, true, true}, {true, false, true}, {false, false, false}}) {
        CAPTURE(c.writeHydro);
        CAPTURE(c.computeEccentricities);
        Parameters param;
        makeEvolutionTestParam(param, N);
        param.output.writeHydro = c.writeHydro;
        param.output.writeJazma = false;
        param.output.writeTmunu = false;
        param.output.computeEccentricities = c.computeEccentricities;
        // an empty output grid: the hydro file has only its header
        param.output.sizeOutput = 0;
        param.output.LOutput = param.lattice.L;
        Lattice lat(&param, N);

        const double a = param.lattice.L / N;
        const double tau = it * param.run.dtau * a;
        for (int pos = 0; pos < N * N; ++pos) {
            Cell *cell = lat.cells[pos];
            cell->setTtautau((e + p) * gamma * gamma - p);
            cell->setTtaux((e + p) * gamma * gamma * v);
            cell->setTxx((e + p) * gamma * gamma * v * v + p);
            cell->setTyy(p);
            cell->setTetaeta(p / (tau * tau));
        }

        int nn[2] = {N, N};
        Evolution evolution(nn);
        evolution.finalFlowMeasurement(&lat, &param, it);

        for (int pos = 0; pos < N * N; ++pos) {
            CAPTURE(pos);
            if (c.expectSolve) {
                CHECK(lat.cells[pos]->getEpsilon() == doctest::Approx(e));
                CHECK(lat.cells[pos]->getutau() == doctest::Approx(gamma));
            } else {
                // the expensive solve is skipped when nothing needs it
                CHECK(lat.cells[pos]->getutau() == 0.);
            }
        }
        std::remove("epsilon-u-Hydro-TauHydro-0.dat");
    }
}

TEST_CASE(
    "Evolution::outputSteps rounds the output times down to distinct steps "
    "before the final one") {
    const double stepLength = 0.01;  // fm/c
    const int itmax = 40;
    CHECK(Evolution::outputSteps({}, stepLength, itmax).empty());
    CHECK(
        Evolution::outputSteps({0.1, 0.2, 0.3}, stepLength, itmax)
        == std::vector<int> {10, 20, 30});
    // unsorted, duplicates after rounding, 0.005 is step 0 and 0.4 and
    // 0.5 are not before the final step
    CHECK(
        Evolution::outputSteps(
            {0.3, 0.104, 0.1, 0.005, 0.4, 0.5}, stepLength, itmax)
        == std::vector<int> {10, 30});
}

TEST_CASE(
    "Evolution::evolveU rotates the links by exp(i omega E) and keeps them in "
    "SU(3)") {
    const int N = 4;
    Parameters param;
    makeEvolutionTestParam(param, N);
    param.coupling.g = 1.3;
    Lattice lat(&param, N);
    // E large enough that the old 1 + X + X^2/2 left SU(3) visibly
    Matrix E(0.);
    E.set(0, 1, complex<double>(0.4, 0.3));
    E.set(1, 0, complex<double>(0.4, -0.3));
    E.set(1, 2, complex<double>(-0.2, 0.5));
    E.set(2, 1, complex<double>(-0.2, -0.5));
    E.set(0, 0, 0.6);
    E.set(1, 1, -0.1);
    E.set(2, 2, -0.5);
    for (int pos = 0; pos < N * N; ++pos) {
        lat.U[pos] = E;
        lat.U2[pos] = -1. * E;
        lat.Ux[pos] = Matrix(1.);
        lat.Uy[pos] = Matrix(1.);
    }
    const double dtau = 0.1, tau = 0.3;
    int nn[2] = {N, N};
    Evolution evolution(nn);
    for (int step = 0; step < 50; step++) {
        evolution.evolveU(&lat, &param, dtau, tau);
    }
    // 50 rotations about the same generator add up
    const double omega =
        param.coupling.g * param.coupling.g * dtau / (tau + dtau / 2.);
    const Matrix expectedX = Matrix::expiHermitian((50. * omega) * E);
    const Matrix expectedY = Matrix::expiHermitian((-50. * omega) * E);
    for (int pos = 0; pos < N * N; ++pos) {
        for (int i = 0; i < 9; i++) {
            CHECK(std::abs(lat.Ux[pos](i) - expectedX(i)) < 1e-12);
            CHECK(std::abs(lat.Uy[pos](i) - expectedY(i)) < 1e-12);
        }
        const Matrix unitarity = Matrix::prodAconjB(lat.Ux[pos], lat.Ux[pos]);
        for (int i = 0; i < 9; i++) {
            CHECK(std::abs(unitarity(i) - ((i % 4 == 0) ? 1. : 0.)) < 1e-13);
        }
        CHECK(std::abs(lat.Ux[pos].det() - 1.) < 1e-13);
    }
}

TEST_CASE(
    "With periodic boundaries, the forward light cone, the evolution and "
    "T^{mu nu} commute with a cyclic shift of the Wilson lines") {
    // A constant color-charge density (useNucleus 0) is an infinite target,
    // whose lattice is periodic: no site is special, so shifting the Wilson
    // lines of both nuclei by s shifts every field and T^{mu nu} by s, also
    // where the stencils cross an edge. With clamped boundaries the edges
    // are special, and the shift changes the result.
    const int N = 6;
    const int steps = 3;
    const int sx = 1;
    const int sy = 4;
    Parameters param;
    makeEvolutionTestParam(param, N);
    param.coupling.g = 1.3;
    Group group;
    std::mt19937 rng(41);
    std::vector<Matrix> VA, VB;
    for (int pos = 0; pos < N * N; ++pos) {
        VA.push_back(randomSU3(rng, 0.5));
        VB.push_back(randomSU3(rng, 0.5));
    }

    // the largest difference between the fields and T^{mu nu} of the
    // shifted Wilson lines at x + s and those of the original ones at x
    auto shiftDeviation = [&](LatticeBoundary boundary) {
        Lattice lat(&param, N, boundary);
        Lattice shifted(&param, N, boundary);
        auto target = [&](int pos) {
            return lat.positionFromXY(
                (lat.xFromPosition(pos) + sx) % N,
                (lat.yFromPosition(pos) + sy) % N);
        };
        for (int pos = 0; pos < N * N; ++pos) {
            lat.U[pos] = VA[pos];
            lat.U2[pos] = VB[pos];
            shifted.U[target(pos)] = VA[pos];
            shifted.U2[target(pos)] = VB[pos];
        }
        int nn[2] = {N, N};
        Evolution evolution(nn);
        const double dtau = param.run.dtau;
        for (Lattice *l : {&lat, &shifted}) {
            ForwardLightCone(&group).initialize(l, &param);
            evolution.evolvePhi(l, &param, dtau, 0.);
            evolution.evolveU(l, &param, dtau, 0.);
            for (int it = 1; it < steps; ++it) {
                evolution.evolvePi(l, &param, dtau, it * dtau);
                evolution.evolveE(l, &param, dtau, it * dtau);
                evolution.evolvePhi(l, &param, dtau, it * dtau);
                evolution.evolveU(l, &param, dtau, it * dtau);
            }
            EnergyMomentumTensor::compute(l, &param, steps);
        }

        double deviation = 0.;
        for (int pos = 0; pos < N * N; ++pos) {
            const int to = target(pos);
            for (auto field :
                 {&Lattice::Ux, &Lattice::Uy, &Lattice::U, &Lattice::U2,
                  &Lattice::Ux2, &Lattice::Uy2}) {
                for (int k = 0; k < 9; ++k) {
                    deviation = std::max(
                        deviation, std::abs(
                                       (shifted.*field)[to].get(k)
                                       - (lat.*field)[pos].get(k)));
                }
            }
            const Cell &cell = *lat.cells[pos];
            const Cell &moved = *shifted.cells[to];
            for (auto component :
                 {&Cell::getEpsilon, &Cell::getTxx, &Cell::getTyy,
                  &Cell::getTxy, &Cell::getTtaux, &Cell::getTtauy,
                  &Cell::getTtaueta, &Cell::getTxeta, &Cell::getTyeta}) {
                const double value = (cell.*component)();
                deviation = std::max(
                    deviation, std::abs((moved.*component)() - value)
                                   / (1. + std::abs(value)));
            }
        }
        return deviation;
    };
    CHECK(shiftDeviation(LatticeBoundary::Periodic) < 1e-10);
    CHECK(shiftDeviation(LatticeBoundary::Clamped) > 1e-3);
}

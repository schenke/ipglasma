// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "CollisionGeometry.h"
#include "Glauber.h"
#include "Init.h"
#include "InputFile.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "PhysConst.h"
#include "Random.h"
#include "doctest.h"
#include "test_helpers.h"

using PhysConst::hbarc;

TEST_CASE(
    "Init::computeEffectiveRapidity passes rapidity through unchanged "
    "when usePseudoRapidity is off") {
    Parameters param;
    makeInitTestParam(param, 4);
    param.colorCharge.usePseudoRapidity = false;
    param.colorCharge.rapidity = 1.5;

    int nn[2] = {4, 4};
    Init init(nn);
    CHECK(init.computeEffectiveRapidity(&param) == doctest::Approx(1.5));
}

TEST_CASE(
    "Init::computeEffectiveRapidity converts pseudorapidity to rapidity "
    "when enabled") {
    Parameters param;
    makeInitTestParam(param, 4);
    param.colorCharge.usePseudoRapidity = true;
    param.colorCharge.rapidity = 1.0;
    param.colorCharge.jacobianMass = 0.14;
    param.collision.sqrtS = 200.0;

    int nn[2] = {4, 4};
    Init init(nn);
    auto rapidityOf = [&](double eta) {
        param.colorCharge.rapidity = eta;
        return init.computeEffectiveRapidity(&param);
    };

    // independent reference: y = atanh(p_z/E) for p_T = P and mass m
    const double m = 0.14;
    const double P = 0.13 + 0.32 * std::pow(200.0 / 1000., 0.115);
    for (double eta : {-2.5, -1.0, 0.3, 1.0, 4.0}) {
        CAPTURE(eta);
        const double pz = P * std::sinh(eta);
        const double E =
            std::sqrt(P * P * std::cosh(eta) * std::cosh(eta) + m * m);
        CHECK(rapidityOf(eta) == doctest::Approx(std::atanh(pz / E)));
        // a massive particle has |y| < |eta|, with the sign of eta
        CHECK(std::abs(rapidityOf(eta)) < std::abs(eta));
        CHECK(rapidityOf(-eta) == doctest::Approx(-rapidityOf(eta)));
    }
    // eta = 0 is y = 0 (issue #66: a misplaced parenthesis gave y != 0)
    CHECK(rapidityOf(0.) == doctest::Approx(0.).epsilon(1e-12));

    // massless particles have y = eta
    param.colorCharge.jacobianMass = 0.;
    CHECK(rapidityOf(1.7) == doctest::Approx(1.7));
}

TEST_CASE(
    "Init::computeCellColorCharge: Q_s^2 at projectileX/targetX with a fixed "
    "x, and at jimwlkInitialX with JIMWLK") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);  // a = 1 fm
    param.colorCharge.QsMuRatio = 1.;
    param.coupling.g = 1.;
    param.colorCharge.projectileX = 1e-3;  // y = ln 10
    param.colorCharge.targetX = 1e-4;      // y = ln 100
    param.jimwlk.initialX = 5e-3;          // y = ln 2
    Lattice lat(&param, N);
    const int pos = 5;
    lat.cells[pos]->setTpA(1.);
    lat.cells[pos]->setTpB(2.);

    int nn[2] = {N, N};
    Init init(nn);
    const std::string fileName = writeLinearQsTable();  // Qs^2 = 2 T + 3 y
    init.readQsTable(fileName);
    std::remove(fileName.c_str());
    // g^2 mu^2 = Qs^2 a^2 / hbarc^2 with QsMuRatio = g = 1
    const double toLattice = 1. / (hbarc * hbarc);

    param.jimwlk.enabled = false;
    param.colorCharge.useFluctuatingX = false;
    init.computeCellColorCharge(&lat, &param, pos, 1., /*rapidity=*/5.);
    CHECK(
        lat.cells[pos]->getg2mu2A()
        == doctest::Approx((2. * 1. + 3. * std::log(10.)) * toLattice));
    CHECK(
        lat.cells[pos]->getg2mu2B()
        == doctest::Approx((2. * 2. + 3. * std::log(100.)) * toLattice));

    param.jimwlk.enabled = true;
    init.computeCellColorCharge(&lat, &param, pos, 1., /*rapidity=*/5.);
    CHECK(
        lat.cells[pos]->getg2mu2A()
        == doctest::Approx((2. * 1. + 3. * std::log(2.)) * toLattice));
    CHECK(
        lat.cells[pos]->getg2mu2B()
        == doctest::Approx((2. * 2. + 3. * std::log(2.)) * toLattice));
}

TEST_CASE(
    "Init::computeWilsonLineMomentumKernel: massless kernel is 1/kt^2, zero "
    "at kt=0") {
    const int N = 4;
    int nn[2] = {N, N};
    Init init(nn);

    std::vector<double> kernel = init.computeWilsonLineMomentumKernel(
        N, N * N, /*m=*/0.0, /*UVdamp=*/0.0);
    REQUIRE(kernel.size() == static_cast<std::size_t>(N * N));

    // pos=0 -> (i=0,j=0): kx=ky=-pi, kt2 = 4*(sin(-pi/2)^2*2) = 8.
    CHECK(kernel[0] == doctest::Approx(1.0 / 8.0));
    // pos = (N/2)*N + N/2 -> (i=j=N/2): kx=ky=0, kt2=0 -> defined as 0.
    CHECK(kernel[latticeIndex(N / 2, N / 2, N)] == doctest::Approx(0.0));
}

TEST_CASE(
    "Init::computeWilsonLineMomentumKernel: massive kernel matches "
    "1/(kt2+m^2)*exp(-sqrt(kt2)*UVdamp)") {
    const int N = 4;
    int nn[2] = {N, N};
    Init init(nn);
    const double m = 0.5;

    std::vector<double> kernel =
        init.computeWilsonLineMomentumKernel(N, N * N, m, /*UVdamp=*/0.0);
    // At kt2=0 (i=j=N/2), the kernel is exactly 1/m^2.
    CHECK(
        kernel[latticeIndex(N / 2, N / 2, N)]
        == doctest::Approx(1.0 / (m * m)));
}

TEST_CASE(
    "Init::computeWilsonLineMomentumKernel: the UV damping "
    "exp(-|kt| UVdamp) applies with and without the mass, and K(0) = 0 "
    "stays for m = 0") {
    const int N = 8;
    int nn[2] = {N, N};
    Init init(nn);
    const double UVdamp = 0.7;
    for (const double m : {0.0, 0.5}) {
        CAPTURE(m);
        const std::vector<double> kernel =
            init.computeWilsonLineMomentumKernel(N, N * N, m, UVdamp);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                const double kx =
                    2. * M_PI * (-0.5 + static_cast<double>(i) / N);
                const double ky =
                    2. * M_PI * (-0.5 + static_cast<double>(j) / N);
                const double kt2 = 4.
                                   * (std::sin(kx / 2.) * std::sin(kx / 2.)
                                      + std::sin(ky / 2.) * std::sin(ky / 2.));
                const double value = kernel[latticeIndex(i, j, N)];
                if (m == 0. && kt2 == 0.) {
                    CHECK(value == 0.);
                } else {
                    CHECK(
                        value
                        == doctest::Approx(
                            std::exp(-std::sqrt(kt2) * UVdamp)
                            / (kt2 + m * m)));
                }
            }
        }
    }
}

TEST_CASE(
    "Init::computeWilsonLineColorChargeScales: g*sqrt(g2mu2*invNy) per "
    "nucleus") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    Lattice lat(&param, N);
    for (int pos = 0; pos < N * N; ++pos) {
        lat.cells[pos]->setg2mu2A(4.0);
        lat.cells[pos]->setg2mu2B(9.0);
    }

    int nn[2] = {N, N};
    Init init(nn);
    const ColorChargeScales scales = init.computeWilsonLineColorChargeScales(
        &lat, N * N, /*g=*/2.0, /*invNy=*/0.25);
    const std::vector<double> &scaleA = scales.projectile;
    const std::vector<double> &scaleB = scales.target;

    REQUIRE(scaleA.size() == static_cast<std::size_t>(N * N));
    REQUIRE(scaleB.size() == static_cast<std::size_t>(N * N));
    for (int pos = 0; pos < N * N; ++pos) {
        CHECK(scaleA[pos] == doctest::Approx(2.0 * std::sqrt(4.0 * 0.25)));
        CHECK(scaleB[pos] == doctest::Approx(2.0 * std::sqrt(9.0 * 0.25)));
    }
}

TEST_CASE(
    "Init::setConstantColorChargeDensity: flat background sets g2mu2A=g2mu2B="
    "(g2mu/g)^2 everywhere, with g2mu = g2muGeV in lattice units") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    // a = 1 fm, so this is g2mu = 6 in lattice units
    param.collision.g2muGeV = 6.0 * PhysConst::hbarc;
    param.coupling.g = 2.0;
    Lattice lat(&param, N);

    int nn[2] = {N, N};
    Init init(nn);
    init.setConstantColorChargeDensity(&lat, &param);

    for (int pos = 0; pos < N * N; ++pos) {
        CHECK(lat.cells[pos]->getg2mu2A() == doctest::Approx(9.0));
        CHECK(lat.cells[pos]->getg2mu2B() == doctest::Approx(9.0));
    }
    CHECK(param.event.success == 1);
}

TEST_CASE(
    "Init::shiftFieldsWithImpactParameter leaves the fields unchanged at "
    "b = 0 and moves them by the lattice shifts otherwise") {
    // L = 16 fm, N = 48: truncating (x + L/2)/a, with x = -L/2 + a ix,
    // used to copy 10 rows from their neighbours at b = 0
    const int N = 48;
    Parameters param;
    makeInitTestParam(param, N);
    param.lattice.L = 16.;
    const double a = param.lattice.L / N;
    int nn[2] = {N, N};
    Init init(nn);
    Lattice lat(&param, N);
    // a different multiple of the identity at every site, none of them 1
    auto fill = [&]() {
        for (int pos = 0; pos < N * N; pos++) {
            lat.U[pos] = Matrix(2. + pos);
            lat.U2[pos] = Matrix(-2. - pos);
        }
    };

    fill();
    param.event.b = 0.;
    param.event.phiRP = 1.;
    init.shiftFieldsWithImpactParameter(&lat, &param);
    int changed = 0;
    for (int pos = 0; pos < N * N; pos++) {
        if (lat.U[pos].get(0).real() != 2. + pos) changed++;
        if (lat.U2[pos].get(0).real() != -2. - pos) changed++;
    }
    CHECK(changed == 0);

    fill();
    param.event.b = 1.3;
    param.event.phiRP = 2.2;
    const NucleusShifts shifts = CollisionGeometry::latticeShifts(1.3, 2.2, a);
    REQUIRE(shifts.xA - shifts.xB != 0);
    REQUIRE(shifts.yA - shifts.yB != 0);
    init.shiftFieldsWithImpactParameter(&lat, &param);
    // the value of site (ix - sx, iy - sy) before the shift, the identity
    // where the field is shifted in from outside the lattice
    auto before = [&](int ix, int iy, int sx, int sy, double sign) {
        const int jx = ix - sx, jy = iy - sy;
        if (jx < 0 || jx >= N || jy < 0 || jy >= N) return 1.;
        return sign * (2. + lat.positionFromXY(jx, jy));
    };
    int wrong = 0;
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            const int pos = lat.positionFromXY(ix, iy);
            if (lat.U[pos].get(0).real()
                != before(ix, iy, shifts.xA, shifts.yA, 1.))
                wrong++;
            if (lat.U2[pos].get(0).real()
                != before(ix, iy, shifts.xB, shifts.yB, -1.))
                wrong++;
        }
    }
    CHECK(wrong == 0);
}

TEST_CASE(
    "Init::computeSmoothNucleusThickness builds centered projectile (A) and "
    "target (B) profiles, and CollisionGeometry::scanOverlap finds their "
    "shifted overlap") {
    // Note: Glauber::interNuPInSP()/interNuTInST() cache their tables in
    // function-local statics on first use, so they hold Cu/Au for every
    // later call in this test binary.
    const int N = 64;
    const double L = 30.;
    const double a = L / N;
    Parameters param;
    makeInitTestParam(param, N);
    param.lattice.L = L;
    param.collision.useNucleus = true;
    param.nucleus.useSmoothNucleus = true;
    param.coupling.g = 1.;
    param.colorCharge.QsMuRatio = 0.643;
    param.event.b = 5.;  // a stale b must not shift the profiles

    Glauber glauber;
    glauber.initGlauber(
        42., /*target=*/"Au", /*projectile=*/"Cu", 0., false, 0., 0., 0., 0.,
        0., 0., false, 0., 0., 0., 1000);
    REQUIRE(glauber.nucleusA1() == 63);   // projectile
    REQUIRE(glauber.nucleusA2() == 197);  // target

    int nn[2] = {N, N};
    Init init(nn);
    Lattice lat(&param, N);
    init.computeSmoothNucleusThickness(&lat, &param, &glauber);

    // each profile integrates to its own mass number (up to the cut tails)
    // and is centered at the origin
    double sumA = 0., sumB = 0., xA = 0., xB = 0.;
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            const int pos = lat.positionFromXY(ix, iy);
            const double x = -L / 2. + a * ix;
            const double TA = lat.cells[pos]->getTpA() * a * a / hbarc / hbarc;
            const double TB = lat.cells[pos]->getTpB() * a * a / hbarc / hbarc;
            sumA += TA;
            sumB += TB;
            xA += x * TA;
            xB += x * TB;
            // g2mu2 is only used for the Qs averages; any positive
            // profile will do for the overlap test below
            lat.cells[pos]->setg2mu2A(lat.cells[pos]->getTpA());
            lat.cells[pos]->setg2mu2B(lat.cells[pos]->getTpB());
        }
    }
    CHECK(sumA == doctest::Approx(63.).epsilon(0.02));
    CHECK(sumB == doctest::Approx(197.).epsilon(0.02));
    CHECK(std::abs(xA / sumA) < a);
    CHECK(std::abs(xB / sumB) < a);

    // a smooth nucleus has no nucleons
    std::vector<ReturnValue> noNucleonsA, noNucleonsB;
    CollisionGeometry geometry(noNucleonsA, noNucleonsB);
    auto overlapCells = [&](double b) {
        return geometry.scanOverlap(&lat, &param, N, a, b, /*phiRP=*/0.).count;
    };
    const int central = overlapCells(0.);
    const int shifted = overlapCells(4.);
    CHECK(central > 0);
    CHECK(shifted > 0);
    CHECK(shifted < central);
    CHECK(overlapCells(40.) == 0);
}

namespace {
// Reads the shipped example input into param, with the values of the keys
// in changes replaced.
void readExampleInput(
    Parameters &param, const std::map<std::string, std::string> &changes) {
    std::ifstream file(std::string(IPGLASMA_SOURCE_DIR) + "/input");
    REQUIRE(file.good());
    std::string line, text;
    while (std::getline(file, line)) {
        std::istringstream tokens(line);
        std::string key;
        tokens >> key;
        const auto change = changes.find(key);
        text += (change == changes.end()) ? line : key + " " + change->second;
        text += "\n";
    }
    std::istringstream in(text);
    REQUIRE(param.readInput(InputFile(in, "test")).empty());
}

// The positions and proton flags (x y proton) of the nucleons in
// NpartList0.dat, without the collided flag, which Gaussian wounding draws
// at random.
std::vector<std::string> nucleonLines() {
    std::ifstream in("NpartList0.dat");
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line);) {
        std::istringstream columns(line);
        std::string x, y, proton;
        if (columns >> x >> y >> proton) {
            lines.push_back(x + " " + y + " " + proton);
        }
    }
    return lines;
}

// A small lattice for colliding light nuclei with Init::sampleCollision().
struct LightCollision {
    Parameters param;
    Glauber glauber;
    Random random;

    LightCollision(
        const std::string &nucleus, double b, int fixedNpart,
        bool writeNpartList) {
        readExampleInput(
            param, {{"projectile", nucleus},
                    {"target", nucleus},
                    {"size", "32"},
                    {"L", "8"},
                    {"bMin", std::to_string(b)},
                    {"bMax", std::to_string(b)},
                    {"rotateReactionPlane", "0"},
                    {"useFixedNpart", std::to_string(fixedNpart)},
                    {"writeNpartList", writeNpartList ? "1" : "0"},
                    {"writeNcollList", "0"},
                    {"writeNgluonEstimators", "0"}});
        param.run.MPIRank = 0;
        param.run.MPISize = 1;
        param.event.eventId = 0;
        glauber.initGlauber(
            param.collision.sigmaNN, nucleus, nucleus, 0., false, 0., 0., 0.,
            0., 0., 0., false, 0., 0., 0., 100);
        random.init_genrand64(11);
        random.gslRandomInit(11);
    }
};

void removeCollisionFiles() {
    std::remove("NpartList0.dat");
    std::remove("usedParameters0.dat");
}
}  // namespace

TEST_CASE(
    "Init::sampleCollision accepts a central p+p collision on the first "
    "attempt") {
    LightCollision collision("p", /*b=*/0., /*fixedNpart=*/0, false);
    Parameters &param = collision.param;
    const int N = param.lattice.size;
    int nn[2] = {N, N};
    Init init(nn);
    init.readQsTable(
        std::string(IPGLASMA_SOURCE_DIR) + "/"
        + param.colorCharge.nucleusQsTableFileName);
    Lattice lat(&param, N);
    CHECK(init.sampleCollision(
        &lat, &param, &collision.random, &collision.glauber,
        /*resampleNuclei=*/true, /*maxConfigurations=*/1,
        /*maxLatticeRejections=*/1));
    CHECK(param.event.success == 1);
    CHECK(param.event.Npart == 2);
    CHECK(param.event.b == 0.);
    removeCollisionFiles();
}

TEST_CASE(
    "Init::sampleCollision stops after the given numbers of attempts and of "
    "lattice rejections") {
    SUBCASE("p+p without overlap: every attempt is rejected on the lattice") {
        // the protons are 7 fm apart, beyond twice the wounding distance
        LightCollision collision("p", /*b=*/7., /*fixedNpart=*/0, false);
        Parameters &param = collision.param;
        const int N = param.lattice.size;
        int nn[2] = {N, N};
        Init init(nn);
        init.readQsTable(
            std::string(IPGLASMA_SOURCE_DIR) + "/"
            + param.colorCharge.nucleusQsTableFileName);
        Lattice lat(&param, N);
        CHECK_FALSE(init.sampleCollision(
            &lat, &param, &collision.random, &collision.glauber, true,
            /*maxConfigurations=*/100, /*maxLatticeRejections=*/3));
        CHECK(param.event.success == 0);
    }
    SUBCASE("an Npart that p+p cannot reach is rejected before the lattice") {
        LightCollision collision("p", /*b=*/0., /*fixedNpart=*/3, false);
        Parameters &param = collision.param;
        const int N = param.lattice.size;
        int nn[2] = {N, N};
        Init init(nn);
        init.readQsTable(
            std::string(IPGLASMA_SOURCE_DIR) + "/"
            + param.colorCharge.nucleusQsTableFileName);
        Lattice lat(&param, N);
        CHECK_FALSE(init.sampleCollision(
            &lat, &param, &collision.random, &collision.glauber, true,
            /*maxConfigurations=*/20, /*maxLatticeRejections=*/1));
        CHECK(param.event.success == 0);
    }
    removeCollisionFiles();
}

TEST_CASE(
    "Init::sampleCollision draws new nuclei for every attempt, or keeps "
    "them for read Wilson lines") {
    // d+d cannot have Npart = 5, so every attempt is rejected after writing
    // the nucleon positions (at b = 0, unshifted) to NpartList0.dat
    LightCollision collision("d", /*b=*/0., /*fixedNpart=*/5, true);
    Parameters &param = collision.param;
    const int N = param.lattice.size;
    int nn[2] = {N, N};
    Init init(nn);
    Lattice lat(&param, N);
    auto attempt = [&](bool resampleNuclei) {
        CHECK_FALSE(init.sampleCollision(
            &lat, &param, &collision.random, &collision.glauber, resampleNuclei,
            /*maxConfigurations=*/1,
            /*maxLatticeRejections=*/1));
        return nucleonLines();
    };
    const std::vector<std::string> first = attempt(true);
    REQUIRE(first.size() == 4);
    const std::vector<std::string> redrawn = attempt(true);
    CHECK(redrawn.size() == 4);
    CHECK(redrawn != first);
    const std::vector<std::string> kept = attempt(false);
    CHECK(kept == redrawn);
    removeCollisionFiles();
}

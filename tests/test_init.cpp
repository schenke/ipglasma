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
    // with removeZeroMode, only K(0) changes, to 0
    std::vector<double> removed = init.computeWilsonLineMomentumKernel(
        N, N * N, m, /*UVdamp=*/0.0, /*removeZeroMode=*/true);
    for (int pos = 0; pos < N * N; pos++) {
        CAPTURE(pos);
        CHECK(
            removed[pos]
            == (pos == latticeIndex(N / 2, N / 2, N) ? 0. : kernel[pos]));
    }
}

TEST_CASE(
    "Init::computeWilsonLineMomentumKernel: the UV damping "
    "exp(-|kt| UVdamp) applies with and without the mass, and K(0) = 0 "
    "stays for m = 0 and with removeZeroMode") {
    const int N = 8;
    int nn[2] = {N, N};
    Init init(nn);
    const double UVdamp = 0.7;
    for (const bool removeZeroMode : {false, true}) {
        for (const double m : {0.0, 0.5}) {
            CAPTURE(m);
            CAPTURE(removeZeroMode);
            const std::vector<double> kernel =
                init.computeWilsonLineMomentumKernel(
                    N, N * N, m, UVdamp, removeZeroMode);
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    const double kx =
                        2. * M_PI * (-0.5 + static_cast<double>(i) / N);
                    const double ky =
                        2. * M_PI * (-0.5 + static_cast<double>(j) / N);
                    const double kt2 =
                        4.
                        * (std::sin(kx / 2.) * std::sin(kx / 2.)
                           + std::sin(ky / 2.) * std::sin(ky / 2.));
                    const double value = kernel[latticeIndex(i, j, N)];
                    if ((m == 0. || removeZeroMode) && kt2 == 0.) {
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
    "Init::computeSmoothNucleusThickness gives a smooth deuteron the size of "
    "its nucleon density, half of the Hulthen neutron-proton separation") {
    const int N = 128;
    const double L = 20.;
    const double a = L / N;
    Parameters param;
    makeInitTestParam(param, N);
    param.lattice.L = L;
    param.collision.useNucleus = true;
    param.nucleus.useSmoothNucleus = true;

    Glauber glauber;
    glauber.initGlauber(
        42., /*target=*/"Au", /*projectile=*/"d", 0., false, 0., 0., 0., 0., 0.,
        0., false, 0., 0., 0., 1000);
    REQUIRE(glauber.nucleusA1() == 2);
    int nn[2] = {N, N};
    Init init(nn);
    Lattice lat(&param, N);
    init.computeSmoothNucleusThickness(&lat, &param, &glauber);

    double sum = 0., sumR2 = 0.;
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            const double x = -L / 2. + a * ix;
            const double y = -L / 2. + a * iy;
            const double T = lat.cells[lat.positionFromXY(ix, iy)]->getTpA() * a
                             * a / hbarc / hbarc;
            sum += T;
            sumR2 += (x * x + y * y) * T;
        }
    }
    CHECK(sum == doctest::Approx(2.).epsilon(0.02));
    const double rms = std::sqrt(sumR2 / sum);

    // The Hulthen thickness itself is that of the neutron-proton separation
    // r. Evaluated with the same cut of the tails (T < 0.001) on a lattice
    // with twice the spacing and size, it must be exactly twice as wide as
    // the nucleons, which sit at r/2 from the center.
    double sepSum = 0., sepSumR2 = 0.;
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            const double x = 2. * (-L / 2. + a * ix);
            const double y = 2. * (-L / 2. + a * iy);
            double T = glauber.interNuPInSP(std::sqrt(x * x + y * y));
            if (T < 0.001) T = 0.;
            sepSum += T;
            sepSumR2 += (x * x + y * y) * T;
        }
    }
    CHECK(rms == doctest::Approx(0.5 * std::sqrt(sepSumR2 / sepSum)));

    // the shape is the Hulthen one: along x, T(s)/T(0) is the projected
    // separation density at 2s, (int dz rho(sqrt(4s^2 + z^2)))/(same at 0)
    // with rho(r) = (e^{-ar} - e^{-br})^2/r^2 (a Woods-Saxon such as Cu
    // would give 0.86 instead of 0.21 at s = 1 fm)
    auto projectedHulthen = [](double u) {
        const double ha = 0.228, hb = 1.18;
        auto rho = [&](double r) {
            if (r < 1e-9) return (hb - ha) * (hb - ha);
            const double f = (std::exp(-ha * r) - std::exp(-hb * r)) / r;
            return f * f;
        };
        // Simpson's rule on z in [0, 60] fm
        const int n = 6000;
        const double h = 60. / n;
        double integral = rho(u) + rho(std::sqrt(u * u + 3600.));
        for (int i = 1; i < n; i++) {
            const double z = i * h;
            integral += (i % 2 ? 4. : 2.) * rho(std::sqrt(u * u + z * z));
        }
        return integral * h / 3.;
    };
    const double centralT =
        lat.cells[lat.positionFromXY(N / 2, N / 2)]->getTpA();
    for (const int k : {3, 6, 13}) {
        const double sk = a * k;  // 0.47, 0.94 and 2.03 fm
        CAPTURE(sk);
        CHECK(
            lat.cells[lat.positionFromXY(N / 2 + k, N / 2)]->getTpA() / centralT
            == doctest::Approx(projectedHulthen(2. * sk) / projectedHulthen(0.))
                   .epsilon(0.02));
    }

    // Analytically, the separation density (e^{-ar} - e^{-br})^2/r^2,
    // a = 0.228/fm, b = 1.18/fm, has <r^2> = (4.07 fm)^2, so a nucleon has
    // the transverse <s^2> = 2/3 <r^2>/4, rms 1.66 fm; the cut tails make it
    // about 8% smaller (the separation density itself would give 3.32 fm)
    const double ha = 0.228, hb = 1.18;
    auto moment = [&](int n) {
        return std::tgamma(n + 1.)
               * (std::pow(2. * ha, -n - 1.) - 2. * std::pow(ha + hb, -n - 1.)
                  + std::pow(2. * hb, -n - 1.));
    };
    const double expectedRms = std::sqrt(2. / 3. * moment(2) / moment(0) / 4.);
    CHECK(rms < expectedRms);
    CHECK(rms > 0.9 * expectedRms);
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

TEST_CASE(
    "Init::init: with colorChargeSeed the color charges come from their own "
    "stream, seeded with it and the event id (issue #127)") {
    const int N = 4;
    // the projectile's Wilson lines of one event with constant color charges
    auto wilsonLines = [&](long long seed, long long colorChargeSeed,
                           int eventId) {
        Parameters param;
        makeInitTestParam(param, N);
        param.collision.useNucleus = false;
        param.collision.g2muGeV = 1.;
        param.coupling.g = 1.;
        param.colorCharge.Ny = 2;
        param.subnucleon.m = 0.2;
        param.random.seed = seed;
        param.random.colorChargeSeed = colorChargeSeed;
        param.event.eventId = eventId;
        Lattice lat(&param, N);
        int nn[2] = {N, N};
        Init init(nn);
        Glauber glauber;
        Random random;
        random.init_genrand64(static_cast<unsigned long long>(seed));
        init.init(
            &lat, &param, &random, &glauber,
            InitializationMethod::SampleColorCharges);
        return lat.U;
    };
    auto same = [](const std::vector<Matrix> &a, const std::vector<Matrix> &b) {
        for (std::size_t i = 0; i < a.size(); ++i) {
            if (!matricesClose(a[i], b[i], 1e-15)) return false;
        }
        return true;
    };
    // the main stream (seed) no longer matters, colorChargeSeed and the
    // event do
    CHECK(same(wilsonLines(1, 5, 0), wilsonLines(2, 5, 0)));
    CHECK_FALSE(same(wilsonLines(1, 5, 0), wilsonLines(1, 6, 0)));
    CHECK_FALSE(same(wilsonLines(1, 5, 0), wilsonLines(1, 5, 1)));
    // -1: the color charges come from the stream of the run
    CHECK_FALSE(same(wilsonLines(1, -1, 0), wilsonLines(2, -1, 0)));
    CHECK_FALSE(same(wilsonLines(1, -1, 0), wilsonLines(1, 5, 0)));
    CHECK(Init::colorChargeRandomSeed(5, 3) == 5ULL + (3ULL << 32));
}

TEST_CASE(
    "Init::largeXQs2Factor: the extrapolation of Q_s^2 to x > 0.01 with "
    "useFluctuatingX 1 is 1 at x = 0.01, falls with x and is 0 for x >= 1") {
    CHECK(Init::largeXQs2Factor(0.01) == doctest::Approx(1.));
    CHECK(
        Init::largeXQs2Factor(0.1)
        == doctest::Approx(std::pow(0.9 / 0.99, 5.6) * std::pow(0.1, 0.2)));
    CHECK(Init::largeXQs2Factor(0.5) < Init::largeXQs2Factor(0.1));
    CHECK(Init::largeXQs2Factor(0.999) > 0.);
    // no gluons at x >= 1: zero, not NaN, so that the cell gets no color
    // charge and the fluctuating-x iteration ends there
    CHECK(Init::largeXQs2Factor(1.) == 0.);
    CHECK(Init::largeXQs2Factor(1.5) == 0.);
    CHECK(Init::largeXQs2Factor(40.) == 0.);
}

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "CollisionGeometry.h"
#include "Glauber.h"
#include "Init.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "PhysConst.h"
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

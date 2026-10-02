#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

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
    "Init::computeEffectiveRapidities passes rapidity through unchanged "
    "when usePseudoRapidity is off") {
    Parameters param;
    makeInitTestParam(param, 4);
    param.colorCharge.usePseudoRapidity = false;
    param.colorCharge.rapidityA = 1.5;
    param.colorCharge.rapidityB = -0.8;

    int nn[2] = {4, 4};
    Init init(nn);
    double rapidityA = 0., rapidityB = 0.;
    init.computeEffectiveRapidities(&param, rapidityA, rapidityB);

    CHECK(rapidityA == doctest::Approx(1.5));
    CHECK(rapidityB == doctest::Approx(-0.8));
}

TEST_CASE(
    "Init::computeEffectiveRapidities converts pseudorapidity to rapidity "
    "when enabled") {
    Parameters param;
    makeInitTestParam(param, 4);
    param.colorCharge.usePseudoRapidity = true;
    param.colorCharge.rapidityA = 1.0;
    param.colorCharge.rapidityB = 1.0;
    param.colorCharge.jacobianMass = 0.14;
    param.collision.sqrtS = 200.0;

    int nn[2] = {4, 4};
    Init init(nn);
    double rapidityA = 0., rapidityB = 0.;
    init.computeEffectiveRapidities(&param, rapidityA, rapidityB);

    // Same pseudorapidity input on both sides must give the same rapidity
    // output (the conversion formula is symmetric in A vs B), and the
    // conversion must actually have changed the value (it wouldn't at
    // pseudorapidity 0, but 1.0 is not a fixed point).
    CHECK(rapidityA == doctest::Approx(rapidityB));
    CHECK(rapidityA != doctest::Approx(1.0));
    CHECK(std::isfinite(rapidityA));
}

TEST_CASE(
    "Init::computeAndSetRunningAlphaS respects nFlavors/LambdaQCD instead "
    "of assuming 3 flavors and LambdaQCD=0.2, and uses the same muZero/c "
    "as Evolution::computeRunningCouplingGfactor()") {
    Parameters param;
    makeInitTestParam(param, 4);
    param.coupling.runningCoupling = true;
    param.coupling.runWithKt = false;
    param.coupling.runWithQs = 1;  // average Qs
    param.coupling.runningCouplingQsFactor = 0.5;
    param.event.averageQsAvg = 1.3;
    param.coupling.mu0 = 0.3;
    param.coupling.c = 0.2;

    int nn[2] = {4, 4};
    Init init(nn);

    param.coupling.nFlavors = 3;
    param.coupling.LambdaQCD = 0.2;
    init.computeAndSetRunningAlphaS(&param);
    CHECK(param.event.alphas == doctest::Approx(0.5922901359561648));

    param.coupling.nFlavors = 4;
    param.coupling.LambdaQCD = 0.25;
    init.computeAndSetRunningAlphaS(&param);
    CHECK(param.event.alphas == doctest::Approx(0.789051392062008));
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
    std::vector<double> scaleA, scaleB;
    init.computeWilsonLineColorChargeScales(
        &lat, N * N, /*g=*/2.0, /*invNy=*/0.25, scaleA, scaleB);

    REQUIRE(scaleA.size() == static_cast<std::size_t>(N * N));
    REQUIRE(scaleB.size() == static_cast<std::size_t>(N * N));
    for (int pos = 0; pos < N * N; ++pos) {
        CHECK(scaleA[pos] == doctest::Approx(2.0 * std::sqrt(4.0 * 0.25)));
        CHECK(scaleB[pos] == doctest::Approx(2.0 * std::sqrt(9.0 * 0.25)));
    }
}

TEST_CASE(
    "Init::setConstantColorChargeDensity: flat background sets g2mu2A=g2mu2B="
    "(g2mu/g)^2 everywhere") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    param.collision.useGaussian = false;
    param.collision.g2mu = 6.0;
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
    "Init::setConstantColorChargeDensity: Gaussian background peaks at the "
    "lattice center with the documented envelope") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    param.collision.useGaussian = true;
    param.collision.g2mu = 6.0;
    param.coupling.g = 2.0;
    Lattice lat(&param, N);

    int nn[2] = {N, N};
    Init init(nn);
    init.setConstantColorChargeDensity(&lat, &param);

    // Center cell: x = ix*L/N - L/2 = 0 at ix = N/2 (same for y).
    const double sigmax = 0.35, sigmay = 0.5;
    const double envelope = 1.0 / (2.0 * M_PI * sigmax * sigmay);
    const int centerPos = latticeIndex(N / 2, N / 2, N);
    const double expected = envelope * 6.0 * 6.0 / 2.0 / 2.0;
    CHECK(lat.cells[centerPos]->getg2mu2A() == doctest::Approx(expected));
    CHECK(lat.cells[centerPos]->getg2mu2B() == doctest::Approx(expected));
}

TEST_CASE(
    "Init::computeNcollList (hard-sphere mode) marks a p+p pair collided and "
    "counts it once") {
    Parameters param;
    param.collision.sigmaNN = 4.2;
    param.collision.gaussianWounding = false;
    param.event.eventId = 0;

    Glauber glauber;
    Random random;
    int nn[2] = {4, 4};
    Init init(nn);
    // A1=A2=1 ("p"): sampleTA places both nucleons at the origin (the
    // random polarization rotation leaves them there).
    glauber.initGlauber(
        4.2, "p", "p", /*inb=*/0.0, /*setWSDeformParams=*/false, 0., 0., 0., 0.,
        0., 0., /*forceDminFlag=*/false, 0., 0., 0., /*imax=*/1000);
    param.collision.nucleiToAverage = 1;
    init.sampleTA(&param, &random, &glauber);

    const double d2 = param.collision.sigmaNN / (M_PI * 10.);  // in fm^2
    int Ncoll = 0;
    init.computeNcollList(&param, d2, /*b=*/0.0, /*phiRP=*/0.0, Ncoll);

    CHECK(Ncoll == 1);
    std::remove("NcollList0.dat");
}

TEST_CASE(
    "Init::computeSmoothNucleusThickness builds centered projectile (A) and "
    "target (B) profiles, and scanCollisionGeometry finds their shifted "
    "overlap") {
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

    auto overlapCells = [&](double b) {
        double averageQs = 0., averageQs2 = 0., averageQs2Avg = 0.,
               averageQs2min = 0., averageQs2min2 = 0., Tpp = 0.;
        int count = 0;
        init.scanCollisionGeometry(
            &lat, &param, N, a, b, /*phiRP=*/0., averageQs, averageQs2,
            averageQs2Avg, averageQs2min, averageQs2min2, Tpp, count);
        return count;
    };
    const int central = overlapCells(0.);
    const int shifted = overlapCells(4.);
    CHECK(central > 0);
    CHECK(shifted > 0);
    CHECK(shifted < central);
    CHECK(overlapCells(40.) == 0);
}

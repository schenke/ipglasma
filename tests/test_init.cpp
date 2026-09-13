#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "Glauber.h"
#include "Init.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "PhysConst.h"
#include "doctest.h"

using PhysConst::hbarc;

namespace {
// Matches test_lattice.cpp's makeLatticeParam: Parameters holds a
// PrettyOstream member, which holds a non-copyable/non-movable
// std::ostringstream, so this takes an out-parameter rather than
// returning by value.
void makeInitTestParam(Parameters &param, int size) {
    param.setSize(size);
    param.setL(static_cast<double>(size));  // a = 1 fm
    param.setMPIRank(0);
    param.setEventId(0);
    param.setSeed(0);
    param.setMPISize(1);
    param.setRapidityA(0.0);
    param.setRapidityB(0.0);
}

// A deterministic, position-dependent (not merely diagonal) matrix, so
// tests exercise genuine matrix multiplication/conjugation rather than
// something a transposition or index-swap bug could accidentally pass.
Matrix makeTestMatrix(int seed) {
    Matrix m(Matrix::noInit);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            m.set(
                i, j,
                std::complex<double>(
                    0.1 * seed + 0.01 * (i + 1), 0.05 * seed - 0.02 * (j + 1)));
        }
    }
    return m;
}

bool matricesClose(const Matrix &a, const Matrix &b, double tol) {
    for (int k = 0; k < 9; ++k) {
        if (std::abs(a.get(k) - b.get(k)) >= tol) return false;
    }
    return true;
}
}  // namespace

TEST_CASE(
    "Init::computeEffectiveRapidities passes rapidity through unchanged "
    "when usePseudoRapidity is off") {
    Parameters param;
    makeInitTestParam(param, 4);
    param.setUsePseudoRapidity(0);
    param.setRapidityA(1.5);
    param.setRapidityB(-0.8);

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
    param.setUsePseudoRapidity(1);
    param.setRapidityA(1.0);
    param.setRapidityB(1.0);
    param.setJacobianm(0.14);
    param.setRoots(200.0);

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
    "Init::computeWilsonLineMomentumKernel: massless kernel is 1/kt^2, zero "
    "at kt=0") {
    const int N = 4;
    int nn[2] = {N, N};
    Init init(nn);

    std::vector<double> kernel =
        init.computeWilsonLineMomentumKernel(N, N * N, /*m=*/0.0, /*UVdamp=*/0.0);
    REQUIRE(kernel.size() == static_cast<std::size_t>(N * N));

    // pos=0 -> (i=0,j=0): kx=ky=-pi, kt2 = 4*(sin(-pi/2)^2*2) = 8.
    CHECK(kernel[0] == doctest::Approx(1.0 / 8.0));
    // pos = (N/2)*N + N/2 -> (i=j=N/2): kx=ky=0, kt2=0 -> defined as 0.
    CHECK(kernel[(N / 2) * N + (N / 2)] == doctest::Approx(0.0));
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
    CHECK(kernel[(N / 2) * N + (N / 2)] == doctest::Approx(1.0 / (m * m)));
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
    param.setUseGaussian(0);
    param.setg2mu(6.0);
    param.setg(2.0);
    Lattice lat(&param, N);

    int nn[2] = {N, N};
    Init init(nn);
    init.setConstantColorChargeDensity(&lat, &param);

    for (int pos = 0; pos < N * N; ++pos) {
        CHECK(lat.cells[pos]->getg2mu2A() == doctest::Approx(9.0));
        CHECK(lat.cells[pos]->getg2mu2B() == doctest::Approx(9.0));
    }
    CHECK(param.getSuccess() == 1);
}

TEST_CASE(
    "Init::setConstantColorChargeDensity: Gaussian background peaks at the "
    "lattice center with the documented envelope") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    param.setUseGaussian(1);
    param.setg2mu(6.0);
    param.setg(2.0);
    Lattice lat(&param, N);

    int nn[2] = {N, N};
    Init init(nn);
    init.setConstantColorChargeDensity(&lat, &param);

    // Center cell: x = ix*L/N - L/2 = 0 at ix = N/2 (same for y).
    const double sigmax = 0.35, sigmay = 0.5;
    const double envelope = 1.0 / (2.0 * M_PI * sigmax * sigmay);
    const int centerPos = (N / 2) * N + (N / 2);
    const double expected = envelope * 6.0 * 6.0 / 2.0 / 2.0;
    CHECK(lat.cells[centerPos]->getg2mu2A() == doctest::Approx(expected));
    CHECK(lat.cells[centerPos]->getg2mu2B() == doctest::Approx(expected));
}

TEST_CASE(
    "Init::computeNucleonThicknessAtCell: single-Gaussian branch at zero "
    "separation reduces to 1/(2*pi*BG)") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    param.setUseConstituentQuarkProton(0);
    param.setBG(1.0);

    std::vector<ReturnValue> nucleus(1);
    nucleus[0].x = 0.0;
    nucleus[0].y = 0.0;
    nucleus[0].phi = 0.0;
    nucleus[0].collided = 0;
    nucleus[0].proton = true;

    std::vector<std::vector<double>> xq(1), yq(1), BGq(1);
    std::vector<std::vector<double>> gauss(1, std::vector<double>(1, 1.0));

    int nn[2] = {N, N};
    Init init(nn);
    const double Tp = init.computeNucleonThicknessAtCell(
        &param, nucleus, xq, yq, BGq, gauss, /*x=*/0.0, /*y=*/0.0, /*xi=*/0.0,
        /*nucleiInAverage=*/1.0);

    CHECK(Tp == doctest::Approx(1.0 / (2.0 * M_PI)));
}

TEST_CASE(
    "Init::computeNucleonThicknessAtCell: constituent-quark branch with one "
    "quark at zero separation matches the single-Gaussian formula") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    param.setUseConstituentQuarkProton(1);

    std::vector<ReturnValue> nucleus(1);
    nucleus[0].x = 0.0;
    nucleus[0].y = 0.0;
    nucleus[0].phi = 0.0;
    nucleus[0].collided = 0;
    nucleus[0].proton = true;

    std::vector<std::vector<double>> xq(1, std::vector<double>(1, 0.0));
    std::vector<std::vector<double>> yq(1, std::vector<double>(1, 0.0));
    std::vector<std::vector<double>> BGq(1, std::vector<double>(1, 1.0));
    std::vector<std::vector<double>> gauss(1, std::vector<double>(1, 1.0));

    int nn[2] = {N, N};
    Init init(nn);
    const double Tp = init.computeNucleonThicknessAtCell(
        &param, nucleus, xq, yq, BGq, gauss, /*x=*/0.0, /*y=*/0.0, /*xi=*/0.0,
        /*nucleiInAverage=*/1.0);

    CHECK(Tp == doctest::Approx(1.0 / (2.0 * M_PI)));
}

TEST_CASE(
    "Init::readWilsonLineText/readWilsonLineBinary round-trip a synthetic "
    "file into lat->U") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    param.setb(3.0);

    int nn[2] = {N, N};
    Init init(nn);

    // Text format: N*N lines of "i j Re0 Im0 ... Re8 Im8" (dummies i,j are
    // discarded by the reader; loop order alone fixes position).
    const std::string textPath = "ipglasma_test_init_wilson.txt";
    {
        std::ofstream out(textPath);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                out << i << " " << j;
                for (int k = 0; k < 9; ++k) {
                    out << " " << (1000. * i + 100. * j + k) << " "
                        << -(1000. * i + 100. * j + k);
                }
                out << "\n";
            }
        }
    }

    Lattice lat(&param, N);
    init.readWilsonLineText(textPath, &param, NucleusRole::Projectile, lat.U);
    std::remove(textPath.c_str());

    // a=1 (L=N), b=3: isProjectile shifts x by -b/2=-1.5, so
    // ix = floor(i - 1.5); i=0,1 -> negative -> skipped, i=2 -> ix=0,
    // i=3 -> ix=1.
    for (int j = 0; j < N; ++j) {
        const int pos2 = 0 * N + j;  // from i=2
        CHECK(lat.U[pos2].get(0).real() == doctest::Approx(2000. + j * 100.));
        const int pos3 = 1 * N + j;  // from i=3
        CHECK(lat.U[pos3].get(0).real() == doctest::Approx(3000. + j * 100.));
    }

    // Binary format: header (N, Nc, L, a, dummy) then N*N*9 (re,im) pairs in
    // ix-outer/iy-inner block order.
    const std::string binPath = "ipglasma_test_init_wilson_bin";
    {
        std::ofstream out(binPath, std::ios::out | std::ios::binary);
        int Nc = 3;
        double L = param.getL();
        double a = L / N;
        double dummy = 0.;
        out.write(reinterpret_cast<const char *>(&N), sizeof(int));
        out.write(reinterpret_cast<const char *>(&Nc), sizeof(int));
        out.write(reinterpret_cast<const char *>(&L), sizeof(double));
        out.write(reinterpret_cast<const char *>(&a), sizeof(double));
        out.write(reinterpret_cast<const char *>(&dummy), sizeof(double));
        for (int ix = 0; ix < N; ++ix) {
            for (int iy = 0; iy < N; ++iy) {
                for (int row = 0; row < 3; ++row) {
                    for (int col = 0; col < 3; ++col) {
                        double re = 1000. * ix + 100. * iy + 10. * row + col;
                        double im = -re;
                        out.write(reinterpret_cast<const char *>(&re), sizeof(double));
                        out.write(reinterpret_cast<const char *>(&im), sizeof(double));
                    }
                }
            }
        }
    }

    Lattice lat2(&param, N);
    init.readWilsonLineBinary(binPath, &param, NucleusRole::Target, lat2.U2);
    std::remove(binPath.c_str());

    // isProjectile=false: ix = round(ixRaw + 1.5). ixRaw=0 -> ix=round(1.5)=2
    // (round-half-to-even or away-from-zero both give 2 here).
    for (int iy = 0; iy < N; ++iy) {
        const int pos = 2 * N + iy;
        CHECK(lat2.U2[pos].get(0).real() == doctest::Approx(0. + 10. * 0 + iy * 100.));
    }
}

TEST_CASE("Init::sanitizeForwardLightconeU replaces a NaN U/U2 with identity") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    Lattice lat(&param, N);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    lat.U[5].set(0, 0, std::complex<double>(nan, 0.0));
    lat.U2[7].set(4, std::complex<double>(0.0, nan));

    int nn[2] = {N, N};
    Init init(nn);
    init.sanitizeForwardLightconeU(&lat, N * N);

    const Matrix identity(1.0);
    CHECK(matricesClose(lat.U[5], identity, 1e-14));
    CHECK(matricesClose(lat.U2[7], identity, 1e-14));
    // An untouched cell must remain the identity it started as.
    CHECK(matricesClose(lat.U[0], identity, 1e-14));
}

TEST_CASE(
    "Init::computeForwardLightconeLinksTeam computes Ux1/Uy1/Ux2/Uy2 as "
    "U * conjg(U at the +x/+y neighbor)") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    Lattice lat(&param, N);
    for (int pos = 0; pos < N * N; ++pos) {
        lat.U[pos] = makeTestMatrix(pos + 1);
        lat.U2[pos] = makeTestMatrix(pos + 101);
    }

    int nn[2] = {N, N};
    Init init(nn);
    Init::ForwardLightconeLinkScratch scratch;
    init.computeForwardLightconeLinksTeam(&lat, N * N, scratch);

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
    "Init::computeForwardLightconePlaquetteTeam wires the four links in the "
    "documented order") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    Lattice lat(&param, N);
    for (int pos = 0; pos < N * N; ++pos) {
        lat.Ux[pos] = makeTestMatrix(pos + 1);
        lat.Uy[pos] = makeTestMatrix(pos + 51);
    }

    int nn[2] = {N, N};
    Init init(nn);
    Init::ForwardLightconePlaquetteScratch scratch;
    init.computeForwardLightconePlaquetteTeam(&lat, N * N, scratch);

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
    "Init::resetForwardLightconeFieldsTeam zeroes U/U2/Uy2 and resets Ux1 to "
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

    int nn[2] = {N, N};
    Init init(nn);
    init.resetForwardLightconeFieldsTeam(&lat, N * N);

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
    "Init::computeNcollList (hard-sphere mode) marks a p+p pair collided and "
    "counts it once") {
    Parameters param;
    param.setSigmaNN(4.2);
    param.setGaussianWounding(0);
    param.setEventId(0);

    Glauber glauber;
    Random random;
    int nn[2] = {4, 4};
    Init init(nn);
    // A1=A2=1 ("p"): sampleTAWoodsSaxon places both nucleons at the origin
    // deterministically, with no random draws.
    glauber.initGlauber(
        4.2, "p", "p", /*inb=*/0.0, /*setWSDeformParams=*/false, 0., 0., 0.,
        0., 0., 0., /*forceDminFlag=*/false, 0., 0., 0., /*imax=*/1000);
    param.setAverageOverNuclei(1);
    init.sampleTAWoodsSaxon(&param, &random, &glauber);

    const double d2 = param.getSigmaNN() / (M_PI * 10.);  // in fm^2
    int Ncoll = 0;
    init.computeNcollList(&param, d2, /*b=*/0.0, /*phiRP=*/0.0, Ncoll);

    CHECK(Ncoll == 1);
    std::remove("NcollList0.dat");
}

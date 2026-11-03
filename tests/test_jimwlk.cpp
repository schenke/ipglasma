// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include "Group.h"
#include "Init.h"
#include "JIMWLK.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "PhysConst.h"
#include "Random.h"
#include "WilsonLineIO.h"
#include "doctest.h"

namespace {
void makeJimwlkTestParam(Parameters &param, int size) {
    param.lattice.size = size;
    param.lattice.L = static_cast<double>(size);
    param.run.MPIRank = 0;
    param.event.eventId = 0;
    param.random.seed = 0;
    param.run.MPISize = 1;
    param.colorCharge.projectileX = 0.01;
    param.colorCharge.targetX = 0.01;
    param.jimwlk.mu0 = 0.2;
    param.jimwlk.LambdaQCD = 0.2;
    param.jimwlk.c = 0.2;
    param.coupling.nFlavors = 3;
    param.jimwlk.mass = 0.0;    // skips the Bessel mass-regulator branch
    param.jimwlk.alphaS = 0.3;  // fixed coupling, skips running-coupling
    param.jimwlk.Ds = 0.001;
}
}  // namespace

TEST_CASE(
    "JIMWLK::evolutionStep runs without crashing (regression test: "
    "VxsiVx_/VxsiVy_ used to be allocated only conditionally on a "
    "since-removed input flag, but evolutionStep needs them regardless)") {
    const int N = 8;
    Parameters param;
    makeJimwlkTestParam(param, N);

    Group group;
    Random random;
    random.init_genrand64(42ULL);
    Lattice lat(&param, N);

    JIMWLK jimwlk(param, &group, &lat, &random);
    jimwlk.evolutionStep(NucleusRole::Projectile);
    jimwlk.evolutionStep(NucleusRole::Target);

    for (int pos = 0; pos < N * N; ++pos) {
        for (int k = 0; k < 9; ++k) {
            CHECK(std::isfinite(lat.U[pos].get(k).real()));
            CHECK(std::isfinite(lat.U[pos].get(k).imag()));
            CHECK(std::isfinite(lat.U2[pos].get(k).real()));
            CHECK(std::isfinite(lat.U2[pos].get(k).imag()));
        }
    }
}

TEST_CASE(
    "JIMWLK::evolution doesn't crash when steps_1/steps_2 come out under 10 "
    "(regression test: printSteps = steps_1/10 used to be 0 in that case, "
    "making \"ids % printSteps\" an integer modulo-by-zero)") {
    const int N = 8;
    Parameters param;
    makeJimwlkTestParam(param, N);
    param.jimwlk.saveSnapshots = 0;
    param.jimwlk.xSnapshotList = std::vector<double>();
    // Fixed coupling (jimwlk.alphaS > 0): steps_1 = as*log(x0/x_proj) /
    // (pi^2*ds) + 0.5. These values give steps_1 = steps_2 = 1.
    param.jimwlk.initialX = 0.01;
    param.colorCharge.projectileX = 0.008;
    param.colorCharge.targetX = 0.008;

    Group group;
    Random random;
    random.init_genrand64(42ULL);
    Lattice lat(&param, N);

    JIMWLK jimwlk(param, &group, &lat, &random);
    jimwlk.evolution();

    for (int pos = 0; pos < N * N; ++pos) {
        for (int k = 0; k < 9; ++k) {
            CHECK(std::isfinite(lat.U[pos].get(k).real()));
            CHECK(std::isfinite(lat.U2[pos].get(k).real()));
        }
    }
}

TEST_CASE(
    "JIMWLK::getAlphas: fixed coupling returns 1 so alpha_s is not counted "
    "twice (it is already in the step count)") {
    const int N = 8;
    Parameters param;
    makeJimwlkTestParam(param, N);
    param.jimwlk.alphaS = 0.3;

    Group group;
    Random random;
    random.init_genrand64(42ULL);
    Lattice lat(&param, N);
    JIMWLK jimwlk(param, &group, &lat, &random);

    CHECK(jimwlk.getAlphas(0.1, 0.15) == doctest::Approx(1.0));
    CHECK(jimwlk.getAlphas(-0.3, 0.2) == doctest::Approx(1.0));
}

TEST_CASE(
    "JIMWLK::getAlphas: running-coupling branch matches an independently "
    "computed reference value, and is sensitive to nFlavors/jimwlkC") {
    const int N = 8;
    Parameters param;
    makeJimwlkTestParam(param, N);
    param.jimwlk.mu0 = 0.28;
    param.jimwlk.LambdaQCD = 0.04;
    param.jimwlk.alphaS = 0.0;  // forces the running-coupling branch

    Group group;
    Random random;
    random.init_genrand64(42ULL);
    Lattice lat(&param, N);

    const double x = 0.1;
    const double y = 0.15;

    // JIMWLK stores Parameters by reference, so getAlphas() always reads
    // whatever param currently holds -- read each value right after
    // setting it, before mutating param again.
    JIMWLK jimwlk(param, &group, &lat, &random);

    const double alphasDefault = jimwlk.getAlphas(x, y);
    CHECK(alphasDefault == doctest::Approx(0.3482995062039298));

    param.coupling.nFlavors = 4;
    const double alphasNf4 = jimwlk.getAlphas(x, y);
    CHECK(alphasNf4 == doctest::Approx(0.37616346670024414));
    CHECK(alphasNf4 > alphasDefault);
    param.coupling.nFlavors = 3;

    param.jimwlk.c = 0.25;
    const double alphasC025 = jimwlk.getAlphas(x, y);
    CHECK(alphasC025 == doctest::Approx(0.3453365551834896));
}

TEST_CASE(
    "JIMWLK::snapshotSteps saves each x at the closest step and skips x "
    "outside the evolution") {
    // x_k = 0.01 e^{-0.1 k}, k = 0..20
    const double x0 = 0.01;
    const double dlogx = 0.1;
    const int steps = 20;
    const std::vector<double> requested = {
        x0,                    // before the first step
        x0 * std::exp(-0.53),  // closer to step 5 than to 6
        x0 * std::exp(-0.57),  // closer to step 6
        x0 * std::exp(-2.04),  // half a step past the end at most: step 20
        x0 * std::exp(0.04),   // just above x0: step 0
        x0 * std::exp(-2.2),   // beyond the end
        x0 * std::exp(0.2),    // above x0
        0.};
    CHECK(
        JIMWLK::snapshotSteps(requested, x0, dlogx, steps)
        == std::vector<int> {0, 5, 6, 20, 0, -1, -1, -1});
}

TEST_CASE(
    "JIMWLK::evolution writes the snapshots between the initial and the "
    "final x, but not a second copy of the initial or final Wilson lines "
    "(issue #120)") {
    const int N = 8;
    Parameters param;
    makeJimwlkTestParam(param, N);
    const std::filesystem::path dir =
        std::filesystem::temp_directory_path() / "ipglasma_test_snapshots";
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    param.wilsonLines.wilsonLinePath = dir.string();
    param.wilsonLines.writeWilsonLines = 2;
    param.run.eventsPerRank = 1;
    param.jimwlk.enabled = true;  // the initial x is jimwlkInitialX
    param.jimwlk.saveSnapshots = true;
    // 7 steps from 0.01 to 0.008 for both nuclei; 0.009 lies in between,
    // and 0.01 and 0.008 are the initial and the final x, whose Wilson
    // lines Init::setV() and main write under the same names
    param.jimwlk.initialX = 0.01;
    param.colorCharge.projectileX = 0.008;
    param.colorCharge.targetX = 0.008;
    param.jimwlk.xSnapshotList = {0.01, 0.009, 0.008};

    Group group;
    Random random;
    random.init_genrand64(42ULL);
    Lattice lat(&param, N);
    JIMWLK jimwlk(param, &group, &lat, &random);
    jimwlk.evolution();

    for (const NucleusRole nucleus :
         {NucleusRole::Projectile, NucleusRole::Target}) {
        CHECK(std::filesystem::exists(
            WilsonLineIO::fileName(&param, 0.009, nucleus)));
        CHECK_FALSE(std::filesystem::exists(
            WilsonLineIO::fileName(&param, 0.01, nucleus)));
        CHECK_FALSE(std::filesystem::exists(
            WilsonLineIO::fileName(&param, 0.008, nucleus)));
    }
    std::filesystem::remove_all(dir);
}

TEST_CASE(
    "JIMWLK: with colorChargeSeed the noise comes from a stream of its own, "
    "seeded by colorChargeSeed and the event, not from that of the run") {
    const int N = 8;
    // one step from the same nontrivial Wilson lines (the step leaves unit
    // ones unchanged): the Wilson line of the projectile at one site, and
    // the next number of the run's stream after the step
    struct Step {
        std::complex<double> wilsonLine;
        unsigned long long nextOfRun;
    };
    auto evolved = [&](unsigned long long runSeed, long long colorChargeSeed,
                       int eventId) {
        Parameters param;
        makeJimwlkTestParam(param, N);
        param.random.colorChargeSeed = colorChargeSeed;
        param.event.eventId = eventId;
        Group group;
        Random random;
        random.init_genrand64(runSeed);
        Lattice lat(&param, N);
        for (int pos = 0; pos < N * N; ++pos) {
            std::vector<double> q(8);
            for (int a = 0; a < 8; ++a) q[a] = std::sin(1.3 * pos + 0.7 * a);
            lat.U[pos] = Matrix::fromAlgebraExponent(q);
        }
        JIMWLK jimwlk(param, &group, &lat, &random);
        jimwlk.evolutionStep(NucleusRole::Projectile);
        return Step {lat.U[3].get(1), random.genrand64_int64()};
    };
    // the next number of the stream of runSeed after count Gaussian ones
    auto nextAfter = [](unsigned long long runSeed, std::size_t count) {
        Random random;
        random.init_genrand64(runSeed);
        std::vector<double> gaussians(count);
        std::vector<double> scratch;
        random.gaussBulk(gaussians.data(), count, scratch);
        return random.genrand64_int64();
    };
    // the FFT plans can differ in the last bits between two runs
    auto same = [](const Step &a, const Step &b) {
        return std::abs(a.wilsonLine - b.wilsonLine) < 1e-12;
    };
    // -1: the noise comes from the stream of the run, which the step
    // advances by its 2 (Nc^2 - 1) N^2 Gaussian numbers
    CHECK_FALSE(same(evolved(1, -1, 0), evolved(2, -1, 0)));
    CHECK(
        evolved(1, -1, 0).nextOfRun
        == nextAfter(
            1, static_cast<std::size_t>(2 * PhysConst::Nc2m1 * N * N)));
    // colorChargeSeed: the seed of the run no longer matters, the
    // colorChargeSeed and the event do, and the step leaves the run's stream
    // alone
    CHECK(same(evolved(1, 5, 0), evolved(2, 5, 0)));
    CHECK_FALSE(same(evolved(1, 5, 0), evolved(1, 6, 0)));
    CHECK_FALSE(same(evolved(1, 5, 0), evolved(1, 5, 1)));
    CHECK(evolved(1, 5, 0).nextOfRun == nextAfter(1, 0));
    // exactly the noise of a run whose stream is seeded with
    // noiseRandomSeed(), a stream other than that of the color charges
    CHECK(
        same(evolved(1, 5, 3), evolved(JIMWLK::noiseRandomSeed(5, 3), -1, 3)));
    CHECK(
        JIMWLK::noiseRandomSeed(5, 3)
        == splitmix64(Init::colorChargeRandomSeed(5, 3)));
    CHECK(JIMWLK::noiseRandomSeed(5, 3) != Init::colorChargeRandomSeed(5, 3));
}

TEST_CASE(
    "JIMWLK::getAlphasMomentum: the noise coupling alpha_s(k) of the lattice "
    "momentum matches an independently computed reference value") {
    const int N = 8;
    Parameters param;
    makeJimwlkTestParam(param, N);  // L = 8 fm, so a = 1 fm
    param.jimwlk.mu0 = 0.28;
    param.jimwlk.LambdaQCD = 0.04;
    param.jimwlk.alphaS = 0.0;
    param.jimwlk.runningCoupling = 2;

    Group group;
    Random random;
    random.init_genrand64(42ULL);
    Lattice lat(&param, N);
    JIMWLK jimwlk(param, &group, &lat, &random);

    // k~^2 = (2/a)^2 [sin^2(pi x) + sin^2(pi y)]
    CHECK(
        jimwlk.getAlphasMomentum(0.25, -0.125)
        == doctest::Approx(0.33305093593417584));
    // alpha_s(-k) = alpha_s(k), so the filtered noise stays real
    CHECK(
        jimwlk.getAlphasMomentum(-0.25, 0.125)
        == doctest::Approx(jimwlk.getAlphasMomentum(0.25, -0.125)));
    // the infrared value, set by mu0
    CHECK(
        jimwlk.getAlphasMomentum(0., 0.)
        == doctest::Approx(0.35876872379572905));
    param.jimwlk.c = 0.25;
    CHECK(
        jimwlk.getAlphasMomentum(0.25, -0.125)
        == doctest::Approx(0.3308474501922769));
}

TEST_CASE(
    "JIMWLK: with a constant alpha_s the noise and the square-root running "
    "couplings give the same evolution") {
    // With mu0 far above every momentum on the lattice, alpha_s(r) and
    // alpha_s(k) both freeze to 4 pi / (beta ln(mu0^2/Lambda^2)). The noise
    // coupling then only multiplies the noise by sqrt(alpha_s), and the
    // square-root coupling the kernel, so the two steps agree for the
    // same random numbers.
    const int N = 8;
    std::vector<Matrix> start;
    std::vector<Matrix> result[2];
    for (const int prescription : {1, 2}) {
        Parameters param;
        makeJimwlkTestParam(param, N);
        param.jimwlk.mu0 = 100.;
        param.jimwlk.LambdaQCD = 0.2;
        param.jimwlk.alphaS = 0.0;
        param.jimwlk.Ds = 0.005;
        param.jimwlk.runningCoupling = prescription;

        Group group;
        Random random;
        random.init_genrand64(42ULL);
        Lattice lat(&param, N);
        // the identity is a fixed point of the evolution, so start from
        // the same nontrivial SU(3) matrices exp(i w^a t^a) for both
        Random initial;
        initial.init_genrand64(7ULL);
        for (int pos = 0; pos < N * N; ++pos) {
            Matrix generator(0.);
            for (int a = 0; a < 8; a++) {
                generator += initial.gauss() * group.getT(a);
            }
            lat.U[pos] = (std::complex<double>(0., 1.) * generator).expm();
        }
        start = lat.U;
        JIMWLK jimwlk(param, &group, &lat, &random);
        for (int step = 0; step < 3; step++) {
            jimwlk.evolutionStep(NucleusRole::Projectile);
        }
        result[prescription - 1] = lat.U;
    }
    double maxDiff = 0.;
    double maxChange = 0.;
    for (int pos = 0; pos < N * N; ++pos) {
        for (int k = 0; k < 9; ++k) {
            maxDiff = std::max(
                maxDiff,
                std::abs(result[0][pos].get(k) - result[1][pos].get(k)));
            maxChange = std::max(
                maxChange, std::abs(result[1][pos].get(k) - start[pos].get(k)));
        }
    }
    CHECK(maxChange > 1e-3);  // the Wilson lines did evolve
    CHECK(maxDiff < 1e-10);
}

TEST_CASE(
    "JIMWLK: the noise of the noise running coupling has the correlator of "
    "Eq. (24) of Lappi:2012vw on the lattice") {
    // eta is linear in xi, eta_x = sum_z M_xz xi_z, so for white noise
    // <eta_x eta_y> = sum_z M_xz M_yz. On the lattice (noise in units of
    // the lattice spacing) Eq. (24) is (1/N^2) sum_k alpha_s(k) e^{ik(x-y)}
    // over the N^2 momenta, with alpha_s(k) of Eq. (39) there.
    for (const int N : {8, 16}) {
        CAPTURE(N);
        Parameters param;
        makeJimwlkTestParam(param, N);
        // the setup of Fig. 6 there, in units of the lattice size L
        const double L = param.lattice.L * PhysConst::invHbarc;  // 1/GeV
        param.jimwlk.mu0 = 15. / L;
        param.jimwlk.LambdaQCD = 6. / L;
        param.jimwlk.alphaS = 0.0;
        param.jimwlk.runningCoupling = 2;

        Group group;
        Random random;
        random.init_genrand64(42ULL);
        Lattice lat(&param, N);
        JIMWLK jimwlk(param, &group, &lat, &random);

        // the infrared value quoted in the paper
        CHECK(
            jimwlk.getAlphasMomentum(0., 0.)
            == doctest::Approx(0.7619).epsilon(1e-4));

        std::vector<std::vector<double>> columns;  // columns[z][x] = M_xz
        for (int z = 0; z < N * N; ++z) {
            std::vector<double> xi(N * N, 0.);
            xi[z] = 1.;
            columns.push_back(jimwlk.correlatedNoise(xi));
        }

        // alpha_s(k) written as in Eq. (39), of the lattice momentum
        const double a = L / N;
        const double beta = 11. - 2. * param.coupling.nFlavors / 3.;
        const double c = param.jimwlk.c;
        auto alphas = [&](int nx, int ny) {
            const double k2 = 4. / (a * a)
                              * (std::pow(std::sin(M_PI * nx / N), 2)
                                 + std::pow(std::sin(M_PI * ny / N), 2));
            const double Lambda2 =
                param.jimwlk.LambdaQCD * param.jimwlk.LambdaQCD;
            return 4. * M_PI
                   / (beta
                      * std::log(
                          std::pow(
                              std::pow(
                                  param.jimwlk.mu0 * param.jimwlk.mu0 / Lambda2,
                                  1. / c)
                                  + std::pow(k2 / Lambda2, 1. / c),
                              c)));
        };
        double maxDiff = 0.;
        for (int x = 0; x < N * N; ++x) {
            for (int y = 0; y < N * N; ++y) {
                double correlator = 0.;
                for (int z = 0; z < N * N; ++z) {
                    correlator += columns[z][x] * columns[z][y];
                }
                const int dx = x / N - y / N;
                const int dy = x % N - y % N;
                double expected = 0.;
                for (int nx = -N / 2; nx < N / 2; ++nx) {
                    for (int ny = -N / 2; ny < N / 2; ++ny) {
                        expected +=
                            alphas(nx, ny)
                            * std::cos(2. * M_PI * (nx * dx + ny * dy) / N);
                    }
                }
                expected /= N * N;
                maxDiff = std::max(maxDiff, std::abs(correlator - expected));
            }
        }
        CHECK(maxDiff < 1e-12);
    }
}

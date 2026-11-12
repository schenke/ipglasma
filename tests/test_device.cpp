// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

// The device layer (src/device/) against the CPU code: the FFT, the JIMWLK
// step and the Yang-Mills leapfrog, for every backend this executable was
// built with, in double and in single precision. The host backend is
// always there; CUDA and OpenCL devices are tested when the build has them
// and the machine has a device (otherwise the case is skipped with a
// message).

#include <math.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "Evolution.h"
#include "FFT.h"
#include "Group.h"
#include "JIMWLK.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "Random.h"
#include "device/DeviceBackend.h"
#include "device/DeviceDrivers.h"
#include "doctest.h"
#include "test_helpers.h"

namespace {

// the device kernels on the host, for direct tests of their helpers
namespace kernels {
#define IPG_REAL double
#define IPG_FN static inline
#define IPG_GLOBAL
#include "device/ipg_kernels.h"
}  // namespace kernels

struct DeviceCase {
    std::string label;
    std::unique_ptr<ipgdev::Backend> backend;
    // the agreement with the CPU code expected of one step: rounding for
    // double, the single-precision rounding of O(1) numbers otherwise
    double tolerance;
};

// every backend and precision that can be created on this machine
std::vector<DeviceCase> deviceCases() {
    std::vector<DeviceCase> cases;
    for (const std::string &name : ipgdev::availableBackends()) {
        for (int precision = 1; precision <= 2; precision++) {
            ipgdev::Config config;
            config.backend = name;
            config.precision = precision;
            const std::string label =
                name + (precision == 1 ? "/single" : "/double");
            try {
                DeviceCase deviceCase;
                deviceCase.label = label;
                deviceCase.backend = ipgdev::createBackend(config);
                deviceCase.tolerance = (precision == 1) ? 2e-5 : 1e-11;
                cases.push_back(std::move(deviceCase));
            } catch (const std::exception &e) {
                std::cout << "[test_device] skipping " << label << ": "
                          << e.what() << std::endl;
            }
        }
    }
    return cases;
}

// IPGLASMA_TEST_VERBOSE=1 prints the differences the checks compare
void report(const std::string &what, double difference) {
    if (std::getenv("IPGLASMA_TEST_VERBOSE") != nullptr) {
        std::cout << "[test_device] " << what << ": " << difference
                  << std::endl;
    }
}

double maxDifference(
    const std::vector<Matrix> &a, const std::vector<Matrix> &b) {
    double difference = 0.;
    for (std::size_t i = 0; i < a.size(); i++) {
        for (int k = 0; k < 9; k++) {
            difference =
                std::max(difference, std::abs(a[i].get(k) - b[i].get(k)));
        }
    }
    return difference;
}

// the largest element of U^dagger U - 1 and of det U - 1 over a field
double maxNonUnitarity(const std::vector<Matrix> &field) {
    double deviation = 0.;
    for (const Matrix &u : field) {
        Matrix product = u;
        product.conjg();
        product = product * u;
        for (int k = 0; k < 9; k++) {
            const double unit = (k % 4 == 0) ? 1. : 0.;
            deviation = std::max(deviation, std::abs(product.get(k) - unit));
        }
        deviation = std::max(deviation, std::abs(u.det() - 1.));
    }
    return deviation;
}

// a Hermitian traceless matrix with elements of size `width`
Matrix randomAlgebraElement(std::mt19937 &engine, Group &group, double width) {
    std::normal_distribution<double> gauss(0., width);
    Matrix m(0.);
    for (int a = 0; a < 8; a++) m += gauss(engine) * group.getT(a);
    return m;
}

void makeDeviceTestParam(Parameters &param, int size) {
    param.lattice.size = size;
    param.lattice.L = 0.1 * size;
    param.run.MPIRank = 0;
    param.run.MPISize = 1;
    param.event.eventId = 0;
    param.random.seed = 0;
    param.colorCharge.projectileX = 0.01;
    param.colorCharge.targetX = 0.01;
    param.coupling.g = 1.0;
    param.coupling.nFlavors = 3;
    param.run.dtau = 0.1;
    param.jimwlk.mu0 = 0.28;
    param.jimwlk.LambdaQCD = 0.04;
    param.jimwlk.c = 0.2;
    param.jimwlk.mass = 0.4;
    param.jimwlk.alphaS = 0.;
    param.jimwlk.runningCoupling = 1;
    param.jimwlk.Ds = 0.003;
}

}  // namespace

TEST_CASE("fftRadices factorizes the lattice sizes of the device FFT") {
    auto product = [](const std::vector<int> &radices) {
        int result = 1;
        for (const int r : radices) result *= r;
        return result;
    };
    for (const int n : {2,  4,  8,  16, 32, 64, 128, 256, 512, 1024, 6,
                        10, 12, 14, 20, 26, 30, 100, 340, 600, 800,  1536}) {
        const std::vector<int> radices = ipgdev::fftRadices(n);
        CHECK(product(radices) == n);
        for (const int r : radices) CHECK((r >= 2 && r <= 32));
    }
    CHECK(ipgdev::fftRadices(800) == std::vector<int>({8, 4, 5, 5}));
    CHECK(ipgdev::fftRadices(1024) == std::vector<int>({4, 4, 4, 4, 4}));
    // a prime factor above 32 is not supported
    CHECK(ipgdev::fftRadices(2 * 37).empty());
    CHECK(ipgdev::fftRadices(4 * 41 * 3).empty());
}

TEST_CASE("the device FFT is the centered FFT of FFT::fftnArray") {
    for (DeviceCase &device : deviceCases()) {
        CAPTURE(device.label);
        // sizes with every radix: 8, 4, 2, 3, 5, 7, 11, 13
        for (const int n : {8, 12, 20, 16, 30, 32, 14, 22, 26, 40}) {
            CAPTURE(n);
            const int nn[2] = {n, n};
            const int ntot = n * n;
            const int mDim = 3;
            FFT fft(nn);
            std::mt19937 engine(1234 + n);
            std::normal_distribution<double> gauss(0., 1.);

            std::vector<std::complex<double>> data(ntot * mDim);
            std::vector<std::complex<double>> reference(ntot * mDim);
            std::vector<std::complex<double> *> dataPtr(ntot), refPtr(ntot);
            for (int pos = 0; pos < ntot; pos++) {
                dataPtr[pos] = &data[pos * mDim];
                refPtr[pos] = &reference[pos * mDim];
            }
            for (auto &value : data) {
                value = std::complex<double>(gauss(engine), gauss(engine));
            }

            ipgdev::Backend &backend = *device.backend;
            const ipgdev::Buffer buffer = backend.createBuffer(2 * ntot * mDim);
            std::vector<std::complex<double>> result(ntot * mDim);
            // the FFT sums n^2 numbers of size 1
            const double tolerance = device.tolerance * n;
            for (const int isign : {1, -1}) {
                CAPTURE(isign);
                fft.fftnArray(dataPtr.data(), refPtr.data(), nn, isign, mDim);
                backend.upload(
                    buffer, reinterpret_cast<const double *>(data.data()),
                    2 * ntot * mDim);
                backend.centeredFft(buffer, n, mDim, isign);
                backend.download(
                    buffer, reinterpret_cast<double *>(result.data()),
                    2 * ntot * mDim);
                double difference = 0.;
                double size = 0.;
                for (int i = 0; i < ntot * mDim; i++) {
                    difference = std::max(
                        difference, std::abs(result[i] - reference[i]));
                    size = std::max(size, std::abs(reference[i]));
                }
                report(
                    device.label + " FFT n=" + std::to_string(n)
                        + " isign=" + std::to_string(isign) + " relative",
                    difference / size);
                CHECK(size > 0.1);
                CHECK(difference < tolerance * (isign == 1 ? 1. : 1. / n));
            }
            backend.releaseBuffer(buffer);
        }
    }
}

TEST_CASE("the JIMWLK step on a device agrees with JIMWLK::evolutionStep") {
    struct Coupling {
        const char *name;
        double alphaS;
        int runningCoupling;
        double mass;
    };
    const std::vector<Coupling> couplings = {
        {"fixed coupling, no mass", 0.3, 1, 0.},
        {"square-root running coupling, mass", 0., 1, 0.4},
        {"noise running coupling, mass", 0., 2, 0.4},
    };
    for (DeviceCase &device : deviceCases()) {
        CAPTURE(device.label);
        for (const Coupling &coupling : couplings) {
            CAPTURE(coupling.name);
            // 2^2 * 5: radices 4 and 5
            const int N = 20;
            Parameters param;
            makeDeviceTestParam(param, N);
            param.jimwlk.alphaS = coupling.alphaS;
            param.jimwlk.runningCoupling = coupling.runningCoupling;
            param.jimwlk.mass = coupling.mass;

            Group group;
            Lattice cpuLattice(&param, N);
            Lattice deviceLattice(&param, N);
            std::mt19937 engine(99);
            for (int pos = 0; pos < N * N; pos++) {
                cpuLattice.U[pos] = randomSU3(engine, 1.0);
                deviceLattice.U[pos] = cpuLattice.U[pos];
            }
            // the same random numbers for both
            Random cpuRandom, deviceRandom;
            cpuRandom.init_genrand64(4711ULL);
            deviceRandom.init_genrand64(4711ULL);
            JIMWLK cpuJimwlk(param, &group, &cpuLattice, &cpuRandom);
            JIMWLK deviceJimwlk(param, &group, &deviceLattice, &deviceRandom);

            std::unique_ptr<ipgdev::JimwlkDevice> deviceEvolution =
                deviceJimwlk.makeDevice(*device.backend);
            deviceEvolution->upload(deviceLattice.U);
            const std::vector<Matrix> initial = cpuLattice.U;
            const int steps = 5;
            for (int step = 0; step < steps; step++) {
                cpuJimwlk.evolutionStep(NucleusRole::Projectile);
                deviceJimwlk.evolutionStepOnDevice(*deviceEvolution);
            }
            deviceEvolution->download(deviceLattice.U);

            // the steps moved the Wilson lines ...
            CHECK(maxDifference(cpuLattice.U, initial) > 1e-2);
            // ... to the same place
            report(
                device.label + " JIMWLK 5 steps, " + coupling.name,
                maxDifference(cpuLattice.U, deviceLattice.U));
            CHECK(
                maxDifference(cpuLattice.U, deviceLattice.U)
                < steps * device.tolerance);
            // the Wilson lines come back unitary to the rounding of the
            // host, also from single precision (the driver projects them)
            report(
                device.label + " JIMWLK non-unitarity",
                maxNonUnitarity(deviceLattice.U));
            CHECK(maxNonUnitarity(deviceLattice.U) < 1e-13);
            // both consumed the same random numbers
            CHECK(
                cpuRandom.genrand64_int64() == deviceRandom.genrand64_int64());
        }
    }
}

TEST_CASE("the Yang-Mills leapfrog on a device agrees with Evolution") {
    for (DeviceCase &device : deviceCases()) {
        CAPTURE(device.label);
        // N = 96 is a launch over thousands of sites: an OpenCL driver for
        // CPUs then forms large work groups, whose private variables once
        // overflowed the stack (the backend now sets the group size)
        struct Setup {
            LatticeBoundary boundary;
            int N;
        };
        for (const Setup &setup :
             {Setup {LatticeBoundary::Clamped, 12},
              Setup {LatticeBoundary::Periodic, 12},
              Setup {LatticeBoundary::Clamped, 96}}) {
            const LatticeBoundary boundary = setup.boundary;
            const int N = setup.N;
            CAPTURE(N);
            CAPTURE(boundary == LatticeBoundary::Periodic);
            Parameters param;
            makeDeviceTestParam(param, N);
            param.coupling.g = 1.3;
            const int nn[2] = {N, N};

            Group group;
            Lattice cpuLattice(&param, N, boundary);
            Lattice deviceLattice(&param, N, boundary);
            std::mt19937 engine(2026);
            for (int pos = 0; pos < N * N; pos++) {
                cpuLattice.Ux[pos] = randomSU3(engine, 0.7);
                cpuLattice.Uy[pos] = randomSU3(engine, 0.7);
                cpuLattice.U[pos] = randomAlgebraElement(engine, group, 0.3);
                cpuLattice.U2[pos] = randomAlgebraElement(engine, group, 0.3);
                cpuLattice.Ux2[pos] = randomAlgebraElement(engine, group, 0.3);
                cpuLattice.Uy2[pos] = randomAlgebraElement(engine, group, 0.3);
            }
            const std::vector<Matrix> initialE = cpuLattice.U;
            const std::vector<Matrix> initialUx = cpuLattice.Ux;
            deviceLattice.U = cpuLattice.U;
            deviceLattice.U2 = cpuLattice.U2;
            deviceLattice.Ux = cpuLattice.Ux;
            deviceLattice.Uy = cpuLattice.Uy;
            deviceLattice.Ux2 = cpuLattice.Ux2;
            deviceLattice.Uy2 = cpuLattice.Uy2;

            Evolution evolution(nn);
            ipgdev::YangMillsDevice deviceEvolution(
                *device.backend, N, boundary == LatticeBoundary::Periodic,
                param.coupling.g);
            deviceEvolution.upload(deviceLattice);

            // the sequence of Evolution::run(): the initial half step of
            // the coordinates, full steps, a measurement half step from a
            // backup, and the final half step of the momenta
            const double dtau = param.run.dtau;
            evolution.evolvePhi(&cpuLattice, &param, dtau, 0.);
            evolution.evolveU(&cpuLattice, &param, dtau, 0.);
            deviceEvolution.evolvePhi(dtau, 0.);
            deviceEvolution.evolveU(dtau, 0.);
            const int steps = 6;
            for (int it = 1; it < steps; it++) {
                const double tau = it * dtau;
                evolution.evolvePi(&cpuLattice, &param, dtau, tau);
                evolution.evolveE(&cpuLattice, &param, dtau, tau);
                evolution.evolvePhi(&cpuLattice, &param, dtau, tau);
                evolution.evolveU(&cpuLattice, &param, dtau, tau);
                deviceEvolution.step(dtau, tau, true);
            }
            deviceEvolution.backupMomenta();
            deviceEvolution.step(dtau / 2., steps * dtau, false);
            deviceEvolution.restoreMomenta();
            evolution.evolvePi(&cpuLattice, &param, dtau / 2., steps * dtau);
            evolution.evolveE(&cpuLattice, &param, dtau / 2., steps * dtau);
            deviceEvolution.step(dtau / 2., steps * dtau, false);
            deviceEvolution.download(deviceLattice);

            CHECK(maxDifference(cpuLattice.U, initialE) > 1e-2);
            CHECK(maxDifference(cpuLattice.Ux, initialUx) > 1e-2);
            const double tolerance = steps * device.tolerance;
            report(
                device.label + " Yang-Mills 6 steps, E1",
                maxDifference(cpuLattice.U, deviceLattice.U));
            report(
                device.label + " Yang-Mills 6 steps, Ux",
                maxDifference(cpuLattice.Ux, deviceLattice.Ux));
            CHECK(maxDifference(cpuLattice.U, deviceLattice.U) < tolerance);
            CHECK(maxDifference(cpuLattice.U2, deviceLattice.U2) < tolerance);
            CHECK(maxDifference(cpuLattice.Ux, deviceLattice.Ux) < tolerance);
            CHECK(maxDifference(cpuLattice.Uy, deviceLattice.Uy) < tolerance);
            CHECK(maxDifference(cpuLattice.Ux2, deviceLattice.Ux2) < tolerance);
            CHECK(maxDifference(cpuLattice.Uy2, deviceLattice.Uy2) < tolerance);
            // the links come back unitary, see the JIMWLK test
            CHECK(maxNonUnitarity(deviceLattice.Ux) < 1e-13);
            CHECK(maxNonUnitarity(deviceLattice.Uy) < 1e-13);
        }
    }
}

TEST_CASE("configFromEnvironment reads the device of the run") {
    auto withEnvironment = [](const char *device, const char *precision,
                              const char *stages) {
        if (device != nullptr) {
            setenv("IPGLASMA_DEVICE", device, 1);
        } else {
            unsetenv("IPGLASMA_DEVICE");
        }
        if (precision != nullptr) {
            setenv("IPGLASMA_DEVICE_PRECISION", precision, 1);
        } else {
            unsetenv("IPGLASMA_DEVICE_PRECISION");
        }
        if (stages != nullptr) {
            setenv("IPGLASMA_DEVICE_STAGES", stages, 1);
        } else {
            unsetenv("IPGLASMA_DEVICE_STAGES");
        }
    };
    withEnvironment(nullptr, nullptr, nullptr);
    CHECK_FALSE(ipgdev::configFromEnvironment().enabled());
    withEnvironment("cpu", nullptr, nullptr);
    CHECK_FALSE(ipgdev::configFromEnvironment().enabled());

    withEnvironment("CUDA", nullptr, nullptr);
    ipgdev::Config config = ipgdev::configFromEnvironment();
    CHECK(config.backend == "cuda");
    CHECK(config.precision == 0);
    CHECK(config.jimwlk);
    CHECK(config.yangMills);

    withEnvironment("opencl", "single", "jimwlk");
    config = ipgdev::configFromEnvironment();
    CHECK(config.backend == "opencl");
    CHECK(config.precision == 1);
    CHECK(config.jimwlk);
    CHECK_FALSE(config.yangMills);

    withEnvironment("host", "double", "yangmills");
    config = ipgdev::configFromEnvironment();
    CHECK(config.precision == 2);
    CHECK_FALSE(config.jimwlk);
    CHECK(config.yangMills);

    withEnvironment("metal", nullptr, nullptr);
    CHECK_THROWS(ipgdev::configFromEnvironment());
    withEnvironment("host", "quad", nullptr);
    CHECK_THROWS(ipgdev::configFromEnvironment());
    withEnvironment("host", nullptr, "tmunu");
    CHECK_THROWS(ipgdev::configFromEnvironment());

    // the noise of the JIMWLK evolution
    withEnvironment("host", nullptr, nullptr);
    CHECK_FALSE(ipgdev::configFromEnvironment().noiseOnDevice);
    setenv("IPGLASMA_DEVICE_NOISE", "device", 1);
    CHECK(ipgdev::configFromEnvironment().noiseOnDevice);
    setenv("IPGLASMA_DEVICE_NOISE", "host", 1);
    CHECK_FALSE(ipgdev::configFromEnvironment().noiseOnDevice);
    setenv("IPGLASMA_DEVICE_NOISE", "gpu", 1);
    CHECK_THROWS(ipgdev::configFromEnvironment());
    unsetenv("IPGLASMA_DEVICE_NOISE");
    withEnvironment(nullptr, nullptr, nullptr);
}

TEST_CASE("ipg_threefry2x32 reproduces the known answers of Random123") {
    // kat_vectors of the Random123 library, threefry2x32 with 20 rounds:
    // counter, key, result
    struct Vector {
        unsigned int c0, c1, k0, k1, r0, r1;
    };
    const Vector vectors[] = {
        {0x00000000u, 0x00000000u, 0x00000000u, 0x00000000u, 0x6b200159u,
         0x99ba4efeu},
        {0xffffffffu, 0xffffffffu, 0xffffffffu, 0xffffffffu, 0x1cb996fcu,
         0xbb002be7u},
        {0x243f6a88u, 0x85a308d3u, 0x13198a2eu, 0x03707344u, 0xc4923a9cu,
         0x483df7a0u},
    };
    for (const Vector &v : vectors) {
        const kernels::ipg_u32_pair result =
            kernels::ipg_threefry2x32(v.c0, v.c1, v.k0, v.k1);
        CHECK(result.a == v.r0);
        CHECK(result.b == v.r1);
    }
}

TEST_CASE("the device noise is Gaussian, independent and reproducible") {
    std::vector<double> hostDoubleNoise;
    for (DeviceCase &device : deviceCases()) {
        CAPTURE(device.label);
        const int N = 64;
        const std::size_t sites = static_cast<std::size_t>(N) * N;
        Parameters param;
        makeDeviceTestParam(param, N);
        Group group;
        Lattice lattice(&param, N);
        for (int pos = 0; pos < N * N; pos++) lattice.U[pos] = Matrix(1.);
        Random random;
        random.init_genrand64(7ULL);
        JIMWLK jimwlk(param, &group, &lattice, &random);
        std::unique_ptr<ipgdev::JimwlkDevice> evolution =
            jimwlk.makeDevice(*device.backend);
        evolution->upload(lattice.U);

        const unsigned int key0 = 0x9E3779B9u;
        const unsigned int key1 = 0x7F4A7C15u;
        std::vector<double> noise, again, otherStep, otherKey;
        evolution->stepWithDeviceNoise(key0, key1, 3);
        evolution->downloadNoise(noise);
        evolution->stepWithDeviceNoise(key0, key1, 3);
        evolution->downloadNoise(again);
        evolution->stepWithDeviceNoise(key0, key1, 4);
        evolution->downloadNoise(otherStep);
        evolution->stepWithDeviceNoise(key0 + 1, key1, 3);
        evolution->downloadNoise(otherKey);
        REQUIRE(noise.size() == 16 * sites);

        // the same (key, step) gives the same numbers, another step or key
        // other ones
        CHECK(noise == again);
        const std::size_t n = noise.size();
        double sum = 0., sum2 = 0., sum3 = 0., sum4 = 0.;
        double withStep = 0., withKey = 0., withNext = 0., largest = 0.;
        std::size_t beyondTwo = 0;
        for (std::size_t k = 0; k < n; k++) {
            const double g = noise[k];
            sum += g;
            sum2 += g * g;
            sum3 += g * g * g;
            sum4 += g * g * g * g;
            withStep += g * otherStep[k];
            withKey += g * otherKey[k];
            // the neighbor in the buffer: the partner of a Box-Muller pair,
            // the next color or the next site
            withNext += g * noise[(k + 1) % n];
            largest = std::max(largest, std::abs(g));
            if (std::abs(g) > 2.) beyondTwo++;
        }
        // n = 65536 numbers: the errors of the moments are 1/sqrt(n) = 0.004
        // times sqrt(1), sqrt(2), sqrt(15), sqrt(96); 5 standard deviations
        const double sigma = 1. / std::sqrt(static_cast<double>(n));
        CHECK(std::abs(sum / n) < 5. * sigma);
        CHECK(std::abs(sum2 / n - 1.) < 5. * std::sqrt(2.) * sigma);
        CHECK(std::abs(sum3 / n) < 5. * std::sqrt(15.) * sigma);
        CHECK(std::abs(sum4 / n - 3.) < 5. * std::sqrt(96.) * sigma);
        CHECK(std::abs(withStep / n) < 5. * sigma);
        CHECK(std::abs(withKey / n) < 5. * sigma);
        CHECK(std::abs(withNext / n) < 5. * sigma);
        // P(|g| > 2) = 0.0455003
        const double fraction = static_cast<double>(beyondTwo) / n;
        CHECK(
            std::abs(fraction - 0.0455003)
            < 5. * std::sqrt(0.0455 * 0.9545) * sigma);
        CHECK(largest > 3.5);
        CHECK(largest < 7.);

        // every backend draws the same numbers (up to its rounding): those
        // of the host in double precision
        if (hostDoubleNoise.empty()) {
            std::unique_ptr<ipgdev::Backend> host =
                ipgdev::createHostBackendDouble();
            std::unique_ptr<ipgdev::JimwlkDevice> hostEvolution =
                jimwlk.makeDevice(*host);
            hostEvolution->upload(lattice.U);
            hostEvolution->stepWithDeviceNoise(key0, key1, 3);
            hostEvolution->downloadNoise(hostDoubleNoise);
        }
        {
            double difference = 0.;
            for (std::size_t k = 0; k < n; k++) {
                difference = std::max(
                    difference, std::abs(noise[k] - hostDoubleNoise[k]));
            }
            report(
                device.label + " device noise against host/double", difference);
            // single precision: the rounding of ln and of the angle
            CHECK(difference < device.tolerance);
        }
    }
}

TEST_CASE(
    "JIMWLK with noise drawn on the device evolves like JIMWLK with the "
    "noise of the host") {
    // The two kinds of noise are different random numbers, so the runs
    // agree only on average. The Wilson lines start as a smooth field,
    // whose neighbors are almost equal: Re tr(V^dagger(x) V(x + e_x)) / 3
    // is close to 1. The evolution decorrelates neighbors and lowers the
    // site average of this dipole. The kernel is long ranged, so the sites
    // of one run are not independent samples: the test evolves several
    // realizations with each kind of noise and compares the means of the
    // dipole with the error that the spread of the realizations gives.
    for (DeviceCase &device : deviceCases()) {
        CAPTURE(device.label);
        const int N = 24;
        const int realizations = 16;
        Parameters param;
        makeDeviceTestParam(param, N);
        param.jimwlk.alphaS = 0.3;
        param.jimwlk.mass = 0.;
        param.jimwlk.Ds = 0.01;
        Group group;

        auto dipole = [&](const Lattice &lattice) {
            double sum = 0.;
            for (int ix = 0; ix < N; ix++) {
                for (int iy = 0; iy < N; iy++) {
                    const Matrix &V = lattice.U[ix * N + iy];
                    const Matrix &W = lattice.U[((ix + 1) % N) * N + iy];
                    sum += Matrix::prodAconjB(V, W).trace().real() / 3.;
                }
            }
            return sum / (static_cast<double>(N) * N);
        };
        double initial = 0.;
        // the mean of the dipole after the evolution and the error of this
        // mean, over the realizations
        auto evolve = [&](bool deviceNoise, double &mean, double &error) {
            double sum = 0., sum2 = 0.;
            for (int realization = 0; realization < realizations;
                 realization++) {
                Lattice lattice(&param, N);
                for (int ix = 0; ix < N; ix++) {
                    for (int iy = 0; iy < N; iy++) {
                        std::vector<double> q(8, 0.);
                        q[0] = 0.8 * std::cos(2. * M_PI * ix / N);
                        q[3] = 0.8 * std::sin(2. * M_PI * iy / N);
                        q[6] = 0.5 * std::cos(2. * M_PI * (ix + iy) / N);
                        lattice.U[ix * N + iy] = Matrix::fromAlgebraExponent(q);
                    }
                }
                initial = dipole(lattice);
                Random random;
                random.init_genrand64(1000ULL + realization);
                JIMWLK jimwlk(param, &group, &lattice, &random);
                std::unique_ptr<ipgdev::JimwlkDevice> evolution =
                    jimwlk.makeDevice(*device.backend);
                evolution->upload(lattice.U);
                unsigned int key0 = 0, key1 = 0;
                jimwlk.drawDeviceNoiseKey(key0, key1);
                for (unsigned int step = 0; step < 20; step++) {
                    if (deviceNoise) {
                        evolution->stepWithDeviceNoise(key0, key1, step);
                    } else {
                        jimwlk.evolutionStepOnDevice(*evolution);
                    }
                }
                evolution->download(lattice.U);
                const double value = dipole(lattice);
                sum += value;
                sum2 += value * value;
            }
            mean = sum / realizations;
            error = std::sqrt(
                (sum2 / realizations - mean * mean) / (realizations - 1));
        };
        double hostMean, hostError, deviceMean, deviceError;
        evolve(false, hostMean, hostError);
        evolve(true, deviceMean, deviceError);
        report(device.label + " neighbor dipole before", initial);
        report(device.label + " neighbor dipole, host noise", hostMean);
        report(device.label + "   its error", hostError);
        report(device.label + " neighbor dipole, device noise", deviceMean);
        report(device.label + "   its error", deviceError);
        // the evolution lowered the dipole well beyond its error ...
        CHECK(initial > 0.98);
        CHECK(initial - hostMean > 10. * hostError);
        CHECK(initial - deviceMean > 10. * deviceError);
        // ... by the same amount for both kinds of noise
        CHECK(
            std::abs(hostMean - deviceMean)
            < 4.
                  * std::sqrt(
                      hostError * hostError + deviceError * deviceError));
    }
}

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <gsl/gsl_sf_gamma.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "Random.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#include "doctest.h"

TEST_CASE("Random: gaussBulk reproduces the exact gauss() stream") {
    const std::size_t n = 2001;  // odd, to exercise the leftover-pair path
    const unsigned long long seed = 12345ULL;

    Random scalarRng;
    scalarRng.init_genrand64(seed);
    std::vector<double> scalarStream(n);
    for (std::size_t i = 0; i < n; ++i) scalarStream[i] = scalarRng.gauss();
    // gauss() uses the Box-Muller transform, which produces values in
    // pairs and caches the second one for the next call. n is odd, so at
    // this point scalarRng has one cached value left over; draw it now so
    // the comparison below also covers gaussBulk() leaving the same value
    // cached for its own next scalar draw (an implementation that returns
    // the right n values but clears or overwrites that cache would
    // otherwise still pass this test).
    const double scalarNext = scalarRng.gauss();

    Random bulkRng;
    bulkRng.init_genrand64(seed);
    std::vector<double> bulkStream(n);
    std::vector<double> scratch;
    bulkRng.gaussBulk(bulkStream.data(), n, scratch);

    for (std::size_t i = 0; i < n; ++i) {
        CHECK(bulkStream[i] == scalarStream[i]);
    }
    CHECK(bulkRng.gauss() == scalarNext);
}

TEST_CASE(
    "Random: gaussBulk reproduces the scalar stream for any number of "
    "threads and sizes, mixed with scalar draws (issue #121)") {
    // the acceptance scan runs in parallel, with a prefix sum over the
    // threads; sizes below the number of threads leave some of them idle
    const std::vector<std::size_t> sizes = {1, 2, 3, 7, 1000, 524289, 4};
    const unsigned long long seed = 777ULL;
    Random scalarRng;
    scalarRng.init_genrand64(seed);
    std::vector<double> expected;
    for (const std::size_t n : sizes) {
        for (std::size_t i = 0; i < n; ++i) {
            expected.push_back(scalarRng.gauss());
        }
        expected.push_back(scalarRng.genrand64_real2());
    }
#ifdef _OPENMP
    const int maxThreads = omp_get_max_threads();
    for (const int threads : {1, 3, 8}) {
        omp_set_num_threads(threads);
#else
    for (const int threads : {1}) {
#endif
        CAPTURE(threads);
        Random bulkRng;
        bulkRng.init_genrand64(seed);
        std::vector<double> scratch, values, stream;
        for (const std::size_t n : sizes) {
            values.resize(n);
            bulkRng.gaussBulk(values.data(), n, scratch);
            stream.insert(stream.end(), values.begin(), values.end());
            stream.push_back(bulkRng.genrand64_real2());
        }
        REQUIRE(stream.size() == expected.size());
        std::size_t mismatches = 0;
        for (std::size_t i = 0; i < stream.size(); ++i) {
            if (stream[i] != expected[i]) mismatches++;
        }
        CHECK(mismatches == 0);
    }
#ifdef _OPENMP
    omp_set_num_threads(maxThreads);
#endif
}

TEST_CASE("Random: gaussBulk output has ~zero mean and ~unit variance") {
    const std::size_t n = 20000;
    Random rng;
    rng.init_genrand64(42ULL);

    std::vector<double> out(n);
    std::vector<double> scratch;
    rng.gaussBulk(out.data(), n, scratch);

    double sum = 0.0;
    for (double v : out) sum += v;
    const double mean = sum / static_cast<double>(n);

    double sqSum = 0.0;
    for (double v : out) sqSum += (v - mean) * (v - mean);
    const double variance = sqSum / static_cast<double>(n);

    // Generous tolerances (a few standard errors) so this isn't flaky for a
    // fixed seed while still catching a badly broken generator.
    CHECK(std::abs(mean) < 0.05);
    CHECK(std::abs(variance - 1.0) < 0.1);
}

TEST_CASE("Random: init_genrand64 makes genrand64_int64/int63 deterministic") {
    Random a, b;
    a.init_genrand64(777ULL);
    b.init_genrand64(777ULL);
    for (int i = 0; i < 100; ++i) {
        CHECK(a.genrand64_int64() == b.genrand64_int64());
    }

    a.init_genrand64(777ULL);
    b.init_genrand64(777ULL);
    for (int i = 0; i < 100; ++i) {
        CHECK(a.genrand64_int63() == b.genrand64_int63());
    }
}

TEST_CASE(
    "Random: genrand64_int63 is genrand64_int64 shifted right by one (bottom "
    "bit dropped, sign bit clear)") {
    Random rng;
    rng.init_genrand64(1ULL);
    Random rng2;
    rng2.init_genrand64(1ULL);
    for (int i = 0; i < 50; ++i) {
        const unsigned long long full = rng.genrand64_int64();
        const long long half = rng2.genrand64_int63();
        CHECK(half >= 0);
        CHECK(static_cast<unsigned long long>(half) == (full >> 1));
    }
}

TEST_CASE(
    "Random: genrand64_real1/real2/real3 stay within their documented ranges") {
    Random rng;
    rng.init_genrand64(2024ULL);
    for (int i = 0; i < 5000; ++i) {
        const double r1 = rng.genrand64_real1();  // [0, 1]
        CHECK(r1 >= 0.0);
        CHECK(r1 <= 1.0);
        const double r2 = rng.genrand64_real2();  // [0, 1)
        CHECK(r2 >= 0.0);
        CHECK(r2 < 1.0);
        const double r3 = rng.genrand64_real3();  // (0, 1)
        CHECK(r3 > 0.0);
        CHECK(r3 < 1.0);
    }
}

TEST_CASE(
    "Random: gslRandomInit makes poisson() deterministic and ~correct mean") {
    Random a, b;
    a.gslRandomInit(2025ULL);
    b.gslRandomInit(2025ULL);

    const double mean = 5.0;
    const int n = 5000;
    long sum = 0;
    for (int i = 0; i < n; ++i) {
        const int drawA = a.poisson(mean);
        const int drawB = b.poisson(mean);
        CHECK(drawA == drawB);  // same seed -> same stream
        CHECK(drawA >= 0);
        sum += drawA;
    }
    const double empiricalMean = static_cast<double>(sum) / n;
    CHECK(std::abs(empiricalMean - mean) < 0.3);  // generous, fixed-seed check
}

namespace {
/// Kolmogorov-Smirnov statistic sqrt(n) D of \p samples against \p cdf,
/// below 1.95 with probability 0.999 if they follow it. With \p stride
/// > 1 the CDF is only evaluated at every stride-th sorted sample, which
/// underestimates D by at most about 2 stride/n.
template <typename Cdf>
double kolmogorovSmirnov(
    std::vector<double> samples, Cdf cdf, std::size_t stride = 1) {
    std::sort(samples.begin(), samples.end());
    const double n = static_cast<double>(samples.size());
    double d = 0.;
    for (std::size_t i = 0; i < samples.size(); i += stride) {
        const double F = cdf(samples[i]);
        d = std::max({d, F - i / n, (i + 1) / n - F});
    }
    return std::sqrt(n) * d;
}
}  // namespace

TEST_CASE("Random: gammaVariate follows the Gamma distribution") {
    Random rng;
    rng.init_genrand64(101ULL);
    const int n = 200000;
    for (const double shape : {0.5, 1., 4. / 3., 3., 11.}) {
        CAPTURE(shape);
        std::vector<double> samples(n);
        double sum = 0., sum2 = 0.;
        bool positive = true;
        for (double &x : samples) {
            x = rng.gammaVariate(shape);
            if (!(x > 0.)) positive = false;
            sum += x;
            sum2 += x * x;
        }
        CHECK(positive);
        // mean and variance are both shape; the fourth central moment
        // is 3 shape^2 + 6 shape
        const double mean = sum / n;
        const double var = sum2 / n - mean * mean;
        CHECK(std::abs(mean - shape) < 5. * std::sqrt(shape / n));
        CHECK(
            std::abs(var - shape)
            < 5. * std::sqrt((2. * shape * shape + 6. * shape) / n));
        CHECK(
            kolmogorovSmirnov(
                samples,
                [shape](double x) { return gsl_sf_gamma_inc_P(shape, x); })
            < 1.95);
    }
}

TEST_CASE(
    "Random: sampleGammaInc draws x with the density omega Q(1/omega, x)") {
    // With s = 1/omega the CDF is
    // int_0^x Q(s, t) dt / s = x Q(s, x)/s + P(s + 1, x), the mean
    // (s + 1)/2 and the second moment (s + 1)(s + 2)/3. A table of 4000
    // bins whose left edges were returned (issue #104) gave
    // sqrt(n) D = 5-15 here.
    Random rng;
    rng.init_genrand64(99ULL);
    const int n = 1000000;
    for (const double omega : {0.3, 0.5, 2., 3.}) {
        CAPTURE(omega);
        const double s = 1. / omega;
        std::vector<double> samples(n);
        double sum = 0.;
        bool positive = true;
        for (double &x : samples) {
            x = rng.sampleGammaInc(omega);
            if (!(x > 0.)) positive = false;
            sum += x;
        }
        CHECK(positive);
        const double mean = (s + 1.) / 2.;
        const double var = (s + 1.) * (s + 2.) / 3. - mean * mean;
        CHECK(std::abs(sum / n - mean) < 5. * std::sqrt(var / n));
        CHECK(
            kolmogorovSmirnov(
                samples,
                [s](double x) {
                    return x * gsl_sf_gamma_inc_Q(s, x) / s
                           + gsl_sf_gamma_inc_P(s + 1., x);
                },
                100)
            < 1.95);
    }
}

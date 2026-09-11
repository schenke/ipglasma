#include <cmath>
#include <vector>

#include "Random.h"
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

TEST_CASE("Random: setGammaIncCDF/sampleGammaInc stay within [0, xmax]") {
    Random rng;
    rng.init_genrand64(99ULL);

    const double omega = 1.0;
    rng.setGammaIncCDF(omega);
    const double xmax = std::max(5.0, 5.0 / omega);

    bool sawNonzero = false;
    for (int i = 0; i < 2000; ++i) {
        const double x = rng.sampleGammaInc();
        CHECK(x >= 0.0);
        CHECK(x <= xmax);
        if (x > 0.0) sawNonzero = true;
    }
    CHECK(sawNonzero);
}

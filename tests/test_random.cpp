#include "doctest.h"

#include <cmath>
#include <vector>

#include "Random.h"

TEST_CASE("Random: gaussBulk reproduces the exact gauss() stream") {
    const std::size_t n = 2001;  // odd, to exercise the leftover-pair path
    const unsigned long long seed = 12345ULL;

    Random scalarRng;
    scalarRng.init_genrand64(seed);
    std::vector<double> scalarStream(n);
    for (std::size_t i = 0; i < n; ++i) scalarStream[i] = scalarRng.gauss();

    Random bulkRng;
    bulkRng.init_genrand64(seed);
    std::vector<double> bulkStream(n);
    std::vector<double> scratch;
    bulkRng.gaussBulk(bulkStream.data(), n, scratch);

    for (std::size_t i = 0; i < n; ++i) {
        CHECK(bulkStream[i] == scalarStream[i]);
    }
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

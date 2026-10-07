#include <cmath>
#include <vector>

#include "Group.h"
#include "JIMWLK.h"
#include "Lattice.h"
#include "Parameters.h"
#include "Random.h"
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
    "computed reference value, and is sensitive to nFlavors/c_jimwlk") {
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

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
    param.colorCharge.RapidityA = 0.0;
    param.colorCharge.RapidityB = 0.0;
    param.jimwlk.mu0_jimwlk = 0.2;
    param.jimwlk.Lambda_QCD_jimwlk = 0.2;
    param.jimwlk.c_jimwlk = 0.2;
    param.coupling.nFlavors = 3;
    param.jimwlk.m_jimwlk = 0.0;       // skips the Bessel mass-regulator branch
    param.jimwlk.alphas_jimwlk = 0.3;  // fixed coupling, skips running-coupling
    param.jimwlk.Ds_jimwlk = 0.001;
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
    // Fixed coupling (jimwlk.alphas_jimwlk > 0): steps_1 = as*log(x0/x_proj) /
    // (pi^2*ds) + 0.5. These values give steps_1 = steps_2 = 1.
    param.jimwlk.jimwlk_ic_x = 0.01;
    param.jimwlk.x_projectile_jimwlk = 0.008;
    param.jimwlk.x_target_jimwlk = 0.008;

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
    param.jimwlk.alphas_jimwlk = 0.3;

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
    param.jimwlk.mu0_jimwlk = 0.28;
    param.jimwlk.Lambda_QCD_jimwlk = 0.04;
    param.jimwlk.alphas_jimwlk = 0.0;  // forces the running-coupling branch

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

    param.jimwlk.c_jimwlk = 0.25;
    const double alphasC025 = jimwlk.getAlphas(x, y);
    CHECK(alphasC025 == doctest::Approx(0.3453365551834896));
}

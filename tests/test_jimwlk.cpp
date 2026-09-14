#include <cmath>

#include "Group.h"
#include "JIMWLK.h"
#include "Lattice.h"
#include "Parameters.h"
#include "Random.h"
#include "doctest.h"

namespace {
void makeJimwlkTestParam(Parameters &param, int size) {
    param.setSize(size);
    param.setL(static_cast<double>(size));
    param.setMPIRank(0);
    param.setEventId(0);
    param.setSeed(0);
    param.setMPISize(1);
    param.setRapidityA(0.0);
    param.setRapidityB(0.0);
    param.setMu0_jimwlk(0.2);
    param.setLambdaQCD_jimwlk(0.2);
    param.setm_jimwlk(0.0);       // skips the Bessel mass-regulator branch
    param.setJimwlk_alphas(0.3);  // fixed coupling, skips running-coupling
    param.setDs_jimwlk(0.001);
}
}  // namespace

TEST_CASE(
    "JIMWLK::evolutionStep runs with simpleLangevin=0 without crashing "
    "(regression test: VxsiVx_/VxsiVy_ used to be allocated only when "
    "simpleLangevin was true, but evolutionStep needs them regardless)") {
    const int N = 8;
    Parameters param;
    makeJimwlkTestParam(param, N);
    param.setSimpleLangevin(0);

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

TEST_CASE("JIMWLK::evolutionStep runs with simpleLangevin=1 without crashing") {
    const int N = 8;
    Parameters param;
    makeJimwlkTestParam(param, N);
    param.setSimpleLangevin(1);

    Group group;
    Random random;
    random.init_genrand64(42ULL);
    Lattice lat(&param, N);

    JIMWLK jimwlk(param, &group, &lat, &random);
    jimwlk.evolutionStep(NucleusRole::Projectile);

    for (int pos = 0; pos < N * N; ++pos) {
        for (int k = 0; k < 9; ++k) {
            CHECK(std::isfinite(lat.U[pos].get(k).real()));
            CHECK(std::isfinite(lat.U[pos].get(k).imag()));
        }
    }
}

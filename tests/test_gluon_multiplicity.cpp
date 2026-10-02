#include <cstdio>
#include <fstream>
#include <random>
#include <string>

#include "GluonMultiplicity.h"
#include "Group.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "doctest.h"
#include "test_helpers.h"

namespace {
const int N = 8;
const int it = 4;

// a = 1 fm and maxTime = it * a * dtau, so that it is the final step
void makeMultiplicityTestParam(Parameters &param) {
    makeLatticeParam(param, N);
    param.lattice.L = N;
    param.coupling.g = 1.;
    param.run.dtau = 0.1;
    param.evolution.maxTime = 0.4;
    param.evolution.inverseQsForMaxTime = false;
    param.coupling.runningCoupling = false;
    param.output.writeOutputs = 0;
}

void setVacuum(Lattice &lat) {
    for (int pos = 0; pos < N * N; pos++) {
        lat.Ux[pos] = Matrix(1.);
        lat.Uy[pos] = Matrix(1.);
        lat.U[pos] = Matrix(0.);
        lat.U2[pos] = Matrix(0.);
        lat.Ux2[pos] = Matrix(0.);
        lat.Uy2[pos] = Matrix(0.);
    }
}

void removeOutput() {
    std::remove("NpartdNdy-t0.4-0.dat");
    std::remove("gluonMultiplicity0.json");
}

bool fileExists(const std::string &name) { return std::ifstream(name).good(); }
}  // namespace

TEST_CASE("GluonMultiplicity::compute reports no collision in the vacuum") {
    Parameters param;
    makeMultiplicityTestParam(param);
    Lattice lat(&param, N);
    setVacuum(lat);
    Group group;
    int nn[2] = {N, N};
    GluonMultiplicity multiplicity(nn);
    removeOutput();
    param.event.success = 1;

    CHECK(multiplicity.compute(&lat, &group, &param, it) == 0);
    CHECK_FALSE(fileExists("NpartdNdy-t0.4-0.dat"));
    removeOutput();
}

TEST_CASE(
    "GluonMultiplicity::compute measures a positive multiplicity for "
    "electric fields and writes the event files") {
    Parameters param;
    makeMultiplicityTestParam(param);
    Lattice lat(&param, N);
    setVacuum(lat);
    Group group;
    std::mt19937 rng(5);
    std::normal_distribution<double> gauss(0., 0.3);
    for (int pos = 0; pos < N * N; pos++) {
        for (Matrix *field : {&lat.U[pos], &lat.U2[pos]}) {
            for (int a = 0; a < 8; a++) *field += gauss(rng) * group.getT(a);
        }
    }
    int nn[2] = {N, N};
    GluonMultiplicity multiplicity(nn);
    removeOutput();
    param.event.success = 0;

    REQUIRE(multiplicity.compute(&lat, &group, &param, it) == 1);
    CHECK(param.event.success == 1);
    std::ifstream in("NpartdNdy-t0.4-0.dat");
    REQUIRE(in.good());
    double Npart, dNdy, Tpp, b, dEdy;
    in >> Npart >> dNdy >> Tpp >> b >> dEdy;
    CHECK(dNdy > 0.);
    CHECK(dEdy > 0.);
    CHECK(fileExists("gluonMultiplicity0.json"));
    removeOutput();
}

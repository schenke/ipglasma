#include <algorithm>
#include <cstdio>
#include <fstream>
#include <random>
#include <sstream>
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
    param.output.writeHadronSpectrum = false;
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
    // the 14 columns documented in OUTPUT.md, with placeholders in 7-9
    std::string seed, na1, na2, na3;
    in >> seed >> na1 >> na2 >> na3;
    CHECK(na1 == "N/A");
    CHECK(na2 == "N/A");
    CHECK(na3 == "N/A");
    int rest = 0;
    for (std::string token; in >> token;) rest++;
    CHECK(rest == 5);

    // the JSON file has the documented keys and four arrays of 100 bins
    std::ifstream json("gluonMultiplicity0.json");
    REQUIRE(json.good());
    std::stringstream text;
    text << json.rdbuf();
    const std::string content = text.str();
    for (const char *key :
         {"\"format\": \"ipglasma-gluon-target\"", "\"dN\":", "\"dE_GeV\":",
          "\"tau_fm\":", "\"rapidity_variable\": \"y\""}) {
        CAPTURE(key);
        CHECK(content.find(key) != std::string::npos);
    }
    for (const char *array :
         {"\"kt_GeV\": [", "\"dN_d2k_GeV_minus2\": [",
          "\"dE_d2k_GeV_minus1\": [", "\"lattice_bin_counts\": ["}) {
        CAPTURE(array);
        const std::size_t start = content.find(array);
        REQUIRE(start != std::string::npos);
        const std::size_t end = content.find(']', start);
        const std::string values = content.substr(start, end - start);
        CHECK(std::count(values.begin(), values.end(), ',') == 99);
    }
    removeOutput();
}

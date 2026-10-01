#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "Cell.h"
#include "Evolution.h"
#include "Lattice.h"
#include "Parameters.h"
#include "doctest.h"

namespace {
// Matches test_lattice.cpp's makeLatticeParam: Parameters holds a
// PrettyOstream member, which holds a non-copyable/non-movable
// std::ostringstream, so this takes an out-parameter rather than
// returning by value.
void makeEvolutionTestParam(Parameters &param, int size) {
    param.setSize(size);
    param.setL(static_cast<double>(size));  // a = 1 fm
    param.setMPIRank(0);
    param.setEventId(0);
    param.setSeed(0);
    param.setMPISize(1);
    param.setRapidityA(0.0);
    param.setRapidityB(0.0);
    param.setg(1.0);
    param.setdtau(0.1);
    param.setRunningCoupling(0);  // gfactor == 1 everywhere; see below
}

}  // namespace

TEST_CASE(
    "Evolution::eccentricity(doAniso=1): the unrotated Txx-Tyy/Txx+Tyy ratio "
    "(num2/den2) is the same at every sampled angle") {
    const int N = 8;
    Parameters param;
    makeEvolutionTestParam(param, N);
    param.setQsmuRatio(0.643);
    param.setQsmuRatioB(0.643);
    Lattice lat(&param, N);

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            const int pos = ix * N + iy;
            const double x = ix - N / 2.0 + 0.3;
            const double y = iy - N / 2.0 - 0.2;
            const double r2 = x * x + y * y + 1.0;
            lat.cells[pos]->setEpsilon(1.0 / r2);
            lat.cells[pos]->setutau(1.0);
            lat.cells[pos]->setTxx(1.5 / r2 + 0.1 * x);
            lat.cells[pos]->setTyy(0.9 / r2 - 0.05 * y);
            lat.cells[pos]->setTxy(0.2 * x * y / r2);
            lat.cells[pos]->setux(0.01 * x);
            lat.cells[pos]->setuy(0.02 * y - 0.005 * x);
            lat.cells[pos]->setg2mu2A(0.5 / r2);
            lat.cells[pos]->setg2mu2B(0.4 / r2);
        }
    }

    // Independently compute the (Psi-independent) unrotated ratio the same
    // way computeRotatedAnisotropy does, as this file's own ground truth.
    double num2 = 0., den2 = 0.;
    for (int pos = 0; pos < N * N; ++pos) {
        num2 += lat.cells[pos]->getTxx() - lat.cells[pos]->getTyy();
        den2 += lat.cells[pos]->getTxx() + lat.cells[pos]->getTyy();
    }
    const double expectedRatio = num2 / den2;

    int nn[2] = {N, N};
    Evolution evo(nn);
    std::remove("anisotropy0.dat");
    std::remove("eccentricities0.dat");
    evo.eccentricity(&lat, &param, /*it=*/5, /*cutoff=*/0.0, /*doAniso=*/1);

    std::ifstream in("anisotropy0.dat");
    REQUIRE(in.good());
    std::string line;
    std::getline(in, line);  // "Psi2=..." header
    std::getline(in, line);  // "PsiU=..." header

    int dataLines = 0;
    double tau, ratio, ratio2;
    // The line's last field is "angle=<value>" with no space before the
    // number, i.e. one whitespace-delimited token; its value isn't needed
    // here (only that ratio2, the Psi-independent num2/den2, holds).
    std::string angleToken;
    while (in >> tau >> ratio >> ratio2 >> angleToken) {
        CHECK(ratio2 == doctest::Approx(expectedRatio));
        ++dataLines;
    }
    CHECK(dataLines == 10);  // Psi = PsiU + k*pi/8 for k = 0..9

    std::remove("anisotropy0.dat");
    std::remove("eccentricities0.dat");
}

TEST_CASE(
    "Evolution::finalFlowMeasurement fills epsilon and u^tau for "
    "eccentricity() also with writeEpsilonUHydro off") {
    // Every cell holds the T^{mu nu} of an ideal fluid with energy density
    // e, pressure p = e/3 and velocity v along x, so the flow solve must
    // return epsilon = e and u^tau = gamma.
    const int N = 4;
    const int it = 5;
    const double e = 1.0;
    const double p = e / 3.;
    const double v = 0.2;
    const double gamma = 1. / std::sqrt(1. - v * v);

    struct Case {
        int writeEpsilonUHydro;
        int computeGluonMultiplicity;
        bool expectSolve;
    };
    for (const Case &c :
         std::vector<Case> {{0, 1, true}, {1, 0, true}, {0, 0, false}}) {
        CAPTURE(c.writeEpsilonUHydro);
        CAPTURE(c.computeGluonMultiplicity);
        Parameters param;
        makeEvolutionTestParam(param, N);
        param.setWriteOutputs(0);  // no output files
        param.setWriteEpsilonUHydro(c.writeEpsilonUHydro);
        param.setComputeGluonMultiplicity(c.computeGluonMultiplicity);
        Lattice lat(&param, N);

        const double a = param.getL() / N;
        const double tau = it * param.getdtau() * a;
        for (int pos = 0; pos < N * N; ++pos) {
            Cell *cell = lat.cells[pos];
            cell->setTtautau((e + p) * gamma * gamma - p);
            cell->setTtaux((e + p) * gamma * gamma * v);
            cell->setTxx((e + p) * gamma * gamma * v * v + p);
            cell->setTyy(p);
            cell->setTetaeta(p / (tau * tau));
        }

        int nn[2] = {N, N};
        Evolution evolution(nn);
        evolution.finalFlowMeasurement(&lat, &param, it);

        for (int pos = 0; pos < N * N; ++pos) {
            CAPTURE(pos);
            if (c.expectSolve) {
                CHECK(lat.cells[pos]->getEpsilon() == doctest::Approx(e));
                CHECK(lat.cells[pos]->getutau() == doctest::Approx(gamma));
            } else {
                // the expensive solve is skipped when nothing needs it
                CHECK(lat.cells[pos]->getutau() == 0.);
            }
        }
    }
}

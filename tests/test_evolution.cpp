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
    param.lattice.size = size;
    param.lattice.L = static_cast<double>(size);  // a = 1 fm
    param.run.MPIRank = 0;
    param.event.eventId = 0;
    param.random.seed = 0;
    param.run.MPISize = 1;
    param.colorCharge.rapidityA = 0.0;
    param.colorCharge.rapidityB = 0.0;
    param.coupling.g = 1.0;
    param.run.dtau = 0.1;
    param.coupling.runningCoupling =
        false;  // gfactor == 1 everywhere; see below
}

}  // namespace

TEST_CASE(
    "Evolution::finalFlowMeasurement fills epsilon and u^tau for "
    "Eccentricity::compute() also with writeEpsilonUHydro off") {
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
        param.output.writeOutputs = 0;  // no output files
        param.output.writeEpsilonUHydro = c.writeEpsilonUHydro;
        param.output.computeGluonMultiplicity = c.computeGluonMultiplicity;
        Lattice lat(&param, N);

        const double a = param.lattice.L / N;
        const double tau = it * param.run.dtau * a;
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

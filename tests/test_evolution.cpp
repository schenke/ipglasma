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
    param.colorCharge.projectileX = 0.01;
    param.colorCharge.targetX = 0.01;
    param.coupling.g = 1.0;
    param.run.dtau = 0.1;
    param.coupling.runningCoupling =
        false;  // gfactor == 1 everywhere; see below
}

}  // namespace

TEST_CASE(
    "Evolution::finalFlowMeasurement fills epsilon and u^tau for "
    "Eccentricity::compute() also without hydro output") {
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
        bool writeHydro;
        bool computeEccentricities;
        bool expectSolve;
    };
    for (const Case &c : std::vector<Case> {
             {false, true, true}, {true, false, true}, {false, false, false}}) {
        CAPTURE(c.writeHydro);
        CAPTURE(c.computeEccentricities);
        Parameters param;
        makeEvolutionTestParam(param, N);
        param.output.writeHydro = c.writeHydro;
        param.output.writeJazma = false;
        param.output.writeTmunu = false;
        param.output.computeEccentricities = c.computeEccentricities;
        // an empty output grid: the hydro file has only its header
        param.output.sizeOutput = 0;
        param.output.LOutput = param.lattice.L;
        param.output.etaSizeOutput = 0;
        param.output.dEtaOutput = 0.;
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
        std::remove("epsilon-u-Hydro-TauHydro-0.dat");
    }
}

TEST_CASE(
    "Evolution::outputSteps rounds the output times down to distinct steps "
    "before the final one") {
    const double stepLength = 0.01;  // fm/c
    const int itmax = 40;
    CHECK(Evolution::outputSteps({}, stepLength, itmax).empty());
    CHECK(
        Evolution::outputSteps({0.1, 0.2, 0.3}, stepLength, itmax)
        == std::vector<int> {10, 20, 30});
    // unsorted, duplicates after rounding, 0.005 is step 0 and 0.4 and
    // 0.5 are not before the final step
    CHECK(
        Evolution::outputSteps(
            {0.3, 0.104, 0.1, 0.005, 0.4, 0.5}, stepLength, itmax)
        == std::vector<int> {10, 30});
}

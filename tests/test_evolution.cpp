#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Evolution.h"
#include "Lattice.h"
#include "Parameters.h"
#include "PhysConst.h"
#include "doctest.h"

using PhysConst::hbarc;

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

struct EpsilonPlotRow {
    double x, y, value;
};

std::vector<EpsilonPlotRow> readEpsilonPlot(const std::string &path) {
    std::vector<EpsilonPlotRow> rows;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        EpsilonPlotRow row;
        if (iss >> row.x >> row.y >> row.value) rows.push_back(row);
    }
    return rows;
}
}  // namespace

TEST_CASE(
    "Evolution::writeEpsilonInitialPlot: with runningCoupling off, output is "
    "exactly hbarc*|epsilon| at each cell's (x, y)") {
    const int N = 4;
    Parameters param;
    makeEvolutionTestParam(param, N);
    Lattice lat(&param, N);

    // One positive, one negative (exercises the abs()), rest left at the
    // Lattice default.
    lat.cells[0]->setEpsilon(2.5);
    lat.cells[1]->setEpsilon(-1.25);

    int nn[2] = {N, N};
    Evolution evo(nn);
    std::remove("epsilonInitialPlot0.dat");
    evo.writeEpsilonInitialPlot(&lat, &param);

    std::vector<EpsilonPlotRow> rows =
        readEpsilonPlot("epsilonInitialPlot0.dat");
    std::remove("epsilonInitialPlot0.dat");
    REQUIRE(rows.size() == static_cast<std::size_t>(N * N));

    const double a = param.getL() / N;
    const double L = param.getL();
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            const int pos = ix * N + iy;
            const double expectedX = -L / 2. + a * ix;
            const double expectedY = -L / 2. + a * iy;
            const double expectedValue =
                hbarc * std::abs(lat.cells[pos]->getEpsilon());
            CHECK(rows[pos].x == doctest::Approx(expectedX));
            CHECK(rows[pos].y == doctest::Approx(expectedY));
            CHECK(rows[pos].value == doctest::Approx(expectedValue));
        }
    }
}

TEST_CASE(
    "Evolution::writeEpsilonIntermediatePlot matches writeEpsilonInitialPlot "
    "when runningCoupling is off (both reduce to gfactor=1)") {
    const int N = 4;
    Parameters param;
    makeEvolutionTestParam(param, N);
    Lattice lat(&param, N);
    for (int pos = 0; pos < N * N; ++pos) {
        lat.cells[pos]->setEpsilon(0.1 * (pos + 1));
    }

    int nn[2] = {N, N};
    Evolution evo(nn);
    std::remove("epsilonInitialPlot0.dat");
    std::remove("epsilonIntermediatePlot0.dat");
    evo.writeEpsilonInitialPlot(&lat, &param);
    evo.writeEpsilonIntermediatePlot(&lat, &param);

    std::vector<EpsilonPlotRow> initial =
        readEpsilonPlot("epsilonInitialPlot0.dat");
    std::vector<EpsilonPlotRow> intermediate =
        readEpsilonPlot("epsilonIntermediatePlot0.dat");
    std::remove("epsilonInitialPlot0.dat");
    std::remove("epsilonIntermediatePlot0.dat");

    REQUIRE(initial.size() == intermediate.size());
    for (std::size_t i = 0; i < initial.size(); ++i) {
        CHECK(initial[i].x == doctest::Approx(intermediate[i].x));
        CHECK(initial[i].y == doctest::Approx(intermediate[i].y));
        CHECK(initial[i].value == doctest::Approx(intermediate[i].value));
    }
}

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

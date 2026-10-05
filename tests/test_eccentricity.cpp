#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Cell.h"
#include "Eccentricity.h"
#include "Lattice.h"
#include "Parameters.h"
#include "doctest.h"

namespace {
// Parameters holds a PrettyOstream member, which holds a
// non-copyable/non-movable std::ostringstream, so this takes an
// out-parameter rather than returning by value.
void makeEccentricityTestParam(Parameters &param, int size) {
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
    param.coupling.runningCoupling = false;  // gfactor == 1 everywhere
}
}  // namespace

TEST_CASE(
    "Eccentricity::compute: an elliptic Gaussian energy density has "
    "epsilon_2 = (sigma_y^2 - sigma_x^2) / (sigma_x^2 + sigma_y^2) along y "
    "and no odd harmonics") {
    const int N = 48;
    Parameters param;
    makeEccentricityTestParam(param, N);
    // a = 0.5 fm, the Gaussian centered on x = y = 0; the moments are summed
    // over the asymmetric window ix, iy = 2..N-3, so the tails must be
    // negligible there for the odd harmonics to vanish
    param.lattice.L = 24.;
    const double a = param.lattice.L / N;
    const double sigmaX = 1., sigmaY = 1.5;
    Lattice lat(&param, N);
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            const double x = -param.lattice.L / 2. + a * ix;
            const double y = -param.lattice.L / 2. + a * iy;
            const int pos = lat.positionFromXY(ix, iy);
            lat.cells[pos]->setEpsilon(std::exp(
                -x * x / (2. * sigmaX * sigmaX)
                - y * y / (2. * sigmaY * sigmaY)));
            lat.cells[pos]->setutau(1.0);
        }
    }

    std::remove("eccentricities0.dat");
    const int it = 3;
    Eccentricity::compute(&lat, &param, it, /*cutoff=*/0.0, /*doAniso=*/0);

    std::ifstream in("eccentricities0.dat");
    REQUIRE(in.good());
    double tau, eps[7], psi[7];
    in >> tau;
    for (int n = 1; n <= 6; n++) in >> eps[n] >> psi[n];
    REQUIRE(in.good());
    CHECK(tau == doctest::Approx(it * a * param.run.dtau));
    const double sx2 = sigmaX * sigmaX, sy2 = sigmaY * sigmaY;
    CHECK(eps[2] == doctest::Approx((sy2 - sx2) / (sx2 + sy2)).epsilon(1e-4));
    CHECK(psi[2] == doctest::Approx(M_PI / 2.));
    CHECK(eps[1] < 1e-6);
    CHECK(eps[3] < 1e-6);
    CHECK(eps[5] < 1e-6);
    // the 23 columns documented in OUTPUT.md
    std::ifstream again("eccentricities0.dat");
    std::string line;
    std::getline(again, line);
    std::istringstream columns(line);
    int count = 0;
    for (std::string token; columns >> token;) count++;
    CHECK(count == 23);
    again.close();
    std::remove("eccentricities0.dat");
}

TEST_CASE(
    "Eccentricity::compute(doAniso=1): the unrotated Txx-Tyy/Txx+Tyy ratio "
    "(num2/den2) is the same at every sampled angle") {
    const int N = 8;
    Parameters param;
    makeEccentricityTestParam(param, N);
    param.colorCharge.QsMuRatio = 0.643;
    Lattice lat(&param, N);

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            const int pos = lat.positionFromXY(ix, iy);
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

    std::remove("anisotropy0.dat");
    std::remove("eccentricities0.dat");
    Eccentricity::compute(
        &lat, &param, /*it=*/5, /*cutoff=*/0.0, /*doAniso=*/1);

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
    "Eccentricity::compute leaves out cells below the cutoff, given in "
    "GeV/fm^3") {
    const int N = 16;
    Parameters param;
    makeEccentricityTestParam(param, N);  // a = 1 fm
    Lattice lat(&param, N);
    // epsilon = 1 fm^-4 (0.197 GeV/fm^3) in a 4 x 6 block, 0.1 fm^-4
    // (0.0197 GeV/fm^3) elsewhere
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            const bool inBlock = ix >= 6 && ix < 10 && iy >= 5 && iy < 11;
            const int pos = lat.positionFromXY(ix, iy);
            lat.cells[pos]->setEpsilon(inBlock ? 1. : 0.1);
            lat.cells[pos]->setutau(1.0);
        }
    }
    auto area = [&](double cutoff) {
        std::remove("eccentricities0.dat");
        Eccentricity::compute(&lat, &param, /*it=*/3, cutoff, /*doAniso=*/0);
        std::ifstream in("eccentricities0.dat");
        std::string line;
        std::getline(in, line);
        in.close();
        std::remove("eccentricities0.dat");
        std::istringstream columns(line);
        std::vector<double> values;
        for (double value; columns >> value;) values.push_back(value);
        REQUIRE(values.size() == 23);
        CHECK(values[13] == doctest::Approx(cutoff));  // column 14
        return values[19];  // column 20: area of the cells kept [fm^2]
    };
    CHECK(area(0.) == doctest::Approx(N * N));
    CHECK(area(0.1) == doctest::Approx(4 * 6));  // between the two values
    CHECK(area(0.01) == doctest::Approx(N * N));
}

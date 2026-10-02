#include <algorithm>
#include <cmath>

#include "Lattice.h"
#include "Parameters.h"
#include "PhysConst.h"
#include "RunningCoupling.h"
#include "doctest.h"
#include "test_helpers.h"

TEST_CASE(
    "computeAlphaS: matches an independently computed reference value at "
    "the default 3-flavor, Lambda_QCD=0.2 point") {
    CHECK(
        computeAlphaS(
            /*muZero=*/0.3, /*c=*/0.2, /*lambdaQCD=*/0.2,
            /*nFlavors=*/3, /*scale=*/1.0)
        == doctest::Approx(0.43377345548158003));
    CHECK(
        computeAlphaS(
            /*muZero=*/0.5, /*c=*/0.35, /*lambdaQCD=*/0.2,
            /*nFlavors=*/3, /*scale=*/2.5)
        == doctest::Approx(0.2764060975087722));
}

TEST_CASE("computeAlphaS: increasing nFlavors increases alpha_s") {
    const double alphasNf3 =
        computeAlphaS(0.3, 0.2, /*lambdaQCD=*/0.2, /*nFlavors=*/3, 1.0);
    const double alphasNf4 =
        computeAlphaS(0.3, 0.2, /*lambdaQCD=*/0.2, /*nFlavors=*/4, 1.0);

    CHECK(alphasNf3 == doctest::Approx(0.43377345548158003));
    CHECK(alphasNf4 == doctest::Approx(0.4684753319201063));
    CHECK(alphasNf4 > alphasNf3);
}

TEST_CASE("computeAlphaS: increasing Lambda_QCD increases alpha_s") {
    const double alphasDefault =
        computeAlphaS(0.3, 0.2, /*lambdaQCD=*/0.2, /*nFlavors=*/3, 1.0);
    const double alphasLargerLambda =
        computeAlphaS(0.3, 0.2, /*lambdaQCD=*/0.25, /*nFlavors=*/3, 1.0);

    CHECK(alphasDefault == doctest::Approx(0.43377345548158003));
    CHECK(alphasLargerLambda == doctest::Approx(0.5035953568090804));
    CHECK(alphasLargerLambda > alphasDefault);
}

TEST_CASE("computeRunningCouplingGfactorFromScale: equals g^2/(4*pi*alpha_s)") {
    const double g = 1.3;
    const double muZero = 0.3;
    const double c = 0.2;
    const double lambdaQCD = 0.2;
    const int nFlavors = 3;
    const double scale = 1.0;

    const double alphas = computeAlphaS(muZero, c, lambdaQCD, nFlavors, scale);
    const double expected = g * g / (4. * M_PI * alphas);

    CHECK(
        computeRunningCouplingGfactorFromScale(
            g, muZero, c, lambdaQCD, nFlavors, scale)
        == doctest::Approx(expected));
}

TEST_CASE(
    "computeRunningCouplingGfactor: 1 without running coupling, otherwise "
    "the factor at the selected event-averaged or local Q_s") {
    const int N = 4;
    Parameters param;
    makeLatticeParam(param, N);
    param.coupling.g = 1.3;
    param.coupling.LambdaQCD = 0.2;
    param.coupling.nFlavors = 3;
    param.coupling.runningCouplingQsFactor = 0.5;
    param.colorCharge.QsMuRatio = 0.6;
    param.event.averageQsmin = 0.8;
    param.event.averageQsAvg = 1.1;
    param.event.averageQs = 1.7;
    Lattice lat(&param, N);
    for (int pos = 0; pos < N * N; pos++) {
        lat.cells[pos]->setg2mu2A(0.02 * (pos + 1));
        lat.cells[pos]->setg2mu2B(0.03 * (N * N - pos));
    }
    const double a = 0.25, c = 0.2, muZero = 0.3, g = param.coupling.g;
    auto factorAt = [&](double Qs) {
        return computeRunningCouplingGfactorFromScale(
            g, muZero, c, 0.2, 3, 0.5 * Qs);
    };
    const int interior = lat.positionFromXY(1, 2);

    param.coupling.runningCoupling = false;
    CHECK(
        computeRunningCouplingGfactor(
            &lat, &param, interior, N, a, g, c, muZero)
        == 1.);

    param.coupling.runningCoupling = true;
    param.coupling.runWithLocalQs = false;
    const double averages[3] = {0.8, 1.1, 1.7};  // runWithQs 0 / 1 / 2
    for (int choice = 0; choice < 3; choice++) {
        CAPTURE(choice);
        param.coupling.runWithQs = choice;
        CHECK(
            computeRunningCouplingGfactor(
                &lat, &param, interior, N, a, g, c, muZero)
            == doctest::Approx(factorAt(averages[choice])));
    }

    param.coupling.runWithLocalQs = true;
    const double g2mu2A = lat.cells[interior]->getg2mu2A();
    const double g2mu2B = lat.cells[interior]->getg2mu2B();
    auto localQs = [&](double g2mu2) {
        return std::sqrt(
            g2mu2 * 0.6 * 0.6 / a / a * PhysConst::hbarc * PhysConst::hbarc * g
            * g);
    };
    const double locals[3] = {
        localQs(std::min(g2mu2A, g2mu2B)), localQs((g2mu2A + g2mu2B) / 2.),
        localQs(std::max(g2mu2A, g2mu2B))};
    for (int choice = 0; choice < 3; choice++) {
        CAPTURE(choice);
        param.coupling.runWithQs = choice;
        CHECK(
            computeRunningCouplingGfactor(
                &lat, &param, interior, N, a, g, c, muZero)
            == doctest::Approx(factorAt(locals[choice])));
        // the outermost ring runs with Q_s = 0
        CHECK(
            computeRunningCouplingGfactor(&lat, &param, 0, N, a, g, c, muZero)
            == doctest::Approx(factorAt(0.)));
    }
}

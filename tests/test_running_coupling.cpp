#include "RunningCoupling.h"
#include "doctest.h"

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

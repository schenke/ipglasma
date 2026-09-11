#include "doctest.h"

#include "Glauber.h"

namespace {
Nucleus lookup(const std::string &name) {
    Nucleus nucleus{};
    Glauber glauber;
    glauber.findNucleusData(
        &nucleus, name, /*setWSDeformParams=*/false, 0., 0., 0., 0., 0., 0.,
        /*forceDminFlag=*/false, 0., 0., 0.);
    return nucleus;
}
}  // namespace

TEST_CASE("Glauber::findNucleusData: known species have the right A/Z") {
    Nucleus pb = lookup("Pb");
    CHECK(pb.A == 208);
    CHECK(pb.Z == 82);
    CHECK(pb.anumFunc == 3);  // "3Fermi"

    Nucleus au = lookup("Au");
    CHECK(au.A == 197);
    CHECK(au.Z == 79);

    Nucleus deuteron = lookup("d");
    CHECK(deuteron.A == 2);
    CHECK(deuteron.Z == 1);
    CHECK(deuteron.anumFunc == 8);  // "Hulthen"

    Nucleus proton = lookup("p");
    CHECK(proton.A == 1);
    CHECK(proton.Z == 1);
}

TEST_CASE("Glauber::findNucleusData: setWSDeformParams overrides beta2/R_WS") {
    Nucleus nucleus{};
    Glauber glauber;
    glauber.findNucleusData(
        &nucleus, "Pb", /*setWSDeformParams=*/true, 7.5, 0.6, 0.11, 0.0, 0.0,
        0.0, /*forceDminFlag=*/false, 0., 0., 0.);

    CHECK(nucleus.A == 208);  // species-intrinsic data still set
    CHECK(nucleus.R_WS == doctest::Approx(7.5));
    CHECK(nucleus.a_WS == doctest::Approx(0.6));
    CHECK(nucleus.beta2 == doctest::Approx(0.11));
}

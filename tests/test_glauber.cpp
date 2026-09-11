#include "doctest.h"

#include <string>

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

namespace {
// calcRho() solves for nucleus.rho_WS such that the corresponding anum*()
// normalization integral equals nucleus.A exactly (that is its whole
// purpose): each anum*() is linear in rho_WS, so
//     newRho = A * oldRho / anum*(oldRho)
// makes anum*(newRho) == A by construction, regardless of what placeholder
// rho_WS calcRho started from. This exercises calcRho, the integral()
// machinery, and each density profile's anum*/anum*Int pair together,
// without needing a hand-derived reference number.
double anumForProfile(Glauber &glauber, Nucleus &nucleus) {
    switch (nucleus.anumFunc) {
        case 1: return glauber.anum2HO();
        case 2: return glauber.anum3Gauss(nucleus.R_WS);
        case 3: return glauber.anum3Fermi(nucleus.R_WS);
        case 8: return glauber.anumHulthen();
        default: return -1.0;  // unreachable for the species used below
    }
}
}  // namespace

TEST_CASE("Glauber::calcRho: solved rho_WS reproduces A for every density profile") {
    // One species per density-profile branch calcRho dispatches on.
    for (const char *name : {"Pb", "S", "C", "d"}) {
        Nucleus nucleus{};
        Glauber glauber;
        glauber.findNucleusData(
            &nucleus, name, /*setWSDeformParams=*/false, 0., 0., 0., 0., 0.,
            0., /*forceDminFlag=*/false, 0., 0., 0.);

        glauber.calcRho(&nucleus);

        CHECK(anumForProfile(glauber, nucleus) == doctest::Approx(nucleus.A));
    }
}

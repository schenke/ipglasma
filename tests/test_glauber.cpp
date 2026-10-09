// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <cmath>
#include <string>

#include "Glauber.h"
#include "doctest.h"

namespace {
Nucleus lookup(const std::string &name) {
    Nucleus nucleus {};
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

TEST_CASE(
    "Glauber::findNucleusData: U and Xe have built-in deformations, "
    "independent of the beta2 argument") {
    Glauber glauber;
    struct Case {
        const char *name;
        double beta2, beta4;
    };
    for (const Case &c :
         {Case {"U", 0.28, 0.093}, Case {"Xe", 0.162, -0.003}}) {
        CAPTURE(c.name);
        Nucleus nucleus {};
        glauber.findNucleusData(
            &nucleus, c.name, /*setWSDeformParams=*/false, 0., 0., 0.5, 0., 0.,
            0., /*forceDminFlag=*/false, 0., 0., 0.);
        CHECK(nucleus.beta2 == doctest::Approx(c.beta2));
        CHECK(nucleus.beta4 == doctest::Approx(c.beta4));

        glauber.findNucleusData(
            &nucleus, c.name, /*setWSDeformParams=*/true, 6., 0.5, 0.2, 0.,
            0.05, 0., /*forceDminFlag=*/false, 0., 0., 0.);
        CHECK(nucleus.beta2 == doctest::Approx(0.2));
        CHECK(nucleus.beta4 == doctest::Approx(0.05));
    }
}

TEST_CASE("Glauber::findNucleusData: setWSDeformParams overrides beta2/R_WS") {
    Nucleus nucleus {};
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
        case 1:
            return glauber.anum2HO();
        case 2:
            return glauber.anum3Gauss(nucleus.R_WS);
        case 3:
            return glauber.anum3Fermi(nucleus.R_WS);
        case 8:
            return glauber.anumHulthen();
        default:
            return -1.0;  // unreachable for the species used below
    }
}
}  // namespace

TEST_CASE(
    "Glauber::calcRho: solved rho_WS reproduces A for every density profile") {
    // One species per density-profile branch calcRho dispatches on.
    for (const char *name : {"Pb", "S", "C", "d"}) {
        Nucleus nucleus {};
        Glauber glauber;
        glauber.findNucleusData(
            &nucleus, name, /*setWSDeformParams=*/false, 0., 0., 0., 0., 0., 0.,
            /*forceDminFlag=*/false, 0., 0., 0.);

        glauber.calcRho(&nucleus);

        CHECK(anumForProfile(glauber, nucleus) == doctest::Approx(nucleus.A));
    }
}

TEST_CASE(
    "Glauber::nuInS: the transverse integral of T(s) reproduces A for the "
    "Woods-Saxon and Hulthen profiles") {
    const double sigmaNN = 42.;  // mb
    for (const char *name : {"Pb", "d"}) {
        CAPTURE(name);
        Glauber glauber;
        // initGlauber sets sigma_NN; calcRho selects the nucleus nuInS uses
        glauber.initGlauber(
            sigmaNN, name, name, 0., /*setWSDeformParams=*/false, 0., 0., 0.,
            0., 0., 0., /*forceDminFlag=*/false, 0., 0., 0., 100);
        Nucleus nucleus {};
        glauber.findNucleusData(
            &nucleus, name, /*setWSDeformParams=*/false, 0., 0., 0., 0., 0., 0.,
            /*forceDminFlag=*/false, 0., 0., 0.);
        glauber.calcRho(&nucleus);
        const double A = nucleus.A;

        // nuInS includes a factor sigma_NN [fm^2]; integrate 2 pi s T(s)
        // with the trapezoidal rule out to where T is negligible.
        const double sigmaFm2 = sigmaNN * 0.1;
        const double sMax = 40.;
        const int nSteps = 4000;
        const double ds = sMax / nSteps;
        double integral = 0.;
        for (int i = 1; i <= nSteps; i++) {
            const double s = i * ds;
            const double weight = (i == nSteps) ? 0.5 : 1.0;
            integral += weight * 2. * M_PI * s * glauber.nuInS(s) * ds;
        }
        CHECK(integral / sigmaFm2 == doctest::Approx(A).epsilon(1e-3));
    }
}

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <array>
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
        double beta2, beta4, gamma;
    };
    for (const Case &c :
         {Case {"U", 0.28, 0.093, 0.},
          Case {"Xe", 0.207, 0., 26.93 * M_PI / 180.}}) {
        CAPTURE(c.name);
        Nucleus nucleus {};
        glauber.findNucleusData(
            &nucleus, c.name, /*setWSDeformParams=*/false, 0., 0., 0.5, 0., 0.,
            0., /*forceDminFlag=*/false, 0., 0., 0.);
        CHECK(nucleus.beta2 == doctest::Approx(c.beta2));
        CHECK(nucleus.beta4 == doctest::Approx(c.beta4));
        CHECK(nucleus.gamma == doctest::Approx(c.gamma));

        glauber.findNucleusData(
            &nucleus, c.name, /*setWSDeformParams=*/true, 6., 0.5, 0.2, 0.,
            0.05, 0., /*forceDminFlag=*/false, 0., 0., 0.);
        CHECK(nucleus.beta2 == doctest::Approx(0.2));
        CHECK(nucleus.beta4 == doctest::Approx(0.05));
        CHECK(nucleus.gamma == 0.);
    }
}

TEST_CASE("Glauber::findNucleusData: the built-in parameters decided in #177") {
    struct Case {
        const char *name;
        int A, Z;
        double R, w, a, beta2, beta3, beta4, gammaDegrees;
    };
    // 12C oscillator of De Jager et al. 1974; muonic 196Pt of Fricke et al.
    // 1995 (skin thickness 2.30 fm); Xe triaxial (Bally et al. 2022); the
    // others from De Vries et al. 1987 and FRDM
    const double aPt = 2.30 / (4. * std::log(3.));
    for (const Case &c : {
             Case {"C", 12, 6, 2.46, 1.247, 1.649, 0., 0., 0., 0.},
             Case {"O", 16, 8, 2.608, -0.051, 0.513, -0.01, -0.258, -0.122, 0.},
             Case {"Ne", 20, 10, 2.805, 0., 0.571, 0., 0., 0., 0.},
             Case {"Ar", 40, 18, 3.53, 0., 0.542, 0., 0., 0., 0.},
             Case {"Cu", 63, 29, 4.163, 0., 0.606, 0.162, 0., -0.006, 0.},
             Case {"Pt", 196, 78, 6.5472, 0., aPt, 0., 0., 0., 0.},
             Case {"Au", 197, 79, 6.37, 0., 0.535, -0.13, 0., -0.03, 0.},
             Case {"U", 238, 92, 6.874, 0., 0.556, 0.28, 0., 0.093, 0.},
             Case {"Xe", 129, 54, 5.601, 0., 0.492, 0.207, 0., 0., 26.93},
         }) {
        CAPTURE(c.name);
        const Nucleus nucleus = lookup(c.name);
        CHECK(nucleus.A == c.A);
        CHECK(nucleus.Z == c.Z);
        CHECK(nucleus.R_WS == doctest::Approx(c.R));
        CHECK(nucleus.w_WS == doctest::Approx(c.w));
        CHECK(nucleus.a_WS == doctest::Approx(c.a));
        CHECK(nucleus.beta2 == doctest::Approx(c.beta2));
        CHECK(nucleus.beta3 == doctest::Approx(c.beta3));
        CHECK(nucleus.beta4 == doctest::Approx(c.beta4));
        CHECK(nucleus.gamma == doctest::Approx(c.gammaDegrees * M_PI / 180.));
    }
    CHECK(aPt == doctest::Approx(0.52339).epsilon(1e-5));
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

TEST_CASE(
    "Glauber::findNucleusData: setWSDeformParams makes nuclei with A > 2 "
    "Woods-Saxon nuclei and keeps the proton and the deuteron") {
    Glauber glauber;
    struct Case {
        const char *name;
        bool overridden;
    };
    for (const Case &c :
         {Case {"p", false}, Case {"d", false}, Case {"He4", true},
          Case {"C", true}}) {
        CAPTURE(c.name);
        Nucleus builtIn {};
        glauber.findNucleusData(
            &builtIn, c.name, /*setWSDeformParams=*/false, 0., 0., 0., 0., 0.,
            0., /*forceDminFlag=*/false, 0., 0., 0.);
        Nucleus nucleus {};
        glauber.findNucleusData(
            &nucleus, c.name, /*setWSDeformParams=*/true, 6.6, 0.55, 0.1, 0.,
            0., 0., /*forceDminFlag=*/false, 0., 0.1, 0.);
        if (c.overridden) {
            CHECK(nucleus.densityFunc == 3);
            CHECK(nucleus.R_WS == doctest::Approx(6.6));
            CHECK(nucleus.a_WS == doctest::Approx(0.55));
            CHECK(nucleus.w_WS == 0.);
            CHECK(nucleus.beta2 == doctest::Approx(0.1));
            CHECK(nucleus.dR_np == doctest::Approx(0.1));
        } else {
            CHECK(nucleus.densityFunc == builtIn.densityFunc);
            CHECK(nucleus.R_WS == builtIn.R_WS);
            CHECK(nucleus.a_WS == builtIn.a_WS);
            CHECK(nucleus.w_WS == builtIn.w_WS);
            CHECK(nucleus.beta2 == builtIn.beta2);
            CHECK(nucleus.dR_np == 0.);
        }
    }
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

TEST_CASE(
    "Glauber::interNuPInSP, interNuTInST: the thickness tables shared by all "
    "instances follow the Woods-Saxon parameters and sigma_NN") {
    // a posterior file can change these from event to event, each event
    // with its own Glauber
    auto thickness = [](double sigmaNN, double radius, double diffuseness) {
        Glauber glauber;
        glauber.initGlauber(
            sigmaNN, "Pb", "Pb", 0., /*setWSDeformParams=*/true, radius,
            diffuseness, 0., 0., 0., 0., /*forceDminFlag=*/false, 0., 0., 0.,
            100);
        return std::array<double, 2> {
            glauber.interNuPInSP(5.), glauber.interNuTInST(5.)};
    };
    const std::array<double, 2> first = thickness(67., 6.62, 0.546);
    for (const std::array<double, 2> &other :
         {thickness(67., 6.8, 0.546), thickness(67., 6.62, 0.6),
          thickness(50., 6.62, 0.546)}) {
        CHECK(other[0] != first[0]);
        CHECK(other[1] != first[1]);
    }
    // the same values again after the tables were rebuilt for the others
    CHECK(thickness(67., 6.62, 0.546) == first);
}

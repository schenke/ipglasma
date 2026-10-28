// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <gsl/gsl_sf_bessel.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "GluonMultiplicity.h"
#include "Group.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "doctest.h"
#include "test_helpers.h"

namespace {
const int N = 8;
const int it = 4;

// a = 1 fm and maxTime = it * a * dtau, so that it is the final step
void makeMultiplicityTestParam(Parameters &param) {
    makeLatticeParam(param, N);
    param.lattice.L = N;
    param.coupling.g = 1.;
    param.run.dtau = 0.1;
    param.evolution.maxTime = 0.4;
    param.evolution.inverseQsForMaxTime = false;
    param.coupling.runningCoupling = false;
    param.output.writeHadronSpectrum = false;
}

void setVacuum(Lattice &lat) {
    for (int pos = 0; pos < N * N; pos++) {
        lat.Ux[pos] = Matrix(1.);
        lat.Uy[pos] = Matrix(1.);
        lat.U[pos] = Matrix(0.);
        lat.U2[pos] = Matrix(0.);
        lat.Ux2[pos] = Matrix(0.);
        lat.Uy2[pos] = Matrix(0.);
    }
}

void removeOutput() {
    std::remove("NpartdNdy-t0.4-0.dat");
    std::remove("gluonMultiplicity0.json");
}

bool fileExists(const std::string &name) { return std::ifstream(name).good(); }
}  // namespace

TEST_CASE("GluonMultiplicity::compute reports no collision in the vacuum") {
    Parameters param;
    makeMultiplicityTestParam(param);
    Lattice lat(&param, N);
    setVacuum(lat);
    Group group;
    int nn[2] = {N, N};
    GluonMultiplicity multiplicity(nn);
    removeOutput();
    param.event.success = 1;

    CHECK(multiplicity.compute(&lat, &group, &param, it) == 0);
    CHECK_FALSE(fileExists("NpartdNdy-t0.4-0.dat"));
    removeOutput();
}

TEST_CASE(
    "GluonMultiplicity::compute measures a positive multiplicity for "
    "electric fields and writes the event files") {
    Parameters param;
    makeMultiplicityTestParam(param);
    Lattice lat(&param, N);
    setVacuum(lat);
    Group group;
    std::mt19937 rng(5);
    std::normal_distribution<double> gauss(0., 0.3);
    for (int pos = 0; pos < N * N; pos++) {
        for (Matrix *field : {&lat.U[pos], &lat.U2[pos]}) {
            for (int a = 0; a < 8; a++) *field += gauss(rng) * group.getT(a);
        }
    }
    int nn[2] = {N, N};
    GluonMultiplicity multiplicity(nn);
    removeOutput();
    param.event.success = 0;

    REQUIRE(multiplicity.compute(&lat, &group, &param, it) == 1);
    CHECK(param.event.success == 1);
    std::ifstream in("NpartdNdy-t0.4-0.dat");
    REQUIRE(in.good());
    double Npart, dNdy, Tpp, b, dEdy;
    in >> Npart >> dNdy >> Tpp >> b >> dEdy;
    CHECK(dNdy > 0.);
    CHECK(dEdy > 0.);
    // the 14 columns documented in OUTPUT.md: the estimator in 7 (the
    // default, the projection), placeholders in 8 and 9
    std::string seed, estimator, na2, na3;
    in >> seed >> estimator >> na2 >> na3;
    CHECK(estimator == "1");
    CHECK(na2 == "N/A");
    CHECK(na3 == "N/A");
    int rest = 0;
    for (std::string token; in >> token;) rest++;
    CHECK(rest == 5);

    // the JSON file has the documented keys and four arrays of 100 bins
    std::ifstream json("gluonMultiplicity0.json");
    REQUIRE(json.good());
    std::stringstream text;
    text << json.rdbuf();
    const std::string content = text.str();
    for (const char *key :
         {"\"format\": \"ipglasma-gluon-target\"", "\"dN\":", "\"dE_GeV\":",
          "\"tau_fm\":", "\"rapidity_variable\": \"y\"",
          "\"estimator\": \"bjorken_projection\"",
          "\"dN_electric_estimator\":"}) {
        CAPTURE(key);
        CHECK(content.find(key) != std::string::npos);
    }
    for (const char *array :
         {"\"kt_GeV\": [", "\"dN_d2k_GeV_minus2\": [",
          "\"dE_d2k_GeV_minus1\": [", "\"lattice_bin_counts\": ["}) {
        CAPTURE(array);
        const std::size_t start = content.find(array);
        REQUIRE(start != std::string::npos);
        const std::size_t end = content.find(']', start);
        const std::string values = content.substr(start, end - start);
        CHECK(std::count(values.begin(), values.end(), ',') == 99);
    }
    removeOutput();
}

TEST_CASE(
    "GluonMultiplicity::compute normalizes the electric and pi terms at the "
    "same time tau_n = it dtau") {
    // The same field once as E^x and once as pi: the occupation numbers
    // differ only by the factors g^2/tau of the E^i term and tau of the pi
    // term, so dN_E/dN_pi = g^2/tau^2. After Evolution::run()'s final half
    // step both momenta are at tau_n = it dtau (with tau_{n-1/2}, the ratio
    // would be 1/0.35^2 = 8.16 instead of 1/0.4^2 = 6.25 here).
    auto multiplicityOf = [](bool asPi) {
        Parameters param;
        makeMultiplicityTestParam(param);
        // the electric estimator, 2|E|^2/omega, whose two terms are tested
        param.output.gluonSpectrumEstimator = 0;
        Lattice lat(&param, N);
        setVacuum(lat);
        Group group;
        std::mt19937 rng(11);
        std::normal_distribution<double> gauss(0., 0.3);
        for (int pos = 0; pos < N * N; pos++) {
            Matrix field(0.);
            for (int a = 0; a < 8; a++) field += gauss(rng) * group.getT(a);
            (asPi ? lat.Ux2[pos] : lat.U[pos]) = field;
        }
        int nn[2] = {N, N};
        GluonMultiplicity multiplicity(nn);
        removeOutput();
        REQUIRE(multiplicity.compute(&lat, &group, &param, it) == 1);
        std::ifstream in("NpartdNdy-t0.4-0.dat");
        REQUIRE(in.good());
        double Npart, dNdy, Tpp, b, dEdy;
        in >> Npart >> dNdy >> Tpp >> b >> dEdy;
        in.close();
        removeOutput();
        return std::make_pair(dNdy, dEdy);
    };
    const auto electric = multiplicityOf(false);
    const auto pi = multiplicityOf(true);
    Parameters param;
    makeMultiplicityTestParam(param);
    const double tauN = it * param.run.dtau;
    const double g2 = param.coupling.g * param.coupling.g;
    CHECK(
        electric.first / pi.first
        == doctest::Approx(g2 / (tauN * tauN)).epsilon(1e-4));
    CHECK(
        electric.second / pi.second
        == doctest::Approx(g2 / (tauN * tauN)).epsilon(1e-4));
}

TEST_CASE(
    "GluonMultiplicity::bjorkenProjectionCoefficients: the regular free modes "
    "hold 1/pi gluons "
    "at every k tau") {
    // the boost-invariant modes regular at tau = 0: transverse
    // q = J0(k tau), p = tau dq/dtau = -x J1(x); longitudinal
    // q = tau J1(k tau), p = (1/tau) dq/dtau = k J0(x). The projection gives
    // 1/pi independent of tau, while 2|E|^2/omega (|p|^2/(k tau) and
    // tau |p|^2/k) gives pi x J1(x)^2 and pi x J0(x)^2 times that.
    const double k = 1.7;
    for (const double x : {0.25, 0.5, 1., 1.5, 2., 3., 4., 6., 10.}) {
        CAPTURE(x);
        const double tau = x / k;
        std::complex<double> cq, cp;
        GluonMultiplicity::bjorkenProjectionCoefficients(k, tau, false, cq, cp);
        const double qT = gsl_sf_bessel_J0(x), pT = -x * gsl_sf_bessel_J1(x);
        CHECK(std::norm(cq * qT + cp * pT) == doctest::Approx(1. / M_PI));
        GluonMultiplicity::bjorkenProjectionCoefficients(k, tau, true, cq, cp);
        const double qL = tau * gsl_sf_bessel_J1(x);
        const double pL = k * gsl_sf_bessel_J0(x);
        CHECK(std::norm(cq * qL + cp * pL) == doctest::Approx(1. / M_PI));
        // the electric estimator against the exact number
        CHECK(
            M_PI * pT * pT / (k * tau)
            == doctest::Approx(M_PI * x * std::pow(gsl_sf_bessel_J1(x), 2)));
        CHECK(
            M_PI * tau * pL * pL / k
            == doctest::Approx(M_PI * x * std::pow(gsl_sf_bessel_J0(x), 2)));
    }
    // two values of the table of #140: the longitudinal mode at x = 2 and
    // the transverse one at x = 4 are counted 0.315 and 0.055 times
    CHECK(
        M_PI * 2. * std::pow(gsl_sf_bessel_J0(2.), 2)
        == doctest::Approx(0.3150).epsilon(1e-3));
    CHECK(
        M_PI * 4. * std::pow(gsl_sf_bessel_J1(4.), 2)
        == doctest::Approx(0.0547).epsilon(1e-2));
}

TEST_CASE(
    "bjorkenProjectionCoefficients: k = 0 or tau = 0, where Y_0 diverges, "
    "throw") {
    std::complex<double> cq, cp;
    CHECK_THROWS_AS(
        GluonMultiplicity::bjorkenProjectionCoefficients(1., 0., false, cq, cp),
        std::invalid_argument);
    CHECK_THROWS_AS(
        GluonMultiplicity::bjorkenProjectionCoefficients(0., 1., true, cq, cp),
        std::invalid_argument);
}

TEST_CASE(
    "GluonMultiplicity::compute with the projection counts 1/pi gluons per "
    "unit amplitude in free plane waves") {
    // A transverse plane wave in the links and E^i and a longitudinal one in
    // phi and pi, each the free mode regular at tau = 0 at tau_n = it dtau,
    // in one color direction: the projection gives N^2 |a|^2/pi gluons in
    // each of the modes k and -k of a wave of amplitude a. A wrong sign of
    // the Coulomb momentum c_j = 1 - e^{-i k_j} or a wrong momentum label of
    // the FFT would move the transverse wave into the longitudinal part.
    Parameters param;
    makeMultiplicityTestParam(param);
    param.output.gluonSpectrumEstimator = 1;
    Lattice lat(&param, N);
    setVacuum(lat);
    Group group;
    const Matrix &T = group.getT(0);
    const double g = param.coupling.g;
    const double tau = it * param.run.dtau;
    const std::complex<double> I(0., 1.);
    const double a = 1e-3, b = 2e-3;
    // transverse: the mode (2, 1), A = t a e^{ik.x} + c.c.,
    // t = (-c_y, c_x)/k~
    const double kx = 2. * M_PI * 2. / N, ky = 2. * M_PI * 1. / N;
    const std::complex<double> cx = 1. - std::exp(-I * kx);
    const std::complex<double> cy = 1. - std::exp(-I * ky);
    const double kt = std::sqrt(std::norm(cx) + std::norm(cy));
    const double x = kt * tau;
    // longitudinal: the mode (1, 3)
    const double lx = 2. * M_PI * 1. / N, ly = 2. * M_PI * 3. / N;
    const double lt =
        2.
        * std::sqrt(
            std::pow(std::sin(lx / 2.), 2) + std::pow(std::sin(ly / 2.), 2));
    const double xl = lt * tau;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            const int pos = lat.positionFromXY(i, j);
            const std::complex<double> wave =
                a * std::exp(I * (kx * i + ky * j));
            const double Ax = 2. * std::real(-cy / kt * wave);
            const double Ay = 2. * std::real(cx / kt * wave);
            // q = J_0(x), p = g E = -x J_1(x)
            lat.Ux[pos] =
                Matrix::expiHermitian(g * Ax * gsl_sf_bessel_J0(x) * T);
            lat.Uy[pos] =
                Matrix::expiHermitian(g * Ay * gsl_sf_bessel_J0(x) * T);
            lat.U[pos] = (-x * gsl_sf_bessel_J1(x) * Ax / g) * T;
            lat.U2[pos] = (-x * gsl_sf_bessel_J1(x) * Ay / g) * T;
            // q = phi = tau J_1(x), p = pi = k~ J_0(x)
            const double l = 2. * b * std::cos(lx * i + ly * j);
            lat.Uy2[pos] = (l * tau * gsl_sf_bessel_J1(xl)) * T;
            lat.Ux2[pos] = (l * lt * gsl_sf_bessel_J0(xl)) * T;
        }
    }
    int nn[2] = {N, N};
    GluonMultiplicity multiplicity(nn);
    removeOutput();
    REQUIRE(multiplicity.compute(&lat, &group, &param, it) == 1);
    std::ifstream in("NpartdNdy-t0.4-0.dat");
    REQUIRE(in.good());
    double Npart, dNdy;
    in >> Npart >> dNdy;
    in.close();
    removeOutput();
    CHECK(
        dNdy
        == doctest::Approx(2. * N * N * (a * a + b * b) / M_PI).epsilon(1e-6));
}

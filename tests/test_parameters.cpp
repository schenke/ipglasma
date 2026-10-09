// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "Glauber.h"
#include "Parameters.h"
#include "doctest.h"

namespace {
// Parameters' default constructor does not initialize every member, so every
// field validationErrors() reads must be set explicitly.
//
// Takes an out-parameter rather than returning by value: Parameters holds a
// PrettyOstream member, which holds a non-copyable std::ostringstream.
//
// Checks of a single value run while reading the input file and are tested
// in test_input_file.cpp; these tests cover the checks that combine several
// parameters.
void makeValidBaseline(Parameters &param) {
    param.collision.useNucleus = true;
    param.wilsonLines.writeWilsonLines = 2;
    param.wilsonLines.wilsonLinePath = ".";
    param.jimwlk.saveSnapshots = 0;
    param.coupling.runningCoupling = false;
    param.coupling.mu0 = 0.3;
    param.coupling.LambdaQCD = 0.2;
    param.jimwlk.enabled = 0;
    param.jimwlk.alphaS = 0.3;  // fixed JIMWLK coupling
    param.jimwlk.mu0 = 0.28;
    param.jimwlk.LambdaQCD = 0.04;
}
}  // namespace

TEST_CASE("Parameters::validationErrors: accepts a normal configuration") {
    Parameters param;
    makeValidBaseline(param);
    CHECK(param.validationErrors().empty());
}

TEST_CASE(
    "Parameters::validationErrors: rejects jimwlkSaveSnapshots without "
    "writeWilsonLines") {
    Parameters param;
    makeValidBaseline(param);
    param.jimwlk.saveSnapshots = 1;
    param.wilsonLines.writeWilsonLines = 0;
    CHECK(param.validationErrors().size() == 1);

    param.jimwlk.saveSnapshots = 0;
    CHECK(param.validationErrors().empty());
}

TEST_CASE(
    "Parameters::validationErrors: rejects a missing Wilson-line directory "
    "only when Wilson lines are written") {
    Parameters param;
    makeValidBaseline(param);
    param.wilsonLines.wilsonLinePath = "this_directory_does_not_exist";
    CHECK(param.validationErrors().size() == 1);
    param.wilsonLines.writeWilsonLines = 0;
    CHECK(param.validationErrors().empty());
}

TEST_CASE(
    "Parameters::validationErrors: rejects LambdaQCD >= mu0 only with "
    "running coupling") {
    Parameters param;
    makeValidBaseline(param);
    param.coupling.mu0 = 0.2;
    param.coupling.LambdaQCD = 0.2;  // equal: log argument is 0 at the boundary
    CHECK(param.validationErrors().empty());

    param.coupling.runningCoupling = true;
    CHECK(param.validationErrors().size() == 1);
    param.coupling.LambdaQCD = 0.3;  // larger: log argument is negative
    CHECK(param.validationErrors().size() == 1);
    param.coupling.LambdaQCD = 0.1;
    CHECK(param.validationErrors().empty());
}

TEST_CASE(
    "Parameters::validationErrors: rejects jimwlkLambdaQCD >= jimwlkMu0 "
    "only with the JIMWLK running coupling") {
    Parameters param;
    makeValidBaseline(param);
    param.jimwlk.LambdaQCD = 0.3;             // >= jimwlkMu0
    CHECK(param.validationErrors().empty());  // JIMWLK off
    param.jimwlk.enabled = 1;
    CHECK(param.validationErrors().empty());  // fixed JIMWLK coupling
    param.jimwlk.alphaS = 0.;
    CHECK(param.validationErrors().size() == 1);
}

TEST_CASE(
    "Parameters::validationErrors: rejects averaging over nuclei with a "
    "proton") {
    Parameters param;
    makeValidBaseline(param);
    param.collision.projectile = "p";
    param.collision.target = "Pb";
    param.collision.nucleiToAverage = 2;
    CHECK(param.validationErrors().size() == 1);

    param.collision.projectile = "Pb";
    CHECK(param.validationErrors().empty());
    param.collision.target = "p";
    CHECK(param.validationErrors().size() == 1);
    param.collision.nucleiToAverage = 1;
    CHECK(param.validationErrors().empty());
}

TEST_CASE("Parameters::validationErrors: output times must be before maxTime") {
    Parameters param;
    makeValidBaseline(param);
    param.evolution.maxTime = 0.4;
    param.evolution.inverseQsForMaxTime = false;
    param.output.outputTimes = {0.1, 0.3};
    CHECK(param.validationErrors().empty());
    param.output.outputTimes = {0.1, 0.4, 0.5};
    CHECK(param.validationErrors().size() == 2);
    // with 1/<Qs> the final time is not known yet
    param.evolution.inverseQsForMaxTime = true;
    param.collision.useNucleus = true;
    param.wilsonLines.readInitialWilsonLines = 0;
    CHECK(param.validationErrors().empty());
}

TEST_CASE("Parameters::validationErrors: reports every failed check") {
    Parameters param;
    makeValidBaseline(param);
    param.jimwlk.saveSnapshots = 1;
    param.wilsonLines.writeWilsonLines = 0;
    param.coupling.runningCoupling = true;
    param.coupling.LambdaQCD = 0.5;
    CHECK(param.validationErrors().size() == 2);
}

TEST_CASE(
    "Parameters::validationErrors: inverseQsForMaxTime and running coupling "
    "need the event-averaged Qs") {
    for (int setting = 0; setting < 2; ++setting) {
        CAPTURE(setting);
        Parameters param;
        makeValidBaseline(param);
        if (setting == 0) {
            param.evolution.inverseQsForMaxTime = true;
        } else {
            param.coupling.runningCoupling = true;
            param.coupling.mu0 = 0.3;
            param.coupling.LambdaQCD = 0.2;
        }
        CHECK(param.validationErrors().empty());
        param.collision.useNucleus = false;
        CHECK(param.validationErrors().size() == 1);
        // read Wilson lines come with their nuclei's geometry
        param.collision.useNucleus = true;
        param.wilsonLines.readInitialWilsonLines = 2;
        CHECK(param.validationErrors().empty());
    }
}

TEST_CASE(
    "Parameters::validationErrors: JIMWLK evolves sampled Wilson lines from "
    "jimwlkInitialX only to smaller x") {
    Parameters param;
    makeValidBaseline(param);
    param.jimwlk.enabled = 1;
    param.jimwlk.initialX = 0.01;
    param.colorCharge.projectileX = 1e-3;
    param.colorCharge.targetX = 0.01;  // no evolution
    CHECK(param.validationErrors().empty());
    param.colorCharge.targetX = 0.02;
    std::vector<std::string> errors = param.validationErrors();
    REQUIRE(errors.size() == 1);
    CHECK(
        errors[0].find("must not be larger than jimwlkInitialX (0.01)")
        != std::string::npos);
    param.colorCharge.targetX = 1e-3;
    param.colorCharge.projectileX = 0.02;
    errors = param.validationErrors();
    REQUIRE(errors.size() == 1);
    CHECK(
        errors[0].find("must not be larger than jimwlkInitialX")
        != std::string::npos);
}

TEST_CASE(
    "Parameters::validationErrors: JIMWLK evolves read Wilson lines from "
    "readWilsonLinesX only to smaller x") {
    Parameters param;
    makeValidBaseline(param);
    param.jimwlk.enabled = 1;
    param.jimwlk.initialX = 0.01;
    param.colorCharge.projectileX = 1e-3;
    param.colorCharge.targetX = 2e-3;
    param.wilsonLines.readInitialWilsonLines = 2;
    param.wilsonLines.readX = 5e-3;
    CHECK(param.validationErrors().empty());
    param.wilsonLines.readX = 1.5e-3;  // below targetX
    std::vector<std::string> errors = param.validationErrors();
    REQUIRE(errors.size() == 1);
    CHECK(
        errors[0].find("must not be larger than readWilsonLinesX (0.0015)")
        != std::string::npos);
    // the read x, not jimwlkInitialX, is where JIMWLK starts
    param.wilsonLines.readX = 5e-3;
    param.jimwlk.initialX = 1e-3;
    CHECK(param.validationErrors().empty());
    param.wilsonLines.readX = 0.;  // the initial lines, at jimwlkInitialX
    CHECK(param.validationErrors().size() == 1);
    param.jimwlk.initialX = 0.01;
    CHECK(param.validationErrors().empty());
    param.jimwlk.enabled = 0;
    param.wilsonLines.readX = 1e-4;
    CHECK(param.validationErrors().empty());
}

TEST_CASE(
    "Parameters::validationErrors: posterior parameter sets require hot-spot "
    "nucleons") {
    Parameters param;
    makeValidBaseline(param);
    for (int type : {1, 2, 4}) {
        CAPTURE(type);
        param.subnucleon.subNucleonParamType = type;
        param.subnucleon.nucleonModel = "hotspots";
        CHECK(param.validationErrors().empty());
        param.subnucleon.nucleonModel = "gaussian";
        CHECK(param.validationErrors().size() == 1);
    }
}

namespace {
class TempCsvFile {
  public:
    explicit TempCsvFile(const std::string &contents) {
        path_ = "ipglasma_test_parameters_tmp_posterior.csv";
        std::ofstream out(path_);
        out << contents;
    }
    ~TempCsvFile() { std::remove(path_.c_str()); }
    const std::string &path() const { return path_; }

  private:
    std::string path_;
};
}  // namespace

TEST_CASE(
    "Parameters::loadPosteriorParameterSetsFromFile parses a CSV, skipping the "
    "header") {
    TempCsvFile file(
        "m,BG,BGq,smearingWidth,NqBase,QsmuRatio,dqmin\n"
        "0.4,3.3,0.3,0.6,3,0.643,0.2\n"
        "0.5,3.1,0.25,0.55,4,0.6,0.3\n");
    Parameters param;

    std::vector<std::vector<float>> parsed;
    CHECK(param.loadPosteriorParameterSetsFromFile(file.path(), parsed) == "");

    REQUIRE(parsed.size() == 2);
    REQUIRE(parsed[0].size() == 7);
    CHECK(parsed[0][0] == doctest::Approx(0.4));
    CHECK(parsed[0][1] == doctest::Approx(3.3));
    CHECK(parsed[1][4] == doctest::Approx(4.0));
    CHECK(parsed[1][6] == doctest::Approx(0.3));
}

TEST_CASE(
    "Parameters::validationErrors: rejects an x above the Q_s table only "
    "where Q_s^2 is looked up at it") {
    Parameters param;
    makeValidBaseline(param);
    param.wilsonLines.readInitialWilsonLines = 0;
    param.colorCharge.useFluctuatingX = false;
    param.colorCharge.projectileX = 0.01;
    param.colorCharge.targetX = 0.02;
    CHECK(param.validationErrors().size() == 1);
    param.colorCharge.projectileX = 0.02;
    param.colorCharge.targetX = 0.01;
    CHECK(param.validationErrors().size() == 1);

    // a fluctuating x does not use projectileX and targetX
    param.colorCharge.useFluctuatingX = true;
    CHECK(param.validationErrors().empty());
    // with JIMWLK, Q_s^2 is looked up at jimwlkInitialX
    param.colorCharge.useFluctuatingX = false;
    param.jimwlk.enabled = 1;
    param.jimwlk.initialX = 0.01;
    param.colorCharge.projectileX = 1e-3;
    param.colorCharge.targetX = 1e-3;
    CHECK(param.validationErrors().empty());
    param.jimwlk.initialX = 0.02;
    CHECK(param.validationErrors().size() == 1);
    param.jimwlk.enabled = 0;
    param.colorCharge.projectileX = 0.02;
    // no color charges are sampled
    param.wilsonLines.readInitialWilsonLines = 2;
    CHECK(param.validationErrors().empty());
    param.wilsonLines.readInitialWilsonLines = 0;
    param.collision.useNucleus = false;
    CHECK(param.validationErrors().empty());
}

TEST_CASE(
    "Parameters::loadPosteriorParameterSetsFromFile reports a missing file "
    "and a value that is not a number") {
    Parameters param;
    std::vector<std::vector<float>> parsed;
    CHECK(
        param.loadPosteriorParameterSetsFromFile(
            "no_such_posterior_file.csv", parsed)
        == "cannot open the posterior parameter file "
           "no_such_posterior_file.csv");

    TempCsvFile file("m,BG\n0.4,3.3\n0.5,x\n");
    const std::string problem =
        param.loadPosteriorParameterSetsFromFile(file.path(), parsed);
    CHECK(problem == file.path() + ":3: x is not a number");
}

TEST_CASE(
    "Parameters::initialX: jimwlkInitialX with JIMWLK, none with a "
    "fluctuating x, otherwise projectileX/targetX") {
    Parameters param;
    param.jimwlk.enabled = false;
    param.colorCharge.useFluctuatingX = false;
    param.colorCharge.projectileX = 1e-3;
    param.colorCharge.targetX = 2e-3;
    param.jimwlk.initialX = 0.005;
    CHECK(param.initialX(NucleusRole::Projectile) == 1e-3);
    CHECK(param.initialX(NucleusRole::Target) == 2e-3);
    param.colorCharge.useFluctuatingX = true;
    CHECK(param.initialX(NucleusRole::Projectile) < 0.);
    CHECK(param.initialX(NucleusRole::Target) < 0.);
    param.colorCharge.useFluctuatingX = false;
    param.jimwlk.enabled = true;
    CHECK(param.initialX(NucleusRole::Projectile) == 0.005);
    CHECK(param.initialX(NucleusRole::Target) == 0.005);
}

TEST_CASE(
    "Parameters::validationErrors: rejects JIMWLK together with a "
    "fluctuating x") {
    Parameters param;
    makeValidBaseline(param);
    param.jimwlk.enabled = 1;
    param.jimwlk.initialX = 0.01;
    param.colorCharge.useFluctuatingX = false;
    CHECK(param.validationErrors().empty());
    param.colorCharge.useFluctuatingX = true;
    const std::vector<std::string> errors = param.validationErrors();
    REQUIRE(errors.size() == 1);
    CHECK(
        errors[0].find("useJIMWLK = 1 and useFluctuatingX = 1 are mutually "
                       "exclusive")
        != std::string::npos);
    param.jimwlk.enabled = 0;
    CHECK(param.validationErrors().empty());
}

TEST_CASE(
    "Parameters::xBeforeJimwlk: readWilsonLinesX if Wilson lines are read "
    "and it is set, otherwise the initial x") {
    Parameters param;
    param.jimwlk.enabled = false;
    param.colorCharge.useFluctuatingX = false;
    param.colorCharge.projectileX = 1e-3;
    param.wilsonLines.readInitialWilsonLines = 2;
    param.wilsonLines.readX = 0.;
    CHECK(param.xBeforeJimwlk(NucleusRole::Projectile) == 1e-3);
    param.wilsonLines.readX = 2e-4;
    CHECK(param.xBeforeJimwlk(NucleusRole::Projectile) == 2e-4);
    CHECK(param.xBeforeJimwlk(NucleusRole::Target) == 2e-4);
    // not read: the initial x, jimwlkInitialX where JIMWLK starts
    param.wilsonLines.readInitialWilsonLines = 0;
    param.jimwlk.enabled = true;
    param.jimwlk.initialX = 0.01;
    CHECK(param.xBeforeJimwlk(NucleusRole::Projectile) == 0.01);
}

TEST_CASE("Parameters::validationErrors: bMin must not be larger than bMax") {
    Parameters param;
    makeValidBaseline(param);
    param.collision.bMin = 3.;
    param.collision.bMax = 2.;
    const std::vector<std::string> errors = param.validationErrors();
    REQUIRE(errors.size() == 1);
    CHECK(errors[0] == "bMin (3) must not be larger than bMax (2)");
    param.collision.bMax = 3.;
    CHECK(param.validationErrors().empty());
    // b is not sampled without nuclei
    param.collision.bMax = 2.;
    param.collision.useNucleus = false;
    CHECK(param.validationErrors().empty());
}

TEST_CASE(
    "Parameters::validationErrors: lightNucleusOption must select a "
    "configuration file when the files are read") {
    CHECK(speciesMassNumber("Ar") == 40);
    CHECK(speciesMassNumber("p") == 1);
    CHECK(speciesMassNumber("Xx") == 0);

    Parameters param;
    makeValidBaseline(param);
    param.collision.projectile = "Ar";
    param.collision.target = "Pb";
    param.nucleus.nucleonPositionsFromFile = true;
    param.nucleus.lightNucleusOption = 5;
    std::vector<std::string> errors = param.validationErrors();
    REQUIRE(errors.size() == 1);
    CHECK(errors[0]
          == "projectile: lightNucleusOption 5 has no configuration file for "
             "A = 40; allowed values: 0 4");

    // the same species on both sides is reported once
    param.collision.target = "Ar";
    errors = param.validationErrors();
    REQUIRE(errors.size() == 1);
    CHECK(errors[0].rfind("projectile and target: ", 0) == 0);

    // an option that exists, or files that are not read, are fine
    param.nucleus.lightNucleusOption = 4;
    CHECK(param.validationErrors().empty());
    param.nucleus.lightNucleusOption = 5;
    param.wilsonLines.readInitialWilsonLines = 2;
    CHECK(param.validationErrors().empty());
    param.wilsonLines.readInitialWilsonLines = 0;
    param.nucleus.nucleonPositionsFromFile = false;
    CHECK(param.validationErrors().empty());
}

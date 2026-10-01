#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

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
        param.collision.useNucleus = true;
        param.wilsonLines.readInitialWilsonLines = 2;
        CHECK(param.validationErrors().size() == 1);
    }
}

TEST_CASE(
    "Parameters::validationErrors: posterior types 2 and 4 need Nq 0 or 3") {
    Parameters param;
    makeValidBaseline(param);
    for (int type : {2, 4}) {
        CAPTURE(type);
        param.subnucleon.subNucleonParamType = type;
        for (double nq : {0., 3.}) {
            param.subnucleon.Nq = nq;
            CHECK(param.validationErrors().empty());
        }
        param.subnucleon.Nq = 5.;
        CHECK(param.validationErrors().size() == 1);
    }
    param.subnucleon.subNucleonParamType = 1;  // variable Nq
    CHECK(param.validationErrors().empty());
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
    param.loadPosteriorParameterSetsFromFile(file.path(), parsed);

    REQUIRE(parsed.size() == 2);
    REQUIRE(parsed[0].size() == 7);
    CHECK(parsed[0][0] == doctest::Approx(0.4));
    CHECK(parsed[0][1] == doctest::Approx(3.3));
    CHECK(parsed[1][4] == doctest::Approx(4.0));
    CHECK(parsed[1][6] == doctest::Approx(0.3));
}

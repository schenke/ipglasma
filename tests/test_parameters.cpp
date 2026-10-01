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
    param.wilsonLines.writeWilsonLines = 2;
    param.wilsonLines.wilsonLinePath = ".";
    param.jimwlk.saveSnapshots = 0;
    param.coupling.runningCoupling = false;
    param.coupling.muZero = 0.3;
    param.coupling.LambdaQCD = 0.2;
    param.jimwlk.useJIMWLK = 0;
    param.jimwlk.alphas_jimwlk = 0.3;  // fixed JIMWLK coupling
    param.jimwlk.mu0_jimwlk = 0.28;
    param.jimwlk.Lambda_QCD_jimwlk = 0.04;
}
}  // namespace

TEST_CASE("Parameters::validationErrors: accepts a normal configuration") {
    Parameters param;
    makeValidBaseline(param);
    CHECK(param.validationErrors().empty());
    CHECK(param.ValidParameters() == true);
}

TEST_CASE(
    "Parameters::validationErrors: rejects saveSnapshots without "
    "writeWilsonLines") {
    Parameters param;
    makeValidBaseline(param);
    param.jimwlk.saveSnapshots = 1;
    param.wilsonLines.writeWilsonLines = 0;
    CHECK(param.validationErrors().size() == 1);
    CHECK(param.ValidParameters() == false);

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
    "Parameters::validationErrors: rejects LambdaQCD >= muZero only with "
    "running coupling") {
    Parameters param;
    makeValidBaseline(param);
    param.coupling.muZero = 0.2;
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
    "Parameters::validationErrors: rejects Lambda_QCD_jimwlk >= mu0_jimwlk "
    "only with the JIMWLK running coupling") {
    Parameters param;
    makeValidBaseline(param);
    param.jimwlk.Lambda_QCD_jimwlk = 0.3;     // >= mu0_jimwlk
    CHECK(param.validationErrors().empty());  // JIMWLK off
    param.jimwlk.useJIMWLK = 1;
    CHECK(param.validationErrors().empty());  // fixed JIMWLK coupling
    param.jimwlk.alphas_jimwlk = 0.;
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

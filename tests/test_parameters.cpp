#include "doctest.h"

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "Parameters.h"

namespace {
// Parameters' default constructor does not initialize every member (no
// default member initializers), so every field ValidParameters() reads
// must be set explicitly before calling it, or the test would depend on
// indeterminate values.
//
// Takes an out-parameter rather than returning by value: Parameters now
// holds a PrettyOstream member (for ValidParameters()'s own error
// messages), and PrettyOstream holds a std::ostringstream, which is not
// copyable or movable -- exactly like Lattice, which already explicitly
// deletes its copy/move for the same reason.
void makeValidBaseline(Parameters &param) {
    param.setSize(256);
    param.setWriteWilsonLines(2);
    param.setSaveSnapshots(0);
}
}  // namespace

TEST_CASE("Parameters::ValidParameters: accepts a normal configuration") {
    Parameters param;
    makeValidBaseline(param);
    CHECK(param.ValidParameters() == true);
}

TEST_CASE("Parameters::ValidParameters: rejects a non-positive lattice size") {
    Parameters param;
    makeValidBaseline(param);
    param.setSize(0);
    CHECK(param.ValidParameters() == false);
}

TEST_CASE("Parameters::ValidParameters: rejects an invalid Wilson-line format") {
    Parameters param;
    makeValidBaseline(param);
    param.setWriteWilsonLines(3);  // only 0 (off), 1 (text), 2 (binary) are valid
    CHECK(param.ValidParameters() == false);
}

TEST_CASE("Parameters::ValidParameters: writeWilsonLines=0 is valid by itself") {
    Parameters param;
    makeValidBaseline(param);
    param.setWriteWilsonLines(0);
    CHECK(param.ValidParameters() == true);
}

TEST_CASE(
    "Parameters::ValidParameters: rejects saveSnapshots without "
    "writeWilsonLines") {
    Parameters param;
    makeValidBaseline(param);
    param.setWriteWilsonLines(0);
    param.setSaveSnapshots(1);
    CHECK(param.ValidParameters() == false);
}

TEST_CASE("Parameters: int-to-bool coercing setters treat any nonzero as true") {
    // setSaveSnapshots/setForceDmin/setComputeGluonMultiplicity/... all
    // share the same "x == 0 -> false, else -> true" pattern; this checks
    // a representative sample rather than every one of them individually.
    Parameters param;

    param.setSaveSnapshots(0);
    CHECK(param.getSaveSnapshots() == false);
    param.setSaveSnapshots(1);
    CHECK(param.getSaveSnapshots() == true);
    param.setSaveSnapshots(2);
    CHECK(param.getSaveSnapshots() == true);
    param.setSaveSnapshots(-1);
    CHECK(param.getSaveSnapshots() == true);

    param.setForceDmin(0);
    CHECK(param.getForceDmin() == false);
    param.setForceDmin(5);
    CHECK(param.getForceDmin() == true);

    param.setComputeGluonMultiplicity(0);
    CHECK(param.getComputeGluonMultiplicity() == false);
    param.setComputeGluonMultiplicity(-3);
    CHECK(param.getComputeGluonMultiplicity() == true);
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

TEST_CASE("Parameters::loadPosteriorParameterSetsFromFile parses a CSV, skipping the header") {
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

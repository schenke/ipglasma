#include "doctest.h"

#include "Parameters.h"

namespace {
// Parameters' default constructor does not initialize every member (no
// default member initializers), so every field ValidParameters() reads
// must be set explicitly before calling it, or the test would depend on
// indeterminate values.
Parameters makeValidBaseline() {
    Parameters param;
    param.setSize(256);
    param.setWriteWilsonLines(2);
    param.setSaveSnapshots(0);
    return param;
}
}  // namespace

TEST_CASE("Parameters::ValidParameters: accepts a normal configuration") {
    Parameters param = makeValidBaseline();
    CHECK(param.ValidParameters() == true);
}

TEST_CASE("Parameters::ValidParameters: rejects a non-positive lattice size") {
    Parameters param = makeValidBaseline();
    param.setSize(0);
    CHECK(param.ValidParameters() == false);
}

TEST_CASE("Parameters::ValidParameters: rejects an invalid Wilson-line format") {
    Parameters param = makeValidBaseline();
    param.setWriteWilsonLines(3);  // only 0 (off), 1 (text), 2 (binary) are valid
    CHECK(param.ValidParameters() == false);
}

TEST_CASE("Parameters::ValidParameters: writeWilsonLines=0 is valid by itself") {
    Parameters param = makeValidBaseline();
    param.setWriteWilsonLines(0);
    CHECK(param.ValidParameters() == true);
}

TEST_CASE(
    "Parameters::ValidParameters: rejects saveSnapshots without "
    "writeWilsonLines") {
    Parameters param = makeValidBaseline();
    param.setWriteWilsonLines(0);
    param.setSaveSnapshots(1);
    CHECK(param.ValidParameters() == false);
}

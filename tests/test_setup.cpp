#include "doctest.h"

#include <cstdio>
#include <fstream>
#include <string>

#include "Setup.h"

namespace {
// A small IP-Glasma-style "key value" input file, written to a fixed path
// in the current working directory (ctest runs the test binary with a
// writable cwd) and cleaned up at the end of each test case.
//
// Note: Setup::stringFind()/iFind() (the non-"Optional" variants) call
// exit(1) on a missing key or missing file -- that would kill the whole
// test binary, not just fail an assertion, so those paths are
// deliberately not exercised here. Only the success path (key present)
// is covered for those two; the Optional variants' "key absent, file
// present" fallback-to-default path is safe to test and is covered below.
class TempInputFile {
  public:
    explicit TempInputFile(const std::string &contents) {
        path_ = "ipglasma_test_setup_tmp_input.txt";
        std::ofstream out(path_);
        out << contents;
    }
    ~TempInputFile() { std::remove(path_.c_str()); }
    const std::string &path() const { return path_; }

  private:
    std::string path_;
};
}  // namespace

TEST_CASE("Setup: stringFind/iFind read present keys") {
    TempInputFile file(
        "foo 42\n"
        "name hello\n"
        "EndOfFile\n");
    Setup setup;

    CHECK(setup.stringFind(file.path(), "name") == "hello");
    CHECK(setup.iFind(file.path(), "foo") == 42);
}

TEST_CASE("Setup: stringFindOptional/iFindOptional fall back to default") {
    TempInputFile file(
        "foo 42\n"
        "name hello\n"
        "EndOfFile\n");
    Setup setup;

    CHECK(setup.stringFindOptional(file.path(), "name", "fallback") == "hello");
    CHECK(
        setup.stringFindOptional(file.path(), "missingKey", "fallback")
        == "fallback");
    CHECK(setup.iFindOptional(file.path(), "foo", -1) == 42);
    CHECK(setup.iFindOptional(file.path(), "missingKey", -1) == -1);
}

TEST_CASE("Setup: isFile reflects whether the path exists") {
    TempInputFile file("EndOfFile\n");
    Setup setup;

    CHECK(setup.isFile(file.path()) == 1);
    CHECK(setup.isFile("ipglasma_test_setup_definitely_missing.txt") == 0);
}

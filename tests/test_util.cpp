#include "doctest.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

#include "Util.h"

namespace {
// Util::stringFind/dFind/iFind (the non-"Optional" variants) call exit(1)
// on a missing key or missing file -- same as Setup's equivalents -- so
// only the success path (key present) is exercised here.
class TempInputFile {
  public:
    explicit TempInputFile(const std::string &contents) {
        path_ = "ipglasma_test_util_tmp_input.txt";
        std::ofstream out(path_);
        out << contents;
    }
    ~TempInputFile() { std::remove(path_.c_str()); }
    const std::string &path() const { return path_; }

  private:
    std::string path_;
};
}  // namespace

TEST_CASE("Util::vector_malloc returns a zero-initialized array") {
    double *v = Util::vector_malloc(10);
    for (int i = 0; i < 10; ++i) CHECK(v[i] == 0.0);
    Util::vector_free(v);
}

TEST_CASE("Util::char_malloc returns an empty string") {
    char *c = Util::char_malloc(16);
    CHECK(std::strlen(c) == 0);
    Util::char_free(c);
}

TEST_CASE("Util::isFile reflects whether the path exists") {
    TempInputFile file("EndOfData\n");
    CHECK(Util::isFile(file.path()) == 1);
    CHECK(Util::isFile("ipglasma_test_util_definitely_missing.txt") == 0);
}

TEST_CASE("Util::stringFind/dFind/iFind read present keys") {
    TempInputFile file(
        "name hello\n"
        "count 7\n"
        "pi 3.14159\n"
        "EndOfData\n");

    CHECK(Util::stringFind(file.path(), "name") == "hello");
    CHECK(Util::iFind(file.path(), "count") == 7);
    CHECK(Util::dFind(file.path(), "pi") == doctest::Approx(3.14159));
}

#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

#include "Util.h"
#include "doctest.h"

namespace {
// Util::stringFind/dFind/iFind (the non-"Optional" variants) exit(1) on an
// empty file name or an absent key, same as Setup's equivalents. A missing
// but non-empty file name instead falls back to reading an empty "input"
// file it creates in the current directory -- which then almost always
// hits the "key not found" exit(1) path too, since that fallback file has
// nothing in it. Only the success path (key present in an existing file)
// is exercised here; the exit(1) paths would need a subprocess-based
// death test, which this suite doesn't have infrastructure for.
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

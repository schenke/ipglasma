#include "doctest.h"

#include <cstdlib>
#include <string>

#include "PrettyOstream.h"

TEST_CASE("PrettyOstream::getMemoryUsage returns a plausible \"<number> MB\"") {
    PrettyOstream messager;
    std::string usage = messager.getMemoryUsage();

    REQUIRE(usage.size() > 3);
    CHECK(usage.substr(usage.size() - 3) == " MB");

    const double value = std::atof(usage.c_str());
    CHECK(value >= 0.0);
    CHECK(value < 1e7);  // sanity bound, not a tight one
}

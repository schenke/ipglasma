// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <cstdlib>
#include <string>

#include "PrettyOstream.h"
#include "doctest.h"

TEST_CASE("PrettyOstream::getMemoryUsage returns a plausible \"<number> MB\"") {
    PrettyOstream messager;
    std::string usage = messager.getMemoryUsage();

    REQUIRE(usage.size() > 3);
    CHECK(usage.substr(usage.size() - 3) == " MB");

    const double value = std::atof(usage.c_str());
    CHECK(value >= 0.0);
    CHECK(value < 1e7);  // sanity bound, not a tight one
}

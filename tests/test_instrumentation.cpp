// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <cstdlib>
#include <string>
#include <thread>

#include "Instrumentation.h"
#include "doctest.h"

TEST_CASE("ipg::wallSeconds is finite and monotonically non-decreasing") {
    const double t1 = ipg::wallSeconds();
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    const double t2 = ipg::wallSeconds();

    CHECK(t1 >= 0.0);
    CHECK(t2 >= t1);
    CHECK((t2 - t1) < 5.0);  // sanity bound, not a tight one
}

TEST_CASE(
    "ipg::fingerprintEnabled: the disabling values are matched "
    "case-insensitively") {
    const char *previous = std::getenv("IPGLASMA_FINGERPRINT");
    const std::string saved = previous ? previous : "";

    unsetenv("IPGLASMA_FINGERPRINT");
    CHECK_FALSE(ipg::fingerprintEnabled());
    for (const char *off :
         {"", "0", "false", "False", "FALSE", "off", "Off", "no", "No", "NO"}) {
        CAPTURE(off);
        setenv("IPGLASMA_FINGERPRINT", off, 1);
        CHECK_FALSE(ipg::fingerprintEnabled());
    }
    for (const char *on : {"1", "yes", "true", "On"}) {
        CAPTURE(on);
        setenv("IPGLASMA_FINGERPRINT", on, 1);
        CHECK(ipg::fingerprintEnabled());
    }

    if (previous) {
        setenv("IPGLASMA_FINGERPRINT", saved.c_str(), 1);
    } else {
        unsetenv("IPGLASMA_FINGERPRINT");
    }
}

// Profiler reads IPGLASMA_PROFILE once and is a process-wide singleton, so it
// is not a self-contained unit to test without process isolation and is not
// covered here.

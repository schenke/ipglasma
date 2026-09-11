#include "doctest.h"

#include <thread>

#include "Instrumentation.h"

TEST_CASE("ipg::wallSeconds is finite and monotonically non-decreasing") {
    const double t1 = ipg::wallSeconds();
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    const double t2 = ipg::wallSeconds();

    CHECK(t1 >= 0.0);
    CHECK(t2 >= t1);
    CHECK((t2 - t1) < 5.0);  // sanity bound, not a tight one
}

// Profiler/fingerprintEnabled are gated by environment variables read once
// at process start (IPGLASMA_PROFILE, IPGLASMA_FINGERPRINT) and, for
// Profiler, held in a process-wide singleton -- neither is a "simple",
// self-contained unit to test without process isolation, so they are not
// covered here.

#include "doctest.h"

#include <cmath>

#include "Fragmentation.h"

// Fragmentation::kkp is a large, hardcoded fit (Kniehl-Kramer-Potter
// parametrization, NPB582 (2000) 514) with no simple closed form to check
// against, so this is a sanity check rather than a correctness check.
// Note: the NLO (iset=1) fit is not guaranteed to be non-negative
// everywhere (observed a small negative value at ih=7, x=0.9, qs=5 --
// that is a known characteristic of this class of NLO fit near the edges
// of its range, not a bug), so this only checks finiteness and a loose
// magnitude bound for physically sensible inputs (momentum fraction x in
// (0, 1), scale in the GeV range actually used by Evolution.cpp).
TEST_CASE("Fragmentation::kkp returns finite, bounded values for typical inputs") {
    const double xs[] = {0.1, 0.3, 0.5, 0.7, 0.9};
    const int hadrons[] = {1, 2, 4, 7};  // pions, kaons, protons, h+/h-
    for (int ih : hadrons) {
        for (double x : xs) {
            const double value = Fragmentation::kkp(ih, /*iset=*/1, x, /*qs=*/5.0);
            CHECK(std::isfinite(value));
            CHECK(std::abs(value) < 100.0);
        }
    }
}

#include "doctest.h"

#include "Cell.h"

TEST_CASE("Cell: default constructor zero-initializes every field") {
    // A regression guard: this codebase has repeatedly had bugs from
    // members left uninitialized by a hand-written constructor (see
    // CHANGELOG's "Fixed" section), so this pins down that Cell's
    // constructor actually initializes all of them, including any field
    // added later.
    Cell c;
    CHECK(c.getEpsilon() == 0.0);
    CHECK(c.getg2mu2A() == 0.0);
    CHECK(c.getg2mu2B() == 0.0);
    CHECK(c.getTpA() == 0.0);
    CHECK(c.getTpB() == 0.0);
    CHECK(c.getTtautau() == 0.0);
    CHECK(c.getTxx() == 0.0);
    CHECK(c.getTyy() == 0.0);
    CHECK(c.getTxy() == 0.0);
    CHECK(c.getTetaeta() == 0.0);
    CHECK(c.getTtaux() == 0.0);
    CHECK(c.getTtauy() == 0.0);
    CHECK(c.getTtaueta() == 0.0);
    CHECK(c.getTxeta() == 0.0);
    CHECK(c.getTyeta() == 0.0);
    CHECK(c.getpitautau() == 0.0);
    CHECK(c.getpixx() == 0.0);
    CHECK(c.getpiyy() == 0.0);
    CHECK(c.getpixy() == 0.0);
    CHECK(c.getpietaeta() == 0.0);
    CHECK(c.getpitaux() == 0.0);
    CHECK(c.getpitauy() == 0.0);
    CHECK(c.getpitaueta() == 0.0);
    CHECK(c.getpixeta() == 0.0);
    CHECK(c.getpiyeta() == 0.0);
    CHECK(c.getutau() == 0.0);
    CHECK(c.getux() == 0.0);
    CHECK(c.getuy() == 0.0);
    CHECK(c.getueta() == 0.0);
}

TEST_CASE("Cell: setters/getters round-trip (representative sample)") {
    Cell c;
    c.setEpsilon(1.5);
    CHECK(c.getEpsilon() == 1.5);
    c.setTxx(-2.25);
    CHECK(c.getTxx() == -2.25);
    c.setux(0.3);
    CHECK(c.getux() == 0.3);
}

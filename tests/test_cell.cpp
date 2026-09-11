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

TEST_CASE("Cell: every setter/getter pair round-trips") {
    Cell c;
    c.setEpsilon(1.5);
    CHECK(c.getEpsilon() == 1.5);
    c.setg2mu2A(2.5);
    CHECK(c.getg2mu2A() == 2.5);
    c.setg2mu2B(3.5);
    CHECK(c.getg2mu2B() == 3.5);
    c.setTpA(4.5);
    CHECK(c.getTpA() == 4.5);
    c.setTpB(5.5);
    CHECK(c.getTpB() == 5.5);

    c.setTtautau(1.1);
    CHECK(c.getTtautau() == 1.1);
    c.setTxx(-2.25);
    CHECK(c.getTxx() == -2.25);
    c.setTyy(2.2);
    CHECK(c.getTyy() == 2.2);
    c.setTxy(3.3);
    CHECK(c.getTxy() == 3.3);
    c.setTetaeta(4.4);
    CHECK(c.getTetaeta() == 4.4);
    c.setTtaux(5.5);
    CHECK(c.getTtaux() == 5.5);
    c.setTtauy(6.6);
    CHECK(c.getTtauy() == 6.6);
    c.setTtaueta(7.7);
    CHECK(c.getTtaueta() == 7.7);
    c.setTxeta(8.8);
    CHECK(c.getTxeta() == 8.8);
    c.setTyeta(9.9);
    CHECK(c.getTyeta() == 9.9);

    c.setpitautau(1.01);
    CHECK(c.getpitautau() == 1.01);
    c.setpixx(1.02);
    CHECK(c.getpixx() == 1.02);
    c.setpiyy(1.03);
    CHECK(c.getpiyy() == 1.03);
    c.setpixy(1.04);
    CHECK(c.getpixy() == 1.04);
    c.setpietaeta(1.05);
    CHECK(c.getpietaeta() == 1.05);
    c.setpitaux(1.06);
    CHECK(c.getpitaux() == 1.06);
    c.setpitauy(1.07);
    CHECK(c.getpitauy() == 1.07);
    c.setpitaueta(1.08);
    CHECK(c.getpitaueta() == 1.08);
    c.setpixeta(1.09);
    CHECK(c.getpixeta() == 1.09);
    c.setpiyeta(1.10);
    CHECK(c.getpiyeta() == 1.10);

    c.setutau(0.1);
    CHECK(c.getutau() == 0.1);
    c.setux(0.3);
    CHECK(c.getux() == 0.3);
    c.setuy(0.4);
    CHECK(c.getuy() == 0.4);
    c.setueta(0.5);
    CHECK(c.getueta() == 0.5);
}

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <complex>
#include <cstdio>
#include <fstream>
#include <string>

#include "Glauber.h"  // for NucleusRole
#include "Lattice.h"
#include "Parameters.h"
#include "doctest.h"
#include "test_helpers.h"

TEST_CASE(
    "Lattice: constructor allocates fields as identity and sizes them "
    "correctly") {
    const int length = 4;
    Parameters param;
    makeLatticeParam(param, length);
    Lattice lat(&param, length);

    CHECK(lat.getSize() == length * length);
    CHECK(lat.U.size() == static_cast<std::size_t>(length * length));
    CHECK(lat.cells.size() == static_cast<std::size_t>(length * length));

    const Matrix identity(1.0);
    for (int pos = 0; pos < length * length; ++pos) {
        for (int k = 0; k < 9; ++k) {
            CHECK(std::abs(lat.U[pos].get(k) - identity.get(k)) < 1e-14);
            CHECK(std::abs(lat.Uy2[pos].get(k) - identity.get(k)) < 1e-14);
        }
    }
}

TEST_CASE(
    "Lattice: neighbor index arrays wrap at the boundary (clamped, not "
    "periodic)") {
    const int length = 4;
    Parameters param;
    makeLatticeParam(param, length);
    Lattice lat(&param, length);

    // Site (0, 0): both minus-neighbors clamp back to the site itself.
    CHECK(lat.posmX[0] == 0);
    CHECK(lat.posmY[0] == 0);
    // Site (0, 0): plus-neighbors move one step in x or y.
    CHECK(lat.pospX[0] == length);  // (1, 0)
    CHECK(lat.pospY[0] == 1);       // (0, 1)

    // Site (length-1, length-1): both plus-neighbors clamp back to itself.
    const int last = length * length - 1;
    CHECK(lat.pospX[last] == last);
    CHECK(lat.pospY[last] == last);
}

TEST_CASE("Lattice::writeSU3Matrices writes non-empty Phi/Pi files") {
    const int length = 4;
    Parameters param;
    makeLatticeParam(param, length);
    Lattice lat(&param, length);

    const std::string prefix = "ipglasma_test_lattice_su3_";
    lat.writeSU3Matrices(prefix, &param);

    // Phi suffix = eventId + 2*seed*MPISize = 0; Pi suffix =
    // eventId + (1 + 2*seed)*MPISize = 1, with eventId=seed=0, MPISize=1.
    for (const char *name : {"Phi-0.txt", "Pi-1.txt"}) {
        const std::string path = prefix + name;
        std::ifstream in(path);
        REQUIRE(in.good());
        std::string firstLine;
        std::getline(in, firstLine);
        CHECK(!firstLine.empty());
        in.close();
        std::remove(path.c_str());
    }
}

TEST_CASE(
    "BufferLattice: allocates buffer1/buffer2 as identity, sized correctly") {
    const int length = 4;
    BufferLattice buf(length);
    CHECK(buf.buffer1.size() == static_cast<std::size_t>(length * length));
    CHECK(buf.buffer2.size() == static_cast<std::size_t>(length * length));

    const Matrix identity(1.0);
    for (int k = 0; k < 9; ++k) {
        CHECK(std::abs(buf.buffer1[0].get(k) - identity.get(k)) < 1e-14);
    }
}

TEST_CASE(
    "Lattice::positionFromXY/xFromPosition/yFromPosition are consistent "
    "with each other, with latticeIndex() and with the neighbor tables") {
    const int length = 5;
    Parameters param;
    makeLatticeParam(param, length);
    Lattice lat(&param, length);
    CHECK(lat.getLength() == length);

    for (int ix = 0; ix < length; ++ix) {
        for (int iy = 0; iy < length; ++iy) {
            const int pos = lat.positionFromXY(ix, iy);
            CHECK(pos == latticeIndex(ix, iy, length));
            CHECK(lat.xFromPosition(pos) == ix);
            CHECK(lat.yFromPosition(pos) == iy);
            CHECK(latticeX(pos, length) == ix);
            CHECK(latticeY(pos, length) == iy);
            // interior sites: the neighbor tables step one site in x / y
            if (ix > 0 && ix < length - 1 && iy > 0 && iy < length - 1) {
                CHECK(lat.pospX[pos] == lat.positionFromXY(ix + 1, iy));
                CHECK(lat.posmX[pos] == lat.positionFromXY(ix - 1, iy));
                CHECK(lat.pospY[pos] == lat.positionFromXY(ix, iy + 1));
                CHECK(lat.posmY[pos] == lat.positionFromXY(ix, iy - 1));
            }
        }
    }
}

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "NuclearQsTable.h"
#include "doctest.h"
#include "test_helpers.h"

namespace {
// The Bjorken x of the tabulated rapidity y.
double xAt(double y) { return 0.01 * std::exp(-y); }
}  // namespace

TEST_CASE("NuclearQsTable interpolates bilinearly inside the table") {
    const std::string fileName = writeLinearQsTable();
    NuclearQsTable table;
    table.read(fileName);
    std::remove(fileName.c_str());

    for (double T : {0.5, 0.73, 3.21, 20.0}) {
        for (double y : {0., 0.1, 4.6, 10.5}) {
            CAPTURE(T);
            CAPTURE(y);
            CHECK(table.qs2(T, xAt(y)) == doctest::Approx(2. * T + 3. * y));
        }
    }
}

TEST_CASE(
    "NuclearQsTable gives 0 below the tabulated T_p and clamps above it") {
    const std::string fileName = writeLinearQsTable();
    NuclearQsTable table;
    table.read(fileName);
    std::remove(fileName.c_str());

    CHECK(table.qs2(0.4, xAt(2.0)) == 0.);
    const double Tmax = 0.5 + 0.1 * 239;
    CHECK(table.qs2(100., xAt(2.1)) == doctest::Approx(2. * Tmax + 3. * 2.1));
}

TEST_CASE("NuclearQsTable stays inside the table at its upper edges") {
    const std::string fileName = writeLinearQsTable();
    NuclearQsTable table;
    table.read(fileName);
    std::remove(fileName.c_str());

    // the largest tabulated T_p (written as "24.4") and rapidity
    const double Tmax = 24.4;
    const double ymax = 10.75;
    for (double T : {24.35, Tmax}) {
        for (double y : {0., 10.6, ymax}) {
            CAPTURE(T);
            CAPTURE(y);
            CHECK(table.qs2(T, xAt(y)) == doctest::Approx(2. * T + 3. * y));
        }
    }
}

TEST_CASE(
    "NuclearQsTable clamps an x below the table (a rapidity above it) to the "
    "smallest tabulated one") {
    const std::string fileName = writeLinearQsTable();
    NuclearQsTable table;
    table.read(fileName);
    std::remove(fileName.c_str());

    for (double T : {0.73, 24.4, 100.}) {
        for (double y : {10.8, 11.0, 13.5}) {
            CAPTURE(T);
            CAPTURE(y);
            CHECK(table.qs2(T, xAt(y)) == table.qs2(T, xAt(10.75)));
        }
    }
    CHECK(table.qs2(0.73, xAt(12.)) == doctest::Approx(2. * 0.73 + 3. * 10.75));
}

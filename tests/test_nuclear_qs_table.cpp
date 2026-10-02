#include <cstdio>
#include <fstream>
#include <string>

#include "NuclearQsTable.h"
#include "doctest.h"

namespace {
// A table with Qs^2 = 2 T + 3 y, which bilinear interpolation reproduces
// exactly: 240 values T = 0.5 + 0.1 iT (outer) times 44 rapidities
// y = 0.25 iy (inner), one "y T Qs^2" line each.
std::string writeLinearTable() {
    const std::string fileName = "ipglasma_test_qs_table.in";
    std::ofstream out(fileName);
    for (int iT = 0; iT < 240; ++iT) {
        for (int iy = 0; iy < 44; ++iy) {
            const double T = 0.5 + 0.1 * iT;
            const double y = 0.25 * iy;
            out << y << " " << T << " " << 2. * T + 3. * y << "\n";
        }
    }
    return fileName;
}
}  // namespace

TEST_CASE("NuclearQsTable interpolates bilinearly inside the table") {
    const std::string fileName = writeLinearTable();
    NuclearQsTable table;
    table.read(fileName);
    std::remove(fileName.c_str());

    for (double T : {0.5, 0.73, 3.21, 20.0}) {
        for (double y : {0., 0.1, 4.6, 10.5}) {
            CAPTURE(T);
            CAPTURE(y);
            CHECK(table.qs2(T, y) == doctest::Approx(2. * T + 3. * y));
        }
    }
}

TEST_CASE(
    "NuclearQsTable gives 0 below the tabulated T_p and clamps above it") {
    const std::string fileName = writeLinearTable();
    NuclearQsTable table;
    table.read(fileName);
    std::remove(fileName.c_str());

    CHECK(table.qs2(0.4, 2.0) == 0.);
    const double Tmax = 0.5 + 0.1 * 239;
    CHECK(table.qs2(100., 2.1) == doctest::Approx(2. * Tmax + 3. * 2.1));
}

TEST_CASE("NuclearQsTable stays inside the table at its upper edges") {
    const std::string fileName = writeLinearTable();
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
            CHECK(table.qs2(T, y) == doctest::Approx(2. * T + 3. * y));
        }
    }
}

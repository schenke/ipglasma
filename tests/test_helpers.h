// Helpers shared by several unit test files.

#ifndef TESTS_TEST_HELPERS_H_
#define TESTS_TEST_HELPERS_H_

#include <complex>
#include <fstream>
#include <string>

#include "Matrix.h"
#include "Parameters.h"

// Writes a nuclear Q_s table with Qs^2 = 2 T + 3 y, which bilinear
// interpolation reproduces exactly: 240 values T = 0.5 + 0.1 iT (outer)
// times 44 rapidities y = 0.25 iy (inner), one "y T Qs^2" line each.
inline std::string writeLinearQsTable() {
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

// Matches test_lattice.cpp's makeLatticeParam: Parameters holds a
// PrettyOstream member, which holds a non-copyable/non-movable
// std::ostringstream, so this takes an out-parameter rather than
// returning by value.
inline void makeInitTestParam(Parameters &param, int size) {
    param.lattice.size = size;
    param.lattice.L = static_cast<double>(size);  // a = 1 fm
    param.run.MPIRank = 0;
    param.event.eventId = 0;
    param.random.seed = 0;
    param.run.MPISize = 1;
    param.colorCharge.projectileX = 0.01;
    param.colorCharge.targetX = 0.01;
}

// A deterministic, position-dependent (not merely diagonal) matrix, so
// tests exercise genuine matrix multiplication/conjugation rather than
// something a transposition or index-swap bug could accidentally pass.
inline Matrix makeTestMatrix(int seed) {
    Matrix m(Matrix::noInit);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            m.set(
                i, j,
                std::complex<double>(
                    0.1 * seed + 0.01 * (i + 1), 0.05 * seed - 0.02 * (j + 1)));
        }
    }
    return m;
}

inline bool matricesClose(const Matrix &a, const Matrix &b, double tol) {
    for (int k = 0; k < 9; ++k) {
        if (std::abs(a.get(k) - b.get(k)) >= tol) return false;
    }
    return true;
}

// Takes an out-parameter rather than returning by value: Parameters holds
// a PrettyOstream member, which holds a non-copyable/non-movable
// std::ostringstream -- see the identical note in test_parameters.cpp.
inline void makeLatticeParam(Parameters &param, int size) {
    param.lattice.L = 10.0;
    param.run.MPIRank = 0;
    param.lattice.size = size;
    param.event.eventId = 0;
    param.random.seed = 0;
    param.run.MPISize = 1;
    param.colorCharge.projectileX = 0.01;
    param.colorCharge.targetX = 0.01;
    param.wilsonLines.wilsonLinePath = ".";
}

#endif  // TESTS_TEST_HELPERS_H_

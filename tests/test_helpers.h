// Helpers shared by several unit test files.

#ifndef TESTS_TEST_HELPERS_H_
#define TESTS_TEST_HELPERS_H_

#include <complex>

#include "Matrix.h"
#include "Parameters.h"

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
    param.colorCharge.rapidityA = 0.0;
    param.colorCharge.rapidityB = 0.0;
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
    param.colorCharge.rapidityA = 0.0;
    param.colorCharge.rapidityB = 0.0;
    param.wilsonLines.wilsonLinePath = ".";
}

#endif  // TESTS_TEST_HELPERS_H_

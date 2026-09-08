#ifndef SU3_h
#define SU3_h

#include <complex>

#include "Matrix.h"

// Fixed-size SU(3) hot-path helpers. These routines intentionally do not
// provide a generic-N fallback: IP-Glasma is validated to run with Nc == 3
// at input time. Keeping the kernels dimension-free lets the compiler see the
// complete 3x3 operation and avoids Matrix temporaries when only a trace is
// required.
namespace su3 {

using Complex = std::complex<double>;

struct Matrix3 {
    Complex e[9];
};

inline Matrix3 multiply(const Matrix &a, const Matrix &b) {
    const Complex *A = a.data();
    const Complex *B = b.data();
    Matrix3 c;
    c.e[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
    c.e[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
    c.e[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
    c.e[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
    c.e[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
    c.e[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
    c.e[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
    c.e[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
    c.e[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
    return c;
}

inline Matrix3 multiplyABdagger(const Matrix &a, const Matrix &b) {
    const Complex *A = a.data();
    const Complex *B = b.data();
    Matrix3 c;
    c.e[0] = A[0] * std::conj(B[0]) + A[1] * std::conj(B[1])
             + A[2] * std::conj(B[2]);
    c.e[1] = A[0] * std::conj(B[3]) + A[1] * std::conj(B[4])
             + A[2] * std::conj(B[5]);
    c.e[2] = A[0] * std::conj(B[6]) + A[1] * std::conj(B[7])
             + A[2] * std::conj(B[8]);
    c.e[3] = A[3] * std::conj(B[0]) + A[4] * std::conj(B[1])
             + A[5] * std::conj(B[2]);
    c.e[4] = A[3] * std::conj(B[3]) + A[4] * std::conj(B[4])
             + A[5] * std::conj(B[5]);
    c.e[5] = A[3] * std::conj(B[6]) + A[4] * std::conj(B[7])
             + A[5] * std::conj(B[8]);
    c.e[6] = A[6] * std::conj(B[0]) + A[7] * std::conj(B[1])
             + A[8] * std::conj(B[2]);
    c.e[7] = A[6] * std::conj(B[3]) + A[7] * std::conj(B[4])
             + A[8] * std::conj(B[5]);
    c.e[8] = A[6] * std::conj(B[6]) + A[7] * std::conj(B[7])
             + A[8] * std::conj(B[8]);
    return c;
}

inline Matrix3 commutator(const Matrix &a, const Matrix &b) {
    Matrix3 ab = multiply(a, b);
    Matrix3 ba = multiply(b, a);
    Matrix3 c;
    for (int i = 0; i < 9; ++i) c.e[i] = ab.e[i] - ba.e[i];
    return c;
}

inline Complex trace(const Matrix &a) {
    const Complex *A = a.data();
    return A[0] + A[4] + A[8];
}

inline Complex trace(const Matrix3 &a) { return a.e[0] + a.e[4] + a.e[8]; }

inline Complex traceAB(const Matrix &a, const Matrix &b) {
    const Complex *A = a.data();
    const Complex *B = b.data();
    const Complex t0 = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
    const Complex t1 = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
    const Complex t2 = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
    return t0 + t1 + t2;
}

inline Complex traceAB(const Matrix3 &a, const Matrix &b) {
    const Complex *B = b.data();
    const Complex t0 = a.e[0] * B[0] + a.e[1] * B[3] + a.e[2] * B[6];
    const Complex t1 = a.e[3] * B[1] + a.e[4] * B[4] + a.e[5] * B[7];
    const Complex t2 = a.e[6] * B[2] + a.e[7] * B[5] + a.e[8] * B[8];
    return t0 + t1 + t2;
}

inline Complex traceAB(const Matrix3 &a, const Matrix3 &b) {
    const Complex t0 = a.e[0] * b.e[0] + a.e[1] * b.e[3] + a.e[2] * b.e[6];
    const Complex t1 = a.e[3] * b.e[1] + a.e[4] * b.e[4] + a.e[5] * b.e[7];
    const Complex t2 = a.e[6] * b.e[2] + a.e[7] * b.e[5] + a.e[8] * b.e[8];
    return t0 + t1 + t2;
}

inline Complex traceABdagger(const Matrix &a, const Matrix &b) {
    const Complex *A = a.data();
    const Complex *B = b.data();
    const Complex t0 = A[0] * std::conj(B[0]) + A[1] * std::conj(B[1])
                       + A[2] * std::conj(B[2]);
    const Complex t1 = A[3] * std::conj(B[3]) + A[4] * std::conj(B[4])
                       + A[5] * std::conj(B[5]);
    const Complex t2 = A[6] * std::conj(B[6]) + A[7] * std::conj(B[7])
                       + A[8] * std::conj(B[8]);
    return t0 + t1 + t2;
}

inline Complex traceSquare(const Matrix &a) { return traceAB(a, a); }

inline Complex traceDifferenceSquare(const Matrix &a, const Matrix &b) {
    const Complex *A = a.data();
    const Complex *B = b.data();
    const Complex d0 = A[0] - B[0];
    const Complex d1 = A[1] - B[1];
    const Complex d2 = A[2] - B[2];
    const Complex d3 = A[3] - B[3];
    const Complex d4 = A[4] - B[4];
    const Complex d5 = A[5] - B[5];
    const Complex d6 = A[6] - B[6];
    const Complex d7 = A[7] - B[7];
    const Complex d8 = A[8] - B[8];
    const Complex t0 = d0 * d0 + d1 * d3 + d2 * d6;
    const Complex t1 = d3 * d1 + d4 * d4 + d5 * d7;
    const Complex t2 = d6 * d2 + d7 * d5 + d8 * d8;
    return t0 + t1 + t2;
}

inline Complex traceABC(const Matrix &a, const Matrix &b, const Matrix &c) {
    const Complex *A = a.data();
    const Complex *B = b.data();
    const Complex *C = c.data();
    Complex tr = 0.;
    for (int i = 0; i < 3; ++i) {
        const int ai = 3 * i;
        for (int j = 0; j < 3; ++j) {
            const Complex ab =
                A[ai] * B[j] + A[ai + 1] * B[3 + j] + A[ai + 2] * B[6 + j];
            tr += ab * C[3 * j + i];
        }
    }
    return tr;
}

inline Complex traceABCD(
    const Matrix &a, const Matrix &b, const Matrix &c, const Matrix &d) {
    const Complex *A = a.data();
    const Complex *B = b.data();
    const Complex *C = c.data();
    const Complex *D = d.data();
    Complex tr = 0.;
    for (int i = 0; i < 3; ++i) {
        const int ai = 3 * i;
        for (int j = 0; j < 3; ++j) {
            const Complex ab =
                A[ai] * B[j] + A[ai + 1] * B[3 + j] + A[ai + 2] * B[6 + j];
            const int cj = 3 * j;
            const Complex cd =
                C[cj] * D[i] + C[cj + 1] * D[3 + i] + C[cj + 2] * D[6 + i];
            tr += ab * cd;
        }
    }
    return tr;
}

}  // namespace su3

#endif

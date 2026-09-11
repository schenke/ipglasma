#include "doctest.h"

#include "Matrix.h"
#include "SU3.h"

namespace {
Matrix fromMatrix3(const su3::Matrix3 &m) {
    Matrix out(Matrix::noInit);
    for (int i = 0; i < 9; ++i) out.set(i, m.e[i]);
    return out;
}

su3::Matrix3 toMatrix3(const Matrix &m) {
    su3::Matrix3 out;
    for (int i = 0; i < 9; ++i) out.e[i] = m.get(i);
    return out;
}

bool closeTo(std::complex<double> a, std::complex<double> b, double tol) {
    return std::abs(a - b) < tol;
}

bool matricesClose(const Matrix &a, const Matrix &b, double tol) {
    for (int i = 0; i < 9; ++i) {
        if (!closeTo(a.get(i), b.get(i), tol)) return false;
    }
    return true;
}

// Two fixed, arbitrary matrices reused across the cross-checks below.
Matrix makeA() {
    Matrix a(Matrix::noInit);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            a.set(i, j, std::complex<double>(0.1 * (i + 1), 0.05 * (j + 1)));
        }
    }
    return a;
}

Matrix makeB() {
    Matrix b(Matrix::noInit);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            b.set(i, j, std::complex<double>(0.2 * (j + 1), -0.03 * (i + 1)));
        }
    }
    return b;
}
}  // namespace

TEST_CASE("su3::multiply matches Matrix's operator*") {
    Matrix a = makeA(), b = makeB();
    Matrix product = fromMatrix3(su3::multiply(a, b));
    CHECK(matricesClose(product, a * b, 1e-13));
}

TEST_CASE("su3::multiplyABdagger matches a * conjg(b)") {
    Matrix a = makeA(), b = makeB();
    Matrix bDagger = b;
    bDagger.conjg();
    Matrix product = fromMatrix3(su3::multiplyABdagger(a, b));
    CHECK(matricesClose(product, a * bDagger, 1e-13));
}

TEST_CASE("su3::commutator matches a*b - b*a") {
    Matrix a = makeA(), b = makeB();
    Matrix comm = fromMatrix3(su3::commutator(a, b));
    CHECK(matricesClose(comm, a * b - b * a, 1e-13));
}

TEST_CASE("su3::trace matches Matrix::trace, for both Matrix and Matrix3") {
    Matrix a = makeA();
    CHECK(closeTo(su3::trace(a), a.trace(), 1e-13));

    su3::Matrix3 a3 = su3::multiply(a, a);
    CHECK(closeTo(su3::trace(a3), fromMatrix3(a3).trace(), 1e-13));
}

TEST_CASE("su3::traceAB overloads match (a*b).trace()") {
    Matrix a = makeA(), b = makeB();
    const std::complex<double> expected = (a * b).trace();

    CHECK(closeTo(su3::traceAB(a, b), expected, 1e-13));

    const su3::Matrix3 a3 = toMatrix3(a);
    CHECK(closeTo(su3::traceAB(a3, b), expected, 1e-13));

    const su3::Matrix3 b3 = toMatrix3(b);
    CHECK(closeTo(su3::traceAB(a3, b3), expected, 1e-13));
}

TEST_CASE("su3::traceABdagger matches (a*conjg(b)).trace()") {
    Matrix a = makeA(), b = makeB();
    Matrix bDagger = b;
    bDagger.conjg();
    CHECK(closeTo(su3::traceABdagger(a, b), (a * bDagger).trace(), 1e-13));
}

TEST_CASE("su3::traceSquare matches (a*a).trace()") {
    Matrix a = makeA();
    CHECK(closeTo(su3::traceSquare(a), (a * a).trace(), 1e-13));
}

TEST_CASE("su3::traceDifferenceSquare matches ((a-b)*(a-b)).trace()") {
    Matrix a = makeA(), b = makeB();
    Matrix diff = a - b;
    CHECK(closeTo(
        su3::traceDifferenceSquare(a, b), (diff * diff).trace(), 1e-13));
}

TEST_CASE("su3::traceABC matches (a*b*c).trace()") {
    Matrix a = makeA(), b = makeB(), c = makeA() + makeB();
    CHECK(closeTo(su3::traceABC(a, b, c), (a * b * c).trace(), 1e-12));
}

TEST_CASE("su3::traceABCD matches (a*b*c*d).trace()") {
    Matrix a = makeA(), b = makeB(), c = makeA() + makeB(), d = makeB() - makeA();
    CHECK(closeTo(su3::traceABCD(a, b, c, d), (a * b * c * d).trace(), 1e-12));
}

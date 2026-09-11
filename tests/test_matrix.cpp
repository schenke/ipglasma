#include "doctest.h"

#include "Matrix.h"

namespace {
bool closeTo(std::complex<double> a, std::complex<double> b, double tol) {
    return std::abs(a - b) < tol;
}

bool matricesClose(Matrix a, Matrix b, double tol) {
    for (int i = 0; i < a.getNDim(); ++i) {
        for (int j = 0; j < a.getNDim(); ++j) {
            if (!closeTo(a.get(i, j), b.get(i, j), tol)) return false;
        }
    }
    return true;
}
}  // namespace

TEST_CASE("Matrix: identity properties") {
    Matrix I(1.);
    CHECK(closeTo(I.trace(), 3.0, 1e-14));
    CHECK(closeTo(I.det(), 1.0, 1e-14));
    CHECK(I.oneNorm() == doctest::Approx(1.0));
    CHECK(I.frobeniusNorm() == doctest::Approx(std::sqrt(3.0)));
}

TEST_CASE("Matrix: default constructor is the zero matrix") {
    Matrix Z;
    CHECK(closeTo(Z.trace(), 0.0, 1e-14));
    CHECK(Z.frobeniusNorm() == doctest::Approx(0.0));
}

TEST_CASE("Matrix: +=/-=/*=//= operators") {
    Matrix I(1.);
    Matrix M = I;
    M += I;
    CHECK(matricesClose(M, Matrix(2.), 1e-14));
    M -= I;
    CHECK(matricesClose(M, I, 1e-14));
    M *= std::complex<double>(2.0, 0.0);
    CHECK(matricesClose(M, Matrix(2.), 1e-14));
    M /= std::complex<double>(2.0, 0.0);
    CHECK(matricesClose(M, I, 1e-14));
}

TEST_CASE("Matrix: inv() inverts a concrete invertible matrix") {
    // A fixed, arbitrary invertible 3x3 matrix (det != 0).
    Matrix M(Matrix::noInit);
    M.set(0, 0, std::complex<double>(2.0, 0.0));
    M.set(0, 1, std::complex<double>(0.0, 1.0));
    M.set(0, 2, std::complex<double>(0.0, 0.0));
    M.set(1, 0, std::complex<double>(0.0, -1.0));
    M.set(1, 1, std::complex<double>(3.0, 0.0));
    M.set(1, 2, std::complex<double>(1.0, 0.0));
    M.set(2, 0, std::complex<double>(1.0, 0.0));
    M.set(2, 1, std::complex<double>(0.0, 0.0));
    M.set(2, 2, std::complex<double>(4.0, 0.0));

    Matrix Minv = M;
    Minv.inv();

    Matrix product = M * Minv;
    CHECK(matricesClose(product, Matrix(1.), 1e-10));
}

TEST_CASE("Matrix: logm(expm(A)) recovers A for a small generator") {
    // A small, fixed, non-trivial matrix -- well within the convergence
    // radius of both the Pade exponential (expm) and the inverse
    // scaling-and-squaring logarithm (logm), so the round trip should be
    // accurate to close to machine precision.
    Matrix A(Matrix::noInit);
    A.set(0, 0, std::complex<double>(0.0, 0.05));
    A.set(0, 1, std::complex<double>(0.03, -0.02));
    A.set(0, 2, std::complex<double>(0.0, 0.0));
    A.set(1, 0, std::complex<double>(-0.03, -0.02));
    A.set(1, 1, std::complex<double>(0.0, -0.05));
    A.set(1, 2, std::complex<double>(0.01, 0.0));
    A.set(2, 0, std::complex<double>(0.0, 0.0));
    A.set(2, 1, std::complex<double>(-0.01, 0.0));
    A.set(2, 2, std::complex<double>(0.0, 0.0));

    Matrix X = A;
    X.expm();  // X = exp(A)
    X.logm();  // X = log(exp(A)), should recover A

    CHECK(matricesClose(X, A, 1e-8));
}

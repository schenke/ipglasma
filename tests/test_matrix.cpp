// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <algorithm>
#include <cmath>
#include <complex>
#include <random>
#include <vector>

#include "Group.h"
#include "Matrix.h"
#include "doctest.h"

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

TEST_CASE("Matrix: logmPade computes log(I+A) directly for small A") {
    Matrix A(Matrix::noInit);
    A.set(0, 0, std::complex<double>(0.01, 0.0));
    A.set(0, 1, std::complex<double>(0.0, 0.005));
    A.set(0, 2, std::complex<double>(0.0, 0.0));
    A.set(1, 0, std::complex<double>(0.0, -0.005));
    A.set(1, 1, std::complex<double>(-0.01, 0.0));
    A.set(1, 2, std::complex<double>(0.002, 0.0));
    A.set(2, 0, std::complex<double>(0.0, 0.0));
    A.set(2, 1, std::complex<double>(0.002, 0.0));
    A.set(2, 2, std::complex<double>(0.0, 0.0));

    Matrix L = A;
    L.logmPade(20);  // L should now hold log(I + A)

    Matrix reconstructed = L;
    reconstructed.expm();  // exp(log(I + A)) should recover I + A

    CHECK(matricesClose(reconstructed, Matrix(1.) + A, 1e-10));
}

TEST_CASE("Matrix: sqrtm computes the principal square root") {
    Matrix I(1.);
    Matrix sqrtI = I;
    sqrtI.sqrtm();
    CHECK(matricesClose(sqrtI, I, 1e-10));

    // A concrete diagonal positive-definite matrix: sqrt is unambiguous.
    Matrix M(Matrix::noInit);
    M.set(0, 0, std::complex<double>(4.0, 0.0));
    M.set(0, 1, std::complex<double>(0.0, 0.0));
    M.set(0, 2, std::complex<double>(0.0, 0.0));
    M.set(1, 0, std::complex<double>(0.0, 0.0));
    M.set(1, 1, std::complex<double>(9.0, 0.0));
    M.set(1, 2, std::complex<double>(0.0, 0.0));
    M.set(2, 0, std::complex<double>(0.0, 0.0));
    M.set(2, 1, std::complex<double>(0.0, 0.0));
    M.set(2, 2, std::complex<double>(16.0, 0.0));

    Matrix S = M;
    S.sqrtm();

    Matrix expected(Matrix::noInit);
    expected.set(0, 0, std::complex<double>(2.0, 0.0));
    expected.set(0, 1, std::complex<double>(0.0, 0.0));
    expected.set(0, 2, std::complex<double>(0.0, 0.0));
    expected.set(1, 0, std::complex<double>(0.0, 0.0));
    expected.set(1, 1, std::complex<double>(3.0, 0.0));
    expected.set(1, 2, std::complex<double>(0.0, 0.0));
    expected.set(2, 0, std::complex<double>(0.0, 0.0));
    expected.set(2, 1, std::complex<double>(0.0, 0.0));
    expected.set(2, 2, std::complex<double>(4.0, 0.0));

    CHECK(matricesClose(S, expected, 1e-8));
    CHECK(matricesClose(S * S, M, 1e-8));
}

TEST_CASE("Matrix: setRe/setIm preserve the other component") {
    Matrix M(Matrix::noInit);
    M.set(0, 1, std::complex<double>(1.0, 2.0));

    M.setRe(0, 1, 5.0);
    CHECK(closeTo(M.get(0, 1), std::complex<double>(5.0, 2.0), 1e-14));
    CHECK(M.getRe(1) == doctest::Approx(5.0));  // flat index of (0, 1) is 1
    CHECK(M.getIm(1) == doctest::Approx(2.0));

    M.setIm(0, 1, 7.0);
    CHECK(closeTo(M.get(0, 1), std::complex<double>(5.0, 7.0), 1e-14));

    // Flat-index overloads address the same underlying element as (i, j).
    M.setRe(1, 9.0);
    CHECK(M.getRe(1) == doctest::Approx(9.0));
    M.setIm(1, 11.0);
    CHECK(M.getIm(1) == doctest::Approx(11.0));
}

TEST_CASE("Matrix: square() is half the squared Frobenius norm") {
    Matrix I(1.);
    CHECK(I.square() == doctest::Approx(1.5));  // 0.5 * 3

    Matrix M(2.);
    CHECK(M.square() == doctest::Approx(6.0));  // 0.5 * (3 * 2^2)
}

TEST_CASE("Matrix: conjg() computes the conjugate transpose") {
    Matrix M(Matrix::noInit);
    M.set(0, 0, std::complex<double>(1.0, 0.0));
    M.set(0, 1, std::complex<double>(2.0, 3.0));
    M.set(0, 2, std::complex<double>(0.0, 0.0));
    M.set(1, 0, std::complex<double>(2.0, -3.0));
    M.set(1, 1, std::complex<double>(4.0, 0.0));
    M.set(1, 2, std::complex<double>(0.0, 1.0));
    M.set(2, 0, std::complex<double>(0.0, 0.0));
    M.set(2, 1, std::complex<double>(0.0, -1.0));
    M.set(2, 2, std::complex<double>(5.0, 0.0));

    // M as built above is already Hermitian: conjg() should be a no-op.
    Matrix hermitianCopy = M;
    hermitianCopy.conjg();
    CHECK(matricesClose(hermitianCopy, M, 1e-14));

    // For a non-Hermitian matrix, conjg() must match the manual transpose.
    Matrix N(Matrix::noInit);
    N.set(0, 0, std::complex<double>(1.0, 1.0));
    N.set(0, 1, std::complex<double>(2.0, 0.0));
    N.set(0, 2, std::complex<double>(0.0, 0.0));
    N.set(1, 0, std::complex<double>(0.0, 0.0));
    N.set(1, 1, std::complex<double>(3.0, -1.0));
    N.set(1, 2, std::complex<double>(1.0, 0.0));
    N.set(2, 0, std::complex<double>(0.0, 2.0));
    N.set(2, 1, std::complex<double>(0.0, 0.0));
    N.set(2, 2, std::complex<double>(1.0, 0.0));

    Matrix adjoint = N;
    adjoint.conjg();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            CHECK(closeTo(adjoint.get(i, j), std::conj(N.get(j, i)), 1e-14));
        }
    }
}

TEST_CASE("Matrix: prodABconj/prodAconjB match a*conjg(b)/conjg(a)*b") {
    Matrix a(Matrix::noInit), b(Matrix::noInit);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            a.set(i, j, std::complex<double>(0.1 * (i + 1), 0.05 * (j + 1)));
            b.set(i, j, std::complex<double>(0.2 * (j + 1), -0.03 * (i + 1)));
        }
    }

    Matrix bDagger = b;
    bDagger.conjg();
    Matrix expectedABconj = a * bDagger;

    Matrix actualABconj(Matrix::noInit);
    actualABconj = actualABconj.prodABconj(a, b);
    CHECK(matricesClose(actualABconj, expectedABconj, 1e-13));

    Matrix aDagger = a;
    aDagger.conjg();
    Matrix expectedAconjB = aDagger * b;

    Matrix actualAconjB(Matrix::noInit);
    actualAconjB = actualAconjB.prodAconjB(a, b);
    CHECK(matricesClose(actualAconjB, expectedAconjB, 1e-13));
}

TEST_CASE("Matrix: traceOfProductOfMatrix matches (a*b).trace()") {
    Matrix a(Matrix::noInit), b(Matrix::noInit);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            a.set(i, j, std::complex<double>(0.3 * (i + 1), -0.1 * j));
            b.set(i, j, std::complex<double>(-0.2 * j, 0.15 * (i + 1)));
        }
    }
    Matrix dummy;
    CHECK(closeTo(dummy.traceOfProductOfMatrix(a, b), (a * b).trace(), 1e-13));
}

TEST_CASE(
    "Matrix::fromAlgebraExponent matches expm(i*sum Q^a t^a) with Group's "
    "generators") {
    Group group;
    const std::vector<double> Q = {0.10,  -0.07, 0.15,  0.02,
                                   -0.12, 0.08,  -0.05, 0.09};

    Matrix H;  // H = sum_a Q[a] * t_a (traceless Hermitian)
    for (int a = 0; a < 8; ++a) H += Q[a] * group.getT(a);

    Matrix reference = std::complex<double>(0.0, 1.0) * H;
    reference.expm();  // exp(i H)

    CHECK(matricesClose(Matrix::fromAlgebraExponent(Q), reference, 1e-9));
    // the same matrix as expiHermitian of H
    CHECK(matricesClose(
        Matrix::fromAlgebraExponent(Q), Matrix::expiHermitian(H), 1e-15));
}

TEST_CASE("Matrix: MatrixToString formats elements in row-major order") {
    // A non-symmetric matrix with a distinct value at every entry, so a
    // transposed (column-major) order would not go undetected.
    Matrix M(Matrix::noInit);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            M.set(i, j, std::complex<double>(10 * i + j, 100 + 10 * i + j));
        }
    }
    // Row-major: (0,0),(0,1),(0,2), (1,0),(1,1),(1,2), (2,0),(2,1),(2,2).
    CHECK(
        M.MatrixToString()
        == "0 100 1 101 2 102 10 110 11 111 12 112 20 120 21 121 22 122");
}

TEST_CASE("Matrix: getNDim/getNN report the fixed 3x3 SU(3) shape") {
    Matrix M;
    CHECK(M.getNDim() == 3);
    CHECK(M.getNN() == 9);
}

TEST_CASE("Matrix: operator() matches get()") {
    Matrix M(Matrix::noInit);
    M.set(1, 2, std::complex<double>(0.7, -0.3));
    CHECK(closeTo(M(1, 2), M.get(1, 2), 1e-14));
    CHECK(closeTo(M(5), M.get(5), 1e-14));  // flat-index overload
}

TEST_CASE("Matrix: free +/- operators match the in-place versions") {
    Matrix a(2.), b(1.);
    Matrix sum = a + b;
    Matrix diff = a - b;
    CHECK(matricesClose(sum, Matrix(3.), 1e-14));
    CHECK(matricesClose(diff, Matrix(1.), 1e-14));
}

TEST_CASE(
    "Matrix::fromAlgebraExponent: exp(i Q^a t^a) is unitary, also for "
    "degenerate eigenvalues, and Q = 0 gives the identity") {
    // a generic exponent: unitary with determinant 1
    std::vector<double> Q = {0.3, -1.2, 0.7, 2.1, -0.4, 0.9, 1.5, -2.2};
    const Matrix U = Matrix::fromAlgebraExponent(Q);
    Matrix Udag = U;
    Udag.conjg();
    const Matrix product = U * Udag;
    const Matrix identity(1.);
    for (int k = 0; k < 9; k++) {
        CHECK(std::abs(product.get(k) - identity.get(k)) < 1e-12);
    }
    CHECK(std::abs(U.det() - std::complex<double>(1., 0.)) < 1e-12);

    // Q = 0 is exactly the identity
    const Matrix vacuum =
        Matrix::fromAlgebraExponent(std::vector<double>(8, 0.));
    for (int k = 0; k < 9; k++) CHECK(vacuum.get(k) == identity.get(k));
}

namespace {
// exp(iQ) from its power series, with scaling and squaring: the reference
// for Matrix::expiHermitian(), independent of its closed form.
Matrix seriesExpi(const Matrix &Q) {
    double norm = 0.;
    for (int i = 0; i < 9; i++) norm += std::norm(Q(i));
    int squarings = 0;
    while (std::sqrt(norm) / std::pow(2., squarings) > 0.25) squarings++;
    const Matrix X = complex<double>(0., std::pow(0.5, squarings)) * Q;
    Matrix result(1.), term(1.);
    for (int k = 1; k <= 20; k++) {
        term = (1. / k) * (term * X);
        result += term;
    }
    for (int i = 0; i < squarings; i++) result = result * result;
    return result;
}

double maxDifference(const Matrix &a, const Matrix &b) {
    double largest = 0.;
    for (int i = 0; i < 9; i++)
        largest = std::max(largest, std::abs(a(i) - b(i)));
    return largest;
}

// max |U^dagger U - 1| and |det U - 1|
double distanceFromSU3(const Matrix &U) {
    const Matrix product = Matrix::prodAconjB(U, U);
    return std::max(maxDifference(product, Matrix(1.)), std::abs(U.det() - 1.));
}

// A random traceless Hermitian matrix with Frobenius norm `norm`.
Matrix randomTracelessHermitian(std::mt19937 &engine, double norm) {
    std::normal_distribution<double> gauss;
    Matrix A;
    for (int i = 0; i < 9; i++)
        A.set(i, complex<double>(gauss(engine), gauss(engine)));
    Matrix Ad = A;
    Ad.conjg();
    Matrix Q = 0.5 * (A + Ad);
    const complex<double> traceThird = Q.trace() / 3.;
    Q = Q - traceThird * Matrix(1.);
    double frobenius = 0.;
    for (int i = 0; i < 9; i++) frobenius += std::norm(Q(i));
    return (norm / std::sqrt(frobenius)) * Q;
}
}  // namespace

TEST_CASE(
    "Matrix::expiHermitian gives exp(iQ) in SU(3), also for equal "
    "eigenvalues and tiny Q") {
    std::mt19937 engine(7);
    std::vector<Matrix> cases;
    for (const double norm : {1e-13, 1e-7, 1e-3, 0.1, 0.5, 2., 5.}) {
        for (int i = 0; i < 20; i++) {
            cases.push_back(randomTracelessHermitian(engine, norm));
        }
    }
    // equal eigenvalues (lambda_8 and -lambda_8: theta = 0 and pi in the
    // closed form), det Q = 0 (lambda_3), nearly equal eigenvalues
    for (const double s : {1e-6, 0.3, 2.}) {
        Matrix lambda8(0.);
        lambda8.set(0, 0, s);
        lambda8.set(1, 1, s);
        lambda8.set(2, 2, -2. * s);
        cases.push_back(lambda8);
        cases.push_back(-1. * lambda8);
        cases.push_back(lambda8 + randomTracelessHermitian(engine, 1e-9 * s));
        Matrix lambda3(0.);
        lambda3.set(0, 0, s);
        lambda3.set(1, 1, -s);
        cases.push_back(lambda3);
    }
    for (const Matrix &Q : cases) {
        const Matrix U = Matrix::expiHermitian(Q);
        CAPTURE(Q(0));
        CAPTURE(Q(1));
        CHECK(maxDifference(U, seriesExpi(Q)) < 1e-12);
        CHECK(distanceFromSU3(U) < 1e-13);
    }
    // no field: exactly the identity
    CHECK(maxDifference(Matrix::expiHermitian(Matrix(0.)), Matrix(1.)) == 0.);
}

TEST_CASE(
    "Matrix::fromAlgebraExponent is exact along +-t^8 and for nearly "
    "degenerate eigenvalues (issue #100)") {
    // the former closed form returned the identity for Q = +-s t^8 (two
    // equal eigenvalues) and lost precision as 1/splitting^2 close to them
    Group group;
    auto error = [&](const std::vector<double> &Q) {
        Matrix H;
        for (int a = 0; a < 8; ++a) H += Q[a] * group.getT(a);
        return maxDifference(Matrix::fromAlgebraExponent(Q), seriesExpi(H));
    };
    for (const double s : {0.5, 1., -0.5, -1.}) {
        CAPTURE(s);
        std::vector<double> Q(8, 0.);
        Q[7] = s;
        CHECK(error(Q) < 1e-13);
        CHECK(maxDifference(Matrix::fromAlgebraExponent(Q), Matrix(1.)) > 0.1);
    }
    for (const double splitting : {1e-1, 1e-3, 1e-5, 1e-7}) {
        CAPTURE(splitting);
        std::vector<double> Q(8, 0.);
        Q[7] = 1.;
        Q[2] = splitting;
        CHECK(error(Q) < 1e-13);
    }
}

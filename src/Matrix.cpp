// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include "Matrix.h"

constexpr Matrix::NoInitTag Matrix::noInit;

#include <gsl/gsl_integration.h>  // include gsl for Gauss-Legendre nodes and weights for log Pade

#include <algorithm>
#include <sstream>
#include <vector>

#include "PhysConst.h"
#include "PrettyOstream.h"

using PhysConst::Nc2m1;
using PhysConst::smallEps;

static_assert(
    sizeof(Matrix) == 9 * sizeof(std::complex<double>),
    "Matrix must remain an exact contiguous complex<double>[9]");

namespace {
/**
 * Analytic inverse of a general \f$3\times3\f$ matrix via its cofactor
 * (adjugate) matrix. Shared by Matrix::inv() (inverting `*this`) and
 * Matrix::expm() (inverting the Pade denominator \f$Q\f$), which both
 * need this exact formula.
 * \param[in] Q Matrix to invert.
 * \return \f$Q^{-1}\f$.
 */
Matrix cofactorInverse(const Matrix &Q) {
    Matrix H2(Matrix::noInit);
    H2.set(0, 0, (Q(1, 1) * Q(2, 2) - Q(1, 2) * Q(2, 1)));
    H2.set(0, 1, (Q(0, 2) * Q(2, 1) - Q(0, 1) * Q(2, 2)));
    H2.set(0, 2, (Q(0, 1) * Q(1, 2) - Q(0, 2) * Q(1, 1)));
    H2.set(1, 0, (Q(1, 2) * Q(2, 0) - Q(1, 0) * Q(2, 2)));
    H2.set(1, 1, (Q(0, 0) * Q(2, 2) - Q(0, 2) * Q(2, 0)));
    H2.set(1, 2, (Q(0, 2) * Q(1, 0) - Q(0, 0) * Q(1, 2)));
    H2.set(2, 0, (Q(1, 0) * Q(2, 1) - Q(1, 1) * Q(2, 0)));
    H2.set(2, 1, (Q(0, 1) * Q(2, 0) - Q(0, 0) * Q(2, 1)));
    H2.set(2, 2, (Q(0, 0) * Q(1, 1) - Q(0, 1) * Q(1, 0)));
    H2 *= 1. / Q.det();
    return H2;
}
}  // namespace

Matrix::Matrix() {
    for (int i = 0; i < 9; ++i) e_[i] = complex<double>(0.0, 0.0);
}

Matrix::Matrix(double a) {
    for (int i = 0; i < 9; ++i) e_[i] = complex<double>(0.0, 0.0);
    e_[0] = complex<double>(a, 0.0);
    e_[4] = complex<double>(a, 0.0);
    e_[8] = complex<double>(a, 0.0);
}

Matrix::Matrix(NoInitTag) {}

// operators:

Matrix operator*(const Matrix &a, const Matrix &b) {
    Matrix c(Matrix::noInit);
    const complex<double> *A = a.data();
    const complex<double> *B = b.data();
    complex<double> *C = c.data();
    C[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
    C[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
    C[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
    C[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
    C[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
    C[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
    C[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
    C[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
    C[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
    return c;
}

//-
Matrix operator-(const Matrix &a, const Matrix &b) {
    Matrix aa(Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) aa.set(i, a(i) - b(i));
    return aa;
}

//+
Matrix operator+(const Matrix &a, const Matrix &b) {
    Matrix aa(Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) aa.set(i, a(i) + b(i));
    return aa;
}

//* multiply by a real scalar
Matrix operator*(const Matrix &a, const double s) {
    Matrix aa(Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) {
        aa.set(i, a(i) * s);
    }
    return aa;
}
Matrix operator*(const double s, const Matrix &a) {
    Matrix aa(Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) {
        aa.set(i, a(i) * s);
    }
    return aa;
}

//* multiply by a complex number
Matrix operator*(const complex<double> s, const Matrix &a) {
    Matrix aa(Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) {
        aa.set(i, a(i) * s);
    }
    return aa;
}

// Division by a real scalar. Not declared in Matrix.h and never called
// from within this file either, so this is currently unreachable from
// anywhere else in the codebase (a free function needs a visible
// declaration to be called from another translation unit) -- found
// while documenting this file.
Matrix operator/(const Matrix &a, const double s) {
    Matrix aa(Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) aa.set(i, a(i) / s);
    return aa;
}

Matrix &Matrix::conjg() {
    const complex<double> a01 = e_[1];
    const complex<double> a02 = e_[2];
    const complex<double> a12 = e_[5];
    e_[0] = conj(e_[0]);
    e_[4] = conj(e_[4]);
    e_[8] = conj(e_[8]);
    e_[1] = conj(e_[3]);
    e_[2] = conj(e_[6]);
    e_[5] = conj(e_[7]);
    e_[3] = conj(a01);
    e_[6] = conj(a02);
    e_[7] = conj(a12);
    return *this;
}

Matrix Matrix::prodABconj(const Matrix &a, const Matrix &b) {
    Matrix c(Matrix::noInit);
    c.set(
        0, 0,
        a(0, 0) * conj(b(0, 0)) + a(0, 1) * conj(b(0, 1))
            + a(0, 2) * conj(b(0, 2)));
    c.set(
        0, 1,
        a(0, 0) * conj(b(1, 0)) + a(0, 1) * conj(b(1, 1))
            + a(0, 2) * conj(b(1, 2)));
    c.set(
        0, 2,
        a(0, 0) * conj(b(2, 0)) + a(0, 1) * conj(b(2, 1))
            + a(0, 2) * conj(b(2, 2)));
    c.set(
        1, 0,
        a(1, 0) * conj(b(0, 0)) + a(1, 1) * conj(b(0, 1))
            + a(1, 2) * conj(b(0, 2)));
    c.set(
        1, 1,
        a(1, 0) * conj(b(1, 0)) + a(1, 1) * conj(b(1, 1))
            + a(1, 2) * conj(b(1, 2)));
    c.set(
        1, 2,
        a(1, 0) * conj(b(2, 0)) + a(1, 1) * conj(b(2, 1))
            + a(1, 2) * conj(b(2, 2)));
    c.set(
        2, 0,
        a(2, 0) * conj(b(0, 0)) + a(2, 1) * conj(b(0, 1))
            + a(2, 2) * conj(b(0, 2)));
    c.set(
        2, 1,
        a(2, 0) * conj(b(1, 0)) + a(2, 1) * conj(b(1, 1))
            + a(2, 2) * conj(b(1, 2)));
    c.set(
        2, 2,
        a(2, 0) * conj(b(2, 0)) + a(2, 1) * conj(b(2, 1))
            + a(2, 2) * conj(b(2, 2)));
    return c;
}

Matrix Matrix::prodAconjB(const Matrix &a, const Matrix &b) {
    Matrix c(Matrix::noInit);
    c.set(
        0, 0,
        conj(a(0, 0)) * b(0, 0) + conj(a(1, 0)) * b(1, 0)
            + conj(a(2, 0)) * b(2, 0));
    c.set(
        0, 1,
        conj(a(0, 0)) * b(0, 1) + conj(a(1, 0)) * b(1, 1)
            + conj(a(2, 0)) * b(2, 1));
    c.set(
        0, 2,
        conj(a(0, 0)) * b(0, 2) + conj(a(1, 0)) * b(1, 2)
            + conj(a(2, 0)) * b(2, 2));
    c.set(
        1, 0,
        conj(a(0, 1)) * b(0, 0) + conj(a(1, 1)) * b(1, 0)
            + conj(a(2, 1)) * b(2, 0));
    c.set(
        1, 1,
        conj(a(0, 1)) * b(0, 1) + conj(a(1, 1)) * b(1, 1)
            + conj(a(2, 1)) * b(2, 1));
    c.set(
        1, 2,
        conj(a(0, 1)) * b(0, 2) + conj(a(1, 1)) * b(1, 2)
            + conj(a(2, 1)) * b(2, 2));
    c.set(
        2, 0,
        conj(a(0, 2)) * b(0, 0) + conj(a(1, 2)) * b(1, 0)
            + conj(a(2, 2)) * b(2, 0));
    c.set(
        2, 1,
        conj(a(0, 2)) * b(0, 1) + conj(a(1, 2)) * b(1, 1)
            + conj(a(2, 2)) * b(2, 1));
    c.set(
        2, 2,
        conj(a(0, 2)) * b(0, 2) + conj(a(1, 2)) * b(1, 2)
            + conj(a(2, 2)) * b(2, 2));
    return c;
}

// matrix exponential e^iQ of a traceless Hermitian 3x3 matrix Q, as
// f0 + f1 Q + f2 Q^2 (Morningstar and Peardon)
Matrix Matrix::expiHermitian(const Matrix &Q) {
    const Matrix Q2 = Q * Q;
    // c1 = tr(Q^2)/2, c0 = det Q (real for Hermitian Q)
    double c1 = 0.;
    for (int i = 0; i < 9; i++) c1 += std::norm(Q.e_[i]);
    c1 *= 0.5;
    const complex<double> I(0., 1.);
    if (c1 < 1e-24) {
        // |Q| < 1e-12: the next term, Q^3/6, is below 1e-36
        return Matrix(1.) + I * Q - 0.5 * Q2;
    }
    const double c0 = Q.det().real();
    const double c1Over3 = c1 / 3.;
    const double c0max = 2. * c1Over3 * std::sqrt(c1Over3);
    // evaluate for |c0|: theta in [0, pi/2], so 9u^2 - w^2 >= 2 c1
    const double ratio = std::min(std::abs(c0) / c0max, 1.);
    const double thetaOverThree = std::acos(ratio) / 3.;
    const double u = std::sqrt(c1Over3) * std::cos(thetaOverThree);
    const double w = std::sqrt(c1) * std::sin(thetaOverThree);
    const double w2 = w * w;
    // xi0 = sin(w)/w, expanded for small w
    const double xi0 = (w <= 0.05)
                           ? 1. - w2 / 6. * (1. - w2 / 20. * (1. - w2 / 42.))
                           : std::sin(w) / w;
    const double cosw = std::cos(w);
    const complex<double> exp2iu(std::cos(2. * u), std::sin(2. * u));
    const complex<double> expmiu(std::cos(u), -std::sin(u));
    const double u2 = u * u;
    const double den = 9. * u2 - w2;
    complex<double> f0 =
        ((u2 - w2) * exp2iu
         + expmiu * (8. * u2 * cosw + 2. * I * u * (3. * u2 + w2) * xi0))
        / den;
    complex<double> f1 =
        (2. * u * exp2iu - expmiu * (2. * u * cosw - I * (3. * u2 - w2) * xi0))
        / den;
    complex<double> f2 = (exp2iu - expmiu * (cosw + 3. * I * u * xi0)) / den;
    if (c0 < 0.) {
        // f_j(-c0) = (-1)^j f_j(c0)^*
        f0 = std::conj(f0);
        f1 = -std::conj(f1);
        f2 = std::conj(f2);
    }
    Matrix result(noInit);
    for (int i = 0; i < 9; i++) result.e_[i] = f1 * Q.e_[i] + f2 * Q2.e_[i];
    result.e_[0] += f0;
    result.e_[4] += f0;
    result.e_[8] += f0;
    return result;
}

// matrix exponential using Pade approximant
// t is a scalar that multiplies the matrix (default: t=1) and p is the order in
// the Pade approximant (default: p=6)
Matrix &Matrix::expm(double t, const int p) {
    const int n = this->getNDim();
    const Matrix I(1.);
    Matrix U, H2, P, Q;
    double norm = 0.0;
    // Calculate Pade coefficients
    if (p < 6) {
        PrettyOstream messager;
        messager.error("[Matrix::expm]: p should be at least 6. Exiting.");
        exit(1);
    }
    // hard coded values for speed
    std::vector<double> c(p + 1, 0);
    c[0] = 1.;
    c[1] = 0.5;
    c[2] = 0.1136363636;
    c[3] = 0.01515151515;
    c[4] = 0.001262626263;
    c[5] = 6.313131313e-05;
    c[6] = 1.503126503e-06;
    if (p > 6) {
        for (int i = 6; i < p; ++i) {
            c[i + 1] = c[i] * ((p - i) / ((i + 1.0) * (2.0 * p - i)));
        }
    }
    // Calculate the infinty norm of e, which is defined as the largest row sum
    // of a matrix
    for (int i = 0; i < n; ++i) {
        double temp = 0.0;
        for (int j = 0; j < n; j++) temp += abs((*this)(i, j));
        norm = t * std::max<double>(norm, temp);
    }
    // If norm = 0, and all H elements are not nan or infinity but zero,
    // then U should be identity.
    if (norm == 0.0) {
        bool all_H_are_zero = 1;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((*this)(i, j) != 0.0) {
                    all_H_are_zero = 0;
                }
            }
        }
        if (all_H_are_zero) {
            *this = I;
            return *this;
        } else {
            //	    Some error happens, H has elements which are NaN or
            // infinity.
            std::ostringstream errorMsg;
            errorMsg << "[Matrix::expm]: input matrix has a zero norm but "
                        "non-zero elements (likely NaN or Inf). Input: "
                     << *this;
            PrettyOstream messager;
            messager.error(errorMsg.str());
            exit(1);
        }
    }

    // Scaling, seek s such that || e*2^(-s) || < 1/2, and set scale = 2^(-s)
    int s = 0;
    double scale = 1.0;
    if (norm > 0.5) {
        s = std::max<int>(0, static_cast<int>((log(norm) / log(2.0) + 2.0)));
        scale /= double(pow(2.0, s));
        U = (scale * t)
            * (*this);  // Here U is used as temp value due to that H is const
    } else
        U = t * (*this);

    // Horner evaluation of the irreducible fraction.
    // Initialize P (numerator) and Q (denominator)
    H2 = U * U;
    Q = c[p] * I;
    P = c[p - 1] * I;
    int odd = 1;

    for (int k = p - 1; k > 0; --k) {
        if (odd == 1) {
            Q = Q * H2 + (c[k - 1] * I);
        } else {
            P = P * H2 + (c[k - 1] * I);
        }
        odd = 1 - odd;
    }
    if (odd == 1) {
        Q = Q * U;
    } else {
        P = P * U;
    }

    Q -= P;

    // Invert Q (SU(3) only):
    H2 = cofactorInverse(Q);

    if (odd == 1) {
        U = -1. * ((2.0 * H2 * P) + I);
    } else {
        U = (2.0 * H2 * P) + I;
    }

    // square
    for (int i = 0; i < s; ++i) U = U * U;

    *this = U;

    return *this;
}

complex<double> Matrix::det() const {
    return e_[0] * e_[4] * e_[8] + e_[1] * e_[5] * e_[6] + e_[2] * e_[3] * e_[7]
           - e_[2] * e_[4] * e_[6] - e_[1] * e_[3] * e_[8]
           - e_[5] * e_[7] * e_[0];
}

complex<double> Matrix::trace() const { return e_[0] + e_[4] + e_[8]; }

complex<double> Matrix::traceOfProductOfMatrix(Matrix &a, Matrix &b) const {
    return a(0) * b(0) + a(1) * b(3) + a(2) * b(6) + a(3) * b(1) + a(4) * b(4)
           + a(5) * b(7) + a(6) * b(2) + a(7) * b(5) + a(8) * b(8);
}

std::string Matrix::MatrixToString() {
    std::stringstream output;
    output.precision(15);
    for (int i = 0; i < 9; i++) {
        if (i > 0) output << " ";
        output << e_[i].real() << " " << e_[i].imag();
    }
    return output.str();
}

double Matrix::frobeniusNorm() const {
    int n = this->getNDim();
    double norm = 0.;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            norm += abs((*this)(i, j)) * abs((*this)(i, j));
        }
    }

    norm = sqrt(norm);

    return norm;
}

double Matrix::oneNorm() const {
    int n = this->getNDim();
    double maxColSum = 0.0;

    for (int j = 0; j < n; j++) {
        double colSum = 0.0;

        for (int i = 0; i < n; i++) {
            colSum += abs((*this)(i, j));
        }

        maxColSum = std::max(maxColSum, colSum);
    }

    return maxColSum;
}

Matrix &Matrix::inv() {
    *this = cofactorInverse(*this);
    return *this;
}

// Pade approximant of log(I+A) (I is unit matrix). good for A\sim I
Matrix &Matrix::logmPade(const int m) {
    Matrix S(0.);
    Matrix A;
    A = *this;
    Matrix I(1.);
    Matrix D;     // denominator
    Matrix invD;  // denominator
    double xi;
    double wi;
    gsl_integration_glfixed_table *table;
    table = gsl_integration_glfixed_table_alloc(m);

    for (int i = 0; i < m; i++) {
        gsl_integration_glfixed_point(0., 1., i, &xi, &wi, table);
        D = I + xi * A;
        // compute inverse of D:
        invD = D;
        invD.inv();
        S = S + wi * (A * invD);
    }

    gsl_integration_glfixed_table_free(table);

    *this = S;

    return *this;
}

// Matrix square root by product from Denman-Beavers (DB) iteration.
// computes principal square root X of the matrix A using the product form
// of the Denman-Beavers iteration. The matrix M tends to I.
// scale specifies scaling: 0, no scaling. 1, determinant scaling (default)
// maxit is the number of iterations.
// Adabted from The Matrix Function Toolbox by Nick Higham (MATLAB code)
Matrix &Matrix::sqrtm(const int scale) {
    const int n = this->getNDim();
    int sc = scale;
    double eps = 1e-2;
    double tol = sqrt(static_cast<double>(n)) * smallEps / 2.;
    double g;
    double Mres;
    double reldiff;
    Matrix X;
    Matrix Xold;
    Matrix M;
    Matrix invM;
    Matrix I(1.);
    Matrix Mr;
    Matrix XmXo;

    X = *this;
    M = *this;

    int maxit = 25;  // maximal number of iterations

    for (int k = 0; k < maxit; k++) {
        if (sc == 1) {
            g = pow(abs(M.det()), -1. / (2. * n));
            X = g * X;
            M = g * g * M;
        }

        Xold = X;
        invM = M;
        invM.inv();

        X = X * (I + invM) / 2.;
        M = 0.5 * (I + (M + invM) / 2.);

        Mr = M - I;
        Mres = Mr.frobeniusNorm();

        XmXo = X - Xold;

        reldiff = XmXo.frobeniusNorm() / X.frobeniusNorm();
        if (reldiff < eps) sc = 0;  // switch to no scaling

        if (Mres <= tol) break;
    }

    *this = X;

    return *this;
}

// matrix logarithm using Pade approximant (inverse scaling and squaring)
// A.H. Al-Mohy and N.J. Higham, Improved Inverse Scaling and Squaring
// algorithms for the matrix logarithm, MIMS eprint 2011.83 this is not the
// improved one, just standard.
Matrix &Matrix::logm() {
    const Matrix I(1.);
    Matrix X;
    Matrix L;
    int k, p, itk;
    double normdiff;
    int j1, j2 = 0.;
    Matrix M;
    int m;

    double xvals[16] = {
        1.586970738772063e-005, 2.313807884242979e-003, 1.938179313533253e-002,
        6.209171588994762e-002, 1.276404810806775e-001, 2.060962623452836e-001,
        2.879093714241194e-001, 3.666532675959788e-001, 4.389227326152340e-001,
        5.034050432047666e-001, 5.600071293013720e-001, 6.092525642521717e-001,
        6.519202543720032e-001, 6.888477797186464e-001, 7.208340678820352e-001,
        7.485977242539218e-001};

    X = *this;
    k = 0;
    p = 0;
    itk = 5;

    while (1) {
        M = X - I;
        normdiff = M.oneNorm();

        if (normdiff <= xvals[15]) {
            p = p + 1;
            // set j1
            for (int i = 0; i < 16; i++) {
                if (normdiff <= xvals[i]) {
                    j1 = i;
                    break;
                }
            }

            // set j2
            for (int i = 0; i < 16; i++) {
                if (normdiff / 2. <= xvals[i]) {
                    j2 = i;
                    break;
                }
            }

            if ((2 * static_cast<double>(j1 - j2) / 3. < itk) || (p == 2)) {
                m = j1;
                break;  // break while loop
            }
        }

        X.sqrtm();  // take the square root

        k = k + 1;

    }  // while(1) loop

    L = X - I;
    L.logmPade(m);

    X = pow(2., k) * L;

    *this = X;

    return *this;
}

Matrix Matrix::fromAlgebraExponent(const std::vector<double> &Q) {
    // the traceless Hermitian Q^a t^a, with t^a = lambda^a/2 (Gell-Mann, the
    // generators of Group)
    const complex<double> I(0., 1.);
    const double invSqrt3 = 1. / std::sqrt(3.);
    Matrix H(noInit);
    H.set(0, 0, 0.5 * (Q[2] + invSqrt3 * Q[7]));
    H.set(1, 1, 0.5 * (-Q[2] + invSqrt3 * Q[7]));
    H.set(2, 2, -invSqrt3 * Q[7]);
    H.set(0, 1, 0.5 * (Q[0] - I * Q[1]));
    H.set(1, 0, 0.5 * (Q[0] + I * Q[1]));
    H.set(0, 2, 0.5 * (Q[3] - I * Q[4]));
    H.set(2, 0, 0.5 * (Q[3] + I * Q[4]));
    H.set(1, 2, 0.5 * (Q[5] - I * Q[6]));
    H.set(2, 1, 0.5 * (Q[5] + I * Q[6]));
    return expiHermitian(H);
}

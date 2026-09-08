#include "Matrix.h"

constexpr Matrix::NoInitTag Matrix::noInit;

#include <gsl/gsl_integration.h>  // include gsl for Gauss-Legendre nodes and weights for log Pade

#include <iostream>
#include <sstream>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::vector;

static_assert(
    sizeof(Matrix) == 9 * sizeof(std::complex<double>),
    "Matrix must remain an exact contiguous complex<double>[9]");

namespace {

inline void requireSU3Dimension(int n) {
    if (n != 3) {
        std::cerr << "Error: fixed Matrix is SU(3)-only; requested " << n << "x"
                  << n << " matrix. Exiting." << std::endl;
        std::exit(1);
    }
}

}  // namespace

Matrix::Matrix() {
    for (int i = 0; i < 9; ++i) e[i] = complex<double>(0.0, 0.0);
}

Matrix::Matrix(int n) {
    requireSU3Dimension(n);
    for (int i = 0; i < 9; ++i) e[i] = complex<double>(0.0, 0.0);
}

Matrix::Matrix(int n, double a) {
    requireSU3Dimension(n);
    for (int i = 0; i < 9; ++i) e[i] = complex<double>(0.0, 0.0);
    e[0] = complex<double>(a, 0.0);
    e[4] = complex<double>(a, 0.0);
    e[8] = complex<double>(a, 0.0);
}

Matrix::Matrix(int n, NoInitTag) { requireSU3Dimension(n); }

// MaxTr version of reunitarization
void Matrix::reu2() {
    Matrix A1(3, 0.);
    Matrix A2(3, 0.);
    Matrix A3(3, 0.);

    Matrix G(3);
    Matrix E(3);

    for (int i = 0; i < 10; i++) {
        E = *this;
        complex<double> N1 = sqrt(
            (conj(e[0]) + e[4]) * conj(conj(e[0]) + e[4])
            + (conj(e[3]) - e[1]) * conj(conj(e[3]) - e[1]));
        complex<double> N2 = sqrt(
            (conj(e[0]) + e[8]) * conj(conj(e[0]) + e[8])
            + (conj(e[6]) - e[2]) * conj(conj(e[6]) - e[2]));
        complex<double> N3 = sqrt(
            (conj(e[4]) + e[8]) * conj(conj(e[4]) + e[8])
            + (conj(e[7]) - e[5]) * conj(conj(e[7]) - e[5]));

        G = (1. / N1) * E;
        A1.set(0, 0, conj(G(0)) + G(4));
        A1.set(0, 1, -G(1) + conj(G(3)));
        A1.set(0, 2, 0.);
        A1.set(1, 0, conj(G(1)) - G(3));
        A1.set(1, 1, G(0) + conj(G(4)));
        A1.set(1, 2, 0.);
        A1.set(2, 0, 0.);
        A1.set(2, 1, 0.);
        A1.set(2, 2, 1.);

        G = (1. / N2) * E;
        A2.set(0, 0, conj(G(0)) + G(8));
        A2.set(0, 1, 0.);
        A2.set(0, 2, -G(2) + conj(G(6)));
        A2.set(1, 0, 0.);
        A2.set(1, 1, 1.);
        A2.set(1, 2, 0.);
        A2.set(2, 0, conj(G(2)) - G(6));
        A2.set(2, 1, 0.);
        A2.set(2, 2, G(0) + conj(G(8)));

        G = (1. / N3) * E;
        A3.set(0, 0, 1.);
        A3.set(0, 1, 0.);
        A3.set(0, 2, 0.);
        A3.set(1, 0, 0.);
        A3.set(1, 1, conj(G(4)) + G(8));
        A3.set(1, 2, -G(5) + conj(G(7)));
        A3.set(2, 0, 0.);
        A3.set(2, 1, conj(G(5)) - G(7));
        A3.set(2, 2, G(4) + conj(G(8)));

        *this = A1 * A2 * A3;
    }
}

// operators:

Matrix operator*(const Matrix &a, const Matrix &b) {
    Matrix c(3, Matrix::noInit);
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
    Matrix aa(a.getNDim(), Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) aa.set(i, a(i) - b(i));
    return aa;
}

//+
Matrix operator+(const Matrix &a, const Matrix &b) {
    Matrix aa(a.getNDim(), Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) aa.set(i, a(i) + b(i));
    return aa;
}

//* multiply by a real scalar
Matrix operator*(const Matrix &a, const double s) {
    Matrix aa(a.getNDim(), Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) {
        aa.set(i, a(i) * s);
    }
    return aa;
}
Matrix operator*(const double s, const Matrix &a) {
    Matrix aa(a.getNDim(), Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) {
        aa.set(i, a(i) * s);
    }
    return aa;
}

//* multiply by a complex number
Matrix operator*(const complex<double> s, const Matrix &a) {
    Matrix aa(a.getNDim(), Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) {
        aa.set(i, a(i) * s);
    }
    return aa;
}

// / division by scalar
Matrix operator/(const Matrix &a, const double s) {
    Matrix aa(a.getNDim(), Matrix::noInit);
    for (int i = 0; i < a.getNN(); i++) aa.set(i, a(i) / s);
    return aa;
}

Matrix &Matrix::conjg() {
    const complex<double> a01 = e[1];
    const complex<double> a02 = e[2];
    const complex<double> a12 = e[5];
    e[0] = conj(e[0]);
    e[4] = conj(e[4]);
    e[8] = conj(e[8]);
    e[1] = conj(e[3]);
    e[2] = conj(e[6]);
    e[5] = conj(e[7]);
    e[3] = conj(a01);
    e[6] = conj(a02);
    e[7] = conj(a12);
    return *this;
}

Matrix Matrix::prodABconj(const Matrix &a, const Matrix &b) {
    Matrix c(3, Matrix::noInit);
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
    Matrix c(3, Matrix::noInit);
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

Matrix &Matrix::imag() {
    Matrix dagger = *this;
    dagger.conjg();
    *this -= dagger;
    return *this;
}

// matrix exponential e^iQ of traceless Hermitian matrices, using coefficients
// Q^a of generators t^a as argument. Dimension is Nc
void Matrix::expmCoeff(const double *Q, complex<double> result[9]) const {
    const int Nc2m1 = 8;
    double sqrt3 = sqrt(3.);
    complex<double> f0, f1, f2, iu, u0, ua[8];
    double c0 = 0., c0max, u, w, xi0, den, thetaOverThree;

    c0 = sqrt3 * (Q[0] * Q[0] * Q[7] + Q[1] * Q[1] * Q[7] + Q[2] * Q[2] * Q[7]);
    c0 -= Q[7] * Q[7] * Q[7] / sqrt3;
    c0 -= (sqrt3 / 2.)
          * (Q[3] * Q[3] * Q[7] + Q[4] * Q[4] * Q[7] + Q[5] * Q[5] * Q[7]
             + Q[6] * Q[6] * Q[7]);
    c0 += 3.
          * (Q[0] * Q[3] * Q[5] + Q[0] * Q[4] * Q[6] + Q[1] * Q[4] * Q[5]
             - Q[1] * Q[3] * Q[6]);
    c0 += 1.5
          * (Q[2] * Q[3] * Q[3] + Q[2] * Q[4] * Q[4] - Q[2] * Q[5] * Q[5]
             - Q[2] * Q[6] * Q[6]);

    c0 /= 12.;

    double c1 = 0.;
    for (int a = 0; a < Nc2m1; a++) {
        c1 += Q[a] * Q[a];
    }
    c1 *= 0.25;

    const double c1Over3 = c1 / 3.;
    const double sqrtC1Over3 = sqrt(c1Over3);
    c0max = std::max(1e-15, 2. * c1Over3 * sqrtC1Over3);

    thetaOverThree = acos(c0 / c0max) / 3.;

    u = sqrtC1Over3 * cos(thetaOverThree);
    w = sqrt(c1) * sin(thetaOverThree);

    xi0 = sin(w) / w;

    den = 9. * u * u - w * w;

    iu = complex<double>(0, 1) * u;

    double cosw = cos(w);
    // iu is purely imaginary.  Construct the two phase factors directly
    // instead of routing them through the general complex exponential.
    const double sinu = sin(u);
    const double cosu = cos(u);
    complex<double> exp2iu(cos(2. * u), sin(2. * u));
    complex<double> expmiu(cosu, -sinu);

    f0 = (u * u - w * w) * exp2iu
         + expmiu * (8. * u * u * cosw + 2. * iu * xi0 * (3. * u * u + w * w));
    f0 /= den;

    f1 = 2. * u * exp2iu
         - expmiu
               * (2. * u * cosw
                  - complex<double>(0., 1.) * (3. * u * u - w * w) * xi0);
    f1 /= den;

    f2 = exp2iu - expmiu * (cosw + 3. * iu * xi0);
    f2 /= den;

    u0 = f0 + 2. / 3. * c1 * f2;
    // if (std::isnan(real(u0))) {
    //     for (int i = 0; i < Nc2m1; i++) {
    //         cout << Q[i] << " ";
    //     }
    //     cout << endl;
    //     cout << "c0 = " << c0 << endl;
    //     cout << "c0max=" << c0max << endl;
    //     cout << "c0/c0max=" << c0 / c0max << endl;
    //     cout << "thetaOverThree=" << thetaOverThree << endl;
    //     cout << "u = " << u << endl;
    //     cout << "w = " << w << endl;
    //     cout << "f0=" << f0 << endl;
    //     cout << "f2=" << f2 << endl;
    //     cout << "c1=" << c1 << endl;
    //     exit(1);
    // }

    // The historical code divided f1 by (0.5*f2), added the real
    // quadratic SU(3) coefficients, and then multiplied the whole result by
    // (0.5*f2) again.  Algebraically cancel that complex division:
    //   (f1/(0.5*f2) * Q + q2) * (0.5*f2)
    //       = f1 * Q + q2 * (0.5*f2).
    // Besides being cheaper, this avoids an unnecessary division when f2 is
    // very small.
    const complex<double> halfF2 = 0.5 * f2;

    ua[0] = Q[3] * Q[5] + Q[4] * Q[6] + 2. / sqrt3 * Q[0] * Q[7];
    ua[1] = 2. * Q[1] * Q[7] / sqrt3 - Q[3] * Q[6] + Q[4] * Q[5];
    ua[2] = 2. * Q[2] * Q[7] / sqrt3 + 0.5 * Q[3] * Q[3] + 0.5 * Q[4] * Q[4]
            - 0.5 * Q[5] * Q[5] - 0.5 * Q[6] * Q[6];
    ua[3] = -1. / sqrt3 * Q[3] * Q[7] + Q[0] * Q[5] - Q[1] * Q[6] + Q[2] * Q[3];
    ua[4] = -1. / sqrt3 * Q[4] * Q[7] + Q[0] * Q[6] + Q[1] * Q[5] + Q[2] * Q[4];
    ua[5] = -1. / sqrt3 * Q[5] * Q[7] + Q[0] * Q[3] + Q[1] * Q[4] - Q[2] * Q[5];
    ua[6] = -1. / sqrt3 * Q[6] * Q[7] + Q[0] * Q[4] - Q[1] * Q[3] - Q[2] * Q[6];
    ua[7] = (Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2] - Q[7] * Q[7]
             - 0.5 * Q[3] * Q[3] - 0.5 * Q[4] * Q[4] - 0.5 * Q[5] * Q[5]
             - 0.5 * Q[6] * Q[6])
            / sqrt3;

    result[0] = u0;
    for (int i = 0; i < 8; i++) {
        result[i + 1] = f1 * Q[i] + halfF2 * ua[i];
    }

    // Check potential NaNs
    for (int i = 0; i < 9; i++) {
        if (std::isnan(result[i].real()) or std::isnan(result[i].imag())) {
            // Sometimes in the very low density region we may encounter
            // (numerically) 0/0 situations In that case, set coefficient to 0,
            // so this contributes only a unit matrix (=vacuum contribution)
            result[i] = 0;
        }
    }
}

vector<complex<double>> Matrix::expmCoeff(std::vector<double> &Q, int Nc) {
    requireSU3Dimension(Nc);
    complex<double> coeff[9];
    expmCoeff(Q.data(), coeff);
    return vector<complex<double>>(coeff, coeff + 9);
}

// matrix exponential using Pade approximant
// t is a scalar that multiplies the matrix (default: t=1) and p is the order in
// the Pade approximant (default: p=6)
Matrix &Matrix::expm(double t, const int p) {
    const int n = this->getNDim();
    const Matrix I(n, 1.);
    Matrix U(n), H2(n), P(n), Q(n);
    double norm = 0.0;
    // Calculate Pade coefficients
    if (p < 6) {
        cout << "Matrix::expm: p should be at least 6. Exiting." << endl;
        exit(0);
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
            cerr << "Null input error in the template expm_pad.\n";
            cout << "Null INPUT : " << *this << "\n";
            exit(0);
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
        U = *this;

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
    H2.set(0, 0, (Q(1, 1) * Q(2, 2) - Q(1, 2) * Q(2, 1)));
    H2.set(0, 1, (Q(0, 2) * Q(2, 1) - Q(0, 1) * Q(2, 2)));
    H2.set(0, 2, (Q(0, 1) * Q(1, 2) - Q(0, 2) * Q(1, 1)));
    H2.set(1, 0, (Q(1, 2) * Q(2, 0) - Q(1, 0) * Q(2, 2)));
    H2.set(1, 1, (Q(0, 0) * Q(2, 2) - Q(0, 2) * Q(2, 0)));
    H2.set(1, 2, (Q(0, 2) * Q(1, 0) - Q(0, 0) * Q(1, 2)));
    H2.set(2, 0, (Q(1, 0) * Q(2, 1) - Q(1, 1) * Q(2, 0)));
    H2.set(2, 1, (Q(0, 1) * Q(2, 0) - Q(0, 0) * Q(2, 1)));
    H2.set(2, 2, (Q(0, 0) * Q(1, 1) - Q(0, 1) * Q(1, 0)));
    H2 *= 1.
          / (Q(0, 0) * Q(1, 1) * Q(2, 2) + Q(0, 1) * Q(1, 2) * Q(2, 0)
             + Q(0, 2) * Q(1, 0) * Q(2, 1) - Q(0, 2) * Q(1, 1) * Q(2, 0)
             - Q(0, 1) * Q(1, 0) * Q(2, 2) - Q(1, 2) * Q(2, 1) * Q(0, 0));

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

complex<double> Matrix::det() {
    return e[0] * e[4] * e[8] + e[1] * e[5] * e[6] + e[2] * e[3] * e[7]
           - e[2] * e[4] * e[6] - e[1] * e[3] * e[8] - e[5] * e[7] * e[0];
}

complex<double> Matrix::trace() const { return e[0] + e[4] + e[8]; }

complex<double> Matrix::traceOfProdcutOfMatrix(Matrix &M1, Matrix &M2) const {
    return M1(0) * M2(0) + M1(1) * M2(3) + M1(2) * M2(6) + M1(3) * M2(1)
           + M1(4) * M2(4) + M1(5) * M2(7) + M1(6) * M2(2) + M1(7) * M2(5)
           + M1(8) * M2(8);
}

std::string Matrix::MatrixToString() {
    std::stringstream output;
    output.precision(15);
    output << e[0].real() << " " << e[0].imag() << " " << e[3].real() << " "
           << e[3].imag() << " " << e[6].real() << " " << e[6].imag() << " "
           << e[1].real() << " " << e[1].imag() << " " << e[4].real() << " "
           << e[4].imag() << " " << e[7].real() << " " << e[7].imag() << " "
           << e[2].real() << " " << e[2].imag() << " " << e[5].real() << " "
           << e[5].imag() << " " << e[8].real() << " " << e[8].imag();
    return output.str();
}

double Matrix::FrobeniusNorm() {
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

double Matrix::OneNorm() {
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
    Matrix Q = *this;
    Matrix H2(3, Matrix::noInit);
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
    *this = H2;
    return *this;
}

// Pade approximant of log(I+A) (I is unit matrix). good for A\sim I
Matrix &Matrix::logm_pade(const int m) {
    const int n = this->getNDim();
    Matrix S(n, 0.);
    Matrix A(n);
    A = *this;
    Matrix I(n, 1.);
    Matrix D(n);     // denominator
    Matrix invD(n);  // denominator
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
    double tol = sqrt(static_cast<double>(n)) * 1e-16 / 2.;
    double g;
    double Mres;
    double reldiff;
    Matrix X(n);
    Matrix Xold(n);
    Matrix M(n);
    Matrix invM(n);
    Matrix I(n, 1.);
    Matrix Mr(n);
    Matrix XmXo(n);

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
        Mres = Mr.FrobeniusNorm();

        XmXo = X - Xold;

        reldiff = XmXo.FrobeniusNorm() / X.FrobeniusNorm();
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
    const int n = this->getNDim();
    const Matrix I(n, 1.);
    Matrix X(n);
    Matrix L(n);
    int k, p, itk;
    double normdiff;
    int j1, j2 = 0.;
    Matrix M(n);
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
        normdiff = M.OneNorm();

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
    L.logm_pade(m);

    X = pow(2., k) * L;

    *this = X;

    return *this;
}

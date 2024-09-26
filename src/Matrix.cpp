#include "Matrix.h"

#include <gsl/gsl_integration.h>  // include gsl for Gauss-Legendre nodes and weights for log Pade

#include <iostream>
#include <sstream>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::stringstream;

// constructor if just dimension is given
Matrix::Matrix(int n) {
    ndim = n;
    nn = ndim * ndim;
    // e = new complex<double> [nn];
    e.resize(nn);
    for (int i = 0; i < nn; i++) e[i] = complex<double>(0.0, 0.0);
}

// constructor if value for a and dimensions are given (a is the real value on
// the diagonal)
Matrix::Matrix(int n, double a) {
    ndim = n;
    nn = ndim * ndim;
    // e = new complex<double> [nn];
    e.resize(nn);
    // if(e==0)
    //  {
    //    cout << "Matrix: cannot allocate memory (matrix:) e= " << e << endl;
    //    abort();
    //  }

    for (int i = 0; i < nn; i++) e[i] = complex<double>(0.0, 0.0);
    for (int i = 0; i < ndim; i++) e[i * ndim + i] = complex<double>(a, 0.0);
}

// MaxTr version of reunitarization
void Matrix::reu2() {
    Matrix A1(ndim, 0.);
    Matrix A2(ndim, 0.);
    Matrix A3(ndim, 0.);

    Matrix G(ndim);
    Matrix E(ndim);

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
    int n = a.getNDim();
    Matrix c(n);
    if (n == 2) {
        c.set(0, 0, a(0, 0) * b(0, 0) + a(0, 1) * b(1, 0));
        c.set(0, 1, a(0, 0) * b(0, 1) + a(0, 1) * b(1, 1));
        c.set(1, 0, a(1, 0) * b(0, 0) + a(1, 1) * b(1, 0));
        c.set(1, 1, a(1, 0) * b(0, 1) + a(1, 1) * b(1, 1));
        return c;
    } else if (n == 3) {
        c.set(0, 0, a(0, 0) * b(0, 0) + a(0, 1) * b(1, 0) + a(0, 2) * b(2, 0));
        c.set(0, 1, a(0, 0) * b(0, 1) + a(0, 1) * b(1, 1) + a(0, 2) * b(2, 1));
        c.set(0, 2, a(0, 0) * b(0, 2) + a(0, 1) * b(1, 2) + a(0, 2) * b(2, 2));
        c.set(1, 0, a(1, 0) * b(0, 0) + a(1, 1) * b(1, 0) + a(1, 2) * b(2, 0));
        c.set(1, 1, a(1, 0) * b(0, 1) + a(1, 1) * b(1, 1) + a(1, 2) * b(2, 1));
        c.set(1, 2, a(1, 0) * b(0, 2) + a(1, 1) * b(1, 2) + a(1, 2) * b(2, 2));
        c.set(2, 0, a(2, 0) * b(0, 0) + a(2, 1) * b(1, 0) + a(2, 2) * b(2, 0));
        c.set(2, 1, a(2, 0) * b(0, 1) + a(2, 1) * b(1, 1) + a(2, 2) * b(2, 1));
        c.set(2, 2, a(2, 0) * b(0, 2) + a(2, 1) * b(1, 2) + a(2, 2) * b(2, 2));
        return c;
    } else if (n == 8) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                c.set(
                    i, j,
                    a(i, 0) * b(0, j) + a(i, 1) * b(1, j) + a(i, 2) * b(2, j)
                        + a(i, 3) * b(3, j) + a(i, 4) * b(4, j)
                        + a(i, 5) * b(5, j) + a(i, 6) * b(6, j)
                        + a(i, 7) * b(7, j));
            }
        }
        return c;
    } else {
        cout << "[Matrix::operator *]: Matrix product is only defined for 2x2, "
                "3x3, and 8x8 matrixes. You gave me "
             << n << "x" << n << ". Exiting." << endl;
        return c;
    }
}

//-
Matrix operator-(const Matrix &a, const Matrix &b) {
    Matrix aa(a.getNDim());
    for (int i = 0; i < a.getNN(); i++) aa.set(i, a(i) - b(i));
    return aa;
}

//+
Matrix operator+(const Matrix &a, const Matrix &b) {
    Matrix aa(a.getNDim());
    for (int i = 0; i < a.getNN(); i++) aa.set(i, a(i) + b(i));
    return aa;
}

//* multiply by a real scalar
Matrix operator*(const Matrix &a, const double s) {
    Matrix aa(a.getNDim());
    for (int i = 0; i < a.getNN(); i++) {
        aa.set(i, a(i) * s);
    }
    return aa;
}
Matrix operator*(const double s, const Matrix &a) {
    Matrix aa(a.getNDim());
    for (int i = 0; i < a.getNN(); i++) {
        aa.set(i, a(i) * s);
    }
    return aa;
}

//* multiply by a complex number
Matrix operator*(const complex<double> s, const Matrix &a) {
    Matrix aa(a.getNDim());
    for (int i = 0; i < a.getNN(); i++) {
        aa.set(i, a(i) * s);
    }
    return aa;
}

// / division by scalar
Matrix operator/(const Matrix &a, const double s) {
    Matrix aa(a.getNDim());
    for (int i = 0; i < a.getNN(); i++) aa.set(i, a(i) / s);
    return aa;
}

Matrix &Matrix::conjg() {
    //  complex<double>* temp;
    // temp = new complex<double>[nn];
    complex<double> temp[9];
    // save half the matrix
    for (int i = 0; i < ndim; i++)
        for (int j = i; j < ndim; j++) {
            temp[i * ndim + j] = conj(e[i * ndim + j]);
        }
    // transpose half the matrix
    for (int i = 0; i < ndim; i++)
        for (int j = 0; j < i; j++) {
            e[j * ndim + i] = conj(e[i * ndim + j]);
        }
    // transpose the other half using the saved values
    for (int i = 0; i < ndim; i++)
        for (int j = i; j < ndim; j++) {
            e[j * ndim + i] = temp[i * ndim + j];
        }
    //  delete temp;
    return *this;
}

Matrix Matrix::prodABconj(const Matrix &a, const Matrix &b) {
    Matrix c(3);
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
    Matrix c(3);
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
    if (ndim == 2) {
        complex<double> e0 = e[0];
        complex<double> e1 = e[1];
        complex<double> e2 = e[2];
        complex<double> e3 = e[3];
        e[0] -= conj(e0);
        e[1] -= conj(e2);
        e[2] -= conj(e1);
        e[3] -= conj(e3);
    } else {
        cerr << " (Matrix::) invalid dimension nn= " << nn << endl;
        exit(1);
    }
    return *this;
}

// matrix exponential e^iQ of traceless Hermitian matrices, using coefficients
// Q^a of generators t^a as argument. Dimension is Nc
vector<complex<double>> Matrix::expmCoeff(std::vector<double> &Q, int Nc) {
    int Nc2m1 = Nc * Nc - 1;
    vector<complex<double>> result;
    result.reserve(9);
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

    c0max = std::max(1e-15, 2. * pow(c1 / 3., 1.5));

    thetaOverThree = acos(std::max(-1., std::min(1., c0 / c0max))) / 3.;

    u = sqrt(c1 / 3.) * cos(thetaOverThree);
    w = sqrt(c1) * sin(thetaOverThree);

    xi0 = sin(w) / w;

    den = 9. * u * u - w * w;

    iu = complex<double>(0, 1) * u;

    double cosw = cos(w);
    complex<double> exp2iu = exp(2. * iu);
    complex<double> expmiu = exp(-iu);

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

    f1 /= (0.5 * f2);  // will multiply everything by 0.5 f2 again later

    for (int i = 0; i < 8; i++) {
        ua[i] = f1 * Q[i];
    }

    ua[0] += (Q[3] * Q[5] + Q[4] * Q[6] + 2. / sqrt3 * Q[0] * Q[7]);
    ua[1] += (2. * Q[1] * Q[7] / sqrt3 - Q[3] * Q[6] + Q[4] * Q[5]);
    ua[2] +=
        (2. * Q[2] * Q[7] / sqrt3 + 0.5 * Q[3] * Q[3] + 0.5 * Q[4] * Q[4]
         - 0.5 * Q[5] * Q[5] - 0.5 * Q[6] * Q[6]);
    ua[3] +=
        (-1. / sqrt3 * Q[3] * Q[7] + Q[0] * Q[5] - Q[1] * Q[6] + Q[2] * Q[3]);
    ua[4] +=
        (-1. / sqrt3 * Q[4] * Q[7] + Q[0] * Q[6] + Q[1] * Q[5] + Q[2] * Q[4]);
    ua[5] +=
        (-1. / sqrt3 * Q[5] * Q[7] + Q[0] * Q[3] + Q[1] * Q[4] - Q[2] * Q[5]);
    ua[6] +=
        (-1. / sqrt3 * Q[6] * Q[7] + Q[0] * Q[4] - Q[1] * Q[3] - Q[2] * Q[6]);
    ua[7] += (Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2] - Q[7] * Q[7]
              - 0.5 * Q[3] * Q[3] - 0.5 * Q[4] * Q[4] - 0.5 * Q[5] * Q[5]
              - 0.5 * Q[6] * Q[6])
             / sqrt3;

    result.push_back(u0);

    for (int i = 0; i < 8; i++) {
        result.push_back(ua[i] * 0.5 * f2);
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
    return result;
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

    // Invert Q:
    if (n == 2) {
        H2.set(0, 0, Q(1, 1));
        H2.set(0, 1, -Q(0, 1));
        H2.set(1, 0, -Q(1, 0));
        H2.set(1, 1, Q(0, 0));
        H2 *= 1. / (Q(0, 0) * Q(1, 1) - Q(0, 1) * Q(1, 0));  // divide by det(H)
    } else if (n == 3) {
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
        // divide by det(H)
    } else {
        cerr << "My matrix exponential works only for up to 3x3 matrices. If "
                "you "
                "need more, include MatrixExp.h and the boost library."
             << endl;
    }

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
    int n = this->getNDim();
    Matrix Q(n);
    Q = *this;
    complex<double> det;

    if (n == 2) {
        det = Q(0, 0) * Q(1, 1) - Q(0, 1) * Q(1, 0);
    } else if (n == 3) {
        det = Q(0, 0) * Q(1, 1) * Q(2, 2) + Q(0, 1) * Q(1, 2) * Q(2, 0)
              + Q(0, 2) * Q(1, 0) * Q(2, 1) - Q(0, 2) * Q(1, 1) * Q(2, 0)
              - Q(0, 1) * Q(1, 0) * Q(2, 2) - Q(1, 2) * Q(2, 1) * Q(0, 0);
    }

    return det;
}

complex<double> Matrix::trace() {
    int n = this->getNDim();
    Matrix Q(n);
    Q = *this;
    complex<double> trace;

    if (n == 2) {
        trace = Q(0, 0) + Q(1, 1);
    } else if (n == 3) {
        trace = Q(0, 0) + Q(1, 1) + Q(2, 2);
    }

    return trace;
}

complex<double> Matrix::traceOfProdcutOfMatrix(Matrix &M1, Matrix &M2) const {
    // this function returns the trace of the product of two matrices
    complex<double> trace = 0.;
    int n = M1.getNDim();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            trace += M1(i, j) * M2(j, i);
        }
    }
    return (trace);
}

string Matrix::MatrixToString() {
    int n = this->getNDim();
    Matrix Q(n);
    Q = *this;
    stringstream output;
    output.precision(15);

    if (n == 3) {
        output << Q(0, 0).real() << " " << Q(0, 0).imag() << " "
               << Q(1, 0).real() << " " << Q(1, 0).imag() << " "
               << Q(2, 0).real() << " " << Q(2, 0).imag() << " "
               << Q(0, 1).real() << " " << Q(0, 1).imag() << " "
               << Q(1, 1).real() << " " << Q(1, 1).imag() << " "
               << Q(2, 1).real() << " " << Q(2, 1).imag() << " "
               << Q(0, 2).real() << " " << Q(0, 2).imag() << " "
               << Q(1, 2).real() << " " << Q(1, 2).imag() << " "
               << Q(2, 2).real() << " " << Q(2, 2).imag();
    } else if (n == 2) {
    }

    return output.str();
}

double Matrix::FrobeniusNorm() {
    int n = this->getNDim();
    Matrix Q(n);
    Q = *this;
    double norm = 0.;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            norm = abs(Q(i, j)) * abs(Q(i, j));
        }
    }

    norm = sqrt(norm);

    return norm;
}

double Matrix::OneNorm() {
    int n = this->getNDim();
    Matrix Q(n);
    Q = *this;
    double norm[3];
    double onenorm;

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            norm[j] = abs(Q(i, j));
        }
    }

    onenorm = std::max(norm[0], norm[1]);
    onenorm = std::max(onenorm, norm[2]);

    return onenorm;
}

Matrix &Matrix::inv() {
    int n = this->getNDim();
    Matrix H2(n);
    Matrix Q(n);
    Q = *this;

    if (n == 2) {
        H2.set(0, 0, Q(1, 1));
        H2.set(0, 1, -Q(0, 1));
        H2.set(1, 0, -Q(1, 0));
        H2.set(1, 1, Q(0, 0));
        H2 *= 1. / (Q(0, 0) * Q(1, 1) - Q(0, 1) * Q(1, 0));  // divide by det(H)
    } else if (n == 3) {
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
                 - Q(0, 1) * Q(1, 0) * Q(2, 2)
                 - Q(1, 2) * Q(2, 1) * Q(0, 0));  // divide by det(H)
    }

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

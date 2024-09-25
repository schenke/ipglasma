#ifndef Matrix_h
#define Matrix_h

#include <gsl/gsl_integration.h>  // include gsl for Gauss-Legendre nodes and weights for log Pade

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "Spinor.h"

using namespace std;

class Matrix {
  private:
    int ndim;
    int nn;
    std::vector<complex<double>> e;

  public:
    // constructor(s)
    Matrix() = default;
    Matrix(int n);
    Matrix(int n, double a);

    // destructor
    ~Matrix() {}

    Matrix &inv();
    Matrix &logm_pade(const int m);
    Matrix &sqrtm(const int scale = 1);
    double OneNorm();
    double FrobeniusNorm();

    Matrix &logm();

    void setRe(int i, double a) { e[i] = complex<double>(a, e[i].imag()); };
    void setRe(int i, int j, double a) {
        e[j + ndim * i] = complex<double>(a, e[j + ndim * i].imag());
    };
    void setIm(int i, double a) { e[i] = complex<double>(e[i].real(), a); };
    void setIm(int i, int j, double a) {
        e[j + ndim * i] = complex<double>(e[j + ndim * i].real(), a);
    };

    void set(int i, complex<double> a) { e[i] = a; };
    void set(int i, int j, complex<double> a) { e[j + ndim * i] = a; };

    complex<double> get(int i) { return e[i]; };
    complex<double> get(int i, int j) { return e[j + ndim * i]; };

    double getRe(int i) { return e[i].real(); };
    double getIm(int i) { return e[i].imag(); };

    int getNDim() const { return ndim; }
    int getNN() const { return nn; }

    string MatrixToString();

    Matrix &expm(double t = 1.0, const int p = 6);

    // Matrix exponential of traceless hermitian
    // matrix using coefficients of t^a as input
    vector<complex<double>> expmCoeff(std::vector<double> &Q, int n);

    complex<double> det();
    complex<double> trace();

    void reu() {
        if (ndim == 3) {
            Spinor e1(ndim);
            Spinor e2(ndim);
            Spinor e3(ndim);

            Spinor a1(ndim, e[0], e[1], e[2]);
            Spinor a2(ndim, e[3], e[4], e[5]);

            e1 = a1.normalize();
            e2 = a2.GramSchmidt(e1);
            e3 = (e1 % e2).normalize();

            e[0] = e1(0);
            e[1] = e1(1);
            e[2] = e1(2);
            e[3] = e2(0);
            e[4] = e2(1);
            e[5] = e2(2);
            e[6] = e3(0);
            e[7] = e3(1);
            e[8] = e3(2);
        } else if (ndim == 2) {
            Spinor e1(ndim);
            Spinor e2(ndim);
            Spinor a1(ndim, e[0], e[1]);
            Spinor a2(ndim, e[2], e[3]);

            // cout << "a1=" << a1 << endl << endl;
            // cout << "a2=" << a2 << endl << endl;

            e1 = a1.normalize();
            e2 = a2.GramSchmidt(e1);

            // cout << "e1=" << e1 << endl << endl;
            // cout << "e2=" << e2 << endl << endl;

            e[0] = e1(0);
            e[1] = e1(1);
            e[2] = e2(0);
            e[3] = e2(1);
        }
    }

    void reu2();

    // operators:

    //()
    std::complex<double> operator()(const int i) const { return e[i]; }
    std::complex<double> operator()(const int i, const int j) const {
        return e[j + ndim * i];
    }
    // std::complex<double> operator () (const int i, const int j) { return
    // e[j+ndim*i]; }

    //=
    // const Matrix& operator = (const Matrix& p)  {
    //    nn = p.getNN();
    //    if (&p != this ) {
    //        for (int i=0; i<nn; i++) {
    //            e[i] = p.e[i];
    //        }
    //    }
    //    return *this;
    //}

    //==
    bool operator==(const Matrix &p) const {
        for (int i = 0; i < nn; i++)
            if (e[i] != p.e[i]) return false;
        return true;
    }

    //!=
    bool operator!=(const Matrix &p) const {
        for (int i = 0; i < nn; i++)
            if (e[i] != p.e[i]) return true;
        return false;
    }

    //+=
    Matrix &operator+=(const Matrix &a) {
        for (int i = 0; i < nn; i++) e[i] += a.e[i];
        return *this;
    }

    //-=
    Matrix &operator-=(const Matrix &a) {
        for (int i = 0; i < nn; i++) e[i] -= a.e[i];
        return *this;
    }

    //*=
    Matrix &operator*=(const complex<double> a) {
        for (int i = 0; i < nn; i++) e[i] *= a;
        return *this;
    }

    // /=
    Matrix &operator/=(const complex<double> a) {
        for (int i = 0; i < nn; i++) e[i] /= a;
        return *this;
    }

    double square() const {
        double tr = 0.0;
        for (int i = 0; i < nn; i++) {
            tr += e[i].real() * e[i].real() + e[i].imag() * e[i].imag();
        }
        return 0.5 * tr;
    }

    Matrix &imag();

    Matrix &conjg();
    Matrix prodABconj(const Matrix &a, const Matrix &b);
    Matrix prodAconjB(const Matrix &a, const Matrix &b);

    complex<double> traceOfProdcutOfMatrix(Matrix &M1, Matrix &M2) const;

    //<<
    friend ostream &operator<<(ostream &os, const Matrix &p) {
        for (int i = 0; i < p.getNDim(); i++) {
            for (int j = 0; j < p.getNDim(); j++) os << p(i, j);
            if (i < p.getNDim() - 1) os << std::endl;
        }
        return os;
    }
};

Matrix operator+(const Matrix &a, const Matrix &b);
Matrix operator-(const Matrix &a, const Matrix &b);
Matrix operator-(const Matrix &a);
Matrix operator/(const Matrix &a, const Matrix &b);

Matrix operator*(const double a, const Matrix &b);
Matrix operator*(const std::complex<double> a, const Matrix &b);
Matrix operator*(const Matrix &a, const double b);
Matrix operator*(const Matrix &a, const Matrix &b);

#endif

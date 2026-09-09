#ifndef Matrix_h
#define Matrix_h

#include <complex>
#include <string>
#include <vector>

#include "Spinor.h"

using std::complex;
using std::ostream;

// Fundamental-color matrix used throughout IP-Glasma.
//
// IP-Glasma is SU(3)-only, so this type is deliberately a fixed 3x3 matrix.
// The object contains exactly nine std::complex<double> values and no shape,
// heap, or self-pointer metadata.  This makes std::vector<Matrix> a genuinely
// contiguous complex<double>[9] lattice field with a 144-byte stride.
class Matrix {
  private:
    static constexpr int kN = 3;
    static constexpr int kNN = 9;
    complex<double> e[kNN];

  public:
    struct NoInitTag {};
    static constexpr NoInitTag noInit {};

    Matrix();
    explicit Matrix(int n);
    Matrix(int n, double a);
    Matrix(int n, NoInitTag);

    Matrix(const Matrix &) = default;
    Matrix &operator=(const Matrix &) = default;
    ~Matrix() = default;

    complex<double> *data() { return e; }
    const complex<double> *data() const { return e; }

    Matrix &inv();
    Matrix &logm_pade(const int m);
    Matrix &sqrtm(const int scale = 1);
    double OneNorm();
    double FrobeniusNorm();

    Matrix &logm();

    void setRe(int i, double a) { e[i] = complex<double>(a, e[i].imag()); }
    void setRe(int i, int j, double a) {
        e[j + kN * i] = complex<double>(a, e[j + kN * i].imag());
    }
    void setIm(int i, double a) { e[i] = complex<double>(e[i].real(), a); }
    void setIm(int i, int j, double a) {
        e[j + kN * i] = complex<double>(e[j + kN * i].real(), a);
    }

    void set(int i, complex<double> a) { e[i] = a; }
    void set(int i, int j, complex<double> a) { e[j + kN * i] = a; }

    complex<double> get(int i) const { return e[i]; }
    complex<double> get(int i, int j) const { return e[j + kN * i]; }

    double getRe(int i) const { return e[i].real(); }
    double getIm(int i) const { return e[i].imag(); }

    int getNDim() const { return kN; }
    int getNN() const { return kNN; }

    std::string MatrixToString();

    Matrix &expm(double t = 1.0, const int p = 6);

    // Matrix exponential of traceless hermitian matrix using coefficients of
    // the eight SU(3) fundamental generators as input.
    // Allocation-free SU(3) exponential coefficients for hot paths.
    void expmCoeff(const double *Q, complex<double> out[9]) const;
    std::vector<complex<double>> expmCoeff(std::vector<double> &Q, int n);

    complex<double> det();
    complex<double> trace() const;

    void reu() {
        Spinor e1(kN);
        Spinor e2(kN);
        Spinor e3(kN);

        Spinor a1(kN, e[0], e[1], e[2]);
        Spinor a2(kN, e[3], e[4], e[5]);

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
    }

    void reu2();

    std::complex<double> operator()(const int i) const { return e[i]; }
    std::complex<double> operator()(const int i, const int j) const {
        return e[j + kN * i];
    }

    bool operator==(const Matrix &p) const {
        for (int i = 0; i < kNN; ++i)
            if (e[i] != p.e[i]) return false;
        return true;
    }

    bool operator!=(const Matrix &p) const {
        for (int i = 0; i < kNN; ++i)
            if (e[i] != p.e[i]) return true;
        return false;
    }

    Matrix &operator+=(const Matrix &a) {
        for (int i = 0; i < kNN; ++i) e[i] += a.e[i];
        return *this;
    }

    Matrix &operator-=(const Matrix &a) {
        for (int i = 0; i < kNN; ++i) e[i] -= a.e[i];
        return *this;
    }

    Matrix &operator*=(const complex<double> a) {
        for (int i = 0; i < kNN; ++i) e[i] *= a;
        return *this;
    }

    Matrix &operator/=(const complex<double> a) {
        for (int i = 0; i < kNN; ++i) e[i] /= a;
        return *this;
    }

    double square() const {
        double tr = 0.0;
        for (int i = 0; i < kNN; ++i) {
            tr += e[i].real() * e[i].real() + e[i].imag() * e[i].imag();
        }
        return 0.5 * tr;
    }

    Matrix &imag();

    Matrix &conjg();
    Matrix prodABconj(const Matrix &a, const Matrix &b);
    Matrix prodAconjB(const Matrix &a, const Matrix &b);

    complex<double> traceOfProdcutOfMatrix(Matrix &M1, Matrix &M2) const;

    friend ostream &operator<<(ostream &os, const Matrix &p) {
        for (int i = 0; i < kN; ++i) {
            for (int j = 0; j < kN; ++j) os << p(i, j);
            if (i < kN - 1) os << std::endl;
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

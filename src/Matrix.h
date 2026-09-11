#ifndef SRC_MATRIX_H_
#define SRC_MATRIX_H_

#include <complex>
#include <string>
#include <vector>

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
    complex<double> e_[kNN];

  public:
    struct NoInitTag {};
    static constexpr NoInitTag noInit {};

    Matrix();
    explicit Matrix(double a);
    explicit Matrix(NoInitTag);

    Matrix(const Matrix &) = default;
    Matrix &operator=(const Matrix &) = default;
    ~Matrix() = default;

    complex<double> *data() { return e_; }
    const complex<double> *data() const { return e_; }

    Matrix &inv();
    Matrix &logmPade(const int m);
    Matrix &sqrtm(const int scale = 1);
    double oneNorm();
    double frobeniusNorm();

    Matrix &logm();

    void setRe(int i, double a) { e_[i] = complex<double>(a, e_[i].imag()); }
    void setRe(int i, int j, double a) {
        e_[j + kN * i] = complex<double>(a, e_[j + kN * i].imag());
    }
    void setIm(int i, double a) { e_[i] = complex<double>(e_[i].real(), a); }
    void setIm(int i, int j, double a) {
        e_[j + kN * i] = complex<double>(e_[j + kN * i].real(), a);
    }

    void set(int i, complex<double> a) { e_[i] = a; }
    void set(int i, int j, complex<double> a) { e_[j + kN * i] = a; }

    complex<double> get(int i) const { return e_[i]; }
    complex<double> get(int i, int j) const { return e_[j + kN * i]; }

    double getRe(int i) const { return e_[i].real(); }
    double getIm(int i) const { return e_[i].imag(); }

    int getNDim() const { return kN; }
    int getNN() const { return kNN; }

    std::string MatrixToString();

    Matrix &expm(double t = 1.0, const int p = 6);

    // Matrix exponential of traceless hermitian matrix using coefficients of
    // the eight SU(3) fundamental generators as input.
    // Allocation-free SU(3) exponential coefficients for hot paths.
    void expmCoeff(const double *Q, complex<double> out[9]) const;

    complex<double> det();
    complex<double> trace() const;

    std::complex<double> operator()(const int i) const { return e_[i]; }
    std::complex<double> operator()(const int i, const int j) const {
        return e_[j + kN * i];
    }

    Matrix &operator+=(const Matrix &a) {
        for (int i = 0; i < kNN; ++i) e_[i] += a.e_[i];
        return *this;
    }

    Matrix &operator-=(const Matrix &a) {
        for (int i = 0; i < kNN; ++i) e_[i] -= a.e_[i];
        return *this;
    }

    Matrix &operator*=(const complex<double> a) {
        for (int i = 0; i < kNN; ++i) e_[i] *= a;
        return *this;
    }

    Matrix &operator/=(const complex<double> a) {
        for (int i = 0; i < kNN; ++i) e_[i] /= a;
        return *this;
    }

    double square() const {
        double tr = 0.0;
        for (int i = 0; i < kNN; ++i) {
            tr += e_[i].real() * e_[i].real() + e_[i].imag() * e_[i].imag();
        }
        return 0.5 * tr;
    }

    Matrix &conjg();
    Matrix prodABconj(const Matrix &a, const Matrix &b);
    Matrix prodAconjB(const Matrix &a, const Matrix &b);

    complex<double> traceOfProdcutOfMatrix(Matrix &a, Matrix &b) const;

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

Matrix operator*(const double a, const Matrix &b);
Matrix operator*(const std::complex<double> a, const Matrix &b);
Matrix operator*(const Matrix &a, const double b);
Matrix operator*(const Matrix &a, const Matrix &b);

#endif  // SRC_MATRIX_H_

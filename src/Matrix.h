#ifndef SRC_MATRIX_H_
#define SRC_MATRIX_H_

#include <complex>
#include <string>
#include <vector>

using std::complex;
using std::ostream;

/**
 * Fundamental-color \f$3\times3\f$ complex matrix used throughout
 * IP-Glasma.
 *
 * IP-Glasma is SU(3)-only, so this type is deliberately a fixed
 * \f$3\times3\f$ matrix. The object contains exactly nine
 * `std::complex<double>` values and no shape, heap, or self-pointer
 * metadata. This makes `std::vector<Matrix>` a genuinely contiguous
 * `complex<double>[9]` lattice field with a 144-byte stride.
 */
class Matrix {
  private:
    /// Matrix side length (fixed at 3; see the class-level note on why
    /// this codebase doesn't generalize to other \f$N_c\f$).
    static constexpr int kN = 3;
    /// Total element count, `kN*kN`.
    static constexpr int kNN = 9;
    /// The 9 matrix elements, row-major (`e_[kN*row+col]`).
    complex<double> e_[kNN];

  public:
    /// Tag type selecting the uninitialized constructor.
    struct NoInitTag {};
    /// Tag value selecting the uninitialized constructor, for hot
    /// paths that immediately overwrite every element anyway.
    static constexpr NoInitTag noInit {};

    /**
     * Constructs the zero matrix.
     */
    Matrix();
    /**
     * Constructs \p a times the identity matrix.
     * \param[in] a Diagonal value.
     */
    explicit Matrix(double a);
    /**
     * Constructs a Matrix with uninitialized elements, for callers that
     * will immediately overwrite every element via set()/data() anyway
     * and want to skip the zero-initialization cost.
     */
    explicit Matrix(NoInitTag);

    Matrix(const Matrix &) = default;
    Matrix &operator=(const Matrix &) = default;
    ~Matrix() = default;

    /**
     * Returns a mutable pointer to the 9 underlying elements, row-major.
     * \return Pointer to `e_[0]`.
     */
    complex<double> *data() { return e_; }
    /**
     * Returns a read-only pointer to the 9 underlying elements,
     * row-major.
     * \return Pointer to `e_[0]`.
     */
    const complex<double> *data() const { return e_; }

    /**
     * Inverts this matrix in place via its cofactor (adjugate) matrix.
     * \return `*this`, now \f$A^{-1}\f$.
     */
    Matrix &inv();
    /**
     * Replaces this matrix with \f$\log(I+A)\f$ (\f$A\f$ = this
     * matrix), via an \p m-point Gauss-Legendre quadrature of the
     * integral representation
     * \f$\log(I+A) = \int_0^1 A(I+xA)^{-1}\,dx\f$. Good for \f$A\f$
     * close to \f$0\f$ (i.e. this matrix close to \f$I\f$); used by
     * logm() on an already square-root-reduced matrix.
     * \param[in] m Number of quadrature points.
     * \return `*this`, now \f$\log(I+A)\f$.
     */
    Matrix &logmPade(const int m);
    /**
     * Replaces this matrix with its principal square root, via the
     * product form of the Denman-Beavers iteration (up to 25
     * iterations, or until the residual \f$\|M-I\|_F\f$ falls below a
     * dimension-scaled tolerance), adapted from Nick Higham's Matrix
     * Function Toolbox.
     * \param[in] scale `1` to use determinant scaling for faster
     * initial convergence (switched off automatically once the
     * iteration is already converging well), `0` for no scaling.
     * \return `*this`, now \f$A^{1/2}\f$.
     */
    Matrix &sqrtm(const int scale = 1);
    /**
     * Computes this matrix's induced 1-norm (maximum absolute column
     * sum).
     * \return \f$\|A\|_1 = \max_j \sum_i |A_{ij}|\f$.
     */
    double oneNorm() const;
    /**
     * Computes this matrix's Frobenius norm.
     * \return \f$\|A\|_F = \sqrt{\sum_{ij} |A_{ij}|^2}\f$.
     */
    double frobeniusNorm() const;

    /**
     * Replaces this matrix with its principal matrix logarithm, via
     * the (standard, not the "improved") inverse scaling-and-squaring
     * algorithm of Al-Mohy & Higham (MIMS eprint 2011.83): repeatedly
     * takes square roots via sqrtm() until \f$\|X-I\|_1\f$ is small
     * enough for a table-selected Pade order to converge well, then
     * applies logmPade() and rescales by the number of square roots
     * taken.
     * \return `*this`, now \f$\log(A)\f$.
     */
    Matrix &logm();

    /**
     * Sets one element's real part, by flat index, keeping its
     * imaginary part.
     * \param[in] i Flat index, `0`-`8`.
     * \param[in] a New real part.
     */
    void setRe(int i, double a) { e_[i] = complex<double>(a, e_[i].imag()); }
    /**
     * Sets one element's real part, by row/column, keeping its
     * imaginary part.
     * \param[in] i Row index, `0`-`2`.
     * \param[in] j Column index, `0`-`2`.
     * \param[in] a New real part.
     */
    void setRe(int i, int j, double a) {
        e_[j + kN * i] = complex<double>(a, e_[j + kN * i].imag());
    }
    /**
     * Sets one element's imaginary part, by flat index, keeping its
     * real part.
     * \param[in] i Flat index, `0`-`8`.
     * \param[in] a New imaginary part.
     */
    void setIm(int i, double a) { e_[i] = complex<double>(e_[i].real(), a); }
    /**
     * Sets one element's imaginary part, by row/column, keeping its
     * real part.
     * \param[in] i Row index, `0`-`2`.
     * \param[in] j Column index, `0`-`2`.
     * \param[in] a New imaginary part.
     */
    void setIm(int i, int j, double a) {
        e_[j + kN * i] = complex<double>(e_[j + kN * i].real(), a);
    }

    /**
     * Sets one element, by flat index.
     * \param[in] i Flat index, `0`-`8`.
     * \param[in] a New value.
     */
    void set(int i, complex<double> a) { e_[i] = a; }
    /**
     * Sets one element, by row/column.
     * \param[in] i Row index, `0`-`2`.
     * \param[in] j Column index, `0`-`2`.
     * \param[in] a New value.
     */
    void set(int i, int j, complex<double> a) { e_[j + kN * i] = a; }

    /**
     * Returns one element, by flat index.
     * \param[in] i Flat index, `0`-`8`.
     * \return The element's value.
     */
    complex<double> get(int i) const { return e_[i]; }
    /**
     * Returns one element, by row/column.
     * \param[in] i Row index, `0`-`2`.
     * \param[in] j Column index, `0`-`2`.
     * \return The element's value.
     */
    complex<double> get(int i, int j) const { return e_[j + kN * i]; }

    /**
     * Returns one element's real part, by flat index.
     * \param[in] i Flat index, `0`-`8`.
     * \return The element's real part.
     */
    double getRe(int i) const { return e_[i].real(); }
    /**
     * Returns one element's imaginary part, by flat index.
     * \param[in] i Flat index, `0`-`8`.
     * \return The element's imaginary part.
     */
    double getIm(int i) const { return e_[i].imag(); }

    /**
     * Returns this matrix's side length.
     * \return `3`.
     */
    int getNDim() const { return kN; }
    /**
     * Returns this matrix's total element count.
     * \return `9`.
     */
    int getNN() const { return kNN; }

    /**
     * Serializes this matrix's 9 elements as a single-line,
     * space-separated string, column-major (`Re(0,0) Im(0,0) Re(1,0)
     * Im(1,0) Re(2,0) Im(2,0) Re(0,1) ...`), at 15-digit precision --
     * the format writeSU3Matrices()/writeWilsonLines()'s text output
     * uses.
     * \return The serialized string.
     */
    std::string MatrixToString();

    /**
     * Replaces this matrix with \f$\exp(tA)\f$ (\f$A\f$ = this
     * matrix), via a scaled-and-squared Pade approximant: finds a
     * scale \f$2^{-s}\f$ making \f$\|t A \cdot 2^{-s}\|_\infty <
     * 1/2\f$, evaluates the order-\p p Pade approximant there via
     * Horner's method (inverting the denominator through
     * cofactorInverse()), then squares the result \f$s\f$ times.
     * \param[in] t Scalar multiplying this matrix before
     * exponentiating.
     * \param[in] p Pade approximant order; must be at least 6 (exits
     * with an error otherwise). Orders above 6 use a generated
     * coefficient recurrence rather than the hardcoded low-order
     * values.
     * \return `*this`, now \f$\exp(tA)\f$.
     */
    Matrix &expm(double t = 1.0, const int p = 6);

    /**
     * Computes the coefficients of \f$\exp(iQ)\f$'s expansion in the
     * identity and the eight SU(3) fundamental generators, for
     * \f$Q = Q^a t^a\f$ a traceless Hermitian matrix given by its eight
     * real generator coefficients -- exactly, via the matrix's
     * characteristic-polynomial (Cayley-Hamilton) solution, rather
     * than a series/Pade approximation. Allocation-free, for hot paths
     * (e.g. GaugeFix.cpp's \c expGaugeRotationSU3, which combines
     * these coefficients into the actual exponential matrix). Any
     * coefficient that comes out NaN (a 0/0 case that can occur in the
     * very-low-density region) is replaced by `0`, contributing only
     * the identity/vacuum piece.
     * \param[in] Q The eight real coefficients \f$Q^1,\ldots,Q^8\f$ of
     * \f$Q\f$ in the generator basis (see Group).
     * \param[out] out The nine coefficients: `out[0]` for the identity,
     * `out[1..8]` for \f$t^1,\ldots,t^8\f$.
     */
    void expmCoeff(const double *Q, complex<double> out[9]) const;

    /**
     * Computes this matrix's determinant via the explicit
     * \f$3\times3\f$ cofactor expansion.
     * \return \f$\det(A)\f$.
     */
    complex<double> det() const;
    /**
     * Computes this matrix's trace.
     * \return \f$\mathrm{Tr}(A)\f$.
     */
    complex<double> trace() const;

    /**
     * Read-only element access, by flat index.
     * \param[in] i Flat index, `0`-`8`.
     * \return The element's value.
     */
    std::complex<double> operator()(const int i) const { return e_[i]; }
    /**
     * Read-only element access, by row/column.
     * \param[in] i Row index, `0`-`2`.
     * \param[in] j Column index, `0`-`2`.
     * \return The element's value.
     */
    std::complex<double> operator()(const int i, const int j) const {
        return e_[j + kN * i];
    }

    /**
     * Adds another matrix to this one, elementwise, in place.
     * \param[in] a Matrix to add.
     * \return `*this`.
     */
    Matrix &operator+=(const Matrix &a) {
        for (int i = 0; i < kNN; ++i) e_[i] += a.e_[i];
        return *this;
    }

    /**
     * Subtracts another matrix from this one, elementwise, in place.
     * \param[in] a Matrix to subtract.
     * \return `*this`.
     */
    Matrix &operator-=(const Matrix &a) {
        for (int i = 0; i < kNN; ++i) e_[i] -= a.e_[i];
        return *this;
    }

    /**
     * Scales this matrix by a complex scalar, in place.
     * \param[in] a Scalar factor.
     * \return `*this`.
     */
    Matrix &operator*=(const complex<double> a) {
        for (int i = 0; i < kNN; ++i) e_[i] *= a;
        return *this;
    }

    /**
     * Divides this matrix by a complex scalar, in place.
     * \param[in] a Scalar divisor.
     * \return `*this`.
     */
    Matrix &operator/=(const complex<double> a) {
        for (int i = 0; i < kNN; ++i) e_[i] /= a;
        return *this;
    }

    /**
     * Computes half the sum of squared magnitudes of every element,
     * equal to \f$\mathrm{Tr}(A^\dagger A)/2\f$ for this matrix
     * \f$A\f$.
     * \return \f$\tfrac{1}{2}\sum_{ij}|A_{ij}|^2\f$.
     */
    double square() const {
        double tr = 0.0;
        for (int i = 0; i < kNN; ++i) {
            tr += e_[i].real() * e_[i].real() + e_[i].imag() * e_[i].imag();
        }
        return 0.5 * tr;
    }

    /**
     * Conjugate-transposes this matrix in place.
     * \return `*this`, now \f$A^\dagger\f$.
     */
    Matrix &conjg();
    /**
     * Computes a matrix product with a Hermitian-conjugated right
     * factor.
     * \param[in] a Left factor.
     * \param[in] b Factor to conjugate-transpose before multiplying.
     * \return \f$AB^\dagger\f$.
     */
    static Matrix prodABconj(const Matrix &a, const Matrix &b);
    /**
     * Computes a matrix product with a Hermitian-conjugated left
     * factor.
     * \param[in] a Factor to conjugate-transpose before multiplying.
     * \param[in] b Right factor.
     * \return \f$A^\dagger B\f$.
     */
    static Matrix prodAconjB(const Matrix &a, const Matrix &b);

    /**
     * Computes \f$\mathrm{Tr}(ab)\f$ without materializing the product.
     *
     * \note Despite being a member function, this does not read
     * `*this` at all -- it depends only on \p a and \p b. Call it on
     * any Matrix instance.
     * \param[in] a Left factor.
     * \param[in] b Right factor.
     * \return \f$\mathrm{Tr}(ab)\f$.
     */
    complex<double> traceOfProductOfMatrix(Matrix &a, Matrix &b) const;

    friend ostream &operator<<(ostream &os, const Matrix &p) {
        for (int i = 0; i < kN; ++i) {
            for (int j = 0; j < kN; ++j) os << p(i, j);
            if (i < kN - 1) os << std::endl;
        }
        return os;
    }
};

/**
 * Matrix sum \f$C = A+B\f$.
 * \param[in] a First matrix.
 * \param[in] b Second matrix.
 * \return \f$A+B\f$.
 */
Matrix operator+(const Matrix &a, const Matrix &b);
/**
 * Matrix difference \f$C = A-B\f$.
 * \param[in] a First matrix.
 * \param[in] b Second matrix.
 * \return \f$A-B\f$.
 */
Matrix operator-(const Matrix &a, const Matrix &b);

/**
 * Scales a matrix by a real scalar (scalar on the left).
 * \param[in] a Scalar factor.
 * \param[in] b Matrix to scale.
 * \return \f$aB\f$.
 */
Matrix operator*(const double a, const Matrix &b);
/**
 * Scales a matrix by a complex scalar (scalar on the left).
 * \param[in] a Scalar factor.
 * \param[in] b Matrix to scale.
 * \return \f$aB\f$.
 */
Matrix operator*(const std::complex<double> a, const Matrix &b);
/**
 * Scales a matrix by a real scalar (scalar on the right).
 * \param[in] a Matrix to scale.
 * \param[in] b Scalar factor.
 * \return \f$Ab\f$.
 */
Matrix operator*(const Matrix &a, const double b);
/**
 * Matrix product \f$C = AB\f$.
 * \param[in] a Left factor.
 * \param[in] b Right factor.
 * \return \f$AB\f$.
 */
Matrix operator*(const Matrix &a, const Matrix &b);

#endif  // SRC_MATRIX_H_

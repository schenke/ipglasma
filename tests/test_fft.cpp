#include "doctest.h"

#include <cmath>
#include <complex>
#include <vector>

#include "FFT.h"
#include "Matrix.h"

TEST_CASE("FFT::fftnComplexArray: forward+backward recovers the input") {
    const int nn[2] = {8, 8};
    const int ntot = nn[0] * nn[1];
    FFT fft(nn);

    std::vector<std::complex<double>> original(ntot);
    for (int i = 0; i < ntot; ++i) {
        original[i] = std::complex<double>(
            std::sin(0.3 * i + 1.0), std::cos(0.7 * i - 0.4));
    }
    std::vector<std::complex<double>> transformed(ntot);
    std::vector<std::complex<double>> recovered(ntot);

    std::complex<double> *dataIn = original.data();
    std::complex<double> *dataOut = transformed.data();
    std::complex<double> *dataBack = recovered.data();

    fft.fftnComplexArray(&dataIn, &dataOut, nn, /*isign=*/1, /*mDim=*/1);
    fft.fftnComplexArray(&dataOut, &dataBack, nn, /*isign=*/-1, /*mDim=*/1);

    for (int i = 0; i < ntot; ++i) {
        CHECK(std::abs(recovered[i] - original[i]) < 1e-10);
    }
}

TEST_CASE("FFT::fftnArray: forward+backward recovers the input (position-major)") {
    // fftnArray uses the opposite layout from fftnComplexArray: data[pos]
    // is itself a contiguous mDim-length array (one per lattice site),
    // matching how JIMWLK's noise arrays are laid out.
    const int nn[2] = {8, 8};
    const int ntot = nn[0] * nn[1];
    const int mDim = 2;
    FFT fft(nn);

    std::vector<std::complex<double>> originalFlat(ntot * mDim);
    std::vector<std::complex<double>> transformedFlat(ntot * mDim);
    std::vector<std::complex<double>> recoveredFlat(ntot * mDim);
    std::vector<std::complex<double> *> originalPtrs(ntot);
    std::vector<std::complex<double> *> transformedPtrs(ntot);
    std::vector<std::complex<double> *> recoveredPtrs(ntot);
    for (int pos = 0; pos < ntot; ++pos) {
        originalPtrs[pos] = &originalFlat[pos * mDim];
        transformedPtrs[pos] = &transformedFlat[pos * mDim];
        recoveredPtrs[pos] = &recoveredFlat[pos * mDim];
        for (int k = 0; k < mDim; ++k) {
            originalPtrs[pos][k] = std::complex<double>(
                std::sin(0.2 * pos + k + 0.5), std::cos(0.4 * pos - k));
        }
    }

    fft.fftnArray(originalPtrs.data(), transformedPtrs.data(), nn, 1, mDim);
    fft.fftnArray(transformedPtrs.data(), recoveredPtrs.data(), nn, -1, mDim);

    for (int pos = 0; pos < ntot; ++pos) {
        for (int k = 0; k < mDim; ++k) {
            CHECK(std::abs(recoveredPtrs[pos][k] - originalPtrs[pos][k]) < 1e-10);
        }
    }
}

TEST_CASE("FFT::fftn<Matrix>: forward+backward recovers the input") {
    const int nn[2] = {8, 8};
    const int ntot = nn[0] * nn[1];
    FFT fft(nn);

    std::vector<Matrix> original(ntot, Matrix(Matrix::noInit));
    std::vector<Matrix> transformed(ntot, Matrix(Matrix::noInit));
    std::vector<Matrix> recovered(ntot, Matrix(Matrix::noInit));
    for (int pos = 0; pos < ntot; ++pos) {
        for (int k = 0; k < 9; ++k) {
            original[pos].set(
                k, std::complex<double>(
                       std::sin(0.1 * pos + k), std::cos(0.2 * pos - k)));
        }
    }

    std::vector<Matrix *> originalPtrs(ntot), transformedPtrs(ntot),
        recoveredPtrs(ntot);
    for (int pos = 0; pos < ntot; ++pos) {
        originalPtrs[pos] = &original[pos];
        transformedPtrs[pos] = &transformed[pos];
        recoveredPtrs[pos] = &recovered[pos];
    }

    fft.fftn(originalPtrs.data(), transformedPtrs.data(), nn, 1);
    fft.fftn(transformedPtrs.data(), recoveredPtrs.data(), nn, -1);

    for (int pos = 0; pos < ntot; ++pos) {
        for (int k = 0; k < 9; ++k) {
            CHECK(std::abs(recovered[pos].get(k) - original[pos].get(k)) < 1e-10);
        }
    }
}

TEST_CASE("FFT::fftnVector: forward+backward recovers the input") {
    const int nn[2] = {8, 8};
    const int ntot = nn[0] * nn[1];
    const int mDim = 3;
    FFT fft(nn);

    std::vector<std::vector<std::complex<double>>> original(
        ntot, std::vector<std::complex<double>>(mDim));
    std::vector<std::vector<std::complex<double>>> transformed(
        ntot, std::vector<std::complex<double>>(mDim));
    std::vector<std::vector<std::complex<double>>> recovered(
        ntot, std::vector<std::complex<double>>(mDim));
    for (int pos = 0; pos < ntot; ++pos) {
        for (int k = 0; k < mDim; ++k) {
            original[pos][k] = std::complex<double>(
                std::sin(0.15 * pos + k), std::cos(0.25 * pos - k));
        }
    }

    std::vector<std::vector<std::complex<double>> *> originalPtrs(ntot),
        transformedPtrs(ntot), recoveredPtrs(ntot);
    for (int pos = 0; pos < ntot; ++pos) {
        originalPtrs[pos] = &original[pos];
        transformedPtrs[pos] = &transformed[pos];
        recoveredPtrs[pos] = &recovered[pos];
    }

    fft.fftnVector(originalPtrs.data(), transformedPtrs.data(), nn, 1);
    fft.fftnVector(transformedPtrs.data(), recoveredPtrs.data(), nn, -1);

    for (int pos = 0; pos < ntot; ++pos) {
        for (int k = 0; k < mDim; ++k) {
            CHECK(std::abs(recovered[pos][k] - original[pos][k]) < 1e-10);
        }
    }
}

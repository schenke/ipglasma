#include "doctest.h"

#include <cmath>
#include <complex>
#include <vector>

#include "FFT.h"

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

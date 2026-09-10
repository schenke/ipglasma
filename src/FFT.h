// FFT.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_FFT_H_
#define SRC_FFT_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fftw3.h>

#include <unistd.h>

#include <algorithm>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

using std::complex;
using std::vector;

template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
    // assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(
        a.begin(), a.end(), b.begin(), std::back_inserter(result),
        std::plus<T>());
    return result;
}

template <typename T>
std::vector<T> operator*(
    const std::vector<T> &a, const std::complex<double> &b) {
    std::vector<T> result;
    result = a;
    int size = a.size();
    for (int i = 0; i < size; i++) {
        result.at(i) = b * a.at(i);
    }
    return result;
}

template <typename T>
std::vector<T> operator/(const std::vector<T> &a, const double b) {
    std::vector<T> result;
    result = a;
    int size = a.size();
    for (int i = 0; i < size; i++) {
        result.at(i) = a.at(i) / b;
    }
    return result;
}

// FFTW planner flag. FFTW_MEASURE benchmarks candidate algorithms at plan
// time by actually timing them, which is sensitive to whatever else the
// machine is doing at that moment; on separate invocations of the same
// binary, for the same problem size, this can and does pick different
// algorithms with different floating-point rounding, making event output
// non-reproducible at fixed seed even with no other source of randomness.
// Building with -DIPGLASMA_DETERMINISTIC_FFT=ON (see CMakeLists.txt) fixes
// this via an on-disk wisdom cache (see the FFT constructor below): the
// first run in a given directory measures and records a plan, and every run
// after that imports it and reuses it verbatim instead of measuring again,
// making the FFT step -- and hence the whole event -- bit-reproducible. Off
// by default, since it makes the first run in a directory take as long as
// FFTW_MEASURE always does today, in exchange for every later run there
// being both faster to plan and reproducible.
#define IPG_FFTW_PLAN_FLAG FFTW_MEASURE

class FFT {
  private:
    fftw_complex *input, *output;
    fftw_complex *inputMany, *outputMany;
    fftw_plan p_, pback_;

  public:
    // Largest number of planes any batched transform (fftn's 9 Matrix
    // components, fftnArray's up to 2*(Nc^2-1)=16 noise components) packs
    // into inputMany/outputMany at once.
    static constexpr int kMaxBatchDim = 16;

    // Constructor.
    FFT(const int nn[]) {
        //      if(fftw_init_threads()==0)
        //  cerr << "Error initializing multi-threaded fftw." << endl;
        // fftw_plan_with_nthreads(omp_get_max_threads());
        input =
            (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nn[0] * nn[1]);
        output =
            (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nn[0] * nn[1]);
#ifdef IPGLASMA_DETERMINISTIC_FFT
        // Pin the plan across runs: import any wisdom already on disk before
        // planning, so a previously chosen algorithm is reused verbatim
        // instead of FFTW deciding again (and possibly differently -- see
        // the comment on IPG_FFTW_PLAN_FLAG above).
        const char *wisdomFile = "ipglasma_fftw_wisdom.dat";
        fftw_import_wisdom_from_filename(wisdomFile);
#endif
        p_ = fftw_plan_dft_2d(
            nn[0], nn[1], input, output, FFTW_FORWARD, IPG_FFTW_PLAN_FLAG);
        pback_ = fftw_plan_dft_2d(
            nn[0], nn[1], input, output, FFTW_BACKWARD, IPG_FFTW_PLAN_FLAG);
#ifdef IPGLASMA_DETERMINISTIC_FFT
        {
            // Persist whatever wisdom now exists (including any plan just
            // created above) so future runs -- and other concurrently
            // starting MPI ranks -- converge onto the same plan. Write to a
            // per-process temp file and rename into place atomically so a
            // concurrent writer can never observe a partially written file.
            std::string tmpFile =
                std::string(wisdomFile) + ".tmp." + std::to_string(getpid());
            if (fftw_export_wisdom_to_filename(tmpFile.c_str())) {
                std::rename(tmpFile.c_str(), wisdomFile);
            } else {
                std::remove(tmpFile.c_str());
            }
        }
#endif
        inputMany = (fftw_complex *)fftw_malloc(
            sizeof(fftw_complex) * nn[0] * nn[1] * kMaxBatchDim);
        outputMany = (fftw_complex *)fftw_malloc(
            sizeof(fftw_complex) * nn[0] * nn[1] * kMaxBatchDim);
    };
    // Destructor
    ~FFT() {
        fftw_destroy_plan(p_);
        fftw_destroy_plan(pback_);
        fftw_free(input);
        fftw_free(output);
        fftw_free(inputMany);
        fftw_free(outputMany);
        // fftw_cleanup_threads();
    };
    void fftnVector(
        vector<complex<double>> **data, vector<complex<double>> **outdata,
        const int nn[], const int isign);
    void fftnArray(
        complex<double> **data, complex<double> **outdata, const int nn[],
        const int isign, const int mDim);

    template <class T>
    void fftn(T **data, T **outdata, const int nn[], const int isign);

    void fftnComplex(
        complex<double> *data, complex<double> *outdata, const int nn[],
        const int isign);
};

#endif  // SRC_FFT_H_

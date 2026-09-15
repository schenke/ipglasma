// FFT.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_FFT_H_
#define SRC_FFT_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fftw3.h>
#include <unistd.h>

#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "PrettyOstream.h"

using std::complex;
using std::vector;

/**
 * \def IPG_FFTW_PLAN_FLAG
 * FFTW planner flag used for every plan FFT creates.
 *
 * FFTW_MEASURE benchmarks candidate algorithms at plan time by actually
 * timing them, which is sensitive to whatever else the machine is doing
 * at that moment; on separate invocations of the same binary, for the
 * same problem size, this can and does pick different algorithms with
 * different floating-point rounding, making event output
 * non-reproducible at fixed seed even with no other source of
 * randomness. Building with `-DIPGLASMA_DETERMINISTIC_FFT=ON` (see
 * CMakeLists.txt) fixes this via an on-disk wisdom cache (see the FFT
 * constructor below): the first run in a given directory measures and
 * records a plan, and every run after that imports it and reuses it
 * verbatim instead of measuring again, making the FFT step -- and hence
 * the whole event -- bit-reproducible. Off by default, since it makes
 * the first run in a directory take as long as FFTW_MEASURE always does
 * today, in exchange for every later run being both faster to plan and
 * reproducible.
 */
#define IPG_FFTW_PLAN_FLAG FFTW_MEASURE

/**
 * Thin RAII wrapper around FFTW plans/scratch buffers for the 2D
 * complex-to-complex transforms used throughout the classical
 * Yang-Mills evolution, JIMWLK, and Coulomb gauge fixing.
 *
 * One FFT instance owns exactly one pair of forward/backward plans,
 * sized for a fixed transverse lattice (the `nn` passed to the
 * constructor), plus scratch buffers reused by every transform call.
 * All transform methods take an `isign` of `+1` (forward) or `-1`
 * (inverse, normalized by \f$1/N_0 N_1\f$ so a forward+inverse round
 * trip is the identity). Every method shares the same underlying
 * "shift-by-half-a-lattice" convention: input/output are indexed as a
 * function of \f$-N/2\f$ to \f$N/2\f$, not \f$0\f$ to \f$N\f$, matching
 * how the rest of the lattice is addressed.
 */
class FFT {
  private:
    /// Log sink for error/warning messages (e.g. an oversized \p mDim).
    PrettyOstream messager_;
    /// Single-plane scratch buffers used by fftnVector(), sized
    /// `nn[0]*nn[1]`.
    fftw_complex *input, *output;
    /// Batched, multi-plane scratch buffers used by fftn()/fftnArray()/
    /// fftnComplexArray(), sized `nn[0]*nn[1]*kMaxBatchDim`.
    fftw_complex *inputMany, *outputMany;
    /// Forward (\c p_) and backward (\c pback_) FFTW plans for this
    /// instance's fixed lattice size.
    fftw_plan p_, pback_;

  public:
    /// Largest number of planes any batched transform (fftn()'s 9
    /// Matrix components, fftnArray()'s up to \f$2(N_c^2-1)=16\f$
    /// noise components) packs into inputMany/outputMany at once.
    static constexpr int kMaxBatchDim = 16;

    /**
     * Allocates FFTW buffers and creates the forward/backward plans for
     * a fixed transverse lattice size.
     *
     * When built with `-DIPGLASMA_DETERMINISTIC_FFT=ON`, this also
     * imports any FFTW wisdom already on disk before planning (so a
     * previously chosen algorithm is reused verbatim instead of FFTW
     * possibly deciding differently -- see \ref IPG_FFTW_PLAN_FLAG) and
     * persists whatever wisdom exists afterwards, written to a
     * per-process temp file and renamed into place atomically so a
     * concurrent writer (e.g. another MPI rank starting at the same
     * time) can never observe a partially written file.
     * \param[in] nn Two-element array `{N0, N1}`, the transverse
     * lattice dimensions every transform on this instance will use.
     */
    FFT(const int nn[]) {
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
    /**
     * Destroys the plans and frees every FFTW buffer this instance
     * owns.
     */
    ~FFT() {
        fftw_destroy_plan(p_);
        fftw_destroy_plan(pback_);
        fftw_free(input);
        fftw_free(output);
        fftw_free(inputMany);
        fftw_free(outputMany);
    };

    // Owns raw FFTW buffers/plans freed in the destructor; default copies
    // would double-free them, so disable copying (nothing needs it).
    FFT(const FFT &) = delete;
    FFT &operator=(const FFT &) = delete;

    /**
     * Transforms one 2-component-per-site vector field in place.
     *
     * Used by JIMWLK::initializeK() to transform the two-component
     * transverse momentum kernel into position space once, at
     * construction. Handles the historical four-quadrant shift
     * explicitly (see fftn() for the equivalent, faster
     * checkerboard-sign approach used by every other transform here).
     * \param[in] data Array of `nn[0]*nn[1]` pointers, one per lattice
     * site, each to a same-length `vector<complex<double>>`; every
     * vector index is transformed independently.
     * \param[out] outdata Array of `nn[0]*nn[1]` pointers, one per
     * lattice site, receiving the transformed values; may alias \p
     * data for an in-place transform.
     * \param[in] nn Two-element array `{N0, N1}`, the transverse
     * lattice dimensions (must match what this instance was
     * constructed with).
     * \param[in] isign `+1` for a forward transform, `-1` for an
     * inverse transform (normalized by \f$1/N_0 N_1\f$).
     */
    void fftnVector(
        vector<complex<double>> **data, vector<complex<double>> **outdata,
        const int nn[], const int isign);
    /**
     * Transforms up to \p mDim independent single-complex-value planes
     * at once, batched and parallelized over planes.
     *
     * Used by JIMWLK::evolutionStep() for its noise fields (`xi2_`,
     * `CKxi_`). Each of the \p mDim planes is packed into its own slice
     * of the shared \c inputMany/outputMany scratch and executed
     * concurrently via FFTW's thread-safe new-array interface, using
     * the same checkerboard-sign trick as fftn() (a half-lattice input/
     * output shift is equivalent to multiplying by \f$(-1)^{i+j}\f$,
     * plus a global sign on the output) instead of an explicit
     * four-quadrant index remap.
     * \param[in] data Array of `nn[0]*nn[1]` pointers, one per lattice
     * site, each to an \p mDim-length array of values (one value per
     * plane at that site).
     * \param[out] outdata Array of `nn[0]*nn[1]` pointers, one per
     * lattice site, receiving the transformed values; may alias \p
     * data for an in-place transform.
     * \param[in] nn Two-element array `{N0, N1}`, the transverse
     * lattice dimensions (must match what this instance was
     * constructed with).
     * \param[in] isign `+1` for a forward transform, `-1` for an
     * inverse transform (normalized by \f$1/N_0 N_1\f$).
     * \param[in] mDim Number of planes to transform; exits with an
     * error if greater than kMaxBatchDim.
     */
    void fftnArray(
        complex<double> **data, complex<double> **outdata, const int nn[],
        const int isign, const int mDim);

    /**
     * Transforms every component of a per-site object of class \p T at
     * once, batched and parallelized over components.
     *
     * \p T must provide `int getNDim() const` (its side length, e.g. 3
     * for an SU(3) Matrix) and `complex<double> *data()` (its
     * `getNDim()^2` components, contiguous); \f$getNDim()^2\f$ planes
     * are packed into the shared \c inputMany/outputMany scratch and
     * executed concurrently via FFTW's thread-safe new-array interface,
     * using the same checkerboard-sign trick as fftnArray(). Used with
     * `T = Matrix` for every SU(3)-matrix-field FFT in the classical
     * Yang-Mills evolution (Evolution.cpp), JIMWLK (`VxsiVx_`/
     * `VxsiVy_`), and Coulomb gauge fixing (GaugeFix.cpp).
     * \param[in] data Array of `nn[0]*nn[1]` pointers, one per lattice
     * site, each to an object of class \p T.
     * \param[out] outdata Array of `nn[0]*nn[1]` pointers, one per
     * lattice site, receiving the transformed values; may alias \p
     * data for an in-place transform.
     * \param[in] nn Two-element array `{N0, N1}`, the transverse
     * lattice dimensions (must match what this instance was
     * constructed with).
     * \param[in] isign `+1` for a forward transform, `-1` for an
     * inverse transform (normalized by \f$1/N_0 N_1\f$).
     */
    template <class T>
    void fftn(T **data, T **outdata, const int nn[], const int isign);

    /**
     * Transforms \p mDim independent 2D planes at once, each plane
     * already a standalone contiguous `nn[0]*nn[1]`-length array.
     *
     * Used by Init.cpp's Wilson-line construction to transform the
     * `Nc^2-1` color-charge Fourier coefficients when solving the
     * 2D Poisson equation for the classical field (forward transform,
     * apply the momentum kernel, inverse transform back). Batched and
     * parallelized identically to fftnArray()/fftn(); unlike
     * fftnArray() (whose planes are gathered from a per-site,
     * pos-major array), packing here is a per-plane copy, since each
     * plane is already contiguous.
     * \param[in] data Array of \p mDim pointers, each to a standalone
     * contiguous `nn[0]*nn[1]`-length plane.
     * \param[out] outdata Array of \p mDim pointers, each to a
     * standalone contiguous `nn[0]*nn[1]`-length plane receiving the
     * transformed values; may alias \p data for an in-place transform.
     * \param[in] nn Two-element array `{N0, N1}`, the transverse
     * lattice dimensions (must match what this instance was
     * constructed with).
     * \param[in] isign `+1` for a forward transform, `-1` for an
     * inverse transform (normalized by \f$1/N_0 N_1\f$).
     * \param[in] mDim Number of planes to transform; exits with an
     * error if greater than kMaxBatchDim.
     */
    void fftnComplexArray(
        complex<double> **data, complex<double> **outdata, const int nn[],
        const int isign, const int mDim);
};

#endif  // SRC_FFT_H_

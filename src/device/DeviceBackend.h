// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#ifndef SRC_DEVICE_DEVICEBACKEND_H_
#define SRC_DEVICE_DEVICEBACKEND_H_

#include <cstddef>
#include <initializer_list>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

/**
 * The device layer: the JIMWLK evolution and the classical Yang-Mills
 * leapfrog on a GPU.
 *
 * The physics is in the kernels of ipg_kernels.h, one source for all
 * backends. A Backend only owns buffers of real numbers on its device,
 * copies them from and to the host, and runs a kernel over an index range.
 * The drivers (DeviceJimwlk, DeviceYangMills) are written against this
 * interface, so they are the same code for every backend:
 *
 *   - `host`   the kernels on the CPU (with OpenMP), always built: the
 *              reference of the device layer, in double or single
 *              precision,
 *   - `cuda`   NVIDIA GPUs (CudaBackend.cu, `-DIPGLASMA_CUDA=ON`), with
 *              cuFFT,
 *   - `opencl` OpenCL devices, among them the GPUs of Apple M processors
 *              (OpenCLBackend.cpp, `-DIPGLASMA_OPENCL=ON`).
 *
 * The backend of a run is chosen with the environment variable
 * `IPGLASMA_DEVICE` (see deviceConfigFromEnvironment()); without it the
 * code runs as before, on the CPU, and this layer is not used.
 */
namespace ipgdev {

/// The kernels of ipg_kernels.h, in the order of its IPG_KERNEL_LIST.
enum Kernel {
    kCopy = 0,
    kRealToComplex,
    kComplexToReal,
    kScaleBySite,
    kCheckerScale,
    kFftPass,
    kGaussNoise,
    kJimwlkKernelNoise,
    kJimwlkRotate,
    kJimwlkKernelMatrix,
    kJimwlkUpdate,
    kCymPi,
    kCymE,
    kCymPhi,
    kCymU,
    kGfDivergence,
    kGfUpdate,
    kGfTransform,
    kWilsonCharge,
    kWilsonExponent,
    kNumKernels
};

/// The name of a kernel in ipg_kernels.h, e.g. "ipg_k_copy".
const char *kernelName(Kernel kernel);

/// The parameters of a kernel launch (ipg_params of ipg_kernels.h, with
/// the real numbers in double precision; a backend converts them).
struct Params {
    int n = 0;
    int ncomp = 1;
    int i0 = 0;
    int i1 = 0;
    int i2 = 0;
    int i3 = 0;
    int i4 = 0;
    int i5 = 0;
    double r0 = 0.;
    double r1 = 0.;
    double r2 = 0.;
    double r3 = 0.;
};

/// Handle of a buffer of real numbers on the device.
using Buffer = int;

/// The radices of the passes of the device FFT of length \p n (8, 4 and 2
/// for the powers of two, then the odd primes). Empty if \p n has a prime
/// factor above \p maxRadix, for which the device FFT is not available.
std::vector<int> fftRadices(int n, int maxRadix = 32);

class Backend {
  public:
    virtual ~Backend() = default;

    /// `host`, `cuda` or `opencl`.
    virtual std::string name() const = 0;
    /// The device and the precision, for the log.
    virtual std::string description() const = 0;
    /// Whether the device computes in single precision.
    virtual bool singlePrecision() const = 0;

    /// A buffer of \p nReals real numbers, zero initialized.
    virtual Buffer createBuffer(std::size_t nReals) = 0;
    virtual void releaseBuffer(Buffer buffer) = 0;
    virtual std::size_t bufferSize(Buffer buffer) const = 0;
    /// Copies \p nReals doubles of the host to the start of the buffer
    /// (converted to the precision of the device).
    virtual void upload(
        Buffer buffer, const double *source, std::size_t nReals) = 0;
    /// Copies the first \p nReals reals of the buffer to the host.
    virtual void download(
        Buffer buffer, double *destination, std::size_t nReals) = 0;
    /// Exchanges the contents of two buffers of the same size (the handles
    /// keep their meaning, no data is copied).
    virtual void swapBuffers(Buffer a, Buffer b) = 0;

    /// Runs \p kernel for the indices 0 .. \p count - 1 with up to six
    /// buffers b0 .. b5.
    virtual void launch(
        Kernel kernel, std::size_t count, const Params &params,
        std::initializer_list<Buffer> buffers) = 0;

    /// Waits until all launches have finished (for timings).
    virtual void finish() {}

    /**
     * The centered two-dimensional FFT of FFT::fftn() (FFT.cpp), in place,
     * of a field with \p ncomp complex components per site on an
     * \p n x \p n lattice: forward for \p isign = 1, backward and
     * normalized by 1/n^2 for -1. The default is the mixed-radix FFT of
     * the kernel ipg_k_fft_pass; a backend with an FFT library overrides
     * it.
     */
    virtual void centeredFft(Buffer field, int n, int ncomp, int isign);

    /// Whether centeredFft() can transform a lattice of side \p n.
    virtual bool fftAvailable(int n) const {
        return n % 2 == 0 && !fftRadices(n).empty();
    }

  protected:
    /// A scratch buffer of the size of \p like, kept for later calls.
    Buffer scratchLike(Buffer like);
    /// The table (cos, sin)(2 pi t / n), t < n, on the device.
    Buffer twiddles(int n);
    /// Releases the scratch buffers and tables; a backend calls it in its
    /// destructor, while its device is still alive.
    void releaseShared();

  private:
    std::map<std::size_t, Buffer> scratch_;
    std::map<int, Buffer> twiddles_;
};

/// The choice of a device for a run.
struct Config {
    /// Empty: no device, the CPU code. Otherwise `host`, `cuda`, `opencl`.
    std::string backend;
    /// 0: the default of the device (double where it has it), 1: single
    /// precision, 2: double precision (an error where the device has none).
    int precision = 0;
    /// Whether the JIMWLK evolution / the Yang-Mills evolution run on the
    /// device.
    bool jimwlk = true;
    bool yangMills = true;
    /// Whether the Wilson lines of the nuclei (Init::setV()) and the
    /// Coulomb gauge fixing (GaugeFix::fftChi()) run on the device.
    bool wilsonLines = true;
    bool gaugeFix = true;
    /// Whether every kernel launch, FFT and copy is timed
    /// (IPGLASMA_DEVICE_PROFILE): the device then waits after each of
    /// them, so the run gets slower, and the times go to the profile of
    /// the run (IPGLASMA_PROFILE) as device.kernel.<name>, device.fft,
    /// device.upload and device.download.
    bool profileKernels = false;
    /// Where the Gaussian numbers of the color charges and the noise of
    /// the JIMWLK evolution are drawn: on the
    /// host from the random numbers of the run, exactly as in a CPU run
    /// (false), or on the device with a counter-based generator (true):
    /// faster, statistically equivalent, but not the same event as the
    /// CPU run with this seed.
    bool noiseOnDevice = false;
    bool enabled() const { return !backend.empty(); }
};

/**
 * Reads the device of the run from the environment:
 *   - `IPGLASMA_DEVICE` = `cuda`, `opencl`, `host`; unset, empty, `cpu`,
 *     `none` or `off` for the CPU code,
 *   - `IPGLASMA_DEVICE_PRECISION` = `single` or `double` (default: double
 *     if the device has it, e.g. not on the GPUs of Apple M processors),
 *   - `IPGLASMA_DEVICE_PROFILE` = 1 times every launch (see
 *     Config::profileKernels),
 *   - `IPGLASMA_DEVICE_STAGES` = a comma separated list of `wilsonlines`,
 *     `jimwlk`, `yangmills` and `gaugefix` (default: all),
 *   - `IPGLASMA_DEVICE_NOISE` = `host` (default) or `device`, see
 *     Config::noiseOnDevice.
 * Throws std::runtime_error for other values.
 */
Config configFromEnvironment();

/// The backends this executable was built with.
std::vector<std::string> availableBackends();

/**
 * Creates the backend of \p config. Throws std::runtime_error if it is not
 * built in, has no device, or the device lacks the requested precision.
 */
std::unique_ptr<Backend> createBackend(const Config &config);

/**
 * The backend of the run: created once from configFromEnvironment() and
 * shared by the events. nullptr without a device.
 */
Backend *runBackend();
/// The configuration runBackend() was created from.
const Config &runConfig();
/// Replaces the backend of the run (for tests); an empty Config::backend
/// removes it.
void setRunBackend(const Config &config);

// factories of the backends, defined where they are built
std::unique_ptr<Backend> createHostBackendDouble();
std::unique_ptr<Backend> createHostBackendSingle();
#ifdef IPGLASMA_CUDA
std::unique_ptr<Backend> createCudaBackend(int precision);
#endif
#ifdef IPGLASMA_OPENCL
std::unique_ptr<Backend> createOpenCLBackend(int precision);
#endif

}  // namespace ipgdev

#endif  // SRC_DEVICE_DEVICEBACKEND_H_

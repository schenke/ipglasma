// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include "DeviceBackend.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "../Instrumentation.h"

namespace ipgdev {

std::vector<int> fftRadices(int n, int maxRadix) {
    std::vector<int> radices;
    if (n < 1) return radices;
    int rest = n;
    int twos = 0;
    while (rest % 2 == 0) {
        rest /= 2;
        twos++;
    }
    // 2^twos as radices 8 and 4 where possible: fewer passes over the
    // memory for the same number of operations per element
    if (twos % 2 == 1 && twos >= 3) {
        radices.push_back(8);
        twos -= 3;
    }
    while (twos >= 2) {
        radices.push_back(4);
        twos -= 2;
    }
    if (twos == 1) radices.push_back(2);
    for (int prime = 3; rest > 1; prime += 2) {
        while (rest % prime == 0) {
            if (prime > maxRadix) return std::vector<int>();
            radices.push_back(prime);
            rest /= prime;
        }
        if (prime > maxRadix && rest > 1) return std::vector<int>();
    }
    return radices;
}

Buffer Backend::scratchLike(Buffer like) {
    const std::size_t size = bufferSize(like);
    auto found = scratch_.find(size);
    if (found != scratch_.end()) return found->second;
    const Buffer buffer = createBuffer(size);
    scratch_[size] = buffer;
    return buffer;
}

Buffer Backend::twiddles(int n) {
    auto found = twiddles_.find(n);
    if (found != twiddles_.end()) return found->second;
    std::vector<double> table(2 * static_cast<std::size_t>(n));
    for (int t = 0; t < n; t++) {
        const double angle = 2. * M_PI * t / n;
        table[2 * t] = std::cos(angle);
        table[2 * t + 1] = std::sin(angle);
    }
    // exact values where they are known, so that the transform of a real
    // symmetric input has no rounding in them
    table[0] = 1.;
    table[1] = 0.;
    if (n % 2 == 0) {
        table[n] = -1.;
        table[n + 1] = 0.;
    }
    if (n % 4 == 0) {
        table[n / 2] = 0.;
        table[n / 2 + 1] = 1.;
        table[3 * n / 2] = 0.;
        table[3 * n / 2 + 1] = -1.;
    }
    const Buffer buffer = createBuffer(table.size());
    upload(buffer, table.data(), table.size());
    twiddles_[n] = buffer;
    return buffer;
}

void Backend::releaseShared() {
    for (const auto &entry : scratch_) releaseBuffer(entry.second);
    for (const auto &entry : twiddles_) releaseBuffer(entry.second);
    scratch_.clear();
    twiddles_.clear();
}

namespace {

// the order of ipgdev::Kernel
const char *const kKernelNames[kNumKernels] = {
    "ipg_k_copy",
    "ipg_k_real_to_complex",
    "ipg_k_complex_to_real",
    "ipg_k_scale_by_site",
    "ipg_k_checker_scale",
    "ipg_k_fft_pass",
    "ipg_k_gauss_noise",
    "ipg_k_jimwlk_kernel_noise",
    "ipg_k_jimwlk_rotate",
    "ipg_k_jimwlk_kernel_matrix",
    "ipg_k_jimwlk_update",
    "ipg_k_cym_pi",
    "ipg_k_cym_e",
    "ipg_k_cym_phi",
    "ipg_k_cym_u",
    "ipg_k_gf_divergence",
    "ipg_k_gf_update",
    "ipg_k_gf_transform",
    "ipg_k_wilson_charge",
    "ipg_k_wilson_exponent",
};
static_assert(
    sizeof(kKernelNames) / sizeof(kKernelNames[0]) == kNumKernels,
    "one name per kernel");

}  // namespace

const char *kernelName(Kernel kernel) {
    const int k = static_cast<int>(kernel);
    return (k >= 0 && k < kNumKernels) ? kKernelNames[k] : "unknown";
}

void Backend::centeredFft(Buffer field, int n, int ncomp, int isign) {
    const std::vector<int> radices = fftRadices(n);
    if (radices.empty()) {
        std::ostringstream message;
        message << "the device FFT needs a lattice size with prime factors "
                   "up to 32 only; size = "
                << n << " has a larger one";
        throw std::runtime_error(message.str());
    }
    const std::size_t elements =
        static_cast<std::size_t>(n) * n * static_cast<std::size_t>(ncomp);
    const Buffer table = twiddles(n);
    const Buffer scratch = scratchLike(field);

    // S F S = G M F M (FFT.cpp): the checkerboard signs M go into the first
    // and the last pass, the global sign G and the normalization into the
    // last
    const double globalSign = (((n / 2 + n / 2) & 1) != 0) ? -1. : 1.;
    const double norm = (isign == -1) ? 1. / (static_cast<double>(n) * n) : 1.;
    const bool checker = (n % 2 == 0);
    if (!checker && n > 1) {
        throw std::runtime_error(
            "the centered device FFT needs an even lattice size");
    }

    const int passes = 2 * static_cast<int>(radices.size());
    int pass = 0;
    Buffer source = field;
    Buffer destination = scratch;
    for (int axis = 0; axis < 2; axis++) {
        int s = 1;
        for (const int radix : radices) {
            Params p;
            p.n = n;
            p.ncomp = ncomp;
            p.i0 = axis;
            p.i1 = radix;
            p.i2 = s;
            p.i3 = (isign == 1) ? 1 : -1;
            p.i4 = 0;
            p.r0 = 1.;
            if (pass == 0 && checker) p.i4 |= 1;
            if (pass == passes - 1) {
                if (checker) p.i4 |= 2;
                p.r0 = globalSign * norm;
            }
            // one element per butterfly of radix `radix`
            launch(
                kFftPass, elements / static_cast<std::size_t>(radix), p,
                {source, destination, table});
            std::swap(source, destination);
            s *= radix;
            pass++;
        }
    }
    // the result is in `source`; an even number of passes leaves it in
    // `field`
    if (source != field) swapBuffers(field, scratch);
}

namespace {

std::string lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return text;
}

std::string environment(const char *name) {
    const char *value = std::getenv(name);
    return value == nullptr ? std::string() : lower(std::string(value));
}

struct RunDevice {
    bool initialized = false;
    Config config;
    std::unique_ptr<Backend> backend;
};

RunDevice &runDevice() {
    static RunDevice device;
    return device;
}

}  // namespace

Config configFromEnvironment() {
    Config config;
    const std::string device = environment("IPGLASMA_DEVICE");
    if (device.empty() || device == "cpu" || device == "none" || device == "off"
        || device == "0") {
        return config;
    }
    if (device != "host" && device != "cuda" && device != "opencl") {
        throw std::runtime_error(
            "IPGLASMA_DEVICE must be cuda, opencl, host or cpu, not '" + device
            + "'");
    }
    config.backend = device;

    const std::string precision = environment("IPGLASMA_DEVICE_PRECISION");
    if (precision == "single" || precision == "float") {
        config.precision = 1;
    } else if (precision == "double") {
        config.precision = 2;
    } else if (!precision.empty()) {
        throw std::runtime_error(
            "IPGLASMA_DEVICE_PRECISION must be single or double, not '"
            + precision + "'");
    }

    const std::string stages = environment("IPGLASMA_DEVICE_STAGES");
    if (!stages.empty()) {
        config.jimwlk = false;
        config.yangMills = false;
        config.wilsonLines = false;
        config.gaugeFix = false;
        std::stringstream list(stages);
        std::string stage;
        while (std::getline(list, stage, ',')) {
            if (stage == "jimwlk") {
                config.jimwlk = true;
            } else if (stage == "yangmills" || stage == "cym") {
                config.yangMills = true;
            } else if (stage == "wilsonlines") {
                config.wilsonLines = true;
            } else if (stage == "gaugefix") {
                config.gaugeFix = true;
            } else if (!stage.empty()) {
                throw std::runtime_error(
                    "IPGLASMA_DEVICE_STAGES is a comma separated list of "
                    "wilsonlines, jimwlk, yangmills and gaugefix, not '"
                    + stage + "'");
            }
        }
    }

    const std::string profile = lower(environment("IPGLASMA_DEVICE_PROFILE"));
    config.profileKernels = !profile.empty() && profile != "0"
                            && profile != "false" && profile != "off"
                            && profile != "no";

    const std::string noise = environment("IPGLASMA_DEVICE_NOISE");
    if (noise == "device") {
        config.noiseOnDevice = true;
    } else if (!noise.empty() && noise != "host") {
        throw std::runtime_error(
            "IPGLASMA_DEVICE_NOISE must be host or device, not '" + noise
            + "'");
    }
    return config;
}

std::vector<std::string> availableBackends() {
    std::vector<std::string> backends;
    backends.reserve(3);
    backends.emplace_back("host");
#ifdef IPGLASMA_CUDA
    backends.emplace_back("cuda");
#endif
#ifdef IPGLASMA_OPENCL
    backends.emplace_back("opencl");
#endif
    return backends;
}

namespace {

/**
 * A backend that times every launch, FFT and copy of another
 * (IPGLASMA_DEVICE_PROFILE): it waits for the device before and after
 * each, so that a time is the work of that call alone, and adds it to the
 * profile of the run. The waiting removes the overlap of host and device,
 * so the run is slower than without.
 */
class ProfilingBackend : public Backend {
  public:
    explicit ProfilingBackend(std::unique_ptr<Backend> inner)
        : inner_(std::move(inner)) {}

    std::string name() const override { return inner_->name(); }
    std::string description() const override {
        return inner_->description() + ", every launch timed";
    }
    bool singlePrecision() const override { return inner_->singlePrecision(); }
    Buffer createBuffer(std::size_t nReals) override {
        return inner_->createBuffer(nReals);
    }
    void releaseBuffer(Buffer buffer) override {
        inner_->releaseBuffer(buffer);
    }
    std::size_t bufferSize(Buffer buffer) const override {
        return inner_->bufferSize(buffer);
    }
    void upload(Buffer buffer, const double *source, std::size_t nReals)
        override {
        Timer timer(*inner_, "device.upload");
        inner_->upload(buffer, source, nReals);
    }
    void download(Buffer buffer, double *destination, std::size_t nReals)
        override {
        Timer timer(*inner_, "device.download");
        inner_->download(buffer, destination, nReals);
    }
    void swapBuffers(Buffer a, Buffer b) override {
        inner_->swapBuffers(a, b);
    }
    void launch(
        Kernel kernel, std::size_t count, const Params &params,
        std::initializer_list<Buffer> buffers) override {
        // "ipg_k_" left out
        Timer timer(
            *inner_, std::string("device.kernel.") + (kernelName(kernel) + 6));
        inner_->launch(kernel, count, params, buffers);
    }
    void finish() override { inner_->finish(); }
    void centeredFft(Buffer field, int n, int ncomp, int isign) override {
        Timer timer(*inner_, "device.fft");
        inner_->centeredFft(field, n, ncomp, isign);
    }
    bool fftAvailable(int n) const override { return inner_->fftAvailable(n); }

  private:
    /// Times its scope, from an idle device to an idle device.
    class Timer {
      public:
        Timer(Backend &backend, std::string phase)
            : backend_(backend), phase_(std::move(phase)) {
            backend_.finish();
            start_ = std::chrono::steady_clock::now();
        }
        ~Timer() {
            backend_.finish();
            const std::chrono::duration<double> elapsed =
                std::chrono::steady_clock::now() - start_;
            ipg::Profiler::instance().add(phase_, elapsed.count());
        }
        Timer(const Timer &) = delete;
        Timer &operator=(const Timer &) = delete;

      private:
        Backend &backend_;
        const std::string phase_;
        std::chrono::steady_clock::time_point start_;
    };

    std::unique_ptr<Backend> inner_;
};

std::unique_ptr<Backend> createPlainBackend(const Config &config);

}  // namespace

std::unique_ptr<Backend> createBackend(const Config &config) {
    std::unique_ptr<Backend> backend = createPlainBackend(config);
    if (config.profileKernels) {
        backend.reset(new ProfilingBackend(std::move(backend)));
    }
    return backend;
}

namespace {

std::unique_ptr<Backend> createPlainBackend(const Config &config) {
    if (config.backend == "host") {
        return (config.precision == 1) ? createHostBackendSingle()
                                       : createHostBackendDouble();
    }
    if (config.backend == "cuda") {
#ifdef IPGLASMA_CUDA
        return createCudaBackend(config.precision);
#else
        throw std::runtime_error(
            "IPGLASMA_DEVICE=cuda, but this executable was built without "
            "CUDA (configure with -DIPGLASMA_CUDA=ON)");
#endif
    }
    if (config.backend == "opencl") {
#ifdef IPGLASMA_OPENCL
        return createOpenCLBackend(config.precision);
#else
        throw std::runtime_error(
            "IPGLASMA_DEVICE=opencl, but this executable was built without "
            "OpenCL (configure with -DIPGLASMA_OPENCL=ON)");
#endif
    }
    throw std::runtime_error("unknown device backend '" + config.backend + "'");
}

}  // namespace

Backend *runBackend() {
    RunDevice &device = runDevice();
    if (!device.initialized) {
        device.initialized = true;
        device.config = configFromEnvironment();
        if (device.config.enabled()) {
            device.backend = createBackend(device.config);
        }
    }
    return device.backend.get();
}

const Config &runConfig() {
    runBackend();
    return runDevice().config;
}

void setRunBackend(const Config &config) {
    RunDevice &device = runDevice();
    device.backend.reset();
    device.initialized = true;
    device.config = config;
    if (config.enabled()) device.backend = createBackend(config);
}

}  // namespace ipgdev

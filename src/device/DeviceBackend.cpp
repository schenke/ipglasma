// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include "DeviceBackend.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

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
            launch(kFftPass, elements, p, {source, destination, table});
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
        std::stringstream list(stages);
        std::string stage;
        while (std::getline(list, stage, ',')) {
            if (stage == "jimwlk") {
                config.jimwlk = true;
            } else if (stage == "yangmills" || stage == "cym") {
                config.yangMills = true;
            } else if (!stage.empty()) {
                throw std::runtime_error(
                    "IPGLASMA_DEVICE_STAGES is a comma separated list of "
                    "jimwlk and yangmills, not '"
                    + stage + "'");
            }
        }
    }

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

std::unique_ptr<Backend> createBackend(const Config &config) {
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

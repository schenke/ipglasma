// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

// The host backend: the kernels of ipg_kernels.h on the CPU. Included by
// HostBackendDouble.cpp and HostBackendSingle.cpp, which define
// IPG_HOST_REAL and IPG_HOST_FACTORY; not a header for other use.

#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "DeviceBackend.h"

namespace {

#define IPG_REAL IPG_HOST_REAL
#define IPG_FN static inline
#define IPG_GLOBAL
#define IPG_SQRT(x) std::sqrt(x)
#define IPG_SIN(x) std::sin(x)
#define IPG_COS(x) std::cos(x)
#define IPG_ACOS(x) std::acos(x)
#define IPG_LOG(x) std::log(x)
#define IPG_FABS(x) std::fabs(x)
// every kernel has the same arguments and uses some of them
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include "ipg_kernels.h"
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

class HostBackend : public ipgdev::Backend {
  public:
    ~HostBackend() override { releaseShared(); }

    std::string name() const override { return "host"; }
    std::string description() const override {
        return std::string("host CPU, ")
               + (singlePrecision() ? "single" : "double") + " precision";
    }
    bool singlePrecision() const override {
        return sizeof(ipg_real) == sizeof(float);
    }

    ipgdev::Buffer createBuffer(std::size_t nReals) override {
        for (std::size_t slot = 0; slot < buffers_.size(); slot++) {
            if (!used_[slot]) {
                buffers_[slot].assign(nReals, ipg_real(0));
                used_[slot] = true;
                return static_cast<ipgdev::Buffer>(slot);
            }
        }
        buffers_.emplace_back(nReals, ipg_real(0));
        used_.push_back(true);
        return static_cast<ipgdev::Buffer>(buffers_.size() - 1);
    }
    void releaseBuffer(ipgdev::Buffer buffer) override {
        std::vector<ipg_real>().swap(storage(buffer));
        used_[static_cast<std::size_t>(buffer)] = false;
    }
    std::size_t bufferSize(ipgdev::Buffer buffer) const override {
        return buffers_.at(static_cast<std::size_t>(buffer)).size();
    }
    void upload(ipgdev::Buffer buffer, const double *source, std::size_t nReals)
        override {
        std::vector<ipg_real> &data = storage(buffer);
        checkSize(data, nReals);
#pragma omp parallel for
        for (std::size_t k = 0; k < nReals; k++) {
            data[k] = static_cast<ipg_real>(source[k]);
        }
    }
    void download(
        ipgdev::Buffer buffer, double *destination,
        std::size_t nReals) override {
        const std::vector<ipg_real> &data = storage(buffer);
        checkSize(data, nReals);
#pragma omp parallel for
        for (std::size_t k = 0; k < nReals; k++) {
            destination[k] = static_cast<double>(data[k]);
        }
    }
    void swapBuffers(ipgdev::Buffer a, ipgdev::Buffer b) override {
        storage(a).swap(storage(b));
    }

    void launch(
        ipgdev::Kernel kernel, std::size_t count, const ipgdev::Params &params,
        std::initializer_list<ipgdev::Buffer> buffers) override {
        ipg_params p;
        p.n = params.n;
        p.ncomp = params.ncomp;
        p.i0 = params.i0;
        p.i1 = params.i1;
        p.i2 = params.i2;
        p.i3 = params.i3;
        p.i4 = params.i4;
        p.i5 = params.i5;
        p.r0 = static_cast<ipg_real>(params.r0);
        p.r1 = static_cast<ipg_real>(params.r1);
        p.r2 = static_cast<ipg_real>(params.r2);
        p.r3 = static_cast<ipg_real>(params.r3);
        ipg_real *b[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
        int slot = 0;
        for (const ipgdev::Buffer buffer : buffers) {
            if (slot >= 6) throw std::runtime_error("too many kernel buffers");
            b[slot++] = storage(buffer).data();
        }
        const int n = static_cast<int>(count);
        switch (kernel) {
#define IPG_HOST_CASE(id, function)                             \
    case ipgdev::id: {                                          \
        _Pragma("omp parallel for") for (int i = 0; i < n; i++) \
            function(i, p, b[0], b[1], b[2], b[3], b[4], b[5]); \
    } break;
            IPG_HOST_CASE(kCopy, ipg_k_copy)
            IPG_HOST_CASE(kRealToComplex, ipg_k_real_to_complex)
            IPG_HOST_CASE(kComplexToReal, ipg_k_complex_to_real)
            IPG_HOST_CASE(kScaleBySite, ipg_k_scale_by_site)
            IPG_HOST_CASE(kCheckerScale, ipg_k_checker_scale)
            IPG_HOST_CASE(kFftPass, ipg_k_fft_pass)
            IPG_HOST_CASE(kGaussNoise, ipg_k_gauss_noise)
            IPG_HOST_CASE(kJimwlkKernelNoise, ipg_k_jimwlk_kernel_noise)
            IPG_HOST_CASE(kJimwlkRotate, ipg_k_jimwlk_rotate)
            IPG_HOST_CASE(kJimwlkKernelMatrix, ipg_k_jimwlk_kernel_matrix)
            IPG_HOST_CASE(kJimwlkUpdate, ipg_k_jimwlk_update)
            IPG_HOST_CASE(kCymPi, ipg_k_cym_pi)
            IPG_HOST_CASE(kCymE, ipg_k_cym_e)
            IPG_HOST_CASE(kCymPhi, ipg_k_cym_phi)
            IPG_HOST_CASE(kCymU, ipg_k_cym_u)
            IPG_HOST_CASE(kGfDivergence, ipg_k_gf_divergence)
            IPG_HOST_CASE(kPartialSum, ipg_k_partial_sum)
            IPG_HOST_CASE(kGfUpdate, ipg_k_gf_update)
            IPG_HOST_CASE(kGfTransform, ipg_k_gf_transform)
            IPG_HOST_CASE(kWilsonCharge, ipg_k_wilson_charge)
            IPG_HOST_CASE(kWilsonExponent, ipg_k_wilson_exponent)
#undef IPG_HOST_CASE
            default:
                throw std::runtime_error("unknown device kernel");
        }
    }

  private:
    std::vector<ipg_real> &storage(ipgdev::Buffer buffer) {
        return buffers_.at(static_cast<std::size_t>(buffer));
    }
    const std::vector<ipg_real> &storage(ipgdev::Buffer buffer) const {
        return buffers_.at(static_cast<std::size_t>(buffer));
    }
    static void checkSize(const std::vector<ipg_real> &data, std::size_t n) {
        if (n > data.size()) {
            throw std::runtime_error("device buffer too small for the copy");
        }
    }

    std::vector<std::vector<ipg_real>> buffers_;
    std::vector<bool> used_;
};

}  // namespace

namespace ipgdev {
std::unique_ptr<Backend> IPG_HOST_FACTORY() {
    return std::unique_ptr<Backend>(new HostBackend());
}
}  // namespace ipgdev

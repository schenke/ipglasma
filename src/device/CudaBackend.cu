// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

// The CUDA backend of the device layer: the kernels of ipg_kernels.h as
// CUDA kernels, and cuFFT for the Fourier transforms. Double precision by
// default; IPGLASMA_DEVICE_PRECISION=single for single precision.
//
// Environment:
//   IPGLASMA_CUDA_DEVICE  number of the GPU (default: the rank of the
//                         process on its node, ipgdev::nodeLocalRank(),
//                         modulo the number of GPUs; 0 without MPI)
//   IPGLASMA_CUDA_FFT     `cufft` (default) or `kernel`, the mixed-radix
//                         FFT of ipg_kernels.h (slower; to cross-check
//                         cuFFT)

#include <cuda_runtime.h>
#include <cufft.h>

#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "DeviceBackend.h"

namespace {

#define IPG_FN __device__ inline
#define IPG_GLOBAL

namespace cuda_double {
#define IPG_REAL double
#include "CudaPrecisionImpl.h"
#undef IPG_REAL
}  // namespace cuda_double

namespace cuda_single {
#define IPG_REAL float
#include "CudaPrecisionImpl.h"
#undef IPG_REAL
}  // namespace cuda_single

void check(cudaError_t status, const char *what) {
    if (status != cudaSuccess) {
        std::ostringstream message;
        message << "CUDA: " << what << " failed: "
                << cudaGetErrorString(status);
        throw std::runtime_error(message.str());
    }
}

/// check() for the calls that wait for the device (copies,
/// synchronization): kernels run asynchronously, so an error of an
/// earlier kernel only shows up here.
void checkAfterWait(cudaError_t status, const char *what) {
    if (status != cudaSuccess) {
        std::ostringstream message;
        message << "CUDA: " << what << " failed: " << cudaGetErrorString(status)
                << " (this can be the error of an earlier kernel; "
                   "IPGLASMA_DEVICE_PROFILE=1 waits after every kernel and "
                   "names the one that fails)";
        throw std::runtime_error(message.str());
    }
}

void checkFft(cufftResult status, const char *what) {
    if (status != CUFFT_SUCCESS) {
        std::ostringstream message;
        message << "cuFFT: " << what << " failed with error "
                << static_cast<int>(status);
        throw std::runtime_error(message.str());
    }
}

class CudaBackend : public ipgdev::Backend {
  public:
    explicit CudaBackend(int precision) : single_(precision == 1) {
        int count = 0;
        check(cudaGetDeviceCount(&count), "cudaGetDeviceCount");
        // by default every MPI rank on a node gets its own GPU (with as
        // many ranks as GPUs); with one GPU per rank, e.g. from Slurm, the
        // rank sees only that one
        const int localRank = ipgdev::nodeLocalRank();
        int device = (localRank >= 0 && count > 0) ? localRank % count : 0;
        if (const char *request = std::getenv("IPGLASMA_CUDA_DEVICE")) {
            device = std::atoi(request);
        }
        if (device < 0 || device >= count) {
            std::ostringstream message;
            message << "CUDA: no GPU number " << device << " (" << count
                    << " found)";
            throw std::runtime_error(message.str());
        }
        check(cudaSetDevice(device), "cudaSetDevice");
        for (Staging &staging : staging_) {
            check(
                cudaEventCreateWithFlags(&staging.done, cudaEventDisableTiming),
                "cudaEventCreateWithFlags");
        }
        cudaDeviceProp properties;
        check(
            cudaGetDeviceProperties(&properties, device),
            "cudaGetDeviceProperties");
        deviceName_ = properties.name;
        if (count > 1) {
            deviceName_ += " (GPU " + std::to_string(device) + " of "
                           + std::to_string(count) + ")";
        }
        if (const char *fft = std::getenv("IPGLASMA_CUDA_FFT")) {
            const std::string choice(fft);
            if (choice == "kernel") {
                useCufft_ = false;
            } else if (choice != "cufft" && !choice.empty()) {
                throw std::runtime_error(
                    "IPGLASMA_CUDA_FFT must be cufft or kernel, not '" + choice
                    + "'");
            }
        }
    }

    ~CudaBackend() override {
        releaseShared();
        // errors are ignored here (e.g. a runtime that is shutting down)
        cudaDeviceSynchronize();
        for (Staging &staging : staging_) {
            if (staging.data != nullptr) cudaFreeHost(staging.data);
            if (staging.done != nullptr) cudaEventDestroy(staging.done);
        }
        for (const auto &plan : plans_) cufftDestroy(plan.second);
        for (void *buffer : buffers_) {
            if (buffer != nullptr) cudaFree(buffer);
        }
    }

    std::string name() const override { return "cuda"; }
    std::string description() const override {
        return deviceName_ + ", " + (single_ ? "single" : "double")
               + " precision, " + (useCufft_ ? "cuFFT" : "kernel FFT");
    }
    bool singlePrecision() const override { return single_; }

    ipgdev::Buffer createBuffer(std::size_t nReals) override {
        const std::size_t bytes = (nReals > 0 ? nReals : 1) * realSize();
        void *buffer = nullptr;
        check(cudaMalloc(&buffer, bytes), "cudaMalloc");
        check(cudaMemset(buffer, 0, bytes), "cudaMemset");
        for (std::size_t slot = 0; slot < buffers_.size(); slot++) {
            if (buffers_[slot] == nullptr) {
                buffers_[slot] = buffer;
                sizes_[slot] = nReals;
                return static_cast<ipgdev::Buffer>(slot);
            }
        }
        buffers_.push_back(buffer);
        sizes_.push_back(nReals);
        return static_cast<ipgdev::Buffer>(buffers_.size() - 1);
    }
    // Never throws, since the destructors of the drivers and of the backend
    // call it: a buffer released before is left alone, the error of a CUDA
    // runtime that is already shutting down (at program exit) is ignored,
    // and other errors are only reported.
    void releaseBuffer(ipgdev::Buffer buffer) override {
        const std::size_t index = static_cast<std::size_t>(buffer);
        if (buffer < 0 || index >= buffers_.size()) return;
        void *pointer = buffers_[index];
        buffers_[index] = nullptr;
        sizes_[index] = 0;
        if (pointer == nullptr) return;
        const cudaError_t status = cudaFree(pointer);
        if (status != cudaSuccess && status != cudaErrorCudartUnloading) {
            std::cerr << "[Warning] CUDA: cudaFree failed: "
                      << cudaGetErrorString(status) << std::endl;
        }
    }
    std::size_t bufferSize(ipgdev::Buffer buffer) const override {
        return sizes_.at(static_cast<std::size_t>(buffer));
    }
    // The copy to the device is asynchronous: the numbers are copied (in
    // single precision converted) with OpenMP threads into one of two
    // pinned staging buffers, and cudaMemcpyAsync sends them from there,
    // so that the host returns at once, e.g. to draw the noise of the next
    // JIMWLK step while the device copies and computes this one. A staging
    // buffer is filled again only after its previous copy has finished
    // (its event); the caller may reuse source at once, as before.
    void upload(
        ipgdev::Buffer buffer, const double *source,
        std::size_t nReals) override {
        checkSize(buffer, nReals);
        if (nReals == 0) return;
        const std::size_t bytes = nReals * realSize();
        Staging &staging = staging_[nextStaging_];
        nextStaging_ = 1 - nextStaging_;
        check(cudaEventSynchronize(staging.done), "cudaEventSynchronize");
        if (staging.bytes < bytes) {
            if (staging.data != nullptr) cudaFreeHost(staging.data);
            staging.data = nullptr;
            staging.bytes = 0;
            check(cudaMallocHost(&staging.data, bytes), "cudaMallocHost");
            staging.bytes = bytes;
        }
        if (single_) {
            ipgdev::parallelCopy(
                source, static_cast<float *>(staging.data), nReals);
        } else {
            ipgdev::parallelCopy(
                source, static_cast<double *>(staging.data), nReals);
        }
        check(
            cudaMemcpyAsync(
                memory(buffer), staging.data, bytes, cudaMemcpyHostToDevice, 0),
            "cudaMemcpyAsync to the device");
        check(cudaEventRecord(staging.done, 0), "cudaEventRecord");
    }
    void download(
        ipgdev::Buffer buffer, double *destination,
        std::size_t nReals) override {
        checkSize(buffer, nReals);
        void *data = destination;
        if (single_) {
            convert_.resize(nReals);
            data = convert_.data();
        }
        checkAfterWait(
            cudaMemcpy(
                data, memory(buffer), nReals * realSize(),
                cudaMemcpyDeviceToHost),
            "cudaMemcpy from the device");
        if (single_) {
            ipgdev::parallelCopy(convert_.data(), destination, nReals);
        }
    }
    void swapBuffers(ipgdev::Buffer a, ipgdev::Buffer b) override {
        std::swap(
            buffers_.at(static_cast<std::size_t>(a)),
            buffers_.at(static_cast<std::size_t>(b)));
    }

    void launch(
        ipgdev::Kernel kernel, std::size_t count, const ipgdev::Params &params,
        std::initializer_list<ipgdev::Buffer> buffers) override {
        if (count == 0) return;
        void *arguments[6] = {nullptr, nullptr, nullptr,
                              nullptr, nullptr, nullptr};
        int slot = 0;
        for (const ipgdev::Buffer buffer : buffers) {
            if (slot >= 6) throw std::runtime_error("too many kernel buffers");
            arguments[slot++] = memory(buffer);
        }
        if (single_) {
            cuda_single::launchKernel(
                kernel, static_cast<int>(count), params, arguments);
        } else {
            cuda_double::launchKernel(
                kernel, static_cast<int>(count), params, arguments);
        }
        check(cudaGetLastError(), "a kernel launch");
    }

    void finish() override {
        checkAfterWait(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
    }

    void releaseCaches() override {
        releaseShared();
        for (const auto &plan : plans_) cufftDestroy(plan.second);
        plans_.clear();
    }

    bool fftAvailable(int n) const override {
        // cuFFT transforms any size; the centering needs an even one
        return useCufft_ ? (n % 2 == 0) : Backend::fftAvailable(n);
    }

    void centeredFft(
        ipgdev::Buffer field, int n, int ncomp, int isign) override {
        if (!useCufft_) {
            Backend::centeredFft(field, n, ncomp, isign);
            return;
        }
        if (n % 2 != 0) {
            throw std::runtime_error(
                "the centered device FFT needs an even lattice size");
        }
        const long long elements = static_cast<long long>(n) * n * ncomp;
        // S F S = G M F M (FFT.cpp): checkerboard signs M before and after
        // the transform, with the global sign G and the normalization
        const double globalSign = (((n / 2 + n / 2) & 1) != 0) ? -1. : 1.;
        const double scaleAfter =
            globalSign
            * ((isign == -1) ? 1. / (static_cast<double>(n) * n) : 1.);
        // the fields are planar (ipg_kernels.h), the layout of cuFFT's
        // fast path: a batch of ncomp planes, stride 1, transformed in
        // place, with the signs before and the signs and scale after
        const cufftHandle plan = planFor(n, ncomp);
        const int direction = (isign == 1) ? CUFFT_FORWARD : CUFFT_INVERSE;
        ipgdev::Params p;
        p.n = n;
        p.ncomp = ncomp;
        p.i0 = 1;
        p.r0 = 1.;
        launch(
            ipgdev::kCheckerScale, static_cast<std::size_t>(elements), p,
            {field});
        if (single_) {
            cufftComplex *data = static_cast<cufftComplex *>(memory(field));
            checkFft(cufftExecC2C(plan, data, data, direction), "cufftExecC2C");
        } else {
            cufftDoubleComplex *data =
                static_cast<cufftDoubleComplex *>(memory(field));
            checkFft(cufftExecZ2Z(plan, data, data, direction), "cufftExecZ2Z");
        }
        p.r0 = scaleAfter;
        launch(
            ipgdev::kCheckerScale, static_cast<std::size_t>(elements), p,
            {field});
    }

  private:
    // The plan for the ncomp planes of a field (ipg_kernels.h):
    // a batch of ncomp transforms with stride 1 and distance n * n.
    cufftHandle planFor(int n, int ncomp) {
        const std::pair<int, int> key(n, ncomp);
        const auto found = plans_.find(key);
        if (found != plans_.end()) return found->second;
        cufftHandle plan;
        int dims[2] = {n, n};
        int embed[2] = {n, n};
        checkFft(
            cufftPlanMany(
                &plan, 2, dims, embed, 1, n * n, embed, 1, n * n,
                single_ ? CUFFT_C2C : CUFFT_Z2Z, ncomp),
            "cufftPlanMany");
        plans_[key] = plan;
        return plan;
    }
    std::size_t realSize() const {
        return single_ ? sizeof(float) : sizeof(double);
    }
    void *memory(ipgdev::Buffer buffer) const {
        void *pointer = buffers_.at(static_cast<std::size_t>(buffer));
        if (pointer == nullptr) {
            throw std::runtime_error("CUDA: use of a released buffer");
        }
        return pointer;
    }
    void checkSize(ipgdev::Buffer buffer, std::size_t nReals) const {
        if (nReals > bufferSize(buffer)) {
            throw std::runtime_error("device buffer too small for the copy");
        }
    }

    bool single_;
    bool useCufft_ = true;
    std::string deviceName_;
    std::vector<void *> buffers_;
    std::vector<std::size_t> sizes_;
    std::vector<float> convert_;
    /// A pinned host buffer of the asynchronous uploads, with the event of
    /// its last copy.
    struct Staging {
        void *data = nullptr;
        std::size_t bytes = 0;
        cudaEvent_t done = nullptr;
    };
    Staging staging_[2];
    int nextStaging_ = 0;
    std::map<std::pair<int, int>, cufftHandle> plans_;
};

}  // namespace

namespace ipgdev {
std::unique_ptr<Backend> createCudaBackend(int precision) {
    return std::unique_ptr<Backend>(new CudaBackend(precision));
}
}  // namespace ipgdev

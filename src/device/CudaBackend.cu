// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

// The CUDA backend of the device layer: the kernels of ipg_kernels.h as
// CUDA kernels, and cuFFT for the Fourier transforms. Double precision by
// default; IPGLASMA_DEVICE_PRECISION=single for single precision.
//
// Environment:
//   IPGLASMA_CUDA_DEVICE  number of the GPU (default 0)
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

// The layout of cuFFT's fast path. The device layer keeps the ncomp
// complex components of a site next to each other (component c of site s
// at c + ncomp * s), which cuFFT can only transform with stride ncomp, its
// slow path (about 3 times slower in double and 6-9 times in single
// precision at 800^2 on an RTX 4070). centeredFft() therefore copies a
// field into ncomp contiguous planes (component c of site s at
// c * sites + s), transforms them with stride 1, and copies them back. The
// checkerboard signs and the scale of the centered FFT (FFT.cpp) are
// applied in these copies, so they need no passes of their own.

/// field (interleaved) -> planes, times scale * (-1)^(x + y); i over the
/// complex elements of the field.
template <class Complex, class Real>
__global__ void interleavedToPlanes(
    const Complex *__restrict__ field, Complex *__restrict__ planes, int n,
    int ncomp, long long elements, Real scale) {
    const long long i =
        static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= elements) return;
    const long long sites = static_cast<long long>(n) * n;
    const long long site = i / ncomp;
    const int c = static_cast<int>(i % ncomp);
    const int x = static_cast<int>(site / n);
    const int y = static_cast<int>(site % n);
    const Real s = (((x + y) & 1) != 0) ? -scale : scale;
    Complex v = field[i];
    v.x *= s;
    v.y *= s;
    planes[c * sites + site] = v;
}

/// planes -> field (interleaved), times scale * (-1)^(x + y).
template <class Complex, class Real>
__global__ void planesToInterleaved(
    const Complex *__restrict__ planes, Complex *__restrict__ field, int n,
    int ncomp, long long elements, Real scale) {
    const long long i =
        static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= elements) return;
    const long long sites = static_cast<long long>(n) * n;
    const long long site = i / ncomp;
    const int c = static_cast<int>(i % ncomp);
    const int x = static_cast<int>(site / n);
    const int y = static_cast<int>(site % n);
    const Real s = (((x + y) & 1) != 0) ? -scale : scale;
    Complex v = planes[c * sites + site];
    v.x *= s;
    v.y *= s;
    field[i] = v;
}

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
        int device = 0;
        if (const char *request = std::getenv("IPGLASMA_CUDA_DEVICE")) {
            device = std::atoi(request);
        }
        int count = 0;
        check(cudaGetDeviceCount(&count), "cudaGetDeviceCount");
        if (device < 0 || device >= count) {
            std::ostringstream message;
            message << "CUDA: no GPU number " << device << " (" << count
                    << " found)";
            throw std::runtime_error(message.str());
        }
        check(cudaSetDevice(device), "cudaSetDevice");
        cudaDeviceProp properties;
        check(
            cudaGetDeviceProperties(&properties, device),
            "cudaGetDeviceProperties");
        deviceName_ = properties.name;
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
    void upload(
        ipgdev::Buffer buffer, const double *source,
        std::size_t nReals) override {
        checkSize(buffer, nReals);
        const void *data = source;
        if (single_) {
            convert_.resize(nReals);
            for (std::size_t k = 0; k < nReals; k++) {
                convert_[k] = static_cast<float>(source[k]);
            }
            data = convert_.data();
        }
        checkAfterWait(
            cudaMemcpy(
                memory(buffer), data, nReals * realSize(),
                cudaMemcpyHostToDevice),
            "cudaMemcpy to the device");
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
            for (std::size_t k = 0; k < nReals; k++) {
                destination[k] = static_cast<double>(convert_[k]);
            }
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
        // the transform, with the global sign G and the normalization, here
        // in the copies into and out of the planes of cuFFT's fast path
        const double globalSign = (((n / 2 + n / 2) & 1) != 0) ? -1. : 1.;
        const double scaleAfter =
            globalSign
            * ((isign == -1) ? 1. / (static_cast<double>(n) * n) : 1.);
        const ipgdev::Buffer planes = scratchLike(field);
        const cufftHandle plan = planFor(n, ncomp);
        const int direction = (isign == 1) ? CUFFT_FORWARD : CUFFT_INVERSE;
        const int threads = 256;
        const unsigned int blocks =
            static_cast<unsigned int>((elements + threads - 1) / threads);
        if (single_) {
            cufftComplex *data = static_cast<cufftComplex *>(memory(field));
            cufftComplex *planar = static_cast<cufftComplex *>(memory(planes));
            interleavedToPlanes<<<blocks, threads>>>(
                data, planar, n, ncomp, elements, 1.f);
            check(cudaGetLastError(), "a kernel launch");
            checkFft(
                cufftExecC2C(plan, planar, planar, direction), "cufftExecC2C");
            planesToInterleaved<<<blocks, threads>>>(
                planar, data, n, ncomp, elements,
                static_cast<float>(scaleAfter));
        } else {
            cufftDoubleComplex *data =
                static_cast<cufftDoubleComplex *>(memory(field));
            cufftDoubleComplex *planar =
                static_cast<cufftDoubleComplex *>(memory(planes));
            interleavedToPlanes<<<blocks, threads>>>(
                data, planar, n, ncomp, elements, 1.);
            check(cudaGetLastError(), "a kernel launch");
            checkFft(
                cufftExecZ2Z(plan, planar, planar, direction), "cufftExecZ2Z");
            planesToInterleaved<<<blocks, threads>>>(
                planar, data, n, ncomp, elements, scaleAfter);
        }
        check(cudaGetLastError(), "a kernel launch");
    }

  private:
    // The plan for the ncomp planes of a field (see interleavedToPlanes):
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
    std::map<std::pair<int, int>, cufftHandle> plans_;
};

}  // namespace

namespace ipgdev {
std::unique_ptr<Backend> createCudaBackend(int precision) {
    return std::unique_ptr<Backend>(new CudaBackend(precision));
}
}  // namespace ipgdev

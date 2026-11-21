// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

// The kernels of ipg_kernels.h as CUDA kernels for one precision. Included
// twice by CudaBackend.cu, in the namespaces cuda_double and cuda_single,
// with IPG_REAL defined; not a header for other use.

#undef SRC_DEVICE_IPG_KERNELS_H_
#include "ipg_kernels.h"

// One __global__ function per kernel: thread i of the launch runs index i.
#define IPG_CUDA_KERNEL(id, name)                                          \
    __global__ void name##_cuda(                                           \
        int count, ipg_params p, ipg_real *b0, ipg_real *b1, ipg_real *b2, \
        ipg_real *b3, ipg_real *b4, ipg_real *b5) {                        \
        const int i = (int)(blockIdx.x * blockDim.x + threadIdx.x);        \
        if (i < count) name(i, p, b0, b1, b2, b3, b4, b5);                 \
    }
IPG_KERNEL_LIST(IPG_CUDA_KERNEL)
#undef IPG_CUDA_KERNEL

inline void launchKernel(
    ipgdev::Kernel kernel, int count, const ipgdev::Params &params,
    void *const buffers[6]) {
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
    ipg_real *b[6];
    for (int k = 0; k < 6; k++) b[k] = static_cast<ipg_real *>(buffers[k]);
    const unsigned int block = 256;
    const unsigned int grid =
        (static_cast<unsigned int>(count) + block - 1) / block;
    switch (kernel) {
#define IPG_CUDA_CASE(id, function)                        \
    case ipgdev::id:                                       \
        function##_cuda<<<grid, block>>>(                  \
            count, p, b[0], b[1], b[2], b[3], b[4], b[5]); \
        break;
        IPG_KERNEL_LIST(IPG_CUDA_CASE)
#undef IPG_CUDA_CASE
        default:
            throw std::runtime_error("unknown device kernel");
    }
}

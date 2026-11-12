/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
 * This file is part of IP-Glasma; the license text is in LICENSE.
 *
 * Device kernels of IP-Glasma: the per-site (or per-element) bodies of the
 * JIMWLK step, of the classical Yang-Mills leapfrog and of a mixed-radix
 * FFT, in the common subset of C99, C++, CUDA C++ and OpenCL C 1.2, so
 * that the same source runs
 *
 *   - on the host (HostBackend.cpp, the reference of the device layer),
 *   - on NVIDIA GPUs with CUDA (CudaBackend.cu),
 *   - on OpenCL devices, among them the GPUs of Apple M processors
 *     (OpenCLBackend.cpp).
 *
 * Rules of this file (they are what makes it portable; OpenCL C 1.2 is the
 * narrowest of the four languages):
 *   - no C++ features, no recursion, no function pointers, no variable
 *     length arrays, no static data;
 *   - helper functions take and return structs by value; only a kernel and
 *     the load/store helpers touch a buffer, through IPG_GLOBAL pointers;
 *   - every floating-point constant goes through IPG_C(), every math
 *     function through an IPG_ macro, so that single precision does not
 *     use a double anywhere (Apple GPUs have no double precision);
 *   - every kernel has the same signature, (index, parameters, six
 *     buffers), so that a backend needs one launch path for all of them.
 *
 * The includer defines, before including this file:
 *   IPG_REAL      float or double
 *   IPG_FN        the qualifiers of a device function
 *   IPG_GLOBAL    the address space of a buffer pointer (OpenCL: __global)
 *   IPG_SQRT ...  the math functions for IPG_REAL (defaults below)
 *
 * Layout of the data: a complex number is two consecutive reals (re, im);
 * an SU(3) matrix is 9 complex numbers in row-major order, as in Matrix; a
 * field has its components of one site next to each other, the sites in
 * the order pos = ix * N + iy of LatticeIndex.h.
 */

#ifndef SRC_DEVICE_IPG_KERNELS_H_
#define SRC_DEVICE_IPG_KERNELS_H_

#ifndef IPG_REAL
#define IPG_REAL double
#endif
#ifndef IPG_FN
#define IPG_FN static inline
#endif
#ifndef IPG_GLOBAL
#define IPG_GLOBAL
#endif
#ifndef IPG_SQRT
#define IPG_SQRT(x) sqrt(x)
#define IPG_SIN(x) sin(x)
#define IPG_COS(x) cos(x)
#define IPG_ACOS(x) acos(x)
#define IPG_LOG(x) log(x)
#define IPG_FABS(x) fabs(x)
#endif

typedef IPG_REAL ipg_real;
#define IPG_C(x) ((ipg_real)(x))
/* 32 bits in all four languages */
typedef unsigned int ipg_u32;

/* The parameters of a kernel launch. Only int and ipg_real members, the
 * ints first and an even number of them, so that the host and the device
 * agree on the layout without padding. */
typedef struct {
    int n;     /* lattice side N */
    int ncomp; /* complex components per site of the field of the launch */
    int i0;
    int i1;
    int i2;
    int i3;
    int i4;
    int i5;
    ipg_real r0;
    ipg_real r1;
    ipg_real r2;
    ipg_real r3;
} ipg_params;

typedef struct {
    ipg_real e[18];
} ipg_su3;

/* ------------------------------------------------------------------ */
/* SU(3) algebra on values                                             */
/* ------------------------------------------------------------------ */

IPG_FN ipg_su3 ipg_load(IPG_GLOBAL const ipg_real *b, int idx) {
    ipg_su3 m;
    int k;
    for (k = 0; k < 18; k++) m.e[k] = b[18 * idx + k];
    return m;
}

IPG_FN void ipg_store(IPG_GLOBAL ipg_real *b, int idx, ipg_su3 m) {
    int k;
    for (k = 0; k < 18; k++) b[18 * idx + k] = m.e[k];
}

IPG_FN ipg_su3 ipg_zero(void) {
    ipg_su3 m;
    int k;
    for (k = 0; k < 18; k++) m.e[k] = IPG_C(0.);
    return m;
}

IPG_FN ipg_su3 ipg_identity(void) {
    ipg_su3 m = ipg_zero();
    m.e[0] = IPG_C(1.);
    m.e[8] = IPG_C(1.);
    m.e[16] = IPG_C(1.);
    return m;
}

/* a b */
IPG_FN ipg_su3 ipg_mul(ipg_su3 a, ipg_su3 b) {
    ipg_su3 c;
    int i, j, k;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            ipg_real re = IPG_C(0.);
            ipg_real im = IPG_C(0.);
            for (k = 0; k < 3; k++) {
                const ipg_real ar = a.e[2 * (3 * i + k)];
                const ipg_real ai = a.e[2 * (3 * i + k) + 1];
                const ipg_real br = b.e[2 * (3 * k + j)];
                const ipg_real bi = b.e[2 * (3 * k + j) + 1];
                re += ar * br - ai * bi;
                im += ar * bi + ai * br;
            }
            c.e[2 * (3 * i + j)] = re;
            c.e[2 * (3 * i + j) + 1] = im;
        }
    }
    return c;
}

/* a b^dagger (Matrix::prodABconj) */
IPG_FN ipg_su3 ipg_mul_adj(ipg_su3 a, ipg_su3 b) {
    ipg_su3 c;
    int i, j, k;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            ipg_real re = IPG_C(0.);
            ipg_real im = IPG_C(0.);
            for (k = 0; k < 3; k++) {
                const ipg_real ar = a.e[2 * (3 * i + k)];
                const ipg_real ai = a.e[2 * (3 * i + k) + 1];
                const ipg_real br = b.e[2 * (3 * j + k)];
                const ipg_real bi = -b.e[2 * (3 * j + k) + 1];
                re += ar * br - ai * bi;
                im += ar * bi + ai * br;
            }
            c.e[2 * (3 * i + j)] = re;
            c.e[2 * (3 * i + j) + 1] = im;
        }
    }
    return c;
}

/* a^dagger b (Matrix::prodAconjB) */
IPG_FN ipg_su3 ipg_adj_mul(ipg_su3 a, ipg_su3 b) {
    ipg_su3 c;
    int i, j, k;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            ipg_real re = IPG_C(0.);
            ipg_real im = IPG_C(0.);
            for (k = 0; k < 3; k++) {
                const ipg_real ar = a.e[2 * (3 * k + i)];
                const ipg_real ai = -a.e[2 * (3 * k + i) + 1];
                const ipg_real br = b.e[2 * (3 * k + j)];
                const ipg_real bi = b.e[2 * (3 * k + j) + 1];
                re += ar * br - ai * bi;
                im += ar * bi + ai * br;
            }
            c.e[2 * (3 * i + j)] = re;
            c.e[2 * (3 * i + j) + 1] = im;
        }
    }
    return c;
}

IPG_FN ipg_su3 ipg_add(ipg_su3 a, ipg_su3 b) {
    ipg_su3 c;
    int k;
    for (k = 0; k < 18; k++) c.e[k] = a.e[k] + b.e[k];
    return c;
}

IPG_FN ipg_su3 ipg_sub(ipg_su3 a, ipg_su3 b) {
    ipg_su3 c;
    int k;
    for (k = 0; k < 18; k++) c.e[k] = a.e[k] - b.e[k];
    return c;
}

IPG_FN ipg_su3 ipg_scale(ipg_su3 a, ipg_real s) {
    ipg_su3 c;
    int k;
    for (k = 0; k < 18; k++) c.e[k] = s * a.e[k];
    return c;
}

/* a^dagger */
IPG_FN ipg_su3 ipg_adj(ipg_su3 a) {
    ipg_su3 c;
    int i, j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            c.e[2 * (3 * i + j)] = a.e[2 * (3 * j + i)];
            c.e[2 * (3 * i + j) + 1] = -a.e[2 * (3 * j + i) + 1];
        }
    }
    return c;
}

/* a + s * identity */
IPG_FN ipg_su3 ipg_add_diag(ipg_su3 a, ipg_real s) {
    a.e[0] += s;
    a.e[8] += s;
    a.e[16] += s;
    return a;
}

/* the inverse from the cofactors, as cofactorInverse() of Matrix.cpp */
IPG_FN ipg_su3 ipg_inv(ipg_su3 q) {
    /* cofactor (i, j) = q[r1][c1] q[r2][c2] - q[r1][c2] q[r2][c1], with the
     * index pairs of cofactorInverse() */
    ipg_su3 h;
    ipg_real detRe = IPG_C(0.);
    ipg_real detIm = IPG_C(0.);
    ipg_real inv;
    int i, j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            /* adjugate(i, j) = cofactor(j, i): rows of q without j, columns
             * without i, in cyclic order, which carries the sign */
            const int r1 = (j + 1) % 3;
            const int r2 = (j + 2) % 3;
            const int c1 = (i + 1) % 3;
            const int c2 = (i + 2) % 3;
            const ipg_real ar = q.e[2 * (3 * r1 + c1)];
            const ipg_real ai = q.e[2 * (3 * r1 + c1) + 1];
            const ipg_real br = q.e[2 * (3 * r2 + c2)];
            const ipg_real bi = q.e[2 * (3 * r2 + c2) + 1];
            const ipg_real cr = q.e[2 * (3 * r1 + c2)];
            const ipg_real ci = q.e[2 * (3 * r1 + c2) + 1];
            const ipg_real dr = q.e[2 * (3 * r2 + c1)];
            const ipg_real di = q.e[2 * (3 * r2 + c1) + 1];
            h.e[2 * (3 * i + j)] = (ar * br - ai * bi) - (cr * dr - ci * di);
            h.e[2 * (3 * i + j) + 1] =
                (ar * bi + ai * br) - (cr * di + ci * dr);
        }
    }
    /* det = sum_j q[0][j] * adjugate[j][0] */
    for (j = 0; j < 3; j++) {
        const ipg_real qr = q.e[2 * j];
        const ipg_real qi = q.e[2 * j + 1];
        const ipg_real hr = h.e[2 * (3 * j)];
        const ipg_real hi = h.e[2 * (3 * j) + 1];
        detRe += qr * hr - qi * hi;
        detIm += qr * hi + qi * hr;
    }
    inv = IPG_C(1.) / (detRe * detRe + detIm * detIm);
    for (i = 0; i < 9; i++) {
        const ipg_real hr = h.e[2 * i];
        const ipg_real hi = h.e[2 * i + 1];
        h.e[2 * i] = (hr * detRe + hi * detIm) * inv;
        h.e[2 * i + 1] = (hi * detRe - hr * detIm) * inv;
    }
    return h;
}

/* exp(a) by a (6,6) Pade approximant with scaling and squaring: the
 * algorithm of Matrix::expm(), step by step */
IPG_FN ipg_su3 ipg_expm(ipg_su3 a) {
    ipg_real c[7];
    ipg_real norm = IPG_C(0.);
    ipg_su3 u, h2, pm, qm;
    int i, j, k, s, odd;

    c[0] = IPG_C(1.);
    for (k = 1; k <= 6; k++) {
        c[k] = c[k - 1]
               * ((IPG_C(7.) - (ipg_real)k)
                  / ((ipg_real)k * (IPG_C(13.) - (ipg_real)k)));
    }
    for (i = 0; i < 3; i++) {
        ipg_real row = IPG_C(0.);
        for (j = 0; j < 3; j++) {
            const ipg_real re = a.e[2 * (3 * i + j)];
            const ipg_real im = a.e[2 * (3 * i + j) + 1];
            row += IPG_SQRT(re * re + im * im);
        }
        if (row > norm) norm = row;
    }
    if (!(norm > IPG_C(0.))) return ipg_identity();

    s = 0;
    u = a;
    if (norm > IPG_C(0.5)) {
        ipg_real scale = IPG_C(1.);
        s = (int)(IPG_LOG(norm) / IPG_LOG(IPG_C(2.)) + IPG_C(2.));
        if (s < 0) s = 0;
        for (i = 0; i < s; i++) scale *= IPG_C(0.5);
        u = ipg_scale(a, scale);
    }

    h2 = ipg_mul(u, u);
    qm = ipg_scale(ipg_identity(), c[6]);
    pm = ipg_scale(ipg_identity(), c[5]);
    odd = 1;
    for (k = 5; k > 0; k--) {
        if (odd == 1) {
            qm = ipg_add_diag(ipg_mul(qm, h2), c[k - 1]);
        } else {
            pm = ipg_add_diag(ipg_mul(pm, h2), c[k - 1]);
        }
        odd = 1 - odd;
    }
    if (odd == 1) {
        qm = ipg_mul(qm, u);
    } else {
        pm = ipg_mul(pm, u);
    }
    qm = ipg_sub(qm, pm);
    h2 = ipg_inv(qm);
    u = ipg_add_diag(ipg_scale(ipg_mul(h2, pm), IPG_C(2.)), IPG_C(1.));
    if (odd == 1) u = ipg_scale(u, IPG_C(-1.));
    for (i = 0; i < s; i++) u = ipg_mul(u, u);
    return u;
}

/* exp(i Q) of a Hermitian traceless Q from the Cayley-Hamilton theorem: the
 * algorithm of Matrix::expiHermitian(), step by step */
IPG_FN ipg_su3 ipg_expi_hermitian(ipg_su3 q) {
    const ipg_su3 q2 = ipg_mul(q, q);
    ipg_su3 result;
    ipg_real c1 = IPG_C(0.);
    ipg_real c0, c1Over3, c0max, ratio, thetaOverThree, u, w, w2, xi0, cosw;
    ipg_real e2r, e2i, emr, emi, u2, den;
    ipg_real f0r, f0i, f1r, f1i, f2r, f2i, tr, ti;
    int k;

    for (k = 0; k < 18; k++) c1 += q.e[k] * q.e[k];
    c1 *= IPG_C(0.5);
    if (c1 < IPG_C(1e-24)) {
        /* 1 + i Q - Q^2 / 2 */
        for (k = 0; k < 9; k++) {
            result.e[2 * k] = -q.e[2 * k + 1] - IPG_C(0.5) * q2.e[2 * k];
            result.e[2 * k + 1] = q.e[2 * k] - IPG_C(0.5) * q2.e[2 * k + 1];
        }
        return ipg_add_diag(result, IPG_C(1.));
    }
    /* c0 = Re det Q */
    {
        ipg_real dr = IPG_C(0.);
        int t;
        for (t = 0; t < 3; t++) {
            /* + q[0][t] q[1][t+1] q[2][t+2] - q[0][t] q[1][t+2] q[2][t+1] */
            const int t1 = (t + 1) % 3;
            const int t2 = (t + 2) % 3;
            const ipg_real ar = q.e[2 * t];
            const ipg_real ai = q.e[2 * t + 1];
            const ipg_real br = q.e[2 * (3 + t1)];
            const ipg_real bi = q.e[2 * (3 + t1) + 1];
            const ipg_real cr = q.e[2 * (6 + t2)];
            const ipg_real ci = q.e[2 * (6 + t2) + 1];
            const ipg_real dr2 = q.e[2 * (3 + t2)];
            const ipg_real di2 = q.e[2 * (3 + t2) + 1];
            const ipg_real er = q.e[2 * (6 + t1)];
            const ipg_real ei = q.e[2 * (6 + t1) + 1];
            const ipg_real abr = ar * br - ai * bi;
            const ipg_real abi = ar * bi + ai * br;
            const ipg_real adr = ar * dr2 - ai * di2;
            const ipg_real adi = ar * di2 + ai * dr2;
            dr += (abr * cr - abi * ci) - (adr * er - adi * ei);
        }
        c0 = dr;
    }
    c1Over3 = c1 / IPG_C(3.);
    c0max = IPG_C(2.) * c1Over3 * IPG_SQRT(c1Over3);
    ratio = IPG_FABS(c0) / c0max;
    if (ratio > IPG_C(1.)) ratio = IPG_C(1.);
    thetaOverThree = IPG_ACOS(ratio) / IPG_C(3.);
    u = IPG_SQRT(c1Over3) * IPG_COS(thetaOverThree);
    w = IPG_SQRT(c1) * IPG_SIN(thetaOverThree);
    w2 = w * w;
    xi0 = (w <= IPG_C(0.05))
              ? IPG_C(1.)
                    - w2 / IPG_C(6.)
                          * (IPG_C(1.)
                             - w2 / IPG_C(20.) * (IPG_C(1.) - w2 / IPG_C(42.)))
              : IPG_SIN(w) / w;
    cosw = IPG_COS(w);
    e2r = IPG_COS(IPG_C(2.) * u);
    e2i = IPG_SIN(IPG_C(2.) * u);
    emr = IPG_COS(u);
    emi = -IPG_SIN(u);
    u2 = u * u;
    den = IPG_C(9.) * u2 - w2;
    /* f0 = ((u2 - w2) e2 + em (8 u2 cosw + 2 i u (3 u2 + w2) xi0)) / den */
    tr = IPG_C(8.) * u2 * cosw;
    ti = IPG_C(2.) * u * (IPG_C(3.) * u2 + w2) * xi0;
    f0r = ((u2 - w2) * e2r + (emr * tr - emi * ti)) / den;
    f0i = ((u2 - w2) * e2i + (emr * ti + emi * tr)) / den;
    /* f1 = (2 u e2 - em (2 u cosw - i (3 u2 - w2) xi0)) / den */
    tr = IPG_C(2.) * u * cosw;
    ti = -(IPG_C(3.) * u2 - w2) * xi0;
    f1r = (IPG_C(2.) * u * e2r - (emr * tr - emi * ti)) / den;
    f1i = (IPG_C(2.) * u * e2i - (emr * ti + emi * tr)) / den;
    /* f2 = (e2 - em (cosw + 3 i u xi0)) / den */
    tr = cosw;
    ti = IPG_C(3.) * u * xi0;
    f2r = (e2r - (emr * tr - emi * ti)) / den;
    f2i = (e2i - (emr * ti + emi * tr)) / den;
    if (c0 < IPG_C(0.)) {
        /* f0 = conj(f0), f1 = -conj(f1), f2 = conj(f2) */
        f0i = -f0i;
        f1r = -f1r;
        f2i = -f2i;
    }
    for (k = 0; k < 9; k++) {
        const ipg_real qr = q.e[2 * k];
        const ipg_real qi = q.e[2 * k + 1];
        const ipg_real pr = q2.e[2 * k];
        const ipg_real pi = q2.e[2 * k + 1];
        result.e[2 * k] = (f1r * qr - f1i * qi) + (f2r * pr - f2i * pi);
        result.e[2 * k + 1] = (f1r * qi + f1i * qr) + (f2r * pi + f2i * pr);
    }
    result.e[0] += f0r;
    result.e[1] += f0i;
    result.e[8] += f0r;
    result.e[9] += f0i;
    result.e[16] += f0r;
    result.e[17] += f0i;
    return result;
}

/* ------------------------------------------------------------------ */
/* Lattice neighbors                                                   */
/* ------------------------------------------------------------------ */

/* The coordinate i + d on a lattice of side n: clamped at the edges, or
 * periodic (Lattice::boundedCoordinate()). d is -1, 0 or 1. */
IPG_FN int ipg_coord(int i, int d, int n, int periodic) {
    int j = i + d;
    if (periodic != 0) {
        if (j < 0) j += n;
        if (j >= n) j -= n;
    } else {
        if (j < 0) j = 0;
        if (j >= n) j = n - 1;
    }
    return j;
}

/* ------------------------------------------------------------------ */
/* Random numbers on the device                                        */
/* ------------------------------------------------------------------ */

typedef struct {
    ipg_u32 a;
    ipg_u32 b;
} ipg_u32_pair;

/* Threefry-2x32 with 20 rounds (Salmon, Moraes, Dror and Shaw, "Parallel
 * random numbers: as easy as 1, 2, 3", SC'11): a counter-based generator,
 * so that every element of a launch computes its own random numbers from
 * (counter, key) without any shared state, and with 32-bit additions,
 * rotations and exclusive-ors only, which every device has. */
IPG_FN ipg_u32_pair
ipg_threefry2x32(ipg_u32 c0, ipg_u32 c1, ipg_u32 k0, ipg_u32 k1) {
    ipg_u32 ks[3];
    ipg_u32 x0, x1;
    ipg_u32_pair out;
    int r;
    ks[0] = k0;
    ks[1] = k1;
    ks[2] = 0x1BD11BDAu ^ k0 ^ k1;
    x0 = c0 + ks[0];
    x1 = c1 + ks[1];
    for (r = 0; r < 20; r++) {
        int rot;
        switch (r % 8) {
            case 0:
                rot = 13;
                break;
            case 1:
                rot = 15;
                break;
            case 2:
                rot = 26;
                break;
            case 3:
                rot = 6;
                break;
            case 4:
                rot = 17;
                break;
            case 5:
                rot = 29;
                break;
            case 6:
                rot = 16;
                break;
            default:
                rot = 24;
                break;
        }
        x0 += x1;
        x1 = (x1 << rot) | (x1 >> (32 - rot));
        x1 ^= x0;
        if (r % 4 == 3) {
            const int inject = r / 4 + 1;
            x0 += ks[inject % 3];
            x1 += ks[(inject + 1) % 3] + (ipg_u32)inject;
        }
    }
    out.a = x0;
    out.b = x1;
    return out;
}

/* ------------------------------------------------------------------ */
/* Kernels. Common signature:                                          */
/*   (i, p, b0, b1, b2, b3, b4, b5)                                    */
/* i runs over the sites, or over the elements named in the comment.   */
/* ------------------------------------------------------------------ */

#define IPG_KERNEL_ARGS                                                    \
    int i, ipg_params p, IPG_GLOBAL ipg_real *b0, IPG_GLOBAL ipg_real *b1, \
        IPG_GLOBAL ipg_real *b2, IPG_GLOBAL ipg_real *b3,                  \
        IPG_GLOBAL ipg_real *b4, IPG_GLOBAL ipg_real *b5

/* b1[i] = b0[i]; i over reals */
IPG_FN void ipg_k_copy(IPG_KERNEL_ARGS) { b1[i] = b0[i]; }

/* complex b1[i] = (real b0[i], 0); i over the complex elements */
IPG_FN void ipg_k_real_to_complex(IPG_KERNEL_ARGS) {
    b1[2 * i] = b0[i];
    b1[2 * i + 1] = IPG_C(0.);
}

/* real b1[i] = Re(complex b0[i]); i over the complex elements */
IPG_FN void ipg_k_complex_to_real(IPG_KERNEL_ARGS) { b1[i] = b0[2 * i]; }

/* complex b0[site][c] *= real b1[site]; i over the complex elements,
 * p.ncomp components per site */
IPG_FN void ipg_k_scale_by_site(IPG_KERNEL_ARGS) {
    const ipg_real s = b1[i / p.ncomp];
    b0[2 * i] *= s;
    b0[2 * i + 1] *= s;
}

/* complex b0[i] *= r0 * (-1)^(ix + iy) if i0, else r0; i over the complex
 * elements, p.ncomp components per site. The checkerboard signs turn an
 * FFT in natural order into the centered FFT of FFT.cpp (see there). */
IPG_FN void ipg_k_checker_scale(IPG_KERNEL_ARGS) {
    ipg_real s = p.r0;
    if (p.i0 != 0) {
        const int site = i / p.ncomp;
        if ((((site / p.n) + (site % p.n)) & 1) != 0) s = -s;
    }
    b0[2 * i] *= s;
    b0[2 * i + 1] *= s;
}

/* One pass of a mixed-radix Stockham FFT along one axis of a field with
 * p.ncomp complex components per site; i over the complex elements.
 *   b0 source, b1 destination (a different buffer), b2 the table
 *   (cos, sin)(2 pi t / N), t = 0 .. N-1
 *   i0 axis: 0 along iy, 1 along ix
 *   i1 radix r of the pass
 *   i2 s, the product of the radices of the earlier passes of the axis
 *   i3 +1 forward (exp(-i...), FFTW_FORWARD), -1 backward
 *   i4 bit 0: multiply the source by (-1)^(ix+iy); bit 1: the destination
 *   r0 factor of the destination (the normalization and global sign)
 * With n = N / s the length of the subsequences and m = n / r, the
 * destination t = q + s (r pp + k) (q < s, k < r, pp < m) is
 *   sum_j source[q + s (pp + j m)] exp(-/+ 2 pi i k (pp + j m) / n).
 * After the passes of all radices the axis is transformed, in natural
 * order. */
IPG_FN void ipg_k_fft_pass(IPG_KERNEL_ARGS) {
    const int n = p.n;
    const int comp = i % p.ncomp;
    const int site = i / p.ncomp;
    const int ix = site / n;
    const int iy = site % n;
    const int t = (p.i0 == 0) ? iy : ix;
    const int other = (p.i0 == 0) ? ix : iy;
    const int r = p.i1;
    const int s = p.i2;
    const int m = (n / s) / r;
    const int q = t % s;
    const int k = (t / s) % r;
    const int pp = (t / s) / r;
    ipg_real re = IPG_C(0.);
    ipg_real im = IPG_C(0.);
    ipg_real f = p.r0;
    int j;
    for (j = 0; j < r; j++) {
        const int l = pp + j * m;
        const int tj = q + s * l;
        const int srcSite = (p.i0 == 0) ? (ix * n + tj) : (tj * n + iy);
        const int src = srcSite * p.ncomp + comp;
        const int tw = (k * l * s) % n;
        const ipg_real wr = b2[2 * tw];
        const ipg_real wi = (p.i3 > 0) ? -b2[2 * tw + 1] : b2[2 * tw + 1];
        ipg_real xr = b0[2 * src];
        ipg_real xi = b0[2 * src + 1];
        if ((p.i4 & 1) != 0 && ((tj + other) & 1) != 0) {
            xr = -xr;
            xi = -xi;
        }
        re += xr * wr - xi * wi;
        im += xr * wi + xi * wr;
    }
    if ((p.i4 & 2) != 0 && ((ix + iy) & 1) != 0) f = -f;
    b1[2 * i] = f * re;
    b1[2 * i + 1] = f * im;
}

/* ---- JIMWLK (JIMWLK::evolutionStep) ---- */

/* Two independent Gaussian random numbers of unit variance, b0[2 i] and
 * b0[2 i + 1], from the Threefry counter (i, i2) with the key (i0, i1), by
 * the Box-Muller transformation; i over the pairs. The first uniform
 * number is (w + 1/2) / 2^32, in (0, 1], which keeps the resolution of the
 * 32-bit word near zero also in single precision, so that the Gaussian
 * tail reaches 6.7 sigma. */
IPG_FN void ipg_k_gauss_noise(IPG_KERNEL_ARGS) {
    const ipg_u32_pair w = ipg_threefry2x32(
        (ipg_u32)i, (ipg_u32)p.i2, (ipg_u32)p.i0, (ipg_u32)p.i1);
    const ipg_real u1 =
        ((ipg_real)w.a + IPG_C(0.5)) * IPG_C(2.3283064365386962890625e-10);
    const ipg_real u2 =
        ((ipg_real)w.b + IPG_C(0.5)) * IPG_C(2.3283064365386962890625e-10);
    const ipg_real radius = IPG_SQRT(IPG_C(-2.) * IPG_LOG(u1));
    const ipg_real angle = IPG_C(6.283185307179586476925286766559) * u2;
    b0[2 * i] = radius * IPG_COS(angle);
    b0[2 * i + 1] = radius * IPG_SIN(angle);
}

/* CK[site][a] = K_x xi_x^a + K_y xi_y^a in momentum space.
 *   b0 kernel K [site][2] complex, b1 xi [site][16] complex,
 *   b2 CK [site][8] complex */
IPG_FN void ipg_k_jimwlk_kernel_noise(IPG_KERNEL_ARGS) {
    const ipg_real kxr = b0[4 * i];
    const ipg_real kxi = b0[4 * i + 1];
    const ipg_real kyr = b0[4 * i + 2];
    const ipg_real kyi = b0[4 * i + 3];
    int a;
    for (a = 0; a < 8; a++) {
        const ipg_real xr = b1[32 * i + 2 * a];
        const ipg_real xi = b1[32 * i + 2 * a + 1];
        const ipg_real yr = b1[32 * i + 2 * (a + 8)];
        const ipg_real yi = b1[32 * i + 2 * (a + 8) + 1];
        b2[16 * i + 2 * a] = (kxr * xr - kxi * xi) + (kyr * yr - kyi * yi);
        b2[16 * i + 2 * a + 1] = (kxr * xi + kxi * xr) + (kyr * yi + kyi * yr);
    }
}

/* V_d[site] = U (sum_a xi_d^a t_a) U^dagger for d = x, y.
 *   b0 Wilson lines U, b1 real noise [site][16], b2 generators [8] su3,
 *   b3 V [site][2] su3 */
IPG_FN void ipg_k_jimwlk_rotate(IPG_KERNEL_ARGS) {
    const ipg_su3 U = ipg_load(b0, i);
    int d, a, k;
    for (d = 0; d < 2; d++) {
        ipg_su3 A = ipg_zero();
        for (a = 0; a < 8; a++) {
            const ipg_real c = b1[16 * i + 8 * d + a];
            for (k = 0; k < 18; k++) A.e[k] += c * b2[18 * a + k];
        }
        ipg_store(b3, 2 * i + d, ipg_mul(U, ipg_mul_adj(A, U)));
    }
}

/* W[site] = K_x V_x + K_y V_y in momentum space.
 *   b0 kernel K, b1 V [site][2] su3, b2 W [site] su3 */
IPG_FN void ipg_k_jimwlk_kernel_matrix(IPG_KERNEL_ARGS) {
    const ipg_real kxr = b0[4 * i];
    const ipg_real kxi = b0[4 * i + 1];
    const ipg_real kyr = b0[4 * i + 2];
    const ipg_real kyi = b0[4 * i + 3];
    int k;
    for (k = 0; k < 9; k++) {
        const ipg_real xr = b1[36 * i + 2 * k];
        const ipg_real xi = b1[36 * i + 2 * k + 1];
        const ipg_real yr = b1[36 * i + 18 + 2 * k];
        const ipg_real yi = b1[36 * i + 18 + 2 * k + 1];
        b2[18 * i + 2 * k] = (kxr * xr - kxi * xi) + (kyr * yr - kyi * yi);
        b2[18 * i + 2 * k + 1] = (kxr * xi + kxi * xr) + (kyr * yi + kyi * yr);
    }
}

/* U <- exp(-i sqrt(ds) W) U exp(i sqrt(ds) sum_a Re(CK^a) t_a).
 *   b0 Wilson lines U, b1 W [site] su3, b2 CK [site][8] complex,
 *   b3 generators; r0 = sqrt(ds) */
IPG_FN void ipg_k_jimwlk_update(IPG_KERNEL_ARGS) {
    const ipg_su3 W = ipg_load(b1, i);
    ipg_su3 left, right, A;
    int a, k;
    /* -i sqrt(ds) (re + i im) = sqrt(ds) (im - i re) */
    for (k = 0; k < 9; k++) {
        left.e[2 * k] = p.r0 * W.e[2 * k + 1];
        left.e[2 * k + 1] = -p.r0 * W.e[2 * k];
    }
    A = ipg_zero();
    for (a = 0; a < 8; a++) {
        const ipg_real c = b2[16 * i + 2 * a];
        for (k = 0; k < 18; k++) A.e[k] += c * b3[18 * a + k];
    }
    /* i sqrt(ds) (re + i im) = sqrt(ds) (-im + i re) */
    for (k = 0; k < 9; k++) {
        right.e[2 * k] = -p.r0 * A.e[2 * k + 1];
        right.e[2 * k + 1] = p.r0 * A.e[2 * k];
    }
    ipg_store(
        b0, i,
        ipg_mul(ipg_mul(ipg_expm(left), ipg_load(b0, i)), ipg_expm(right)));
}

/* ---- classical Yang-Mills leapfrog (Evolution.cpp) ----
 * Fields: E1 = Lattice::U, E2 = U2, links Ux, Uy, pi = Ux2, phi = Uy2.
 * i0 = 1 for a periodic lattice, 0 for clamped edges. */

/* pi += r0 (sum of the covariant neighbors of phi - 4 phi); r0 = dtau/tau.
 *   b0 Ux, b1 Uy, b2 phi, b3 pi */
IPG_FN void ipg_k_cym_pi(IPG_KERNEL_ARGS) {
    const int n = p.n;
    const int ix = i / n;
    const int iy = i % n;
    const int pX = ipg_coord(ix, 1, n, p.i0) * n + iy;
    const int mX = ipg_coord(ix, -1, n, p.i0) * n + iy;
    const int pY = ix * n + ipg_coord(iy, 1, n, p.i0);
    const int mY = ix * n + ipg_coord(iy, -1, n, p.i0);
    const ipg_su3 Ux = ipg_load(b0, i);
    const ipg_su3 Uy = ipg_load(b1, i);
    const ipg_su3 phi = ipg_load(b2, i);
    const ipg_su3 UxXm1 = ipg_load(b0, mX);
    const ipg_su3 UyYm1 = ipg_load(b1, mY);
    const ipg_su3 phiX = ipg_mul(Ux, ipg_mul_adj(ipg_load(b2, pX), Ux));
    const ipg_su3 phiY = ipg_mul(Uy, ipg_mul_adj(ipg_load(b2, pY), Uy));
    const ipg_su3 phimX = ipg_mul(ipg_adj_mul(UxXm1, ipg_load(b2, mX)), UxXm1);
    const ipg_su3 phimY = ipg_mul(ipg_adj_mul(UyYm1, ipg_load(b2, mY)), UyYm1);
    ipg_su3 pi = ipg_load(b3, i);
    int k;
    for (k = 0; k < 18; k++) {
        pi.e[k] += p.r0
                   * (phiX.e[k] + phimX.e[k] + phiY.e[k] + phimY.e[k]
                      - IPG_C(4.) * phi.e[k]);
    }
    ipg_store(b3, i, pi);
}

/* addEForceSU3() of Evolution.cpp: E += i cp (traceless anti-Hermitian
 * part of a + bSign b) + i cc [phiN, phi], then E is made traceless */
IPG_FN ipg_su3 ipg_add_e_force(
    ipg_su3 E, ipg_su3 a, ipg_su3 b, ipg_real bSign, ipg_su3 phiN, ipg_su3 phi,
    ipg_real cp, ipg_real cc) {
    const ipg_su3 comm = ipg_sub(ipg_mul(phiN, phi), ipg_mul(phi, phiN));
    ipg_su3 mm;
    ipg_real traceThirdIm, etr, eti;
    int row, col, k;
    for (k = 0; k < 18; k++) mm.e[k] = a.e[k] + bSign * b.e[k];
    /* (m_ii - conj(m_ii)) summed over i, over 3: purely imaginary */
    traceThirdIm = IPG_C(2.) * (mm.e[1] + mm.e[9] + mm.e[17]) / IPG_C(3.);
    for (row = 0; row < 3; row++) {
        for (col = 0; col < 3; col++) {
            const int idx = 3 * row + col;
            const int tidx = 3 * col + row;
            /* plaq = m_ij - conj(m_ji) */
            const ipg_real pr = mm.e[2 * idx] - mm.e[2 * tidx];
            ipg_real pim = mm.e[2 * idx + 1] + mm.e[2 * tidx + 1];
            if (row == col) pim -= traceThirdIm;
            /* i cp (pr + i pim) + i cc (cr + i ci) */
            E.e[2 * idx] += -cp * pim - cc * comm.e[2 * idx + 1];
            E.e[2 * idx + 1] += cp * pr + cc * comm.e[2 * idx];
        }
    }
    etr = (E.e[0] + E.e[8] + E.e[16]) / IPG_C(3.);
    eti = (E.e[1] + E.e[9] + E.e[17]) / IPG_C(3.);
    E.e[0] -= etr;
    E.e[1] -= eti;
    E.e[8] -= etr;
    E.e[9] -= eti;
    E.e[16] -= etr;
    E.e[17] -= eti;
    return E;
}

/* The update of the electric fields E1, E2 from the plaquettes and phi.
 *   b0 Ux, b1 Uy, b2 phi, b3 E1, b4 E2;
 *   r0 = tau dtau / (2 g^2), r1 = dtau / tau */
IPG_FN void ipg_k_cym_e(IPG_KERNEL_ARGS) {
    const int n = p.n;
    const int ix = i / n;
    const int iy = i % n;
    const int xp = ipg_coord(ix, 1, n, p.i0);
    const int xm = ipg_coord(ix, -1, n, p.i0);
    const int yp = ipg_coord(iy, 1, n, p.i0);
    const int ym = ipg_coord(iy, -1, n, p.i0);
    const int pX = xp * n + iy;
    const int mX = xm * n + iy;
    const int pY = ix * n + yp;
    const int mY = ix * n + ym;
    const int mXpY = xm * n + yp;
    const int pXmY = xp * n + ym;
    const ipg_su3 Ux = ipg_load(b0, i);
    const ipg_su3 Uy = ipg_load(b1, i);
    const ipg_su3 phi = ipg_load(b2, i);
    /* U12 = Ux(x) Uy(x+ex) Ux(x+ey)^dagger Uy(x)^dagger */
    const ipg_su3 U12 = ipg_mul(
        ipg_mul(Ux, ipg_load(b1, pX)),
        ipg_mul_adj(ipg_adj(ipg_load(b0, pY)), Uy));
    /* U1m2 = Ux(x) Uy(x+ex-ey)^dagger Ux(x-ey)^dagger Uy(x-ey) */
    const ipg_su3 U1m2 = ipg_mul(
        ipg_mul_adj(Ux, ipg_load(b1, pXmY)),
        ipg_adj_mul(ipg_load(b0, mY), ipg_load(b1, mY)));
    /* U2m1 = Uy(x) Ux(x-ex+ey)^dagger Uy(x-ex)^dagger Ux(x-ex) */
    const ipg_su3 U2m1 = ipg_mul(
        ipg_mul_adj(Uy, ipg_load(b0, mXpY)),
        ipg_adj_mul(ipg_load(b1, mX), ipg_load(b0, mX)));
    ipg_su3 phiN = ipg_mul(Ux, ipg_mul_adj(ipg_load(b2, pX), Ux));
    ipg_store(
        b3, i,
        ipg_add_e_force(
            ipg_load(b3, i), U12, U1m2, IPG_C(1.), phiN, phi, p.r0, p.r1));
    phiN = ipg_mul(Uy, ipg_mul_adj(ipg_load(b2, pY), Uy));
    ipg_store(
        b4, i,
        ipg_add_e_force(
            ipg_load(b4, i), U2m1, U12, IPG_C(-1.), phiN, phi, p.r0, p.r1));
}

/* phi += r0 pi; r0 = (tau + dtau/2) dtau. b0 phi, b1 pi; i over reals */
IPG_FN void ipg_k_cym_phi(IPG_KERNEL_ARGS) { b0[i] += p.r0 * b1[i]; }

/* Ux <- exp(i r0 E1) Ux, Uy <- exp(i r0 E2) Uy; r0 = g^2 dtau/(tau+dtau/2).
 *   b0 Ux, b1 Uy, b2 E1, b3 E2 */
IPG_FN void ipg_k_cym_u(IPG_KERNEL_ARGS) {
    ipg_store(
        b0, i,
        ipg_mul(
            ipg_expi_hermitian(ipg_scale(ipg_load(b2, i), p.r0)),
            ipg_load(b0, i)));
    ipg_store(
        b1, i,
        ipg_mul(
            ipg_expi_hermitian(ipg_scale(ipg_load(b3, i), p.r0)),
            ipg_load(b1, i)));
}

/* The kernels, for the launch tables of the backends: X(name). The order
 * defines the kernel ids (DeviceBackend.h). */
#define IPG_KERNEL_LIST(X)        \
    X(ipg_k_copy)                 \
    X(ipg_k_real_to_complex)      \
    X(ipg_k_complex_to_real)      \
    X(ipg_k_scale_by_site)        \
    X(ipg_k_checker_scale)        \
    X(ipg_k_fft_pass)             \
    X(ipg_k_gauss_noise)          \
    X(ipg_k_jimwlk_kernel_noise)  \
    X(ipg_k_jimwlk_rotate)        \
    X(ipg_k_jimwlk_kernel_matrix) \
    X(ipg_k_jimwlk_update)        \
    X(ipg_k_cym_pi)               \
    X(ipg_k_cym_e)                \
    X(ipg_k_cym_phi)              \
    X(ipg_k_cym_u)

#endif /* SRC_DEVICE_IPG_KERNELS_H_ */

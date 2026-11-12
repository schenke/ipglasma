// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include "DeviceDrivers.h"

#include <cmath>
#include <complex>
#include <stdexcept>

#include "../Instrumentation.h"

namespace ipgdev {

namespace {

// a Matrix is 9 complex<double>, i.e. 18 doubles (re, im), which is the
// layout of the device
static_assert(
    sizeof(Matrix) == 18 * sizeof(double),
    "Matrix must be 9 complex<double> without padding");

const double *reals(const std::vector<Matrix> &field) {
    return reinterpret_cast<const double *>(field.data());
}
double *reals(std::vector<Matrix> &field) {
    return reinterpret_cast<double *>(field.data());
}

void checkSites(const std::vector<Matrix> &field, std::size_t sites) {
    if (field.size() != sites) {
        throw std::runtime_error(
            "device driver: a field does not have the size of the lattice");
    }
}

void checkIndexRange(std::size_t sites) {
    // the kernels index with int; the largest index is 2 * 18 * sites (the
    // reals of two matrices per site)
    if (sites * 36 >= 2147483647ULL) {
        throw std::runtime_error(
            "device driver: the lattice is too large for the 32-bit indices "
            "of the device kernels");
    }
}

/**
 * Projects the nearly unitary 3x3 matrices of \p field onto SU(3), in
 * double precision on the host: the unitary factor of the polar
 * decomposition (the closest unitary matrix), by the Newton-Schulz
 * iteration X <- X (3 - X^dagger X) / 2, with the phase of the determinant
 * removed.
 *
 * A device in single precision rounds every product of group elements at
 * 6e-8, and the rounding is not unitary: after the about 300 steps of the
 * JIMWLK evolution X^dagger X differs from 1 by several 1e-6. The host
 * code takes its matrices to be unitary, and the magnetic energy,
 * Nc - Re Tr of a plaquette, is linear in a non-unitary error of the links
 * but quadratic in a unitary one. Without the projection the vacuum
 * around the nuclei gets an energy density of the order of this error,
 * which grows relative to the fields as a^-4 with the lattice size. The
 * unitary part of the rounding error stays; it is a small field.
 */
void projectOnSU3(std::vector<Matrix> &field) {
    typedef std::complex<double> cplx;
    const long sites = static_cast<long>(field.size());
    cplx *all = reinterpret_cast<cplx *>(field.data());
#pragma omp parallel for
    for (long site = 0; site < sites; site++) {
        cplx *x = all + 9 * site;
        for (int iteration = 0; iteration < 4; iteration++) {
            // h = 3 - X^dagger X
            cplx h[9];
            double deviation = 0.;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    cplx sum = 0.;
                    for (int k = 0; k < 3; k++) {
                        sum += std::conj(x[3 * k + i]) * x[3 * k + j];
                    }
                    const cplx unit = (i == j) ? 1. : 0.;
                    deviation += std::norm(sum - unit);
                    h[3 * i + j] = 3. * unit - sum;
                }
            }
            if (deviation < 1e-30) break;
            cplx next[9];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    cplx sum = 0.;
                    for (int k = 0; k < 3; k++) {
                        sum += x[3 * i + k] * h[3 * k + j];
                    }
                    next[3 * i + j] = 0.5 * sum;
                }
            }
            for (int k = 0; k < 9; k++) x[k] = next[k];
        }
        const cplx det = x[0] * (x[4] * x[8] - x[5] * x[7])
                         - x[1] * (x[3] * x[8] - x[5] * x[6])
                         + x[2] * (x[3] * x[7] - x[4] * x[6]);
        const cplx phase = std::polar(1., -std::arg(det) / 3.);
        for (int k = 0; k < 9; k++) x[k] *= phase;
    }
}

}  // namespace

JimwlkDevice::JimwlkDevice(
    Backend &backend, int n, const double *generators, const double *kernel,
    const double *sqrtAlphasNoise, double ds)
    : backend_(backend),
      n_(n),
      sites_(static_cast<std::size_t>(n) * n),
      sqrtDs_(std::sqrt(ds)),
      noiseCoupling_(sqrtAlphasNoise != nullptr) {
    checkIndexRange(sites_);
    if (!backend_.fftAvailable(n)) {
        throw std::runtime_error(
            "device JIMWLK: the FFT of the device is not available for this "
            "lattice size (it needs an even size with prime factors up to "
            "16)");
    }
    wilson_ = backend_.createBuffer(18 * sites_);
    generators_ = backend_.createBuffer(8 * 18);
    kernel_ = backend_.createBuffer(4 * sites_);
    noise_ = backend_.createBuffer(16 * sites_);
    xi_ = backend_.createBuffer(32 * sites_);
    ck_ = backend_.createBuffer(16 * sites_);
    v_ = backend_.createBuffer(36 * sites_);
    w_ = backend_.createBuffer(18 * sites_);
    alphas_ = -1;
    eta_ = -1;
    backend_.upload(generators_, generators, 8 * 18);
    backend_.upload(kernel_, kernel, 4 * sites_);
    if (noiseCoupling_) {
        alphas_ = backend_.createBuffer(sites_);
        eta_ = backend_.createBuffer(32 * sites_);
        backend_.upload(alphas_, sqrtAlphasNoise, sites_);
    }
}

JimwlkDevice::~JimwlkDevice() {
    for (const Buffer buffer :
         {wilson_, generators_, kernel_, noise_, xi_, ck_, v_, w_, alphas_,
          eta_}) {
        if (buffer >= 0) backend_.releaseBuffer(buffer);
    }
}

void JimwlkDevice::upload(const std::vector<Matrix> &wilsonLines) {
    checkSites(wilsonLines, sites_);
    backend_.upload(wilson_, reals(wilsonLines), 18 * sites_);
}

void JimwlkDevice::download(std::vector<Matrix> &wilsonLines) {
    checkSites(wilsonLines, sites_);
    backend_.download(wilson_, reals(wilsonLines), 18 * sites_);
    // the rounding of single precision is not unitary, see projectOnSU3()
    if (backend_.singlePrecision()) projectOnSU3(wilsonLines);
}

void JimwlkDevice::step(const double *noise) {
    backend_.upload(noise_, noise, 16 * sites_);
    evolveWithNoise();
}

void JimwlkDevice::stepWithDeviceNoise(
    unsigned int key0, unsigned int key1, unsigned int stepNumber) {
    Params p;
    p.n = n_;
    // bit patterns; the kernel reads them as unsigned again
    p.i0 = static_cast<int>(key0);
    p.i1 = static_cast<int>(key1);
    p.i2 = static_cast<int>(stepNumber);
    // 16 numbers per site, two per element of the launch
    backend_.launch(kGaussNoise, 8 * sites_, p, {noise_});
    evolveWithNoise();
}

void JimwlkDevice::downloadNoise(std::vector<double> &noise) {
    noise.resize(16 * sites_);
    backend_.download(noise_, noise.data(), 16 * sites_);
}

void JimwlkDevice::evolveWithNoise() {
    Params p;
    p.n = n_;

    // the noise xi and its Fourier transform
    backend_.launch(kRealToComplex, 16 * sites_, p, {noise_, xi_});
    backend_.centeredFft(xi_, n_, 16, 1);

    if (noiseCoupling_) {
        // eta = F^{-1}(sqrt(alpha_s(k)) F(xi)) replaces xi in both factors
        // (JIMWLK::applyNoiseCoupling())
        p.ncomp = 16;
        backend_.launch(kScaleBySite, 16 * sites_, p, {xi_, alphas_});
        backend_.launch(kCopy, 32 * sites_, p, {xi_, eta_});
        backend_.centeredFft(eta_, n_, 16, -1);
        backend_.launch(kComplexToReal, 16 * sites_, p, {eta_, noise_});
    }

    // C(K_i, xi_i^a) = F^{-1}(F(K_x) F(xi_x^a) + F(K_y) F(xi_y^a))
    backend_.launch(kJimwlkKernelNoise, sites_, p, {kernel_, xi_, ck_});
    backend_.centeredFft(ck_, n_, 8, -1);

    // V xi V^dagger, convolved with K
    backend_.launch(
        kJimwlkRotate, sites_, p, {wilson_, noise_, generators_, v_});
    backend_.centeredFft(v_, n_, 18, 1);
    backend_.launch(kJimwlkKernelMatrix, sites_, p, {kernel_, v_, w_});
    backend_.centeredFft(w_, n_, 9, -1);

    // the new Wilson lines
    p.r0 = sqrtDs_;
    backend_.launch(kJimwlkUpdate, sites_, p, {wilson_, w_, ck_, generators_});
}

YangMillsDevice::YangMillsDevice(
    Backend &backend, int n, bool periodic, double g)
    : backend_(backend),
      n_(n),
      sites_(static_cast<std::size_t>(n) * n),
      periodic_(periodic ? 1 : 0),
      g_(g) {
    checkIndexRange(sites_);
    for (Buffer *buffer :
         {&e1_, &e2_, &ux_, &uy_, &pi_, &phi_, &e1Backup_, &e2Backup_,
          &piBackup_}) {
        *buffer = backend_.createBuffer(18 * sites_);
    }
}

YangMillsDevice::~YangMillsDevice() {
    for (const Buffer buffer :
         {e1_, e2_, ux_, uy_, pi_, phi_, e1Backup_, e2Backup_, piBackup_}) {
        backend_.releaseBuffer(buffer);
    }
}

void YangMillsDevice::upload(const Lattice &lat) {
    checkSites(lat.U, sites_);
    backend_.upload(e1_, reals(lat.U), 18 * sites_);
    backend_.upload(e2_, reals(lat.U2), 18 * sites_);
    backend_.upload(ux_, reals(lat.Ux), 18 * sites_);
    backend_.upload(uy_, reals(lat.Uy), 18 * sites_);
    backend_.upload(pi_, reals(lat.Ux2), 18 * sites_);
    backend_.upload(phi_, reals(lat.Uy2), 18 * sites_);
}

void YangMillsDevice::download(Lattice &lat) {
    downloadMomenta(lat);
    backend_.download(ux_, reals(lat.Ux), 18 * sites_);
    backend_.download(uy_, reals(lat.Uy), 18 * sites_);
    backend_.download(phi_, reals(lat.Uy2), 18 * sites_);
    // the links, for the plaquettes of the observables; the fields on the
    // device are not changed
    if (backend_.singlePrecision()) {
        projectOnSU3(lat.Ux);
        projectOnSU3(lat.Uy);
    }
}

void YangMillsDevice::downloadMomenta(Lattice &lat) {
    checkSites(lat.U, sites_);
    backend_.download(e1_, reals(lat.U), 18 * sites_);
    backend_.download(e2_, reals(lat.U2), 18 * sites_);
    backend_.download(pi_, reals(lat.Ux2), 18 * sites_);
}

void YangMillsDevice::evolvePhi(double dtau, double tau) {
    Params p;
    p.n = n_;
    p.r0 = (tau + dtau / 2.) * dtau;
    backend_.launch(kCymPhi, 18 * sites_, p, {phi_, pi_});
}

void YangMillsDevice::evolveU(double dtau, double tau) {
    Params p;
    p.n = n_;
    p.r0 = g_ * g_ * dtau / (tau + dtau / 2.);
    backend_.launch(kCymU, sites_, p, {ux_, uy_, e1_, e2_});
}

void YangMillsDevice::evolvePi(double dtau, double tau) {
    Params p;
    p.n = n_;
    p.i0 = periodic_;
    p.r0 = dtau / tau;
    backend_.launch(kCymPi, sites_, p, {ux_, uy_, phi_, pi_});
}

void YangMillsDevice::evolveE(double dtau, double tau) {
    Params p;
    p.n = n_;
    p.i0 = periodic_;
    p.r0 = tau * dtau / (2. * g_ * g_);
    p.r1 = dtau / tau;
    backend_.launch(kCymE, sites_, p, {ux_, uy_, phi_, e1_, e2_});
}

void YangMillsDevice::step(double dtau, double tau, bool updateCoordinates) {
    evolvePi(dtau, tau);
    evolveE(dtau, tau);
    if (updateCoordinates) {
        evolvePhi(dtau, tau);
        evolveU(dtau, tau);
    }
}

void YangMillsDevice::backupMomenta() {
    Params p;
    p.n = n_;
    backend_.launch(kCopy, 18 * sites_, p, {e1_, e1Backup_});
    backend_.launch(kCopy, 18 * sites_, p, {e2_, e2Backup_});
    backend_.launch(kCopy, 18 * sites_, p, {pi_, piBackup_});
}

void YangMillsDevice::restoreMomenta() {
    Params p;
    p.n = n_;
    backend_.launch(kCopy, 18 * sites_, p, {e1Backup_, e1_});
    backend_.launch(kCopy, 18 * sites_, p, {e2Backup_, e2_});
    backend_.launch(kCopy, 18 * sites_, p, {piBackup_, pi_});
}

}  // namespace ipgdev

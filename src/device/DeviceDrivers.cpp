// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include "DeviceDrivers.h"

#include <algorithm>
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

}  // namespace

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
            "lattice size (the FFT of the kernels needs an even size with "
            "prime factors up to 32, cuFFT an even size)");
    }
    wilson_ = DeviceBuffer(backend_, 18 * sites_);
    generators_ = DeviceBuffer(backend_, 8 * 18);
    kernel_ = DeviceBuffer(backend_, 4 * sites_);
    noise_ = DeviceBuffer(backend_, 16 * sites_);
    zxi_ = DeviceBuffer(backend_, 16 * sites_);
    zv_ = DeviceBuffer(backend_, 16 * sites_);
    ck_ = DeviceBuffer(backend_, 8 * sites_);
    w_ = DeviceBuffer(backend_, 8 * sites_);
    backend_.upload(generators_, generators, 8 * 18);
    // K with K(-k) = K(k)^*, so that the convolutions of real fields are
    // real (ipg_k_jimwlk_combine): (K(k) + K(-k)^*) / 2; the CPU code takes
    // the real part of its convolutions instead, the same up to rounding
    std::vector<double> symmetric(4 * sites_);
    for (int ix = 0; ix < n; ix++) {
        for (int iy = 0; iy < n; iy++) {
            const std::size_t k = static_cast<std::size_t>(ix) * n + iy;
            const std::size_t m =
                static_cast<std::size_t>((n - ix) % n) * n + (n - iy) % n;
            for (int d = 0; d < 2; d++) {
                symmetric[4 * k + 2 * d] =
                    0.5 * (kernel[4 * k + 2 * d] + kernel[4 * m + 2 * d]);
                symmetric[4 * k + 2 * d + 1] =
                    0.5
                    * (kernel[4 * k + 2 * d + 1] - kernel[4 * m + 2 * d + 1]);
            }
        }
    }
    backend_.upload(kernel_, symmetric.data(), 4 * sites_);
    if (noiseCoupling_) {
        alphas_ = DeviceBuffer(backend_, sites_);
        eta_ = DeviceBuffer(backend_, 16 * sites_);
        backend_.upload(alphas_, sqrtAlphasNoise, sites_);
    }
}

JimwlkDevice::~JimwlkDevice() { backend_.releaseCaches(); }

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

    // the transform of the noise, two real fields per complex one:
    // Z^a = F(xi_x^a + i xi_y^a)
    backend_.launch(kJimwlkPackNoise, sites_, p, {noise_, zxi_});
    backend_.centeredFft(zxi_, n_, 8, 1);

    if (noiseCoupling_) {
        // eta = F^{-1}(sqrt(alpha_s(k)) F(xi)) replaces xi in both factors
        // (JIMWLK::applyNoiseCoupling()); sqrt(alpha_s(k)) is real and
        // even, so the pairs stay real
        p.ncomp = 8;
        backend_.launch(kScaleBySite, 8 * sites_, p, {zxi_, alphas_});
        backend_.launch(kCopy, 16 * sites_, p, {zxi_, eta_});
        backend_.centeredFft(eta_, n_, 8, -1);
        backend_.launch(kJimwlkUnpackNoise, sites_, p, {eta_, noise_});
    }

    // CK^a = K_x * xi_x^a + K_y * xi_y^a (convolutions), real
    backend_.launch(kJimwlkCombine, sites_, p, {kernel_, zxi_, ck_});
    backend_.centeredFft(ck_, n_, 4, -1);

    // the coefficients of V xi V^dagger, convolved with K: W^a
    backend_.launch(
        kJimwlkRotate, sites_, p, {wilson_, noise_, generators_, zv_});
    backend_.centeredFft(zv_, n_, 8, 1);
    backend_.launch(kJimwlkCombine, sites_, p, {kernel_, zv_, w_});
    backend_.centeredFft(w_, n_, 4, -1);

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
    for (DeviceBuffer *buffer :
         {&e1_, &e2_, &ux_, &uy_, &pi_, &phi_, &e1Backup_, &e2Backup_,
          &piBackup_}) {
        *buffer = DeviceBuffer(backend_, 18 * sites_);
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

// ---------------------------------------------------------------------

GaugeFixDevice::GaugeFixDevice(Backend &backend, int n)
    : backend_(backend), n_(n), sites_(static_cast<std::size_t>(n) * n) {
    checkIndexRange(sites_);
    if (!backend_.fftAvailable(n)) {
        throw std::runtime_error(
            "device gauge fixing: the FFT of the device is not available "
            "for this lattice size");
    }
    ux0_ = DeviceBuffer(backend_, 18 * sites_);
    uy0_ = DeviceBuffer(backend_, 18 * sites_);
    ux_ = DeviceBuffer(backend_, 18 * sites_);
    uy_ = DeviceBuffer(backend_, 18 * sites_);
    g_ = DeviceBuffer(backend_, 18 * sites_);
    chi_ = DeviceBuffer(backend_, 8 * sites_);
    residual_ = DeviceBuffer(backend_, sites_);
    // the residual is summed on the device into partialCount_ partial
    // sums, which the host adds up: a small download per iteration
    partialCount_ = std::min<std::size_t>(sites_, 1024);
    partial_ = DeviceBuffer(backend_, partialCount_);
    partialSums_.resize(partialCount_);
    poisson_ = DeviceBuffer(backend_, sites_);

    // the factor of the Poisson step of GaugeFix::fftChi()
    std::vector<double> factor(sites_);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            const double kx = std::sin(
                M_PI
                * (-0.5 + static_cast<double>(i) / static_cast<double>(n)));
            const double ky = std::sin(
                M_PI
                * (-0.5 + static_cast<double>(j) / static_cast<double>(n)));
            const double kt2 = 4. * (kx * kx + ky * ky);
            factor[static_cast<std::size_t>(i) * n + j] =
                -1.5 * (1. / (kt2 + 1e-9));
        }
    }
    // The zero mode of chi is the sum of the divergence over the lattice,
    // zero up to rounding, and GaugeFix::fftChi() multiplies it by 1.5e9.
    // In double precision that is a global rotation of 1e-10 per
    // iteration; in single precision it would be of order 0.1, harmless
    // for the gauge condition but pointless: leave the zero mode out.
    if (backend_.singlePrecision()) {
        factor[static_cast<std::size_t>(n / 2) * n + n / 2] = 0.;
    }
    backend_.upload(poisson_, factor.data(), sites_);
}

GaugeFixDevice::~GaugeFixDevice() { backend_.releaseCaches(); }

void GaugeFixDevice::upload(
    const std::vector<Matrix> &ux, const std::vector<Matrix> &uy) {
    checkSites(ux, sites_);
    checkSites(uy, sites_);
    backend_.upload(ux0_, reals(ux), 18 * sites_);
    backend_.upload(uy0_, reals(uy), 18 * sites_);
    Params p;
    p.n = n_;
    backend_.launch(kCopy, 18 * sites_, p, {ux0_, ux_});
    backend_.launch(kCopy, 18 * sites_, p, {uy0_, uy_});
    const std::vector<Matrix> one(sites_, Matrix(1.));
    backend_.upload(g_, reals(one), 18 * sites_);
}

double GaugeFixDevice::divergence() {
    Params p;
    p.n = n_;
    backend_.launch(kGfDivergence, sites_, p, {ux_, uy_, chi_, residual_});
    // in another order than on the host: the residual differs from that of
    // the host in the last digits; the host checks it again after the
    // device (GaugeFix::fftChi())
    p.i0 = static_cast<int>(sites_);
    p.i1 = static_cast<int>(partialCount_);
    backend_.launch(kPartialSum, partialCount_, p, {residual_, partial_});
    backend_.download(partial_, partialSums_.data(), partialCount_);
    double sum = 0.;
    for (const double partial : partialSums_) sum += partial;
    return sum / static_cast<double>(sites_);
}

void GaugeFixDevice::relax() {
    Params p;
    p.n = n_;
    // chi packed as 4 complex components (ipg_gf_unpack)
    backend_.centeredFft(chi_, n_, 4, 1);
    p.ncomp = 4;
    backend_.launch(kScaleBySite, 4 * sites_, p, {chi_, poisson_});
    backend_.centeredFft(chi_, n_, 4, -1);
    backend_.launch(kGfUpdate, sites_, p, {chi_, g_});
    backend_.launch(kGfTransform, sites_, p, {g_, ux0_, uy0_, ux_, uy_});
}

void GaugeFixDevice::download(std::vector<Matrix> &g) {
    g.resize(sites_);
    backend_.download(g_, reals(g), 18 * sites_);
    // Also in double precision: the zero mode of chi, rounding times
    // 1.5e9, is not exactly Hermitian and traceless, which leaves G with a
    // phase and a non-unitary part of 1e-9.
    projectOnSU3(g);
}

// ---------------------------------------------------------------------

WilsonLineDevice::WilsonLineDevice(
    Backend &backend, int n, const double *momentumKernel)
    : backend_(backend), n_(n), sites_(static_cast<std::size_t>(n) * n) {
    checkIndexRange(sites_);
    if (!backend_.fftAvailable(n)) {
        throw std::runtime_error(
            "device Wilson lines: the FFT of the device is not available "
            "for this lattice size");
    }
    wilson_ = DeviceBuffer(backend_, 18 * sites_);
    kernel_ = DeviceBuffer(backend_, sites_);
    scale_ = DeviceBuffer(backend_, sites_);
    gauss_ = DeviceBuffer(backend_, 8 * sites_);
    rho_ = DeviceBuffer(backend_, 8 * sites_);
    backend_.upload(kernel_, momentumKernel, sites_);
}

WilsonLineDevice::~WilsonLineDevice() { backend_.releaseCaches(); }

void WilsonLineDevice::upload(
    const std::vector<Matrix> &wilsonLines, const double *scale) {
    checkSites(wilsonLines, sites_);
    backend_.upload(wilson_, reals(wilsonLines), 18 * sites_);
    backend_.upload(scale_, scale, sites_);
}

void WilsonLineDevice::download(std::vector<Matrix> &wilsonLines) {
    checkSites(wilsonLines, sites_);
    backend_.download(wilson_, reals(wilsonLines), 18 * sites_);
    if (backend_.singlePrecision()) projectOnSU3(wilsonLines);
}

void WilsonLineDevice::sheet(const double *gauss) {
    backend_.upload(gauss_, gauss, 8 * sites_);
    multiplySheet();
}

void WilsonLineDevice::sheetWithDeviceNoise(
    unsigned int key0, unsigned int key1, unsigned int sheetNumber) {
    Params p;
    p.n = n_;
    p.i0 = static_cast<int>(key0);
    p.i1 = static_cast<int>(key1);
    p.i2 = static_cast<int>(sheetNumber);
    // 8 numbers per site, two per element of the launch
    backend_.launch(kGaussNoise, 4 * sites_, p, {gauss_});
    multiplySheet();
}

void WilsonLineDevice::multiplySheet() {
    Params p;
    p.n = n_;
    // the 8 real charges of a site as 4 complex components
    backend_.launch(kWilsonCharge, 8 * sites_, p, {gauss_, scale_, rho_});
    backend_.centeredFft(rho_, n_, 4, 1);
    p.ncomp = 4;
    backend_.launch(kScaleBySite, 4 * sites_, p, {rho_, kernel_});
    backend_.centeredFft(rho_, n_, 4, -1);
    backend_.launch(kWilsonExponent, sites_, p, {rho_, wilson_});
}

}  // namespace ipgdev

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include "JIMWLK.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Init.h"
#include "Instrumentation.h"
#include "RunningCoupling.h"
#include "WilsonLineIO.h"

using PhysConst::invHbarc;
using PhysConst::isClose;
using PhysConst::Nc;
using PhysConst::Nc2m1;
using PhysConst::smallEps;

JIMWLK::JIMWLK(Parameters &param, Group *group, Lattice *lat, Random *random)
    : param_(param),
      Ngrid_(param.lattice.size),
      Ncells_(param.lattice.size * param.lattice.size) {
    nn_[0] = param_.lattice.size;
    nn_[1] = param_.lattice.size;

    fft_ptr_ = std::make_shared<FFT>(nn_);

    group_ptr_ = group;
    // with colorChargeSeed, the configurations of a fixed geometry differ
    // also in the noise of the evolution, which then comes from a stream of
    // its own, like the color charges (issue #127); without it, from the
    // stream of the run
    if (param_.random.colorChargeSeed >= 0) {
        noiseRandom_.init_genrand64(noiseRandomSeed(
            param_.random.colorChargeSeed, param_.event.eventId));
        random_ptr_ = &noiseRandom_;
    } else {
        random_ptr_ = random;
    }
    lat_ptr_ = lat;

    initializeK();
    initializeNoise();
    // evolutionStep() dereferences VxsiVx_/VxsiVy_ unconditionally, so these
    // must always be allocated -- a since-removed input flag used to gate
    // this allocation conditionally (matching the pre-refactor new[]-based
    // version), leaving an empty view here whenever that flag was off, so
    // the first evolutionStep() call indexed past an empty vector.
    VxsiVxStorage_.assign(Ncells_, Matrix(0.));
    VxsiVyStorage_.assign(Ncells_, Matrix(0.));
    VxsiVx_.resize(Ncells_);
    VxsiVy_.resize(Ncells_);
    for (int i = 0; i < Ncells_; i++) {
        VxsiVx_[i] = &VxsiVxStorage_[i];
        VxsiVy_[i] = &VxsiVyStorage_[i];
    }
}

void JIMWLK::initializeK() {
    if (initializedK_) {
        return;
    }
    // sized to its final length up front (always x,y) so it never has to
    // grow/reallocate via push_back; entries default to (0,0)
    K_storage_.assign(Ncells_, std::vector<std::complex<double> >(2));
    K_.resize(Ncells_);
    for (int i = 0; i < Ncells_; i++) {
        K_[i] = &K_storage_[i];
    }

    // with the noise running coupling alpha_s(k) is in the noise instead
    const bool alphasInKernel = param_.jimwlk.runningCoupling != 2
                                || !isClose(param_.jimwlk.alphaS, 0.);
    if (!alphasInKernel) {
        sqrtAlphasNoise_.resize(Ncells_);
#pragma omp parallel for
        for (int pos = 0; pos < Ncells_; pos++) {
            // same (shifted) site <-> momentum mapping as the FFT output
            const double kx =
                (latticeX(pos, Ngrid_) - static_cast<double>(Ngrid_) / 2.)
                / Ngrid_;
            const double ky =
                (latticeY(pos, Ngrid_) - static_cast<double>(Ngrid_) / 2.)
                / Ngrid_;
            sqrtAlphasNoise_[pos] = sqrt(getAlphasMomentum(kx, ky));
        }
    }

    // set once here rather than on every getMassRegulator() call below
    gsl_set_error_handler_off();

#pragma omp parallel for
    for (int pos = 0; pos < Ncells_; pos++) {
        double x = latticeX(pos, Ngrid_) - static_cast<double>(Ngrid_) / 2.;
        double y = latticeY(pos, Ngrid_) - static_cast<double>(Ngrid_) / 2.;
        x /= Ngrid_;
        y /= Ngrid_;
        double r2 = x * x + y * y;
        if (r2 < smallEps) {
            continue;  // K_[pos] is already (0, 0)
        }
        double mass_regulator = getMassRegulator(x, y);
        double alphas_sqroot = alphasInKernel ? sqrt(getAlphas(x, y)) : 1.;
        // discretization without singularities
        double tmpk1 = cos(M_PI * y) * sin(2. * M_PI * x) / (2. * M_PI);
        double tmpk2 = cos(M_PI * x) * sin(2. * M_PI * y) / (2. * M_PI);

        // Regulate long distance tails, does nothing if m=0
        double sin_x = sin(M_PI * x) / M_PI;
        double sin_y = sin(M_PI * y) / M_PI;
        double denom = sin_x * sin_x + sin_y * sin_y;
        double factor = alphas_sqroot * mass_regulator / Ngrid_ / denom;

        (*K_[pos])[0] = tmpk1 * factor;
        (*K_[pos])[1] = tmpk2 * factor;
    }
    fft_ptr_->fftnVector(K_.data(), K_.data(), nn_, 1);
    initializedK_ = true;
}

unsigned long long JIMWLK::noiseRandomSeed(
    long long colorChargeSeed, int eventId) {
    return splitmix64(Init::colorChargeRandomSeed(colorChargeSeed, eventId));
}

double JIMWLK::getMassRegulator(const double x, const double y) const {
    // if m suppresses long distance tails,
    // K is multiplied by this, which is m*r*K_1(m*r)
    double mass_regulator = 1.0;
    double m = param_.jimwlk.mass;
    if (m < smallEps) {
        return mass_regulator;
    }
    double length = param_.lattice.L;

    // Lattice units
    // Here x is [-N/2, N/2]
    double lat_x = sin(M_PI * x) / (M_PI);
    double lat_y = sin(M_PI * y) / (M_PI);
    // lat_x and lat_y are now in [-1/2,1/2] as x/nn[0] is in [-N/2, N/2]
    double lat_r = sqrt(lat_x * lat_x + lat_y * lat_y) * Ngrid_;
    // lat_r now tells how many lattice units the distance is
    double a = length / Ngrid_;
    double lat_m = m * a * invHbarc;
    double bessel_argument = lat_m * lat_r;
    // use gsl bessel function to be compatible with the AppleClang compiler
    // (assumes the GSL error handler was already disabled by the caller,
    // since this runs once per cell from initializeK())
    gsl_sf_result bes;
    int status = gsl_sf_bessel_K1_e(bessel_argument, &bes);
    if (status != GSL_SUCCESS) {
        mass_regulator = 0.0;
    } else {
        mass_regulator = bessel_argument * bes.val;
    }
    return mass_regulator;
}

double JIMWLK::getAlphas(const double x, const double y) const {
    // Fixed coupling: alpha_s is already absorbed into the evolution step
    // count in evolution() (ds = alpha_s dy / pi^2), so the kernel must not
    // carry it again.
    if (!isClose(param_.jimwlk.alphaS, 0.)) {
        return 1.0;
    }

    const double length = param_.lattice.L;
    double phys_x = x * length;  // in fm
    double phys_y = y * length;
    double phys_r2 = phys_x * phys_x + phys_y * phys_y;

    // at the scale 2/r, in GeV
    return computeAlphaS(
        param_.jimwlk.mu0, param_.jimwlk.c, param_.jimwlk.LambdaQCD,
        param_.coupling.nFlavors, 2. / (sqrt(phys_r2) * invHbarc));
}

double JIMWLK::getAlphasMomentum(const double kx, const double ky) const {
    const double a = param_.lattice.L / Ngrid_;  // in fm
    const double sin_x = sin(M_PI * kx);
    const double sin_y = sin(M_PI * ky);
    // lattice momentum squared, in GeV^2
    const double k2 =
        4. * (sin_x * sin_x + sin_y * sin_y) / (a * a * invHbarc * invHbarc);

    return computeAlphaS(
        param_.jimwlk.mu0, param_.jimwlk.c, param_.jimwlk.LambdaQCD,
        param_.coupling.nFlavors, sqrt(k2));
}

void JIMWLK::initializeNoise() {
    if (initializedNoise_) {
        return;
    }
    // one contiguous allocation per array instead of Ncells_ separate
    // ones, so consecutive cells are adjacent in memory
    xi_data_.assign(Ncells_ * 2 * Nc2m1, std::complex<double>());
    xi2_data_.assign(Ncells_ * 2 * Nc2m1, std::complex<double>());
    CKxi_data_.assign(Ncells_ * Nc2m1, std::complex<double>());

    xi_.resize(Ncells_);
    xi2_.resize(Ncells_);
    CKxi_.resize(Ncells_);
    for (int i = 0; i < Ncells_; i++) {
        xi_[i] = xi_data_.data() + i * 2 * Nc2m1;
        xi2_[i] = xi2_data_.data() + i * 2 * Nc2m1;
        CKxi_[i] = CKxi_data_.data() + i * Nc2m1;
    }
    initializedNoise_ = true;
}

void JIMWLK::evolution() {
    IPG_PROFILE_SCOPE("jimwlk.evolution");
    initializeNoise();

    // Calculate evolution steps for different nuclei, from the x of the
    // read Wilson lines, otherwise from jimwlkInitialX
    const double x0A = param_.xBeforeJimwlk(NucleusRole::Projectile);
    const double x0B = param_.xBeforeJimwlk(NucleusRole::Target);
    double ds = param_.jimwlk.Ds;
    bool saveSnapshots = param_.jimwlk.saveSnapshots;
    std::vector<double> xSnapshotList = param_.jimwlk.xSnapshotList;
    double dlogx = M_PI * M_PI * ds;
    int steps_1 = 0;
    int steps_2 = 0;
    double as = param_.jimwlk.alphaS;
    if (!isClose(as, 0.)) {
        // Fixed coupling
        steps_1 = static_cast<int>(
            as * std::log(x0A / param_.colorCharge.projectileX)
                / (M_PI * M_PI * ds)
            + 0.5);
        steps_2 = static_cast<int>(
            as * std::log(x0B / param_.colorCharge.targetX) / (M_PI * M_PI * ds)
            + 0.5);
        dlogx = M_PI * M_PI * ds / as;
    } else {
        // Running coupling
        steps_1 = static_cast<int>(
            std::log(x0A / param_.colorCharge.projectileX) / (M_PI * M_PI * ds)
            + 0.5);
        steps_2 = static_cast<int>(
            std::log(x0B / param_.colorCharge.targetX) / (M_PI * M_PI * ds)
            + 0.5);
    }

    runEvolutionLoop(
        NucleusRole::Projectile, steps_1, x0A, dlogx, saveSnapshots,
        xSnapshotList);
    runEvolutionLoop(
        NucleusRole::Target, steps_2, x0B, dlogx, saveSnapshots, xSnapshotList);
}

void JIMWLK::runEvolutionLoop(
    NucleusRole nucleus, int steps, double x0, double dlogx, bool saveSnapshots,
    const std::vector<double> &xSnapshotList) {
    const std::string label =
        (nucleus == NucleusRole::Projectile) ? "projectile" : "target";
    messager_ << "[JIMWLK::evolution]: Evolving " << label
              << ", evolution steps " << steps;
    messager_.flush("info");

    std::vector<int> snapshotStep;
    if (saveSnapshots) {
        snapshotStep = snapshotSteps(xSnapshotList, x0, dlogx, steps);
        // the Wilson lines before and after the evolution have their own
        // files anyway (written by Init::setV() and main, or read), so a
        // snapshot under one of these names would only write them again
        const double finalX = (nucleus == NucleusRole::Projectile)
                                  ? param_.colorCharge.projectileX
                                  : param_.colorCharge.targetX;
        const std::string initialFile =
            WilsonLineIO::fileName(&param_, x0, nucleus);
        const std::string finalFile =
            WilsonLineIO::fileName(&param_, finalX, nucleus);
        for (std::size_t i = 0; i < xSnapshotList.size(); i++) {
            if (snapshotStep[i] < 0) {
                messager_ << "[JIMWLK::evolution]: x = " << xSnapshotList[i]
                          << " of jimwlkXSnapshotList is not between " << x0
                          << " and " << x0 * exp(-steps * dlogx) << " (the "
                          << label
                          << " evolution); no snapshot is saved for it.";
                messager_.flush("warning");
                continue;
            }
            const std::string file =
                WilsonLineIO::fileName(&param_, xSnapshotList[i], nucleus);
            if ((snapshotStep[i] == 0 && file == initialFile)
                || (snapshotStep[i] == steps && file == finalFile)) {
                messager_ << "[JIMWLK::evolution]: x = " << xSnapshotList[i]
                          << " of jimwlkXSnapshotList is the "
                          << (snapshotStep[i] == 0 ? "initial" : "final")
                          << " x of the " << label
                          << " evolution, whose Wilson lines are in " << file
                          << " anyway; no separate snapshot is saved.";
                messager_.flush("info");
                snapshotStep[i] = -1;
            }
        }
    }
    // writes the snapshots due after `step` steps, named with their
    // requested x
    auto saveSnapshotsAfter = [&](int step) {
        for (std::size_t i = 0; i < snapshotStep.size(); i++) {
            if (snapshotStep[i] == step) {
                WilsonLineIO().write(
                    lat_ptr_, &param_, nucleus, xSnapshotList[i]);
            }
        }
    };

    saveSnapshotsAfter(0);

    // on a device (IPGLASMA_DEVICE), the Wilson lines of the nucleus stay
    // there during the evolution and come back for the snapshots and at
    // the end
    std::vector<Matrix> &wilsonLines =
        (nucleus == NucleusRole::Projectile) ? lat_ptr_->U : lat_ptr_->U2;
    std::unique_ptr<ipgdev::JimwlkDevice> device;
    if (steps > 0 && ipgdev::runBackend() != nullptr
        && ipgdev::runConfig().jimwlk) {
        ipgdev::Backend &backend = *ipgdev::runBackend();
        messager_ << "[JIMWLK::evolution]: on the device " << backend.name()
                  << " (" << backend.description() << ")";
        messager_.flush("info");
        device = makeDevice(backend);
        device->upload(wilsonLines);
    }
    // IPGLASMA_DEVICE_NOISE=device: the noise comes from the generator of
    // the device, with a key from the random numbers of the run
    const bool deviceNoise = device && ipgdev::runConfig().noiseOnDevice;
    unsigned int noiseKey0 = 0;
    unsigned int noiseKey1 = 0;
    if (deviceNoise) {
        drawDeviceNoiseKey(noiseKey0, noiseKey1);
        messager_ << "[JIMWLK::evolution]: the noise is drawn on the device "
                     "(IPGLASMA_DEVICE_NOISE=device): this is not the event "
                     "a CPU run gives for this seed.";
        messager_.flush("info");
    }
    auto snapshotDueAfter = [&](int step) {
        return std::find(snapshotStep.begin(), snapshotStep.end(), step)
               != snapshotStep.end();
    };

    for (int ids = 0; ids < steps; ids++) {
        // steps (from user-configurable JIMWLK parameters, no lower bound
        // enforced) can be under 10, making this 0; guard against the
        // resulting integer modulo-by-zero below.
        int printSteps = std::max(1, steps / 10);
        if (ids % printSteps == 0) {
            messager_ << "[JIMWLK::evolution]: Step " << ids;
            messager_.flush("info");
        }
        if (device) {
            if (deviceNoise) {
                IPG_PROFILE_SCOPE("jimwlk.step");
                IPG_PROFILE_SCOPE("jimwlk.device_step");
                device->stepWithDeviceNoise(
                    noiseKey0, noiseKey1, static_cast<unsigned int>(ids));
            } else {
                evolutionStepOnDevice(*device);
            }
            if (ids + 1 == steps || snapshotDueAfter(ids + 1)) {
                IPG_PROFILE_SCOPE("jimwlk.device_download");
                device->download(wilsonLines);
            }
        } else {
            evolutionStep(nucleus);
        }
        saveSnapshotsAfter(ids + 1);
    }
    messager_ << "[JIMWLK::evolution]: Done.";
    messager_.flush("info");
}

std::vector<int> JIMWLK::snapshotSteps(
    const std::vector<double> &xSnapshotList, double x0, double dlogx,
    int steps) {
    std::vector<int> result;
    for (const double x : xSnapshotList) {
        // the number of steps to reach x exactly
        const double exact = std::log(x0 / x) / dlogx;
        if (!(exact >= -0.5 && exact <= steps + 0.5)) {
            result.push_back(-1);
        } else {
            const int step = static_cast<int>(std::floor(exact + 0.5));
            result.push_back(std::min(std::max(step, 0), steps));
        }
    }
    return result;
}

void JIMWLK::applyNoiseCoupling() {
    if (sqrtAlphasNoise_.empty()) {
        return;
    }
    IPG_PROFILE_SCOPE("jimwlk.noise_coupling");
    // eta = F^{-1}(sqrt(alpha_s(k)) F(xi)), whose correlator is the Fourier
    // transform of alpha_s(k) (Eq. (24) of Lappi:2012vw); xi then holds
    // F(eta) and xi2 eta
#pragma omp parallel for
    for (int i = 0; i < Ncells_; i++) {
        for (int n = 0; n < 2 * Nc2m1; n++) {
            xi_[i][n] *= sqrtAlphasNoise_[i];
        }
    }
    fft_ptr_->fftnArray(xi_.data(), xi2_.data(), nn_, -1, 2 * Nc2m1);
    // alpha_s(k) = alpha_s(-k), so eta is real up to rounding
#pragma omp parallel for
    for (int i = 0; i < Ncells_; i++) {
        for (int n = 0; n < 2 * Nc2m1; n++) {
            xi2_[i][n] = std::complex<double>(xi2_[i][n].real(), 0.);
        }
    }
}

std::vector<double> JIMWLK::correlatedNoise(const std::vector<double> &xi) {
    for (int i = 0; i < Ncells_; i++) {
        for (int n = 0; n < 2 * Nc2m1; n++) {
            xi2_[i][n] = std::complex<double>(xi[i], 0.);
        }
    }
    fft_ptr_->fftnArray(xi2_.data(), xi_.data(), nn_, 1, 2 * Nc2m1);
    applyNoiseCoupling();
    std::vector<double> eta(Ncells_);
    for (int i = 0; i < Ncells_; i++) {
        eta[i] = xi2_[i][0].real();
    }
    return eta;
}

std::unique_ptr<ipgdev::JimwlkDevice> JIMWLK::makeDevice(
    ipgdev::Backend &backend) {
    std::vector<double> generators(8 * 18);
    for (int a = 0; a < Nc2m1; a++) {
        const Matrix &Ta = group_ptr_->getT(a);
        for (int k = 0; k < 9; k++) {
            generators[18 * a + 2 * k] = Ta(k).real();
            generators[18 * a + 2 * k + 1] = Ta(k).imag();
        }
    }
    std::vector<double> kernel(4 * static_cast<std::size_t>(Ncells_));
    for (int i = 0; i < Ncells_; i++) {
        for (int d = 0; d < 2; d++) {
            kernel[4 * static_cast<std::size_t>(i) + 2 * d] =
                (*K_[i])[d].real();
            kernel[4 * static_cast<std::size_t>(i) + 2 * d + 1] =
                (*K_[i])[d].imag();
        }
    }
    return std::unique_ptr<ipgdev::JimwlkDevice>(new ipgdev::JimwlkDevice(
        backend, Ngrid_, generators.data(), kernel.data(),
        sqrtAlphasNoise_.empty() ? nullptr : sqrtAlphasNoise_.data(),
        param_.jimwlk.Ds));
}

void JIMWLK::drawDeviceNoiseKey(unsigned int &key0, unsigned int &key1) {
    const unsigned long long word = random_ptr_->genrand64_int64();
    key0 = static_cast<unsigned int>(word & 0xFFFFFFFFULL);
    key1 = static_cast<unsigned int>(word >> 32);
}

void JIMWLK::evolutionStepOnDevice(ipgdev::JimwlkDevice &device) {
    IPG_PROFILE_SCOPE("jimwlk.step");
    // the same draw as in evolutionStep(), so that the random numbers of
    // the run are consumed identically
    const std::size_t count = static_cast<std::size_t>(Ncells_) * 2 * Nc2m1;
    {
        IPG_PROFILE_SCOPE("jimwlk.gauss_bulk");
        gaussNoise_.resize(count);
        random_ptr_->gaussBulk(gaussNoise_.data(), count, gaussNoiseScratch_);
    }
    IPG_PROFILE_SCOPE("jimwlk.device_step");
    device.step(gaussNoise_.data());
}

void JIMWLK::evolutionStep(NucleusRole nucleus) {
    IPG_PROFILE_SCOPE("jimwlk.step");
    const bool evolveProjectile = (nucleus == NucleusRole::Projectile);
    std::vector<Matrix> &wilsonLines =
        evolveProjectile ? lat_ptr_->U : lat_ptr_->U2;

    const complex<double> I(0., 1.);
    const double ds_sqrt = std::sqrt(param_.jimwlk.Ds);
    const complex<double> negI_dssqrt = -I * ds_sqrt;
    const complex<double> posI_dssqrt = I * ds_sqrt;

    // generate random Gaussian noise in every cell for Nc^2-1 color
    // components and 2 spatial components x and y.
    // gaussBulk reproduces the exact same value stream, in the exact same
    // order, as calling random_ptr_->gauss() this many times in a row (see
    // its own comments in Random.cpp), so this is a drop-in replacement for
    // what used to be a serial scalar loop -- the scatter into xi2_ (which,
    // unlike the draw itself, touches no shared RNG state) can then run in
    // parallel.
    {
        IPG_PROFILE_SCOPE("jimwlk.gauss_bulk");
        const std::size_t count = static_cast<std::size_t>(Ncells_) * 2 * Nc2m1;
        gaussNoise_.resize(count);
        random_ptr_->gaussBulk(gaussNoise_.data(), count, gaussNoiseScratch_);
#pragma omp parallel for
        for (int i = 0; i < Ncells_; i++) {
            const double *cellNoise =
                gaussNoise_.data() + static_cast<std::size_t>(i) * 2 * Nc2m1;
            for (int n = 0; n < 2 * Nc2m1; n++) {
                xi2_[i][n] = std::complex<double>(cellNoise[n], 0.);
            }
        }
    }

    // Fourier transform the noise in xi2 into xi, leaving xi2 unchanged
    fft_ptr_->fftnArray(xi2_.data(), xi_.data(), nn_, 1, 2 * Nc2m1);

    // with the noise running coupling, the correlated noise eta replaces
    // xi in both arrays, i.e. in both factors of the step below
    applyNoiseCoupling();

    // now compute C(K_i,xi_i^a) = F^{-1}(F(K_i)F(xi_i^a))
    //                           = F^{-1}(F(K_x)F(xi_x^a)+F(K_y)F(xi_y^a))
#pragma omp parallel for
    for (int i = 0; i < Ncells_; i++) {
        for (int n = 0; n < Nc2m1; n++) {
            CKxi_[i][n] =
                (*K_[i])[0] * xi_[i][n] + (*K_[i])[1] * xi_[i][n + Nc2m1];
            // product of x components + product of y components
        }
    }

    // now CKxi contains C(K_i,xi_i^a) - it is a vector with a components
    fft_ptr_->fftnArray(CKxi_.data(), CKxi_.data(), nn_, -1, Nc2m1);

#pragma omp parallel for
    for (int i = 0; i < Ncells_; i++) {
        *VxsiVx_[i] = zero_;
        *VxsiVy_[i] = zero_;
        const Matrix &U = wilsonLines[i];
        for (int a = 0; a < Nc2m1; a++) {
            Matrix UTa = U * group_ptr_->getT(a);
            // UTa * U^dagger, folding in the conjugate transpose of U
            // instead of building it explicitly (Nc=3 specialization)
            Matrix UTUconj = UTa.prodABconj(UTa, U);
            *VxsiVx_[i] += xi2_[i][a] * UTUconj;
            *VxsiVy_[i] += xi2_[i][a + Nc2m1] * UTUconj;
        }
    }

    // FFT V xi V
    fft_ptr_->fftn(VxsiVx_.data(), VxsiVx_.data(), nn_, 1);
    fft_ptr_->fftn(VxsiVy_.data(), VxsiVy_.data(), nn_, 1);

#pragma omp parallel for
    for (int i = 0; i < Ncells_; i++) {
        *VxsiVx_[i] = (*K_[i])[0] * (*VxsiVx_[i]) + (*K_[i])[1] * (*VxsiVy_[i]);
    }

    // FFT back
    fft_ptr_->fftn(VxsiVx_.data(), VxsiVx_.data(), nn_, -1);

    // Evolve Matrix
#pragma omp parallel for
    for (int i = 0; i < Ncells_; i++) {
        Matrix left = negI_dssqrt * (*VxsiVx_[i]);
        Matrix right(0.);

        for (int a = 0; a < Nc2m1; a++) {
            const Matrix &Ta = group_ptr_->getT(a);
            const double c = real(CKxi_[i][a]);
            for (int idx = 0; idx < Ta.getNN(); idx++) {
                right.set(idx, right(idx) + c * Ta(idx));
            }
        }
        right *= posI_dssqrt;
        wilsonLines[i] = left.expm() * wilsonLines[i] * right.expm();
    }
}

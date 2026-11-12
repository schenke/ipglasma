// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#ifndef SRC_DEVICE_DEVICEDRIVERS_H_
#define SRC_DEVICE_DEVICEDRIVERS_H_

#include <vector>

#include "../Lattice.h"
#include "../Matrix.h"
#include "DeviceBackend.h"

namespace ipgdev {

/**
 * One nucleus' JIMWLK evolution on a device: the steps of
 * JIMWLK::evolutionStep() with the kernels of ipg_kernels.h.
 *
 * The momentum-space kernel K, the generators and the coupling of the
 * noise are computed on the host, by the code the CPU evolution uses, and
 * copied to the device once. The Gaussian noise of every step is drawn on
 * the host as well, from the random numbers of the run, so that a device
 * run consumes them exactly like a CPU run and the two can be compared
 * site by site.
 */
class JimwlkDevice {
  public:
    /**
     * \param[in] backend The device.
     * \param[in] n Lattice side.
     * \param[in] generators The 8 generators, 18 reals each.
     * \param[in] kernel K in momentum space, [site][x, y] complex (4 reals
     * per site), as JIMWLK::initializeK() leaves it.
     * \param[in] sqrtAlphasNoise sqrt(alpha_s(k)) per site for the noise
     * running coupling (`jimwlkRunningCoupling 2`), or nullptr.
     * \param[in] ds The step `jimwlkDs`.
     */
    JimwlkDevice(
        Backend &backend, int n, const double *generators, const double *kernel,
        const double *sqrtAlphasNoise, double ds);
    ~JimwlkDevice();
    JimwlkDevice(const JimwlkDevice &) = delete;
    JimwlkDevice &operator=(const JimwlkDevice &) = delete;

    /// Copies the Wilson lines of a nucleus to the device.
    void upload(const std::vector<Matrix> &wilsonLines);
    /// Copies them back.
    void download(std::vector<Matrix> &wilsonLines);
    /// One step with the noise \p noise, [site][2 x 8] real numbers in the
    /// order of the CPU code (x components first).
    void step(const double *noise);
    /**
     * One step with noise drawn on the device (the kernel
     * ipg_k_gauss_noise, a counter-based generator): independent Gaussian
     * numbers of unit variance for every site, direction and color, as in
     * step(), without the copy from the host.
     * \param[in] key0, key1 The key of the generator; the caller draws it
     * from the random numbers of the run, once per nucleus.
     * \param[in] stepNumber The number of the step, the second word of
     * the counter: every step of a key needs another one.
     */
    void stepWithDeviceNoise(
        unsigned int key0, unsigned int key1, unsigned int stepNumber);
    /// The noise of the last step, [site][16] (for tests).
    void downloadNoise(std::vector<double> &noise);

  private:
    Backend &backend_;
    const int n_;
    const std::size_t sites_;
    const double sqrtDs_;
    const bool noiseCoupling_;
    Buffer wilson_, generators_, kernel_, alphas_;
    Buffer noise_, xi_, eta_, ck_, v_, w_;

    /// The step once the noise is in noise_.
    void evolveWithNoise();
};

/**
 * The classical Yang-Mills leapfrog on a device: the updates of
 * Evolution.cpp (evolvePi, evolveE, evolvePhi, evolveU) with the kernels
 * of ipg_kernels.h. The observables stay on the host: the fields are
 * copied back where Evolution::run() measures.
 */
class YangMillsDevice {
  public:
    YangMillsDevice(Backend &backend, int n, bool periodic, double g);
    ~YangMillsDevice();
    YangMillsDevice(const YangMillsDevice &) = delete;
    YangMillsDevice &operator=(const YangMillsDevice &) = delete;

    /// Copies E1, E2 (Lattice::U, U2), the links (Ux, Uy), pi (Ux2) and
    /// phi (Uy2) to the device.
    void upload(const Lattice &lat);
    /// Copies all six fields back.
    void download(Lattice &lat);
    /// Copies only E1, E2 and pi back.
    void downloadMomenta(Lattice &lat);

    void evolvePhi(double dtau, double tau);
    void evolveU(double dtau, double tau);
    void evolvePi(double dtau, double tau);
    void evolveE(double dtau, double tau);
    /// evolveStepPersistent() of Evolution.cpp.
    void step(double dtau, double tau, bool updateCoordinates);

    /// Keeps a copy of E1, E2 and pi on the device ...
    void backupMomenta();
    /// ... and restores them from it.
    void restoreMomenta();

  private:
    Backend &backend_;
    const int n_;
    const std::size_t sites_;
    const int periodic_;
    const double g_;
    Buffer e1_, e2_, ux_, uy_, pi_, phi_;
    Buffer e1Backup_, e2Backup_, piBackup_;
};

}  // namespace ipgdev

#endif  // SRC_DEVICE_DEVICEDRIVERS_H_

#include "JIMWLK.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <vector>

#include "Instrumentation.h"

using PhysConst::invHbarc;
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
    random_ptr_ = random;
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

    double mu0 = param_.jimwlk.mu0;
    double Lambda2 = param_.jimwlk.LambdaQCD * param_.jimwlk.LambdaQCD;

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
        double alphas_sqroot = sqrt(getAlphas(x, y));
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
    if (param_.jimwlk.alphaS > 1e-10) {
        return 1.0;
    }

    const double c = param_.jimwlk.c;
    const int Nf = param_.coupling.nFlavors;
    const double length = param_.lattice.L;
    const double mu0 = param_.jimwlk.mu0;
    const double Lambda2 = param_.jimwlk.LambdaQCD * param_.jimwlk.LambdaQCD;
    double phys_x = x * length;  // in fm
    double phys_y = y * length;
    double phys_r2 = phys_x * phys_x + phys_y * phys_y;

    // Alphas in physical units! Lambda2 is lambda_QCD^2 in GeV
    double alphas =
        4. * M_PI
        / ((11. * Nc - 2. * Nf) / 3. * c
           * log((
               pow(mu0 * mu0 / Lambda2, 1. / c)
               + pow(4. / (phys_r2 * Lambda2 * invHbarc * invHbarc), 1. / c))));
    return alphas;
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
    initializeNoise();

    // Calculate evolution steps for different nuclei
    double x0 = param_.jimwlk.initialX;
    double ds = param_.jimwlk.Ds;
    bool saveSnapshots = param_.jimwlk.saveSnapshots;
    std::vector<double> xSnapshotList = param_.jimwlk.xSnapshotList;
    double dlogx = M_PI * M_PI * ds;
    int steps_1 = 0;
    int steps_2 = 0;
    double as = param_.jimwlk.alphaS;
    if (as > 1e-10) {
        // Fixed coupling
        steps_1 = static_cast<int>(
            as * std::log(x0 / param_.jimwlk.xProjectile) / (M_PI * M_PI * ds)
            + 0.5);
        steps_2 = static_cast<int>(
            as * std::log(x0 / param_.jimwlk.xTarget) / (M_PI * M_PI * ds)
            + 0.5);
        dlogx = M_PI * M_PI * ds / as;
    } else {
        // Running coupling
        steps_1 = static_cast<int>(
            std::log(x0 / param_.jimwlk.xProjectile) / (M_PI * M_PI * ds)
            + 0.5);
        steps_2 = static_cast<int>(
            std::log(x0 / param_.jimwlk.xTarget) / (M_PI * M_PI * ds) + 0.5);
    }

    runEvolutionLoop(
        NucleusRole::Projectile, steps_1, x0, dlogx, saveSnapshots,
        xSnapshotList);
    runEvolutionLoop(
        NucleusRole::Target, steps_2, x0, dlogx, saveSnapshots, xSnapshotList);
}

void JIMWLK::runEvolutionLoop(
    NucleusRole nucleus, int steps, double x0, double dlogx, bool saveSnapshots,
    const std::vector<double> &xSnapshotList) {
    const std::string label =
        (nucleus == NucleusRole::Projectile) ? "projectile" : "target";
    messager_ << "[JIMWLK::evolution]: Evolving " << label
              << ", evolution steps " << steps;
    messager_.flush("info");
    unsigned int iSnapshot = 0;
    for (int ids = 0; ids < steps; ids++) {
        // steps (from user-configurable JIMWLK parameters, no lower bound
        // enforced) can be under 10, making this 0; guard against the
        // resulting integer modulo-by-zero below.
        int printSteps = std::max(1, steps / 10);
        if (ids % printSteps == 0) {
            messager_ << "[JIMWLK::evolution]: Step " << ids;
            messager_.flush("info");
        }
        double xLoc = x0 * exp(-ids * dlogx);
        evolutionStep(nucleus);
        if (saveSnapshots) {
            if (iSnapshot < xSnapshotList.size()) {
                if (xLoc > xSnapshotList[iSnapshot]
                    && xLoc * exp(-dlogx) < xSnapshotList[iSnapshot]) {
                    lat_ptr_->writeWilsonLines(&param_, nucleus, xLoc);
                    iSnapshot++;
                }
            }
        }
    }
    messager_ << "[JIMWLK::evolution]: Done.";
    messager_.flush("info");
}

void JIMWLK::evolutionStep(NucleusRole nucleus) {
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

    // the local xi now contains the Fourier transform of xi,
    // while the original xi is stored in the array xi2
    fft_ptr_->fftnArray(xi2_.data(), xi_.data(), nn_, 1, 2 * Nc2m1);

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

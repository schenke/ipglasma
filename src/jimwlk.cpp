#include "jimwlk.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <vector>

JIMWLK::JIMWLK(Parameters &param, Group *group, Lattice *lat, Random *random)
    : param_(param),
      Nc_(param.getNc()),
      Nc2m1_(param.getNc() * param.getNc() - 1),
      Ngrid_(param.getSize()),
      Ncells_(param.getSize() * param.getSize()) {
    nn_[0] = param_.getSize();
    nn_[1] = param_.getSize();

    fft_ptr_ = std::make_shared<FFT>(nn_);

    group_ptr_ = group;
    random_ptr_ = random;
    lat_ptr_ = lat;

    initializeK();
    initializeNoise();
    if (param_.getSimpleLangevin()) {
        VxsiVx_ = new Matrix *[Ncells_];
        VxsiVy_ = new Matrix *[Ncells_];
        for (int i = 0; i < Ncells_; i++) {
            VxsiVx_[i] = new Matrix(Nc_, 0);
            VxsiVy_[i] = new Matrix(Nc_, 0);
        }
    }
}

JIMWLK::~JIMWLK() {
    if (initializedK_) {
        for (int i = 0; i < Ncells_; i++) {
            delete K_[i];
        }
        delete[] K_;
    }

    if (initializedNoise_) {
        delete[] xi_data_;
        delete[] xi2_data_;
        delete[] CKxi_data_;
        delete[] xi_;
        delete[] xi2_;
        delete[] CKxi_;
    }

    if (param_.getSimpleLangevin()) {
        for (int i = 0; i < Ncells_; i++) {
            delete VxsiVx_[i];
            delete VxsiVy_[i];
        }
        delete[] VxsiVx_;
        delete[] VxsiVy_;
    }
}

void JIMWLK::initializeK() {
    if (initializedK_) {
        return;
    }
    K_ = new std::vector<std::complex<double> > *[Ncells_];
    for (int i = 0; i < Ncells_; i++) {
        // sized to its final length up front (always x,y) so it never
        // has to grow/reallocate via push_back; entries default to (0,0)
        K_[i] = new std::vector<std::complex<double> >(2);
    }

    double mu0 = param_.getMu0_jimwlk();
    double Lambda2 =
        param_.getLambdaQCD_jimwlk() * param_.getLambdaQCD_jimwlk();

    // set once here rather than on every getMassRegulator() call below
    gsl_set_error_handler_off();

#pragma omp parallel for
    for (int pos = 0; pos < Ncells_; pos++) {
        double x = pos / Ngrid_ - static_cast<double>(Ngrid_) / 2.;
        double y = pos % Ngrid_ - static_cast<double>(Ngrid_) / 2.;
        x /= Ngrid_;
        y /= Ngrid_;
        double r2 = x * x + y * y;
        if (r2 < 1e-16) {
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
    fft_ptr_->fftnVector(K_, K_, nn_, 1);
    initializedK_ = true;
}

double JIMWLK::getMassRegulator(const double x, const double y) const {
    // if m suppresses long distance tails,
    // K is multiplied by this, which is m*r*K_1(m*r)
    double mass_regulator = 1.0;
    double m = param_.getm_jimwlk();
    if (m < 1e-16) {
        return mass_regulator;
    }
    double length = param_.getL();

    // Lattice units
    // Here x is [-N/2, N/2]
    double lat_x = sin(M_PI * x) / (M_PI);
    double lat_y = sin(M_PI * y) / (M_PI);
    // lat_x and lat_y are now in [-1/2,1/2] as x/nn[0] is in [-N/2, N/2]
    double lat_r = sqrt(lat_x * lat_x + lat_y * lat_y) * Ngrid_;
    // lat_r now tells how many lattice units the distance is
    double a = length / Ngrid_;
    double lat_m = m * a * fmgev;
    double bessel_argument = lat_m * lat_r;
    // double bes = std::cyl_bessel_k(1, bessel_argument);
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
    double alphas = 1.0;
    if (param_.getJimwlk_alphas() > 1e-10) {
        return alphas;
    }

    const double c = 0.2;
    const int Nf = 3;
    const double length = param_.getL();
    const double mu0 = param_.getMu0_jimwlk();
    const double Lambda2 =
        param_.getLambdaQCD_jimwlk() * param_.getLambdaQCD_jimwlk();
    double phys_x = x * length;  // in fm
    double phys_y = y * length;
    double phys_r2 = phys_x * phys_x + phys_y * phys_y;

    // Alphas in physical units! Lambda2 is lambda_QCD^2 in GeV
    alphas = 4. * M_PI
             / ((11. * Nc_ - 2. * Nf) / 3. * c
                * log(
                    (pow(mu0 * mu0 / Lambda2, 1. / c)
                     + pow(4. / (phys_r2 * Lambda2 * fmgev * fmgev), 1. / c))));
    return alphas;
}

void JIMWLK::initializeNoise() {
    if (initializedNoise_) {
        return;
    }
    // one contiguous allocation per array instead of Ncells_ separate
    // ones, so consecutive cells are adjacent in memory
    xi_data_ = new std::complex<double>[Ncells_ * 2 * Nc2m1_];
    xi2_data_ = new std::complex<double>[Ncells_ * 2 * Nc2m1_];
    CKxi_data_ = new std::complex<double>[Ncells_ * Nc2m1_];

    xi_ = new std::complex<double> *[Ncells_];
    xi2_ = new std::complex<double> *[Ncells_];
    CKxi_ = new std::complex<double> *[Ncells_];
    for (int i = 0; i < Ncells_; i++) {
        xi_[i] = xi_data_ + i * 2 * Nc2m1_;
        xi2_[i] = xi2_data_ + i * 2 * Nc2m1_;
        CKxi_[i] = CKxi_data_ + i * Nc2m1_;
    }
    initializedNoise_ = true;
}

void JIMWLK::evolution() {
    initializeNoise();

    // Calculate evolution steps for different nuclei
    double x0 = param_.getJimwlk_x0();
    double ds = param_.getDs_jimwlk();
    bool saveSnapshots = param_.getSaveSnapshots();
    std::vector<double> xSnapshotList = param_.getxSnapshotList();
    double dlogx = M_PI * M_PI * ds;
    int steps_1 = 0;
    int steps_2 = 0;
    double as = param_.getJimwlk_alphas();
    if (as > 1e-10) {
        // Fixed coupling
        steps_1 = static_cast<int>(
            as * std::log(x0 / param_.getJimwlk_x_projectile())
                / (M_PI * M_PI * ds)
            + 0.5);
        steps_2 = static_cast<int>(
            as * std::log(x0 / param_.getJimwlk_x_target()) / (M_PI * M_PI * ds)
            + 0.5);
        dlogx = M_PI * M_PI * ds / as;
    } else {
        // Running coupling
        steps_1 = static_cast<int>(
            std::log(x0 / param_.getJimwlk_x_projectile()) / (M_PI * M_PI * ds)
            + 0.5);
        steps_2 = static_cast<int>(
            std::log(x0 / param_.getJimwlk_x_target()) / (M_PI * M_PI * ds)
            + 0.5);
    }

    unsigned int iSnapshot = 0;
    std::cout << "Evolving projectile, evolution steps " << steps_1
              << std::endl;
    for (int ids = 0; ids < steps_1; ids++) {
        int printSteps = steps_1 / 10;
        if (ids % printSteps == 0) {
            std::cout << "Step " << ids << std::endl;
        }
        double xLoc = x0 * exp(-ids * dlogx);
        evolutionStep(NucleusRole::Projectile);
        if (saveSnapshots) {
            if (iSnapshot < xSnapshotList.size()) {
                if (xLoc > xSnapshotList[iSnapshot]
                    && xLoc * exp(-dlogx) < xSnapshotList[iSnapshot]) {
                    std::stringstream ss;
                    ss << "JIMWLKSnapshot_x_" << xLoc << "_";
                    lat_ptr_->WriteWilsonLines(ss.str(), &param_, 1);
                    iSnapshot++;
                }
            }
        }
    }
    std::cout << "Done." << std::endl;

    std::cout << "Evolving target, evolution steps " << steps_2 << std::endl;
    iSnapshot = 0;
    for (int ids = 0; ids < steps_2; ids++) {
        int printSteps = steps_2 / 10;
        if (ids % printSteps == 0) {
            std::cout << "Step " << ids << std::endl;
        }
        double xLoc = x0 * exp(-ids * dlogx);
        evolutionStep(NucleusRole::Target);
        if (saveSnapshots) {
            if (iSnapshot < xSnapshotList.size()) {
                if (xLoc > xSnapshotList[iSnapshot]
                    && xLoc * exp(-dlogx) < xSnapshotList[iSnapshot]) {
                    std::stringstream ss;
                    ss << "JIMWLKSnapshot_x_" << xLoc << "_";
                    lat_ptr_->WriteWilsonLines(ss.str(), &param_, 2);
                    iSnapshot++;
                }
            }
        }
    }
    std::cout << "Done." << std::endl;
}

void JIMWLK::evolutionStep(NucleusRole nucleus) {
    const bool evolveProjectile = (nucleus == NucleusRole::Projectile);
    std::vector<Matrix> &wilsonLines =
        evolveProjectile ? lat_ptr_->U : lat_ptr_->U2;

    const complex<double> I(0., 1.);
    const double ds_sqrt = std::sqrt(param_.getDs_jimwlk());
    const complex<double> negI_dssqrt = -I * ds_sqrt;
    const complex<double> posI_dssqrt = I * ds_sqrt;

    // generate random Gaussian noise in every cell for Nc^2-1 color
    // components and 2 spatial components x and y
    // (kept serial: random_ptr_->Gauss() mutates shared RNG state)
    for (int i = 0; i < Ncells_; i++) {
        for (int n = 0; n < 2 * Nc2m1_; n++) {
            xi2_[i][n] = std::complex<double>(random_ptr_->Gauss(), 0.);
        }
    }

    // the local xi now contains the Fourier transform of xi,
    // while the original xi is stored in the array xi2
    fft_ptr_->fftnArray(xi2_, xi_, nn_, 1, 2 * Nc2m1_);

    // now compute C(K_i,xi_i^a) = F^{-1}(F(K_i)F(xi_i^a))
    //                           = F^{-1}(F(K_x)F(xi_x^a)+F(K_y)F(xi_y^a))
#pragma omp parallel for
    for (int i = 0; i < Ncells_; i++) {
        for (int n = 0; n < Nc2m1_; n++) {
            CKxi_[i][n] =
                (*K_[i])[0] * xi_[i][n] + (*K_[i])[1] * xi_[i][n + Nc2m1_];
            // product of x components + product of y components
        }
    }

    // now CKxi contains C(K_i,xi_i^a) - it is a vector with a components
    fft_ptr_->fftnArray(CKxi_, CKxi_, nn_, -1, Nc2m1_);

#pragma omp parallel for
    for (int i = 0; i < Ncells_; i++) {
        *VxsiVx_[i] = zero_;
        *VxsiVy_[i] = zero_;
        const Matrix &U = wilsonLines[i];
        for (int a = 0; a < Nc2m1_; a++) {
            Matrix UTa = U * group_ptr_->getT(a);
            // UTa * U^dagger, folding in the conjugate transpose of U
            // instead of building it explicitly (Nc=3 specialization)
            Matrix UTUconj = UTa.prodABconj(UTa, U);
            *VxsiVx_[i] += xi2_[i][a] * UTUconj;
            *VxsiVy_[i] += xi2_[i][a + Nc2m1_] * UTUconj;
        }
    }

    // FFT V xi V
    fft_ptr_->fftn(VxsiVx_, VxsiVx_, nn_, 1);
    fft_ptr_->fftn(VxsiVy_, VxsiVy_, nn_, 1);

#pragma omp parallel for
    for (int i = 0; i < Ncells_; i++) {
        *VxsiVx_[i] = (*K_[i])[0] * (*VxsiVx_[i]) + (*K_[i])[1] * (*VxsiVy_[i]);
    }

    // FFT back
    fft_ptr_->fftn(VxsiVx_, VxsiVx_, nn_, -1);

    // Evolve Matrix
#pragma omp parallel for
    for (int i = 0; i < Ncells_; i++) {
        Matrix left = negI_dssqrt * (*VxsiVx_[i]);
        Matrix right(Nc_, 0.);

        for (int a = 0; a < Nc2m1_; a++) {
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

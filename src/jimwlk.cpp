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
    evolution();
}

JIMWLK::~JIMWLK() {
    if (initializedK_) {
        for (int i = 0; i < Ncells_; i++) {
            delete K_[i];
        }
        delete[] K_;
    }

    if (initializedNoise_) {
        for (int i = 0; i < Ncells_; i++) {
            delete xi_[i];
            delete xi2_[i];
            delete CKxi_[i];
        }
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
        K_[i] = new std::vector<std::complex<double> >;
    }

    double mu0 = param_.getMu0_jimwlk();
    double Lambda2 = param_.getLambdaQCD_jimwlk() * param_.getLambdaQCD_jimwlk();

    for (int pos = 0; pos < Ncells_; pos++) {
        double x = pos / Ngrid_ - static_cast<double>(Ngrid_) / 2.;
        double y = pos % Ngrid_ - static_cast<double>(Ngrid_) / 2.;
        x /= Ngrid_;
        y /= Ngrid_;
        double r2 = x * x + y * y;
        if (r2 < 1e-16) {
            K_[pos]->push_back(0.);
            K_[pos]->push_back(0.);
            continue;
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
        tmpk1 *= factor;
        tmpk2 *= factor;

        K_[pos]->push_back(tmpk1);
        K_[pos]->push_back(tmpk2);
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
    gsl_sf_result bes;
    gsl_set_error_handler_off();  // so that when we have too large
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
    if (param_.getRunningCoupling_jimwlk() == 0) {
        return alphas;
    }

    const double c = 0.2;
    const int Nf = 3;
    const double length = param_.getL();
    const double mu0 = param_.getMu0_jimwlk();
    const double Lambda2 = param_.getLambdaQCD_jimwlk() * param_.getLambdaQCD_jimwlk();
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
    xi_ = new std::complex<double> *[Ncells_];
    xi2_ = new std::complex<double> *[Ncells_];
    CKxi_ = new std::complex<double> *[Ncells_];
    for (int i = 0; i < Ncells_; i++) {
        xi_[i] = new std::complex<double>[2 * Nc2m1_];
        xi2_[i] = new std::complex<double>[2 * Nc2m1_];
        CKxi_[i] = new std::complex<double>[Nc2m1_];
    }
    initializedNoise_ = true;
}

void JIMWLK::evolution() {
    initializeNoise();
    const int steps = param_.getSteps_jimwlk();
    std::cout << "Beginning evolution ..." << std::endl;
    for (int ids = 0; ids < steps; ids++) {
        if (ids % 100 == 0) {
            std::cout << "Step " << ids << std::endl;
        }
        evolutionStep();
        evolutionStep2();
    }
    std::cout << "Done." << std::endl;
}

void JIMWLK::evolutionStep() {
    const complex<double> I(0., 1.);
    const double ds_sqrt = std::sqrt(param_.getDs_jimwlk());

    // generate random Gaussian noise in every cell for Nc^2-1 color
    // components and 2 spatial components x and y
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
    for (int i = 0; i < Ncells_; i++) {
        for (int n = 0; n < Nc2m1_; n++) {
            CKxi_[i][n] = (*K_[i]).at(0) * xi_[i][n]
                          + (*K_[i]).at(1) * xi_[i][n + Nc2m1_];
            // product of x components + product of y components
        }
    }

    // now CKxi contains C(K_i,xi_i^a) - it is a vector with a components
    fft_ptr_->fftnArray(CKxi_, CKxi_, nn_, -1, Nc2m1_);

    for (int i = 0; i < Ncells_; i++) {
        *VxsiVx_[i] = zero_;
        *VxsiVy_[i] = zero_;
        for (int a = 0; a < Nc2m1_; a++) {
            Matrix Uconj = lat_ptr_->cells[i]->getU();
            Uconj.conjg();
            *VxsiVx_[i] = (*VxsiVx_[i])
                          + xi2_[i][a] * lat_ptr_->cells[i]->getU()
                                * group_ptr_->getT(a) * Uconj;
            *VxsiVy_[i] = (*VxsiVy_[i])
                          + xi2_[i][a + Nc2m1_] * lat_ptr_->cells[i]->getU()
                                * group_ptr_->getT(a) * Uconj;
        }
    }

    // FFT V xi V
    fft_ptr_->fftn(VxsiVx_, VxsiVx_, nn_, 1);
    fft_ptr_->fftn(VxsiVy_, VxsiVy_, nn_, 1);

    for (int i = 0; i < Ncells_; i++) {
        *VxsiVx_[i] = (*K_[i])[0] * (*VxsiVx_[i]) + (*K_[i])[1] * (*VxsiVy_[i]);
    }

    // FFT back
    fft_ptr_->fftn(VxsiVx_, VxsiVx_, nn_, -1);

    // Evolve Matrix
    for (int i = 0; i < Ncells_; i++) {
        Matrix left(Nc_, 0.);
        left = -I * ds_sqrt * (*VxsiVx_[i]);
        Matrix right(Nc_, 0.);

        for (int a = 0; a < Nc2m1_; a++) {
            right = right + real(CKxi_[i][a]) * group_ptr_->getT(a);
        }
        right = I * ds_sqrt * right;
        lat_ptr_->cells[i]->setU(
            left.expm() * lat_ptr_->cells[i]->getU() * right.expm());
    }
}

void JIMWLK::evolutionStep2() {
    const complex<double> I(0., 1.);
    const double ds_sqrt = std::sqrt(param_.getDs_jimwlk());

    // generate random Gaussian noise in every cell for Nc^2-1 color
    // components and 2 spatial components x and y
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
    for (int i = 0; i < Ncells_; i++) {
        for (int n = 0; n < Nc2m1_; n++) {
            CKxi_[i][n] = (*K_[i]).at(0) * xi_[i][n]
                          + (*K_[i]).at(1) * xi_[i][n + Nc2m1_];
            // product of x components + product of y components
        }
    }

    // now CKxi contains C(K_i,xi_i^a) - it is a vector with a components
    fft_ptr_->fftnArray(CKxi_, CKxi_, nn_, -1, Nc2m1_);

    for (int i = 0; i < Ncells_; i++) {
        *VxsiVx_[i] = zero_;
        *VxsiVy_[i] = zero_;
        for (int a = 0; a < Nc2m1_; a++) {
            Matrix Uconj = lat_ptr_->cells[i]->getU2();
            Uconj.conjg();
            *VxsiVx_[i] = (*VxsiVx_[i])
                          + xi2_[i][a] * lat_ptr_->cells[i]->getU2()
                                * group_ptr_->getT(a) * Uconj;
            *VxsiVy_[i] = (*VxsiVy_[i])
                          + xi2_[i][a + Nc2m1_] * lat_ptr_->cells[i]->getU2()
                                * group_ptr_->getT(a) * Uconj;
        }
    }

    // FFT V xi V
    fft_ptr_->fftn(VxsiVx_, VxsiVx_, nn_, 1);
    fft_ptr_->fftn(VxsiVy_, VxsiVy_, nn_, 1);

    for (int i = 0; i < Ncells_; i++) {
        *VxsiVx_[i] = (*K_[i])[0] * (*VxsiVx_[i]) + (*K_[i])[1] * (*VxsiVy_[i]);
    }

    // FFT back
    fft_ptr_->fftn(VxsiVx_, VxsiVx_, nn_, -1);

    // Evolve Matrix
    for (int i = 0; i < Ncells_; i++) {
        Matrix left(Nc_, 0.);
        left = -I * ds_sqrt * (*VxsiVx_[i]);
        Matrix right(Nc_, 0.);

        for (int a = 0; a < Nc2m1_; a++) {
            right = right + real(CKxi_[i][a]) * group_ptr_->getT(a);
        }
        right = I * ds_sqrt * right;
        lat_ptr_->cells[i]->setU2(
            left.expm() * lat_ptr_->cells[i]->getU2() * right.expm());
    }
}


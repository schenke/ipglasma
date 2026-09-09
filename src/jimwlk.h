#ifndef SRC_JIMWLK_H_
#define SRC_JIMWLK_H_

#include <complex>
#include <memory>
#include <vector>

#include "FFT.h"
#include "Group.h"
#include "Lattice.h"
#include "Parameters.h"
#include "Random.h"

enum class NucleusRole { Projectile, Target };

class JIMWLK {
  private:
    Parameters &param_;
    std::shared_ptr<FFT> fft_ptr_;
    int nn_[2];

    const double fmgev = 5.068;

    const int Nc_;
    const int Nc2m1_;
    const int Ngrid_;
    const int Ncells_;

    // could change to smart pointers later
    Group *group_ptr_;
    Random *random_ptr_;
    Lattice *lat_ptr_;

    bool initializedK_ = false;
    bool initializedNoise_ = false;
    std::vector<std::complex<double> > **K_;  // data type matches FFT.h

    // xi_/xi2_/CKxi_ are pointer-per-cell views into one contiguous
    // backing buffer each (xi_data_/xi2_data_/CKxi_data_), so the hot
    // per-cell loops in evolutionStep() get cache-friendly, sequential
    // access instead of chasing Ncells_ independent heap allocations.
    std::complex<double> **xi_;    // noise
    std::complex<double> **xi2_;   // noise
    std::complex<double> **CKxi_;  // noise
    std::complex<double> *xi_data_;
    std::complex<double> *xi2_data_;
    std::complex<double> *CKxi_data_;

    Matrix **VxsiVx_;
    Matrix **VxsiVy_;
    Matrix zero_ = Matrix(Nc_, 0);

  public:
    JIMWLK() = delete;
    JIMWLK(Parameters &param, Group *group, Lattice *lat, Random *random);
    ~JIMWLK();

    void initializeK();
    double getMassRegulator(const double x, const double y) const;
    double getAlphas(const double x, const double y) const;
    void initializeNoise();

    void evolution();
    void evolutionStep(NucleusRole nucleus);
};

#endif  // SRC_JIMWLK_H_

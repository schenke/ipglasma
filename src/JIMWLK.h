#ifndef SRC_JIMWLK_H_
#define SRC_JIMWLK_H_

#include <complex>
#include <memory>
#include <vector>

#include "FFT.h"
#include "Glauber.h"
#include "Group.h"
#include "Lattice.h"
#include "Parameters.h"
#include "PrettyOstream.h"
#include "Random.h"

class JIMWLK {
  private:
    Parameters &param_;
    std::shared_ptr<FFT> fft_ptr_;
    int nn_[2];
    PrettyOstream messager_;

    const double fmgev = 5.068;

    static constexpr int Nc_ = 3;
    static constexpr int Nc2m1_ = Nc_ * Nc_ - 1;
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
    Matrix zero_ = Matrix(0.);

    // Persistent scratch for the per-step bulk noise draw in
    // evolutionStep(), reused across steps to avoid reallocating.
    std::vector<double> gaussNoise_;
    std::vector<double> gaussNoiseScratch_;

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

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

    // K_storage_ owns the per-cell 2-vectors; K_ is a pointer-per-cell view
    // over it (data type matches FFT.h's T** interface).
    std::vector<std::vector<std::complex<double> > > K_storage_;
    std::vector<std::vector<std::complex<double> > *> K_;

    // xi_/xi2_/CKxi_ are pointer-per-cell views into one contiguous
    // backing buffer each (xi_data_/xi2_data_/CKxi_data_), so the hot
    // per-cell loops in evolutionStep() get cache-friendly, sequential
    // access instead of chasing Ncells_ independent heap allocations.
    std::vector<std::complex<double> *> xi_;    // noise
    std::vector<std::complex<double> *> xi2_;   // noise
    std::vector<std::complex<double> *> CKxi_;  // noise
    std::vector<std::complex<double> > xi_data_;
    std::vector<std::complex<double> > xi2_data_;
    std::vector<std::complex<double> > CKxi_data_;

    // VxsiVx_/VxsiVy_ likewise own their Ncells_ Matrix storage via
    // *Storage_ and expose a pointer-per-cell view for FFT::fftn's T**
    // interface.
    std::vector<Matrix> VxsiVxStorage_;
    std::vector<Matrix> VxsiVyStorage_;
    std::vector<Matrix *> VxsiVx_;
    std::vector<Matrix *> VxsiVy_;
    Matrix zero_ = Matrix(0.);

    // Persistent scratch for the per-step bulk noise draw in
    // evolutionStep(), reused across steps to avoid reallocating.
    std::vector<double> gaussNoise_;
    std::vector<double> gaussNoiseScratch_;

  public:
    JIMWLK() = delete;
    JIMWLK(Parameters &param, Group *group, Lattice *lat, Random *random);
    ~JIMWLK() = default;

    // K_/xi_/xi2_/CKxi_/VxsiVx_/VxsiVy_ are pointer-per-cell views aliasing
    // this object's own *_storage_/*_data_ vectors (see their declarations
    // below). A compiler-generated copy would deep-copy the storage but
    // shallow-copy the views, leaving the copy's views pointing into the
    // original's memory instead of its own -- so disable copying (nothing
    // needs it; JIMWLK is only ever stack-constructed once in main.cpp).
    JIMWLK(const JIMWLK &) = delete;
    JIMWLK &operator=(const JIMWLK &) = delete;

    void initializeK();
    double getMassRegulator(const double x, const double y) const;
    double getAlphas(const double x, const double y) const;
    void initializeNoise();

    void evolution();
    void evolutionStep(NucleusRole nucleus);

  private:
    // Shared by evolution()'s projectile and target passes: runs `steps`
    // Langevin steps for `nucleus`, logging progress and writing any
    // snapshots due in xSnapshotList along the way.
    void runEvolutionLoop(
        NucleusRole nucleus, int steps, double x0, double dlogx,
        bool saveSnapshots, const std::vector<double> &xSnapshotList);
};

#endif  // SRC_JIMWLK_H_

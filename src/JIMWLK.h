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
#include "PhysConst.h"
#include "PrettyOstream.h"
#include "Random.h"

/**
 * JIMWLK small-x evolution of the projectile's and target's Wilson
 * lines, run before the classical Yang-Mills evolution.
 *
 * Each evolutionStep() advances one nucleus' Wilson-line field by one
 * step of the JIMWLK Langevin equation: real Gaussian noise is drawn
 * per site/color, convolved (via FFT) with a regularized momentum-space
 * kernel \c K_ built once by initializeK(), gauge-covariantly rotated
 * by the current Wilson line, convolved with \c K_ a second time, and
 * exponentiated back onto the Wilson line. evolution() repeats this for
 * as many steps as the requested rapidity range requires, for both
 * nuclei in turn.
 */
class JIMWLK {
  private:
    /// Simulation parameters.
    Parameters &param_;
    /// FFT instance sized for this lattice, used for every noise/kernel
    /// transform.
    std::shared_ptr<FFT> fft_ptr_;
    /// Lattice dimensions `{Ngrid_, Ngrid_}`, in FFT::fftn()'s `nn[]`
    /// format.
    int nn_[2];
    /// Log sink for progress messages.
    PrettyOstream messager_;

    /// Number of colors; fixed at 3.
    static constexpr int Nc_ = 3;
    /// SU(3) adjoint dimension, \f$N_c^2-1=8\f$.
    static constexpr int Nc2m1_ = Nc_ * Nc_ - 1;
    /// Lattice side length.
    const int Ngrid_;
    /// Total number of lattice sites, `Ngrid_*Ngrid_`.
    const int Ncells_;

    // could change to smart pointers later
    /// Non-owning pointer to the shared Group instance (the SU(3)
    /// generators).
    Group *group_ptr_;
    /// Non-owning pointer to the shared Random instance.
    Random *random_ptr_;
    /// Non-owning pointer to the Lattice whose Wilson lines
    /// (`lat_ptr_->U`/`U2`) this instance evolves.
    Lattice *lat_ptr_;

    /// Whether initializeK() has already built \c K_.
    bool initializedK_ = false;
    /// Whether initializeNoise() has already allocated the noise
    /// buffers.
    bool initializedNoise_ = false;

    // K_storage_ owns the per-cell 2-vectors; K_ is a pointer-per-cell view
    // over it (data type matches FFT.h's T** interface).
    /// Owns the per-cell 2-vectors \c K_ points into.
    std::vector<std::vector<std::complex<double> > > K_storage_;
    /// Regularized momentum-space evolution kernel, one 2-component
    /// (\f$x\f$, \f$y\f$) vector per site, built once by initializeK()
    /// and reused by every evolutionStep() call.
    std::vector<std::vector<std::complex<double> > *> K_;

    // xi_/xi2_/CKxi_ are pointer-per-cell views into one contiguous
    // backing buffer each (xi_data_/xi2_data_/CKxi_data_), so the hot
    // per-cell loops in evolutionStep() get cache-friendly, sequential
    // access instead of chasing Ncells_ independent heap allocations.
    /// Momentum-space Gaussian noise, one \f$2(N_c^2-1)\f$-component
    /// array (\f$x\f$/\f$y\f$ spatial components times 8 colors) per
    /// site, filled by evolutionStep()'s forward FFT of \c xi2_.
    std::vector<std::complex<double> *> xi_;  // noise
    /// Position-space Gaussian noise, freshly drawn every
    /// evolutionStep() call.
    std::vector<std::complex<double> *> xi2_;  // noise
    /// The kernel-noise contraction \f$C(K,\xi)^a = K_x\xi_x^a +
    /// K_y\xi_y^a\f$, one 8-component (color) array per site.
    std::vector<std::complex<double> *> CKxi_;  // noise
    /// Backing storage for \c xi_'s per-cell views.
    std::vector<std::complex<double> > xi_data_;
    /// Backing storage for \c xi2_'s per-cell views.
    std::vector<std::complex<double> > xi2_data_;
    /// Backing storage for \c CKxi_'s per-cell views.
    std::vector<std::complex<double> > CKxi_data_;

    // VxsiVx_/VxsiVy_ likewise own their Ncells_ Matrix storage via
    // *Storage_ and expose a pointer-per-cell view for FFT::fftn's T**
    // interface.
    /// Backing storage for \c VxsiVx_'s per-cell views.
    std::vector<Matrix> VxsiVxStorage_;
    /// Backing storage for \c VxsiVy_'s per-cell views.
    std::vector<Matrix> VxsiVyStorage_;
    /// The gauge-covariant noise combination
    /// \f$\sum_a \xi_x^a\, U t^a U^\dagger\f$ (\f$x\f$ component),
    /// reused as scratch across evolutionStep()'s FFT convolution
    /// steps.
    std::vector<Matrix *> VxsiVx_;
    /// Same as \c VxsiVx_, \f$y\f$ component.
    std::vector<Matrix *> VxsiVy_;
    /// Reusable zero matrix, to reset \c VxsiVx_/\c VxsiVy_ at the
    /// start of each evolutionStep() without reallocating.
    Matrix zero_ = Matrix(0.);

    // Persistent scratch for the per-step bulk noise draw in
    // evolutionStep(), reused across steps to avoid reallocating.
    /// Flat buffer evolutionStep() draws this step's Gaussian noise
    /// into via Random::gaussBulk(), before scattering it into \c xi2_.
    std::vector<double> gaussNoise_;
    /// Scratch buffer Random::gaussBulk() reuses internally across
    /// calls.
    std::vector<double> gaussNoiseScratch_;

  public:
    JIMWLK() = delete;
    /**
     * Constructs a JIMWLK for a given lattice: allocates the FFT
     * instance, builds the momentum-space kernel and noise buffers, and
     * allocates the gauge-covariant-noise scratch fields (unconditionally,
     * since evolutionStep() needs them regardless of
     * `param.getSimpleLangevin()`).
     * \param[in] param Simulation parameters; `getSize()` sets the
     * lattice dimensions.
     * \param[in] group Non-owning pointer to the shared Group instance;
     * must outlive this JIMWLK.
     * \param[in] lat Non-owning pointer to the Lattice to evolve; must
     * outlive this JIMWLK.
     * \param[in] random Non-owning pointer to the shared Random
     * instance; must outlive this JIMWLK.
     */
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

    /**
     * Builds the regularized momentum-space evolution kernel \c K_
     * (idempotent: a no-op if already built). Disables the GSL error
     * handler as a side effect, since getMassRegulator() relies on
     * that to handle a Bessel-function evaluation failure gracefully.
     */
    void initializeK();
    /**
     * Evaluates the long-distance mass regulator
     * \f$m r K_1(m r)\f$ (via the modified Bessel function
     * \f$K_1\f$) at one momentum-space point, used by initializeK() to
     * build \c K_.
     * \param[in] x Momentum-space \f$x\f$ coordinate, in units where
     * the lattice spans \f$[-1/2, 1/2]\f$.
     * \param[in] y Momentum-space \f$y\f$ coordinate, same units.
     * \return `1.0` if `param_.getm_jimwlk()` is `0` (regulator
     * disabled); otherwise \f$m r K_1(m r)\f$, or `0.0` if the Bessel
     * function evaluation fails.
     */
    double getMassRegulator(const double x, const double y) const;
    /**
     * Evaluates the running coupling \f$\alpha_s\f$ at one
     * momentum-space point, used by initializeK() to build \c K_.
     * \param[in] x Momentum-space \f$x\f$ coordinate, in units where
     * the lattice spans \f$[-1/2, 1/2]\f$.
     * \param[in] y Momentum-space \f$y\f$ coordinate, same units.
     * \return `param_.getJimwlk_alphas()` directly if positive (fixed
     * coupling); otherwise a one-loop running-coupling evaluation from
     * the dipole size implied by \f$(x,y)\f$.
     */
    double getAlphas(const double x, const double y) const;
    /**
     * Allocates the per-step Gaussian-noise buffers \c xi_/\c xi2_/\c
     * CKxi_ (idempotent: a no-op if already allocated).
     */
    void initializeNoise();

    /**
     * Runs the full JIMWLK evolution: computes how many Langevin steps
     * each nucleus needs to reach its requested \f$x\f$ (fixed- or
     * running-coupling formula, per `param_.getJimwlk_alphas()`), then
     * evolves the projectile and target in turn via runEvolutionLoop().
     */
    void evolution();
    /**
     * Advances one nucleus' Wilson-line field by one JIMWLK Langevin
     * step: draws Gaussian noise per site/color (via
     * Random::gaussBulk()), FFT-convolves it with \c K_, builds the
     * gauge-covariant combination \f$\sum_a \xi_i^a\, U t^a
     * U^\dagger\f$, FFT-convolves that with \c K_ again, and updates
     * the Wilson line as \f$U \to e^{L} U e^{R}\f$ for the resulting
     * left/right generator combinations \f$L\f$, \f$R\f$.
     * \param[in] nucleus Which Wilson-line field (`lat_ptr_->U` for
     * Projectile, `lat_ptr_->U2` for Target) to evolve.
     */
    void evolutionStep(NucleusRole nucleus);

  private:
    /**
     * Shared by evolution()'s projectile and target passes: runs \p
     * steps Langevin steps for \p nucleus, logging progress and writing
     * any snapshots due in \p xSnapshotList along the way.
     * \param[in] nucleus Which nucleus to evolve.
     * \param[in] steps Number of Langevin steps to run.
     * \param[in] x0 Starting Bjorken \f$x\f$ for this evolution.
     * \param[in] dlogx Logarithmic \f$x\f$ step size per Langevin step.
     * \param[in] saveSnapshots Whether to write a Wilson-line snapshot
     * whenever \f$x\f$ crosses one of \p xSnapshotList's values.
     * \param[in] xSnapshotList Sorted \f$x\f$ values to snapshot at, if
     * \p saveSnapshots.
     */
    void runEvolutionLoop(
        NucleusRole nucleus, int steps, double x0, double dlogx,
        bool saveSnapshots, const std::vector<double> &xSnapshotList);
};

#endif  // SRC_JIMWLK_H_

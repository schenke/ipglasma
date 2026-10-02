// Init.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_INIT_H_
#define SRC_INIT_H_

#include <cstdint>
#include <memory>

#include "FFT.h"
#include "Glauber.h"
#include "Group.h"
#include "Lattice.h"
#include "Matrix.h"
#include "NuclearQsTable.h"
#include "NucleonModel.h"
#include "NucleusSampler.h"
#include "Parameters.h"
#include "PrettyOstream.h"
#include "Random.h"

/**
 * Selects how Init::init() populates the initial Wilson lines.
 */
enum class InitializationMethod {
    /// Sample nucleon positions and color charges, then construct the
    /// Wilson lines from scratch (Init::sampleTA()/setColorChargeDensity()/
    /// setV()).
    SampleColorCharges,
    /// Read previously written Wilson lines from a plain-text file
    /// (WilsonLineIO::read() with format `1`).
    ReadWlineText,
    /// Read previously written Wilson lines from a binary file
    /// (WilsonLineIO::read() with format `2`).
    ReadWlineBinary
};

/**
 * Builds the initial classical Yang-Mills fields: samples nucleon
 * positions and color-charge densities, constructs each nucleus'
 * Wilson line by solving the classical field equations in the presence
 * of its color source, then matches the two nuclei's fields across the
 * forward light cone to get the post-collision gauge fields and
 * electric field that seed the subsequent time evolution.
 */
class Init {
  private:
    /// FFT instance used for the Wilson-line Poisson solve (setV()).
    FFT fft_;
    /// The nuclear \f$Q_s^2(T_p, y)\f$ table, read by init().
    NuclearQsTable qsTable_;

    /// Samples the nucleon positions (\c nucleusA_/\c nucleusB_).
    NucleusSampler nucleusSampler_;

    /// This event's sampled projectile nucleon positions.
    std::vector<ReturnValue> nucleusA_;
    /// This event's sampled target nucleon positions.
    std::vector<ReturnValue> nucleusB_;
    /// Sampled transverse structure of each nucleon of nucleus A, in the
    /// order of \c nucleusA_ (see sampleNucleonProfiles()).
    std::vector<std::unique_ptr<NucleonProfile>> profilesA_;
    /// Same as \c profilesA_, for nucleus B.
    std::vector<std::unique_ptr<NucleonProfile>> profilesB_;

    /// Log sink for progress/warning/error messages.
    PrettyOstream messager_;
    /// Non-owning pointer to the shared Random instance, set by init().
    Random *random_ptr_;

    /// Reusable identity matrix.
    Matrix one_;

  public:
    /**
     * Constructs an Init for a lattice of the given dimensions.
     * \param[in] nn Two-element array `{N, N}`, the transverse lattice
     * dimensions, forwarded to the owned FFT instance.
     */
    explicit Init(const int nn[]) : fft_(nn), one_(1.) {};

    /**
     * Destroys this Init (nothing to release).
     */
    ~Init() {};

    /**
     * Top-level initialization entry point: reads the \f$Q_s^2\f$ table
     * and any pre-tabulated nucleon configurations, then either reads
     * the Wilson lines from disk or samples nucleon positions/color
     * charges and constructs them from scratch, depending on \p
     * init_method.
     * \param[in,out] lat Lattice to populate.
     * \param[in] param Simulation parameters.
     * \param[in] random Non-owning pointer to the shared Random
     * instance; stored in \c random_ptr_.
     * \param[in] glauber Configured Glauber instance providing nuclear
     * geometry.
     * \param[in] init_method Selects how the Wilson lines are obtained
     * (see InitializationMethod).
     */
    void init(
        Lattice *lat, Parameters *param, Random *random, Glauber *glauber,
        InitializationMethod init_method);
    /**
     * Shifts the projectile's and target's Wilson-line fields by
     * \f$\mp b/2\f$ along the impact-parameter direction so they sit at
     * their correct separated positions on the shared lattice, filling
     * with the identity outside the lattice bounds.
     * \param[in,out] lat Lattice whose `U`/`U2` fields are shifted in
     * place.
     * \param[in] param Simulation parameters; `event.b`/`event.phiRP`
     * give the impact parameter and reaction-plane angle.
     */
    void shiftFieldsWithImpactParameter(Lattice *lat, Parameters *param);
    /**
     * Samples this event's impact parameter \f$b\f$ (linearly or
     * uniformly distributed between `collision.bMin`/`collision.bMax`, or `0`
     * for the constant-color-charge-density case) and reaction-plane angle, and
     * resets every nucleon's `.collided` flag to `0`.
     * \param[in,out] param Simulation parameters; `event.b`/`event.phiRP`
     * store the sampled values.
     */
    void sampleImpactParameter(Parameters *param);
    /**
     * Samples this event's nucleon positions of both nuclei into
     * \c nucleusA_/\c nucleusB_ (NucleusSampler::sample()). Both nuclei
     * are centered at the origin.
     * \param[in,out] param Simulation parameters.
     * \param[in,out] random Random-number source.
     * \param[in] glauber Configured Glauber instance providing nuclear
     * geometry.
     */
    void sampleTA(Parameters *param, Random *random, Glauber *glauber);
    /**
     * Sets \f$g^2\mu_A^2\f$/\f$g^2\mu_B^2\f$ at one cell from its
     * already-accumulated \f$T_p^A\f$/\f$T_p^B\f$, via NuclearQsTable::qs2()
     * and (if enabled) the fluctuating-\f$x\f$ iterative solve. Called
     * from setColorChargeDensity()'s per-cell loop
     * \param[in,out] lat Lattice to read \f$T_p\f$ from and write
     * \f$g^2\mu^2\f$ into.
     * \param[in] param Simulation parameters.
     * \param[in] ipos Flat cell index.
     * \param[in] a Lattice spacing [fm].
     * \param[in] rapidityA Effective rapidity for the projectile (see
     * computeEffectiveRapidities()).
     * \param[in] rapidityB Effective rapidity for the target.
     */
    void computeCellColorCharge(
        Lattice *lat, Parameters *param, int ipos, double a, double rapidityA,
        double rapidityB);
    /**
     * computeCellColorCharge()'s `useFluctuatingX` iterative solve
     * for one nucleus's \f$g^2\mu^2\f$ at this cell: iterates the
     * self-consistent local rapidity/\f$Q_s\f$ relation (following the
     * \f$x\f$-dependent \f$Q_s\f$ suppression of \cite Rezaeian:2012ji
     * Eq. (17)) until the rapidity estimate converges to within
     * \f$10^{-3}\f$.
     * \param[in] param Simulation parameters.
     * \param[in] a Lattice spacing [fm].
     * \param[in] rapidity Effective rapidity to start the iteration
     * from.
     * \param[in] Tp Nuclear thickness \f$T_p\f$ at this cell.
     * \param[in] qsmuRatio \f$Q_s/\mu\f$ ratio for this nucleus.
     * \param[in] ySign `+1` for nucleus A, `-1` for nucleus B -- the
     * only sign difference between the two originally-copy-pasted
     * solves.
     * \return The converged \f$g^2\mu^2\f$ [lattice units] (`0` if the
     * iteration ever produces \f$Q_s=0\f$ or NaN).
     */
    double computeFluctuatingXG2mu2(
        Parameters *param, double a, double rapidity, double Tp,
        double qsmuRatio, double ySign);
    /**
     * Sets the color-charge density (\f$g^2\mu_A^2\f$/\f$g^2\mu_B^2\f$)
     * over the whole lattice: either the constant-background case
     * (`!useNucleus`, via setConstantColorChargeDensity()) or the
     * full nucleus pipeline -- sample proton-anisotropy angles and
     * constituent-quark geometry, compute each nucleus' thickness
     * function (smooth or nucleon-summed), then set
     * \f$g^2\mu^2\f$ per cell via computeCellColorCharge().
     * \param[in,out] lat Lattice to populate.
     * \param[in] param Simulation parameters.
     * \param[in,out] random Random-number source.
     * \param[in] glauber Configured Glauber instance.
     */
    void setColorChargeDensity(
        Lattice *lat, Parameters *param, Random *random, Glauber *glauber);
    /**
     * Converts \p param's input rapidity to true rapidity when it's
     * flagged as pseudorapidity, otherwise passes it through unchanged.
     * \param[in] param Simulation parameters; `colorCharge.usePseudoRapidity`
     * selects the conversion, `colorCharge.jacobianMass`/`collision.sqrtS`
     * parameterize it.
     * \param[out] rapidityA Projectile's effective rapidity.
     * \param[out] rapidityB Target's effective rapidity.
     */
    void computeEffectiveRapidities(
        Parameters *param, double &rapidityA, double &rapidityB);
    /**
     * setColorChargeDensity()'s `!useNucleus` (constant \f$g^2\mu\f$
     * background) branch: sets every cell's \f$g^2\mu_A^2\f$/
     * \f$g^2\mu_B^2\f$ to `param->collision.g2mu`'s value, optionally
     * modulated by a fixed Gaussian envelope (`useGaussian`); marks
     * the event a success.
     * \param[in,out] lat Lattice to populate.
     * \param[in] param Simulation parameters.
     */
    void setConstantColorChargeDensity(Lattice *lat, Parameters *param);
    /**
     * Samples each nucleon's proton-anisotropy angle \f$\phi\f$
     * (uniform in \f$[0, 2\pi)\f$ if `param->subnucleon.protonAnisotropy` is
     * nonzero, `0` otherwise) for both nuclei.
     * \param[in] param Simulation parameters.
     * \param[in,out] random Random-number source.
     */
    void sampleNucleonAnisotropyAngles(Parameters *param, Random *random);
    /**
     * Samples the transverse structure of every nucleon of both nuclei
     * (\c profilesA_/\c profilesB_) with the NucleonModel selected by
     * `subnucleon.nucleonModel`.
     * \param[in] param Simulation parameters.
     * \param[in,out] random Random-number source.
     */
    void sampleNucleonProfiles(Parameters *param, Random *random);
    /**
     * Sums the thickness of all nucleons of one nucleus at a lattice
     * position.
     * \param[in] profiles The nucleons' profiles (\c profilesA_ or
     * \c profilesB_).
     * \param[in] x Transverse \f$x\f$ position [fm].
     * \param[in] y Transverse \f$y\f$ position [fm].
     * \param[in] nucleiInAverage Number of nuclei averaged over
     * (`collision.nucleiToAverage`); each \f$T_p\f$ is divided by it.
     * \return \f$\sum T_p\f$ [GeV\f$^2\f$].
     */
    double computeNucleonThicknessAtCell(
        const std::vector<std::unique_ptr<NucleonProfile>> &profiles, double x,
        double y, double nucleiInAverage) const;
    /**
     * setColorChargeDensity()'s `useSmoothNucleus` branch: sets
     * \f$T_p^A\f$/\f$T_p^B\f$ from the smooth (undeformed) Woods-Saxon
     * thickness functions (Glauber::interNuPInSP() for the projectile
     * A, Glauber::interNuTInST() for the target B), centered at the
     * origin and normalized to each nucleus' mass number. The impact
     * parameter is applied later by shiftFieldsWithImpactParameter().
     * \param[in,out] lat Lattice to populate.
     * \param[in] param Simulation parameters.
     * \param[in] glauber Configured Glauber instance.
     */
    void computeSmoothNucleusThickness(
        Lattice *lat, Parameters *param, Glauber *glauber);
    /**
     * setColorChargeDensity()'s default branch: sets \f$T_p^A\f$/\f$T_p^B\f$
     * by summing each sampled nucleon's (constituent-quark or
     * single-Gaussian) thickness, via computeNucleonThicknessAtCell().
     * \param[in,out] lat Lattice to populate.
     * \param[in] param Simulation parameters.
     * \param[in] nucleiInAverage Number of nuclei being averaged over
     * (`param->collision.nucleiToAverage`), used to normalize the
     * sum.
     */
    void computeThicknessFromNucleons(
        Lattice *lat, Parameters *param, double nucleiInAverage);
    /**
     * Determines \f$N_{\text{part}}\f$/\f$N_{\text{coll}}\f$ from the
     * (already-sampled) nucleon positions, writes
     * `NcollList*.dat`/`NpartList*.dat`, and sets
     * `param->event.Npart = `.
     * \param[in] param Simulation parameters.
     * \param[out] Npart Number of participants.
     * \param[out] Ncoll Number of binary collisions.
     * \return `false` (having called `param->event.success = 0`) if
     * `useFixedNpart` is set and this event's \f$N_{\text{part}}\f$
     * doesn't match, signaling the caller to abort and resample;
     * `true` otherwise.
     */
    bool determineNpartAndNcoll(Parameters *param, int &Npart, int &Ncoll);
    /**
     * determineNpartAndNcoll()'s binary-collision pair loop: writes
     * `NcollList<id>.dat` and marks each colliding nucleon pair's
     * `.collided`, using either a hard-sphere (\f$d_{ij}^2 <\f$ \p d2)
     * or Gaussian-profile wounding criterion depending on
     * `param->collision.gaussianWounding`.
     * \param[in] param Simulation parameters.
     * \param[in] d2 Squared wounding distance
     * (\f$\sigma_{NN}/(10\pi)\f$) [fm\f$^2\f$].
     * \param[in] b Impact parameter [fm].
     * \param[in] phiRP Reaction-plane angle [rad].
     * \param[in,out] Ncoll Incremented for each colliding pair found.
     */
    void computeNcollList(
        Parameters *param, double d2, double b, double phiRP, int &Ncoll);
    /**
     * Sets \p param's running-coupling \f$\alpha_s\f$ from whichever
     * \f$Q_s\f$ choice `param->coupling.runWithQs` selects, or a
     * fixed value if running coupling is disabled or \f$\alpha_s\f$ runs with
     * \f$k_T\f$ instead (handled per-cell elsewhere via
     * `computeRunningCouplingGfactor`, which shares computeAlphaS()
     * with this function).
     * \param[in,out] param Simulation parameters;
     * `event.alphas` stores the result.
     */
    void computeAndSetRunningAlphaS(Parameters *param);
    /**
     * Computes and logs this event's collision-geometry summary
     * (\f$N_{\text{part}}\f$, \f$N_{\text{coll}}\f$, \f$T_{pp}\f$,
     * average \f$Q_s\f$, \f$\alpha_s\f$), writes the
     * `usedParameters*.dat`/`NgluonEstimators*.dat` files, and marks
     * the event a success or failure (e.g. no overlap region, no
     * physical \f$Q_s\f$, or \f$Q_{s,\min}^2 S_T\f$ below
     * `param->colorCharge.minimumQs2ST`).
     * \param[in] lat Lattice to read color-charge densities from.
     * \param[in,out] param Simulation parameters.
     */
    void computeCollisionGeometryQuantities(Lattice *lat, Parameters *param);
    /**
     * Scans the full lattice, accumulating the \f$Q_s\f$/\f$T_{pp}\f$
     * collision-geometry averages computeCollisionGeometryQuantities()
     * reports and stores (only over cells within the wounding distance
     * of at least one collided nucleon from each nucleus).
     * \param[in] lat Lattice to read color-charge densities from.
     * \param[in] param Simulation parameters.
     * \param[in] N Lattice side length.
     * \param[in] a Lattice spacing [fm].
     * \param[in] b Impact parameter [fm].
     * \param[in] phiRP Reaction-plane angle [rad].
     * \param[out] averageQs Running sum of \f$Q_s\f$ (max of the two
     * nuclei at each cell).
     * \param[out] averageQs2 Running sum of \f$Q_s^2\f$ (max).
     * \param[out] averageQs2Avg Running sum of \f$Q_s^2\f$ (average of
     * the two nuclei).
     * \param[out] averageQs2min Running sum of \f$Q_s^2\f$ (min).
     * \param[out] averageQs2min2 Running sum of \f$Q_s^2\f$ (min),
     * accumulated over every lattice cell rather than just the
     * overlap region.
     * \param[out] Tpp Running sum of \f$T_p^A T_p^B\f$ [fm\f$^{-2}\f$].
     * \param[out] count Number of cells included in the overlap-region
     * averages.
     */
    void scanCollisionGeometry(
        Lattice *lat, Parameters *param, int N, double a, double b,
        double phiRP, double &averageQs, double &averageQs2,
        double &averageQs2Avg, double &averageQs2min, double &averageQs2min2,
        double &Tpp, int &count);
    /**
     * Logs computeCollisionGeometryQuantities()'s
     * \f$N_{\text{part}}\f$/\f$N_{\text{coll}}\f$/\f$T_{pp}\f$/\f$Q_s\f$/
     * \f$\alpha_s\f$ summary.
     * \param[in] param Simulation parameters.
     * \param[in] Npart Number of participants.
     * \param[in] Ncoll Number of binary collisions.
     * \param[in] Tpp \f$T_{pp}\f$ [fm\f$^{-2}\f$].
     * \param[in] a Lattice spacing [fm].
     * \param[in] averageQs2 Average \f$Q_s^2\f$ (max).
     * \param[in] averageQs2Avg Average \f$Q_s^2\f$ (average).
     * \param[in] averageQs2min Average \f$Q_s^2\f$ (min).
     * \param[in] averageQs2min2 Average \f$Q_s^2\f$ (min), whole
     * lattice.
     * \param[in] count Number of cells the overlap-region averages were
     * computed over.
     */
    void logCollisionGeometryQuantities(
        Parameters *param, int Npart, int Ncoll, double Tpp, double a,
        double averageQs2, double averageQs2Avg, double averageQs2min,
        double averageQs2min2, int count);
    /**
     * Appends this event's `usedParameters<id>.dat` entry (called only
     * when computeCollisionGeometryQuantities() marks the event a
     * success).
     * \param[in] param Simulation parameters.
     * \param[in] phiRP Reaction-plane angle [rad].
     * \param[in] Npart Number of participants.
     * \param[in] Ncoll Number of binary collisions.
     */
    void writeUsedParametersFile(
        Parameters *param, double phiRP, int Npart, int Ncoll);
    /**
     * Writes this event's `NgluonEstimators<id>.dat` file (rough
     * gluon-multiplicity estimator inputs: \f$Q_s^2 S_T\f$ combinations
     * for the min/average/max \f$Q_s\f$ choices).
     * \param[in] param Simulation parameters.
     * \param[in] a Lattice spacing [fm].
     * \param[in] averageQs2 Average \f$Q_s^2\f$ (max).
     * \param[in] averageQs2Avg Average \f$Q_s^2\f$ (average).
     * \param[in] averageQs2min2 Average \f$Q_s^2\f$ (min), whole
     * lattice.
     * \param[in] count Number of cells the overlap-region averages were
     * computed over.
     */
    void writeNgluonEstimatorsFile(
        Parameters *param, double a, double averageQs2, double averageQs2Avg,
        double averageQs2min2, int count);
    /**
     * Constructs both nuclei's Wilson lines by solving the classical
     * color-source Poisson problem in momentum space, one longitudinal
     * sheet (of `param->colorCharge.Ny`) at a time: sample Gaussian
     * color-charge fluctuations, FFT to momentum space, apply the lattice
     * Poisson/UV-damping kernel (computeWilsonLineMomentumKernel()),
     * inverse FFT, exponentiate into an incremental SU(3) rotation
     * (Matrix::fromAlgebraExponent()), and left-multiply onto the running
     * Wilson line. Optionally writes ML training data (`writeOutputs==5`)
     * and/or an initial Wilson-line snapshot.
     * \param[in,out] lat Lattice whose `U`/`U2` are set.
     * \param[in] param Simulation parameters.
     * \param[in,out] random Random-number source.
     */
    void setV(Lattice *lat, Parameters *param, Random *random);
    /**
     * setV()'s lattice Poisson/UV kernel: depends only on transverse
     * momentum and run parameters, so it's computed once and reused for
     * every longitudinal sheet of both nuclei.
     * \param[in] N Lattice side length.
     * \param[in] sites Total number of lattice sites (`N*N`).
     * \param[in] m Infrared mass regulator [lattice units].
     * \param[in] UVdamp UV damping parameter [lattice units].
     * \return The per-site kernel value, length \p sites.
     */
    std::vector<double> computeWilsonLineMomentumKernel(
        int N, int sites, double m, double UVdamp);
    /**
     * setV()'s per-site \f$\sqrt{g^2\mu^2/N_y}\f$ scale, cached once per
     * nucleus instead of being recomputed in every longitudinal sheet.
     * \param[in] lat Lattice to read \f$g^2\mu^2\f$ from.
     * \param[in] sites Total number of lattice sites.
     * \param[in] g Coupling \f$g\f$.
     * \param[in] invNy \f$1/N_y\f$.
     * \param[out] colorChargeScaleA Filled with the projectile's
     * per-site scale, length \p sites.
     * \param[out] colorChargeScaleB Filled with the target's per-site
     * scale, length \p sites.
     */
    void computeWilsonLineColorChargeScales(
        Lattice *lat, int sites, double g, double invNy,
        std::vector<double> &colorChargeScaleA,
        std::vector<double> &colorChargeScaleB);
};

#endif  // SRC_INIT_H_

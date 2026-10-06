// Init.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_INIT_H_
#define SRC_INIT_H_

#include <memory>
#include <vector>

#include "CollisionGeometry.h"
#include "FFT.h"
#include "Glauber.h"
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

/// The rapidities of the two nuclei, see Init::computeEffectiveRapidities().
struct Rapidities {
    /// Rapidity of the projectile (nucleus A).
    double projectile;
    /// Rapidity of the target (nucleus B).
    double target;
};

/// Per-site color-charge scales of the two nuclei, see
/// Init::computeWilsonLineColorChargeScales().
struct ColorChargeScales {
    /// The projectile's scale per site.
    std::vector<double> projectile;
    /// The target's scale per site.
    std::vector<double> target;
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

    /// This event's sampled projectile nucleon positions.
    std::vector<ReturnValue> nucleusA_;
    /// This event's sampled target nucleon positions.
    std::vector<ReturnValue> nucleusB_;
    /// Sampled transverse structure of each nucleon of nucleus A, in the
    /// order of \c nucleusA_ (see sampleNucleonProfiles()).
    std::vector<std::unique_ptr<NucleonProfile>> profilesA_;
    /// Same as \c profilesA_, for nucleus B.
    std::vector<std::unique_ptr<NucleonProfile>> profilesB_;
    /// Samples the nucleon positions (\c nucleusA_/\c nucleusB_).
    NucleusSampler nucleusSampler_;
    /// Impact parameter, wounded nucleons and overlap averages of
    /// \c nucleusA_/\c nucleusB_ (declared after them, which it refers to).
    CollisionGeometry collisionGeometry_ {nucleusA_, nucleusB_};

    /// Log sink for progress/warning/error messages.
    PrettyOstream messager_;
    /// Non-owning pointer to the shared Random instance, set by init().
    Random *random_ptr_ = nullptr;

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

    /// collisionGeometry_ refers to nucleusA_/nucleusB_ of this instance,
    /// so a copy would refer to the original's; nothing needs copying.
    Init(const Init &) = delete;
    /// \copydoc Init(const Init &)
    Init &operator=(const Init &) = delete;

    /**
     * Top-level initialization entry point: reads the \f$Q_s^2\f$ table
     * and any pre-tabulated nucleon configurations, then either reads
     * the Wilson lines and (with `useNucleus 1`) their nuclei's geometry
     * from disk (readGeometry()) or samples nucleon positions/color
     * charges and constructs them from scratch, depending on \p
     * init_method. Either way, the impact parameter is sampled afterwards
     * by main's collision-geometry loop.
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
     * \param[in] param Simulation parameters;
     * `param->event.b`/`param->event.phiRP` give the impact parameter and
     * reaction-plane angle.
     */
    void shiftFieldsWithImpactParameter(Lattice *lat, Parameters *param);
    /**
     * Samples this event's impact parameter and reaction-plane angle
     * (CollisionGeometry::sampleImpactParameter()).
     * \param[in,out] param Simulation parameters;
     * `param->event.b`/`param->event.phiRP` store the sampled values.
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
     * \param[in] param Simulation parameters;
     * `param->colorCharge.usePseudoRapidity` selects the conversion,
     * `param->colorCharge.jacobianMass`/`param->collision.sqrtS` parameterize
     * it.
     * \return The effective rapidities of the projectile and the target.
     */
    Rapidities computeEffectiveRapidities(Parameters *param);
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
     * `param->subnucleon.nucleonModel`.
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
     * (`nucleiToAverage`); each \f$T_p\f$ is divided by it.
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
     * Computes this event's collision geometry and decides whether it is
     * accepted (CollisionGeometry::computeQuantities()).
     * \param[in] lat Lattice to read color-charge densities from.
     * \param[in,out] param Simulation parameters.
     */
    void computeCollisionGeometryQuantities(Lattice *lat, Parameters *param);
    /**
     * Constructs both nuclei's Wilson lines by solving the classical
     * color-source Poisson problem in momentum space, one longitudinal
     * sheet (of `param->colorCharge.Ny`) at a time: sample Gaussian
     * color-charge fluctuations, FFT to momentum space, apply the lattice
     * Poisson/UV-damping kernel (computeWilsonLineMomentumKernel()),
     * inverse FFT, exponentiate into an incremental SU(3) rotation
     * (Matrix::fromAlgebraExponent()), and left-multiply onto the running
     * Wilson line. Optionally writes the binary initial Wilson-line
     * snapshot (`param->output.writeWilsonLineSnapshot`) and/or the
     * Wilson-line files (`param->wilsonLines.writeWilsonLines`).
     * \param[in,out] lat Lattice whose `U`/`U2` are set.
     * \param[in] param Simulation parameters.
     * \param[in,out] random Random-number source.
     */
    void setV(Lattice *lat, Parameters *param, Random *random);
    /**
     * Writes both nuclei's geometry files (WilsonLineIO::writeGeometry())
     * when the event writes Wilson lines (`writeWilsonLines` 1 or 2) of
     * nuclei (`useNucleus 1`), so that readGeometry() can read them back.
     * \param[in] lat Lattice holding the color-charge densities and
     * thicknesses.
     * \param[in] param Simulation parameters.
     */
    void writeGeometry(Lattice *lat, Parameters *param);
    /**
     * Reads both nuclei's geometry files (WilsonLineIO::readGeometry()) for
     * Wilson lines read from disk: sets the color-charge densities and
     * thicknesses on \p lat, the nucleon lists, and
     * `param->colorCharge.QsMuRatio` to the value the files were written
     * with. Exits with an error if the two files disagree on it.
     * \param[in,out] lat Lattice to set.
     * \param[in,out] param Simulation parameters.
     */
    void readGeometry(Lattice *lat, Parameters *param);
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
     * \return The projectile's and the target's per-site scale, length
     * \p sites each.
     */
    ColorChargeScales computeWilsonLineColorChargeScales(
        Lattice *lat, int sites, double g, double invNy);
};

#endif  // SRC_INIT_H_

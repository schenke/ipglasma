// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#ifndef SRC_INIT_H_
#define SRC_INIT_H_

#include <memory>
#include <string>
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
     * The seed of the separate generator of the color charges with
     * `colorChargeSeed` (issue #127): \p colorChargeSeed in the lower and
     * \p eventId in the upper 32 bits, so that every event of a run, and
     * every `colorChargeSeed` below \f$2^{32}\f$, gets its own stream.
     * \param[in] colorChargeSeed The input parameter (not `-1`).
     * \param[in] eventId The event id.
     * \return The seed for Random::init_genrand64().
     */
    static unsigned long long colorChargeRandomSeed(
        long long colorChargeSeed, int eventId) {
        return static_cast<unsigned long long>(colorChargeSeed)
               + (static_cast<unsigned long long>(eventId) << 32);
    }

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
     * init_method. With `runEvolution 1` and `useNucleus 1`, the collision
     * is accepted first (sampleCollision()): the sampled nuclei are drawn
     * together with the impact parameter, before their Wilson lines are
     * built; read nuclei keep their configuration and only the impact
     * parameter is redrawn. Exits with an error if no collision is
     * accepted within the limits of sampleCollision().
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
     * How well the lattice resolves the color-charge density: the median of
     * \f$g^2\mu a\f$ over the given cells, weighted with \f$(g^2\mu)^2\f$
     * like the gluon multiplicity. The multiplicity converges only for
     * values below about 1, the energy more slowly.
     * \param[in] g2mu2 The stored \f$g^2\mu^2\f$ of the cells
     * (Cell::getg2mu2A(), Cell::getg2mu2B()), in lattice units and divided
     * by \f$g^2\f$, so that \f$g^2\mu a = g \sqrt{g^2\mu^2}\f$.
     * \param[in] g The coupling of the stored values.
     * \return The weighted median of \f$g^2\mu a\f$, or 0 if no cell has a
     * color charge.
     */
    static double medianG2muA(std::vector<double> g2mu2, double g);
    /**
     * Prints the median \f$g^2\mu a\f$ of both nuclei of the event
     * (medianG2muA()) and, if it is above 1, a warning that the lattice
     * does not resolve the color-charge density. Does nothing if no cell
     * has a color charge.
     * \param[in] lat Lattice whose cells hold \f$g^2\mu^2\f$ of both nuclei.
     * \param[in] param Simulation parameters (size, L, g).
     */
    void reportColorChargeResolution(Lattice *lat, Parameters *param);
    /**
     * Reads the nuclear \f$Q_s\f$ table that setColorChargeDensity() and
     * computeCellColorCharge() interpolate (NuclearQsTable::read()).
     * \param[in] fileName Path of the table file.
     */
    void readQsTable(const std::string &fileName);
    /**
     * Shifts the projectile's and target's Wilson-line fields by whole
     * lattice sites, CollisionGeometry::latticeShifts() (about
     * \f$\pm b/2\f$ along the impact-parameter direction), so they sit at
     * their separated positions on the shared lattice, filling with the
     * identity where a field is shifted in from outside the lattice. With
     * \f$b = 0\f$ the fields are unchanged.
     * \param[in,out] lat Lattice whose `U`/`U2` fields are shifted in
     * place.
     * \param[in] param Simulation parameters;
     * `param->event.b`/`param->event.phiRP` give the impact parameter and
     * reaction-plane angle.
     */
    void shiftFieldsWithImpactParameter(Lattice *lat, Parameters *param);
    /**
     * Draws collisions until one is accepted: each attempt samples the
     * nuclei (sampleTA(), if \p resampleNuclei), the impact parameter and
     * reaction-plane angle (CollisionGeometry::sampleImpactParameter()) and
     * the wounded nucleons (CollisionGeometry::determineNpartAndNcoll());
     * only an attempt with \f$N_{\text{part}} \ge 2\f$ (and the
     * requested `useFixedNpart`) goes on to the color-charge densities
     * (setColorChargeDensity(), if \p resampleNuclei) and the lattice
     * criteria (CollisionGeometry::computeLatticeQuantities()). Every
     * configuration thus enters with its probability to collide, as in a
     * Monte Carlo Glauber sample.
     * \param[in,out] lat Lattice whose color-charge densities and
     * thicknesses are set (if \p resampleNuclei) and read.
     * \param[in,out] param Simulation parameters; `param->event` holds the
     * accepted collision, with `b`/`phiRP` of the lattice vector
     * separating the nuclei (CollisionGeometry::sampleImpactParameter()).
     * \param[in,out] random Random-number source.
     * \param[in] glauber Configured Glauber instance (only used if \p
     * resampleNuclei).
     * \param[in] resampleNuclei `true` to draw new nuclei for every
     * attempt, `false` to keep the nuclei and color-charge densities on
     * \p lat (Wilson lines read from files).
     * \param[in] maxConfigurations Largest number of attempts.
     * \param[in] maxLatticeRejections Largest number of attempts that
     * reach the lattice criteria and are rejected by them (each costs a
     * scan of the lattice and, with \p resampleNuclei, new color-charge
     * densities).
     * \return `true` if a collision is accepted, `false` (after an error
     * message naming the limit) if a limit is reached.
     */
    bool sampleCollision(
        Lattice *lat, Parameters *param, Random *random, Glauber *glauber,
        bool resampleNuclei, int maxConfigurations = 100000,
        int maxLatticeRejections = 1000);
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
     * and (with `useFluctuatingX 1`) the fluctuating-\f$x\f$ iterative
     * solve. Called from setColorChargeDensity()'s per-cell loop
     * \param[in,out] lat Lattice to read \f$T_p\f$ from and write
     * \f$g^2\mu^2\f$ into.
     * \param[in] param Simulation parameters.
     * \param[in] ipos Flat cell index.
     * \param[in] a Lattice spacing [fm].
     * \param[in] rapidity With `useFluctuatingX 1`, the rapidity the
     * fluctuating-\f$x\f$ iteration starts from (see
     * computeEffectiveRapidity()); not used with a fixed \f$x\f$, which
     * is Parameters::initialX().
     * \return `false` if the fluctuating-\f$x\f$ iteration did not converge
     * for one of the nuclei, `true` otherwise.
     */
    bool computeCellColorCharge(
        Lattice *lat, Parameters *param, int ipos, double a, double rapidity);
    /**
     * computeCellColorCharge()'s `useFluctuatingX` iterative solve
     * for one nucleus's \f$g^2\mu^2\f$ at this cell: iterates the
     * self-consistent local rapidity/\f$Q_s\f$ relation (following the
     * \f$x\f$-dependent \f$Q_s\f$ suppression of \cite Rezaeian:2012ji
     * Eq. (17)) until the rapidity estimate converges to within
     * \f$10^{-3}\f$, for at most 1000 steps.
     * \param[in] param Simulation parameters.
     * \param[in] a Lattice spacing [fm].
     * \param[in] rapidity Effective rapidity to start the iteration
     * from.
     * \param[in] Tp Nuclear thickness \f$T_p\f$ at this cell.
     * \param[in] qsmuRatio \f$Q_s/\mu\f$ ratio for this nucleus.
     * \param[in] ySign `+1` for nucleus A, `-1` for nucleus B -- the
     * only sign difference between the two originally-copy-pasted
     * solves.
     * \param[out] converged `false` if the iteration did not converge in
     * 1000 steps.
     * \return The converged \f$g^2\mu^2\f$ [lattice units] (`0` if the
     * iteration ever produces \f$Q_s=0\f$ or NaN).
     */
    double computeFluctuatingXG2mu2(
        Parameters *param, double a, double rapidity, double Tp,
        double qsmuRatio, double ySign, bool &converged);
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
     * Converts `param->colorCharge.rapidity` to true rapidity when it's
     * flagged as pseudorapidity, otherwise passes it through unchanged.
     * \param[in] param Simulation parameters;
     * `param->colorCharge.usePseudoRapidity` selects the conversion,
     * `param->colorCharge.jacobianMass`/`param->collision.sqrtS` parameterize
     * it.
     * \return The effective rapidity.
     */
    double computeEffectiveRapidity(Parameters *param);
    /**
     * setColorChargeDensity()'s `!useNucleus` (constant \f$g^2\mu\f$
     * background) branch: sets every cell's \f$g^2\mu_A^2\f$/
     * \f$g^2\mu_B^2\f$ to \f$(g^2\mu/g)^2\f$, with \f$g^2\mu\f$ =
     * `param->collision.g2muGeV` converted to lattice units
     * (\f$\times a/\hbar c\f$); marks the event a success.
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
     * \f$T_p^A\f$/\f$T_p^B\f$ from the smooth (undeformed) thickness
     * functions of the density profiles (Glauber::interNuPInSP() for the
     * projectile A, Glauber::interNuTInST() for the target B), centered
     * at the origin and normalized to each nucleus' mass number. For the
     * deuteron, whose Hulthen density is that of the neutron-proton
     * separation \f$r\f$, the thickness is taken at \f$2r\f$, since its
     * nucleons sit at \f$r/2\f$ from the center. The impact
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
     * Constructs both nuclei's Wilson lines by solving the classical
     * color-source Poisson problem in momentum space, one longitudinal
     * sheet (of `param->colorCharge.Ny`) at a time: sample Gaussian
     * color-charge fluctuations, FFT to momentum space, apply the lattice
     * Poisson/UV-damping kernel (computeWilsonLineMomentumKernel()),
     * inverse FFT, exponentiate into an incremental SU(3) rotation
     * (Matrix::fromAlgebraExponent()), and left-multiply onto the running
     * Wilson line. Optionally writes the binary initial Wilson-line
     * snapshot (`param->output.writeWilsonLineSnapshot`), the initial
     * Wilson-line files (`param->wilsonLines.writeWilsonLines`) and the
     * nuclei's geometry files (writeGeometry()).
     * \param[in,out] lat Lattice whose `U`/`U2` are set.
     * \param[in] param Simulation parameters.
     * \param[in,out] random Random-number source.
     */
    void setV(Lattice *lat, Parameters *param, Random *random);
    /**
     * Writes both nuclei's geometry files (WilsonLineIO::writeGeometry())
     * when the event writes Wilson lines (`writeWilsonLines` 1 or 2) of
     * nuclei (`useNucleus 1`) and `writeWilsonLineGeometry` is set, so that
     * readGeometry() can read them back.
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
     * setV()'s lattice Poisson/UV kernel
     * \f$K(\tilde k) = e^{-|\tilde k|\,\mathrm{UVdamp}} /
     * (\tilde k^2 + m^2)\f$, with the lattice momentum
     * \f$\tilde k^2 = 4[\sin^2(k_x/2) + \sin^2(k_y/2)]\f$. For
     * \f$m = 0\f$ the zero mode is removed, \f$K(0) = 0\f$, and with
     * \p removeZeroMode also for \f$m > 0\f$; the UV damping applies for
     * any \f$m\f$. The kernel depends only on transverse
     * momentum and run parameters, so it's computed once and reused for
     * every longitudinal sheet of both nuclei.
     * \param[in] N Lattice side length.
     * \param[in] sites Total number of lattice sites (`N*N`).
     * \param[in] m Infrared mass regulator [lattice units].
     * \param[in] UVdamp UV damping parameter [lattice units].
     * \param[in] removeZeroMode Whether \f$K(0) = 0\f$ also for \p m > 0
     * (`removeZeroMode`); for \p m = 0 it always is.
     * \return The per-site kernel value, length \p sites.
     */
    std::vector<double> computeWilsonLineMomentumKernel(
        int N, int sites, double m, double UVdamp, bool removeZeroMode = false);
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

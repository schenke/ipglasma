// Init.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_INIT_H_
#define SRC_INIT_H_

#include <cstdint>

#include "FFT.h"
#include "Glauber.h"
#include "Group.h"
#include "Lattice.h"
#include "Matrix.h"
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
    /// (Init::readVFromFile() with format `1`).
    ReadWlineText,
    /// Read previously written Wilson lines from a binary file
    /// (Init::readVFromFile() with format `2`).
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
    /// Number of tabulated rapidity points in \c Qs2Nuclear_/readNuclearQs().
    int const static iymaxNuc_ = 44;  // for the Tp-y table

    /// Number of tabulated \f$T_p\f$ points in \c Qs2Nuclear_/\c Tlist_
    /// (updated in Sep 2026 to a 10x extended \f$T_p\f$ range).
    int const static iTpmax_ = 240;

    /// Rapidity spacing [in the same units as \c Qs2Nuclear_'s tabulated
    /// \f$y\f$] between consecutive entries of \c Qs2Nuclear_.
    double const deltaYNuc_ = 0.25;  // for the new table
    /// FFT instance used for the Wilson-line Poisson solve
    /// (setV()) and the forward-lightcone momentum-space kernel.
    FFT fft_;
    /// Tabulated \f$Q_s^2(T_p, y)\f$ values read by readNuclearQs(),
    /// indexed `[iT][iy]`; used by getNuclearQs2()'s bilinear
    /// interpolation.
    double Qs2Nuclear_[iTpmax_][iymaxNuc_];
    /// Tabulated \f$T_p\f$ abscissas corresponding to \c
    /// Qs2Nuclear_'s first index.
    double Tlist_[iTpmax_];

    /// Pre-tabulated nucleon configurations for the projectile, loaded
    /// from file by readInNucleusConfigs() (empty if sampled fresh
    /// instead).
    std::vector<vector<float>> nucleonPosArrA_;
    /// Same as \c nucleonPosArrA_, for the target.
    std::vector<vector<float>> nucleonPosArrB_;

    /// This event's sampled projectile nucleon positions.
    std::vector<ReturnValue> nucleusA_;
    /// This event's sampled target nucleon positions.
    std::vector<ReturnValue> nucleusB_;

    /// Log sink for progress/warning/error messages.
    PrettyOstream messager_;

    /// Number of colors; fixed at 3.
    static constexpr int Nc_ = 3;
    /// SU(3) adjoint dimension, \f$N_c^2-1=8\f$.
    static constexpr int Nc2m1_ = Nc_ * Nc_ - 1;
    /// Non-owning pointer to the shared Group instance, set by init().
    Group *group_ptr_;
    /// Non-owning pointer to the shared Random instance, set by init().
    Random *random_ptr_;

    /// Reusable identity matrix.
    Matrix one_;
    /// Sampled constituent-quark ("hot spot") transverse positions and
    /// widths for the projectile (\c xq1_/\c yq1_/\c BGq1_, one inner
    /// vector per nucleon) and target (\c xq2_/\c yq2_/\c BGq2_); \c
    /// gauss1_/\c gauss2_ hold each nucleon's (or, if substructure is
    /// off, each constituent quark's) \f$Q_s\f$-normalization factor
    /// from sampleQsNormalization(). Populated by
    /// sampleConstituentQuarkGeometry(), consumed by
    /// computeNucleonThicknessAtCell().
    vector<vector<double>> xq1_, xq2_, yq1_, yq2_, BGq1_, BGq2_, gauss1_,
        gauss2_;

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
     * \param[in] group Non-owning pointer to the shared Group instance;
     * stored in \c group_ptr_.
     * \param[in] param Simulation parameters.
     * \param[in] random Non-owning pointer to the shared Random
     * instance; stored in \c random_ptr_.
     * \param[in] glauber Configured Glauber instance providing nuclear
     * geometry.
     * \param[in] init_method Selects how the Wilson lines are obtained
     * (see InitializationMethod).
     */
    void init(
        Lattice *lat, Group *group, Parameters *param, Random *random,
        Glauber *glauber, InitializationMethod init_method);
    /**
     * Shifts the projectile's and target's Wilson-line fields by
     * \f$\mp b/2\f$ along the impact-parameter direction so they sit at
     * their correct separated positions on the shared lattice, filling
     * with the identity outside the lattice bounds.
     * \param[in,out] lat Lattice whose \c U/\c U2 fields are shifted in
     * place.
     * \param[in] param Simulation parameters; `getb()`/`getPhiRP()`
     * give the impact parameter and reaction-plane angle.
     */
    void shiftFieldsWithImpactParameter(Lattice *lat, Parameters *param);
    /**
     * Matches the projectile's and target's Wilson lines across the
     * forward light cone to obtain the post-collision gauge links (\c
     * Ux/\c Uy), electric field (\c U/\c U2, reused as scratch and then
     * as the actual electric-field components), and momentum
     * \f$\pi\f$ (\c Ux2) that seed the subsequent classical Yang-Mills
     * evolution. Runs its steps (see below) in sequence inside one
     * shared `#pragma omp parallel` region; each step has its own
     * `#pragma omp for`, so the implicit barrier at the end of each
     * keeps them correctly ordered.
     * \param[in,out] lat Lattice whose fields are transformed in place.
     * \param[in] param Simulation parameters.
     */
    void initializeForwardLightCone(Lattice *lat, Parameters *param);
    /**
     * Replaces any NaN \c U/\c U2 (left over from a failed
     * forward-lightcone solve at a previous stage) with the identity.
     * \param[in,out] lat Lattice to sanitize.
     * \param[in] N2 Total number of lattice sites.
     */
    void sanitizeForwardLightconeU(Lattice *lat, int N2);
    /**
     * Scratch matrices for computeForwardLightconeLinksTeam(), reused
     * across cells to avoid reallocating.
     */
    struct ForwardLightconeLinkScratch {
        /// Conjugate-transposed neighbor link in \f$x\f$.
        Matrix UDx;
        /// Conjugate-transposed neighbor link in \f$y\f$.
        Matrix UDy;
    };
    /**
     * Computes the pre-collision forward/backward gauge links \c
     * Ux1/\c Uy1 (from \c U) and \c Ux2/\c Uy2 (from \c U2) needed by
     * the forward-lightcone matching.
     * \param[in,out] lat Lattice to read \c U/\c U2 from and write \c
     * Ux1/\c Uy1/\c Ux2/\c Uy2 into.
     * \param[in] N2 Total number of lattice sites.
     * \param[in,out] scratch Thread-local scratch storage.
     */
    void computeForwardLightconeLinksTeam(
        Lattice *lat, int N2, ForwardLightconeLinkScratch &scratch);
    /**
     * Scratch matrices for computeForwardLightconeUxUyTeam(), reused
     * across cells to avoid reallocating.
     */
    struct ForwardLightconeUScratch {
        /// Newton-solved matching result for this cell/direction.
        Matrix temp2;
        /// Projectile-side link in \f$x\f$, passed to
        /// findUInForwardLightcone().
        Matrix UDx1;
        /// Target-side link in \f$x\f$, passed to
        /// findUInForwardLightcone().
        Matrix UDx2;
        /// Projectile-side link in \f$y\f$, passed to
        /// findUInForwardLightcone().
        Matrix UDy1;
        /// Target-side link in \f$y\f$, passed to
        /// findUInForwardLightcone().
        Matrix UDy2;
    };
    /**
     * Solves for the post-collision gauge links \c Ux/\c Uy from \c
     * Ux1/\c Ux2 and \c Uy1/\c Uy2 via findUInForwardLightcone(),
     * logging a warning at each cell where the Newton solve didn't
     * converge.
     * \param[in,out] lat Lattice to read \c Ux1/\c Uy1/\c Ux2/\c Uy2
     * from and write \c Ux/\c Uy into.
     * \param[in] param Simulation parameters; only used for the
     * warning message's cell coordinates and `getRandomSeed()`/
     * `getEventId()` (to seed findUInForwardLightcone()'s deterministic
     * retry stream).
     * \param[in] N2 Total number of lattice sites.
     * \param[in,out] scratch Thread-local scratch storage.
     */
    void computeForwardLightconeUxUyTeam(
        Lattice *lat, Parameters *param, int N2,
        ForwardLightconeUScratch &scratch);
    /**
     * Scratch matrices for computeForwardLightconeElectricFieldTeam(),
     * reused across cells to avoid reallocating.
     */
    struct ForwardLightconeElectricFieldScratch {
        Matrix temp2;
        Matrix Ux1mUx2;
        Matrix UDx1;
        Matrix UDx2;
        Matrix UDx1mUDx2;
        Matrix Ux;
        Matrix UDx;
        Matrix Uy1mUy2;
        Matrix UDy1;
        Matrix UDy2;
        Matrix UDy1mUDy2;
        Matrix Uy;
        Matrix UDy;
    };
    /**
     * Computes the initial electric field's contribution from one
     * direction (\p neighborX/\p neighborY select minus-shifted or
     * plus-shifted neighbors), written into \p outputField. Called
     * once with `(posmX, posmY, lat->U)` and once with `(pospX, pospY,
     * lat->U2)` -- previously two copy-pasted loops.
     * \param[in] lat Lattice to read \c Ux1/\c Uy1/\c Ux2/\c Uy2/\c
     * Ux/\c Uy from.
     * \param[in] N2 Total number of lattice sites.
     * \param[in] neighborX Neighbor-index table to use in \f$x\f$
     * (`lat->posmX` or `lat->pospX`).
     * \param[in] neighborY Neighbor-index table to use in \f$y\f$
     * (`lat->posmY` or `lat->pospY`).
     * \param[in,out] outputField Field this contribution is written
     * into (`lat->U` or `lat->U2`, reused as scratch here ahead of
     * resetForwardLightconeFieldsTeam() zeroing them).
     * \param[in,out] scratch Thread-local scratch storage.
     */
    void computeForwardLightconeElectricFieldTeam(
        Lattice *lat, int N2, const std::vector<int> &neighborX,
        const std::vector<int> &neighborY, std::vector<Matrix> &outputField,
        ForwardLightconeElectricFieldScratch &scratch);
    /**
     * Scratch matrices for computeForwardLightconePlaquetteTeam(),
     * reused across cells to avoid reallocating.
     */
    struct ForwardLightconePlaquetteScratch {
        Matrix UDx;
        Matrix UDy;
        Matrix Uplaq;
    };
    /**
     * Computes the spatial plaquette from \c Ux/\c Uy into \c
     * lat->Uy1 (reused as scratch here, ahead of
     * computeForwardLightconePiTeam()/resetForwardLightconeFieldsTeam()
     * repurposing it further).
     * \param[in,out] lat Lattice to read \c Ux/\c Uy from and write \c
     * Uy1 into.
     * \param[in] N2 Total number of lattice sites.
     * \param[in,out] scratch Thread-local scratch storage.
     */
    void computeForwardLightconePlaquetteTeam(
        Lattice *lat, int N2, ForwardLightconePlaquetteScratch &scratch);
    /**
     * Sets \c lat->Ux2 to the initial momentum \f$\pi\f$ (\f$E^\eta\f$)
     * in lattice units, from \c lat->U (the electric field just
     * computed by computeForwardLightconeElectricFieldTeam()).
     * \param[in,out] lat Lattice to read \c U from and write \c Ux2
     * into.
     * \param[in] param Simulation parameters; `getg()` is used.
     * \param[in] N2 Total number of lattice sites.
     */
    void computeForwardLightconePiTeam(Lattice *lat, Parameters *param, int N2);
    /**
     * Zeroes \c lat->U/\c U2/\c Uy2 and resets \c lat->Ux1 to the
     * identity, now that this event's forward-lightcone fields have
     * been consumed by the steps above.
     * \param[in,out] lat Lattice to reset.
     * \param[in] N2 Total number of lattice sites.
     */
    void resetForwardLightconeFieldsTeam(Lattice *lat, int N2);
    /**
     * Samples this event's impact parameter \f$b\f$ (linearly or
     * uniformly distributed between `getbmin()`/`getbmax()`, or `0` for
     * the constant-color-charge-density case) and reaction-plane angle,
     * and resets every nucleon's `.collided` flag to `0`.
     * \param[in,out] param Simulation parameters; `setb()`/`setPhiRP()`
     * store the sampled values.
     */
    void sampleImpactParameter(Parameters *param);
    /**
     * Samples nucleon positions for both nuclei (via
     * sampleTAWoodsSaxon()/sampleTAFromConfigFiles()/
     * sampleTAFromAlvioliFiles(), depending on
     * `param->getNucleonPositionsFromFile()`) and applies each
     * nucleus' global polarization rotation. Both nuclei are centered
     * at the origin.
     * \param[in,out] param Simulation parameters.
     * \param[in,out] random Random-number source.
     * \param[in] glauber Configured Glauber instance providing nuclear
     * geometry.
     */
    void sampleTA(Parameters *param, Random *random, Glauber *glauber);
    /**
     * Samples nucleon positions via rejection sampling from each
     * nucleus' Woods-Saxon (or deformed Woods-Saxon) thickness
     * function, special-casing a single proton (\f$A=1\f$, placed at
     * the origin) and the deuteron (\f$A=2\f$, via
     * Glauber::sampleTARejection() and its neutron-proton distance
     * convention).
     * \param[in,out] param Simulation parameters; exits with an error
     * if `getAverageOverNuclei() > 1` and either nucleus is a proton.
     * \param[in,out] random Random-number source.
     * \param[in] glauber Configured Glauber instance providing nuclear
     * geometry.
     */
    void sampleTAWoodsSaxon(
        Parameters *param, Random *random, Glauber *glauber);
    /**
     * Samples nucleon positions by drawing a random pre-tabulated
     * configuration from \c nucleonPosArrA_/\c nucleonPosArrB_,
     * falling back to sampleTAWoodsSaxon()-style generation for
     * whichever nucleus has no configurations loaded.
     * \param[in,out] random Random-number source.
     * \param[in] glauber Configured Glauber instance providing nuclear
     * geometry.
     */
    void sampleTAFromConfigFiles(Random *random, Glauber *glauber);
    /**
     * Samples nucleon positions from Alvioli's correlated Pb-208
     * configuration files, via readOneAlvioliNucleus().
     * \param[in,out] random Random-number source.
     * \param[in] glauber Configured Glauber instance providing nuclear
     * geometry; exits with an error unless both nuclei are Pb-208, or
     * the projectile is a proton and the target is Pb-208.
     */
    void sampleTAFromAlvioliFiles(Random *random, Glauber *glauber);
    /**
     * Reads one nucleus's worth of nucleon positions from a randomly
     * selected Alvioli correlated-Pb-208 configuration file, appending
     * them to \p nucleus. \p label (`"A"` or `"B"`) is used only for
     * log messages.
     * \param[in,out] random Random-number source.
     * \param[in] nucleonCount Number of nucleon positions to read
     * (`1` places a single nucleon at the origin instead of reading).
     * \param[in] label Nucleus label for log messages.
     * \param[out] nucleus Appended with the read (or, for
     * `nucleonCount == 1`, synthesized) positions.
     */
    void readOneAlvioliNucleus(
        Random *random, int nucleonCount, const std::string &label,
        std::vector<ReturnValue> &nucleus);
    /**
     * Applies sampleTA()'s global nucleus rotation for one polarization
     * flag; called once each for the projectile and target.
     * \param[in,out] random Random-number source.
     * \param[in] polarizationFlag `0` random 3D rotation
     * (rotateNucleus3D()), `1` longitudinal polarization (random
     * azimuthal rotation only), `2` transverse polarization (rotates
     * \f$J\f$ to the \f$+y\f$ axis).
     * \param[in,out] nucleus Nucleon positions to rotate in place.
     */
    void applyPolarizationRotation(
        Random *random, int polarizationFlag,
        std::vector<ReturnValue> &nucleus);
    /**
     * Reads the \f$Q_s^2(T_p, y)\f$ lookup table (\c Qs2Nuclear_/\c
     * Tlist_) from `param->getNucleusQsTableFileName()`.
     * \param[in] param Simulation parameters; exits with an error if
     * the file doesn't exist or ends prematurely.
     */
    void readNuclearQs(Parameters *param);
    /**
     * Solves the dense linear system \f$J x = F\f$ via GSL LU
     * decomposition, sized for the SU(3) adjoint dimension
     * (\f$8\times8\f$); used by findUInForwardLightcone()'s Newton
     * iteration.
     * \param[in] Jab Row-major \f$8\times8\f$ Jacobian.
     * \param[in] Fa Length-8 right-hand side.
     * \param[out] xvec Filled with the length-8 solution.
     */
    void solveAxb(double *Jab, double *Fa, std::vector<double> &xvec);

    /**
     * Bilinearly interpolates the tabulated \f$Q_s^2(T_p, y)\f$ table
     * (\c Qs2Nuclear_/\c Tlist_).
     * \param[in] T Nuclear thickness \f$T_p\f$ to interpolate at.
     * \param[in] y Rapidity to interpolate at; exits with an error if
     * above the tabulated range.
     * \return Interpolated \f$Q_s^2\f$; `0` if \p T is below the
     * tabulated range, clamped to the maximal tabulated \f$T_p\f$ (with
     * a warning) if above it.
     */
    double getNuclearQs2(double T, double y);
    /**
     * Sets \f$g^2\mu_A^2\f$/\f$g^2\mu_B^2\f$ at one cell from its
     * already-accumulated \f$T_p^A\f$/\f$T_p^B\f$, via getNuclearQs2()
     * and (if enabled) the fluctuating-\f$x\f$ iterative solve. Called
     * from setColorChargeDensity()'s per-cell loop; a no-op at cells
     * further than `param->getRmax()` from both nuclei (unless a smooth
     * nucleus or JIMWLK is in use, in which case every cell is set).
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
     * computeCellColorCharge()'s `useFluctuatingx==1` iterative solve
     * for one nucleus's \f$g^2\mu^2\f$ at this cell: iterates the
     * self-consistent local rapidity/\f$Q_s\f$ relation (following the
     * \f$x\f$-dependent \f$Q_s\f$ suppression of arXiv:1212.2974
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
     * (`useNucleus==0`, via setConstantColorChargeDensity()) or the
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
     * \param[in] param Simulation parameters; `getUsePseudoRapidity()`
     * selects the conversion, `getJacobianm()`/`getRoots()` parameterize
     * it.
     * \param[out] rapidityA Projectile's effective rapidity.
     * \param[out] rapidityB Target's effective rapidity.
     */
    void computeEffectiveRapidities(
        Parameters *param, double &rapidityA, double &rapidityB);
    /**
     * setColorChargeDensity()'s `useNucleus==0` (constant \f$g^2\mu\f$
     * background) branch: sets every cell's \f$g^2\mu_A^2\f$/
     * \f$g^2\mu_B^2\f$ to `param->getg2mu()`'s value, optionally
     * modulated by a fixed Gaussian envelope (`useGaussian==1`); marks
     * the event a success.
     * \param[in,out] lat Lattice to populate.
     * \param[in] param Simulation parameters.
     */
    void setConstantColorChargeDensity(Lattice *lat, Parameters *param);
    /**
     * Samples each nucleon's proton-anisotropy angle \f$\phi\f$
     * (uniform in \f$[0, 2\pi)\f$ if `param->getProtonAnisotropy()` is
     * nonzero, `0` otherwise) for both nuclei.
     * \param[in] param Simulation parameters.
     * \param[in,out] random Random-number source.
     */
    void sampleNucleonAnisotropyAngles(Parameters *param, Random *random);
    /**
     * Samples constituent-quark positions/widths (\c xq1_/\c yq1_/\c
     * BGq1_ etc., via samplePartonPositions()) and each nucleon's
     * \f$Q_s\f$-normalization factor (via sampleQsNormalization()), for
     * both nuclei.
     * \param[in] param Simulation parameters;
     * `getUseConstituentQuarkProton()` selects whether substructure is
     * sampled at all.
     * \param[in,out] random Random-number source.
     */
    void sampleConstituentQuarkGeometry(Parameters *param, Random *random);
    /**
     * setColorChargeDensity()'s `useSmoothNucleus==1` branch: sets
     * \f$T_p^A\f$/\f$T_p^B\f$ from the smooth (undeformed) Woods-Saxon
     * thickness functions (Glauber::interNuTInST()/interNuPInSP()),
     * normalized to each nucleus' mass number; marks the event a
     * success.
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
     * (`param->getAverageOverNuclei()`), used to normalize the sum.
     */
    void computeThicknessFromNucleons(
        Lattice *lat, Parameters *param, double nucleiInAverage);
    /**
     * computeThicknessFromNucleons()'s per-cell, per-nucleus \f$T_p\f$
     * sum (constituent-quark or single-Gaussian, depending on
     * `param->getUseConstituentQuarkProton()`); called once per
     * nucleus.
     * \param[in] param Simulation parameters.
     * \param[in] nucleus Nucleon positions to sum over.
     * \param[in] xq Per-nucleon constituent-quark \f$x\f$ offsets (\c
     * xq1_/\c xq2_); ignored if substructure is off.
     * \param[in] yq Per-nucleon constituent-quark \f$y\f$ offsets.
     * \param[in] BGq Per-nucleon constituent-quark widths.
     * \param[in] gauss Per-nucleon (or per-constituent-quark)
     * \f$Q_s\f$-normalization factors.
     * \param[in] x Cell's \f$x\f$ position [fm].
     * \param[in] y Cell's \f$y\f$ position [fm].
     * \param[in] xi Proton thickness-function anisotropy
     * (`param->getProtonAnisotropy()`); only used in the
     * single-Gaussian branch.
     * \param[in] nucleiInAverage Number of nuclei being averaged over,
     * used to normalize the sum.
     * \return This nucleus' contribution to \f$T_p\f$ at `(x, y)`
     * [GeV\f$^2\f$].
     */
    double computeNucleonThicknessAtCell(
        Parameters *param, const std::vector<ReturnValue> &nucleus,
        const vector<vector<double>> &xq, const vector<vector<double>> &yq,
        const vector<vector<double>> &BGq, const vector<vector<double>> &gauss,
        double x, double y, double xi, double nucleiInAverage);
    /**
     * Determines \f$N_{\text{part}}\f$/\f$N_{\text{coll}}\f$ from the
     * (already-sampled) nucleon positions, writes
     * `NcollList*.dat`/`NpartList*.dat`, and sets
     * `param->setNpart()`.
     * \param[in] param Simulation parameters.
     * \param[out] Npart Number of participants.
     * \param[out] Ncoll Number of binary collisions.
     * \return `false` (having called `param->setSuccess(0)`) if
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
     * `param->getGaussianWounding()`.
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
     * \f$Q_s\f$ choice `param->getRunWithQs()` selects, or a fixed
     * value if running coupling is disabled or \f$\alpha_s\f$ runs
     * with \f$k_T\f$ instead (handled per-cell elsewhere via
     * `computeRunningCouplingGfactor`).
     * \param[in,out] param Simulation parameters;
     * `setalphas()` stores the result.
     */
    void computeAndSetRunningAlphaS(Parameters *param);
    /**
     * Computes and logs this event's collision-geometry summary
     * (\f$N_{\text{part}}\f$, \f$N_{\text{coll}}\f$, \f$T_{pp}\f$,
     * average \f$Q_s\f$, \f$\alpha_s\f$), writes the
     * `usedParameters*.dat`/`NgluonEstimators*.dat` files, and marks
     * the event a success or failure (e.g. no overlap region, no
     * physical \f$Q_s\f$, or \f$Q_{s,\min}^2 S_T\f$ below
     * `param->getMinimumQs2ST()`).
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
     * sheet (of `param->getNy()`) at a time: sample Gaussian color-charge
     * fluctuations, FFT to momentum space, apply the lattice
     * Poisson/UV-damping kernel (computeWilsonLineMomentumKernel()),
     * inverse FFT, exponentiate into an incremental SU(3) rotation
     * (getUfromExponent()), and left-multiply onto the running Wilson
     * line. Optionally writes ML training data (`writeOutputs==5`)
     * and/or an initial Wilson-line snapshot.
     * \param[in,out] lat Lattice whose \c U/\c U2 are set.
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
    /**
     * Reads both nuclei's Wilson lines from disk, using
     * Lattice::generateWilsonLineDataFileName() to build each file
     * name, in text or binary format depending on \p format.
     * \param[in,out] lat Lattice whose \c U/\c U2 are set.
     * \param[in] param Simulation parameters.
     * \param[in] format `1` for plain text, `2` for binary; exits with
     * an error for any other value
     * (`Lattice::IsValidWilsonLineDataFormat()`).
     * \param[in] x If non-negative, embedded in the generated file
     * names (see Lattice::generateWilsonLineDataFileName()); the
     * default `-1` omits it, matching a file written without an
     * explicit \f$x\f$ (e.g. by Init::setV() outside JIMWLK/fluctuating-
     * \f$x\f$ runs).
     */
    void readVFromFile(
        Lattice *lat, Parameters *param, int format, double x = -1);
    /**
     * readVFromFile()'s `format==1` branch: reads one nucleus' Wilson
     * line from a plain-text file, shifting it by \f$\mp b/2\f$ along
     * the impact-parameter direction and dropping any resulting
     * out-of-bounds column (on the low side for the projectile, the
     * high side for the target).
     * \param[in] fileName Path to read from; exits with an error if it
     * doesn't exist.
     * \param[in] param Simulation parameters.
     * \param[in] role Which nucleus this file belongs to (selects the
     * sign of the \f$b/2\f$ shift and which side out-of-bounds columns
     * are dropped on).
     * \param[out] U Wilson-line field to fill (`lat->U` or `lat->U2`).
     */
    void readWilsonLineText(
        const std::string &fileName, Parameters *param, NucleusRole role,
        std::vector<Matrix> &U);
    /**
     * readVFromFile()'s `format==2` branch: reads one nucleus' Wilson
     * line from a binary file (the format Lattice::writeWilsonLines()'s
     * binary mode writes), shifting it by \f$\mp b/2\f$ along the
     * impact-parameter direction.
     * \param[in] fileName Path to read from; exits with an error if it
     * doesn't exist, or if the file's lattice size/physical length
     * don't match \p param.
     * \param[in] param Simulation parameters.
     * \param[in] role Which nucleus this file belongs to (selects the
     * sign of the \f$b/2\f$ shift).
     * \param[out] U Wilson-line field to fill (`lat->U` or `lat->U2`).
     */
    void readWilsonLineBinary(
        const std::string &fileName, Parameters *param, NucleusRole role,
        std::vector<Matrix> &U);

    /**
     * Assembles the SU(3) matrix \f$\exp(i\sum_a Q_a t_a)\f$ from its
     * eight real generator coefficients, via Matrix::expmCoeff().
     * \param[in] Q The eight coefficients \f$Q_1,\ldots,Q_8\f$.
     * \return The assembled matrix (the identity if the exponential's
     * identity-coefficient comes out numerically zero).
     */
    Matrix getUfromExponent(std::vector<double> &Q);
    /**
     * Solves the forward-lightcone matching condition
     * \f$U_1+U_2 = U_{\text{sol}} U_1 U_2 + (U_1 U_2)^\dagger
     * U_{\text{sol}}^\dagger\f$ for \f$U_{\text{sol}}\f$ via Newton's
     * method in \f$U_{\text{sol}}\f$'s eight SU(3) generator
     * coefficients (computeForwardLightconeResidual()/
     * computeForwardLightconeJacobian()/solveAxb()). If the iteration
     * diverges or fails to converge within 2000 iterations, restarts
     * from a fresh initial guess drawn from a deterministic,
     * seed-independent pseudo-random stream (so restarts are
     * reproducible without touching the shared Random state), up to
     * 200 restarts; returns the best (lowest-residual) estimate found
     * across all restarts either way.
     * \param[in] U1 Projectile-side link.
     * \param[in] U2 Target-side link.
     * \param[out] Usol The solved (or best-effort) matching matrix.
     * \param[in] retrySeed Seed for the deterministic restart stream
     * (see \c forwardLightconeRetrySeed() in Init.cpp), unique per
     * cell/direction/event/run so restarts are reproducible but
     * uncorrelated across cells.
     * \return `true` if the iteration converged (residual below
     * \f$10^{-6}\f$); `false` if it exhausted all restarts without
     * converging (a warning is logged in that case).
     */
    bool findUInForwardLightcone(
        Matrix &U1, Matrix &U2, Matrix &Usol, std::uint64_t retrySeed);
    /**
     * findUInForwardLightcone()'s per-iteration residual: computes
     * \f$F_a\f$ (the quantity the Newton iteration drives to zero) into
     * \p Fa.
     * \param[in] U1pU2 \f$U_1+U_2\f$.
     * \param[in] U1pU2dagger \f$(U_1+U_2)^\dagger\f$.
     * \param[in] Usol Current estimate of the matching matrix.
     * \param[in] Usoldagger \f$U_{\text{sol}}^\dagger\f$.
     * \param[in] traceCache Per-generator trace terms that don't depend
     * on \p Usol, precomputed once by findUInForwardLightcone().
     * \param[out] Fa Filled with the eight residual components.
     * \return \f$F_{\text{zero}} = \sum_a |F_a|\f$, the convergence
     * criterion.
     */
    double computeForwardLightconeResidual(
        const Matrix &U1pU2, const Matrix &U1pU2dagger, const Matrix &Usol,
        const Matrix &Usoldagger,
        const std::vector<complex<double>> &traceCache, double *Fa);
    /**
     * Computes the Jacobian \f$dF/d\alpha\f$ into \p Jab, via a
     * numerical finite difference, falling back to an analytical
     * approximation if the numerical one comes out singular. \p alpha
     * is left unchanged on return (each finite-difference step adds
     * then subtracts its own \f$d\alpha_{b_i}\f$).
     * \param[in] U0 \f$U_1 U_2\f$.
     * \param[in] U1pU2 \f$U_1+U_2\f$.
     * \param[in] Usoldagger \f$U_{\text{sol}}^\dagger\f$ at the current
     * estimate.
     * \param[in,out] MtempArr Scratch, reused across calls; filled with
     * \f$t_a (U_1+U_2)\f$ for the analytical fallback.
     * \param[in,out] alpha Current Newton estimate; perturbed and
     * restored in place during the finite-difference computation.
     * \param[out] Jab Filled with the row-major \f$8\times8\f$
     * Jacobian.
     */
    void computeForwardLightconeJacobian(
        const Matrix &U0, const Matrix &U1pU2, const Matrix &Usoldagger,
        std::vector<Matrix> &MtempArr, std::vector<double> &alpha, double *Jab);

    /**
     * Loads pre-tabulated binary nucleon configurations for light
     * nuclei (deuteron through Pb-208) into \p nucleonPosArr, selecting
     * the file by \p nucleusA/\p lightNucleusOption (and, for the
     * deuteron, \p polarizationFlag/\p polJz). A no-op if \p
     * nucleonPosArr is already populated, or if \p nucleusA isn't one
     * of the supported species (in which case nucleon positions are
     * sampled fresh instead, elsewhere).
     * \param[in] nucleusA Mass number of the nucleus to load
     * configurations for.
     * \param[in] lightNucleusOption Selects which configuration variant
     * to load for species with more than one available (Woods-Saxon,
     * variational Monte Carlo, alpha clusters, PGCM, NLEFT).
     * \param[in] polarizationFlag Deuteron only: `0` samples a random
     * polarization file per call.
     * \param[in] polJz Deuteron only: selects the \f$J_z=\pm1\f$
     * configuration file if `|polJz| == 1`.
     * \param[out] nucleonPosArr Filled with one row per configuration
     * read from the file.
     * \param[in] param Simulation parameters;
     * `getNuclearConfigurationsPath()` gives the directory to read
     * from.
     */
    void readInNucleusConfigs(
        const int nucleusA, const int lightNucleusOption,
        const int polarizationFlag, const double polJz,
        vector<vector<float>> &nucleonPosArr, Parameters *param);
    /**
     * Dispatches to the appropriate nucleus-configuration generator
     * based on which deformation parameters are nonzero: plain
     * Woods-Saxon if none are, otherwise one of the deformed variants
     * (selected by whether \p gamma is nonzero -- triaxial vs. axially
     * symmetric -- and whether \p forceDminFlag requests a strictly
     * enforced minimum inter-nucleon distance).
     * \param[in,out] random Random-number source.
     * \param[in] A Mass number.
     * \param[in] Z Atomic number.
     * \param[in] a_WS Woods-Saxon surface diffuseness [fm].
     * \param[in] R_WS Woods-Saxon half-density radius [fm].
     * \param[in] beta2 Quadrupole deformation.
     * \param[in] beta3 Octupole deformation.
     * \param[in] beta4 Hexadecapole deformation.
     * \param[in] gamma Triaxiality angle [rad].
     * \param[in] forceDminFlag Whether to strictly enforce \p d_min via
     * full resampling.
     * \param[in] d_min Minimum inter-nucleon distance [fm].
     * \param[in] dR_np Neutron-skin radius offset [fm].
     * \param[in] da_np Neutron-skin diffuseness offset [fm].
     * \param[out] nucleus Appended with the generated positions.
     */
    void generateNucleusConfiguration(
        Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
        double beta3, double beta4, double gamma, bool forceDminFlag,
        double d_min, double dR_np, double da_np,
        std::vector<ReturnValue> &nucleus);
    /**
     * Generates an undeformed (spherically symmetric) Woods-Saxon
     * nucleon configuration: samples each nucleon's radius via
     * sampleRFromWoodsSaxon() (protons and neutrons from
     * possibly-different \p a_WS/\p R_WS, via \p dR_np/\p da_np), then
     * places them at random angles subject to a best-effort (up to 100
     * retries) minimum-distance rejection, and recenters the result.
     * \param[in,out] random Random-number source.
     * \param[in] A Mass number.
     * \param[in] Z Atomic number.
     * \param[in] a_WS Proton surface diffuseness [fm].
     * \param[in] R_WS Proton half-density radius [fm].
     * \param[in] d_min Minimum inter-nucleon distance [fm].
     * \param[in] dR_np Neutron-skin radius offset [fm].
     * \param[in] da_np Neutron-skin diffuseness offset [fm].
     * \param[out] nucleus Appended with the generated positions.
     */
    void generateNucleusConfigurationWithWoodsSaxon(
        Random *random, int A, int Z, double a_WS, double R_WS, double d_min,
        double dR_np, double da_np, std::vector<ReturnValue> &nucleus);
    /**
     * Generates an axially symmetric (\f$\gamma=0\f$) deformed
     * Woods-Saxon nucleon configuration: samples each nucleon's
     * `(r, cos\theta)` jointly via
     * sampleRAndCosthetaFromDeformedWoodsSaxon(), then \f$\phi\f$
     * uniformly, subject to the same best-effort minimum-distance
     * rejection as generateNucleusConfigurationWithWoodsSaxon().
     * \param[in,out] random Random-number source.
     * \param[in] A Mass number.
     * \param[in] Z Atomic number.
     * \param[in] a_WS Proton surface diffuseness [fm].
     * \param[in] R_WS Proton half-density radius [fm].
     * \param[in] beta2 Quadrupole deformation.
     * \param[in] beta3 Octupole deformation.
     * \param[in] beta4 Hexadecapole deformation.
     * \param[in] d_min Minimum inter-nucleon distance [fm].
     * \param[in] dR_np Neutron-skin radius offset [fm].
     * \param[in] da_np Neutron-skin diffuseness offset [fm].
     * \param[out] nucleus Appended with the generated positions.
     */
    void generateNucleusConfigurationWithDeformedWoodsSaxon(
        Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
        double beta3, double beta4, double d_min, double dR_np, double da_np,
        std::vector<ReturnValue> &nucleus);
    /**
     * generateNucleusConfiguration()'s `forceDminFlag` variant: samples
     * each nucleon's full `(r, \theta, \phi)` jointly against the
     * (possibly triaxial, \p gamma-dependent) deformed Woods-Saxon
     * surface, and strictly enforces \p d_min by fully resampling any
     * candidate position closer than \p d_min to an already-placed
     * nucleon (no retry cap).
     * \param[in,out] random Random-number source.
     * \param[in] A Mass number.
     * \param[in] Z Atomic number.
     * \param[in] a_WS Proton surface diffuseness [fm].
     * \param[in] R_WS Proton half-density radius [fm].
     * \param[in] beta2 Quadrupole deformation.
     * \param[in] beta3 Octupole deformation.
     * \param[in] beta4 Hexadecapole deformation.
     * \param[in] gamma Triaxiality angle [rad].
     * \param[in] d_min Minimum inter-nucleon distance [fm], strictly
     * enforced.
     * \param[in] dR_np Neutron-skin radius offset [fm].
     * \param[in] da_np Neutron-skin diffuseness offset [fm].
     * \param[out] nucleus Appended with the generated positions.
     */
    void generateNucleusConfigurationWithDeformedWoodsSaxonForceDmin(
        Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
        double beta3, double beta4, double gamma, double d_min, double dR_np,
        double da_np, std::vector<ReturnValue> &nucleus);
    /**
     * generateNucleusConfiguration()'s triaxial (\f$\gamma\neq0\f$),
     * non-forced-\f$d_{\min}\f$ variant: like
     * generateNucleusConfigurationWithDeformedWoodsSaxonForceDmin()'s
     * per-nucleon `(r, \theta, \phi)` sampling against the triaxial
     * surface, but without any minimum-distance rejection.
     * \param[in,out] random Random-number source.
     * \param[in] A Mass number.
     * \param[in] Z Atomic number.
     * \param[in] a_WS Proton surface diffuseness [fm].
     * \param[in] R_WS Proton half-density radius [fm].
     * \param[in] beta2 Quadrupole deformation.
     * \param[in] beta3 Octupole deformation.
     * \param[in] beta4 Hexadecapole deformation.
     * \param[in] gamma Triaxiality angle [rad].
     * \param[in] dR_np Neutron-skin radius offset [fm].
     * \param[in] da_np Neutron-skin diffuseness offset [fm].
     * \param[out] nucleus Appended with the generated positions.
     */
    void generateNucleusConfigurationWithDeformedWoodsSaxon2(
        Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
        double beta3, double beta4, double gamma, double dR_np, double da_np,
        std::vector<ReturnValue> &nucleus);
    /**
     * Samples one nucleon's radius from an undeformed Woods-Saxon
     * distribution via rejection sampling: draws \f$r\f$ with density
     * \f$\propto r^2\f$ (the correct 3D volume-element weighting, via
     * inverse-cube-root of a uniform draw) up to a generous cutoff, and
     * accepts it with probability given by fermiDistribution().
     * \param[in,out] random Random-number source.
     * \param[in] a_WS Surface diffuseness [fm].
     * \param[in] R_WS Half-density radius [fm].
     * \return The sampled radius [fm].
     */
    double sampleRFromWoodsSaxon(
        Random *random, double a_WS, double R_WS) const;
    /**
     * Samples one nucleon's radius and polar angle jointly from an
     * axially symmetric (\f$\gamma=0\f$) deformed Woods-Saxon
     * distribution: draws `(r, cos\theta)` uniformly over their
     * respective ranges and accepts with probability given by
     * fermiDistribution() evaluated at the angle-dependent surface
     * radius \f$R(\theta) = R_{WS}(1+\beta_2 Y_{20}+\beta_3
     * Y_{30}+\beta_4 Y_{40})\f$.
     * \param[in,out] random Random-number source.
     * \param[in] a_WS Surface diffuseness [fm].
     * \param[in] R_WS Half-density radius [fm].
     * \param[in] beta2 Quadrupole deformation.
     * \param[in] beta3 Octupole deformation.
     * \param[in] beta4 Hexadecapole deformation.
     * \param[out] r Sampled radius [fm].
     * \param[out] costheta Sampled \f$\cos\theta\f$.
     */
    void sampleRAndCosthetaFromDeformedWoodsSaxon(
        Random *random, double a_WS, double R_WS, double beta2, double beta3,
        double beta4, double &r, double &costheta) const;
    /**
     * Evaluates the Woods-Saxon (Fermi) density profile.
     * \param[in] r Radius [fm].
     * \param[in] R_WS Half-density radius [fm].
     * \param[in] a_WS Surface diffuseness [fm].
     * \return \f$1/(1+\exp((r-R_{WS})/a_{WS}))\f$.
     */
    double fermiDistribution(double r, double R_WS, double a_WS) const;
    /**
     * Evaluates the axially symmetric (\f$m=0\f$) real spherical
     * harmonic \f$Y_{l0}(\theta)\f$ for \f$l=2\f$, `3`, or `4`.
     * \param[in] l Degree; must be `2`, `3`, or `4` (returns `0`
     * otherwise).
     * \param[in] ct \f$\cos\theta\f$.
     * \return \f$Y_{l0}(\theta)\f$.
     */
    double sphericalHarmonics(int l, double ct) const;
    /**
     * Evaluates the real spherical harmonic \f$Y_{22}(\theta,\phi)\f$,
     * used for triaxial (\f$\gamma\neq0\f$) nuclear deformation.
     * \param[in] ct \f$\cos\theta\f$.
     * \param[in] phi Azimuthal angle \f$\phi\f$ [rad].
     * \return \f$Y_{22}(\theta,\phi)\f$.
     */
    double sphericalHarmonicsY22(double ct, double phi) const;
    /**
     * Shifts a set of coordinates so their center of mass sits at the
     * origin.
     * \param[in,out] x \f$x\f$ coordinates, shifted in place.
     * \param[in,out] y \f$y\f$ coordinates, shifted in place.
     * \param[in,out] z \f$z\f$ coordinates, shifted in place.
     */
    void recenterNucleus(
        std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);
    /**
     * Shifts a nucleus' nucleon positions so their center of mass sits
     * at the origin.
     * \param[in,out] nucleus Nucleon positions, shifted in place.
     */
    void recenterNucleus(std::vector<ReturnValue> &nucleus);
    /**
     * Randomly assigns \p Z of a nucleus' nucleons to be protons (the
     * rest neutrons), via a Fisher-Yates shuffle of the nucleon list
     * followed by labeling the first \p Z entries.
     * \param[in,out] random Random-number source.
     * \param[in,out] nucleus Nucleon positions; shuffled in place, with
     * `.proton` set on every entry.
     * \param[in] Z Number of protons (`std::abs(Z)` is used, in case a
     * negative value is ever passed).
     */
    void assignProtons(
        Random *random, std::vector<ReturnValue> &nucleus, const int Z);
    /**
     * Rotates a nucleus' nucleon positions by a fixed azimuthal angle
     * \p phi_global and polar angle \p theta_global (used for
     * longitudinal and transverse polarization, where the rotation
     * axis/angle is determined rather than random).
     * \param[in] phi_global Azimuthal rotation angle [rad].
     * \param[in] theta_global Polar rotation angle [rad].
     * \param[in,out] nucleus Nucleon positions, rotated in place.
     */
    void rotateNucleus(
        double phi_global, double theta_global,
        std::vector<ReturnValue> &nucleus);
    /**
     * Rotates a nucleus' nucleon positions by a uniformly random
     * 3D (Euler-angle) rotation, as needed for an unpolarized
     * (including triaxially deformed) nucleus.
     * \param[in,out] random Random-number source.
     * \param[in,out] nucleus Nucleon positions, rotated in place.
     */
    void rotateNucleus3D(Random *random, std::vector<ReturnValue> &nucleus);

    /**
     * Samples one nucleon's constituent-quark ("hot spot")
     * substructure: the number of quarks (sampleNumberOfPartons()),
     * their radial distances (3D Gaussian if `omega==1`, otherwise a
     * gamma-distribution-based 2D radial profile via
     * Random::sampleGammaInc()), and their angular placement subject to
     * a minimum-distance (`param->getDqmin()`) rejection criterion, with
     * an optional recentering of the constituent-quark center of mass.
     * \param[in] param Simulation parameters.
     * \param[in,out] random Random-number source.
     * \param[out] x_array Filled with each constituent quark's \f$x\f$
     * offset [fm].
     * \param[out] y_array Filled with each constituent quark's \f$y\f$
     * offset [fm].
     * \param[out] z_array Filled with each constituent quark's \f$z\f$
     * offset [fm] (`0` when `omega != 1`, since that branch assumes a
     * 2D transverse profile).
     * \param[out] BGq_array Filled with each constituent quark's width
     * (the same log-normally sampled value for every quark in this
     * nucleon).
     */
    void samplePartonPositions(
        Parameters *param, Random *random, std::vector<double> &x_array,
        std::vector<double> &y_array, std::vector<double> &z_array,
        std::vector<double> &BGq_array);

    /**
     * Draws one sample from a log-normal distribution with the given
     * mean and variance (not the underlying normal's \f$\mu\f$/\f$\sigma\f$).
     * \param[in,out] random Random-number source.
     * \param[in] mean Desired mean of the log-normal distribution.
     * \param[in] variance Desired variance of the log-normal
     * distribution.
     * \return One log-normally distributed sample.
     */
    double sampleLogNormalDistribution(
        Random *random, const double mean, const double variance);

    /**
     * Samples each of \p Nq constituent quarks' (or, if substructure is
     * off, the single nucleon's) \f$Q_s\f$-normalization factor: `1`
     * for every quark if `param->getSmearQs()` is off, otherwise an
     * independent log-normal draw (mean 1) per quark.
     * \param[in,out] random Random-number source.
     * \param[in] param Simulation parameters.
     * \param[in] Nq Number of quarks (or `1` if substructure is off).
     * \param[out] gauss_array Filled with \p Nq normalization factors.
     */
    void sampleQsNormalization(
        Random *random, Parameters *param, const int Nq,
        std::vector<double> &gauss_array);
    /**
     * Samples the number of constituent quarks for one nucleon: the
     * integer part of `param->getNqBase()`, rounded up with probability
     * equal to its fractional part, plus a Poisson-distributed
     * fluctuation (`param->getNqFluc()`'s mean).
     * \param[in,out] random Random-number source.
     * \param[in] param Simulation parameters.
     * \return The sampled quark count, at least `1`.
     */
    int sampleNumberOfPartons(Random *random, Parameters *param);
};

#endif  // SRC_INIT_H_

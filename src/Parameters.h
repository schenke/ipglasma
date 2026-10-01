// Parameters.h is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.

#ifndef SRC_PARAMETERS_H_
#define SRC_PARAMETERS_H_

#include <iosfwd>
#include <string>
#include <vector>

#include "PrettyOstream.h"

class InputFile;

/// Lattice geometry.
struct LatticeParameters {
    /// Lattice side length (recommended to be \f$2^n\f$).
    int size = 0;
    /// Physical lattice size [fm].
    double L = 0.;
};

/// Classical Yang-Mills evolution.
struct EvolutionParameters {
    /// Run mode: `1` run the evolution, `2` analysis with files from
    /// disk.
    int mode = 0;
    /// Maximal evolution time [fm/c].
    double maxtime = 0.;
    /// Whether to use \f$1/Q_s\f$ as the maximal evolution time (`1`)
    /// or the manually entered \c maxtime (`0`).
    int inverseQsForMaxTime = 0;
};

/// Random-number seeding.
struct RandomParameters {
    /// Random seed added to the current time to generate the full seed
    /// (or the full seed itself, depending on \c useTimeForSeed).
    unsigned long long int seed = 0;
    /// Whether to read random seeds from a file (`1`); overrides \c
    /// useTimeForSeed if set.
    int useSeedList = 0;
    /// Whether to use the system time to generate a random seed (`1`)
    /// or not (`0`).
    int useTimeForSeed = 0;
};

/// Collision system, impact parameter and nucleon wounding.
struct CollisionParameters {
    /// Projectile nucleus' species name.
    std::string Projectile;
    /// Target nucleus' species name (see Glauber::findNucleusData()).
    std::string Target;
    /// Center-of-mass energy \f$\sqrt{s}\f$ of the collision [GeV].
    double roots = 0.;
    /// Inelastic nucleon-nucleon cross section [mb].
    double SigmaNN = 0.;
    /// Minimum impact parameter [fm] to sample from.
    double bmin = 0.;
    /// Maximum impact parameter [fm] to sample to.
    double bmax = 0.;
    /// Whether to sample the impact parameter from a linear (`1`) or
    /// uniform (`0`) distribution.
    int samplebFromLinearDistribution = 0;
    /// Whether to randomly rotate the event's reaction plane.
    bool rotateReactionPlane = false;
    /// Whether to use nuclei with finite geometry (`1`) or a constant
    /// \f$g^2\mu\f$ distribution over the lattice (`0`).
    int useNucleus = 0;
    /// Whether to use a Gaussian profile on top of the constant
    /// background (`1`).
    int useGaussian = 0;
    /// \f$g^2\mu\f$ [lattice units], used for the constant
    /// (`useNucleus=0`) color-charge-density mode.
    double g2mu = 0.;
    /// If `0`, don't demand a given \f$N_{\text{part}}\f$; if `>1`,
    /// resample the initial configuration until this
    /// \f$N_{\text{part}}\f$ is reached.
    int useFixedNpart = 0;
    /// Number of nuclei to average over, for a smoother thickness
    /// distribution.
    int averageOverThisManyNuclei = 0;
    /// Whether to use a hard-sphere profile (`0`) or Gaussian cross
    /// section (`1`) to decide whether a nucleon is wounded.
    int gaussianWounding = 0;
};

/// Nucleon positions: configuration files, polarization and the (deformed)
/// Woods-Saxon distribution.
struct NucleusParameters {
    /// Whether to sample nucleon positions (`0`) or read them from a
    /// file (`1`).
    int nucleonPositionsFromFile = 0;
    /// Path to the nuclear configuration files (used when \c
    /// nucleonPositionsFromFile is `1`).
    std::string nuclearConfigurationsPath;
    /// Light-nucleus (carbon, oxygen) sampling method: `1` Woods-Saxon,
    /// `2` variational Monte Carlo, `3` alpha clusters.
    int lightNucleusOption = 0;
    /// Projectile polarization: `0` unpolarized, `1` longitudinally
    /// polarized, `2` transversely polarized.
    int polariztionProjectile = 0;
    /// Same as \c polariztionProjectile, for the target.
    int polariztionTarget = 0;
    /// Projectile's \f$J_z\f$ polarization.
    double polarizationProjectileJz = 0.;
    /// Target's \f$J_z\f$ polarization.
    double polarizationTargetJz = 0.;
    /// Whether to use a smooth Woods-Saxon distribution for a heavy
    /// nucleus (`1`) instead of sampling discrete nucleons.
    int useSmoothNucleus = 0;
    /// Whether `R_WS`/`a_WS`/`beta2`/`beta3`/`beta4`/\c
    /// gamma override a nucleus species' built-in deformation
    /// parameters.
    bool setWSDeformParams = false;
    /// Woods-Saxon half-density radius [fm] override.
    double R_WS = 0.;
    /// Woods-Saxon surface diffuseness [fm] override.
    double a_WS = 0.;
    /// Quadrupole deformation parameter \f$\beta_2\f$ override (e.g.
    /// to test sensitivity in Uranium).
    double beta2 = 0.;
    /// Octupole deformation parameter \f$\beta_3\f$ override.
    double beta3 = 0.;
    /// Hexadecapole deformation parameter \f$\beta_4\f$ override.
    double beta4 = 0.;
    /// Triaxiality angle \f$\gamma\f$ override.
    double gamma = 0.;
    /// Neutron-skin radius offset [fm] between proton and neutron
    /// density profiles.
    double dR_np = 0.;
    /// Neutron-skin diffuseness offset [fm] between proton and neutron
    /// density profiles.
    double da_np = 0.;
    /// Whether to enforce \c d_min as a minimum inter-nucleon
    /// distance when sampling nucleon positions.
    bool force_dmin_flag = false;
    /// Minimum inter-nucleon distance [fm], enforced when \c
    /// force_dmin_flag is set.
    double d_min = 0.;
};

/// Nucleon substructure (constituent quarks) and its Bayesian-posterior
/// parameter sets. With SubNucleonParamType > 0, m, BG, BGq, smearingWidth,
/// NqBase, QsmuRatio and dqMin are overwritten from the posterior set every
/// event.
struct SubnucleonParameters {
    /// Infrared mass regulator [GeV] cutting off the Coulomb tail;
    /// should be of order \f$\Lambda_{QCD}=0.2\f$ GeV.
    double m = 0.;
    /// Width [GeV\f$^{-2}\f$] of the Gaussian describing the proton's
    /// shape, \f$T \sim e^{-b^2/(2B)}\f$.
    double BG = 0.;
    /// Mean width [GeV\f$^{-2}\f$] of the Gaussian describing one
    /// constituent quark's ("hot spot") shape.
    double BGq = 0.;
    /// Variance [GeV\f$^{-4}\f$] of the constituent-quark Gaussian
    /// width.
    double BGqVar = 0.;
    /// Minimum distance [fm] between valence (constituent) quarks.
    double dqMin = 0.;
    /// Gamma-distribution shape parameter for constituent-quark
    /// radial-position sampling (see Random::setGammaIncCDF()); `1`
    /// reduces to plain 3D Gaussian sampling.
    double omega = 0.;
    /// If `>0`, use a proton made of this many constituent quarks
    /// ("hot spots").
    int useConstituentQuarkProton = 0;
    /// Base number of constituent quarks (posterior-fit parameter; see
    /// setParamsWithPosteriorParameterSet()).
    double NqBase = 0.;
    /// Fluctuation in the number of constituent quarks.
    double NqFluc = 0.;
    /// Whether to shift the constituent-quark center of mass to the
    /// origin after sampling hot-spot positions (`1`).
    int shiftConstituentQuarkProtonOrigin = 0;
    /// Anisotropy \f$\xi\f$ of the proton thickness function,
    /// \f$T \propto \exp[-(x^2+\xi y^2)/2B]/(2\pi B\sqrt{\xi})\f$ (an
    /// initial test parameter).
    double protonAnisotropy = 0.;
    /// Selects which Bayesian-posterior-fit table
    /// setParamsWithPosteriorParameterSet() draws from: `1` (variable
    /// \f$N_q\f$) or `2`/`4` (fixed
    /// \f$N_q=3\f$); see loadPosteriorParameterSets().
    int SubNucleonParamType = 0;
    /// Index into the posterior parameter set selected by \c
    /// SubNucleonParamType, modulo the table's row count.
    int SubNucleonParamSet = 0;
    /// Whether to smear \f$Q_s\f$ using a Poisson distribution around
    /// its mean at every transverse position (`1`) or not (`0`).
    int smearQs = 0;
    /// Width of the Gaussian smearing around the mean \f$g^2\mu^2\f$
    /// (parameter \f$\sigma\f$ in Eq. (23) of \cite Mantysaari:2016jaz).
    double smearingWidth = 0.;
    /// UV damping parameter.
    double UVdamp = 0.;
};

/// Saturation scale, color charges and rapidity/x of the nuclei.
struct ColorChargeParameters {
    /// Ratio between \f$Q_s\f$ and \f$\mu\f$ for nucleus A:
    /// \f$Q_s = \text{QsmuRatio} \cdot g^2\mu\f$.
    double QsmuRatio = 0.;
    /// File name of the table giving \f$Q_s^2\f$ as a function of
    /// rapidity \f$Y\f$ and \f$Q_s^2(Y=0)\f$.
    std::string NucleusQsTableFileName;
    /// If `>0`, exclude events with \f$Q_{s,\min}^2 S_T <\f$ this
    /// value (used to trigger on high-multiplicity events).
    int minimumQs2ST = 0;
    /// Longitudinal "resolution" (see Lappi, Eur. Phys. J. C55, 285).
    int Ny = 0;
    /// Rapidity used to pick Bjorken \f$x\f$ from IP-Sat for the
    /// projectile.
    double RapidityA = 0.;
    /// Same as \c RapidityA, for the target.
    double RapidityB = 0.;
    /// Whether `RapidityA`/`RapidityB` hold pseudorapidity instead
    /// of rapidity (`1`), applying the corresponding Jacobian.
    int usePseudoRapidity = 0;
    /// Mass term [GeV] in the Jacobian converting rapidity \f$y\f$ to
    /// pseudorapidity \f$\eta\f$.
    double Jacobianm = 0.;
    /// Whether Bjorken \f$x\f$ should fluctuate as the local \f$Q_s\f$
    /// (`1`, \f$x = Q_s\beta/\sqrt{s}\f$) or always use the input
    /// file's rapidity value (`0`).
    int useFluctuatingx = 0;
    /// Factor \f$\beta\f$ in \f$x = Q_s\beta/\sqrt{s}\f$; only used
    /// when \c useFluctuatingx is set.
    double xFromThisFactorTimesQs = 0.;

    /// Mean of the projectile's and target's rapidity.
    double rapidity() const { return (RapidityA + RapidityB) / 2.; }
};

/// Fixed and running coupling.
struct CouplingParameters {
    /// Coupling \f$g\f$ needed where \f$g^2\mu\f$ does not scale out
    /// (e.g. the Wilson-line exponential).
    double g = 0.;
    /// Whether \f$\alpha_s\f$ should run: `0` fixed, `1` running.
    int runningCoupling = 0;
    /// \f$\mu_0\f$ in the running-coupling formula (keeps it infrared
    /// finite).
    double muZero = 0.;
    /// Controls how smooth the running-coupling cutoff is.
    double c = 0.;
    /// Number of active quark flavors \f$N_f\f$ in the one-loop QCD
    /// beta-function coefficient \f$\beta_0=(11 N_c-2N_f)/3\f$, shared
    /// by the classical-evolution/hydro-output running-coupling formula
    /// (see computeAlphaS()) and JIMWLK::getAlphas().
    int nFlavors = 0;
    /// \f$\Lambda_{QCD}\f$ [GeV] in the running-coupling formula (see
    /// computeAlphaS()); distinct from \c Lambda_QCD_jimwlk, which
    /// scales the separate JIMWLK small-x evolution coupling.
    double LambdaQCD = 0.;
    /// Whether \f$\alpha_s\f$ should run with the maximum (`2`),
    /// average (`1`), or minimum (`0`) of \f$Q_s\f$ from nuclei A and
    /// B.
    int runWith0Min1Avg2MaxQs = 0;
    /// Factor multiplying \f$Q_s\f$ under the log in the running
    /// \f$\alpha_s\f$ formula.
    double runWithThisFactorTimesQs = 0.;
    /// Whether \f$\alpha_s\f$ should run with the local \f$Q_s\f$ from
    /// nuclei A and B (`1`) or the average (`0`); both still use \c
    /// runWith0Min1Avg2MaxQs's max/average/min choice.
    int runWithLocalQs = 0;
    /// Whether \f$\alpha_s\f$ should run with \f$k_T\f$ (`1`) instead;
    /// overrides any \c runWith0Min1Avg2MaxQs-based running if set.
    int runWithkt = 0;
};

/// Observables and output files.
struct OutputParameters {
    /// Whether to compute the gluon multiplicity spectrum (requires
    /// GaugeFix::fftChi()'s Coulomb-gauge fixing).
    bool computeGluonMultiplicity = false;
    /// Whether to read the gluon spectrum \f$dN/d^2k_T\f$ from file and
    /// compute the integrated rate from it (`1`).
    int readMultFromFile = 0;
    /// Whether to write large output files such as hydro input data
    /// (`1`) or not (`0`); see README.md's "writeOutputs" bit values.
    int writeOutputs = 0;
    /// Whether to run the flow-velocity/hydro-output calculation (`1`)
    /// or write only \f$T^{\mu\nu}\f$ at measurement times (`0`).
    int writeEpsilonUHydro = 0;
    /// Whether to write \f$T^{\mu\nu}\f$ as compact binary `.ipgt`
    /// (`1`) or formatted text `.dat` (`0`).
    int writeTmunuBinary = 0;
    /// Whether to collect all output files into one HDF5 file (`1`) or
    /// not (`0`).
    int writeOutputsToHDF5 = 0;
    /// Physical lattice size for the output grid [fm].
    double LOutput = 0.;
    /// Lattice side length for hydro/Tmunu output; must be `<= size`.
    int sizeOutput = 0;
    /// Lattice side length in rapidity for the output data.
    int etaSizeOutput = 0;
    /// Output-grid step size in rapidity.
    double detaOutput = 0.;
};

/// Writing and reading Wilson lines.
struct WilsonLineParameters {
    /// Whether to write generated Wilson lines (before any evolution)
    /// as text (`1`), binary (`2`), or not at all (`0`).
    int writeWilsonLines = 0;
    /// Path to the directory where generated Wilson lines are written
    /// to, or existing ones read from (see
    /// Lattice::generateWilsonLineDataFileName()).
    std::string wilsonLinePath;
    /// Whether to generate initial Wilson lines (`0`), or read them
    /// from plain text (`1`) or binary (`2`).
    int readInitialWilsonLines = 0;
};

/// JIMWLK small-x evolution.
struct JimwlkParameters {
    /// Whether to run JIMWLK evolution before the classical Yang-Mills
    /// stage.
    bool useJIMWLK = false;
    /// Regulator [GeV] in JIMWLK's running \f$\alpha_s(r)\f$ (Eq. (22)
    /// of \cite Mantysaari:2022sux).
    double mu0_jimwlk = 0.;
    /// \f$\Lambda_{QCD}\f$ [GeV] in JIMWLK's running \f$\alpha_s(r)\f$
    /// (Eq. (22) of \cite Mantysaari:2022sux).
    double Lambda_QCD_jimwlk = 0.;
    /// Cutoff smoothness parameter in JIMWLK's running \f$\alpha_s(r)\f$
    /// (Eq. (22) of \cite Mantysaari:2022sux); distinct from \c c, the
    /// analogous parameter in the classical-evolution/hydro-output
    /// running-coupling formula (see computeAlphaS()). \c nFlavors is
    /// shared between the two formulas instead, since \f$N_f\f$ is the
    /// same physical quantity in both.
    double c_jimwlk = 0.;
    /// Infrared regulator [GeV] in the JIMWLK kernel (see Eq. (21) of
    /// \cite Mantysaari:2022sux).
    double m_jimwlk = 0.;
    /// JIMWLK coupling: `0` running coupling, a positive value fixed
    /// coupling.
    double alphas_jimwlk = 0.;
    /// JIMWLK evolution step size (recommended `0.005` with running
    /// coupling, `0.0005` with fixed coupling).
    double Ds_jimwlk = 0.;
    /// Bjorken \f$x\f$ at the initial condition of the JIMWLK
    /// evolution.
    double jimwlk_ic_x = 0.;
    /// Bjorken \f$x\f$ the projectile (nucleus A) is evolved to.
    double x_projectile_jimwlk = 0.;
    /// Bjorken \f$x\f$ the target (nucleus B) is evolved to.
    double x_target_jimwlk = 0.;
    /// Whether to save Wilson-line snapshots at the \f$x\f$ values in
    /// \c xSnapshotList during JIMWLK evolution.
    bool saveSnapshots = false;
    /// Bjorken \f$x\f$ values to save a JIMWLK snapshot at, if \c
    /// saveSnapshots.
    std::vector<double> xSnapshotList;
};

/// Set by the program for the whole run (not read from the input file).
struct RunState {
    /// This process' MPI rank.
    int MPIRank = 0;
    /// Total number of MPI ranks.
    int MPISize = 0;
    /// Evolution time step [lattice units].
    double dtau = 0.;
    /// The random seed actually used this run (so the event can be
    /// reproduced), as opposed to \c seed's configured input.
    unsigned long long int randomSeed = 0;
};

/// Per-event values computed by the program (not read from the input file).
struct EventState {
    /// Current event's identifier (used in output file names).
    int eventId = 0;
    /// Whether a collision happened this event: `0` no collision
    /// (restart), `1` collision happened.
    int success = 0;
    /// Impact parameter [fm] actually used this event.
    double b = 0.;
    /// Reaction-plane angle [rad].
    double phiRP = 0.;
    /// Number of participants.
    int Npart = 0;
    /// Convolution of two nuclear thickness functions,
    /// \f$T_{pp}(b_T) = \sum \delta^2x_T\, T_p(x_T) T_p(x_T-b_T)\f$,
    /// used to weight different impact parameters.
    double Tpp = 0.;
    /// Area [fm\f$^2\f$] of the initial interaction region.
    double area = 0.;
    /// Initial event-plane angle \f$\Psi_2\f$ (geometric/spatial).
    double psi = 0.;
    /// Average \f$Q_s\f$ (maximum of nuclei A and B), used as the
    /// running-coupling scale when \c runWith0Min1Avg2MaxQs selects the
    /// maximum.
    double averageQs = 0.;
    /// Average \f$Q_s\f$ (average of nuclei A and B), used when \c
    /// runWith0Min1Avg2MaxQs selects the average.
    double averageQsAvg = 0.;
    /// Average \f$Q_s\f$ (minimum of nuclei A and B), used when \c
    /// runWith0Min1Avg2MaxQs selects the minimum.
    double averageQsmin = 0.;
    /// \f$\alpha_s\f$ computed at the scale set by the chosen average
    /// \f$Q_s\f$.
    double alphas = 0.;
    /// Same as \c QsmuRatio, for nucleus B.
    double QsmuRatioB = 0.;
};

/**
 * All simulation parameters, grouped by topic: the input-file parameters
 * (read by readInput(), see ParameterTable.cpp for the list) and the
 * values the program sets during a run (\c run) and an event (\c event).
 *
 * The fields are plain public members, e.g. `param->lattice.size` or
 * `param->jimwlk.mu0_jimwlk`; input parameters are named like their
 * input-file key.
 */
class Parameters {
  private:
    /// Log sink for error/info messages.
    PrettyOstream messager_;

    /// Loaded rows of `tables/posterior.csv` (variable \f$N_q\f$),
    /// each `{m, BG, BGq, smearingWidth, NqBase, QsmuRatio, dq_min}`.
    std::vector<std::vector<float>> posteriorParamSets_;

    /// Loaded rows of `tables/posterior_Nq3.csv` or
    /// `tables/posterior5020_Nq3.csv` (fixed \f$N_q=3\f$), each
    /// `{m, BG, BGq, smearingWidth, QsmuRatio, dq_min}`.
    std::vector<std::vector<float>> posteriorParamSetsNq3_;

  public:
    /// Lattice geometry.
    LatticeParameters lattice;
    /// Classical Yang-Mills evolution.
    EvolutionParameters evolution;
    /// Random-number seeding.
    RandomParameters random;
    /// Collision system, impact parameter and nucleon wounding.
    CollisionParameters collision;
    /// Nucleon positions: configuration files, polarization and the (deformed)
    /// Woods-Saxon distribution.
    NucleusParameters nucleus;
    /// Nucleon substructure (constituent quarks) and its Bayesian-posterior
    /// parameter sets. With SubNucleonParamType > 0, m, BG, BGq, smearingWidth,
    /// NqBase, QsmuRatio and dqMin are overwritten from the posterior set every
    /// event.
    SubnucleonParameters subnucleon;
    /// Saturation scale, color charges and rapidity/x of the nuclei.
    ColorChargeParameters colorCharge;
    /// Fixed and running coupling.
    CouplingParameters coupling;
    /// Observables and output files.
    OutputParameters output;
    /// Writing and reading Wilson lines.
    WilsonLineParameters wilsonLines;
    /// JIMWLK small-x evolution.
    JimwlkParameters jimwlk;
    /// Set by the program for the whole run (not read from the input file).
    RunState run;
    /// Per-event values computed by the program (not read from the input file).
    EventState event;

    /**
     * Loads a posterior-fit parameter table from a CSV file (skipping
     * its header row) into \p ParamSet.
     * \param[in] posteriorFileName Path to the CSV file; exits with an
     * error if it can't be opened.
     * \param[out] ParamSet Appended with one row per CSV data line,
     * each a vector of the comma-separated values parsed as `float`.
     */
    void loadPosteriorParameterSetsFromFile(
        std::string posteriorFileName,
        std::vector<std::vector<float>> &ParamSet);
    /**
     * Loads the posterior-fit parameter table selected by \p itype into
     * \c posteriorParamSets_ or \c posteriorParamSetsNq3_.
     * \param[in] itype `1` loads `tables/posterior.csv` (variable
     * \f$N_q\f$); `2` loads `tables/posterior_Nq3.csv`; `4` loads
     * `tables/posterior5020_Nq3.csv` (both fixed \f$N_q=3\f$); any
     * other value is a no-op.
     */
    void loadPosteriorParameterSets(const int itype);
    /**
     * Applies one row of a loaded posterior-fit table to this
     * instance's
     * `subnucleon.m`/`subnucleon.BG`/`subnucleon.BGq`/`subnucleon.smearingWidth`/`subnucleon.NqBase`/\c
     * colorCharge.QsmuRatio/\c subnucleon.dqMin.
     * \param[in] itype `1` uses \c posteriorParamSets_ (also sets \c
     * subnucleon.NqBase from the table); `2` or `4` use \c
     * posteriorParamSetsNq3_ (fixing \c subnucleon.NqBase to `3`); any other
     * value is a no-op.
     * \param[in] iset Row index, taken modulo the selected table's row
     * count.
     */
    void setParamsWithPosteriorParameterSet(const int itype, int iset);

    /**
     * Sets every input parameter from \p input, using the parameter
     * table in ParameterTable.cpp: parses and checks each value, applies
     * the defaults of optional parameters, and then sets the derived
     * values (time step, NqBase, posterior parameter sets). Nothing
     * derived is set if there are errors.
     * \param[in] input The parsed input file.
     * \return Every problem found (including those from reading the
     * file): unknown keys, missing required keys, malformed values and
     * values failing their checks. Empty on success.
     */
    std::vector<std::string> readInput(const InputFile &input);
    /**
     * Writes every input parameter's current value, one `key value` line
     * each in input-file syntax and table order, so the output can be
     * used as an input file again.
     * \param[out] out Stream to write to.
     */
    void writeInputParameters(std::ostream &out) const;
    /**
     * Checks that combine several parameters (e.g. snapshots require
     * writing Wilson lines, LambdaQCD < muZero with running coupling).
     * Checks of a single value run in readInput().
     * \return One message per failed check; empty if all pass.
     */
    std::vector<std::string> validationErrors() const;
    /**
     * Logs every validationErrors() message as an error.
     * \return `true` if there are none.
     */
    bool ValidParameters();
};
#endif  // SRC_PARAMETERS_H_

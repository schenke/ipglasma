// Parameters.h is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.

#ifndef SRC_PARAMETERS_H_
#define SRC_PARAMETERS_H_

#include <string>
#include <vector>

#include "PrettyOstream.h"

/**
 * Every user-configurable and derived simulation parameter, read from
 * the input file (see README.md's "Input parameters" section) by
 * main.cpp/Setup and threaded through to every stage of the pipeline.
 *
 * Almost every field is a plain value plus a `getX()`/`setX()` pair;
 * only the handful of methods below with real bodies (the posterior-
 * parameter-set loaders and ValidParameters()) do anything beyond
 * storing/returning a value.
 */
class Parameters {
  private:
    /// Log sink for error/info messages.
    PrettyOstream messager_;
    /// Selects which Bayesian-posterior-fit table
    /// setParamsWithPosteriorParameterSet() draws from: `1` (variable
    /// \f$N_q\f$) or `2`/`4` (fixed
    /// \f$N_q=3\f$); see loadPosteriorParameterSets().
    int subNucleonParamType_;
    /// Index into the posterior parameter set selected by \c
    /// subNucleonParamType_, modulo the table's row count.
    int subNucleonParamSet_;
    /// Loaded rows of `tables/posterior.csv` (variable \f$N_q\f$),
    /// each `{m, BG, BGq, smearingWidth, NqBase, QsmuRatio, dq_min}`.
    std::vector<std::vector<float>> posteriorParamSets_;
    /// Loaded rows of `tables/posterior_Nq3.csv` or
    /// `tables/posterior5020_Nq3.csv` (fixed \f$N_q=3\f$), each
    /// `{m, BG, BGq, smearingWidth, QsmuRatio, dq_min}`.
    std::vector<std::vector<float>> posteriorParamSetsNq3_;

    /// Lattice side length (recommended to be \f$2^n\f$).
    int size_;
    /// Lattice side length for hydro/Tmunu output; must be `<= size_`.
    int sizeOutput_;
    /// Lattice side length in rapidity for the output data.
    int etaSizeOutput_;
    /// Output-grid step size in rapidity.
    double detaOutput_;
    /// Whether \f$\alpha_s\f$ should run: `0` fixed, `1` running.
    int runningCoupling_;
    /// Whether to use the system time to generate a random seed (`1`)
    /// or not (`0`).
    int useTimeForSeed_;
    /// Whether to read random seeds from a file (`1`); overrides \c
    /// useTimeForSeed_ if set.
    int useSeedList_;
    /// Random seed added to the current time to generate the full seed
    /// (or the full seed itself, depending on \c useTimeForSeed_).
    unsigned long long int seed_;
    /// Longitudinal "resolution" (see Lappi, Eur. Phys. J. C55, 285).
    int Ny_;
    /// \f$g^2\mu\f$ [lattice units], used for the constant
    /// (`useNucleus_=0`) color-charge-density mode.
    double g2mu_;
    /// Run mode: `1` run the evolution, `2` analysis with files from
    /// disk.
    int mode_;
    /// Whether \f$\alpha_s\f$ should run with the maximum (`2`),
    /// average (`1`), or minimum (`0`) of \f$Q_s\f$ from nuclei A and
    /// B.
    int runWithQs_;
    /// Whether \f$\alpha_s\f$ should run with \f$k_T\f$ (`1`) instead;
    /// overrides any \c runWithQs_-based running if set.
    int runWithkt_;
    /// Whether \f$\alpha_s\f$ should run with the local \f$Q_s\f$ from
    /// nuclei A and B (`1`) or the average (`0`); both still use \c
    /// runWithQs_'s max/average/min choice.
    int runWithLocalQs_;
    /// Factor multiplying \f$Q_s\f$ under the log in the running
    /// \f$\alpha_s\f$ formula.
    double runWithThisFactorTimesQs_;
    /// Coupling \f$g\f$ needed where \f$g^2\mu\f$ does not scale out
    /// (e.g. the Wilson-line exponential).
    double g_;
    /// Infrared mass regulator [GeV] cutting off the Coulomb tail;
    /// should be of order \f$\Lambda_{QCD}=0.2\f$ GeV.
    double m_;
    /// Mass term [GeV] in the Jacobian converting rapidity \f$y\f$ to
    /// pseudorapidity \f$\eta\f$.
    double Jacobianm_;
    /// Ratio between \f$Q_s\f$ and \f$\mu\f$ for nucleus A:
    /// \f$Q_s = \text{QsmuRatio} \cdot g^2\mu\f$.
    double QsmuRatio_;
    /// Same as \c QsmuRatio_, for nucleus B.
    double QsmuRatioB_;
    /// Rapidity used to pick Bjorken \f$x\f$ from IP-Sat for the
    /// projectile.
    double rapidityA_;
    /// Same as \c rapidityA_, for the target.
    double rapidityB_;
    /// Whether \c rapidityA_/\c rapidityB_ hold pseudorapidity instead
    /// of rapidity (`1`), applying the corresponding Jacobian.
    int usePseudoRapidity_;
    /// Average \f$Q_s\f$ (maximum of nuclei A and B), used as the
    /// running-coupling scale when \c runWithQs_ selects the maximum.
    double averageQs_;
    /// Average \f$Q_s\f$ (average of nuclei A and B), used when \c
    /// runWithQs_ selects the average.
    double averageQsAvg_;
    /// Average \f$Q_s\f$ (minimum of nuclei A and B), used when \c
    /// runWithQs_ selects the minimum.
    double averageQsmin_;
    /// \f$\alpha_s\f$ computed at the scale set by the chosen average
    /// \f$Q_s\f$.
    double alphas_;
    /// Whether to write large output files such as hydro input data
    /// (`1`) or not (`0`); see README.md's "writeOutputs" bit values.
    int writeOutputs_;
    /// Whether to run the flow-velocity/hydro-output calculation (`1`)
    /// or write only \f$T^{\mu\nu}\f$ at measurement times (`0`).
    int writeEpsilonUHydro_;
    /// Whether to write \f$T^{\mu\nu}\f$ as compact binary `.ipgt`
    /// (`1`) or formatted text `.dat` (`0`).
    int writeTmunuBinary_;
    /// Whether to collect all output files into one HDF5 file (`1`) or
    /// not (`0`).
    int writeOutputsToHDF5_;
    /// Whether to write generated Wilson lines (before any evolution)
    /// as text (`1`), binary (`2`), or not at all (`0`).
    int writeWilsonLines_;
    /// Path to the directory where generated Wilson lines are written
    /// to, or existing ones read from (see
    /// Lattice::generateWilsonLineDataFileName()).
    std::string wilsonLinePath_;
    /// Whether to generate initial Wilson lines (`0`), or read them
    /// from plain text (`1`) or binary (`2`).
    int readInitialWilsonLines_;
    /// The random seed actually used this run (so the event can be
    /// reproduced), as opposed to \c seed_'s configured input.
    unsigned long long int randomSeed_;
    /// File name of the table giving \f$Q_s^2\f$ as a function of
    /// rapidity \f$Y\f$ and \f$Q_s^2(Y=0)\f$.
    std::string nucleusQsTableFileName_;
    /// Width [GeV\f$^{-2}\f$] of the Gaussian describing the proton's
    /// shape, \f$T \sim e^{-b^2/(2B)}\f$.
    double BG_;
    /// Mean width [GeV\f$^{-2}\f$] of the Gaussian describing one
    /// constituent quark's ("hot spot") shape.
    double BGq_;
    /// Variance [GeV\f$^{-4}\f$] of the constituent-quark Gaussian
    /// width.
    double BGqVar_;
    /// Gamma-distribution shape parameter for constituent-quark
    /// radial-position sampling (see Random::setGammaIncCDF()); `1`
    /// reduces to plain 3D Gaussian sampling.
    double omega_;
    /// Minimum distance [fm] between valence (constituent) quarks.
    double dq_min_;
    /// \f$\mu_0\f$ in the running-coupling formula (keeps it infrared
    /// finite).
    double muZero_;
    /// Controls how smooth the running-coupling cutoff is.
    double c_;
    /// Center-of-mass energy \f$\sqrt{s}\f$ of the collision [GeV].
    double roots_;
    /// Whether Bjorken \f$x\f$ should fluctuate as the local \f$Q_s\f$
    /// (`1`, \f$x = Q_s\beta/\sqrt{s}\f$) or always use the input
    /// file's rapidity value (`0`).
    int useFluctuatingx_;
    /// Factor \f$\beta\f$ in \f$x = Q_s\beta/\sqrt{s}\f$; only used
    /// when \c useFluctuatingx_ is set.
    double xFromThisFactorTimesQs_;
    /// Convolution of two nuclear thickness functions,
    /// \f$T_{pp}(b_T) = \sum \delta^2x_T\, T_p(x_T) T_p(x_T-b_T)\f$,
    /// used to weight different impact parameters.
    double Tpp_;
    /// Whether to use \f$1/Q_s\f$ as the maximal evolution time (`1`)
    /// or the manually entered \c maxtime_ (`0`).
    int inverseQsForMaxTime_;
    /// Area [fm\f$^2\f$] of the initial interaction region.
    double area_;
    /// Initial event-plane angle \f$\Psi_2\f$ (geometric/spatial).
    double psi_;
    // Glauber parameters:
    /// Inelastic nucleon-nucleon cross section [mb].
    double sigmaNN_;
    /// Impact parameter [fm] actually used this event.
    double b_;
    /// Minimum impact parameter [fm] to sample from.
    double bmin_;
    /// Maximum impact parameter [fm] to sample to.
    double bmax_;
    /// Reaction-plane angle [rad].
    double phiRP_;
    /// Whether to sample the impact parameter from a linear (`1`) or
    /// uniform (`0`) distribution.
    int linearb_;
    /// Target nucleus' species name (see Glauber::findNucleusData()).
    std::string target_;
    /// Projectile nucleus' species name.
    std::string projectile_;
    /// Physical lattice size [fm].
    double L_;
    /// Physical lattice size for the output grid [fm].
    double LOutput_;
    /// Whether to use nuclei with finite geometry (`1`) or a constant
    /// \f$g^2\mu\f$ distribution over the lattice (`0`).
    int useNucleus_;
    /// Light-nucleus (carbon, oxygen) sampling method: `1` Woods-Saxon,
    /// `2` variational Monte Carlo, `3` alpha clusters.
    int lightNucleusOption_;

    /// Projectile polarization: `0` unpolarized, `1` longitudinally
    /// polarized, `2` transversely polarized.
    int polarizationFlagProjectile_;
    /// Same as \c polarizationFlagProjectile_, for the target.
    int polarizationFlagTarget_;
    /// Projectile's \f$J_z\f$ polarization.
    double polJzProjectile_;
    /// Target's \f$J_z\f$ polarization.
    double polJzTarget_;

    /// Whether to use a Gaussian profile on top of the constant
    /// background (`1`).
    int useGaussian_;
    /// Evolution time step [lattice units].
    double dtau_;
    /// Maximal evolution time [fm/c].
    double maxtime_;
    /// Number of participants.
    int Npart_;
    /// Number of nuclei to average over, for a smoother thickness
    /// distribution.
    int averageOverNuclei_;
    /// Whether to sample nucleon positions (`0`) or read them from a
    /// file (`1`) or Alvioli's correlated Pb-208 files (`2`).
    int nucleonPositionsFromFile_;
    /// Path to the nuclear configuration files (used when \c
    /// nucleonPositionsFromFile_ is `1`).
    std::string nuclearConfigurationsPath_;
    /// If `0`, don't demand a given \f$N_{\text{part}}\f$; if `>1`,
    /// resample the initial configuration until this
    /// \f$N_{\text{part}}\f$ is reached.
    int useFixedNpart_;
    /// Whether to smear \f$Q_s\f$ using a Poisson distribution around
    /// its mean at every transverse position (`1`) or not (`0`).
    int smearQs_;
    /// Width of the Gaussian smearing around the mean \f$g^2\mu^2\f$
    /// (parameter \f$\sigma\f$ in Eq. (23) of arXiv:1607.01711).
    double smearingWidth_;
    /// Whether to use a hard-sphere profile (`0`) or Gaussian cross
    /// section (`1`) to decide whether a nucleon is wounded.
    int gaussianWounding_;
    /// This process' MPI rank.
    int MPIrank_;
    /// Total number of MPI ranks.
    int MPIsize_;
    /// Current event's identifier (used in output file names).
    int event_id_;
    /// Whether a collision happened this event: `0` no collision
    /// (restart), `1` collision happened.
    int success_;
    /// Whether to read the gluon spectrum \f$dN/d^2k_T\f$ from file and
    /// compute the integrated rate from it (`1`).
    int readMultFromFile_;
    /// Radius [fm] at which the per-nucleon thickness distribution is
    /// cut off.
    double rmax_;
    /// Anisotropy \f$\xi\f$ of the proton thickness function,
    /// \f$T \propto \exp[-(x^2+\xi y^2)/2B]/(2\pi B\sqrt{\xi})\f$ (an
    /// initial test parameter).
    double protonAnisotropy_;
    /// If `>0`, use a proton made of this many constituent quarks
    /// ("hot spots").
    int useConstituentQuarkProton_;
    /// Base number of constituent quarks (posterior-fit parameter; see
    /// setParamsWithPosteriorParameterSet()).
    double NqBase_;
    /// Fluctuation in the number of constituent quarks.
    double NqFluc_;
    /// Whether to use a smooth Woods-Saxon distribution for a heavy
    /// nucleus (`1`) instead of sampling discrete nucleons.
    int useSmoothNucleus_;
    /// Whether to shift the constituent-quark center of mass to the
    /// origin after sampling hot-spot positions (`1`).
    int shiftConstituentQuarkProtonOrigin_;
    /// UV damping parameter.
    double UVdamp_;
    /// If `>0`, exclude events with \f$Q_{s,\min}^2 S_T <\f$ this
    /// value (used to trigger on high-multiplicity events).
    int minimumQs2ST_;
    /// Woods-Saxon half-density radius [fm] override.
    double R_WS_;
    /// Woods-Saxon surface diffuseness [fm] override.
    double a_WS_;
    /// Quadrupole deformation parameter \f$\beta_2\f$ override (e.g.
    /// to test sensitivity in Uranium).
    double beta2_;
    /// Octupole deformation parameter \f$\beta_3\f$ override.
    double beta3_;
    /// Hexadecapole deformation parameter \f$\beta_4\f$ override.
    double beta4_;
    /// Triaxiality angle \f$\gamma\f$ override.
    double gamma_;
    /// Minimum inter-nucleon distance [fm], enforced when \c
    /// forceDminFlag_ is set.
    double d_min_;
    /// Whether \c R_WS_/\c a_WS_/\c beta2_/\c beta3_/\c beta4_/\c
    /// gamma_ override a nucleus species' built-in deformation
    /// parameters.
    bool setWSDeformParams_;
    /// Whether to enforce \c d_min_ as a minimum inter-nucleon
    /// distance when sampling nucleon positions.
    bool forceDminFlag_;
    /// Neutron-skin radius offset [fm] between proton and neutron
    /// density profiles.
    double WSdR_np_;
    /// Neutron-skin diffuseness offset [fm] between proton and neutron
    /// density profiles.
    double WSda_np_;

    /// Whether to randomly rotate the event's reaction plane.
    bool rotateReactionPlane_;

    /// Whether to compute the gluon multiplicity spectrum (requires
    /// GaugeFix::fftChi()'s Coulomb-gauge fixing).
    bool computeGluonMultiplicity_;

    /// Whether to run JIMWLK evolution before the classical Yang-Mills
    /// stage.
    bool useJIMWLK_;
    /// Whether JIMWLK uses the simple Langevin discretization of
    /// arXiv:1212.4825.
    bool simpleLangevin_;

    /// JIMWLK coupling: `0` running coupling, a positive value fixed
    /// coupling.
    double jimwlk_alphas_;
    /// Infrared regulator [GeV] in the JIMWLK kernel (see Eq. (21) of
    /// arXiv:2207.03712).
    double m_jimwlk_;
    /// Regulator [GeV] in JIMWLK's running \f$\alpha_s(r)\f$ (Eq. (22)
    /// of arXiv:2207.03712).
    double mu0_jimwlk_;
    /// \f$\Lambda_{QCD}\f$ [GeV] in JIMWLK's running \f$\alpha_s(r)\f$
    /// (Eq. (22) of arXiv:2207.03712).
    double LambdaQCD_jimwlk_;
    /// JIMWLK evolution step size (recommended `0.005` with running
    /// coupling, `0.0005` with fixed coupling).
    double ds_jimwlk_;
    /// Bjorken \f$x\f$ at the initial condition of the JIMWLK
    /// evolution.
    double x0_jimwlk_;

    /// Bjorken \f$x\f$ the projectile (nucleus A) is evolved to.
    double jimwlk_x1_;
    /// Bjorken \f$x\f$ the target (nucleus B) is evolved to.
    double jimwlk_x2_;
    /// Whether to save Wilson-line snapshots at the \f$x\f$ values in
    /// \c xSnapshotList_ during JIMWLK evolution.
    bool saveSnapshots_;
    /// Bjorken \f$x\f$ values to save a JIMWLK snapshot at, if \c
    /// saveSnapshots_.
    std::vector<double> xSnapshotList_;

  public:
    /**
     * Constructs a Parameters with every field default-initialized
     * (i.e. left indeterminate for scalar types); callers populate it
     * via Setup/main.cpp before use.
     */
    Parameters() {};

    // functions to access the private variables:
    /**
     * Sets which Bayesian-posterior-fit table type is in use.
     * \param[in] paramType New value.
     */
    void setSubNucleonParamType(int paramType) {
        subNucleonParamType_ = paramType;
    }
    /**
     * Returns which Bayesian-posterior-fit table type is in use.
     * \return The stored value.
     */
    int getSubNucleonParamType() const { return (subNucleonParamType_); }
    /**
     * Sets the posterior parameter set index.
     * \param[in] paramSet New value.
     */
    void setSubNucleonParamSet(int paramSet) { subNucleonParamSet_ = paramSet; }
    /**
     * Returns the posterior parameter set index.
     * \return The stored value.
     */
    int getSubNucleonParamSet() const { return (subNucleonParamSet_); }
    /**
     * Sets the random seed.
     * \param[in] x New seed value.
     */
    void setSeed(unsigned long long int x) { seed_ = x; }
    /**
     * Returns the configured random seed.
     * \return The stored seed value.
     */
    unsigned long long int getSeed() const { return seed_; }
    /**
     * Sets the longitudinal resolution \c Ny_.
     * \param[in] x New value.
     */
    void setNy(int x) { Ny_ = x; }
    /**
     * Returns the longitudinal resolution.
     * \return The stored value.
     */
    int getNy() const { return Ny_; }
    /**
     * Sets the lattice side length.
     * \param[in] x New side length.
     */
    void setSize(int x) { size_ = x; }
    /**
     * Returns the lattice side length.
     * \return The stored side length.
     */
    int getSize() const { return size_; }
    /**
     * Sets the proton thickness-function anisotropy \f$\xi\f$.
     * \param[in] x New value.
     */
    void setProtonAnisotropy(double x) { protonAnisotropy_ = x; }
    /**
     * Returns the proton thickness-function anisotropy.
     * \return The stored value.
     */
    double getProtonAnisotropy() const { return protonAnisotropy_; }

    /**
     * Sets the output lattice side length.
     * \param[in] x New value.
     */
    void setSizeOutput(int x) { sizeOutput_ = x; }
    /**
     * Returns the output lattice side length.
     * \return The stored value.
     */
    int getSizeOutput() const { return sizeOutput_; }
    /**
     * Sets the output lattice side length in rapidity.
     * \param[in] x New value.
     */
    void setEtaSizeOutput(int x) { etaSizeOutput_ = x; }
    /**
     * Returns the output lattice side length in rapidity.
     * \return The stored value.
     */
    int getEtaSizeOutput() const { return etaSizeOutput_; }
    /**
     * Sets the output-grid step size in rapidity.
     * \param[in] x New value.
     */
    void setDetaOutput(double x) { detaOutput_ = x; }
    /**
     * Returns the output-grid step size in rapidity.
     * \return The stored value.
     */
    double getDetaOutput() const { return detaOutput_; }

    /**
     * Sets the number of nuclei to average over.
     * \param[in] x New value.
     */
    void setAverageOverNuclei(int x) { averageOverNuclei_ = x; }
    /**
     * Returns the number of nuclei to average over.
     * \return The stored value.
     */
    int getAverageOverNuclei() const { return averageOverNuclei_; }
    /**
     * Sets \f$g^2\mu\f$ for the constant color-charge-density mode.
     * \param[in] x New value [lattice units].
     */
    void setg2mu(double x) { g2mu_ = x; }
    /**
     * Returns \f$g^2\mu\f$ for the constant color-charge-density mode.
     * \return The stored value [lattice units].
     */
    double getg2mu() const { return g2mu_; }
    /**
     * Sets the run mode.
     * \param[in] x New value (`1` evolution, `2` analysis).
     */
    void setMode(int x) { mode_ = x; };
    /**
     * Returns the run mode.
     * \return The stored value.
     */
    int getMode() const { return mode_; }
    /**
     * Sets whether \f$\alpha_s\f$ runs.
     * \param[in] x New value (`0` fixed, `1` running).
     */
    void setRunningCoupling(int x) { runningCoupling_ = x; };
    /**
     * Returns whether \f$\alpha_s\f$ runs.
     * \return The stored value.
     */
    int getRunningCoupling() const { return runningCoupling_; }
    /**
     * Sets the coupling \f$g\f$.
     * \param[in] x New value.
     */
    void setg(double x) { g_ = x; }
    /**
     * Returns the coupling \f$g\f$.
     * \return The stored value.
     */
    double getg() const { return g_; }
    /**
     * Sets the inelastic nucleon-nucleon cross section.
     * \param[in] x New value [mb].
     */
    void setSigmaNN(double x) { sigmaNN_ = x; }
    /**
     * Returns the inelastic nucleon-nucleon cross section.
     * \return The stored value [mb].
     */
    double getSigmaNN() const { return sigmaNN_; }
    /**
     * Sets the impact parameter used this event.
     * \param[in] x New value [fm].
     */
    void setb(double x) { b_ = x; }
    /**
     * Returns the impact parameter used this event.
     * \return The stored value [fm].
     */
    double getb() const { return b_; }
    /**
     * Sets the reaction-plane angle.
     * \param[in] x New value [rad].
     */
    void setPhiRP(double x) { phiRP_ = x; }
    /**
     * Returns the reaction-plane angle.
     * \return The stored value [rad].
     */
    double getPhiRP() const { return phiRP_; }
    /**
     * Sets the minimum impact parameter to sample from.
     * \param[in] x New value [fm].
     */
    void setbmin(double x) { bmin_ = x; }
    /**
     * Returns the minimum impact parameter to sample from.
     * \return The stored value [fm].
     */
    double getbmin() const { return bmin_; }
    /**
     * Sets the maximum impact parameter to sample to.
     * \param[in] x New value [fm].
     */
    void setbmax(double x) { bmax_ = x; }
    /**
     * Returns the maximum impact parameter to sample to.
     * \return The stored value [fm].
     */
    double getbmax() const { return bmax_; }
    /**
     * Sets the target nucleus' species name.
     * \param[in] x New species name.
     */
    void setTarget(std::string x) { target_ = x; }
    /**
     * Returns the target nucleus' species name.
     * \return The stored species name.
     */
    std::string getTarget() const { return target_; }
    /**
     * Sets the projectile nucleus' species name.
     * \param[in] x New species name.
     */
    void setProjectile(std::string x) { projectile_ = x; }
    /**
     * Returns the projectile nucleus' species name.
     * \return The stored species name.
     */
    std::string getProjectile() const { return projectile_; }
    /**
     * Sets the physical lattice size.
     * \param[in] x New value [fm].
     */
    void setL(double x) { L_ = x; }
    /**
     * Returns the physical lattice size.
     * \return The stored value [fm].
     */
    double getL() const { return L_; }
    /**
     * Sets the physical output-grid size.
     * \param[in] x New value [fm].
     */
    void setLOutput(double x) { LOutput_ = x; }
    /**
     * Returns the physical output-grid size.
     * \return The stored value [fm].
     */
    double getLOutput() const { return LOutput_; }
    /**
     * Sets the infrared mass regulator.
     * \param[in] x New value [GeV].
     */
    void setm(double x) { m_ = x; }
    /**
     * Returns the infrared mass regulator.
     * \return The stored value [GeV].
     */
    double getm() const { return m_; }
    /**
     * Sets the rapidity-to-pseudorapidity Jacobian mass term.
     * \param[in] x New value [GeV].
     */
    void setJacobianm(double x) { Jacobianm_ = x; }
    /**
     * Returns the rapidity-to-pseudorapidity Jacobian mass term.
     * \return The stored value [GeV].
     */
    double getJacobianm() const { return Jacobianm_; }
    /**
     * Sets the \f$Q_s/\mu\f$ ratio for nucleus A.
     * \param[in] x New value.
     */
    void setQsmuRatio(double x) { QsmuRatio_ = x; }
    /**
     * Returns the \f$Q_s/\mu\f$ ratio for nucleus A.
     * \return The stored value.
     */
    double getQsmuRatio() const { return QsmuRatio_; }
    /**
     * Sets the \f$Q_s/\mu\f$ ratio for nucleus B.
     * \param[in] x New value.
     */
    void setQsmuRatioB(double x) { QsmuRatioB_ = x; }
    /**
     * Returns the \f$Q_s/\mu\f$ ratio for nucleus B.
     * \return The stored value.
     */
    double getQsmuRatioB() const { return QsmuRatioB_; }
    /**
     * Sets the projectile's rapidity.
     * \param[in] x New value.
     */
    void setRapidityA(double x) { rapidityA_ = x; }
    /**
     * Returns the projectile's rapidity.
     * \return The stored value.
     */
    double getRapidityA() const { return rapidityA_; }
    /**
     * Sets the target's rapidity.
     * \param[in] x New value.
     */
    void setRapidityB(double x) { rapidityB_ = x; }
    /**
     * Returns the target's rapidity.
     * \return The stored value.
     */
    double getRapidityB() const { return rapidityB_; }
    /**
     * Returns the mean of the projectile's and target's rapidity.
     * \return \f$(\text{rapidityA} + \text{rapidityB})/2\f$.
     */
    double getRapidity() const { return (rapidityA_ + rapidityB_) / 2.; }
    /**
     * Sets the maximal evolution time.
     * \param[in] x New value [fm/c].
     */
    void setMaxtime(double x) { maxtime_ = x; }
    /**
     * Returns the maximal evolution time.
     * \return The stored value [fm/c].
     */
    double getMaxtime() const { return maxtime_; }
    /**
     * Sets the evolution time step.
     * \param[in] x New value [lattice units].
     */
    void setdtau(double x) { dtau_ = x; }
    /**
     * Returns the evolution time step.
     * \return The stored value [lattice units].
     */
    double getdtau() const { return dtau_; }
    /**
     * Sets the number of participants.
     * \param[in] x New value.
     */
    void setNpart(int x) { Npart_ = x; };
    /**
     * Returns the number of participants.
     * \return The stored value.
     */
    int getNpart() const { return Npart_; }
    /**
     * Sets the average (maximum) \f$Q_s\f$.
     * \param[in] x New value.
     */
    void setAverageQs(double x) { averageQs_ = x; }
    /**
     * Returns the average (maximum) \f$Q_s\f$.
     * \return The stored value.
     */
    double getAverageQs() const { return averageQs_; }
    /**
     * Sets the average (mean) \f$Q_s\f$.
     * \param[in] x New value.
     */
    void setAverageQsAvg(double x) { averageQsAvg_ = x; }
    /**
     * Returns the average (mean) \f$Q_s\f$.
     * \return The stored value.
     */
    double getAverageQsAvg() const { return averageQsAvg_; }
    /**
     * Sets the average (minimum) \f$Q_s\f$.
     * \param[in] x New value.
     */
    void setAverageQsmin(double x) { averageQsmin_ = x; }
    /**
     * Returns the average (minimum) \f$Q_s\f$.
     * \return The stored value.
     */
    double getAverageQsmin() const { return averageQsmin_; }
    /**
     * Sets \f$\alpha_s\f$ computed at the chosen average \f$Q_s\f$.
     * \param[in] x New value.
     */
    void setalphas(double x) { alphas_ = x; }
    /**
     * Returns \f$\alpha_s\f$ computed at the chosen average
     * \f$Q_s\f$.
     * \return The stored value.
     */
    double getalphas() const { return alphas_; }
    /**
     * Sets the random seed actually used this run.
     * \param[in] x New value.
     */
    void setRandomSeed(unsigned long long int x) { randomSeed_ = x; };
    /**
     * Returns the random seed actually used this run.
     * \return The stored value.
     */
    unsigned long long int getRandomSeed() const { return randomSeed_; }
    /**
     * Sets whether to use the system time to generate a seed.
     * \param[in] x New value (`1` use system time, `0` don't).
     */
    void setUseTimeForSeed(int x) { useTimeForSeed_ = x; };
    /**
     * Returns whether to use the system time to generate a seed.
     * \return The stored value.
     */
    int getUseTimeForSeed() const { return useTimeForSeed_; }
    /**
     * Sets whether to read random seeds from a file.
     * \param[in] x New value.
     */
    void setUseSeedList(int x) { useSeedList_ = x; };
    /**
     * Returns whether to read random seeds from a file.
     * \return The stored value.
     */
    int getUseSeedList() const { return useSeedList_; }
    /**
     * Sets the \f$Q_s^2(Y)\f$ lookup table's file name.
     * \param[in] x New file name.
     */
    void setNucleusQsTableFileName(std::string x) {
        nucleusQsTableFileName_ = x;
    }
    /**
     * Returns the \f$Q_s^2(Y)\f$ lookup table's file name.
     * \return The stored file name.
     */
    std::string getNucleusQsTableFileName() const {
        return nucleusQsTableFileName_;
    }
    /**
     * Sets the proton-width Gaussian parameter \c BG_.
     * \param[in] x New value [GeV\f$^{-2}\f$].
     */
    void setBG(double x) { BG_ = x; }
    /**
     * Returns the proton-width Gaussian parameter \c BG_.
     * \return The stored value [GeV\f$^{-2}\f$].
     */
    double getBG() const { return BG_; }
    /**
     * Sets the constituent-quark mean width \c BGq_.
     * \param[in] x New value [GeV\f$^{-2}\f$].
     */
    void setBGq(double x) { BGq_ = x; }
    /**
     * Returns the constituent-quark mean width \c BGq_.
     * \return The stored value [GeV\f$^{-2}\f$].
     */
    double getBGq() const { return BGq_; }
    /**
     * Sets the constituent-quark width variance.
     * \param[in] BGqVar New value [GeV\f$^{-4}\f$].
     */
    void setBGqVar(double BGqVar) { BGqVar_ = BGqVar; }
    /**
     * Returns the constituent-quark width variance.
     * \return The stored value [GeV\f$^{-4}\f$].
     */
    double getBGqVar() const { return BGqVar_; }
    /**
     * Sets the gamma-distribution shape parameter \c omega_.
     * \param[in] x New value.
     */
    void setOmega(double x) { omega_ = x; }
    /**
     * Returns the gamma-distribution shape parameter \c omega_.
     * \return The stored value.
     */
    double getOmega() const { return omega_; }
    /**
     * Sets the minimum inter-quark distance.
     * \param[in] dq_min New value [fm].
     */
    void setDqmin(double dq_min) { dq_min_ = dq_min; }
    /**
     * Returns the minimum inter-quark distance.
     * \return The stored value [fm].
     */
    double getDqmin() const { return dq_min_; }
    /**
     * Sets \f$\mu_0\f$ in the running-coupling formula.
     * \param[in] x New value.
     */
    void setMuZero(double x) { muZero_ = x; }
    /**
     * Returns \f$\mu_0\f$ in the running-coupling formula.
     * \return The stored value.
     */
    double getMuZero() const { return muZero_; }
    /**
     * Sets the running-coupling cutoff smoothness parameter.
     * \param[in] x New value.
     */
    void setc(double x) { c_ = x; }
    /**
     * Returns the running-coupling cutoff smoothness parameter.
     * \return The stored value.
     */
    double getc() const { return c_; }
    /**
     * Sets the center-of-mass energy.
     * \param[in] x New value [GeV].
     */
    void setRoots(double x) { roots_ = x; }
    /**
     * Returns the center-of-mass energy.
     * \return The stored value [GeV].
     */
    double getRoots() const { return roots_; }
    /**
     * Sets whether Bjorken \f$x\f$ fluctuates with the local
     * \f$Q_s\f$.
     * \param[in] x New value.
     */
    void setUseFluctuatingx(int x) { useFluctuatingx_ = x; }
    /**
     * Returns whether Bjorken \f$x\f$ fluctuates with the local
     * \f$Q_s\f$.
     * \return The stored value.
     */
    int getUseFluctuatingx() const { return useFluctuatingx_; }
    /**
     * Sets the factor multiplying \f$Q_s\f$ under the running
     * \f$\alpha_s\f$ log.
     * \param[in] x New value.
     */
    void setRunWithThisFactorTimesQs(double x) {
        runWithThisFactorTimesQs_ = x;
    };
    /**
     * Returns the factor multiplying \f$Q_s\f$ under the running
     * \f$\alpha_s\f$ log.
     * \return The stored value.
     */
    double getRunWithThisFactorTimesQs() const {
        return runWithThisFactorTimesQs_;
    }
    /**
     * Sets the factor \f$\beta\f$ in \f$x = Q_s\beta/\sqrt{s}\f$.
     * \param[in] x New value.
     */
    void setxFromThisFactorTimesQs(double x) { xFromThisFactorTimesQs_ = x; };
    /**
     * Returns the factor \f$\beta\f$ in \f$x = Q_s\beta/\sqrt{s}\f$.
     * \return The stored value.
     */
    double getxFromThisFactorTimesQs() const { return xFromThisFactorTimesQs_; }
    /**
     * Sets the nuclear-thickness convolution \f$T_{pp}\f$.
     * \param[in] x New value.
     */
    void setTpp(double x) { Tpp_ = x; }
    /**
     * Returns the nuclear-thickness convolution \f$T_{pp}\f$.
     * \return The stored value.
     */
    double getTpp() const { return Tpp_; }
    /**
     * Sets the fixed-\f$N_{\text{part}}\f$ resampling threshold.
     * \param[in] x New value.
     */
    void setUseFixedNpart(int x) { useFixedNpart_ = x; }
    /**
     * Returns the fixed-\f$N_{\text{part}}\f$ resampling threshold.
     * \return The stored value.
     */
    int getUseFixedNpart() const { return useFixedNpart_; }
    /**
     * Sets the initial interaction region's area.
     * \param[in] x New value [fm\f$^2\f$].
     */
    void setArea(double x) { area_ = x; }
    /**
     * Returns the initial interaction region's area.
     * \return The stored value [fm\f$^2\f$].
     */
    double getArea() const { return area_; }
    /**
     * Sets the initial event-plane angle \f$\Psi_2\f$.
     * \param[in] x New value.
     */
    void setPsi(double x) { psi_ = x; }
    /**
     * Returns the initial event-plane angle \f$\Psi_2\f$.
     * \return The stored value.
     */
    double getPsi() const { return psi_; }
    /**
     * Sets the \f$Q_s\f$-smearing Gaussian width.
     * \param[in] x New value.
     */
    void setSmearingWidth(double x) { smearingWidth_ = x; }
    /**
     * Returns the \f$Q_s\f$-smearing Gaussian width.
     * \return The stored value.
     */
    double getSmearingWidth() const { return smearingWidth_; }
    /**
     * Sets this process' MPI rank.
     * \param[in] x New value.
     */
    void setMPIRank(int x) { MPIrank_ = x; }
    /**
     * Returns this process' MPI rank.
     * \return The stored value.
     */
    int getMPIRank() const { return MPIrank_; }
    /**
     * Sets the current event's identifier.
     * \param[in] x New value.
     */
    void setEventId(int x) { event_id_ = x; }
    /**
     * Returns the current event's identifier.
     * \return The stored value.
     */
    int getEventId() const { return event_id_; }
    /**
     * Sets the total number of MPI ranks.
     * \param[in] x New value.
     */
    void setMPISize(int x) { MPIsize_ = x; }
    /**
     * Returns the total number of MPI ranks.
     * \return The stored value.
     */
    int getMPISize() const { return MPIsize_; }
    /**
     * Sets whether a collision happened this event.
     * \param[in] x New value (`0` no collision, `1` collision).
     */
    void setSuccess(int x) { success_ = x; }
    /**
     * Returns whether a collision happened this event.
     * \return The stored value.
     */
    int getSuccess() const { return success_; }
    /**
     * Sets the per-nucleon thickness-distribution cutoff radius.
     * \param[in] x New value [fm].
     */
    void setRmax(double x) { rmax_ = x; }
    /**
     * Returns the per-nucleon thickness-distribution cutoff radius.
     * \return The stored value [fm].
     */
    double getRmax() const { return rmax_; }
    /**
     * Sets the UV damping parameter.
     * \param[in] x New value.
     */
    void setUVdamp(double x) { UVdamp_ = x; }
    /**
     * Returns the UV damping parameter.
     * \return The stored value.
     */
    double getUVdamp() const { return UVdamp_; }
    /**
     * Sets whether \c R_WS_/\c a_WS_/\c beta2_/\c beta3_/\c beta4_/\c
     * gamma_ override a species' built-in deformation parameters.
     * \param[in] x Non-zero to enable the override.
     */
    void setSetWSDeformParams(int x) {
        if (x == 0)
            setWSDeformParams_ = false;
        else
            setWSDeformParams_ = true;
    }
    /**
     * Returns whether the Woods-Saxon deformation-parameter override
     * is enabled.
     * \return The stored value.
     */
    bool getSetWSDeformParams() const { return setWSDeformParams_; }
    /**
     * Sets the Woods-Saxon half-density radius override.
     * \param[in] x New value [fm].
     */
    void setR_WS(double x) { R_WS_ = x; }
    /**
     * Returns the Woods-Saxon half-density radius override.
     * \return The stored value [fm].
     */
    double getR_WS() const { return R_WS_; }
    /**
     * Sets the Woods-Saxon surface diffuseness override.
     * \param[in] x New value [fm].
     */
    void setA_WS(double x) { a_WS_ = x; }
    /**
     * Returns the Woods-Saxon surface diffuseness override.
     * \return The stored value [fm].
     */
    double getA_WS() const { return a_WS_; }
    /**
     * Sets the quadrupole deformation override \f$\beta_2\f$.
     * \param[in] x New value.
     */
    void setBeta2(double x) { beta2_ = x; }
    /**
     * Returns the quadrupole deformation override \f$\beta_2\f$.
     * \return The stored value.
     */
    double getBeta2() const { return beta2_; }
    /**
     * Sets the octupole deformation override \f$\beta_3\f$.
     * \param[in] x New value.
     */
    void setBeta3(double x) { beta3_ = x; }
    /**
     * Returns the octupole deformation override \f$\beta_3\f$.
     * \return The stored value.
     */
    double getBeta3() const { return beta3_; }
    /**
     * Sets the hexadecapole deformation override \f$\beta_4\f$.
     * \param[in] x New value.
     */
    void setBeta4(double x) { beta4_ = x; }
    /**
     * Returns the hexadecapole deformation override \f$\beta_4\f$.
     * \return The stored value.
     */
    double getBeta4() const { return beta4_; }
    /**
     * Sets the triaxiality angle override \f$\gamma\f$.
     * \param[in] x New value [rad].
     */
    void setGamma(double x) { gamma_ = x; }
    /**
     * Returns the triaxiality angle override \f$\gamma\f$.
     * \return The stored value [rad].
     */
    double getGamma() const { return gamma_; }
    /**
     * Sets the minimum inter-nucleon distance.
     * \param[in] x New value [fm].
     */
    void setDmin(double x) { d_min_ = x; }
    /**
     * Returns the minimum inter-nucleon distance.
     * \return The stored value [fm].
     */
    double getDmin() const { return d_min_; }
    /**
     * Sets whether to enforce a minimum inter-nucleon distance.
     * \param[in] x Non-zero to enable enforcement.
     */
    void setForceDmin(int x) {
        if (x == 0)
            forceDminFlag_ = false;
        else
            forceDminFlag_ = true;
    }
    /**
     * Returns whether a minimum inter-nucleon distance is enforced.
     * \return The stored value.
     */
    bool getForceDmin() const { return (forceDminFlag_); }
    /**
     * Sets the neutron-skin radius offset.
     * \param[in] dR_np New value [fm].
     */
    void setWSdR_np(double dR_np) { WSdR_np_ = dR_np; }
    /**
     * Returns the neutron-skin radius offset.
     * \return The stored value [fm].
     */
    double getWSdR_np() const { return WSdR_np_; }
    /**
     * Sets the neutron-skin diffuseness offset.
     * \param[in] da_np New value [fm].
     */
    void setWSda_np(double da_np) { WSda_np_ = da_np; }
    /**
     * Returns the neutron-skin diffuseness offset.
     * \return The stored value [fm].
     */
    double getWSda_np() const { return WSda_np_; }

    /**
     * Sets whether to randomly rotate the reaction plane.
     * \param[in] iflag Non-zero to enable the rotation.
     */
    void setRotateReactionPlane(int iflag) {
        if (iflag == 0) {
            rotateReactionPlane_ = false;
        } else {
            rotateReactionPlane_ = true;
        }
    }
    /**
     * Returns whether the reaction plane is randomly rotated.
     * \return The stored value.
     */
    bool getRotateReactionPlane() const { return rotateReactionPlane_; }

    // switches:
    /**
     * Sets whether to use nuclei with finite geometry.
     * \param[in] x New value (`1` finite nucleus, `0` constant
     * \f$g^2\mu\f$).
     */
    void setUseNucleus(int x) { useNucleus_ = x; };
    /**
     * Returns whether nuclei with finite geometry are used.
     * \return The stored value.
     */
    int getUseNucleus() const { return useNucleus_; }
    /**
     * Sets whether to use a Gaussian profile on top of the constant
     * background.
     * \param[in] x New value.
     */
    void setUseGaussian(int x) { useGaussian_ = x; };
    /**
     * Returns whether a Gaussian profile on top of the constant
     * background is used.
     * \return The stored value.
     */
    int getUseGaussian() const { return useGaussian_; }
    /**
     * Sets the light-nucleus sampling method.
     * \param[in] x New value (`1` Woods-Saxon, `2` variational MC, `3`
     * alpha clusters).
     */
    void setlightNucleusOption(int x) { lightNucleusOption_ = x; };
    /**
     * Returns the light-nucleus sampling method.
     * \return The stored value.
     */
    int getlightNucleusOption() const { return lightNucleusOption_; }
    /**
     * Sets the projectile's polarization flag.
     * \param[in] x New value (`0` unpolarized, `1` longitudinal, `2`
     * transverse).
     */
    void setPolarizationProjectile(int x) { polarizationFlagProjectile_ = x; };
    /**
     * Returns the projectile's polarization flag.
     * \return The stored value.
     */
    int getPolarizationProjectile() const {
        return polarizationFlagProjectile_;
    }
    /**
     * Sets the target's polarization flag.
     * \param[in] x New value (`0` unpolarized, `1` longitudinal, `2`
     * transverse).
     */
    void setPolarizationTarget(int x) { polarizationFlagTarget_ = x; };
    /**
     * Returns the target's polarization flag.
     * \return The stored value.
     */
    int getPolarizationTarget() const { return polarizationFlagTarget_; }
    /**
     * Sets the projectile's \f$J_z\f$ polarization.
     * \param[in] x New value.
     */
    void setPolarizationProjectileJz(double x) { polJzProjectile_ = x; };
    /**
     * Returns the projectile's \f$J_z\f$ polarization.
     * \return The stored value.
     */
    double getPolarizationProjectileJz() const { return polJzProjectile_; }
    /**
     * Sets the target's \f$J_z\f$ polarization.
     * \param[in] x New value.
     */
    void setPolarizationTargetJz(double x) { polJzTarget_ = x; };
    /**
     * Returns the target's \f$J_z\f$ polarization.
     * \return The stored value.
     */
    double getPolarizationTargetJz() const { return polJzTarget_; }
    /**
     * Sets which \f$Q_s\f$ average \f$\alpha_s\f$ runs with.
     * \param[in] x New value (`0` min, `1` average, `2` max).
     */
    void setRunWithQs(int x) { runWithQs_ = x; };
    /**
     * Returns which \f$Q_s\f$ average \f$\alpha_s\f$ runs with.
     * \return The stored value.
     */
    int getRunWithQs() const { return runWithQs_; }
    /**
     * Sets whether \f$\alpha_s\f$ runs with \f$k_T\f$.
     * \param[in] x New value.
     */
    void setRunWithkt(int x) { runWithkt_ = x; };
    /**
     * Returns whether \f$\alpha_s\f$ runs with \f$k_T\f$.
     * \return The stored value.
     */
    int getRunWithkt() const { return runWithkt_; }
    /**
     * Sets whether \f$\alpha_s\f$ runs with the local \f$Q_s\f$.
     * \param[in] x New value.
     */
    void setRunWithLocalQs(int x) { runWithLocalQs_ = x; };
    /**
     * Returns whether \f$\alpha_s\f$ runs with the local \f$Q_s\f$.
     * \return The stored value.
     */
    int getRunWithLocalQs() const { return runWithLocalQs_; }
    /**
     * Sets whether to sample the impact parameter linearly.
     * \param[in] x New value.
     */
    void setLinearb(int x) { linearb_ = x; };
    /**
     * Returns whether the impact parameter is sampled linearly.
     * \return The stored value.
     */
    int getLinearb() const { return linearb_; }
    /**
     * Sets the output-file bitmask.
     * \param[in] x New value (see README.md's "writeOutputs" bit
     * values).
     */
    void setWriteOutputs(int x) { writeOutputs_ = x; };
    /**
     * Returns the output-file bitmask.
     * \return The stored value.
     */
    int getWriteOutputs() const { return writeOutputs_; }
    /**
     * Sets whether to run the flow-velocity/hydro-output calculation.
     * \param[in] x New value.
     */
    void setWriteEpsilonUHydro(int x) { writeEpsilonUHydro_ = x; };
    /**
     * Returns whether the flow-velocity/hydro-output calculation runs.
     * \return The stored value.
     */
    int getWriteEpsilonUHydro() const { return writeEpsilonUHydro_; }
    /**
     * Sets whether \f$T^{\mu\nu}\f$ is written as binary or text.
     * \param[in] x New value (`1` binary, `0` text).
     */
    void setWriteTmunuBinary(int x) { writeTmunuBinary_ = x; };
    /**
     * Returns whether \f$T^{\mu\nu}\f$ is written as binary or text.
     * \return The stored value.
     */
    int getWriteTmunuBinary() const { return writeTmunuBinary_; }
    /**
     * Sets whether to collect outputs into one HDF5 file.
     * \param[in] x New value.
     */
    void setWriteOutputsToHDF5(int x) { writeOutputsToHDF5_ = x; };
    /**
     * Returns whether outputs are collected into one HDF5 file.
     * \return The stored value.
     */
    int getWriteOutputsToHDF5() const { return writeOutputsToHDF5_; }
    /**
     * Sets the Wilson-line output format.
     * \param[in] x New value (`0` none, `1` text, `2` binary).
     */
    void setWriteWilsonLines(int x) { writeWilsonLines_ = x; }
    /**
     * Returns the Wilson-line output format.
     * \return The stored value.
     */
    int getWriteWilsonLines() const { return writeWilsonLines_; }
    /**
     * Sets the directory Wilson lines are written to or read from.
     * \param[in] x New path.
     */
    void setWilsonLinePath(std::string x) { wilsonLinePath_ = x; }
    /**
     * Returns the directory Wilson lines are written to or read from.
     * \return The stored path.
     */
    std::string getWilsonLinePath() const { return wilsonLinePath_; }
    /**
     * Sets whether/how to read initial Wilson lines from disk.
     * \param[in] x New value (`0` generate, `1` text, `2` binary).
     */
    void setReadInitialWilsonLines(int x) { readInitialWilsonLines_ = x; }
    /**
     * Returns whether/how initial Wilson lines are read from disk.
     * \return The stored value.
     */
    int getReadInitialWilsonLines() const { return readInitialWilsonLines_; }
    /**
     * Sets how nucleon positions are obtained.
     * \param[in] x New value (`0` sample, `1` file, `2` Alvioli
     * Pb-208 files).
     */
    void setNucleonPositionsFromFile(int x) { nucleonPositionsFromFile_ = x; }
    /**
     * Returns how nucleon positions are obtained.
     * \return The stored value.
     */
    int getNucleonPositionsFromFile() const {
        return nucleonPositionsFromFile_;
    }
    /**
     * Sets the nuclear configuration files' path.
     * \param[in] x New path.
     */
    void setNuclearConfigurationsPath(std::string x) {
        nuclearConfigurationsPath_ = x;
    }
    /**
     * Returns the nuclear configuration files' path.
     * \return The stored path.
     */
    std::string getNuclearConfigurationsPath() const {
        return nuclearConfigurationsPath_;
    }
    /**
     * Sets whether to use \f$1/Q_s\f$ as the maximal evolution time.
     * \param[in] x New value.
     */
    void setInverseQsForMaxTime(int x) { inverseQsForMaxTime_ = x; };
    /**
     * Returns whether \f$1/Q_s\f$ is used as the maximal evolution
     * time.
     * \return The stored value.
     */
    int getInverseQsForMaxTime() const { return inverseQsForMaxTime_; }
    /**
     * Sets whether to smear \f$Q_s\f$ with a Poisson distribution.
     * \param[in] x New value.
     */
    void setSmearQs(int x) { smearQs_ = x; }
    /**
     * Returns whether \f$Q_s\f$ is smeared with a Poisson
     * distribution.
     * \return The stored value.
     */
    int getSmearQs() const { return smearQs_; }
    /**
     * Sets whether to read the gluon spectrum from file.
     * \param[in] x New value.
     */
    void setReadMultFromFile(int x) { readMultFromFile_ = x; }
    /**
     * Returns whether the gluon spectrum is read from file.
     * \return The stored value.
     */
    int getReadMultFromFile() const { return readMultFromFile_; }
    /**
     * Sets whether wounding uses a Gaussian cross section.
     * \param[in] x New value (`0` hard sphere, `1` Gaussian).
     */
    void setGaussianWounding(int x) { gaussianWounding_ = x; }
    /**
     * Returns whether wounding uses a Gaussian cross section.
     * \return The stored value.
     */
    int getGaussianWounding() const { return gaussianWounding_; }
    /**
     * Sets whether \c rapidityA_/\c rapidityB_ hold pseudorapidity.
     * \param[in] x New value.
     */
    void setUsePseudoRapidity(int x) { usePseudoRapidity_ = x; }
    /**
     * Returns whether \c rapidityA_/\c rapidityB_ hold pseudorapidity.
     * \return The stored value.
     */
    int getUsePseudoRapidity() const { return usePseudoRapidity_; }
    /**
     * Sets the number of constituent quarks ("hot spots") per proton.
     * \param[in] x New value; `0` disables substructure.
     */
    void setUseConstituentQuarkProton(int x) { useConstituentQuarkProton_ = x; }
    /**
     * Returns the number of constituent quarks per proton.
     * \return The stored value.
     */
    int getUseConstituentQuarkProton() const {
        return useConstituentQuarkProton_;
    }
    /**
     * Sets the base number of constituent quarks.
     * \param[in] NqBase New value.
     */
    void setNqBase(double NqBase) { NqBase_ = NqBase; }
    /**
     * Returns the base number of constituent quarks.
     * \return The stored value.
     */
    double getNqBase() const { return NqBase_; }
    /**
     * Sets the constituent-quark-count fluctuation.
     * \param[in] NqFluc New value.
     */
    void setNqFluc(double NqFluc) { NqFluc_ = NqFluc; }
    /**
     * Returns the constituent-quark-count fluctuation.
     * \return The stored value.
     */
    double getNqFluc() const { return NqFluc_; }
    /**
     * Sets whether to use a smooth Woods-Saxon distribution for a
     * heavy nucleus.
     * \param[in] x New value.
     */
    void setUseSmoothNucleus(int x) { useSmoothNucleus_ = x; }
    /**
     * Returns whether a smooth Woods-Saxon distribution is used for a
     * heavy nucleus.
     * \return The stored value.
     */
    int getUseSmoothNucleus() const { return useSmoothNucleus_; }
    /**
     * Sets whether to shift the constituent-quark center of mass to
     * the origin.
     * \param[in] x New value.
     */
    void setShiftConstituentQuarkProtonOrigin(int x) {
        shiftConstituentQuarkProtonOrigin_ = x;
    }
    /**
     * Returns whether the constituent-quark center of mass is shifted
     * to the origin.
     * \return The stored value.
     */
    int getShiftConstituentQuarkProtonOrigin() const {
        return shiftConstituentQuarkProtonOrigin_;
    }
    /**
     * Sets the minimum \f$Q_{s,\min}^2 S_T\f$ trigger threshold.
     * \param[in] x New value; `0` disables the trigger.
     */
    void setMinimumQs2ST(int x) { minimumQs2ST_ = x; }
    /**
     * Returns the minimum \f$Q_{s,\min}^2 S_T\f$ trigger threshold.
     * \return The stored value.
     */
    int getMinimumQs2ST() const { return minimumQs2ST_; }

    /**
     * Sets whether to compute the gluon multiplicity spectrum.
     * \param[in] x Non-zero to enable.
     */
    void setComputeGluonMultiplicity(int x) {
        if (x == 0) {
            computeGluonMultiplicity_ = false;
        } else {
            computeGluonMultiplicity_ = true;
        }
    }
    /**
     * Returns whether the gluon multiplicity spectrum is computed.
     * \return The stored value.
     */
    bool getComputeGluonMultiplicity() const {
        return computeGluonMultiplicity_;
    }

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
     * instance's \c m_/\c BG_/\c BGq_/\c smearingWidth_/\c NqBase_/\c
     * QsmuRatio_/\c dq_min_.
     * \param[in] itype `1` uses \c posteriorParamSets_ (also sets \c
     * NqBase_ from the table); `2` or `4` use \c
     * posteriorParamSetsNq3_ (fixing \c NqBase_ to `3`); any other
     * value is a no-op.
     * \param[in] iset Row index, taken modulo the selected table's row
     * count.
     */
    void setParamsWithPosteriorParameterSet(const int itype, int iset);

    // JIMWLK functions
    /**
     * Sets JIMWLK's infrared regulator.
     * \param[in] x New value [GeV].
     */
    void setm_jimwlk(double x) { m_jimwlk_ = x; };
    /**
     * Returns JIMWLK's infrared regulator.
     * \return The stored value [GeV].
     */
    double getm_jimwlk() const { return m_jimwlk_; }
    /**
     * Sets JIMWLK's running-\f$\alpha_s\f$ regulator \f$\mu_0\f$.
     * \param[in] x New value [GeV].
     */
    void setMu0_jimwlk(double x) { mu0_jimwlk_ = x; }
    /**
     * Returns JIMWLK's running-\f$\alpha_s\f$ regulator \f$\mu_0\f$.
     * \return The stored value [GeV].
     */
    double getMu0_jimwlk() const { return mu0_jimwlk_; }
    /**
     * Sets whether JIMWLK uses the simple Langevin discretization.
     * \param[in] x New value; `0` disables it.
     */
    void setSimpleLangevin(int x) {
        if (x == 0) {
            simpleLangevin_ = false;
        } else {
            simpleLangevin_ = x;
        }
    }
    /**
     * Returns whether JIMWLK uses the simple Langevin discretization.
     * \return The stored value.
     */
    bool getSimpleLangevin() const { return simpleLangevin_; }
    /**
     * Sets JIMWLK's \f$\Lambda_{QCD}\f$.
     * \param[in] x New value [GeV].
     */
    void setLambdaQCD_jimwlk(double x) { LambdaQCD_jimwlk_ = x; }
    /**
     * Returns JIMWLK's \f$\Lambda_{QCD}\f$.
     * \return The stored value [GeV].
     */
    double getLambdaQCD_jimwlk() const { return LambdaQCD_jimwlk_; }
    /**
     * Sets the Bjorken \f$x\f$ the projectile is evolved to.
     * \param[in] x New value.
     */
    void setJimwlk_x_projectile(double x) { jimwlk_x1_ = x; }
    /**
     * Returns the Bjorken \f$x\f$ the projectile is evolved to.
     * \return The stored value.
     */
    double getJimwlk_x_projectile() const { return jimwlk_x1_; }
    /**
     * Sets the Bjorken \f$x\f$ the target is evolved to.
     * \param[in] x New value.
     */
    void setJimwlk_x_target(double x) { jimwlk_x2_ = x; }
    /**
     * Returns the Bjorken \f$x\f$ the target is evolved to.
     * \return The stored value.
     */
    double getJimwlk_x_target() const { return jimwlk_x2_; }
    /**
     * Sets the JIMWLK evolution step size.
     * \param[in] x New value.
     */
    void setDs_jimwlk(double x) { ds_jimwlk_ = x; }
    /**
     * Returns the JIMWLK evolution step size.
     * \return The stored value.
     */
    double getDs_jimwlk() const { return ds_jimwlk_; }
    /**
     * Sets JIMWLK's coupling.
     * \param[in] as New value; `0` selects running coupling.
     */
    void setJimwlk_alphas(double as) { jimwlk_alphas_ = as; }
    /**
     * Returns JIMWLK's coupling.
     * \return The stored value.
     */
    double getJimwlk_alphas() const { return jimwlk_alphas_; }
    /**
     * Sets JIMWLK's initial Bjorken \f$x\f$.
     * \param[in] x New value.
     */
    void setJimwlk_x0(double x) { x0_jimwlk_ = x; }
    /**
     * Returns JIMWLK's initial Bjorken \f$x\f$.
     * \return The stored value.
     */
    double getJimwlk_x0() const { return x0_jimwlk_; }
    /**
     * Returns whether JIMWLK evolution is enabled.
     * \return The stored value.
     */
    bool getUseJIMWLK() const { return useJIMWLK_; }
    /**
     * Sets whether JIMWLK evolution is enabled.
     * \param[in] x Non-zero to enable.
     */
    void setUseJIMWLK(int x) {
        if (x == 0) {
            useJIMWLK_ = false;
        } else {
            useJIMWLK_ = true;
        }
    }
    /**
     * Sets the JIMWLK snapshot \f$x\f$ list.
     * \param[in] xList New list of Bjorken \f$x\f$ values.
     */
    void setxSnapshotList(std::vector<double> xList) { xSnapshotList_ = xList; }
    /**
     * Returns the JIMWLK snapshot \f$x\f$ list.
     * \return The stored list.
     */
    std::vector<double> getxSnapshotList() const { return xSnapshotList_; }
    /**
     * Sets whether to save Wilson-line snapshots during JIMWLK
     * evolution.
     * \param[in] x Non-zero to enable.
     */
    void setSaveSnapshots(int x) {
        if (x == 0)
            saveSnapshots_ = false;
        else
            saveSnapshots_ = true;
    }
    /**
     * Returns whether Wilson-line snapshots are saved during JIMWLK
     * evolution.
     * \return The stored value.
     */
    bool getSaveSnapshots() const { return saveSnapshots_; }

    /**
     * Checks whether the current parameter set is internally
     * consistent enough to run (e.g. a positive lattice size, a valid
     * Wilson-line data format, snapshots only requested when Wilson
     * lines are actually written), logging an error and returning
     * `false` on the first problem found.
     * \return `true` if every check passes, `false` otherwise.
     */
    bool ValidParameters();
};
#endif  // SRC_PARAMETERS_H_

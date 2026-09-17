#ifndef SRC_GLAUBER_H_
#define SRC_GLAUBER_H_

#include <string>
#include <vector>

#include "PrettyOstream.h"
#include "Random.h"

/// Adaptive-quadrature convergence tolerance used by integral()/qnc7().
#define TOL (1.0e-6)
/// Small-number regularizer used throughout the density-profile
/// integrands to keep \f$\xi=0\f$/\f$\xi=1\f$ away from a
/// \f$\log(0)\f$ or division-by-zero singularity.
#define TINY (1.0e-10)
/// Maximum number of integrand evaluations qnc7() will perform before
/// giving up on reaching \c TOL and returning its best estimate so far.
#define LIMIT 10000

/**
 * Which of the two colliding nuclei a lattice/nucleon-sampling operation
 * applies to.
 */
enum class NucleusRole {
    /// The first (\f$A\f$) nucleus.
    Projectile,
    /// The second (\f$B\f$) nucleus.
    Target,
};

/**
 * Selects which integrand evaluateIntegrand()/integral()/qnc7() sample.
 *
 * `NuInt2HO`/`NuInt3Gauss`/`NuInt3Fermi`/`NuIntHulthen` mirror
 * `Nucleus::densityFunc`'s values (1/2/3/8) and are used for nuInS()'s
 * normalization integral, whose integrand depends on the active
 * nucleus's density profile; the other four are fixed,
 * density-profile-specific integrands used by exactly one caller each
 * (anum3Fermi()/anum3Gauss()/anum2HO()'s own normalization, and tAB()'s
 * overlap integral) and never vary at runtime.
 */
enum class IntegrandId {
    /// nuInt2HO(): 2-parameter harmonic-oscillator density integrand.
    NuInt2HO = 1,
    /// nuInt3Gauss(): 3-parameter Gaussian density integrand.
    NuInt3Gauss = 2,
    /// nuInt3Fermi(): 3-parameter Fermi (Woods-Saxon) density
    /// integrand.
    NuInt3Fermi = 3,
    /// anum3FermiInt(): 3-parameter Fermi normalization integrand.
    Anum3FermiInt = 4,
    /// anum3GaussInt(): 3-parameter Gaussian normalization integrand.
    Anum3GaussInt = 5,
    /// anum2HOInt(): 2-parameter harmonic-oscillator normalization
    /// integrand.
    Anum2HOInt = 6,
    /// oLSIntegrand(): nuclear overlap function \f$T_{AB}\f$ integrand.
    OLSIntegrand = 7,
    /// nuIntHulthen(): Hulthen (deuteron) density integrand.
    NuIntHulthen = 8,
};

/**
 * One sampled nucleon's transverse position, longitudinal angle, and
 * species/collision bookkeeping.
 */
struct ReturnValue {
    /// Transverse \f$x\f$ position [fm], relative to the nucleus
    /// center.
    double x;
    /// Transverse \f$y\f$ position [fm], relative to the nucleus
    /// center.
    double y;
    /// Longitudinal \f$z\f$ position [fm], relative to the nucleus
    /// center (0 for a purely transverse sampling such as
    /// Glauber::sampleTARejection()).
    double z;
    /// Anisotropy/orientation angle [rad]; not set by every producer
    /// (see Glauber::sampleTARejection()) -- callers that need it set
    /// for every nucleon do so afterward (e.g.
    /// Init::sampleNucleonAnisotropyAngles()).
    double phi;
    /// `1` if this nucleon has undergone at least one binary collision,
    /// `0` otherwise. Initialized to `0` by every producer (both
    /// Glauber::sampleTARejection() and Init.cpp's other position
    /// samplers); later set to `1` by Init::computeNcollList() for
    /// nucleons found to collide, and reset to `0` again at the start
    /// of each new event by Init::sampleImpactParameter().
    int collided;
    /// `true` if this nucleon is a proton, `false` if a neutron.
    bool proton;
};

/**
 * One nucleus species' Woods-Saxon/density-profile parameters, as
 * resolved by Glauber::findNucleusData().
 */
struct Nucleus {
    /// Species name (e.g. `"Pb"`, `"p"`, `"d"`), as passed to
    /// Glauber::findNucleusData().
    std::string name;
    /// Mass number (nucleon count).
    int A;
    /// Atomic number (proton count).
    int Z;
    /// Density-function selector used by Glauber::calcRho() to pick
    /// which \c anumX() normalization to call; same encoding as \c
    /// densityFunc (1=2HO, 2=3Gauss, 3=3Fermi, 8=Hulthen).
    int anumFunc;
    /// Set by Glauber::findNucleusData() to the same value as \c
    /// anumFunc, but never read anywhere in the codebase.
    int anumFuncIntegrand;
    /// Density-function selector used by Glauber::nuInS() (via
    /// IntegrandId) to pick which \c nuIntX() integrand to sample; same
    /// encoding as \c anumFunc (1=2HO, 2=3Gauss, 3=3Fermi, 8=Hulthen).
    int densityFunc;
    /// Woods-Saxon/profile-shape parameter \f$w\f$ (meaning depends on
    /// \c densityFunc: e.g. the deformation-adjacent weight for 3Fermi/
    /// 3Gauss, or the second Hulthen decay constant \f$b\f$)
    /// [dimensionless or 1/fm, profile-dependent].
    double w_WS;
    /// Woods-Saxon surface diffuseness \f$a\f$ [fm].
    double a_WS;
    /// Woods-Saxon half-density radius \f$R\f$ [fm].
    double R_WS;
    /// Central density, normalized by Glauber::calcRho() so the
    /// profile integrates to \c A [1/fm^3].
    double rho_WS;
    /// Quadrupole deformation parameter \f$\beta_2\f$
    /// [dimensionless].
    double beta2;
    /// Octupole deformation parameter \f$\beta_3\f$ [dimensionless].
    double beta3;
    /// Hexadecapole deformation parameter \f$\beta_4\f$
    /// [dimensionless].
    double beta4;
    /// Triaxiality angle \f$\gamma\f$ [rad].
    double gamma;
    /// Whether a minimum inter-nucleon distance (\c d_min) is enforced
    /// when sampling this nucleus' nucleon positions.
    bool forceDminFlag;
    /// Minimum inter-nucleon distance enforced when \c forceDminFlag is
    /// set [fm].
    double d_min;
    /// Neutron-skin radius offset \f$\Delta R_{np}\f$ applied between
    /// proton and neutron density profiles [fm].
    double dR_np;
    /// Neutron-skin diffuseness offset \f$\Delta a_{np}\f$ applied
    /// between proton and neutron density profiles [fm].
    double da_np;
};

/**
 * Collision-geometry configuration for one projectile/target pair, as
 * assembled by Glauber::initGlauber().
 */
struct Data {
    /// Inelastic nucleon-nucleon cross section \f$\sigma_{NN}\f$ [fm^2]
    /// (converted from the input file's mb via Glauber::initGlauber()).
    double sigmaNN;
    /// Target nucleus' resolved profile parameters.
    Nucleus target;
    /// Projectile nucleus' resolved profile parameters.
    Nucleus projectile;
    /// Upper integration cutoff [fm] for Glauber::tAB()'s overlap
    /// integral and the interpolation tables built by
    /// Glauber::interNuPInSP()/interNuTInST().
    double sCutoff;
    /// Number of grid points used by Glauber::interNuPInSP()/
    /// interNuTInST()'s interpolation tables; set from
    /// Glauber::initGlauber()'s \c imax. Unrelated to qnc7(), whose
    /// evaluation-count bound is the fixed \c LIMIT macro.
    int interMax;
};

/**
 * Glauber-model nuclear geometry: resolves a nucleus name to its
 * Woods-Saxon/density-profile parameters, provides the corresponding
 * thickness and overlap functions (\f$T_A\f$, \f$T_{AB}\f$), and
 * samples individual nucleon transverse positions from them.
 *
 * One Glauber instance is constructed fresh per event in main.cpp and
 * configured once via initGlauber() (which resolves both the
 * projectile and target species via findNucleusData()); every other
 * method then operates on that instance's resolved Data.
 */
class Glauber {
  private:
    /// Log sink for progress/warning/error messages.
    PrettyOstream messager_;
    /// Scratch state passed to the currently-active \c anumXInt()
    /// integrand via calcRho()/anum3Fermi()/anum3Gauss()/anum2HO(): the
    /// Woods-Saxon radius in units of \c a_WS.
    double AnumR_;
    /// Scratch state passed to the currently-active \c nuIntX()
    /// integrand via nuInS(): the transverse offset \f$s\f$ [fm] at
    /// which the thickness function is being evaluated.
    double NuInS_S_;
    /// The nucleus currently being integrated over by calcRho() and the
    /// \c anumX()/\c nuIntX() family (whichever of \c
    /// glauberData_.projectile/target was passed to calcRho() last).
    Nucleus *Nuc_WS_;
    /// This instance's resolved projectile/target configuration, set by
    /// initGlauber().
    Data glauberData_;
    /// Impact parameter \f$b\f$ [fm], set by initGlauber(); used by
    /// oLSIntegrand()'s overlap-function angular average.
    double b_;  // impact parameter
    /// Projectile mass number, cached from \c glauberData_.projectile.A
    /// by initGlauber() for nucleusA1().
    int currentA1_;
    /// Target mass number, cached from \c glauberData_.target.A by
    /// initGlauber() for nucleusA2().
    int currentA2_;
    /// Projectile atomic number, cached from \c
    /// glauberData_.projectile.Z by initGlauber() for nucleusZ1().
    int currentZ1_;
    /// Target atomic number, cached from \c glauberData_.target.Z by
    /// initGlauber() for nucleusZ2().
    int currentZ2_;

  public:
    /**
     * Constructs a Glauber with no configuration; call initGlauber()
     * before using any other method.
     */
    Glauber() {};
    /**
     * Removes the scratch file `tmp.dat` some callers may have written
     * (harmless if it doesn't exist).
     */
    ~Glauber() { remove("tmp.dat"); }

    /**
     * Returns the projectile's mass number.
     * \return \f$A\f$ of the projectile nucleus.
     */
    int nucleusA1() const { return currentA1_; }
    /**
     * Returns the target's mass number.
     * \return \f$A\f$ of the target nucleus.
     */
    int nucleusA2() const { return currentA2_; }
    /**
     * Returns the projectile's atomic number.
     * \return \f$Z\f$ of the projectile nucleus.
     */
    int nucleusZ1() const { return currentZ1_; }
    /**
     * Returns the target's atomic number.
     * \return \f$Z\f$ of the target nucleus.
     */
    int nucleusZ2() const { return currentZ2_; }
    /**
     * Returns this instance's resolved projectile/target configuration.
     * \return Const reference to the Data set up by initGlauber().
     */
    const Data &getGlauberData() const { return glauberData_; }
    /**
     * Resolves a nucleus name to its Woods-Saxon/density-profile
     * parameters (mass/atomic number, radius, diffuseness, deformation,
     * density-function selector) from the built-in table of supported
     * species, applying any explicit overrides requested.
     * \param[out] nucleus Filled with the resolved parameters.
     * \param[in] name Species name (e.g. `"Pb"`, `"Au"`, `"p"`, `"d"`);
     * exits with an error if not one of the supported species (see
     * Glauber.cpp's \c kNucleusTemplates table).
     * \param[in] setWSDeformParams If `true`, `R_WS`/`a_WS`/\p
     * beta2/`beta3`/`beta4`/`gamma` override the species' built-in
     * values instead of being ignored.
     * \param[in] R_WS Woods-Saxon half-density radius override [fm],
     * used only if \p setWSDeformParams.
     * \param[in] a_WS Woods-Saxon surface diffuseness override [fm],
     * used only if \p setWSDeformParams.
     * \param[in] beta2 Quadrupole deformation override
     * [dimensionless], used only if \p setWSDeformParams (also used
     * unconditionally for the handful of species whose built-in \c
     * beta2 is a sentinel meaning "use the caller's value").
     * \param[in] beta3 Octupole deformation override [dimensionless],
     * used only if \p setWSDeformParams.
     * \param[in] beta4 Hexadecapole deformation override
     * [dimensionless], used only if \p setWSDeformParams.
     * \param[in] gamma Triaxiality angle override [rad], used only if
     * \p setWSDeformParams.
     * \param[in] forceDminFlag Whether to enforce a minimum
     * inter-nucleon distance when sampling this nucleus (applied
     * unconditionally, regardless of \p setWSDeformParams).
     * \param[in] d_min Minimum inter-nucleon distance [fm], used only
     * if \p forceDminFlag.
     * \param[in] dR_np Neutron-skin radius offset [fm] (applied
     * unconditionally).
     * \param[in] da_np Neutron-skin diffuseness offset [fm] (applied
     * unconditionally).
     */
    void findNucleusData(
        Nucleus *nucleus, std::string name, bool setWSDeformParams, double R_WS,
        double a_WS, double beta2, double beta3, double beta4, double gamma,
        bool forceDminFlag, double d_min, double dR_np, double da_np);
    /**
     * Logs this instance's \c sigmaNN/interMax/sCutoff at debug level.
     *
     * \note No call site currently exists anywhere in this codebase.
     */
    void printGlauberData();
    /**
     * Logs one nucleus' name/A/Z/w_WS/a_WS/R_WS at info level.
     * \param[in] nucleus Nucleus whose parameters to log.
     *
     * \note No call site currently exists anywhere in this codebase.
     */
    void printNucleusData(Nucleus *nucleus);
    /**
     * Finds the index of the first of four evenly-spaced tabulated
     * points to use for a cubic interpolation at \p x, clamped so all
     * four points stay in range.
     * \param[in] x Point to interpolate at.
     * \param[in] Vx Evenly-spaced table abscissas, length \p ymax + 1.
     * \param[in] ymax Index of \p Vx's last entry.
     * \return Index of the first of the four interpolation points, in
     * `[0, ymax-3]`.
     */
    int linearFindXorg(double x, double *Vx, int ymax);
    /**
     * Cubic-polynomial interpolation of a tabulated function at \p x,
     * using the four points starting at \p x_org.
     * \param[in] x Point to interpolate at.
     * \param[in] Vx Evenly-spaced table abscissas.
     * \param[in] Vy Table ordinates, same length as \p Vx.
     * \param[in] h Spacing between consecutive \p Vx entries.
     * \param[in] x_org Index of the first of the four interpolation
     * points (see linearFindXorg()).
     * \return The interpolated value.
     */
    double fourPtInterpolate(
        double x, double *Vx, double *Vy, double h, int x_org);
    /**
     * Computes the cubic polynomial coefficients (in powers of
     * \f$x - Vx[\text{x\_org}]\f$) that pass through the four
     * consecutive table points starting at \p x_org.
     * \param[out] a Cubic coefficient.
     * \param[out] b Quadratic coefficient.
     * \param[out] c Linear coefficient.
     * \param[out] d Constant coefficient (equals `Vy[x_org]`).
     * \param[in] Vy Table ordinates.
     * \param[in] h Spacing between consecutive table abscissas.
     * \param[in] x_org Index of the first of the four points.
     */
    void makeCoeff(
        double *a, double *b, double *c, double *d, double *Vy, double h,
        int x_org);
    /**
     * Interpolates a tabulated function at \p x via linearFindXorg() +
     * fourPtInterpolate(), assuming \p Vx is evenly spaced.
     * \param[in] x Point to interpolate at; exits with an error if
     * outside `[Vx[0], Vx[ymax]]`.
     * \param[in] Vx Evenly-spaced table abscissas.
     * \param[in] Vy Table ordinates, same length as \p Vx.
     * \param[in] ymax Index of \p Vx's last entry.
     * \return The interpolated value.
     */
    double vInterpolate(double x, double *Vx, double *Vy, int ymax);
    /**
     * Builds an evenly-spaced table of abscissas.
     * \param[in] down First abscissa.
     * \param[in] up Last abscissa.
     * \param[in] maxi_num Number of intervals; the table has \p
     * maxi_num + 1 entries.
     * \return Table `{down, down+dx, ..., up}` with `dx = (up-down)/`
     * \p maxi_num.
     */
    std::vector<double> makeVx(double down, double up, int maxi_num);
    /**
     * Builds the corresponding table of ordinates for makeVx()'s
     * abscissas, by evaluating nuInS() (the currently-configured
     * nucleus' thickness-like normalization integral) at each one.
     * \param[in] vx Abscissa table from makeVx().
     * \param[in] maxi_num Number of intervals; \p vx and the returned
     * table both have \p maxi_num + 1 entries.
     * \return Table of `nuInS(vx[i])` for each `i`.
     */
    std::vector<double> makeVy(const double *vx, int maxi_num);
    /**
     * Reads a tabulated \f$V_x\f$ array from a text file (format:
     * whitespace-separated tokens up to and including an `EndOfData`
     * marker, then \p maxi_num + 1 pairs of numbers, keeping only the
     * first of each pair).
     * \param[in] file_name Path to the file; exits with an error if it
     * can't be opened, is empty/malformed, is missing its `EndOfData`
     * marker, or has fewer than \p maxi_num + 1 entries.
     * \param[in] maxi_num Number of entries to read minus one; the
     * returned table has \p maxi_num + 1 entries.
     * \param[in] quiet If `1`, logs a "reading in" message before
     * opening the file.
     * \return The table's first-of-each-pair values.
     *
     * \note No call site currently exists anywhere in this codebase.
     */
    std::vector<double> readInVx(char *file_name, int maxi_num, int quiet);
    /**
     * Same file format as readInVx(), but keeps the second of each
     * pair instead of the first.
     * \param[in] file_name Path to the file; exits with an error under
     * the same conditions as readInVx().
     * \param[in] maxi_num Number of entries to read minus one; the
     * returned table has \p maxi_num + 1 entries.
     * \param[in] quiet If `1`, logs a "reading in" message before
     * opening the file.
     * \return The table's second-of-each-pair values.
     *
     * \note No call site currently exists anywhere in this codebase.
     */
    std::vector<double> readInVy(char *file_name, int maxi_num, int quiet);

    /**
     * Projectile thickness function \f$T_A(s)\f$ at transverse offset
     * \p s, memoized in a lookup table built on the first call.
     *
     * \warning The lookup table is a function-static local, shared by
     * every call to this method regardless of which Glauber instance
     * makes it: it is built once, from whichever projectile is
     * configured on the *first* call made anywhere in the process, and
     * never rebuilt afterward. This is safe today because a run's
     * projectile species and Woods-Saxon parameters never change after
     * the first event, but would silently return stale results if that
     * ever changed.
     * \param[in] s Transverse offset [fm]; returns `0.0` for a
     * single-nucleon (\f$A=1\f$) projectile, or for \p s beyond the
     * table's upper bound (\f$2 \times\f$ \c sCutoff).
     * \return \f$T_A(s)\f$ [1/fm^2] (clamped to `0.0` if the
     * interpolation would go negative).
     */
    double interNuPInSP(double s);
    /**
     * Target thickness function \f$T_B(s)\f$ at transverse offset \p
     * s; see interNuPInSP() (same caching caveat, same table
     * construction, mirrored for the target nucleus).
     * \param[in] s Transverse offset [fm].
     * \return \f$T_B(s)\f$ [1/fm^2].
     */
    double interNuTInST(double s);
    /**
     * Normalizes a nucleus' central density \c rho_WS so its density
     * profile integrates to \p nucleus's mass number \c A, and records
     * \p nucleus as the one the \c anumX()/\c nuIntX() integrand family
     * (and hence nuInS()) currently operates on.
     * \param[in,out] nucleus Nucleus whose \c rho_WS is normalized in
     * place, using its \c anumFunc to pick which \c anumX() to call
     * (defaulting to \c anum3Fermi() for an unrecognized value).
     */
    void calcRho(Nucleus *nucleus);
    /**
     * Thickness-like normalization integral
     * \f$\int_0^\infty dz\, \rho(\sqrt{s^2+z^2})\f$ for the
     * currently-\link calcRho() configured\endlink nucleus at
     * transverse offset \p s, used to build the interNuPInSP()/
     * interNuTInST() lookup tables.
     * \param[in] s Transverse offset [fm].
     * \return The integral's value, via integral() dispatching on the
     * nucleus' \c densityFunc.
     */
    double nuInS(double s);

    /**
     * Normalization integral \f$\int d^3r\,\rho_{\text{3Fermi}}(r)\f$
     * for the 3-parameter Fermi (Woods-Saxon) profile, used by
     * calcRho() to fix \c rho_WS.
     * \param[in] R_WS Woods-Saxon half-density radius [fm].
     * \return The (unnormalized) integral's value.
     */
    double anum3Fermi(double R_WS);
    /**
     * Integrand for anum3Fermi(), in the substitution
     * \f$\xi = e^{-r/a}\f$.
     * \param[in] xi Integration variable in `(0, 1)`, clamped away from
     * the endpoints by \c TINY.
     * \return The integrand's value at \p xi.
     */
    double anum3FermiInt(double xi);
    /**
     * 3-parameter Fermi (Woods-Saxon) density integrand for nuInS(), in
     * the substitution \f$\xi = e^{-z/a}\f$ at fixed transverse offset
     * \link NuInS_S_ s\endlink.
     * \param[in] xi Integration variable in `(0, 1)`, clamped away from
     * the endpoints by \c TINY.
     * \return The integrand's value at \p xi.
     */
    double nuInt3Fermi(double xi);
    /**
     * Normalization integral \f$\int d^3r\,\rho_{\text{3Gauss}}(r)\f$
     * for the 3-parameter Gaussian profile, used by calcRho() to fix
     * \c rho_WS.
     * \param[in] R_WS Gaussian-profile radius parameter [fm].
     * \return The (unnormalized) integral's value.
     */
    double anum3Gauss(double R_WS);
    /**
     * Integrand for anum3Gauss(), in the substitution
     * \f$\xi = e^{-r^2/a^2}\f$.
     * \param[in] xi Integration variable in `(0, 1)`, clamped away from
     * the endpoints by \c TINY.
     * \return The integrand's value at \p xi.
     */
    double anum3GaussInt(double xi);
    /**
     * 3-parameter Gaussian density integrand for nuInS(), in the
     * substitution \f$\xi = e^{-z^2/a^2}\f$ at fixed transverse offset
     * \link NuInS_S_ s\endlink.
     * \param[in] xi Integration variable in `(0, 1)`, clamped away from
     * the endpoints by \c TINY.
     * \return The integrand's value at \p xi.
     */
    double nuInt3Gauss(double xi);
    /**
     * Normalization integral \f$\int d^3r\,\rho_{\text{2HO}}(r)\f$ for
     * the 2-parameter harmonic-oscillator profile, used by calcRho()
     * to fix \c rho_WS.
     * \return The (unnormalized) integral's value.
     */
    double anum2HO();
    /**
     * Integrand for anum2HO(), in the substitution
     * \f$\xi = e^{-r^2}\f$.
     * \param[in] xi Integration variable in `(0, 1)`, clamped away from
     * the endpoints by \c TINY.
     * \return The integrand's value at \p xi.
     */
    double anum2HOInt(double xi);
    /**
     * 2-parameter harmonic-oscillator density integrand for nuInS(), in
     * the substitution \f$\xi = e^{-z^2}\f$ at fixed transverse offset
     * \link NuInS_S_ s\endlink.
     * \param[in] xi Integration variable in `(0, 1)`, clamped away from
     * the endpoints by \c TINY.
     * \return The integrand's value at \p xi.
     */
    double nuInt2HO(double xi);
    /**
     * Normalization integral for the Hulthen (deuteron) profile
     * \f$\rho \propto (e^{-ar} - e^{-br})^2/r^2\f$, evaluated in
     * closed form; used by calcRho() to fix \c rho_WS.
     * \return The (unnormalized) integral's value, times \c rho_WS
     * (calcRho() divides it back out).
     */
    double anumHulthen();
    /**
     * Hulthen (deuteron) density integrand for nuInS(), in the
     * substitution \f$\xi = e^{-za}\f$ at fixed transverse offset
     * \link NuInS_S_ s\endlink.
     * \param[in] xi Integration variable in `(0, 1)`, clamped away from
     * the endpoints by \c TINY.
     * \return The integrand's value at \p xi.
     */
    double nuIntHulthen(double xi);

    /**
     * Adaptive Newton-Cotes (7-point) quadrature of the density
     * profile integrand selected by \p id, over `[down, up]`.
     * \param[in] id Which integrand to sample (see evaluateIntegrand()).
     * \param[in] down Lower integration bound.
     * \param[in] up Upper integration bound.
     * \param[in] tol Relative convergence tolerance for qnc7()'s
     * adaptive subdivision.
     * \param[out] count Number of integrand evaluations used; capped at
     * \c LIMIT by qnc7().
     * \return The integral's estimated value (`0.0` if `down == up`).
     */
    double integral(
        IntegrandId id, double down, double up, double tol, int *count);
    /**
     * Recursive adaptive step of the 7-point Newton-Cotes quadrature
     * used by integral(): estimates the integral over the left and
     * right halves of `[down, down+12*dx]` and recurses into whichever
     * half hasn't converged to \p tol yet.
     * \param[in] id Which integrand to sample (see evaluateIntegrand()).
     * \param[in] tol Relative convergence tolerance for this level's
     * halves.
     * \param[in] down Lower bound of this call's interval.
     * \param[in] dx One-sixth of this call's interval width.
     * \param[in] f_of Integrand values at the 7 equally-spaced points
     * across this call's interval.
     * \param[in] pre_sum This interval's estimate from the parent call,
     * used to judge convergence.
     * \param[in] area Running estimate of the total absolute area
     * integrated so far, used to scale the convergence test.
     * \param[in,out] count Number of integrand evaluations used so
     * far; recursion stops once this reaches \c LIMIT regardless of
     * convergence.
     * \return This interval's integral estimate.
     */
    double qnc7(
        IntegrandId id, double tol, double down, double dx, double *f_of,
        double pre_sum, double area, int *count);
    /**
     * Dispatches to the density-profile integrand selected by \p id.
     * \param[in] id Which integrand to evaluate.
     * \param[in] xi Integration variable (meaning depends on \p id; see
     * the individual \c nuIntX()/\c anumXInt()/oLSIntegrand() methods).
     * \return The selected integrand's value at \p xi.
     */
    double evaluateIntegrand(IntegrandId id, double xi);
    /**
     * Integrand for tAB()'s nuclear overlap function: the projectile
     * thickness at transverse offset \p s, times a 20-point angular
     * average of the target thickness at
     * \f$\sqrt{s^2+b^2+2sb\cos\theta}\f$ (\f$b\f$ = \link b_
     * impact parameter\endlink).
     * \param[in] s Transverse offset [fm].
     * \return The overlap integrand's value at \p s.
     */
    double oLSIntegrand(double s);
    /**
     * Nuclear overlap function \f$T_{AB}(b)\f$ at this instance's
     * configured impact parameter, integrated via oLSIntegrand().
     * \return \f$T_{AB}(b)\f$ [1/fm^2] (dimensionless number of binary
     * collisions per unit \f$\sigma_{NN}\f$).
     */
    double tAB();
    /**
     * Configures this instance for one projectile/target pair: resolves
     * both species via findNucleusData() and records the collision
     * geometry (cross section, impact parameter, interpolation-table
     * resolution) in \c glauberData_. Must be called before any other
     * method.
     * \param[in] sigmaNN Inelastic nucleon-nucleon cross section [mb];
     * converted to fm^2 internally.
     * \param[in] target Target species name (see findNucleusData()).
     * \param[in] projectile Projectile species name (see
     * findNucleusData()).
     * \param[in] inb Impact parameter \f$b\f$ [fm], stored for
     * oLSIntegrand()'s use.
     * \param[in] setWSDeformParams Forwarded to findNucleusData() for
     * both species.
     * \param[in] R_WS Forwarded to findNucleusData() for both species.
     * \param[in] a_WS Forwarded to findNucleusData() for both species.
     * \param[in] beta2 Forwarded to findNucleusData() for both species.
     * \param[in] beta3 Forwarded to findNucleusData() for both species.
     * \param[in] beta4 Forwarded to findNucleusData() for both species.
     * \param[in] gamma Forwarded to findNucleusData() for both species.
     * \param[in] forceDminFlag Forwarded to findNucleusData() for both
     * species.
     * \param[in] d_min Forwarded to findNucleusData() for both species.
     * \param[in] dR_np Forwarded to findNucleusData() for both species.
     * \param[in] da_np Forwarded to findNucleusData() for both species.
     * \param[in] imax Interpolation-table resolution, stored as \c
     * glauberData_.interMax (\c sCutoff is fixed at 12 fm).
     */
    void initGlauber(
        double sigmaNN, std::string target, std::string projectile, double inb,
        bool setWSDeformParams, double R_WS, double a_WS, double beta2,
        double beta3, double beta4, double gamma, bool forceDminFlag,
        double d_min, double dR_np, double da_np, int imax);
    /**
     * Envelope-area helper for sampleTARejection()'s rejection sampling:
     * the (unnormalized) area under its Gaussian-like envelope function
     * out to radius \p x.
     * \param[in] x Upper radius [fm].
     * \param[in] A Overall envelope-height scale (from \p A's inelastic
     * cross section, see sampleTARejection()).
     * \return The envelope's area out to \p x.
     */
    double areaTA(double x, double A);
    /**
     * Samples one nucleon's transverse position from the projectile's
     * or target's thickness function \f$T_A(r)\f$/\f$T_B(r)\f$ via
     * rejection sampling against a Gaussian-like envelope, logging a
     * warning if the true distribution ever exceeds that envelope.
     * \param[in] random Random-number source.
     * \param[in] nucleus Which nucleus' thickness function to sample
     * from.
     * \return A ReturnValue with `x`/`y` set to the sampled position
     * and \c collided set to `0`; \c z, \c phi and \c proton are left
     * unset (see ReturnValue::phi for why that's safe today).
     */
    ReturnValue sampleTARejection(Random *random, NucleusRole nucleus);
};
#endif  // SRC_GLAUBER_H_

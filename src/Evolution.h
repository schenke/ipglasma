// Evolution.h is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_EVOLUTION_H_
#define SRC_EVOLUTION_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include "FFT.h"
#include "Group.h"
#include "Lattice.h"
#include "Parameters.h"
#include "PrettyOstream.h"

/**
 * Advances the classical Yang-Mills fields built by Init via a leapfrog
 * time evolution: the coordinate fields (\f$\phi\f$, the transverse
 * gauge links `Ux`/`Uy`) live at integer-\f$\tau\f$ steps and are
 * updated by evolvePhi()/evolveU(), while the momentum fields (the
 * electric fields `U`/`U2`, and \f$\pi\f$) live at half-integer steps
 * and are updated by evolvePi()/evolveE(). Also computes the
 * energy-momentum tensor (tmunu()), derived observables (eccentricity(),
 * u()), and an optional final-time gluon spectrum/multiplicity estimate
 * (multiplicity()) via Coulomb-gauge fixing and FFT.
 */
class Evolution {
  private:
    /// Owned FFT instance, sized to the lattice at construction; used by
    /// multiplicity() to Fourier-transform the gauge-fixed fields.
    FFT *fft_;
    /// \f$n(k_T)\f$ spectrum read by readNkt() from a previous run's
    /// `multiplicity<id>.dat`; zero-initialized so a short input file
    /// (fewer than 100 rows) leaves untouched entries at a defined `0`
    /// rather than uninitialized garbage.
    double nIn_[100] = {};
    /// Log sink for progress/warning/error messages.
    PrettyOstream messager_;

  public:
    /**
     * Constructs an Evolution for a lattice of the given dimensions.
     * \param[in] nn Two-element array `{N, N}`, the transverse lattice
     * dimensions, forwarded to the owned FFT instance.
     */
    Evolution(const int nn[]) { fft_ = new FFT(nn); }

    /**
     * Destroys this Evolution, freeing the owned FFT instance.
     */
    ~Evolution() { delete fft_; }

    /// Owns \c fft_, a raw pointer freed in the destructor; default
    /// copies would double-free it, so disable copying (nothing needs
    /// it).
    Evolution(const Evolution &) = delete;
    /// \copydoc Evolution(const Evolution &)
    Evolution &operator=(const Evolution &) = delete;

    /**
     * Top-level evolution driver: after an initial coordinate half-step
     * (evolvePhi()/evolveU() from \f$\tau=0\f$ to \f$\tau=d\tau\f$,
     * since the momenta at \f$\tau=d\tau/2\f$ are taken equal to their
     * \f$\tau=0\f$ values), leapfrog-advances the fields
     * (evolveStepPersistent(), an anonymous-namespace helper in
     * Evolution.cpp that runs one shared `#pragma omp parallel` team
     * through evolvePi()/evolveE()/evolvePhi()/evolveU() in order) up to
     * `param->getMaxtime()` (or, if `getInverseQsForMaxTime()`, up to
     * \f$1/Q_s\f$). At the final time step (and, if `getWriteOutputs()
     * == 5`, at four additional fixed intermediate times), temporarily
     * recenters the momenta from \f$\tau_{n-1/2}\f$ to \f$\tau_n\f$ to
     * measure tmunu() and either u() or `MyEigen::writeTmunu4D()`
     * (depending on `getWriteEpsilonUHydro()`), then restores the
     * unmodified momenta so the measurement cannot perturb the
     * trajectory. At the very end, runs checkGaussLaw(), and -- if
     * `getComputeGluonMultiplicity()` -- eccentricity() and
     * multiplicity(), stopping early if multiplicity() reports no
     * collision.
     * \param[in,out] lat Lattice to evolve in place.
     * \param[in] group Group instance, forwarded to multiplicity().
     * \param[in,out] param Simulation parameters.
     */
    void run(Lattice *lat, Group *group, Parameters *param);
    /**
     * `writeOutputs==3` diagnostic: writes `epsilonInitialPlot<id>.dat`,
     * a per-cell map of the (running-coupling-weighted) energy density
     * at the initial time.
     * \note Computes its running-coupling factor via only the "3
     * flavors, local \f$Q_s\f$" formula, unconditionally -- unlike
     * writeEpsilonIntermediatePlot() (which shares the
     * anonymous-namespace `computeRunningCouplingGfactor` helper and so
     * respects `getRunWithLocalQs()`). This looks like a pre-existing
     * inconsistency between the two writers, preserved exactly as found
     * rather than silently changed.
     * \param[in] lat Lattice to read the energy density from.
     * \param[in] param Simulation parameters.
     */
    void writeEpsilonInitialPlot(Lattice *lat, Parameters *param);
    /**
     * `writeOutputs==3` diagnostic: writes `epsilonIntermediatePlot<id>
     * .dat`, the same per-cell energy-density map as
     * writeEpsilonInitialPlot() but at \f$\tau=\f$`itmax/2`, using the
     * shared `computeRunningCouplingGfactor` helper for its
     * running-coupling factor.
     * \param[in] lat Lattice to read the energy density from.
     * \param[in] param Simulation parameters.
     */
    void writeEpsilonIntermediatePlot(Lattice *lat, Parameters *param);
    /**
     * Leapfrog coordinate update for the transverse gauge links: rotates
     * `Ux`/`Uy` by \f$\exp(i g^2 d\tau/(\tau+d\tau/2)\,U)\f$ (a
     * second-order Padé approximant of the exponential), using the
     * electric fields `U`/`U2` as the generator. Runs
     * `evolveUTeam` (an anonymous-namespace helper in Evolution.cpp)
     * inside its own `#pragma omp parallel` region.
     * \param[in,out] lat Lattice whose `Ux`/`Uy` are updated in place.
     * \param[in] param Simulation parameters.
     * \param[in] dtau Time step [lattice units].
     * \param[in] tau Current proper time [lattice units], evaluated at
     * the half-step \f$\tau+d\tau/2\f$ inside the update.
     */
    void evolveU(Lattice *lat, Parameters *param, double dtau, double tau);
    /**
     * Leapfrog coordinate update for the scalar field: \f$\phi \mathrel{+}=
     * (\tau+d\tau/2)\,d\tau\,\pi\f$. Runs `evolvePhiTeam` (an
     * anonymous-namespace helper in Evolution.cpp) inside its own
     * `#pragma omp parallel` region.
     * \param[in,out] lat Lattice whose \c Uy2 (\f$\phi\f$) is updated in
     * place.
     * \param[in] param Simulation parameters.
     * \param[in] dtau Time step [lattice units].
     * \param[in] tau Current proper time [lattice units].
     */
    void evolvePhi(Lattice *lat, Parameters *param, double dtau, double tau);
    /**
     * Leapfrog momentum update for \f$\pi\f$: adds
     * \f$(d\tau/\tau)\f$ times the covariant discrete Laplacian of
     * \f$\phi\f$ (parallel-transported via `Ux`/`Uy` to its four
     * neighbors). Runs `evolvePiTeam` (an anonymous-namespace helper in
     * Evolution.cpp) inside its own `#pragma omp parallel` region.
     * \param[in,out] lat Lattice whose \c Ux2 (\f$\pi\f$) is updated in
     * place.
     * \param[in] param Simulation parameters.
     * \param[in] dtau Time step [lattice units].
     * \param[in] tau Current proper time [lattice units].
     */
    void evolvePi(Lattice *lat, Parameters *param, double dtau, double tau);
    /**
     * Leapfrog momentum update for the electric fields `U`/`U2`: adds
     * the traceless anti-Hermitian plaquette force (from the four
     * spatial plaquettes touching each link) and the \f$\phi\f$-\f$\pi\f$
     * commutator force. Runs `evolveETeam` (an anonymous-namespace
     * helper in Evolution.cpp, using the shared `addEForceSU3` kernel)
     * inside its own `#pragma omp parallel` region.
     * \param[in,out] lat Lattice whose `U`/`U2` are updated in place.
     * \param[in] param Simulation parameters.
     * \param[in] dtau Time step [lattice units].
     * \param[in] tau Current proper time [lattice units].
     */
    void evolveE(Lattice *lat, Parameters *param, double dtau, double tau);
    /**
     * Diagnostic: computes, at every cell, the discrete Gauss-law
     * constraint violation
     * \f$U_x(x-\hat x)^\dagger E_1(x-\hat x) U_x(x-\hat x) - E_1(x) +
     * U_y(x-\hat y)^\dagger E_2(x-\hat y) U_y(x-\hat y) - E_2(x) -
     * i[\phi,\pi]\f$ (evaluated using the momenta at their current
     * \f$\tau-d\tau/2\f$), and logs the largest value found.
     * \param[in] lat Lattice to read `U`/`U2`/`Ux`/`Uy`/`Uy2`/\c
     * Ux2 from.
     * \param[in] param Simulation parameters.
     */
    void checkGaussLaw(Lattice *lat, Parameters *param);
    /**
     * Computes the spatial eccentricities \f$\varepsilon_1,\ldots,
     * \varepsilon_6\f$ and their event-plane angles
     * \f$\Psi_1,\ldots,\Psi_6\f$ from the energy-density-weighted
     * (running-coupling-corrected, `getEpsilon()*getutau()`-weighted,
     * cells below \p cutoff excluded) spatial moments, first recentering
     * on the energy-weighted centroid. If \p doAniso is `0`, appends one
     * row to `eccentricities<id>.dat`. If \p doAniso is `1`, instead
     * appends to `anisotropy<id>.dat` the \f$T^{xx}-T^{yy}\f$ spatial
     * anisotropy (via the anonymous-namespace `computeRotatedAnisotropy`
     * helper) resampled at ten angles around the flow-velocity event
     * plane \f$\Psi_U\f$ (computed from `getux()`/`getuy()`).
     * \param[in] lat Lattice to read the energy density (and, for
     * \p doAniso==1, `Txx`/`Txy`/`Tyy`/flow velocity) from.
     * \param[in] param Simulation parameters.
     * \param[in] it Current time step index, used for the output file's
     * time column and (on `it==1`) to seed `param->setPsi()`.
     * \param[in] cutoff Energy-density cutoff below which a cell is
     * excluded from the weighted averages [\f$\Lambda_{QCD}^4\f$-like
     * units, i.e. roughly 1/fm\f$^4\f$].
     * \param[in] doAniso Selects which output file/quantity is written
     * (`0`: eccentricities, `1`: rotated-tensor anisotropy).
     */
    void eccentricity(
        Lattice *lat, Parameters *param, int it, double cutoff, int doAniso);
    /**
     * Computes every component of the energy-momentum tensor
     * \f$T^{\mu\nu}\f$ at every cell (the spatial plaquette, the
     * diagonal electric and magnetic/gradient contributions, and the
     * six off-diagonal components), storing the results into
     * `lat->cells`. Runs its constituent anonymous-namespace team
     * helpers (`tmunuPlaquetteTeam`, `tmunuDiagonalElectricTeam`,
     * `tmunuDiagonalMagneticTeam`, `tmunuNormalizeDiagonalTeam`,
     * `tmunuOffDiagonalTeam`) in sequence inside one shared
     * `#pragma omp parallel` region.
     * \param[in,out] lat Lattice to read the fields from and write
     * \f$T^{\mu\nu}\f$ into (via `lat->cells`).
     * \param[in] param Simulation parameters.
     * \param[in] it Current time step index (\f$\tau=it\cdot d\tau\f$),
     * used to convert the momentum fields' lattice normalization to
     * physical units.
     */
    void tmunu(Lattice *lat, Parameters *param, int it);
    /**
     * Writes a binary snapshot (`evolvedFields<id>_it<it>.ipgf`) of the
     * six evolved fields (\f$\phi\f$, \f$\pi\f$, \c U (\f$E_1\f$), \c U2
     * (\f$E_2\f$), \c Ux, \c Uy) as full \f$3\times3\f$ complex
     * matrices per cell, little-endian single-precision, preceded by an
     * 8-byte magic string, an 8-byte metadata length, and a JSON
     * metadata header describing the layout/units/staggering.
     * \param[in] lat Lattice to read the fields from.
     * \param[in] param Simulation parameters.
     * \param[in] it Current time step index, used for the file name and
     * the metadata's time/tau fields.
     */
    void writeEvolvedFields(Lattice *lat, Parameters *param, int it);
    /**
     * Computes the hydrodynamic flow-velocity field \f$u^\mu\f$ by
     * delegating to `MyEigen::flowVelocity4D`.
     * \param[in,out] lat Lattice to read \f$T^{\mu\nu}\f$ from and write
     * the flow velocity into.
     * \param[in] param Simulation parameters.
     * \param[in] it Current time step index.
     * \param[in] finalFlag Forwarded to `MyEigen::flowVelocity4D`
     * (selects final- vs. intermediate-time output naming/behavior).
     */
    void u(Lattice *lat, Parameters *param, int it, bool finalFlag);
    /**
     * Fixes transverse Coulomb gauge (`GaugeFix::fftChi`), then computes
     * the azimuthally averaged gluon transverse-momentum spectrum by
     * Fourier-transforming the gauge-fixed electric fields \c U, \c U2,
     * and \f$\pi\f$ (via the anonymous-namespace `prepareSpectrumField`/
     * `accumulateGluonSpectrum` helpers) and binning \f$dN/d^2k_T\f$ and
     * \f$dE/d^2k_T\f$. Integrates these into \f$dN/dy\f$ (or
     * \f$dN/d\eta\f$) and \f$dE/dy\f$, with separate cut sums above
     * \f$k_T>3\f$ and \f$6\f$ GeV. At the final time step, optionally
     * hadronizes the spectrum (hadronizeAndWriteMultiplicity(), if
     * `getWriteOutputs()==3`) and writes `NpartdNdy-t<t>-<id>.dat` and
     * the `gluonMultiplicity<id>.json` target
     * (writeGluonMultiplicityTarget()).
     * \param[in,out] lat Lattice to read the fields from (gauge-fixed in
     * place by `GaugeFix::fftChi`).
     * \param[in] group Group instance, forwarded to `GaugeFix::fftChi`.
     * \param[in] param Simulation parameters.
     * \param[in] it Current time step index.
     * \return `1` on success (`param->setSuccess(1)` is also called);
     * `0` if no collision was found (\f$dN/dy=0\f$), signaling the
     * caller (run()) to stop this event so it can be resampled with a
     * new random seed.
     */
    int multiplicity(Lattice *lat, Group *group, Parameters *param, int it);
    /**
     * multiplicity()'s `writeOutputs==3` hadronization step: convolves
     * the binned gluon spectrum \p n with the KKP fragmentation function
     * (`Fragmentation::kkp`), integrated over the fragmentation variable
     * \f$z\f$ via a GSL cubic spline, to get a hadron \f$p_T\f$ spectrum,
     * writing `multiplicityHadrons<id>.dat`. \p Nhgsl (length
     * `hbins+1`) is scratch/output space owned by the caller.
     * \param[in] param Simulation parameters.
     * \param[in] a Lattice spacing [fm].
     * \param[in] dkt Momentum-bin width [lattice units].
     * \param[in] bins Number of gluon-spectrum bins in \p n.
     * \param[in] n Binned, azimuthally averaged gluon spectrum
     * \f$dN/d^2k_T\f$ [lattice units], length \p bins.
     * \param[out] Nhgsl Filled with the hadron \f$dN/dp_T\f$ spectrum,
     * length `hbins+1`.
     * \param[in] hbins Number of hadron-spectrum bins.
     */
    void hadronizeAndWriteMultiplicity(
        Parameters *param, double a, double dkt, int bins, const double *n,
        double *Nhgsl, int hbins);

    /**
     * Writes `gluonMultiplicity<id>.json`, a structured summary of this
     * event's gluon spectrum and integrated multiplicity/energy
     * (including the cut sums above 3 and 6 GeV and the binned spectrum
     * itself), intended as a downstream analysis/ML training target.
     * \param[in] param Simulation parameters.
     * \param[in] it Current time step index.
     * \param[in] a Lattice spacing [fm].
     * \param[in] dtau Time step [lattice units], used to convert \p it
     * to a physical time.
     * \param[in] dNPrimary \f$dN/dy\f$ (or \f$d\eta\f$) from the direct
     * per-bin sum.
     * \param[in] dNBinned \f$dN/dy\f$ from the alternate binned-weight
     * integration, included as a cross-check.
     * \param[in] dEPrimary \f$dE/dy\f$ from the direct per-bin sum
     * [GeV].
     * \param[in] dEBinned \f$dE/dy\f$ from the alternate binned-weight
     * integration [GeV], included as a cross-check.
     * \param[in] dNCut3 \f$dN/dy\f$ restricted to \f$k_T>3\f$ GeV.
     * \param[in] dECut3 \f$dE/dy\f$ restricted to \f$k_T>3\f$ GeV [GeV].
     * \param[in] dNCut6 \f$dN/dy\f$ restricted to \f$k_T>6\f$ GeV.
     * \param[in] dECut6 \f$dE/dy\f$ restricted to \f$k_T>6\f$ GeV [GeV].
     * \param[in] spectrumN Binned \f$dN/d^2k_T\f$ [lattice units],
     * length \p bins.
     * \param[in] spectrumE Binned \f$dE/d^2k_T\f$ [lattice units],
     * length \p bins.
     * \param[in] spectrumCounts Number of lattice momentum modes
     * averaged into each bin, length \p bins.
     * \param[in] bins Number of spectrum bins.
     * \param[in] dkt Momentum-bin width [lattice units].
     */
    void writeGluonMultiplicityTarget(
        Parameters *param, int it, double a, double dtau, double dNPrimary,
        double dNBinned, double dEPrimary, double dEBinned, double dNCut3,
        double dECut3, double dNCut6, double dECut6, const double *spectrumN,
        const double *spectrumE, const int *spectrumCounts, int bins,
        double dkt);
    /**
     * Standalone post-processing utility: reads a previous run's
     * `multiplicity<id>.dat` (into \c nIn_) and `NpartdNdy<id>.dat`,
     * recomputes \f$dN/d\eta\f$ from the pseudorapidity Jacobian
     * (`param->getJacobianm()`/`getRoots()`), and writes
     * `NpartdNdy-mod.dat`. Not part of the normal run() flow; terminates
     * the process (`exit(1)`) unconditionally when done, and also exits
     * early if either input file is missing.
     * \param[in] param Simulation parameters.
     */
    void readNkt(Parameters *param);
};

#endif  // SRC_EVOLUTION_H_

// GluonMultiplicity.h is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_GLUONMULTIPLICITY_H_
#define SRC_GLUONMULTIPLICITY_H_

#include "FFT.h"
#include "Group.h"
#include "Lattice.h"
#include "Parameters.h"
#include "PrettyOstream.h"

/**
 * The gluon transverse-momentum spectrum and multiplicity of the evolved
 * classical fields (in transverse Coulomb gauge), with the optional
 * hadronization and the per-event output files.
 */
class GluonMultiplicity {
  public:
    /**
     * Constructs the measurement for a lattice of the given dimensions.
     * \param[in] nn Two-element array `{N, N}`, the transverse lattice
     * dimensions, forwarded to the owned FFT instance.
     */
    explicit GluonMultiplicity(const int nn[]) : fft_(nn) {}

    /**
     * Fixes transverse Coulomb gauge (`GaugeFix::fftChi`), then computes
     * the azimuthally averaged gluon transverse-momentum spectrum by
     * Fourier-transforming the gauge-fixed electric fields \c U, \c U2,
     * and \f$\pi\f$ (via the anonymous-namespace `prepareSpectrumField`/
     * `accumulateGluonSpectrum` helpers) and binning \f$dN/d^2k_T\f$ and
     * \f$dE/d^2k_T\f$. Integrates these into \f$dN/dy\f$ (or
     * \f$dN/d\eta\f$) and \f$dE/dy\f$, with separate cut sums above
     * \f$k_T>3\f$ and \f$6\f$ GeV. At the final time step, optionally
     * hadronizes the spectrum (hadronizeAndWrite(), if
     * `param->output.writeHadronSpectrum`) and writes `NpartdNdy-t<t>-<id>.dat`
     * and
     * the `gluonMultiplicity<id>.json` target
     * (writeTarget()).
     * The file is described in \ref md_OUTPUT "OUTPUT.md".
     * \param[in,out] lat Lattice to read the fields from (gauge-fixed in
     * place by `GaugeFix::fftChi`).
     * \param[in] group Group instance, forwarded to `GaugeFix::fftChi`.
     * \param[in] param Simulation parameters.
     * \param[in] it Current time step index.
     * \return `1` on success (`param->event.success = 1` is also called);
     * `0` if no collision was found (\f$dN/dy=0\f$), signaling the
     * caller (Evolution::run()) to stop this event so it can be resampled with
     * a new random seed.
     */
    int compute(Lattice *lat, Group *group, Parameters *param, int it);
    /**
     * Standalone post-processing utility: reads a previous run's
     * `multiplicity<id>.dat` and `NpartdNdy<id>.dat`,
     * recomputes \f$dN/d\eta\f$ from the pseudorapidity Jacobian
     * (`param->colorCharge.jacobianMass`/`param->collision.sqrtS`), and writes
     * `NpartdNdy-mod.dat`. Not called by the program (no version of it
     * writes these input files under these names); kept for post-processing
     * such files. Terminates the process (`exit(1)`) unconditionally when
     * done, and also exits early if either input file is missing.
     * \param[in] param Simulation parameters.
     */
    static void readNkt(Parameters *param);

  private:
    /**
     * compute()'s `writeHadronSpectrum` hadronization step: convolves
     * the binned gluon spectrum \p n with the KKP fragmentation function
     * (`Fragmentation::kkp`), integrated over the fragmentation variable
     * \f$z\f$ via a GSL cubic spline, to get a hadron \f$p_T\f$ spectrum,
     * writing `multiplicityHadrons<id>.dat`. \p Nhgsl (length
     * `hbins+1`) is scratch/output space owned by the caller.
     * The file is described in \ref md_OUTPUT "OUTPUT.md".
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
    void hadronizeAndWrite(
        Parameters *param, double a, double dkt, int bins, const double *n,
        double *Nhgsl, int hbins);

    /**
     * Writes `gluonMultiplicity<id>.json`, a structured summary of this
     * event's gluon spectrum and integrated multiplicity/energy
     * (including the cut sums above 3 and 6 GeV and the binned spectrum
     * itself), intended as a downstream analysis/ML training target.
     * The file is described in \ref md_OUTPUT "OUTPUT.md".
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
    void writeTarget(
        Parameters *param, int it, double a, double dtau, double dNPrimary,
        double dNBinned, double dEPrimary, double dEBinned, double dNCut3,
        double dECut3, double dNCut6, double dECut6, const double *spectrumN,
        const double *spectrumE, const int *spectrumCounts, int bins,
        double dkt);
    /// FFT used to Fourier-transform the gauge-fixed fields (and by the
    /// Coulomb gauge fixing).
    FFT fft_;
    /// Log sink for progress/warning/error messages.
    PrettyOstream messager_;
};

#endif  // SRC_GLUONMULTIPLICITY_H_

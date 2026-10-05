// GluonMultiplicity.cpp is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "GluonMultiplicity.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Fragmentation.h"
#include "GaugeFix.h"
#include "Instrumentation.h"
#include "LatticeIndex.h"
#include "PhysConst.h"
#include "RunningCoupling.h"

using Fragmentation::kkp;
using PhysConst::hbarc;
using PhysConst::m_kaon;
using PhysConst::m_pion;
using PhysConst::m_proton;

using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;

namespace {

/**
 * Records the elapsed wall time since \p started under \p phase in
 * the global profiler and resets \p started to the current time, so the next
 * phase's elapsed time is measured from here.
 * \param[in] phase Profiler phase name to add the elapsed time to.
 * \param[in,out] started Wall-clock time the phase began; updated to
 * now on return.
 */
void addPhaseAndRestart(const char *phase, double &started) {
    const double now = ipg::wallSeconds();
    ipg::Profiler::instance().add(phase, now - started);
    started = now;
}

/**
 * Fills \p E1[pos] from \p sourceField[pos] (one of `lat->U`/`U2`/
 * `Ux2`), scaling by the square root of the local running-coupling
 * g-factor (computeRunningCouplingGfactor()) unless \f$\alpha_s\f$
 * runs with \f$k_T\f$ (in which case the \f$k_T\f$-dependent factor
 * is applied later, per-mode, in accumulateGluonSpectrum() instead).
 * \param[in] lat Lattice to read \f$g^2\mu_A^2\f$/\f$g^2\mu_B^2\f$
 * from, forwarded to computeRunningCouplingGfactor().
 * \param[in] param Simulation parameters.
 * \param[in] N Lattice side length.
 * \param[in] a Lattice spacing [fm].
 * \param[in] g Coupling \f$g\f$.
 * \param[in] c Running-coupling shape parameter.
 * \param[in] muZero \f$\mu_0\f$ in the running-coupling formula.
 * \param[in] sourceField Field to copy from (`lat->U`, `lat->U2`, or
 * `lat->Ux2`).
 * \param[out] E1 Filled with the (optionally rescaled) field, as a
 * pointer-per-cell view ready for `FFT::fftn`.
 */
void prepareSpectrumField(
    Lattice *lat, Parameters *param, int N, double a, double g, double c,
    double muZero, const std::vector<Matrix> &sourceField,
    std::vector<Matrix *> &E1) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int pos = lat->positionFromXY(i, j);
            double gfactor = computeRunningCouplingGfactor(
                lat, param, pos, N, a, g, c, muZero);
            if (!param->coupling.runWithKt) {
                *E1[pos] = sourceField[pos] * sqrt(gfactor);
            } else {
                *E1[pos] = sourceField[pos];
            }
        }
    }
}

/**
 * Accumulates \p E1's (already FFT'd) momentum-space spectrum into
 * `dNdeta`/`dEdeta` and the `n`/`E`/`n2` \f$k_T\f$ bins. Called
 * once each for the \c U, \c U2, and \c Ux2 (\f$\pi\f$) fields by
 * GluonMultiplicity::compute().
 * \param[in] param Simulation parameters.
 * \param[in] N Lattice side length.
 * \param[in] it Current time step index.
 * \param[in] dtau Time step [lattice units].
 * \param[in] g Coupling \f$g\f$.
 * \param[in] a Lattice spacing [fm].
 * \param[in] c Running-coupling shape parameter.
 * \param[in] muZero \f$\mu_0\f$ in the running-coupling formula.
 * \param[in] dkt Momentum-bin width [lattice units].
 * \param[in] bins Number of \f$k_T\f$ bins in `n`/`E`/`n2`/
 * \p counter.
 * \param[in] E1 The field's FFT'd momentum-space values, as a
 * pointer-per-cell view.
 * \param[in] useElectricNormalization Selects \c nkt's electric-field
 * normalization (`true`, \f$g^2/((it-0.5)d\tau)\f$) vs. the
 * \f$\pi\f$-field normalization (`false`, \f$(it-0.5)d\tau\f$, no
 * \f$g^2\f$).
 * \param[in] accumulateCounter Whether to also increment \p counter;
 * only one of the three per-field passes needs to, since all three
 * share the same \f$k_T\f$ grid.
 * \param[in,out] dNdeta Running sum of \f$dN/dy\f$ (or \f$d\eta\f$),
 * incremented by this field's contribution.
 * \param[in,out] dEdeta Running sum of \f$dE/dy\f$ (or \f$d\eta\f$),
 * incremented by this field's contribution.
 * \param[in,out] n Binned \f$dN/d^2k_T\f$, incremented by this
 * field's contribution, length \p bins.
 * \param[in,out] E Binned \f$dE/d^2k_T\f$, incremented by this
 * field's contribution, length \p bins.
 * \param[in,out] n2 Alternate binned \f$dN/d^2k_T\f$ normalization
 * (used for a cross-check), incremented by this field's contribution,
 * length \p bins.
 * \param[in,out] counter Number of lattice momentum modes falling
 * into each bin, incremented only if \p accumulateCounter, length
 * \p bins.
 */
void accumulateGluonSpectrum(
    Parameters *param, int N, int it, double dtau, double g, double a, double c,
    double muZero, double dkt, int bins, const std::vector<Matrix *> &E1,
    bool useElectricNormalization, bool accumulateCounter, double &dNdeta,
    double &dEdeta, double *n, double *E, double *n2, int *counter) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double nkt = 0.;
            int pos = latticeIndex(i, j, N);
            int npos = latticeIndex(N - i, N - j, N);

            double kx =
                2. * M_PI
                * (-0.5 + static_cast<double>(i) / static_cast<double>(N));
            double ky =
                2. * M_PI
                * (-0.5 + static_cast<double>(j) / static_cast<double>(N));
            double kt2 =
                4.
                * (sin(kx / 2.) * sin(kx / 2.) + sin(ky / 2.) * sin(ky / 2.));
            double omega2 =
                4.
                * (sin(kx / 2.) * sin(kx / 2.)
                   + sin(ky / 2.) * sin(ky / 2.));  // lattice dispersion
                                                    // relation (this is
                                                    // omega squared)

            // i=0 or j=0 have no negative k_T value available
            if (i != 0 && j != 0) {
                if (omega2 != 0) {
                    if (useElectricNormalization) {
                        nkt =
                            2. / sqrt(omega2) / static_cast<double>(N * N)
                            * (g * g / ((it - 0.5) * dtau)
                               * ((((*E1[pos]) * (*E1[npos])).trace()).real()));
                    } else {
                        nkt =
                            2. / sqrt(omega2) / static_cast<double>(N * N)
                            * (((it - 0.5) * dtau)
                               * ((((*E1[pos]) * (*E1[npos])).trace()).real()));
                    }
                    if (param->coupling.runWithKt) {
                        nkt *= computeRunningCouplingGfactorFromScale(
                            g, muZero, c, param->coupling.LambdaQCD,
                            param->coupling.nFlavors,
                            param->coupling.runningCouplingQsFactor * sqrt(kt2)
                                * hbarc / a);
                    }
                }

                dNdeta += nkt;
                dEdeta += nkt * sqrt(omega2) * hbarc / a;

                for (int ik = 0; ik < bins; ik++) {
                    if (abs(sqrt(kt2)) > ik * dkt
                        && abs(sqrt(kt2)) <= (ik + 1) * dkt) {
                        n[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2) * 2 * M_PI
                                 * sqrt(kt2) * dkt * N * N / M_PI / M_PI / 2.
                                 / 2.;
                        E[ik] += sqrt(omega2) * hbarc / a * nkt / dkt / 2 / M_PI
                                 / sqrt(kt2) * 2 * M_PI * sqrt(kt2) * dkt * N
                                 * N / M_PI / M_PI / 2. / 2.;
                        n2[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2);
                        // dividing by bin size; bin is dkt times Jacobian
                        // k(=ik*dkt) times 2Pi in phi times the correct
                        // number of counts for an infinite lattice: area in
                        // bin divided by total area
                        if (accumulateCounter) {
                            counter[ik] += 1;  // number of entries in n[ik]
                        }
                    }
                }
            }
        }
    }
}

/**
 * GluonMultiplicity::compute()'s per-bin \f$dN/dy\f$, \f$dE/dy\f$
 * weight: the phase-space factor \f$(ik+0.5)\,dk_T^2\,2\pi\f$, times
 * a Jacobian ratio when the rapidity input is actually a
 * pseudorapidity (the same factor previously computed identically
 * three times -- unconditionally, and again inside the \f$k_T>3\f$
 * and \f$k_T>6\f$ GeV cuts -- for both `usePseudoRapidity` branches).
 * \param[in] param Simulation parameters.
 * \param[in] m Jacobian mass term [GeV] (`param->colorCharge.jacobianMass`).
 * \param[in] ik Bin index.
 * \param[in] dkt Momentum-bin width [lattice units].
 * \param[in] a Lattice spacing [fm].
 * \return The per-bin weight.
 */
double computeMultiplicityBinWeight(
    Parameters *param, double m, int ik, double dkt, double a) {
    const double base = (ik + 0.5) * dkt * dkt * 2. * M_PI;
    if (!param->colorCharge.usePseudoRapidity) {
        return base;
    }
    return base * cosh(param->colorCharge.rapidity())
           / (sqrt(
               pow(cosh(param->colorCharge.rapidity()), 2.)
               + m * m
                     / (((ik + 0.5) * dkt / a * hbarc)
                        * ((ik + 0.5) * dkt / a * hbarc))));
}

}  // namespace

int GluonMultiplicity::compute(
    Lattice *lat, Group *group, Parameters *param, int it) {
    IPG_PROFILE_SCOPE("observables.gluon_multiplicity");
    int N = param->lattice.size;
    int npos, pos;
    double L = param->lattice.L;
    double a = L / N;  // lattice spacing in fm
    double kx, ky, kt2, omega2;
    double g = param->coupling.g;
    int nn[2];
    nn[0] = N;
    nn[1] = N;
    double dtau = param->run.dtau;
    double nkt;
    const int bins = 100;
    double n[bins];   // k_T array
    double E[bins];   // k_T array
    double n2[bins];  // k_T array
    int counter[bins];
    double dkt = 2.83 / static_cast<double>(bins);
    double dNdeta = 0.;
    double dNdeta2 = 0.;
    double dNdetaCut = 0.;
    double dNdetaCut2 = 0.;
    double dEdetaCut = 0.;
    double dEdetaCut2 = 0.;
    double dEdeta = 0.;
    double dEdeta2 = 0.;

    stringstream strNpartdNdy_name;
    strNpartdNdy_name << "NpartdNdy-t" << it * dtau * a << "-"
                      << param->event.eventId << ".dat";
    string NpartdNdy_name;
    NpartdNdy_name = strNpartdNdy_name.str();
    messager_ << "[GluonMultiplicity::compute]: Measuring multiplicity ... ";
    messager_.flush("info");

    // fix transverse Coulomb gauge
    GaugeFix gaugefix;

    double maxtime;
    if (param->evolution.inverseQsForMaxTime) {
        maxtime = 1. / param->event.averageQs * hbarc;
        messager_ << "[GluonMultiplicity::compute]: maximal evolution time = "
                  << maxtime << " fm";
        messager_.flush("info");
    } else {
        maxtime = param->evolution.maxTime;  // maxtime is in fm
    }

    int itmax = static_cast<int>(floor(maxtime / (a * dtau) + 1e-10));

    double multiplicityPhaseStart = ipg::wallSeconds();
    gaugefix.fftChi(&fft_, lat, group, param, 4000);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.gauge_fix", multiplicityPhaseStart);
    // gauge is fixed

    // E1Storage owns the N*N scratch matrices; E1 is a pointer-per-cell
    // view over it for FFT::fftn's T** interface.
    std::vector<Matrix> E1Storage(N * N, Matrix(0.));
    std::vector<Matrix *> E1(N * N);
    for (int i = 0; i < N * N; i++) {
        E1[i] = &E1Storage[i];
    }
    addPhaseAndRestart(
        "observables.gluon_multiplicity.allocate", multiplicityPhaseStart);

    double c = param->coupling.c;
    double muZero = param->coupling.mu0;

    prepareSpectrumField(lat, param, N, a, g, c, muZero, lat->U, E1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.prepare_E1", multiplicityPhaseStart);

    // do Fourier transforms
    fft_.fftn(E1.data(), E1.data(), nn, 1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.fft_E1", multiplicityPhaseStart);

    for (int ik = 0; ik < bins; ik++) {
        n[ik] = 0.;
        E[ik] = 0.;
        n2[ik] = 0.;
        counter[ik] = 0;
    }

    const int hbins = 2000;
    double Nhgsl[hbins + 1];

    addPhaseAndRestart(
        "observables.gluon_multiplicity.setup_bins", multiplicityPhaseStart);

    accumulateGluonSpectrum(
        param, N, it, dtau, g, a, c, muZero, dkt, bins, E1, true, true, dNdeta,
        dEdeta, n, E, n2, counter);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.spectrum_E1", multiplicityPhaseStart);

    /// -------- 2 ---------

    prepareSpectrumField(lat, param, N, a, g, c, muZero, lat->U2, E1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.prepare_E2", multiplicityPhaseStart);

    fft_.fftn(E1.data(), E1.data(), nn, 1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.fft_E2", multiplicityPhaseStart);

    accumulateGluonSpectrum(
        param, N, it, dtau, g, a, c, muZero, dkt, bins, E1, true, false, dNdeta,
        dEdeta, n, E, n2, counter);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.spectrum_E2", multiplicityPhaseStart);

    /// ------3 --------

    prepareSpectrumField(lat, param, N, a, g, c, muZero, lat->Ux2, E1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.prepare_pi", multiplicityPhaseStart);

    // do Fourier transforms
    fft_.fftn(E1.data(), E1.data(), nn, 1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.fft_pi", multiplicityPhaseStart);

    accumulateGluonSpectrum(
        param, N, it, dtau, g, a, c, muZero, dkt, bins, E1, false, false,
        dNdeta, dEdeta, n, E, n2, counter);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.spectrum_pi", multiplicityPhaseStart);

    double m, P;
    m = param->colorCharge.jacobianMass;                           // in GeV
    P = 0.13 + 0.32 * pow(param->collision.sqrtS / 1000., 0.115);  // in GeV

    for (int ik = 0; ik < bins; ik++) {
        if (counter[ik] > 0) {
            n[ik] = n[ik] / static_cast<double>(counter[ik]);
            E[ik] = E[ik] / static_cast<double>(counter[ik]);
            // integrate, gives a ik*dkt*2pi*dkt
            const double weight =
                computeMultiplicityBinWeight(param, m, ik, dkt, a);
            dNdeta2 += n[ik] * weight;
            dEdeta2 += E[ik] * weight;
            if (ik * dkt / a * hbarc > 3.)  //
            {
                dNdetaCut += n[ik] * weight;
                dEdetaCut += E[ik] * weight;
            }
            if (ik * dkt / a * hbarc > 6.)  // large cut
            {
                dNdetaCut2 += n[ik] * weight;
                dEdetaCut2 += E[ik] * weight;
            }
        }
    }

    addPhaseAndRestart(
        "observables.gluon_multiplicity.bin_postprocess",
        multiplicityPhaseStart);

    // compute hadrons using fragmentation function
    if (it == itmax && param->output.writeHadronSpectrum) {
        hadronizeAndWrite(param, a, dkt, bins, n, Nhgsl, hbins);
        multiplicityPhaseStart = ipg::wallSeconds();
    }

    if (!param->colorCharge.usePseudoRapidity && param->run.MPIRank == 0) {
        messager_ << "[GluonMultiplicity::compute]: dN/dy 1 = " << dNdeta
                  << ", dE/dy 1 = " << dEdeta;
        messager_.flush("info");
        messager_ << "[GluonMultiplicity::compute]: dN/dy 2 = " << dNdeta2
                  << ", dE/dy 2 = " << dEdeta2;
        messager_.flush("info");
        messager_ << "[GluonMultiplicity::compute]: gluon <p_T> = "
                  << dEdeta / dNdeta;
        messager_.flush("info");
    } else if (param->colorCharge.usePseudoRapidity) {
        m = param->colorCharge.jacobianMass;                           // in GeV
        P = 0.13 + 0.32 * pow(param->collision.sqrtS / 1000., 0.115);  // in GeV
        dNdeta *= cosh(param->colorCharge.rapidity())
                  / (sqrt(
                      pow(cosh(param->colorCharge.rapidity()), 2.)
                      + m * m / (P * P)));
        dEdeta *= cosh(param->colorCharge.rapidity())
                  / (sqrt(
                      pow(cosh(param->colorCharge.rapidity()), 2.)
                      + m * m / (P * P)));

        if (param->run.MPIRank == 0) {
            messager_ << "[GluonMultiplicity::compute]: dN/deta 1 = " << dNdeta
                      << ", dE/deta 1 = " << dEdeta;
            messager_.flush("info");
            messager_ << "[GluonMultiplicity::compute]: dN/deta 2 = " << dNdeta2
                      << ", dE/deta 2 = " << dEdeta2;
            messager_.flush("info");
            messager_ << "[GluonMultiplicity::compute]: dN/deta_cut 1 = "
                      << dNdetaCut;
            messager_.flush("info");
            messager_ << "[GluonMultiplicity::compute]: dN/deta_cut 2 = "
                      << dNdetaCut2;
            messager_.flush("info");
            messager_ << "[GluonMultiplicity::compute]: gluon <p_T> = "
                      << dEdeta / dNdeta;
            messager_.flush("info");
        }
    }

    addPhaseAndRestart(
        "observables.gluon_multiplicity.report", multiplicityPhaseStart);

    if (dNdeta == 0.) {
        messager_ << "[GluonMultiplicity::compute]: No collision happened on "
                     "rank "
                  << param->run.MPIRank
                  << ". Restarting with new random number...";
        messager_.flush("warning");
        addPhaseAndRestart(
            "observables.gluon_multiplicity.cleanup", multiplicityPhaseStart);
        return 0;
    }

    if (it == itmax) {
        ofstream foutNN(NpartdNdy_name.c_str(), std::ios::out);
        foutNN << param->event.Npart << " " << dNdeta << " " << param->event.Tpp
               << " " << param->event.b << " " << dEdeta << " "
               << param->run.randomSeed << " "
               << "N/A"
               << " "
               << "N/A"
               << " "
               << "N/A"
               << " " << dNdetaCut << " " << dEdetaCut << " " << dNdetaCut2
               << " " << dEdetaCut2 << " "
               << computeRunningCouplingGfactorFromScale(
                      g, muZero, c, param->coupling.LambdaQCD,
                      param->coupling.nFlavors,
                      param->coupling.runningCouplingQsFactor
                          * param->event.averageQs)
               << endl;
        foutNN.close();
        writeTarget(
            param, it, a, dtau, dNdeta, dNdeta2, dEdeta, dEdeta2, dNdetaCut,
            dEdetaCut, dNdetaCut2, dEdetaCut2, n, E, counter, bins, dkt);
    }
    addPhaseAndRestart("output.gluon_multiplicity", multiplicityPhaseStart);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.cleanup", multiplicityPhaseStart);

    messager_ << "[GluonMultiplicity::compute]:  done.";
    messager_.flush("info");
    param->event.success = 1;
    return 1;
}

void GluonMultiplicity::hadronizeAndWrite(
    Parameters *param, double a, double dkt, int bins, const double *n,
    double *Nhgsl, int hbins) {
    const double hadronizationStart = ipg::wallSeconds();
    messager_ << "[GluonMultiplicity::compute]:  Hadronizing ... ";
    messager_.flush("info");
    double z, frac;
    double mypt, kt, Ng;
    int ik;
    const int steps = 6000;
    double dz = 0.95 / static_cast<double>(steps);
    double zValues[steps + 1];
    double zintegrand[steps + 1];
    gsl_interp_accel *zacc = gsl_interp_accel_alloc();
    gsl_spline *zspline = gsl_spline_alloc(gsl_interp_cspline, steps + 1);

    for (int ih = 0; ih <= hbins; ih++) {
        mypt = ih * (20. / static_cast<double>(hbins));

        for (int iz = 0; iz <= steps; iz++) {
            z = 0.05 + iz * dz;
            zValues[iz] = z;

            kt = mypt / z;

            ik = static_cast<int>(
                floor(kt * a / hbarc / dkt - 0.5 + 0.00000001));

            frac = (kt - (ik + 0.5) * dkt / a * hbarc) / (dkt / a * hbarc);

            if (ik + 1 < bins && ik >= 0)
                Ng = ((1. - frac) * n[ik] + frac * n[ik + 1]) * a / hbarc * a
                     / hbarc;  // to make dN/d^2k_T fo k_T in GeV
            else
                Ng = 0.;

            if (!param->colorCharge.usePseudoRapidity) {
                zintegrand[iz] = 1. / (z * z) * Ng * kkp(7, 1, z, kt);
            } else {
                zintegrand[iz] =
                    1. / (z * z) * Ng * 2.
                    * (kkp(1, 1, z, kt) * cosh(param->colorCharge.rapidity())
                           / (sqrt(
                               pow(cosh(param->colorCharge.rapidity()), 2.)
                               + m_pion * m_pion / (mypt * mypt)))
                       + kkp(2, 1, z, kt) * cosh(param->colorCharge.rapidity())
                             / (sqrt(
                                 pow(cosh(param->colorCharge.rapidity()), 2.)
                                 + m_kaon * m_kaon / (mypt * mypt)))
                       + kkp(4, 1, z, kt) * cosh(param->colorCharge.rapidity())
                             / (sqrt(
                                 pow(cosh(param->colorCharge.rapidity()), 2.)
                                 + m_proton * m_proton / (mypt * mypt))));
            }
        }

        zValues[steps] = 1.;  // set exactly 1

        gsl_spline_init(zspline, zValues, zintegrand, steps + 1);
        Nhgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);
    }

    gsl_spline_free(zspline);
    gsl_interp_accel_free(zacc);

    stringstream strmultHad_name;
    strmultHad_name << "multiplicityHadrons" << param->event.eventId << ".dat";
    string multHad_name;
    multHad_name = strmultHad_name.str();

    ofstream foutdNdpt(multHad_name.c_str(), std::ios::out);
    for (int ih = 0; ih <= hbins; ih++) {
        if (ih % 10 == 0)
            foutdNdpt << ih * 20. / static_cast<double>(hbins) << " "
                      << Nhgsl[ih] << " " << 0. << " " << 0. << " "
                      << param->event.Tpp << " " << param->event.b
                      << endl;  // leaving out the L and H ones for now
    }
    foutdNdpt.close();

    messager_ << "[GluonMultiplicity::compute]:  done.";
    messager_.flush("info");

    ipg::Profiler::instance().add(
        "observables.gluon_multiplicity.hadronization",
        ipg::wallSeconds() - hadronizationStart);
}

void GluonMultiplicity::writeTarget(
    Parameters *param, int it, double a, double dtau, double dNPrimary,
    double dNBinned, double dEPrimary, double dEBinned, double dNCut3,
    double dECut3, double dNCut6, double dECut6, const double *spectrumN,
    const double *spectrumE, const int *spectrumCounts, int bins, double dkt) {
    IPG_PROFILE_SCOPE("output.gluon_target");
    stringstream filename;
    filename << "gluonMultiplicity" << param->event.eventId << ".json";
    ofstream output(filename.str().c_str(), std::ios::out | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "could not open gluon-multiplicity target " + filename.str());
    }

    const char *rapidityVariable =
        (!param->colorCharge.usePseudoRapidity) ? "y" : "eta";
    const double meanKt = (dNPrimary != 0.0) ? dEPrimary / dNPrimary : 0.0;
    const double spectrumUnitFactor = (a / hbarc) * (a / hbarc);

    output << std::setprecision(17) << "{\n"
           << "  \"format\": \"ipglasma-gluon-target\",\n"
           << "  \"version\": 1,\n"
           << "  \"event_id\": " << param->event.eventId << ",\n"
           << "  \"step\": " << it << ",\n"
           << "  \"tau_fm\": " << static_cast<double>(it) * dtau * a << ",\n"
           << "  \"rapidity_variable\": \"" << rapidityVariable << "\",\n"
           << "  \"dN\": " << dNPrimary << ",\n"
           << "  \"dN_binned_check\": " << dNBinned << ",\n"
           << "  \"dE_GeV\": " << dEPrimary << ",\n"
           << "  \"dE_binned_check_GeV\": " << dEBinned << ",\n"
           << "  \"mean_kT_GeV\": " << meanKt << ",\n"
           << "  \"dN_kT_gt_3_GeV\": " << dNCut3 << ",\n"
           << "  \"dE_kT_gt_3_GeV\": " << dECut3 << ",\n"
           << "  \"dN_kT_gt_6_GeV\": " << dNCut6 << ",\n"
           << "  \"dE_kT_gt_6_GeV\": " << dECut6 << ",\n"
           << "  \"Npart\": " << param->event.Npart << ",\n"
           << "  \"Tpp\": " << param->event.Tpp << ",\n"
           << "  \"impact_parameter_fm\": " << param->event.b << ",\n"
           << "  \"random_seed\": " << param->run.randomSeed << ",\n"
           << "  \"spectrum_definition\": \"azimuthally averaged Coulomb-gauge "
              "gluon spectrum used by GluonMultiplicity::compute\",\n"
           << "  \"kt_GeV\": [";

    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << (static_cast<double>(ik) + 0.5) * dkt / a * hbarc;
    }
    output << "],\n  \"dN_d2k_GeV_minus2\": [";
    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << spectrumN[ik] * spectrumUnitFactor;
    }
    output << "],\n  \"dE_d2k_GeV_minus1\": [";
    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << spectrumE[ik] * spectrumUnitFactor;
    }
    output << "],\n  \"lattice_bin_counts\": [";
    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << spectrumCounts[ik];
    }
    output << "]\n}\n";
    output.close();

    if (!output) {
        throw std::runtime_error(
            "failed while writing gluon-multiplicity target " + filename.str());
    }
    messager_ << "[GluonMultiplicity::writeTarget]: Wrote gluon "
                 "target dN/d"
              << rapidityVariable << "=" << dNPrimary << " to "
              << filename.str();
    messager_.flush("info");
}

void GluonMultiplicity::readNkt(Parameters *param) {
    // static: no FFT needed; a local log sink and spectrum buffer
    PrettyOstream messager;
    double nIn[100] = {};
    messager << "[GluonMultiplicity::readNkt]: Reading n(k_T) from file ";
    messager.flush("info");
    string Npart, dummy;
    string kt, nkt, Tpp, b;
    double dkt = 0.;
    double dNdeta = 0.;

    // open file

    ifstream fin;
    stringstream strmult_name;
    strmult_name << "multiplicity" << param->event.eventId << ".dat";
    string mult_name;
    mult_name = strmult_name.str();
    fin.open(mult_name.c_str());
    messager << "[GluonMultiplicity::readNkt]: File " << mult_name.c_str()
             << " ... ";
    messager.flush("info");

    // open file

    ifstream fin2;
    stringstream strmult_name2;
    strmult_name2 << "NpartdNdy" << param->event.eventId << ".dat";
    string mult_name2;
    mult_name2 = strmult_name2.str();
    fin2.open(mult_name2.c_str());
    messager << "[GluonMultiplicity::readNkt]: File " << mult_name2.c_str()
             << " ... ";
    messager.flush("info");

    // read file

    if (fin) {
        for (int ikt = 0; ikt < 100; ikt++) {
            if (!fin.eof()) {
                fin >> dummy;
                fin >> kt;
                fin >> nkt;
                nIn[ikt] = atof(nkt.c_str());
                fin >> dummy >> Tpp >> b >> Npart;
                if (ikt == 0) dkt = atof(kt.c_str());
                if (ikt == 1) dkt = dkt - atof(kt.c_str());
            }
            messager << "[GluonMultiplicity::readNkt]: " << nIn[ikt];
            messager.flush("info");
        }
        fin.close();
        messager << "[GluonMultiplicity::readNkt]:  done.";
        messager.flush("info");
    } else {
        messager << "[GluonMultiplicity::readNkt]: File " << mult_name.c_str()
                 << " does not exist. Exiting.";
        messager.flush("error");
        exit(1);
    }

    if (fin2) {
        if (!fin2.eof()) {
            fin2 >> Npart;
            fin2 >> nkt;
            fin2 >> Tpp;
            fin2 >> b;

            dNdeta = atof(nkt.c_str());
        }
        fin2.close();
        messager << "[GluonMultiplicity::readNkt]:  done.";
        messager.flush("info");
    } else {
        messager << "[GluonMultiplicity::readNkt]: File " << mult_name2.c_str()
                 << " does not exist. Exiting.";
        messager.flush("error");
        exit(1);
    }

    double m, P;
    m = param->colorCharge.jacobianMass;                           // in GeV
    P = 0.13 + 0.32 * pow(param->collision.sqrtS / 1000., 0.115);  // in GeV
    double dNdeta2;
    dNdeta2 = 0.;

    for (int ik = 0; ik < 100; ik++) {
        if (!param->colorCharge.usePseudoRapidity) {
            dNdeta2 += nIn[ik] * (ik + 0.5) * dkt * dkt * 2.
                       * M_PI;  // integrate, gives a ik*dkt*2pi*dkt
        } else {
            dNdeta2 +=
                nIn[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI
                * cosh(param->colorCharge.rapidity())
                / (sqrt(
                    pow(cosh(param->colorCharge.rapidity()), 2.)
                    + m * m / (((ik + 0.5) * dkt) * ((ik + 0.5) * dkt))));
        }
    }

    dNdeta *=
        cosh(param->colorCharge.rapidity())
        / (sqrt(pow(cosh(param->colorCharge.rapidity()), 2.) + m * m / P / P));

    ofstream foutNN("NpartdNdy-mod.dat", std::ios::out);
    foutNN << Npart << " " << dNdeta << " " << dNdeta2 << " "
           << atof(Tpp.c_str()) << " " << atof(b.c_str()) << endl;
    foutNN.close();

    exit(1);
}

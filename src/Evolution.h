// Init.h is part of the IP-Glasma solver.
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

class Evolution {
  private:
    FFT *fft_;
    double nIn_[100];  // k_T array
    PrettyOstream messager_;

  public:
    // Constructor
    Evolution(const int nn[]) { fft_ = new FFT(nn); }

    ~Evolution() { delete fft_; }

    void run(Lattice *lat, Group *group, Parameters *param);
    // run()'s writeOutputs==3 epsilon-plot writers.
    // NOTE: this one computes its running-coupling factor via only the
    // "local Qs" formula, unconditionally -- unlike
    // writeEpsilonIntermediatePlot and computeRunningCouplingGfactor, it
    // never checks getRunWithLocalQs(). Preserved exactly as found; this
    // looks like a pre-existing inconsistency, not something to silently
    // change while extracting it.
    void writeEpsilonInitialPlot(Lattice *lat, Parameters *param);
    // Uses the same formula as computeRunningCouplingGfactor (an
    // anonymous-namespace free function further up in Evolution.cpp).
    void writeEpsilonIntermediatePlot(Lattice *lat, Parameters *param);
    void evolveU(Lattice *lat, Parameters *param, double dtau, double tau);
    void evolvePhi(Lattice *lat, Parameters *param, double dtau, double tau);
    void evolvePi(Lattice *lat, Parameters *param, double dtau, double tau);
    void evolveE(Lattice *lat, Parameters *param, double dtau, double tau);
    void checkGaussLaw(Lattice *lat, Parameters *param);
    void eccentricity(
        Lattice *lat, Parameters *param, int it, double cutoff, int doAniso);
    void tmunu(Lattice *lat, Parameters *param, int it);
    void writeEvolvedFields(Lattice *lat, Parameters *param, int it);
    void u(Lattice *lat, Parameters *param, int it, bool finalFlag);
    int multiplicity(Lattice *lat, Group *group, Parameters *param, int it);
    // multiplicity()'s writeOutputs==3 hadronization step: convolves the
    // binned gluon spectrum n[] with the KKP fragmentation function to get
    // a hadron p_T spectrum, writing multiplicityHadrons<id>.dat. Nhgsl
    // (length hbins+1) is scratch/output space owned by the caller.
    void hadronizeAndWriteMultiplicity(
        Parameters *param, double a, double dkt, int bins, const double *n,
        double *Nhgsl, int hbins);

    void writeGluonMultiplicityTarget(
        Parameters *param, int it, double a, double dtau, double dNPrimary,
        double dNBinned, double dEPrimary, double dEBinned, double dNCut3,
        double dECut3, double dNCut6, double dECut6, const double *spectrumN,
        const double *spectrumE, const int *spectrumCounts, int bins,
        double dkt);
    void readNkt(Parameters *param);
};

#endif  // SRC_EVOLUTION_H_

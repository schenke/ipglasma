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

    void writeGluonMultiplicityTarget(
        Parameters *param, int it, double a, double dtau, double dNPrimary,
        double dNBinned, double dEPrimary, double dEBinned, double dNCut3,
        double dECut3, double dNCut6, double dECut6, const double *spectrumN,
        const double *spectrumE, const int *spectrumCounts, int bins,
        double dkt);
    void readNkt(Parameters *param);
};

#endif  // SRC_EVOLUTION_H_

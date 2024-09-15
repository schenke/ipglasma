// Init.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef Evolution_H
#define Evolution_H

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "FFT.h"
#include "GaugeFix.h"
#include "Glauber.h"
#include "Group.h"
#include "Lattice.h"
#include "Matrix.h"
#include "MyEigen.h"
#include "Parameters.h"
#include "Random.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

class Evolution {
private:
  FFT *fft;
  double nIn[100]; // k_T array

public:
  // Constructor
  Evolution(const int nn[]) { fft = new FFT(nn); }

  ~Evolution() { delete fft; }

  void run(Lattice *lat, BufferLattice *bufferlat, Group *group,
           Parameters *param);
  void evolveU(Lattice *lat, BufferLattice *bufferlat, Parameters *param,
               double dtau, double tau);
  void evolveUfast(Lattice *lat, Group *group, Parameters *param, double dtau,
                   double tau);
  void evolvePhi(Lattice *lat, BufferLattice *bufferlat, Parameters *param,
                 double dtau, double tau);
  void evolvePi(Lattice *lat, BufferLattice *bufferlat, Parameters *param,
                double dtau, double tau);
  void evolveE(Lattice *lat, BufferLattice *bufferlat, Parameters *param,
               double dtau, double tau);
  void checkGaussLaw(Lattice *lat, Parameters *param);
  void eccentricity(Lattice *lat, Parameters *param, int it, double cutoff,
                    int doAniso);
  void Tmunu(Lattice *lat, Parameters *param, int it);
  void u(Lattice *lat, Parameters *param, int it, bool finalFlag);
  int multiplicity(Lattice *lat, Group *group, Parameters *param, int it);
  int multiplicitynkxky(Lattice *lat, Group *group, Parameters *param, int it);
  int correlations(Lattice *lat, Group *group, Parameters *param, int it);
  void anisotropy(Lattice *lat, Parameters *param, int it);
  void readNkt(Parameters *param);
};

#endif // Evolution_H

// Init.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef Evolution_H
#define Evolution_H

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <unistd.h>

#include "Lattice.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Random.h"
#include "Group.h"
#include "FFT.h"
#include "Glauber.h"
#include "GaugeFix.h"
#include "MyEigen.h"
#include "Fragmentation.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>


using namespace std;

class Evolution {

 private:
  FFT *fft;
  Fragmentation  *frag;
  
  double nIn[100]; //k_T array


 public:
  
  // Constructor.
  Evolution(const int nn[]) 
    {
      fft = new FFT(nn);
      frag = new Fragmentation();
    };
  
  ~Evolution() 
    { 
      delete fft;
      delete frag;
    };
  
  void run(Lattice* lat, BufferLattice* bufferlat, Group* group, Parameters *param);
  void evolveU(Lattice* lat, BufferLattice *bufferlat, Group* group, Parameters *param, double dtau, double tau);
  void evolveUfast(Lattice* lat, Group* group, Parameters *param, double dtau, double tau);
  void evolvePhi(Lattice* lat, BufferLattice *bufferlat, Group* group, Parameters *param, double dtau, double tau);
  void evolvePi(Lattice* lat, BufferLattice * bufferlat, Group* group, Parameters *param, double dtau, double tau);
  void evolveE(Lattice* lat, BufferLattice *bufferlat, Group* group, Parameters *param, double dtau, double tau);
  void checkGaussLaw(Lattice* lat, Group* group, Parameters *param, double dtau, double tau);
  void eccentricity(Lattice *lat, Group *group, Parameters *param, int it, double cutoff, int doAniso);
  void Tmunu(Lattice *lat, Group *group, Parameters *param, int it);
  void u(Lattice *lat, Group *group, Parameters *param, int it);
  int multiplicity(Lattice *lat, Group *group, Parameters *param, int it);
  int multiplicitynkxky(Lattice *lat, Group *group, Parameters *param, int it);
  int correlations(Lattice *lat, Group *group, Parameters *param, int it);
  int correlationsColor(Lattice *lat, Group *group, Parameters *param, int it);
  void anisotropy(Lattice *lat, Group *group, Parameters *param, int it);
  void readNkt(Parameters *param);
  void twoPointFunctionInK(Parameters *param, Lattice *lat, int ids);
};

#endif // Evolution_H

// Init.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef Init_H
#define Init_H

#include <complex>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FFT.h"
#include "GaugeFix.h"
#include "Glauber.h"
#include "Group.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "Random.h"
#include "gsl/gsl_linalg.h"
#include "pretty_ostream.h"

class Init {

private:
  int const static iymaxNuc = 44; // for the Tp-y table

  int const static iTpmax = 200; // updated in March 2019 to a larger T_A range

  double const deltaYNuc = 0.25; // for the new table

  FFT fft;
  //  Matrix** A;
  //  Glauber *glauber;
  double Qs2Nuclear[iTpmax][iymaxNuc];
  double Tlist[iTpmax];

  double As[1];

  std::vector<ReturnValue>
      nucleusA; // list of x and y coordinates of nucleons in nucleus A
  std::vector<ReturnValue>
      nucleusB; // list of x and y coordinates of nucleons in nucleus B

  pretty_ostream messager;

public:
  // Constructor.
  Init(const int nn[]) : fft(nn){};

  ~Init(){};

  void init(Lattice *lat, Group *group, Parameters *param, Random *random,
            Glauber *glauber, int READFROMFILE);
  void sampleTA(Parameters *param, Random *random, Glauber *glauber);
  void readNuclearQs(Parameters *param);
  std::vector<complex<double>> solveAxb(Parameters *param, complex<double> *A,
                                        complex<double> *b);
  double getNuclearQs2(double Qs2atZeroY, double y);
  void setColorChargeDensity(Lattice *lat, Parameters *param, Random *random,
                             Glauber *glauber);
  void setV(Lattice *lat, Group *group, Parameters *param, Random *random);
  void readV(Lattice *lat, Parameters *param);
  // void eccentricity(Lattice *lat, Group *group, Parameters *param, Random
  // *random, Glauber *glauber);
  void multiplicity(Lattice *lat, Parameters *param);

  void generate_nucleus_configuration(Random *random, int A, int Z, double a_WS,
                                      double R_WS, double beta2, double beta4,
                                      std::vector<ReturnValue> *nucleus);
  void generate_nucleus_configuration_with_woods_saxon(
      Random *random, int A, int Z, double a_WS, double R_WS,
      std::vector<ReturnValue> *nucleus);
  void generate_nucleus_configuration_with_deformed_woods_saxon(
      Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
      double beta4, std::vector<ReturnValue> *nucleus);
  double sample_r_from_woods_saxon(Random *random, double a_WS,
                                   double R_WS) const;
  void sample_r_and_costheta_from_deformed_woods_saxon(Random *random,
                                                       double a_WS, double R_WS,
                                                       double beta2,
                                                       double beta4, double &r,
                                                       double &costheta) const;
  double fermi_distribution(double r, double R_WS, double a_WS) const;
  double spherical_harmonics(int l, double ct) const;
  void recenter_nucleus(std::vector<double> &x, std::vector<double> &y,
                        std::vector<double> &z);
  void rotate_nucleus(double phi, double theta, std::vector<double> &x,
                      std::vector<double> &y, std::vector<double> &z);

  void samplePartonPositions(Parameters *param, Random *random,
                             std::vector<double> &x_array,
                             std::vector<double> &y_array,
                             std::vector<double> &z_array,
                             std::vector<double> &BGq_array);

  double sampleLogNormalDistribution(Random *random, const double mean,
                                     const double variance);

  void sampleQsNormalization(Random *random, Parameters *param,
                             const int Nq,
                             std::vector<double> &gauss_array);
};

#endif  // Init_H

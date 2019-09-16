// Init.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef Init_H
#define Init_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <limits> 
#include <ctime>

#include "mpi.h"
#include "Lattice.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Random.h"
#include "Group.h"
#include "FFT.h"
#include "Glauber.h"
#include "GaugeFix.h"
#include "gsl/gsl_linalg.h"
#include "pretty_ostream.h"

using namespace std;

class Init {

 private:
  int const static iymaxNuc = 44; // for the Tp-y table

  int const static iTpmax = 200; // updated in March 2019 to a larger T_A range
 
  double const deltaYNuc = 0.25; // for the new table
  
  FFT *fft;
  //  Matrix** A;
  //  Glauber *glauber;
  double Qs2Nuclear[iTpmax][iymaxNuc];
  double Tlist[iTpmax];
  
  double As[1];
  
  vector<ReturnValue> nucleusA;  // list of x and y coordinates of nucleons in nucleus A      
  vector<ReturnValue> nucleusB;  // list of x and y coordinates of nucleons in nucleus B 

    pretty_ostream messager;

 public:
  
  // Constructor.

  Init(const int nn[]) 
    {
      //  random = new Random;
      fft = new FFT(nn);
      //      glauber = new Glauber;
    };
  
  ~Init() 
    { 
      delete fft; 
      //delete glauber;
      //delete random;
    };
  
  void init(Lattice *lat, Group *group, Parameters *param, Random *random, Glauber* glauber, int READFROMFILE);
  void sampleTA(Parameters *param, Random *random, Glauber* glauber);
  void readNuclearQs(Parameters *param);
  vector <complex<double> > solveAxb(Parameters *param, complex<double>* A, complex<double>* b);
  double getNuclearQs2(Parameters *param, Random *random, double Qs2atZeroY, double y);
  void setColorChargeDensity(Lattice *lat, Parameters *param, Random *random, Glauber *glauber);
  void setV(Lattice *lat, Group* group, Parameters *param, Random* random, Glauber *glauber);
  void readV(Lattice *lat, Group* group, Parameters *param);
  // void eccentricity(Lattice *lat, Group *group, Parameters *param, Random *random, Glauber *glauber);
  void multiplicity(Lattice *lat, Group *group, Parameters *param, Random *random, Glauber *glauber);

    void generate_nucleus_configuration(
            Parameters *param, Random *random,
            int A, int Z, double a_WS, double R_WS, double beta2, double beta4,
            std::vector<ReturnValue> *nucleus);
    void generate_nucleus_configuration_with_woods_saxon(
            Parameters *param, Random *random,
            int A, int Z, double a_WS, double R_WS,
            std::vector<ReturnValue> *nucleus);
    void generate_nucleus_configuration_with_deformed_woods_saxon(
           Parameters *param, Random *random,
           int A, int Z, double a_WS, double R_WS, double beta2, double beta4,
           std::vector<ReturnValue> *nucleus);
    double sample_r_from_woods_saxon(Random *random, double a_WS, double R_WS) const;
    void sample_r_and_costheta_from_deformed_woods_saxon(
            Random *random, double a_WS, double R_WS, double beta2, double beta4,
            double &r, double &costheta) const;
    double fermi_distribution(double r, double R_WS, double a_WS) const;
    double spherical_harmonics(int l, double ct) const;
    void recenter_nucleus(std::vector<double> &x, std::vector<double> &y,
                          std::vector<double> &z);
    void rotate_nucleus(double phi, double theta,
                        std::vector<double> &x, std::vector<double> &y,
                        std::vector<double> &z);

};

#endif // Init_H

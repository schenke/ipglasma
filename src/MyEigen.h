// MyEigen.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef MyEigen_H
#define MyEigen_H

#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_eigen.h"
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"

using namespace std;

class MyEigen {
private:
public:
  // Constructor.
  MyEigen(){};

  ~MyEigen(){};
  void test();
  void flowVelocity4D(Lattice *lat, Parameters *param, int it, bool finalFlag);
};

#endif // MyEigen_H

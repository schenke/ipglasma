// MyEigen.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef MyEigen_H
#define MyEigen_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>

#include "gsl/gsl_eigen.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"

#include "Lattice.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Group.h"

using namespace std;

class MyEigen {

 private:
 public:
  
  // Constructor.
  MyEigen() 
    {
    };
  
  ~MyEigen() 
    { 
    };
  
  void test();
  void flowVelocity4D(Lattice *lat, Group *group, Parameters *param, int it);


};

#endif // MyEigen_H

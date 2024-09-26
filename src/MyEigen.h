// MyEigen.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef MyEigen_H
#define MyEigen_H

#include "Lattice.h"
#include "Parameters.h"

class MyEigen {
  private:
  public:
    // Constructor.
    MyEigen() {};

    ~MyEigen() {};
    void test();
    void flowVelocity4D(
        Lattice *lat, Parameters *param, int it, bool finalFlag);
};

#endif  // MyEigen_H

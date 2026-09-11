// MyEigen.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_MYEIGEN_H_
#define SRC_MYEIGEN_H_

#include "Lattice.h"
#include "Parameters.h"
#include "PrettyOstream.h"

class MyEigen {
  private:
    PrettyOstream messager_;

  public:
    // Constructor.
    MyEigen() {};

    ~MyEigen() {};
    void flowVelocity4D(
        Lattice *lat, Parameters *param, int it, bool finalFlag);
    void writeTmunu4D(Lattice *lat, Parameters *param, int it);

  private:
    void flowVelocity4DImpl(
        Lattice *lat, Parameters *param, int it, bool finalFlag,
        bool tmunuOnly);
};

#endif  // SRC_MYEIGEN_H_

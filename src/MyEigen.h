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

    // Writes the hydro-flow text output (epsilon-u-Hydro*.dat): per-cell
    // interpolated energy density, u^mu, and pi^munu on the output grid.
    // Returns the total energy Etot integrated over the grid, or 0 if this
    // output is disabled (param->getWriteOutputs() % 2 != 1) -- writeJazma
    // needs Etot even when this writer itself is off.
    double writeHydroText(
        Lattice *lat, Parameters *param, int it, bool finalFlag,
        bool tmunuOnly, int N, double L, double a, double dtau,
        double gfactor, int hx, int hy, int heta, double hL, double deta,
        double ha, double tau0);

    // Writes the raw/binary Tmunu output (Tmunu-t*.dat or *.ipgt):
    // per-cell interpolated T^munu components on the output grid.
    void writeRawTmunu(
        Lattice *lat, Parameters *param, int it, int N, double L, double a,
        double dtau, double gfactor, int hx, int hy, int heta, double hL,
        double deta, double ha, double tau0);

    // Writes the Jazma output (Jazma-Hydro-t*.dat): per-cell g2mu2A*g2mu2B,
    // normalized so its grid integral matches Etot (from writeHydroText).
    void writeJazma(
        Lattice *lat, Parameters *param, int it, double Etot, int N,
        double L, double a, double dtau, int hx, int hy, int heta, double hL,
        double deta, double ha);
};

#endif  // SRC_MYEIGEN_H_

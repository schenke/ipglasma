// GaugeFix.h is part of the CYM solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef GaugeFix_H
#define GaugeFix_H

#include "FFT.h"
#include "Group.h"
#include "Lattice.h"
#include "Parameters.h"

class GaugeFix {
  private:
  public:
    // Constructor.
    GaugeFix() {};

    // Destructor.
    ~GaugeFix() {};

    void gaugeTransform(Lattice *lat, Parameters *param, int i, int j);
    void FFTChi(
        FFT *fft, Lattice *lat, Group *group, Parameters *param, int steps);
};

#endif  // GaugeFix_H

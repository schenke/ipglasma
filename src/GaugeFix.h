// GaugeFix.h is part of the CYM solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_GAUGEFIX_H_
#define SRC_GAUGEFIX_H_

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

    void FFTChi(
        FFT *fft, Lattice *lat, Group *group, Parameters *param, int steps);
};

#endif  // SRC_GAUGEFIX_H_

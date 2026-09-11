// GaugeFix.h is part of the CYM solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_GAUGEFIX_H_
#define SRC_GAUGEFIX_H_

#include "FFT.h"
#include "Group.h"
#include "Lattice.h"
#include "Parameters.h"
#include "PrettyOstream.h"

class GaugeFix {
  private:
    PrettyOstream messager_;

  public:
    // Constructor.
    GaugeFix() {};

    // Destructor.
    ~GaugeFix() {};

    void fftChi(
        FFT *fft, Lattice *lat, Group *group, Parameters *param, int steps);
};

#endif  // SRC_GAUGEFIX_H_

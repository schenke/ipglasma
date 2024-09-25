// GaugeFix.h is part of the CYM solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef GaugeFix_H
#define GaugeFix_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "FFT.h"
#include "Glauber.h"
#include "Group.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "Random.h"
#include "Spinor.h"

using namespace std;

class GaugeFix {
  private:
    // FFT *fft;
    // Random *random;

  public:
    // Constructor.
    GaugeFix() {
        //  fft = new FFT(nn);
        // random = new Random();
    };

    // Destructor.
    ~GaugeFix() {
        // delete fft;
        // delete random;
    };

    void gaugeTransform(Lattice *lat, Parameters *param, int i, int j);
    void FFTChi(
        FFT *fft, Lattice *lat, Group *group, Parameters *param, int steps);
};

#endif  // GaugeFix_H

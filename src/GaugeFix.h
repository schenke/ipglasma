// GaugeFix.h is part of the CYM solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef GaugeFix_H
#define GaugeFix_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>

#include "Lattice.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Spinor.h"
#include "Random.h"
#include "Group.h"
#include "FFT.h"
#include "Glauber.h"

using namespace std;

class GaugeFix {
 private:
    //FFT *fft;
    //Random *random;

 public:
    // Constructor.
    GaugeFix() {
        //  fft = new FFT(nn);
        //random = new Random();
    };

    // Destructor.
    ~GaugeFix() {
        //delete fft;
        //delete random;
    };

    void gaugeTransform(Lattice* lat, Parameters *param, int i, int j);
    void FFTChi(FFT* fft, Lattice* lat, Group* group, Parameters *param,
                int steps);
};

#endif // GaugeFix_H

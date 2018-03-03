// Fragmentation.h is part of the IP-Glasma solver.
// Copyright (C) 2013 Bjoern Schenke.

#ifndef Fragmentation_H
#define Fragmentation_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>

using namespace std;

class Fragmentation {

public:

  // Constructor.
  Fragmentation() {}

  double kkp(int ih, int iset, double x, double qs);
 
};

#endif // Fragmentation_H

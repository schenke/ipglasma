#ifndef Lattice_h
#define Lattice_h

#include <vector>

#include "Cell.h"
#include "Parameters.h"

// The Lattice class is a level higher than the Cell class
// It takes care of the overall structure of the lattice and the arrangement of
// individual cells "cells" is an array of pointers to individual cells of the
// lattice. The values of the quantities in a cell can be modified or retrieved
// by the public functions in both lattice and cells.

class Lattice {
  private:
    int size;  // the total number of cells (length*length)
    int Nc;  // the number of colors in SU(Nc): Determines the dimension of the
             // used matrices

  public:
    // constructor
    Lattice(Parameters *param, int N, int length);
    // destructor
    ~Lattice();

    // functions to access values within individual cells
    int getSize() { return size; };

    std::vector<Cell *> cells;  // the actual array of cells, the "lattice".
                                // cells is an array of pointers to cell objects

    std::vector<int> posmX;
    std::vector<int> pospX;
    std::vector<int> posmY;
    std::vector<int> pospY;
};

class BufferLattice {
  private:
    int size;  // the total number of cells (length*length)
    int Nc;  // the number of colors in SU(Nc): Determines the dimension of the
             // used matrices

  public:
    // constructor
    BufferLattice(int N, int length);
    // destructor
    ~BufferLattice();

    std::vector<SmallCell *> cells;
    // the actual array of cells, the "lattice".
    // cells is an array of pointers to cell objects
};

#endif

#include "Lattice.h"

// constructor
Lattice::Lattice(Parameters *param, int N, int length) {
  Nc = N;
  size = length * length;
  double a = param->getL() / static_cast<double>(length);

  cout << "Allocating square lattice of size " << length << "x" << length
       << " with a=" << a << " fm ...";

  // initialize the array of cells
  for (int i = 0; i < size; i++) {
    Cell *cell;
    cell = new Cell(Nc);
    cells.push_back(cell);
  }

  for (int i = 0; i < length; i++) {
    for (int j = 0; j < length; j++) {
      // pos = i*length+j;
      pospX.push_back((std::min(length - 1, i + 1)) * length + j);
      pospY.push_back(i * length + std::min(length - 1, j + 1));

      posmX.push_back((std::max(0, i - 1)) * length + j);
      posmY.push_back(i * length + std::max(0, j - 1));

      posmXpY.push_back((std::max(0, i - 1)) * length
                        + std::min(length - 1, j + 1));

      pospXmY.push_back((std::min(length - 1, i + 1)) * length
                        + std::max(0, j - 1));
    }
  }
  cout << " done on rank " << param->getMPIRank() << "." << endl;
}

Lattice::~Lattice() {
  for (int i = 0; i < size; i++)
    delete cells[i];
  cells.clear();
}

// constructor
BufferLattice::BufferLattice(int N, int length) {
  Nc = N;
  size = length * length;

  for (int i = 0; i < size; i++) {
    SmallCell *cell;
    cell = new SmallCell(Nc);
    cells.push_back(cell);
  }
}

BufferLattice::~BufferLattice() {
  for (int i = 0; i < size; i++)
    delete cells[i];
  cells.clear();
}

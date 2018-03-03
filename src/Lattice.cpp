#include "Lattice.h"

//constructor
Lattice::Lattice(Parameters *param, int N, int length)
{
  Nc = N;
  size = length*length;
  double a = param->getL()/static_cast<double>(length);

  cout << "Allocating square lattice of size " << length << "x" << length << " with a=" << a << " fm ...";

  // initialize the array of cells
  //  cells = new Cell*[size];
  for(int i=0; i<size; i++)
    {
      Cell* cell;
      cell = new Cell(Nc);
      cells.push_back(cell);
    }

//   // x (i) is the outer loop, y (j) the inner
//   for (int i=0; i<length; i++)
//     for (int j=0; j<length; j++)
//       {
// 	int pos = i*length+j;
// 	cells[pos]->setParity((i+j)%2); // set parity (like a checker board)
//       }

  cout << " done on rank " << param->getMPIRank() << "." << endl;
}

Lattice::~Lattice()
{
  for(int i=0; i<size; i++)
    delete cells[i];
  cells.clear();
  //delete[] cells;
}


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

  for (int i=0; i<length; i++)
    {
      for (int j=0; j<length; j++)
	{
          // pos = i*length+j;
          pospX.push_back(((i+1)%length)*length+j);
          pospY.push_back(i*length+(j+1)%length);
          
          if(i>0)
	    posmX.push_back((i-1)*length+j);
	  else
	    posmX.push_back((length-1)*length+j);
          if(j>0)
	    posmY.push_back(i*length+(j-1));
	  else
	    posmY.push_back(i*length+(length-1));

          if(i>0)
	    posmXpY.push_back((i-1)*length+(j+1)%length);
	  else
	    posmXpY.push_back((length-1)*length+(j+1)%length);
	  
	  if(j>0)
	    pospXmY.push_back(((i+1)%length)*length+j-1);
	  else
	    pospXmY.push_back(((i+1)%length)*length+length-1);
        }
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
}

//constructor
BufferLattice::BufferLattice(Parameters *param, int N, int length)
{
  Nc = N;
  size = length*length;

  for(int i=0; i<size; i++)
    {
      SmallCell* cell;
      cell = new SmallCell(Nc);
      cells.push_back(cell);
    }
}

BufferLattice::~BufferLattice()
{
  for(int i=0; i<size; i++)
    delete cells[i];
  cells.clear();
}


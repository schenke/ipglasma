#include "Cell.h"

Cell::Cell(int N)
{
  Nc = N;
  g2mu2A = new double;
  g2mu2B = new double;
  TpA = new double;
  TpB = new double;
  U = new Matrix(Nc,1.);
  U2 = new Matrix(Nc,1.);
  Ux = new Matrix(Nc,1.);
  Uy = new Matrix(Nc,1.);
  Ux1 = new Matrix(Nc,1.);
  Uy1 = new Matrix(Nc,1.);
  Ux2 = new Matrix(Nc,1.);
  Uy2 = new Matrix(Nc,1.);
 }

Cell::~Cell()
{
  delete g2mu2A;
  delete g2mu2B;
  delete TpA;
  delete TpB;
  delete U;
  delete U2;
  delete Ux;
  delete Uy;
  delete Ux1;
  delete Ux2;
  delete Uy1;
  delete Uy2;
}

SmallCell::SmallCell(int N)
{
  Nc = N;
  buffer1 = new Matrix(Nc,1.);
  buffer2 = new Matrix(Nc,1.);
}

SmallCell::~SmallCell()
{
  delete buffer1;
  delete buffer2;
}

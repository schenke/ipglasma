#ifndef Group_h
#define Group_h

#include <complex>
#include <iostream>
#include <cstdlib>
#include "Matrix.h"

using namespace std;

class Group
{
 private:
  Matrix** t;  // generators of the group 
  Matrix** tA; // adjoint representation of generators of the group 
  int Nc; //number of colors
  
 public:
  
  //constructor(s)
  Group(int N);
  ~Group();

  Matrix& getT(int i) const { return *t[i]; };
  Matrix& getTA(int i) const { return *tA[i]; };

};
#endif


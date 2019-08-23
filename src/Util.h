#ifndef PI
#define PI (3.14159265358979324)
#endif

#ifndef hbarc
#define hbarc (0.197326938)
#endif

#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <string.h>

using namespace std;

class Util{
 public:
  double *vector_malloc(int ); 
  char *char_malloc(int );
  void char_free(char *);
  void vector_free(double *);
  int IsFile(string );
  string StringFind(string file_name, string st);
  double DFind(string file_name, string st);
  int IFind(string file_name, string st);
};
#endif

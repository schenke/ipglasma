#ifndef UTIL_H
#define UTIL_H

#include <string>

namespace Util {
double *vector_malloc(int);
char *char_malloc(int);
void char_free(char *);
void vector_free(double *);
int IsFile(std::string);
std::string stringFind(std::string file_name, std::string st);
double DFind(std::string file_name, std::string st);
int IFind(std::string file_name, std::string st);
} // namespace Util

#endif

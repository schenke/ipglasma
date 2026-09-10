#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_

#include <string>

namespace Util {
double *vector_malloc(int);
char *char_malloc(int);
void char_free(char *);
void vector_free(double *);
int isFile(std::string);
std::string stringFind(std::string file_name, std::string st);
double dFind(std::string file_name, std::string st);
int iFind(std::string file_name, std::string st);
}  // namespace Util

#endif  // SRC_UTIL_H_

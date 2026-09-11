// Setup.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_SETUP_H_
#define SRC_SETUP_H_

#include <string>
#include <vector>

#include "PrettyOstream.h"

class Setup {
  private:
    PrettyOstream messager_;

  public:
    // Constructor.
    Setup() {}

    std::string stringFind(std::string file_name, std::string st);
    std::string stringFindOptional(
        std::string file_name, std::string st, std::string defaultValue);
    int iFind(std::string file_name, std::string st);
    int iFindOptional(std::string file_name, std::string st, int defaultValue);
    unsigned long long int uLLIFind(std::string file_name, std::string st);
    double dFind(std::string file_name, std::string st);
    int isFile(std::string file_name);
    std::vector<double> listFind(std::string file_name, std::string st);
};

#endif  // SRC_SETUP_H_

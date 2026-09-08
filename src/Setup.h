// Setup.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef Setup_H
#define Setup_H

#include <string>
#include <vector>

class Setup {
  public:
    // Constructor.
    Setup() {}

    std::string StringFind(std::string file_name, std::string st);
    int IFind(std::string file_name, std::string st);
    int IFindOptional(std::string file_name, std::string st, int defaultValue);
    unsigned long long int ULLIFind(std::string file_name, std::string st);
    double DFind(std::string file_name, std::string st);
    int IsFile(std::string file_name);
    std::vector<double> ListFind(std::string file_name, std::string st);
};

#endif  // Setup_H

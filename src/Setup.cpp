// Setup.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Setup.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using std::ifstream;
using std::string;
using std::stringstream;

//**************************************************************************
// Setup class.

//**************************************************************************
// Parameter I/O

// reads a string
string Setup::stringFind(string file_name, string st) {
    string inputname = file_name;
    string tmpfilename;
    string str = st;

    string s;
    string xstr;

    tmpfilename = "input";

    int ind;
    // Check every call rather than gating on a static "already checked"
    // flag: with that flag, a call for a missing file_name after an
    // earlier call already succeeded for a different file_name would skip
    // this whole check, then hang forever reading from a stream that
    // never opened (see the while loop below).
    if (!isFile(file_name)) {
        messager_ << "[Setup::stringFind]: The input file named " << file_name
                  << " is absent. Exiting.";
        messager_.flush("error");
        exit(1);
    }

    ifstream input(inputname.c_str());

    input >> s;

    ind = 0;
    while (s.compare("EndOfFile") != 0) {
        input >> xstr;
        if (s.compare(str) == 0) {
            ind++;
            input.close();
            return xstr;
        } /* if right, return */
        s.clear();
        input >> s;
    } /* while */

    input.close();

    if (ind == 0) {
        messager_ << "[Setup::stringFind]: " << str << " not found in "
                  << inputname << ". Create a complete input file.";
        messager_.flush("error");
        exit(1);
    }
    return "";
} /* stringFind */

std::vector<double> Setup::listFind(string file_name, string st) {
    std::vector<double> varlist;
    if (!isFile(file_name)) {
        messager_ << "[Setup::listFind]: The input file named " << file_name
                  << " is absent. Exiting.";
        messager_.flush("error");
        exit(1);
    }
    ifstream input(file_name.c_str());

    string line;
    while (std::getline(input, line)) {
        if (line.find(st) != string::npos) {
            std::stringstream lineStream(line);
            std::string cell;
            lineStream >> cell;
            while (std::getline(lineStream, cell, ',')) {
                varlist.push_back(std::stod(cell));
            }
        }
    }
    return varlist;
}

// reads a double using stringfind:
double Setup::dFind(string file_name, string st) {
    string s, s2;
    double x;
    stringstream stm;
    s = stringFind(file_name, st);
    stm << s;
    s2 = stm.str();
    x = ::atof(s2.c_str());
    // x << stm;
    return x;
} /* dFind */

// reads an integer using stringfind:
int Setup::iFind(string file_name, string st) {
    double f;
    f = dFind(file_name, st);

    // return (int)(f + 0.5);
    return static_cast<int>(f);
} /* iFind */

int Setup::iFindOptional(string file_name, string st, int defaultValue) {
    ifstream input(file_name.c_str());
    if (!input.is_open()) {
        messager_ << "[Setup::iFindOptional]: The input file named "
                  << file_name << " is absent. Exiting.";
        messager_.flush("error");
        exit(1);
    }

    string key;
    string value;
    while (input >> key) {
        if (key == "EndOfFile") break;
        if (!(input >> value)) break;
        if (key == st) {
            return static_cast<int>(::atof(value.c_str()));
        }
    }
    return defaultValue;
}

string Setup::stringFindOptional(
    string file_name, string st, string defaultValue) {
    ifstream input(file_name.c_str());
    if (!input.is_open()) {
        messager_ << "[Setup::stringFindOptional]: The input file named "
                  << file_name << " is absent. Exiting.";
        messager_.flush("error");
        exit(1);
    }

    string key;
    string value;
    while (input >> key) {
        if (key == "EndOfFile") break;
        if (!(input >> value)) break;
        if (key == st) {
            return value;
        }
    }
    return defaultValue;
}

// reads an integer using stringfind:
unsigned long long int Setup::uLLIFind(string file_name, string st) {
    double f;
    f = dFind(file_name, st);

    return (unsigned long long int)(f + 0.5);
} /* iFind */

int Setup::isFile(string file_name) {
    FILE *temp;

    if ((temp = fopen(file_name.c_str(), "r")) == NULL)
        return 0;
    else {
        fclose(temp);
        return 1;
    }
} /* isFile */

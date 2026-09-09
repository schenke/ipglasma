// Setup.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Setup.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using std::cerr;
using std::cout;
using std::endl;
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
    static int flag = 0;
    if (flag == 0) {
        if (!isFile(file_name)) {
            cerr << "The input file named " << file_name
                 << " is absent. Exiting." << endl;
            exit(1);
        }
        flag = 1;
    } /* if flag == 0 */

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
        cerr << str << " not found in " << inputname << endl;
        cout << "Create a complete input file." << endl;
        // return xstr;
        exit(1);
    }
    return (0);
} /* stringFind */

std::vector<double> Setup::listFind(string fileName, string paramName) {
    std::vector<double> varlist;
    if (!isFile(fileName)) {
        cerr << "The input file named " << fileName << " is absent. Exiting."
             << endl;
        exit(1);
    }
    ifstream input(fileName.c_str());

    string line;
    while (std::getline(input, line)) {
        if (line.find(paramName) != string::npos) {
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
    // cout << "ccheck1" << endl;
    string s, s2;
    double x;
    stringstream stm;
    s = stringFind(file_name, st);
    // cout << "ccheck2" << endl;
    stm << s;
    s2 = stm.str();
    // cout << "ccheck3" << endl;
    x = ::atof(s2.c_str());
    // x << stm;
    // cout << "ccheck4" << endl;
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
        cerr << "The input file named " << file_name << " is absent. Exiting."
             << endl;
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

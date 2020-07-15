// Setup.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Setup.h"

//**************************************************************************
// Setup class.

//**************************************************************************
// Parameter I/O

// reads a string
string Setup::StringFind(string file_name, string st) {
  string inputname = file_name;
  string tmpfilename;
  string str = st;

  string s;
  string xstr;

  tmpfilename = "input";

  int ind;
  static int flag = 0;
  if (flag == 0) {
    if (!IsFile(file_name)) {
      cerr << "The input file named " << file_name << " is absent. Exiting."
           << endl;
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
} /* StringFind */

// reads a double using stringfind:
double Setup::DFind(string file_name, string st) {
  // cout << "ccheck1" << endl;
  string s, s2;
  double x;
  stringstream stm;
  s = StringFind(file_name, st);
  // cout << "ccheck2" << endl;
  stm << s;
  s2 = stm.str();
  // cout << "ccheck3" << endl;
  x = ::atof(s2.c_str());
  // x << stm;
  // cout << "ccheck4" << endl;
  return x;
} /* DFind */

// reads an integer using stringfind:
int Setup::IFind(string file_name, string st) {
  double f;
  f = DFind(file_name, st);

  return (int)(f + 0.5);
} /* IFind */

// reads an integer using stringfind:
unsigned long long int Setup::ULLIFind(string file_name, string st) {
  double f;
  f = DFind(file_name, st);

  return (unsigned long long int)(f + 0.5);
} /* IFind */

int Setup::IsFile(string file_name) {
  FILE *temp;

  if ((temp = fopen(file_name.c_str(), "r")) == NULL)
    return 0;
  else {
    fclose(temp);
    return 1;
  }
} /* IsFile */

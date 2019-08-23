#include "Util.h"
#include <iostream>

double *Util::vector_malloc(int n1)
{
 double *d1_ptr;
 int i;

 /* pointer to the n1 array */
 d1_ptr = new double[n1];
 for(i=0; i<n1; i++) d1_ptr[i] = 0.0; 
    
 return d1_ptr;
}


void Util::vector_free(double *vec)
{
  delete [] vec;
}


char *Util::char_malloc(int n1)
{
    char *char_ptr;

    /* pointer to the n1 array */
    char_ptr = (char *) malloc (sizeof(char)*n1);
    //char_ptr = new char[n1];

    strcpy(char_ptr, "");
return char_ptr;
}


void Util::char_free(char *vec)
{
  free(vec);
}

int Util::IsFile(string file_name)
{
 FILE *temp;

 if( (temp = fopen(file_name.c_str(),"r")) == NULL) return 0;
 else 
  {
   fclose(temp);
   return 1;
  }
}/* IsFile */


string Util::StringFind(string file_name, string st)
{
  string inputname = file_name;
  string tmpfilename;
  string str = st;

  string s;
  string xstr;
 
  tmpfilename = "input";

  int ind;
  static int flag = 0;
  if(flag == 0)
    {
      if(!IsFile(file_name))
	{
	  cerr << "The file named " << file_name << " is absent." << endl;
	  if(file_name == "") 
	    {
	      cerr << "No input file name specified." << endl;
	      exit(1);
	    }
	  else 
	    {
	      cout << "Creating " << tmpfilename << "..." << endl;
	    }
	  ofstream tmp_file(tmpfilename.c_str());
	  tmp_file << "EndOfData" << endl;
	  tmp_file.close();
	}/* if isfile */
      flag = 1;
    }/* if flag == 0 */
  
    ifstream input(inputname.c_str());
    
    input >> s;

    ind = 0;
    while(s.compare("EndOfData") != 0)
      {
	input >> xstr;
	if(s.compare(str) == 0)
	  {
	    ind++;
	    input.close();
	    return xstr;
	  }/* if right, return */
	s.clear();
	input >> s;
      }/* while */

    input.close();

    if(ind == 0)
      {
	cerr << str << " not found in " << inputname << endl; 
	cout << "Create an input file." << endl;
	exit(1);
	return xstr;
      }    
    return(0);
}/* StringFind */

double Util::DFind(string file_name, string st)
{
  string s;
  double x;
  stringstream stm;
  s = StringFind(file_name, st);
  stm << s;
  stm >> x;
  return x;
}/* DFind */


int Util::IFind(string file_name, string st)
{
 double f;
 f = DFind(file_name, st);

 return (int) (f + 0.5);
}/* IFind */


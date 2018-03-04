// Init.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Init.h"

//**************************************************************************
// Init class.

vector <complex<double> > Init::solveAxb(Parameters *param, complex<double>* A, complex<double>* b_in)
{
  int Nc = param->getNc();
  int Nc2m1 = Nc*Nc-1;

  vector <complex<double> > xvec;
  xvec.reserve(Nc2m1); 
  //  xvec = new complex<double>[Nc2m1];
  //complex<double> xvec[8];
  
  double a_data[128];

#pragma omp prallel   
  for(int i=0; i<64; i++)
    {
      a_data[2*i] = real(A[i]); 
      a_data[2*i+1] = imag(A[i]);
    }

  double b_data[2*Nc2m1];

  for(int i=0; i<Nc2m1; i++)
    {
      b_data[2*i] = real(b_in[i]); 
      b_data[2*i+1] = imag(b_in[i]);
    }

  gsl_matrix_complex_view m = gsl_matrix_complex_view_array (a_data, Nc2m1, Nc2m1);
  gsl_vector_complex_view c = gsl_vector_complex_view_array (b_data, Nc2m1);
  gsl_vector_complex *x = gsl_vector_complex_alloc (Nc2m1);

  //  gsl_vector_complex_fprintf (stdout, &c.vector, "%g");

  int s;
  gsl_permutation * p = gsl_permutation_alloc (Nc2m1);
  gsl_linalg_complex_LU_decomp (&m.matrix, p, &s);
  gsl_linalg_complex_LU_solve (&m.matrix, p, &c.vector, x);
  //  printf ("x = \n");
  //gsl_vector_complex_fprintf (stdout, x, "%g");
  gsl_permutation_free(p);
  
  if(Nc == 3)
    {
      xvec.push_back(complex<double>(GSL_REAL(gsl_vector_complex_get(x,0)),GSL_IMAG(gsl_vector_complex_get(x,0))));
      xvec.push_back(complex<double>(GSL_REAL(gsl_vector_complex_get(x,1)),GSL_IMAG(gsl_vector_complex_get(x,1))));
      xvec.push_back(complex<double>(GSL_REAL(gsl_vector_complex_get(x,2)),GSL_IMAG(gsl_vector_complex_get(x,2))));
      xvec.push_back(complex<double>(GSL_REAL(gsl_vector_complex_get(x,3)),GSL_IMAG(gsl_vector_complex_get(x,3))));
      xvec.push_back(complex<double>(GSL_REAL(gsl_vector_complex_get(x,4)),GSL_IMAG(gsl_vector_complex_get(x,4))));
      xvec.push_back(complex<double>(GSL_REAL(gsl_vector_complex_get(x,5)),GSL_IMAG(gsl_vector_complex_get(x,5))));
      xvec.push_back(complex<double>(GSL_REAL(gsl_vector_complex_get(x,6)),GSL_IMAG(gsl_vector_complex_get(x,6))));
      xvec.push_back(complex<double>(GSL_REAL(gsl_vector_complex_get(x,7)),GSL_IMAG(gsl_vector_complex_get(x,7))));
    }
  else if (Nc == 2)
    {
      xvec.push_back(complex<double>(GSL_REAL(gsl_vector_complex_get(x,0)),GSL_IMAG(gsl_vector_complex_get(x,0))));
      xvec.push_back(complex<double>(GSL_REAL(gsl_vector_complex_get(x,1)),GSL_IMAG(gsl_vector_complex_get(x,1))));
      xvec.push_back(complex<double>(GSL_REAL(gsl_vector_complex_get(x,2)),GSL_IMAG(gsl_vector_complex_get(x,2))));
    }
  gsl_vector_complex_free(x);
   
  return(xvec);
  // delete [] xvec;
}


void Init::sampleTA(Parameters *param, Random* random, Glauber* glauber)
{ 
  ReturnValue rv, rv2;
  cout << "Sampling nucleon positions ... ";

  if(param->getNucleonPositionsFromFile()==0)
    {
      int A1,A2;
      A1 = static_cast<int>(glauber->nucleusA1())*param->getAverageOverNuclei(); // projectile
      A2 = static_cast<int>(glauber->nucleusA2())*param->getAverageOverNuclei(); // target
      if((glauber->nucleusA1()==1 || glauber->nucleusA2()==1) && param->getAverageOverNuclei()>1)
	{
	  cerr << "Averaging not supported for collisions involving protons ... Exiting." << endl;
	  exit(1);
	}
      
      if(A1==1)
	{
	  rv.x=0.;
	  rv.y=0;
	  rv.collided=0;
	  nucleusA.push_back(rv);  
	}   
      else if(A1==2) // deuteron
	{

// 	  //testing
// 	  double r, scale=10.;
// 	  double rbins[100];
// 	  int count=0;
// 	  for (int b=0; b<100; b++)
// 	    {
// 	      rbins[b]=0;
// 	    }

// 	  for (int i=0; i<100000; i++)
// 	    {
// 	      rv = glauber->SampleTARejection(random,1);
// 	      r = sqrt((rv.x)*(rv.x)+(rv.y)*(rv.y));
// 	      for (int b=0; b<100; b++)
// 		{
// 		  if (r > scale/100.*(b) && r < scale/100.*(b+1)) 
// 		    {
// 		      rbins[b]+=100./scale;
// 		      count++;
// 		    }
// 		}
// 	    }
// 	  for (int b=0; b<100; b++)
// 	    {
// 	      cout << b*scale/100. << " " << rbins[b]/count <<  endl;
// 	    }

	  
	  rv = glauber->SampleTARejection(random,1);
	  param->setRnp(sqrt(rv.x*rv.x+rv.y*rv.y));
	  // we sample the neutron proton distance, so distance to the center needs to be divided by 2
	  rv.x = rv.x/2.;
	  rv.y = rv.y/2.;
	  nucleusA.push_back(rv);
	 
	  // other nucleon is 180 degrees rotated:	  
	  rv.x = -rv.x;
	  rv.y = -rv.y;
	  rv.collided=0;
	  nucleusA.push_back(rv);

	}   
      else if(A1==3) // He3
	{
	  //sample the position in the file
	  ifstream fin;
	  fin.open("he3_plaintext.dat"); 
	     
	  double dummy;
	  double ran2 = random->genrand64_real3();   // sample the position in the file uniformly (13699 events in file)
	  int nucleusNumber = static_cast<int>(ran2*13699);

	  cout << "using nucleus Number = " << nucleusNumber << endl;
	  
	  // go to the correct line in the file
	  fin.seekg(std::ios::beg);
	  for(int i=0; i < nucleusNumber; ++i)
	    {
	      fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	    }
	  // am now at the correct line in the file
	  
	  // start reading one nucleus (3 positions)
	  int A=0;

	  while(A<glauber->nucleusA1())
	    {
	      if(!fin.eof())
		{  
		  fin >> rv.x;
		  fin >> rv.y;
		  fin >> dummy; // don't care about z direction
		  rv.collided=0;
		  nucleusA.push_back(rv);
		  A++;
		  cout << "A=" << A << ", x=" << rv.x << ", y=" << rv.y << endl;
		}
	    }
	  
	  	  
	  fin.close();
	  
	  param->setA1FromFile(A);
	
	}   
      else
	{
	  for (int i = 0; i < A1; i++) // get all nucleon coordinates
	    {
	      rv = glauber->SampleTARejection(random,1);
	      nucleusA.push_back(rv);
	    }
	}
    
      if(A2==1)
	{
	  rv2.x=0.;
	  rv2.y=0;
	  rv2.collided=0;
	  nucleusB.push_back(rv2);
	}   
      else if(A2==2) // deuteron
	{
	  rv = glauber->SampleTARejection(random,2);
	  // we sample the neutron proton distance, so distance to the center needs to be divided by 2
	  param->setRnp(sqrt(rv.x*rv.x+rv.y*rv.y));

	  rv.x = rv.x/2.;
	  rv.y = rv.y/2.;
	  nucleusB.push_back(rv);
	 
	  // other nucleon is 180 degrees rotated:	  
	  rv.x = -rv.x;
	  rv.y = -rv.y;
	  rv.collided=0;
	  nucleusB.push_back(rv);

	}   
    else if(A2==3) // He3
	{
	  //sample the position in the file
	  ifstream fin;
	  fin.open("he3_plaintext.dat"); 
	     
	  double dummy;
	  double ran2 = random->genrand64_real3();   // sample the position in the file uniformly (13699 events in file)
	  int nucleusNumber = static_cast<int>(ran2*13699);

	  cout << "using nucleus Number = " << nucleusNumber << endl;
	  
	  // go to the correct line in the file
	  fin.seekg(std::ios::beg);
	  for(int i=0; i < nucleusNumber; ++i)
	    {
	      fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	    }
	  // am now at the correct line in the file
	  
	  // start reading one nucleus (3 positions)
	  int A=0;

	  while(A<glauber->nucleusA2())
	    {
	      if(!fin.eof())
		{  
		  fin >> rv.x;
		  fin >> rv.y;
		  fin >> dummy; // don't care about z direction
		  rv.collided=0;
		  nucleusB.push_back(rv);
		  A++;
		  cout << "A=" << A << ", x=" << rv.x << ", y=" << rv.y << endl;
		}
	    }
	  
	  	  
	  fin.close();
	  
	  param->setA2FromFile(A);

	  //	  cout << glauber->nucleusA1() << " " <<glauber->nucleusA1() << endl;
	
	}   
      else
	{
	  for (int i = 0; i < A2; i++) // get all nucleon coordinates
	    {
	      rv2 = glauber->SampleTARejection(random,2);
	      nucleusB.push_back(rv2);
	    }
	}
      cout << "done. " << endl;
    }
  else if (param->getNucleonPositionsFromFile()==1)
    {
      ifstream fin;
      fin.open("nucleus1.dat"); 
      cout << "Reading nucleon positions from files 'nucleus1.dat' and 'nucleus2.dat' ... " << endl;
      int A=0;
      int A2=0;
      if(fin)
	{
	  while(!fin.eof())
	    {  
	      fin >> rv.x;
	      fin >> rv.y;
	      rv.collided=0;
	      nucleusA.push_back(rv);
	      A++;
	    }
	}
    
      
      fin.close();
      fin.open("nucleus2.dat"); 
      
      if(fin)
	{
	  while(!fin.eof())
	    {  
	      fin >> rv.x;
	      fin >> rv.y;
	      rv.collided=0;
	      nucleusB.push_back(rv);
	      A2++;
	    }
	}
      
      A=A-1;
      A2=A2-1;

      cout << "A1 (from file) = " << A << endl;
      cout << "A2 (from file) = " << A2 << endl;

      param->setA1FromFile(A);
      param->setA2FromFile(A2);
      
      fin.close();
 
      cout << " ... done." << endl;
    }
  else if (param->getNucleonPositionsFromFile()==2) // Read in Alvioli's nucleon positions including correlations
    {
      if (glauber->nucleusA1()!=208 && glauber->nucleusA2()!=208)
	{
	  cerr << "[Init.cpp]: The option 'getNucleonPositionsFromFile == 2' only works for either both nuclei Pb-208 or Projectile p and Target Pb-208. Exiting." << endl;
	  exit(1);
	}

      cout << endl << "Retrieving nuclei from " << endl;
      
      //generate the file name
      double ran = random->genrand64_real3();      //sample the file name uniformly
      int fileNumber = static_cast<int>(ran*10+1);
	  
	  stringstream str_file;
	  str_file.str("");
	  if(fileNumber < 10)
	    str_file << "/global/homes/s/schenke/Alvioli-Pb208/pb208-0";
	  else
	    str_file << "/global/homes/s/schenke/Alvioli-Pb208/pb208-";
	  str_file << fileNumber;
	  str_file << ".dat";
	  string fileName = str_file.str();  
	  
	  //open the file
	  ifstream fin;
	  fin.open(fileName.c_str());
	  if (!fin)
	    {
	      cerr << "File " << fileName << " not found. Trying alternative location:" << endl;
	      str_file.str("");
	      if(fileNumber < 10)
		str_file << "./Alvioli-Pb208/pb208-0";
	      else
		str_file << "./Alvioli-Pb208/pb208-";
	      str_file << fileNumber;
	      str_file << ".dat";
	      fileName = str_file.str();  
	      fin.open(fileName.c_str());
	    }
	  
	  if (!fin)
	    {
	      cerr << "File " << fileName << " not found. Exiting." << endl;
	      exit(1);
	    }
	  
	  
	  cout << "Reading nucleon positions for nuceus A from file " << fileName << " ... " << endl;
	  
	  
	  //sample the position in the file
	  double ran2 = random->genrand64_real3();   // sample the position in the file uniformly (10,000 events per file)
	  int nucleusNumber = static_cast<int>(ran2*10000);
	  cout << "Nucleus Number = " << nucleusNumber << endl;
	  
	  int A=0;
	  int A2=0;
	  double dummy;
	  
	  // go to the correct line in the file
	  fin.seekg(std::ios::beg);
	  for(int i=0; i < (nucleusNumber)*glauber->nucleusA1(); ++i)
	    {
	      fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	    }
	  // am now at the correct line in the file
	  
	  // start reading one nucleus (208 positions)

	  if(glauber->nucleusA1()==1)
	    {
	      rv.x = 0;
	      rv.y = 0;
	      rv.collided=0;
	      nucleusA.push_back(rv);
	      A = 1;
	    }
	  else
	    {
	      while(A<glauber->nucleusA1())
		if(!fin.eof())
		  {  
		    fin >> rv.x;
		    fin >> rv.y;
		    fin >> dummy; // don't care about z direction
		    fin >> dummy; // don't care about isospin
		    rv.collided=0;
		    nucleusA.push_back(rv);
		    A++;
		    //		    cout << "A=" << A << "/" << glauber->nucleusA1()<<endl;
		    //cout << rv.x << " " << rv.y << endl;
		  }
	      
	    }  
	  fin.close();
	    
	  // do the second nucleus (Target)
	  
	  
	  //generate the file name
	  ran = random->genrand64_real3();      //sample the file name uniformly
	  fileNumber = static_cast<int>(ran*10+1);
	  
	  str_file.str("");
	  if(fileNumber < 10)
	    str_file << "/global/homes/s/schenke/Alvioli-Pb208/pb208-0";
	  else
	    str_file << "/global/homes/s/schenke/Alvioli-Pb208/pb208-";
	  str_file << fileNumber;
	  str_file << ".dat";
	  fileName = str_file.str();  
	  
	  //open the file
	  fin.open(fileName.c_str());
	  if (!fin)
	    {
	      cerr << "File " << fileName << " not found. Trying alternative location:" << endl;
	      str_file.str("");
	      if(fileNumber < 10)
		str_file << "./Alvioli-Pb208/pb208-0";
	      else
		str_file << "./Alvioli-Pb208/pb208-";
	      str_file << fileNumber;
	      str_file << ".dat";
	      fileName = str_file.str();  
	      fin.open(fileName.c_str());
	    }
	  
	  if (!fin)
	    {
	      cerr << "File " << fileName << " not found. Exiting." << endl;
	      exit(1);
	    }
	  
	  
	  cout << "Reading nucleon positions for nucleus B from file " << fileName << " ... " << endl;
	  
	  
	  //sample the position in the file
	  ran2 = random->genrand64_real3();   // sample the position in the file uniformly (10,000 events per file)
	  nucleusNumber = static_cast<int>(ran2*10000);
	  cout << "Nucleus Number = " << nucleusNumber << endl;
	  
	  // go to the correct line in the file
	  fin.seekg(std::ios::beg);
	  for(int i=0; i < (nucleusNumber)*glauber->nucleusA1(); ++i)
	    {
	      fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	    }
	  // am now at the correct line in the file
	  
	  // start reading one nucleus (208 positions)
	  if(glauber->nucleusA2()==1)
	    {
	      rv.x = 0;
	      rv.y = 0;
	      rv.collided=0;
	      nucleusB.push_back(rv);
	      A2 = 1;
	    }
	  else
	    {
	      while(A2<glauber->nucleusA2())
		if(!fin.eof())
		  {  
		    fin >> rv.x;
		    fin >> rv.y;
		    fin >> dummy; // don't care about z direction
		    fin >> dummy; // don't care about isospin
		    rv.collided=0;
		    nucleusB.push_back(rv);
		    A2++;
		    //		cout << "A2=" << A2 << "/" << glauber->nucleusA2()<<endl;
		    //cout << rv.x << " " << rv.y << endl;
		  }
	    }
	  
	  fin.close();
	  
	  param->setA1FromFile(A);
	  param->setA2FromFile(A2);
	  
	  cout << " ... done." << endl;
	}
      else 
	{
	  cerr << "NucleonPositionsFromFile can be 0 (sample nucleons) or 1 or 2 (read from files) - you chose " << 
	    param->getNucleonPositionsFromFile() << ". Exiting." << endl;
	  exit(1);
	}
    }  
  
  
  void Init::readNuclearQs(Parameters *param)
  {
  cout << "Reading Q_s(sum(T_p),y) from file ";
  // steps in qs0 and Y in the file
  // double y[iymaxNuc];
  // double qs0[ibmax];
  string dummy;
  string T, Qs;
  // open file
  ifstream fin;
  fin.open((param->getNucleusQsTableFileName()).c_str()); 

  cout << param->getNucleusQsTableFileName() << " ... " ;

  if(fin)
    {
      for (int iT=0; iT<iTpmax; iT++)
	{
	  for (int iy=0; iy<iymaxNuc; iy++)
	    {
	      if (!fin.eof())
		{  
		  fin >> dummy;
		  fin >> T;
		  Tlist[iT]=atof(T.c_str());
		  fin >> Qs;
		  Qs2Nuclear[iT][iy]=atof(Qs.c_str());
		  //cout << iT << " " << iy <<  " " << T << " " << iy*deltaYNuc << " " << Qs2Nuclear[iT][iy] << endl;
		}
	      else 
		{
		  cerr << " End of file reached prematurely. Did the file change? Exiting." << endl;
		  exit(1);
		}
	    }
	}
      fin.close();
      cout << " done." << endl;
    }
  else
    {
      cout << "[Init.cpp:readNuclearQs]: File qs2_Adj_Y_qs20_IPSat.dat does not exist. Exiting." << endl;
      exit(1);
    }
}  

// Q_s as a function of \sum T_p and y (new in this version of the code - v1.2 and up)
double Init::getNuclearQs2(Parameters *param, Random* random, double T, double y)
{
  double value, fracy, fracT, QsYdown, QsYup;
  int posb, posy, check=0;
  posy = static_cast<int>(floor(y/deltaYNuc+0.0000001));

  if (y>iymaxNuc*deltaYNuc)
    {
      cout << " [Init:getNuclearQs2]:ERROR: y out of range. Maximum y value is " << iymaxNuc*deltaYNuc << ", you used " << y << ". Exiting." << endl;
      exit(1);
    }

  if ( T > Qs2Nuclear[iTpmax-1][iymaxNuc-1] )
    {
      cout << "T=" << T << ", maximal T in table=" << Tlist[iTpmax-1] << endl;
      cout << " [Init:getNuclearQs2]:ERROR: out of range. Exiting." << endl;
      exit(1);
    }

  if ( T < Tlist[0] )
    {
      check = 1;
      return 0.;
    }
  
  for(int iT=0; iT<iTpmax; iT++)
    {
      if (T>=Tlist[iT] && T<Tlist[iT+1])
	{
	  fracT = (T-Tlist[iT])/(Tlist[iT+1]-Tlist[iT]);
	  fracy = (y-static_cast<double>(posy)*deltaYNuc)/deltaYNuc;
	 
	  QsYdown = (fracT)*(Qs2Nuclear[iT+1][posy])+(1.-fracT)*(Qs2Nuclear[iT][posy]);
	  QsYup = (fracT)*(Qs2Nuclear[iT+1][posy+1])+(1.-fracT)*(Qs2Nuclear[iT][posy+1]);

	  //	  cout << posy << endl;
	  //cout << Qs2Nuclear[iT+1][posy] << " " <<  QsYdown << " " << QsYup << endl;

	  value = (fracy*QsYup+(1.-fracy)*QsYdown);//*hbarc*hbarc;
		  
//   	  cout << "T=" << T << ", lowT=" << Tlist[iT] << ", highT=" << Tlist[iT+1] << endl;
//   	  cout << "y=" << y << ", lowy=" << (posy)*deltaYNuc << ", highy=" << (posy+1)*deltaYNuc << endl;
// 	  cout << "fracy=" << fracy << endl;
// 	  cout << "Qs^2=" << value << endl;
//   	  cout << "Qs2Nuclear[iT][posy]=" << Qs2Nuclear[iT][posy] << endl;
// 	  cout << "Qs2Nuclear[iT][posy+1]=" << Qs2Nuclear[iT][posy+1] << endl;
// 	  cout << "Qs2Nuclear[iT+1][posy]=" << Qs2Nuclear[iT+1][posy] << endl;
// 	  cout << "Qs2Nuclear[iT+1][posy+1]=" << Qs2Nuclear[iT+1][posy+1] << endl;

	  check++;
	  continue;
	}
    }

   if (check!=1)
    {
      cout << check << ": T=" << T << endl ;
      // cerr << " [Init:getNuclearQs2]:ERROR: something went wrong in determining the value of Qs^2. Exiting." << endl;
      // exit(1);
      cerr << " [Init:getNuclearQs2]:ERROR: something went wrong in determining the value of Qs^2. Using maximal T_p" << endl;
      value = Tlist[iTpmax-1];
    }
 
  return value;
}  


// set g^2\mu^2 as the sum of the individual nucleons' g^2\mu^2, using Q_s(b,y) prop tp g^mu(b,y)
// also compute N_part using Glauber
void Init::setColorChargeDensity(Lattice *lat, Parameters *param, Random *random, Glauber *glauber)
{
  int pos,posA,posB;
  int N = param->getSize();
  int A1, A2;
  int check=0;
  double xVal;
  if(param->getNucleonPositionsFromFile()==0)
    {
      A1 = static_cast<int>(glauber->nucleusA1())*param->getAverageOverNuclei();
      A2 = static_cast<int>(glauber->nucleusA2())*param->getAverageOverNuclei();
    }
  else
    {
      A1 = param->getA1FromFile();
      A2 = param->getA2FromFile();
    }

  int Npart = 0;
  int Ncoll = 0;
  double g2mu2A, g2mu2B;
  double b = param->getb();
  double x, xm;
  double y, ym;
  double r;
  double L = param->getL();
  double rapidity;
  double P,m;
  if(param->getUsePseudoRapidity()==0)
    rapidity = param->getRapidity();
  else
    {
      // when using pseudorapidity as input convert to rapidity here. later include Jacobian in multiplicity and energy
      cout << "Using pseudorapidity " << param->getRapidity() << endl;
      m=param->getJacobianm(); // in GeV
      P=0.13+0.32*pow(param->getRoots()/1000.,0.115); //in GeV
      rapidity = 0.5 * log(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(P*P))+sinh(param->getRapidity())
			   / (sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(P*P))-sinh(param->getRapidity())));
      cout << "Corresponds to rapidity " << rapidity << endl;
    }
  
  double yIn = rapidity;//param->getRapidity();
  double a = L/N; // lattice spacing in fm
  double dx, dy, dij;
  double d2 = param->getSigmaNN()/(PI*10.);          // in fm^2
  double averageQs = 0.;
  double averageQs2 = 0.;
  double averageQs2Avg = 0.;
  double averageQs2min = 0.;
  double averageQs2min2 = 0.;
  int count = 0;
  int count2 = 0;
  double nucleiInAverage;
  nucleiInAverage = static_cast<double>(param->getAverageOverNuclei());



  //testing distributions
//   double binsNBD[30];
//   double binsGauss[30];
//   double nbd;
//   double scale=2.;
//   int bins=30;
//   int events=2000;


  double gaussA[A1][3];
  double gaussB[A2][3];

  for (int i = 0; i<A1; i++) 
    {
      for (int iq = 0; iq<3; iq++) 
	{
	  gaussA[i][iq]=1.;
	}
    }    
  for (int i = 0; i<A2; i++) 
    {
      for (int iq = 0; iq<3; iq++) 
	{
	  gaussB[i][iq]=1.;
	}    
    }
  
  // let the log fluctuate
  if(param->getSmearQs() == 1)
    {
      if (A1>0) 
	{
	  for (int i = 0; i<A1; i++) 
	    {
	      for (int iq = 0; iq<3; iq++) 
		{
		  gaussA[i][iq] = (exp(random->Gauss(0,param->getSmearingWidth())))/1.13; // dividing by 1.13 restores the same mean Q_s 
		  //cout << i << " " << iq << " " << gaussA[i][iq] << endl; 
		  //	  if (gaussA[i]<0)
		  //  gaussA[i]=0.;
		}
	    }
	}
      if (A2>0) 
	{
	  for (int i = 0; i<A2; i++) 
	    {
	      for (int iq = 0; iq<3; iq++) 
		{
		  gaussB[i][iq] = (exp(random->Gauss(0,param->getSmearingWidth())))/1.13;
		  
		  //	      cout << i << " " << iq << " " << gaussB[i][iq] << endl; 
		  //if (gaussB[i]<0)
		  //  gaussB[i]=0.;
		}
	    }
	}
    }

//   if(param->getSmearQs() == 1)
//     {
//       // set to different ratios for nucleus A and B
//       gauss = random->Gauss(1,param->getSmearingWidth());
//       cout << "QsmuRatio before=" << param->getQsmuRatio() << endl;
//       if (gauss<0)
// 	gauss=0.;
//       param->setQsmuRatio(1./(1./param->getQsmuRatio())*gauss);

//       gauss = random->Gauss(1,param->getSmearingWidth());
//       if (gauss<0)
// 	gauss=0.;
//       param->setQsmuRatioB(1./(1./param->getQsmuRatio())*gauss);
      
// //       ofstream fout1("NBD.dat",ios::out); 
// //       ofstream fout2("Gauss.dat",ios::out); 
      
// //       for(int ib=0; ib<30; ib++)
// // 	{
// // 	  binsNBD[ib]=0.;
// // 	  binsGauss[ib]=0.;
// // 	}

// //       for (int l=0; l<events; l++)
// // 	{
// // 	  nbd = random->NBD(100,40)/100.;
// // 	  gauss = random->Gauss(1,0.18);
	  
// // 	  for(int ib=0; ib<bins; ib++)
// // 	    {
// // 	      if (nbd >= ib*(scale/static_cast<double>(bins)) && nbd < (ib+1)*(scale/static_cast<double>(bins)))
// // 		binsNBD[ib]+=1/static_cast<double>(events);
// // 	      if (gauss >= ib*(scale/static_cast<double>(bins)) && gauss < (ib+1)*(scale/static_cast<double>(bins)))
// // 		binsGauss[ib]+=1/static_cast<double>(events);
// // 	    }
// // 	}

// //       for(int ib=0; ib<bins; ib++)
// // 	{
// // 	  fout1 << (ib+0.5)*(scale/static_cast<double>(bins)) << " " << binsNBD[ib] << endl;
// // 	  fout2 << (ib+0.5)*(scale/static_cast<double>(bins)) << " " << binsGauss[ib] << endl;
// // 	}


//       cout << "QsmuRatio for nucleus A smeared with width " << param->getSmearingWidth() << " is " << param->getQsmuRatio() << endl;
//       cout << "QsmuRatio for nucleus B smeared with width " << param->getSmearingWidth() << " is " << param->getQsmuRatioB() << endl;
//     }
//   else
  
  param->setQsmuRatioB(param->getQsmuRatio());

  if(param->getUseNucleus() == 0)
    {
      for(int ix=0; ix<N; ix++) // loop over all positions
	{
	  for(int iy=0; iy<N; iy++)
	    {
	      pos = ix*N+iy;
	      lat->cells[pos]->setg2mu2A(param->getg2mu()*param->getg2mu()/param->getg()/param->getg());
	      lat->cells[pos]->setg2mu2B(param->getg2mu()*param->getg2mu()/param->getg()/param->getg());
	    }
	}
      param->setSuccess(1);
      cout << "constant color charge density set" << endl;
      return;
    }
  

  for(int ix=0; ix<N; ix++) // loop over all positions
    {
      for(int iy=0; iy<N; iy++)
	{
	  pos = ix*N+iy;
	  lat->cells[pos]->setg2mu2A(0.);
	  lat->cells[pos]->setg2mu2B(0.);
	}
    }
  
  // compute N_part
  // positions are shifted here. not later as in previous versions. bshift below (in init(..)) is zero.

  if(A1 == 1 && A2 > 1) 
    {
      for (int i = 0; i<A2; i++) 
	{
	  nucleusB.at(i).x=nucleusB.at(i).x+b;
	}   
    }
  else if(A2 == 1 && A1 > 1) 
    {
      for (int i = 0; i<A1; i++) 
	{
	  nucleusA.at(i).x=nucleusA.at(i).x-b;
	}   
    }
  else
    {
      for (int i = 0; i<A1; i++) // shift the nuclei's position by -b/2 or +b/2 respectively
	{
	  nucleusA.at(i).x=nucleusA.at(i).x-b/2.;
	}   
      
      for (int i = 0; i<A2; i++) // shift the nuclei's position by -b/2 or +b/2 respectively
	{
	  nucleusB.at(i).x=nucleusB.at(i).x+b/2.;
	}   
    }

  
  double maxT=0;

  double bp2,T,BG;
  BG = param->getBG();
  double BGq = param->getBGq(); // quark size in GeV^-2
  double xi = param->getProtonAnisotropy();
  double phi;

  if(xi!=0.)
    {
      for (int i = 0; i<A1; i++) 
	{
	  nucleusA.at(i).phi = 2*M_PI*random->genrand64_real2();
	}      
      
      for (int i = 0; i<A2; i++) 
	{
	  nucleusB.at(i).phi = 2*M_PI*random->genrand64_real2();
	}      
    }

  //  cout << "BG=" << BG << endl;


  double xq[A1][param->getUseConstituentQuarkProton()], xq2[A2][param->getUseConstituentQuarkProton()];
  double yq[A1][param->getUseConstituentQuarkProton()], yq2[A2][param->getUseConstituentQuarkProton()];
  double avgxq=0.;
  double avgyq=0.;

  if(param->getUseConstituentQuarkProton()>0)
    {
      for (int i=0; i<A1; i++)
	{
	  avgxq=0.;
	  avgyq=0.;
	  for (int iq=0; iq<param->getUseConstituentQuarkProton(); iq++)
	    {
	      xq[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
	      yq[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
	    }
	  for (int iq=0; iq<3; iq++)
	    {
	      avgxq += xq[i][iq];
	      avgyq += yq[i][iq];
	    }
	  for (int iq=0; iq<3; iq++)
	    {
	      xq[i][iq] -= avgxq/3.;
	      yq[i][iq] -= avgyq/3.;
	      //	      cout << xq[i][iq] << " " << yq[i][iq] << endl;
	    }
	  
	  
	  // avgyq=0.;
	  // for (int iq=0; iq<3; iq++)
	  //   {
	  //     avgxq += xq[i][iq];
	  //     avgyq += yq[i][iq];
	  //   }
	  // cout << avgyq << endl;
	}
    }




      // stringstream strq_name;
      // strq_name << "qPos-" << param->getMPIRank() + param->getSeed()*16 << ".txt";
      // string q_name;
      // q_name = strq_name.str();
      // ofstream fout(q_name.c_str(),ios::out); 
      
      // fout << " x[fm]     y[fm] " << endl;
      
      // for (int i=0; i<A1; i++)
      // 	{
      // 	  for (int iq=0; iq<3; iq++)
      // 	    {
      // 	      fout << xq[i][iq] << " " << yq[i][iq] << endl;
      // 	    }
      // 	}
  //    }
  
  
  
  if(param->getUseConstituentQuarkProton()>0)
    {
      for (int i=0; i<A2; i++)
	{
	  avgxq=0.;
	  avgyq=0.;
	  for (int iq=0; iq<param->getUseConstituentQuarkProton(); iq++)
	    {
	      xq2[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
	      yq2[i][iq] = sqrt(BG*hbarc*hbarc)*random->Gauss();
	    }
	  
	  for (int iq=0; iq<3; iq++)
	    {
	      avgxq += xq2[i][iq];
	      avgyq += yq2[i][iq];
	    }
	  for (int iq=0; iq<3; iq++)
	    {
	      xq2[i][iq] -= avgxq/3.;
	      yq2[i][iq] -= avgyq/3.;
	    }
	  // avgyq=0.;
	  // for (int iq=0; iq<3; iq++)
	  //   {
	  //     avgxq += xq2[i][iq];
	  //     avgyq += yq2[i][iq];
	  //   }
	  // cout << avgyq << endl;
	}

      // stringstream strq_name;
      // strq_name << "qPos-" << param->getMPISize() + param->getMPIRank() + param->getSeed()*16 << ".txt";
      // string q_name;

      // q_name = strq_name.str();
      // ofstream fout(q_name.c_str(),ios::out); 
      
      // fout << " x[fm]     y[fm] " << endl;
      
      // for (int i=0; i<A2; i++)
      // 	{
      // 	  for (int iq=0; iq<3; iq++)
      // 	    {
      // 	      fout << xq2[i][iq] << " " << yq2[i][iq] << endl;
      // 	    }
      // 	}
    }


  //add all T_p's (new in version 1.2)
  for(int ix=0; ix<N; ix++) // loop over all positions
    {
      x = -L/2.+a*ix;
      for(int iy=0; iy<N; iy++)
	{
	  y = -L/2.+a*iy;
	   
	  pos = ix*N+iy;
	
	  // nucleus A 
	  lat->cells[pos]->setTpA(0.);
	  for (int i = 0; i<A1; i++) 
	    {
	      xm = nucleusA.at(i).x;
	      ym = nucleusA.at(i).y;
	      
	      if(param->getUseConstituentQuarkProton()>0)
		{
		  T = 0.;
		  for (int iq=0; iq<param->getUseConstituentQuarkProton(); iq++)
		    {
		      bp2 = (xm+xq[i][iq]-x)*(xm+xq[i][iq]-x)+(ym+yq[i][iq]-y)*(ym+yq[i][iq]-y);
		      bp2 /= hbarc*hbarc;

		      T += exp(-bp2/(2.*BGq))/(2.*PI*BGq)/(double(param->getUseConstituentQuarkProton()))*gaussA[i][iq]; // I removed the 2/3 here to make it a bit bigger
		      //	      cout << "A " << i << " " << iq << " " << gaussA[i][iq] << endl;
		    }
		}
	      else
		{
		  phi = nucleusA.at(i).phi;

		  bp2 = (xm-x)*(xm-x)+(ym-y)*(ym-y) + xi*pow((xm-x)*cos(phi) + (ym-y)*sin(phi),2.);
		  bp2 /= hbarc*hbarc;     	  
		  
		  T = sqrt(1+xi)*exp(-bp2/(2.*BG))/(2.*PI*BG)*gaussA[i][0]; // T_p in this cell for the current nucleon
		}

	      lat->cells[pos]->setTpA(lat->cells[pos]->getTpA()+T/nucleiInAverage); // add up all T_p
	      
	      maxT=max(lat->cells[pos]->getTpA()+T,maxT);

	    }
	  
	  // nucleus B 
	  lat->cells[pos]->setTpB(0.);
	  for (int i = 0; i<A2; i++) 
	    {
	      xm = nucleusB.at(i).x;
	      ym = nucleusB.at(i).y;

	      if(param->getUseConstituentQuarkProton()>0)
		{
		  T = 0.;
		  for (int iq=0; iq<param->getUseConstituentQuarkProton(); iq++)
		    {
		      bp2 = (xm+xq2[i][iq]-x)*(xm+xq2[i][iq]-x)+(ym+yq2[i][iq]-y)*(ym+yq2[i][iq]-y);
		      bp2 /= hbarc*hbarc;

		      T += exp(-bp2/(2.*BGq))/(2.*PI*BGq)/double(param->getUseConstituentQuarkProton())*gaussB[i][iq];
		      //	      cout << "B " << i << " " << iq << " " << gaussA[i][iq] << endl;

		    }
		}
	      else
		{
		  phi = nucleusB.at(i).phi;
	      
		  bp2 = (xm-x)*(xm-x)+(ym-y)*(ym-y) + xi*pow((xm-x)*cos(phi) + (ym-y)*sin(phi),2.);
		  bp2 /= hbarc*hbarc;
		  
		  T = sqrt(1+xi)*exp(-bp2/(2.*BG))/(2.*PI*BG)*gaussB[i][0]; // T_p in this cell for the current nucleon
		}
	      
	      lat->cells[pos]->setTpB(lat->cells[pos]->getTpB()+T/nucleiInAverage); // add up all T_p
	    
	      maxT=max(lat->cells[pos]->getTpB()+T,maxT);
	      
	    }
	}
    }

  //  cout << "maximal used T=" << maxT << endl;

  stringstream strNcoll_name;
  strNcoll_name << "NcollList" << param->getMPIRank() << ".dat";
  string Ncoll_name;  Ncoll_name = strNcoll_name.str();

 ofstream foutNcoll(Ncoll_name.c_str(),ios::out); 
 

  if (param->getGaussianWounding() == 0)
    {
      for (int i = 0; i<A1; i++) 
	{
	  for (int j = 0 ; j<A2 ;j++) 
	    {
	      dx = nucleusB.at(j).x-nucleusA.at(i).x;
	      dy = nucleusB.at(j).y-nucleusA.at(i).y;
	      dij = dx*dx+dy*dy;
	      if (dij < d2) 
		{
		  foutNcoll << (nucleusB.at(j).x+nucleusA.at(i).x)/2. << " " << (nucleusB.at(j).y+nucleusA.at(i).y)/2. << endl;
		  Ncoll++;
		  nucleusB.at(j).collided=1;
		  nucleusA.at(i).collided=1;
		}
	    }
	}
    }
  else
    {
      double p;
      double G=0.92;
      double ran;
      
      for (int i = 0; i<A1; i++) 
	{
	  for (int j = 0 ; j<A2 ;j++) 
	    {
	      dx = nucleusB.at(j).x-nucleusA.at(i).x;
	      dy = nucleusB.at(j).y-nucleusA.at(i).y;
	      dij = dx*dx+dy*dy;
	      
	      p = G * exp(-G*dij/d2); // Gaussian profile 
	      
	      ran = random->genrand64_real1();
	      
	      if (ran < p) 
		{
		  foutNcoll << (nucleusB.at(j).x+nucleusA.at(i).x)/2. << " " << (nucleusB.at(j).y+nucleusA.at(i).y)/2. << endl;
		  Ncoll++;
		  nucleusB.at(j).collided=1;
		  nucleusA.at(i).collided=1;
		}
	    }
	}
    }

  foutNcoll.close();
  

//   for (int i = 0; i<A1; i++) 
//     {
//       for (int j = 0 ; j<A2 ;j++) 
// 	{
// 	  dx = nucleusB.at(j).x-nucleusA.at(i).x;
// 	  dy = nucleusB.at(j).y-nucleusA.at(i).y;
// 	  dij = dx*dx+dy*dy;
// 	  if (dij < d2) 
// 	    {
// 	      nucleusB.at(j).collided=1;
// 	      nucleusA.at(i).collided=1;
// 	    }
// 	}
//     }

  // in p+p assume that they collided in any case
  if ( A1 == 1 && A2 == 1 )
    {
      nucleusB.at(0).collided=1;
      nucleusA.at(0).collided=1;
    }
  
  // stringstream strgmuA_name;
  // strgmuA_name << "gmuA" << param->getMPIRank() << ".dat";
  // string gmuA_name;
  // gmuA_name = strgmuA_name.str();
	  
  // ofstream fout(gmuA_name.c_str(),ios::out); 
  double outvalue;
  double alphas;
  double Ydeviation = 10000;
  double QsA, QsB, distanceA, distanceB;

  QsA = 1;
  QsB = 1;

  Npart = 0;

  for (int i = 0; i<A1; i++) 
    {
      if (nucleusA.at(i).collided==1)
	Npart++;
    }
  
  for (int i = 0; i<A2; i++) 
    {
      if (nucleusB.at(i).collided==1)
	Npart++;
    }
  
//   if ( Npart == 0 && param->getUseFixedNpart()==0)
//     {
//       cout << "no collision happened. Exiting." << endl;
//       exit(1);
//     }
  
  param->setNpart(Npart);

  if(param->getUseFixedNpart()!=0 && Npart!=param->getUseFixedNpart())
    {
      cout << "current Npart = " << Npart << endl;
      return;
    }

  // get Q_s^2 (and from that g^2mu^2) for a given \sum T_p and Y
  for(int ix=0; ix<N; ix++) // loop over all positions
    {
      x = -L/2.+a*ix;
      for(int iy=0; iy<N; iy++)
	{
	  check = 0;
	  y = -L/2.+a*iy;
	  Ydeviation = 10000;
	  pos = ix*N+iy;
	  
	  // this version removes noise outside the interaction region 
	  // by checking whether we are inside a wounded nucleon
	  // using 2 times the nucleon size as a cutoff
// 	  for (int i = 0; i<A1; i++) 
// 	    {
// 	      xm = nucleusA.at(i).x;
// 	      ym = nucleusA.at(i).y;
// 	      r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
	      
// 	      if(r<sqrt(2.*0.1*param->getSigmaNN()/PI) && nucleusA.at(i).collided==1)
// 		{
// 		  check=1;
// 		}
// 	    }
	  
// 	  for (int i = 0; i<A2; i++) 
// 	    {
// 	      xm = nucleusB.at(i).x;
// 	      ym = nucleusB.at(i).y;
// 	      r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
	      
// 	      if(r<sqrt(2.*0.1*param->getSigmaNN()/PI) && nucleusB.at(i).collided==1 && check==1)
// 		check=2;
// 	    }
	  

// 	  // cut at a radius of ~1.8 fm
// 	  for (int i = 0; i<A1; i++) 
// 	    {
// 	      if(lat->cells[pos]->getTpA() > 0.0019 && nucleusA.at(i).collided==1)
// 		{
// 		  check=1;
// 		}
// 	    }
	  
// 	  for (int i = 0; i<A2; i++) 
// 	    {
// 	      if(lat->cells[pos]->getTpB() > 0.0019 && nucleusB.at(i).collided==1 && check==1)
// 		check=2;
// 	    }


// 	  // cut at a radius of ~1.8 fm
// 	  if(lat->cells[pos]->getTpA() > 0.0019 )
// 	    {
// 	      check=1;
// 	    }
	  
// 	  if(lat->cells[pos]->getTpB() > 0.0019 && check==1)
// 	    check=2;
	  


//	  cut proton at a radius of rmax [fm] (about twice the gluonic radius to be generous)	  
	  if(log(2*M_PI*BG*lat->cells[pos]->getTpA())<0.)
	    distanceA = sqrt(-2.*BG*log(2*M_PI*BG*lat->cells[pos]->getTpA()))*hbarc;
	  else distanceA=0.;

	  if(log(2*M_PI*BG*lat->cells[pos]->getTpB())<0.)
	    distanceB = sqrt(-2.*BG*log(2*M_PI*BG*lat->cells[pos]->getTpB()))*hbarc;
	  else distanceB=0.;

	  // if(distanceA>0.1)
	  //   cout << lat->cells[pos]->getTpA()<< " " << distanceA << " " << param->getRmax() << endl;
	  
	  if(distanceA < param->getRmax())
	    {
	      check=1;
	    }
	 
	  if(distanceB < param->getRmax() && check==1)
	    {		
	      check=2;
	    }


	  double exponent=5.6; // see 1212.2974 Eq. (17)
	  if(check==2)
	    {
	      if ( param->getUseFluctuatingx() == 1)
		{
		  // iterative loops here to determine the fluctuating Y
		  // _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
		  while (abs(Ydeviation) > 0.001) 
		    {
		      if(rapidity>=0)
			QsA = sqrt(getNuclearQs2(param, random, lat->cells[pos]->getTpA(), abs(rapidity)));
		      else 
			{
			  xVal = QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn);
			  //cout << " QsA=" << QsA << ", param->getRoots()=" << param->getRoots() << ", exp(yIn)=" << exp(yIn) << endl; 
			  if(xVal==0)
			    QsA=0.;
			  else
			    QsA = sqrt(getNuclearQs2(param, random, lat->cells[pos]->getTpA(), 0.))*
			      sqrt(pow((1-xVal)/(1-0.01),exponent)*pow((0.01/xVal),0.2));
			  //cout << "xVal=" << xVal << endl;
			  //cout << "QsA=" << QsA << endl;
			}
		      if(QsA == 0)
			{
			  Ydeviation = 0;
			  lat->cells[pos]->setg2mu2A(0.);
			}
		      else
			{
			  // nucleus A 
			  lat->cells[pos]->setg2mu2A(QsA*QsA/param->getQsmuRatio()/param->getQsmuRatio()
						     *a*a/hbarc/hbarc/param->getg()); // lattice units? check
			  
			  
			  Ydeviation = rapidity - log(0.01/(QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn)));
			  rapidity = log(0.01/(QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn)));
			}	    
		    }
		  if(lat->cells[pos]->getg2mu2A()!=lat->cells[pos]->getg2mu2A())
		    {
		      lat->cells[pos]->setg2mu2A(0.);
		    }
		  
		  // if(ix==N/2 && iy==N/2)
		  //   {
		  //     if(QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn)<=0.01)
		  // 	cout  << "rapidity_A=" << rapidity << endl;
		  //     else
		  // 	cout  << "rapidity_A=" << log(0.01/xVal) << endl;
		  // 	cout  << "xVal=" << xVal << endl;
		  //     cout  << "Q_sA=" << QsA << endl;
		  //   }
		   
		  

		  Ydeviation = 10000;
		  
		  while (abs(Ydeviation) > 0.001) 
		    {
		      if(rapidity>=0)
			QsB = sqrt(getNuclearQs2(param, random, lat->cells[pos]->getTpB(), abs(rapidity)));
		      else
			{
			  xVal = QsB*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(-yIn);
			  if(xVal==0)
			    QsB=0.;
			  else
			    QsB = sqrt(getNuclearQs2(param, random, lat->cells[pos]->getTpB(), 0.))*
			    sqrt(pow((1-xVal)/(1-0.01),exponent)*pow((0.01/xVal),0.2));
			}
		      if(QsB == 0)
			{
			  Ydeviation = 0;
			  lat->cells[pos]->setg2mu2B(0.);
			}
		      else
			{
			  // nucleus B 
			  lat->cells[pos]->setg2mu2B(QsB*QsB/param->getQsmuRatioB()/param->getQsmuRatioB()
						     *a*a/hbarc/hbarc/param->getg()/param->getg()); 
			  
			  
			  //      cout << ix << " " << iy << endl;
			  // 		      cout << " QsB = " << QsB << " GeV" << endl;
			  // 		      cout << " x= " << (QsB/2./param->getRoots()) << endl;
			  
			  Ydeviation = rapidity - log(0.01/(QsB*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(-yIn)));
			  rapidity = log(0.01/(QsB*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(-yIn)));
			}
		    }   
		  if(lat->cells[pos]->getg2mu2B()!=lat->cells[pos]->getg2mu2B())
		    {
		      lat->cells[pos]->setg2mu2B(0.);
		    }
		// if(ix==N/2 && iy==N/2)
		//   {
		//     if(QsB*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(-yIn)<=0.01)
		//       cout  << "rapidity_B=" << rapidity << endl;
		//     else
		//       cout  << "rapidity_B=" << log(0.01/xVal) << endl;
		//     cout  << "Q_sB=" << QsB << endl;
		//   }
	   
		  
		  // _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
		  // end iterative loops here
		}
	      else
		{
		  // nucleus A 
		  lat->cells[pos]->setg2mu2A(getNuclearQs2(param, random, lat->cells[pos]->getTpA(), rapidity)/param->getQsmuRatio()/param->getQsmuRatio()
					     *a*a/hbarc/hbarc/param->getg()/param->getg()); // lattice units? check
		  
		  // nucleus B 
		  lat->cells[pos]->setg2mu2B(getNuclearQs2(param, random, lat->cells[pos]->getTpB(), rapidity)/param->getQsmuRatioB()/param->getQsmuRatioB()
					     *a*a/hbarc/hbarc/param->getg()/param->getg()); 
		 
		}
	    }
	}
    }
  
  
  // output gmu 
  count=0;
  count2=0;
  double Tpp=0.;
      
  for(int ix=0; ix<N; ix++) // loop over all positions
    {
      for(int iy=0; iy<N; iy++)
	{
	  check = 0;
	  pos = ix*N+iy;
	  x = -L/2.+a*ix;
	  y = -L/2.+a*iy;
	  //  outvalue = sqrt(lat->cells[pos]->getg2mu2B())/a*hbarc; // in GeV
	  outvalue = lat->cells[pos]->getg2mu2A();
	  
	  // posA = static_cast<int>(floor((x-b/2.+L/2.)/a+0.00000001))*N+iy;
	  // posB = static_cast<int>(floor((x+b/2.+L/2.)/a+0.00000001))*N+iy;
	  
	  posA = pos;
	  posB = pos;
	  
	  if(posA>0 && posA<(N-1)*N+N-1)
	    {
	      g2mu2A = lat->cells[posA]->getg2mu2A();
	    }
	  else 
	    g2mu2A = 0;
	  
	  if(posB>0 && posB<(N-1)*N+N-1)
	    {
	      g2mu2B = lat->cells[posB]->getg2mu2B();
	    }
	  else
	    g2mu2B = 0;
	  
	  if(g2mu2B>=g2mu2A) 
	    {
	      averageQs2min2 += g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg();
	    }
	  else
	    {
	      averageQs2min2 += g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg();
	    }
	  
	  // if( param->getWriteOutputs() == 1 )
	  //   fout << x << " " << y << " " << " " << outvalue << endl;
	  
	  for (int i = 0; i<A1; i++) 
	    {
	      xm = nucleusA.at(i).x;
	      ym = nucleusA.at(i).y;
	      r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
	      
	      if(r<sqrt(0.1*param->getSigmaNN()/PI) && nucleusA.at(i).collided==1)
		{
		  check=1;
		}
	    }
	  
	  for (int i = 0; i<A2; i++) 
	    {
	      xm = nucleusB.at(i).x;
	      ym = nucleusB.at(i).y;
	      r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
	      
	      if(r<sqrt(0.1*param->getSigmaNN()/PI) && nucleusB.at(i).collided==1 && check==1)
		check=2;
	    }
	  
	  if(check==2)
	    {
	      if(g2mu2B>g2mu2A) 
		{
		  averageQs += sqrt(g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg());
		  averageQs2 += g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg();
		  averageQs2min += g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg();
		}
	      else
		{
		  averageQs += sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg());
		  averageQs2 += g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg();
		  averageQs2min += g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()/a/a*hbarc*hbarc*param->getg()*param->getg();
		}
	      averageQs2Avg += (g2mu2B*param->getQsmuRatioB()*param->getQsmuRatioB()+g2mu2A*param->getQsmuRatio()*param->getQsmuRatio())
		/2./a/a*hbarc*hbarc*param->getg()*param->getg();
	      count++;
	    }
	  //else
	  //	    fout << x << " " << y << " " << 0. << endl;
	  
	  // compute T_pp

	  Tpp += lat->cells[pos]->getTpB()*lat->cells[pos]->getTpA()*a*a/hbarc/hbarc/hbarc/hbarc; // now this quantity is in fm^-2
	  // remember: Tp is in GeV^2
	}
      //      fout << endl;
    }
  //  fout.close();
  
  averageQs/=static_cast<double>(count);
  averageQs2/=static_cast<double>(count);
  averageQs2Avg/=static_cast<double>(count);
  averageQs2min/=static_cast<double>(count);
  
  param->setAverageQs(sqrt(averageQs2));
  param->setAverageQsAvg(sqrt(averageQs2Avg));
  param->setAverageQsmin(sqrt(averageQs2min));
  
  param->setTpp(Tpp);

  cout << "N_part=" << Npart << endl;
  cout << "N_coll=" << Ncoll << endl;
  cout << "T_pp(" << param->getb() << " fm) = " << Tpp << " 1/fm^2" << endl;
  cout << "Q_s^2(max) S_T = " << averageQs2*a*a/hbarc/hbarc*static_cast<double>(count) << endl;
  cout << "Q_s^2(avg) S_T = " << averageQs2Avg*a*a/hbarc/hbarc*static_cast<double>(count) << endl;
  cout << "Q_s^2(min) S_T = " << averageQs2min2*a*a/hbarc/hbarc << endl;
  
  cout << "Area = " << a*a*count << " fm^2" << endl;
  
  cout << "Average Qs(max) = " << param->getAverageQs() << " GeV" << endl;
  cout << "Average Qs(avg) = " << param->getAverageQsAvg() << " GeV" << endl;
  cout << "Average Qs(min) = " << param->getAverageQsmin() << " GeV" << endl;
  
  cout << "resulting Y(Qs(max)*" << param->getxFromThisFactorTimesQs() << ") = " 
  << log(0.01/(param->getAverageQs()*param->getxFromThisFactorTimesQs()/param->getRoots())) << endl;
  cout << "resulting Y(Qs(avg)*" << param->getxFromThisFactorTimesQs() << ") = "  << log(0.01/(param->getAverageQsAvg()*param->getxFromThisFactorTimesQs()/param->getRoots())) << endl;
  cout << "resulting Y(Qs(min)*" << param->getxFromThisFactorTimesQs() << ") =  " << log(0.01/(param->getAverageQsmin()*param->getxFromThisFactorTimesQs()/param->getRoots())) << endl;
   
  cout << "Color charge densities for nucleus A and B set. " << endl;
  
 
  // stringstream strQs2ST_name;
  // strQs2ST_name << "Qs2ST" << param->getMPIRank() << ".dat";
  // string Qs2ST_name;
  // Qs2ST_name = strQs2ST_name.str();
  // ofstream foutQ(Qs2ST_name.c_str(),ios::out);
  // foutQ << b << " " << Npart << " " << averageQs2min2*a*a/hbarc/hbarc << " " << a*a*count << endl;
  // foutQ.close();
     


  if(param->getRunningCoupling() && param->getRunWithkt()==0)
    {
      if(param->getRunWithQs()==2)
	{
	  cout << "running with " << param->getRunWithThisFactorTimesQs() << " Q_s(max)" << endl;
	  alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2)); // 3 flavors
	  cout << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " Qs_max)=" << alphas << endl;
	}
      else if(param->getRunWithQs()==0)
	{
	  cout << "running with " << param->getRunWithThisFactorTimesQs() << " Q_s(min)" << endl;
	  alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2)); // 3 flavors
	  cout << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " Qs_min)=" << alphas << endl;
	}
      else if(param->getRunWithQs()==1)
	{
	  cout << "running with " << param->getRunWithThisFactorTimesQs() << " <Q_s>" << endl;
	  alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2)); // 3 flavors
	  cout << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " <Qs>)=" << alphas << endl;
	}
    }
  else if(param->getRunningCoupling() && param->getRunWithkt()==1)
    {
      cout << "Multiplicity with running alpha_s(k_T)" << endl;
    }
  else
    {
      cout << "Using fixed alpha_s" << endl;
      alphas = param->getg()*param->getg()/4./PI;
    }
 
  if(param->getAverageQs() > 0 && param->getAverageQsAvg()>0 && averageQs2>0  && param->getAverageQsmin()>0 && averageQs2Avg>0 && alphas>0 && Npart>=2)
    param->setSuccess(1);
 
  param->setalphas(alphas);
  
  stringstream strup_name;
  strup_name << "usedParameters" << param->getMPIRank() << ".dat";
  string up_name;
  up_name = strup_name.str();

  ofstream fout1(up_name.c_str(),ios::app);
  fout1 << " " << endl;
  fout1 << " Output by setColorChargeDensity in Init.cpp: " << endl;
  fout1 << " " << endl;
  fout1 << "b = " << b << " fm" << endl;
  fout1 << "Npart = " << Npart << endl;
  fout1 << "Ncoll = " << Ncoll << endl;
  if(param->getRunningCoupling())
    {
      if(param->getRunWithQs()==2)
	fout1 << "<Q_s>(max) = " << param->getAverageQs() << endl;
      else if(param->getRunWithQs()==1)
	fout1 << "<Q_s>(avg) = " << param->getAverageQsAvg() << endl;
      else if(param->getRunWithQs()==0)
	fout1 << "<Q_s>(min) = " << param->getAverageQsmin() << endl;
      fout1 << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " <Q_s>) = " << param->getalphas() << endl;
    }
  else
    fout1 << "using fixed coupling alpha_s=" << param->getalphas() << endl;
  //  fout1 << "Q_s ~ x^-" << param->getxExponent() << endl;
  fout1.close();

 //  // output gmu 
 //  ofstream foutB("gmuB.dat",ios::out); 
 //  for(int ix=0; ix<N; ix++) // loop over all positions
 //    {
 //      for(int iy=0; iy<N; iy++)
 // 	{
 // 	  pos = ix*N+iy;
 // 	  outvalue = sqrt(lat->cells[pos]->getg2mu2B());
 // 	  foutB << ix << " " << iy << " " << outvalue << endl;
 // 	}
 //      foutB << endl;
 //    }
 //  foutB.close();
 // // output gmu 
 //  ofstream foutA("gmuA.dat",ios::out); 
 //  for(int ix=0; ix<N; ix++) // loop over all positions
 //    {
 //      for(int iy=0; iy<N; iy++)
 // 	{
 // 	  pos = ix*N+iy;
 // 	  outvalue = sqrt(lat->cells[pos]->getg2mu2A());
 // 	  foutA << ix << " " << iy << " " << outvalue << endl;
 // 	}
 //      foutA << endl;
 //    }
 //  foutA.close();
}



void Init::setV(Lattice *lat, Group* group, Parameters *param, Random* random, Glauber *glauber)
{
  int pos;
  int N = param->getSize();
  int Ny=param->getNy();
  int Nc = param->getNc();
  int Nc2m1 = Nc*Nc-1;
  int nn[2];
  nn[0]=N;
  nn[1]=N; 

  double corr;
  double g2mu;
  double L = param->getL();
  double x, y;
  double a = L/N; // lattice spacing in fm
  double kt2, kx, ky;
  double m = param->getm(); //GeV
  m=m*a/hbarc;
  double temp3;
  Matrix** rhoA;
  Matrix** rhoB;
  Matrix** AA;
  Matrix** AB;
  Matrix temp(Nc,1.);
  Matrix temp2(Nc,0.);
  Matrix tempNew(Nc,0.);
  Matrix Udag(Nc);
  Matrix zero(Nc,0.);
  Matrix one(Nc,1.);
 
  //  rhoA = new Matrix*[N*N];
  //rhoB = new Matrix*[N*N];

  //AA = new Matrix*[N*N];
  //AB = new Matrix*[N*N];


  


  complex<double>** rhoACoeff;
  complex<double>** rhoBCoeff;

  rhoACoeff = new complex<double>*[Nc2m1];
  rhoBCoeff = new complex<double>*[Nc2m1];

  for(int i=0; i<Nc2m1; i++)
    {
      rhoACoeff[i] = new complex<double>[N*N];
      rhoBCoeff[i] = new complex<double>[N*N];
    }
  

  // loop over longitudinal direction
  for(int k=0; k<Ny; k++)
    {

      //  for(int i=0; i<N*N; i++)
      //	{
	  // rhoA[i] = new Matrix(Nc,0.);
      	  //rhoB[i] = new Matrix(Nc,0.);
	  //	  AA[i] = new Matrix(Nc,0.);
      	  //AB[i] = new Matrix(Nc,0.);
      //	}
   

     // compute \rho
      for (int i=0; i<N; i++)
	{
	  for (int j=0; j<N; j++)
	    {
	      pos = i*N+j;
	      for(int n=0; n<Nc2m1; n++)
		{
		  g2mu = param->getg()*sqrt(lat->cells[pos]->getg2mu2A()/static_cast<double>(Ny));
		  rhoACoeff[n][pos] = g2mu*random->Gauss();
		  //	  *rhoA[pos]+=rhoACoeff[n][pos]*group->getT(n);
		  g2mu = param->getg()*sqrt(lat->cells[pos]->getg2mu2B()/static_cast<double>(Ny));
		  rhoBCoeff[n][pos] = g2mu*random->Gauss();
		  //*rhoB[pos]+=rhoBCoeff[n][pos]*group->getT(n);
		}
	    }
	}
   

      // --------


      //// output rho 
      //      if(k==0)
      // 	{
      	  // ofstream foutr("RhoOne.txt",ios::out); 
      	  // for(int ix=0; ix<N; ix++)
      	  //   {
      	  //     for(int iy=0; iy<N; iy++) // loop over all positions
      	  // 	{
      	  // 	  pos = ix*N+iy;
      	  // 	  foutr << ix << " " << iy << " " << setprecision(17) << rhoACoeff[0][pos] << " " << rhoACoeff[1][pos] << " " << rhoACoeff[2][pos] << " " << rhoACoeff[3][pos] << " " << rhoACoeff[4][pos] << " " << rhoACoeff[5][pos] << " " << rhoACoeff[6][pos] << " " << rhoACoeff[7][pos] <<  endl;
      	  // 	}
      	  //   }
      	  // foutr.close();
      	  // ofstream foutr2("RhoTwo.txt",ios::out); 
      	  // for(int ix=0; ix<N; ix++)
      	  //   {
      	  //     for(int iy=0; iy<N; iy++) // loop over all positions
      	  // 	{
      	  // 	  pos = ix*N+iy;
      	  // 	  foutr2 << ix << " " << iy << " " << setprecision(17) << rhoBCoeff[0][pos] << " " << rhoBCoeff[1][pos] << " " << rhoBCoeff[2][pos] << endl;
      	  // 	}
      	  //   }
      	  // foutr2.close();
	  //      	  cout << "wrote rho's into file." << endl;
      
      //	}

      // delete [] rhoACoeff;
      //delete [] rhoBCoeff;

	  
      // for (int i=0; i<N; i++)
      // 	{
      // 	  for (int j=0; j<N; j++)
      // 	    {
      // 	      pos = i*N+j;
      // 	      if (j==N/2)
      // 		{
      // 		  cout << "k=" << k << ", i=" << i << ", rhoB=" << *rhoB[pos] << endl;
      // 		}
      // 	    }
      // 	}
      
      // --------


      // Fourier transform rho
      // fft->fftn(rhoA,AA,nn,2,1);
      // fft->fftn(rhoB,AB,nn,2,1);

  
      // begin new
      for(int n=0; n<Nc2m1; n++)
	{
	  fft->fftnComplex(rhoACoeff[n],rhoACoeff[n],nn,2,1);
	  fft->fftnComplex(rhoBCoeff[n],rhoBCoeff[n],nn,2,1);
	}
      // end new
 
      // compute A^+
      for (int i=0; i<N; i++)
	{
	  for (int j=0; j<N; j++)
	    {
	      pos = i*N+j;
	      kx = 2.*param->getPi()*(-0.5+static_cast<double>(i)/static_cast<double>(N));
	      ky = 2.*param->getPi()*(-0.5+static_cast<double>(j)/static_cast<double>(N));
	      kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice momentum
	      if(m==0)
		{
		  if(kt2!=0)
		    {

		      // begin new
		      for(int n=0; n<Nc2m1; n++)
			{
			  rhoACoeff[n][pos] =  rhoACoeff[n][pos]*(1./(kt2));
			  rhoBCoeff[n][pos] =  rhoBCoeff[n][pos]*(1./(kt2));
			}
		      // end new
		      //*AA[pos] = *AA[pos]*(1./(kt2)); // rho contains A to save memory // check a's, is dimensionless as it should be
		      //*AB[pos] = *AB[pos]*(1./(kt2)); // rho contains A to save memory
		    }
		  else
		    {
		      // begin new
		      for(int n=0; n<Nc2m1; n++)
			{
			  rhoACoeff[n][pos] = 0.;
			  rhoBCoeff[n][pos] = 0.;
			}
		      // end new
		      //*AA[pos] = zero; // rho contains A to save memory // check a's, is dimensionless as it should be
		      //*AB[pos] = zero; // rho contains A to save memory
		    }
		}
	      else
		{
		  // begin new
		  for(int n=0; n<Nc2m1; n++)
		    {
		      rhoACoeff[n][pos] =  rhoACoeff[n][pos]*(1./(kt2+m*m));
		      rhoBCoeff[n][pos] =  rhoBCoeff[n][pos]*(1./(kt2+m*m));
	
		    }
		  // end new
		  //		  *AA[pos] = *AA[pos]*(1./(kt2+m*m)); // rho contains A to save memory // check a's, is dimensionless as it should be
		  //*AB[pos] = *AB[pos]*(1./(kt2+m*m)); // rho contains A to save memory
		}
	    }
	}
      
      // Fourier transform back A^+
      //fft->fftn(AA,AA,nn,2,-1);
      //fft->fftn(AB,AB,nn,2,-1);

      // begin new
      for(int n=0; n<Nc2m1; n++)
	{
	  fft->fftnComplex(rhoACoeff[n],rhoACoeff[n],nn,2,-1);
	  fft->fftnComplex(rhoBCoeff[n],rhoBCoeff[n],nn,2,-1);
	}
      // end new

      

      // --------

      // // output phi
      //      if(k==0)
      // 	{
      	  // ofstream foutph("PhiOne.txt",ios::out); 
      	  // for(int ix=0; ix<N; ix++)
      	  //   {
      	  //     for(int iy=0; iy<N; iy++) // loop over all positions
      	  // 	{
      	  // 	  pos = ix*N+iy;
	  // 	  //      		  foutph << ix << " " << iy << " " <<  setprecision(17) << 2.* ((*AA[pos])*group->getT(0)).trace().real() << " "
	  // 	  //	 << 2.* ((*AA[pos])*group->getT(1)).trace().real() << " " << 2.* ((*AA[pos])*group->getT(2)).trace().real() << endl;
      	  // 	  foutph << ix << " " << iy << " " <<  setprecision(17) 
	  // 		 << (*AA[pos]).MatrixToString() << endl;
		  
      	  // 	}
      	  //   }
   


   // 	  foutph.close();
      // 	  ofstream foutph2("PhiTwo.txt",ios::out); 
      // 	  for(int ix=0; ix<N; ix++)
      // 	    {
//       // 	      for(int iy=0; iy<N; iy++) // loop over all positions
//       // 		{
//       // 		  pos = ix*N+iy;
//       // 		  foutph2 << ix << " " << iy << " " << setprecision(17) << 2.* ((*AB[pos])*group->getT(0)).trace().real() << " "
//       // 			  << 2.* ((*AB[pos])*group->getT(1)).trace().real() << " " << 2.* ((*AB[pos])*group->getT(2)).trace().real() << endl;

//       // 		}
//       // 	    }
//       // 	  foutph2.close();
//       	}

//       // --------

      
// //       double distanceA, distanceB;
// //       double BG=param->getBG();
// //       int check=0;

//       // clean up noise where it shouldn't be - far away from interaction region: causes problems with convergence when solving for U3
//  //      for (int i=0; i<N; i++)
// // 	{
// // 	  for (int j=0; j<N; j++)
// // 	    {
// // 	      pos = i*N+j;
	     
	   
// // 	      if(lat->cells[pos]->getTpA() < 0.0019)
// // 		{
// // 		  check=1;
// // 		}
	      
// // 	      if(lat->cells[pos]->getTpB() < 0.0019 && check==1)
// // 		{		
// // 		  check=2;
// // 		}
	  
// // 	      //if(*rhoB[pos]==zero && *rhoA[pos]==zero && check==2)
// // 	      if( check==2)
// // 		{
// // 		  *AA[pos] = zero;
// // 		  *AB[pos] = zero;
// // 		}
// // 	    }
      //	}      


      
      
      // begin new
      // for (int i=0; i<N; i++)
      // 	{
      // 	  for (int j=0; j<N; j++)
      // 	    {
      // 	      pos = i*N+j;
      // 	      *rhoA[pos]=zero;
      // 	      *rhoB[pos]=zero;
      // 	    }
      // 	}


      // for (int i=0; i<N; i++)
      // 	{
      // 	  for (int j=0; j<N; j++)
      // 	    {
      // 	      pos = i*N+j;
      // 	      for(int n=0; n<Nc2m1; n++)
      // 		{
      // 		  *rhoA[pos]+=rhoACoeff[n][pos]*group->getT(n);
      // 		  *rhoB[pos]+=rhoBCoeff[n][pos]*group->getT(n);
      // 		}
      // 	    }
      // 	}
    

      // for (int i=0; i<N; i++)
      // 	{
      // 	  for (int j=0; j<N; j++)
      // 	    {
      // 	      pos = i*N+j;
	      
      // 	      cout << "new: \n" <<*rhoA[pos] << " \n old:\n" << *AA[pos] << endl;
      // 	    }
      // 	}
      
      //      exit(1);
      // end new
  
      double in[8];
      vector <complex<double> > U;
      // compute U

      // new method
      //  clock_t start;
      //double duration;
      //start = clock();
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*N+j;
	      
	      for (int a=0; a<Nc2m1; a++)
		{
		  in[a] = -(rhoACoeff[a][pos]).real(); // expmCoeff wil calculate exp(i in[a]t[a]), so just multiply by -1 (not -i)
		}
	      
	      U = temp2.expmCoeff(in, Nc);

	      // if(U[0].real()!=U[0].real())
	      // 	{
	      // 	  cout << "PROBLEM 1" << " " << g2mu << endl;
	      // 	  for (int a=0; a<Nc2m1; a++)
	      // 	    {
	      // 	      cout << in[a] << endl; 
	      // 	    }
	      // 	}
	      
	      tempNew = U[0]*one + U[1]*group->getT(0) + U[2]*group->getT(1) + U[3]*group->getT(2) + U[4]*group->getT(3) + 
		U[5]*group->getT(4) + U[6]*group->getT(5) + U[7]*group->getT(6) + U[8]*group->getT(7);

	      temp = tempNew * lat->cells[pos]->getU();
	      // set U
	      lat->cells[pos]->setU(temp);
	      
	      
	      for (int a=0; a<Nc2m1; a++)
		{
		  in[a] = -(rhoBCoeff[a][pos]).real(); // expmCoeff wil calculate exp(i in[a]t[a]), so just multiply by -1 (not -i)
		}
	      
	      U = temp2.expmCoeff(in, Nc);
	      
	      // if(U[0].real()!=U[0].real())
	      // 	{
	      // 	  cout << "PROBLEM" << endl;
	      // 	}
	      

	      tempNew = U[0]*one + U[1]*group->getT(0) + U[2]*group->getT(1) + U[3]*group->getT(2) + U[4]*group->getT(3) + 
		U[5]*group->getT(4) + U[6]*group->getT(5) + U[7]*group->getT(6) + U[8]*group->getT(7);

	      temp = tempNew * lat->cells[pos]->getU2();

	      // set U
	      lat->cells[pos]->setU2(temp);
	    }
	}
      //      duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
      // cout << "duration of new method = " << duration << endl;
     
      // start = clock();
      // // old method
      // for (int i=0; i<nn[0]; i++)
      // 	{
      // 	  for (int j=0; j<nn[1]; j++)
      // 	    {
      // 	      pos = i*N+j;
	      
      // 	      //multiply by -i:
      // 	      temp2=complex<double>(0.,-1.)*(*AA[pos]);
      // 	      temp2.expm();

      // 	      temp = temp2 * lat->cells[pos]->getU();
      // 	      // set U
      // 	      lat->cells[pos]->setU(temp);
	      
      // 	      if(i==0 && j ==0)
      // 		cout << "OLD: " << endl << temp2 << endl;

      // 	      //multiply by -i:
      // 	      temp2=complex<double>(0.,-1.)*(*AB[pos]);
      // 	      temp2.expm();
  
      // 	      temp = temp2 * lat->cells[pos]->getU2();

      // 	      // set U
      // 	      lat->cells[pos]->setU2(temp);
      // 	    }
      // 	}

      // duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
      // cout << "duration of old method = " << duration << endl;
      
      

    }//Ny loop

      for(int ic=0; ic<Nc2m1; ic++)
      	{
       	  delete [] rhoACoeff[ic];
       	  delete [] rhoBCoeff[ic];
       	}
      delete [] rhoACoeff;
      delete [] rhoBCoeff;


   //output of U loop over index (real imag)

  // --------

  
        
  // // output U
  // stringstream strVOne_name;
  // strVOne_name << "V1-" << param->getMPIRank() << ".txt";
  // string VOne_name;
  // VOne_name = strVOne_name.str();

  // ofstream foutU(VOne_name.c_str(),ios::out); 
  // foutU.precision(15);

  // for(int ix=0; ix<N; ix++)
  //   {
  //     for(int iy=0; iy<N; iy++) // loop over all positions
  // 	{
  // 	  pos = ix*N+iy;
  // 	  foutU << ix << " " << iy << " "  << (lat->cells[pos]->getU()).MatrixToString() << endl;
  // 	}
  //     foutU << endl;
  //   }
  // foutU.close();

  // cout<<"wrote " << strVOne_name.str() <<endl;
  
  // stringstream strVTwo_name;
  // strVTwo_name << "V2-" << param->getMPIRank() << ".txt";
  // string VTwo_name;
  // VTwo_name = strVTwo_name.str();

  // ofstream foutU2(VTwo_name.c_str(),ios::out); 
  // foutU2.precision(15);
  // for(int ix=0; ix<N; ix++)
  //   {
  //     for(int iy=0; iy<N; iy++) // loop over all positions
  // 	{
  // 	  pos = ix*N+iy;
  // 	  foutU2 << ix << " " << iy << " "  << (lat->cells[pos]->getU2()).MatrixToString() << endl;
  // 	}
  //     foutU2 << endl;
  //   }
  // foutU2.close();
  
  // cout<<"wrote " << strVTwo_name.str() <<endl;
 
  // --------


  //  delete[] rhoA;
  //delete[] rhoB;
  //delete[] AA;
  //delete[] AB;
  
  cout << " Wilson lines V_A and V_B set on rank " << param->getMPIRank() << ". " << endl; 

  //  MPI::Finalize();
  // exit(1);

  // output correlator 
  //  Udag = lat->cells[N*N/2+N/2]->getU();
  //Udag.conjg();
  
  if(param->getWriteOutputs() == 1)
    {
      // stringstream strVdagV_name;
      // strVdagV_name << "VdagV" << param->getMPIRank() << ".dat";
      // string VdagV_name;
      // VdagV_name = strVdagV_name.str();

      // ofstream fout(VdagV_name.c_str(),ios::out); 
      
      // for (int i=0; i<nn[0]; i++)
      // 	{
      // 	  for (int j=0; j<nn[1]; j++)
      // 	    {
      // 	      pos = i*N+j;
      // 	      x = -L/2.+a*i;
      // 	      y = -L/2.+a*j;
	      
      // 	      corr = real((Udag*lat->cells[pos]->getU()).trace())/static_cast<double>(Nc);
	      
      // 	      //  	      fout << x << ", " << y << ", " << corr << ", " << endl;
      // 	    }
      // 	  //  	  fout << endl;
      // 	}
      
      // fout.close();
    }
}


void Init::readV(Lattice *lat, Group* group, Parameters *param)
{
  int pos;
  int N = param->getSize();
  int Ny=param->getNy();
  int Nc = param->getNc();
  int Nc2m1 = Nc*Nc-1;
  int nn[2];
  nn[0]=N;
  nn[1]=N; 

  Matrix temp(Nc,1.);
 
  double Re[9], Im[9];
  double dummy;
  
  stringstream strVOne_name;
  strVOne_name << "V1Y-0.4.txt";
  string VOne_name;
  VOne_name = strVOne_name.str();
  ifstream finV1(VOne_name.c_str(),ios::in); 

  if (!finV1)
    {
      cerr << "File " << VOne_name << " not found. Exiting." << endl;
      exit(1);
    }	 

  cout << "Reading Wilson line from file " << VOne_name << " ..." << endl;


  // set V for nucleus A
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{

	  finV1 >> dummy >> dummy 
		>> Re[0] >> Im[0] >> Re[1] >> Im[1] >> Re[2] >> Im[2] 
		>> Re[3] >> Im[3] >> Re[4] >> Im[4] >> Re[5] >> Im[5] 
		>> Re[6] >> Im[6] >> Re[7] >> Im[7] >> Re[8] >> Im[8];
	  
	  
	  temp.set(0,0,complex<double>(Re[0],Im[0]));
	  temp.set(0,1,complex<double>(Re[1],Im[1]));
	  temp.set(0,2,complex<double>(Re[2],Im[2]));
	  temp.set(1,0,complex<double>(Re[3],Im[3]));
	  temp.set(1,1,complex<double>(Re[4],Im[4]));
	  temp.set(1,2,complex<double>(Re[5],Im[5]));
	  temp.set(2,0,complex<double>(Re[6],Im[6]));
	  temp.set(2,1,complex<double>(Re[7],Im[7]));
	  temp.set(2,2,complex<double>(Re[8],Im[8]));
	  
	  pos = i*N+j;	  
	  lat->cells[pos]->setU(temp);	  

	}
    }

  finV1.close();
  
  stringstream strVTwo_name;
  strVTwo_name << "V3Y0.4.txt";
  string VTwo_name;
  VTwo_name = strVTwo_name.str();
  ifstream finV2(VTwo_name.c_str(),ios::in); 
 
  if (!finV2)
    {
      cerr << "File " << VTwo_name << " not found. Exiting." << endl;
      exit(1);
    }	 

  cout << "Reading Wilson line from file " << VTwo_name << " ..." << endl;

  // set V for nucleus B
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{

	  finV2 >> dummy >> dummy 
		>> Re[0] >> Im[0] >> Re[1] >> Im[1] >> Re[2] >> Im[2] 
		>> Re[3] >> Im[3] >> Re[4] >> Im[4] >> Re[5] >> Im[5] 
		>> Re[6] >> Im[6] >> Re[7] >> Im[7] >> Re[8] >> Im[8];
	  
	  
	  temp.set(0,0,complex<double>(Re[0],Im[0]));
	  temp.set(0,1,complex<double>(Re[1],Im[1]));
	  temp.set(0,2,complex<double>(Re[2],Im[2]));
	  temp.set(1,0,complex<double>(Re[3],Im[3]));
	  temp.set(1,1,complex<double>(Re[4],Im[4]));
	  temp.set(1,2,complex<double>(Re[5],Im[5]));
	  temp.set(2,0,complex<double>(Re[6],Im[6]));
	  temp.set(2,1,complex<double>(Re[7],Im[7]));
	  temp.set(2,2,complex<double>(Re[8],Im[8]));
	  

	  pos = i*N+j;	  
	  lat->cells[pos]->setU2(temp);
	}
    }
  
  finV2.close();

  cout << " Wilson lines V_A and V_B set on rank " << param->getMPIRank() << ". " << endl; 

}






void Init::init(Lattice *lat, Group *group, Parameters *param, Random *random, Glauber *glauber, int READFROMFILE)
{
  int maxIterations = 100000;
  int N = param->getSize();
  int Ny= param->getNy();
  int Nc = param->getNc();
  int bins = param->getSize();
  int ir;
  int count[bins];
  int Nc2m1 = Nc*Nc-1;
  int nn[2];
  int pos, pos1, pos2, pos3, posx, posy, posxm, posym, posxmym;
  int counts, countMe;
  int checkConvergence;
  int alphaCheck;
  int bShift; // number of cells to be shifted by due to impact parameter
  int posU;

  double dNc = static_cast<double>(Nc);
  double Fold;
  double Fnew;
  double r;
  double L = param->getL();
  double x;
  double y;
  double a = L/N; // lattice spacing in fm
  cout << "Initializing fields ... " << endl;
  param->setRnp(0.);
  
  double b;
  double bmin=param->getbmin();
  double bmax=param->getbmax();
  
  double xb = random->genrand64_real1(); // uniformly distributed random variable                                                               
  
  if(param->getUseNucleus() == 0) // use b=0 fm for the constant g^2 mu case
    {
      param->setSuccess(1);
      b=0.;
      cout << "Setting b=0 for constant color charge density case." << endl;
    }
  else 
    {
      if(param->getLinearb()==1) // use a linear probability distribution for b if we are doing nuclei
	{
	  cout << "Sampling linearly distributed b between " << bmin << " and " << bmax << "fm. Found ";
	  b = sqrt((bmax*bmax-bmin*bmin)*xb+bmin*bmin);
	}
      else // use a uniform distribution instead
	{
	  cout << "Sampling uniformly distributed b between " << bmin << " and " << bmax << "fm. Found ";
	  b = (bmax-bmin)*xb+bmin;
	}
    }

  param->setb(b);
  cout << "b=" << b << " fm." << endl;

  double Qs2G;
  double temp3;
  double g2mu;
  double trATB[N*N];
  double dr=a;
  double rtrAT2[bins];
  double epsilon;
  double m = param->getm(); //GeV
  double lambda;
  double trATA[N*N];
  double avgEps;
  double avgEpsMag;
  double avgEpsEl;

  m=m*a/hbarc;
  //  cout << "m_lat =" << m << endl; 

  complex<double>* M;
  complex<double>* F;
  complex<double>* result;
  complex<double>* alpha;
  complex<double>* alphaSave;
  vector <complex<double> > Dalpha;
  Dalpha.reserve(Nc2m1);

  M = new complex<double>[Nc2m1*Nc2m1];
  F = new complex<double>[Nc2m1];
  result = new complex<double>[Nc2m1];
  alpha = new complex<double>[Nc2m1];
  alphaSave = new complex<double>[Nc2m1];
  //Dalpha = new complex<double>[Nc2m1];

  // read Q_s^2 from file
  if(param->getUseNucleus() == 1)
    {
      readNuclearQs(param);
    }

  Matrix temp(Nc,1.);
  Matrix tempNew(Nc,1.);
  double in[8];
  vector <complex<double> > U;
  //   Matrix U(Nc,1.);
  Matrix temp2(Nc,0.);
  Matrix expAlpha(Nc,0.);
  Matrix expNegAlpha(Nc,0.);
  Matrix Ux(int(Nc),0.);
  Matrix Uy(int(Nc),0.);
  Matrix Ux1(int(Nc),0.);
  Matrix Uy1(int(Nc),0.);
  Matrix Ux2(int(Nc),0.);
  Matrix Uy2(int(Nc),0.);
  Matrix UD(int(Nc),0.);
  Matrix UDx(int(Nc),0.);
  Matrix UDy(int(Nc),0.);
  Matrix UDx1(int(Nc),0.);
  Matrix UDy1(int(Nc),0.);
 
  Matrix Uplaq(int(Nc),0.);
  Matrix Uplaq1(int(Nc),0.);
  Matrix Uplaq2(int(Nc),0.);
  Matrix Uplaq3(int(Nc),0.);
  Matrix Uplaq4(int(Nc),0.);

  Matrix UD2(int(Nc),0.);
  Matrix UDx2(int(Nc),0.);
  Matrix UDy2(int(Nc),0.);

  Matrix Ax(int(Nc),0.);
  Matrix Ay(int(Nc),0.);
  Matrix Ax1(int(Nc),0.);
  Matrix Ay1(int(Nc),0.);
  Matrix AT(int(Nc),0.);

  Matrix Ax2(int(Nc),0.);
  Matrix Ay2(int(Nc),0.);
  Matrix AT2(int(Nc),0.);

  // field strength tensor
  Matrix Fxy(int(Nc),0.);
  Matrix Fyx(int(Nc),0.);

  Matrix AM(int(Nc),0.);
  Matrix AP(int(Nc),0.);

  Matrix AxUpY(int(Nc),0.);
  Matrix AyUpY(int(Nc),0.);

  Matrix Aeta2(int(Nc),0.);
  Matrix Ux1pUx2(int(Nc),0.);
  Matrix UDx1pUDx2(int(Nc),0.);
  Matrix Uy1pUy2(int(Nc),0.);
  Matrix UDy1pUDy2(int(Nc),0.);
  Matrix Ux1mUx2(int(Nc),0.);
  Matrix UDx1mUDx2(int(Nc),0.);
  Matrix Uy1mUy2(int(Nc),0.);
  Matrix UDy1mUDy2(int(Nc),0.);

  Matrix one(Nc,1.);
  Matrix zero(Nc,0.);

  

  nn[0]=N;
  nn[1]=N;

  // sample nucleon positions
  nucleusA.clear();
  nucleusB.clear();

  // to read Wilson lines from file (e.g. after JIMWLK evolution for the 3DGlasma)
  if(READFROMFILE)
    {
      readV(lat, group, param);
      param->setSuccess(1);
    }
  // to generate your own Wilson lines
  else
    {
      if(param->getUseNucleus() == 1)
	sampleTA(param, random, glauber);                           // populate the lists nucleusA and nucleusB with position data of the 
      
      // set color charge densities
      setColorChargeDensity(lat, param, random, glauber);
      
      if(param->getUseNucleus() == 1 && param->getUseFixedNpart()!=0 && param->getNucleonPositionsFromFile()!=1)
	{
	  if(param->getNpart()!=param->getUseFixedNpart())
	    {
	      while(param->getNpart()!=param->getUseFixedNpart())
		{
		  cout << "resampling... desired Npart=" << param->getUseFixedNpart() << endl;
		  nucleusA.clear();
		  nucleusB.clear();
		  
		  xb = random->genrand64_real1(); // uniformly distributed random variable                                                               
		  
		  if(param->getLinearb()==1) // use a linear probability distribution for b if we are doing nuclei
		    {
		      cout << "Sampling linearly distributed b between " << bmin << " and " << bmax << "fm." << endl;
		      b = sqrt((bmax*bmax-bmin*bmin)*xb+bmin*bmin);
		    }
		  else // use a uniform distribution instead
		    {
		      cout << "Sampling uniformly distributed b between " << bmin << " and " << bmax << "fm." << endl;
		      b = (bmax-bmin)*xb+bmin;
		    }
		  
		  param->setb(b);
		  cout << "Using b=" << b << " fm" << endl;
		  
		  sampleTA(param, random, glauber);                           // populate the lists nucleusA and nucleusB with position data of the 
		  setColorChargeDensity(lat, param, random, glauber);
		}
	    }
	  cout << "Using fixed Npart=" << param->getNpart() << endl;
	}
      
      
      if(param->getSuccess()==0)
	{
	  delete M;
	  delete F;
	  delete result;
	  delete alpha;
	  delete alphaSave;
	  //delete Dalpha;
	  cout << "No collision happened on rank " << param->getMPIRank() << ". Restarting with new random number..." << endl;
	  return;
	}
      
     
      
      // sample color charges and find Wilson lines V_A and V_B
      setV(lat, group, param, random, glauber);
      
    }
  
  // output Wilson lines

  // ofstream fout("V.dat",ios::out); 
  // fout << "# Wilson lines. Format: x and y coordinate in [fm], then 3x3 matrix: Re(V_{i,j}) Im(V_{i,j}), i is the row, j the column, j is the inner loop, i.e., the order is Re(V_{0,0}) Im(V_{0,0}) Re(V_{0,1}) Im(V_{0,1}) Re(V_{0,2}) Im(V_{0,2}) Re(V_{1,0}) Im(V_{1,0}) ..." << endl;

  // for (int i=0; i<nn[0]; i++)      //loops over all cells
  //   {
  //     for (int j=0; j<nn[1]; j++)      //loops over all cells
  // 	{
  // 	  pos = i*N+j;
  // 	  x = -L/2.+a*i;
  // 	  y = -L/2.+a*j;

  // 	  fout << x << " " << y << " " 
  // 	       << lat->cells[pos]->getU().getRe(0) << " " << lat->cells[pos]->getU().getIm(0) << " " 
  // 	       << lat->cells[pos]->getU().getRe(1) << " " << lat->cells[pos]->getU().getIm(1) << " "
  // 	       << lat->cells[pos]->getU().getRe(2) << " " << lat->cells[pos]->getU().getIm(2) << " "
  // 	       << lat->cells[pos]->getU().getRe(3) << " " << lat->cells[pos]->getU().getIm(3) << " "
  // 	       << lat->cells[pos]->getU().getRe(4) << " " << lat->cells[pos]->getU().getIm(4) << " "
  // 	       << lat->cells[pos]->getU().getRe(5) << " " << lat->cells[pos]->getU().getIm(5) << " "
  // 	       << lat->cells[pos]->getU().getRe(6) << " " << lat->cells[pos]->getU().getIm(6) << " "
  // 	       << lat->cells[pos]->getU().getRe(7) << " " << lat->cells[pos]->getU().getIm(7) << " "
  // 	       << lat->cells[pos]->getU().getRe(8) << " " << lat->cells[pos]->getU().getIm(8) 
  // 	       << endl;
	  
	  
  // 	}
  //   }

  // fout.close();      


  // compute Ux(3) Uy(3) after the collision
  //#pragma omp parallel for collapse(2)
  for (int i=0; i<nn[0]; i++)      //loops over all cells
    {
      for (int j=0; j<nn[1]; j++)      //loops over all cells
	{

	  pos = i*N+j;
	  pos1 = (i)*N+j;
	  posx = ((i+1)%N)*N+j; //<-- this may be a bottleneck
	  posy = (i)*N+(j+1)%N;
	  
	  // get Ux, Uy:
	  // Ux = U(r) UD(r+ax) (ax meaning one cell in the x direction away)
	  // if(i<N-1)
	  //   UDx=lat->cells[posx]->getU();
	  // else
	  //   UDx=lat->cells[(N-1)*N+j]->getU();

	  if(lat->cells[pos]->getU().trace()!=lat->cells[pos]->getU().trace())
	     {
	       lat->cells[pos]->setU(one);
	     }

	  if(lat->cells[pos1]->getU().trace()!=lat->cells[pos1]->getU().trace())
	     {
	       lat->cells[pos1]->setU(one);
	     }

	  if(lat->cells[posx]->getU().trace()!=lat->cells[posx]->getU().trace())
	     {
	       lat->cells[posx]->setU(one);
	     }

	  if(lat->cells[posy]->getU().trace()!=lat->cells[posy]->getU().trace())
	     {
	       lat->cells[posy]->setU(one);
	     }

	  UDx=lat->cells[posx]->getU();
	  
	  UDx.conjg();

	  lat->cells[pos]->setUx1(lat->cells[pos1]->getU()*UDx);
	  


	  // Uy = U(r) UD(r+ay) 
	  // if(j<N-1)
	  //   UDy=lat->cells[posy]->getU();
	  // else
	  //   UDy=lat->cells[(i)*N+N-1]->getU();
	

	  UDy=lat->cells[posy]->getU();

	  UDy.conjg();
	  
	  lat->cells[pos]->setUy1(lat->cells[pos1]->getU()*UDy);
	    
	  pos1 = (i)*N+j;
	  posx = ((i+1)%N)*N+j;
	  posy = (i)*N+(j+1)%N;
	  
	  // if(i<N-1)
	  //   {
	  //     UDx=lat->cells[posx]->getU2();
	  //   }
	  // else 
	  //   {
	  //     UDx=lat->cells[(N-1)*N+j]->getU2();
	  //   }


	  if(lat->cells[pos]->getU2().trace()!=lat->cells[pos]->getU2().trace())
	     {
	       
	       lat->cells[pos]->setU2(one);
	     }

	  if(lat->cells[pos1]->getU2().trace()!=lat->cells[pos1]->getU2().trace())
	     {
	       
	       lat->cells[pos1]->setU2(one);
	     }

	  if(lat->cells[posx]->getU2().trace()!=lat->cells[posx]->getU2().trace())
	     {
	       
	       lat->cells[posx]->setU2(one);
	     }

	  if(lat->cells[posy]->getU2().trace()!=lat->cells[posy]->getU2().trace())
	     {
	       
	       lat->cells[posy]->setU2(one);
	     }

	  UDx=lat->cells[posx]->getU2();

	  UDx.conjg();

	  //	  Ux2 = lat->cells[pos1]->getU2()*UDx;
	  lat->cells[pos]->setUx2(lat->cells[pos1]->getU2()*UDx);
	  
	  // if(j<N-1)
	  //   UDy=lat->cells[posy]->getU2();
	  // else
	  //   UDy=lat->cells[(i)*N+N-1]->getU2();
	  
          UDy=lat->cells[posy]->getU2();

	  UDy.conjg();	 

	  //  Uy2 = lat->cells[pos1]->getU2()*UDy;
	  lat->cells[pos]->setUy2(lat->cells[pos1]->getU2()*UDy);
	  
	  // -----------------------------------------------------------------
	  // from Ux(1,2) and Uy(1,2) compute Ux(3) and Uy(3):

	  UDx1 = lat->cells[pos]->getUx1();
	  UDx2 = lat->cells[pos]->getUx2();
	  Ux1pUx2 = UDx1+UDx2;

	  UDx1.conjg();
	  UDx2.conjg();
	  UDx1pUDx2 = UDx1+UDx2;

	  UDy1 = lat->cells[pos]->getUy1();
	  UDy2 = lat->cells[pos]->getUy2();
	  Uy1pUy2 = UDy1+UDy2;

	  UDy1.conjg();
	  UDy2.conjg();
	  UDy1pUDy2 = UDy1+UDy2;


	  


	  // do Ux(3) first
	  //initial guess for alpha
	  for (int ai=0; ai<Nc2m1; ai++)
	    {
	      alpha[ai] = 0.;
	    }

	  // ---- set new Ux(3) --------------------------------------------
	  
	  
	  
	  // //set exponent to zero
	  // for(int l=0; l<Nc*Nc; l++) temp2.set(l,complex<double>(0.0,0.0));
	  // //set exponent to alpha_b t^b
	  // for(int ci=0; ci<Nc2m1; ci++)
	  //   {
	  //     temp2 = temp2 + real(alpha[ci])*group->getT(ci); // alpha is real anyways but this removes numerical noise
	  //   }
	  // //set exponent to i alpha_b t^b
	  // for(int nc=0; nc<Nc*Nc; nc++)
	  //   {
	  //     temp3 = temp2.getRe(nc);
	  //     temp2.setRe(nc,-temp2.getIm(nc));
	  //     temp2.setIm(nc,temp3);
	  //   }
	  // temp2.expm();
	 
	  // cout << "OLD: " << endl << temp2 << endl;

	  // lat->cells[pos]->setUx(temp2); 

	  lat->cells[pos]->setUx(one); 
	  
	  


	  // solve for alpha iteratively (U(3)=exp(i alpha_b t^b))
	  for (int ni=0; ni<maxIterations; ni++)
	    {
	      //cout << "iteration " << ni << endl;
	      //set exponential term
	      temp2=lat->cells[pos]->getUx(); // contains exp(i alpha_b t^b)
	      expAlpha=temp2;
	      temp2.conjg();
	      expNegAlpha=temp2; // contains exp(-i alpha_b t^b)
	      
	      //cout << "Ux(3)=" << expAlpha << endl;
	      // compute Jacobian
	      countMe = 0;
	      for(int ai=0; ai<Nc2m1; ai++)
		{
		  for(int bi=0; bi<Nc2m1; bi++)
		    {
		      temp = group->getT(ai)*Ux1pUx2*group->getT(bi)*expNegAlpha+group->getT(ai)*expAlpha*group->getT(bi)*UDx1pUDx2;
		      // -i times trace of temp gives my Jacobian matrix elements:
		      M[countMe] = complex<double>(0.,-1.)*temp.trace();
			//complex<double>(temp.getIm(0)+temp.getIm(4)+temp.getIm(8),-temp.getRe(0)-temp.getRe(4)-temp.getRe(8)); 
		      countMe++;
		    }
		}

	      // compute function F that needs to be zero
	      for(int ai=0; ai<Nc2m1; ai++)
		{
		  temp = group->getT(ai)*(Ux1pUx2-UDx1pUDx2)+group->getT(ai)*Ux1pUx2*expNegAlpha-group->getT(ai)*expAlpha*UDx1pUDx2;
		  // minus trace if temp gives -F_ai
		  F[ai] = (-1.)*temp.trace();
		    //complex<double>(-temp.getRe(0)-temp.getRe(4)-temp.getRe(8),-temp.getIm(0)-temp.getIm(4)-temp.getIm(8));
		}
	      Dalpha = solveAxb(param,M,F);
	     
		
	      // for(int ai=0; ai<Nc2m1; ai++)
	      // 	{
	      // 	  cout << Dalpha[ai] << endl;
	      // 	}
	      
	      Fold = 0.;
	      lambda=1.;
	      for(int ai=0; ai<Nc2m1; ai++)
		{
		  alphaSave[ai] = alpha[ai];
		  Fold += 0.5*(real(F[ai])*real(F[ai])+imag(F[ai])*imag(F[ai]));
		}
	      
	      alphaCheck=0;
	      // solve J_{ab} \Dalpha_b = -F_a and do alpha -> alpha+Dalpha
	      // reject or accept the new alpha:
	      if (Dalpha[0].real()!=Dalpha[0].real())
		{
		  //		  cout<<"setting to one" << endl;
		  alphaCheck=1;
		  lat->cells[pos]->setUx(one); 
		}

	      while(alphaCheck==0)
		{
		  for(int ai=0; ai<Nc2m1; ai++)
		    {
		      alpha[ai] = alphaSave[ai]+lambda*Dalpha[ai];
		      //  cout << "Dalpha=" << Dalpha[ai] << endl;
		    }
		  
		  // ---- set new Ux(3) --------------------------------------------


		  for (int a=0; a<Nc2m1; a++)
		    {
		      in[a] = (alpha[a]).real(); // expmCoeff wil calculate exp(i in[a]t[a])
		    }
		  
		  U = tempNew.expmCoeff(in, Nc);
		  
		  temp2 = U[0]*one + U[1]*group->getT(0) + U[2]*group->getT(1) + U[3]*group->getT(2) + U[4]*group->getT(3) + 
		    U[5]*group->getT(4) + U[6]*group->getT(5) + U[7]*group->getT(6) + U[8]*group->getT(7);
		  
		  lat->cells[pos]->setUx(temp2); 
	      
		  //	  cout << "NEW: " << endl << tempNew << endl;

		  // //set exponent to zero
		  // for(int l=0; l<Nc*Nc; l++) temp2.set(l,complex<double>(0.0,0.0));
		  // //set exponent to alpha_b t^b
		  // for(int ci=0; ci<Nc2m1; ci++)
		  //   {
		  //     temp2 = temp2 + real(alpha[ci])*group->getT(ci); // alpha is real anyways but this removes numerical noise
		  //   }
		
		  // //set exponent to i alpha_b t^b
		  // for(int nc=0; nc<Nc*Nc; nc++)
		  //   {
		  //     temp3 = temp2.getRe(nc);
		  //     temp2.setRe(nc,-temp2.getIm(nc));
		  //     temp2.setIm(nc,temp3);
		  //   }
		  // temp2.expm();
		  //lat->cells[pos]->setUx(temp2); 

		  //	  cout << "OLD: " << endl << temp2 << endl;
		  
		  expAlpha=temp2;
		  temp2.conjg();
		  expNegAlpha=temp2; // contains exp(-i alpha_b t^b)
		  
		  for(int ai=0; ai<Nc2m1; ai++)
		    {
		      temp = group->getT(ai)*(Ux1pUx2-UDx1pUDx2)+group->getT(ai)*Ux1pUx2*expNegAlpha-group->getT(ai)*expAlpha*UDx1pUDx2;
		      // minus trace of temp gives -F_ai
		      F[ai] = (-1.)*temp.trace();
			//complex<double>(-temp.getRe(0)-temp.getRe(4)-temp.getRe(8),-temp.getIm(0)-temp.getIm(4)-temp.getIm(8));
		    }
		  
		  // ---- done: set U(3) ------------------------------------------
		  
		  // quit the misery and try a new start
		  if (lambda==0.1)
		    {
		      for (int ai=0; ai<Nc2m1; ai++)
			{
			  alpha[ai] = 0.1*random->Gauss();
			}
		 
		      for (int a=0; a<Nc2m1; a++)
			{
			  in[a] = (alpha[a]).real(); // expmCoeff wil calculate exp(i in[a]t[a])
			}
		      
		      U = tempNew.expmCoeff(in, Nc);
		      
		      temp2 = U[0]*one + U[1]*group->getT(0) + U[2]*group->getT(1) + U[3]*group->getT(2) + U[4]*group->getT(3) + 
			U[5]*group->getT(4) + U[6]*group->getT(5) + U[7]*group->getT(6) + U[8]*group->getT(7);
		      
		      lat->cells[pos]->setUx(temp2); 
	

		      // //set exponent to zeo
		      // for(int l=0; l<Nc*Nc; l++) temp2.set(l,complex<double>(0,0));
		      // //set exponent to alpha_b t^b
		      // for(int ci=0; ci<Nc2m1; ci++)
		      // 	{
		      // 	  temp2 = temp2 + real(alpha[ci])*group->getT(ci); // alpha is real anyways but this removes numerical noise
		      // 	}
		      // //set exponent to i alpha_b t^b
		      // for(int nc=0; nc<Nc*Nc; nc++)
		      // 	{
		      // 	  temp3 = temp2.getRe(nc);
		      // 	  temp2.setRe(nc,-temp2.getIm(nc));
		      // 	  temp2.setIm(nc,temp3);
		      // 	}
		      // temp2.expm();
		      // lat->cells[pos]->setUx(temp2); 
		      lambda = 1.;
		      //cout <<"break" << endl;
		      break;
		    }
		  
		  Fnew = 0.;
		  for(int ai=0; ai<Nc2m1; ai++)
		    {
		      Fnew += 0.5*(real(F[ai])*real(F[ai])+imag(F[ai])*imag(F[ai]));
		    }
		  		  
		  if(Fnew>Fold-0.00001*(Fnew*2.))
		    {
		      lambda = max(lambda*0.9,0.1);
		    }
		  else
		    {
		      alphaCheck=1;
		    }
		}
	      
	      checkConvergence=1;

	      if(Nc==2 && Fnew<0.00000001)
		checkConvergence=0;

	      if(Nc==3 && Fnew<0.0001)
		checkConvergence=0;
	      
	      if (Dalpha[0].real()!=Dalpha[0].real())
		checkConvergence=0;
	

	      if(checkConvergence==0)
		{
		  break;
		}
	      else if (ni==maxIterations-1)
		{
		  cout << i << " " << j << " result for Ux(3) did not converge for x!" << endl;
		  cout << "last Dalpha = " << endl;
		  for(int ai=0; ai<Nc2m1; ai++)
		    {
		      cout << "Dalpha/alpha=" << Dalpha[ai]/alpha[ai] << endl;
		      cout << "Dalpha=" << Dalpha[ai] << endl;
		      cout << param->getAverageQs() << " " << param->getAverageQsAvg() << " " << param->getAverageQsmin() << endl;
		      cout << param->getb() << " " << endl;
		      cout <<  lat->cells[pos]->getUx1() << " " <<  lat->cells[pos]->getUx2() << endl;
		      cout << lat->cells[pos]->getU() << endl;
		    }
		}
	    }// iteration loop
	  









	  // ------------------------------------------------------------------
	  // now do Uy(3)
	  // initial guess for alpha
 	  for (int ai=0; ai<Nc2m1; ai++)
 	    {
 	      alpha[ai] = 0.;
 	    }

	  countMe = 0;

// 	  for(int ai=0; ai<Nc2m1; ai++)
// 	    {
// 	      for(int bi=0; bi<Nc2m1; bi++)
// 		{
// 		  temp = group->getT(bi)*group->getT(ai)*Uy1pUy2+group->getT(ai)*group->getT(bi)*UDy1pUDy2;
// 		  M[countMe] = temp.trace(); 
// 		    //complex<double>(temp.getRe(0)+temp.getRe(4)+temp.getRe(8),temp.getIm(0)+temp.getIm(4)+temp.getIm(8)); 
// 		  //cout << ni << ": M[" << ai << "," << bi << "]=" << M[countMe] << endl;
// 		  countMe++;
// 		}
// 	    }

// 	  for(int ai=0; ai<Nc2m1; ai++)
// 	    {
// 	      temp = (-2.)*group->getT(ai)*(Uy1pUy2-UDy1pUDy2);
// 	      // factor of i below
// 	      F[ai] = complex<double>(0.,1.)*temp.trace();
// 		complex<double>(-temp.getIm(0)-temp.getIm(4)-temp.getIm(8),temp.getRe(0)+temp.getRe(4)+temp.getRe(8));
// 	    }
	  
// 	  alpha = solveAxb(param,M,F);

	  // ---- set initial Uy(3) --------------------------------------------

	  // //set exponent to zero
	  // for(int l=0; l<Nc*Nc; l++) temp2.set(l,complex<double>(0.0,0.0));
	  // //set exponent to alpha_b t^b
	  // for(int ci=0; ci<Nc2m1; ci++)
	  //   {
	  //     temp2 = temp2 + real(alpha[ci])*group->getT(ci);
	  //   }
	  // //set exponent to i alpha_b t^b
	  // for(int nc=0; nc<Nc*Nc; nc++)
	  //   {
	  //     temp3 = temp2.getRe(nc);
	  //     temp2.setRe(nc,-temp2.getIm(nc));
	  //     temp2.setIm(nc,temp3);
	  //   }
	  // temp2.expm();
	  
	  //  lat->cells[pos]->setUy(temp2); 
	  lat->cells[pos]->setUy(one); 

	  // ---- done: set U(3) ------------------------------------------

	  // solve for alpha iteratively (U(3)=exp(i alpha_b t^b))
	  for (int ni=0; ni<maxIterations; ni++)
	    {
	      //set exponential term
	      temp2=lat->cells[pos]->getUy(); // contains exp(i alpha_b t^b)
	      expAlpha=temp2;
	      temp2.conjg();
	      expNegAlpha=temp2; // contains exp(-i alpha_b t^b)

	      // compute Jacobian
	      countMe = 0;
	      for(int ai=0; ai<Nc2m1; ai++)
		{
		  for(int bi=0; bi<Nc2m1; bi++)
		    {
		   
		      temp = group->getT(ai)*Uy1pUy2*group->getT(bi)*expNegAlpha
			+group->getT(ai)*expAlpha*group->getT(bi)*UDy1pUDy2;
		      // -i times trace of temp gives my Jacobian matrix elements:
		      M[countMe] = complex<double>(0.,-1.)*temp.trace();
		      //complex<double>(temp.getIm(0)+temp.getIm(4)+temp.getIm(8),-temp.getRe(0)-temp.getRe(4)-temp.getRe(8)); 
		      //cout << "M[" << ai << "," << bi << "]=" << M[countMe] << endl;
		      countMe++;
		    }
		}

	      // compute function F that needs to be zero
	      for(int ai=0; ai<Nc2m1; ai++)
		{
		  temp = group->getT(ai)*(Uy1pUy2-UDy1pUDy2)
		    +group->getT(ai)*Uy1pUy2*expNegAlpha
		    -group->getT(ai)*expAlpha*UDy1pUDy2;
		  // minus trace if temp gives -F_ai
		  F[ai] = (-1.)*temp.trace();
	      //   complex<double>(-temp.getRe(0)-temp.getRe(4)-temp.getRe(8),-temp.getIm(0)-temp.getIm(4)-temp.getIm(8));
		}

	      // solve J_{ab} \Dalpha_b = -F_a and do alpha -> alpha+Dalpha
	      Dalpha = solveAxb(param,M,F);
	      
	      Fold = 0.;
	      lambda=1.;
	      for(int ai=0; ai<Nc2m1; ai++)
		{
		  alphaSave[ai] = alpha[ai];
		  Fold += 0.5*(real(F[ai])*real(F[ai])+imag(F[ai])*imag(F[ai]));
		}
	      
	      alphaCheck=0;
	      // solve J_{ab} \Dalpha_b = -F_a and do alpha -> alpha+Dalpha
	      // reject or accept the new alpha:
	      if (Dalpha[0].real()!=Dalpha[0].real())
		{
		  alphaCheck=1;
		  lat->cells[pos]->setUy(one); 
		}
	    
	      while(alphaCheck==0)
		{
		  for(int ai=0; ai<Nc2m1; ai++)
		    {
		      alpha[ai] = alphaSave[ai]+lambda*Dalpha[ai];
		      //  cout << "Dalpha=" << Dalpha[ai] << endl;
		    }
	

// 		  for(int ai=0; ai<Nc2m1; ai++)
// 		    {
// 		      alpha[ai]=alphaSave[ai];
// 		    }
// 		  for(int ai=0; ai<Nc2m1; ai++)
// 		    {
// 		      alpha[ai] = alpha[ai]+lambda*Dalpha[ai];
// 		      //  cout << "Dalpha=" << Dalpha[ai] << endl;
// 		    }
		  
		  // ---- set new Uy(3) --------------------------------------------
	

		  for (int a=0; a<Nc2m1; a++)
		    {
		      in[a] = (alpha[a]).real(); // expmCoeff wil calculate exp(i in[a]t[a])
		    }
		  
		  U = tempNew.expmCoeff(in, Nc);
		  
		  temp2 = U[0]*one + U[1]*group->getT(0) + U[2]*group->getT(1) + U[3]*group->getT(2) + U[4]*group->getT(3) + 
		    U[5]*group->getT(4) + U[6]*group->getT(5) + U[7]*group->getT(6) + U[8]*group->getT(7);
		  
		  lat->cells[pos]->setUy(temp2); 


		  // //set exponent to zero
		  // for(int l=0; l<Nc*Nc; l++) temp2.set(l,complex<double>(0.0,0.0));
		  // //set exponent to alpha_b t^b
		  // for(int ci=0; ci<Nc2m1; ci++)
		  //   {
		  //     temp2 = temp2 + real(alpha[ci])*group->getT(ci); // alpha is real anyways but this removes numerical noise
		  //   }
		  // //set exponent to i alpha_b t^b
		  // for(int nc=0; nc<Nc*Nc; nc++)
		  //   {
		  //     temp3 = temp2.getRe(nc);
		  //     temp2.setRe(nc,-temp2.getIm(nc));
		  //     temp2.setIm(nc,temp3);
		  //   }
		  // temp2.expm();
		  // lat->cells[pos]->setUy(temp2); 
		  
		  expAlpha=temp2;
		  temp2.conjg();
		  expNegAlpha=temp2; // contains exp(-i alpha_b t^b)
		  
		  for(int ai=0; ai<Nc2m1; ai++)
		    {
		      temp = group->getT(ai)*(Uy1pUy2-UDy1pUDy2)
			+group->getT(ai)*Uy1pUy2*expNegAlpha
			-group->getT(ai)*expAlpha*UDy1pUDy2;
		      // minus trace if temp gives -F_ai
		      F[ai] = (-1.)*temp.trace();
			//complex<double>(-temp.getRe(0)-temp.getRe(4)-temp.getRe(8),-temp.getIm(0)-temp.getIm(4)-temp.getIm(8));
		    }
		  
		  //cout << ni << " Ux=" << lat->cells[pos]->getUx() << endl; 
		  //cout << "UxUDx=" << lat->cells[pos]->getUx()*temp2 << endl; 
		  // ---- done: set U(3) ------------------------------------------
		  
		 
		  if (lambda==0.1)
		    {
		      for (int ai=0; ai<Nc2m1; ai++)
			{
			  alpha[ai] = 0.1*random->Gauss();
			}
		  
		      for (int a=0; a<Nc2m1; a++)
			{
			  in[a] = (alpha[a]).real(); // expmCoeff will calculate exp(i in[a]t[a])
			}
		      
		      U = tempNew.expmCoeff(in, Nc);
		      
		      temp2 = U[0]*one + U[1]*group->getT(0) + U[2]*group->getT(1) + U[3]*group->getT(2) + U[4]*group->getT(3) + 
			U[5]*group->getT(4) + U[6]*group->getT(5) + U[7]*group->getT(6) + U[8]*group->getT(7);
		      
		      lat->cells[pos]->setUy(temp2); 


		      // //set exponent to noise
		      // for(int l=0; l<Nc*Nc; l++) temp2.set(l,complex<double>(0,0));
		      // //set exponent to alpha_b t^b
		      // for(int ci=0; ci<Nc2m1; ci++)
		      // 	{
		      // 	  temp2 = temp2 + real(alpha[ci])*group->getT(ci); // alpha is real anyways but this removes numerical noise
		      // 	}
		      // //set exponent to i alpha_b t^b
		      // for(int nc=0; nc<Nc*Nc; nc++)
		      // 	{
		      // 	  temp3 = temp2.getRe(nc);
		      // 	  temp2.setRe(nc,-temp2.getIm(nc));
		      // 	  temp2.setIm(nc,temp3);
		      // 	}
		      // temp2.expm();
		      // lat->cells[pos]->setUy(temp2); 
		      lambda = 1.;
		      //cout <<"break" << endl;
		      break;
		    }
		 

 	  
		  Fnew = 0.;
		  for(int ai=0; ai<Nc2m1; ai++)
		    {
		      Fnew += 0.5*(real(F[ai])*real(F[ai])+imag(F[ai])*imag(F[ai]));
		    }
		  
		  if(Fnew>Fold-0.00001*(Fnew*2.))
		    {
		      lambda = max(lambda*0.9,0.1);
		    }
		  else
		    alphaCheck=1;
		}

	      checkConvergence=1;
	      if(Nc==2 && Fnew<0.00000001)
		checkConvergence-=1;
	      if(Nc==3 && Fnew<0.0001)
		checkConvergence-=1;
	      if (Dalpha[0].real()!=Dalpha[0].real())
		checkConvergence=0;
		     
	      if(checkConvergence==0)
		{
		  //cout << "result converged after " << ni << " iterations" << endl;
		  break;
		}
	      else if (ni==maxIterations-1)
		{
		  cout << i << " " << j << " result for Uy(3) did not converge for y!" << endl;
		  cout << "last Dalpha = " << endl;
		  for(int ai=0; ai<Nc2m1; ai++)
		    {
		      cout << "Dalpha/alpha=" << Dalpha[ai]/alpha[ai] << endl;
		      cout << "Dalpha=" << Dalpha[ai] << endl;
		      cout << param->getAverageQs() << " " << param->getAverageQsAvg() << " " << param->getAverageQsmin() << endl;

		    }
		}
	    }//iteration loop

	

	}
    }//loops over x and y
  
  // ofstream foutU5("Ux.txt",ios::out); 
  // foutU5.precision(15);
  // for(int ix=0; ix<N; ix++)
  //   {
  //     for(int iy=0; iy<N; iy++) // loop over all positions
  // 	{
  // 	  pos = ix*N+iy;
  // 	  foutU5 << ix << " " << iy << " "  << (lat->cells[pos]->getUx()).MatrixToString() << endl;
  // 	}
  //   }
  // foutU5.close();
  
  // ofstream foutU6("Uy.txt",ios::out); 
  // foutU6.precision(15);
  // for(int ix=0; ix<N; ix++)
  //   {
  //     for(int iy=0; iy<N; iy++) // loop over all positions
  // 	{
  // 	  pos = ix*N+iy;
  // 	  foutU6 << ix << " " << iy << " "  << (lat->cells[pos]->getUy()).MatrixToString() << endl;
  // 	}
  //   }
  // foutU6.close();
  

  // ------


 //  stringstream strtU1_name;
 //  strtU1_name << "VOne-" << param->getMPIRank() << ".dat";
 //  string tU1_name;
 //  tU1_name = strtU1_name.str();

 //  ofstream foutU1(tU1_name.c_str(),ios::out); 

 // //  // output U
 //  ofstream foutU("UOne.txt",ios::out); 
 //  for(int ix=0; ix<N; ix++)
 //    {
 //      for(int iy=0; iy<N; iy++) // loop over all positions
 // 	{
 // 	  pos = ix*N+iy;
 // 	  foutU << ix << " " << iy << " " <<  setprecision(10) 
 // 		<< (complex<double>(0.,1.)*(lat->cells[pos]->getUx1()*group->getT(0)).trace()).real() << " "
 // 		<< (complex<double>(0.,1.)*(lat->cells[pos]->getUx1()*group->getT(1)).trace()).real() << " " 
 // 		<< (complex<double>(0.,1.) * (lat->cells[pos]->getUx1()*group->getT(2)).trace()).real() << " " 
 // 		<< (lat->cells[pos]->getUx1().trace()).real() << " " 
 // 		<< (complex<double>(0.,1.)*(lat->cells[pos]->getUy1()*group->getT(0)).trace()).real() << " "
 // 		<< (complex<double>(0.,1.)*(lat->cells[pos]->getUy1()*group->getT(1)).trace()).real() << " " 
 // 		<< (complex<double>(0.,1.) * (lat->cells[pos]->getUy1()*group->getT(2)).trace()).real() << " " 
 // 		<< (lat->cells[pos]->getUy1().trace()).real() << endl;
 // 	}
 //    }
 //  foutU.close();
 //  ofstream foutU2("UTwo.txt",ios::out); 
 //  for(int ix=0; ix<N; ix++)
 //    {
 //      for(int iy=0; iy<N; iy++) // loop over all positions
 // 	{
 // 	  pos = ix*N+iy;
 // 	  foutU2 << ix << " " << iy << " " <<  setprecision(10) 
 // 		 << (complex<double>(0.,1.)*(lat->cells[pos]->getUx2()*group->getT(0)).trace()).real() << " "
 // 		 << (complex<double>(0.,1.)*(lat->cells[pos]->getUx2()*group->getT(1)).trace()).real() << " " 
 // 		 << (complex<double>(0.,1.) * (lat->cells[pos]->getUx2()*group->getT(2)).trace()).real() << " " 
 // 		 << (lat->cells[pos]->getUx2().trace()).real() << " " 
 // 		 << (complex<double>(0.,1.)*(lat->cells[pos]->getUy2()*group->getT(0)).trace()).real() << " "
 // 		 << (complex<double>(0.,1.)*(lat->cells[pos]->getUy2()*group->getT(1)).trace()).real() << " " 
 // 		 << (complex<double>(0.,1.) * (lat->cells[pos]->getUy2()*group->getT(2)).trace()).real() << " " 
 // 		 << (lat->cells[pos]->getUy2().trace()).real() << endl;
 // 	}
 //    }
 //  foutU2.close();




 //  // ------


	  
  double UDU[bins];
  
  for(int i=0; i<bins; i++)
    {
      UDU[i] = 0.;
      count[i]=0;
    }

  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*N+j;
	  UD=lat->cells[pos]->getU();
	  UD.conjg();
	  Ux = lat->cells[N/2*N+N/2]->getU()*UD;
	  x = -L/2.+a*i;
	  y = -L/2.+a*j;
	  r = sqrt(x*x+y*y);
	  ir = static_cast<int>(floor(r/dr+0.000001));
	  count[ir]++;
	  UDU[ir] += (Ux.trace()).real();//Ux.getRe(0)+Ux.getRe(4)+Ux.getRe(8);
     	}
    }

  // // output UDU
  // if( param->getWriteOutputs() == 1)
  //   {
  //     cout << "output UDU" << endl;

  //     stringstream strUDU_name;
  //     strUDU_name << "UDU" << param->getMPIRank() << ".dat";
  //     string UDU_name;
  //     UDU_name = strUDU_name.str();
  //     ofstream foutrUDU(UDU_name.c_str(),ios::out); 
  //     for(int i=0; i<N; i++) // loop over all positions
  // 	{
  // 	  foutrUDU << i*dr << " " << 1./param->getNc()*UDU[i]/count[i] << endl;
  // 	}
  //     foutrUDU.close();
  //   }
  // cout << "done." << endl;


  // compute initial electric field
  // with minus ax, ay
  //#pragma omp parallel for collapse(2)
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*N+j;
	  if(i>0)
	    posxm = (i-1)*N+j;
	  else
	    posxm = (N-1)*N+j;

	  if(j>0)
	    posym = i*N+(j-1);
	  else
	    posym = i*N+(N-1);

	  // x part in sum:
	  Ux1mUx2 = lat->cells[pos]->getUx1()-lat->cells[pos]->getUx2();
	  UDx1 = lat->cells[pos]->getUx1();
	  UDx1.conjg();
	  UDx2 = lat->cells[pos]->getUx2();
	  UDx2.conjg();
	  UDx1mUDx2 = UDx1 - UDx2;
	  
	  Ux = lat->cells[pos]->getUx();
	  UDx = Ux;
	  UDx.conjg();

	  temp2 = Ux1mUx2*UDx - Ux1mUx2 - Ux*UDx1mUDx2 + UDx1mUDx2;
	  
	  Ux1mUx2 = lat->cells[posxm]->getUx1()-lat->cells[posxm]->getUx2();
	  UDx1 = lat->cells[posxm]->getUx1();
	  UDx1.conjg();
	  UDx2 = lat->cells[posxm]->getUx2();
	  UDx2.conjg();
	  UDx1mUDx2 = UDx1 - UDx2;
	  
	  Ux = lat->cells[posxm]->getUx();
	  UDx = Ux;
	  UDx.conjg();
	  
	  temp2 = temp2 - UDx*Ux1mUx2 + Ux1mUx2 + UDx1mUDx2*Ux - UDx1mUDx2; 
	  
	  // y part in sum
	  Uy1mUy2 = lat->cells[pos]->getUy1()-lat->cells[pos]->getUy2();
	  UDy1 = lat->cells[pos]->getUy1();
	  UDy1.conjg();
	  UDy2 = lat->cells[pos]->getUy2();
	  UDy2.conjg();
	  UDy1mUDy2 = UDy1 - UDy2;
	  
	  Uy = lat->cells[pos]->getUy();
	  UDy = Uy;
	  UDy.conjg();

	  // y part of the sum:
	  temp2 = temp2 + Uy1mUy2*UDy - Uy1mUy2 - Uy*UDy1mUDy2 + UDy1mUDy2;

	  Uy1mUy2 = lat->cells[posym]->getUy1()-lat->cells[posym]->getUy2();
	  UDy1 = lat->cells[posym]->getUy1();
	  UDy1.conjg();
	  UDy2 = lat->cells[posym]->getUy2();
	  UDy2.conjg();
	  UDy1mUDy2 = UDy1 - UDy2;
	  
	  Uy = lat->cells[posym]->getUy();
	  UDy = Uy;
	  UDy.conjg();
	  
	  temp2 = temp2 - UDy*Uy1mUy2 + Uy1mUy2 + UDy1mUDy2*Uy - UDy1mUDy2; 
	  
	  lat->cells[pos]->setE1((1./8.)*temp2);

	  //	  if(i==nn[0]-1)
	  //  lat->cells[pos]->setE1(zero);
	    
	}
    }

  // with plus ax, ay
  //#pragma omp parallel for collapse(2)
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*N+j;
	  if(i<N-1)
	    posx = (i+1)*N+j;
	  else 
	    posx = (N-1)*N+j;

	  if(j<N-1)
	    posy = i*N+(j+1);
	  else
	    posy = i*N+N-1;

	  // x part in sum:
	  Ux1mUx2 = lat->cells[pos]->getUx1()-lat->cells[pos]->getUx2();
	  UDx1 = lat->cells[pos]->getUx1();
	  UDx1.conjg();
	  UDx2 = lat->cells[pos]->getUx2();
	  UDx2.conjg();
	  UDx1mUDx2 = UDx1 - UDx2;
	  
	  Ux = lat->cells[pos]->getUx();
	  UDx = Ux;
	  UDx.conjg();

	  temp2 = Ux1mUx2*UDx - Ux1mUx2 - Ux*UDx1mUDx2 + UDx1mUDx2;
	  
	  Ux1mUx2 = lat->cells[posx]->getUx1()-lat->cells[posx]->getUx2();
	  UDx1 = lat->cells[posx]->getUx1();
	  UDx1.conjg();
	  UDx2 = lat->cells[posx]->getUx2();
	  UDx2.conjg();
	  UDx1mUDx2 = UDx1 - UDx2;
	  
	  Ux = lat->cells[posx]->getUx();
	  UDx = Ux;
	  UDx.conjg();
	  
	  temp2 = temp2 - UDx*Ux1mUx2 + Ux1mUx2 + UDx1mUDx2*Ux - UDx1mUDx2; 
	  
	  // y part in sum
	  Uy1mUy2 = lat->cells[pos]->getUy1()-lat->cells[pos]->getUy2();
	  UDy1 = lat->cells[pos]->getUy1();
	  UDy1.conjg();
	  UDy2 = lat->cells[pos]->getUy2();
	  UDy2.conjg();
	  UDy1mUDy2 = UDy1 - UDy2;
	  
	  Uy = lat->cells[pos]->getUy();
	  UDy = Uy;
	  UDy.conjg();

	  // y part of the sum:
	  temp2 = temp2 + Uy1mUy2*UDy - Uy1mUy2 - Uy*UDy1mUDy2 + UDy1mUDy2;

	  Uy1mUy2 = lat->cells[posy]->getUy1()-lat->cells[posy]->getUy2();
	  UDy1 = lat->cells[posy]->getUy1();
	  UDy1.conjg();
	  UDy2 = lat->cells[posy]->getUy2();
	  UDy2.conjg();
	  UDy1mUDy2 = UDy1 - UDy2;
	  
	  Uy = lat->cells[posy]->getUy();
	  UDy = Uy;
	  UDy.conjg();
	  
	  temp2 = temp2 - UDy*Uy1mUy2 + Uy1mUy2 + UDy1mUDy2*Uy - UDy1mUDy2; 
	  
	  lat->cells[pos]->setE2((1./8.)*temp2);
	  
	}
    }

  // compute the plaquette
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*N+j;
	  posx = ((i+1)%N)*N+j;
	  posy = i*N+(j+1)%N;

	  if(i<N-1 && j<N-1)
	    {
	      UDx = lat->cells[posy]->getUx();
	      UDy = lat->cells[pos]->getUy();
	      UDx.conjg();
	      UDy.conjg();

	      Uplaq = lat->cells[pos]->getUx()*(lat->cells[posx]->getUy()*(UDx*UDy));
	    }
	  else
	    {
	      UDx = lat->cells[pos]->getUx();
	      UDx.conjg();
	      Uplaq = lat->cells[pos]->getUx()*UDx; //indentity
	    }

	  lat->cells[pos]->setUplaq(Uplaq);
	}
    }

    
    
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
       	  pos = i*N+j;
	  
	  //	  AM = (lat->cells[pos]->getAetaM());//+lat->cells[pos]->getAetaP());
	  //      AP = (lat->cells[pos]->getAetaP());//+lat->cells[pos]->getAetaP());

	  AM = (lat->cells[pos]->getE1());//+lat->cells[pos]->getAetaP());
	  AP = (lat->cells[pos]->getE2());//+lat->cells[pos]->getAetaP());

	  // this is pi in lattice units as needed for the evolution. (later, the a^4 gives the right units for the energy density
	  lat->cells[pos]->setpi(complex<double>(0.,-2./param->getg())*(AM)); // factor -2 because I have A^eta (note the 1/8 before) but want \pi (E^z). Soren sagt, da muss ein minus hin... check .
	  // lat->cells[pos]->setpi(complex<double>(0.,-1./param->getg())*(AM+AP)); // factor -2 because I have A^eta (note the 1/8 before) but want \pi (E^z). Soren sagt, da muss ein minus hin... check .
	  lat->cells[pos]->setE1(zero);
	  lat->cells[pos]->setE2(zero);
	  lat->cells[pos]->setphi(zero);
	  lat->cells[pos]->setUx1(one); // reset the Ux1 to be used for other purposes later
	}
    }

  

 // // output U
 //  ofstream foutU3("U.txt",ios::out); 
 //  for(int ix=0; ix<N; ix++)
 //    {
 //      for(int iy=0; iy<N; iy++) // loop over all positions
 // 	{
 // 	  pos = ix*N+iy;
 // 	  foutU3 << ix << " " << iy << " " << (lat->cells[pos]->getUx()).MatrixToString() << " " 
 // 		 << (lat->cells[pos]->getUy()).MatrixToString() << endl;
 
 // 		// << (complex<double>(0.,1.)*(lat->cells[pos]->getUx()*group->getT(0)).trace()).real() << " "
 // 		// << (complex<double>(0.,1.)*(lat->cells[pos]->getUx()*group->getT(1)).trace()).real() << " " 
 // 		// << (complex<double>(0.,1.)*(lat->cells[pos]->getUx()*group->getT(2)).trace()).real() << " " 
 // 		// << (complex<double>(0.,1.)*(lat->cells[pos]->getUx()*group->getT(3)).trace()).real() << " " 
 // 		// << (complex<double>(0.,1.)*(lat->cells[pos]->getUx()*group->getT(4)).trace()).real() << " " 
 // 		// << (complex<double>(0.,1.)*(lat->cells[pos]->getUx()*group->getT(5)).trace()).real() << " " 
 // 		// << (complex<double>(0.,1.)*(lat->cells[pos]->getUx()*group->getT(6)).trace()).real() << " " 
 // 		// << (complex<double>(0.,1.)*(lat->cells[pos]->getUx()*group->getT(7)).trace()).real() << " " 
 // 		// << (lat->cells[pos]->getUx().trace()).real() << " " 
 // 		// << (complex<double>(0.,1.)*(lat->cells[pos]->getUy()*group->getT(0)).trace()).real() << " "
 // 		// << (complex<double>(0.,1.)*(lat->cells[pos]->getUy()*group->getT(1)).trace()).real() << " " 
 // 		// << (complex<double>(0.,1.)*(lat->cells[pos]->getUy()*group->getT(2)).trace()).real() << " " 
 // 		// << (lat->cells[pos]->getUy().trace()).real() << endl;
 // 	}
 //    }
 //  foutU3.close();

 //  ofstream foutU4("pi.txt",ios::out); 
 //  for(int ix=0; ix<N; ix++)
 //    {
 //      for(int iy=0; iy<N; iy++) // loop over all positions
 // 	{
 // 	  pos = ix*N+iy;
 // 	  foutU4 << ix << " " << iy << " " 
 // 		 << 2.*((lat->cells[pos]->getpi()*group->getT(0)).trace()).real() << " "
 // 		 << 2.*((lat->cells[pos]->getpi()*group->getT(1)).trace()).real() << " " 
 // 		 << 2.*((lat->cells[pos]->getpi()*group->getT(2)).trace()).real() << " " 
 // 		 << 2.*((lat->cells[pos]->getpi()*group->getT(3)).trace()).real() << " " 
 // 		 << 2.*((lat->cells[pos]->getpi()*group->getT(4)).trace()).real() << " " 
 // 		 << 2.*((lat->cells[pos]->getpi()*group->getT(5)).trace()).real() << " " 
 // 		 << 2.*((lat->cells[pos]->getpi()*group->getT(6)).trace()).real() << " " 
 // 		 << 2.*((lat->cells[pos]->getpi()*group->getT(7)).trace()).real() << endl;
 // 	}
 //    }
 //  foutU4.close();


  // // output E
  // ofstream foutE("E.txt",ios::out); 
  // for(int ix=0; ix<N; ix++)
  //   {
  //     for(int iy=0; iy<N; iy++) // loop over all positions
  // 	{
  // 	  pos = ix*N+iy;
  // 	  foutE << ix << " " << iy << " " <<  setprecision(10) 
  // 		<< 2.*((lat->cells[pos]->getpi()*group->getT(0)).trace()).real() << " "
  // 		<< 2.*((lat->cells[pos]->getpi()*group->getT(1)).trace()).real() << " " 
  // 		<< 2.*((lat->cells[pos]->getpi()*group->getT(2)).trace()).real() << endl;
  // 	}
  //   }
  // foutE.close();


  delete M;
  delete F;
  delete result;
  delete alpha;
  delete alphaSave;
  //  delete Dalpha;
  
  // done. 
  
  // -----------------------------------------------------------------------------
  // finish
  // -----------------------------------------------------------------------------
}

void Init::multiplicity(Lattice *lat, Group *group, Parameters *param, Random *random, Glauber *glauber)
{
  int N = param->getSize();
  int pos;
  double epsilonSum=0.;
  double L = param->getL();
  double a = L/N; // lattice spacing in fm
     
  for(int ix=0; ix<N; ix++) 
    {
      for(int iy=0; iy<N; iy++)
	{
	  pos = ix*N+iy;
	  epsilonSum += a*a*lat->cells[pos]->getEpsilon()*hbarc;
	}
    }
  stringstream strtE_name;
  strtE_name << "totalEnergy" << param->getMPIRank() << ".dat";
  string tE_name;
  tE_name = strtE_name.str();

  ofstream fout(tE_name.c_str(),ios::out); 
  fout << epsilonSum << endl;
  fout.close();      
}


// void Init::eccentricity(Lattice *lat, Group *group, Parameters *param, Random *random, Glauber *glauber)
// {
//   int N = param->getSize();
//   int pos, posNew;
//   double rA, phiA, x, y;
//   double L = param->getL();
//   double a = L/N; // lattice spacing in fm
//   double eccentricity1, eccentricity2, eccentricity3, eccentricity4, eccentricity5, eccentricity6;
//   double avcos, avsin, avcos1, avsin1, avcos3, avsin3, avrSq, avr1, avr3, avcos4, avsin4, avr4, avcos5, avsin5, avr5, avcos6, avsin6, avr6;
//   double Pi;
//   Pi = param->getPi();
//   double Psi1, Psi2, Psi3, Psi4, Psi5, Psi6;


//   avrSq=0.;
//   avr3=0.;

//   double avx = 0.;
//   double avy = 0.;
//   double toteps = 0.;
//   int xshift;
//   int yshift;
  
//   // first shift to the center
//   for(int ix=0; ix<N; ix++) 
//     {
//       x = -L/2.+a*ix;
//       for(int iy=0; iy<N; iy++)
// 	{
// 	  pos = ix*N+iy;
// 	  y = -L/2.+a*iy;
// 	  avx += x*lat->cells[pos]->getEpsilon();
// 	  avy += y*lat->cells[pos]->getEpsilon();
// 	  toteps += lat->cells[pos]->getEpsilon();
// 	}
//     }
      
//   avx/=toteps;
//   avy/=toteps;
  
//   cout << "avx=" << avx << endl;
//   cout << "avy=" << avy << endl;
  
//   xshift = static_cast<int>(floor(avx/a+0.00000000001));
//   yshift = static_cast<int>(floor(avy/a+0.00000000001));
  
//   cout << "xshift=" << xshift << endl;
//   cout << "yshift=" << yshift << endl;
  
//   cout << "a xshift=" << a*xshift << endl;
//   cout << "a yshift=" << a*yshift << endl;
  
  
//   avcos1 = 0.;
//   avsin1 = 0.;
//   avcos = 0.;
//   avsin = 0.;
//   avcos3 = 0.;
//   avsin3 = 0.;
//   avcos4 = 0.;
//   avsin4 = 0.;
//   avcos5 = 0.;
//   avsin5 = 0.;
//   avcos6 = 0.;
//   avsin6 = 0.;
//   avr1=0.;
//   avrSq=0.;
//   avr3=0.;
//   avr4=0.;
//   avr5=0.;
//   avr6=0.;

   
//   for(int ix=2; ix<N-2; ix++) 
// 	{
// 	  x = -L/2.+a*ix-avx;
// 	  for(int iy=2; iy<N-2; iy++)
// 	    {
// 	      pos = ix*N+iy;
// 	      y = -L/2.+a*iy-avy;
// 	      if (x>=0)
// 		{
// 		  phiA = atan(y/x);
// 		  if (x==0) 
// 		    {
// 		      if (y>=0) phiA=PI/2.;
// 		      else if (y<0) phiA=3.*PI/2.;
// 		    }
// 		}
// 	      else
// 		{
// 		  phiA = atan(y/x)+PI;
// 		}

// 	      if (lat->cells[pos]->getEpsilon()<10.) // in fm^{-4}
// 		lat->cells[pos]->setEpsilon(0.);

// 	      rA = sqrt( x*x + y*y );
// 	      avr1 += rA*(lat->cells[pos]->getEpsilon());
// 	      avrSq += rA*rA*(lat->cells[pos]->getEpsilon()); // compute average r^2
// 	      avr3 += rA*rA*rA*(lat->cells[pos]->getEpsilon());
// 	      avr4 += rA*rA*rA*rA*(lat->cells[pos]->getEpsilon());
// 	      avr5 += rA*rA*rA*rA*rA*(lat->cells[pos]->getEpsilon());
// 	      avr6 += rA*rA*rA*rA*rA*rA*(lat->cells[pos]->getEpsilon());
	      
// 	      avcos1 += rA*cos(1.*phiA)*(lat->cells[pos]->getEpsilon());
// 	      avsin1 += rA*sin(1.*phiA)*(lat->cells[pos]->getEpsilon());
// 	      avcos  += rA*rA*cos(2.*phiA)*(lat->cells[pos]->getEpsilon());
// 	      avsin  += rA*rA*sin(2.*phiA)*(lat->cells[pos]->getEpsilon());
// 	      avcos3 += rA*rA*rA*cos(3.*phiA)*(lat->cells[pos]->getEpsilon());
// 	      avsin3 += rA*rA*rA*sin(3.*phiA)*(lat->cells[pos]->getEpsilon());
// 	      avcos4 += rA*rA*rA*rA*cos(4.*phiA)*(lat->cells[pos]->getEpsilon());
// 	      avsin4 += rA*rA*rA*rA*sin(4.*phiA)*(lat->cells[pos]->getEpsilon());
// 	      avcos5 += rA*rA*rA*rA*rA*cos(5.*phiA)*(lat->cells[pos]->getEpsilon());
// 	      avsin5 += rA*rA*rA*rA*rA*sin(5.*phiA)*(lat->cells[pos]->getEpsilon());
// 	      avcos6 += rA*rA*rA*rA*rA*rA*cos(6.*phiA)*(lat->cells[pos]->getEpsilon());
// 	      avsin6 += rA*rA*rA*rA*rA*rA*sin(6.*phiA)*(lat->cells[pos]->getEpsilon());
// 	    }
// 	}
      
//       // compute and print eccentricity and angles:
//       Psi1 = (atan(avsin1/avcos1)+PI)/1.;
//       Psi2 = (atan(avsin/avcos))/2.+PI/2.;
//       Psi3 = (atan(avsin3/avcos3)+PI)/3.;
//       Psi4 = (atan(avsin4/avcos4)+PI)/4.;
//       Psi5 = (atan(avsin5/avcos5)+PI)/5.;
//       Psi6 = (atan(avsin6/avcos6)+PI)/6.;
//       eccentricity1 = sqrt(avcos1*avcos1+avsin1*avsin1)/avr1;
//       eccentricity2 = sqrt(avcos*avcos+avsin*avsin)/avrSq;
//       eccentricity3 = sqrt(avcos3*avcos3+avsin3*avsin3)/avr3;
//       eccentricity4 = sqrt(avcos4*avcos4+avsin4*avsin4)/avr4;
//       eccentricity5 = sqrt(avcos5*avcos5+avsin5*avsin5)/avr5;
//       eccentricity6 = sqrt(avcos6*avcos6+avsin6*avsin6)/avr6;
      
//       cout << endl;
//       cout << "ecc1=" << eccentricity1 << endl;
//       cout << "ecc2=" << eccentricity2 << endl;
//       cout << "ecc3=" << eccentricity3 << endl;
//       cout << "ecc4=" << eccentricity4 << endl;
//       cout << "ecc5=" << eccentricity5 << endl;
//       cout << "ecc6=" << eccentricity6 << endl;
//       cout << endl;
  
//       double avx2 = avx;
//       double avy2 = avy;
//       avx=0.;
//       avy=0.;
//       toteps=0.;

//       for(int ix=0; ix<N; ix++) 
// 	{
// 	  x = -L/2.+a*ix-avx2;
// 	  for(int iy=0; iy<N; iy++)
// 	    {
// 	      pos = ix*N+iy;
// 	      y = -L/2.+a*iy-avy2;
// 	      avx += x*lat->cells[pos]->getEpsilon();
// 	      avy += y*lat->cells[pos]->getEpsilon();
// 	      toteps += lat->cells[pos]->getEpsilon();
// 	    }
// 	}
 
//       avx/=toteps;
//       avy/=toteps;
      
//       cout << "new avx=" << avx << endl;
//       cout << "new avy=" << avy << endl;
      
//       ofstream foutEcc("eccentricities.dat",ios::app); 
//       foutEcc << "0 " << eccentricity1 << " " << Psi1 << " " << eccentricity2 << " " << Psi2 << " " 
// 	      << eccentricity3 << " " << Psi3 << " " << eccentricity4 << " " << Psi4 
// 	      << " " << eccentricity5 << " " << Psi5 << " " << eccentricity6 << " " << Psi6 << endl;
//       foutEcc.close();          
// }

  


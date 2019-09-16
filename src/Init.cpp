// Init.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Init.h"
#include<algorithm>
#include <utility>
#include "Phys_consts.h"

using PhysConst::hbarc;

//**************************************************************************
// Init class.

vector <complex<double> > Init::solveAxb(Parameters *param, complex<double>* A, complex<double>* b_in)
{
  int Nc = param->getNc();
  int Nc2m1 = Nc*Nc-1;

  vector <complex<double> > xvec;
  xvec.reserve(Nc2m1); 
  double a_data[128];

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

  int s;
  gsl_permutation * p = gsl_permutation_alloc (Nc2m1);
  gsl_linalg_complex_LU_decomp (&m.matrix, p, &s);
  gsl_linalg_complex_LU_solve (&m.matrix, p, &c.vector, x);
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
}


void Init::sampleTA(Parameters *param, Random* random, Glauber* glauber)
{ 
  ReturnValue rv, rv2;
  cout << "Sampling nucleon positions ... ";

  if(param->getNucleonPositionsFromFile()==0)
    {
      int A1,A2,Z1,Z2;
      A1 = static_cast<int>(glauber->nucleusA1())*param->getAverageOverNuclei(); // projectile
      A2 = static_cast<int>(glauber->nucleusA2())*param->getAverageOverNuclei(); // target
      Z1 = static_cast<int>(glauber->nucleusZ1())*param->getAverageOverNuclei(); // projectile
      Z2 = static_cast<int>(glauber->nucleusZ2())*param->getAverageOverNuclei(); // target

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
          rv.proton=1;
	  nucleusA.push_back(rv);  
	}   
      else if(A1==2) // deuteron
	{
	  rv = glauber->SampleTARejection(random,1);
	  param->setRnp(sqrt(rv.x*rv.x+rv.y*rv.y));
	  // we sample the neutron proton distance, so distance to the center needs to be divided by 2
	  rv.x = rv.x/2.;
	  rv.y = rv.y/2.;
          rv.proton = 1;
	  rv.collided=0;
	  nucleusA.push_back(rv);
	  // other nucleon is 180 degrees rotated:	  
	  rv.x = -rv.x;
	  rv.y = -rv.y;
	  rv.collided=0;
	  nucleusA.push_back(rv);

	}   
      else if(A1==3) // He3
	{
          int rank = param->getMPIRank();
          int size;
          MPI_Comm_size (MPI_COMM_WORLD, &size);
          double package[6]; 
 
          if (rank==0)
            {
              //sample the position in the file
              ifstream fin;
              
              // stringstream strhe3_name;
              // strhe3_name << "he3_plaintext-" << param->getMPIRank()%10 << ".dat";
              // string he3_name;
              // he3_name = strhe3_name.str();
              // cout << "reading from file " << he3_name << "." << endl;
              // fin.open(he3_name); 
              
              fin.open("he3_plaintext.dat"); 
              
              double dummy;
              
              for (int event=0;event<size;event++)
                {
                  double ran2 = random->genrand64_real3();   // sample the position in the file uniformly (13699 events in file)
                  int nucleusNumber = static_cast<int>(ran2*13699);
                  
                  cout << "using nucleus Number = " << nucleusNumber << endl;
                  
                  // go to the correct line in the file
                  if(fin)
                    {
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
                              if (A==2) 
                                rv.proton=0;
                              else 
                                rv.proton=1;
                              if (event==0)
                                nucleusA.push_back(rv);
                              else 
                                {
                                  if(A==0)
                                    {
                                      package[0] = rv.x;
                                      package[1] = rv.y;
                                    }
                                  else if (A==1)
                                    {
                                      package[2] = rv.x;
                                      package[3] = rv.y;
                                    }
                                  else if (A==2)
                                    {
                                      package[4] = rv.x;
                                      package[5] = rv.y;
                                    }
                                }
                              A++;
                              cout << "A=" << A << ", x=" << rv.x << ", y=" << rv.y << endl;
                            }
                        }
                      if(event>0)
                        MPI_Send(package,6,MPI_DOUBLE, event, 1, MPI_COMM_WORLD);
                      
                      param->setA1FromFile(A);
                    }
                  else
                    {
                      cerr << " file he3_plaintext.dat not found. exiting." << endl;
                      exit(1);
                    }
                }
              fin.close();
                     
            }
          else
            {
              MPI_Recv(package,6,MPI_DOUBLE,0,1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              rv.x = package[0];
              rv.y = package[1];
              rv.proton = 1;
              nucleusA.push_back(rv);
              cout << "A=1, x=" << rv.x << ", y=" << rv.y << endl;
              rv.x = package[2];
              rv.y = package[3];
              rv.proton = 0;
              nucleusA.push_back(rv);
              cout << "A=2, x=" << rv.x << ", y=" << rv.y << endl;
              rv.x = package[4];
              rv.y = package[5];
              rv.proton = 0;
              nucleusA.push_back(rv);
              cout << "A=3, x=" << rv.x << ", y=" << rv.y << endl;
              param->setA1FromFile(3);
            }
        }
      else
	{
          generate_nucleus_configuration(random, A1, Z1,
                                         glauber->GlauberData.Projectile.a_WS,
                                         glauber->GlauberData.Projectile.R_WS,
                                         glauber->GlauberData.Projectile.beta2,
                                         glauber->GlauberData.Projectile.beta4,
                                         &nucleusA);
	}
    
      if(A2==1)
	{
	  rv2.x=0.;
	  rv2.y=0;
	  rv2.collided=0;
          rv2.proton=1;
	  nucleusB.push_back(rv2);
	}   
      else if(A2==2) // deuteron
	{
	  rv = glauber->SampleTARejection(random,2);
	  // we sample the neutron proton distance, so distance to the center needs to be divided by 2
	  param->setRnp(sqrt(rv.x*rv.x+rv.y*rv.y));

	  rv.x = rv.x/2.;
	  rv.y = rv.y/2.;
          rv2.proton=1;
	  rv.collided=0;
	  nucleusB.push_back(rv);
	 
	  // other nucleon is 180 degrees rotated:	  
	  rv.x = -rv.x;
	  rv.y = -rv.y;
          rv2.proton=0;
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
         if(fin)
            {
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
                      if (A==2) 
                        rv.proton=0;
                      else 
                        rv.proton=1;
                      nucleusB.push_back(rv);
                      A++;
                      cout << "A=" << A << ", x=" << rv.x << ", y=" << rv.y << endl;
                    }
                }
              
	  	  
              fin.close();
            
              param->setA2FromFile(A);
            }
         else
           {
             cerr << " file he3_plaintext.dat not found. exiting." << endl;
             exit(1);
           }
	}   
      else
	{
          generate_nucleus_configuration(random, A2, Z2,
                                         glauber->GlauberData.Target.a_WS,
                                         glauber->GlauberData.Target.R_WS,
                                         glauber->GlauberData.Target.beta2,
                                         glauber->GlauberData.Target.beta4,
                                         &nucleusB);
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
  // steps in qs0 and Y in the file
  // double y[iymaxNuc];
  // double qs0[ibmax];
  string dummy;
  string T, Qs;
  // open file

  cout << param->getNucleusQsTableFileName() << " ... " ;

      cout << "Reading Q_s(sum(T_p),y) from file ";
      ifstream fin;
      fin.open((param->getNucleusQsTableFileName()).c_str()); 
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
          cout << "[Init.cpp:readNuclearQs]: File " << param->getNucleusQsTableFileName() << " does not exist. Exiting." << endl;
          exit(1);
        }
}  


// void Init::readNuclearQs(Parameters *param)
// {
//   int rank = param->getMPIRank();
//   int size;
//   MPI_Comm_size (MPI_COMM_WORLD, &size);
//   if(rank==0)
//     {
//       double package[iTpmax*iymaxNuc];
//       // steps in qs0 and Y in the file
//       string dummy;
//       string T, Qs;
//       // open file
//       ifstream fin;
//       fin.open((param->getNucleusQsTableFileName()).c_str()); 
      
//       cout << param->getNucleusQsTableFileName() << " ... " ;
      
//       cout << "Reading Q_s(sum(T_p),y) from file ";
//       if(fin)
//         {
//           for (int iT=0; iT<iTpmax; iT++)
//             {
//               for (int iy=0; iy<iymaxNuc; iy++)
//                 {
//                   if (!fin.eof())
//                     {  
//                       fin >> dummy;
//                       fin >> T;
//                       Tlist[iT]=atof(T.c_str());
//                       fin >> Qs;
//                       Qs2Nuclear[iT][iy]=atof(Qs.c_str());
//                       package[iT*iymaxNuc+iy] = Qs2Nuclear[iT][iy];
//                     }
//                   else 
//                     {
//                       cerr << " End of file reached prematurely. Did the file change? Exiting." << endl;
//                       exit(1);
//                     }
//                 }
//             }
//           fin.close();
//           for (int target=1; target<size; target++)
//             {
//               MPI::COMM_WORLD.Send(package,iTpmax*iymaxNuc,MPI::DOUBLE,target,target);
//               MPI::COMM_WORLD.Send(Tlist,iTpmax,MPI::DOUBLE,target,target+size);
//             }
//           cout << " done." << endl;
//         }
//       else
//         {
//           cout << "[Init.cpp:readNuclearQs]: File " << param->getNucleusQsTableFileName() << " does not exist. Exiting." << endl;
//           exit(1);
//         }
//     }
//   else
//     {
//       double package[iTpmax*iymaxNuc];
//       MPI::COMM_WORLD.Recv(package,iTpmax*iymaxNuc,MPI::DOUBLE,0,rank);
//       MPI::COMM_WORLD.Recv(Tlist,iTpmax,MPI::DOUBLE,0,rank+size);
//       for (int iT=0; iT<iTpmax; iT++)
//         {
//           for (int iy=0; iy<iymaxNuc; iy++)
//             {
//               Qs2Nuclear[iT][iy]= package[iT*iymaxNuc+iy];
//             }
//         }
//     }
// }

// Q_s as a function of \sum T_p and y (new in this version of the code - v1.2 and up)
double Init::getNuclearQs2(Parameters *param, Random* random, double T, double y)
  {
  double value, fracy, fracT, QsYdown, QsYup;
  int posy, check=0;
  fracy=0.;
  posy = static_cast<int>(floor(y/deltaYNuc+0.0000001));

  if (y>iymaxNuc*deltaYNuc)
    {
      cout << " [Init:getNuclearQs2]:ERROR: y out of range. Maximum y value is " << iymaxNuc*deltaYNuc << ", you used " << y << ". Exiting." << endl;
      exit(1);
    }

  //  if ( T > Qs2Nuclear[iTpmax-1][iymaxNuc-1] )
  if ( T > Tlist[iTpmax-1] )
    {
      cerr << "T=" << T << ", maximal T in table=" << Tlist[iTpmax-1] << endl;
      cerr << " [Init:getNuclearQs2]:WARNING: out of range. Using maximal T in table." << endl;
      check = 1;
      fracy = (y-static_cast<double>(posy)*deltaYNuc)/deltaYNuc;
      QsYdown = (Qs2Nuclear[iTpmax-1][posy]);
      QsYup =   (Qs2Nuclear[iTpmax-1][posy+1]);
      value = (fracy*QsYup+(1.-fracy)*QsYdown);//*hbarc*hbarc;      
      return value;
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
	  value = (fracy*QsYup+(1.-fracy)*QsYdown);//*hbarc*hbarc;
		  
	  check++;
	  continue;
	}
    }

   if (check!=1)
    {
      cout << check << ": T=" << T << endl ;
      cerr << " [Init:getNuclearQs2]:ERROR: something went wrong in determining the value of Qs^2. Using maximal T_p" << endl;
      value = (fracy*Qs2Nuclear[iTpmax-1][posy+1]+(1.-fracy)*Qs2Nuclear[iTpmax-1][posy]);
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
  // int check=0;
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
  double r;
  double L = param->getL();
  double P,m;
  double rapidity;
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
  double d2 = param->getSigmaNN()/(M_PI*10.);          // in fm^2
  double averageQs = 0.;
  double averageQs2 = 0.;
  double averageQs2Avg = 0.;
  double averageQs2min = 0.;
  double averageQs2min2 = 0.;
  int count = 0;
  double nucleiInAverage;
  nucleiInAverage = static_cast<double>(param->getAverageOverNuclei());
 
  // Arrays to store Q_s fluctuations
  // Make sure that array size is always at least 1 
  int len_quark_array=param->getUseConstituentQuarkProton();
  if (len_quark_array==0) len_quark_array=1;
  double gaussA[A1][len_quark_array];
  double gaussB[A2][len_quark_array];

  for (int i = 0; i<A1; i++) 
    {
      if (param->getUseConstituentQuarkProton()>0)
      {
        for (int iq = 0; iq<param->getUseConstituentQuarkProton(); iq++) 
	  {
	    gaussA[i][iq]=1.;
	  }
      }
      else
        gaussA[i][0]=1.; 
    }    
  for (int i = 0; i<A2; i++) 
    {
      if (param->getUseConstituentQuarkProton()>0)
      {
        for (int iq = 0; iq<param->getUseConstituentQuarkProton(); iq++) 
	  {
	    gaussB[i][iq]=1.;
	  }
      }
      else
        gaussB[i][0]=1.; 
    }
  
  // let the log fluctuate
  if(param->getSmearQs() == 1)
    {
      if (A1>0) 
	{
	  for (int i = 0; i<A1; i++) 
	    {
              // Note: len_quark_array is the number of constituent quarks if useConstituentQuarkProton>0, and 1, if 
              // fluctuations are not included. This way Q_s fluctuations can be implemented for each nucleon even
              // if useConstituentQuarkProton=0
	      for (int iq = 0; iq<len_quark_array; iq++) 
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
	      for (int iq = 0; iq<len_quark_array; iq++) 
		{
		  gaussB[i][iq] = (exp(random->Gauss(0,param->getSmearingWidth())))/1.13;
		  
		  //	      cout << i << " " << iq << " " << gaussB[i][iq] << endl; 
		  //if (gaussB[i]<0)
		  //  gaussB[i]=0.;
		}
	    }
	}
    }

  param->setQsmuRatioB(param->getQsmuRatio());

  if(param->getUseNucleus() == 0)
    {
      for(int ix=0; ix<N; ix++) // loop over all positions
	{
	  for(int iy=0; iy<N; iy++)
	    {
              int localpos = ix*N+iy;
	      lat->cells[localpos]->setg2mu2A(param->getg2mu()*param->getg2mu()/param->getg()/param->getg());
	      lat->cells[localpos]->setg2mu2B(param->getg2mu()*param->getg2mu()/param->getg()/param->getg());
	    }
	}
      param->setSuccess(1);
      cout << "constant color charge density set" << endl;
      return;
    }
  
#pragma omp parallel for
  for(int ix=0; ix<N; ix++) // loop over all positions
    {
      for(int iy=0; iy<N; iy++)
	{
	  int localpos = ix*N+iy;
	  lat->cells[localpos]->setg2mu2A(0.);
	  lat->cells[localpos]->setg2mu2B(0.);
	}
    }
  
  // compute N_part
  // positions are shifted here. not later as in previous versions. bshift below (in init(..)) is zero.

  if(A1 < 4 && A2 > 1) 
    {
      for (int i = 0; i<A2; i++) 
	{
	  nucleusB.at(i).x=nucleusB.at(i).x+b;
	}   
    }
  else if(A2 < 4 && A1 > 1) 
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

  
  double BG;
  BG = param->getBG();
  double BGq = param->getBGq(); // quark size in GeV^-2
  double xi = param->getProtonAnisotropy();

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
  else
    {
      for (int i = 0; i<A1; i++) 
	{
	  nucleusA.at(i).phi = 0.;
	}      
      
      for (int i = 0; i<A2; i++) 
	{
	  nucleusB.at(i).phi = 0.;
	}      
    }

  //  cout << "BG=" << BG << endl;


  double xq[A1][len_quark_array], xq2[A2][len_quark_array];
  double yq[A1][len_quark_array], yq2[A2][len_quark_array];
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
	  for (int iq=0; iq<param->getUseConstituentQuarkProton(); iq++)
	    {
	      avgxq += xq[i][iq];
	      avgyq += yq[i][iq];
	    }
          // Move center of mass to the origin
          // Note that 1607.01711 this is not done, so parameters quoted in
          // that paper can't be used if this is done
	  for (int iq=0; iq<param->getUseConstituentQuarkProton(); iq++)
	    {
              if (param->getShiftConstituentQuarkProtonOrigin())
              {
                xq[i][iq] -= avgxq/double(param->getUseConstituentQuarkProton());
                yq[i][iq] -= avgyq/double(param->getUseConstituentQuarkProton());
              }
	    }
	}
    }
  
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
	  
	  for (int iq=0; iq<param->getUseConstituentQuarkProton(); iq++)
	    {
	      avgxq += xq2[i][iq];
	      avgyq += yq2[i][iq];
	    }
          // Move center of mass to the origin, see comment above
	  for (int iq=0; iq<param->getUseConstituentQuarkProton(); iq++)
	    {
              if (param->getShiftConstituentQuarkProtonOrigin())
              {
                xq2[i][iq] -= avgxq/double(param->getUseConstituentQuarkProton());
                yq2[i][iq] -= avgyq/double(param->getUseConstituentQuarkProton());
              }
	    }
	}
    }


  // test what a smmoth Woods-Saxon would give
  if(param->getUseSmoothNucleus()==1)
    {
      cout << "Using smooth nucleus for test purposes. Does not include deformation." << endl;
      Npart = 2; //avoid break below 
      double xA, xB;
      double y;
      double T;
      double localpos;
      double normA = 0.;
      double normB = 0.;
      double b = param->getb();
      for(int ix=0; ix<N; ix++) // loop over all positions
        {
          xA = -L/2.+a*ix-b/2.;
          xB = -L/2.+a*ix+b/2.;
          for(int iy=0; iy<N; iy++)
            {
              y = -L/2.+a*iy;
              
              localpos = ix*N+iy;
              
              // nucleus A 
              r = sqrt(xA*xA+y*y);
              T = glauber->InterNuTInST(r);
              lat->cells[localpos]->setTpA(T);
              
              normA+=T*a*a;

              // nucleus B
              r = sqrt(xB*xB+y*y);
              T = glauber->InterNuPInSP(r);
              lat->cells[localpos]->setTpB(T);

              normB+=T*a*a;

              //remove potential stuff outside the interaction region
              if (lat->cells[localpos]->getTpA() < 0.001 || lat->cells[localpos]->getTpB() < 0.001)
                {
                  lat->cells[localpos]->setTpA(0.);
                  lat->cells[localpos]->setTpB(0.);
                }

            }
        }
      for(int ix=0; ix<N; ix++) // loop over all positions
        {
          for(int iy=0; iy<N; iy++)
            {
              localpos = ix*N+iy;
              lat->cells[localpos]->setTpA(lat->cells[localpos]->getTpA()/normA*glauber->nucleusA1()*hbarc*hbarc);
              lat->cells[localpos]->setTpB(lat->cells[localpos]->getTpB()/normB*glauber->nucleusA2()*hbarc*hbarc);
            }
        }

      // double normTest=0.;
      // for(int ix=0; ix<N; ix++) // loop over all positions
      //   {
      //     for(int iy=0; iy<N; iy++)
      //       {
      //         localpos = ix*N+iy;
      //         normTest+=lat->cells[localpos]->getTpA()*a*a;
      //       }
      //   }
      // cout << "normTest=" << normTest << endl;
      param->setSuccess(1);
    }
  else
    {
      //add all T_p's (new in version 1.2)
#pragma omp parallel
      {
        double x, xm;
        double y, ym;
        int localpos;
        double bp2,T, phi;
        
#pragma omp for   
        for(int ix=0; ix<N; ix++) // loop over all positions
          {
            x = -L/2.+a*ix;
            for(int iy=0; iy<N; iy++)
              {
                y = -L/2.+a*iy;
                
                localpos = ix*N+iy;
                
                // nucleus A 
                lat->cells[localpos]->setTpA(0.);
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
                            
                            T += exp(-bp2/(2.*BGq))/(2.*M_PI*BGq)/(double(param->getUseConstituentQuarkProton()))*gaussA[i][iq]; // I removed the 2/3 here to make it a bit bigger
                          }
                      }
                    else
                      {
                        phi = nucleusA.at(i).phi;
                        
                        bp2 = (xm-x)*(xm-x)+(ym-y)*(ym-y) + xi*pow((xm-x)*cos(phi) + (ym-y)*sin(phi),2.);
                        bp2 /= hbarc*hbarc;     	  
                        T = sqrt(1+xi)*exp(-bp2/(2.*BG))/(2.*M_PI*BG)*gaussA[i][0]; // T_p in this cell for the current nucleon
                      }
                    lat->cells[localpos]->setTpA(lat->cells[localpos]->getTpA()+T/nucleiInAverage); // add up all T_p
                  }
                
                
                // nucleus B 
                lat->cells[localpos]->setTpB(0.);
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
                            
                            T += exp(-bp2/(2.*BGq))/(2.*M_PI*BGq)/double(param->getUseConstituentQuarkProton())*gaussB[i][iq];
                          }
                      }
                    else
                      {
                        phi = nucleusB.at(i).phi;
                        
                        bp2 = (xm-x)*(xm-x)+(ym-y)*(ym-y) + xi*pow((xm-x)*cos(phi) + (ym-y)*sin(phi),2.);
                        bp2 /= hbarc*hbarc;
                        
                        T = sqrt(1+xi)*exp(-bp2/(2.*BG))/(2.*M_PI*BG)*gaussB[i][0]; // T_p in this cell for the current nucleon
                      }
                    
                    lat->cells[localpos]->setTpB(lat->cells[localpos]->getTpB()+T/nucleiInAverage); // add up all T_p	      
                  }
              }
          }
      }
    }
 
  if(param->getUseSmoothNucleus()==0)
    {
      stringstream strNcoll_name;
      strNcoll_name << "NcollList" << param->getEventId() << ".dat";
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
      
      
      stringstream strNpart_name;
      strNpart_name << "NpartList" << param->getEventId() << ".dat";
      string Npart_name;  Npart_name = strNpart_name.str();
      
      ofstream foutNpart(Npart_name.c_str(),ios::out); 
      
      for (int i = 0; i<A1; i++) 
        {
          foutNpart << nucleusA.at(i).x << " " << nucleusA.at(i).y << " " << nucleusA.at(i).proton << " " << nucleusA.at(i).collided << endl;
        }
      foutNpart << endl;
      for (int i = 0; i<A2; i++) 
        {
          foutNpart << nucleusB.at(i).x << " " << nucleusB.at(i).y << " " << nucleusB.at(i).proton << " " << nucleusB.at(i).collided << endl;
        }
      foutNpart.close();
      
      // in p+p assume that they collided in any case
      if ( A1 == 1 && A2 == 1 )
        {
          nucleusB.at(0).collided=1;
          nucleusA.at(0).collided=1;
        }
      
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
      
      param->setNpart(Npart);
      
      if(param->getUseFixedNpart()!=0 && Npart!=param->getUseFixedNpart())
        {
          cout << "current Npart = " << Npart << endl;
          return;
        }
    }
  // get Q_s^2 (and from that g^2mu^2) for a given \sum T_p and Y
#pragma omp parallel
  {
    //    double x;
    //    double y;
    int localpos;
    double Ydeviation = 10000;
    double QsA, QsB, distanceA, distanceB;
    double xVal;
    double localrapidity=rapidity;
    int check;
    QsA = 1;
    QsB = 1;
     
#pragma omp for   
  for(int ix=0; ix<N; ix++) // loop over all positions
    {
      // x = -L/2.+a*ix;
      for(int iy=0; iy<N; iy++)
	{
          Ydeviation = 10000;
	  check = 0;
          // y = -L/2.+a*iy;
	  localpos = ix*N+iy;
	  

          if(param->getUseSmoothNucleus()==1)
            check=2;
          else
            {

              //	  cut proton at a radius of rmax [fm] (about twice the gluonic radius to be generous)	  
              
              if(log(2*M_PI*BG*lat->cells[localpos]->getTpA())<0.)
                distanceA = sqrt(-2.*BG*log(2*M_PI*BG*lat->cells[localpos]->getTpA()))*hbarc;
              else distanceA=0.;
              
              if(log(2*M_PI*BG*lat->cells[localpos]->getTpB())<0.)
                distanceB = sqrt(-2.*BG*log(2*M_PI*BG*lat->cells[localpos]->getTpB()))*hbarc;
              else distanceB=0.;
              
              if(distanceA < param->getRmax())
                {
                  check=1;
                }
              
              if(distanceB < param->getRmax() && check==1)
                {		
                  check=2;
                }
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
		      if(localrapidity>=0)
                        {
                          QsA = sqrt(getNuclearQs2(param, random, lat->cells[localpos]->getTpA(), abs(localrapidity)));
                        }
		      else 
			{
			  xVal = QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn);
			  if(xVal==0)
			    QsA=0.;
			  else
			    QsA = sqrt(getNuclearQs2(param, random, lat->cells[localpos]->getTpA(), 0.))*
			      sqrt(pow((1-xVal)/(1-0.01),exponent)*pow((0.01/xVal),0.2));
			}
		      if(QsA == 0)
			{
			  Ydeviation = 0;
			  lat->cells[localpos]->setg2mu2A(0.);
			}
		      else
			{
			  // nucleus A 
			  lat->cells[localpos]->setg2mu2A(QsA*QsA/param->getQsmuRatio()/param->getQsmuRatio()
						     *a*a/hbarc/hbarc/param->getg()/param->getg()); // lattice units? check
			  
			  
			  Ydeviation = localrapidity - log(0.01/(QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn)));
			  localrapidity = log(0.01/(QsA*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(yIn)));
			}
		    }
		  if(lat->cells[localpos]->getg2mu2A()!=lat->cells[localpos]->getg2mu2A())
		    {
		      lat->cells[localpos]->setg2mu2A(0.);
		    }
		  
                  localrapidity=rapidity;
		  Ydeviation = 10000;
		  while (abs(Ydeviation) > 0.001) 
		    {
		      if(localrapidity>=0)
			QsB = sqrt(getNuclearQs2(param, random, lat->cells[localpos]->getTpB(), abs(localrapidity)));
		      else
			{
			  xVal = QsB*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(-yIn);
			  if(xVal==0)
			    QsB=0.;
			  else
			    QsB = sqrt(getNuclearQs2(param, random, lat->cells[localpos]->getTpB(), 0.))*
			    sqrt(pow((1-xVal)/(1-0.01),exponent)*pow((0.01/xVal),0.2));
			}
		      if(QsB == 0)
			{
			  Ydeviation = 0;
			  lat->cells[localpos]->setg2mu2B(0.);
			}
		      else
			{
			  // nucleus B 
			  lat->cells[localpos]->setg2mu2B(QsB*QsB/param->getQsmuRatioB()/param->getQsmuRatioB()
						     *a*a/hbarc/hbarc/param->getg()/param->getg()); 			  
			  Ydeviation = localrapidity - log(0.01/(QsB*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(-yIn)));
			  localrapidity = log(0.01/(QsB*param->getxFromThisFactorTimesQs()/param->getRoots()*exp(-yIn)));
			}
		    }   
		  if(lat->cells[localpos]->getg2mu2B()!=lat->cells[localpos]->getg2mu2B())
		    {
		      lat->cells[localpos]->setg2mu2B(0.);
		    }
		  // _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
		  // end iterative loops here
		}
	      else
		{
		  // nucleus A 
		  lat->cells[localpos]->setg2mu2A(getNuclearQs2(param, random, lat->cells[localpos]->getTpA(), localrapidity)/param->getQsmuRatio()/param->getQsmuRatio()
					     *a*a/hbarc/hbarc/param->getg()/param->getg()); // lattice units? check
		  
		  // nucleus B 
		  lat->cells[localpos]->setg2mu2B(getNuclearQs2(param, random, lat->cells[localpos]->getTpB(), localrapidity)/param->getQsmuRatioB()/param->getQsmuRatioB()
					     *a*a/hbarc/hbarc/param->getg()/param->getg()); 
		 
		}
	    }
	}
    }
  }
  

 
  count=0;
  double Tpp=0.;
  double x,xm, y,ym;
  double alphas=0.;
  int check=0;
  for(int ix=0; ix<N; ix++) // loop over all positions
    {
      for(int iy=0; iy<N; iy++)
	{
	  check = 0;
	  pos = ix*N+iy;
	  x = -L/2.+a*ix;
	  y = -L/2.+a*iy;
          //	  outvalue = lat->cells[pos]->getg2mu2A();
	  
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
	  
	  for (int i = 0; i<A1; i++) 
	    {
	      xm = nucleusA.at(i).x;
	      ym = nucleusA.at(i).y;
	      r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
	      
	      if(r<sqrt(0.1*param->getSigmaNN()/M_PI) && nucleusA.at(i).collided==1)
		{
		  check=1;
		}
	    }
	  
	  for (int i = 0; i<A2; i++) 
	    {
	      xm = nucleusB.at(i).x;
	      ym = nucleusB.at(i).y;
	      r = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym));
	      
	      if(r<sqrt(0.1*param->getSigmaNN()/M_PI) && nucleusB.at(i).collided==1 && check==1)
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
	  // compute T_pp
	  Tpp += lat->cells[pos]->getTpB()*lat->cells[pos]->getTpA()*a*a/hbarc/hbarc/hbarc/hbarc; // now this quantity is in fm^-2
	  // remember: Tp is in GeV^2
	}
    }
  
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
  
  if(param->getRunningCoupling() && param->getRunWithkt()==0)
    {
      if(param->getRunWithQs()==2)
	{
	  cout << "running with " << param->getRunWithThisFactorTimesQs() << " Q_s(max)" << endl;
	  alphas = 12.*M_PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2)); // 3 flavors
	  cout << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " Qs_max)=" << alphas << endl;
	}
      else if(param->getRunWithQs()==0)
	{
	  cout << "running with " << param->getRunWithThisFactorTimesQs() << " Q_s(min)" << endl;
	  alphas = 12.*M_PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2)); // 3 flavors
	  cout << "alpha_s(" << param->getRunWithThisFactorTimesQs() << " Qs_min)=" << alphas << endl;
	}
      else if(param->getRunWithQs()==1)
	{
	  cout << "running with " << param->getRunWithThisFactorTimesQs() << " <Q_s>" << endl;
	  alphas = 12.*M_PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2)); // 3 flavors
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
      alphas = param->getg()*param->getg()/4./M_PI;
    }
 
  if(param->getAverageQs() > 0 && param->getAverageQsAvg()>0 && averageQs2>0  && param->getAverageQsmin()>0 && averageQs2Avg>0 && alphas>0 && Npart>=2)
    param->setSuccess(1);
 

  param->setalphas(alphas);
  
  stringstream strup_name;
  strup_name << "usedParameters" << param->getEventId() << ".dat";
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
  fout1.close();
}



void Init::setV(Lattice *lat, Group* group, Parameters *param, Random* random, Glauber *glauber)
{
  cout << "Setting Wilson lines ..." << endl;
  const int N = param->getSize();
  const int Ny=param->getNy();
  const int Nc = param->getNc();
  const int Nc2m1 = Nc*Nc-1;
  const int nn[2] = {N,N};
  const double L = param->getL();
  const double a = L/N; // lattice spacing in fm
  const double m = param->getm()*a/hbarc;
  const Matrix one(Nc,1.);
  double UVdamp = param->getUVdamp(); //GeV^-1
  UVdamp = UVdamp/a*hbarc;
  complex<double>** rhoACoeff;
  rhoACoeff = new complex<double>*[Nc2m1];
  for(int i=0; i<Nc2m1; i++)
    {
      rhoACoeff[i] = new complex<double>[N*N];
    }
 
  // loop over longitudinal direction
  for(int k=0; k<Ny; k++)
    {
      double g2muA;
      for (int pos=0; pos<N*N; pos++)
        {
          for(int n=0; n<Nc2m1; n++)
            {
              g2muA = param->getg()*sqrt(lat->cells[pos]->getg2mu2A()/static_cast<double>(Ny));
              rhoACoeff[n][pos] = g2muA*random->Gauss();
            }
        }
        
        for(int n=0; n<Nc2m1; n++)
          {
            fft->fftnComplex(rhoACoeff[n],rhoACoeff[n],nn,2,1);
          }
        
        // compute A^+
#pragma omp parallel for
        for (int i=0; i<N; i++)
          {
            for (int j=0; j<N; j++)
              {
                double kt2, kx, ky;
                int localpos = i*N+j; 
                kx = 2.*M_PI*(-0.5+static_cast<double>(i)/static_cast<double>(N));
                ky = 2.*M_PI*(-0.5+static_cast<double>(j)/static_cast<double>(N));
                kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice momentum
                if(m==0)
                  {
                    if(kt2!=0)
                      {
                        for(int n=0; n<Nc2m1; n++)
                          {
                            rhoACoeff[n][localpos] =  rhoACoeff[n][localpos]*(1./(kt2));
                          }
                      }
                    else
                      {
                        for(int n=0; n<Nc2m1; n++)
                          {
                            rhoACoeff[n][localpos] = 0.;
                          }
                      }
                  }
                else
                  {
                    for(int n=0; n<Nc2m1; n++)
                      {
                        rhoACoeff[n][localpos] *= (1./(kt2+m*m))*exp(-sqrt(kt2)*UVdamp);
                      }
                  }
              }
          }
        
        // Fourier transform back A^+
        for(int n=0; n<Nc2m1; n++)
          {
            fft->fftnComplex(rhoACoeff[n],rhoACoeff[n],nn,2,-1);
          }
        // compute U
  
#pragma omp parallel
        {      
          double in[8];
          vector <complex<double> > U;
          Matrix temp(Nc,1.);
          Matrix temp2(Nc,0.);
          Matrix tempNew(Nc,0.);
        
          #pragma omp for
             for (int pos=0; pos<N*N; pos++)
               {
                 for (int a=0; a<Nc2m1; a++)
                   {
                     in[a] = -(rhoACoeff[a][pos]).real(); // expmCoeff wil calculate exp(i in[a]t[a]), so just multiply by -1 (not -i)
                   }
                 
                 U = temp2.expmCoeff(in, Nc);
                 
                 tempNew = U[0]*one + U[1]*group->getT(0) + U[2]*group->getT(1) + U[3]*group->getT(2) + U[4]*group->getT(3) + 
                   U[5]*group->getT(4) + U[6]*group->getT(5) + U[7]*group->getT(6) + U[8]*group->getT(7);
                 
                 temp = tempNew * lat->cells[pos]->getU();
                 // set U
                 lat->cells[pos]->setU(temp);
                 
               }
        }

    }//Ny loop
  
  // loop over longitudinal direction
  for(int k=0; k<Ny; k++)
    {
      double g2muB;
      for (int pos=0; pos<N*N; pos++)
        {
          for(int n=0; n<Nc2m1; n++)
            {
              g2muB = param->getg()*sqrt(lat->cells[pos]->getg2mu2B()/static_cast<double>(Ny));
              rhoACoeff[n][pos] = g2muB*random->Gauss();
            }
        }
        
        for(int n=0; n<Nc2m1; n++)
          {
            fft->fftnComplex(rhoACoeff[n],rhoACoeff[n],nn,2,1);
          }
        
        // compute A^+
#pragma omp parallel for
        for (int i=0; i<N; i++)
          {
            for (int j=0; j<N; j++)
              {
                double kt2, kx, ky;
                int localpos = i*N+j; 
                kx = 2.*M_PI*(-0.5+static_cast<double>(i)/static_cast<double>(N));
                ky = 2.*M_PI*(-0.5+static_cast<double>(j)/static_cast<double>(N));
                kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice momentum
                if(m==0)
                  {
                    if(kt2!=0)
                      {
                        for(int n=0; n<Nc2m1; n++)
                          {
                            rhoACoeff[n][localpos] =  rhoACoeff[n][localpos]*(1./(kt2));
                          }
                      }
                    else
                      {
                        for(int n=0; n<Nc2m1; n++)
                          {
                            rhoACoeff[n][localpos] = 0.;
                          }
                      }
                  }
                else
                  {
                    for(int n=0; n<Nc2m1; n++)
                      {
                        rhoACoeff[n][localpos] *= (1./(kt2+m*m))*exp(-sqrt(kt2)*UVdamp);
                      }
                  }
              }
          }
        
        // Fourier transform back A^+
        for(int n=0; n<Nc2m1; n++)
          {
            fft->fftnComplex(rhoACoeff[n],rhoACoeff[n],nn,2,-1);
          }
        // compute U
  
#pragma omp parallel
        {      
          double in[8];
          vector <complex<double> > U;
          Matrix temp(Nc,1.);
          Matrix temp2(Nc,0.);
          Matrix tempNew(Nc,0.);
        
          #pragma omp for
             for (int pos=0; pos<N*N; pos++)
               {
      
                 for (int a=0; a<Nc2m1; a++)
                   {
                     in[a] = -(rhoACoeff[a][pos]).real(); // expmCoeff wil calculate exp(i in[a]t[a]), so just multiply by -1 (not -i)
                   }
                 
                 U = temp2.expmCoeff(in, Nc);
                 
                 
                 tempNew = U[0]*one + U[1]*group->getT(0) + U[2]*group->getT(1) + U[3]*group->getT(2) + U[4]*group->getT(3) + 
                   U[5]*group->getT(4) + U[6]*group->getT(5) + U[7]*group->getT(6) + U[8]*group->getT(7);
                 
                 temp = tempNew * lat->cells[pos]->getU2();
                 // set U
                 lat->cells[pos]->setU2(temp);
               }
        }

    }//Ny loop

  
  // --------
    for(int ic=0; ic<Nc2m1; ic++)
      {
        delete [] rhoACoeff[ic];
      }
    delete [] rhoACoeff;
  
  
        
  // // output U
  if (param->getWriteInitialWilsonLines())
  {
   stringstream strVOne_name;
   //strVOne_name << "V1-" << param->getMPIRank() << ".txt";
   strVOne_name << "V-" <<  param->getEventId() + 2*param->getSeed()*param->getMPISize() << ".txt";
   string VOne_name;
   VOne_name = strVOne_name.str();

   ofstream foutU(VOne_name.c_str(),ios::out); 
   foutU.precision(15);

   for(int ix=0; ix<N; ix++)
     {
       for(int iy=0; iy<N; iy++) // loop over all positions
   	{
   	  int pos = ix*N+iy;
   	  foutU << ix << " " << iy << " "  << (lat->cells[pos]->getU()).MatrixToString() << endl;
   	}
       foutU << endl;
     }
   foutU.close();

   cout<<"wrote " << strVOne_name.str() <<endl;
  
   stringstream strVTwo_name;
   // strVTwo_name << "V2-" << param->getMPIRank() << ".txt";
   strVTwo_name << "V-" <<  param->getEventId() + (1+2*param->getSeed())*param->getMPISize() << ".txt";
   string VTwo_name;
   VTwo_name = strVTwo_name.str();

   ofstream foutU2(VTwo_name.c_str(),ios::out); 
   foutU2.precision(15);
   for(int ix=0; ix<N; ix++)
     {
       for(int iy=0; iy<N; iy++) // loop over all positions
   	{
   	  int pos = ix*N+iy;
   	  foutU2 << ix << " " << iy << " "  << (lat->cells[pos]->getU2()).MatrixToString() << endl;
   	}
       foutU2 << endl;
     }
   foutU2.close();
  
   cout<<"wrote " << strVTwo_name.str() <<endl;
  } 
  // --------


  cout << " Wilson lines V_A and V_B set on rank " << param->getMPIRank() << ". " << endl; 
}


void Init::readV(Lattice *lat, Group* group, Parameters *param)
{
  int pos;
  int N = param->getSize();
  int Nc = param->getNc();
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
  const int maxIterations = 100000;
  const int N = param->getSize();
  const int Nc = param->getNc();
  const int Nc2m1 = Nc*Nc-1;
  const double bmin=param->getbmin();
  const double bmax=param->getbmax();
  const Matrix one(Nc,1.);
  const Matrix zero(Nc,0.);

  messager.info("Initializing fields ... ");
  param->setRnp(0.);
  
  double b; 
  double xb = random->genrand64_real1(); // uniformly distributed random variable                                                               
  
  if(param->getUseNucleus() == 0) // use b=0 fm for the constant g^2 mu case
    {
      param->setSuccess(1);
      b=0.;
      messager << "Setting b=0 for constant color charge density case.";
        messager.flush("info");
    }
  else 
    {
      if(param->getLinearb()==1) // use a linear probability distribution for b if we are doing nuclei
	{
	  messager << "Sampling linearly distributed b between " << bmin << " and " << bmax << "fm. Found ";
	  b = sqrt((bmax*bmax-bmin*bmin)*xb+bmin*bmin);
	}
      else // use a uniform distribution instead
	{
	  messager << "Sampling uniformly distributed b between " << bmin << " and " << bmax << "fm. Found ";
	  b = (bmax-bmin)*xb+bmin;
	}
    }

  param->setb(b);
  messager << "b=" << b << " fm.";
  messager.flush("info");

 
  // read Q_s^2 from file
  if(param->getUseNucleus() == 1)
    {
      readNuclearQs(param);
    }
  
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

  messager.info("Finding fields in forward lightcone...");

#pragma omp parallel
  {
    int countMe;
    int checkConvergence;
    int alphaCheck;

    double Fold;
    double Fnew;
    Fnew = 0.;
    double lambda;

    complex<double>* M         = new complex<double>[Nc2m1*Nc2m1];
    complex<double>* F         = new complex<double>[Nc2m1];
    complex<double>* result    = new complex<double>[Nc2m1];
    complex<double>* alpha     = new complex<double>[Nc2m1];
    complex<double>* alphaSave = new complex<double>[Nc2m1];

    vector <complex<double> > Dalpha;
    Dalpha.reserve(Nc2m1);

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

    
    // compute Ux(3) Uy(3) after the collision
#pragma omp for 
    for (int pos=0; pos<N*N; pos++)      //loops over all cells
      {
        if(lat->cells[pos]->getU().trace()!=lat->cells[pos]->getU().trace())
          {
            lat->cells[pos]->setU(one);
          }
        
        if(lat->cells[pos]->getU2().trace()!=lat->cells[pos]->getU2().trace())
          {
            lat->cells[pos]->setU2(one);
          }
      }        
  
#pragma omp for 
    for (int pos=0; pos<N*N; pos++)      //loops over all cells
      {
        UDx=lat->cells[lat->pospX[pos]]->getU();
        UDx.conjg();
        lat->cells[pos]->setUx1(lat->cells[pos]->getU()*UDx);
        
        UDy=lat->cells[lat->pospY[pos]]->getU();
        UDy.conjg();  
        lat->cells[pos]->setUy1(lat->cells[pos]->getU()*UDy);
        
        UDx=lat->cells[lat->pospX[pos]]->getU2();
        UDx.conjg();
        lat->cells[pos]->setUx2(lat->cells[pos]->getU2()*UDx);
        
        UDy=lat->cells[lat->pospY[pos]]->getU2();
        UDy.conjg();	 
        lat->cells[pos]->setUy2(lat->cells[pos]->getU2()*UDy);
      }      
        // -----------------------------------------------------------------
        // from Ux(1,2) and Uy(1,2) compute Ux(3) and Uy(3):

#pragma omp for 
    for (int pos=0; pos<N*N; pos++)      //loops over all cells
      {
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
        
        lat->cells[pos]->setUx(one); 
        int ni=0;
        // solve for alpha iteratively (U(3)=exp(i alpha_b t^b))
        checkConvergence=1;
        while (ni<maxIterations && checkConvergence)
          {
            ni++;
            //set exponential term
            temp2=lat->cells[pos]->getUx(); // contains exp(i alpha_b t^b)
            expAlpha=temp2;
            temp2.conjg();
            expNegAlpha=temp2; // contains exp(-i alpha_b t^b)
            
            // compute Jacobian
            countMe = 0;
            for(int ai=0; ai<Nc2m1; ai++)
              {
                for(int bi=0; bi<Nc2m1; bi++)
                  {
                    temp = group->getT(ai)*Ux1pUx2*group->getT(bi)*expNegAlpha+group->getT(ai)*expAlpha*group->getT(bi)*UDx1pUDx2;
                    // -i times trace of temp gives my Jacobian matrix elements:
                    M[countMe] = complex<double>(0.,-1.)*temp.trace();
                    countMe++;
                  }
              }
            
            // compute function F that needs to be zero
            for(int ai=0; ai<Nc2m1; ai++)
              {
                temp = group->getT(ai)*(Ux1pUx2-UDx1pUDx2)+group->getT(ai)*Ux1pUx2*expNegAlpha-group->getT(ai)*expAlpha*UDx1pUDx2;
                // minus trace if temp gives -F_ai
                F[ai] = (-1.)*temp.trace();
              }
            Dalpha = solveAxb(param,M,F);
            
            Fold = 0.;
            lambda=1.;
            
#pragma omp simd reduction(+:Fold)
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
                lat->cells[pos]->setUx(one); 
              }
            
            while(alphaCheck==0)
              {
                for(int ai=0; ai<Nc2m1; ai++)
                  {
                    alpha[ai] = alphaSave[ai]+lambda*Dalpha[ai];
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
                
                expAlpha=temp2;
                temp2.conjg();
                expNegAlpha=temp2; // contains exp(-i alpha_b t^b)
                
                for(int ai=0; ai<Nc2m1; ai++)
                  {
                    temp = group->getT(ai)*(Ux1pUx2-UDx1pUDx2)+group->getT(ai)*Ux1pUx2*expNegAlpha-group->getT(ai)*expAlpha*UDx1pUDx2;
                    // minus trace of temp gives -F_ai
                    F[ai] = (-1.)*temp.trace();
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
                    
                    lambda = 1.;
                    alphaCheck=1;
                  }
                
                Fnew = 0.;
#pragma omp simd reduction(+:Fnew)
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
            
            if(Nc==2 && Fnew<0.00000001)
              checkConvergence=0;
            
            if(Nc==3 && Fnew<0.0001)
              checkConvergence=0;
            
            if (Dalpha[0].real()!=Dalpha[0].real())
              checkConvergence=0;
            
            else if (ni==maxIterations-1)
              {
                cout << pos << " result for Ux(3) did not converge for x!" << endl;
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
        
        lat->cells[pos]->setUy(one); 
        
        // ---- done: set U(3) ------------------------------------------
        checkConvergence=1;
        ni=0;
        // solve for alpha iteratively (U(3)=exp(i alpha_b t^b))
        while (ni<maxIterations && checkConvergence)
          {
            ni++;
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
              }
            
            // solve J_{ab} \Dalpha_b = -F_a and do alpha -> alpha+Dalpha
            Dalpha = solveAxb(param,M,F);
	    
            Fold = 0.;
            lambda=1.;
#pragma omp simd reduction(+:Fold)
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
                  }
                
                // ---- set new Uy(3) --------------------------------------------	
                
                for (int a=0; a<Nc2m1; a++)
                  {
                    in[a] = (alpha[a]).real(); // expmCoeff wil calculate exp(i in[a]t[a])
                  }
                
                U = tempNew.expmCoeff(in, Nc);
		
                temp2 = U[0]*one + U[1]*group->getT(0) + U[2]*group->getT(1) + U[3]*group->getT(2) + U[4]*group->getT(3) + 
                  U[5]*group->getT(4) + U[6]*group->getT(5) + U[7]*group->getT(6) + U[8]*group->getT(7);
                
                lat->cells[pos]->setUy(temp2); 
		
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
                  }
                
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
                    
                    lambda = 1.;
                    alphaCheck=1;
                  }
                
                Fnew = 0.;
#pragma omp simd reduction(+:Fnew)
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
            
            if(Nc==2 && Fnew<0.00000001)
              checkConvergence=0;
            if(Nc==3 && Fnew<0.0001)
              checkConvergence=0;
            if (Dalpha[0].real()!=Dalpha[0].real())
              checkConvergence=0;
            
            else if (ni==maxIterations-1)
              {
                cout << pos << " result for Uy(3) did not converge for y!" << endl;
                cout << "last Dalpha = " << endl;
                for(int ai=0; ai<Nc2m1; ai++)
                  {
                    cout << "Dalpha/alpha=" << Dalpha[ai]/alpha[ai] << endl;
                    cout << "Dalpha=" << Dalpha[ai] << endl;
                    cout << param->getAverageQs() << " " << param->getAverageQsAvg() << " " << param->getAverageQsmin() << endl;
                    
                  }
              }
          }//iteration loop
      }//loop over pos

// compute initial electric field
// with minus ax, ay
#pragma omp for
  for (int pos=0; pos<N*N; pos++)
    {
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
	  
      Ux1mUx2 = lat->cells[lat->posmX[pos]]->getUx1()-lat->cells[lat->posmX[pos]]->getUx2();
      UDx1 = lat->cells[lat->posmX[pos]]->getUx1();
      UDx1.conjg();
      UDx2 = lat->cells[lat->posmX[pos]]->getUx2();
      UDx2.conjg();
      UDx1mUDx2 = UDx1 - UDx2;
      
      Ux = lat->cells[lat->posmX[pos]]->getUx();
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
      
      Uy1mUy2 = lat->cells[lat->posmY[pos]]->getUy1()-lat->cells[lat->posmY[pos]]->getUy2();
      UDy1 = lat->cells[lat->posmY[pos]]->getUy1();
      UDy1.conjg();
      UDy2 = lat->cells[lat->posmY[pos]]->getUy2();
      UDy2.conjg();
      UDy1mUDy2 = UDy1 - UDy2;
      
      Uy = lat->cells[lat->posmY[pos]]->getUy();
      UDy = Uy;
      UDy.conjg();
      
      temp2 = temp2 - UDy*Uy1mUy2 + Uy1mUy2 + UDy1mUDy2*Uy - UDy1mUDy2; 
      
      lat->cells[pos]->setE1((1./8.)*temp2);
    }
  
  // with plus ax, ay
#pragma omp for
  for (int pos=0; pos<N*N; pos++)
    {
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
      
      Ux1mUx2 = lat->cells[lat->pospX[pos]]->getUx1()-lat->cells[lat->pospX[pos]]->getUx2();
      UDx1 = lat->cells[lat->pospX[pos]]->getUx1();
      UDx1.conjg();
      UDx2 = lat->cells[lat->pospX[pos]]->getUx2();
      UDx2.conjg();
      UDx1mUDx2 = UDx1 - UDx2;
      
      Ux = lat->cells[lat->pospX[pos]]->getUx();
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
      
      Uy1mUy2 = lat->cells[lat->pospY[pos]]->getUy1()-lat->cells[lat->pospY[pos]]->getUy2();
      UDy1 = lat->cells[lat->pospY[pos]]->getUy1();
      UDy1.conjg();
      UDy2 = lat->cells[lat->pospY[pos]]->getUy2();
      UDy2.conjg();
      UDy1mUDy2 = UDy1 - UDy2;
      
      Uy = lat->cells[lat->pospY[pos]]->getUy();
      UDy = Uy;
      UDy.conjg();
      
      temp2 = temp2 - UDy*Uy1mUy2 + Uy1mUy2 + UDy1mUDy2*Uy - UDy1mUDy2; 
      
      lat->cells[pos]->setE2((1./8.)*temp2);
      
    }
  // compute the plaquette
  #pragma omp for
  for (int pos=0; pos<N*N; pos++)
    {
      UDx = lat->cells[lat->pospY[pos]]->getUx();
      UDy = lat->cells[pos]->getUy();
      UDx.conjg();
      UDy.conjg();
      
      Uplaq = lat->cells[pos]->getUx()*(lat->cells[lat->pospX[pos]]->getUy()*(UDx*UDy));
      lat->cells[pos]->setUplaq(Uplaq);
    }

#pragma omp for
  for (int pos=0; pos<N*N; pos++)
    {
	  AM = (lat->cells[pos]->getE1());//+lat->cells[pos]->getAetaP());
	  AP = (lat->cells[pos]->getE2());//+lat->cells[pos]->getAetaP());
	  // this is pi in lattice units as needed for the evolution. (later, the a^4 gives the right units for the energy density
	  lat->cells[pos]->setpi(complex<double>(0.,-2./param->getg())*(AM)); // factor -2 because I have A^eta (note the 1/8 before) but want \pi (E^z).
	  // lat->cells[pos]->setpi(complex<double>(0.,-1./param->getg())*(AM+AP)); // factor -2 because I have A^eta (note the 1/8 before) but want \pi (E^z). 
    }

#pragma omp for  
  for (int pos=0; pos<N*N; pos++)
    {
      lat->cells[pos]->setE1(zero);
      lat->cells[pos]->setE2(zero);
      lat->cells[pos]->setphi(zero);
      lat->cells[pos]->setUx1(one); // reset the Ux1 to be used for other purposes later
    }

  delete[] M;
  delete[] F;
  delete[] result;
  delete[] alpha;
  delete[] alphaSave;

  }

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
  strtE_name << "totalEnergy" << param->getEventId() << ".dat";
  string tE_name;
  tE_name = strtE_name.str();

  ofstream fout(tE_name.c_str(),ios::out); 
  fout << epsilonSum << endl;
  fout.close();      
}

void Init::generate_nucleus_configuration(
                Random *random,
                int A, int Z, double a_WS, double R_WS, double beta2, double beta4,
                std::vector<ReturnValue> *nucleus) {
    if (std::abs(beta2) < 1e-15 && std::abs(beta4) < 1e-15) {
        generate_nucleus_configuration_with_woods_saxon(
                                                        random, A, Z, a_WS, R_WS, nucleus);
    } else {
        generate_nucleus_configuration_with_deformed_woods_saxon(
                                                                 random, A, Z, a_WS, R_WS, beta2, beta4, nucleus);
    }
}

void Init::generate_nucleus_configuration_with_woods_saxon(
                                        Random *random,
                                        int A, int Z, double a_WS, double R_WS,
                                        std::vector<ReturnValue> *nucleus) {
    std::vector<double> r_array(A, 0.);
    for (int i = 0; i < A; i++) {
        r_array[i] = sample_r_from_woods_saxon(random, a_WS, R_WS);
    }
    std::sort(r_array.begin(), r_array.end());

    std::vector<double> x_array(A, 0.), y_array(A, 0.), z_array(A, 0.);
    const double d_min    = 0.9;
    const double d_min_sq = d_min*d_min;
    for (unsigned int i = 0; i < r_array.size(); i++) {
        double r_i = r_array[i];
        int reject_flag = 0;
        int iter = 0;
        double x_i, y_i, z_i;
        do {
            iter++;
            reject_flag = 0;
            double phi    = 2.*M_PI*random->genrand64_real3();
            double theta  = acos(1. - 2.*random->genrand64_real3());
            x_i = r_i*sin(theta)*cos(phi);
            y_i = r_i*sin(theta)*sin(phi);
            z_i = r_i*cos(theta);
            for (int j = i - 1; j >= 0; j--) {
                if ((r_i - r_array[j])*(r_i - r_array[j]) > d_min_sq) break;
                double dsq = (  (x_i - x_array[j])*(x_i - x_array[j])
                              + (y_i - y_array[j])*(y_i - y_array[j])
                              + (z_i - z_array[j])*(z_i - z_array[j]));
                if (dsq < d_min_sq) {
                    reject_flag = 1;
                    break;
                }
            }
        } while (reject_flag == 1 && iter < 100);
        //if (iter == 100) {
        //    cout << "[Warning] can not find configuration : "
        //         << "r[i] = " << r_i << ", r[i-1] = " << r_array[i-1] << endl;
        //}
        x_array[i] = x_i;
        y_array[i] = y_i;
        z_array[i] = z_i;
    }

    recenter_nucleus(x_array, y_array, z_array);

    for (unsigned int i = 0; i < r_array.size(); i++) {
        ReturnValue rv;
        rv.x = x_array[i];
        rv.y = y_array[i];
        rv.phi = atan2(y_array[i], x_array[i]);
        rv.collided = 0;
        nucleus->push_back(rv);
    }

    std::random_shuffle ( nucleus->begin(), nucleus->end() );
    
    for (unsigned int i = 0; i < r_array.size(); i++) 
      {
        if(i<abs(Z))
          nucleus->at(i).proton = 1;
        else
          nucleus->at(i).proton = 0;
      }

}


double Init::sample_r_from_woods_saxon(Random *random, double a_WS, double R_WS) const {
    double rmaxCut = R_WS + 10.*a_WS;
    double r = 0.;
    do {
        r = rmaxCut*pow(random->genrand64_real3(), 1.0/3.0);
    } while (random->genrand64_real3() > fermi_distribution(r, R_WS, a_WS));
    return(r);
}


double Init::fermi_distribution(double r, double R_WS, double a_WS) const {
    double f = 1./(1. + exp((r - R_WS)/a_WS));
    return (f);
}


void Init::generate_nucleus_configuration_with_deformed_woods_saxon(
                Random *random,
                int A, int Z, double a_WS, double R_WS, double beta2, double beta4,
                std::vector<ReturnValue> *nucleus) {
    std::vector<double> r_array(A, 0.);
    std::vector<double> costheta_array(A, 0.);
    std::vector<std::pair<double, double>> pair_array;
    for (int i = 0; i < A; i++) {
        sample_r_and_costheta_from_deformed_woods_saxon(
                random, a_WS, R_WS, beta2, beta4,
                r_array[i], costheta_array[i]);
        pair_array.push_back(std::make_pair(r_array[i], costheta_array[i]));
    }
    //std::sort(r_array.begin(), r_array.end());
    std::sort(pair_array.begin(), pair_array.end());

    std::vector<double> x_array(A, 0.), y_array(A, 0.), z_array(A, 0.);
    const double d_min    = 0.9;
    const double d_min_sq = d_min*d_min;
    for (unsigned int i = 0; i < abs(A); i++) {
        //const double r_i     = r_array[i];
        //const double theta_i = acos(costheta_array[i]);
        const double r_i     = pair_array[i].first;
        const double theta_i = acos(pair_array[i].second);
        int reject_flag = 0;
        int iter = 0;
        double x_i, y_i, z_i;
        do {
            iter++;
            reject_flag = 0;
            double phi    = 2.*M_PI*random->genrand64_real3();
            x_i = r_i*sin(theta_i)*cos(phi);
            y_i = r_i*sin(theta_i)*sin(phi);
            z_i = r_i*cos(theta_i);
            for (int j = i - 1; j >= 0; j--) {
                //if ((r_i - r_array[j])*(r_i - r_array[j]) > d_min_sq) break;
                if ((r_i - pair_array[j].first)*(r_i - pair_array[j].first) > d_min_sq) break;
                double dsq = (  (x_i - x_array[j])*(x_i - x_array[j])
                              + (y_i - y_array[j])*(y_i - y_array[j])
                              + (z_i - z_array[j])*(z_i - z_array[j]));
                if (dsq < d_min_sq) {
                    reject_flag = 1;
                    break;
                }
            }
        } while (reject_flag == 1 && iter < 100);
        x_array[i] = x_i;
        y_array[i] = y_i;
        z_array[i] = z_i;
    }
    recenter_nucleus(x_array, y_array, z_array);

    double phi   = 2.*M_PI*random->genrand64_real3();
    double theta = acos(1. - 2.*random->genrand64_real3());
    rotate_nucleus(phi, theta, x_array, y_array, z_array);


    for (unsigned int i = 0; i < r_array.size(); i++) {
        ReturnValue rv;
        rv.x = x_array[i];
        rv.y = y_array[i];
        rv.phi = atan2(y_array[i], x_array[i]);
        rv.collided = 0;
        nucleus->push_back(rv);
    }

    std::random_shuffle ( nucleus->begin(), nucleus->end() );
    
    for (unsigned int i = 0; i < r_array.size(); i++) 
      {
        if(i<abs(Z))
          nucleus->at(i).proton = 1;
        else
          nucleus->at(i).proton = 0;
      }

}


void Init::sample_r_and_costheta_from_deformed_woods_saxon(
        Random *random, double a_WS, double R_WS, double beta2, double beta4,
        double &r, double &costheta) const {
    double rmaxCut = R_WS + 10.*a_WS;
    double R_WS_theta = R_WS;
    do {
        r = rmaxCut*pow(random->genrand64_real3(), 1.0/3.0);
        costheta = 1.0 - 2.0*random->genrand64_real3();
        auto y20 = spherical_harmonics(2, costheta);
        auto y40 = spherical_harmonics(4, costheta);
        R_WS_theta = R_WS*(1.0 + beta2*y20 + beta4*y40);
    } while (random->genrand64_real3()
             > fermi_distribution(r, R_WS_theta, a_WS));
}


double Init::spherical_harmonics(int l, double ct) const {
    // Currently assuming m=0 and available for Y_{20} and Y_{40}
    // "ct" is cos(theta)
    double ylm = 0.0;
    if (l == 2) {
        ylm = 3.0*ct*ct-1.0;
        ylm *= 0.31539156525252005;  // pow(5.0/16.0/M_PI,0.5);
    } else if (l == 4) {
        ylm  = 35.0*ct*ct*ct*ct;
        ylm -= 30.0*ct*ct;
        ylm += 3.0;
        ylm *= 0.10578554691520431;  // 3.0/16.0/pow(M_PI,0.5);
    }
    return(ylm);
}


void Init::recenter_nucleus(std::vector<double> &x, std::vector<double> &y,
                            std::vector<double> &z) {
    // compute the center of mass position and shift it to (0, 0, 0)
    double meanx = 0., meany = 0., meanz = 0.;
    for (unsigned int i = 0; i < x.size(); i++) {
        meanx += x[i];
        meany += y[i];
        meanz += z[i];
    }
      
    meanx /= static_cast<double>(x.size());
    meany /= static_cast<double>(y.size());
    meanz /= static_cast<double>(z.size());

    for (unsigned int i = 0; i < x.size(); i++) {
        x[i] -= meanx;
        y[i] -= meany;
        z[i] -= meanz;
    }
}


void Init::rotate_nucleus(double phi, double theta,
                          std::vector<double> &x, std::vector<double> &y,
                          std::vector<double> &z) {
    auto cth  = cos(theta);
    auto sth  = sin(theta);
    auto cphi = cos(phi);
    auto sphi = sin(phi);
    for (unsigned int i = 0; i < x.size(); i++) {
        auto x_new = cth*cphi*x[i] - sphi*y[i] + sth*cphi*z[i];
        auto y_new = cth*sphi*x[i] + cphi*y[i] + sth*sphi*z[i];
        auto z_new = -sth    *x[i] + 0.  *y[i] + cth     *z[i];
        x[i] = x_new; y[i] = y_new; z[i] = z_new;
    }
}

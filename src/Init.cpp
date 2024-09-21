// Init.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Init.h"
#include "Phys_consts.h"
#include <algorithm>
#include <utility>

using namespace std;
using PhysConst::hbarc;

//**************************************************************************
// Init class.

vector<double> Init::solveAxb(double *Jab, double *Fa) {
  vector<double> xvec(Nc2m1_, 0.);

  gsl_matrix_complex_view m =
      gsl_matrix_complex_view_array(Jab, Nc2m1_, Nc2m1_);
  gsl_vector_complex_view c = gsl_vector_complex_view_array(Fa, Nc2m1_);
  gsl_vector_complex *x = gsl_vector_complex_alloc(Nc2m1_);

  int s;
  gsl_permutation *p = gsl_permutation_alloc(Nc2m1_);
  gsl_linalg_complex_LU_decomp(&m.matrix, p, &s);
  gsl_linalg_complex_LU_solve(&m.matrix, p, &c.vector, x);
  gsl_permutation_free(p);

  for (int i = 0; i < Nc2m1_; i++) {
    xvec[i] = GSL_REAL(gsl_vector_complex_get(x, i));
  }

  gsl_vector_complex_free(x);
  return (xvec);
}


void Init::sampleTA(Parameters *param, Random *random, Glauber *glauber) {
  ReturnValue rv, rv2;
  messager.info("Sampling nucleon positions ... ");

  if (param->getNucleonPositionsFromFile() == 0) {
    int A1, A2, Z1, Z2;
    A1 = static_cast<int>(glauber->nucleusA1()) *
         param->getAverageOverNuclei(); // projectile
    A2 = static_cast<int>(glauber->nucleusA2()) *
         param->getAverageOverNuclei(); // target
    Z1 = static_cast<int>(glauber->nucleusZ1()) *
         param->getAverageOverNuclei(); // projectile
    Z2 = static_cast<int>(glauber->nucleusZ2()) *
         param->getAverageOverNuclei(); // target

    if (param->getAverageOverNuclei() > 1) {
        if ((glauber->nucleusA1() == 1 || glauber->nucleusA2() == 1)) {
            cerr << "Averaging not supported for collisions involving protons "
                    "... Exiting."
                 << std::endl;
            exit(1);
        }
    }

    if (A1 == 1) {
      rv.x = 0.;
      rv.y = 0;
      rv.z = 0;
      rv.collided = 0;
      rv.proton = 1;
      nucleusA_.push_back(rv);
    } else if (A1 == 2) {
      // deuteron
      rv = glauber->SampleTARejection(random, 1);
      param->setRnp(sqrt(rv.x * rv.x + rv.y * rv.y));
      // we sample the neutron proton distance, so distance to the center needs
      // to be divided by 2
      rv.x = rv.x / 2.;
      rv.y = rv.y / 2.;
      rv.z = 0.;
      rv.proton = 1;
      rv.collided = 0;
      nucleusA_.push_back(rv);
      // other nucleon is 180 degrees rotated:
      rv.x = -rv.x;
      rv.y = -rv.y;
      rv.z = -rv.z;
      rv.collided = 0;
      nucleusA_.push_back(rv);
    } else {
      generate_nucleus_configuration(
          random, A1, Z1,
          glauber->GlauberData.Projectile.a_WS,
          glauber->GlauberData.Projectile.R_WS,
          glauber->GlauberData.Projectile.beta2,
          glauber->GlauberData.Projectile.beta3,
          glauber->GlauberData.Projectile.beta4,
          glauber->GlauberData.Projectile.gamma,
          glauber->GlauberData.Projectile.forceDminFlag,
          glauber->GlauberData.Projectile.d_min,
          glauber->GlauberData.Projectile.dR_np,
          glauber->GlauberData.Projectile.da_np,
          nucleusA_);
    }

    if (A2 == 1) {
      rv2.x = 0.;
      rv2.y = 0;
      rv2.z = 0;
      rv2.collided = 0;
      rv2.proton = 1;
      nucleusB_.push_back(rv2);
    } else if (A2 == 2) {
      // deuteron
      rv = glauber->SampleTARejection(random, 2);
      // we sample the neutron proton distance, so distance to the center needs
      // to be divided by 2
      param->setRnp(sqrt(rv.x * rv.x + rv.y * rv.y));

      rv.x = rv.x / 2.;
      rv.y = rv.y / 2.;
      rv.z = 0.;
      rv.proton = 1;
      rv.collided = 0;
      nucleusB_.push_back(rv);

      // other nucleon is 180 degrees rotated:
      rv.x = -rv.x;
      rv.y = -rv.y;
      rv.z = -rv.z;
      rv.proton = 0;
      rv.collided = 0;
      nucleusB_.push_back(rv);
    } else {
      generate_nucleus_configuration(
          random, A2, Z2,
          glauber->GlauberData.Target.a_WS,
          glauber->GlauberData.Target.R_WS,
          glauber->GlauberData.Target.beta2,
          glauber->GlauberData.Target.beta3,
          glauber->GlauberData.Target.beta4,
          glauber->GlauberData.Target.gamma,
          glauber->GlauberData.Target.forceDminFlag,
          glauber->GlauberData.Target.d_min,
          glauber->GlauberData.Target.dR_np,
          glauber->GlauberData.Target.da_np,
          nucleusB_);
    }
  } else if (param->getNucleonPositionsFromFile() == 1) {
      if (nucleonPosArrA_.size() > 0) {
          double ran2 =random->genrand64_real3();
          int nucleusNumber = static_cast<int>(ran2 * nucleonPosArrA_.size());
          std::cout << "using nucleus Number = " << nucleusNumber << std::endl;
          for (int iA = 0; iA < glauber->nucleusA1(); iA++) {
              rv.x = nucleonPosArrA_[nucleusNumber][3*iA];
              rv.y = nucleonPosArrA_[nucleusNumber][3*iA + 1];
              rv.z = nucleonPosArrA_[nucleusNumber][3*iA + 2];
              rv.collided = 0;
              if (iA % 2 == 0) {
                rv.proton = 0;
              } else {
                rv.proton = 1;
              }
              nucleusA_.push_back(rv);
          }
          assignProtons(nucleusA_, glauber->nucleusZ1());
          recenter_nucleus(nucleusA_);
      } else {
          // no configurations, sample with Woods-Saxon
          messager << "configuration file for A = " << glauber->nucleusA1()
                   << " is not available, generate the nucleus configuration "
                   << "using Woods-Saxon distribution instead.";
          messager.flush("info");
          generate_nucleus_configuration(
              random, glauber->nucleusA1(), glauber->nucleusZ1(),
              glauber->GlauberData.Projectile.a_WS,
              glauber->GlauberData.Projectile.R_WS,
              glauber->GlauberData.Projectile.beta2,
              glauber->GlauberData.Projectile.beta3,
              glauber->GlauberData.Projectile.beta4,
              glauber->GlauberData.Projectile.gamma,
              glauber->GlauberData.Projectile.forceDminFlag,
              glauber->GlauberData.Projectile.d_min,
              glauber->GlauberData.Projectile.dR_np,
              glauber->GlauberData.Projectile.da_np,
              nucleusA_);
      }

      if (nucleonPosArrB_.size() > 0) {
          double ran2 =random->genrand64_real3();
          int nucleusNumber = static_cast<int>(ran2 * nucleonPosArrB_.size());
          std::cout << "using nucleus Number = " << nucleusNumber << std::endl;
          for (int iA = 0; iA < glauber->nucleusA2(); iA++) {
              rv.x = nucleonPosArrB_[nucleusNumber][3*iA];
              rv.y = nucleonPosArrB_[nucleusNumber][3*iA + 1];
              rv.z = nucleonPosArrB_[nucleusNumber][3*iA + 2];
              rv.collided = 0;
              if (iA % 2 == 0) {
                rv.proton = 0;
              } else {
                rv.proton = 1;
              }
              nucleusB_.push_back(rv);
          }
          assignProtons(nucleusB_, glauber->nucleusZ2());
          recenter_nucleus(nucleusB_);
      } else {
          // no configurations, sample with Woods-Saxon
          messager << "configuration file for A = " << glauber->nucleusA2()
                   << " is not available, generate the nucleus configuration "
                   << "using Woods-Saxon distribution instead.";
          messager.flush("info");
          generate_nucleus_configuration(
              random, glauber->nucleusA2(), glauber->nucleusZ2(),
              glauber->GlauberData.Target.a_WS,
              glauber->GlauberData.Target.R_WS,
              glauber->GlauberData.Target.beta2,
              glauber->GlauberData.Target.beta3,
              glauber->GlauberData.Target.beta4,
              glauber->GlauberData.Target.gamma,
              glauber->GlauberData.Target.forceDminFlag,
              glauber->GlauberData.Target.d_min,
              glauber->GlauberData.Target.dR_np,
              glauber->GlauberData.Target.da_np,
              nucleusB_);
      }
  } else if (param->getNucleonPositionsFromFile() == 2) {
    // Read in Alvioli's nucleon positions including correlations
    if (glauber->nucleusA1() != 208 && glauber->nucleusA2() != 208) {
      cerr << "[Init.cpp]: The option 'getNucleonPositionsFromFile == 2' only "
              "works for either both nuclei Pb-208 or Projectile p and Target "
              "Pb-208. Exiting."
           << std::endl;
      exit(1);
    }

    std::cout << "Retrieving nuclei from " << std::endl;

    // generate the file name
    double ran = random->genrand64_real3(); // sample the file name uniformly
    int fileNumber = static_cast<int>(ran * 10 + 1);

    stringstream str_file;
    str_file.str("");
    if (fileNumber < 10) {
      str_file << "/global/homes/s/schenke/Alvioli-Pb208/pb208-0";
    } else {
      str_file << "/global/homes/s/schenke/Alvioli-Pb208/pb208-";
    }
    str_file << fileNumber;
    str_file << ".dat";
    string fileName = str_file.str();

    // open the file
    ifstream fin;
    fin.open(fileName.c_str());
    if (!fin) {
      cerr << "File " << fileName
           << " not found. Trying alternative location:" << endl;
      str_file.str("");
      if (fileNumber < 10)
        str_file << "./Alvioli-Pb208/pb208-0";
      else
        str_file << "./Alvioli-Pb208/pb208-";
      str_file << fileNumber;
      str_file << ".dat";
      fileName = str_file.str();
      fin.open(fileName.c_str());
    }

    if (!fin) {
      cerr << "File " << fileName << " not found. Exiting." << endl;
      exit(1);
    }

    cout << "Reading nucleon positions for nuceus A from file " << fileName
         << " ... " << endl;

    // sample the position in the file
    // sample the position in the file uniformly (10,000 events per file)
    double ran2 = random->genrand64_real3();
    int nucleusNumber = static_cast<int>(ran2 * 10000);
    cout << "Nucleus Number = " << nucleusNumber << endl;

    int A = 0;
    int A2 = 0;
    double dummy;

    // go to the correct line in the file
    fin.seekg(std::ios::beg);
    for (int i = 0; i < (nucleusNumber)*glauber->nucleusA1(); ++i) {
      fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    // am now at the correct line in the file

    // start reading one nucleus (208 positions)
    if (glauber->nucleusA1() == 1) {
      rv.x = 0;
      rv.y = 0;
      rv.z = 0;
      rv.collided = 0;
      nucleusA_.push_back(rv);
      A = 1;
    } else {
      while (A < glauber->nucleusA1()) {
        if (!fin.eof()) {
          fin >> rv.x;
          fin >> rv.y;
          fin >> rv.z;
          fin >> dummy; // don't care about isospin
          rv.collided = 0;
          nucleusA_.push_back(rv);
          A++;
          //        cout << "A=" << A << "/" <<
          //glauber->nucleusA1()<<endl; cout << rv.x << " " << rv.y << endl;
        }
      }
    }
    fin.close();
    // do the second nucleus (Target)

    // generate the file name
    ran = random->genrand64_real3(); // sample the file name uniformly
    fileNumber = static_cast<int>(ran * 10 + 1);

    str_file.str("");
    if (fileNumber < 10)
      str_file << "/global/homes/s/schenke/Alvioli-Pb208/pb208-0";
    else
      str_file << "/global/homes/s/schenke/Alvioli-Pb208/pb208-";
    str_file << fileNumber;
    str_file << ".dat";
    fileName = str_file.str();

    // open the file
    fin.open(fileName.c_str());
    if (!fin) {
      cerr << "File " << fileName
           << " not found. Trying alternative location:" << endl;
      str_file.str("");
      if (fileNumber < 10)
        str_file << "./Alvioli-Pb208/pb208-0";
      else
        str_file << "./Alvioli-Pb208/pb208-";
      str_file << fileNumber;
      str_file << ".dat";
      fileName = str_file.str();
      fin.open(fileName.c_str());
    }

    if (!fin) {
      cerr << "File " << fileName << " not found. Exiting." << endl;
      exit(1);
    }

    cout << "Reading nucleon positions for nucleus B from file " << fileName
         << " ... " << endl;

    // sample the position in the file
    ran2 = random->genrand64_real3(); // sample the position in the file
                                      // uniformly (10,000 events per file)
    nucleusNumber = static_cast<int>(ran2 * 10000);
    cout << "Nucleus Number = " << nucleusNumber << endl;

    // go to the correct line in the file
    fin.seekg(std::ios::beg);
    for (int i = 0; i < (nucleusNumber)*glauber->nucleusA1(); ++i) {
      fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    // am now at the correct line in the file

    // start reading one nucleus (208 positions)
    if (glauber->nucleusA2() == 1) {
      rv.x = 0;
      rv.y = 0;
      rv.z = 0;
      rv.collided = 0;
      nucleusB_.push_back(rv);
      A2 = 1;
    } else {
      while (A2 < glauber->nucleusA2()) {
        if (!fin.eof()) {
          fin >> rv.x;
          fin >> rv.y;
          fin >> rv.z; // don't care about z direction
          fin >> dummy; // don't care about isospin
          rv.collided = 0;
          nucleusB_.push_back(rv);
          A2++;
          //    cout << "A2=" << A2 << "/" <<
          //glauber->nucleusA2()<<endl; cout << rv.x << " " << rv.y << endl;
        }
      }
    }

    fin.close();

    param->setA1FromFile(A);
    param->setA2FromFile(A2);
  } else {
    cerr << "NucleonPositionsFromFile can be 0 (sample nucleons) or 1 or 2 "
            "(read from files) - you chose "
         << param->getNucleonPositionsFromFile() << ". Exiting."
         << std::endl;
    exit(1);
  }

  // global rotation of the nucleus
  //rotate_nucleus(random, nucleusA_);
  //rotate_nucleus(random, nucleusB_);
  rotate_nucleus_3D(random, nucleusA_);
  rotate_nucleus_3D(random, nucleusB_);
}


void Init::readNuclearQs(Parameters *param) {
  // steps in qs0 and Y in the file
  // double y[iymaxNuc];
  // double qs0[ibmax];
  string dummy;
  string T, Qs;
  // open file

  std::cout << "Reading Q_s(sum(T_p),y) from file ";

  std::cout << param->getNucleusQsTableFileName() << " ... " << std::endl;

  ifstream fin;
  fin.open((param->getNucleusQsTableFileName()).c_str());
  if (fin) {
    for (int iT = 0; iT < iTpmax; iT++) {
      for (int iy = 0; iy < iymaxNuc; iy++) {
        if (!fin.eof()) {
          fin >> dummy;
          fin >> T;
          Tlist[iT] = atof(T.c_str());
          fin >> Qs;
          Qs2Nuclear[iT][iy] = atof(Qs.c_str());
        } else {
          cerr << " End of file reached prematurely. Did the file change? "
                  "Exiting."
               << std::endl;
          exit(1);
        }
      }
    }
    fin.close();
  } else {
    std::cout << "[Init.cpp:readNuclearQs]: File "
         << param->getNucleusQsTableFileName() << " does not exist. Exiting."
         << std::endl;
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
//                       cerr << " End of file reached prematurely. Did the file
//                       change? Exiting." << endl; exit(1);
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
//           cout << "[Init.cpp:readNuclearQs]: File " <<
//           param->getNucleusQsTableFileName() << " does not exist. Exiting."
//           << endl; exit(1);
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


void Init::readInNucleusConfigs(const int nucleusA,
                                const int lightNucleusOption,
                                vector< vector<float> > &nucleonPosArr) {
    if (nucleonPosArr.size() > 0) return;
    std::string path = "nucleusConfigurations/";
    std::string fileName;
    bool readFlag = true;
    if (nucleusA == 3) {
        fileName = "He3.bin.in";
    } else if (nucleusA == 12) {
        if (lightNucleusOption == 2) {
            fileName = "C12_VMC.bin.in";
        } else if (lightNucleusOption == 3) {
            fileName = "C12_alphaCluster.bin.in";
        }
    } else if (nucleusA == 16) {
        if (lightNucleusOption == 2) {
            fileName = "O16_VMC.bin.in";
        } else if (lightNucleusOption == 3) {
            fileName = "O16_alphaCluster.bin.in";
        } else if (lightNucleusOption == 4) {
            fileName = "O16_PGCM.bin.in";
        } else if (lightNucleusOption == 5) {
            fileName = "O16_NLEFT.bin.in";
        }
    } else if (nucleusA == 20) {
        fileName = "Ne20_PGCM.bin.in";
    } else if (nucleusA == 40) {
        fileName = "Ar40_VMC.bin.in";
    } else {
        readFlag = false;
    }

    if (!readFlag) return;

    fileName = path + fileName;
    messager << "read in nucleus configurations from " << fileName;
    messager.flush("info");
    std::ifstream inFile(fileName, std::ios::binary);
    while (true) {
        vector<float> tempPos;
        for (int i = 0; i < nucleusA; i++) {
            for (int j = 0; j < 3; j++) {
                float temp;
                inFile.read(reinterpret_cast<char*>(&temp), sizeof(float));
                tempPos.push_back(temp);
            }
        }
        if (inFile.eof())
            break;
        nucleonPosArr.push_back(tempPos);
    }
    inFile.close();
    messager << "read in " << nucleonPosArr.size() << " configurations.";
    messager.flush("info");
}


void Init::samplePartonPositions(Parameters *param, Random *random,
                                 vector<double> &x_array,
                                 vector<double> &y_array,
                                 vector<double> &z_array,
                                 vector<double> &BGq_array) {
    const double sqrtBG = sqrt(param->getBG())*hbarc;    // fm
    const double BGqMean = param->getBGq();
    const double BGqVar = param->getBGqVar();
    const double BGq = (
        0.09 + sampleLogNormalDistribution(random, BGqMean - 0.09, BGqVar));
    const double QsSmearWidth = param->getSmearingWidth();
    const int Nq = sampleNumberOfPartons(random, param);
    const double dq_min = param->getDqmin();             // fm
    const double dq_min_sq = dq_min*dq_min;

    vector<double> r_array(Nq, 0.);
    BGq_array.resize(Nq, BGq);
    for (int iq = 0; iq < Nq; iq++) {
        double xq = sqrtBG*random->Gauss();
        double yq = sqrtBG*random->Gauss();
        double zq = sqrtBG*random->Gauss();
        r_array[iq] = sqrt(xq*xq + yq*yq + zq*zq);
    }
    std::sort(r_array.begin(), r_array.end());

    x_array.resize(Nq, 0.);
    y_array.resize(Nq, 0.);
    z_array.resize(Nq, 0.);
    for (unsigned int i = 0; i < r_array.size(); i++) {
        double r_i = r_array[i];
        int reject_flag = 0;
        int iter = 0;
        double x_i, y_i, z_i;
        do {
            iter++;
            reject_flag  = 0;
            double phi   = 2.*M_PI*random->genrand64_real2();
            double theta = acos(1. - 2.*random->genrand64_real2());
            x_i = r_i*sin(theta)*cos(phi);
            y_i = r_i*sin(theta)*sin(phi);
            z_i = r_i*cos(theta);
            for (int j = i - 1; j >= 0; j--) {
                if ((r_i - r_array[j])*(r_i - r_array[j]) > dq_min_sq) break;
                double dsq = (  (x_i - x_array[j])*(x_i - x_array[j])
                              + (y_i - y_array[j])*(y_i - y_array[j])
                              + (z_i - z_array[j])*(z_i - z_array[j]));
                if (dsq < dq_min_sq) {
                    reject_flag = 1;
                    break;
                }
            }
        } while (reject_flag == 1 && iter < 100);
        x_array[i] = x_i;
        y_array[i] = y_i;
        z_array[i] = z_i;
    }
    double avgxq = 0.;
    double avgyq = 0.;
    double avgzq = 0.;
    if (param->getShiftConstituentQuarkProtonOrigin()) {
        for (int iq = 0; iq < Nq; iq++) {
            avgxq += x_array[iq];
            avgyq += y_array[iq];
            avgzq += z_array[iq];
        }
        avgxq /= static_cast<double>(Nq);
        avgyq /= static_cast<double>(Nq);
        avgzq /= static_cast<double>(Nq);
        for (int iq = 0; iq < Nq; iq++) {
            x_array[iq] -= avgxq;
            y_array[iq] -= avgyq;
            z_array[iq] -= avgzq;
        }
    }
}


// Q_s as a function of \sum T_p and y (new in this version of the code - v1.2
// and up)
double Init::getNuclearQs2(double T, double y) {
  double value, fracy, fracT, QsYdown, QsYup;
  int posy, check = 0;
  fracy = 0.;
  posy = static_cast<int>(floor(y / deltaYNuc + 0.0000001));

  if (y > iymaxNuc * deltaYNuc) {
    cout << " [Init:getNuclearQs2]:ERROR: y out of range. Maximum y value is "
         << iymaxNuc * deltaYNuc << ", you used " << y << ". Exiting." << endl;
    exit(1);
  }

  //  if ( T > Qs2Nuclear[iTpmax-1][iymaxNuc-1] )
  if (T > Tlist[iTpmax - 1]) {
    cerr << "T=" << T << ", maximal T in table=" << Tlist[iTpmax - 1] << endl;
    cerr << " [Init:getNuclearQs2]:WARNING: out of range. Using maximal T in "
            "table."
         << endl;
    check = 1;
    fracy = (y - static_cast<double>(posy) * deltaYNuc) / deltaYNuc;
    QsYdown = (Qs2Nuclear[iTpmax - 1][posy]);
    QsYup = (Qs2Nuclear[iTpmax - 1][posy + 1]);
    value = (fracy * QsYup + (1. - fracy) * QsYdown); //*hbarc*hbarc;
    return value;
  }

  if (T < Tlist[0]) {
    check = 1;
    return 0.;
  }

  for (int iT = 0; iT < iTpmax; iT++) {
    if (T >= Tlist[iT] && T < Tlist[iT + 1]) {
      fracT = (T - Tlist[iT]) / (Tlist[iT + 1] - Tlist[iT]);
      fracy = (y - static_cast<double>(posy) * deltaYNuc) / deltaYNuc;

      QsYdown = (fracT) * (Qs2Nuclear[iT + 1][posy]) +
                (1. - fracT) * (Qs2Nuclear[iT][posy]);
      QsYup = (fracT) * (Qs2Nuclear[iT + 1][posy + 1]) +
              (1. - fracT) * (Qs2Nuclear[iT][posy + 1]);
      value = (fracy * QsYup + (1. - fracy) * QsYdown); //*hbarc*hbarc;

      check++;
      continue;
    }
  }

  if (check != 1) {
    cout << check << ": T=" << T << endl;
    cerr << " [Init:getNuclearQs2]:ERROR: something went wrong in determining "
            "the value of Qs^2. Using maximal T_p"
         << endl;
    value = (fracy * Qs2Nuclear[iTpmax - 1][posy + 1] +
             (1. - fracy) * Qs2Nuclear[iTpmax - 1][posy]);
  }

  return value;
}

// set g^2\mu^2 as the sum of the individual nucleons' g^2\mu^2, using Q_s(b,y)
// prop tp g^mu(b,y) also compute N_part using Glauber
void Init::setColorChargeDensity(Lattice *lat, Parameters *param,
                                 Random *random, Glauber *glauber) {
  std::cout << "set color charge density ..." << std::endl;
  int pos, posA, posB;
  int N = param->getSize();
  const int A1 = nucleusA_.size();
  const int A2 = nucleusB_.size();
  // int check=0;
  //if (param->getNucleonPositionsFromFile() == 2) {
  //  A1 = param->getA1FromFile();
  //  A2 = param->getA2FromFile();
  //} else {
  //  A1 = static_cast<int>(glauber->nucleusA1()) * param->getAverageOverNuclei();
  //  A2 = static_cast<int>(glauber->nucleusA2()) * param->getAverageOverNuclei();
  //}

  int Npart = 0;
  int Ncoll = 0;
  double g2mu2A, g2mu2B;
  const double impact_b = param->getb();
  double r;
  double L = param->getL();
  double P, m;
  double rapidity;
  if (param->getUsePseudoRapidity() == 0)
    rapidity = param->getRapidity();
  else {
    // when using pseudorapidity as input convert to rapidity here. later
    // include Jacobian in multiplicity and energy
    cout << "Using pseudorapidity " << param->getRapidity() << endl;
    m = param->getJacobianm();                               // in GeV
    P = 0.13 + 0.32 * pow(param->getRoots() / 1000., 0.115); // in GeV
    rapidity =
        0.5 *
        log(sqrt(pow(cosh(param->getRapidity()), 2.) + m * m / (P * P)) +
            sinh(param->getRapidity()) /
                (sqrt(pow(cosh(param->getRapidity()), 2.) + m * m / (P * P)) -
                 sinh(param->getRapidity())));
    cout << "Corresponds to rapidity " << rapidity << endl;
  }

  double yIn = rapidity; // param->getRapidity();
  double a = L / N;      // lattice spacing in fm
  double dx, dy, dij;
  double d2 = param->getSigmaNN() / (M_PI * 10.); // in fm^2
  double averageQs = 0.;
  double averageQs2 = 0.;
  double averageQs2Avg = 0.;
  double averageQs2min = 0.;
  double averageQs2min2 = 0.;
  int count = 0;
  double nucleiInAverage;
  nucleiInAverage = static_cast<double>(param->getAverageOverNuclei());

  param->setQsmuRatioB(param->getQsmuRatio());

  if (param->getUseNucleus() == 0) {
    if (param->getUseGaussian() == 1) {
      double sigmax = 0.35;
      double sigmay = 0.5;
      for (int ix = 0; ix < N; ix++) // loop over all positions
      {
        double x = ix * L / double(N) - L / 2.;
        for (int iy = 0; iy < N; iy++) {
          double y = iy * L / double(N) - L / 2.;
          int localpos = ix * N + iy;
          double envelope = exp(-(x * x / (2. * sigmax * sigmax) +
                                  y * y / (2. * sigmay * sigmay))) /
                            (2. * M_PI * sigmax * sigmay);
          lat->cells[localpos]->setg2mu2A(envelope * param->getg2mu() *
                                          param->getg2mu() / param->getg() /
                                          param->getg());
          lat->cells[localpos]->setg2mu2B(envelope * param->getg2mu() *
                                          param->getg2mu() / param->getg() /
                                          param->getg());
        }
      }
    } else {
      for (int ix = 0; ix < N; ix++) // loop over all positions
      {
        for (int iy = 0; iy < N; iy++) {
          int localpos = ix * N + iy;
          lat->cells[localpos]->setg2mu2A(param->getg2mu() * param->getg2mu() /
                                          param->getg() / param->getg());
          lat->cells[localpos]->setg2mu2B(param->getg2mu() * param->getg2mu() /
                                          param->getg() / param->getg());
        }
      }
    }
    param->setSuccess(1);
    cout << "constant color charge density set" << endl;
    return;
  }

#pragma omp parallel for
  for (int ix = 0; ix < N; ix++) // loop over all positions
  {
    for (int iy = 0; iy < N; iy++) {
      int localpos = ix * N + iy;
      lat->cells[localpos]->setg2mu2A(0.);
      lat->cells[localpos]->setg2mu2B(0.);
    }
  }

  // compute N_part
  // positions are shifted here. not later as in previous versions. bshift below
  // (in init(..)) is zero.

  double phiRP = 0.;
  if (param->getRotateReactionPlane()) {
    phiRP = 2 * M_PI * random->genrand64_real2();
  }
  for (unsigned int i = 0; i < nucleusA_.size(); i++) {
    // shift the nuclei's position by -b/2 or +b/2 respectively
    nucleusA_.at(i).x -= impact_b/2.*cos(phiRP);
    nucleusA_.at(i).y -= impact_b/2.*sin(phiRP);
  }
  for (unsigned int i = 0; i < nucleusB_.size(); i++) {
    // shift the nuclei's position by -b/2 or +b/2 respectively
    nucleusB_.at(i).x += impact_b/2.*cos(phiRP);
    nucleusB_.at(i).y += impact_b/2.*sin(phiRP);
  }

  double xi = param->getProtonAnisotropy();

  if (xi != 0.) {
    for (int i = 0; i < A1; i++) {
      nucleusA_.at(i).phi = 2 * M_PI * random->genrand64_real2();
    }

    for (int i = 0; i < A2; i++) {
      nucleusB_.at(i).phi = 2 * M_PI * random->genrand64_real2();
    }
  } else {
    for (int i = 0; i < A1; i++) {
      nucleusA_.at(i).phi = 0.;
    }

    for (int i = 0; i < A2; i++) {
      nucleusB_.at(i).phi = 0.;
    }
  }

  const int NqFlag = param->getUseConstituentQuarkProton();
  vector< vector<double> > xq1, xq2, yq1, yq2, BGq1, BGq2, gauss1, gauss2;
  vector<double> x_array, y_array, z_array, BGq_array, gauss_array;

  for (int i = 0; i < A1; i++) {
    x_array.clear();
    if (NqFlag > 0) {
      samplePartonPositions(param, random, x_array, y_array, z_array,
                            BGq_array);
      // if (param->getShiftConstituentQuarkProtonOrigin())
      // Move center of mass to the origin
      // Note that 1607.01711 this is not done, so parameters quoted in
      // that paper can't be used if this is done
      xq1.push_back(x_array);
      yq1.push_back(y_array);
      BGq1.push_back(BGq_array);
    }
    int Npartons = std::max(1, static_cast<int>(x_array.size()));
    sampleQsNormalization(random, param, Npartons, gauss_array);
    gauss1.push_back(gauss_array);
  }

  for (int i = 0; i < A2; i++) {
    x_array.clear();
    if (NqFlag > 0) {
      samplePartonPositions(param, random, x_array, y_array, z_array,
                            BGq_array);
      xq2.push_back(x_array);
      yq2.push_back(y_array);
      BGq2.push_back(BGq_array);
    }
    int Npartons = std::max(1, static_cast<int>(x_array.size()));
    sampleQsNormalization(random, param, Npartons, gauss_array);
    gauss2.push_back(gauss_array);
  }

  // test what a smmoth Woods-Saxon would give
  if (param->getUseSmoothNucleus() == 1) {
    cout << "Using smooth nucleus for test purposes. Does not include "
            "deformation."
         << endl;
    Npart = 2; // avoid break below
    double xA, xB;
    double y;
    double T;
    double localpos;
    double normA = 0.;
    double normB = 0.;
    double bb = param->getb();
    for (int ix = 0; ix < N; ix++) // loop over all positions
    {
      xA = -L / 2. + a * ix - bb / 2.;
      xB = -L / 2. + a * ix + bb / 2.;
      for (int iy = 0; iy < N; iy++) {
        y = -L / 2. + a * iy;

        localpos = ix * N + iy;

        // nucleus A
        r = sqrt(xA * xA + y * y);
        T = glauber->InterNuTInST(r);
        lat->cells[localpos]->setTpA(T);

        normA += T * a * a;

        // nucleus B
        r = sqrt(xB * xB + y * y);
        T = glauber->InterNuPInSP(r);
        lat->cells[localpos]->setTpB(T);

        normB += T * a * a;

        // remove potential stuff outside the interaction region
        if (lat->cells[localpos]->getTpA() < 0.001 ||
            lat->cells[localpos]->getTpB() < 0.001) {
          lat->cells[localpos]->setTpA(0.);
          lat->cells[localpos]->setTpB(0.);
        }
      }
    }
    for (int ix = 0; ix < N; ix++) // loop over all positions
    {
      for (int iy = 0; iy < N; iy++) {
        localpos = ix * N + iy;
        lat->cells[localpos]->setTpA(lat->cells[localpos]->getTpA() / normA *
                                     glauber->nucleusA1() * hbarc * hbarc);
        lat->cells[localpos]->setTpB(lat->cells[localpos]->getTpB() / normB *
                                     glauber->nucleusA2() * hbarc * hbarc);
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
  } else {
    // add all T_p's (new in version 1.2)
#pragma omp parallel
    {
      double x, xm;
      double y, ym;
      int localpos;
      double bp2, T, phi;

#pragma omp for
      for (int ix = 0; ix < N; ix++) // loop over all positions
      {
        x = -L / 2. + a * ix;
        for (int iy = 0; iy < N; iy++) {
          y = -L / 2. + a * iy;

          localpos = ix * N + iy;

          // nucleus A
          lat->cells[localpos]->setTpA(0.);
          for (int i = 0; i < A1; i++) {
            xm = nucleusA_.at(i).x;
            ym = nucleusA_.at(i).y;

            if (param->getUseConstituentQuarkProton() > 0) {
              T = 0.;
              for (unsigned int iq = 0; iq < xq1[i].size(); iq++) {
                bp2 = (xm + xq1[i][iq] - x) * (xm + xq1[i][iq] - x) +
                      (ym + yq1[i][iq] - y) * (ym + yq1[i][iq] - y);
                bp2 /= hbarc * hbarc;

                T += exp(-bp2 / (2. * BGq1[i][iq])) / (2. * M_PI * BGq1[i][iq]) /
                     (static_cast<double>(xq1[i].size())) *
                     gauss1[i][iq]; // I removed the 2/3 here to make it a bit
                                    // bigger
              }
            } else {
              const double BG = param->getBG();
              phi = nucleusA_.at(i).phi;

              bp2 = (xm - x) * (xm - x) + (ym - y) * (ym - y) +
                    xi * pow((xm - x) * cos(phi) + (ym - y) * sin(phi), 2.);
              bp2 /= hbarc * hbarc;
              T = sqrt(1 + xi) * exp(-bp2 / (2. * BG)) / (2. * M_PI * BG) *
                  gauss1[i][0]; // T_p in this cell for the current nucleon
            }
            lat->cells[localpos]->setTpA(lat->cells[localpos]->getTpA() +
                                         T / nucleiInAverage); // add up all T_p
          }

          // nucleus B
          lat->cells[localpos]->setTpB(0.);
          for (int i = 0; i < A2; i++) {
            xm = nucleusB_.at(i).x;
            ym = nucleusB_.at(i).y;

            if (param->getUseConstituentQuarkProton() > 0) {
              T = 0.;
              for (unsigned int iq = 0; iq < xq2[i].size(); iq++) {
                bp2 = (xm + xq2[i][iq] - x) * (xm + xq2[i][iq] - x) +
                      (ym + yq2[i][iq] - y) * (ym + yq2[i][iq] - y);
                bp2 /= hbarc * hbarc;

                T += exp(-bp2 / (2. * BGq2[i][iq])) / (2. * M_PI * BGq2[i][iq]) /
                     (static_cast<double>(xq2[i].size())) *
                     gauss2[i][iq];
              }
            } else {
              const double BG = param->getBG();
              phi = nucleusB_.at(i).phi;

              bp2 = (xm - x) * (xm - x) + (ym - y) * (ym - y) +
                    xi * pow((xm - x) * cos(phi) + (ym - y) * sin(phi), 2.);
              bp2 /= hbarc * hbarc;

              T = sqrt(1 + xi) * exp(-bp2 / (2. * BG)) / (2. * M_PI * BG) *
                  gauss2[i][0]; // T_p in this cell for the current nucleon
            }

            lat->cells[localpos]->setTpB(lat->cells[localpos]->getTpB() +
                                         T / nucleiInAverage); // add up all T_p
          }
        }
      }
    }
  }

  if (param->getUseSmoothNucleus() == 0) {
    stringstream strNcoll_name;
    strNcoll_name << "NcollList" << param->getEventId() << ".dat";
    string Ncoll_name;
    Ncoll_name = strNcoll_name.str();

    ofstream foutNcoll(Ncoll_name.c_str(), ios::out);

    if (param->getGaussianWounding() == 0) {
      for (int i = 0; i < A1; i++) {
        for (int j = 0; j < A2; j++) {
          dx = nucleusB_.at(j).x - nucleusA_.at(i).x;
          dy = nucleusB_.at(j).y - nucleusA_.at(i).y;
          dij = dx * dx + dy * dy;
          if (dij < d2) {
            foutNcoll << (nucleusB_.at(j).x + nucleusA_.at(i).x) / 2. << " "
                      << (nucleusB_.at(j).y + nucleusA_.at(i).y) / 2. << endl;
            Ncoll++;
            nucleusB_.at(j).collided = 1;
            nucleusA_.at(i).collided = 1;
          }
        }
      }
    } else {
      double p;
      double G = 0.92;
      double ran;

      for (int i = 0; i < A1; i++) {
        for (int j = 0; j < A2; j++) {
          dx = nucleusB_.at(j).x - nucleusA_.at(i).x;
          dy = nucleusB_.at(j).y - nucleusA_.at(i).y;
          dij = dx * dx + dy * dy;

          p = G * exp(-G * dij / d2); // Gaussian profile

          ran = random->genrand64_real1();

          if (ran < p) {
            foutNcoll << (nucleusB_.at(j).x + nucleusA_.at(i).x) / 2. << " "
                      << (nucleusB_.at(j).y + nucleusA_.at(i).y) / 2. << endl;
            Ncoll++;
            nucleusB_.at(j).collided = 1;
            nucleusA_.at(i).collided = 1;
          }
        }
      }
    }

    foutNcoll.close();

    stringstream strNpart_name;
    strNpart_name << "NpartList" << param->getEventId() << ".dat";
    string Npart_name;
    Npart_name = strNpart_name.str();

    ofstream foutNpart(Npart_name.c_str(), ios::out);

    for (int i = 0; i < A1; i++) {
      foutNpart << nucleusA_.at(i).x << " " << nucleusA_.at(i).y << " "
                << nucleusA_.at(i).proton << " " << nucleusA_.at(i).collided
                << endl;
    }
    foutNpart << endl;
    for (int i = 0; i < A2; i++) {
      foutNpart << nucleusB_.at(i).x << " " << nucleusB_.at(i).y << " "
                << nucleusB_.at(i).proton << " " << nucleusB_.at(i).collided
                << endl;
    }
    foutNpart.close();

    // in p+p assume that they collided in any case
    if (A1 == 1 && A2 == 1) {
      nucleusB_.at(0).collided = 1;
      nucleusA_.at(0).collided = 1;
    }

    Npart = 0;

    for (int i = 0; i < A1; i++) {
      if (nucleusA_.at(i).collided == 1)
        Npart++;
    }

    for (int i = 0; i < A2; i++) {
      if (nucleusB_.at(i).collided == 1)
        Npart++;
    }

    param->setNpart(Npart);

    if (param->getUseFixedNpart() != 0 && Npart != param->getUseFixedNpart()) {
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
    double localrapidity = rapidity;
    int check;
    QsA = 1;
    QsB = 1;

#pragma omp for
    for (int ix = 0; ix < N; ix++) // loop over all positions
    {
      // x = -L/2.+a*ix;
      for (int iy = 0; iy < N; iy++) {
        Ydeviation = 10000;
        check = 0;
        // y = -L/2.+a*iy;
        localpos = ix * N + iy;

        if (param->getUseSmoothNucleus() == 1)
          check = 2;
        else {
          const double BG = param->getBG();

          //    cut proton at a radius of rmax [fm] (about twice the gluonic
          //radius to be generous)

          if (log(2 * M_PI * BG * lat->cells[localpos]->getTpA()) < 0.)
            {
              if (isinf(log(2 * M_PI * BG * lat->cells[localpos]->getTpA()))==1)
                distanceA= param->getRmax()+1.;
              else 
                distanceA =
                  sqrt(-2. * BG *
                       log(2 * M_PI * BG * lat->cells[localpos]->getTpA())) *
                  hbarc;
                //cout << log(2 * M_PI * BG * lat->cells[localpos]->getTpA()) << endl;
            }
          else
            distanceA = 0.;

          if (log(2 * M_PI * BG * lat->cells[localpos]->getTpB()) < 0.){
            if (isinf(log(2 * M_PI * BG * lat->cells[localpos]->getTpB()))==1)
              distanceB= param->getRmax()+1.;
              else 
                distanceB =
                  sqrt(-2. * BG *
                       log(2 * M_PI * BG * lat->cells[localpos]->getTpB())) *
                  hbarc;
          }
          else
            distanceB = 0.;

          if (distanceA < param->getRmax()) {
            check = 1;
          }
          // else
          //   cout << "large: " << distanceA << " " << ix << " " << iy << endl;

          if (distanceB < param->getRmax() && check == 1) {
            check = 2;
          }
        }

        double exponent = 5.6; // see 1212.2974 Eq. (17)
        if (check == 2) {
          if (param->getUseFluctuatingx() == 1) {
            // iterative loops here to determine the fluctuating Y
            // _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
            while (abs(Ydeviation) > 0.001) {
              if (localrapidity >= 0) {
                QsA = sqrt(getNuclearQs2(lat->cells[localpos]->getTpA(),
                                         abs(localrapidity)));
              } else {
                xVal = QsA * param->getxFromThisFactorTimesQs() /
                       param->getRoots() * exp(yIn);
                if (xVal == 0)
                  QsA = 0.;
                else
                  QsA =
                      sqrt(getNuclearQs2(lat->cells[localpos]->getTpA(), 0.)) *
                      sqrt(pow((1 - xVal) / (1 - 0.01), exponent) *
                           pow((0.01 / xVal), 0.2));
              }
              if (QsA == 0) {
                Ydeviation = 0;
                lat->cells[localpos]->setg2mu2A(0.);
              } else {
                // nucleus A
                lat->cells[localpos]->setg2mu2A(
                    QsA * QsA / param->getQsmuRatio() / param->getQsmuRatio() *
                    a * a / hbarc / hbarc / param->getg() /
                    param->getg()); // lattice units? check

                Ydeviation =
                    localrapidity -
                    log(0.01 / (QsA * param->getxFromThisFactorTimesQs() /
                                param->getRoots() * exp(yIn)));
                localrapidity =
                    log(0.01 / (QsA * param->getxFromThisFactorTimesQs() /
                                param->getRoots() * exp(yIn)));
              }
            }
            if (lat->cells[localpos]->getg2mu2A() !=
                lat->cells[localpos]->getg2mu2A()) {
              lat->cells[localpos]->setg2mu2A(0.);
            }

            localrapidity = rapidity;
            Ydeviation = 10000;
            while (abs(Ydeviation) > 0.001) {
              if (localrapidity >= 0)
                QsB = sqrt(getNuclearQs2(lat->cells[localpos]->getTpB(),
                                         abs(localrapidity)));
              else {
                xVal = QsB * param->getxFromThisFactorTimesQs() /
                       param->getRoots() * exp(-yIn);
                if (xVal == 0)
                  QsB = 0.;
                else
                  QsB =
                      sqrt(getNuclearQs2(lat->cells[localpos]->getTpB(), 0.)) *
                      sqrt(pow((1 - xVal) / (1 - 0.01), exponent) *
                           pow((0.01 / xVal), 0.2));
              }
              if (QsB == 0) {
                Ydeviation = 0;
                lat->cells[localpos]->setg2mu2B(0.);
              } else {
                // nucleus B
                lat->cells[localpos]->setg2mu2B(
                    QsB * QsB / param->getQsmuRatioB() /
                    param->getQsmuRatioB() * a * a / hbarc / hbarc /
                    param->getg() / param->getg());
                Ydeviation =
                    localrapidity -
                    log(0.01 / (QsB * param->getxFromThisFactorTimesQs() /
                                param->getRoots() * exp(-yIn)));
                localrapidity =
                    log(0.01 / (QsB * param->getxFromThisFactorTimesQs() /
                                param->getRoots() * exp(-yIn)));
              }
            }
            if (lat->cells[localpos]->getg2mu2B() !=
                lat->cells[localpos]->getg2mu2B()) {
              lat->cells[localpos]->setg2mu2B(0.);
            }
            // _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
            // end iterative loops here
          } else {
            // nucleus A
            lat->cells[localpos]->setg2mu2A(
                getNuclearQs2(lat->cells[localpos]->getTpA(), localrapidity) /
                param->getQsmuRatio() / param->getQsmuRatio() * a * a / hbarc /
                hbarc / param->getg() / param->getg()); // lattice units? check

            // nucleus B
            lat->cells[localpos]->setg2mu2B(
                getNuclearQs2(lat->cells[localpos]->getTpB(), localrapidity) /
                param->getQsmuRatioB() / param->getQsmuRatioB() * a * a /
                hbarc / hbarc / param->getg() / param->getg());
          }
        }
      }
    }
  }

  count = 0;
  double Tpp = 0.;
  double x, xm, y, ym;
  double alphas = 0.;
  int check = 0;
  for (int ix = 0; ix < N; ix++) // loop over all positions
  {
    for (int iy = 0; iy < N; iy++) {
      check = 0;
      pos = ix * N + iy;
      x = -L / 2. + a * ix;
      y = -L / 2. + a * iy;
      //    outvalue = lat->cells[pos]->getg2mu2A();

      posA = pos;
      posB = pos;

      if (posA > 0 && posA < (N - 1) * N + N - 1) {
        g2mu2A = lat->cells[posA]->getg2mu2A();
      } else
        g2mu2A = 0;

      if (posB > 0 && posB < (N - 1) * N + N - 1) {
        g2mu2B = lat->cells[posB]->getg2mu2B();
      } else
        g2mu2B = 0;

      if (g2mu2B >= g2mu2A) {
        averageQs2min2 += g2mu2A * param->getQsmuRatio() *
                          param->getQsmuRatio() / a / a * hbarc * hbarc *
                          param->getg() * param->getg();
      } else {
        averageQs2min2 += g2mu2B * param->getQsmuRatioB() *
                          param->getQsmuRatioB() / a / a * hbarc * hbarc *
                          param->getg() * param->getg();
      }

      for (int i = 0; i < A1; i++) {
        xm = nucleusA_.at(i).x;
        ym = nucleusA_.at(i).y;
        r = sqrt((x - xm) * (x - xm) + (y - ym) * (y - ym));

        if (r < sqrt(0.1 * param->getSigmaNN() / M_PI) &&
            nucleusA_.at(i).collided == 1) {
          check = 1;
        }
      }

      for (int i = 0; i < A2; i++) {
        xm = nucleusB_.at(i).x;
        ym = nucleusB_.at(i).y;
        r = sqrt((x - xm) * (x - xm) + (y - ym) * (y - ym));

        if (r < sqrt(0.1 * param->getSigmaNN() / M_PI) &&
            nucleusB_.at(i).collided == 1 && check == 1)
          check = 2;
      }

      if (check == 2) {
        if (g2mu2B > g2mu2A) {
          averageQs +=
              sqrt(g2mu2B * param->getQsmuRatioB() * param->getQsmuRatioB() /
                   a / a * hbarc * hbarc * param->getg() * param->getg());
          averageQs2 += g2mu2B * param->getQsmuRatioB() *
                        param->getQsmuRatioB() / a / a * hbarc * hbarc *
                        param->getg() * param->getg();
          averageQs2min += g2mu2A * param->getQsmuRatio() *
                           param->getQsmuRatio() / a / a * hbarc * hbarc *
                           param->getg() * param->getg();
        } else {
          averageQs +=
              sqrt(g2mu2A * param->getQsmuRatio() * param->getQsmuRatio() / a /
                   a * hbarc * hbarc * param->getg() * param->getg());
          averageQs2 += g2mu2A * param->getQsmuRatio() * param->getQsmuRatio() /
                        a / a * hbarc * hbarc * param->getg() * param->getg();
          averageQs2min += g2mu2B * param->getQsmuRatioB() *
                           param->getQsmuRatioB() / a / a * hbarc * hbarc *
                           param->getg() * param->getg();
        }
        averageQs2Avg +=
            (g2mu2B * param->getQsmuRatioB() * param->getQsmuRatioB() +
             g2mu2A * param->getQsmuRatio() * param->getQsmuRatio()) /
            2. / a / a * hbarc * hbarc * param->getg() * param->getg();
        count++;
      }
      // compute T_pp
      Tpp += lat->cells[pos]->getTpB() * lat->cells[pos]->getTpA() * a * a /
             hbarc / hbarc / hbarc / hbarc; // now this quantity is in fm^-2
      // remember: Tp is in GeV^2
    }
  }

  averageQs /= static_cast<double>(count);
  averageQs2 /= static_cast<double>(count);
  averageQs2Avg /= static_cast<double>(count);
  averageQs2min /= static_cast<double>(count);

  param->setAverageQs(sqrt(averageQs2));
  param->setAverageQsAvg(sqrt(averageQs2Avg));
  param->setAverageQsmin(sqrt(averageQs2min));

  param->setTpp(Tpp);

  messager << "N_part=" << Npart;
  messager.flush("info");
  messager << "N_coll=" << Ncoll;
  messager.flush("info");
  cout << "T_pp(" << param->getb() << " fm) = " << Tpp << " 1/fm^2" << endl;
  cout << "Q_s^2(max) S_T = "
       << averageQs2 * a * a / hbarc / hbarc * static_cast<double>(count)
       << endl;
  cout << "Q_s^2(avg) S_T = "
       << averageQs2Avg * a * a / hbarc / hbarc * static_cast<double>(count)
       << endl;
  cout << "Q_s^2(min) S_T = " << averageQs2min * a * a / hbarc / hbarc * static_cast<double>(count) << endl;
  cout << "Q_s^2(min) S_T = " << averageQs2min2 * a * a / hbarc / hbarc << endl;

  cout << "Area = " << a * a * count << " fm^2" << endl;

  cout << "Average Qs(max) = " << param->getAverageQs() << " GeV" << endl;
  cout << "Average Qs(avg) = " << param->getAverageQsAvg() << " GeV" << endl;
  cout << "Average Qs(min) = " << param->getAverageQsmin() << " GeV" << endl;

  cout << "resulting Y(Qs(max)*" << param->getxFromThisFactorTimesQs() << ") = "
       << log(0.01 / (param->getAverageQs() *
                      param->getxFromThisFactorTimesQs() / param->getRoots()))
       << endl;
  cout << "resulting Y(Qs(avg)*" << param->getxFromThisFactorTimesQs() << ") = "
       << log(0.01 / (param->getAverageQsAvg() *
                      param->getxFromThisFactorTimesQs() / param->getRoots()))
       << endl;
  cout << "resulting Y(Qs(min)*" << param->getxFromThisFactorTimesQs()
       << ") =  "
       << log(0.01 / (param->getAverageQsmin() *
                      param->getxFromThisFactorTimesQs() / param->getRoots()))
       << endl;

  messager.info("Color charge densities for nucleus A and B set. ");

  if (param->getRunningCoupling() && param->getRunWithkt() == 0) {
    if (param->getRunWithQs() == 2) {
      cout << "running with " << param->getRunWithThisFactorTimesQs()
           << " Q_s(max)" << endl;
      alphas = 12. * M_PI /
               ((27.) * 2. *
                log(param->getRunWithThisFactorTimesQs() *
                    param->getAverageQs() / 0.2)); // 3 flavors
      cout << "alpha_s(" << param->getRunWithThisFactorTimesQs()
           << " Qs_max)=" << alphas << endl;
    } else if (param->getRunWithQs() == 0) {
      cout << "running with " << param->getRunWithThisFactorTimesQs()
           << " Q_s(min)" << endl;
      alphas = 12. * M_PI /
               ((27.) * 2. *
                log(param->getRunWithThisFactorTimesQs() *
                    param->getAverageQsmin() / 0.2)); // 3 flavors
      cout << "alpha_s(" << param->getRunWithThisFactorTimesQs()
           << " Qs_min)=" << alphas << endl;
    } else if (param->getRunWithQs() == 1) {
      cout << "running with " << param->getRunWithThisFactorTimesQs()
           << " <Q_s>" << endl;
      alphas = 12. * M_PI /
               ((27.) * 2. *
                log(param->getRunWithThisFactorTimesQs() *
                    param->getAverageQsAvg() / 0.2)); // 3 flavors
      cout << "alpha_s(" << param->getRunWithThisFactorTimesQs()
           << " <Qs>)=" << alphas << endl;
    }
  } else if (param->getRunningCoupling() && param->getRunWithkt() == 1) {
    cout << "Multiplicity with running alpha_s(k_T)" << endl;
  } else {
    cout << "Using fixed alpha_s" << endl;
    alphas = param->getg() * param->getg() / 4. / M_PI;
  }

  if (param->getAverageQs() > 0 && param->getAverageQsAvg() > 0 &&
      averageQs2 > 0 && param->getAverageQsmin() > 0 && averageQs2Avg > 0 &&
      alphas > 0 && Npart >= 2 && averageQs2min2 * a * a / hbarc / hbarc > param->getMinimumQs2ST())
    {
      param->setSuccess(1);
      stringstream strup_name;
      strup_name << "usedParameters" << param->getEventId() << ".dat";
      string up_name;
      up_name = strup_name.str();

      ofstream fout1(up_name.c_str(), ios::app);
      fout1 << " " << endl;
      fout1 << " Output by setColorChargeDensity in Init.cpp: " << endl;
      fout1 << " " << endl;
      fout1 << "b = " << impact_b << " fm" << endl;
      fout1 << "phiRP = " << phiRP << endl;
      fout1 << "Npart = " << Npart << endl;
      fout1 << "Ncoll = " << Ncoll << endl;
      if (param->getRunningCoupling()) {
        if (param->getRunWithQs() == 2)
          fout1 << "<Q_s>(max) = " << param->getAverageQs() << endl;
        else if (param->getRunWithQs() == 1)
          fout1 << "<Q_s>(avg) = " << param->getAverageQsAvg() << endl;
        else if (param->getRunWithQs() == 0)
          fout1 << "<Q_s>(min) = " << param->getAverageQsmin() << endl;
        fout1 << "alpha_s(" << param->getRunWithThisFactorTimesQs()
              << " <Q_s>) = " << param->getalphas() << endl;
      } else
        fout1 << "using fixed coupling alpha_s=" << param->getalphas() << endl;
      fout1.close();
    }
  if ( averageQs2min2 * a * a / hbarc / hbarc < param->getMinimumQs2ST()) 
    cout << " **** Rejected event - Qsmin^2 S_T=" << averageQs2min2 * a * a / hbarc / hbarc << " too small ( < " << param->getMinimumQs2ST() << ")." << endl;


  param->setalphas(alphas);

  stringstream strNEst_name;
  strNEst_name << "NgluonEstimators" << param->getEventId() << ".dat";
  string NEst_name;
  NEst_name = strNEst_name.str();

  ofstream foutNEst(NEst_name.c_str(), ios::out);

  foutNEst << "#Q_s^2(min) S_T  " <<  "Q_s^2(avg) S_T  " << "Q_s^2(max) S_T " << " Q_s^2(min) S_T Log^2( Q_s^2(max) / Q_s^2(min))  " << endl;
  foutNEst << averageQs2min2 * a * a / hbarc / hbarc  <<  "         " << averageQs2Avg * a * a / hbarc / hbarc * static_cast<double>(count) << "         " << averageQs2 * a * a / hbarc / hbarc * static_cast<double>(count) << "         " << averageQs2min2 * a * a / hbarc / hbarc * pow(log(averageQs2* static_cast<double>(count)/averageQs2min2),2.) << endl;

  foutNEst.close();
}

void Init::setV(Lattice *lat, Parameters *param, Random *random) {
  messager.info("Setting Wilson lines ...");
  const int N = param->getSize();
  const int Ny = param->getNy();
  const int nn[2] = {N, N};
  const double L = param->getL();
  const double a = L / N; // lattice spacing in fm
  const double m = param->getm() * a / hbarc;
  double UVdamp = param->getUVdamp(); // GeV^-1
  UVdamp = UVdamp / a * hbarc;
  complex<double> **rhoACoeff;
  rhoACoeff = new complex<double> *[Nc2m1_];
  for (int i = 0; i < Nc2m1_; i++) {
    rhoACoeff[i] = new complex<double>[N * N];
  }

  // loop over longitudinal direction
  for (int k = 0; k < Ny; k++) {
    double g2muA;
    for (int pos = 0; pos < N * N; pos++) {
      for (int n = 0; n < Nc2m1_; n++) {
        g2muA = param->getg() *
                sqrt(lat->cells[pos]->getg2mu2A() / static_cast<double>(Ny));
        rhoACoeff[n][pos] = g2muA * random->Gauss();
      }
    }

    for (int n = 0; n < Nc2m1_; n++) {
      fft.fftnComplex(rhoACoeff[n], rhoACoeff[n], nn, 1);
    }

    // compute A^+
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        double kt2, kx, ky;
        int localpos = i * N + j;
        kx = 2. * M_PI *
             (-0.5 + static_cast<double>(i) / static_cast<double>(N));
        ky = 2. * M_PI *
             (-0.5 + static_cast<double>(j) / static_cast<double>(N));
        kt2 = 4. * (sin(kx / 2.) * sin(kx / 2.) +
                    sin(ky / 2.) * sin(ky / 2.)); // lattice momentum
        if (m == 0) {
          if (kt2 != 0) {
            for (int n = 0; n < Nc2m1_; n++) {
              rhoACoeff[n][localpos] = rhoACoeff[n][localpos] * (1. / (kt2));
            }
          } else {
            for (int n = 0; n < Nc2m1_; n++) {
              rhoACoeff[n][localpos] = 0.;
            }
          }
        } else {
          for (int n = 0; n < Nc2m1_; n++) {
            rhoACoeff[n][localpos] *=
                (1. / (kt2 + m * m)) * exp(-sqrt(kt2) * UVdamp);
          }
        }
      }
    }

    // Fourier transform back A^+
    for (int n = 0; n < Nc2m1_; n++) {
      fft.fftnComplex(rhoACoeff[n], rhoACoeff[n], nn, -1);
    }
    // compute U

#pragma omp parallel
    {
      std::vector<double> in(Nc2m1_, 0.);
      Matrix temp(Nc_, 1.);
      Matrix tempNew(Nc_, 0.);

#pragma omp for
      for (int pos = 0; pos < N * N; pos++) {
        for (int aa = 0; aa < Nc2m1_; aa++) {
          // expmCoeff will calculate exp(i in[a]t[a]),
          // so just multiply by -1 (not -i)
          in[aa] = -(rhoACoeff[aa][pos]).real();
        }
        tempNew = getUfromExponent(in);
        temp = tempNew * lat->cells[pos]->getU();
        // set U
        lat->cells[pos]->setU(temp);
      }
    }

  } // Ny loop

  // loop over longitudinal direction
  for (int k = 0; k < Ny; k++) {
    double g2muB;
    for (int pos = 0; pos < N * N; pos++) {
      for (int n = 0; n < Nc2m1_; n++) {
        g2muB = param->getg() *
                sqrt(lat->cells[pos]->getg2mu2B() / static_cast<double>(Ny));
        rhoACoeff[n][pos] = g2muB * random->Gauss();
      }
    }

    for (int n = 0; n < Nc2m1_; n++) {
      fft.fftnComplex(rhoACoeff[n], rhoACoeff[n], nn, 1);
    }

    // compute A^+
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        double kt2, kx, ky;
        int localpos = i * N + j;
        kx = 2. * M_PI *
             (-0.5 + static_cast<double>(i) / static_cast<double>(N));
        ky = 2. * M_PI *
             (-0.5 + static_cast<double>(j) / static_cast<double>(N));
        kt2 = 4. * (sin(kx / 2.) * sin(kx / 2.) +
                    sin(ky / 2.) * sin(ky / 2.)); // lattice momentum
        if (m == 0) {
          if (kt2 != 0) {
            for (int n = 0; n < Nc2m1_; n++) {
              rhoACoeff[n][localpos] = rhoACoeff[n][localpos] * (1. / (kt2));
            }
          } else {
            for (int n = 0; n < Nc2m1_; n++) {
              rhoACoeff[n][localpos] = 0.;
            }
          }
        } else {
          for (int n = 0; n < Nc2m1_; n++) {
            rhoACoeff[n][localpos] *=
                (1. / (kt2 + m * m)) * exp(-sqrt(kt2) * UVdamp);
          }
        }
      }
    }

    // Fourier transform back A^+
    for (int n = 0; n < Nc2m1_; n++) {
      fft.fftnComplex(rhoACoeff[n], rhoACoeff[n], nn, -1);
    }
    // compute U

    //#pragma omp parallel
    {
      std::vector<double> in(Nc2m1_, 0.);
      Matrix temp(Nc_, 1.);
      Matrix tempNew(Nc_, 0.);

#pragma omp for
      for (int pos = 0; pos < N * N; pos++) {

        for (int aa = 0; aa < Nc2m1_; aa++) {
          // expmCoeff will calculate exp(i in[a]t[a]), so
          // just multiply by -1 (not -i)
          in[aa] = -(rhoACoeff[aa][pos]).real();
        }
        tempNew = getUfromExponent(in);
        temp = tempNew * lat->cells[pos]->getU2();

        // set U
        lat->cells[pos]->setU2(temp);
      }
    }

  } // Ny loop

  // --------
  for (int ic = 0; ic < Nc2m1_; ic++) {
    delete[] rhoACoeff[ic];
  }
  delete[] rhoACoeff;


  // // output U
  if (param->getWriteInitialWilsonLines() > 0) {
      if (std::abs(param->getb()) > 1e-5)
          messager.warning("Writing Wilson lines with a non-zero impact parameter b!");

    stringstream strVOne_name;
    // strVOne_name << "V1-" << param->getMPIRank() << ".txt";
    strVOne_name << "V-"
                 << param->getEventId() +
                        2 * param->getSeed() * param->getMPISize();
    if (param->getWriteInitialWilsonLines() == 1)
      strVOne_name << ".txt";
    string VOne_name;
    VOne_name = strVOne_name.str();

    stringstream strVTwo_name;
    // strVTwo_name << "V2-" << param->getMPIRank() << ".txt";
    strVTwo_name << "V-"
                 << param->getEventId() +
                        (1 + 2 * param->getSeed()) * param->getMPISize();
    if (param->getWriteInitialWilsonLines() == 1)
      strVTwo_name << ".txt";
    string VTwo_name;
    VTwo_name = strVTwo_name.str();
    // Output in text
    if (param->getWriteInitialWilsonLines() == 1) {
      ofstream foutU(VOne_name.c_str(), ios::out);
      foutU.precision(15);

      for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) // loop over all positions
        {
          int pos = ix * N + iy;
          foutU << ix << " " << iy << " "
                << (lat->cells[pos]->getU()).MatrixToString() << endl;
        }
        foutU << endl;
      }
      foutU.close();

      cout << "wrote " << strVOne_name.str() << endl;

      ofstream foutU2(VTwo_name.c_str(), ios::out);
      foutU2.precision(15);
      for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) // loop over all positions
        {
          int pos = ix * N + iy;
          foutU2 << ix << " " << iy << " "
                 << (lat->cells[pos]->getU2()).MatrixToString() << endl;
        }
        foutU2 << endl;
      }
      foutU2.close();

      cout << "wrote " << strVTwo_name.str() << endl;
    } // end output in text
    else if (param->getWriteInitialWilsonLines() == 2) {
      std::ofstream Outfile1, Outfile2;
      Outfile1.open(VOne_name.c_str(), ios::out | ios::binary);
      Outfile2.open(VTwo_name.c_str(), ios::out | ios::binary);

      double temp = param->getRapidity();

      // print header ------------- //
      Outfile1.write((char *)&N, sizeof(int));
      Outfile1.write((char *)&Nc_, sizeof(int));
      Outfile1.write((char *)&L, sizeof(double));
      Outfile1.write((char *)&a, sizeof(double));
      Outfile1.write((char *)&temp, sizeof(double));

      Outfile2.write((char *)&N, sizeof(int));
      Outfile2.write((char *)&Nc_, sizeof(int));
      Outfile2.write((char *)&L, sizeof(double));
      Outfile2.write((char *)&a, sizeof(double));
      Outfile2.write((char *)&temp, sizeof(double));
      //

      double *val1 = new double[2];
      double *val2 = new double[2];

      for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
          for (int a1 = 0; a1 < 3; a1++) {
            for (int b = 0; b < 3; b++) {
              int indx = N * iy + ix;
              val1[0] = (lat->cells[indx]->getU()).getRe(a1 * Nc_ + b);
              val1[1] = (lat->cells[indx]->getU()).getIm(a1 * Nc_ + b);
              val2[0] = (lat->cells[indx]->getU2()).getRe(a1 * Nc_ + b);
              val2[1] = (lat->cells[indx]->getU2()).getIm(a1 * Nc_ + b);

              Outfile1.write((char *)val1, 2 * sizeof(double));
              Outfile2.write((char *)val2, 2 * sizeof(double));
            }
          }
        }
      }

      if (Outfile1.good() == false || Outfile2.good() == false) {
        std::cerr << "#CRTICAL ERROR -- BINARY OUTPUT OF VECTOR CURRENTS FAILED"
                  << std::endl;
        exit(1);
      }

      delete[] val1;
      delete[] val2;

      Outfile1.close();
      Outfile2.close();
      cout << "wrote " << strVOne_name.str() << " and " << strVTwo_name.str()
           << endl;
    } // end binary output
    else {
      std::cerr << "# Unknwon option param->getWriteInitialWilsonLines()=="
                << param->getWriteInitialWilsonLines() << std::endl;
      exit(1);
    }
  }
  // --------

  messager << " Wilson lines V_A and V_B set on rank " << param->getMPIRank()
           << ". ";
  messager.flush("info");
}

void Init::readV(Lattice *lat, Parameters *param, int format) {
  // format 1 = plain text, 2 = binary

  if (format > 2 or format < 1)
  {
    messager << "Unknown format " << format << " when reading the initial Wilson lines, supported formats: 1,2";
    messager.flush("info");
    exit(1);
  }

  stringstream strVOne_name;
  // strVOne_name << "V1-" << param->getMPIRank() << ".txt";
  strVOne_name << "V-"
               << param->getEventId() +
                      2 * param->getSeed() * param->getMPISize();
  if (format == 1)
    strVOne_name << ".txt";
  string VOne_name;
  VOne_name = strVOne_name.str();

  stringstream strVTwo_name;
  // strVTwo_name << "V2-" << param->getMPIRank() << ".txt";
  strVTwo_name << "V-"
                 << param->getEventId() +
                        (1 + 2 * param->getSeed()) * param->getMPISize();
  if (format == 1)
    strVTwo_name << ".txt";
  string VTwo_name;
  VTwo_name = strVTwo_name.str();

  messager << "Reading Wilson lines from files " << VOne_name << " and " << VTwo_name;
  messager.flush("info");

  if (format == 1)
  {
    int N = param->getSize();

    double L = param->getL();
    double a = L/static_cast<double>(N);

    int nn[2];
    nn[0] = N;
    nn[1] = N;

    Matrix temp(Nc_, 1.);

    double Re[9], Im[9];
    double dummy;

    ifstream finV1(VOne_name.c_str(), ios::in);

    if (!finV1) {
      messager << "File " << VOne_name << " not found. Exiting.";
      messager.flush("info");
      exit(1);
    }

    messager << "Reading Wilson line from file " << VOne_name << " ..." ;

    // set V for nucleus A
    for (int i = 0; i < nn[0]; i++) {
      for (int j = 0; j < nn[1]; j++) {

        finV1 >> dummy >> dummy >> Re[0] >> Im[0] >> Re[1] >> Im[1] >> Re[2] >>
            Im[2] >> Re[3] >> Im[3] >> Re[4] >> Im[4] >> Re[5] >> Im[5] >>
            Re[6] >> Im[6] >> Re[7] >> Im[7] >> Re[8] >> Im[8];

        temp.set(0, 0, complex<double>(Re[0], Im[0]));
        temp.set(0, 1, complex<double>(Re[1], Im[1]));
        temp.set(0, 2, complex<double>(Re[2], Im[2]));
        temp.set(1, 0, complex<double>(Re[3], Im[3]));
        temp.set(1, 1, complex<double>(Re[4], Im[4]));
        temp.set(1, 2, complex<double>(Re[5], Im[5]));
        temp.set(2, 0, complex<double>(Re[6], Im[6]));
        temp.set(2, 1, complex<double>(Re[7], Im[7]));
        temp.set(2, 2, complex<double>(Re[8], Im[8]));

        double bb = param->getb();
        a = L/static_cast<double>(N);

        double xtemp = a * i - bb / 2.;
        int ix = xtemp / a;

        if (ix<0)
          continue;

        int pos = ix * N + j;
        lat->cells[pos]->setU(temp);
      }
    }

    finV1.close();


    ifstream finV2(VTwo_name.c_str(), ios::in);

    if (!finV2) {
      cerr << "File " << VTwo_name << " not found. Exiting." << endl;
      exit(1);
    }

    cout << "Reading Wilson line from file " << VTwo_name << " ..." << endl;

    // set V for nucleus B
    for (int i = 0; i < nn[0]; i++) {
      for (int j = 0; j < nn[1]; j++) {

        finV2 >> dummy >> dummy >> Re[0] >> Im[0] >> Re[1] >> Im[1] >> Re[2] >>
            Im[2] >> Re[3] >> Im[3] >> Re[4] >> Im[4] >> Re[5] >> Im[5] >>
            Re[6] >> Im[6] >> Re[7] >> Im[7] >> Re[8] >> Im[8];

        temp.set(0, 0, complex<double>(Re[0], Im[0]));
        temp.set(0, 1, complex<double>(Re[1], Im[1]));
        temp.set(0, 2, complex<double>(Re[2], Im[2]));
        temp.set(1, 0, complex<double>(Re[3], Im[3]));
        temp.set(1, 1, complex<double>(Re[4], Im[4]));
        temp.set(1, 2, complex<double>(Re[5], Im[5]));
        temp.set(2, 0, complex<double>(Re[6], Im[6]));
        temp.set(2, 1, complex<double>(Re[7], Im[7]));
        temp.set(2, 2, complex<double>(Re[8], Im[8]));

        double bb = param->getb();
        a = L/static_cast<double>(N);

        double xtemp = a * i + bb / 2.;
        int ix = xtemp / a;

        if (ix>=N)
          continue;

        int pos = ix * N + j;
        lat->cells[pos]->setU2(temp);
      }
    }

    finV2.close();
  }
  else if (format == 2)
  {

    std::ifstream InStream;
    InStream.precision(15);
    InStream.open(VOne_name.c_str(), std::ios::in | std::ios::binary);
    int N;
    double L, a, temp;

    Matrix tempM(Nc_, 1.);

    if (!InStream.good())
    {
        messager << "File " << VOne_name.c_str() << " does not exist!";
        messager.flush("info");
        exit(1);
    }

    if(InStream.is_open())
    {
        // READING IN PARAMETERS
        InStream.read(reinterpret_cast<char*>(&N), sizeof(int));
        InStream.read(reinterpret_cast<char*>(&Nc_), sizeof(int));
        InStream.read(reinterpret_cast<char*>(&L), sizeof(double));
        InStream.read(reinterpret_cast<char*>(&a), sizeof(double));
        InStream.read(reinterpret_cast<char*>(&temp), sizeof(double));

        if(N != param->getSize())
        {
          messager << "# ERROR wrong lattice size, data is " << N
          << " but you have specified " << param->getSize() ;
           exit(0);
        }
      if(std::abs(L - param->getL()) > 1e-5)
      {
        messager << "# ERROR grid length, data has " << L
        << " but you have specified " << param->getL();
         exit(0);
      }


        // READING ACTUAL DATA
        double ValueBuffer;
        int INPUT_CTR=0;
        double re,im;
        re = 0.;
        im = 0.;

        while( InStream.read(reinterpret_cast<char*>(&ValueBuffer), sizeof(double)))
        {
            if(INPUT_CTR%2==0)              //this is the real part
            {
                 re=ValueBuffer;
            }
            else                            // this is the imaginary part, write then to variable //
            {
                im=ValueBuffer;
                int TEMPINDX=((INPUT_CTR-1)/2);
                int PositionIndx = TEMPINDX / 9;

                // shift here by half an impact parameter
                int iy = PositionIndx / N;
                int ixIn = PositionIndx - N*iy;

                double bb = param->getb();
                a = L/static_cast<double>(N);

                double xtemp = a * ixIn - bb / 2.;

                int ix = round(xtemp / a);

                // cout << ixIn << " " << ix << endl;

                int MatrixIndx=TEMPINDX - PositionIndx*9;
                int j=MatrixIndx/3;
                int k=MatrixIndx-j*3;

                int indx = N*ix + iy;
                if (indx >= N*N || indx < 0)
              {
                if (bb==0){
                  messager << "Warning: datafile " << VOne_name << " has an element " << indx << " (iy=" << iy << ", ix=" << ix << "), but the grid is N="<< N << ". Element is (" << re <<" + " << im << "i), skipping it";
                  messager.flush("info");
                }
                INPUT_CTR++;
                continue;

              }
                lat->cells[indx]->getU().set(j,k, complex<double> (re,im));

            }
            INPUT_CTR++;
        }

      InStream.close();

      std::ifstream InStream2;
      InStream2.precision(15);
      InStream2.open(VTwo_name.c_str(), std::ios::in | std::ios::binary);
      if (!InStream2.good())
      {
          messager << "File " << VTwo_name.c_str() << " does not exist!";
          messager.flush("info");
          exit(1);
      }

      INPUT_CTR=0;
      if(InStream2.is_open())
      {
          // READING IN PARAMETERS
          InStream2.read(reinterpret_cast<char*>(&N), sizeof(int));
          InStream2.read(reinterpret_cast<char*>(&Nc_), sizeof(int));
          InStream2.read(reinterpret_cast<char*>(&L), sizeof(double));
          InStream2.read(reinterpret_cast<char*>(&a), sizeof(double));
          InStream2.read(reinterpret_cast<char*>(&temp), sizeof(double));


          if(N != param->getSize())
          {
            messager << "# ERROR wrong lattice size, data is " << N
            << " but you have specified " << param->getSize();
            messager.flush("info");
             exit(0);
          }
        if(std::abs(L - param->getL()) > 1e-5)
        {
          messager << "# ERROR grid length, dataas " << L
          << " but you have specified " << param->getL() ;
          messager.flush("info");
           exit(0);
        }


          // READING ACTUAL DATA
          while( InStream2.read(reinterpret_cast<char*>(&ValueBuffer), sizeof(double)))
          {
              if(INPUT_CTR%2==0)              //this is the real part
              {
                   re=ValueBuffer;
              }
              else                            // this is the imaginary part, write then to variable //
              {
                  im=ValueBuffer;

                  int TEMPINDX=((INPUT_CTR-1)/2);
                  int PositionIndx = TEMPINDX / 9;

                  // shift here by half an impact parameter
                  int iy = PositionIndx / N;
                  int ixIn = PositionIndx - N*iy;

                  double bb = param->getb();
                  a = L/static_cast<double>(N);

                  double xtemp = a * ixIn + bb / 2.;

                  int ix = round(xtemp / a);

                  //                  if (ixIn != ix)
                  //  cout << ixIn << " " << ix << endl;
                  int MatrixIndx=TEMPINDX - PositionIndx*9;
                  int j=MatrixIndx/3;
                  int k=MatrixIndx-j*3;

                  int indx = N*ix + iy;

                  if (indx >= N*N || indx < 0)
                  {
                    if (bb==0){
                      messager << "Warning: datafile " << VTwo_name << " has an element " << indx << " (iy=" << iy << ", ix=" << ix << "), but the grid is N="<< N << ". Element is (" << re <<" + " << im << "i), skipping it";
                      messager.flush("info");
                    }
                    INPUT_CTR++;
                    continue;
                  }
                  lat->cells[indx]->getU2().set(j,k, complex<double> (re,im));
               // if (indx > 65000) cout << "Save ok" << endl;

              }
              INPUT_CTR++;
        }
        InStream2.close();
      }
    }
  }

  messager << " Wilson lines V_A and V_B set on rank " << param->getMPIRank()
           << ". ";
  messager.flush("info");
}

void Init::init(Lattice *lat, Group *group, Parameters *param, Random *random,
                Glauber *glauber, int READFROMFILE) {
  const int maxIterations = 100000;
  const int N = param->getSize();
  Nc_ = param->getNc();
  Nc2m1_ = Nc_ * Nc_ - 1;
  group_ptr_ = group;
  random_ptr_ = random;
  one_ = Matrix(Nc_, 1.);
  const double bmin = param->getbmin();
  const double bmax = param->getbmax();
  const Matrix zero(Nc_, 0.);


/////test
//    std::vector<double> in1(Nc2m1_, 0);
//    std::vector<double> in2(Nc2m1_, 0);
//    for (int ii = 0; ii < 8; ii++) {
//        in1[ii] = 10*random->genrand64_real1();
//        in2[ii] = 10*random->genrand64_real1();
//    }
//    Matrix U1 = getUfromExponent(in1);
//    Matrix U2 = getUfromExponent(in2);
//    Matrix USol;
//    findUInForwardLightcone(U1, U2, USol);
//    for (int ii = 0; ii < Nc_; ii++) {
//        for (int jj = 0; jj < Nc_; jj++) {
//            cout << U1(ii, jj) << " "
//                 << U2(ii, jj) << " "
//                 << USol(ii, jj) << endl;
//        }
//    }
//    exit(0);
/////

  messager.info("Initializing fields ... ");
  param->setRnp(0.);

  double b = 0.;
  double xb = random->genrand64_real1();  // uniformly distributed random variable

  if (param->getUseNucleus() == 0) {
    // use b=0 fm for the constant g^2 mu case
    param->setSuccess(1);
    b = 0.;
    messager << "Setting b=0 for constant color charge density case.";
    messager.flush("info");
  } else {
    if (param->getLinearb() == 1) {
      // use a linear probability distribution for b if we are doing nuclei
      messager << "Sampling linearly distributed b between " << bmin << " and "
               << bmax << "fm. Found ";
      b = sqrt((bmax * bmax - bmin * bmin) * xb + bmin * bmin);
    } else {
      // use a uniform distribution instead
      messager << "Sampling uniformly distributed b between " << bmin << " and "
               << bmax << "fm. Found ";
      b = (bmax - bmin) * xb + bmin;
    }
  }

  param->setb(b);
  messager << "b=" << b << " fm.";
  messager.flush("info");

  // read Q_s^2 from file
  if (param->getUseNucleus() == 1) {
    readNuclearQs(param);
  }

  readInNucleusConfigs(static_cast<int>(glauber->nucleusA1()),
                       param->getlightNucleusOption(),
                       nucleonPosArrA_);
  readInNucleusConfigs(static_cast<int>(glauber->nucleusA2()),
                       param->getlightNucleusOption(),
                       nucleonPosArrB_);

  // sample nucleon positions
  nucleusA_.clear();
  nucleusB_.clear();

  // to read Wilson lines from file (e.g. after JIMWLK evolution for the
  // 3DGlasma)
  if (READFROMFILE > 0) {
    readV(lat, param, READFROMFILE);
    param->setSuccess(1);
  } else {
    // to generate your own Wilson lines
    if (param->getUseNucleus() == 1) {
      sampleTA(param, random, glauber); // populate the lists nucleusA_ and
                                        // nucleusB_ with position data of the
    }

    // set color charge densities
    setColorChargeDensity(lat, param, random, glauber);

    // for enforcing a specific Npart:
    if (param->getUseNucleus() == 1 && param->getUseFixedNpart() != 0) {
      if (param->getNpart() != param->getUseFixedNpart()) {
        while (param->getNpart() != param->getUseFixedNpart()) {
          cout << "resampling... desired Npart=" << param->getUseFixedNpart()
               << endl;
          nucleusA_.clear();
          nucleusB_.clear();

          xb = random
                   ->genrand64_real1(); // uniformly distributed random variable

          if (param->getLinearb() == 1) // use a linear probability distribution
                                        // for b if we are doing nuclei
          {
            cout << "Sampling linearly distributed b between " << bmin
                 << " and " << bmax << "fm." << endl;
            b = sqrt((bmax * bmax - bmin * bmin) * xb + bmin * bmin);
          } else // use a uniform distribution instead
          {
            cout << "Sampling uniformly distributed b between " << bmin
                 << " and " << bmax << "fm." << endl;
            b = (bmax - bmin) * xb + bmin;
          }

          param->setb(b);
          cout << "Using b=" << b << " fm" << endl;

          // populate the lists nucleusA_ and nucleusB_ with position data
          sampleTA(param, random, glauber);

          setColorChargeDensity(lat, param, random, glauber);
        }
      }
      cout << "Using fixed Npart=" << param->getNpart() << endl;
    }

    if (param->getSuccess() == 0) {
      cout << "No collision happened on rank " << param->getMPIRank()
           << ". Restarting with new random number..." << endl;
      return;
    }
    // sample color charges and find Wilson lines V_A and V_B
    setV(lat, param, random);
  }

  // output Wilson lines (used also for the proton plots)
  double L = param->getL();
  double a = L / N; // lattice spacing in fm
  double x,y;

  messager.info("Finding fields in forward lightcone...");

#pragma omp parallel
  {
    Matrix temp2(Nc_, 0.);
    Matrix Ux(Nc_, 0.);
    Matrix Uy(Nc_, 0.);
    Matrix Ux1(Nc_, 0.);
    Matrix Uy1(Nc_, 0.);
    Matrix Ux2(Nc_, 0.);
    Matrix Uy2(Nc_, 0.);
    Matrix UD(Nc_, 0.);
    Matrix UDx(Nc_, 0.);
    Matrix UDy(Nc_, 0.);
    Matrix UDx1(Nc_, 0.);
    Matrix UDy1(Nc_, 0.);

    Matrix Uplaq(Nc_, 0.);
    Matrix Uplaq1(Nc_, 0.);
    Matrix Uplaq2(Nc_, 0.);
    Matrix Uplaq3(Nc_, 0.);
    Matrix Uplaq4(Nc_, 0.);

    Matrix UD2(Nc_, 0.);
    Matrix UDx2(Nc_, 0.);
    Matrix UDy2(Nc_, 0.);

    Matrix Ax(Nc_, 0.);
    Matrix Ay(Nc_, 0.);
    Matrix Ax1(Nc_, 0.);
    Matrix Ay1(Nc_, 0.);
    Matrix AT(Nc_, 0.);

    Matrix Ax2(Nc_, 0.);
    Matrix Ay2(Nc_, 0.);
    Matrix AT2(Nc_, 0.);

    // field strength tensor
    Matrix Fxy(Nc_, 0.);
    Matrix Fyx(Nc_, 0.);

    Matrix AM(Nc_, 0.);
    Matrix AP(Nc_, 0.);

    Matrix AxUpY(Nc_, 0.);
    Matrix AyUpY(Nc_, 0.);

    Matrix Aeta2(Nc_, 0.);
    Matrix Ux1pUx2(Nc_, 0.);
    Matrix UDx1pUDx2(Nc_, 0.);
    Matrix Uy1pUy2(Nc_, 0.);
    Matrix UDy1pUDy2(Nc_, 0.);
    Matrix Ux1mUx2(Nc_, 0.);
    Matrix UDx1mUDx2(Nc_, 0.);
    Matrix Uy1mUy2(Nc_, 0.);
    Matrix UDy1mUDy2(Nc_, 0.);

    // compute Ux(3) Uy(3) after the collision

    //    ofstream fout("test1", ios::out);
 #pragma omp for
    for (int pos = 0; pos < N * N; pos++) // loops over all cells
    {
      auto checkU = lat->cells[pos]->getU().trace();
      if (checkU != checkU) {
        lat->cells[pos]->setU(one_);
      }

      checkU = lat->cells[pos]->getU2().trace();
      if (checkU != checkU) {
        lat->cells[pos]->setU2(one_);
      }

      ////check - remove later
      //      fout << lat->cells[pos]->getU() << endl;
      //fout << lat->cells[pos]->getU2() << endl;

    }

    //fout.close();

#pragma omp for
    for (int pos = 0; pos < N * N; pos++) // loops over all cells
    {
      UDx = lat->cells[lat->pospX[pos]]->getU();
      UDx.conjg();
      lat->cells[pos]->setUx1(lat->cells[pos]->getU() * UDx);

      UDy = lat->cells[lat->pospY[pos]]->getU();
      UDy.conjg();
      lat->cells[pos]->setUy1(lat->cells[pos]->getU() * UDy);

      UDx = lat->cells[lat->pospX[pos]]->getU2();
      UDx.conjg();
      lat->cells[pos]->setUx2(lat->cells[pos]->getU2() * UDx);

      UDy = lat->cells[lat->pospY[pos]]->getU2();
      UDy.conjg();
      lat->cells[pos]->setUy2(lat->cells[pos]->getU2() * UDy);
    }
    // -----------------------------------------------------------------
    // from Ux(1,2) and Uy(1,2) compute Ux(3) and Uy(3):

#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
      // loops over all cells
      UDx1 = lat->cells[pos]->getUx1();
      UDx2 = lat->cells[pos]->getUx2();
      bool status = findUInForwardLightcone(UDx1, UDx2, temp2);
      lat->cells[pos]->setUx(temp2);
      if (!status) {
        cout << "pos x = " << pos/512 << " y = " << pos%512 << endl;
      }

      UDy1 = lat->cells[pos]->getUy1();
      UDy2 = lat->cells[pos]->getUy2();
      status = findUInForwardLightcone(UDy1, UDy2, temp2);
      lat->cells[pos]->setUy(temp2);
      if (!status) {
        cout << "pos x = " << pos/512 << " y = " << pos%512 << endl;
      }
    }

// compute initial electric field
// with minus ax, ay
#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
      // x part in sum:
      Ux1mUx2 = lat->cells[pos]->getUx1() - lat->cells[pos]->getUx2();
      UDx1 = lat->cells[pos]->getUx1();
      UDx1.conjg();
      UDx2 = lat->cells[pos]->getUx2();
      UDx2.conjg();
      UDx1mUDx2 = UDx1 - UDx2;

      Ux = lat->cells[pos]->getUx();
      UDx = Ux;
      UDx.conjg();

      temp2 = Ux1mUx2 * UDx - Ux1mUx2 - Ux * UDx1mUDx2 + UDx1mUDx2;

      Ux1mUx2 = lat->cells[lat->posmX[pos]]->getUx1() -
                lat->cells[lat->posmX[pos]]->getUx2();
      UDx1 = lat->cells[lat->posmX[pos]]->getUx1();
      UDx1.conjg();
      UDx2 = lat->cells[lat->posmX[pos]]->getUx2();
      UDx2.conjg();
      UDx1mUDx2 = UDx1 - UDx2;

      Ux = lat->cells[lat->posmX[pos]]->getUx();
      UDx = Ux;
      UDx.conjg();

      temp2 = temp2 - UDx * Ux1mUx2 + Ux1mUx2 + UDx1mUDx2 * Ux - UDx1mUDx2;

      // y part in sum
      Uy1mUy2 = lat->cells[pos]->getUy1() - lat->cells[pos]->getUy2();
      UDy1 = lat->cells[pos]->getUy1();
      UDy1.conjg();
      UDy2 = lat->cells[pos]->getUy2();
      UDy2.conjg();
      UDy1mUDy2 = UDy1 - UDy2;

      Uy = lat->cells[pos]->getUy();
      UDy = Uy;
      UDy.conjg();

      // y part of the sum:
      temp2 = temp2 + Uy1mUy2 * UDy - Uy1mUy2 - Uy * UDy1mUDy2 + UDy1mUDy2;

      Uy1mUy2 = lat->cells[lat->posmY[pos]]->getUy1() -
                lat->cells[lat->posmY[pos]]->getUy2();
      UDy1 = lat->cells[lat->posmY[pos]]->getUy1();
      UDy1.conjg();
      UDy2 = lat->cells[lat->posmY[pos]]->getUy2();
      UDy2.conjg();
      UDy1mUDy2 = UDy1 - UDy2;

      Uy = lat->cells[lat->posmY[pos]]->getUy();
      UDy = Uy;
      UDy.conjg();

      temp2 = temp2 - UDy * Uy1mUy2 + Uy1mUy2 + UDy1mUDy2 * Uy - UDy1mUDy2;

      lat->cells[pos]->setE1((1. / 8.) * temp2);
    }

    // with plus ax, ay
#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
      // x part in sum:
      Ux1mUx2 = lat->cells[pos]->getUx1() - lat->cells[pos]->getUx2();
      UDx1 = lat->cells[pos]->getUx1();
      UDx1.conjg();
      UDx2 = lat->cells[pos]->getUx2();
      UDx2.conjg();
      UDx1mUDx2 = UDx1 - UDx2;

      Ux = lat->cells[pos]->getUx();
      UDx = Ux;
      UDx.conjg();

      temp2 = Ux1mUx2 * UDx - Ux1mUx2 - Ux * UDx1mUDx2 + UDx1mUDx2;

      Ux1mUx2 = lat->cells[lat->pospX[pos]]->getUx1() -
                lat->cells[lat->pospX[pos]]->getUx2();
      UDx1 = lat->cells[lat->pospX[pos]]->getUx1();
      UDx1.conjg();
      UDx2 = lat->cells[lat->pospX[pos]]->getUx2();
      UDx2.conjg();
      UDx1mUDx2 = UDx1 - UDx2;

      Ux = lat->cells[lat->pospX[pos]]->getUx();
      UDx = Ux;
      UDx.conjg();

      temp2 = temp2 - UDx * Ux1mUx2 + Ux1mUx2 + UDx1mUDx2 * Ux - UDx1mUDx2;

      // y part in sum
      Uy1mUy2 = lat->cells[pos]->getUy1() - lat->cells[pos]->getUy2();
      UDy1 = lat->cells[pos]->getUy1();
      UDy1.conjg();
      UDy2 = lat->cells[pos]->getUy2();
      UDy2.conjg();
      UDy1mUDy2 = UDy1 - UDy2;

      Uy = lat->cells[pos]->getUy();
      UDy = Uy;
      UDy.conjg();

      // y part of the sum:
      temp2 = temp2 + Uy1mUy2 * UDy - Uy1mUy2 - Uy * UDy1mUDy2 + UDy1mUDy2;

      Uy1mUy2 = lat->cells[lat->pospY[pos]]->getUy1() -
                lat->cells[lat->pospY[pos]]->getUy2();
      UDy1 = lat->cells[lat->pospY[pos]]->getUy1();
      UDy1.conjg();
      UDy2 = lat->cells[lat->pospY[pos]]->getUy2();
      UDy2.conjg();
      UDy1mUDy2 = UDy1 - UDy2;

      Uy = lat->cells[lat->pospY[pos]]->getUy();
      UDy = Uy;
      UDy.conjg();

      temp2 = temp2 - UDy * Uy1mUy2 + Uy1mUy2 + UDy1mUDy2 * Uy - UDy1mUDy2;

      lat->cells[pos]->setE2((1. / 8.) * temp2);
    }
// compute the plaquette
#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
      UDx = lat->cells[lat->pospY[pos]]->getUx();
      UDy = lat->cells[pos]->getUy();
      UDx.conjg();
      UDy.conjg();

      Uplaq = lat->cells[pos]->getUx() *
              (lat->cells[lat->pospX[pos]]->getUy() * (UDx * UDy));
      lat->cells[pos]->setUplaq(Uplaq);
    }

#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
      AM = (lat->cells[pos]->getE1()); //+lat->cells[pos]->getAetaP());
      AP = (lat->cells[pos]->getE2()); //+lat->cells[pos]->getAetaP());
      // this is pi in lattice units as needed for the evolution. (later, the
      // a^4 gives the right units for the energy density
      lat->cells[pos]->setpi(complex<double>(0., -2. / param->getg()) *
                             (AM)); // factor -2 because I have A^eta (note the
                                    // 1/8 before) but want \pi (E^z).
      // lat->cells[pos]->setpi(complex<double>(0.,-1./param->getg())*(AM+AP));
      // // factor -2 because I have A^eta (note the 1/8 before) but want \pi
      // (E^z).
    }

#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
      lat->cells[pos]->setE1(zero);
      lat->cells[pos]->setE2(zero);
      lat->cells[pos]->setphi(zero);
      lat->cells[pos]->setUx1(
          one_); // reset the Ux1 to be used for other purposes later
    }

  }

  // -----------------------------------------------------------------------------
  // finish
  // -----------------------------------------------------------------------------
}

void Init::multiplicity(Lattice *lat, Parameters *param) {
  int N = param->getSize();
  int pos;
  double epsilonSum = 0.;
  double L = param->getL();
  double a = L / N; // lattice spacing in fm

  for (int ix = 0; ix < N; ix++) {
    for (int iy = 0; iy < N; iy++) {
      pos = ix * N + iy;
      epsilonSum += a * a * lat->cells[pos]->getEpsilon() * hbarc;
    }
  }
  stringstream strtE_name;
  strtE_name << "totalEnergy" << param->getEventId() << ".dat";
  string tE_name;
  tE_name = strtE_name.str();

  ofstream fout(tE_name.c_str(), ios::out);
  fout << epsilonSum << endl;
  fout.close();
}

void Init::generate_nucleus_configuration(Random *random, int A, int Z,
                                          double a_WS, double R_WS,
                                          double beta2, double beta3,
                                          double beta4, double gamma,
                                          bool force_dmin_flag, double d_min,
                                          double dR_np, double da_np,
                                          std::vector<ReturnValue> &nucleus) {
  if (std::abs(beta2) < 1e-15 && std::abs(beta4) < 1e-15
          && std::abs(beta3) < 1e-15 && std::abs(gamma) < 1e-15) {
    generate_nucleus_configuration_with_woods_saxon(
                            random, A, Z, a_WS, R_WS, d_min, dR_np, da_np, nucleus);
  } else {
    if (force_dmin_flag) {
        generate_nucleus_configuration_with_deformed_woods_saxon_force_dmin(
                random, A, Z, a_WS, R_WS, beta2, beta3, beta4, gamma, d_min,
                dR_np, da_np, nucleus);
    } else {
      if (std::abs(gamma) > 1e-15) {
        generate_nucleus_configuration_with_deformed_woods_saxon2(
                random, A, Z, a_WS, R_WS, beta2, beta3, beta4, gamma,
                dR_np, da_np, nucleus);
      } else {
        generate_nucleus_configuration_with_deformed_woods_saxon(
                random, A, Z, a_WS, R_WS, beta2, beta3, beta4, d_min,
                dR_np, da_np, nucleus);
      }
    }
  }
}


void Init::generate_nucleus_configuration_with_woods_saxon(
    Random *random, int A, int Z, double a_WS, double R_WS, double d_min,
    double dR_np, double da_np, std::vector<ReturnValue> &nucleus) {
  std::vector<double> r_array(A, 0.);
  std::vector<int> idx_array(A, 0);
  for (int i = 0; i < Z; i++) {
    r_array[i] = sample_r_from_woods_saxon(random, a_WS, R_WS);
    idx_array[i] = i;
  }
  for (int i = Z; i < A; i++) {
    r_array[i] = sample_r_from_woods_saxon(random, a_WS + da_np, R_WS + dR_np);
    idx_array[i] = i;
  }
  std::stable_sort(idx_array.begin(), idx_array.end(),
            [&r_array](int i1, int i2) { return r_array[i1] < r_array[i2]; });
  std::stable_sort(r_array.begin(), r_array.end());

  std::vector<double> x_array(A, 0.), y_array(A, 0.), z_array(A, 0.);
  const double d_min_sq = d_min * d_min;
  for (int i = 0; i < A; i++) {
    double r_i = r_array[i];
    int reject_flag = 0;
    int iter = 0;
    double x_i, y_i, z_i;
    do {
      iter++;
      reject_flag = 0;
      double phi = 2. * M_PI * random->genrand64_real3();
      double theta = acos(1. - 2. * random->genrand64_real3());
      x_i = r_i * sin(theta) * cos(phi);
      y_i = r_i * sin(theta) * sin(phi);
      z_i = r_i * cos(theta);
      for (int j = i - 1; j >= 0; j--) {
        if ((r_i - r_array[j]) * (r_i - r_array[j]) > d_min_sq)
          break;
        double dsq = ((x_i - x_array[j]) * (x_i - x_array[j]) +
                      (y_i - y_array[j]) * (y_i - y_array[j]) +
                      (z_i - z_array[j]) * (z_i - z_array[j]));
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

  for (int i = 0; i < A; i++) {
    ReturnValue rv;
    rv.x = x_array[i];
    rv.y = y_array[i];
    rv.z = z_array[i];
    rv.phi = atan2(y_array[i], x_array[i]);
    rv.collided = 0;
    if (idx_array[i] < Z) {
      rv.proton = 1;
    } else {
      rv.proton = 1;
    }
    nucleus.push_back(rv);
  }
}


double Init::sample_r_from_woods_saxon(Random *random, double a_WS,
                                       double R_WS) const {
  double rmaxCut = R_WS + 10. * a_WS;
  double r = 0.;
  do {
    r = rmaxCut * pow(random->genrand64_real3(), 1.0 / 3.0);
  } while (random->genrand64_real3() > fermi_distribution(r, R_WS, a_WS));
  return (r);
}

double Init::fermi_distribution(double r, double R_WS, double a_WS) const {
  double f = 1. / (1. + exp((r - R_WS) / a_WS));
  return (f);
}

void Init::generate_nucleus_configuration_with_deformed_woods_saxon(
    Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
    double beta3, double beta4, double d_min, double dR_np, double da_np,
    std::vector<ReturnValue> &nucleus) {
  std::vector<double> r_array(A, 0.);
  std::vector<double> costheta_array(A, 0.);
  std::vector<int> idx_array(A, 0);
  for (int i = 0; i < Z; i++) {
    sample_r_and_costheta_from_deformed_woods_saxon(
        random, a_WS, R_WS, beta2, beta3, beta4, r_array[i], costheta_array[i]);
    idx_array[i] = i;
  }
  for (int i = Z; i < A; i++) {
    sample_r_and_costheta_from_deformed_woods_saxon(
        random, a_WS + da_np, R_WS + dR_np, beta2, beta3, beta4, r_array[i], costheta_array[i]);
    idx_array[i] = i;
  }
  std::stable_sort(idx_array.begin(), idx_array.end(),
            [&r_array](int i1, int i2) { return r_array[i1] < r_array[i2]; });
  std::stable_sort(r_array.begin(), r_array.end());

  std::vector<double> x_array(A, 0.), y_array(A, 0.), z_array(A, 0.);
  const double d_min_sq = d_min * d_min;
  for (int i = 0; i < A; i++) {
    const double r_i     = r_array[i];
    const double theta_i = acos(costheta_array[idx_array[i]]);
    int reject_flag = 0;
    int iter = 0;
    double x_i, y_i, z_i;
    do {
      iter++;
      reject_flag = 0;
      double phi = 2. * M_PI * random->genrand64_real3();
      x_i = r_i * sin(theta_i) * cos(phi);
      y_i = r_i * sin(theta_i) * sin(phi);
      z_i = r_i * cos(theta_i);
      for (int j = i - 1; j >= 0; j--) {
        if ((r_i - r_array[j])*(r_i - r_array[j]) > d_min_sq) break;
        double dsq = ((x_i - x_array[j]) * (x_i - x_array[j]) +
                      (y_i - y_array[j]) * (y_i - y_array[j]) +
                      (z_i - z_array[j]) * (z_i - z_array[j]));
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

  for (unsigned int i = 0; i < r_array.size(); i++) {
    ReturnValue rv;
    rv.x = x_array[i];
    rv.y = y_array[i];
    rv.z = z_array[i];
    rv.phi = atan2(y_array[i], x_array[i]);
    rv.collided = 0;
    if (idx_array[i] < Z) {
      rv.proton = 1;
    } else {
      rv.proton = 0;
    }
    nucleus.push_back(rv);
  }
}


void Init::generate_nucleus_configuration_with_deformed_woods_saxon_force_dmin(
    Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
    double beta3, double beta4, double gamma, double d_min,
    double dR_np, double da_np, std::vector<ReturnValue> &nucleus) {
  messager << "Sampling nucleon position forcing d_min = " << d_min
           << " fm ...";
  messager.flush("info");
  double rmaxCut = R_WS + dR_np + 10.*(a_WS + da_np);
  double r = 0.;
  double costheta = 0.;
  double phi = 0.;
  double R_WS_theta = 0.;
  const double d_min_sq = d_min * d_min;
  std::vector<double> x_array(A, 0.), y_array(A, 0.), z_array(A, 0.);
  for (int i = 0; i < A; i++) {
    double R_WS_i = R_WS;
    double a_WS_i = a_WS;
    if (i >= Z) {
        R_WS_i = R_WS + dR_np;
        a_WS_i = a_WS + da_np;
    }
    bool reSampleFlag = false;
    double x_i, y_i, z_i;
    do {
      // sample the position of the nucleon i
      do {
        r = rmaxCut*pow(random->genrand64_real3(), 1.0/3.0);
        costheta = 1.0 - 2.0 * random->genrand64_real3();
        phi  = 2.*M_PI*random->genrand64_real3();
        double y20 = spherical_harmonics(2, costheta);
        double y30 = spherical_harmonics(3, costheta);
        double y40 = spherical_harmonics(4, costheta);
        double y22 = spherical_harmonics_Y22(costheta, phi);
        R_WS_theta = R_WS_i*(1.0
                             + beta2*(cos(gamma)*y20 + sin(gamma)*y22)
                             + beta3*y30 + beta4*y40);
      } while (random->genrand64_real3()
               > fermi_distribution(r, R_WS_theta, a_WS_i));
      double sintheta = sqrt(1. - costheta*costheta);
      x_i = r*sintheta*cos(phi);
      y_i = r*sintheta*sin(phi);
      z_i = r*costheta;
      reSampleFlag = false;
      for (int j = i - 1; j >= 0; j--) {
        double r2 = (  (x_i - x_array[j])*(x_i - x_array[j])
                     + (y_i - y_array[j])*(y_i - y_array[j])
                     + (z_i - z_array[j])*(z_i - z_array[j]));
        if (r2 < d_min_sq) {
          reSampleFlag = true;
          break;
        }
      }
    } while (reSampleFlag);
    x_array[i] = x_i;
    y_array[i] = y_i;
    z_array[i] = z_i;
  }

  recenter_nucleus(x_array, y_array, z_array);

  for (int i = 0; i < A; i++) {
    ReturnValue rv;
    rv.x = x_array[i];
    rv.y = y_array[i];
    rv.z = z_array[i];
    rv.phi = atan2(y_array[i], x_array[i]);
    rv.collided = 0;
    if (i < Z) {
      rv.proton = 1;
    } else {
      rv.proton = 0;
    }
    nucleus.push_back(rv);
  }
}


void Init::generate_nucleus_configuration_with_deformed_woods_saxon2(
    Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
    double beta3, double beta4, double gamma, double dR_np, double da_np,
    std::vector<ReturnValue> &nucleus) {
  double rmaxCut = R_WS + dR_np + 10.*(a_WS + da_np);
  double r = 0.;
  double costheta = 0.;
  double phi = 0.;
  double R_WS_theta = 0.;
  std::vector<double> x_array(A, 0.), y_array(A, 0.), z_array(A, 0.);
  for (int i = 0; i < A; i++) {
    double R_WS_i = R_WS;
    double a_WS_i = a_WS;
    if (i >= Z) {
      // neutrons
      R_WS_i = R_WS + dR_np;
      a_WS_i = a_WS + da_np;
    }
    do {
        r = rmaxCut*pow(random->genrand64_real3(), 1.0/3.0);
        costheta = 1.0 - 2.0 * random->genrand64_real3();
        phi  = 2.*M_PI*random->genrand64_real3();
        double y20 = spherical_harmonics(2, costheta);
        double y30 = spherical_harmonics(3, costheta);
        double y40 = spherical_harmonics(4, costheta);
        double y22 = spherical_harmonics_Y22(costheta, phi);
        R_WS_theta = R_WS_i*(1.0
                             + beta2*(cos(gamma)*y20 + sin(gamma)*y22)
                             + beta3*y30 + beta4*y40);
    } while (random->genrand64_real3()
             > fermi_distribution(r, R_WS_theta, a_WS_i));
    double sintheta = sqrt(1. - costheta*costheta);
    x_array[i] = r*sintheta*cos(phi);
    y_array[i] = r*sintheta*sin(phi);
    z_array[i] = r*costheta;
  }

  recenter_nucleus(x_array, y_array, z_array);

  for (int i = 0; i < A; i++) {
    ReturnValue rv;
    rv.x = x_array[i];
    rv.y = y_array[i];
    rv.z = z_array[i];
    rv.phi = atan2(y_array[i], x_array[i]);
    rv.collided = 0;
    if (i < Z) {
      rv.proton = 1;
    } else {
      rv.proton = 0;
    }
    nucleus.push_back(rv);
  }
}


void Init::sample_r_and_costheta_from_deformed_woods_saxon(
    Random *random, double a_WS, double R_WS, double beta2, double beta3,
    double beta4, double &r, double &costheta) const {
  double rmaxCut = R_WS + 10. * a_WS;
  double R_WS_theta = R_WS;
  do {
    r = rmaxCut * pow(random->genrand64_real3(), 1.0 / 3.0);
    costheta = 1.0 - 2.0 * random->genrand64_real3();
    auto y20 = spherical_harmonics(2, costheta);
    auto y30 = spherical_harmonics(3, costheta);
    auto y40 = spherical_harmonics(4, costheta);
    R_WS_theta = R_WS * (1.0 + beta2 * y20 + beta3 * y30 + beta4 * y40);
  } while (random->genrand64_real3() > fermi_distribution(r, R_WS_theta, a_WS));
}

double Init::spherical_harmonics(int l, double ct) const {
  // Currently assuming m=0 and available for Y_{20} and Y_{40}
  // "ct" is cos(theta)
  double ylm = 0.0;
  if (l == 2) {
    ylm = 3.0 * ct * ct - 1.0;
    ylm *= 0.31539156525252005; // pow(5.0/16.0/M_PI,0.5);
  } else if (l == 3) {
    ylm  = 5.0*ct*ct*ct;
    ylm -= 3.0*ct;
    ylm *= 0.3731763325901154;  // pow(7.0/16.0/M_PI,0.5);
  } else if (l == 4) {
    ylm = 35.0 * ct * ct * ct * ct;
    ylm -= 30.0 * ct * ct;
    ylm += 3.0;
    ylm *= 0.10578554691520431; // 3.0/16.0/pow(M_PI,0.5);
  }
  return (ylm);
}

double Init::spherical_harmonics_Y22(double ct, double phi) const {
    // Y2,2
    double ylm = 0.0;
    ylm  = 1.0 - ct*ct;
    ylm *= cos(2.*phi);
    ylm *= 0.5462742152960397;  // sqrt(2*15)/4/pow(2*M_PI,0.5);
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


void Init::recenter_nucleus(std::vector<ReturnValue> &nucleus) {
    // compute the center of mass position and shift it to (0, 0, 0)
    double meanx = 0., meany = 0., meanz = 0.;
    for (auto &n_i: nucleus) {
        meanx += n_i.x;
        meany += n_i.y;
        meanz += n_i.z;
    }

    meanx /= static_cast<double>(nucleus.size());
    meany /= static_cast<double>(nucleus.size());
    meanz /= static_cast<double>(nucleus.size());

    for (auto &n_i: nucleus) {
        n_i.x -= meanx;
        n_i.y -= meany;
        n_i.z -= meanz;
    }
}


void Init::assignProtons(std::vector<ReturnValue> &nucleus, const int Z) {
    // randomly assign Z nucleons to be protons inside the nucleus
    std::random_shuffle(nucleus.begin(), nucleus.end());
    for (unsigned int i = 0; i < nucleus.size(); i++) {
        if (static_cast<int>(i) < std::abs(Z)) {
            nucleus.at(i).proton = 1;
        } else {
            nucleus.at(i).proton = 0;
        }
    }
}


void Init::rotate_nucleus(Random* random, std::vector<ReturnValue> &nucleus) {
  double phi_global = 2. * M_PI * random->genrand64_real3();
  double theta_global = acos(1. - 2. * random->genrand64_real3());
  auto cth = cos(theta_global);
  auto sth = sin(theta_global);
  auto cphi = cos(phi_global);
  auto sphi = sin(phi_global);
  for (auto &n_i: nucleus) {
    auto x_new = cth * cphi * n_i.x - sphi * n_i.y + sth * cphi * n_i.z;
    auto y_new = cth * sphi * n_i.x + cphi * n_i.y + sth * sphi * n_i.z;
    auto z_new = -sth * n_i.x + 0. * n_i.y + cth * n_i.z;
    n_i.x = x_new;
    n_i.y = y_new;
    n_i.z = z_new;
  }
}


void Init::rotate_nucleus_3D(Random* random,
                             std::vector<ReturnValue> &nucleus) {
  // rotate the nucleus with the full three solid angles
  // required for tri-axial deformed nuclei
  // https://en.wikipedia.org/wiki/Euler_angles
  double alpha = 2*M_PI*random->genrand64_real3();
  double beta = acos(1. - 2. * random->genrand64_real3());
  double gamma = 2*M_PI*random->genrand64_real3();
  auto c1 = cos(alpha);
  auto s1 = sin(alpha);
  auto c2 = cos(beta);
  auto s2 = sin(beta);
  auto c3 = cos(gamma);
  auto s3 = sin(gamma);
  for (auto &n_i: nucleus) {
    auto x_new = c2*n_i.x - c3*s2*n_i.y + s2*s3*n_i.z;
    auto y_new = c1*s2*n_i.x + (c1*c2*c3 - s1*s3)*n_i.y + (-c3*s1 - c1*c2*s3)*n_i.z;
    auto z_new = s1*s2*n_i.x + (c1*s3 + c2*c3*s1)*n_i.y + (c1*c3 - c2*s1*s3)*n_i.z;
    n_i.x = x_new;
    n_i.y = y_new;
    n_i.z = z_new;
  }
}


double Init::sampleLogNormalDistribution(Random *random,
                                         const double mean,
                                         const double variance) {
    const double meansq = mean*mean;
    const double mu = log(meansq/sqrt(variance + meansq));
    const double sigma = sqrt(log(variance/meansq + 1.));
    double sampleX = exp(mu + sigma*random->Gauss());
    return(sampleX);
}


void Init::sampleQsNormalization(Random *random,
                                 Parameters *param,
                                 const int Nq,
                                 vector<double> &gauss_array) {
    const double QsSmearWidth = param->getSmearingWidth();
    gauss_array.resize(Nq, 1.);    // default norm = 1
    if (param->getSmearQs() == 1) {
        // introduce a log-normal distribution for Qs normalization
        // dividing by exp(0.5 sigma^2) to ensure the mean is 1
        // the varirance in this case is exp(sigma) - 1 for the log-normal
        // distribution
        for (int iq = 0; iq < Nq; iq++) {
            gauss_array[iq] = (exp(random->Gauss(0, QsSmearWidth))/
                               exp(QsSmearWidth*QsSmearWidth/2.));
        }
    }
}


int Init::sampleNumberOfPartons(Random *random, Parameters *param) {
    double NqBase = param->getNqBase();
    int NqBaseInt = static_cast<int>(NqBase);
    double ran = random->genrand64_real2();
    int Nq = NqBaseInt;
    if (ran < NqBase - NqBaseInt) {
        Nq += 1;
    }
    Nq += random->Poisson(param->getNqFluc());
    return(std::max(1, Nq));
}


bool Init::findUInForwardLightcone(Matrix &U1, Matrix &U2, Matrix &Usol) {
    const int maxIterations = 2000;
    const int maxRetrys = 100;

    Matrix U1pU2 = U1 + U2;
    Matrix U1pU2dagger = U1pU2;
    U1pU2dagger.conjg();

    std::vector<double> alpha(Nc2m1_, 0.);          // solution
    std::vector<double> Dalpha(Nc2m1_, 0.);

    Matrix temp(Nc_, 0.);
    Matrix Mtemp(Nc_, 0.);
    std::vector<Matrix> MtempArr;
    MtempArr.resize(Nc2m1_);
    std::vector<complex<double>> traceCache(Nc2m1_, 0.);
    for (int ai = 0; ai < Nc2m1_; ai++) {
        temp = group_ptr_->getT(ai) * (U1pU2 - U1pU2dagger);
        traceCache[ai] = temp.trace();
        MtempArr[ai] = group_ptr_->getT(ai) * U1pU2;
    }

    // use raw pointers to interface with gsl
    double *Jab = new double [2 * Nc2m1_ * Nc2m1_];
    double *Fa = new double [2 * Nc2m1_];

    double Fzero = 10.;
    double FzeroMin = 1e6;
    Matrix UsolBestEst(Nc_, 1.);

    // set up initial guess
    Usol = U1*U2;
    Matrix Usoldagger = Usol;
    Usoldagger.conjg();

    int iter = 0;
    int nRestart = 0;
    while (Fzero > 1e-8 && iter < maxIterations && nRestart < maxRetrys) {
        iter++;

        Fzero = 0.;
        // compute function F that needs to be zero
        Mtemp = U1pU2 * Usoldagger - Usol * U1pU2dagger;
        for (int ai = 0; ai < Nc2m1_; ai++) {
            //temp = group_ptr_->getT(ai) * (U1pU2 - U1pU2dagger) +
            //       group_ptr_->getT(ai) * U1pU2 * Usoldagger -
            //       group_ptr_->getT(ai) * Usol * U1pU2dagger;
            //temp = group_ptr_->getT(ai) * Mtemp;
            complex<double> traceLoc = Mtemp.traceOfProdcutOfMatrix(
                                                group_ptr_->getT(ai), Mtemp);

            // minus trace if temp gives -F_ai
            auto traceRes = (-1.) * (traceCache[ai] + traceLoc);
            Fa[2 * ai] = real(traceRes);
            Fa[2 * ai + 1] = imag(traceRes);
            Fzero += std::abs(Fa[2 * ai]) + std::abs(Fa[2 * ai + 1]);
        }

        // compute Jacobian
        for (int bi = 0; bi < Nc2m1_; bi++) {
            Mtemp = group_ptr_->getT(bi) * Usoldagger;
            for (int ai = 0; ai < Nc2m1_; ai++) {
                int countMe = ai * Nc2m1_ + bi;
                //temp = group_ptr_->getT(ai) * U1pU2 * group_ptr_->getT(bi) * Usoldagger +
                //       group_ptr_->getT(ai) * Usol * group_ptr_->getT(bi) * U1pU2dagger;
                // -i times trace of temp gives my Jacobian matrix elements:
                //auto traceRes = complex<double>(0., -1.) * temp.trace();

                //temp = MtempArr[ai] * Mtemp;
                complex<double> traceLoc = Mtemp.traceOfProdcutOfMatrix(
                                                    MtempArr[ai], Mtemp);
                auto traceRes = complex<double>(0., -1.) * (2.*real(traceLoc));

                Jab[2 * countMe] = real(traceRes);
                Jab[2 * countMe + 1] = imag(traceRes);
            }
        }

        Dalpha = solveAxb(Jab, Fa);
        for (int ai = 0; ai < Nc2m1_; ai++) {
            alpha[ai] = alpha[ai] + Dalpha[ai];
        }

        Usol = getUfromExponent(alpha);
        Usoldagger = Usol;
        Usoldagger.conjg();

        if (iter == maxIterations) {
            if (Fzero < FzeroMin) {
                FzeroMin = Fzero;
                UsolBestEst = Usol;
            }

            for (int ai = 0; ai < Nc2m1_; ai++) {
                alpha[ai] = nRestart * random_ptr_->genrand64_real1();
            }
            Usol = getUfromExponent(alpha);
            Usoldagger = Usol;
            Usoldagger.conjg();

            nRestart++;
            iter = 0;
        }
    }
    bool success = true;
    if (nRestart == maxRetrys) {
        std::cout << "Did not converge in findUInForwardLightcone, "
                  << "Fzero: " << FzeroMin << std::endl;
        Usol = UsolBestEst;     // return the best estimate
        success = false;
    }
    delete[] Fa;
    delete[] Jab;
    return(success);
}


Matrix Init::getUfromExponent(std::vector<double> &in) {
    Matrix tempM(Nc_, 0.);

    // expmCoeff wil calculate exp(i in[a]t[a])
    auto U = tempM.expmCoeff(in, Nc_);
    if (std::abs(U[0].real()) < 1e-15) {
        tempM = one_;
    } else {
        tempM = (  U[0] * one_
                 + U[1] * group_ptr_->getT(0)
                 + U[2] * group_ptr_->getT(1)
                 + U[3] * group_ptr_->getT(2)
                 + U[4] * group_ptr_->getT(3)
                 + U[5] * group_ptr_->getT(4)
                 + U[6] * group_ptr_->getT(5)
                 + U[7] * group_ptr_->getT(6)
                 + U[8] * group_ptr_->getT(7)
        );
    }
    return(tempM);
}

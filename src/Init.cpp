// Init.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#include "Init.h"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Instrumentation.h"
#include "Phys_consts.h"
#include "gsl/gsl_linalg.h"

using PhysConst::hbarc;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;

namespace {

// SplitMix64 finalizer/counter step.  This is used only to construct a
// deterministic, stateless retry stream for the forward-light-cone solver.
inline std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

inline std::uint64_t forwardLightconeRetrySeed(
    std::uint64_t runSeed, int eventId, int pos, int direction) {
    std::uint64_t key = splitmix64(runSeed);
    key = splitmix64(
        key
        ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(eventId))
           + 0xD1B54A32D192ED03ULL));
    key = splitmix64(
        key
        ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(pos))
           + 0x94D049BB133111EBULL));
    key = splitmix64(
        key
        ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(direction))
           + 0xBF58476D1CE4E5B9ULL));
    return key;
}

inline double deterministicRetryGaussian(
    std::uint64_t seed, std::uint64_t drawIndex) {
    // Build two open-interval uniform doubles from independent SplitMix64
    // counters, then use Box-Muller.  No shared RNG state is touched.
    const std::uint64_t counter = 2ULL * drawIndex;
    const std::uint64_t r1 =
        splitmix64(seed + 0x9E3779B97F4A7C15ULL * (counter + 1ULL));
    const std::uint64_t r2 =
        splitmix64(seed + 0x9E3779B97F4A7C15ULL * (counter + 2ULL));

    constexpr double invTwo53 = 1.0 / 9007199254740992.0;
    const double u1 = (static_cast<double>(r1 >> 11) + 0.5) * invTwo53;
    const double u2 = (static_cast<double>(r2 >> 11) + 0.5) * invTwo53;
    constexpr double twoPi = 6.283185307179586476925286766559;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(twoPi * u2);
}

}  // namespace

//**************************************************************************
// Init class.

void Init::solveAxbComplex(double *Jab, double *Fa, std::vector<double> &xvec) {
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
}

void Init::solveAxb(double *Jab, double *Fa, std::vector<double> &xvec) {
    gsl_matrix_view m = gsl_matrix_view_array(Jab, Nc2m1_, Nc2m1_);
    gsl_vector_view c = gsl_vector_view_array(Fa, Nc2m1_);
    gsl_vector *x = gsl_vector_alloc(Nc2m1_);

    int s;
    gsl_permutation *p = gsl_permutation_alloc(Nc2m1_);
    gsl_linalg_LU_decomp(&m.matrix, p, &s);
    gsl_linalg_LU_solve(&m.matrix, p, &c.vector, x);
    gsl_permutation_free(p);

    for (int i = 0; i < Nc2m1_; i++) {
        xvec[i] = gsl_vector_get(x, i);
    }

    gsl_vector_free(x);
}

// This function samples the nucleon positions inside the projectile and
// target nuclei. Both nuclei are centered at the origin.
void Init::sampleTA(Parameters *param, Random *random, Glauber *glauber) {
    IPG_PROFILE_SCOPE("initialization.sample_nuclei");
    ReturnValue rv, rv2;
    messager.info("Sampling nucleon positions ... ");

    if (param->getNucleonPositionsFromFile() == 0) {
        if (param->getAverageOverNuclei() > 1) {
            if ((glauber->nucleusA1() == 1 || glauber->nucleusA2() == 1)) {
                cerr << "Averaging not supported for collisions involving "
                        "protons "
                        "... Exiting."
                     << std::endl;
                exit(1);
            }
        }

        int A1 = glauber->nucleusA1();  // projectile
        int A2 = glauber->nucleusA2();  // target
        int Z1 = glauber->nucleusZ1();  // projectile
        int Z2 = glauber->nucleusZ2();  // target

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
            // we sample the neutron proton distance, so distance to the center
            // needs to be divided by 2
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
            rv.proton = 0;
            rv.collided = 0;
            nucleusA_.push_back(rv);
        } else {
            generate_nucleus_configuration(
                random, A1, Z1, glauber->GlauberData.Projectile.a_WS,
                glauber->GlauberData.Projectile.R_WS,
                glauber->GlauberData.Projectile.beta2,
                glauber->GlauberData.Projectile.beta3,
                glauber->GlauberData.Projectile.beta4,
                glauber->GlauberData.Projectile.gamma,
                glauber->GlauberData.Projectile.forceDminFlag,
                glauber->GlauberData.Projectile.d_min,
                glauber->GlauberData.Projectile.dR_np,
                glauber->GlauberData.Projectile.da_np, nucleusA_);
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
            // we sample the neutron proton distance, so distance to the center
            // needs to be divided by 2
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
                random, A2, Z2, glauber->GlauberData.Target.a_WS,
                glauber->GlauberData.Target.R_WS,
                glauber->GlauberData.Target.beta2,
                glauber->GlauberData.Target.beta3,
                glauber->GlauberData.Target.beta4,
                glauber->GlauberData.Target.gamma,
                glauber->GlauberData.Target.forceDminFlag,
                glauber->GlauberData.Target.d_min,
                glauber->GlauberData.Target.dR_np,
                glauber->GlauberData.Target.da_np, nucleusB_);
        }
    } else if (param->getNucleonPositionsFromFile() == 1) {
        if (nucleonPosArrA_.size() > 0) {
            double ran2 = random->genrand64_real3();
            int nucleusNumber = static_cast<int>(ran2 * nucleonPosArrA_.size());
            std::cout << "using nucleus Number = " << nucleusNumber
                      << std::endl;
            for (int iA = 0; iA < glauber->nucleusA1(); iA++) {
                rv.x = nucleonPosArrA_[nucleusNumber][3 * iA];
                rv.y = nucleonPosArrA_[nucleusNumber][3 * iA + 1];
                rv.z = nucleonPosArrA_[nucleusNumber][3 * iA + 2];
                rv.collided = 0;
                nucleusA_.push_back(rv);
            }
            assignProtons(random, nucleusA_, glauber->nucleusZ1());
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
                glauber->GlauberData.Projectile.da_np, nucleusA_);
        }

        if (nucleonPosArrB_.size() > 0) {
            double ran2 = random->genrand64_real3();
            int nucleusNumber = static_cast<int>(ran2 * nucleonPosArrB_.size());
            std::cout << "using nucleus Number = " << nucleusNumber
                      << std::endl;
            for (int iA = 0; iA < glauber->nucleusA2(); iA++) {
                rv.x = nucleonPosArrB_[nucleusNumber][3 * iA];
                rv.y = nucleonPosArrB_[nucleusNumber][3 * iA + 1];
                rv.z = nucleonPosArrB_[nucleusNumber][3 * iA + 2];
                rv.collided = 0;
                nucleusB_.push_back(rv);
            }
            assignProtons(random, nucleusB_, glauber->nucleusZ2());
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
                glauber->GlauberData.Target.da_np, nucleusB_);
        }
    } else if (param->getNucleonPositionsFromFile() == 2) {
        // Read in Alvioli's nucleon positions including correlations
        if (glauber->nucleusA1() != 208 && glauber->nucleusA2() != 208) {
            cerr << "[Init.cpp]: The option 'getNucleonPositionsFromFile == 2' "
                    "only "
                    "works for either both nuclei Pb-208 or Projectile p and "
                    "Target "
                    "Pb-208. Exiting."
                 << std::endl;
            exit(1);
        }

        std::cout << "Retrieving nuclei from " << std::endl;

        // generate the file name
        double ran =
            random->genrand64_real3();  // sample the file name uniformly
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

        cout << "Reading nucleon positions for nucleus A from file " << fileName
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
                    fin >> dummy;  // don't care about isospin
                    rv.collided = 0;
                    nucleusA_.push_back(rv);
                    A++;
                    //        cout << "A=" << A << "/" <<
                    // glauber->nucleusA1()<<endl; cout << rv.x << " " << rv.y
                    // << endl;
                }
            }
        }
        fin.close();
        // do the second nucleus (Target)

        // generate the file name
        ran = random->genrand64_real3();  // sample the file name uniformly
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
        ran2 = random->genrand64_real3();  // sample the position in the file
                                           // uniformly (10,000 events per file)
        nucleusNumber = static_cast<int>(ran2 * 10000);
        cout << "Nucleus Number = " << nucleusNumber << endl;

        // go to the correct line in the file
        fin.seekg(std::ios::beg);
        for (int i = 0; i < (nucleusNumber)*glauber->nucleusA2(); ++i) {
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
                    fin >> rv.z;   // don't care about z direction
                    fin >> dummy;  // don't care about isospin
                    rv.collided = 0;
                    nucleusB_.push_back(rv);
                    A2++;
                    //    cout << "A2=" << A2 << "/" <<
                    // glauber->nucleusA2()<<endl; cout << rv.x << " " << rv.y
                    // << endl;
                }
            }
        }

        fin.close();

        // The files provide only coordinates.
        // Assign proton/neutron labels
        assignProtons(random, nucleusA_, glauber->nucleusZ1());
        assignProtons(random, nucleusB_, glauber->nucleusZ2());

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
    if (param->getPolarizationProjectile() == 0) {
        rotate_nucleus_3D(random, nucleusA_);
    } else if (param->getPolarizationProjectile() == 1) {
        // longitudinal polarization only rotates phi randomly
        double phi = 2. * M_PI * random->genrand64_real3();
        double theta = 0;
        rotate_nucleus(phi, theta, nucleusA_);
    } else if (param->getPolarizationProjectile() == 2) {
        // transverse polarization rotates J to +y axis
        double phi = M_PI / 2;
        double theta = M_PI / 2;
        rotate_nucleus(phi, theta, nucleusA_);
    }

    if (param->getPolarizationTarget() == 0) {
        rotate_nucleus_3D(random, nucleusB_);
    } else if (param->getPolarizationTarget() == 1) {
        // longitudinal polarization only rotates phi randomly
        double phi = 2. * M_PI * random->genrand64_real3();
        double theta = 0;
        rotate_nucleus(phi, theta, nucleusB_);
    } else if (param->getPolarizationTarget() == 2) {
        // transverse polarization rotates J to +y axis
        double phi = M_PI / 2;
        double theta = M_PI / 2;
        rotate_nucleus(phi, theta, nucleusB_);
    }
}

void Init::readNuclearQs(Parameters *param) {
    IPG_PROFILE_SCOPE("initialization.read_qs_table");
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
                    cerr << " End of file reached prematurely. Did the file "
                            "change? "
                            "Exiting."
                         << std::endl;
                    exit(1);
                }
            }
        }
        fin.close();
    } else {
        std::cout << "[Init.cpp:readNuclearQs]: File "
                  << param->getNucleusQsTableFileName()
                  << " does not exist. Exiting." << std::endl;
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

void Init::readInNucleusConfigs(
    const int nucleusA, const int lightNucleusOption,
    const int polarizationFlag, const double polJz,
    vector<vector<float>> &nucleonPosArr) {
    if (nucleonPosArr.size() > 0) return;
    std::string path = "nucleusConfigurations/";
    std::string fileName;
    bool readFlag = true;
    if (nucleusA == 2) {
        fileName = "DeuteronPol0Configs.bin.in";
        if ((std::abs(polJz) - 1) < 1e-8)
            fileName = "DeuteronPolpm1Configs.bin.in";
        if (polarizationFlag == 0) {
            auto ran = random_ptr_->genrand64_real1();
            if (ran < 0.3333) {
                fileName = "DeuteronPol0Configs.bin.in";
            } else {
                fileName = "DeuteronPolpm1Configs.bin.in";
            }
        }
    } else if (nucleusA == 3) {
        fileName = "He3.bin.in";
        if (lightNucleusOption == 1) fileName = "triton.bin.in";
    } else if (nucleusA == 4) {
        fileName = "He4.bin.in";
    } else if (nucleusA == 12) {
        fileName = "C12_VMC.bin.in";
        if (lightNucleusOption == 1) fileName = "C12_alphaCluster.bin.in";
    } else if (nucleusA == 16) {
        fileName = "O16_VMC.bin.in";
        if (lightNucleusOption == 1) {
            fileName = "O16_alphaCluster.bin.in";
        } else if (lightNucleusOption == 2) {
            fileName = "O16_PGCM_clustered_dmin0.bin.in";
        } else if (lightNucleusOption == 3) {
            fileName = "O16_PGCM_uniform_dmin0.bin.in";
        } else if (lightNucleusOption == 4) {
            fileName = "O16_NLEFT_dmin0.5fm_positiveweights.bin.in";
        } else if (lightNucleusOption == 5) {
            fileName = "O16_NLEFT_dmin0.5fm_negativeweights.bin.in";
        }
    } else if (nucleusA == 20) {
        fileName = "Ne20_PGCM_clustered_dmin0.bin.in";
        if (lightNucleusOption == 3) {
            fileName = "Ne20_PGCM_uniform_dmin0.bin.in";
        } else if (lightNucleusOption == 4) {
            fileName = "Ne20_NLEFT_dmin0.5fm_positiveweights.bin.in";
        } else if (lightNucleusOption == 5) {
            fileName = "Ne20_NLEFT_dmin0.5fm_negativeweights.bin.in";
        }
    } else if (nucleusA == 22) {
        fileName = "Ne22_NLEFT.bin.in";
    } else if (nucleusA == 40) {
        fileName = "Ar40_VMC.bin.in";
        if (lightNucleusOption == 4) fileName = "Ar40_NLEFT.bin.in";
    } else if (nucleusA == 197) {
        fileName = "Au197.bin.in";
    } else if (nucleusA == 208) {
        fileName = "Pb208.bin.in";
    } else {
        readFlag = false;
    }

    if (!readFlag) return;

    int Nentry = 3;
    if (nucleusA == 197 || nucleusA == 208) {
        Nentry = 4;
    }

    fileName = path + fileName;
    messager << "read in nucleus configurations from " << fileName;
    messager.flush("info");
    std::ifstream inFile(fileName, std::ios::binary);
    while (true) {
        vector<float> tempPos;
        for (int i = 0; i < nucleusA; i++) {
            for (int j = 0; j < Nentry; j++) {
                float temp;
                inFile.read(reinterpret_cast<char *>(&temp), sizeof(float));
                tempPos.push_back(temp);
            }
        }
        if (inFile.eof()) break;
        nucleonPosArr.push_back(tempPos);
    }
    inFile.close();
    messager << "read in " << nucleonPosArr.size() << " configurations.";
    messager.flush("info");
}

void Init::samplePartonPositions(
    Parameters *param, Random *random, vector<double> &x_array,
    vector<double> &y_array, vector<double> &z_array,
    vector<double> &BGq_array) {
    const double sqrtBG = sqrt(param->getBG()) * hbarc;  // fm
    const double BGqMean = param->getBGq();
    const double BGqVar = param->getBGqVar();
    const double BGq =
        (0.09 + sampleLogNormalDistribution(random, BGqMean - 0.09, BGqVar));
    const double QsSmearWidth = param->getSmearingWidth();
    const int Nq = sampleNumberOfPartons(random, param);
    const double dq_min = param->getDqmin();  // fm
    const double dq_min_sq = dq_min * dq_min;
    const double omega = param->getOmega();

    vector<double> r_array(Nq, 0.);
    BGq_array.assign(Nq, BGq);
    for (int iq = 0; iq < Nq; iq++) {
        if (std::abs(omega - 1) < 1e-8) {
            double xq = sqrtBG * random->Gauss();
            double yq = sqrtBG * random->Gauss();
            double zq = sqrtBG * random->Gauss();
            r_array[iq] = sqrt(xq * xq + yq * yq + zq * zq);
        } else {
            double bperp = sqrtBG * sqrt(omega * random->sampleGammaInc());
            r_array[iq] = bperp;  // bperp in 2D (asuume z = 0)
        }
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
            reject_flag = 0;
            double phi = 2. * M_PI * random->genrand64_real2();
            double theta = acos(1. - 2. * random->genrand64_real2());
            if (std::abs(omega - 1) < 1e-8) {
                x_i = r_i * sin(theta) * cos(phi);
                y_i = r_i * sin(theta) * sin(phi);
                z_i = r_i * cos(theta);
            } else {
                x_i = r_i * cos(phi);
                y_i = r_i * sin(phi);
                z_i = 0.;  // assume z=0
            }
            for (int j = i - 1; j >= 0; j--) {
                if ((r_i - r_array[j]) * (r_i - r_array[j]) > dq_min_sq) break;
                double dsq =
                    ((x_i - x_array[j]) * (x_i - x_array[j])
                     + (y_i - y_array[j]) * (y_i - y_array[j])
                     + (z_i - z_array[j]) * (z_i - z_array[j]));
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

// Q_s as a function of \sum T_p and y (new in this version of the code -
// v1.2 and up)
double Init::getNuclearQs2(double T, double y) {
    double value, fracy, fracT, QsYdown, QsYup;
    int posy, check = 0;
    fracy = 0.;
    posy = static_cast<int>(floor(y / deltaYNuc + 0.0000001));

    if (y > iymaxNuc * deltaYNuc) {
        cout << " [Init:getNuclearQs2]:ERROR: y out of range. Maximum y "
                "value "
                "is "
             << iymaxNuc * deltaYNuc << ", you used " << y << ". Exiting."
             << endl;
        exit(1);
    }

    //  if ( T > Qs2Nuclear[iTpmax-1][iymaxNuc-1] )
    if (T > Tlist[iTpmax - 1]) {
        cerr << "T=" << T << ", maximal T in table=" << Tlist[iTpmax - 1]
             << endl;
        cerr << " [Init:getNuclearQs2]:WARNING: out of range. Using "
                "maximal T "
                "in "
                "table."
             << endl;
        check = 1;
        fracy = (y - static_cast<double>(posy) * deltaYNuc) / deltaYNuc;
        QsYdown = (Qs2Nuclear[iTpmax - 1][posy]);
        QsYup = (Qs2Nuclear[iTpmax - 1][posy + 1]);
        value = (fracy * QsYup + (1. - fracy) * QsYdown);  //*hbarc*hbarc;
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

            QsYdown = (fracT) * (Qs2Nuclear[iT + 1][posy])
                      + (1. - fracT) * (Qs2Nuclear[iT][posy]);
            QsYup = (fracT) * (Qs2Nuclear[iT + 1][posy + 1])
                    + (1. - fracT) * (Qs2Nuclear[iT][posy + 1]);
            value = (fracy * QsYup + (1. - fracy) * QsYdown);  //*hbarc*hbarc;

            check++;
            continue;
        }
    }

    if (check != 1) {
        cout << check << ": T=" << T << endl;
        cerr << " [Init:getNuclearQs2]:ERROR: something went wrong in "
                "determining "
                "the value of Qs^2. Using maximal T_p"
             << endl;
        value =
            (fracy * Qs2Nuclear[iTpmax - 1][posy + 1]
             + (1. - fracy) * Qs2Nuclear[iTpmax - 1][posy]);
    }

    return value;
}

// set g^2\mu^2 as the sum of the individual nucleons' g^2\mu^2, using
// Q_s(b,y) prop to g^mu(b,y) Also compute N_part using Glauber If
// param->getwhich_stage() == 2, then here we shift nuclei back to b=0 for
// JIMLWK evolution (to be shifted back to b after JIMLWK)
void Init::setColorChargeDensity(
    Lattice *lat, Parameters *param, Random *random, Glauber *glauber) {
    IPG_PROFILE_SCOPE("initialization.color_charge_density");
    messager.info("set color charge density ...");

    const int N = param->getSize();
    const double L = param->getL();
    const double a = L / N;  // lattice spacing in fm

    int pos, posA, posB;
    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();

    double rapidityA = 0.;
    double rapidityB = 0.;
    if (param->getUsePseudoRapidity() == 0) {
        rapidityA = param->getRapidityA();
        rapidityB = param->getRapidityB();
    } else {
        // when using pseudorapidity as input convert to rapidity here.
        // later include Jacobian in multiplicity and energy
        cout << "Using pseudorapidity " << param->getRapidityA() << ", "
             << param->getRapidityB() << endl;
        double m = param->getJacobianm();  // in GeV
        double P =
            0.13 + 0.32 * pow(param->getRoots() / 1000., 0.115);  // in GeV
        rapidityA =
            0.5
            * log(
                sqrt(pow(cosh(param->getRapidityA()), 2.) + m * m / (P * P))
                + sinh(param->getRapidityA())
                      / (sqrt(
                             pow(cosh(param->getRapidityA()), 2.)
                             + m * m / (P * P))
                         - sinh(param->getRapidityA())));
        rapidityB =
            0.5
            * log(
                sqrt(pow(cosh(param->getRapidityB()), 2.) + m * m / (P * P))
                + sinh(param->getRapidityB())
                      / (sqrt(
                             pow(cosh(param->getRapidityB()), 2.)
                             + m * m / (P * P))
                         - sinh(param->getRapidityB())));
        cout << "Corresponds to rapidity " << rapidityA << ", " << rapidityB
             << endl;
    }

    double nucleiInAverage = static_cast<double>(param->getAverageOverNuclei());

    param->setQsmuRatioB(param->getQsmuRatio());

    if (param->getUseNucleus() == 0) {
        if (param->getUseGaussian() == 1) {
            double sigmax = 0.35;
            double sigmay = 0.5;
            for (int ix = 0; ix < N; ix++)  // loop over all positions
            {
                double x = ix * L / double(N) - L / 2.;
                for (int iy = 0; iy < N; iy++) {
                    double y = iy * L / double(N) - L / 2.;
                    int localpos = ix * N + iy;
                    double envelope = exp(
                                          -(x * x / (2. * sigmax * sigmax)
                                            + y * y / (2. * sigmay * sigmay)))
                                      / (2. * M_PI * sigmax * sigmay);
                    lat->cells[localpos]->setg2mu2A(
                        envelope * param->getg2mu() * param->getg2mu()
                        / param->getg() / param->getg());
                    lat->cells[localpos]->setg2mu2B(
                        envelope * param->getg2mu() * param->getg2mu()
                        / param->getg() / param->getg());
                }
            }
        } else {
            for (int ix = 0; ix < N; ix++)  // loop over all positions
            {
                for (int iy = 0; iy < N; iy++) {
                    int localpos = ix * N + iy;
                    lat->cells[localpos]->setg2mu2A(
                        param->getg2mu() * param->getg2mu() / param->getg()
                        / param->getg());
                    lat->cells[localpos]->setg2mu2B(
                        param->getg2mu() * param->getg2mu() / param->getg()
                        / param->getg());
                }
            }
        }
        param->setSuccess(1);
        cout << "constant color charge density set" << endl;
        return;
    }

#pragma omp parallel for
    for (int ipos = 0; ipos < N * N; ipos++) {
        lat->cells[ipos]->setg2mu2A(0.);
        lat->cells[ipos]->setg2mu2B(0.);
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
    vector<double> x_array, y_array, z_array, BGq_array, gauss_array;
    xq1.clear();
    xq2.clear();
    yq1.clear();
    yq2.clear();
    BGq1.clear();
    BGq2.clear();
    gauss1.clear();
    gauss2.clear();
    for (int i = 0; i < A1; i++) {
        int Npartons = 1;
        if (NqFlag > 0) {
            samplePartonPositions(
                param, random, x_array, y_array, z_array, BGq_array);
            // if (param->getShiftConstituentQuarkProtonOrigin())
            // Move center of mass to the origin
            // Note that 1607.01711 this is not done, so parameters quoted
            // in that paper can't be used if this is done
            xq1.push_back(x_array);
            yq1.push_back(y_array);
            BGq1.push_back(BGq_array);
            Npartons = std::max(1, static_cast<int>(x_array.size()));
        }
        sampleQsNormalization(random, param, Npartons, gauss_array);
        gauss1.push_back(gauss_array);
    }

    for (int i = 0; i < A2; i++) {
        int Npartons = 1;
        if (NqFlag > 0) {
            samplePartonPositions(
                param, random, x_array, y_array, z_array, BGq_array);
            xq2.push_back(x_array);
            yq2.push_back(y_array);
            BGq2.push_back(BGq_array);
            Npartons = std::max(1, static_cast<int>(x_array.size()));
        }
        sampleQsNormalization(random, param, Npartons, gauss_array);
        gauss2.push_back(gauss_array);
    }

    // test what a smooth Woods-Saxon would give
    if (param->getUseSmoothNucleus() == 1) {
        cout << "Using smooth nucleus for test purposes. Does not include "
                "deformation."
             << endl;
        double xA, xB;
        double y;
        double T;
        double localpos;
        double normA = 0.;
        double normB = 0.;
        double bb = param->getb();
        double r = 0.;
        for (int ix = 0; ix < N; ix++)  // loop over all positions
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
                if (lat->cells[localpos]->getTpA() < 0.001
                    || lat->cells[localpos]->getTpB() < 0.001) {
                    lat->cells[localpos]->setTpA(0.);
                    lat->cells[localpos]->setTpB(0.);
                }
            }
        }
        for (int ix = 0; ix < N; ix++)  // loop over all positions
        {
            for (int iy = 0; iy < N; iy++) {
                localpos = ix * N + iy;
                lat->cells[localpos]->setTpA(
                    lat->cells[localpos]->getTpA() / normA
                    * glauber->nucleusA1() * hbarc * hbarc);
                lat->cells[localpos]->setTpB(
                    lat->cells[localpos]->getTpB() / normB
                    * glauber->nucleusA2() * hbarc * hbarc);
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
        // Non-smooth nucleus add all T_p's (new in version 1.2)

#pragma omp parallel for
        for (int ipos = 0; ipos < N * N; ipos++) {
            // loop over all positions
            int iy = ipos % N;
            int ix = ipos / N;
            double x = -L / 2. + a * ix;
            double y = -L / 2. + a * iy;

            // nucleus A
            lat->cells[ipos]->setTpA(0.);
            for (int i = 0; i < A1; i++) {
                double xm = nucleusA_.at(i).x;
                double ym = nucleusA_.at(i).y;

                double T = 0.;
                double bp2 = 0.;
                if (param->getUseConstituentQuarkProton() > 0) {
                    for (unsigned int iq = 0; iq < xq1[i].size(); iq++) {
                        bp2 = (xm + xq1[i][iq] - x) * (xm + xq1[i][iq] - x)
                              + (ym + yq1[i][iq] - y) * (ym + yq1[i][iq] - y);
                        bp2 /= hbarc * hbarc;

                        T += exp(-bp2 / (2. * BGq1[i][iq]))
                             / (2. * M_PI * BGq1[i][iq])
                             / (static_cast<double>(xq1[i].size()))
                             * gauss1[i][iq];  // I removed the 2/3 here
                                               // to make it a bit bigger
                    }
                } else {
                    const double BG = param->getBG();
                    double phi = nucleusA_.at(i).phi;

                    bp2 = (xm - x) * (xm - x) + (ym - y) * (ym - y)
                          + xi
                                * pow(
                                    (xm - x) * cos(phi) + (ym - y) * sin(phi),
                                    2.);
                    bp2 /= hbarc * hbarc;
                    T = sqrt(1 + xi) * exp(-bp2 / (2. * BG)) / (2. * M_PI * BG)
                        * gauss1[i][0];  // T_p in this cell for the
                                         // current nucleon
                }
                lat->cells[ipos]->setTpA(
                    lat->cells[ipos]->getTpA()
                    + T / nucleiInAverage);  // add up all T_p
            }

            // nucleus B
            lat->cells[ipos]->setTpB(0.);
            for (int i = 0; i < A2; i++) {
                double xm = nucleusB_.at(i).x;
                double ym = nucleusB_.at(i).y;

                double T = 0.;
                double bp2 = 0.;
                if (param->getUseConstituentQuarkProton() > 0) {
                    T = 0.;
                    for (unsigned int iq = 0; iq < xq2[i].size(); iq++) {
                        bp2 = (xm + xq2[i][iq] - x) * (xm + xq2[i][iq] - x)
                              + (ym + yq2[i][iq] - y) * (ym + yq2[i][iq] - y);
                        bp2 /= hbarc * hbarc;

                        T += exp(-bp2 / (2. * BGq2[i][iq]))
                             / (2. * M_PI * BGq2[i][iq])
                             / (static_cast<double>(xq2[i].size()))
                             * gauss2[i][iq];
                    }
                } else {
                    const double BG = param->getBG();
                    double phi = nucleusB_.at(i).phi;

                    bp2 = (xm - x) * (xm - x) + (ym - y) * (ym - y)
                          + xi
                                * pow(
                                    (xm - x) * cos(phi) + (ym - y) * sin(phi),
                                    2.);
                    bp2 /= hbarc * hbarc;

                    T = sqrt(1 + xi) * exp(-bp2 / (2. * BG)) / (2. * M_PI * BG)
                        * gauss2[i][0];  // T_p in this cell for the
                                         // current nucleon
                }

                lat->cells[ipos]->setTpB(
                    lat->cells[ipos]->getTpB()
                    + T / nucleiInAverage);  // add up all T_p
            }
        }
    }

// get Q_s^2 (and from that g^2mu^2) for a given \sum T_p and Y
#pragma omp parallel for
    for (int ipos = 0; ipos < N * N; ipos++) {
        // loop over all positions
        int ix = ipos / N;
        int iy = ipos % N;

        double QsA = 1;
        double QsB = 1;
        int check = 0;
        double distanceA = 0;
        double distanceB = 0;

        if (param->getUseSmoothNucleus() == 1) {
            check = 2;
        } else {
            const double BG = param->getBG();

            // cut proton at a radius of rmax [fm] (about twice the
            // gluonic radius to be generous)

            if (log(2 * M_PI * BG * lat->cells[ipos]->getTpA()) < 0.) {
                if (isinf(log(2 * M_PI * BG * lat->cells[ipos]->getTpA())) == 1)
                    distanceA = param->getRmax() + 1.;
                else
                    distanceA =
                        sqrt(
                            -2. * BG
                            * log(2 * M_PI * BG * lat->cells[ipos]->getTpA()))
                        * hbarc;
                // cout << log(2 * M_PI * BG *
                // lat->cells[ipos]->getTpA()) << endl;
            } else {
                distanceA = 0.;
            }

            if (log(2 * M_PI * BG * lat->cells[ipos]->getTpB()) < 0.) {
                if (isinf(log(2 * M_PI * BG * lat->cells[ipos]->getTpB())) == 1)
                    distanceB = param->getRmax() + 1.;
                else
                    distanceB =
                        sqrt(
                            -2. * BG
                            * log(2 * M_PI * BG * lat->cells[ipos]->getTpB()))
                        * hbarc;
            } else
                distanceB = 0.;

            if (distanceA < param->getRmax()) {
                check = 1;
            }

            if (distanceB < param->getRmax() && check == 1) {
                check = 2;
            }
        }

        if (param->getUseJIMWLK() == 1) {
            // always assgin color charge density for whole lattice
            check = 2;
        }

        double exponent = 5.6;  // see 1212.2974 Eq. (17)
        double xVal = 0.;
        if (check == 2) {
            if (param->getUseFluctuatingx() == 1) {
                double localrapidity = rapidityA;
                double yIn = rapidityA;
                double Ydeviation = 10000;
                // iterative loops here to determine the fluctuating Y
                while (std::abs(Ydeviation) > 0.001) {
                    if (localrapidity >= 0) {
                        QsA = sqrt(getNuclearQs2(
                            lat->cells[ipos]->getTpA(), localrapidity));
                    } else {
                        xVal = QsA * param->getxFromThisFactorTimesQs()
                               / param->getRoots() * exp(yIn);
                        if (xVal == 0)
                            QsA = 0.;
                        else
                            QsA = sqrt(getNuclearQs2(
                                      lat->cells[ipos]->getTpA(), 0.))
                                  * sqrt(
                                      pow((1 - xVal) / (1 - 0.01), exponent)
                                      * pow((0.01 / xVal), 0.2));
                    }
                    if (QsA == 0) {
                        Ydeviation = 0;
                        lat->cells[ipos]->setg2mu2A(0.);
                    } else {
                        // nucleus A
                        lat->cells[ipos]->setg2mu2A(
                            QsA * QsA / param->getQsmuRatio()
                            / param->getQsmuRatio() * a * a / hbarc / hbarc
                            / param->getg()
                            / param->getg());  // lattice units? check

                        Ydeviation =
                            localrapidity
                            - log(
                                0.01
                                / (QsA * param->getxFromThisFactorTimesQs()
                                   / param->getRoots() * exp(yIn)));
                        localrapidity =
                            log(0.01
                                / (QsA * param->getxFromThisFactorTimesQs()
                                   / param->getRoots() * exp(yIn)));
                    }
                }
                if (lat->cells[ipos]->getg2mu2A()
                    != lat->cells[ipos]->getg2mu2A()) {
                    lat->cells[ipos]->setg2mu2A(0.);
                }

                localrapidity = rapidityB;
                yIn = rapidityB;
                Ydeviation = 10000;
                while (std::abs(Ydeviation) > 0.001) {
                    if (localrapidity >= 0)
                        QsB = sqrt(getNuclearQs2(
                            lat->cells[ipos]->getTpB(), localrapidity));
                    else {
                        xVal = QsB * param->getxFromThisFactorTimesQs()
                               / param->getRoots() * exp(-yIn);
                        if (xVal == 0)
                            QsB = 0.;
                        else
                            QsB = sqrt(getNuclearQs2(
                                      lat->cells[ipos]->getTpB(), 0.))
                                  * sqrt(
                                      pow((1 - xVal) / (1 - 0.01), exponent)
                                      * pow((0.01 / xVal), 0.2));
                    }
                    if (QsB == 0) {
                        Ydeviation = 0;
                        lat->cells[ipos]->setg2mu2B(0.);
                    } else {
                        // nucleus B
                        lat->cells[ipos]->setg2mu2B(
                            QsB * QsB / param->getQsmuRatioB()
                            / param->getQsmuRatioB() * a * a / hbarc / hbarc
                            / param->getg() / param->getg());
                        Ydeviation =
                            localrapidity
                            - log(
                                0.01
                                / (QsB * param->getxFromThisFactorTimesQs()
                                   / param->getRoots() * exp(-yIn)));
                        localrapidity =
                            log(0.01
                                / (QsB * param->getxFromThisFactorTimesQs()
                                   / param->getRoots() * exp(-yIn)));
                    }
                }
                if (lat->cells[ipos]->getg2mu2B()
                    != lat->cells[ipos]->getg2mu2B()) {
                    lat->cells[ipos]->setg2mu2B(0.);
                }
            } else {
                // nucleus A
                lat->cells[ipos]->setg2mu2A(
                    getNuclearQs2(lat->cells[ipos]->getTpA(), rapidityA)
                    / param->getQsmuRatio() / param->getQsmuRatio() * a * a
                    / hbarc / hbarc / param->getg()
                    / param->getg());  // lattice units? check

                // nucleus B
                lat->cells[ipos]->setg2mu2B(
                    getNuclearQs2(lat->cells[ipos]->getTpB(), rapidityB)
                    / param->getQsmuRatioB() / param->getQsmuRatioB() * a * a
                    / hbarc / hbarc / param->getg() / param->getg());
            }
        }
    }
    messager.info("Color charge densities for nucleus A and B set. ");
}

// This function compute the collision geometry quantities, such as
// Npart, Ncoll, averageQs, etc.
void Init::computeCollisionGeometryQuantities(Lattice *lat, Parameters *param) {
    const double d2 = param->getSigmaNN() / (M_PI * 10.);  // in fm^2
    const double b = param->getb();
    const double phiRP = param->getPhiRP();
    double averageQs = 0.;
    double averageQs2 = 0.;
    double averageQs2Avg = 0.;
    double averageQs2min = 0.;
    double averageQs2min2 = 0.;
    int Npart = 0;
    int Ncoll = 0;

    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();
    const double L = param->getL();
    const int N = param->getSize();
    const double a = L / N;  // lattice spacing in fm

    // Determine Npart, Ncoll. Do this only during the first stage, as in
    // the 2nd stage nuclei are shifted to b=0
    if (param->getUseSmoothNucleus() == 0) {
        stringstream strNcoll_name;
        strNcoll_name << "NcollList" << param->getEventId() << ".dat";
        string Ncoll_name;
        Ncoll_name = strNcoll_name.str();

        ofstream foutNcoll(Ncoll_name.c_str(), std::ios::out);

        if (param->getGaussianWounding() == 0) {
            for (int i = 0; i < A1; i++) {
                for (int j = 0; j < A2; j++) {
                    double dx =
                        (nucleusB_.at(j).x - nucleusA_.at(i).x
                         - b * cos(phiRP));
                    double dy =
                        (nucleusB_.at(j).y - nucleusA_.at(i).y
                         - b * sin(phiRP));
                    double dij = dx * dx + dy * dy;
                    if (dij < d2) {
                        foutNcoll
                            << (nucleusB_.at(j).x + nucleusA_.at(i).x) / 2.
                            << " "
                            << (nucleusB_.at(j).y + nucleusA_.at(i).y) / 2.
                            << endl;
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
                    double dx =
                        (nucleusB_.at(j).x - nucleusA_.at(i).x
                         - b * cos(phiRP));
                    double dy =
                        (nucleusB_.at(j).y - nucleusA_.at(i).y
                         - b * sin(phiRP));
                    double dij = dx * dx + dy * dy;

                    p = G * exp(-G * dij / d2);  // Gaussian profile

                    ran = random_ptr_->genrand64_real1();

                    if (ran < p) {
                        foutNcoll
                            << (nucleusB_.at(j).x + nucleusA_.at(i).x) / 2.
                            << " "
                            << (nucleusB_.at(j).y + nucleusA_.at(i).y) / 2.
                            << endl;
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

        ofstream foutNpart(Npart_name.c_str(), std::ios::out);

        for (int i = 0; i < A1; i++) {
            foutNpart << nucleusA_.at(i).x + b / 2. * cos(phiRP) << " "
                      << nucleusA_.at(i).y + b / 2. * sin(phiRP) << " "
                      << nucleusA_.at(i).proton << " "
                      << nucleusA_.at(i).collided << endl;
        }
        foutNpart << endl;
        for (int i = 0; i < A2; i++) {
            foutNpart << nucleusB_.at(i).x - b / 2. * cos(phiRP) << " "
                      << nucleusB_.at(i).y - b / 2. * sin(phiRP) << " "
                      << nucleusB_.at(i).proton << " "
                      << nucleusB_.at(i).collided << endl;
        }
        foutNpart.close();

        // in p+p assume that they collided in any case
        if (A1 == 1 && A2 == 1) {
            nucleusB_.at(0).collided = 1;
            nucleusA_.at(0).collided = 1;
        }

        Npart = 0;
        for (int i = 0; i < A1; i++) {
            if (nucleusA_.at(i).collided == 1) {
                Npart++;
            }
        }

        for (int i = 0; i < A2; i++) {
            if (nucleusB_.at(i).collided == 1) {
                Npart++;
            }
        }

        param->setNpart(Npart);

        if (param->getUseFixedNpart() != 0
            && Npart != param->getUseFixedNpart()) {
            cout << "current Npart = " << Npart
                 << " != " << param->getUseFixedNpart() << endl;
            param->setSuccess(0);
            return;
        }
    } else {
        // Smooth nucleus
        Npart = 2;
        Ncoll = 2;
        param->setNpart(Npart);
    }

    int count = 0;
    double Tpp = 0.;
    for (int ipos = 0; ipos < N * N; ipos++) {
        // loop over all positions
        int check = 0;
        int ix = ipos / N;
        int iy = ipos % N;
        double x = -L / 2. + a * ix;
        double y = -L / 2. + a * iy;

        double xA = x - b / 2. * cos(phiRP);
        double yA = y - b / 2. * sin(phiRP);
        double xB = x + b / 2. * cos(phiRP);
        double yB = y + b / 2. * sin(phiRP);

        int ixA = static_cast<int>((xA + L / 2.) / a);
        int iyA = static_cast<int>((yA + L / 2.) / a);
        int ixB = static_cast<int>((xB + L / 2.) / a);
        int iyB = static_cast<int>((yB + L / 2.) / a);

        int posA = ixA * N + iyA;
        int posB = ixB * N + iyB;

        double g2mu2A = 0;
        double TpA = 0;
        if (posA > 0 && posA < N * N) {
            g2mu2A = lat->cells[posA]->getg2mu2A();
            TpA = lat->cells[posA]->getTpA();
        }

        double g2mu2B = 0;
        double TpB = 0;
        if (posB > 0 && posB < N * N) {
            g2mu2B = lat->cells[posB]->getg2mu2B();
            TpB = lat->cells[posB]->getTpB();
        }

        if (g2mu2B >= g2mu2A) {
            averageQs2min2 += g2mu2A * param->getQsmuRatio()
                              * param->getQsmuRatio() / a / a * hbarc * hbarc
                              * param->getg() * param->getg();
        } else {
            averageQs2min2 += g2mu2B * param->getQsmuRatioB()
                              * param->getQsmuRatioB() / a / a * hbarc * hbarc
                              * param->getg() * param->getg();
        }

        for (int i = 0; i < A1; i++) {
            double xm = nucleusA_.at(i).x + b / 2. * cos(phiRP);
            double ym = nucleusA_.at(i).y + b / 2. * sin(phiRP);
            double r = sqrt((x - xm) * (x - xm) + (y - ym) * (y - ym));

            if (r < sqrt(0.1 * param->getSigmaNN() / M_PI)
                && nucleusA_.at(i).collided == 1) {
                check = 1;
            }
        }

        for (int i = 0; i < A2; i++) {
            double xm = nucleusB_.at(i).x - b / 2. * cos(phiRP);
            double ym = nucleusB_.at(i).y - b / 2. * sin(phiRP);
            double r = sqrt((x - xm) * (x - xm) + (y - ym) * (y - ym));

            if (r < sqrt(0.1 * param->getSigmaNN() / M_PI)
                && nucleusB_.at(i).collided == 1 && check == 1) {
                check = 2;
            }
        }

        if (check == 2) {
            if (g2mu2B > g2mu2A) {
                averageQs += sqrt(
                    g2mu2B * param->getQsmuRatioB() * param->getQsmuRatioB() / a
                    / a * hbarc * hbarc * param->getg() * param->getg());
                averageQs2 += g2mu2B * param->getQsmuRatioB()
                              * param->getQsmuRatioB() / a / a * hbarc * hbarc
                              * param->getg() * param->getg();
                averageQs2min += g2mu2A * param->getQsmuRatio()
                                 * param->getQsmuRatio() / a / a * hbarc * hbarc
                                 * param->getg() * param->getg();
            } else {
                averageQs += sqrt(
                    g2mu2A * param->getQsmuRatio() * param->getQsmuRatio() / a
                    / a * hbarc * hbarc * param->getg() * param->getg());
                averageQs2 += g2mu2A * param->getQsmuRatio()
                              * param->getQsmuRatio() / a / a * hbarc * hbarc
                              * param->getg() * param->getg();
                averageQs2min += g2mu2B * param->getQsmuRatioB()
                                 * param->getQsmuRatioB() / a / a * hbarc
                                 * hbarc * param->getg() * param->getg();
            }
            averageQs2Avg +=
                (g2mu2B * param->getQsmuRatioB() * param->getQsmuRatioB()
                 + g2mu2A * param->getQsmuRatio() * param->getQsmuRatio())
                / 2. / a / a * hbarc * hbarc * param->getg() * param->getg();
            count++;
        }
        // compute T_pp
        Tpp += TpA * TpB * a * a / hbarc / hbarc / hbarc
               / hbarc;  // now this quantity is in fm^-2
                         // remember: Tp is in GeV^2
    }

    if (count == 0) {
        param->setAverageQs(0.);
        param->setAverageQsAvg(0.);
        param->setAverageQsmin(0.);
        param->setTpp(Tpp);
        param->setSuccess(0);
        cout << "**** Rejected event - no overlap region (count=0)." << endl;
        return;
    }

    averageQs /= static_cast<double>(count) + 1e-16;
    averageQs2 /= static_cast<double>(count) + 1e-16;
    averageQs2Avg /= static_cast<double>(count) + 1e-16;
    averageQs2min /= static_cast<double>(count) + 1e-16;

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
    cout << "Q_s^2(min) S_T = "
         << averageQs2min * a * a / hbarc / hbarc * static_cast<double>(count)
         << endl;
    cout << "Q_s^2(min) S_T = " << averageQs2min2 * a * a / hbarc / hbarc
         << endl;

    cout << "Area = " << a * a * count << " fm^2" << endl;

    cout << "Average Qs(max) = " << param->getAverageQs() << " GeV" << endl;
    cout << "Average Qs(avg) = " << param->getAverageQsAvg() << " GeV" << endl;
    cout << "Average Qs(min) = " << param->getAverageQsmin() << " GeV" << endl;

    cout << "resulting Y(Qs(max)*" << param->getxFromThisFactorTimesQs()
         << ") = "
         << log(0.01
                / (param->getAverageQs() * param->getxFromThisFactorTimesQs()
                   / param->getRoots()))
         << endl;
    cout << "resulting Y(Qs(avg)*" << param->getxFromThisFactorTimesQs()
         << ") = "
         << log(0.01
                / (param->getAverageQsAvg() * param->getxFromThisFactorTimesQs()
                   / param->getRoots()))
         << endl;
    cout << "resulting Y(Qs(min)*" << param->getxFromThisFactorTimesQs()
         << ") =  "
         << log(0.01
                / (param->getAverageQsmin() * param->getxFromThisFactorTimesQs()
                   / param->getRoots()))
         << endl;

    double alphas = 0.;
    if (param->getRunningCoupling() && param->getRunWithkt() == 0) {
        if (param->getRunWithQs() == 2) {
            cout << "running with " << param->getRunWithThisFactorTimesQs()
                 << " Q_s(max)" << endl;
            alphas = 12. * M_PI
                     / ((27.) * 2.
                        * log(
                            param->getRunWithThisFactorTimesQs()
                            * param->getAverageQs() / 0.2));  // 3 flavors
            cout << "alpha_s(" << param->getRunWithThisFactorTimesQs()
                 << " Qs_max)=" << alphas << endl;
        } else if (param->getRunWithQs() == 0) {
            cout << "running with " << param->getRunWithThisFactorTimesQs()
                 << " Q_s(min)" << endl;
            alphas = 12. * M_PI
                     / ((27.) * 2.
                        * log(
                            param->getRunWithThisFactorTimesQs()
                            * param->getAverageQsmin() / 0.2));  // 3 flavors
            cout << "alpha_s(" << param->getRunWithThisFactorTimesQs()
                 << " Qs_min)=" << alphas << endl;
        } else if (param->getRunWithQs() == 1) {
            cout << "running with " << param->getRunWithThisFactorTimesQs()
                 << " <Q_s>" << endl;
            alphas = 12. * M_PI
                     / ((27.) * 2.
                        * log(
                            param->getRunWithThisFactorTimesQs()
                            * param->getAverageQsAvg() / 0.2));  // 3 flavors
            cout << "alpha_s(" << param->getRunWithThisFactorTimesQs()
                 << " <Qs>)=" << alphas << endl;
        }
    } else if (param->getRunningCoupling() && param->getRunWithkt() == 1) {
        cout << "Multiplicity with running alpha_s(k_T)" << endl;
    } else {
        cout << "Using fixed alpha_s" << endl;
        alphas = param->getg() * param->getg() / 4. / M_PI;
    }

    if (param->getAverageQs() > 0 && param->getAverageQsAvg() > 0
        && averageQs2 > 0 && param->getAverageQsmin() > 0 && averageQs2Avg > 0
        && alphas > 0 && Npart >= 2
        && averageQs2min2 * a * a / hbarc / hbarc > param->getMinimumQs2ST()) {
        param->setSuccess(1);

        stringstream strup_name;
        strup_name << "usedParameters" << param->getEventId() << ".dat";
        string up_name;
        up_name = strup_name.str();

        ofstream fout1(up_name.c_str(), std::ios::app);
        fout1 << " " << endl;
        fout1 << " Output by setColorChargeDensity in Init.cpp: " << endl;
        fout1 << " " << endl;
        fout1 << "b = " << param->getb() << " fm" << endl;
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
            fout1 << "using fixed coupling alpha_s=" << param->getalphas()
                  << endl;
        fout1.close();
    } else {
        param->setSuccess(0);
    }
    if (averageQs2min2 * a * a / hbarc / hbarc < param->getMinimumQs2ST()) {
        cout << " **** Rejected event - Qsmin^2 S_T="
             << averageQs2min2 * a * a / hbarc / hbarc << " too small ( < "
             << param->getMinimumQs2ST() << ")." << endl;
    }

    param->setalphas(alphas);

    stringstream strNEst_name;
    strNEst_name << "NgluonEstimators" << param->getEventId() << ".dat";
    string NEst_name;
    NEst_name = strNEst_name.str();

    ofstream foutNEst(NEst_name.c_str(), std::ios::out);

    foutNEst << "#Q_s^2(min) S_T  " << "Q_s^2(avg) S_T  "
             << "Q_s^2(max) S_T "
             << " Q_s^2(min) S_T Log^2( Q_s^2(max) / Q_s^2(min))  " << endl;
    foutNEst << averageQs2min2 * a * a / hbarc / hbarc << "         "
             << averageQs2Avg * a * a / hbarc / hbarc
                    * static_cast<double>(count)
             << "         "
             << averageQs2 * a * a / hbarc / hbarc * static_cast<double>(count)
             << "         "
             << averageQs2min2 * a * a / hbarc / hbarc
                    * pow(
                        log(averageQs2 * static_cast<double>(count)
                            / averageQs2min2),
                        2.)
             << endl;

    foutNEst.close();
}

namespace {
void writeInitialWilsonTrainingData(Lattice *lat, Parameters *param) {
    const int N = param->getSize();
    const int Nc = param->getNc();
    const double L = param->getL();
    const double a = L / static_cast<double>(N);

    // Payload: [beam, real_or_imag, x, y, row, col], C-order.
    constexpr int nBeams = 2;
    const std::size_t matrixElements =
        static_cast<std::size_t>(N) * N * Nc * Nc;
    std::vector<float> payload(
        static_cast<std::size_t>(nBeams) * 2 * matrixElements);

    for (int beam = 0; beam < nBeams; ++beam) {
        const std::size_t realOffset =
            static_cast<std::size_t>(2 * beam) * matrixElements;
        const std::size_t imagOffset = realOffset + matrixElements;
        for (int x = 0; x < N; ++x) {
            for (int y = 0; y < N; ++y) {
                const int pos = x * N + y;
                const Matrix &matrix = (beam == 0) ? lat->U[pos] : lat->U2[pos];
                const std::complex<double> *elements = matrix.data();
                const std::size_t siteOffset =
                    static_cast<std::size_t>(pos) * Nc * Nc;
                for (int row = 0; row < Nc; ++row) {
                    for (int col = 0; col < Nc; ++col) {
                        const std::size_t element =
                            static_cast<std::size_t>(row) * Nc + col;
                        payload[realOffset + siteOffset + element] =
                            static_cast<float>(elements[element].real());
                        payload[imagOffset + siteOffset + element] =
                            static_cast<float>(elements[element].imag());
                    }
                }
            }
        }
    }

    const std::uint16_t endianProbe = 1;
    if (*reinterpret_cast<const unsigned char *>(&endianProbe) != 1) {
        throw std::runtime_error(
            "writeInitialWilsonTrainingData requires a little-endian host");
    }

    std::stringstream metadata;
    metadata << std::setprecision(17)
             << "{\"format\":\"ipglasma-initial-wilson-lines\","
             << "\"version\":1,"
             << "\"dtype\":\"<f4\","
             << "\"shape\":[2,2," << N << "," << N << "," << Nc << "," << Nc
             << "],"
             << "\"axis_order\":[\"beam\",\"complex_part\",\"x\",\"y\","
                "\"row\",\"col\"],"
             << "\"fields\":[\"VA\",\"VB\"],"
             << "\"complex_part\":[\"real\",\"imag\"],"
             << "\"native_site_index\":\"pos=x*N+y\","
             << "\"event_id\":" << param->getEventId() << ","
             << "\"N\":" << N << ","
             << "\"Nc\":" << Nc << ","
             << "\"L_fm\":" << L << ","
             << "\"a_fm\":" << a << ","
             << "\"rapidity\":" << param->getRapidity() << "}";
    const std::string metadataString = metadata.str();

    std::stringstream filename;
    filename << "initialWilsonLines" << param->getEventId() << ".ipgw";
    std::ofstream output(
        filename.str().c_str(),
        std::ios::out | std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "could not open initial-Wilson snapshot " + filename.str());
    }
    const char magic[8] = {'I', 'P', 'G', 'W', 'I', 'L', '1', '\0'};
    const std::uint64_t metadataBytes =
        static_cast<std::uint64_t>(metadataString.size());
    output.write(magic, sizeof(magic));
    output.write(
        reinterpret_cast<const char *>(&metadataBytes), sizeof(metadataBytes));
    output.write(metadataString.data(), metadataString.size());
    output.write(
        reinterpret_cast<const char *>(payload.data()),
        static_cast<std::streamsize>(payload.size() * sizeof(float)));
    output.close();
    if (!output) {
        throw std::runtime_error(
            "failed while writing initial-Wilson snapshot " + filename.str());
    }
    std::cout << "Wrote incoming Wilson lines to " << filename.str()
              << std::endl;
}
}  // namespace

void Init::setV(Lattice *lat, Parameters *param, Random *random) {
    IPG_PROFILE_SCOPE("initialization.wilson_lines");
    messager.info("Setting Wilson lines ...");
    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();

    const double d2 = param->getSigmaNN() / (M_PI * 10.);  // in fm^2
    const int N = param->getSize();
    const int Ny = param->getNy();
    const int sites = N * N;
    const int nn[2] = {N, N};
    const double L = param->getL();
    const double a = L / N;  // lattice spacing in fm
    const double m = param->getm() * a / hbarc;
    const double g = param->getg();
    const double invNy = 1. / static_cast<double>(Ny);
    double UVdamp = param->getUVdamp();  // GeV^-1
    UVdamp = UVdamp / a * hbarc;

    // The lattice Poisson/UV kernel depends only on transverse momentum and
    // run parameters.  The historical implementation recomputed the same
    // sin/sqrt/exp expressions for every longitudinal sheet of both nuclei.
    std::vector<double> momentumKernel(static_cast<std::size_t>(sites));
#pragma omp parallel for
    for (int pos = 0; pos < sites; ++pos) {
        const int i = pos / N;
        const int j = pos - i * N;
        const double kx =
            2. * M_PI
            * (-0.5 + static_cast<double>(i) / static_cast<double>(N));
        const double ky =
            2. * M_PI
            * (-0.5 + static_cast<double>(j) / static_cast<double>(N));
        const double sx = sin(kx / 2.);
        const double sy = sin(ky / 2.);
        const double kt2 = 4. * (sx * sx + sy * sy);

        if (m == 0.) {
            momentumKernel[static_cast<std::size_t>(pos)] =
                (kt2 != 0.) ? 1. / kt2 : 0.;
        } else {
            momentumKernel[static_cast<std::size_t>(pos)] =
                (1. / (kt2 + m * m)) * exp(-sqrt(kt2) * UVdamp);
        }
    }

    complex<double> **rhoACoeff = new complex<double> *[Nc2m1_];
    for (int i = 0; i < Nc2m1_; i++) {
        rhoACoeff[i] = new complex<double>[sites];
    }

    auto applyMomentumKernel = [&]() {
#pragma omp parallel for
        for (int n = 0; n < Nc2m1_; ++n) {
            complex<double> *rho = rhoACoeff[n];
            for (int pos = 0; pos < sites; ++pos) {
                rho[pos] *= momentumKernel[static_cast<std::size_t>(pos)];
            }
        }
    };

    // Reuse the bulk Gaussian buffers for every longitudinal sheet.  The
    // linear ordering matches the historical pos-major/color-minor Gauss()
    // call sequence exactly.
    std::vector<double> gaussianField(
        static_cast<std::size_t>(sites) * static_cast<std::size_t>(Nc2m1_));
    std::vector<double> gaussianScratch;
    gaussianScratch.reserve(5 * ((gaussianField.size() + 1) / 2));

    // g2mu2 and Ny are fixed throughout Wilson-line construction.  Cache the
    // color-independent site scale once for each nucleus instead of repeating
    // the same sqrt in every longitudinal sheet.
    std::vector<double> colorChargeScaleA(static_cast<std::size_t>(sites));
    std::vector<double> colorChargeScaleB(static_cast<std::size_t>(sites));
#pragma omp parallel for
    for (int pos = 0; pos < sites; ++pos) {
        colorChargeScaleA[static_cast<std::size_t>(pos)] =
            g * sqrt(lat->cells[pos]->getg2mu2A() * invNy);
        colorChargeScaleB[static_cast<std::size_t>(pos)] =
            g * sqrt(lat->cells[pos]->getg2mu2B() * invNy);
    }

    auto fillColorCharge = [&](const std::vector<double> &scale) {
        {
            IPG_PROFILE_SCOPE("initialization.wilson_random.gauss");
            random->GaussBulk(
                gaussianField.data(), gaussianField.size(), gaussianScratch);
        }
        {
            IPG_PROFILE_SCOPE("initialization.wilson_random.scale");
#pragma omp parallel for
            for (int pos = 0; pos < sites; ++pos) {
                const double localScale = scale[static_cast<std::size_t>(pos)];
                const std::size_t base = static_cast<std::size_t>(pos)
                                         * static_cast<std::size_t>(Nc2m1_);
                for (int n = 0; n < Nc2m1_; ++n) {
                    rhoACoeff[n][pos] =
                        localScale
                        * gaussianField[base + static_cast<std::size_t>(n)];
                }
            }
        }
    };

    // loop over longitudinal direction for nucleus A
    for (int k = 0; k < Ny; k++) {
        {
            IPG_PROFILE_SCOPE("initialization.wilson_random");
            fillColorCharge(colorChargeScaleA);
        }

        for (int n = 0; n < Nc2m1_; n++) {
            fft.fftnComplex(rhoACoeff[n], rhoACoeff[n], nn, 1);
        }

        {
            IPG_PROFILE_SCOPE("initialization.wilson_poisson");
            applyMomentumKernel();
        }

        for (int n = 0; n < Nc2m1_; n++) {
            fft.fftnComplex(rhoACoeff[n], rhoACoeff[n], nn, -1);
        }

        {
            IPG_PROFILE_SCOPE("initialization.wilson_exponent");
#pragma omp parallel
            {
                std::vector<double> in(Nc2m1_, 0.);
                Matrix temp(Nc_, 1.);
                Matrix tempNew(Nc_, 0.);

#pragma omp for
                for (int pos = 0; pos < sites; pos++) {
                    for (int aa = 0; aa < Nc2m1_; aa++) {
                        // expmCoeff calculates exp(i in[a] t[a]), so multiply
                        // by -1 (not -i).
                        in[aa] = -(rhoACoeff[aa][pos]).real();
                    }
                    tempNew = getUfromExponent(in);
                    temp = tempNew * lat->U[pos];
                    lat->U[pos] = temp;
                }
            }
        }
    }

    // loop over longitudinal direction for nucleus B
    for (int k = 0; k < Ny; k++) {
        {
            IPG_PROFILE_SCOPE("initialization.wilson_random");
            fillColorCharge(colorChargeScaleB);
        }

        for (int n = 0; n < Nc2m1_; n++) {
            fft.fftnComplex(rhoACoeff[n], rhoACoeff[n], nn, 1);
        }

        {
            IPG_PROFILE_SCOPE("initialization.wilson_poisson");
            applyMomentumKernel();
        }

        for (int n = 0; n < Nc2m1_; n++) {
            fft.fftnComplex(rhoACoeff[n], rhoACoeff[n], nn, -1);
        }

        // The old nucleus-B block had its omp parallel directive commented
        // out, leaving this expensive exponential/multiply pass effectively
        // serial.  Match the nucleus-A implementation.
        {
            IPG_PROFILE_SCOPE("initialization.wilson_exponent");
#pragma omp parallel
            {
                std::vector<double> in(Nc2m1_, 0.);
                Matrix temp(Nc_, 1.);
                Matrix tempNew(Nc_, 0.);

#pragma omp for
                for (int pos = 0; pos < sites; pos++) {
                    for (int aa = 0; aa < Nc2m1_; aa++) {
                        in[aa] = -(rhoACoeff[aa][pos]).real();
                    }
                    tempNew = getUfromExponent(in);
                    temp = tempNew * lat->U2[pos];
                    lat->U2[pos] = temp;
                }
            }
        }
    }

    for (int ic = 0; ic < Nc2m1_; ic++) {
        delete[] rhoACoeff[ic];
    }
    delete[] rhoACoeff;
    if (param->getWriteOutputs() == 5) {
        writeInitialWilsonTrainingData(lat, param);
    }

    // output U
    if (param->getWriteWilsonLines() > 0 && param->getSaveSnapshots()) {
        std::stringstream ss;
        ss << "Initial_x_";
        if (param->getUseJIMWLK()) {
            ss << param->getJimwlk_x0() << "_";
        } else {
            ss << "0.001" << "_";
        }
        std::string wilsonfileHeader = ss.str();
        lat->WriteWilsonLines(wilsonfileHeader, param, 1);  // nucleus A
        lat->WriteWilsonLines(wilsonfileHeader, param, 2);  // nucleus B
    }

    messager << " Wilson lines V_A and V_B set on rank " << param->getMPIRank()
             << ". ";
    messager.flush("info");
}

void Init::readVFromFile(Lattice *lat, Parameters *param, int format) {
    IPG_PROFILE_SCOPE("initialization.read_wilson_lines");
    // format 1 = plain text, 2 = binary

    if (format > 2 or format < 1) {
        messager << "Unknown format " << format
                 << " when reading the initial Wilson lines, supported "
                    "formats: 1,2";
        messager.flush("info");
        exit(1);
    }

    stringstream strVOne_name;
    // strVOne_name << "V1-" << param->getMPIRank() << ".txt";
    strVOne_name << "V-"
                 << param->getEventId()
                        + 2 * param->getSeed() * param->getMPISize();
    if (format == 1) strVOne_name << ".txt";
    string VOne_name;
    VOne_name = strVOne_name.str();

    stringstream strVTwo_name;
    // strVTwo_name << "V2-" << param->getMPIRank() << ".txt";
    strVTwo_name << "V-"
                 << param->getEventId()
                        + (1 + 2 * param->getSeed()) * param->getMPISize();
    if (format == 1) strVTwo_name << ".txt";
    string VTwo_name;
    VTwo_name = strVTwo_name.str();

    messager << "Reading Wilson lines from files " << VOne_name << " and "
             << VTwo_name;
    messager.flush("info");

    if (format == 1) {
        int N = param->getSize();

        double L = param->getL();
        double a = L / static_cast<double>(N);

        int nn[2];
        nn[0] = N;
        nn[1] = N;

        Matrix temp(Nc_, 1.);

        double Re[9], Im[9];
        double dummy;

        ifstream finV1(VOne_name.c_str(), std::ios::in);

        if (!finV1) {
            messager << "File " << VOne_name << " not found. Exiting.";
            messager.flush("info");
            exit(1);
        }

        messager << "Reading Wilson line from file " << VOne_name << " ...";

        // set V for nucleus A
        for (int i = 0; i < nn[0]; i++) {
            for (int j = 0; j < nn[1]; j++) {
                finV1 >> dummy >> dummy >> Re[0] >> Im[0] >> Re[1] >> Im[1]
                    >> Re[2] >> Im[2] >> Re[3] >> Im[3] >> Re[4] >> Im[4]
                    >> Re[5] >> Im[5] >> Re[6] >> Im[6] >> Re[7] >> Im[7]
                    >> Re[8] >> Im[8];

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
                a = L / static_cast<double>(N);

                double xtemp = a * i - bb / 2.;
                int ix = xtemp / a;

                if (ix < 0) continue;

                int pos = ix * N + j;
                lat->U[pos] = (temp);
            }
        }

        finV1.close();

        ifstream finV2(VTwo_name.c_str(), std::ios::in);

        if (!finV2) {
            cerr << "File " << VTwo_name << " not found. Exiting." << endl;
            exit(1);
        }

        cout << "Reading Wilson line from file " << VTwo_name << " ..." << endl;

        // set V for nucleus B
        for (int i = 0; i < nn[0]; i++) {
            for (int j = 0; j < nn[1]; j++) {
                finV2 >> dummy >> dummy >> Re[0] >> Im[0] >> Re[1] >> Im[1]
                    >> Re[2] >> Im[2] >> Re[3] >> Im[3] >> Re[4] >> Im[4]
                    >> Re[5] >> Im[5] >> Re[6] >> Im[6] >> Re[7] >> Im[7]
                    >> Re[8] >> Im[8];

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
                a = L / static_cast<double>(N);

                double xtemp = a * i + bb / 2.;
                int ix = xtemp / a;

                if (ix >= N) continue;

                int pos = ix * N + j;
                lat->U2[pos] = (temp);
            }
        }

        finV2.close();
    } else if (format == 2) {
        std::ifstream InStream;
        InStream.precision(15);
        InStream.open(VOne_name.c_str(), std::ios::in | std::ios::binary);
        int N;
        double L, a, temp;

        Matrix tempM(Nc_, 1.);

        if (!InStream.good()) {
            messager << "File " << VOne_name.c_str() << " does not exist!";
            messager.flush("info");
            exit(1);
        }

        if (InStream.is_open()) {
            // READING IN PARAMETERS
            InStream.read(reinterpret_cast<char *>(&N), sizeof(int));
            InStream.read(reinterpret_cast<char *>(&Nc_), sizeof(int));
            InStream.read(reinterpret_cast<char *>(&L), sizeof(double));
            InStream.read(reinterpret_cast<char *>(&a), sizeof(double));
            InStream.read(reinterpret_cast<char *>(&temp), sizeof(double));

            if (N != param->getSize()) {
                messager << "# ERROR wrong lattice size, data is " << N
                         << " but you have specified " << param->getSize();
                exit(0);
            }
            if (std::abs(L - param->getL()) > 1e-5) {
                messager << "# ERROR grid length, data has " << L
                         << " but you have specified " << param->getL();
                exit(0);
            }

            // READING ACTUAL DATA
            double ValueBuffer;
            int INPUT_CTR = 0;
            double re, im;
            re = 0.;
            im = 0.;

            while (InStream.read(
                reinterpret_cast<char *>(&ValueBuffer), sizeof(double))) {
                if (INPUT_CTR % 2 == 0)  // this is the real part
                {
                    re = ValueBuffer;
                } else  // this is the imaginary part, write then to
                        // variable //
                {
                    im = ValueBuffer;
                    int TEMPINDX = ((INPUT_CTR - 1) / 2);
                    int PositionIndx = TEMPINDX / 9;

                    // shift here by half an impact parameter
                    int iy = PositionIndx / N;
                    int ixIn = PositionIndx - N * iy;

                    double bb = param->getb();
                    a = L / static_cast<double>(N);

                    double xtemp = a * ixIn - bb / 2.;

                    int ix = round(xtemp / a);

                    // cout << ixIn << " " << ix << endl;

                    int MatrixIndx = TEMPINDX - PositionIndx * 9;
                    int j = MatrixIndx / 3;
                    int k = MatrixIndx - j * 3;

                    int indx = N * ix + iy;
                    if (indx >= N * N || indx < 0) {
                        if (bb == 0) {
                            messager << "Warning: datafile " << VOne_name
                                     << " has an element " << indx
                                     << " (iy=" << iy << ", ix=" << ix
                                     << "), but the grid is N=" << N
                                     << ". Element is (" << re << " + " << im
                                     << "i), skipping it";
                            messager.flush("info");
                        }
                        INPUT_CTR++;
                        continue;
                    }
                    lat->U[indx].set(j, k, complex<double>(re, im));
                }
                INPUT_CTR++;
            }

            InStream.close();

            std::ifstream InStream2;
            InStream2.precision(15);
            InStream2.open(VTwo_name.c_str(), std::ios::in | std::ios::binary);
            if (!InStream2.good()) {
                messager << "File " << VTwo_name.c_str() << " does not exist!";
                messager.flush("info");
                exit(1);
            }

            INPUT_CTR = 0;
            if (InStream2.is_open()) {
                // READING IN PARAMETERS
                InStream2.read(reinterpret_cast<char *>(&N), sizeof(int));
                InStream2.read(reinterpret_cast<char *>(&Nc_), sizeof(int));
                InStream2.read(reinterpret_cast<char *>(&L), sizeof(double));
                InStream2.read(reinterpret_cast<char *>(&a), sizeof(double));
                InStream2.read(reinterpret_cast<char *>(&temp), sizeof(double));

                if (N != param->getSize()) {
                    messager << "# ERROR wrong lattice size, data is " << N
                             << " but you have specified " << param->getSize();
                    messager.flush("info");
                    exit(0);
                }
                if (std::abs(L - param->getL()) > 1e-5) {
                    messager << "# ERROR grid length, dataas " << L
                             << " but you have specified " << param->getL();
                    messager.flush("info");
                    exit(0);
                }

                // READING ACTUAL DATA
                while (InStream2.read(
                    reinterpret_cast<char *>(&ValueBuffer), sizeof(double))) {
                    if (INPUT_CTR % 2 == 0)  // this is the real part
                    {
                        re = ValueBuffer;
                    } else  // this is the imaginary part, write then to
                            // variable //
                    {
                        im = ValueBuffer;

                        int TEMPINDX = ((INPUT_CTR - 1) / 2);
                        int PositionIndx = TEMPINDX / 9;

                        // shift here by half an impact parameter
                        int iy = PositionIndx / N;
                        int ixIn = PositionIndx - N * iy;

                        double bb = param->getb();
                        a = L / static_cast<double>(N);

                        double xtemp = a * ixIn + bb / 2.;

                        int ix = round(xtemp / a);

                        //                  if (ixIn != ix)
                        //  cout << ixIn << " " << ix << endl;
                        int MatrixIndx = TEMPINDX - PositionIndx * 9;
                        int j = MatrixIndx / 3;
                        int k = MatrixIndx - j * 3;

                        int indx = N * ix + iy;

                        if (indx >= N * N || indx < 0) {
                            if (bb == 0) {
                                messager << "Warning: datafile " << VTwo_name
                                         << " has an element " << indx
                                         << " (iy=" << iy << ", ix=" << ix
                                         << "), but the grid is N=" << N
                                         << ". Element is (" << re << " + "
                                         << im << "i), skipping it";
                                messager.flush("info");
                            }
                            INPUT_CTR++;
                            continue;
                        }
                        lat->U2[indx].set(j, k, complex<double>(re, im));
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

void Init::sampleImpactParameter(Parameters *param) {
    const double bmin = param->getbmin();
    const double bmax = param->getbmax();
    double b = 0.;
    double xb = random_ptr_->genrand64_real1();
    if (param->getUseNucleus() == 0) {
        // use b=0 fm for the constant g^2 mu case
        messager << "Setting b=0 for constant color charge density case.";
        messager.flush("info");
        b = 0;
    } else {
        if (param->getLinearb() == 1) {
            // use a linear probability distribution for b if we are doing
            // nuclei
            messager << "Sampling linearly distributed b between " << bmin
                     << " and " << bmax << "fm. Found ";
            b = sqrt((bmax * bmax - bmin * bmin) * xb + bmin * bmin);
        } else {
            // use a uniform distribution instead
            messager << "Sampling uniformly distributed b between " << bmin
                     << " and " << bmax << "fm. Found ";
            b = (bmax - bmin) * xb + bmin;
        }
    }
    param->setb(b);
    double phiRP = 0.;
    if (param->getRotateReactionPlane()) {
        phiRP = 2 * M_PI * random_ptr_->genrand64_real2();
    }
    param->setPhiRP(phiRP);
    messager << "b = " << b << " fm, phi_RP = " << phiRP;
    messager.flush("info");

    for (unsigned int i = 0; i < nucleusA_.size(); i++) {
        nucleusA_.at(i).collided = 0;
    }
    for (unsigned int i = 0; i < nucleusB_.size(); i++) {
        nucleusB_.at(i).collided = 0;
    }
}

void Init::init(
    Lattice *lat, Group *group, Parameters *param, Random *random,
    Glauber *glauber, Initialization_method init_method) {
    group_ptr_ = group;
    random_ptr_ = random;

    messager.info("Initializing fields ... ");

    if (param->getUseNucleus() == 0) {
        param->setSuccess(1);
    } else {
        readNuclearQs(param);
    }

    readInNucleusConfigs(
        static_cast<int>(glauber->nucleusA1()), param->getlightNucleusOption(),
        param->getPolarizationProjectile(),
        param->getPolarizationProjectileJz(), nucleonPosArrA_);
    readInNucleusConfigs(
        static_cast<int>(glauber->nucleusA2()), param->getlightNucleusOption(),
        param->getPolarizationTarget(), param->getPolarizationTargetJz(),
        nucleonPosArrB_);

    if (init_method == READ_WLINE_BINARY or init_method == READ_WLINE_TEXT) {
        // to read Wilson lines from file
        readVFromFile(lat, param, (init_method == READ_WLINE_BINARY) ? 2 : 1);
        param->setSuccess(1);
    } else {
        // to generate your own Wilson lines
        if (param->getUseNucleus() == 1) {
            nucleusA_.clear();
            nucleusB_.clear();
            // populate the lists nucleusA_ and nucleusB_ with position data
            sampleTA(param, random, glauber);
            setColorChargeDensity(lat, param, random, glauber);
            // sample color charges and find Wilson lines V_A and V_B
            setV(lat, param, random);
        }
    }
}

void Init::shiftFieldsWithImpactParameter(Lattice *lat, Parameters *param) {
    messager.info("Shifting fields with impact parameter...");
    const double b = param->getb();
    const double phiRP = param->getPhiRP();
    messager << "b = " << b << " fm, phi_RP = " << phiRP;
    messager.flush("info");

    const int N = param->getSize();
    BufferLattice lat_tmp(param->getNc(), param->getSize());
    for (int ipos = 0; ipos < N * N; ipos++) {
        lat_tmp.buffer1[ipos] = lat->U[ipos];
        lat_tmp.buffer2[ipos] = lat->U2[ipos];
    }

    const double L = param->getL();
    const double a = L / N;  // lattice spacing in fm
    for (int ipos = 0; ipos < N * N; ipos++) {
        int ix = ipos / N;
        int iy = ipos % N;
        double x = -L / 2. + a * ix;
        double y = -L / 2. + a * iy;

        double xA = x - b / 2. * cos(phiRP);
        double yA = y - b / 2. * sin(phiRP);
        double xB = x + b / 2. * cos(phiRP);
        double yB = y + b / 2. * sin(phiRP);

        int ixA = static_cast<int>((xA + L / 2.) / a);
        int iyA = static_cast<int>((yA + L / 2.) / a);
        int ixB = static_cast<int>((xB + L / 2.) / a);
        int iyB = static_cast<int>((yB + L / 2.) / a);

        if (ixA < 0 || ixA >= N || iyA < 0 || iyA >= N) {
            lat->U[ipos] = one_;
        } else {
            int posA = ixA * N + iyA;
            lat->U[ipos] = lat_tmp.buffer1[posA];
        }
        if (ixB < 0 || ixB >= N || iyB < 0 || iyB >= N) {
            lat->U2[ipos] = one_;
        } else {
            int posB = ixB * N + iyB;
            lat->U2[ipos] = lat_tmp.buffer2[posB];
        }
    }
}

void Init::initializeForwardLightCone(Lattice *lat, Parameters *param) {
    messager.info("Finding fields in forward lightcone...");
    // output Wilson lines (used also for the proton plots)
    double L = param->getL();
    const int N = param->getSize();
    const int N2 = N * N;
    double a = L / N;  // lattice spacing in fm
    double x, y;
#pragma omp parallel
    {
        Matrix temp2(Nc_, 1);
        Matrix Ux(Nc_, 1);
        Matrix Uy(Nc_, 1);
        Matrix UDx(Nc_, 1);
        Matrix UDy(Nc_, 1);
        Matrix UDx1(Nc_, 1);
        Matrix UDy1(Nc_, 1);

        Matrix Uplaq(Nc_, 1);

        Matrix UDx2(Nc_, 1);
        Matrix UDy2(Nc_, 1);

        Matrix Ux1mUx2(Nc_, 1);
        Matrix UDx1mUDx2(Nc_, 1);
        Matrix Uy1mUy2(Nc_, 1);
        Matrix UDy1mUDy2(Nc_, 1);

// compute Ux(3) Uy(3) after the collision
#pragma omp for
        for (int pos = 0; pos < N2; pos++) {
            // loops over all cells
            auto checkU = lat->U[pos].trace();
            if (checkU != checkU) {
                lat->U[pos] = (one_);
            }

            checkU = lat->U2[pos].trace();
            if (checkU != checkU) {
                lat->U2[pos] = (one_);
            }
        }

#pragma omp for
        for (int pos = 0; pos < N2; pos++) {
            // loops over all cells
            UDx = lat->U[lat->pospX[pos]];
            UDx.conjg();
            lat->Ux1[pos] = (lat->U[pos] * UDx);

            UDy = lat->U[lat->pospY[pos]];
            UDy.conjg();
            lat->Uy1[pos] = (lat->U[pos] * UDy);

            UDx = lat->U2[lat->pospX[pos]];
            UDx.conjg();
            lat->Ux2[pos] = (lat->U2[pos] * UDx);

            UDy = lat->U2[lat->pospY[pos]];
            UDy.conjg();
            lat->Uy2[pos] = (lat->U2[pos] * UDy);
        }
        // -----------------------------------------------------------------
        // from Ux(1,2) and Uy(1,2) compute Ux(3) and Uy(3):

#pragma omp for
        for (int pos = 0; pos < N2; pos++) {
            // loops over all cells
            UDx1 = lat->Ux1[pos];
            UDx2 = lat->Ux2[pos];
            // bool status = findUInForwardLightconeBjoern(UDx1, UDx2, temp2);
            const std::uint64_t retrySeedX = forwardLightconeRetrySeed(
                param->getRandomSeed(), param->getEventId(), pos, 0);
            bool status =
                findUInForwardLightconeChun(UDx1, UDx2, temp2, retrySeedX);
            lat->Ux[pos] = (temp2);
            if (!status) {
                cout << "pos x = " << pos / param->getSize()
                     << " y = " << pos % param->getSize() << endl;
            }

            UDy1 = lat->Uy1[pos];
            UDy2 = lat->Uy2[pos];
            // status = findUInForwardLightconeBjoern(UDy1, UDy2, temp2);
            const std::uint64_t retrySeedY = forwardLightconeRetrySeed(
                param->getRandomSeed(), param->getEventId(), pos, 1);
            status = findUInForwardLightconeChun(UDy1, UDy2, temp2, retrySeedY);
            lat->Uy[pos] = (temp2);
            if (!status) {
                cout << "pos x = " << pos / param->getSize()
                     << " y = " << pos % param->getSize() << endl;
            }
        }

// compute initial electric field
// with minus ax, ay
#pragma omp for
        for (int pos = 0; pos < N2; pos++) {
            // x part in sum:
            Ux1mUx2 = lat->Ux1[pos] - lat->Ux2[pos];
            UDx1 = lat->Ux1[pos];
            UDx1.conjg();
            UDx2 = lat->Ux2[pos];
            UDx2.conjg();
            UDx1mUDx2 = UDx1 - UDx2;

            Ux = lat->Ux[pos];
            UDx = Ux;
            UDx.conjg();

            temp2 = Ux1mUx2 * UDx - Ux1mUx2 - Ux * UDx1mUDx2 + UDx1mUDx2;

            Ux1mUx2 = lat->Ux1[lat->posmX[pos]] - lat->Ux2[lat->posmX[pos]];
            UDx1 = lat->Ux1[lat->posmX[pos]];
            UDx1.conjg();
            UDx2 = lat->Ux2[lat->posmX[pos]];
            UDx2.conjg();
            UDx1mUDx2 = UDx1 - UDx2;

            Ux = lat->Ux[lat->posmX[pos]];
            UDx = Ux;
            UDx.conjg();

            temp2 =
                temp2 - UDx * Ux1mUx2 + Ux1mUx2 + UDx1mUDx2 * Ux - UDx1mUDx2;

            // y part in sum
            Uy1mUy2 = lat->Uy1[pos] - lat->Uy2[pos];
            UDy1 = lat->Uy1[pos];
            UDy1.conjg();
            UDy2 = lat->Uy2[pos];
            UDy2.conjg();
            UDy1mUDy2 = UDy1 - UDy2;

            Uy = lat->Uy[pos];
            UDy = Uy;
            UDy.conjg();

            // y part of the sum:
            temp2 =
                temp2 + Uy1mUy2 * UDy - Uy1mUy2 - Uy * UDy1mUDy2 + UDy1mUDy2;

            Uy1mUy2 = lat->Uy1[lat->posmY[pos]] - lat->Uy2[lat->posmY[pos]];
            UDy1 = lat->Uy1[lat->posmY[pos]];
            UDy1.conjg();
            UDy2 = lat->Uy2[lat->posmY[pos]];
            UDy2.conjg();
            UDy1mUDy2 = UDy1 - UDy2;

            Uy = lat->Uy[lat->posmY[pos]];
            UDy = Uy;
            UDy.conjg();

            temp2 =
                temp2 - UDy * Uy1mUy2 + Uy1mUy2 + UDy1mUDy2 * Uy - UDy1mUDy2;

            lat->U[pos] = ((1. / 8.) * temp2);
        }

        // with plus ax, ay
#pragma omp for
        for (int pos = 0; pos < N * N; pos++) {
            // x part in sum:
            Ux1mUx2 = lat->Ux1[pos] - lat->Ux2[pos];
            UDx1 = lat->Ux1[pos];
            UDx1.conjg();
            UDx2 = lat->Ux2[pos];
            UDx2.conjg();
            UDx1mUDx2 = UDx1 - UDx2;

            Ux = lat->Ux[pos];
            UDx = Ux;
            UDx.conjg();

            temp2 = Ux1mUx2 * UDx - Ux1mUx2 - Ux * UDx1mUDx2 + UDx1mUDx2;

            Ux1mUx2 = lat->Ux1[lat->pospX[pos]] - lat->Ux2[lat->pospX[pos]];
            UDx1 = lat->Ux1[lat->pospX[pos]];
            UDx1.conjg();
            UDx2 = lat->Ux2[lat->pospX[pos]];
            UDx2.conjg();
            UDx1mUDx2 = UDx1 - UDx2;

            Ux = lat->Ux[lat->pospX[pos]];
            UDx = Ux;
            UDx.conjg();

            temp2 =
                temp2 - UDx * Ux1mUx2 + Ux1mUx2 + UDx1mUDx2 * Ux - UDx1mUDx2;

            // y part in sum
            Uy1mUy2 = lat->Uy1[pos] - lat->Uy2[pos];
            UDy1 = lat->Uy1[pos];
            UDy1.conjg();
            UDy2 = lat->Uy2[pos];
            UDy2.conjg();
            UDy1mUDy2 = UDy1 - UDy2;

            Uy = lat->Uy[pos];
            UDy = Uy;
            UDy.conjg();

            // y part of the sum:
            temp2 =
                temp2 + Uy1mUy2 * UDy - Uy1mUy2 - Uy * UDy1mUDy2 + UDy1mUDy2;

            Uy1mUy2 = lat->Uy1[lat->pospY[pos]] - lat->Uy2[lat->pospY[pos]];
            UDy1 = lat->Uy1[lat->pospY[pos]];
            UDy1.conjg();
            UDy2 = lat->Uy2[lat->pospY[pos]];
            UDy2.conjg();
            UDy1mUDy2 = UDy1 - UDy2;

            Uy = lat->Uy[lat->pospY[pos]];
            UDy = Uy;
            UDy.conjg();

            temp2 =
                temp2 - UDy * Uy1mUy2 + Uy1mUy2 + UDy1mUDy2 * Uy - UDy1mUDy2;

            lat->U2[pos] = ((1. / 8.) * temp2);
        }
// compute the plaquette
#pragma omp for
        for (int pos = 0; pos < N * N; pos++) {
            UDx = lat->Ux[lat->pospY[pos]];
            UDy = lat->Uy[pos];

            UDx.conjg();
            UDy.conjg();

            Uplaq = lat->Ux[pos] * (lat->Uy[lat->pospX[pos]] * (UDx * UDy));
            lat->Uy1[pos] = (Uplaq);
        }

#pragma omp for
        for (int pos = 0; pos < N * N; pos++) {
            // AM = (lat->U[pos]); //+lat->cells[pos]->getAetaP());
            // AP = (lat->U2[pos]); //+lat->cells[pos]->getAetaP());

            // this is pi in lattice units as needed for the evolution. (later,
            // the a^4 gives the right units for the energy density
            lat->Ux2[pos] =
                (complex<double>(0., -2. / param->getg()) * (lat->U[pos]));
            // factor -2 because I have A^eta (note the 1/8 before)
            // but want \pi (E^z).

            // lat->Ux2[pos] = (complex<double>(0.,-1./param->getg())*(AM+AP));
            // // factor -2 because I have A^eta (note the 1/8 before) but want
            // \pi (E^z).
        }

        const Matrix zero(Nc_, 0.);
#pragma omp for
        for (int pos = 0; pos < N * N; pos++) {
            lat->U[pos] = (zero);
            lat->U2[pos] = (zero);
            lat->Uy2[pos] = (zero);

            // reset the Ux1 to be used for other purposes later
            lat->Ux1[pos] = (one_);
        }
    }  // omp block
}

void Init::multiplicity(Lattice *lat, Parameters *param) {
    int N = param->getSize();
    int pos;
    double epsilonSum = 0.;
    double L = param->getL();
    double a = L / N;  // lattice spacing in fm

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

    ofstream fout(tE_name.c_str(), std::ios::out);
    fout << epsilonSum << endl;
    fout.close();
}

void Init::generate_nucleus_configuration(
    Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
    double beta3, double beta4, double gamma, bool force_dmin_flag,
    double d_min, double dR_np, double da_np,
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
                    random, A, Z, a_WS, R_WS, beta2, beta3, beta4, gamma, dR_np,
                    da_np, nucleus);
            } else {
                generate_nucleus_configuration_with_deformed_woods_saxon(
                    random, A, Z, a_WS, R_WS, beta2, beta3, beta4, d_min, dR_np,
                    da_np, nucleus);
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
        r_array[i] =
            sample_r_from_woods_saxon(random, a_WS + da_np, R_WS + dR_np);
        idx_array[i] = i;
    }
    std::stable_sort(
        idx_array.begin(), idx_array.end(),
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
                if ((r_i - r_array[j]) * (r_i - r_array[j]) > d_min_sq) break;
                double dsq =
                    ((x_i - x_array[j]) * (x_i - x_array[j])
                     + (y_i - y_array[j]) * (y_i - y_array[j])
                     + (z_i - z_array[j]) * (z_i - z_array[j]));
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
            rv.proton = 0;
        }
        nucleus.push_back(rv);
    }
}

double Init::sample_r_from_woods_saxon(
    Random *random, double a_WS, double R_WS) const {
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
            random, a_WS, R_WS, beta2, beta3, beta4, r_array[i],
            costheta_array[i]);
        idx_array[i] = i;
    }
    for (int i = Z; i < A; i++) {
        sample_r_and_costheta_from_deformed_woods_saxon(
            random, a_WS + da_np, R_WS + dR_np, beta2, beta3, beta4, r_array[i],
            costheta_array[i]);
        idx_array[i] = i;
    }
    std::stable_sort(
        idx_array.begin(), idx_array.end(),
        [&r_array](int i1, int i2) { return r_array[i1] < r_array[i2]; });
    std::stable_sort(r_array.begin(), r_array.end());

    std::vector<double> x_array(A, 0.), y_array(A, 0.), z_array(A, 0.);
    const double d_min_sq = d_min * d_min;
    for (int i = 0; i < A; i++) {
        const double r_i = r_array[i];
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
                if ((r_i - r_array[j]) * (r_i - r_array[j]) > d_min_sq) break;
                double dsq =
                    ((x_i - x_array[j]) * (x_i - x_array[j])
                     + (y_i - y_array[j]) * (y_i - y_array[j])
                     + (z_i - z_array[j]) * (z_i - z_array[j]));
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
    double beta3, double beta4, double gamma, double d_min, double dR_np,
    double da_np, std::vector<ReturnValue> &nucleus) {
    messager << "Sampling nucleon position forcing d_min = " << d_min
             << " fm ...";
    messager.flush("info");
    double rmaxCut = R_WS + dR_np + 10. * (a_WS + da_np);
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
                r = rmaxCut * pow(random->genrand64_real3(), 1.0 / 3.0);
                costheta = 1.0 - 2.0 * random->genrand64_real3();
                phi = 2. * M_PI * random->genrand64_real3();
                double y20 = spherical_harmonics(2, costheta);
                double y30 = spherical_harmonics(3, costheta);
                double y40 = spherical_harmonics(4, costheta);
                double y22 = spherical_harmonics_Y22(costheta, phi);
                R_WS_theta =
                    R_WS_i
                    * (1.0 + beta2 * (cos(gamma) * y20 + sin(gamma) * y22)
                       + beta3 * y30 + beta4 * y40);
            } while (random->genrand64_real3()
                     > fermi_distribution(r, R_WS_theta, a_WS_i));
            double sintheta = sqrt(1. - costheta * costheta);
            x_i = r * sintheta * cos(phi);
            y_i = r * sintheta * sin(phi);
            z_i = r * costheta;
            reSampleFlag = false;
            for (int j = i - 1; j >= 0; j--) {
                double r2 =
                    ((x_i - x_array[j]) * (x_i - x_array[j])
                     + (y_i - y_array[j]) * (y_i - y_array[j])
                     + (z_i - z_array[j]) * (z_i - z_array[j]));
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
    double rmaxCut = R_WS + dR_np + 10. * (a_WS + da_np);
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
            r = rmaxCut * pow(random->genrand64_real3(), 1.0 / 3.0);
            costheta = 1.0 - 2.0 * random->genrand64_real3();
            phi = 2. * M_PI * random->genrand64_real3();
            double y20 = spherical_harmonics(2, costheta);
            double y30 = spherical_harmonics(3, costheta);
            double y40 = spherical_harmonics(4, costheta);
            double y22 = spherical_harmonics_Y22(costheta, phi);
            R_WS_theta = R_WS_i
                         * (1.0 + beta2 * (cos(gamma) * y20 + sin(gamma) * y22)
                            + beta3 * y30 + beta4 * y40);
        } while (random->genrand64_real3()
                 > fermi_distribution(r, R_WS_theta, a_WS_i));
        double sintheta = sqrt(1. - costheta * costheta);
        x_array[i] = r * sintheta * cos(phi);
        y_array[i] = r * sintheta * sin(phi);
        z_array[i] = r * costheta;
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
    } while (random->genrand64_real3()
             > fermi_distribution(r, R_WS_theta, a_WS));
}

double Init::spherical_harmonics(int l, double ct) const {
    // Currently assuming m=0 and available for Y_{20} and Y_{40}
    // "ct" is cos(theta)
    double ylm = 0.0;
    if (l == 2) {
        ylm = 3.0 * ct * ct - 1.0;
        ylm *= 0.31539156525252005;  // pow(5.0/16.0/M_PI,0.5);
    } else if (l == 3) {
        ylm = 5.0 * ct * ct * ct;
        ylm -= 3.0 * ct;
        ylm *= 0.3731763325901154;  // pow(7.0/16.0/M_PI,0.5);
    } else if (l == 4) {
        ylm = 35.0 * ct * ct * ct * ct;
        ylm -= 30.0 * ct * ct;
        ylm += 3.0;
        ylm *= 0.10578554691520431;  // 3.0/16.0/pow(M_PI,0.5);
    }
    return (ylm);
}

double Init::spherical_harmonics_Y22(double ct, double phi) const {
    // Y2,2
    double ylm = 0.0;
    ylm = 1.0 - ct * ct;
    ylm *= cos(2. * phi);
    ylm *= 0.5462742152960397;  // sqrt(2*15)/4/pow(2*M_PI,0.5);
    return (ylm);
}

void Init::recenter_nucleus(
    std::vector<double> &x, std::vector<double> &y, std::vector<double> &z) {
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
    for (auto &n_i : nucleus) {
        meanx += n_i.x;
        meany += n_i.y;
        meanz += n_i.z;
    }

    meanx /= static_cast<double>(nucleus.size());
    meany /= static_cast<double>(nucleus.size());
    meanz /= static_cast<double>(nucleus.size());

    for (auto &n_i : nucleus) {
        n_i.x -= meanx;
        n_i.y -= meany;
        n_i.z -= meanz;
    }
}

void Init::assignProtons(
    Random *random, std::vector<ReturnValue> &nucleus, const int Z) {
    // randomly assign Z nucleons to be protons inside the nucleus
    // Fisher–Yates shuffle using the existing Random instance.
    for (int i = static_cast<int>(nucleus.size()) - 1; i > 0; --i) {
        const int j = static_cast<int>(random->genrand64_real2() * (i + 1));
        std::swap(nucleus[i], nucleus[j]);
    }
    for (unsigned int i = 0; i < nucleus.size(); i++) {
        if (static_cast<int>(i) < std::abs(Z)) {
            nucleus.at(i).proton = 1;
        } else {
            nucleus.at(i).proton = 0;
        }
    }
}

void Init::rotate_nucleus(Random *random, std::vector<ReturnValue> &nucleus) {
    double phi_global = 2. * M_PI * random->genrand64_real3();
    double theta_global = acos(1. - 2. * random->genrand64_real3());
    auto cth = cos(theta_global);
    auto sth = sin(theta_global);
    auto cphi = cos(phi_global);
    auto sphi = sin(phi_global);
    for (auto &n_i : nucleus) {
        auto x_new = cth * cphi * n_i.x - sphi * n_i.y + sth * cphi * n_i.z;
        auto y_new = cth * sphi * n_i.x + cphi * n_i.y + sth * sphi * n_i.z;
        auto z_new = -sth * n_i.x + 0. * n_i.y + cth * n_i.z;
        n_i.x = x_new;
        n_i.y = y_new;
        n_i.z = z_new;
    }
}

void Init::rotate_nucleus(
    double phi_global, double theta_global, std::vector<ReturnValue> &nucleus) {
    auto cth = cos(theta_global);
    auto sth = sin(theta_global);
    auto cphi = cos(phi_global);
    auto sphi = sin(phi_global);
    for (auto &n_i : nucleus) {
        auto x_new = cth * cphi * n_i.x - sphi * n_i.y + sth * cphi * n_i.z;
        auto y_new = cth * sphi * n_i.x + cphi * n_i.y + sth * sphi * n_i.z;
        auto z_new = -sth * n_i.x + 0. * n_i.y + cth * n_i.z;
        n_i.x = x_new;
        n_i.y = y_new;
        n_i.z = z_new;
    }
}

void Init::rotate_nucleus_3D(
    Random *random, std::vector<ReturnValue> &nucleus) {
    // rotate the nucleus with the full three solid angles
    // required for tri-axial deformed nuclei
    // https://en.wikipedia.org/wiki/Euler_angles
    double alpha = 2 * M_PI * random->genrand64_real3();
    double beta = acos(1. - 2. * random->genrand64_real3());
    double gamma = 2 * M_PI * random->genrand64_real3();
    auto c1 = cos(alpha);
    auto s1 = sin(alpha);
    auto c2 = cos(beta);
    auto s2 = sin(beta);
    auto c3 = cos(gamma);
    auto s3 = sin(gamma);
    for (auto &n_i : nucleus) {
        auto x_new = c2 * n_i.x - c3 * s2 * n_i.y + s2 * s3 * n_i.z;
        auto y_new = c1 * s2 * n_i.x + (c1 * c2 * c3 - s1 * s3) * n_i.y
                     + (-c3 * s1 - c1 * c2 * s3) * n_i.z;
        auto z_new = s1 * s2 * n_i.x + (c1 * s3 + c2 * c3 * s1) * n_i.y
                     + (c1 * c3 - c2 * s1 * s3) * n_i.z;
        n_i.x = x_new;
        n_i.y = y_new;
        n_i.z = z_new;
    }
}

double Init::sampleLogNormalDistribution(
    Random *random, const double mean, const double variance) {
    const double meansq = mean * mean;
    const double mu = log(meansq / sqrt(variance + meansq));
    const double sigma = sqrt(log(variance / meansq + 1.));
    double sampleX = exp(mu + sigma * random->Gauss());
    return (sampleX);
}

void Init::sampleQsNormalization(
    Random *random, Parameters *param, const int Nq,
    vector<double> &gauss_array) {
    const double QsSmearWidth = param->getSmearingWidth();
    gauss_array.assign(Nq, 1.);  // default norm = 1
    if (param->getSmearQs() == 1) {
        // introduce a log-normal distribution for Qs normalization
        // dividing by exp(0.5 sigma^2) to ensure the mean is 1
        // the varirance in this case is exp(sigma) - 1 for the log-normal
        // distribution
        for (int iq = 0; iq < Nq; iq++) {
            gauss_array[iq] =
                (exp(random->Gauss(0, QsSmearWidth))
                 / exp(QsSmearWidth * QsSmearWidth / 2.));
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
    return (std::max(1, Nq));
}

bool Init::findUInForwardLightconeBjoern(Matrix &U1, Matrix &U2, Matrix &Usol) {
    const int maxIterations = 100000;

    Matrix U1pU2 = U1 + U2;
    Matrix U1pU2dagger = U1pU2;
    U1pU2dagger.conjg();

    Matrix Mtemp(Nc_, 0.);
    std::vector<Matrix> MtempArr;
    MtempArr.resize(Nc2m1_);
    std::vector<complex<double>> traceCache(Nc2m1_, 0.);
    for (int ai = 0; ai < Nc2m1_; ai++) {
        Mtemp = group_ptr_->getT(ai) * (U1pU2 - U1pU2dagger);
        traceCache[ai] = Mtemp.trace();
        MtempArr[ai] = group_ptr_->getT(ai) * U1pU2;
    }

    // use raw pointers to interface with gsl
    double *Jab = new double[Nc2m1_ * Nc2m1_];
    double *Fa = new double[Nc2m1_];

    double Fzero = 10.;
    double Fprev = 0.;

    // set up initial guess
    std::vector<double> alpha(Nc2m1_, 0.);  // solution
    std::vector<double> alphaSave(Nc2m1_, 0.);
    std::vector<double> Dalpha(Nc2m1_, 0.);
    Usol = getUfromExponent(alpha);
    Matrix Usoldagger = Usol;
    Usoldagger.conjg();

    int iter = 0;
    bool converged = false;
    bool alphaGood = false;
    double lambda = 1.;
    while (!converged && iter < maxIterations) {
        iter++;

        // compute function F that needs to be zero
        Fzero = 0.;
        Mtemp = U1pU2 * Usoldagger - Usol * U1pU2dagger;
        for (int ai = 0; ai < Nc2m1_; ai++) {
            complex<double> traceLoc =
                Mtemp.traceOfProdcutOfMatrix(group_ptr_->getT(ai), Mtemp);
            // minus trace if temp gives -F_ai
            auto traceRes = (-1.) * (traceCache[ai] + traceLoc);
            Fa[ai] = imag(traceRes);
            Fzero += 0.5 * Fa[ai] * Fa[ai];
        }

        // compute Jacobian
        for (int bi = 0; bi < Nc2m1_; bi++) {
            Mtemp = group_ptr_->getT(bi) * Usoldagger;
            for (int ai = 0; ai < Nc2m1_; ai++) {
                int countMe = ai * Nc2m1_ + bi;
                complex<double> traceLoc =
                    Mtemp.traceOfProdcutOfMatrix(MtempArr[ai], Mtemp);
                auto traceRes = -2. * real(traceLoc);
                Jab[countMe] = traceRes;
            }
        }

        solveAxb(Jab, Fa, Dalpha);

        for (int ai = 0; ai < Nc2m1_; ai++) {
            alphaSave[ai] = alpha[ai];
        }

        lambda = 1.;
        Fprev = Fzero;
        alphaGood = false;
        while (!alphaGood) {
            for (int ai = 0; ai < Nc2m1_; ai++) {
                alpha[ai] = alphaSave[ai] + lambda * Dalpha[ai];
            }
            Usol = getUfromExponent(alpha);
            Usoldagger = Usol;
            Usoldagger.conjg();

            Fzero = 0.;
            Mtemp = U1pU2 * Usoldagger - Usol * U1pU2dagger;
            for (int ai = 0; ai < Nc2m1_; ai++) {
                complex<double> traceLoc =
                    Mtemp.traceOfProdcutOfMatrix(group_ptr_->getT(ai), Mtemp);
                // minus trace if temp gives -F_ai
                auto traceRes = (-1.) * (traceCache[ai] + traceLoc);
                Fa[ai] = imag(traceRes);
                Fzero += 0.5 * Fa[ai] * Fa[ai];
            }

            if (lambda < 0.1) {
                for (int ai = 0; ai < Nc2m1_; ai++) {
                    alpha[ai] = 0.1 * random_ptr_->Gauss();
                }
                Usol = getUfromExponent(alpha);
                Usoldagger = Usol;
                Usoldagger.conjg();
                lambda = 1.;
                alphaGood = true;
            }

            if (Fzero > Fprev - 0.00001 * (Fzero * 2.)) {
                lambda *= 0.9;
            } else {
                alphaGood = true;
            }
        }

        if (Fzero < 1e-9) {
            converged = true;
        }
    }
    bool success = true;
    if (iter == maxIterations) {
        std::cout << "Did not converge in findUInForwardLightconeBjoern, "
                  << "Fzero: " << Fzero << std::endl;
        success = false;
        Usol = one_;
    }
    delete[] Fa;
    delete[] Jab;
    return (success);
}

bool Init::findUInForwardLightconeChun(
    Matrix &U1, Matrix &U2, Matrix &Usol, std::uint64_t retrySeed) {
    const int maxIterations = 2000;
    const int maxRetrys = 200;

    Matrix U0 = U1 * U2;
    Matrix U1pU2 = U1 + U2;
    Matrix U1pU2dagger = U1pU2;
    U1pU2dagger.conjg();

    Matrix Mtemp(Nc_, 0.);
    std::vector<Matrix> MtempArr;
    MtempArr.resize(Nc2m1_);
    std::vector<complex<double>> traceCache(Nc2m1_, 0.);
    for (int ai = 0; ai < Nc2m1_; ai++) {
        Mtemp = group_ptr_->getT(ai) * (U1pU2 - U1pU2dagger);
        traceCache[ai] = Mtemp.trace();
        MtempArr[ai] = group_ptr_->getT(ai) * U1pU2;
    }

    // use raw pointers to interface with gsl
    double *Jab = new double[Nc2m1_ * Nc2m1_];
    double *Fa = new double[Nc2m1_];

    double Fzero = 10.;
    double FzeroMin = 1e6;
    Matrix UsolBestEst(Nc_, 1.);

    // set up initial guess
    std::vector<double> alpha(Nc2m1_, 0.);  // solution
    std::vector<double> Dalpha(Nc2m1_, 0.);
    Usol = getUfromExponent(alpha) * U0;
    Matrix Usoldagger = Usol;
    Usoldagger.conjg();

    int iter = 0;
    int nRestart = 0;
    std::uint64_t retryDraw = 0;
    auto nextRetryGaussian = [&]() {
        return deterministicRetryGaussian(retrySeed, retryDraw++);
    };
    while (Fzero > 1e-6 && iter < maxIterations && nRestart < maxRetrys) {
        iter++;

        // compute function F that needs to be zero
        Fzero = 0.;
        Mtemp = U1pU2 * Usoldagger - Usol * U1pU2dagger;
        for (int ai = 0; ai < Nc2m1_; ai++) {
            complex<double> traceLoc =
                Mtemp.traceOfProdcutOfMatrix(group_ptr_->getT(ai), Mtemp);
            // minus trace if temp gives -F_ai
            auto traceRes = (-1.) * (traceCache[ai] + traceLoc);
            Fa[ai] = imag(traceRes);
            Fzero += std::abs(Fa[ai]);
        }

        // compute Jacobian
        // numerical formula
        bool JabGood = true;
        for (int bi = 0; bi < Nc2m1_; bi++) {
            double dalpha_bi =
                (std::max(0.001, std::min(10., 0.01 * std::abs(alpha[bi]))));
            alpha[bi] = alpha[bi] + dalpha_bi;
            Mtemp = getUfromExponent(alpha) * U0;
            Mtemp.conjg();
            Mtemp = U1pU2 * (Mtemp - Usoldagger);
            double Mcheck = 0.;
            for (int ai = 0; ai < Nc2m1_; ai++) {
                int countMe = ai * Nc2m1_ + bi;
                complex<double> traceLoc =
                    Mtemp.traceOfProdcutOfMatrix(group_ptr_->getT(ai), Mtemp);
                Jab[countMe] = 2. * imag(traceLoc) / dalpha_bi;
                Mcheck += std::abs(Jab[countMe]);
                if (Mcheck < 1e-15) {
                    // avoid matrix to be singular
                    JabGood = false;
                }
            }
            alpha[bi] = alpha[bi] - dalpha_bi;
        }
        if (!JabGood) {
            // analytical approximated formula
            for (int bi = 0; bi < Nc2m1_; bi++) {
                Mtemp = group_ptr_->getT(bi) * Usoldagger;
                for (int ai = 0; ai < Nc2m1_; ai++) {
                    int countMe = ai * Nc2m1_ + bi;
                    complex<double> traceLoc =
                        Mtemp.traceOfProdcutOfMatrix(MtempArr[ai], Mtemp);
                    auto traceRes = -2. * real(traceLoc);
                    Jab[countMe] = traceRes;
                }
            }
        }

        solveAxb(Jab, Fa, Dalpha);

        bool DalphaCheck = true;
        for (int ai = 0; ai < Nc2m1_; ai++) {
            if (std::abs(Dalpha[ai]) > 1e5) {
                DalphaCheck = false;
                break;
            }
        }
        if (DalphaCheck) {
            for (int ai = 0; ai < Nc2m1_; ai++) {
                alpha[ai] = alpha[ai] + Dalpha[ai];
            }
        } else {
            for (int ai = 0; ai < Nc2m1_; ai++) {
                alpha[ai] = nextRetryGaussian();
            }
        }
        Usol = getUfromExponent(alpha) * U0;
        Usoldagger = Usol;
        Usoldagger.conjg();

        if (Fzero < FzeroMin) {
            FzeroMin = Fzero;
            UsolBestEst = Usol;
        }
        if (iter == maxIterations) {
            for (int ai = 0; ai < Nc2m1_; ai++) {
                alpha[ai] = nextRetryGaussian();
            }
            Usol = getUfromExponent(alpha) * U0;
            Usoldagger = Usol;
            Usoldagger.conjg();

            nRestart++;
            iter = 0;
        }
    }
    bool success = true;
    if (nRestart == maxRetrys) {
        std::cout << "Did not converge in findUInForwardLightconeChun, "
                  << "Fzero: " << FzeroMin << std::endl;
        Usol = UsolBestEst;  // return the best estimate
        success = false;
    }
    delete[] Fa;
    delete[] Jab;
    return (success);
}

Matrix Init::getUfromExponent(std::vector<double> &in) {
    Matrix tempM(Nc_, Matrix::noInit);
    complex<double> U[9];

    // expmCoeff calculates the coefficients of exp(i in[a] t[a]).  Build
    // the 3x3 matrix directly from the fixed SU(3) generators instead of
    // allocating a coefficient vector and materializing eight scaled Matrix
    // temporaries plus the chained sums.
    tempM.expmCoeff(in.data(), U);
    if (std::abs(U[0].real()) < 1e-15) {
        tempM = one_;
    } else {
        const complex<double> I(0., 1.);
        const double invSqrt3 = 1. / std::sqrt(3.);

        tempM.set(0, 0, U[0] + 0.5 * U[3] + 0.5 * invSqrt3 * U[8]);
        tempM.set(1, 1, U[0] - 0.5 * U[3] + 0.5 * invSqrt3 * U[8]);
        tempM.set(2, 2, U[0] - invSqrt3 * U[8]);

        tempM.set(0, 1, 0.5 * (U[1] - I * U[2]));
        tempM.set(1, 0, 0.5 * (U[1] + I * U[2]));
        tempM.set(0, 2, 0.5 * (U[4] - I * U[5]));
        tempM.set(2, 0, 0.5 * (U[4] + I * U[5]));
        tempM.set(1, 2, 0.5 * (U[6] - I * U[7]));
        tempM.set(2, 1, 0.5 * (U[6] + I * U[7]));
    }
    return tempM;
}

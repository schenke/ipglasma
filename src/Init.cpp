// Init.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#include "Init.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Instrumentation.h"
#include "PhysConst.h"
#include "RunningCoupling.h"
#include "WilsonLineIO.h"
#include "gsl/gsl_linalg.h"

using PhysConst::hbarc;
using PhysConst::mbToFm2;
using PhysConst::Nc;
using PhysConst::Nc2m1;
using PhysConst::smallEps;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;

//**************************************************************************
// Init class.

// This function samples the nucleon positions inside the projectile and
// target nuclei. Both nuclei are centered at the origin.
void Init::sampleTA(Parameters *param, Random *random, Glauber *glauber) {
    IPG_PROFILE_SCOPE("initialization.sample_nuclei");
    messager_.info("[Init::sampleTA]: Sampling nucleon positions ... ");

    if (param->nucleus.nucleonPositionsFromFile) {
        sampleTAFromConfigFiles(random, glauber);
    } else {
        sampleTAWoodsSaxon(param, random, glauber);
    }

    // global rotation of the nucleus
    applyPolarizationRotation(
        random, param->nucleus.polarizationProjectile, nucleusA_);
    applyPolarizationRotation(
        random, param->nucleus.polarizationTarget, nucleusB_);
}

void Init::sampleTAWoodsSaxon(
    Parameters *param, Random *random, Glauber *glauber) {
    ReturnValue rv, rv2;
    if (param->collision.nucleiToAverage > 1) {
        if ((glauber->nucleusA1() == 1 || glauber->nucleusA2() == 1)) {
            messager_ << "[Init::sampleTA]: Averaging over nuclei is not "
                         "supported for collisions involving protons. "
                         "Exiting.";
            messager_.flush("error");
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
        rv = glauber->sampleTARejection(random, NucleusRole::Projectile);
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
        generateNucleusConfiguration(
            random, A1, Z1, glauber->getGlauberData().projectile.a_WS,
            glauber->getGlauberData().projectile.R_WS,
            glauber->getGlauberData().projectile.beta2,
            glauber->getGlauberData().projectile.beta3,
            glauber->getGlauberData().projectile.beta4,
            glauber->getGlauberData().projectile.gamma,
            glauber->getGlauberData().projectile.forceDminFlag,
            glauber->getGlauberData().projectile.d_min,
            glauber->getGlauberData().projectile.dR_np,
            glauber->getGlauberData().projectile.da_np, nucleusA_);
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
        rv = glauber->sampleTARejection(random, NucleusRole::Target);
        // we sample the neutron proton distance, so distance to the center
        // needs to be divided by 2

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
        generateNucleusConfiguration(
            random, A2, Z2, glauber->getGlauberData().target.a_WS,
            glauber->getGlauberData().target.R_WS,
            glauber->getGlauberData().target.beta2,
            glauber->getGlauberData().target.beta3,
            glauber->getGlauberData().target.beta4,
            glauber->getGlauberData().target.gamma,
            glauber->getGlauberData().target.forceDminFlag,
            glauber->getGlauberData().target.d_min,
            glauber->getGlauberData().target.dR_np,
            glauber->getGlauberData().target.da_np, nucleusB_);
    }
}

void Init::sampleTAFromConfigFiles(Random *random, Glauber *glauber) {
    ReturnValue rv;
    if (nucleonPosArrA_.size() > 0) {
        double ran2 = random->genrand64_real3();
        int nucleusNumber = static_cast<int>(ran2 * nucleonPosArrA_.size());
        messager_ << "[Init::sampleTA]: using nucleus Number = "
                  << nucleusNumber;
        messager_.flush("info");
        for (int iA = 0; iA < glauber->nucleusA1(); iA++) {
            rv.x = nucleonPosArrA_[nucleusNumber][3 * iA];
            rv.y = nucleonPosArrA_[nucleusNumber][3 * iA + 1];
            rv.z = nucleonPosArrA_[nucleusNumber][3 * iA + 2];
            rv.collided = 0;
            nucleusA_.push_back(rv);
        }
        assignProtons(random, nucleusA_, glauber->nucleusZ1());
        recenterNucleus(nucleusA_);
    } else {
        // no configurations, sample with Woods-Saxon
        messager_ << "[Init::sampleTA]: configuration file for A = "
                  << glauber->nucleusA1()
                  << " is not available, generate the nucleus configuration "
                  << "using Woods-Saxon distribution instead.";
        messager_.flush("info");

        generateNucleusConfiguration(
            random, glauber->nucleusA1(), glauber->nucleusZ1(),
            glauber->getGlauberData().projectile.a_WS,
            glauber->getGlauberData().projectile.R_WS,
            glauber->getGlauberData().projectile.beta2,
            glauber->getGlauberData().projectile.beta3,
            glauber->getGlauberData().projectile.beta4,
            glauber->getGlauberData().projectile.gamma,
            glauber->getGlauberData().projectile.forceDminFlag,
            glauber->getGlauberData().projectile.d_min,
            glauber->getGlauberData().projectile.dR_np,
            glauber->getGlauberData().projectile.da_np, nucleusA_);
    }

    if (nucleonPosArrB_.size() > 0) {
        double ran2 = random->genrand64_real3();
        int nucleusNumber = static_cast<int>(ran2 * nucleonPosArrB_.size());
        messager_ << "[Init::sampleTA]: using nucleus Number = "
                  << nucleusNumber;
        messager_.flush("info");
        for (int iA = 0; iA < glauber->nucleusA2(); iA++) {
            rv.x = nucleonPosArrB_[nucleusNumber][3 * iA];
            rv.y = nucleonPosArrB_[nucleusNumber][3 * iA + 1];
            rv.z = nucleonPosArrB_[nucleusNumber][3 * iA + 2];
            rv.collided = 0;
            nucleusB_.push_back(rv);
        }
        assignProtons(random, nucleusB_, glauber->nucleusZ2());
        recenterNucleus(nucleusB_);
    } else {
        // no configurations, sample with Woods-Saxon
        messager_ << "[Init::sampleTA]: configuration file for A = "
                  << glauber->nucleusA2()
                  << " is not available, generate the nucleus configuration "
                  << "using Woods-Saxon distribution instead.";
        messager_.flush("info");
        generateNucleusConfiguration(
            random, glauber->nucleusA2(), glauber->nucleusZ2(),
            glauber->getGlauberData().target.a_WS,
            glauber->getGlauberData().target.R_WS,
            glauber->getGlauberData().target.beta2,
            glauber->getGlauberData().target.beta3,
            glauber->getGlauberData().target.beta4,
            glauber->getGlauberData().target.gamma,
            glauber->getGlauberData().target.forceDminFlag,
            glauber->getGlauberData().target.d_min,
            glauber->getGlauberData().target.dR_np,
            glauber->getGlauberData().target.da_np, nucleusB_);
    }
}

void Init::applyPolarizationRotation(
    Random *random, int polarizationFlag, std::vector<ReturnValue> &nucleus) {
    if (polarizationFlag == 0) {
        rotateNucleus3D(random, nucleus);
    } else if (polarizationFlag == 1) {
        // longitudinal polarization only rotates phi randomly
        double phi = 2. * M_PI * random->genrand64_real3();
        double theta = 0;
        rotateNucleus(phi, theta, nucleus);
    } else if (polarizationFlag == 2) {
        // transverse polarization rotates J to +y axis
        double phi = M_PI / 2;
        double theta = M_PI / 2;
        rotateNucleus(phi, theta, nucleus);
    }
}

void Init::readNuclearQs(Parameters *param) {
    IPG_PROFILE_SCOPE("initialization.read_qs_table");
    // steps in qs0 and Y in the file
    string dummy;
    string T, Qs;
    // open file

    messager_ << "[Init::readNuclearQs]: Reading Q_s(sum(T_p),y) from file ";
    messager_ << param->colorCharge.nucleusQsTableFileName << " ... ";
    messager_.flush("info");

    ifstream fin;
    fin.open((param->colorCharge.nucleusQsTableFileName).c_str());
    if (fin) {
        for (int iT = 0; iT < iTpmax_; iT++) {
            for (int iy = 0; iy < iymaxNuc_; iy++) {
                if (!fin.eof()) {
                    fin >> dummy;
                    fin >> T;
                    Tlist_[iT] = atof(T.c_str());
                    fin >> Qs;
                    Qs2Nuclear_[iT][iy] = atof(Qs.c_str());
                } else {
                    messager_ << "[Init::readNuclearQs]: End of file reached "
                                 "prematurely -- did the Q_s table file "
                                 "change? Exiting.";
                    messager_.flush("error");
                    exit(1);
                }
            }
        }
        fin.close();
    } else {
        messager_ << "[Init::readNuclearQs]: File "
                  << param->colorCharge.nucleusQsTableFileName
                  << " does not exist. Exiting.";
        messager_.flush("error");
        exit(1);
    }
}

void Init::readInNucleusConfigs(
    const int nucleusA, const int lightNucleusOption,
    const int polarizationFlag, const double polJz,
    vector<vector<float>> &nucleonPosArr, Parameters *param) {
    if (nucleonPosArr.size() > 0) return;
    std::string path = param->nucleus.nuclearConfigurationsPath + "/";
    std::string fileName;
    bool readFlag = true;
    if (nucleusA == 2) {
        fileName = "DeuteronPol0Configs.bin.in";
        if (std::abs(std::abs(polJz) - 1.) < 1e-8)
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
    messager_
        << "[Init::readInNucleusConfigs]: read in nucleus configurations from "
        << fileName;
    messager_.flush("info");
    std::ifstream inFile(fileName, std::ios::binary);
    if (!inFile) {
        messager_ << "[Init::readInNucleusConfigs]: File " << fileName
                  << " not found. Exiting.";
        messager_.flush("error");
        exit(1);
    }
    while (true) {
        vector<float> tempPos;
        for (int i = 0; i < nucleusA; i++) {
            for (int j = 0; j < Nentry; j++) {
                float temp;
                inFile.read(reinterpret_cast<char *>(&temp), sizeof(float));
                // Au197/Pb208 files carry a 4th per-nucleon entry; only
                // (x, y, z) are kept so callers can index with 3 * iA.
                // Protons are assigned afterwards by assignProtons().
                if (j < 3) tempPos.push_back(temp);
            }
        }
        if (inFile.eof()) break;
        nucleonPosArr.push_back(tempPos);
    }
    inFile.close();
    messager_ << "[Init::readInNucleusConfigs]: read in "
              << nucleonPosArr.size() << " configurations.";
    messager_.flush("info");
}

void Init::samplePartonPositions(
    Parameters *param, Random *random, vector<double> &x_array,
    vector<double> &y_array, vector<double> &z_array,
    vector<double> &BGq_array) {
    const double sqrtBG = sqrt(param->subnucleon.BG) * hbarc;  // fm
    const double BGqMean = param->subnucleon.BGq;
    const double BGqVar = param->subnucleon.BGqVar;
    const double BGq =
        (0.09 + sampleLogNormalDistribution(random, BGqMean - 0.09, BGqVar));
    const double QsSmearWidth = param->subnucleon.smearingWidth;
    const int Nq = sampleNumberOfPartons(random, param);
    const double dq_min = param->subnucleon.dqMin;  // fm
    const double dq_min_sq = dq_min * dq_min;
    const double omega = param->subnucleon.omega;

    vector<double> r_array(Nq, 0.);
    BGq_array.assign(Nq, BGq);
    for (int iq = 0; iq < Nq; iq++) {
        if (std::abs(omega - 1) < 1e-8) {
            double xq = sqrtBG * random->gauss();
            double yq = sqrtBG * random->gauss();
            double zq = sqrtBG * random->gauss();
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
    if (param->subnucleon.shiftConstituentQuarkProtonOrigin) {
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
    posy = static_cast<int>(floor(y / deltaYNuc_ + 0.0000001));

    if (y > iymaxNuc_ * deltaYNuc_) {
        // getNuclearQs2() is called from inside an omp parallel for loop
        // (setColorChargeDensity), so use a fresh, stack-local instance
        // rather than sharing messager_, which is not thread-safe.
        PrettyOstream localMessager;
        localMessager << "[Init::getNuclearQs2]: y=" << y
                      << " is above the tabulated range (max y="
                      << iymaxNuc_ * deltaYNuc_ << "). Exiting.";
        localMessager.flush("error");
        exit(1);
    }

    //  if ( T > Qs2Nuclear_[iTpmax_-1][iymaxNuc_-1] )
    if (T > Tlist_[iTpmax_ - 1]) {
        // Local instance for the same omp thread-safety reason as above.
        PrettyOstream localMessager;
        localMessager << "[Init::getNuclearQs2]: T=" << T
                      << " exceeds the tabulated range (max T="
                      << Tlist_[iTpmax_ - 1]
                      << "); clamping to the maximal tabulated T.";
        localMessager.flush("warning");
        check = 1;
        fracy = (y - static_cast<double>(posy) * deltaYNuc_) / deltaYNuc_;
        QsYdown = (Qs2Nuclear_[iTpmax_ - 1][posy]);
        QsYup = (Qs2Nuclear_[iTpmax_ - 1][posy + 1]);
        value = (fracy * QsYup + (1. - fracy) * QsYdown);  //*hbarc*hbarc;
        return value;
    }

    if (T < Tlist_[0]) {
        check = 1;
        return 0.;
    }

    for (int iT = 0; iT < iTpmax_; iT++) {
        if (T >= Tlist_[iT] && T < Tlist_[iT + 1]) {
            fracT = (T - Tlist_[iT]) / (Tlist_[iT + 1] - Tlist_[iT]);
            fracy = (y - static_cast<double>(posy) * deltaYNuc_) / deltaYNuc_;

            QsYdown = (fracT) * (Qs2Nuclear_[iT + 1][posy])
                      + (1. - fracT) * (Qs2Nuclear_[iT][posy]);
            QsYup = (fracT) * (Qs2Nuclear_[iT + 1][posy + 1])
                    + (1. - fracT) * (Qs2Nuclear_[iT][posy + 1]);
            value = (fracy * QsYup + (1. - fracy) * QsYdown);  //*hbarc*hbarc;

            check++;
            continue;
        }
    }

    if (check != 1) {
        // Local instance: same omp thread-safety reason as above.
        PrettyOstream localMessager;
        localMessager << "[Init::getNuclearQs2]: could not uniquely "
                         "determine Qs^2 (check="
                      << check << ", T=" << T
                      << "); falling back to the maximal tabulated T_p.";
        localMessager.flush("warning");
        value =
            (fracy * Qs2Nuclear_[iTpmax_ - 1][posy + 1]
             + (1. - fracy) * Qs2Nuclear_[iTpmax_ - 1][posy]);
    }

    return value;
}

double Init::computeFluctuatingXG2mu2(
    Parameters *param, double a, double rapidity, double Tp, double qsmuRatio,
    double ySign) {
    const double exponent = 5.6;  // see 1212.2974 Eq. (17)
    double Qs = 1.;
    double xVal = 0.;
    double g2mu2 = 0.;
    double localrapidity = rapidity;
    double yIn = rapidity;
    double Ydeviation = 10000;
    // iterative loops here to determine the fluctuating Y
    while (std::abs(Ydeviation) > 0.001) {
        if (localrapidity >= 0) {
            Qs = sqrt(getNuclearQs2(Tp, localrapidity));
        } else {
            xVal = Qs * param->colorCharge.xQsFactor / param->collision.sqrtS
                   * exp(ySign * yIn);
            if (xVal == 0)
                Qs = 0.;
            else
                Qs = sqrt(getNuclearQs2(Tp, 0.))
                     * sqrt(
                         pow((1 - xVal) / (1 - 0.01), exponent)
                         * pow((0.01 / xVal), 0.2));
        }
        if (Qs == 0) {
            Ydeviation = 0;
            g2mu2 = 0.;
        } else {
            g2mu2 = Qs * Qs / qsmuRatio / qsmuRatio * a * a / hbarc / hbarc
                    / param->coupling.g
                    / param->coupling.g;  // lattice units? check

            Ydeviation = localrapidity
                         - log(
                             0.01
                             / (Qs * param->colorCharge.xQsFactor
                                / param->collision.sqrtS * exp(ySign * yIn)));
            localrapidity =
                log(0.01
                    / (Qs * param->colorCharge.xQsFactor
                       / param->collision.sqrtS * exp(ySign * yIn)));
        }
    }
    if (g2mu2 != g2mu2) {
        g2mu2 = 0.;
    }
    return g2mu2;
}

// set g^2\mu^2 as the sum of the individual nucleons' g^2\mu^2, using
// Q_s(b,y) prop to g^mu(b,y) Also compute N_part using Glauber
void Init::computeCellColorCharge(
    Lattice *lat, Parameters *param, int ipos, double a, double rapidityA,
    double rapidityB) {
    if (param->colorCharge.useFluctuatingX) {  // Local Qs dependent x
        lat->cells[ipos]->setg2mu2A(computeFluctuatingXG2mu2(
            param, a, rapidityA, lat->cells[ipos]->getTpA(),
            param->colorCharge.QsMuRatio, 1.));
        lat->cells[ipos]->setg2mu2B(computeFluctuatingXG2mu2(
            param, a, rapidityB, lat->cells[ipos]->getTpB(),
            param->colorCharge.QsMuRatio, -1.));
    } else {  // Fixed x
        // nucleus A
        lat->cells[ipos]->setg2mu2A(
            getNuclearQs2(lat->cells[ipos]->getTpA(), rapidityA)
            / param->colorCharge.QsMuRatio / param->colorCharge.QsMuRatio * a
            * a / hbarc / hbarc / param->coupling.g
            / param->coupling.g);  // lattice units? check

        // nucleus B
        lat->cells[ipos]->setg2mu2B(
            getNuclearQs2(lat->cells[ipos]->getTpB(), rapidityB)
            / param->colorCharge.QsMuRatio / param->colorCharge.QsMuRatio * a
            * a / hbarc / hbarc / param->coupling.g / param->coupling.g);
    }
}

void Init::setColorChargeDensity(
    Lattice *lat, Parameters *param, Random *random, Glauber *glauber) {
    IPG_PROFILE_SCOPE("initialization.color_charge_density");
    messager_.info(
        "[Init::setColorChargeDensity]: set color charge density ...");

    const int N = param->lattice.size;
    const double a = param->lattice.L / N;  // lattice spacing in fm

    double rapidityA = 0.;
    double rapidityB = 0.;
    computeEffectiveRapidities(param, rapidityA, rapidityB);

    double nucleiInAverage =
        static_cast<double>(param->collision.nucleiToAverage);

    if (!param->collision.useNucleus) {
        setConstantColorChargeDensity(lat, param);
        return;
    }

#pragma omp parallel for
    for (int ipos = 0; ipos < N * N; ipos++) {
        lat->cells[ipos]->setg2mu2A(0.);
        lat->cells[ipos]->setg2mu2B(0.);
    }

    sampleNucleonAnisotropyAngles(param, random);
    sampleConstituentQuarkGeometry(param, random);

    // test what a smooth Woods-Saxon would give
    if (param->nucleus.useSmoothNucleus) {
        computeSmoothNucleusThickness(lat, param, glauber);
    } else {
        // Non-smooth nucleus add all T_p's (new in version 1.2)
        computeThicknessFromNucleons(lat, param, nucleiInAverage);
    }

// get Q_s^2 (and from that g^2mu^2) for a given \sum T_p and Y
#pragma omp parallel for
    for (int ipos = 0; ipos < N * N; ipos++) {
        computeCellColorCharge(lat, param, ipos, a, rapidityA, rapidityB);
    }
    messager_.info(
        "[Init::setColorChargeDensity]: Color charge densities for nucleus A "
        "and B set. ");
}

void Init::computeEffectiveRapidities(
    Parameters *param, double &rapidityA, double &rapidityB) {
    if (!param->colorCharge.usePseudoRapidity) {
        rapidityA = param->colorCharge.rapidityA;
        rapidityB = param->colorCharge.rapidityB;
        return;
    }
    // when using pseudorapidity as input convert to rapidity here.
    // later include Jacobian in multiplicity and energy
    messager_ << "[Init::setColorChargeDensity]: Using pseudorapidity "
              << param->colorCharge.rapidityA << ", "
              << param->colorCharge.rapidityB;
    messager_.flush("info");
    double m = param->colorCharge.jacobianMass;  // in GeV
    double P =
        0.13 + 0.32 * pow(param->collision.sqrtS / 1000., 0.115);  // in GeV
    rapidityA =
        0.5
        * log(
            sqrt(pow(cosh(param->colorCharge.rapidityA), 2.) + m * m / (P * P))
            + sinh(param->colorCharge.rapidityA)
                  / (sqrt(
                         pow(cosh(param->colorCharge.rapidityA), 2.)
                         + m * m / (P * P))
                     - sinh(param->colorCharge.rapidityA)));
    rapidityB =
        0.5
        * log(
            sqrt(pow(cosh(param->colorCharge.rapidityB), 2.) + m * m / (P * P))
            + sinh(param->colorCharge.rapidityB)
                  / (sqrt(
                         pow(cosh(param->colorCharge.rapidityB), 2.)
                         + m * m / (P * P))
                     - sinh(param->colorCharge.rapidityB)));
    messager_ << "[Init::setColorChargeDensity]: Corresponds to rapidity "
              << rapidityA << ", " << rapidityB;
    messager_.flush("info");
}

void Init::setConstantColorChargeDensity(Lattice *lat, Parameters *param) {
    const int N = param->lattice.size;
    const double L = param->lattice.L;
    if (param->collision.useGaussian) {
        double sigmax = 0.35;
        double sigmay = 0.5;
        for (int ix = 0; ix < N; ix++)  // loop over all positions
        {
            double x = ix * L / double(N) - L / 2.;
            for (int iy = 0; iy < N; iy++) {
                double y = iy * L / double(N) - L / 2.;
                int localpos = lat->positionFromXY(ix, iy);
                double envelope = exp(
                                      -(x * x / (2. * sigmax * sigmax)
                                        + y * y / (2. * sigmay * sigmay)))
                                  / (2. * M_PI * sigmax * sigmay);
                lat->cells[localpos]->setg2mu2A(
                    envelope * param->collision.g2mu * param->collision.g2mu
                    / param->coupling.g / param->coupling.g);
                lat->cells[localpos]->setg2mu2B(
                    envelope * param->collision.g2mu * param->collision.g2mu
                    / param->coupling.g / param->coupling.g);
            }
        }
    } else {
        for (int ix = 0; ix < N; ix++)  // loop over all positions
        {
            for (int iy = 0; iy < N; iy++) {
                int localpos = lat->positionFromXY(ix, iy);
                lat->cells[localpos]->setg2mu2A(
                    param->collision.g2mu * param->collision.g2mu
                    / param->coupling.g / param->coupling.g);
                lat->cells[localpos]->setg2mu2B(
                    param->collision.g2mu * param->collision.g2mu
                    / param->coupling.g / param->coupling.g);
            }
        }
    }
    param->event.success = 1;
    messager_.info(
        "[Init::setColorChargeDensity]: constant color charge density set");
}

void Init::sampleNucleonAnisotropyAngles(Parameters *param, Random *random) {
    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();
    double xi = param->subnucleon.protonAnisotropy;
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
}

void Init::sampleConstituentQuarkGeometry(Parameters *param, Random *random) {
    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();
    const double NqFlag = param->subnucleon.Nq;
    vector<double> x_array, y_array, z_array, BGq_array, gauss_array;
    xq1_.clear();
    xq2_.clear();
    yq1_.clear();
    yq2_.clear();
    BGq1_.clear();
    BGq2_.clear();
    gauss1_.clear();
    gauss2_.clear();
    for (int i = 0; i < A1; i++) {
        int Npartons = 1;
        if (NqFlag > 0) {
            samplePartonPositions(
                param, random, x_array, y_array, z_array, BGq_array);
            // if (param->subnucleon.shiftConstituentQuarkProtonOrigin)
            // Move center of mass to the origin
            // Note that 1607.01711 this is not done, so parameters quoted
            // in that paper can't be used if this is done
            xq1_.push_back(x_array);
            yq1_.push_back(y_array);
            BGq1_.push_back(BGq_array);
            Npartons = std::max(1, static_cast<int>(x_array.size()));
        }
        sampleQsNormalization(random, param, Npartons, gauss_array);
        gauss1_.push_back(gauss_array);
    }

    for (int i = 0; i < A2; i++) {
        int Npartons = 1;
        if (NqFlag > 0) {
            samplePartonPositions(
                param, random, x_array, y_array, z_array, BGq_array);
            xq2_.push_back(x_array);
            yq2_.push_back(y_array);
            BGq2_.push_back(BGq_array);
            Npartons = std::max(1, static_cast<int>(x_array.size()));
        }
        sampleQsNormalization(random, param, Npartons, gauss_array);
        gauss2_.push_back(gauss_array);
    }
}

void Init::computeSmoothNucleusThickness(
    Lattice *lat, Parameters *param, Glauber *glauber) {
    messager_ << "[Init::setColorChargeDensity]: Using smooth nucleus for "
                 "test purposes. "
                 "Does not "
                 "include "
                 "deformation.";
    messager_.flush("info");
    const int N = param->lattice.size;
    const double L = param->lattice.L;
    const double a = L / N;
    // Both profiles are centered at the origin: the impact parameter is
    // sampled later and applied by shiftFieldsWithImpactParameter(), the
    // same as for nucleon-based nuclei. Nucleus A is the projectile.
    double normA = 0.;
    double normB = 0.;
    for (int ix = 0; ix < N; ix++) {
        const double x = -L / 2. + a * ix;
        for (int iy = 0; iy < N; iy++) {
            const double y = -L / 2. + a * iy;
            const int localpos = lat->positionFromXY(ix, iy);
            const double r = sqrt(x * x + y * y);

            double TA = glauber->interNuPInSP(r);
            double TB = glauber->interNuTInST(r);
            normA += TA * a * a;
            normB += TB * a * a;

            // remove the far tails of each nucleus separately (the overlap
            // region is not known before b is sampled)
            if (TA < 0.001) TA = 0.;
            if (TB < 0.001) TB = 0.;
            lat->cells[localpos]->setTpA(TA);
            lat->cells[localpos]->setTpB(TB);
        }
    }

    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            const int localpos = lat->positionFromXY(ix, iy);
            lat->cells[localpos]->setTpA(
                lat->cells[localpos]->getTpA() / normA * glauber->nucleusA1()
                * hbarc * hbarc);
            lat->cells[localpos]->setTpB(
                lat->cells[localpos]->getTpB() / normB * glauber->nucleusA2()
                * hbarc * hbarc);
        }
    }
}

double Init::computeNucleonThicknessAtCell(
    Parameters *param, const std::vector<ReturnValue> &nucleus,
    const vector<vector<double>> &xq, const vector<vector<double>> &yq,
    const vector<vector<double>> &BGq, const vector<vector<double>> &gauss,
    double x, double y, double xi, double nucleiInAverage) {
    const int A = nucleus.size();
    double Tp = 0.;
    for (int i = 0; i < A; i++) {
        double xm = nucleus.at(i).x;
        double ym = nucleus.at(i).y;

        double T = 0.;
        double bp2 = 0.;
        if (param->subnucleon.Nq > 0) {
            for (unsigned int iq = 0; iq < xq[i].size(); iq++) {
                bp2 = (xm + xq[i][iq] - x) * (xm + xq[i][iq] - x)
                      + (ym + yq[i][iq] - y) * (ym + yq[i][iq] - y);
                bp2 /= hbarc * hbarc;

                T += exp(-bp2 / (2. * BGq[i][iq])) / (2. * M_PI * BGq[i][iq])
                     / (static_cast<double>(xq[i].size())) * gauss[i][iq];
            }
        } else {
            const double BG = param->subnucleon.BG;
            double phi = nucleus.at(i).phi;

            bp2 = (xm - x) * (xm - x) + (ym - y) * (ym - y)
                  + xi * pow((xm - x) * cos(phi) + (ym - y) * sin(phi), 2.);
            bp2 /= hbarc * hbarc;
            T = sqrt(1 + xi) * exp(-bp2 / (2. * BG)) / (2. * M_PI * BG)
                * gauss[i][0];  // T_p in this cell for the
                                // current nucleon
        }
        Tp += T / nucleiInAverage;  // add up all T_p
    }
    return Tp;
}

void Init::computeThicknessFromNucleons(
    Lattice *lat, Parameters *param, double nucleiInAverage) {
    const int N = param->lattice.size;
    const double L = param->lattice.L;
    const double a = L / N;
    const double xi = param->subnucleon.protonAnisotropy;

#pragma omp parallel for
    for (int ipos = 0; ipos < N * N; ipos++) {
        // loop over all positions
        int iy = lat->yFromPosition(ipos);
        int ix = lat->xFromPosition(ipos);
        double x = -L / 2. + a * ix;
        double y = -L / 2. + a * iy;

        lat->cells[ipos]->setTpA(computeNucleonThicknessAtCell(
            param, nucleusA_, xq1_, yq1_, BGq1_, gauss1_, x, y, xi,
            nucleiInAverage));
        lat->cells[ipos]->setTpB(computeNucleonThicknessAtCell(
            param, nucleusB_, xq2_, yq2_, BGq2_, gauss2_, x, y, xi,
            nucleiInAverage));
    }
}

void Init::computeNcollList(
    Parameters *param, double d2, double b, double phiRP, int &Ncoll) {
    stringstream strNcoll_name;
    strNcoll_name << "NcollList" << param->event.eventId << ".dat";
    string Ncoll_name;
    Ncoll_name = strNcoll_name.str();

    ofstream foutNcoll(Ncoll_name.c_str(), std::ios::out);

    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();
    const bool gaussianWounding = param->collision.gaussianWounding;
    const double G = 0.92;
    for (int i = 0; i < A1; i++) {
        for (int j = 0; j < A2; j++) {
            double dx =
                (nucleusB_.at(j).x - nucleusA_.at(i).x - b * cos(phiRP));
            double dy =
                (nucleusB_.at(j).y - nucleusA_.at(i).y - b * sin(phiRP));
            double dij = dx * dx + dy * dy;

            bool collided;
            if (!gaussianWounding) {
                collided = dij < d2;
            } else {
                double p = G * exp(-G * dij / d2);  // Gaussian profile
                double ran = random_ptr_->genrand64_real1();
                collided = ran < p;
            }

            if (collided) {
                foutNcoll << (nucleusB_.at(j).x + nucleusA_.at(i).x) / 2. << " "
                          << (nucleusB_.at(j).y + nucleusA_.at(i).y) / 2.
                          << endl;
                Ncoll++;
                nucleusB_.at(j).collided = 1;
                nucleusA_.at(i).collided = 1;
            }
        }
    }

    foutNcoll.close();
}

// This function compute the collision geometry quantities, such as
// Npart, Ncoll, averageQs, etc.
// Determines Npart/Ncoll from the (already-sampled) nucleon positions in
// nucleusA_/nucleusB_, writes NcollList*.dat/NpartList*.dat, and sets
// param->event.Npart. Returns false (having called param->event.success = 0) if
// useFixedNpart is set and this event's Npart doesn't match, signaling the
// caller to abort and resample.
bool Init::determineNpartAndNcoll(Parameters *param, int &Npart, int &Ncoll) {
    const double d2 = param->collision.sigmaNN * mbToFm2 / M_PI;  // in fm^2
    const double b = param->event.b;
    const double phiRP = param->event.phiRP;
    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();

    // Determine Npart, Ncoll. Do this only during the first stage, as in
    // the 2nd stage nuclei are shifted to b=0
    if (!param->nucleus.useSmoothNucleus) {
        computeNcollList(param, d2, b, phiRP, Ncoll);

        stringstream strNpart_name;
        strNpart_name << "NpartList" << param->event.eventId << ".dat";
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

        param->event.Npart = Npart;

        if (param->collision.useFixedNpart != 0
            && Npart != param->collision.useFixedNpart) {
            messager_ << "[Init::computeCollisionGeometryQuantities]: "
                         "Npart = "
                      << Npart
                      << " does not match the requested fixed "
                         "Npart = "
                      << param->collision.useFixedNpart << "; resampling.";
            messager_.flush("info");
            param->event.success = 0;
            return false;
        }
    } else {
        // Smooth nucleus
        Npart = 2;
        Ncoll = 2;
        param->event.Npart = Npart;
    }
    return true;
}

// Sets param's running-coupling alpha_s from whichever Qs choice
// param->coupling.runWithQs selects (max/min/avg), or a fixed value
// when running coupling is disabled or alpha_s runs with k_T instead (handled
// per-cell elsewhere via computeRunningCouplingGfactor, which shares
// RunningCoupling.h's computeAlphaS() with this function).
void Init::computeAndSetRunningAlphaS(Parameters *param) {
    double alphas = 0.;
    if (param->coupling.runningCoupling && !param->coupling.runWithKt) {
        // Uses the same regularized formula (and the same muZero/c) as
        // Evolution::computeRunningCouplingGfactor()/MyEigen, instead of
        // the unregularized formula this used to hardcode inline -- so
        // this diagnostic/event-acceptance alpha_s always matches the one
        // actually used during evolution, and can no longer go negative
        // or singular at a small average Qs (validationErrors() already
        // guarantees LambdaQCD < muZero whenever running coupling is on).
        if (param->coupling.runWithQs == 2) {
            messager_
                << "[Init::computeCollisionGeometryQuantities]: running with "
                << param->coupling.runningCouplingQsFactor << " Q_s(max)";
            messager_.flush("info");
            alphas = computeAlphaS(
                param->coupling.mu0, param->coupling.c,
                param->coupling.LambdaQCD, param->coupling.nFlavors,
                param->coupling.runningCouplingQsFactor
                    * param->event.averageQs);
            messager_ << "[Init::computeCollisionGeometryQuantities]: alpha_s("
                      << param->coupling.runningCouplingQsFactor
                      << " Qs_max)=" << alphas;
            messager_.flush("info");
        } else if (param->coupling.runWithQs == 0) {
            messager_
                << "[Init::computeCollisionGeometryQuantities]: running with "
                << param->coupling.runningCouplingQsFactor << " Q_s(min)";
            messager_.flush("info");
            alphas = computeAlphaS(
                param->coupling.mu0, param->coupling.c,
                param->coupling.LambdaQCD, param->coupling.nFlavors,
                param->coupling.runningCouplingQsFactor
                    * param->event.averageQsmin);
            messager_ << "[Init::computeCollisionGeometryQuantities]: alpha_s("
                      << param->coupling.runningCouplingQsFactor
                      << " Qs_min)=" << alphas;
            messager_.flush("info");
        } else if (param->coupling.runWithQs == 1) {
            messager_
                << "[Init::computeCollisionGeometryQuantities]: running with "
                << param->coupling.runningCouplingQsFactor << " <Q_s>";
            messager_.flush("info");
            alphas = computeAlphaS(
                param->coupling.mu0, param->coupling.c,
                param->coupling.LambdaQCD, param->coupling.nFlavors,
                param->coupling.runningCouplingQsFactor
                    * param->event.averageQsAvg);
            messager_ << "[Init::computeCollisionGeometryQuantities]: alpha_s("
                      << param->coupling.runningCouplingQsFactor
                      << " <Qs>)=" << alphas;
            messager_.flush("info");
        }
    } else if (param->coupling.runningCoupling && param->coupling.runWithKt) {
        messager_.info(
            "[Init::computeCollisionGeometryQuantities]: Multiplicity with "
            "running alpha_s(k_T)");
    } else {
        messager_.info(
            "[Init::computeCollisionGeometryQuantities]: Using fixed alpha_s");
        alphas = param->coupling.g * param->coupling.g / 4. / M_PI;
    }
    param->event.alphas = alphas;
}

void Init::computeCollisionGeometryQuantities(Lattice *lat, Parameters *param) {
    int Npart = 0;
    int Ncoll = 0;

    const double L = param->lattice.L;
    const int N = param->lattice.size;
    const double a = L / N;  // lattice spacing in fm
    const double b = param->event.b;
    const double phiRP = param->event.phiRP;

    // Determine Npart, Ncoll only during the first stage, as in the 2nd
    // stage nuclei are shifted to b=0.
    if (!determineNpartAndNcoll(param, Npart, Ncoll)) {
        return;
    }

    double averageQs = 0.;
    double averageQs2 = 0.;
    double averageQs2Avg = 0.;
    double averageQs2min = 0.;
    double averageQs2min2 = 0.;
    double Tpp = 0.;
    int count = 0;
    scanCollisionGeometry(
        lat, param, N, a, b, phiRP, averageQs, averageQs2, averageQs2Avg,
        averageQs2min, averageQs2min2, Tpp, count);

    if (count == 0) {
        param->event.averageQs = 0.;
        param->event.averageQsAvg = 0.;
        param->event.averageQsmin = 0.;
        param->event.Tpp = Tpp;
        param->event.success = 0;
        messager_.warning(
            "[Init::computeCollisionGeometryQuantities]: Rejected event -- "
            "no overlap region (count=0).");
        return;
    }

    averageQs /= static_cast<double>(count) + smallEps;
    averageQs2 /= static_cast<double>(count) + smallEps;
    averageQs2Avg /= static_cast<double>(count) + smallEps;
    averageQs2min /= static_cast<double>(count) + smallEps;

    param->event.averageQs = sqrt(averageQs2);
    param->event.averageQsAvg = sqrt(averageQs2Avg);
    param->event.averageQsmin = sqrt(averageQs2min);
    param->event.Tpp = Tpp;

    logCollisionGeometryQuantities(
        param, Npart, Ncoll, Tpp, a, averageQs2, averageQs2Avg, averageQs2min,
        averageQs2min2, count);

    computeAndSetRunningAlphaS(param);

    // With running alpha_s(k_T) the coupling is evaluated per k_T bin in the
    // multiplicity, and computeAndSetRunningAlphaS() leaves alphas at 0.
    const bool alphasOk =
        param->event.alphas > 0
        || (param->coupling.runningCoupling && param->coupling.runWithKt);
    if (param->event.averageQs > 0 && param->event.averageQsAvg > 0
        && averageQs2 > 0 && param->event.averageQsmin > 0 && averageQs2Avg > 0
        && alphasOk && Npart >= 2
        && averageQs2min2 * a * a / hbarc / hbarc
               > param->colorCharge.minimumQs2ST) {
        param->event.success = 1;
        writeUsedParametersFile(param, phiRP, Npart, Ncoll);
    } else {
        param->event.success = 0;
    }
    if (averageQs2min2 * a * a / hbarc / hbarc
        < param->colorCharge.minimumQs2ST) {
        messager_ << "[Init::computeCollisionGeometryQuantities]: Rejected "
                     "event -- Qsmin^2 S_T="
                  << averageQs2min2 * a * a / hbarc / hbarc
                  << " is below "
                     "the "
                     "minimum ("
                  << param->colorCharge.minimumQs2ST << ").";
        messager_.flush("warning");
    }

    writeNgluonEstimatorsFile(
        param, a, averageQs2, averageQs2Avg, averageQs2min2, count);
}

void Init::scanCollisionGeometry(
    Lattice *lat, Parameters *param, int N, double a, double b, double phiRP,
    double &averageQs, double &averageQs2, double &averageQs2Avg,
    double &averageQs2min, double &averageQs2min2, double &Tpp, int &count) {
    const double L = param->lattice.L;
    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();

    count = 0;
    Tpp = 0.;
    for (int ipos = 0; ipos < N * N; ipos++) {
        // loop over all positions
        int check = 0;
        int ix = lat->xFromPosition(ipos);
        int iy = lat->yFromPosition(ipos);
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

        int posA = lat->positionFromXY(ixA, iyA);
        int posB = lat->positionFromXY(ixB, iyB);

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
            averageQs2min2 += g2mu2A * param->colorCharge.QsMuRatio
                              * param->colorCharge.QsMuRatio / a / a * hbarc
                              * hbarc * param->coupling.g * param->coupling.g;
        } else {
            averageQs2min2 += g2mu2B * param->colorCharge.QsMuRatio
                              * param->colorCharge.QsMuRatio / a / a * hbarc
                              * hbarc * param->coupling.g * param->coupling.g;
        }

        for (int i = 0; i < A1; i++) {
            double xm = nucleusA_.at(i).x + b / 2. * cos(phiRP);
            double ym = nucleusA_.at(i).y + b / 2. * sin(phiRP);
            double r = sqrt((x - xm) * (x - xm) + (y - ym) * (y - ym));

            if (r < sqrt(param->collision.sigmaNN * mbToFm2 / M_PI)
                && nucleusA_.at(i).collided == 1) {
                check = 1;
            }
        }

        for (int i = 0; i < A2; i++) {
            double xm = nucleusB_.at(i).x - b / 2. * cos(phiRP);
            double ym = nucleusB_.at(i).y - b / 2. * sin(phiRP);
            double r = sqrt((x - xm) * (x - xm) + (y - ym) * (y - ym));

            if (r < sqrt(param->collision.sigmaNN * mbToFm2 / M_PI)
                && nucleusB_.at(i).collided == 1 && check == 1) {
                check = 2;
            }
        }

        // A smooth nucleus has no nucleons to test against: any cell where
        // both shifted thickness profiles are nonzero is in the overlap.
        if (param->nucleus.useSmoothNucleus && TpA > 0. && TpB > 0.) {
            check = 2;
        }

        if (check == 2) {
            if (g2mu2B > g2mu2A) {
                averageQs += sqrt(
                    g2mu2B * param->colorCharge.QsMuRatio
                    * param->colorCharge.QsMuRatio / a / a * hbarc * hbarc
                    * param->coupling.g * param->coupling.g);
                averageQs2 += g2mu2B * param->colorCharge.QsMuRatio
                              * param->colorCharge.QsMuRatio / a / a * hbarc
                              * hbarc * param->coupling.g * param->coupling.g;
                averageQs2min += g2mu2A * param->colorCharge.QsMuRatio
                                 * param->colorCharge.QsMuRatio / a / a * hbarc
                                 * hbarc * param->coupling.g
                                 * param->coupling.g;
            } else {
                averageQs += sqrt(
                    g2mu2A * param->colorCharge.QsMuRatio
                    * param->colorCharge.QsMuRatio / a / a * hbarc * hbarc
                    * param->coupling.g * param->coupling.g);
                averageQs2 += g2mu2A * param->colorCharge.QsMuRatio
                              * param->colorCharge.QsMuRatio / a / a * hbarc
                              * hbarc * param->coupling.g * param->coupling.g;
                averageQs2min += g2mu2B * param->colorCharge.QsMuRatio
                                 * param->colorCharge.QsMuRatio / a / a * hbarc
                                 * hbarc * param->coupling.g
                                 * param->coupling.g;
            }
            averageQs2Avg += (g2mu2B * param->colorCharge.QsMuRatio
                                  * param->colorCharge.QsMuRatio
                              + g2mu2A * param->colorCharge.QsMuRatio
                                    * param->colorCharge.QsMuRatio)
                             / 2. / a / a * hbarc * hbarc * param->coupling.g
                             * param->coupling.g;
            count++;
        }
        // compute T_pp
        Tpp += TpA * TpB * a * a / hbarc / hbarc / hbarc
               / hbarc;  // now this quantity is in fm^-2
                         // remember: Tp is in GeV^2
    }
}

void Init::logCollisionGeometryQuantities(
    Parameters *param, int Npart, int Ncoll, double Tpp, double a,
    double averageQs2, double averageQs2Avg, double averageQs2min,
    double averageQs2min2, int count) {
    messager_ << "[Init::computeCollisionGeometryQuantities]: N_part=" << Npart;
    messager_.flush("info");
    messager_ << "[Init::computeCollisionGeometryQuantities]: N_coll=" << Ncoll;
    messager_.flush("info");
    messager_ << "[Init::computeCollisionGeometryQuantities]: T_pp("
              << param->event.b << " fm) = " << Tpp << " 1/fm^2";
    messager_.flush("info");
    messager_ << "[Init::computeCollisionGeometryQuantities]: Q_s^2(max) S_T = "
              << averageQs2 * a * a / hbarc / hbarc
                     * static_cast<double>(count);
    messager_.flush("info");
    messager_ << "[Init::computeCollisionGeometryQuantities]: Q_s^2(avg) S_T = "
              << averageQs2Avg * a * a / hbarc / hbarc
                     * static_cast<double>(count);
    messager_.flush("info");
    messager_ << "[Init::computeCollisionGeometryQuantities]: Q_s^2(min) S_T = "
              << averageQs2min * a * a / hbarc / hbarc
                     * static_cast<double>(count);
    messager_.flush("info");
    messager_ << "[Init::computeCollisionGeometryQuantities]: Q_s^2(min) S_T "
                 "(full lattice) = "
              << averageQs2min2 * a * a / hbarc / hbarc;
    messager_.flush("info");

    messager_ << "[Init::computeCollisionGeometryQuantities]: Area = "
              << a * a * count << " fm^2";
    messager_.flush("info");

    messager_
        << "[Init::computeCollisionGeometryQuantities]: Average Qs(max) = "
        << param->event.averageQs << " GeV";
    messager_.flush("info");
    messager_
        << "[Init::computeCollisionGeometryQuantities]: Average Qs(avg) = "
        << param->event.averageQsAvg << " GeV";
    messager_.flush("info");
    messager_
        << "[Init::computeCollisionGeometryQuantities]: Average Qs(min) = "
        << param->event.averageQsmin << " GeV";
    messager_.flush("info");

    messager_
        << "[Init::computeCollisionGeometryQuantities]: resulting Y(Qs(max)*"
        << param->colorCharge.xQsFactor << ") = "
        << log(0.01
               / (param->event.averageQs * param->colorCharge.xQsFactor
                  / param->collision.sqrtS));
    messager_.flush("info");
    messager_
        << "[Init::computeCollisionGeometryQuantities]: resulting Y(Qs(avg)*"
        << param->colorCharge.xQsFactor << ") = "
        << log(0.01
               / (param->event.averageQsAvg * param->colorCharge.xQsFactor
                  / param->collision.sqrtS));
    messager_.flush("info");
    messager_
        << "[Init::computeCollisionGeometryQuantities]: resulting Y(Qs(min)*"
        << param->colorCharge.xQsFactor << ") =  "
        << log(0.01
               / (param->event.averageQsmin * param->colorCharge.xQsFactor
                  / param->collision.sqrtS));
    messager_.flush("info");
}

void Init::writeUsedParametersFile(
    Parameters *param, double phiRP, int Npart, int Ncoll) {
    stringstream strup_name;
    strup_name << "usedParameters" << param->event.eventId << ".dat";
    string up_name;
    up_name = strup_name.str();

    // comment lines, so the file stays a valid input file
    ofstream fout1(up_name.c_str(), std::ios::app);
    fout1 << "# Collision geometry of this event:" << endl;
    fout1 << "# b = " << param->event.b << " fm" << endl;
    fout1 << "# phiRP = " << phiRP << endl;
    fout1 << "# Npart = " << Npart << endl;
    fout1 << "# Ncoll = " << Ncoll << endl;
    if (param->coupling.runningCoupling) {
        if (param->coupling.runWithQs == 2)
            fout1 << "# <Q_s>(max) = " << param->event.averageQs << endl;
        else if (param->coupling.runWithQs == 1)
            fout1 << "# <Q_s>(avg) = " << param->event.averageQsAvg << endl;
        else if (param->coupling.runWithQs == 0)
            fout1 << "# <Q_s>(min) = " << param->event.averageQsmin << endl;
        fout1 << "# alpha_s(" << param->coupling.runningCouplingQsFactor
              << " <Q_s>) = " << param->event.alphas << endl;
    } else
        fout1 << "# using fixed coupling alpha_s=" << param->event.alphas
              << endl;
    fout1.close();
}

void Init::writeNgluonEstimatorsFile(
    Parameters *param, double a, double averageQs2, double averageQs2Avg,
    double averageQs2min2, int count) {
    stringstream strNEst_name;
    strNEst_name << "NgluonEstimators" << param->event.eventId << ".dat";
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

std::vector<double> Init::computeWilsonLineMomentumKernel(
    int N, int sites, double m, double UVdamp) {
    std::vector<double> momentumKernel(static_cast<std::size_t>(sites));
#pragma omp parallel for
    for (int pos = 0; pos < sites; ++pos) {
        const int i = latticeX(pos, N);
        const int j = latticeY(pos, N);
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
    return momentumKernel;
}

void Init::computeWilsonLineColorChargeScales(
    Lattice *lat, int sites, double g, double invNy,
    std::vector<double> &colorChargeScaleA,
    std::vector<double> &colorChargeScaleB) {
    colorChargeScaleA.resize(static_cast<std::size_t>(sites));
    colorChargeScaleB.resize(static_cast<std::size_t>(sites));
#pragma omp parallel for
    for (int pos = 0; pos < sites; ++pos) {
        colorChargeScaleA[static_cast<std::size_t>(pos)] =
            g * sqrt(lat->cells[pos]->getg2mu2A() * invNy);
        colorChargeScaleB[static_cast<std::size_t>(pos)] =
            g * sqrt(lat->cells[pos]->getg2mu2B() * invNy);
    }
}

void Init::setV(Lattice *lat, Parameters *param, Random *random) {
    IPG_PROFILE_SCOPE("initialization.wilson_lines");
    messager_.info("[Init::setV]: Setting Wilson lines ...");
    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();

    const double d2 = param->collision.sigmaNN * mbToFm2 / M_PI;  // in fm^2
    const int N = param->lattice.size;
    const int Ny = param->colorCharge.Ny;
    const int sites = N * N;
    const int nn[2] = {N, N};
    const double L = param->lattice.L;
    const double a = L / N;  // lattice spacing in fm
    const double m = param->subnucleon.m * a / hbarc;
    const double g = param->coupling.g;
    const double invNy = 1. / static_cast<double>(Ny);
    double UVdamp = param->subnucleon.UVDamp;  // GeV^-1
    UVdamp = UVdamp / a * hbarc;

    // The lattice Poisson/UV kernel depends only on transverse momentum and
    // run parameters.  The historical implementation recomputed the same
    // sin/sqrt/exp expressions for every longitudinal sheet of both nuclei.
    std::vector<double> momentumKernel =
        computeWilsonLineMomentumKernel(N, sites, m, UVdamp);

    // rhoACoeffData owns the Nc2m1*sites backing storage; rhoACoeff is a
    // pointer-per-component view over it for FFT::fftnComplexArray's T**
    // interface.
    std::vector<complex<double>> rhoACoeffData(
        static_cast<std::size_t>(Nc2m1) * sites);
    std::vector<complex<double> *> rhoACoeff(Nc2m1);
    for (int i = 0; i < Nc2m1; i++) {
        rhoACoeff[i] = rhoACoeffData.data() + i * sites;
    }

    auto applyMomentumKernel = [&]() {
#pragma omp parallel for
        for (int n = 0; n < Nc2m1; ++n) {
            complex<double> *rho = rhoACoeff[n];
            for (int pos = 0; pos < sites; ++pos) {
                rho[pos] *= momentumKernel[static_cast<std::size_t>(pos)];
            }
        }
    };

    // Reuse the bulk Gaussian buffers for every longitudinal sheet.  The
    // linear ordering matches the historical pos-major/color-minor gauss()
    // call sequence exactly.
    std::vector<double> gaussianField(
        static_cast<std::size_t>(sites) * static_cast<std::size_t>(Nc2m1));
    std::vector<double> gaussianScratch;
    gaussianScratch.reserve(5 * ((gaussianField.size() + 1) / 2));

    // g2mu2 and Ny are fixed throughout Wilson-line construction.  Cache the
    // color-independent site scale once for each nucleus instead of repeating
    // the same sqrt in every longitudinal sheet.
    std::vector<double> colorChargeScaleA;
    std::vector<double> colorChargeScaleB;
    computeWilsonLineColorChargeScales(
        lat, sites, g, invNy, colorChargeScaleA, colorChargeScaleB);

    auto fillColorCharge = [&](const std::vector<double> &scale) {
        {
            IPG_PROFILE_SCOPE("initialization.wilson_random.gauss");
            random->gaussBulk(
                gaussianField.data(), gaussianField.size(), gaussianScratch);
        }
        {
            IPG_PROFILE_SCOPE("initialization.wilson_random.scale");
#pragma omp parallel for
            for (int pos = 0; pos < sites; ++pos) {
                const double localScale = scale[static_cast<std::size_t>(pos)];
                const std::size_t base = static_cast<std::size_t>(pos)
                                         * static_cast<std::size_t>(Nc2m1);
                for (int n = 0; n < Nc2m1; ++n) {
                    rhoACoeff[n][pos] =
                        localScale
                        * gaussianField[base + static_cast<std::size_t>(n)];
                }
            }
        }
    };

    // loop over longitudinal direction, once for each nucleus
    auto evolveNucleusWilsonLine =
        [&](const std::vector<double> &colorChargeScale,
            std::vector<Matrix> &U) {
            for (int k = 0; k < Ny; k++) {
                {
                    IPG_PROFILE_SCOPE("initialization.wilson_random");
                    fillColorCharge(colorChargeScale);
                }

                fft_.fftnComplexArray(
                    rhoACoeff.data(), rhoACoeff.data(), nn, 1, Nc2m1);

                {
                    IPG_PROFILE_SCOPE("initialization.wilson_Poisson");
                    applyMomentumKernel();
                }

                fft_.fftnComplexArray(
                    rhoACoeff.data(), rhoACoeff.data(), nn, -1, Nc2m1);

                {
                    IPG_PROFILE_SCOPE("initialization.wilson_exponent");
#pragma omp parallel
                    {
                        std::vector<double> in(Nc2m1, 0.);
                        Matrix temp(1.);
                        Matrix tempNew(0.);

#pragma omp for
                        for (int pos = 0; pos < sites; pos++) {
                            for (int aa = 0; aa < Nc2m1; aa++) {
                                // expmCoeff calculates exp(i in[a] t[a]), so
                                // multiply by -1 (not -i).
                                in[aa] = -(rhoACoeff[aa][pos]).real();
                            }
                            tempNew = Matrix::fromAlgebraExponent(in);
                            temp = tempNew * U[pos];
                            U[pos] = temp;
                        }
                    }
                }
            }
        };

    evolveNucleusWilsonLine(colorChargeScaleA, lat->U);
    evolveNucleusWilsonLine(colorChargeScaleB, lat->U2);

    if (param->output.writeOutputs == 5) {
        WilsonLineIO().writeTrainingData(lat, param);
    }

    // output U
    if (param->wilsonLines.writeWilsonLines > 0
        && (param->jimwlk.saveSnapshots || !param->jimwlk.enabled)) {
        double x_projectile, x_target;
        if (param->jimwlk.enabled) {
            x_projectile = x_target = param->jimwlk.initialX;
        } else {
            if (param->colorCharge.useFluctuatingX) {
                // Initial condition does not correspond to a fixed x
                x_projectile = x_target = -1;
            } else {
                x_projectile = 0.01 * std::exp(-param->colorCharge.rapidityA);
                x_target = 0.01 * std::exp(-param->colorCharge.rapidityB);
            }
        }
        WilsonLineIO io;
        io.write(lat, param, NucleusRole::Projectile, x_projectile);
        io.write(lat, param, NucleusRole::Target, x_target);
    }

    messager_ << "[Init::setV]: Wilson lines V_A and V_B set on rank "
              << param->run.MPIRank << ". ";
    messager_.flush("info");
}

void Init::sampleImpactParameter(Parameters *param) {
    const double bmin = param->collision.bMin;
    const double bmax = param->collision.bMax;
    double b = 0.;
    double xb = random_ptr_->genrand64_real1();
    if (!param->collision.useNucleus) {
        // use b=0 fm for the constant g^2 mu case. Deferred flush: this
        // message's tag also covers the shared "b = ..." line below, same
        // as the other two branches.
        messager_ << "[Init::sampleImpactParameter]: Setting b=0 for constant "
                     "color charge density case. ";
        b = 0;
    } else {
        if (param->collision.sampleBFromLinearDistribution) {
            // use a linear probability distribution for b if we are doing
            // nuclei
            messager_ << "[Init::sampleImpactParameter]: Sampling linearly "
                         "distributed b between "
                      << bmin << " and " << bmax << "fm. Found ";
            b = sqrt((bmax * bmax - bmin * bmin) * xb + bmin * bmin);
        } else {
            // use a uniform distribution instead
            messager_ << "[Init::sampleImpactParameter]: Sampling uniformly "
                         "distributed b between "
                      << bmin << " and " << bmax << "fm. Found ";
            b = (bmax - bmin) * xb + bmin;
        }
    }
    param->event.b = b;
    double phiRP = 0.;
    if (param->collision.rotateReactionPlane) {
        phiRP = 2 * M_PI * random_ptr_->genrand64_real2();
    }
    param->event.phiRP = phiRP;
    messager_ << "b = " << b << " fm, phi_RP = " << phiRP;
    messager_.flush("info");

    for (unsigned int i = 0; i < nucleusA_.size(); i++) {
        nucleusA_.at(i).collided = 0;
    }
    for (unsigned int i = 0; i < nucleusB_.size(); i++) {
        nucleusB_.at(i).collided = 0;
    }
}

void Init::init(
    Lattice *lat, Parameters *param, Random *random, Glauber *glauber,
    InitializationMethod init_method) {
    random_ptr_ = random;

    messager_.info("[Init::init]: Initializing fields ... ");

    if (!param->collision.useNucleus) {
        // No real collision geometry in the constant-g^2mu case: skip
        // sampleImpactParameter() entirely, but still give b/phi_RP the
        // same defaults it would have produced.
        param->event.b = 0.;
        param->event.phiRP = 0.;
        param->event.success = 1;
    } else {
        readNuclearQs(param);
    }

    // The configuration files are only used by sampleTAFromConfigFiles();
    // don't require them to exist for Woods-Saxon sampling.
    if (param->collision.useNucleus && param->nucleus.nucleonPositionsFromFile
        && init_method == InitializationMethod::SampleColorCharges) {
        readInNucleusConfigs(
            static_cast<int>(glauber->nucleusA1()),
            param->nucleus.lightNucleusOption,
            param->nucleus.polarizationProjectile,
            param->nucleus.polarizationProjectileJz, nucleonPosArrA_, param);
        readInNucleusConfigs(
            static_cast<int>(glauber->nucleusA2()),
            param->nucleus.lightNucleusOption,
            param->nucleus.polarizationTarget,
            param->nucleus.polarizationTargetJz, nucleonPosArrB_, param);
    }

    if (init_method == InitializationMethod::ReadWlineBinary
        or init_method == InitializationMethod::ReadWlineText) {
        // to read Wilson lines from file
        WilsonLineIO().read(
            lat, param,
            (init_method == InitializationMethod::ReadWlineBinary) ? 2 : 1);
        param->event.success = 1;
    } else {
        // to generate your own Wilson lines
        if (param->collision.useNucleus) {
            nucleusA_.clear();
            nucleusB_.clear();
            // populate the lists nucleusA_ and nucleusB_ with position data
            sampleTA(param, random, glauber);
        }
        setColorChargeDensity(lat, param, random, glauber);
        // sample color charges and find Wilson lines V_A and V_B
        setV(lat, param, random);
    }
}

void Init::shiftFieldsWithImpactParameter(Lattice *lat, Parameters *param) {
    messager_.info(
        "[Init::shiftFieldsWithImpactParameter]: Shifting fields with impact "
        "parameter...");
    const double b = param->event.b;
    const double phiRP = param->event.phiRP;
    messager_ << "[Init::shiftFieldsWithImpactParameter]: b = " << b
              << " fm, phi_RP = " << phiRP;
    messager_.flush("info");

    const int N = param->lattice.size;
    BufferLattice lat_tmp(param->lattice.size);
    for (int ipos = 0; ipos < N * N; ipos++) {
        lat_tmp.buffer1[ipos] = lat->U[ipos];
        lat_tmp.buffer2[ipos] = lat->U2[ipos];
    }

    const double L = param->lattice.L;
    const double a = L / N;  // lattice spacing in fm
    for (int ipos = 0; ipos < N * N; ipos++) {
        int ix = lat->xFromPosition(ipos);
        int iy = lat->yFromPosition(ipos);
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
            int posA = lat->positionFromXY(ixA, iyA);
            lat->U[ipos] = lat_tmp.buffer1[posA];
        }
        if (ixB < 0 || ixB >= N || iyB < 0 || iyB >= N) {
            lat->U2[ipos] = one_;
        } else {
            int posB = lat->positionFromXY(ixB, iyB);
            lat->U2[ipos] = lat_tmp.buffer2[posB];
        }
    }
}

void Init::generateNucleusConfiguration(
    Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
    double beta3, double beta4, double gamma, bool forceDminFlag, double d_min,
    double dR_np, double da_np, std::vector<ReturnValue> &nucleus) {
    if (std::abs(beta2) < 1e-15 && std::abs(beta4) < 1e-15
        && std::abs(beta3) < 1e-15 && std::abs(gamma) < 1e-15) {
        generateNucleusConfigurationWithWoodsSaxon(
            random, A, Z, a_WS, R_WS, d_min, dR_np, da_np, nucleus);
    } else {
        if (forceDminFlag) {
            generateNucleusConfigurationWithDeformedWoodsSaxonForceDmin(
                random, A, Z, a_WS, R_WS, beta2, beta3, beta4, gamma, d_min,
                dR_np, da_np, nucleus);
        } else {
            if (std::abs(gamma) > 1e-15) {
                generateNucleusConfigurationWithDeformedWoodsSaxon2(
                    random, A, Z, a_WS, R_WS, beta2, beta3, beta4, gamma, dR_np,
                    da_np, nucleus);
            } else {
                generateNucleusConfigurationWithDeformedWoodsSaxon(
                    random, A, Z, a_WS, R_WS, beta2, beta3, beta4, d_min, dR_np,
                    da_np, nucleus);
            }
        }
    }
}

void Init::generateNucleusConfigurationWithWoodsSaxon(
    Random *random, int A, int Z, double a_WS, double R_WS, double d_min,
    double dR_np, double da_np, std::vector<ReturnValue> &nucleus) {
    std::vector<double> r_array(A, 0.);
    std::vector<int> idx_array(A, 0);
    for (int i = 0; i < Z; i++) {
        r_array[i] = sampleRFromWoodsSaxon(random, a_WS, R_WS);
        idx_array[i] = i;
    }
    for (int i = Z; i < A; i++) {
        r_array[i] = sampleRFromWoodsSaxon(random, a_WS + da_np, R_WS + dR_np);
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

    recenterNucleus(x_array, y_array, z_array);

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

double Init::sampleRFromWoodsSaxon(
    Random *random, double a_WS, double R_WS) const {
    double rmaxCut = R_WS + 10. * a_WS;
    double r = 0.;
    do {
        r = rmaxCut * pow(random->genrand64_real3(), 1.0 / 3.0);
    } while (random->genrand64_real3() > fermiDistribution(r, R_WS, a_WS));
    return (r);
}

double Init::fermiDistribution(double r, double R_WS, double a_WS) const {
    double f = 1. / (1. + exp((r - R_WS) / a_WS));
    return (f);
}

void Init::generateNucleusConfigurationWithDeformedWoodsSaxon(
    Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
    double beta3, double beta4, double d_min, double dR_np, double da_np,
    std::vector<ReturnValue> &nucleus) {
    std::vector<double> r_array(A, 0.);
    std::vector<double> costheta_array(A, 0.);
    std::vector<int> idx_array(A, 0);
    for (int i = 0; i < Z; i++) {
        sampleRAndCosthetaFromDeformedWoodsSaxon(
            random, a_WS, R_WS, beta2, beta3, beta4, r_array[i],
            costheta_array[i]);
        idx_array[i] = i;
    }
    for (int i = Z; i < A; i++) {
        sampleRAndCosthetaFromDeformedWoodsSaxon(
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
    recenterNucleus(x_array, y_array, z_array);

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

void Init::generateNucleusConfigurationWithDeformedWoodsSaxonForceDmin(
    Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
    double beta3, double beta4, double gamma, double d_min, double dR_np,
    double da_np, std::vector<ReturnValue> &nucleus) {
    messager_
        << "[Init::generateNucleusConfigurationWithDeformedWoodsSaxonForceDmin]"
           ": "
           "Sampling nucleon position forcing dMin = "
        << d_min << " fm ...";
    messager_.flush("info");
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
                double y20 = sphericalHarmonics(2, costheta);
                double y30 = sphericalHarmonics(3, costheta);
                double y40 = sphericalHarmonics(4, costheta);
                double y22 = sphericalHarmonicsY22(costheta, phi);
                R_WS_theta =
                    R_WS_i
                    * (1.0 + beta2 * (cos(gamma) * y20 + sin(gamma) * y22)
                       + beta3 * y30 + beta4 * y40);
            } while (random->genrand64_real3()
                     > fermiDistribution(r, R_WS_theta, a_WS_i));
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

    recenterNucleus(x_array, y_array, z_array);

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

void Init::generateNucleusConfigurationWithDeformedWoodsSaxon2(
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
            double y20 = sphericalHarmonics(2, costheta);
            double y30 = sphericalHarmonics(3, costheta);
            double y40 = sphericalHarmonics(4, costheta);
            double y22 = sphericalHarmonicsY22(costheta, phi);
            R_WS_theta = R_WS_i
                         * (1.0 + beta2 * (cos(gamma) * y20 + sin(gamma) * y22)
                            + beta3 * y30 + beta4 * y40);
        } while (random->genrand64_real3()
                 > fermiDistribution(r, R_WS_theta, a_WS_i));
        double sintheta = sqrt(1. - costheta * costheta);
        x_array[i] = r * sintheta * cos(phi);
        y_array[i] = r * sintheta * sin(phi);
        z_array[i] = r * costheta;
    }

    recenterNucleus(x_array, y_array, z_array);

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

void Init::sampleRAndCosthetaFromDeformedWoodsSaxon(
    Random *random, double a_WS, double R_WS, double beta2, double beta3,
    double beta4, double &r, double &costheta) const {
    double rmaxCut = R_WS + 10. * a_WS;
    double R_WS_theta = R_WS;
    do {
        r = rmaxCut * pow(random->genrand64_real3(), 1.0 / 3.0);
        costheta = 1.0 - 2.0 * random->genrand64_real3();
        auto y20 = sphericalHarmonics(2, costheta);
        auto y30 = sphericalHarmonics(3, costheta);
        auto y40 = sphericalHarmonics(4, costheta);
        R_WS_theta = R_WS * (1.0 + beta2 * y20 + beta3 * y30 + beta4 * y40);
    } while (random->genrand64_real3()
             > fermiDistribution(r, R_WS_theta, a_WS));
}

double Init::sphericalHarmonics(int l, double ct) const {
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

double Init::sphericalHarmonicsY22(double ct, double phi) const {
    // Y2,2
    double ylm = 0.0;
    ylm = 1.0 - ct * ct;
    ylm *= cos(2. * phi);
    ylm *= 0.5462742152960397;  // sqrt(2*15)/4/pow(2*M_PI,0.5);
    return (ylm);
}

void Init::recenterNucleus(
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

void Init::recenterNucleus(std::vector<ReturnValue> &nucleus) {
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

void Init::rotateNucleus(
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

void Init::rotateNucleus3D(Random *random, std::vector<ReturnValue> &nucleus) {
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
    double sampleX = exp(mu + sigma * random->gauss());
    return (sampleX);
}

void Init::sampleQsNormalization(
    Random *random, Parameters *param, const int Nq,
    vector<double> &gauss_array) {
    const double QsSmearWidth = param->subnucleon.smearingWidth;
    gauss_array.assign(Nq, 1.);  // default norm = 1
    if (param->subnucleon.smearQs) {
        // introduce a log-normal distribution for Qs normalization
        // dividing by exp(0.5 sigma^2) to ensure the mean is 1
        // the varirance in this case is exp(sigma) - 1 for the log-normal
        // distribution
        for (int iq = 0; iq < Nq; iq++) {
            gauss_array[iq] =
                (exp(random->gauss(0, QsSmearWidth))
                 / exp(QsSmearWidth * QsSmearWidth / 2.));
        }
    }
}

int Init::sampleNumberOfPartons(Random *random, Parameters *param) {
    double NqBase = param->subnucleon.NqBase;
    int NqBaseInt = static_cast<int>(NqBase);
    double ran = random->genrand64_real2();
    int Nq = NqBaseInt;
    if (ran < NqBase - NqBaseInt) {
        Nq += 1;
    }
    Nq += random->poisson(param->subnucleon.NqFluc);
    return (std::max(1, Nq));
}

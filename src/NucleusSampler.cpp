// NucleusSampler.cpp is part of the IP-Glasma solver.

#include "NucleusSampler.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "Instrumentation.h"

void NucleusSampler::readConfigurations(
    Parameters *param, Glauber *glauber, Random *random) {
    // configurations that are already loaded are kept
    if (configsA_.empty()) {
        configsA_ = readConfigurationFile(
            static_cast<int>(glauber->nucleusA1()),
            param->nucleus.lightNucleusOption,
            param->nucleus.polarizationProjectile,
            param->nucleus.polarizationProjectileJz, random, param);
    }
    if (configsB_.empty()) {
        configsB_ = readConfigurationFile(
            static_cast<int>(glauber->nucleusA2()),
            param->nucleus.lightNucleusOption,
            param->nucleus.polarizationTarget,
            param->nucleus.polarizationTargetJz, random, param);
    }
}

// This function samples the nucleon positions inside the projectile and
// target nuclei. Both nuclei are centered at the origin.
Nuclei NucleusSampler::sample(
    Parameters *param, Random *random, Glauber *glauber) {
    IPG_PROFILE_SCOPE("initialization.sample_nuclei");
    messager_.info("[NucleusSampler::sample]: Sampling nucleon positions ... ");
    Nuclei nuclei;
    std::vector<ReturnValue> &nucleusA = nuclei.projectile;
    std::vector<ReturnValue> &nucleusB = nuclei.target;

    if (param->nucleus.nucleonPositionsFromFile) {
        nucleusA = sampleFromConfigurations(
            random, glauber->getGlauberData().projectile, configsA_);
        nucleusB = sampleFromConfigurations(
            random, glauber->getGlauberData().target, configsB_);
    } else {
        if (param->collision.nucleiToAverage > 1
            && (glauber->nucleusA1() == 1 || glauber->nucleusA2() == 1)) {
            messager_ << "[NucleusSampler::sample]: Averaging over nuclei is "
                         "not supported for collisions involving protons. "
                         "Exiting.";
            messager_.flush("error");
            exit(1);
        }
        nucleusA = sampleWoodsSaxon(random, glauber, NucleusRole::Projectile);
        nucleusB = sampleWoodsSaxon(random, glauber, NucleusRole::Target);
    }

    // global rotation of the nucleus
    applyPolarizationRotation(
        random, param->nucleus.polarizationProjectile, nucleusA);
    applyPolarizationRotation(
        random, param->nucleus.polarizationTarget, nucleusB);
    return nuclei;
}

std::vector<ReturnValue> NucleusSampler::sampleWoodsSaxon(
    Random *random, Glauber *glauber, NucleusRole role) {
    std::vector<ReturnValue> nucleus;
    const Nucleus &data = role == NucleusRole::Projectile
                              ? glauber->getGlauberData().projectile
                              : glauber->getGlauberData().target;
    ReturnValue rv;
    if (data.A == 1) {
        rv.x = 0.;
        rv.y = 0;
        rv.z = 0;
        rv.collided = 0;
        rv.proton = 1;
        nucleus.push_back(rv);
    } else if (data.A == 2) {
        // deuteron
        rv = glauber->sampleTARejection(random, role);
        // we sample the neutron proton distance, so distance to the center
        // needs to be divided by 2
        rv.x = rv.x / 2.;
        rv.y = rv.y / 2.;
        rv.z = 0.;
        rv.proton = 1;
        rv.collided = 0;
        nucleus.push_back(rv);
        // other nucleon is 180 degrees rotated:
        rv.x = -rv.x;
        rv.y = -rv.y;
        rv.z = -rv.z;
        rv.proton = 0;
        rv.collided = 0;
        nucleus.push_back(rv);
    } else {
        nucleus = generate(random, data);
    }
    return nucleus;
}

std::vector<ReturnValue> NucleusSampler::sampleFromConfigurations(
    Random *random, const Nucleus &data,
    const std::vector<std::vector<float>> &configs) {
    std::vector<ReturnValue> nucleus;
    ReturnValue rv;
    if (configs.size() > 0) {
        double ran2 = random->genrand64_real3();
        int nucleusNumber = static_cast<int>(ran2 * configs.size());
        messager_ << "[NucleusSampler::sampleFromConfigurations]: using "
                     "nucleus Number = "
                  << nucleusNumber;
        messager_.flush("info");
        for (int iA = 0; iA < data.A; iA++) {
            rv.x = configs[nucleusNumber][3 * iA];
            rv.y = configs[nucleusNumber][3 * iA + 1];
            rv.z = configs[nucleusNumber][3 * iA + 2];
            rv.collided = 0;
            nucleus.push_back(rv);
        }
        assignProtons(random, nucleus, data.Z);
        recenter(nucleus);
    } else {
        // no configurations, sample with Woods-Saxon
        messager_ << "[NucleusSampler::sampleFromConfigurations]: "
                     "configuration file for A = "
                  << data.A
                  << " is not available, generate the nucleus configuration "
                  << "using Woods-Saxon distribution instead.";
        messager_.flush("info");

        nucleus = generate(random, data);
    }
    return nucleus;
}

std::vector<std::vector<float>> NucleusSampler::readConfigurationFile(
    const int nucleusA, const int lightNucleusOption,
    const int polarizationFlag, const double polJz, Random *random,
    Parameters *param) {
    std::vector<std::vector<float>> nucleonPosArr;
    std::string path = param->nucleus.nuclearConfigurationsPath + "/";
    std::string fileName;
    bool readFlag = true;
    if (nucleusA == 2) {
        fileName = "DeuteronPol0Configs.bin.in";
        if (std::abs(std::abs(polJz) - 1.) < 1e-8)
            fileName = "DeuteronPolpm1Configs.bin.in";
        if (polarizationFlag == 0) {
            auto ran = random->genrand64_real1();
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

    if (!readFlag) return nucleonPosArr;

    int Nentry = 3;
    if (nucleusA == 197 || nucleusA == 208) {
        Nentry = 4;
    }

    fileName = path + fileName;
    messager_ << "[NucleusSampler::readConfigurationFile]: read in nucleus "
                 "configurations from "
              << fileName;
    messager_.flush("info");
    std::ifstream inFile(fileName, std::ios::binary);
    if (!inFile) {
        messager_ << "[NucleusSampler::readConfigurationFile]: File "
                  << fileName << " not found. Exiting.";
        messager_.flush("error");
        exit(1);
    }
    while (true) {
        std::vector<float> tempPos;
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
    messager_ << "[NucleusSampler::readConfigurationFile]: read in "
              << nucleonPosArr.size() << " configurations.";
    messager_.flush("info");
    return nucleonPosArr;
}

std::vector<ReturnValue> NucleusSampler::generate(
    Random *random, const Nucleus &data) {
    const double beta2 = data.beta2;
    const double beta3 = data.beta3;
    const double beta4 = data.beta4;
    const double gamma = data.gamma;
    if (std::abs(beta2) < 1e-15 && std::abs(beta4) < 1e-15
        && std::abs(beta3) < 1e-15 && std::abs(gamma) < 1e-15) {
        return generateWoodsSaxon(random, data);
    } else {
        if (data.forceDminFlag) {
            return generateDeformedWoodsSaxonForceDmin(random, data);
        } else {
            if (std::abs(gamma) > 1e-15) {
                return generateTriaxialWoodsSaxon(random, data);
            } else {
                return generateDeformedWoodsSaxon(random, data);
            }
        }
    }
}

std::vector<ReturnValue> NucleusSampler::generateWoodsSaxon(
    Random *random, const Nucleus &data) {
    std::vector<ReturnValue> nucleus;
    const int A = data.A;
    const int Z = data.Z;
    const double a_WS = data.a_WS;
    const double R_WS = data.R_WS;
    const double d_min = data.d_min;
    const double dR_np = data.dR_np;
    const double da_np = data.da_np;
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

    recenter(x_array, y_array, z_array);

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
    return nucleus;
}

std::vector<ReturnValue> NucleusSampler::generateDeformedWoodsSaxon(
    Random *random, const Nucleus &data) {
    std::vector<ReturnValue> nucleus;
    const int A = data.A;
    const int Z = data.Z;
    const double a_WS = data.a_WS;
    const double R_WS = data.R_WS;
    const double beta2 = data.beta2;
    const double beta3 = data.beta3;
    const double beta4 = data.beta4;
    const double d_min = data.d_min;
    const double dR_np = data.dR_np;
    const double da_np = data.da_np;
    std::vector<double> r_array(A, 0.);
    std::vector<double> costheta_array(A, 0.);
    std::vector<int> idx_array(A, 0);
    for (int i = 0; i < Z; i++) {
        const RadiusAndCosTheta position =
            sampleRAndCosthetaFromDeformedWoodsSaxon(
                random, a_WS, R_WS, beta2, beta3, beta4);
        r_array[i] = position.r;
        costheta_array[i] = position.cosTheta;
        idx_array[i] = i;
    }
    for (int i = Z; i < A; i++) {
        const RadiusAndCosTheta position =
            sampleRAndCosthetaFromDeformedWoodsSaxon(
                random, a_WS + da_np, R_WS + dR_np, beta2, beta3, beta4);
        r_array[i] = position.r;
        costheta_array[i] = position.cosTheta;
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
    recenter(x_array, y_array, z_array);

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
    return nucleus;
}

std::vector<ReturnValue> NucleusSampler::generateDeformedWoodsSaxonForceDmin(
    Random *random, const Nucleus &data) {
    std::vector<ReturnValue> nucleus;
    const int A = data.A;
    const int Z = data.Z;
    const double a_WS = data.a_WS;
    const double R_WS = data.R_WS;
    const double beta2 = data.beta2;
    const double beta3 = data.beta3;
    const double beta4 = data.beta4;
    const double gamma = data.gamma;
    const double d_min = data.d_min;
    const double dR_np = data.dR_np;
    const double da_np = data.da_np;
    messager_ << "[NucleusSampler::generateDeformedWoodsSaxonForceDmin]"
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

    recenter(x_array, y_array, z_array);

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
    return nucleus;
}

std::vector<ReturnValue> NucleusSampler::generateTriaxialWoodsSaxon(
    Random *random, const Nucleus &data) {
    std::vector<ReturnValue> nucleus;
    const int A = data.A;
    const int Z = data.Z;
    const double a_WS = data.a_WS;
    const double R_WS = data.R_WS;
    const double beta2 = data.beta2;
    const double beta3 = data.beta3;
    const double beta4 = data.beta4;
    const double gamma = data.gamma;
    const double dR_np = data.dR_np;
    const double da_np = data.da_np;
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

    recenter(x_array, y_array, z_array);

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
    return nucleus;
}

double NucleusSampler::sampleRFromWoodsSaxon(
    Random *random, double a_WS, double R_WS) {
    double rmaxCut = R_WS + 10. * a_WS;
    double r = 0.;
    do {
        r = rmaxCut * pow(random->genrand64_real3(), 1.0 / 3.0);
    } while (random->genrand64_real3() > fermiDistribution(r, R_WS, a_WS));
    return (r);
}

RadiusAndCosTheta NucleusSampler::sampleRAndCosthetaFromDeformedWoodsSaxon(
    Random *random, double a_WS, double R_WS, double beta2, double beta3,
    double beta4) {
    double r = 0.;
    double costheta = 0.;
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
    return {r, costheta};
}

double NucleusSampler::fermiDistribution(double r, double R_WS, double a_WS) {
    double f = 1. / (1. + exp((r - R_WS) / a_WS));
    return (f);
}

double NucleusSampler::sphericalHarmonics(int l, double ct) {
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

double NucleusSampler::sphericalHarmonicsY22(double ct, double phi) {
    // Y2,2
    double ylm = 0.0;
    ylm = 1.0 - ct * ct;
    ylm *= cos(2. * phi);
    ylm *= 0.5462742152960397;  // sqrt(2*15)/4/pow(2*M_PI,0.5);
    return (ylm);
}

void NucleusSampler::recenter(
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

void NucleusSampler::recenter(std::vector<ReturnValue> &nucleus) {
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

void NucleusSampler::assignProtons(
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

void NucleusSampler::rotate(
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

void NucleusSampler::rotateRandomly(
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

void NucleusSampler::applyPolarizationRotation(
    Random *random, int polarizationFlag, std::vector<ReturnValue> &nucleus) {
    if (polarizationFlag == 0) {
        rotateRandomly(random, nucleus);
    } else if (polarizationFlag == 1) {
        // longitudinal polarization only rotates phi randomly
        double phi = 2. * M_PI * random->genrand64_real3();
        double theta = 0;
        rotate(phi, theta, nucleus);
    } else if (polarizationFlag == 2) {
        // transverse polarization rotates J to +y axis
        double phi = M_PI / 2;
        double theta = M_PI / 2;
        rotate(phi, theta, nucleus);
    }
}

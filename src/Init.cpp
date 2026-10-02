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
#include "NuclearQsTable.h"
#include "NucleonModel.h"
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

void Init::sampleTA(Parameters *param, Random *random, Glauber *glauber) {
    nucleusSampler_.sample(param, random, glauber, nucleusA_, nucleusB_);
}

// Q_s as a function of \sum T_p and y (new in this version of the code -
// v1.2 and up)
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
            Qs = sqrt(qsTable_.qs2(Tp, localrapidity));
        } else {
            xVal = Qs * param->colorCharge.xQsFactor / param->collision.sqrtS
                   * exp(ySign * yIn);
            if (xVal == 0)
                Qs = 0.;
            else
                Qs = sqrt(qsTable_.qs2(Tp, 0.))
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
            qsTable_.qs2(lat->cells[ipos]->getTpA(), rapidityA)
            / param->colorCharge.QsMuRatio / param->colorCharge.QsMuRatio * a
            * a / hbarc / hbarc / param->coupling.g
            / param->coupling.g);  // lattice units? check

        // nucleus B
        lat->cells[ipos]->setg2mu2B(
            qsTable_.qs2(lat->cells[ipos]->getTpB(), rapidityB)
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
    sampleNucleonProfiles(param, random);

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

void Init::sampleNucleonProfiles(Parameters *param, Random *random) {
    const std::unique_ptr<NucleonModel> model = NucleonModel::create(*param);
    profilesA_.clear();
    profilesB_.clear();
    for (const ReturnValue &nucleon : nucleusA_) {
        profilesA_.push_back(model->sample(*random, nucleon));
    }
    for (const ReturnValue &nucleon : nucleusB_) {
        profilesB_.push_back(model->sample(*random, nucleon));
    }
}

double Init::computeNucleonThicknessAtCell(
    const std::vector<std::unique_ptr<NucleonProfile>> &profiles, double x,
    double y, double nucleiInAverage) const {
    double Tp = 0.;
    for (const auto &profile : profiles) {
        Tp += profile->thickness(x, y) / nucleiInAverage;  // add up all T_p
    }
    return Tp;
}

void Init::computeThicknessFromNucleons(
    Lattice *lat, Parameters *param, double nucleiInAverage) {
    const int N = param->lattice.size;
    const double L = param->lattice.L;
    const double a = L / N;

#pragma omp parallel for
    for (int ipos = 0; ipos < N * N; ipos++) {
        // loop over all positions
        int iy = lat->yFromPosition(ipos);
        int ix = lat->xFromPosition(ipos);
        double x = -L / 2. + a * ix;
        double y = -L / 2. + a * iy;

        lat->cells[ipos]->setTpA(
            computeNucleonThicknessAtCell(profilesA_, x, y, nucleiInAverage));
        lat->cells[ipos]->setTpB(
            computeNucleonThicknessAtCell(profilesB_, x, y, nucleiInAverage));
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
        qsTable_.read(param->colorCharge.nucleusQsTableFileName);
    }

    // The configuration files are only used with nucleonPositionsFromFile;
    // don't require them to exist for Woods-Saxon sampling.
    if (param->collision.useNucleus && param->nucleus.nucleonPositionsFromFile
        && init_method == InitializationMethod::SampleColorCharges) {
        nucleusSampler_.readConfigurations(param, glauber, random);
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

// CollisionGeometry.cpp is part of the IP-Glasma solver.

#include "CollisionGeometry.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "PhysConst.h"
#include "RunningCoupling.h"

using PhysConst::hbarc;
using PhysConst::mbToFm2;
using PhysConst::smallEps;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;

void CollisionGeometry::sampleImpactParameter(
    Parameters *param, Random *random) {
    const double bmin = param->collision.bMin;
    const double bmax = param->collision.bMax;
    double b = 0.;
    double xb = random->genrand64_real1();
    if (!param->collision.useNucleus) {
        // use b=0 fm for the constant g^2 mu case. Deferred flush: this
        // message's tag also covers the shared "b = ..." line below, same
        // as the other two branches.
        messager_ << "[CollisionGeometry::sampleImpactParameter]: Setting b=0 "
                     "for constant "
                     "color charge density case. ";
        b = 0;
    } else {
        if (param->collision.sampleBFromLinearDistribution) {
            // use a linear probability distribution for b if we are doing
            // nuclei
            messager_ << "[CollisionGeometry::sampleImpactParameter]: Sampling "
                         "linearly "
                         "distributed b between "
                      << bmin << " and " << bmax << "fm. Found ";
            b = sqrt((bmax * bmax - bmin * bmin) * xb + bmin * bmin);
        } else {
            // use a uniform distribution instead
            messager_ << "[CollisionGeometry::sampleImpactParameter]: Sampling "
                         "uniformly "
                         "distributed b between "
                      << bmin << " and " << bmax << "fm. Found ";
            b = (bmax - bmin) * xb + bmin;
        }
    }
    param->event.b = b;
    double phiRP = 0.;
    if (param->collision.rotateReactionPlane) {
        phiRP = 2 * M_PI * random->genrand64_real2();
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

void CollisionGeometry::computeQuantities(
    Lattice *lat, Parameters *param, Random *random) {
    const double L = param->lattice.L;
    const int N = param->lattice.size;
    const double a = L / N;  // lattice spacing in fm
    const double b = param->event.b;
    const double phiRP = param->event.phiRP;

    // Determine Npart, Ncoll only during the first stage, as in the 2nd
    // stage nuclei are shifted to b=0.
    const WoundedNucleons wounded = determineNpartAndNcoll(param, random);
    if (!wounded.accepted) return;
    const int Npart = wounded.Npart;
    const int Ncoll = wounded.Ncoll;

    const OverlapSums sums = scanOverlap(lat, param, N, a, b, phiRP);
    const int count = sums.count;
    const double Tpp = sums.Tpp;
    const double averageQs2min2 = sums.Qs2minFullLattice;

    if (count == 0) {
        param->event.averageQs = 0.;
        param->event.averageQsAvg = 0.;
        param->event.averageQsmin = 0.;
        param->event.Tpp = Tpp;
        param->event.success = 0;
        messager_.warning(
            "[CollisionGeometry::computeQuantities]: Rejected event -- "
            "no overlap region (count=0).");
        return;
    }

    const double averageQs2 =
        sums.Qs2 / (static_cast<double>(count) + smallEps);
    const double averageQs2Avg =
        sums.Qs2Avg / (static_cast<double>(count) + smallEps);
    const double averageQs2min =
        sums.Qs2min / (static_cast<double>(count) + smallEps);

    param->event.averageQs = sqrt(averageQs2);
    param->event.averageQsAvg = sqrt(averageQs2Avg);
    param->event.averageQsmin = sqrt(averageQs2min);
    param->event.Tpp = Tpp;

    logQuantities(
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
        messager_ << "[CollisionGeometry::computeQuantities]: Rejected "
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

// This function compute the collision geometry quantities, such as
// Npart, Ncoll, averageQs, etc.
// Determines Npart/Ncoll from the (already-sampled) nucleon positions in
// nucleusA_/nucleusB_, writes NcollList*.dat/NpartList*.dat, and sets
// param->event.Npart. Not accepted (having set param->event.success = 0) if
// useFixedNpart is set and this event's Npart doesn't match, signaling the
// caller to abort and resample.
WoundedNucleons CollisionGeometry::determineNpartAndNcoll(
    Parameters *param, Random *random) {
    WoundedNucleons wounded;
    const double d2 = param->collision.sigmaNN * mbToFm2 / M_PI;  // in fm^2
    const double b = param->event.b;
    const double phiRP = param->event.phiRP;
    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();

    // Determine Npart, Ncoll. Do this only during the first stage, as in
    // the 2nd stage nuclei are shifted to b=0
    if (!param->nucleus.useSmoothNucleus) {
        wounded.Ncoll = computeNcollList(param, random, d2, b, phiRP);

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

        wounded.Npart = 0;
        for (int i = 0; i < A1; i++) {
            if (nucleusA_.at(i).collided == 1) {
                wounded.Npart++;
            }
        }

        for (int i = 0; i < A2; i++) {
            if (nucleusB_.at(i).collided == 1) {
                wounded.Npart++;
            }
        }

        param->event.Npart = wounded.Npart;

        if (param->collision.useFixedNpart != 0
            && wounded.Npart != param->collision.useFixedNpart) {
            messager_ << "[CollisionGeometry::computeQuantities]: "
                         "Npart = "
                      << wounded.Npart
                      << " does not match the requested fixed "
                         "Npart = "
                      << param->collision.useFixedNpart << "; resampling.";
            messager_.flush("info");
            param->event.success = 0;
            wounded.accepted = false;
            return wounded;
        }
    } else {
        // Smooth nucleus
        wounded.Npart = 2;
        wounded.Ncoll = 2;
        param->event.Npart = wounded.Npart;
    }
    return wounded;
}

int CollisionGeometry::computeNcollList(
    Parameters *param, Random *random, double d2, double b, double phiRP) {
    int Ncoll = 0;
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
                double ran = random->genrand64_real1();
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
    return Ncoll;
}

OverlapSums CollisionGeometry::scanOverlap(
    Lattice *lat, Parameters *param, int N, double a, double b, double phiRP) {
    OverlapSums sums;
    const double L = param->lattice.L;
    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();

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

        // check each coordinate: a flattened index would let an
        // out-of-range coordinate alias a valid site
        auto onLattice = [N](int column, int row) {
            return column >= 0 && column < N && row >= 0 && row < N;
        };

        double g2mu2A = 0;
        double TpA = 0;
        if (onLattice(ixA, iyA)) {
            const int posA = lat->positionFromXY(ixA, iyA);
            g2mu2A = lat->cells[posA]->getg2mu2A();
            TpA = lat->cells[posA]->getTpA();
        }

        double g2mu2B = 0;
        double TpB = 0;
        if (onLattice(ixB, iyB)) {
            const int posB = lat->positionFromXY(ixB, iyB);
            g2mu2B = lat->cells[posB]->getg2mu2B();
            TpB = lat->cells[posB]->getTpB();
        }

        if (g2mu2B >= g2mu2A) {
            sums.Qs2minFullLattice += g2mu2A * param->colorCharge.QsMuRatio
                                      * param->colorCharge.QsMuRatio / a / a
                                      * hbarc * hbarc * param->coupling.g
                                      * param->coupling.g;
        } else {
            sums.Qs2minFullLattice += g2mu2B * param->colorCharge.QsMuRatio
                                      * param->colorCharge.QsMuRatio / a / a
                                      * hbarc * hbarc * param->coupling.g
                                      * param->coupling.g;
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
                sums.Qs += sqrt(
                    g2mu2B * param->colorCharge.QsMuRatio
                    * param->colorCharge.QsMuRatio / a / a * hbarc * hbarc
                    * param->coupling.g * param->coupling.g);
                sums.Qs2 += g2mu2B * param->colorCharge.QsMuRatio
                            * param->colorCharge.QsMuRatio / a / a * hbarc
                            * hbarc * param->coupling.g * param->coupling.g;
                sums.Qs2min += g2mu2A * param->colorCharge.QsMuRatio
                               * param->colorCharge.QsMuRatio / a / a * hbarc
                               * hbarc * param->coupling.g * param->coupling.g;
            } else {
                sums.Qs += sqrt(
                    g2mu2A * param->colorCharge.QsMuRatio
                    * param->colorCharge.QsMuRatio / a / a * hbarc * hbarc
                    * param->coupling.g * param->coupling.g);
                sums.Qs2 += g2mu2A * param->colorCharge.QsMuRatio
                            * param->colorCharge.QsMuRatio / a / a * hbarc
                            * hbarc * param->coupling.g * param->coupling.g;
                sums.Qs2min += g2mu2B * param->colorCharge.QsMuRatio
                               * param->colorCharge.QsMuRatio / a / a * hbarc
                               * hbarc * param->coupling.g * param->coupling.g;
            }
            sums.Qs2Avg += (g2mu2B * param->colorCharge.QsMuRatio
                                * param->colorCharge.QsMuRatio
                            + g2mu2A * param->colorCharge.QsMuRatio
                                  * param->colorCharge.QsMuRatio)
                           / 2. / a / a * hbarc * hbarc * param->coupling.g
                           * param->coupling.g;
            sums.count++;
        }
        // compute T_pp
        sums.Tpp += TpA * TpB * a * a / hbarc / hbarc / hbarc
                    / hbarc;  // now this quantity is in fm^-2
                              // remember: Tp is in GeV^2
    }
    return sums;
}

// Sets param's running-coupling alpha_s from whichever Qs choice
// param->coupling.runWithQs selects (max/min/avg), or a fixed value
// when running coupling is disabled or alpha_s runs with k_T instead (handled
// per-cell elsewhere via computeRunningCouplingGfactor, which shares
// RunningCoupling.h's computeAlphaS() with this function).
void CollisionGeometry::computeAndSetRunningAlphaS(Parameters *param) {
    double alphas = 0.;
    if (param->coupling.runningCoupling && !param->coupling.runWithKt) {
        // Uses the same regularized formula (and the same muZero/c) as
        // computeRunningCouplingGfactor()/MyEigen, instead of
        // the unregularized formula this used to hardcode inline -- so
        // this diagnostic/event-acceptance alpha_s always matches the one
        // actually used during evolution, and can no longer go negative
        // or singular at a small average Qs (validationErrors() already
        // guarantees LambdaQCD < muZero whenever running coupling is on).
        if (param->coupling.runWithQs == 2) {
            messager_ << "[CollisionGeometry::computeQuantities]: running with "
                      << param->coupling.runningCouplingQsFactor << " Q_s(max)";
            messager_.flush("info");
            alphas = computeAlphaS(
                param->coupling.mu0, param->coupling.c,
                param->coupling.LambdaQCD, param->coupling.nFlavors,
                param->coupling.runningCouplingQsFactor
                    * param->event.averageQs);
            messager_ << "[CollisionGeometry::computeQuantities]: alpha_s("
                      << param->coupling.runningCouplingQsFactor
                      << " Qs_max)=" << alphas;
            messager_.flush("info");
        } else if (param->coupling.runWithQs == 0) {
            messager_ << "[CollisionGeometry::computeQuantities]: running with "
                      << param->coupling.runningCouplingQsFactor << " Q_s(min)";
            messager_.flush("info");
            alphas = computeAlphaS(
                param->coupling.mu0, param->coupling.c,
                param->coupling.LambdaQCD, param->coupling.nFlavors,
                param->coupling.runningCouplingQsFactor
                    * param->event.averageQsmin);
            messager_ << "[CollisionGeometry::computeQuantities]: alpha_s("
                      << param->coupling.runningCouplingQsFactor
                      << " Qs_min)=" << alphas;
            messager_.flush("info");
        } else if (param->coupling.runWithQs == 1) {
            messager_ << "[CollisionGeometry::computeQuantities]: running with "
                      << param->coupling.runningCouplingQsFactor << " <Q_s>";
            messager_.flush("info");
            alphas = computeAlphaS(
                param->coupling.mu0, param->coupling.c,
                param->coupling.LambdaQCD, param->coupling.nFlavors,
                param->coupling.runningCouplingQsFactor
                    * param->event.averageQsAvg);
            messager_ << "[CollisionGeometry::computeQuantities]: alpha_s("
                      << param->coupling.runningCouplingQsFactor
                      << " <Qs>)=" << alphas;
            messager_.flush("info");
        }
    } else if (param->coupling.runningCoupling && param->coupling.runWithKt) {
        messager_.info(
            "[CollisionGeometry::computeQuantities]: Multiplicity with "
            "running alpha_s(k_T)");
    } else {
        messager_.info(
            "[CollisionGeometry::computeQuantities]: Using fixed alpha_s");
        alphas = param->coupling.g * param->coupling.g / 4. / M_PI;
    }
    param->event.alphas = alphas;
}

void CollisionGeometry::logQuantities(
    Parameters *param, int Npart, int Ncoll, double Tpp, double a,
    double averageQs2, double averageQs2Avg, double averageQs2min,
    double averageQs2min2, int count) {
    messager_ << "[CollisionGeometry::computeQuantities]: N_part=" << Npart;
    messager_.flush("info");
    messager_ << "[CollisionGeometry::computeQuantities]: N_coll=" << Ncoll;
    messager_.flush("info");
    messager_ << "[CollisionGeometry::computeQuantities]: T_pp("
              << param->event.b << " fm) = " << Tpp << " 1/fm^2";
    messager_.flush("info");
    messager_ << "[CollisionGeometry::computeQuantities]: Q_s^2(max) S_T = "
              << averageQs2 * a * a / hbarc / hbarc
                     * static_cast<double>(count);
    messager_.flush("info");
    messager_ << "[CollisionGeometry::computeQuantities]: Q_s^2(avg) S_T = "
              << averageQs2Avg * a * a / hbarc / hbarc
                     * static_cast<double>(count);
    messager_.flush("info");
    messager_ << "[CollisionGeometry::computeQuantities]: Q_s^2(min) S_T = "
              << averageQs2min * a * a / hbarc / hbarc
                     * static_cast<double>(count);
    messager_.flush("info");
    messager_ << "[CollisionGeometry::computeQuantities]: Q_s^2(min) S_T "
                 "(full lattice) = "
              << averageQs2min2 * a * a / hbarc / hbarc;
    messager_.flush("info");

    messager_ << "[CollisionGeometry::computeQuantities]: Area = "
              << a * a * count << " fm^2";
    messager_.flush("info");

    messager_ << "[CollisionGeometry::computeQuantities]: Average Qs(max) = "
              << param->event.averageQs << " GeV";
    messager_.flush("info");
    messager_ << "[CollisionGeometry::computeQuantities]: Average Qs(avg) = "
              << param->event.averageQsAvg << " GeV";
    messager_.flush("info");
    messager_ << "[CollisionGeometry::computeQuantities]: Average Qs(min) = "
              << param->event.averageQsmin << " GeV";
    messager_.flush("info");

    messager_ << "[CollisionGeometry::computeQuantities]: resulting Y(Qs(max)*"
              << param->colorCharge.xQsFactor << ") = "
              << log(0.01
                     / (param->event.averageQs * param->colorCharge.xQsFactor
                        / param->collision.sqrtS));
    messager_.flush("info");
    messager_ << "[CollisionGeometry::computeQuantities]: resulting Y(Qs(avg)*"
              << param->colorCharge.xQsFactor << ") = "
              << log(0.01
                     / (param->event.averageQsAvg * param->colorCharge.xQsFactor
                        / param->collision.sqrtS));
    messager_.flush("info");
    messager_ << "[CollisionGeometry::computeQuantities]: resulting Y(Qs(min)*"
              << param->colorCharge.xQsFactor << ") =  "
              << log(0.01
                     / (param->event.averageQsmin * param->colorCharge.xQsFactor
                        / param->collision.sqrtS));
    messager_.flush("info");
}

void CollisionGeometry::writeUsedParametersFile(
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

void CollisionGeometry::writeNgluonEstimatorsFile(
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

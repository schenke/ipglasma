// RunningCoupling.cpp is part of the IP-Glasma solver.

#include "RunningCoupling.h"

#include <algorithm>
#include <cmath>

#include "Lattice.h"
#include "Parameters.h"

using PhysConst::hbarc;

double computeRunningCouplingGfactor(
    Lattice *lat, Parameters *param, int pos, int N, double a, double g,
    double c, double muZero) {
    if (!param->coupling.runningCoupling) return 1.;

    const double lambdaQCD = param->coupling.LambdaQCD;
    const int nFlavors = param->coupling.nFlavors;

    if (param->coupling.runWithLocalQs) {
        // run with the local (in transverse plane) coupling
        const bool inBounds =
            lat->xFromPosition(pos) > 0 && lat->xFromPosition(pos) < N - 1
            && lat->yFromPosition(pos) > 0 && lat->yFromPosition(pos) < N - 1;
        const double g2mu2A = inBounds ? lat->cells[pos]->getg2mu2A() : 0;
        const double g2mu2B = inBounds ? lat->cells[pos]->getg2mu2B() : 0;

        double Qs = 0.;
        if (param->coupling.runWithQs == 2) {
            Qs = sqrt(
                std::max(g2mu2A, g2mu2B) * param->colorCharge.QsMuRatio
                * param->colorCharge.QsMuRatio / a / a * hbarc * hbarc
                * param->coupling.g * param->coupling.g);
        } else if (param->coupling.runWithQs == 0) {
            Qs = sqrt(
                std::min(g2mu2A, g2mu2B) * param->colorCharge.QsMuRatio
                * param->colorCharge.QsMuRatio / a / a * hbarc * hbarc
                * param->coupling.g * param->coupling.g);
        } else if (param->coupling.runWithQs == 1) {
            Qs = sqrt(
                (g2mu2A + g2mu2B) / 2. * param->colorCharge.QsMuRatio
                * param->colorCharge.QsMuRatio / a / a * hbarc * hbarc
                * param->coupling.g * param->coupling.g);
        }

        return computeRunningCouplingGfactorFromScale(
            g, muZero, c, lambdaQCD, nFlavors,
            param->coupling.runningCouplingQsFactor * Qs);
    } else {
        double averageQs = 0.;
        if (param->coupling.runWithQs == 0)
            averageQs = param->event.averageQsmin;
        else if (param->coupling.runWithQs == 1)
            averageQs = param->event.averageQsAvg;
        else if (param->coupling.runWithQs == 2)
            averageQs = param->event.averageQs;

        return computeRunningCouplingGfactorFromScale(
            g, muZero, c, lambdaQCD, nFlavors,
            param->coupling.runningCouplingQsFactor * averageQs);
    }
}

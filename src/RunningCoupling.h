#ifndef SRC_RUNNINGCOUPLING_H_
#define SRC_RUNNINGCOUPLING_H_

#include <cmath>

#include "PhysConst.h"

/**
 * Computes the running \f$\alpha_s\f$ at scale \p scale, via a
 * one-loop formula regularized to stay finite as \p scale \f$\to 0\f$:
 * \f[
 * \alpha_s(\text{scale}) = \frac{4\pi}{\beta_0
 * \ln\left[\left(\mu_0/\Lambda_{QCD}\right)^{2/c}
 * + \left(\text{scale}/\Lambda_{QCD}\right)^{2/c}\right]^c},
 * \f]
 * with the one-loop QCD beta-function coefficient
 * \f$\beta_0=(11N_c-2N_f)/3\f$ (\f$N_c\f$ from PhysConst::Nc). Shared
 * by every classical-evolution/hydro-output site that computes a
 * local-coupling factor from some \f$Q_s\f$- or \f$k_T\f$-like scale
 * (Evolution::computeRunningCouplingGfactor() and its callers,
 * MyEigen::flowVelocity4DImpl()); distinct from JIMWLK::getAlphas(),
 * which uses its own formula and \c LambdaQCD_jimwlk_ parameter for
 * the separate small-x evolution coupling.
 * \param[in] muZero \f$\mu_0\f$ [GeV], keeps \f$\alpha_s\f$ infrared
 * finite as \p scale \f$\to 0\f$.
 * \param[in] c Cutoff smoothness parameter.
 * \param[in] lambdaQCD \f$\Lambda_{QCD}\f$ [GeV].
 * \param[in] nFlavors Number of active quark flavors \f$N_f\f$.
 * \param[in] scale Momentum scale [GeV] (e.g. a local or averaged
 * \f$Q_s\f$, or a \f$k_T\f$-derived scale).
 * \return \f$\alpha_s(\text{scale})\f$.
 */
inline double computeAlphaS(
    double muZero, double c, double lambdaQCD, int nFlavors, double scale) {
    const double beta0 = (11. * PhysConst::Nc - 2. * nFlavors) / 3.;
    return 4. * M_PI
           / (beta0
              * log(
                  pow(pow(muZero / lambdaQCD, 2. / c)
                          + pow(scale / lambdaQCD, 2. / c),
                      c)));
}

/**
 * Computes the local-coupling factor \f$g^2/(4\pi\alpha_s(\text{scale}))\f$,
 * via computeAlphaS().
 * \param[in] g Coupling \f$g\f$.
 * \param[in] muZero \f$\mu_0\f$ [GeV]; forwarded to computeAlphaS().
 * \param[in] c Cutoff smoothness parameter; forwarded to
 * computeAlphaS().
 * \param[in] lambdaQCD \f$\Lambda_{QCD}\f$ [GeV]; forwarded to
 * computeAlphaS().
 * \param[in] nFlavors Number of active quark flavors \f$N_f\f$;
 * forwarded to computeAlphaS().
 * \param[in] scale Momentum scale [GeV]; forwarded to computeAlphaS().
 * \return The local-coupling factor.
 */
inline double computeRunningCouplingGfactorFromScale(
    double g, double muZero, double c, double lambdaQCD, int nFlavors,
    double scale) {
    return g * g
           / (4. * M_PI * computeAlphaS(muZero, c, lambdaQCD, nFlavors, scale));
}

#endif  // SRC_RUNNINGCOUPLING_H_

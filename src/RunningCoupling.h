#ifndef SRC_RUNNINGCOUPLING_H_
#define SRC_RUNNINGCOUPLING_H_

#include <cmath>
#include <vector>

#include "PhysConst.h"

class Lattice;
class Parameters;

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
 * (computeRunningCouplingGfactor() and its callers,
 * MyEigen::flowVelocity4DImpl()); distinct from JIMWLK::getAlphas(),
 * which uses its own formula and its own \c jimwlk.c/\c
 * jimwlk.LambdaQCD parameters for the separate small-x evolution
 * coupling, though it shares this \p nFlavors, since \f$N_f\f$ is the
 * same physical quantity in both.
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

/**
 * Computes the local-coupling factor \f$g^2/(4\pi\alpha_s)\f$ at one
 * cell, used to rescale \f$T^{\mu\nu}\f$/\f$\epsilon\f$-derived
 * quantities when running coupling is enabled (\f$\alpha_s\f$ runs
 * with either the local \f$Q_s\f$ at this cell or one of the
 * event-averaged \f$Q_s\f$ choices, per
 * `param->coupling.runWithLocalQs`/`param->coupling.runWithQs`).
 * \param[in] lat Lattice to read \f$g^2\mu_A^2\f$/\f$g^2\mu_B^2\f$
 * from (only used when `param->coupling.runWithLocalQs`).
 * \param[in] param Simulation parameters.
 * \param[in] pos Flat cell index.
 * \param[in] N Lattice side length.
 * \param[in] a Lattice spacing [fm].
 * \param[in] g Coupling \f$g\f$.
 * \param[in] c Running-coupling shape parameter.
 * \param[in] muZero \f$\mu_0\f$ in the running-coupling formula.
 * \return The local-coupling factor; `1` if running coupling is
 * disabled.
 * \see computeRunningCouplingGfactorFromScale(), which this calls
 * with either the local or an event-averaged \f$Q_s\f$ as the scale.
 */
double computeRunningCouplingGfactor(
    Lattice *lat, Parameters *param, int pos, int N, double a, double g,
    double c, double muZero);

/**
 * The event-averaged saturation scale that sets the running coupling.
 * \param[in] param Simulation parameters; `param->coupling.runWithQs`
 * selects the average over the overlap region of the smaller (`0`), the
 * mean (`1`) or the larger (`2`) of the two nuclei's \f$Q_s\f$.
 * \return That \f$\langle Q_s\rangle\f$ [GeV].
 */
double eventAverageQs(const Parameters *param);

/**
 * The running-coupling factor \f$g^2/(4\pi\alpha_s)\f$ of the event,
 * with \f$\alpha_s\f$ at `param->coupling.runningCouplingQsFactor`
 * times eventAverageQs().
 * \param[in] param Simulation parameters.
 * \return The factor; `1` if running coupling is disabled.
 */
double eventRunningCouplingGfactor(const Parameters *param);

/**
 * The factor computeRunningCouplingGfactor() gives in each lattice
 * cell, computed once for all cells: one value, unless \f$\alpha_s\f$
 * runs with the local \f$Q_s\f$.
 */
struct CouplingFactor {
    /// The factor in every cell, used when \ref perCell is empty.
    double uniform = 1.;
    /// The factor of each cell, indexed by its position, with
    /// `runningCoupling 1` and `runWithLocalQs 1`; empty otherwise.
    std::vector<double> perCell;
};

/**
 * Computes the running-coupling factor of every lattice cell, the one
 * the gluon spectrum and the eccentricity weights use.
 * \param[in] lat Lattice to read \f$g^2\mu_A^2\f$/\f$g^2\mu_B^2\f$
 * from (only with `param->coupling.runWithLocalQs`).
 * \param[in] param Simulation parameters.
 * \return The factor per cell with `runningCoupling 1` and
 * `runWithLocalQs 1`, otherwise the uniform eventRunningCouplingGfactor().
 */
CouplingFactor latticeCouplingFactor(Lattice *lat, Parameters *param);

#endif  // SRC_RUNNINGCOUPLING_H_

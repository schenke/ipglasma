#ifndef SRC_PHYSCONST_H_
#define SRC_PHYSCONST_H_

/**
 * Physical constants shared across the codebase.
 */
namespace PhysConst {
/// Generic small-number regularizer used to guard against division by
/// zero or a singular denominator [dimensionless].
const double smallEps = 1e-16;
/// Number of colors; fixed at 3 (see SU3.h/GaugeFix.cpp's comments on
/// why this codebase doesn't generalize to other \f$N_c\f$).
const int Nc = 3;
/// SU(3) adjoint dimension, \f$N_c^2-1=8\f$.
const int Nc2m1 = Nc * Nc - 1;
/// \f$\hbar c\f$, used to convert lattice/natural-unit energy densities
/// (\f$1/\mathrm{fm}^4\f$) to \f$\mathrm{GeV/fm}^3\f$ [GeV fm].
const double hbarc = 0.1973269718;
/// \f$1/\hbar c\f$ [1/(GeV fm)], used by JIMWLK::getMassRegulator()/
/// getAlphas() to make a GeV-times-fm product dimensionless.
const double invHbarc = 1.0 / hbarc;
/// Charged-pion mass [GeV].
const double m_pion = 0.13957;
/// Charged-kaon mass [GeV].
const double m_kaon = 0.493667;
/// Proton mass [GeV].
const double m_proton = 0.938272;
}  // namespace PhysConst

#endif  // SRC_PHYSCONST_H_

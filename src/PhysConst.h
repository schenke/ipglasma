#ifndef SRC_PHYSCONST_H_
#define SRC_PHYSCONST_H_

/**
 * Physical constants shared across the codebase.
 */
namespace PhysConst {
/// Generic small-number regularizer used to guard against division by
/// zero or a singular denominator [dimensionless].
const double smallEps = 1e-16;
/// \f$\hbar c\f$, used to convert lattice/natural-unit energy densities
/// (\f$1/\mathrm{fm}^4\f$) to \f$\mathrm{GeV/fm}^3\f$ [GeV fm].
const double hbarc = 0.1973269718;
/// Charged-pion mass [GeV].
const double m_pion = 0.13957;
/// Charged-kaon mass [GeV].
const double m_kaon = 0.493667;
/// Proton mass [GeV].
const double m_proton = 0.938272;
}  // namespace PhysConst

#endif  // SRC_PHYSCONST_H_

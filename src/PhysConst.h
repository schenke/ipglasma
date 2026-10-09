// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#ifndef SRC_PHYSCONST_H_
#define SRC_PHYSCONST_H_

#include <cmath>

/**
 * Physical constants shared across the codebase.
 */
namespace PhysConst {
/// Generic small-number regularizer used to guard against division by
/// zero or a singular denominator [dimensionless].
const double smallEps = 1e-16;
/// Tolerance within which an input value counts as one of the special
/// values that select a code path, e.g. `omega 1` or `jimwlkAlphaS 0`
/// (see isClose()) [dimensionless].
const double inputTolerance = 1e-8;
/**
 * Whether an input value equals a special value that selects a code
 * path, within inputTolerance.
 * \param[in] value The input value.
 * \param[in] target The special value.
 * \return Whether \f$|\text{value} - \text{target}| <\f$ inputTolerance.
 */
inline bool isClose(double value, double target) {
    return std::abs(value - target) < inputTolerance;
}
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
/// Millibarn-to-fm\f$^2\f$ conversion factor (\f$1\text{ mb} =
/// 0.1\text{ fm}^2\f$), used to convert the input `sigmaNN` (in mb)
/// to the inelastic nucleon-nucleon cross section in fm\f$^2\f$.
const double mbToFm2 = 0.1;
}  // namespace PhysConst

#endif  // SRC_PHYSCONST_H_

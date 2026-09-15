// Fragmentation.h is part of the IP-Glasma solver.
// Copyright (C) 2013 Bjoern Schenke.

#ifndef SRC_FRAGMENTATION_H_
#define SRC_FRAGMENTATION_H_

/**
 * Kniehl-Kramer-Potter (KKP) parton-to-hadron fragmentation functions
 * (B.A. Kniehl, G. Kramer, B. Potter, Nucl. Phys. B582 (2000) 514),
 * used to convert the gluon spectrum produced by the classical
 * Yang-Mills evolution into a hadron spectrum.
 *
 * \note Once a Doxygen bibliography is set up (see CONTRIBUTING.md),
 * this should switch to \c \\iref/\c \\cite instead of a plain-text
 * citation.
 */
namespace Fragmentation {

/**
 * Gluon-to-hadron KKP fragmentation function \f$D_{g\to h}(x, Q)\f$.
 *
 * The classical Yang-Mills fields IP-Glasma evolves are pure gluon
 * fields, so only the gluon row of the full KKP parton-to-hadron
 * matrix is ever physically relevant here; the quark/antiquark
 * channels are computed internally (needed for \p ih's LO/NLO fits and
 * for summing species into \c ih=7) but never returned.
 * \param[in] ih Hadron species selector:
 * - `1`: \f$(\pi^+ + \pi^-)/2\f$
 * - `2`: \f$(K^+ + K^-)/2\f$
 * - `3`: \f$(K^0 + \bar{K}^0)/2\f$
 * - `4`: \f$(p + \bar{p})/2\f$
 * - `5`: \f$\pi^0\f$
 * - `6`: \f$(n + \bar{n})/2\f$
 * - anything else: \f$(h^+ + h^-)\f$, the sum over pions, kaons and
 *   protons
 * \param[in] iset Fit order: `0` for LO, `1` for NLO. Any other value
 * exits with an error.
 * \param[in] x Longitudinal momentum fraction carried by the hadron
 * [dimensionless], must be in \f$(0, 1)\f$ for a physically sensible
 * result.
 * \param[in] qs Fragmentation scale [GeV]; internally clamped to at
 * least \f$Q_0 = \sqrt{2}\f$ GeV, the fit's validity floor.
 * \return \f$D_{g\to h}(x, Q)\f$ [dimensionless], the number density of
 * hadron \p ih per unit \p x from a fragmenting gluon. Not guaranteed
 * non-negative everywhere for the NLO fit, a known characteristic of
 * this class of NLO fit near the edges of its range, not a bug.
 */
double kkp(int ih, int iset, double x, double qs);
}  // namespace Fragmentation

#endif  // SRC_FRAGMENTATION_H_

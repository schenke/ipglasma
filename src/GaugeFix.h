// GaugeFix.h is part of the CYM solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_GAUGEFIX_H_
#define SRC_GAUGEFIX_H_

#include "FFT.h"
#include "Group.h"
#include "Lattice.h"
#include "Parameters.h"
#include "PrettyOstream.h"

/**
 * Fixes the lattice gauge fields to transverse (Coulomb) gauge.
 *
 * Used by Evolution::multiplicity() before extracting the gluon
 * multiplicity spectrum, since a physical gluon number density is only
 * well-defined once the classical Yang-Mills fields are expressed in a
 * fixed gauge.
 */
class GaugeFix {
  private:
    /// Log sink for progress/warning messages.
    PrettyOstream messager_;

  public:
    /**
     * Constructs a GaugeFix with no state of its own.
     */
    GaugeFix() {};

    /**
     * Destroys this GaugeFix (nothing to release).
     */
    ~GaugeFix() {};

    /**
     * Iteratively relaxes the lattice gauge fields towards transverse
     * (Coulomb) gauge, \f$\partial_i A_i = 0\f$, in place.
     *
     * Each iteration: projects the lattice gauge divergence at every
     * site onto a Hermitian traceless SU(3) algebra element
     * \f$\chi(x)\f$; solves the discretized Poisson equation for
     * \f$\chi\f$ in momentum space (via \p fft, dividing by the
     * lattice momentum-squared kernel); exponentiates the result
     * directly in SU(3) to get a local gauge transformation
     * \f$g(x) = \exp(i\chi(x))\f$; and applies \f$g(x)\f$ to every
     * link and adjoint field on the lattice (\c Ux, \c Uy, \c U, \c
     * U2, \c Uy2, \c Ux2). Stops early once the mean residual
     * \f$\langle\mathrm{Tr}(\chi^2)\rangle/3\f$ drops below
     * \f$10^{-9}\f$, or once it stops improving while already below
     * \f$10^{-6}\f$.
     * \param[in] fft FFT instance used to transform \f$\chi\f$ to and
     * from momentum space; must be constructed for the same lattice
     * size as \p lat.
     * \param[in,out] lat Lattice whose gauge fields are relaxed in
     * place.
     * \param[in] group Unused; kept for interface consistency with
     * other lattice-observable methods that do need a Group instance.
     * \param[in] param Simulation parameters; only `getSize()` is
     * used.
     * \param[in] steps Maximum number of relaxation iterations to run
     * before giving up regardless of the residual.
     */
    void fftChi(
        FFT *fft, Lattice *lat, Group *group, Parameters *param, int steps);
};

#endif  // SRC_GAUGEFIX_H_

// MyEigen.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_MYEIGEN_H_
#define SRC_MYEIGEN_H_

#include "Lattice.h"
#include "Parameters.h"
#include "PrettyOstream.h"

/**
 * Extracts the local fluid four-velocity \f$u^\mu\f$, energy density
 * \f$\epsilon\f$, and shear-stress tensor \f$\pi^{\mu\nu}\f$ from the
 * classical energy-momentum tensor \f$T^{\mu\nu}\f$ (already computed
 * by Evolution::tmunu()) via Landau matching, and writes the resulting
 * hydrodynamic initial condition (and/or raw \f$T^{\mu\nu}\f$) to disk.
 */
class MyEigen {
  private:
    /// Log sink for progress/warning messages.
    PrettyOstream messager_;

  public:
    /**
     * Constructs a MyEigen with no state of its own.
     */
    MyEigen() {};

    /**
     * Destroys this MyEigen (nothing to release).
     */
    ~MyEigen() {};
    /**
     * Solves for the local flow velocity/energy density/shear tensor
     * at every site (see flowVelocity4DImpl()) and writes every output
     * `param->getWriteOutputs()` enables.
     * \param[in] lat Lattice to read \f$T^{\mu\nu}\f$ from and write
     * \f$u^\mu\f$/\f$\epsilon\f$/\f$\pi^{\mu\nu}\f$ into.
     * \param[in] param Simulation parameters.
     * \param[in] it Current evolution time step.
     * \param[in] finalFlag Whether this is the final-time hydro
     * measurement (affects the hydro-text output's file name only).
     */
    void flowVelocity4D(
        Lattice *lat, Parameters *param, int it, bool finalFlag);
    /**
     * Writes only the raw \f$T^{\mu\nu}\f$ output, skipping the
     * flow-velocity solve and the hydro-text/Jazma writers (see
     * flowVelocity4DImpl()'s \c tmunuOnly path).
     * \param[in] lat Lattice to read \f$T^{\mu\nu}\f$ from.
     * \param[in] param Simulation parameters.
     * \param[in] it Current evolution time step.
     */
    void writeTmunu4D(Lattice *lat, Parameters *param, int it);

  private:
    /**
     * Shared implementation of flowVelocity4D()/writeTmunu4D(): solves
     * for \f$u^\mu\f$/\f$\epsilon\f$/\f$\pi^{\mu\nu}\f$ at every site
     * (in parallel; skipped entirely if \p tmunuOnly), then writes
     * whichever of the hydro-text, raw-\f$T^{\mu\nu}\f$, and Jazma
     * outputs `param->getWriteOutputs()` enables (only the raw
     * \f$T^{\mu\nu}\f$ writer runs if \p tmunuOnly).
     * \param[in] lat Lattice to read from and write into.
     * \param[in] param Simulation parameters.
     * \param[in] it Current evolution time step.
     * \param[in] finalFlag Whether this is the final-time hydro
     * measurement (affects the hydro-text output's file name only).
     * \param[in] tmunuOnly If `true`, skip the flow-velocity solve and
     * only write the raw \f$T^{\mu\nu}\f$ output.
     */
    void flowVelocity4DImpl(
        Lattice *lat, Parameters *param, int it, bool finalFlag,
        bool tmunuOnly);

    /**
     * Writes the hydro-flow text output (`epsilon-u-Hydro*.dat`):
     * per-cell interpolated energy density, \f$u^\mu\f$, and
     * \f$\pi^{\mu\nu}\f$ on the output grid.
     * \param[in] lat Lattice to read from.
     * \param[in] param Simulation parameters.
     * \param[in] it Current evolution time step.
     * \param[in] finalFlag Selects the output file name
     * (`...-TauHydro-<id>.dat` if `true`, `...-t<time>-<id>.dat`
     * otherwise).
     * \param[in] tmunuOnly Forwarded from flowVelocity4DImpl(); this
     * writer always no-ops if `true` (see \p N below).
     * \param[in] N Lattice side length.
     * \param[in] L Physical lattice size [fm].
     * \param[in] a Lattice spacing [fm].
     * \param[in] dtau Evolution time step [lattice units].
     * \param[in] gfactor Running-coupling rescaling factor (`1` if
     * running coupling is disabled).
     * \param[in] hx Output grid size in \f$x\f$.
     * \param[in] hy Output grid size in \f$y\f$.
     * \param[in] heta Output grid size in \f$\eta\f$.
     * \param[in] hL Physical output grid size [fm].
     * \param[in] deta Output grid spacing in \f$\eta\f$.
     * \param[in] ha Output grid spacing in \f$x\f$/\f$y\f$ [fm].
     * \param[in] tau0 Current proper time [fm/c].
     * \return The total energy \f$E_{\text{tot}}\f$ integrated over the
     * grid, computed whenever either this text output or writeJazma()'s
     * output is enabled (`param->getWriteOutputs()` bit 0 or bit 1) --
     * writeJazma() needs \f$E_{\text{tot}}\f$ even when this writer's
     * own text file is off (e.g. `writeOutputs=2`, "standalone Jazma
     * output" per README.md).
     */
    double writeHydroText(
        Lattice *lat, Parameters *param, int it, bool finalFlag, bool tmunuOnly,
        int N, double L, double a, double dtau, double gfactor, int hx, int hy,
        int heta, double hL, double deta, double ha, double tau0);

    /**
     * Writes the raw/binary Tmunu output (`Tmunu-t*.dat` or `*.ipgt`):
     * per-cell interpolated \f$T^{\mu\nu}\f$ components on the output
     * grid, in text or a little-endian binary format depending on
     * `param->getWriteTmunuBinary()`/`IPGLASMA_BINARY_TMUNU`. A no-op
     * unless `param->getWriteOutputs()`'s bit 2 (value `4`) is set.
     * \param[in] lat Lattice to read from.
     * \param[in] param Simulation parameters.
     * \param[in] it Current evolution time step.
     * \param[in] N Lattice side length.
     * \param[in] L Physical lattice size [fm].
     * \param[in] a Lattice spacing [fm].
     * \param[in] dtau Evolution time step [lattice units].
     * \param[in] gfactor Running-coupling rescaling factor.
     * \param[in] hx Output grid size in \f$x\f$.
     * \param[in] hy Output grid size in \f$y\f$.
     * \param[in] heta Output grid size in \f$\eta\f$.
     * \param[in] hL Physical output grid size [fm].
     * \param[in] deta Output grid spacing in \f$\eta\f$.
     * \param[in] ha Output grid spacing in \f$x\f$/\f$y\f$ [fm].
     * \param[in] tau0 Current proper time [fm/c].
     */
    void writeRawTmunu(
        Lattice *lat, Parameters *param, int it, int N, double L, double a,
        double dtau, double gfactor, int hx, int hy, int heta, double hL,
        double deta, double ha, double tau0);

    /**
     * Writes the Jazma output (`Jazma-Hydro-t*.dat`): per-cell
     * \f$g^2\mu_A^2 g^2\mu_B^2\f$, normalized so its grid integral
     * matches \p Etot. A no-op unless `param->getWriteOutputs()`'s bit
     * 1 (value `2`) is set.
     * \param[in] lat Lattice to read from.
     * \param[in] param Simulation parameters.
     * \param[in] it Current evolution time step.
     * \param[in] Etot Total energy to normalize the output to (from
     * writeHydroText()).
     * \param[in] N Lattice side length.
     * \param[in] L Physical lattice size [fm].
     * \param[in] a Lattice spacing [fm].
     * \param[in] dtau Evolution time step [lattice units].
     * \param[in] hx Output grid size in \f$x\f$.
     * \param[in] hy Output grid size in \f$y\f$.
     * \param[in] heta Output grid size in \f$\eta\f$.
     * \param[in] hL Physical output grid size [fm].
     * \param[in] deta Output grid spacing in \f$\eta\f$.
     * \param[in] ha Output grid spacing in \f$x\f$/\f$y\f$ [fm].
     */
    void writeJazma(
        Lattice *lat, Parameters *param, int it, double Etot, int N, double L,
        double a, double dtau, int hx, int hy, int heta, double hL, double deta,
        double ha);
};

#endif  // SRC_MYEIGEN_H_

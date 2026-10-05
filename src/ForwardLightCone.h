// ForwardLightCone.h is part of the IP-Glasma solver.

#ifndef SRC_FORWARDLIGHTCONE_H_
#define SRC_FORWARDLIGHTCONE_H_

#include <complex>
#include <cstdint>
#include <vector>

#include "Group.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "PrettyOstream.h"

/**
 * The initial condition in the forward light cone: from the projectile's
 * and target's Wilson lines (lat->U / lat->U2, after
 * Init::shiftFieldsWithImpactParameter()), solves the matching
 * condition for the gauge links and computes the electric field and
 * \f$\pi\f$ that start the classical Yang-Mills evolution
 * (initialize()). The per-link solve is findU().
 */
class ForwardLightCone {
  public:
    /**
     * \param[in] group SU(3) generators used by the matching condition.
     */
    explicit ForwardLightCone(Group *group) : group_(group), one_(1.) {}

    /**
     * Matches the projectile's and target's Wilson lines across the
     * forward light cone to obtain the post-collision gauge links (\c
     * Ux/\c Uy), electric field (`U`/`U2`, reused as scratch and then
     * as the actual electric-field components), and momentum
     * \f$\pi\f$ (\c Ux2) that seed the subsequent classical Yang-Mills
     * evolution. Runs its steps (see below) in sequence inside one
     * shared `#pragma omp parallel` region; each step has its own
     * `#pragma omp for`, so the implicit barrier at the end of each
     * keeps them correctly ordered.
     * \param[in,out] lat Lattice whose fields are transformed in place.
     * \param[in] param Simulation parameters.
     */
    void initialize(Lattice *lat, Parameters *param);
    /**
     * Replaces any NaN `U`/`U2` (left over from a failed
     * forward-lightcone solve at a previous stage) with the identity.
     * \param[in,out] lat Lattice to sanitize.
     * \param[in] N2 Total number of lattice sites.
     */
    void sanitizeU(Lattice *lat, int N2);
    /**
     * Scratch matrices for computeLinksTeam(), reused
     * across cells to avoid reallocating.
     */
    struct LinkScratch {
        /// Conjugate-transposed neighbor link in \f$x\f$.
        Matrix UDx;
        /// Conjugate-transposed neighbor link in \f$y\f$.
        Matrix UDy;
    };
    /**
     * Computes the pre-collision forward/backward gauge links \c
     * Ux1/\c Uy1 (from \c U) and `Ux2`/`Uy2` (from \c U2) needed by
     * the forward-lightcone matching.
     * \param[in,out] lat Lattice to read `U`/`U2` from and write \c
     * Ux1/`Uy1`/`Ux2`/`Uy2` into.
     * \param[in] N2 Total number of lattice sites.
     * \param[in,out] scratch Thread-local scratch storage.
     */
    void computeLinksTeam(Lattice *lat, int N2, LinkScratch &scratch);
    /**
     * Scratch matrices for computeUxUyTeam(), reused
     * across cells to avoid reallocating.
     */
    struct UScratch {
        /// Newton-solved matching result for this cell/direction.
        Matrix temp2;
        /// Projectile-side link in \f$x\f$, passed to
        /// findU().
        Matrix UDx1;
        /// Target-side link in \f$x\f$, passed to
        /// findU().
        Matrix UDx2;
        /// Projectile-side link in \f$y\f$, passed to
        /// findU().
        Matrix UDy1;
        /// Target-side link in \f$y\f$, passed to
        /// findU().
        Matrix UDy2;
    };
    /**
     * Solves for the post-collision gauge links `Ux`/`Uy` from \c
     * Ux1/\c Ux2 and `Uy1`/`Uy2` via findU(),
     * logging a warning at each cell where the Newton solve didn't
     * converge.
     * \param[in,out] lat Lattice to read `Ux1`/`Uy1`/`Ux2`/`Uy2`
     * from and write `Ux`/`Uy` into.
     * \param[in] param Simulation parameters; only used for the
     * warning message's cell coordinates and `param->run.randomSeed`/
     * `param->event.eventId` (to seed findU()'s deterministic
     * retry stream).
     * \param[in] N2 Total number of lattice sites.
     * \param[in,out] scratch Thread-local scratch storage.
     */
    void computeUxUyTeam(
        Lattice *lat, Parameters *param, int N2, UScratch &scratch);
    /**
     * Scratch matrices for computeElectricFieldTeam(),
     * reused across cells to avoid reallocating.
     */
    struct ElectricFieldScratch {
        /// General scratch accumulating the field contribution.
        Matrix temp2;
        /// `lat->Ux1[pos] - lat->Ux2[pos]` (or the same at a neighbor
        /// cell): the difference between the pre-collision forward/
        /// backward links in \f$x\f$.
        Matrix Ux1mUx2;
        /// Conjugate-transposed `lat->Ux1` at the current cell.
        Matrix UDx1;
        /// Conjugate-transposed `lat->Ux2` at the current cell.
        Matrix UDx2;
        /// \c UDx1 minus \c UDx2.
        Matrix UDx1mUDx2;
        /// The post-collision link \c Ux at the current cell.
        Matrix Ux;
        /// Conjugate-transposed \c Ux.
        Matrix UDx;
        /// `lat->Uy1[pos] - lat->Uy2[pos]` (or the same at a neighbor
        /// cell): the difference between the pre-collision forward/
        /// backward links in \f$y\f$.
        Matrix Uy1mUy2;
        /// Conjugate-transposed `lat->Uy1` at the current cell.
        Matrix UDy1;
        /// Conjugate-transposed `lat->Uy2` at the current cell.
        Matrix UDy2;
        /// \c UDy1 minus \c UDy2.
        Matrix UDy1mUDy2;
        /// The post-collision link \c Uy at the current cell.
        Matrix Uy;
        /// Conjugate-transposed \c Uy.
        Matrix UDy;
    };
    /**
     * Computes the initial electric field's contribution from one
     * direction (`neighborX`/`neighborY` select minus-shifted or
     * plus-shifted neighbors), written into \p outputField. Called
     * once with `(posmX, posmY, lat->U)` and once with `(pospX, pospY,
     * lat->U2)` -- previously two copy-pasted loops.
     * \param[in] lat Lattice to read `Ux1`/`Uy1`/`Ux2`/`Uy2`/\c
     * Ux/\c Uy from.
     * \param[in] N2 Total number of lattice sites.
     * \param[in] neighborX Neighbor-index table to use in \f$x\f$
     * (`lat->posmX` or `lat->pospX`).
     * \param[in] neighborY Neighbor-index table to use in \f$y\f$
     * (`lat->posmY` or `lat->pospY`).
     * \param[in,out] outputField Field this contribution is written
     * into (`lat->U` or `lat->U2`, reused as scratch here ahead of
     * resetFieldsTeam() zeroing them).
     * \param[in,out] scratch Thread-local scratch storage.
     */
    void computeElectricFieldTeam(
        Lattice *lat, int N2, const std::vector<int> &neighborX,
        const std::vector<int> &neighborY, std::vector<Matrix> &outputField,
        ElectricFieldScratch &scratch);
    /**
     * Scratch matrices for computePlaquetteTeam(),
     * reused across cells to avoid reallocating.
     */
    struct PlaquetteScratch {
        /// Conjugate-transposed neighbor link in \f$x\f$.
        Matrix UDx;
        /// Conjugate-transposed neighbor link in \f$y\f$.
        Matrix UDy;
        /// The computed spatial plaquette.
        Matrix Uplaq;
    };
    /**
     * Computes the spatial plaquette from `Ux`/`Uy` into \c
     * lat->Uy1 (reused as scratch here, ahead of
     * computePiTeam()/resetFieldsTeam()
     * repurposing it further).
     * \param[in,out] lat Lattice to read `Ux`/`Uy` from and write \c
     * Uy1 into.
     * \param[in] N2 Total number of lattice sites.
     * \param[in,out] scratch Thread-local scratch storage.
     */
    void computePlaquetteTeam(Lattice *lat, int N2, PlaquetteScratch &scratch);
    /**
     * Sets \c lat->Ux2 to the initial momentum \f$\pi\f$ (\f$E^\eta\f$)
     * in lattice units, from \c lat->U (the electric field just
     * computed by computeElectricFieldTeam()).
     * \param[in,out] lat Lattice to read \c U from and write \c Ux2
     * into.
     * \param[in] param Simulation parameters; `param->coupling.g` is used.
     * \param[in] N2 Total number of lattice sites.
     */
    void computePiTeam(Lattice *lat, Parameters *param, int N2);
    /**
     * Zeroes \c lat->U/`U2`/`Uy2` and resets \c lat->Ux1 to the
     * identity, now that this event's forward-lightcone fields have
     * been consumed by the steps above.
     * \param[in,out] lat Lattice to reset.
     * \param[in] N2 Total number of lattice sites.
     */
    void resetFieldsTeam(Lattice *lat, int N2);
    /**
     * Solves the dense linear system \f$J x = F\f$ via GSL LU
     * decomposition, sized for the SU(3) adjoint dimension
     * (\f$8\times8\f$); used by findU()'s Newton
     * iteration.
     * \param[in] Jab Row-major \f$8\times8\f$ Jacobian.
     * \param[in] Fa Length-8 right-hand side.
     * \param[out] xvec Filled with the length-8 solution.
     */
    void solveAxb(double *Jab, double *Fa, std::vector<double> &xvec);
    /**
     * Solves the forward-lightcone matching condition
     * \f$U_1+U_2 = U_{\text{sol}} U_1 U_2 + (U_1 U_2)^\dagger
     * U_{\text{sol}}^\dagger\f$ for \f$U_{\text{sol}}\f$ via Newton's
     * method in \f$U_{\text{sol}}\f$'s eight SU(3) generator
     * coefficients (computeResidual()/
     * computeJacobian()/solveAxb()). If the iteration
     * diverges or fails to converge within 2000 iterations, restarts
     * from a fresh initial guess drawn from a deterministic,
     * seed-independent pseudo-random stream (so restarts are
     * reproducible without touching the shared Random state), up to
     * 200 restarts; returns the best (lowest-residual) estimate found
     * across all restarts either way.
     * \param[in] U1 Projectile-side link.
     * \param[in] U2 Target-side link.
     * \param[out] Usol The solved (or best-effort) matching matrix.
     * \param[in] retrySeed Seed for the deterministic restart stream
     * (see \c forwardLightconeRetrySeed() in ForwardLightCone.cpp), unique per
     * cell/direction/event/run so restarts are reproducible but
     * uncorrelated across cells.
     * \return `true` if the iteration converged (residual below
     * \f$10^{-6}\f$); `false` if it exhausted all restarts without
     * converging (a warning is logged in that case).
     */
    bool findU(Matrix &U1, Matrix &U2, Matrix &Usol, std::uint64_t retrySeed);
    /**
     * findU()'s per-iteration residual: computes
     * \f$F_a\f$ (the quantity the Newton iteration drives to zero) into
     * \p Fa.
     * \param[in] U1pU2 \f$U_1+U_2\f$.
     * \param[in] U1pU2dagger \f$(U_1+U_2)^\dagger\f$.
     * \param[in] Usol Current estimate of the matching matrix.
     * \param[in] Usoldagger \f$U_{\text{sol}}^\dagger\f$.
     * \param[in] traceCache Per-generator trace terms that don't depend
     * on \p Usol, precomputed once by findU().
     * \param[out] Fa Filled with the eight residual components.
     * \return \f$F_{\text{zero}} = \sum_a |F_a|\f$, the convergence
     * criterion.
     */
    double computeResidual(
        const Matrix &U1pU2, const Matrix &U1pU2dagger, const Matrix &Usol,
        const Matrix &Usoldagger,
        const std::vector<complex<double>> &traceCache, double *Fa);
    /**
     * Computes the Jacobian \f$dF/d\alpha\f$ into \p Jab, via a
     * numerical finite difference, falling back to an analytical
     * approximation if the numerical one comes out singular. \p alpha
     * is left unchanged on return (each finite-difference step adds
     * then subtracts its own \f$d\alpha_{b_i}\f$).
     * \param[in] U0 \f$U_1 U_2\f$.
     * \param[in] U1pU2 \f$U_1+U_2\f$.
     * \param[in] Usoldagger \f$U_{\text{sol}}^\dagger\f$ at the current
     * estimate.
     * \param[in,out] MtempArr Scratch, reused across calls; filled with
     * \f$t_a (U_1+U_2)\f$ for the analytical fallback.
     * \param[in,out] alpha Current Newton estimate; perturbed and
     * restored in place during the finite-difference computation.
     * \param[out] Jab Filled with the row-major \f$8\times8\f$
     * Jacobian.
     */
    void computeJacobian(
        const Matrix &U0, const Matrix &U1pU2, const Matrix &Usoldagger,
        std::vector<Matrix> &MtempArr, std::vector<double> &alpha, double *Jab);

  private:
    /// SU(3) generators.
    Group *group_;
    /// The identity matrix.
    Matrix one_;
    /// Log sink for progress/warning messages.
    PrettyOstream messager_;
};

#endif  // SRC_FORWARDLIGHTCONE_H_

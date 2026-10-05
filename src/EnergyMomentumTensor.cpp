// EnergyMomentumTensor.cpp is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "EnergyMomentumTensor.h"

#include <algorithm>
#include <cmath>

#include "Instrumentation.h"
#include "Matrix.h"
#include "SU3.h"

namespace {

/**
 * Computes the traceless part of \p lhs minus \p rhs:
 * \f$(\text{lhs}-\text{rhs}) - \text{tr}(\text{lhs}-\text{rhs})/3\f$.
 * Used throughout tmunuOffDiagonalTeam() to turn a pair of four-link
 * chains into the traceless combination its off-diagonal
 * \f$T^{\mu\nu}\f$ components need.
 * \param[in] lhs Minuend matrix.
 * \param[in] rhs Subtrahend matrix.
 * \param[in] one Reusable identity matrix.
 * \return The traceless difference.
 */
inline Matrix makeTmunuTracelessDifference(
    const Matrix &lhs, const Matrix &rhs, const Matrix &one) {
    Matrix out = lhs - rhs;
    out -= (out.trace() / 3.0) * one;
    return out;
}

/// Scratch matrices for tmunuPlaquetteTeam(), reused across cells to
/// avoid reallocating.
struct TmunuPlaquetteScratch {
    /// Conjugate-transposed \f$U_x\f$ at the plaquette's far corner.
    Matrix UDx;
    /// Conjugate-transposed \f$U_y\f$ at this cell.
    Matrix UDy;
    /// The computed spatial plaquette.
    Matrix Uplaq;
};

/**
 * Precomputes the spatial plaquette \f$U_x(x)\,U_y(x+\hat x)\,
 * U_x(x+\hat y)^\dagger\,U_y(x)^\dagger\f$ at every cell into
 * \c lat->Uy1, consumed by tmunuDiagonalMagneticTeam() below. The
 * outermost ring is a nonphysical guard region (\f$T^{\mu\nu}\f$'s
 * stencils need a genuine one-cell neighborhood), so it gets the
 * identity instead of a clamped, gauge-noncovariant plaquette.
 * \param[in,out] lat Lattice to read `Ux`/`Uy` from and write \c Uy1
 * into.
 * \param[in] N Lattice side length.
 * \param[in] one Reusable identity matrix.
 * \param[in,out] scratch Thread-local scratch storage.
 */
void tmunuPlaquetteTeam(
    Lattice *lat, int N, const Matrix &one, TmunuPlaquetteScratch &scratch) {
    int pos, posX, posY;
#pragma omp for
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = lat->positionFromXY(i, j);
            if (i == 0 || j == 0 || i == N - 1 || j == N - 1) {
                lat->Uy1[pos] = one;
                continue;
            }

            posX = lat->pospX[pos];
            posY = lat->pospY[pos];

            scratch.UDx = lat->Ux[posY];
            scratch.UDy = lat->Uy[pos];
            scratch.UDx.conjg();
            scratch.UDy.conjg();

            scratch.Uplaq =
                lat->Ux[pos] * (lat->Uy[posX] * (scratch.UDx * scratch.UDy));
            lat->Uy1[pos] = (scratch.Uplaq);
        }
    }
}

/// Scratch matrices for tmunuDiagonalElectricTeam(), reused across
/// cells to avoid reallocating.
struct TmunuDiagonalElectricScratch {
    /// This cell's electric field \f$E_1\f$ (\c U).
    Matrix E1;
    /// This cell's electric field \f$E_2\f$ (\c U2).
    Matrix E2;
    /// \f$E_1\f$ at the \f$+\hat y\f$ neighbor.
    Matrix E1p;
    /// \f$E_2\f$ at the \f$+\hat x\f$ neighbor.
    Matrix E2p;
    /// This cell's \f$\pi\f$.
    Matrix pi;
    /// \f$\pi\f$ at the \f$+\hat x\f$ neighbor.
    Matrix piX;
    /// \f$\pi\f$ at the \f$+\hat y\f$ neighbor.
    Matrix piY;
    /// \f$\pi\f$ at the \f$(+\hat x,+\hat y)\f$ neighbor.
    Matrix piXY;
};

/**
 * Sets \f$T^{\tau\tau}\f$, \f$T^{xx}\f$, \f$T^{yy}\f$,
 * \f$T^{\eta\eta}\f$'s electric (\f$E\f$, \f$\pi\f$) contribution.
 * Sets each field outright (rather than adding to it) since this runs
 * before tmunuDiagonalMagneticTeam(), which adds the magnetic/gradient
 * contribution on top; zeroes all four at the nonphysical boundary
 * ring.
 * \param[in,out] lat Lattice to read `U`/`U2`/`Ux2` from and write
 * \f$T^{\tau\tau}\f$/\f$T^{xx}\f$/\f$T^{yy}\f$/\f$T^{\eta\eta}\f$ into
 * (via `lat->cells`).
 * \param[in] N Lattice side length.
 * \param[in] it Current time step index, used to convert to physical
 * units.
 * \param[in] dtau Time step [lattice units].
 * \param[in] g Coupling \f$g\f$.
 * \param[in,out] scratch Thread-local scratch storage.
 */
void tmunuDiagonalElectricTeam(
    Lattice *lat, int N, int it, double dtau, double g,
    TmunuDiagonalElectricScratch &scratch) {
    int pos, posX, posY, posXY;
#pragma omp for
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = lat->positionFromXY(i, j);
            if (i == 0 || j == 0 || i == N - 1 || j == N - 1) {
                lat->cells[pos]->setTtautau(0.);
                lat->cells[pos]->setTxx(0.);
                lat->cells[pos]->setTyy(0.);
                lat->cells[pos]->setTetaeta(0.);
                continue;
            }
            posX = lat->pospX[pos];
            posY = lat->pospY[pos];

            posXY = lat->positionFromXY(
                std::min(N - 1, i + 1), std::min(N - 1, j + 1));

            scratch.E1 = lat->U[pos];
            scratch.E2 = lat->U2[pos];
            scratch.E1p = lat->U[posY];
            scratch.E2p = lat->U2[posX];  // shift y value in x direction

            scratch.pi = lat->Ux2[pos];
            scratch.piX = lat->Ux2[posX];
            scratch.piY = lat->Ux2[posY];
            scratch.piXY = lat->Ux2[posXY];

            // These observables only need traces of matrix squares.
            // Computing the complete 3x3 products here used to create 32
            // Matrix temporaries per site (the same eight traces repeated
            // for four tensor components).  Evaluate each SU(3) trace once
            // and reuse it.
            const double e1Sq = su3::traceSquare(scratch.E1).real();
            const double e1pSq = su3::traceSquare(scratch.E1p).real();
            const double e2Sq = su3::traceSquare(scratch.E2).real();
            const double e2pSq = su3::traceSquare(scratch.E2p).real();
            const double piSq = su3::traceSquare(scratch.pi).real();
            const double piXSq = su3::traceSquare(scratch.piX).real();
            const double piYSq = su3::traceSquare(scratch.piY).real();
            const double piXYSq = su3::traceSquare(scratch.piXY).real();

            const double invTau2 = 1. / (it * dtau) / (it * dtau);
            const double electricPrefactor = g * g / (it * dtau) / (it * dtau);
            const double eSum = e1Sq + e1pSq + e2Sq + e2pSq;
            const double piSum = piSq + piXSq + piYSq + piXYSq;

            lat->cells[pos]->setTtautau(
                electricPrefactor * eSum / 2. + piSum / 4.);
            lat->cells[pos]->setTxx(
                electricPrefactor * (-e1Sq - e1pSq + e2Sq + e2pSq) / 2.
                + piSum / 4.);
            lat->cells[pos]->setTyy(
                electricPrefactor * (e1Sq + e1pSq - e2Sq - e2pSq) / 2.
                + piSum / 4.);
            lat->cells[pos]->setTetaeta(
                invTau2 * (electricPrefactor * eSum / 2. - piSum / 4.));
        }
    }
}

/// Scratch matrices for tmunuDiagonalMagneticTeam(), reused across
/// cells to avoid reallocating.
struct TmunuDiagonalMagneticScratch {
    /// This cell's spatial plaquette, from `lat->Uy1`
    /// (tmunuPlaquetteTeam()'s output).
    Matrix Uplaq;
    /// This cell's \f$\phi\f$.
    Matrix phi;
    /// \f$\phi\f$ at the \f$+\hat x\f$ neighbor.
    Matrix phiX;
    /// \f$\phi\f$ at the \f$+\hat y\f$ neighbor.
    Matrix phiY;
    /// \f$\phi\f$ at the \f$(+\hat x,+\hat y)\f$ neighbor.
    Matrix phiXY;
    /// \f$U_x\f$ used to parallel-transport \f$\phi\f$ (reused for two
    /// different base cells within one iteration).
    Matrix Ux;
    /// \f$U_y\f$ used to parallel-transport \f$\phi\f$ (reused for two
    /// different base cells within one iteration).
    Matrix Uy;
    /// Conjugate-transposed \c Ux.
    Matrix UDx;
    /// Conjugate-transposed \c Uy.
    Matrix UDy;
    /// \f$\phi_X\f$ parallel-transported back to this cell,
    /// \f$U_x\phi_X U_x^\dagger\f$.
    Matrix phiTildeX;
    /// \f$\phi_Y\f$ parallel-transported back to this cell,
    /// \f$U_y\phi_Y U_y^\dagger\f$.
    Matrix phiTildeY;
    /// \f$\phi_{XY}\f$ parallel-transported back to the \f$+\hat
    /// y\f$ neighbor via \f$U_x\f$ there.
    Matrix phiTildeXY1;
    /// \f$\phi_{XY}\f$ parallel-transported back to the \f$+\hat
    /// x\f$ neighbor via \f$U_y\f$ there.
    Matrix phiTildeXY2;
};

/**
 * Adds \f$T^{\tau\tau}\f$, \f$T^{xx}\f$, \f$T^{yy}\f$,
 * \f$T^{\eta\eta}\f$'s magnetic (plaquette) and gradient (\f$\phi\f$)
 * contribution on top of whatever tmunuDiagonalElectricTeam() set
 * (`0` at the boundary, the electric part elsewhere); a no-op at the
 * nonphysical boundary ring.
 * \param[in,out] lat Lattice to read `Uy1`/`Uy2`/`Ux`/`Uy` from
 * and add into
 * \f$T^{\tau\tau}\f$/\f$T^{xx}\f$/\f$T^{yy}\f$/\f$T^{\eta\eta}\f$ (via
 * `lat->cells`).
 * \param[in] N Lattice side length.
 * \param[in] it Current time step index, used to convert to physical
 * units.
 * \param[in] dtau Time step [lattice units].
 * \param[in] g Coupling \f$g\f$.
 * \param[in,out] scratch Thread-local scratch storage.
 */
void tmunuDiagonalMagneticTeam(
    Lattice *lat, int N, int it, double dtau, double g,
    TmunuDiagonalMagneticScratch &scratch) {
    int pos, posX, posY, posXY;
#pragma omp for
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = lat->positionFromXY(i, j);
            if (i == 0 || j == 0 || i == N - 1 || j == N - 1) {
                continue;
            }

            posX = lat->pospX[pos];
            posY = lat->pospY[pos];

            posXY = lat->positionFromXY(
                std::min(N - 1, i + 1), std::min(N - 1, j + 1));

            scratch.Uplaq = lat->Uy1[pos];

            scratch.phi = lat->Uy2[pos];
            scratch.phiX = lat->Uy2[posX];
            scratch.phiY = lat->Uy2[posY];
            scratch.phiXY = lat->Uy2[posXY];

            scratch.Ux = lat->Ux[pos];
            scratch.Uy = lat->Uy[pos];
            scratch.UDx = scratch.Ux;
            scratch.UDx.conjg();
            scratch.UDy = scratch.Uy;
            scratch.UDy.conjg();

            scratch.phiTildeX = scratch.Ux * scratch.phiX * scratch.UDx;
            scratch.phiTildeY = scratch.Uy * scratch.phiY * scratch.UDy;

            // same at one up in the other direction
            scratch.Ux = lat->Ux[posY];
            scratch.Uy = lat->Uy[posX];
            scratch.UDx = scratch.Ux;
            scratch.UDx.conjg();
            scratch.UDy = scratch.Uy;
            scratch.UDy.conjg();

            scratch.phiTildeXY1 = scratch.Ux * scratch.phiXY * scratch.UDx;
            scratch.phiTildeXY2 = scratch.Uy * scratch.phiXY * scratch.UDy;

            // The four covariant-gradient square traces are likewise shared
            // by all diagonal tensor components.  Evaluate (A-B)^2 directly
            // in the trace kernel, avoiding both subtraction and product
            // Matrix temporaries.
            const double gradX0 =
                su3::traceDifferenceSquare(scratch.phi, scratch.phiTildeX)
                    .real();
            const double gradX1 =
                su3::traceDifferenceSquare(scratch.phiY, scratch.phiTildeXY1)
                    .real();
            const double gradY0 =
                su3::traceDifferenceSquare(scratch.phi, scratch.phiTildeY)
                    .real();
            const double gradY1 =
                su3::traceDifferenceSquare(scratch.phiX, scratch.phiTildeXY2)
                    .real();
            const double invTau2 = 1. / (it * dtau) / (it * dtau);
            const double gradientPrefactor = 0.5 / (it * dtau) / (it * dtau);
            const double plaquetteEnergy =
                2. / pow(g, 2.) * (3.0 - su3::trace(scratch.Uplaq).real());

            lat->cells[pos]->setTtautau(
                lat->cells[pos]->getTtautau() + plaquetteEnergy
                + gradientPrefactor * (gradX0 + gradX1 + gradY0 + gradY1));

            lat->cells[pos]->setTxx(
                lat->cells[pos]->getTxx() + plaquetteEnergy
                + gradientPrefactor * (gradX0 + gradX1 - gradY0 - gradY1));

            lat->cells[pos]->setTyy(
                lat->cells[pos]->getTyy() + plaquetteEnergy
                + gradientPrefactor * (-gradX0 - gradX1 + gradY0 + gradY1));

            lat->cells[pos]->setTetaeta(
                lat->cells[pos]->getTetaeta()
                + invTau2
                      * (-plaquetteEnergy
                         + gradientPrefactor
                               * (gradX0 + gradX1 + gradY0 + gradY1)));
        }
    }
}

/**
 * Scratch matrices for tmunuOffDiagonalTeam(), reused across cells to
 * avoid reallocating.
 *
 * Naming convention: `U`/`UD` is the transverse gauge link (`x` or
 * `y` direction) or its conjugate transpose; a trailing `pX`/`mX`/
 * `pY`/`mY` (optionally combined, e.g. `pXpY`) selects a neighbor
 * shifted \f$\pm1\f$ cell in that direction, `p2X`/`p2Y` a neighbor
 * shifted \f$+2\f$ cells; a bare name (no suffix) is this cell. `E1`/
 * `E2`/`pi`/`phi` follow the same neighbor-suffix convention for the
 * electric fields, \f$\pi\f$, and \f$\phi\f$.
 */
struct TmunuOffDiagonalScratch {
    /// This cell's \f$U_x\f$.
    Matrix Ux;
    /// This cell's \f$U_y\f$.
    Matrix Uy;
    /// \f$U_x\f$ at the \f$-\hat x\f$ neighbor.
    Matrix UxmX;
    /// \f$U_y\f$ at the \f$-\hat y\f$ neighbor.
    Matrix UymY;
    /// Conjugate-transposed \c Ux.
    Matrix UDx;
    /// Conjugate-transposed \c Uy.
    Matrix UDy;
    /// Conjugate-transposed \f$U_x\f$ at the \f$-\hat x\f$ neighbor.
    Matrix UDxmX;
    /// Conjugate-transposed \f$U_y\f$ at the \f$-\hat y\f$ neighbor.
    Matrix UDymY;
    /// Conjugate-transposed \f$U_x\f$ at the \f$(-\hat x,+\hat
    /// y)\f$ neighbor.
    Matrix UDxmXpY;
    /// Conjugate-transposed \f$U_x\f$ at the \f$(+\hat x,+\hat
    /// y)\f$ neighbor.
    Matrix UDxpXpY;
    /// \f$U_x\f$ at the \f$+\hat x\f$ neighbor.
    Matrix UxpX;
    /// \f$U_x\f$ at the \f$+\hat y\f$ neighbor.
    Matrix UxpY;
    /// Conjugate-transposed \f$U_x\f$ at the \f$+\hat y\f$ neighbor.
    Matrix UDxpY;
    /// \f$U_x\f$ at the \f$(+\hat x,+\hat y)\f$ neighbor.
    Matrix UxpXpY;
    /// Conjugate-transposed \f$U_y\f$ at the \f$(+\hat x,-\hat
    /// y)\f$ neighbor.
    Matrix UDypXmY;
    /// \f$U_y\f$ at the \f$+\hat y\f$ neighbor.
    Matrix UypY;
    /// \f$U_y\f$ at the \f$+\hat x\f$ neighbor.
    Matrix UypX;
    /// Conjugate-transposed \f$U_y\f$ at the \f$+\hat x\f$ neighbor.
    Matrix UDypX;
    /// \f$U_y\f$ at the \f$(+\hat x,+\hat y)\f$ neighbor.
    Matrix UypXpY;
    /// Conjugate-transposed \f$U_y\f$ at the \f$(+\hat x,+\hat
    /// y)\f$ neighbor.
    Matrix UDypXpY;
    /// \f$U_y\f$ at the \f$-\hat x\f$ neighbor.
    Matrix UymX;
    /// \f$U_x\f$ at the \f$(-\hat x,+\hat y)\f$ neighbor.
    Matrix UxmXpY;
    /// \f$U_x\f$ at the \f$-\hat y\f$ neighbor.
    Matrix UxmY;
    /// Conjugate-transposed \f$U_x\f$ at the \f$-\hat y\f$ neighbor.
    Matrix UDxmY;
    /// \f$U_y\f$ at the \f$(+\hat x,-\hat y)\f$ neighbor.
    Matrix UypXmY;
    /// Conjugate-transposed \f$U_y\f$ at the \f$+2\hat x\f$ neighbor.
    Matrix UDyp2X;
    /// \f$U_y\f$ at the \f$+2\hat x\f$ neighbor.
    Matrix Uyp2X;
    /// Conjugate-transposed \f$U_x\f$ at the \f$+\hat x\f$ neighbor.
    Matrix UDxpX;
    /// \f$U_x\f$ at the \f$+2\hat y\f$ neighbor.
    Matrix Uxp2Y;
    /// Conjugate-transposed \f$U_x\f$ at the \f$+2\hat y\f$ neighbor.
    Matrix UDxp2Y;
    /// Conjugate-transposed \f$U_y\f$ at the \f$+\hat y\f$ neighbor.
    Matrix UDypY;
    /// Conjugate-transposed \f$U_y\f$ at the \f$-\hat x\f$ neighbor.
    Matrix UDymX;
    /// This cell's electric field \f$E_1\f$.
    Matrix E1;
    /// This cell's electric field \f$E_2\f$.
    Matrix E2;
    /// \f$E_1\f$ at the \f$+\hat y\f$ neighbor.
    Matrix E1p;
    /// \f$E_2\f$ at the \f$+\hat x\f$ neighbor.
    Matrix E2p;
    /// This cell's \f$\pi\f$.
    Matrix pi;
    /// \f$\pi\f$ at the \f$+\hat x\f$ neighbor.
    Matrix piX;
    /// \f$\pi\f$ at the \f$+\hat y\f$ neighbor.
    Matrix piY;
    /// \f$\pi\f$ at the \f$(+\hat x,+\hat y)\f$ neighbor.
    Matrix piXY;
    /// This cell's \f$\phi\f$.
    Matrix phi;
    /// \f$\phi\f$ at the \f$+\hat x\f$ neighbor.
    Matrix phiX;
    /// \f$\phi\f$ at the \f$+\hat y\f$ neighbor.
    Matrix phiY;
    /// \f$\phi\f$ at the \f$(+\hat x,+\hat y)\f$ neighbor.
    Matrix phiXY;
    /// \f$\phi\f$ at the \f$-\hat x\f$ neighbor.
    Matrix phimX;
    /// \f$\phi\f$ at the \f$-\hat y\f$ neighbor.
    Matrix phimY;
    /// \f$\phi\f$ at the \f$(+2\hat x,+\hat y)\f$ neighbor.
    Matrix phi2XY;
    /// \f$\phi\f$ at the \f$(+\hat x,+2\hat y)\f$ neighbor.
    Matrix phiX2Y;
    /// \f$\phi\f$ at the \f$+2\hat x\f$ neighbor.
    Matrix phi2X;
    /// \f$\phi\f$ at the \f$+2\hat y\f$ neighbor.
    Matrix phi2Y;
    /// \f$\phi\f$ at the \f$(-\hat x,+\hat y)\f$ neighbor.
    Matrix phimXpY;
    /// \f$\phi\f$ at the \f$(+\hat x,-\hat y)\f$ neighbor.
    Matrix phipXmY;
    /// Scratch for one four-link chain product, reused for each of
    /// the `xMinus*`/`yPlus*` combinations below.
    Matrix chainA;
    /// Scratch for the other four-link chain product paired with \c
    /// chainA.
    Matrix chainB;
    /// Traceless difference of the two four-link chains centered at
    /// this cell, contributing to \f$T^{x\eta}\f$/\f$T^{y\eta}\f$.
    Matrix xMinus0;
    /// Same as \c xMinus0, evaluated at the \f$-\hat x\f$-shifted
    /// chain.
    Matrix xMinusM;
    /// Same as \c xMinus0, evaluated at the \f$+\hat x\f$-shifted
    /// chain.
    Matrix xMinusP;
    /// Same as \c xMinus0, evaluated at the transposed
    /// (\f$x\leftrightarrow y\f$) chain ordering.
    Matrix xMinusT;
    /// \f$\text{xMinus0}+\text{xMinusM}\f$.
    Matrix xMinusSum0;
    /// \f$\text{xMinusP}+\text{xMinusT}\f$.
    Matrix xMinusSum1;
    /// The \f$y\f$-oriented analog of \c xMinus0; equal to
    /// \f$-\text{xMinus0}\f$ by the chains' symmetry, so it's obtained
    /// by a sign flip rather than recomputed.
    Matrix yPlus0;
    /// Same as \c yPlus0, evaluated at the \f$-\hat y\f$-shifted
    /// chain.
    Matrix yPlusM;
    /// Same as \c yPlus0, evaluated at the \f$+\hat y\f$-shifted
    /// chain.
    Matrix yPlusP;
    /// Same as \c yPlus0, evaluated at the transposed chain ordering.
    Matrix yPlusT;
    /// \f$\text{yPlus0}+\text{yPlusM}\f$.
    Matrix yPlusSum0;
    /// \f$\text{yPlusP}+\text{yPlusT}\f$.
    Matrix yPlusSum1;
    /// This cell's covariant gradient of \f$\phi\f$ in \f$x\f$,
    /// \f$\phi_X-\phi\f$ (transported).
    Matrix covGradX0;
    /// This cell's covariant gradient of \f$\phi\f$ in \f$y\f$,
    /// \f$\phi_Y-\phi\f$ (transported).
    Matrix covGradY0;
    /// \c covGradX0-like gradient evaluated at the \f$+\hat y\f$
    /// neighbor.
    Matrix gradXAtY;
    /// \c covGradY0-like gradient evaluated at the \f$+\hat x\f$
    /// neighbor.
    Matrix gradYAtX;
    /// \c gradXAtY parallel-transported back to this cell.
    Matrix gradXAtYToPos;
    /// \c gradYAtX parallel-transported back to this cell.
    Matrix gradYAtXToPos;
    /// \f$E_1\f$ at the \f$+\hat y\f$ neighbor, parallel-transported
    /// back to this cell.
    Matrix E1AtYToPos;
    /// \f$E_2\f$ at the \f$+\hat x\f$ neighbor, parallel-transported
    /// back to this cell.
    Matrix E2AtXToPos;
};

/**
 * Computes the six off-diagonal energy-momentum tensor components
 * (\f$T^{\tau x}\f$, \f$T^{\tau y}\f$, \f$T^{\tau\eta}\f$,
 * \f$T^{xy}\f$, \f$T^{x\eta}\f$, \f$T^{y\eta}\f$) at every cell,
 * zeroing all six at the nonphysical boundary ring.
 * \param[in,out] lat Lattice to read the fields from and write the
 * six components into (via `lat->cells`).
 * \param[in] N Lattice side length.
 * \param[in] it Current time step index, used to convert to physical
 * units.
 * \param[in] dtau Time step [lattice units].
 * \param[in] g Coupling \f$g\f$.
 * \param[in] a Lattice spacing [fm].
 * \param[in] one Reusable identity matrix.
 * \param[in,out] scratch Thread-local scratch storage.
 */
void tmunuOffDiagonalTeam(
    Lattice *lat, int N, int it, double dtau, double g, double a,
    const Matrix &one, TmunuOffDiagonalScratch &scratch) {
    int pos, posX, posY, posmX, posmY, posXY, posmXpY, pospXmY, pos2X, pos2Y,
        posX2Y, pos2XY;
#pragma omp for
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = lat->positionFromXY(i, j);
            if (i == 0 || j == 0 || i == N - 1 || j == N - 1) {
                lat->cells[pos]->setTtaux(0.);
                lat->cells[pos]->setTtauy(0.);
                lat->cells[pos]->setTtaueta(0.);
                lat->cells[pos]->setTxy(0.);
                lat->cells[pos]->setTxeta(0.);
                lat->cells[pos]->setTyeta(0.);
                continue;
            }
            posX = lat->pospX[pos];
            posY = lat->pospY[pos];
            posXY = lat->positionFromXY(
                std::min(N - 1, i + 1), std::min(N - 1, j + 1));

            posmX = lat->posmX[pos];
            posmY = lat->posmY[pos];

            posmXpY = lat->posmXpY[pos];
            pospXmY = lat->pospXmY[pos];

            pos2X = lat->positionFromXY(std::min(N - 1, i + 2), j);
            pos2Y = lat->positionFromXY(i, std::min(N - 1, j + 2));

            pos2XY = lat->positionFromXY(
                std::min(N - 1, i + 2), std::min(N - 1, j + 1));
            posX2Y = lat->positionFromXY(
                std::min(N - 1, i + 1), std::min(N - 1, j + 2));

            scratch.E1 = lat->U[pos];
            scratch.E2 = lat->U2[pos];
            scratch.E1p = lat->U[posY];   // shift x value in y direction
            scratch.E2p = lat->U2[posX];  // shift y value in x direction

            scratch.pi = lat->Ux2[pos];
            scratch.piX = lat->Ux2[posX];
            scratch.piY = lat->Ux2[posY];
            scratch.piXY = lat->Ux2[posXY];

            scratch.phi = lat->Uy2[pos];
            scratch.phimX = lat->Uy2[posmX];
            scratch.phiX = lat->Uy2[posX];
            scratch.phimY = lat->Uy2[posmY];
            scratch.phiY = lat->Uy2[posY];
            scratch.phiXY = lat->Uy2[posXY];
            scratch.phimXpY = lat->Uy2[posmXpY];
            scratch.phipXmY = lat->Uy2[pospXmY];
            scratch.phi2X = lat->Uy2[pos2X];
            scratch.phi2XY = lat->Uy2[pos2XY];
            scratch.phi2Y = lat->Uy2[pos2Y];
            scratch.phiX2Y = lat->Uy2[posX2Y];

            scratch.Ux = lat->Ux[pos];
            scratch.UDx = scratch.Ux;
            scratch.UDx.conjg();

            scratch.UxmX = lat->Ux[posmX];
            scratch.UDxmX = lat->Ux[posmX];
            scratch.UDxmX.conjg();
            scratch.UxmXpY = lat->Ux[posmXpY];
            scratch.UDxmXpY = lat->Ux[posmXpY];
            scratch.UDxmXpY.conjg();

            scratch.UxpX = lat->Ux[posX];
            scratch.UxpY = lat->Ux[posY];
            scratch.UDxpX = scratch.UxpX;
            scratch.UDxpX.conjg();
            scratch.UDxpY = scratch.UxpY;
            scratch.UDxpY.conjg();

            scratch.UxpXpY = lat->Ux[posXY];
            scratch.UDxpXpY = lat->Ux[posXY];
            scratch.UDxpXpY.conjg();
            scratch.UxmXpY = lat->Ux[posmXpY];
            scratch.UDxmXpY = scratch.UxmXpY;
            scratch.UDxmXpY.conjg();

            scratch.Uy = lat->Uy[pos];
            scratch.UDy = scratch.Uy;
            scratch.UDy.conjg();

            scratch.UymY = lat->Uy[posmY];
            scratch.UDymY = lat->Uy[posmY];
            scratch.UDymY.conjg();
            scratch.UypXmY = lat->Uy[pospXmY];
            scratch.UDypXmY = lat->Uy[pospXmY];
            scratch.UDypXmY.conjg();

            scratch.UypY = lat->Uy[posY];
            scratch.UypX = lat->Uy[posX];
            scratch.UDypX = scratch.UypX;
            scratch.UDypX.conjg();

            scratch.UDypY = lat->Uy[posY];
            scratch.UDypY.conjg();

            scratch.UDxpX = lat->Ux[posX];
            scratch.UDxpX.conjg();

            scratch.Uyp2X = lat->Uy[pos2X];
            scratch.UDyp2X = scratch.Uyp2X;
            scratch.UDyp2X.conjg();
            scratch.Uxp2Y = lat->Ux[pos2Y];
            scratch.UDxp2Y = scratch.Uxp2Y;
            scratch.UDxp2Y.conjg();

            scratch.UypXpY = lat->Uy[posXY];
            scratch.UDypXpY = lat->Uy[posXY];
            scratch.UDypXpY.conjg();
            scratch.UymX = lat->Uy[posmX];
            scratch.UDymX = scratch.UymX;
            scratch.UDymX.conjg();
            scratch.UxmY = lat->Ux[posmY];
            scratch.UDxmY = scratch.UxmY;
            scratch.UDxmY.conjg();
            scratch.UypXmY = lat->Uy[pospXmY];

            // Cache the repeated four-link magnetic structures once per
            // site. The historical expressions recomputed every four-link
            // chain once for the matrix difference and again for its trace
            // subtraction, then repeated the same structures in
            // Txeta/Tyeta. Preserve the original product ordering, but
            // materialize each traceless difference only once and reuse it
            // below.
            scratch.chainA =
                scratch.Uy * scratch.UxpY * scratch.UDypX * scratch.UDx;
            scratch.chainB =
                scratch.Ux * scratch.UypX * scratch.UDxpY * scratch.UDy;
            scratch.xMinus0 = makeTmunuTracelessDifference(
                scratch.chainA, scratch.chainB, one);

            scratch.chainA =
                scratch.UDxmX * scratch.UymX * scratch.UxmXpY * scratch.UDy;
            scratch.chainB =
                scratch.Uy * scratch.UDxmXpY * scratch.UDymX * scratch.UxmX;
            scratch.xMinusM = makeTmunuTracelessDifference(
                scratch.chainA, scratch.chainB, one);

            scratch.chainA =
                scratch.UypX * scratch.UxpXpY * scratch.UDyp2X * scratch.UDxpX;
            scratch.chainB =
                scratch.UxpX * scratch.Uyp2X * scratch.UDxpXpY * scratch.UDypX;
            scratch.xMinusP = makeTmunuTracelessDifference(
                scratch.chainA, scratch.chainB, one);

            scratch.chainA =
                scratch.UDx * scratch.Uy * scratch.UxpY * scratch.UDypX;
            scratch.chainB =
                scratch.UypX * scratch.UDxpY * scratch.UDy * scratch.Ux;
            scratch.xMinusT = makeTmunuTracelessDifference(
                scratch.chainA, scratch.chainB, one);

            scratch.xMinusSum0 = scratch.xMinus0 + scratch.xMinusM;
            scratch.xMinusSum1 = scratch.xMinusP + scratch.xMinusT;

            // The first y-oriented difference is the opposite orientation
            // of scratch.xMinus0 and can be reused by a sign flip.
            scratch.yPlus0 = (-1.) * scratch.xMinus0;

            scratch.chainA =
                scratch.UDymY * scratch.UxmY * scratch.UypXmY * scratch.UDx;
            scratch.chainB =
                scratch.Ux * scratch.UDypXmY * scratch.UDxmY * scratch.UymY;
            scratch.yPlusM = makeTmunuTracelessDifference(
                scratch.chainA, scratch.chainB, one);

            scratch.chainA =
                scratch.UxpY * scratch.UypXpY * scratch.UDxp2Y * scratch.UDypY;
            scratch.chainB =
                scratch.UypY * scratch.Uxp2Y * scratch.UDypXpY * scratch.UDxpY;
            scratch.yPlusP = makeTmunuTracelessDifference(
                scratch.chainA, scratch.chainB, one);

            scratch.chainA =
                scratch.UDy * scratch.Ux * scratch.UypX * scratch.UDxpY;
            scratch.chainB =
                scratch.UxpY * scratch.UDypX * scratch.UDx * scratch.Uy;
            scratch.yPlusT = makeTmunuTracelessDifference(
                scratch.chainA, scratch.chainB, one);

            scratch.yPlusSum0 = scratch.yPlus0 + scratch.yPlusM;
            scratch.yPlusSum1 = scratch.yPlusP + scratch.yPlusT;

            // Cache covariant scalar gradients shared by Txy, Ttaueta,
            // Txeta, and Tyeta.
            scratch.covGradX0 =
                scratch.Ux * scratch.phiX * scratch.UDx - scratch.phi;
            scratch.covGradY0 =
                scratch.Uy * scratch.phiY * scratch.UDy - scratch.phi;
            scratch.gradXAtY =
                scratch.UxpY * scratch.phiXY * scratch.UDxpY - scratch.phiY;
            scratch.gradYAtX =
                scratch.UypX * scratch.phiXY * scratch.UDypX - scratch.phiX;
            scratch.gradXAtYToPos = scratch.Uy * scratch.gradXAtY * scratch.UDy;
            scratch.gradYAtXToPos = scratch.Ux * scratch.gradYAtX * scratch.UDx;

            scratch.chainA = scratch.E2 * scratch.xMinusSum0
                             + scratch.E2p * scratch.xMinusSum1;
            const complex<double> ttauxPiTrace =
                su3::traceABCD(
                    scratch.pi, scratch.Ux, scratch.phiX, scratch.UDx)
                - su3::traceABCD(
                    scratch.pi, scratch.UDxmX, scratch.phimX, scratch.UxmX)
                + su3::traceABCD(
                    scratch.piY, scratch.UxpY, scratch.phiXY, scratch.UDxpY)
                - su3::traceABCD(
                    scratch.piY, scratch.UDxmXpY, scratch.phimXpY,
                    scratch.UxmXpY)
                + su3::traceABCD(
                    scratch.piX, scratch.UxpX, scratch.phi2X, scratch.UDxpX)
                - su3::traceABCD(
                    scratch.piX, scratch.UDx, scratch.phi, scratch.Ux)
                + su3::traceABCD(
                    scratch.piXY, scratch.UxpXpY, scratch.phi2XY,
                    scratch.UDxpXpY)
                - su3::traceABCD(
                    scratch.piXY, scratch.UDxpY, scratch.phiY, scratch.UxpY);
            lat->cells[pos]->setTtaux(
                +2. / (it * dtau) / 8. * scratch.chainA.trace().imag()
                - 2. / 8. / (it * dtau) * ttauxPiTrace.real());

            scratch.chainA = scratch.E1 * scratch.yPlusSum0
                             + scratch.E1p * scratch.yPlusSum1;
            const complex<double> ttauyPiTrace =
                su3::traceABCD(
                    scratch.pi, scratch.Uy, scratch.phiY, scratch.UDy)
                - su3::traceABCD(
                    scratch.pi, scratch.UDymY, scratch.phimY, scratch.UymY)
                + su3::traceABCD(
                    scratch.piX, scratch.UypX, scratch.phiXY, scratch.UDypX)
                - su3::traceABCD(
                    scratch.piX, scratch.UDypXmY, scratch.phipXmY,
                    scratch.UypXmY)
                + su3::traceABCD(
                    scratch.piY, scratch.UypY, scratch.phi2Y, scratch.UDypY)
                - su3::traceABCD(
                    scratch.piY, scratch.UDy, scratch.phi, scratch.Uy)
                + su3::traceABCD(
                    scratch.piXY, scratch.UypXpY, scratch.phiX2Y,
                    scratch.UDypXpY)
                - su3::traceABCD(
                    scratch.piXY, scratch.UDypX, scratch.phiX, scratch.UypX);
            lat->cells[pos]->setTtauy(
                +2. / (it * dtau) / 8. * scratch.chainA.trace().imag()
                - 2. / 8. / (it * dtau) * ttauyPiTrace.real());

            const complex<double> ttauetaTrace =
                su3::traceAB(scratch.E1, scratch.covGradX0)
                + su3::traceAB(scratch.E1p, scratch.gradXAtY)
                + su3::traceAB(scratch.E2, scratch.covGradY0)
                + su3::traceAB(scratch.E2p, scratch.gradYAtX);
            lat->cells[pos]->setTtaueta(
                g / (it * dtau) / (it * dtau) / (it * dtau)
                * ttauetaTrace.real());

            scratch.E1AtYToPos = scratch.Uy * scratch.E1p * scratch.UDy;
            scratch.E2AtXToPos = scratch.Ux * scratch.E2p * scratch.UDx;
            scratch.chainA =
                -1. / 4. * g * g * (scratch.E1 + scratch.E1AtYToPos)
                    * (scratch.E2 + scratch.E2AtXToPos)
                + 1. / 4.
                      * (scratch.covGradX0 * scratch.covGradY0
                         + scratch.gradXAtYToPos * scratch.covGradY0
                         + scratch.covGradX0 * scratch.gradYAtXToPos
                         + scratch.gradXAtYToPos * scratch.gradYAtXToPos);
            lat->cells[pos]->setTxy(
                2. / (it * dtau) / (it * dtau) * scratch.chainA.trace().real());

            const complex<double> txetaElectricTrace =
                su3::traceAB(scratch.E1, scratch.pi)
                + su3::traceABCD(
                    scratch.E1, scratch.Ux, scratch.piX, scratch.UDx)
                + su3::traceAB(scratch.E1p, scratch.piY)
                + su3::traceABCD(
                    scratch.E1p, scratch.UxpY, scratch.piXY, scratch.UDxpY);
            scratch.chainA = scratch.xMinusSum0 * scratch.covGradY0
                             + scratch.xMinusSum1 * scratch.gradYAtX;
            lat->cells[pos]->setTxeta(
                -2. / (it * dtau) / (it * dtau)
                * (1. / 4. * g * txetaElectricTrace.real()
                   - 1. / 8. / g * scratch.chainA.trace().imag()));

            const complex<double> tyetaElectricTrace =
                su3::traceAB(scratch.E2, scratch.pi)
                + su3::traceABCD(
                    scratch.E2, scratch.Uy, scratch.piY, scratch.UDy)
                + su3::traceAB(scratch.E2p, scratch.piX)
                + su3::traceABCD(
                    scratch.E2p, scratch.UypX, scratch.piXY, scratch.UDypX);
            scratch.chainA = scratch.yPlusSum0 * scratch.covGradX0
                             + scratch.yPlusSum1 * scratch.gradXAtY;
            lat->cells[pos]->setTyeta(
                -2. / (it * dtau) / (it * dtau)
                * (1. / 4. * g * tyetaElectricTrace.real()
                   - 1. / 8. / g * scratch.chainA.trace().imag()));

            lat->cells[pos]->setTtaux(
                lat->cells[pos]->getTtaux() * 1 / pow(a, 4.));
            lat->cells[pos]->setTtauy(
                lat->cells[pos]->getTtauy() * 1 / pow(a, 4.));
            lat->cells[pos]->setTtaueta(
                lat->cells[pos]->getTtaueta() * 1 / pow(a, 5.));
            lat->cells[pos]->setTxy(lat->cells[pos]->getTxy() * 1 / pow(a, 4.));
            lat->cells[pos]->setTxeta(
                lat->cells[pos]->getTxeta() * 1 / pow(a, 5.));
            lat->cells[pos]->setTyeta(
                lat->cells[pos]->getTyeta() * 1 / pow(a, 5.));
        }
    }
}

/**
 * Sets \f$\epsilon = T^{\tau\tau}\f$ (before lattice-unit rescaling),
 * then rescales \f$T^{\tau\tau}\f$/\f$T^{xx}\f$/\f$T^{yy}\f$/
 * \f$T^{\eta\eta}\f$ and \f$\epsilon\f$ from lattice to physical
 * units; zeroes all five at the nonphysical boundary ring.
 * \param[in,out] lat Lattice to read/write
 * \f$\epsilon\f$/\f$T^{\tau\tau}\f$/\f$T^{xx}\f$/\f$T^{yy}\f$/
 * \f$T^{\eta\eta}\f$ in place (via `lat->cells`).
 * \param[in] N Lattice side length.
 * \param[in] a Lattice spacing [fm].
 */
void tmunuNormalizeDiagonalTeam(Lattice *lat, int N, double a) {
#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
        const int i = lat->xFromPosition(pos);
        const int j = pos - i * N;
        if (i == 0 || j == 0 || i == N - 1 || j == N - 1) {
            lat->cells[pos]->setEpsilon(0.);
            lat->cells[pos]->setTtautau(0.);
            lat->cells[pos]->setTxx(0.);
            lat->cells[pos]->setTyy(0.);
            lat->cells[pos]->setTetaeta(0.);
            continue;
        }
        lat->cells[pos]->setEpsilon(
            lat->cells[pos]->getTtautau() * 1 / pow(a, 4.));
        lat->cells[pos]->setTtautau(
            lat->cells[pos]->getTtautau() * 1 / pow(a, 4.));
        lat->cells[pos]->setTxx(lat->cells[pos]->getTxx() * 1 / pow(a, 4.));
        lat->cells[pos]->setTyy(lat->cells[pos]->getTyy() * 1 / pow(a, 4.));
        lat->cells[pos]->setTetaeta(
            lat->cells[pos]->getTetaeta() * 1 / pow(a, 6.));
    }
}

}  // namespace

void EnergyMomentumTensor::compute(Lattice *lat, Parameters *param, int it) {
    IPG_PROFILE_SCOPE("observables.Tmunu");
    int N = param->lattice.size;
    double L = param->lattice.L;
    double a = L / N;  // lattice spacing in fm
    double g = param->coupling.g;
    double dtau = param->run.dtau;
    Matrix one(1.);

#pragma omp parallel
    {
        TmunuPlaquetteScratch plaquetteScratch;
        tmunuPlaquetteTeam(lat, N, one, plaquetteScratch);

        TmunuDiagonalElectricScratch electricScratch;
        tmunuDiagonalElectricTeam(lat, N, it, dtau, g, electricScratch);

        TmunuDiagonalMagneticScratch magneticScratch;
        tmunuDiagonalMagneticTeam(lat, N, it, dtau, g, magneticScratch);

        tmunuNormalizeDiagonalTeam(lat, N, a);

        TmunuOffDiagonalScratch offDiagonalScratch;
        tmunuOffDiagonalTeam(lat, N, it, dtau, g, a, one, offDiagonalScratch);
    }  // omp parallel
}

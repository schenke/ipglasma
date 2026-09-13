// Evolution.cpp is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Evolution.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "Fragmentation.h"
#include "GaugeFix.h"
#include "Instrumentation.h"
#include "MyEigen.h"
#include "PhysConst.h"
#include "SU3.h"

using Fragmentation::kkp;
using PhysConst::hbarc;
using PhysConst::m_kaon;
using PhysConst::m_pion;
using PhysConst::m_proton;

using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;

//**************************************************************************
// Evolution class.

namespace {

struct EvolveUScratch {
    EvolveUScratch() : one(1.) {}

    Matrix E1;
    Matrix E2;
    Matrix temp1;
    Matrix temp2;
    Matrix one;
};

struct EvolvePhiScratch {
    Matrix phi;
    Matrix pi;
};

struct EvolvePiScratch {
    Matrix Ux;
    Matrix Uy;
    Matrix UxXm1;
    Matrix UyYm1;
    Matrix phi;
    Matrix phiX;
    Matrix phiY;
    Matrix phimX;
    Matrix phimY;
    Matrix bracket;
    Matrix pi;
};

struct EvolveEScratch {
    Matrix Ux;
    Matrix Uy;
    Matrix temp1;
    Matrix temp2;
    Matrix En;
    Matrix phi;
    Matrix phiN;
    Matrix U12;
    Matrix U1m2;
    Matrix U2m1;
};

inline void addEForceSU3(
    Matrix &En, const Matrix &a, const Matrix &b, double bSign,
    const Matrix &phiN, const Matrix &phi, const complex<double> coeffPlaq,
    const complex<double> coeffComm) {
    complex<double> *E = En.data();
    const complex<double> *A = a.data();
    const complex<double> *B = b.data();
    const su3::Matrix3 comm = su3::commutator(phiN, phi);

    // The plaquette force is the traceless anti-Hermitian part of
    // M = a + bSign*b.  Form it directly instead of materializing M, M^dagger,
    // an identity-matrix scale, and the associated Matrix temporaries.
    const complex<double> m00 = A[0] + bSign * B[0];
    const complex<double> m11 = A[4] + bSign * B[4];
    const complex<double> m22 = A[8] + bSign * B[8];
    const complex<double> traceThird =
        ((m00 - std::conj(m00)) + (m11 - std::conj(m11))
         + (m22 - std::conj(m22)))
        / 3.0;

    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            const int idx = 3 * row + col;
            const int tidx = 3 * col + row;
            const complex<double> mij = A[idx] + bSign * B[idx];
            const complex<double> mji = A[tidx] + bSign * B[tidx];
            complex<double> plaq = mij - std::conj(mji);
            if (row == col) plaq -= traceThird;
            E[idx] += coeffPlaq * plaq + coeffComm * comm.e[idx];
        }
    }

    // E is constrained to be traceless.  The old code subtracts
    // trace(E)/3 times the identity Matrix; only the diagonal entries change.
    const complex<double> eTraceThird = (E[0] + E[4] + E[8]) / 3.0;
    E[0] -= eTraceThird;
    E[4] -= eTraceThird;
    E[8] -= eTraceThird;
}

void evolveUTeam(
    Lattice *lat, int N, double g, double dtau, double tau,
    EvolveUScratch &scratch) {
    const int n = 2;
    const complex<double> iOmega(0., g * g * dtau / (tau + dtau / 2.));

#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
        scratch.E1 = iOmega * lat->U[pos];

        scratch.temp2 = scratch.one + 1. / (double)n * scratch.E1;
        for (int in = 0; in < n - 1; in++) {
            scratch.temp1 = scratch.E1 * scratch.temp2;
            scratch.temp2 =
                scratch.one + 1. / (double)(n - 1 - in) * scratch.temp1;
        }

        scratch.E1 = scratch.temp2;

        scratch.E2 = iOmega * lat->U2[pos];

        scratch.temp2 = scratch.one + 1. / (double)n * scratch.E2;
        for (int in = 0; in < n - 1; in++) {
            scratch.temp1 = scratch.E2 * scratch.temp2;
            scratch.temp2 =
                scratch.one + 1. / (double)(n - 1 - in) * scratch.temp1;
        }

        scratch.E2 = scratch.temp2;

        lat->Ux[pos] = (scratch.E1 * lat->Ux[pos]);
        lat->Uy[pos] = (scratch.E2 * lat->Uy[pos]);
    }
}

void evolvePhiTeam(
    Lattice *lat, int N, double dtau, double tau, EvolvePhiScratch &scratch) {
#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
        scratch.phi = lat->Uy2[pos];
        scratch.pi = lat->Ux2[pos];

        scratch.phi = scratch.phi + (tau + dtau / 2.) * dtau * scratch.pi;

        lat->Uy2[pos] = (scratch.phi);
    }
}

void evolvePiTeam(
    Lattice *lat, int N, double dtau, double tau, EvolvePiScratch &scratch) {
    const double dtauOverTau = dtau / tau;

#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
        scratch.Ux = lat->Ux[pos];
        scratch.Uy = lat->Uy[pos];
        scratch.pi = lat->Ux2[pos];
        scratch.phi = lat->Uy2[pos];

        scratch.phiX =
            scratch.Ux
            * scratch.Ux.prodABconj(lat->Uy2[lat->pospX[pos]], scratch.Ux);
        scratch.phiY =
            scratch.Uy
            * scratch.Uy.prodABconj(lat->Uy2[lat->pospY[pos]], scratch.Uy);

        scratch.UxXm1 = lat->Ux[lat->posmX[pos]];
        scratch.UyYm1 = lat->Uy[lat->posmY[pos]];

        scratch.phimX =
            scratch.Ux.prodAconjB(scratch.UxXm1, lat->Uy2[lat->posmX[pos]])
            * scratch.UxXm1;
        scratch.phimY =
            scratch.Ux.prodAconjB(scratch.UyYm1, lat->Uy2[lat->posmY[pos]])
            * scratch.UyYm1;

        scratch.bracket = scratch.phiX + scratch.phimX + scratch.phiY
                          + scratch.phimY - 4. * scratch.phi;

        scratch.pi += dtauOverTau * scratch.bracket;

        lat->Ux2[pos] = (scratch.pi);
    }
}

void evolveETeam(
    Lattice *lat, int N, double g, double dtau, double tau,
    EvolveEScratch &scratch) {
    const complex<double> coeffPlaq =
        complex<double>(0., 1.) * tau * dtau / (2. * g * g);
    const complex<double> coeffComm = complex<double>(0., 1.) * dtau / tau;

#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
        const int posmXpY = lat->posmXpY[pos];
        const int pospXmY = lat->pospXmY[pos];

        scratch.En = lat->U[pos];
        scratch.phi = lat->Uy2[pos];
        scratch.phiN = lat->Uy2[lat->pospX[pos]];
        scratch.Ux = lat->Ux[pos];
        scratch.phiN =
            scratch.Ux * scratch.Ux.prodABconj(scratch.phiN, scratch.Ux);

        scratch.Uy = lat->Uy[pos];
        scratch.temp1 = lat->Ux[lat->pospY[pos]];
        scratch.temp1.conjg();
        scratch.U12 = (scratch.Ux * lat->Uy[lat->pospX[pos]])
                      * (scratch.Ux.prodABconj(scratch.temp1, scratch.Uy));

        scratch.temp1 = lat->Ux[lat->posmY[pos]];
        scratch.temp2 = lat->Uy[pospXmY];
        scratch.U1m2 =
            (scratch.Ux.prodABconj(scratch.Ux, scratch.temp2))
            * (scratch.Ux.prodAconjB(scratch.temp1, lat->Uy[lat->posmY[pos]]));

        scratch.temp1 = lat->Uy[lat->posmX[pos]];
        scratch.temp2 = lat->Ux[posmXpY];
        scratch.U2m1 =
            (scratch.Ux.prodABconj(scratch.Uy, scratch.temp2))
            * (scratch.Ux.prodAconjB(scratch.temp1, lat->Ux[lat->posmX[pos]]));

        // U12 + U1m2 - U12^dagger - U1m2^dagger is the
        // anti-Hermitian part of U12 + U1m2.
        addEForceSU3(
            scratch.En, scratch.U12, scratch.U1m2, 1.0, scratch.phiN,
            scratch.phi, coeffPlaq, coeffComm);
        lat->U[pos] = scratch.En;

        scratch.phiN = lat->Uy2[lat->pospY[pos]];
        scratch.phiN =
            scratch.Uy * scratch.Uy.prodABconj(scratch.phiN, scratch.Uy);

        scratch.En = lat->U2[pos];
        // U12^dagger + U2m1 - U12 - U2m1^dagger is the
        // anti-Hermitian part of U2m1 - U12.
        addEForceSU3(
            scratch.En, scratch.U2m1, scratch.U12, -1.0, scratch.phiN,
            scratch.phi, coeffPlaq, coeffComm);
        lat->U2[pos] = scratch.En;
    }
}

void addTeamPhase(const char *phase, double started) {
    ipg::Profiler::instance().add(phase, ipg::wallSeconds() - started);
}

void addPhaseAndRestart(const char *phase, double &started) {
    const double now = ipg::wallSeconds();
    ipg::Profiler::instance().add(phase, now - started);
    started = now;
}

void evolveStepPersistent(
    Lattice *lat, Parameters *param, double dtau, double tau,
    bool updateCoordinates) {
    IPG_PROFILE_SCOPE("evolution.parallel_step");
    const int N = param->getSize();
    const double g = param->getg();
    double phaseStart = 0.0;

#pragma omp parallel shared(phaseStart)
    {
        EvolveUScratch uScratch;
        EvolvePhiScratch phiScratch;
        EvolvePiScratch piScratch;
        EvolveEScratch eScratch;

#pragma omp single
        {
            phaseStart = ipg::wallSeconds();
        }
        evolvePiTeam(lat, N, dtau, tau, piScratch);
#pragma omp single
        {
            addTeamPhase("evolution.evolvePi", phaseStart);
            phaseStart = ipg::wallSeconds();
        }

        evolveETeam(lat, N, g, dtau, tau, eScratch);
#pragma omp single
        {
            addTeamPhase("evolution.evolveE", phaseStart);
            phaseStart = ipg::wallSeconds();
        }

        if (updateCoordinates) {
            evolvePhiTeam(lat, N, dtau, tau, phiScratch);
#pragma omp single
            {
                addTeamPhase("evolution.evolvePhi", phaseStart);
                phaseStart = ipg::wallSeconds();
            }

            evolveUTeam(lat, N, g, dtau, tau, uScratch);
#pragma omp single
            {
                addTeamPhase("evolution.evolveU", phaseStart);
            }
        }
    }
}

inline void makeTmunuTracelessDifference(
    const Matrix &lhs, const Matrix &rhs, const Matrix &one, Matrix &out) {
    out = lhs - rhs;
    out -= (out.trace() / 3.0) * one;
}

struct TmunuPlaquetteScratch {
    Matrix UDx;
    Matrix UDy;
    Matrix Uplaq;
};

// Precomputes the spatial plaquette U_x(x) U_y(x+xhat) U_x(x+yhat)^dagger
// U_y(x)^dagger at every cell into lat->Uy1, consumed by
// tmunuDiagonalMagneticTeam below. The outermost ring is a nonphysical guard
// region (Tmunu's stencils need a genuine one-cell neighborhood), so it gets
// the identity instead of a clamped, gauge-noncovariant plaquette.
void tmunuPlaquetteTeam(
    Lattice *lat, int N, const Matrix &one, TmunuPlaquetteScratch &scratch) {
    int pos, posX, posY;
#pragma omp for
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;
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

struct TmunuDiagonalElectricScratch {
    Matrix E1;
    Matrix E2;
    Matrix E1p;
    Matrix E2p;
    Matrix pi;
    Matrix piX;
    Matrix piY;
    Matrix piXY;
};

// T^tautau, T^xx, T^yy, T^etaeta: electric (E, pi) contribution. Sets each
// field outright (rather than adding to it) since this runs before
// tmunuDiagonalMagneticTeam, which adds the magnetic/gradient contribution
// on top.
void tmunuDiagonalElectricTeam(
    Lattice *lat, int N, int it, double dtau, double g,
    TmunuDiagonalElectricScratch &scratch) {
    int pos, posX, posY, posXY;
#pragma omp for
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;
            if (i == 0 || j == 0 || i == N - 1 || j == N - 1) {
                lat->cells[pos]->setTtautau(0.);
                lat->cells[pos]->setTxx(0.);
                lat->cells[pos]->setTyy(0.);
                lat->cells[pos]->setTetaeta(0.);
                continue;
            }
            posX = lat->pospX[pos];
            posY = lat->pospY[pos];

            posXY = std::min(N - 1, i + 1) * N + std::min(N - 1, j + 1);

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

struct TmunuDiagonalMagneticScratch {
    Matrix Uplaq;
    Matrix phi;
    Matrix phiX;
    Matrix phiY;
    Matrix phiXY;
    Matrix Ux;
    Matrix Uy;
    Matrix UDx;
    Matrix UDy;
    Matrix phiTildeX;
    Matrix phiTildeY;
    Matrix phiTildeXY1;
    Matrix phiTildeXY2;
};

// T^tautau, T^xx, T^yy, T^etaeta: adds the magnetic (plaquette) and gradient
// (phi) contribution on top of whatever tmunuDiagonalElectricTeam set
// (0 at the boundary, the electric part elsewhere).
void tmunuDiagonalMagneticTeam(
    Lattice *lat, int N, int it, double dtau, double g,
    TmunuDiagonalMagneticScratch &scratch) {
    int pos, posX, posY, posXY;
#pragma omp for
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;
            if (i == 0 || j == 0 || i == N - 1 || j == N - 1) {
                continue;
            }

            posX = lat->pospX[pos];
            posY = lat->pospY[pos];

            posXY = std::min(N - 1, i + 1) * N + std::min(N - 1, j + 1);

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

struct TmunuOffDiagonalScratch {
    Matrix Ux;
    Matrix Uy;
    Matrix UxmX;
    Matrix UymY;
    Matrix UDx;
    Matrix UDy;
    Matrix UDxmX;
    Matrix UDymY;
    Matrix UDxmXpY;
    Matrix UDxpXpY;
    Matrix UxpX;
    Matrix UxpY;
    Matrix UDxpY;
    Matrix UxpXpY;
    Matrix UDypXmY;
    Matrix UypY;
    Matrix UypX;
    Matrix UDypX;
    Matrix UypXpY;
    Matrix UDypXpY;
    Matrix UymX;
    Matrix UxmXpY;
    Matrix UxmY;
    Matrix UDxmY;
    Matrix UypXmY;
    Matrix UDyp2X;
    Matrix Uyp2X;
    Matrix UDxpX;
    Matrix Uxp2Y;
    Matrix UDxp2Y;
    Matrix UDypY;
    Matrix UDymX;
    Matrix E1;
    Matrix E2;
    Matrix E1p;
    Matrix E2p;
    Matrix pi;
    Matrix piX;
    Matrix piY;
    Matrix piXY;
    Matrix phi;
    Matrix phiX;
    Matrix phiY;
    Matrix phiXY;
    Matrix phimX;
    Matrix phimY;
    Matrix phi2XY;
    Matrix phiX2Y;
    Matrix phi2X;
    Matrix phi2Y;
    Matrix phimXpY;
    Matrix phipXmY;
    Matrix chainA;
    Matrix chainB;
    Matrix xMinus0;
    Matrix xMinusM;
    Matrix xMinusP;
    Matrix xMinusT;
    Matrix xMinusSum0;
    Matrix xMinusSum1;
    Matrix yPlus0;
    Matrix yPlusM;
    Matrix yPlusP;
    Matrix yPlusT;
    Matrix yPlusSum0;
    Matrix yPlusSum1;
    Matrix covGradX0;
    Matrix covGradY0;
    Matrix gradXAtY;
    Matrix gradYAtX;
    Matrix gradXAtYToPos;
    Matrix gradYAtXToPos;
    Matrix E1AtYToPos;
    Matrix E2AtXToPos;
};

// T^taux, T^tauy, T^taueta, T^xy, T^xeta, T^yeta.
void tmunuOffDiagonalTeam(
    Lattice *lat, int N, int it, double dtau, double g, double a,
    const Matrix &one, TmunuOffDiagonalScratch &scratch) {
    int pos, posX, posY, posmX, posmY, posXY, posmXpY, pospXmY, pos2X, pos2Y,
        posX2Y, pos2XY;
#pragma omp for
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                pos = i * N + j;
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
                posXY = std::min(N - 1, i + 1) * N + std::min(N - 1, j + 1);

                posmX = lat->posmX[pos];
                posmY = lat->posmY[pos];

                posmXpY = lat->posmXpY[pos];
                pospXmY = lat->pospXmY[pos];

                pos2X = std::min(N - 1, i + 2) * N + j;
                pos2Y = i * N + std::min(N - 1, j + 2);

                pos2XY = std::min(N - 1, i + 2) * N + std::min(N - 1, j + 1);
                posX2Y = std::min(N - 1, i + 1) * N + std::min(N - 1, j + 2);

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
                makeTmunuTracelessDifference(
                    scratch.chainA, scratch.chainB, one, scratch.xMinus0);

                scratch.chainA =
                    scratch.UDxmX * scratch.UymX * scratch.UxmXpY * scratch.UDy;
                scratch.chainB =
                    scratch.Uy * scratch.UDxmXpY * scratch.UDymX * scratch.UxmX;
                makeTmunuTracelessDifference(
                    scratch.chainA, scratch.chainB, one, scratch.xMinusM);

                scratch.chainA =
                    scratch.UypX * scratch.UxpXpY * scratch.UDyp2X
                    * scratch.UDxpX;
                scratch.chainB =
                    scratch.UxpX * scratch.Uyp2X * scratch.UDxpXpY
                    * scratch.UDypX;
                makeTmunuTracelessDifference(
                    scratch.chainA, scratch.chainB, one, scratch.xMinusP);

                scratch.chainA =
                    scratch.UDx * scratch.Uy * scratch.UxpY * scratch.UDypX;
                scratch.chainB =
                    scratch.UypX * scratch.UDxpY * scratch.UDy * scratch.Ux;
                makeTmunuTracelessDifference(
                    scratch.chainA, scratch.chainB, one, scratch.xMinusT);

                scratch.xMinusSum0 = scratch.xMinus0 + scratch.xMinusM;
                scratch.xMinusSum1 = scratch.xMinusP + scratch.xMinusT;

                // The first y-oriented difference is the opposite orientation
                // of scratch.xMinus0 and can be reused by a sign flip.
                scratch.yPlus0 = (-1.) * scratch.xMinus0;

                scratch.chainA =
                    scratch.UDymY * scratch.UxmY * scratch.UypXmY * scratch.UDx;
                scratch.chainB =
                    scratch.Ux * scratch.UDypXmY * scratch.UDxmY * scratch.UymY;
                makeTmunuTracelessDifference(
                    scratch.chainA, scratch.chainB, one, scratch.yPlusM);

                scratch.chainA =
                    scratch.UxpY * scratch.UypXpY * scratch.UDxp2Y
                    * scratch.UDypY;
                scratch.chainB =
                    scratch.UypY * scratch.Uxp2Y * scratch.UDypXpY
                    * scratch.UDxpY;
                makeTmunuTracelessDifference(
                    scratch.chainA, scratch.chainB, one, scratch.yPlusP);

                scratch.chainA =
                    scratch.UDy * scratch.Ux * scratch.UypX * scratch.UDxpY;
                scratch.chainB =
                    scratch.UxpY * scratch.UDypX * scratch.UDx * scratch.Uy;
                makeTmunuTracelessDifference(
                    scratch.chainA, scratch.chainB, one, scratch.yPlusT);

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
                scratch.gradXAtYToPos =
                    scratch.Uy * scratch.gradXAtY * scratch.UDy;
                scratch.gradYAtXToPos =
                    scratch.Ux * scratch.gradYAtX * scratch.UDx;

                scratch.chainA =
                    scratch.E2 * scratch.xMinusSum0
                    + scratch.E2p * scratch.xMinusSum1;
                const complex<double> ttauxPiTrace =
                    su3::traceABCD(
                        scratch.pi, scratch.Ux, scratch.phiX, scratch.UDx)
                    - su3::traceABCD(
                        scratch.pi, scratch.UDxmX, scratch.phimX, scratch.UxmX)
                    + su3::traceABCD(
                        scratch.piY, scratch.UxpY, scratch.phiXY,
                        scratch.UDxpY)
                    - su3::traceABCD(
                        scratch.piY, scratch.UDxmXpY, scratch.phimXpY,
                        scratch.UxmXpY)
                    + su3::traceABCD(
                        scratch.piX, scratch.UxpX, scratch.phi2X,
                        scratch.UDxpX)
                    - su3::traceABCD(
                        scratch.piX, scratch.UDx, scratch.phi, scratch.Ux)
                    + su3::traceABCD(
                        scratch.piXY, scratch.UxpXpY, scratch.phi2XY,
                        scratch.UDxpXpY)
                    - su3::traceABCD(
                        scratch.piXY, scratch.UDxpY, scratch.phiY,
                        scratch.UxpY);
                lat->cells[pos]->setTtaux(
                    +2. / (it * dtau) / 8. * scratch.chainA.trace().imag()
                    - 2. / 8. / (it * dtau) * ttauxPiTrace.real());

                scratch.chainA =
                    scratch.E1 * scratch.yPlusSum0
                    + scratch.E1p * scratch.yPlusSum1;
                const complex<double> ttauyPiTrace =
                    su3::traceABCD(
                        scratch.pi, scratch.Uy, scratch.phiY, scratch.UDy)
                    - su3::traceABCD(
                        scratch.pi, scratch.UDymY, scratch.phimY, scratch.UymY)
                    + su3::traceABCD(
                        scratch.piX, scratch.UypX, scratch.phiXY,
                        scratch.UDypX)
                    - su3::traceABCD(
                        scratch.piX, scratch.UDypXmY, scratch.phipXmY,
                        scratch.UypXmY)
                    + su3::traceABCD(
                        scratch.piY, scratch.UypY, scratch.phi2Y,
                        scratch.UDypY)
                    - su3::traceABCD(
                        scratch.piY, scratch.UDy, scratch.phi, scratch.Uy)
                    + su3::traceABCD(
                        scratch.piXY, scratch.UypXpY, scratch.phiX2Y,
                        scratch.UDypXpY)
                    - su3::traceABCD(
                        scratch.piXY, scratch.UDypX, scratch.phiX,
                        scratch.UypX);
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
                    2. / (it * dtau) / (it * dtau)
                    * scratch.chainA.trace().real());

                const complex<double> txetaElectricTrace =
                    su3::traceAB(scratch.E1, scratch.pi)
                    + su3::traceABCD(
                        scratch.E1, scratch.Ux, scratch.piX, scratch.UDx)
                    + su3::traceAB(scratch.E1p, scratch.piY)
                    + su3::traceABCD(
                        scratch.E1p, scratch.UxpY, scratch.piXY,
                        scratch.UDxpY);
                scratch.chainA =
                    scratch.xMinusSum0 * scratch.covGradY0
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
                        scratch.E2p, scratch.UypX, scratch.piXY,
                        scratch.UDypX);
                scratch.chainA =
                    scratch.yPlusSum0 * scratch.covGradX0
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
                lat->cells[pos]->setTxy(
                    lat->cells[pos]->getTxy() * 1 / pow(a, 4.));
                lat->cells[pos]->setTxeta(
                    lat->cells[pos]->getTxeta() * 1 / pow(a, 5.));
                lat->cells[pos]->setTyeta(
                    lat->cells[pos]->getTyeta() * 1 / pow(a, 5.));
            }
        }
}

// epsilon = T^tautau (before lattice-unit rescaling), then rescales
// T^tautau/T^xx/T^yy/T^etaeta and epsilon from lattice to physical units.
void tmunuNormalizeDiagonalTeam(Lattice *lat, int N, double a) {
#pragma omp for
    for (int pos = 0; pos < N * N; pos++) {
        const int i = pos / N;
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

// Computes the local-coupling factor g^2/(4 pi alpha_s) at one cell, used to
// rescale Tmunu/epsilon-derived quantities when running coupling is enabled
// (alpha_s runs with either the local Qs at this cell or one of the
// event-averaged Qs choices, per param->getRunWithLocalQs()/getRunWithQs()).
// Returns 1 when running coupling is disabled.
double computeRunningCouplingGfactor(
    Lattice *lat, Parameters *param, int pos, int N, double a, double g,
    double c, double muZero) {
    double g2mu2A, g2mu2B, alphas = 0., Qs = 0.;
    if (param->getRunningCoupling()) {
        if (pos / N > 0 && pos / N < N - 1 && pos % N > 0
            && pos % N < N - 1) {
            g2mu2A = lat->cells[pos]->getg2mu2A();
        } else
            g2mu2A = 0;

        if (pos / N > 0 && pos / N < N - 1 && pos % N > 0
            && pos % N < N - 1) {
            g2mu2B = lat->cells[pos]->getg2mu2B();
        } else
            g2mu2B = 0;

        if (param->getRunWithQs() == 2) {
            if (g2mu2A > g2mu2B)
                Qs = sqrt(
                    g2mu2A * param->getQsmuRatio()
                    * param->getQsmuRatio() / a / a * hbarc * hbarc
                    * param->getg() * param->getg());
            else
                Qs = sqrt(
                    g2mu2B * param->getQsmuRatio()
                    * param->getQsmuRatio() / a / a * hbarc * hbarc
                    * param->getg() * param->getg());
        } else if (param->getRunWithQs() == 0) {
            if (g2mu2A < g2mu2B)
                Qs = sqrt(
                    g2mu2A * param->getQsmuRatio()
                    * param->getQsmuRatio() / a / a * hbarc * hbarc
                    * param->getg() * param->getg());
            else
                Qs = sqrt(
                    g2mu2B * param->getQsmuRatio()
                    * param->getQsmuRatio() / a / a * hbarc * hbarc
                    * param->getg() * param->getg());
        } else if (param->getRunWithQs() == 1) {
            Qs = sqrt(
                (g2mu2A + g2mu2B) / 2. * param->getQsmuRatio()
                * param->getQsmuRatio() / a / a * hbarc * hbarc
                * param->getg() * param->getg());
        }

        if (param->getRunWithLocalQs() == 1) {
            // 3 flavors
            alphas = 4. * M_PI
                     / (9.
                        * log(pow(
                            pow(muZero / 0.2, 2. / c)
                                + pow(
                                    param->getRunWithThisFactorTimesQs()
                                        * Qs / 0.2,
                                    2. / c),
                            c)));
            return g * g / (4. * M_PI * alphas);
            // run with the local (in transverse plane) coupling
        } else {
            if (param->getRunWithQs() == 0)
                alphas =
                    4. * M_PI
                    / (9.
                       * log(pow(
                           pow(muZero / 0.2, 2. / c)
                               + pow(
                                   param->getRunWithThisFactorTimesQs()
                                       * param->getAverageQsmin() / 0.2,
                                   2. / c),
                           c)));
            else if (param->getRunWithQs() == 1)
                alphas =
                    4. * M_PI
                    / (9.
                       * log(pow(
                           pow(muZero / 0.2, 2. / c)
                               + pow(
                                   param->getRunWithThisFactorTimesQs()
                                       * param->getAverageQsAvg() / 0.2,
                                   2. / c),
                           c)));
            else if (param->getRunWithQs() == 2)
                alphas =
                    4. * M_PI
                    / (9.
                       * log(pow(
                           pow(muZero / 0.2, 2. / c)
                               + pow(
                                   param->getRunWithThisFactorTimesQs()
                                       * param->getAverageQs() / 0.2,
                                   2. / c),
                           c)));

            return g * g / (4. * M_PI * alphas);
        }
    } else
        return 1.;
}

// Fills E1[pos] from sourceField[pos] (one of lat->U/U2/Ux2), scaling by
// sqrt of the local running-coupling gfactor unless alpha_s runs with k_T
// (in which case the k_T-dependent factor is applied later, per-mode, in
// accumulateGluonSpectrum instead).
void prepareSpectrumField(
    Lattice *lat, Parameters *param, int N, double a, double g, double c,
    double muZero, const std::vector<Matrix> &sourceField,
    std::vector<Matrix *> &E1) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int pos = i * N + j;
            double gfactor = computeRunningCouplingGfactor(
                lat, param, pos, N, a, g, c, muZero);
            if (param->getRunWithkt() == 0) {
                *E1[pos] = sourceField[pos] * sqrt(gfactor);
            } else {
                *E1[pos] = sourceField[pos];
            }
        }
    }
}

// Accumulates E1's (already FFT'd) momentum-space spectrum into
// dNdeta/dEdeta and the n/E/n2 k_T bins. useElectricNormalization selects
// nkt's electric-field (E1/E2 passes: g^2/((it-0.5)dtau)) vs pi-field
// ((it-0.5)dtau, no g^2) normalization. accumulateCounter records bin
// occupancy into counter[]; only one of the three spectrum passes needs to,
// since all three share the same k_T grid.
void accumulateGluonSpectrum(
    Parameters *param, int N, int it, double dtau, double g, double a,
    double c, double muZero, double dkt, int bins,
    const std::vector<Matrix *> &E1, bool useElectricNormalization,
    bool accumulateCounter, double &dNdeta, double &dEdeta, double *n,
    double *E, double *n2, int *counter) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double nkt = 0.;
            int pos = i * N + j;
            int npos = (N - i) * N + (N - j);

            double kx = 2. * M_PI
                        * (-0.5
                           + static_cast<double>(i) / static_cast<double>(N));
            double ky = 2. * M_PI
                        * (-0.5
                           + static_cast<double>(j) / static_cast<double>(N));
            double kt2 = 4.
                         * (sin(kx / 2.) * sin(kx / 2.)
                            + sin(ky / 2.) * sin(ky / 2.));
            double omega2 = 4.
                            * (sin(kx / 2.) * sin(kx / 2.)
                               + sin(ky / 2.)
                                     * sin(ky / 2.));  // lattice dispersion
                                                        // relation (this is
                                                        // omega squared)

            // i=0 or j=0 have no negative k_T value available
            if (i != 0 && j != 0) {
                if (omega2 != 0) {
                    if (useElectricNormalization) {
                        nkt = 2. / sqrt(omega2) / static_cast<double>(N * N)
                              * (g * g / ((it - 0.5) * dtau)
                                 * ((((*E1[pos]) * (*E1[npos])).trace())
                                        .real()));
                    } else {
                        nkt = 2. / sqrt(omega2) / static_cast<double>(N * N)
                              * (((it - 0.5) * dtau)
                                 * ((((*E1[pos]) * (*E1[npos])).trace())
                                        .real()));
                    }
                    if (param->getRunWithkt() == 1) {
                        nkt *=
                            g * g
                            / (4. * M_PI * 4. * M_PI
                               / (9.
                                  * log(pow(
                                      pow(muZero / 0.2, 2. / c)
                                          + pow(
                                              param->getRunWithThisFactorTimesQs()
                                                  * sqrt(kt2) * hbarc / a / 0.2,
                                              2. / c),
                                      c))));
                    }
                }

                dNdeta += nkt;
                dEdeta += nkt * sqrt(omega2) * hbarc / a;

                for (int ik = 0; ik < bins; ik++) {
                    if (abs(sqrt(kt2)) > ik * dkt
                        && abs(sqrt(kt2)) <= (ik + 1) * dkt) {
                        n[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2) * 2 * M_PI
                                 * sqrt(kt2) * dkt * N * N / M_PI / M_PI / 2.
                                 / 2.;
                        E[ik] += sqrt(omega2) * hbarc / a * nkt / dkt / 2
                                 / M_PI / sqrt(kt2) * 2 * M_PI * sqrt(kt2)
                                 * dkt * N * N / M_PI / M_PI / 2. / 2.;
                        n2[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2);
                        // dividing by bin size; bin is dkt times Jacobian
                        // k(=ik*dkt) times 2Pi in phi times the correct
                        // number of counts for an infinite lattice: area in
                        // bin divided by total area
                        if (accumulateCounter) {
                            counter[ik] += 1;  // number of entries in n[ik]
                        }
                    }
                }
            }
        }
    }
}

}  // namespace

void Evolution::evolveU(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    IPG_PROFILE_SCOPE("evolution.evolveU");
    const int N = param->getSize();
    const double g = param->getg();

#pragma omp parallel
    {
        EvolveUScratch scratch;
        evolveUTeam(lat, N, g, dtau, tau, scratch);
    }
}

void Evolution::evolvePhi(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    IPG_PROFILE_SCOPE("evolution.evolvePhi");
    const int N = param->getSize();

#pragma omp parallel
    {
        EvolvePhiScratch scratch;
        evolvePhiTeam(lat, N, dtau, tau, scratch);
    }
}

void Evolution::evolvePi(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    IPG_PROFILE_SCOPE("evolution.evolvePi");
    const int N = param->getSize();

#pragma omp parallel
    {
        EvolvePiScratch scratch;
        evolvePiTeam(lat, N, dtau, tau, scratch);
    }
}

void Evolution::evolveE(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    IPG_PROFILE_SCOPE("evolution.evolveE");
    const int N = param->getSize();
    const double g = param->getg();

#pragma omp parallel
    {
        EvolveEScratch scratch;
        evolveETeam(lat, N, g, dtau, tau, scratch);
    }
}

void Evolution::checkGaussLaw(Lattice *lat, Parameters *param) {
    IPG_PROFILE_SCOPE("diagnostics.gauss_law");
    const int N = param->getSize();

    Matrix Ux;
    Matrix UxXm1;
    Matrix UxYm1;

    Matrix Uy;
    Matrix UyXm1;
    Matrix UyYm1;

    Matrix UxDag;
    Matrix UxXm1Dag;
    Matrix UxYm1Dag;

    Matrix UyDag;
    Matrix UyXm1Dag;
    Matrix UyYm1Dag;

    Matrix E1;
    Matrix E2;
    Matrix E1mX;
    Matrix E2mY;
    Matrix phi;
    Matrix pi;

    Matrix Gauss;
    double largest = 0;

    for (int pos = 0; pos < N * N; pos++) {
        // retrieve current Ux and Uy
        Ux = lat->Ux[pos];
        Uy = lat->Uy[pos];
        UxDag = Ux;
        UxDag.conjg();
        UyDag = Uy;
        UyDag.conjg();

        UxXm1 = lat->Ux[lat->posmX[pos]];
        UxYm1 = lat->Ux[lat->posmY[pos]];
        UxXm1Dag = UxXm1;
        UxXm1Dag.conjg();
        UxYm1Dag = UxYm1;
        UxYm1Dag.conjg();

        UyXm1 = lat->Uy[lat->posmX[pos]];
        UyYm1 = lat->Uy[lat->posmY[pos]];
        UyXm1Dag = UyXm1;
        UyXm1Dag.conjg();
        UyYm1Dag = UyYm1;
        UyYm1Dag.conjg();

        // retrieve current E1 and E2 (that's the one defined at tau-dtau/2)
        E1 = lat->U[pos];
        E2 = lat->U2[pos];
        E1mX = lat->U[lat->posmX[pos]];
        E2mY = lat->U2[lat->posmY[pos]];
        // retrieve current phi (at time tau) at this x_T
        phi = lat->Uy2[pos];
        // retrieve current pi
        pi = lat->Ux2[pos];

        Gauss = UxXm1Dag * E1mX * UxXm1 - E1 + UyYm1Dag * E2mY * UyYm1 - E2
                - complex<double>(0., 1.) * (phi * pi - pi * phi);

        if (Gauss.square() > largest) largest = Gauss.square();
    }
    messager_ << "[Evolution::checkGaussLaw]: Gauss violation=" << largest;
    messager_.flush("info");
}

void Evolution::writeEvolvedFields(Lattice *lat, Parameters *param, int it) {
    IPG_PROFILE_SCOPE("output.evolved_fields");
    const int N = param->getSize();
    constexpr int Nc = 3;
    const double a = param->getL() / static_cast<double>(N);
    const double dtau = param->getdtau();
    const double tauLattice = static_cast<double>(it) * dtau;
    const double tauFm = a * tauLattice;
    const double momentumTauFm =
        (it == 0) ? 0.0 : a * (static_cast<double>(it) - 0.5) * dtau;

    // The payload layout is [field, real_or_imag, x, y, row, col], C-order.
    // Full matrices are stored so that no color information is discarded.
    constexpr int nFields = 6;
    const std::size_t matrixElements =
        static_cast<std::size_t>(N) * N * Nc * Nc;
    const std::size_t payloadElements =
        static_cast<std::size_t>(nFields) * 2 * matrixElements;
    std::vector<float> payload(payloadElements);

    auto matrixAt = [lat](const int field, const int pos) -> const Matrix & {
        switch (field) {
            case 0:
                return lat->Uy2[pos];
            case 1:
                return lat->Ux2[pos];
            case 2:
                return lat->U[pos];
            case 3:
                return lat->U2[pos];
            case 4:
                return lat->Ux[pos];
            case 5:
                return lat->Uy[pos];
            default:
                throw std::runtime_error("invalid evolved-field index");
        }
    };

    for (int field = 0; field < nFields; ++field) {
        const std::size_t realOffset =
            static_cast<std::size_t>(2 * field) * matrixElements;
        const std::size_t imagOffset = realOffset + matrixElements;
        for (int x = 0; x < N; ++x) {
            for (int y = 0; y < N; ++y) {
                const int pos = x * N + y;
                const Matrix &matrix = matrixAt(field, pos);
                const std::complex<double> *elements = matrix.data();
                const std::size_t siteOffset =
                    static_cast<std::size_t>(pos) * Nc * Nc;
                for (int row = 0; row < Nc; ++row) {
                    for (int col = 0; col < Nc; ++col) {
                        const std::size_t element =
                            static_cast<std::size_t>(row) * Nc + col;
                        payload[realOffset + siteOffset + element] =
                            static_cast<float>(elements[element].real());
                        payload[imagOffset + siteOffset + element] =
                            static_cast<float>(elements[element].imag());
                    }
                }
            }
        }
    }

    // This binary format is explicitly little-endian. IP-Glasma production
    // platforms are normally little-endian; fail loudly rather than emit an
    // ambiguous file on another architecture.
    const std::uint16_t endianProbe = 1;
    if (*reinterpret_cast<const unsigned char *>(&endianProbe) != 1) {
        throw std::runtime_error(
            "writeEvolvedFields currently requires a little-endian host");
    }

    std::stringstream metadata;
    metadata << std::setprecision(17)
             << "{\"format\":\"ipglasma-evolved-fields\","
             << "\"version\":1,"
             << "\"dtype\":\"<f4\","
             << "\"shape\":[6,2," << N << "," << N << "," << Nc << "," << Nc
             << "],"
             << "\"axis_order\":[\"field\",\"complex_part\",\"x\",\"y\","
                "\"row\",\"col\"],"
             << "\"fields\":[\"phi\",\"pi\",\"E1\",\"E2\",\"Ux\","
                "\"Uy\"],"
             << "\"complex_part\":[\"real\",\"imag\"],"
             << "\"native_site_index\":\"pos=x*N+y\","
             << "\"event_id\":" << param->getEventId() << ","
             << "\"step\":" << it << ","
             << "\"tau_lattice\":" << tauLattice << ","
             << "\"tau_fm\":" << tauFm << ","
             << "\"momentum_tau_fm\":" << momentumTauFm << ","
             << "\"a_fm\":" << a << ","
             << "\"dtau_lattice\":" << dtau << ","
             << "\"staggering\":\"Ux,Uy,phi at tau; E1,E2,pi at tau-dtau/2 "
                "for step>0; all variables are the initialized tau=0+ values "
                "for step=0\"}";
    const std::string metadataString = metadata.str();

    stringstream filename;
    filename << "evolvedFields" << param->getEventId() << "_it" << std::setw(8)
             << std::setfill('0') << it << ".ipgf";

    ofstream output(
        filename.str().c_str(),
        std::ios::out | std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "could not open evolved-field snapshot " + filename.str());
    }

    const char magic[8] = {'I', 'P', 'G', 'F', 'L', 'D', '1', '\0'};
    const std::uint64_t metadataBytes =
        static_cast<std::uint64_t>(metadataString.size());
    output.write(magic, sizeof(magic));
    output.write(
        reinterpret_cast<const char *>(&metadataBytes), sizeof(metadataBytes));
    output.write(metadataString.data(), metadataString.size());
    output.write(
        reinterpret_cast<const char *>(payload.data()),
        static_cast<std::streamsize>(payload.size() * sizeof(float)));
    output.close();

    if (!output) {
        throw std::runtime_error(
            "failed while writing evolved-field snapshot " + filename.str());
    }
    messager_ << "[Evolution::writeEvolvedFields]: Wrote evolved fields "
                 "at tau="
              << tauFm << " fm/c to " << filename.str();
    messager_.flush("info");
}

void Evolution::writeGluonMultiplicityTarget(
    Parameters *param, int it, double a, double dtau, double dNPrimary,
    double dNBinned, double dEPrimary, double dEBinned, double dNCut3,
    double dECut3, double dNCut6, double dECut6, const double *spectrumN,
    const double *spectrumE, const int *spectrumCounts, int bins, double dkt) {
    IPG_PROFILE_SCOPE("output.gluon_target");
    stringstream filename;
    filename << "gluonMultiplicity" << param->getEventId() << ".json";
    ofstream output(filename.str().c_str(), std::ios::out | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "could not open gluon-multiplicity target " + filename.str());
    }

    const char *rapidityVariable =
        (param->getUsePseudoRapidity() == 0) ? "y" : "eta";
    const double meanKt = (dNPrimary != 0.0) ? dEPrimary / dNPrimary : 0.0;
    const double spectrumUnitFactor = (a / hbarc) * (a / hbarc);

    output << std::setprecision(17) << "{\n"
           << "  \"format\": \"ipglasma-gluon-target\",\n"
           << "  \"version\": 1,\n"
           << "  \"event_id\": " << param->getEventId() << ",\n"
           << "  \"step\": " << it << ",\n"
           << "  \"tau_fm\": " << static_cast<double>(it) * dtau * a << ",\n"
           << "  \"rapidity_variable\": \"" << rapidityVariable << "\",\n"
           << "  \"dN\": " << dNPrimary << ",\n"
           << "  \"dN_binned_check\": " << dNBinned << ",\n"
           << "  \"dE_GeV\": " << dEPrimary << ",\n"
           << "  \"dE_binned_check_GeV\": " << dEBinned << ",\n"
           << "  \"mean_kT_GeV\": " << meanKt << ",\n"
           << "  \"dN_kT_gt_3_GeV\": " << dNCut3 << ",\n"
           << "  \"dE_kT_gt_3_GeV\": " << dECut3 << ",\n"
           << "  \"dN_kT_gt_6_GeV\": " << dNCut6 << ",\n"
           << "  \"dE_kT_gt_6_GeV\": " << dECut6 << ",\n"
           << "  \"Npart\": " << param->getNpart() << ",\n"
           << "  \"Tpp\": " << param->getTpp() << ",\n"
           << "  \"impact_parameter_fm\": " << param->getb() << ",\n"
           << "  \"random_seed\": " << param->getRandomSeed() << ",\n"
           << "  \"spectrum_definition\": \"azimuthally averaged Coulomb-gauge "
              "gluon spectrum used by Evolution::multiplicity\",\n"
           << "  \"kt_GeV\": [";

    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << (static_cast<double>(ik) + 0.5) * dkt / a * hbarc;
    }
    output << "],\n  \"dN_d2k_GeV_minus2\": [";
    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << spectrumN[ik] * spectrumUnitFactor;
    }
    output << "],\n  \"dE_d2k_GeV_minus1\": [";
    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << spectrumE[ik] * spectrumUnitFactor;
    }
    output << "],\n  \"lattice_bin_counts\": [";
    for (int ik = 0; ik < bins; ++ik) {
        if (ik != 0) output << ",";
        output << spectrumCounts[ik];
    }
    output << "]\n}\n";
    output.close();

    if (!output) {
        throw std::runtime_error(
            "failed while writing gluon-multiplicity target " + filename.str());
    }
    messager_ << "[Evolution::writeGluonMultiplicityTarget]: Wrote gluon "
                 "target dN/d"
              << rapidityVariable << "=" << dNPrimary << " to "
              << filename.str();
    messager_.flush("info");
}

void Evolution::run(Lattice *lat, Group *group, Parameters *param) {
    IPG_PROFILE_SCOPE("evolution.total");
    int pos;
    int N = param->getSize();
    double g = param->getg();
    double L = param->getL();
    double a = L / N;  // lattice spacing in fm
    double x, y;

    double alphas = 0.;
    double gfactor;
    double Qs = 0., g2mu2A, g2mu2B;
    double muZero = param->getMuZero();
    double c = param->getc();

    // do the first half step of the momenta (E1,E2,pi)
    // for now I use the \tau=0 value at \tau=d\tau/2.
    double dtau = param->getdtau();  // dtau is in lattice units

    double maxtime = param->getMaxtime();  // maxtime is in fm
    if (param->getInverseQsForMaxTime() == 1) {
        maxtime = 1. / param->getAverageQs() * hbarc;
        messager_ << "[Evolution::run]: maximal evolution time = " << maxtime
                  << " fm";
        messager_.flush("info");
    }

    // E and Pi at tau=dtau/2 are equal to the initial ones (at tau=0)
    // now evolve phi and U to time tau=dtau.
    {
        IPG_PROFILE_SCOPE("evolution.initial_coordinate_half_step");
        evolvePhi(lat, param, dtau, 0.);
        evolveU(lat, param, dtau, 0.);
    }

    int itmax = static_cast<int>(maxtime / (a * dtau) + 0.00000000001);
    int it0 = static_cast<int>(0.1 / (a * dtau) + 0.0000000001);
    int it1 = static_cast<int>(0.2 / (a * dtau) + 0.0000000001);
    int it2 = static_cast<int>(0.3 / (a * dtau) + 0.0000000001);
    int it3 = static_cast<int>(0.4 / (a * dtau) + 0.0000000001);

    // Tmunu is defined at the integer coordinate time tau_n, while the
    // leapfrog momenta E1, E2, and pi live at tau_{n-1/2}.  Keep reusable
    // backups so measurements can temporarily center the momenta without
    // changing the actual evolution trajectory or making it output-dependent.
    const std::size_t latticeSites = static_cast<std::size_t>(N) * N;
    std::vector<Matrix> tmunuE1Backup(latticeSites);
    std::vector<Matrix> tmunuE2Backup(latticeSites);
    std::vector<Matrix> tmunuPiBackup(latticeSites);

    messager_ << "[Evolution::run]: Starting evolution: num of time steps="
              << itmax;
    messager_.flush("info");
    if ((param->getWriteOutputs() == 5)) {
        messager_ << "[Evolution::run]: Measuring at times " << it0 * a * dtau
                  << ", " << it1 * a * dtau << ", " << it2 * a * dtau << ", "
                  << it3 * a * dtau << ", " << itmax * a * dtau << ". ";
        messager_.flush("info");
    }
    messager_ << "[Evolution::run]:  a = " << a;
    messager_.flush("info");
    messager_ << "[Evolution::run]:  dtau = " << dtau;
    messager_.flush("info");
    messager_ << "[Evolution::run]:  it0 = " << it0;
    messager_.flush("info");

    // do evolution
    for (int it = 1; it <= itmax; it++) {
        const bool finalTmunuMeasurement = (it == itmax);
        const bool intermediateTmunuMeasurement =
            (param->getWriteOutputs() == 5)
            && (it == it0 || it == it1 || it == it2 || it == it3);
        const bool measureTmunu =
            finalTmunuMeasurement || intermediateTmunuMeasurement;

        if (measureTmunu) {
            std::copy(lat->U.begin(), lat->U.end(), tmunuE1Backup.begin());
            std::copy(lat->U2.begin(), lat->U2.end(), tmunuE2Backup.begin());
            std::copy(lat->Ux2.begin(), lat->Ux2.end(), tmunuPiBackup.begin());

            // Temporarily move E1, E2, and pi from tau_{n-1/2} to tau_n.
            // Coordinates are not advanced.  The original momenta are restored
            // after output, so enabling Tmunu measurements cannot change the
            // subsequent leapfrog trajectory.
            evolveStepPersistent(lat, param, dtau / 2., it * dtau, false);
        }

        if (finalTmunuMeasurement) {
            tmunu(lat, param, it);
            // Hydro flow fields are optional. Tmunu output remains available
            // through the lightweight writer when the expensive eigen solve is
            // disabled.
            if (param->getWriteEpsilonUHydro() != 0) {
                u(lat, param, it, true);
            } else {
                MyEigen myeigen;
                myeigen.writeTmunu4D(lat, param, it);
            }
        }

        if (intermediateTmunuMeasurement) {
            tmunu(lat, param, it);
            //  Preserve the historical intermediate-time finalFlag=false path
            //  when hydro output is enabled.
            if (param->getWriteEpsilonUHydro() != 0) {
                u(lat, param, it, false);
            } else {
                MyEigen myeigen;
                myeigen.writeTmunu4D(lat, param, it);
            }
        }

        if (measureTmunu) {
            std::copy(
                tmunuE1Backup.begin(), tmunuE1Backup.end(), lat->U.begin());
            std::copy(
                tmunuE2Backup.begin(), tmunuE2Backup.end(), lat->U2.begin());
            std::copy(
                tmunuPiBackup.begin(), tmunuPiBackup.end(), lat->Ux2.begin());
        }

        if (it % 10 == 1) {
            messager_ << "[Evolution::run]: Evolving to time " << it * a * dtau
                      << " fm/c";
            messager_.flush("info");
        }

        // Keep one OpenMP team alive across all leapfrog kernels in this
        // time step. The worksharing loops retain their implicit barriers,
        // preserving the Pi -> E -> phi -> U update order.
        if (it < itmax) {
            evolveStepPersistent(lat, param, dtau, (it)*dtau, true);
        } else {
            evolveStepPersistent(lat, param, dtau / 2., (it)*dtau, false);
        }

        if (it == 1 && param->getWriteOutputs() == 3) {
            IPG_PROFILE_SCOPE("output.epsilon_initial_text");
            stringstream streI_name;
            streI_name << "epsilonInitialPlot" << param->getEventId() << ".dat";
            string eI_name;
            eI_name = streI_name.str();

            ofstream foutEps(eI_name.c_str(), std::ios::out);
            for (int ix = 0; ix < N; ix++)  // loop over all positions
            {
                for (int iy = 0; iy < N; iy++) {
                    pos = ix * N + iy;
                    x = -L / 2. + a * ix;
                    y = -L / 2. + a * iy;

                    if (param->getRunningCoupling()) {
                        if (pos / N > 0 && pos / N < N - 1 && pos % N > 0
                            && pos % N < N - 1) {
                            g2mu2A = lat->cells[pos]->getg2mu2A();
                        } else
                            g2mu2A = 0;

                        if (pos / N > 0 && pos / N < N - 1 && pos % N > 0
                            && pos % N < N - 1) {
                            g2mu2B = lat->cells[pos]->getg2mu2B();
                        } else
                            g2mu2B = 0;

                        if (param->getRunWithQs() == 2) {
                            if (g2mu2A > g2mu2B)
                                Qs = sqrt(
                                    g2mu2A * param->getQsmuRatio()
                                    * param->getQsmuRatio() / a / a * hbarc
                                    * hbarc * param->getg() * param->getg());
                            else
                                Qs = sqrt(
                                    g2mu2B * param->getQsmuRatio()
                                    * param->getQsmuRatio() / a / a * hbarc
                                    * hbarc * param->getg() * param->getg());
                        } else if (param->getRunWithQs() == 0) {
                            if (g2mu2A < g2mu2B)
                                Qs = sqrt(
                                    g2mu2A * param->getQsmuRatio()
                                    * param->getQsmuRatio() / a / a * hbarc
                                    * hbarc * param->getg() * param->getg());
                            else
                                Qs = sqrt(
                                    g2mu2B * param->getQsmuRatio()
                                    * param->getQsmuRatio() / a / a * hbarc
                                    * hbarc * param->getg() * param->getg());
                        } else if (param->getRunWithQs() == 1) {
                            Qs = sqrt(
                                (g2mu2A + g2mu2B) / 2. * param->getQsmuRatio()
                                * param->getQsmuRatio() / a / a * hbarc * hbarc
                                * param->getg() * param->getg());
                        }

                        // 3 flavors
                        alphas =
                            4. * M_PI
                            / (9.
                               * log(pow(
                                   pow(muZero / 0.2, 2. / c)
                                       + pow(
                                           param->getRunWithThisFactorTimesQs()
                                               * Qs / 0.2,
                                           2. / c),
                                   c)));
                        gfactor = g * g / (4. * M_PI * alphas);
                        // run with the local (in transverse plane) coupling
                    } else
                        gfactor = 1.;

                    foutEps
                        << x << " " << y << " "
                        << hbarc * gfactor * abs(lat->cells[pos]->getEpsilon())
                        << endl;
                    // abs just to get rid of negative 10^(-17) numbers at edge
                }
                foutEps << endl;
            }
            foutEps.close();
        }

        if (it == itmax / 2 && param->getWriteOutputs() == 3) {
            IPG_PROFILE_SCOPE("output.epsilon_intermediate_text");
            stringstream streInt_name;
            streInt_name << "epsilonIntermediatePlot" << param->getEventId()
                         << ".dat";
            string eInt_name;
            eInt_name = streInt_name.str();

            ofstream foutEps2(eInt_name.c_str(), std::ios::out);
            for (int ix = 0; ix < N; ix++)  // loop over all positions
            {
                for (int iy = 0; iy < N; iy++) {
                    pos = ix * N + iy;
                    x = -L / 2. + a * ix;
                    y = -L / 2. + a * iy;

                    if (param->getRunningCoupling()) {
                        if (pos / N > 0 && pos / N < N - 1 && pos % N > 0
                            && pos % N < N - 1) {
                            g2mu2A = lat->cells[pos]->getg2mu2A();
                        } else
                            g2mu2A = 0;

                        if (pos / N > 0 && pos / N < N - 1 && pos % N > 0
                            && pos % N < N - 1) {
                            g2mu2B = lat->cells[pos]->getg2mu2B();
                        } else
                            g2mu2B = 0;

                        if (param->getRunWithQs() == 2) {
                            if (g2mu2A > g2mu2B)
                                Qs = sqrt(
                                    g2mu2A * param->getQsmuRatio()
                                    * param->getQsmuRatio() / a / a * hbarc
                                    * hbarc * param->getg() * param->getg());
                            else
                                Qs = sqrt(
                                    g2mu2B * param->getQsmuRatio()
                                    * param->getQsmuRatio() / a / a * hbarc
                                    * hbarc * param->getg() * param->getg());
                        } else if (param->getRunWithQs() == 0) {
                            if (g2mu2A < g2mu2B)
                                Qs = sqrt(
                                    g2mu2A * param->getQsmuRatio()
                                    * param->getQsmuRatio() / a / a * hbarc
                                    * hbarc * param->getg() * param->getg());
                            else
                                Qs = sqrt(
                                    g2mu2B * param->getQsmuRatio()
                                    * param->getQsmuRatio() / a / a * hbarc
                                    * hbarc * param->getg() * param->getg());
                        } else if (param->getRunWithQs() == 1) {
                            Qs = sqrt(
                                (g2mu2A + g2mu2B) / 2. * param->getQsmuRatio()
                                * param->getQsmuRatio() / a / a * hbarc * hbarc
                                * param->getg() * param->getg());
                        }

                        if (param->getRunWithLocalQs() == 1) {
                            // 3 flavors
                            alphas =
                                4. * M_PI
                                / (9.
                                   * log(pow(
                                       pow(muZero / 0.2, 2. / c)
                                           + pow(
                                               param->getRunWithThisFactorTimesQs()
                                                   * Qs / 0.2,
                                               2. / c),
                                       c)));
                            gfactor = g * g / (4. * M_PI * alphas);
                            // run with the local (in transverse plane) coupling
                        } else {
                            if (param->getRunWithQs() == 0)
                                alphas =
                                    4. * M_PI
                                    / (9.
                                       * log(pow(
                                           pow(muZero / 0.2, 2. / c)
                                               + pow(
                                                   param->getRunWithThisFactorTimesQs()
                                                       * param
                                                             ->getAverageQsmin()
                                                       / 0.2,
                                                   2. / c),
                                           c)));
                            else if (param->getRunWithQs() == 1)
                                alphas =
                                    4. * M_PI
                                    / (9.
                                       * log(pow(
                                           pow(muZero / 0.2, 2. / c)
                                               + pow(
                                                   param->getRunWithThisFactorTimesQs()
                                                       * param
                                                             ->getAverageQsAvg()
                                                       / 0.2,
                                                   2. / c),
                                           c)));
                            else if (param->getRunWithQs() == 2)
                                alphas =
                                    4. * M_PI
                                    / (9.
                                       * log(pow(
                                           pow(muZero / 0.2, 2. / c)
                                               + pow(
                                                   param->getRunWithThisFactorTimesQs()
                                                       * param->getAverageQs()
                                                       / 0.2,
                                                   2. / c),
                                           c)));

                            gfactor = g * g / (4. * M_PI * alphas);
                        }
                    } else
                        gfactor = 1.;

                    foutEps2
                        << x << " " << y << " "
                        << hbarc * gfactor * abs(lat->cells[pos]->getEpsilon())
                        << endl;
                    // abs just to get rid of negative 10^(-17) numbers at edge
                }
                foutEps2 << endl;
            }
            foutEps2.close();
        }

        if (it == itmax) {
            checkGaussLaw(lat, param);
        }

        int success = 1;
        if (param->getComputeGluonMultiplicity()) {
            if (it == itmax) {
                eccentricity(lat, param, it, 0.0, 0);
                // eccentricity(lat, param, it, 0.1, 0);
                // eccentricity(lat, param, it, 1., 0);
                // eccentricity(lat, param, it, 10., 0);

                success = multiplicity(lat, group, param, it);
            }
        }

        if (success == 0) break;
    }
}

void Evolution::tmunu(Lattice *lat, Parameters *param, int it) {
    IPG_PROFILE_SCOPE("observables.Tmunu");
    int N = param->getSize();
    double L = param->getL();
    double a = L / N;  // lattice spacing in fm
    double g = param->getg();
    double dtau = param->getdtau();
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

void Evolution::u(Lattice *lat, Parameters *param, int it, bool finalFlag) {
    IPG_PROFILE_SCOPE("observables.flow_velocity");
    MyEigen myeigen;
    myeigen.flowVelocity4D(lat, param, it, finalFlag);
}

void Evolution::eccentricity(
    Lattice *lat, Parameters *param, int it, double cutoff, int doAniso) {
    IPG_PROFILE_SCOPE("observables.eccentricity");
    stringstream strecc_name;
    strecc_name << "eccentricities" << param->getEventId() << ".dat";
    string ecc_name;
    ecc_name = strecc_name.str();

    // cutoff on energy density is 'cutoff' times Lambda_QCD^4
    int N = param->getSize();
    int pos;
    double rA, phiA, x, y;
    double L = param->getL();
    double a = L / N;  // lattice spacing in fm
    double eccentricity1, eccentricity2, eccentricity3, eccentricity4,
        eccentricity5, eccentricity6;
    double avcos, avsin, avcos1, avsin1, avcos3, avsin3, avrSq, avxSq, avySq,
        avr1, avr3, avcos4, avsin4, avr4, avcos5, avsin5, avr5, avcos6, avsin6,
        avr6;
    double Rbar;
    double Psi1, Psi2, Psi3, Psi4, Psi5, Psi6;
    double maxEps = 0;
    double g = param->getg();

    double g2mu2A, g2mu2B, gfactor, alphas = 0., Qs = 0.;
    double c = param->getc();
    double muZero = param->getMuZero();

    double weight;

    double area = 0.;
    double avgeden = 0.;
    int sum = 0;

    avrSq = 0.;
    avr3 = 0.;

    double avx = 0.;
    double avy = 0.;
    double toteps = 0.;
    int xshift;
    int yshift;
    double maxX = 0.;
    double maxY = 0.;

    double smallestX = 0.;
    double smallestY = 0.;
    double avgQs2AQs2B = 0.;

    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            pos = ix * N + iy;
            maxEps = std::max(lat->cells[pos]->getEpsilon(), maxEps);
        }
    }

    // first shift to the center
    for (int ix = 0; ix < N; ix++) {
        x = -L / 2. + a * ix;
        for (int iy = 0; iy < N; iy++) {
            y = -L / 2. + a * iy;
            pos = ix * N + iy;

            gfactor = computeRunningCouplingGfactor(
                lat, param, pos, N, a, g, c, muZero);

            if (lat->cells[pos]->getEpsilon() * gfactor
                < cutoff)  // this is 1/fm^4, so Lambda_QCD^{-4} (because
                           // \Lambda_QCD is roughly 1/fm)
            {
                weight = 0.;
            } else {
                weight =
                    (lat->cells[pos]->getEpsilon() * lat->cells[pos]->getutau()
                     * gfactor);
                area += a * a;
                sum += 1;
                avgeden += lat->cells[pos]->getEpsilon() * hbarc
                           * gfactor;  // GeV/fm^3
                g2mu2A = lat->cells[pos]->getg2mu2A();
                g2mu2B = lat->cells[pos]->getg2mu2B();
                avgQs2AQs2B += g2mu2A * param->getQsmuRatio()
                               * param->getQsmuRatio() * g2mu2B
                               * param->getQsmuRatioB() * param->getQsmuRatioB()
                               / a / a / a / a;
            }
            avx += x * weight;
            avy += y * weight;
            toteps += weight;
        }
    }

    avx /= toteps;
    avy /= toteps;
    avgeden /= double(sum);
    avgQs2AQs2B /= double(sum);
    param->setArea(area);

    xshift = static_cast<int>(floor(avx / a + 0.00000000001));
    yshift = static_cast<int>(floor(avy / a + 0.00000000001));

    avcos1 = 0.;
    avsin1 = 0.;
    avcos = 0.;
    avsin = 0.;
    avcos3 = 0.;
    avsin3 = 0.;
    avcos4 = 0.;
    avsin4 = 0.;
    avcos5 = 0.;
    avsin5 = 0.;
    avcos6 = 0.;
    avsin6 = 0.;
    avr1 = 0.;
    avrSq = 0.;
    avxSq = 0.;
    avySq = 0.;
    avr3 = 0.;
    avr4 = 0.;
    avr5 = 0.;
    avr6 = 0.;

    for (int ix = 2; ix < N - 2; ix++) {
        x = -L / 2. + a * ix - avx;
        for (int iy = 2; iy < N - 2; iy++) {
            pos = ix * N + iy;
            y = -L / 2. + a * iy - avy;
            if (x >= 0) {
                phiA = atan(y / x);
                if (x == 0) {
                    if (y >= 0)
                        phiA = M_PI / 2.;
                    else if (y < 0)
                        phiA = 3. * M_PI / 2.;
                }
            } else {
                phiA = atan(y / x) + M_PI;
            }

            gfactor = computeRunningCouplingGfactor(
                lat, param, pos, N, a, g, c, muZero);

            if (lat->cells[pos]->getEpsilon() * gfactor
                < cutoff)  // this is 1/fm^4, so Lambda_QCD^{-4}
            {
                weight = 0.;
            } else {
                weight =
                    (lat->cells[pos]->getEpsilon() * lat->cells[pos]->getutau()
                     * gfactor);
            }

            rA = sqrt(x * x + y * y);
            avr1 += rA * rA * rA * (weight);
            avrSq += rA * rA * (weight);  // compute average r^2
            avr3 += rA * rA * rA * (weight);
            avr4 += rA * rA * rA * rA * (weight);
            avr5 += rA * rA * rA * rA * rA * (weight);
            avr6 += rA * rA * rA * rA * rA * rA * (weight);

            avcos1 += rA * rA * rA * cos(phiA) * (weight);
            avsin1 += rA * rA * rA * sin(phiA) * (weight);
            avcos += rA * rA * cos(2. * phiA) * (weight);
            avsin += rA * rA * sin(2. * phiA) * (weight);
            avcos3 += rA * rA * rA * cos(3. * phiA) * (weight);
            avsin3 += rA * rA * rA * sin(3. * phiA) * (weight);
            avcos4 += rA * rA * rA * rA * cos(4. * phiA) * (weight);
            avsin4 += rA * rA * rA * rA * sin(4. * phiA) * (weight);
            avcos5 += rA * rA * rA * rA * rA * cos(5. * phiA) * (weight);
            avsin5 += rA * rA * rA * rA * rA * sin(5. * phiA) * (weight);
            avcos6 += rA * rA * rA * rA * rA * rA * cos(6. * phiA) * (weight);
            avsin6 += rA * rA * rA * rA * rA * rA * sin(6. * phiA) * (weight);

            if (weight > cutoff && iy == N / 2 + yshift) {
                maxX = x;
            }
            if (weight > cutoff && ix == N / 2 + xshift) {
                maxY = y;
            }

            if (weight < cutoff && iy == N / 2 + yshift && ix > N / 2 + xshift
                && smallestX == 0) {
                smallestX = x;
            }
            if (weight < cutoff && ix == N / 2 + xshift && iy > N / 2 + yshift
                && smallestY == 0) {
                smallestY = y;
            }
        }
    }

    // compute and print eccentricity and angles:
    Psi1 = (atan(avsin1 / avcos1) + M_PI) / 1.;
    Psi2 = (atan(avsin / avcos) + M_PI) / 2.;
    Psi3 = (atan(avsin3 / avcos3) + M_PI) / 3.;
    Psi4 = (atan(avsin4 / avcos4) + M_PI) / 4.;
    Psi5 = (atan(avsin5 / avcos5) + M_PI) / 5.;
    Psi6 = (atan(avsin6 / avcos6) + M_PI) / 6.;
    eccentricity1 = sqrt(avcos1 * avcos1 + avsin1 * avsin1) / avr1;
    eccentricity2 = sqrt(avcos * avcos + avsin * avsin) / avrSq;
    eccentricity3 = sqrt(avcos3 * avcos3 + avsin3 * avsin3) / avr3;
    eccentricity4 = sqrt(avcos4 * avcos4 + avsin4 * avsin4) / avr4;
    eccentricity5 = sqrt(avcos5 * avcos5 + avsin5 * avsin5) / avr5;
    eccentricity6 = sqrt(avcos6 * avcos6 + avsin6 * avsin6) / avr6;

    double avx2 = avx;
    double avy2 = avy;
    avx = 0.;
    avy = 0.;
    toteps = 0.;

    for (int ix = 0; ix < N; ix++) {
        x = -L / 2. + a * ix - avx2;
        for (int iy = 0; iy < N; iy++) {
            pos = ix * N + iy;

            gfactor = computeRunningCouplingGfactor(
                lat, param, pos, N, a, g, c, muZero);

            if (lat->cells[pos]->getEpsilon() * gfactor
                < cutoff)  // this is 1/fm^4, so Lambda_QCD^{-4}
            {
                weight = 0.;
            } else {
                weight =
                    (lat->cells[pos]->getEpsilon() * lat->cells[pos]->getutau()
                     * gfactor);
            }

            y = -L / 2. + a * iy - avy2;
            avx += x * weight;
            avy += y * weight;
            avxSq += x * x * weight;
            avySq += y * y * weight;
            toteps += weight;
        }
    }
    avx /= toteps;
    avy /= toteps;
    avxSq /= toteps;
    avySq /= toteps;
    avrSq /= toteps;
    Rbar = 1. / sqrt(1. / avxSq + 1. / avySq);
    if (it == 1) param->setPsi(Psi2);

    if (doAniso == 0) {
        ofstream foutEcc(ecc_name.c_str(), std::ios::app);
        foutEcc << it * a * param->getdtau() << " " << eccentricity1 << " "
                << Psi1 << " " << eccentricity2 << " " << Psi2 << " "
                << eccentricity3 << " " << Psi3 << " " << eccentricity4 << " "
                << Psi4 << " " << eccentricity5 << " " << Psi5 << " "
                << eccentricity6 << " " << Psi6 << " " << cutoff << " "
                << sqrt(avrSq) << " " << maxX << " " << maxY << " "
                << param->getb() << " " << param->getTpp() << " "
                << param->getArea() << " " << Rbar << " " << avgeden << " "
                << avgQs2AQs2B * hbarc << endl;
        foutEcc.close();
    }

    if (doAniso == 1) {
        stringstream straniso_name;
        straniso_name << "anisotropy" << param->getEventId() << ".dat";
        string aniso_name;
        aniso_name = straniso_name.str();

        ofstream foutAniso(aniso_name.c_str(), std::ios::app);

        double TxxRot, TyyRot;
        double ux, uy, PsiU;
        double num = 0., den = 0.;
        double unum = 0., uden = 0.;
        double num2 = 0., den2 = 0.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = (ix)*N + (iy);
                ux = lat->cells[pos]->getux();
                uy = lat->cells[pos]->getuy();
                unum += sqrt(ux * ux + uy * uy) * sin(2. * atan2(uy, ux));
                uden += sqrt(ux * ux + uy * uy) * cos(2. * atan2(uy, ux));
            }
        }

        PsiU = atan2(unum, uden) / 2.;

        num = 0.;
        den = 0.;
        num2 = 0.;
        den2 = 0.;
        double Psi = PsiU;  // param->getPsi();//-Pi/2.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = (ix)*N + (iy);

                TxxRot = cos(Psi)
                             * (cos(Psi) * lat->cells[pos]->getTxx()
                                - sin(Psi) * lat->cells[pos]->getTxy())
                         - sin(Psi)
                               * (cos(Psi) * lat->cells[pos]->getTxy()
                                  - sin(Psi) * lat->cells[pos]->getTyy());
                TyyRot = sin(Psi)
                             * (sin(Psi) * lat->cells[pos]->getTxx()
                                + cos(Psi) * lat->cells[pos]->getTxy())
                         + cos(Psi)
                               * (sin(Psi) * lat->cells[pos]->getTxy()
                                  + cos(Psi) * lat->cells[pos]->getTyy());

                num2 += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
                den2 += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();

                num += TxxRot - TyyRot;
                den += TxxRot + TyyRot;
            }
        }

        foutAniso << "Psi2=" << Psi2 << ", cos(Psi2)=" << cos(Psi2)
                  << ", sin(Psi2)=" << sin(Psi2) << endl;
        foutAniso << "PsiU=" << PsiU << ", cos(PsiU)=" << cos(PsiU)
                  << ", sin(PsiU)=" << sin(PsiU) << endl;
        foutAniso << it * a * param->getdtau() << " " << num / den << " "
                  << num2 / den2 << " angle=" << Psi << endl;

        num = 0.;
        den = 0.;
        num2 = 0.;
        den2 = 0.;
        Psi = PsiU + M_PI / 8.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = (ix)*N + (iy);

                TxxRot = cos(Psi)
                             * (cos(Psi) * lat->cells[pos]->getTxx()
                                - sin(Psi) * lat->cells[pos]->getTxy())
                         - sin(Psi)
                               * (cos(Psi) * lat->cells[pos]->getTxy()
                                  - sin(Psi) * lat->cells[pos]->getTyy());
                TyyRot = sin(Psi)
                             * (sin(Psi) * lat->cells[pos]->getTxx()
                                + cos(Psi) * lat->cells[pos]->getTxy())
                         + cos(Psi)
                               * (sin(Psi) * lat->cells[pos]->getTxy()
                                  + cos(Psi) * lat->cells[pos]->getTyy());

                num2 += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
                den2 += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();

                num += TxxRot - TyyRot;
                den += TxxRot + TyyRot;
            }
        }

        foutAniso << it * a * param->getdtau() << " " << num / den << " "
                  << num2 / den2 << " angle=" << Psi << endl;

        num = 0.;
        den = 0.;
        num2 = 0.;
        den2 = 0.;
        Psi = PsiU + M_PI / 4.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = (ix)*N + (iy);

                TxxRot = cos(Psi)
                             * (cos(Psi) * lat->cells[pos]->getTxx()
                                - sin(Psi) * lat->cells[pos]->getTxy())
                         - sin(Psi)
                               * (cos(Psi) * lat->cells[pos]->getTxy()
                                  - sin(Psi) * lat->cells[pos]->getTyy());
                TyyRot = sin(Psi)
                             * (sin(Psi) * lat->cells[pos]->getTxx()
                                + cos(Psi) * lat->cells[pos]->getTxy())
                         + cos(Psi)
                               * (sin(Psi) * lat->cells[pos]->getTxy()
                                  + cos(Psi) * lat->cells[pos]->getTyy());

                num2 += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
                den2 += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();

                num += TxxRot - TyyRot;
                den += TxxRot + TyyRot;
            }
        }

        foutAniso << it * a * param->getdtau() << " " << num / den << " "
                  << num2 / den2 << " angle=" << Psi << endl;

        num = 0.;
        den = 0.;
        num2 = 0.;
        den2 = 0.;
        Psi = PsiU + 3. * M_PI / 8.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = (ix)*N + (iy);

                TxxRot = cos(Psi)
                             * (cos(Psi) * lat->cells[pos]->getTxx()
                                - sin(Psi) * lat->cells[pos]->getTxy())
                         - sin(Psi)
                               * (cos(Psi) * lat->cells[pos]->getTxy()
                                  - sin(Psi) * lat->cells[pos]->getTyy());
                TyyRot = sin(Psi)
                             * (sin(Psi) * lat->cells[pos]->getTxx()
                                + cos(Psi) * lat->cells[pos]->getTxy())
                         + cos(Psi)
                               * (sin(Psi) * lat->cells[pos]->getTxy()
                                  + cos(Psi) * lat->cells[pos]->getTyy());

                num2 += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
                den2 += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();

                num += TxxRot - TyyRot;
                den += TxxRot + TyyRot;
            }
        }

        foutAniso << it * a * param->getdtau() << " " << num / den << " "
                  << num2 / den2 << " angle=" << Psi << endl;

        num = 0.;
        den = 0.;
        num2 = 0.;
        den2 = 0.;
        Psi = PsiU + M_PI / 2.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = (ix)*N + (iy);

                TxxRot = cos(Psi)
                             * (cos(Psi) * lat->cells[pos]->getTxx()
                                - sin(Psi) * lat->cells[pos]->getTxy())
                         - sin(Psi)
                               * (cos(Psi) * lat->cells[pos]->getTxy()
                                  - sin(Psi) * lat->cells[pos]->getTyy());
                TyyRot = sin(Psi)
                             * (sin(Psi) * lat->cells[pos]->getTxx()
                                + cos(Psi) * lat->cells[pos]->getTxy())
                         + cos(Psi)
                               * (sin(Psi) * lat->cells[pos]->getTxy()
                                  + cos(Psi) * lat->cells[pos]->getTyy());

                num2 += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
                den2 += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();

                num += TxxRot - TyyRot;
                den += TxxRot + TyyRot;
            }
        }

        foutAniso << it * a * param->getdtau() << " " << num / den << " "
                  << num2 / den2 << " angle=" << Psi << endl;

        num = 0.;
        den = 0.;
        num2 = 0.;
        den2 = 0.;
        Psi = PsiU + 5. * M_PI / 8.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = (ix)*N + (iy);

                TxxRot = cos(Psi)
                             * (cos(Psi) * lat->cells[pos]->getTxx()
                                - sin(Psi) * lat->cells[pos]->getTxy())
                         - sin(Psi)
                               * (cos(Psi) * lat->cells[pos]->getTxy()
                                  - sin(Psi) * lat->cells[pos]->getTyy());
                TyyRot = sin(Psi)
                             * (sin(Psi) * lat->cells[pos]->getTxx()
                                + cos(Psi) * lat->cells[pos]->getTxy())
                         + cos(Psi)
                               * (sin(Psi) * lat->cells[pos]->getTxy()
                                  + cos(Psi) * lat->cells[pos]->getTyy());

                num2 += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
                den2 += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();

                num += TxxRot - TyyRot;
                den += TxxRot + TyyRot;
            }
        }

        foutAniso << it * a * param->getdtau() << " " << num / den << " "
                  << num2 / den2 << " angle=" << Psi << endl;

        num = 0.;
        den = 0.;
        num2 = 0.;
        den2 = 0.;
        Psi = PsiU + 3. * M_PI / 4.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = (ix)*N + (iy);

                TxxRot = cos(Psi)
                             * (cos(Psi) * lat->cells[pos]->getTxx()
                                - sin(Psi) * lat->cells[pos]->getTxy())
                         - sin(Psi)
                               * (cos(Psi) * lat->cells[pos]->getTxy()
                                  - sin(Psi) * lat->cells[pos]->getTyy());
                TyyRot = sin(Psi)
                             * (sin(Psi) * lat->cells[pos]->getTxx()
                                + cos(Psi) * lat->cells[pos]->getTxy())
                         + cos(Psi)
                               * (sin(Psi) * lat->cells[pos]->getTxy()
                                  + cos(Psi) * lat->cells[pos]->getTyy());

                num2 += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
                den2 += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();

                num += TxxRot - TyyRot;
                den += TxxRot + TyyRot;
            }
        }

        foutAniso << it * a * param->getdtau() << " " << num / den << " "
                  << num2 / den2 << " angle=" << Psi << endl;

        num = 0.;
        den = 0.;
        num2 = 0.;
        den2 = 0.;
        Psi = PsiU + 7. * M_PI / 8.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = (ix)*N + (iy);

                TxxRot = cos(Psi)
                             * (cos(Psi) * lat->cells[pos]->getTxx()
                                - sin(Psi) * lat->cells[pos]->getTxy())
                         - sin(Psi)
                               * (cos(Psi) * lat->cells[pos]->getTxy()
                                  - sin(Psi) * lat->cells[pos]->getTyy());
                TyyRot = sin(Psi)
                             * (sin(Psi) * lat->cells[pos]->getTxx()
                                + cos(Psi) * lat->cells[pos]->getTxy())
                         + cos(Psi)
                               * (sin(Psi) * lat->cells[pos]->getTxy()
                                  + cos(Psi) * lat->cells[pos]->getTyy());

                num2 += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
                den2 += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();

                num += TxxRot - TyyRot;
                den += TxxRot + TyyRot;
            }
        }

        foutAniso << it * a * param->getdtau() << " " << num / den << " "
                  << num2 / den2 << " angle=" << Psi << endl;

        num = 0.;
        den = 0.;
        num2 = 0.;
        den2 = 0.;
        Psi = PsiU + M_PI;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = (ix)*N + (iy);

                TxxRot = cos(Psi)
                             * (cos(Psi) * lat->cells[pos]->getTxx()
                                - sin(Psi) * lat->cells[pos]->getTxy())
                         - sin(Psi)
                               * (cos(Psi) * lat->cells[pos]->getTxy()
                                  - sin(Psi) * lat->cells[pos]->getTyy());
                TyyRot = sin(Psi)
                             * (sin(Psi) * lat->cells[pos]->getTxx()
                                + cos(Psi) * lat->cells[pos]->getTxy())
                         + cos(Psi)
                               * (sin(Psi) * lat->cells[pos]->getTxy()
                                  + cos(Psi) * lat->cells[pos]->getTyy());

                num2 += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
                den2 += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();

                num += TxxRot - TyyRot;
                den += TxxRot + TyyRot;
            }
        }

        foutAniso << it * a * param->getdtau() << " " << num / den << " "
                  << num2 / den2 << " angle=" << Psi << endl;

        num = 0.;
        den = 0.;
        num2 = 0.;
        den2 = 0.;
        Psi = PsiU + 9. * M_PI / 8.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = (ix)*N + (iy);

                TxxRot = cos(Psi)
                             * (cos(Psi) * lat->cells[pos]->getTxx()
                                - sin(Psi) * lat->cells[pos]->getTxy())
                         - sin(Psi)
                               * (cos(Psi) * lat->cells[pos]->getTxy()
                                  - sin(Psi) * lat->cells[pos]->getTyy());
                TyyRot = sin(Psi)
                             * (sin(Psi) * lat->cells[pos]->getTxx()
                                + cos(Psi) * lat->cells[pos]->getTxy())
                         + cos(Psi)
                               * (sin(Psi) * lat->cells[pos]->getTxy()
                                  + cos(Psi) * lat->cells[pos]->getTyy());

                num2 += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
                den2 += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();

                num += TxxRot - TyyRot;
                den += TxxRot + TyyRot;
            }
        }

        foutAniso << it * a * param->getdtau() << " " << num / den << " "
                  << num2 / den2 << " angle=" << Psi << endl;

        foutAniso.close();
    }
}

void Evolution::readNkt(Parameters *param) {
    messager_ << "[Evolution::readNkt]: Reading n(k_T) from file ";
    messager_.flush("info");
    string Npart, dummy;
    string kt, nkt, Tpp, b;
    double dkt = 0.;
    double dNdeta = 0.;

    // open file

    ifstream fin;
    stringstream strmult_name;
    strmult_name << "multiplicity" << param->getEventId() << ".dat";
    string mult_name;
    mult_name = strmult_name.str();
    fin.open(mult_name.c_str());
    messager_ << "[Evolution::readNkt]: File " << mult_name.c_str() << " ... ";
    messager_.flush("info");

    // open file

    ifstream fin2;
    stringstream strmult_name2;
    strmult_name2 << "NpartdNdy" << param->getEventId() << ".dat";
    string mult_name2;
    mult_name2 = strmult_name2.str();
    fin2.open(mult_name2.c_str());
    messager_ << "[Evolution::readNkt]: File " << mult_name2.c_str() << " ... ";
    messager_.flush("info");

    // read file

    if (fin) {
        for (int ikt = 0; ikt < 100; ikt++) {
            if (!fin.eof()) {
                fin >> dummy;
                fin >> kt;
                fin >> nkt;
                nIn_[ikt] = atof(nkt.c_str());
                fin >> dummy >> Tpp >> b >> Npart;
                if (ikt == 0) dkt = atof(kt.c_str());
                if (ikt == 1) dkt = dkt - atof(kt.c_str());
            }
            messager_ << "[Evolution::readNkt]: " << nIn_[ikt];
            messager_.flush("info");
        }
        fin.close();
        messager_ << "[Evolution::readNkt]:  done.";
        messager_.flush("info");
    } else {
        messager_ << "[Evolution::readNkt]: File " << mult_name.c_str()
                  << " does not exist. Exiting.";
        messager_.flush("error");
        exit(1);
    }

    if (fin2) {
        if (!fin2.eof()) {
            fin2 >> Npart;
            fin2 >> nkt;
            fin2 >> Tpp;
            fin2 >> b;

            dNdeta = atof(nkt.c_str());
        }
        fin2.close();
        messager_ << "[Evolution::readNkt]:  done.";
        messager_.flush("info");
    } else {
        messager_ << "[Evolution::readNkt]: File " << mult_name2.c_str()
                  << " does not exist. Exiting.";
        messager_.flush("error");
        exit(1);
    }

    double m, P;
    m = param->getJacobianm();                                // in GeV
    P = 0.13 + 0.32 * pow(param->getRoots() / 1000., 0.115);  // in GeV
    double dNdeta2;
    dNdeta2 = 0.;

    for (int ik = 0; ik < 100; ik++) {
        if (param->getUsePseudoRapidity() == 0) {
            dNdeta2 += nIn_[ik] * (ik + 0.5) * dkt * dkt * 2.
                       * M_PI;  // integrate, gives a ik*dkt*2pi*dkt
        } else {
            dNdeta2 +=
                nIn_[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI
                * cosh(param->getRapidity())
                / (sqrt(
                    pow(cosh(param->getRapidity()), 2.)
                    + m * m / (((ik + 0.5) * dkt) * ((ik + 0.5) * dkt))));
        }
    }

    dNdeta *= cosh(param->getRapidity())
              / (sqrt(pow(cosh(param->getRapidity()), 2.) + m * m / P / P));

    ofstream foutNN("NpartdNdy-mod.dat", std::ios::out);
    foutNN << Npart << " " << dNdeta << " " << dNdeta2 << " "
           << atof(Tpp.c_str()) << " " << atof(b.c_str()) << endl;
    foutNN.close();

    exit(1);
}

int Evolution::multiplicity(
    Lattice *lat, Group *group, Parameters *param, int it) {
    IPG_PROFILE_SCOPE("observables.gluon_multiplicity");
    int N = param->getSize();
    int npos, pos;
    double L = param->getL();
    double a = L / N;  // lattice spacing in fm
    double kx, ky, kt2, omega2;
    double g = param->getg();
    int nn[2];
    nn[0] = N;
    nn[1] = N;
    double dtau = param->getdtau();
    double nkt;
    const int bins = 100;
    double n[bins];   // k_T array
    double E[bins];   // k_T array
    double n2[bins];  // k_T array
    int counter[bins];
    double dkt = 2.83 / static_cast<double>(bins);
    double dNdeta = 0.;
    double dNdeta2 = 0.;
    double dNdetaCut = 0.;
    double dNdetaCut2 = 0.;
    double dEdetaCut = 0.;
    double dEdetaCut2 = 0.;
    double dEdeta = 0.;
    double dEdeta2 = 0.;

    stringstream strNpartdNdy_name;
    strNpartdNdy_name << "NpartdNdy-t" << it * dtau * a << "-"
                      << param->getEventId() << ".dat";
    string NpartdNdy_name;
    NpartdNdy_name = strNpartdNdy_name.str();
    messager_ << "[Evolution::multiplicity]: Measuring multiplicity ... ";
    messager_.flush("info");

    // fix transverse Coulomb gauge
    GaugeFix gaugefix;

    double maxtime;
    if (param->getInverseQsForMaxTime() == 1) {
        maxtime = 1. / param->getAverageQs() * hbarc;
        messager_ << "[Evolution::multiplicity]: maximal evolution time = "
                  << maxtime << " fm";
        messager_.flush("info");
    } else {
        maxtime = param->getMaxtime();  // maxtime is in fm
    }

    int itmax = static_cast<int>(floor(maxtime / (a * dtau) + 1e-10));

    double multiplicityPhaseStart = ipg::wallSeconds();
    gaugefix.fftChi(fft_, lat, group, param, 4000);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.gauge_fix", multiplicityPhaseStart);
    // gauge is fixed

    // E1Storage owns the N*N scratch matrices; E1 is a pointer-per-cell
    // view over it for FFT::fftn's T** interface.
    std::vector<Matrix> E1Storage(N * N, Matrix(0.));
    std::vector<Matrix *> E1(N * N);
    for (int i = 0; i < N * N; i++) {
        E1[i] = &E1Storage[i];
    }
    addPhaseAndRestart(
        "observables.gluon_multiplicity.allocate", multiplicityPhaseStart);

    double c = param->getc();
    double muZero = param->getMuZero();

    prepareSpectrumField(lat, param, N, a, g, c, muZero, lat->U, E1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.prepare_E1", multiplicityPhaseStart);

    // do Fourier transforms
    fft_->fftn(E1.data(), E1.data(), nn, 1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.fft_E1", multiplicityPhaseStart);

    for (int ik = 0; ik < bins; ik++) {
        n[ik] = 0.;
        E[ik] = 0.;
        n2[ik] = 0.;
        counter[ik] = 0;
    }

    const int hbins = 2000;
    double Nhgsl[hbins + 1];
    double Ng;

    addPhaseAndRestart(
        "observables.gluon_multiplicity.setup_bins", multiplicityPhaseStart);

    accumulateGluonSpectrum(
        param, N, it, dtau, g, a, c, muZero, dkt, bins, E1, true, true,
        dNdeta, dEdeta, n, E, n2, counter);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.spectrum_E1", multiplicityPhaseStart);

    /// -------- 2 ---------

    prepareSpectrumField(lat, param, N, a, g, c, muZero, lat->U2, E1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.prepare_E2", multiplicityPhaseStart);

    fft_->fftn(E1.data(), E1.data(), nn, 1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.fft_E2", multiplicityPhaseStart);

    accumulateGluonSpectrum(
        param, N, it, dtau, g, a, c, muZero, dkt, bins, E1, true, false,
        dNdeta, dEdeta, n, E, n2, counter);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.spectrum_E2", multiplicityPhaseStart);

    /// ------3 --------

    prepareSpectrumField(lat, param, N, a, g, c, muZero, lat->Ux2, E1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.prepare_pi", multiplicityPhaseStart);

    // do Fourier transforms
    fft_->fftn(E1.data(), E1.data(), nn, 1);
    addPhaseAndRestart(
        "observables.gluon_multiplicity.fft_pi", multiplicityPhaseStart);

    accumulateGluonSpectrum(
        param, N, it, dtau, g, a, c, muZero, dkt, bins, E1, false, false,
        dNdeta, dEdeta, n, E, n2, counter);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.spectrum_pi", multiplicityPhaseStart);

    double m, P;
    m = param->getJacobianm();                                // in GeV
    P = 0.13 + 0.32 * pow(param->getRoots() / 1000., 0.115);  // in GeV

    for (int ik = 0; ik < bins; ik++) {
        if (counter[ik] > 0) {
            n[ik] = n[ik] / static_cast<double>(counter[ik]);
            E[ik] = E[ik] / static_cast<double>(counter[ik]);
            if (param->getUsePseudoRapidity() == 0) {
                dNdeta2 += n[ik] * (ik + 0.5) * dkt * dkt * 2.
                           * M_PI;  // integrate, gives a ik*dkt*2pi*dkt
                dEdeta2 += E[ik] * (ik + 0.5) * dkt * dkt * 2.
                           * M_PI;  // integrate, gives a ik*dkt*2pi*dkt
                if (ik * dkt / a * hbarc > 3.)  //
                {
                    dNdetaCut += n[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI;
                    dEdetaCut += E[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI;
                }
                if (ik * dkt / a * hbarc > 6.)  // large cut
                {
                    dNdetaCut2 += n[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI;
                    dEdetaCut2 += E[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI;
                }
            } else {
                dNdeta2 += n[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI
                           * cosh(param->getRapidity())
                           / (sqrt(
                               pow(cosh(param->getRapidity()), 2.)
                               + m * m
                                     / (((ik + 0.5) * dkt / a * hbarc)
                                        * ((ik + 0.5) * dkt / a
                                           * hbarc))));  // integrate, gives a
                                                         // ik*dkt*2pi*dkt
                dEdeta2 += E[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI
                           * cosh(param->getRapidity())
                           / (sqrt(
                               pow(cosh(param->getRapidity()), 2.)
                               + m * m
                                     / (((ik + 0.5) * dkt / a * hbarc)
                                        * ((ik + 0.5) * dkt / a
                                           * hbarc))));  // integrate, gives a
                                                         // ik*dkt*2pi*dkt

                if (ik * dkt / a * hbarc > 3.)  //
                {
                    dNdetaCut +=
                        n[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI
                        * cosh(param->getRapidity())
                        / (sqrt(
                            pow(cosh(param->getRapidity()), 2.)
                            + m * m
                                  / (((ik + 0.5) * dkt / a * hbarc)
                                     * ((ik + 0.5) * dkt / a * hbarc))));
                    dEdetaCut +=
                        E[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI
                        * cosh(param->getRapidity())
                        / (sqrt(
                            pow(cosh(param->getRapidity()), 2.)
                            + m * m
                                  / (((ik + 0.5) * dkt / a * hbarc)
                                     * ((ik + 0.5) * dkt / a * hbarc))));
                }
                if (ik * dkt / a * hbarc > 6.)  // large cut
                {
                    dNdetaCut2 +=
                        n[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI
                        * cosh(param->getRapidity())
                        / (sqrt(
                            pow(cosh(param->getRapidity()), 2.)
                            + m * m
                                  / (((ik + 0.5) * dkt / a * hbarc)
                                     * ((ik + 0.5) * dkt / a * hbarc))));
                    dEdetaCut2 +=
                        E[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI
                        * cosh(param->getRapidity())
                        / (sqrt(
                            pow(cosh(param->getRapidity()), 2.)
                            + m * m
                                  / (((ik + 0.5) * dkt / a * hbarc)
                                     * ((ik + 0.5) * dkt / a * hbarc))));
                }
            }
        }
    }

    addPhaseAndRestart(
        "observables.gluon_multiplicity.bin_postprocess",
        multiplicityPhaseStart);

    // compute hadrons using fragmentation function
    if (it == itmax && param->getWriteOutputs() == 3) {
        const double hadronizationStart = ipg::wallSeconds();
        messager_ << "[Evolution::multiplicity]:  Hadronizing ... ";
        messager_.flush("info");
        double z, frac;
        double mypt, kt;
        int ik;
        const int steps = 6000;
        double dz = 0.95 / static_cast<double>(steps);
        double zValues[steps + 1];
        double zintegrand[steps + 1];
        gsl_interp_accel *zacc = gsl_interp_accel_alloc();
        gsl_spline *zspline = gsl_spline_alloc(gsl_interp_cspline, steps + 1);

        for (int ih = 0; ih <= hbins; ih++) {
            mypt = ih * (20. / static_cast<double>(hbins));

            for (int iz = 0; iz <= steps; iz++) {
                z = 0.05 + iz * dz;
                zValues[iz] = z;

                kt = mypt / z;

                ik = static_cast<int>(
                    floor(kt * a / hbarc / dkt - 0.5 + 0.00000001));

                frac = (kt - (ik + 0.5) * dkt / a * hbarc) / (dkt / a * hbarc);

                if (ik + 1 < bins && ik >= 0)
                    Ng = ((1. - frac) * n[ik] + frac * n[ik + 1]) * a / hbarc
                         * a / hbarc;  // to make dN/d^2k_T fo k_T in GeV
                else
                    Ng = 0.;

                if (param->getUsePseudoRapidity() == 0) {
                    zintegrand[iz] = 1. / (z * z) * Ng * kkp(7, 1, z, kt);
                } else {
                    zintegrand[iz] =
                        1. / (z * z) * Ng * 2.
                        * (kkp(1, 1, z, kt) * cosh(param->getRapidity())
                               / (sqrt(
                                   pow(cosh(param->getRapidity()), 2.)
                                   + m_pion * m_pion / (mypt * mypt)))
                           + kkp(2, 1, z, kt) * cosh(param->getRapidity())
                                 / (sqrt(
                                     pow(cosh(param->getRapidity()), 2.)
                                     + m_kaon * m_kaon / (mypt * mypt)))
                           + kkp(4, 1, z, kt) * cosh(param->getRapidity())
                                 / (sqrt(
                                     pow(cosh(param->getRapidity()), 2.)
                                     + m_proton * m_proton / (mypt * mypt))));
                }
            }

            zValues[steps] = 1.;  // set exactly 1

            gsl_spline_init(zspline, zValues, zintegrand, steps + 1);
            Nhgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);
        }

        gsl_spline_free(zspline);
        gsl_interp_accel_free(zacc);

        stringstream strmultHad_name;
        strmultHad_name << "multiplicityHadrons" << param->getEventId()
                        << ".dat";
        string multHad_name;
        multHad_name = strmultHad_name.str();

        ofstream foutdNdpt(multHad_name.c_str(), std::ios::out);
        for (int ih = 0; ih <= hbins; ih++) {
            if (ih % 10 == 0)
                foutdNdpt << ih * 20. / static_cast<double>(hbins) << " "
                          << Nhgsl[ih] << " " << 0. << " " << 0. << " "
                          << param->getTpp() << " " << param->getb()
                          << endl;  // leaving out the L and H ones for now
        }
        foutdNdpt.close();

        messager_ << "[Evolution::multiplicity]:  done.";
        messager_.flush("info");

        ipg::Profiler::instance().add(
            "observables.gluon_multiplicity.hadronization",
            ipg::wallSeconds() - hadronizationStart);
        multiplicityPhaseStart = ipg::wallSeconds();
    }

    if (param->getUsePseudoRapidity() == 0 && param->getMPIRank() == 0) {
        messager_ << "[Evolution::multiplicity]: dN/dy 1 = " << dNdeta
                  << ", dE/dy 1 = " << dEdeta;
        messager_.flush("info");
        messager_ << "[Evolution::multiplicity]: dN/dy 2 = " << dNdeta2
                  << ", dE/dy 2 = " << dEdeta2;
        messager_.flush("info");
        messager_ << "[Evolution::multiplicity]: gluon <p_T> = "
                  << dEdeta / dNdeta;
        messager_.flush("info");
    } else if (param->getUsePseudoRapidity() == 1) {
        m = param->getJacobianm();                                // in GeV
        P = 0.13 + 0.32 * pow(param->getRoots() / 1000., 0.115);  // in GeV
        dNdeta *=
            cosh(param->getRapidity())
            / (sqrt(pow(cosh(param->getRapidity()), 2.) + m * m / (P * P)));
        dEdeta *=
            cosh(param->getRapidity())
            / (sqrt(pow(cosh(param->getRapidity()), 2.) + m * m / (P * P)));

        if (param->getMPIRank() == 0) {
            messager_ << "[Evolution::multiplicity]: dN/deta 1 = " << dNdeta
                      << ", dE/deta 1 = " << dEdeta;
            messager_.flush("info");
            messager_ << "[Evolution::multiplicity]: dN/deta 2 = " << dNdeta2
                      << ", dE/deta 2 = " << dEdeta2;
            messager_.flush("info");
            messager_ << "[Evolution::multiplicity]: dN/deta_cut 1 = "
                      << dNdetaCut;
            messager_.flush("info");
            messager_ << "[Evolution::multiplicity]: dN/deta_cut 2 = "
                      << dNdetaCut2;
            messager_.flush("info");
            messager_ << "[Evolution::multiplicity]: gluon <p_T> = "
                      << dEdeta / dNdeta;
            messager_.flush("info");
        }
    }

    addPhaseAndRestart(
        "observables.gluon_multiplicity.report", multiplicityPhaseStart);

    if (dNdeta == 0.) {
        messager_ << "[Evolution::multiplicity]: No collision happened on "
                     "rank "
                  << param->getMPIRank()
                  << ". Restarting with new random number...";
        messager_.flush("warning");
        addPhaseAndRestart(
            "observables.gluon_multiplicity.cleanup", multiplicityPhaseStart);
        return 0;
    }

    if (it == itmax) {
        ofstream foutNN(NpartdNdy_name.c_str(), std::ios::out);
        foutNN << param->getNpart() << " " << dNdeta << " " << param->getTpp()
               << " " << param->getb() << " " << dEdeta << " "
               << param->getRandomSeed() << " "
               << "N/A"
               << " "
               << "N/A"
               << " "
               << "N/A"
               << " " << dNdetaCut << " " << dEdetaCut << " " << dNdetaCut2
               << " " << dEdetaCut2 << " "
               << g * g
                      / (4. * M_PI * 4. * M_PI
                         / (9.
                            * log(
                                pow(pow(muZero / 0.2, 2. / c)
                                        + pow(
                                            param->getRunWithThisFactorTimesQs()
                                                * param->getAverageQs() / 0.2,
                                            2. / c),
                                    c))))
               << endl;
        foutNN.close();
        writeGluonMultiplicityTarget(
            param, it, a, dtau, dNdeta, dNdeta2, dEdeta, dEdeta2, dNdetaCut,
            dEdetaCut, dNdetaCut2, dEdetaCut2, n, E, counter, bins, dkt);
    }
    addPhaseAndRestart("output.gluon_multiplicity", multiplicityPhaseStart);

    addPhaseAndRestart(
        "observables.gluon_multiplicity.cleanup", multiplicityPhaseStart);

    messager_ << "[Evolution::multiplicity]:  done.";
    messager_.flush("info");
    param->setSuccess(1);
    return 1;
}

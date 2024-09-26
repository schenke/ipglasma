// GaugeFix.cpp is part of the CYM evolution.
// Copyright (C) 2012 Bjoern Schenke.

#include "GaugeFix.h"

#include <complex>
#include <iostream>

#include "Matrix.h"

using std::cout;
using std::endl;

//**************************************************************************
// GaugeFix class.

void GaugeFix::FFTChi(
    FFT *fft, Lattice *lat, Group *group, Parameters *param, int steps) {
    const int N = param->getSize();
    int nn[2];
    nn[0] = N;
    nn[1] = N;
    int pos, posmX, posmY;
    int Nc = param->getNc();
    int Nc2m1 = Nc * Nc - 1;

    Matrix one(Nc, 1.);
    Matrix g(Nc), gdag(Nc);
    Matrix divA(Nc);
    Matrix UDx(Nc), UDy(Nc), UDxMx(Nc), UDyMy(Nc), Ux(Nc), Uy(Nc), UxMx(Nc),
        UyMy(Nc);

    const int max_gfiter = steps;

    Matrix zero(Nc, 0.);
    double gresidual_prev = 10000.;
    double gresidual = 0.;

    Matrix **chi;
    chi = new Matrix *[N * N];

    for (int i = 0; i < N * N; i++) {
        chi[i] = new Matrix(Nc, 0.);
    }

    cout << "gauge fixing" << endl;

    for (int gfiter = 0; gfiter < max_gfiter; gfiter++) {
        gresidual = 0.;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                pos = i * N + j;

                // use periodic boundary conditions to have fast convergence
                if (i == 0) {
                    posmX = (N - 1) * N + j;
                } else {
                    posmX = (i - 1) * N + j;
                }
                if (j == 0) {
                    posmY = i * N + N - 1;
                } else {
                    posmY = i * N + j - 1;
                }

                Ux = UDx = lat->cells[pos]->getUx();
                Uy = UDy = lat->cells[pos]->getUy();
                UxMx = UDxMx = lat->cells[posmX]->getUx();
                UyMy = UDyMy = lat->cells[posmY]->getUy();
                UDx.conjg();
                UDy.conjg();
                UDxMx.conjg();
                UDyMy.conjg();

                divA = (Ux - UxMx + Uy - UyMy);

                g = zero;
                for (int ig = 0; ig < Nc2m1; ig++) {
                    g = g
                        + ((divA)*group->getT(ig)).trace().imag()
                              * group->getT(ig);
                }

                *chi[pos] = gdag = g;
                gdag.conjg();

                gresidual +=
                    ((gdag * g).trace()).real() / static_cast<double>(Nc);

            }  // i loop
        }  // j loop

        gresidual /= N * N;

        if (gfiter % 10 == 0) {
            cout << gfiter << " " << gresidual << endl;
            gresidual_prev = gresidual;
        }

        if (gresidual < 1e-9) {
            break;
        }

        if (gresidual > gresidual_prev && gresidual < 1e-6) {
            // make sure it is progressively converging
            // otherwise, break the loop with less accuracy
            break;
        }

        fft->fftn(chi, chi, nn, 1);

#pragma omp parallel for
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double kx, ky, kt2;
                int localpos = i * N + j;
                kx = sin(
                    M_PI
                    * (-0.5 + static_cast<double>(i) / static_cast<double>(N)));
                ky = sin(
                    M_PI
                    * (-0.5 + static_cast<double>(j) / static_cast<double>(N)));
                kt2 = 4. * (kx * kx + ky * ky);  // lattice momentum squared
                *chi[localpos] = -1.5 * (1. / (kt2 + 1e-9)) * (*chi[localpos]);
            }
        }

        fft->fftn(chi, chi, nn, -1);

#pragma omp parallel
        {
            Matrix localg(Nc);
#pragma omp for
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    int localpos = i * N + j;
                    // exponentiate
                    localg = complex<double>(0, 1.) * (*chi[localpos]);
                    localg.expm();
                    // reunitarize
                    localg.reu();

                    if (localg(2) != localg(2)) {
                        cout << "problem at " << i << " " << j
                             << " with g=" << localg << endl;
                        localg = one;
                    }

                    lat->cells[localpos]->setg(localg);
                }
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                gaugeTransform(lat, param, i, j);
            }
        }
    }  // gfiter loop

    for (int i = 0; i < N * N; i++) {
        delete chi[i];
    }
    delete[] chi;
}

void GaugeFix::gaugeTransform(Lattice *lat, Parameters *param, int i, int j) {
    int N = param->getSize();
    int pos, posmX, posmY;
    int Nc = param->getNc();
    Matrix g(Nc), gdag(Nc);

    pos = i * N + j;

    // use periodic boundary conditions to have fast convergence
    if (i == 0) {
        posmX = (N - 1) * N + j;
    } else {
        posmX = (i - 1) * N + j;
    }
    if (j == 0) {
        posmY = i * N + N - 1;
    } else {
        posmY = i * N + j - 1;
    }

    g = gdag = lat->cells[pos]->getg();
    gdag.conjg();

    // gauge transform Ux and Uy
    lat->cells[pos]->setUx(g * lat->cells[pos]->getUx());
    lat->cells[pos]->setUy(g * lat->cells[pos]->getUy());

    lat->cells[posmX]->setUx(lat->cells[posmX]->getUx() * gdag);
    lat->cells[posmY]->setUy(lat->cells[posmY]->getUy() * gdag);

    // gauge transform Ex and Ey
    lat->cells[pos]->setE1(g * lat->cells[pos]->getE1() * gdag);
    lat->cells[pos]->setE2(g * lat->cells[pos]->getE2() * gdag);

    // gauge transform phi and pi
    lat->cells[pos]->setphi(g * lat->cells[pos]->getphi() * gdag);
    lat->cells[pos]->setpi(g * lat->cells[pos]->getpi() * gdag);
}

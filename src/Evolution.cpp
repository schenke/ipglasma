// Evolution.cpp is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Evolution.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "Fragmentation.h"
#include "Phys_consts.h"

using Fragmentation::kkp;
using PhysConst::hbarc;
using PhysConst::m_kaon;
using PhysConst::m_pion;
using PhysConst::m_proton;

using std::stringstream;

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;

//**************************************************************************
// Evolution class.

void Evolution::evolveU(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    // tau is the current time. The time argument of E^i is tau+dtau/2
    // we evolve to tau+dtau
    const int Nc = param->getNc();
    const int N = param->getSize();
    const double g = param->getg();

    const int n = 2;
    const Matrix one(Nc, 1.);

#pragma omp parallel
    {
        Matrix E1(Nc);
        Matrix E2(Nc);

        Matrix temp1(Nc);
        Matrix temp2(Nc);

#pragma omp for
        for (int pos = 0; pos < N * N; pos++) {
            // retrieve current E1 and E2 (that's the one defined at half a time
            // step in the future (from tau))
            E1 = complex<double>(0., g * g * dtau / (tau + dtau / 2.))
                 * lat->cells[pos]->getE1();
            // E1.expm(); // E1 now contains the exponential of i g^2
            // dtau/(tau+dtau/2)*E1

            temp2 = one + 1. / (double)n * E1;
            for (int in = 0; in < n - 1; in++) {
                temp1 = E1 * temp2;
                temp2 = one + 1. / (double)(n - 1 - in) * temp1;
            }

            E1 = temp2;

            E2 = complex<double>(0., g * g * dtau / (tau + dtau / 2.))
                 * lat->cells[pos]->getE2();
            // E2.expm(); // E2 now contains the exponential of i g^2
            // dtau/(tau+dtau/2)*E2

            temp2 = one + 1. / (double)n * E2;
            for (int in = 0; in < n - 1; in++) {
                temp1 = E2 * temp2;
                temp2 = one + 1. / (double)(n - 1 - in) * temp1;
            }

            E2 = temp2;

            lat->cells[pos]->setUx(E1 * lat->cells[pos]->getUx());
            lat->cells[pos]->setUy(E2 * lat->cells[pos]->getUy());
        }
    }
}

void Evolution::evolvePhi(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    // tau is the current time. The time argument of pi is tau+dtau/2
    // we evolve to tau+dtau
    const int Nc = param->getNc();
    const int N = param->getSize();

#pragma omp parallel
    {
        Matrix phi(Nc);
        Matrix pi(Nc);

#pragma omp for
        for (int pos = 0; pos < N * N; pos++) {
            // retrieve current phi (at time tau)
            phi = lat->cells[pos]->getphi();
            // retrieve current pi (at time tau+dtau/2)
            pi = lat->cells[pos]->getpi();

            phi = phi + (tau + dtau / 2.) * dtau * pi;

            // set the new phi (at time tau+dtau)
            lat->cells[pos]->setphi(phi);
        }
    }
}

void Evolution::evolvePi(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    const int Nc = param->getNc();
    const int N = param->getSize();

#pragma omp parallel
    {
        Matrix Ux(Nc);
        Matrix Uy(Nc);
        Matrix UxXm1(Nc);
        Matrix UyYm1(Nc);

        Matrix phi(Nc);
        Matrix phiX(Nc);   // this is \tilde{phi}_x
        Matrix phiY(Nc);   // this is \tilde{phi}_y
        Matrix phimX(Nc);  // this is \tilde{-phi}_x
        Matrix phimY(Nc);  // this is \tilde{-phi}_y

        // this will hold [phiX+phimX-2*phi+phiY+phimY-2*phi]
        Matrix bracket(Nc);
        Matrix pi(Nc);

#pragma omp for
        for (int pos = 0; pos < N * N; pos++) {
            // retrieve current Ux and Uy and compute conjugates
            Ux = lat->cells[pos]->getUx();
            Uy = lat->cells[pos]->getUy();
            // retrieve current pi (at time tau-dtau/2)
            pi = lat->cells[pos]->getpi();
            // retrieve current phi (at time tau) at this x_T
            phi = lat->cells[pos]->getphi();

            // retrieve current phi (at time tau) at x_T+1
            // parallel transport:
            phiX =
                Ux * Ux.prodABconj(lat->cells[lat->pospX[pos]]->getphi(), Ux);
            phiY =
                Uy * Uy.prodABconj(lat->cells[lat->pospY[pos]]->getphi(), Uy);

            // phi_{-x} should be defined as UxD*phimX*Ux with the Ux and UxD
            // reversed from the phi_{+x} case retrieve current phi (at time
            // tau) at x_T-1 parallel transport:
            UxXm1 = lat->cells[lat->posmX[pos]]->getUx();
            UyYm1 = lat->cells[lat->posmY[pos]]->getUy();

            phimX = Ux.prodAconjB(UxXm1, lat->cells[lat->posmX[pos]]->getphi())
                    * UxXm1;
            phimY = Ux.prodAconjB(UyYm1, lat->cells[lat->posmY[pos]]->getphi())
                    * UyYm1;

            // sum over both directions is included here
            bracket = phiX + phimX + phiY + phimY - 4. * phi;

            pi += dtau / (tau)*bracket;  // divide by \tau because this is
                                         // computing pi(tau+dtau/2) from
                                         // pi(tau-dtau/2) and phi(tau)

            // set the new pi (at time tau+dtau/2)
            lat->cells[pos]->setpi(pi);
        }
    }
}

void Evolution::evolveE(
    Lattice *lat, Parameters *param, double dtau, double tau) {
    const int Nc = param->getNc();
    const int N = param->getSize();
    const double g = param->getg();

#pragma omp parallel
    {
        Matrix Ux(Nc);
        Matrix Uy(Nc);
        Matrix temp1(Nc);  // can contain p or m
        Matrix temp2(Nc);
        Matrix temp3(Nc);
        Matrix En(Nc);
        Matrix phi(Nc);
        Matrix phiN(Nc);  // this is \tilde{phi}_x OR \tilde{phi}_y

        // plaquettes:
        Matrix U12(Nc);
        Matrix U1m2(Nc);
        Matrix U12Dag(Nc);  // equals (U21)
        Matrix U2m1(Nc);
        complex<double> trace;
        Matrix one(Nc, 1.);

#pragma omp for
        for (int pos = 0; pos < N * N; pos++) {
            int i = pos / N;
            int j = pos % N;
            int posmXpY = std::max(0, i - 1) * N + std::min(N - 1, j + 1);
            int pospXmY = std::min(N - 1, i + 1) * N + std::max(0, j - 1);

            // retrieve current E1 and E2 (that's the one defined at tau-dtau/2)
            En = lat->cells[pos]->getE1();
            // retrieve current phi (at time tau) at this x_T
            phi = lat->cells[pos]->getphi();
            // retrieve current phi (at time tau) at x_T+1
            phiN = lat->cells[lat->pospX[pos]]->getphi();
            // parallel transport:
            // retrieve current Ux and Uy
            Ux = lat->cells[pos]->getUx();
            phiN = Ux * Ux.prodABconj(phiN, Ux);

            // compute plaquettes:
            Uy = lat->cells[pos]->getUy();
            temp1 = lat->cells[lat->pospY[pos]]->getUx();  // UxYp1Dag
            temp1.conjg();
            U12 = (Ux * lat->cells[lat->pospX[pos]]->getUy())
                  * (Ux.prodABconj(temp1, Uy));

            temp1 = lat->cells[lat->posmY[pos]]->getUx();  // UxYm1Dag
            temp2 = lat->cells[pospXmY]->getUy();          // UyXp1Ym1Dag
            U1m2 =
                (Ux.prodABconj(Ux, temp2))
                * (Ux.prodAconjB(temp1, lat->cells[lat->posmY[pos]]->getUy()));

            temp1 = lat->cells[lat->posmX[pos]]->getUy();  // UyXm1Dag
            temp2 = lat->cells[posmXpY]->getUx();          // UxXm1Yp1Dag
            U2m1 =
                (Ux.prodABconj(Uy, temp2))
                * (Ux.prodAconjB(temp1, lat->cells[lat->posmX[pos]]->getUx()));

            U12Dag = U12;
            U12Dag.conjg();

            // do E1 update:

            temp3 = U1m2;
            temp3.conjg();

            temp1 = U12;
            temp1 += U1m2;
            temp1 -= U12Dag;
            temp1 -= temp3;

            trace = temp1.trace();
            temp1 -= trace / static_cast<double>(Nc) * one;

            temp2 = phiN * phi - phi * phiN;

            En += complex<double>(0., 1.) * tau * dtau / (2. * g * g) * temp1
                  + complex<double>(0., 1.) * dtau / tau * temp2;

            trace = En.trace();
            En -= (trace / static_cast<double>(Nc)) * one;
            lat->cells[pos]->setE1(En);

            // do E2 update:

            temp3 = U2m1;
            temp3.conjg();

            temp1 = U12Dag;
            temp1 += U2m1;
            temp1 -= U12;
            temp1 -= temp3;
            trace = temp1.trace();
            temp1 -= (trace / static_cast<double>(Nc)) * one;

            phiN = lat->cells[lat->pospY[pos]]->getphi();
            phiN = Uy * Uy.prodABconj(phiN, Uy);

            temp2 = phiN * phi - phi * phiN;

            En = lat->cells[pos]->getE2();
            En += complex<double>(0., 1.) * tau * dtau / (2. * g * g) * temp1
                  + complex<double>(0., 1.) * dtau / tau * temp2;

            trace = En.trace();
            En -= (trace / static_cast<double>(Nc)) * one;
            lat->cells[pos]->setE2(En);
        }
    }
}

void Evolution::checkGaussLaw(Lattice *lat, Parameters *param) {
    const int Nc = param->getNc();
    const int N = param->getSize();

    Matrix Ux(Nc);
    Matrix UxXm1(Nc);
    Matrix UxYm1(Nc);

    Matrix Uy(Nc);
    Matrix UyXm1(Nc);
    Matrix UyYm1(Nc);

    Matrix UxDag(Nc);
    Matrix UxXm1Dag(Nc);
    Matrix UxYm1Dag(Nc);

    Matrix UyDag(Nc);
    Matrix UyXm1Dag(Nc);
    Matrix UyYm1Dag(Nc);

    Matrix E1(Nc);
    Matrix E2(Nc);
    Matrix E1mX(Nc);
    Matrix E2mY(Nc);
    Matrix phi(Nc);
    Matrix pi(Nc);

    Matrix Gauss(Nc);
    double largest = 0;

    for (int pos = 0; pos < N * N; pos++) {
        // retrieve current Ux and Uy
        Ux = lat->cells[pos]->getUx();
        Uy = lat->cells[pos]->getUy();
        UxDag = Ux;
        UxDag.conjg();
        UyDag = Uy;
        UyDag.conjg();

        UxXm1 = lat->cells[lat->posmX[pos]]->getUx();
        UxYm1 = lat->cells[lat->posmY[pos]]->getUx();
        UxXm1Dag = UxXm1;
        UxXm1Dag.conjg();
        UxYm1Dag = UxYm1;
        UxYm1Dag.conjg();

        UyXm1 = lat->cells[lat->posmX[pos]]->getUy();
        UyYm1 = lat->cells[lat->posmY[pos]]->getUy();
        UyXm1Dag = UyXm1;
        UyXm1Dag.conjg();
        UyYm1Dag = UyYm1;
        UyYm1Dag.conjg();

        // retrieve current E1 and E2 (that's the one defined at tau-dtau/2)
        E1 = lat->cells[pos]->getE1();
        E2 = lat->cells[pos]->getE2();
        E1mX = lat->cells[lat->posmX[pos]]->getE1();
        E2mY = lat->cells[lat->posmY[pos]]->getE2();
        // retrieve current phi (at time tau) at this x_T
        phi = lat->cells[pos]->getphi();
        // retrieve current pi
        pi = lat->cells[pos]->getpi();

        Gauss = UxXm1Dag * E1mX * UxXm1 - E1 + UyYm1Dag * E2mY * UyYm1 - E2
                - complex<double>(0., 1.) * (phi * pi - pi * phi);

        if (Gauss.square() > largest) largest = Gauss.square();
    }
    cout << "Gauss violation=" << largest << endl;
}

void Evolution::run(Lattice *lat, Group *group, Parameters *param) {
    int Nc = param->getNc();
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
        cout << "maximal evolution time = " << maxtime << " fm" << endl;
    }

    // E and Pi at tau=dtau/2 are equal to the initial ones (at tau=0)
    // now evolve phi and U to time tau=dtau.
    evolvePhi(lat, param, dtau, 0.);
    evolveU(lat, param, dtau, 0.);

    int itmax = static_cast<int>(maxtime / (a * dtau) + 0.1);
    // int it0   = static_cast<int>(0.1/(a*dtau) + 0.1);
    // int it1   = static_cast<int>(0.2/(a*dtau) + 0.1);
    // int it2   = static_cast<int>(0.4/(a*dtau) + 0.1);
    // int it3   = static_cast<int>(0.6/(a*dtau) + 0.1);

    cout << "Starting evolution" << endl;
    cout << "itmax=" << itmax << endl;

    // do evolution
    for (int it = 1; it <= itmax; it++) {
        // if (it == 1 || it == it0 || it == it1 || it == it2
        //     || it == it3 || it == itmax) {
        if (it == itmax) {
            Tmunu(lat, param, it);
            // computes flow velocity and correct energy density
            u(lat, param, it, true);
        }

        if (it % 10 == 0) {
            cout << "Evolving to time " << it * a * dtau << " fm/c" << endl;
        }

        // evolve from time tau-dtau/2 to tau+dtau/2
        if (it < itmax) {
            evolvePi(lat, param, dtau, (it)*dtau);
            // the last argument is the current time tau.

            evolveE(lat, param, dtau, (it)*dtau);

            // evolve from time tau to tau+dtau
            evolvePhi(lat, param, dtau, (it)*dtau);
            evolveU(lat, param, dtau, (it)*dtau);
        } else if (it == itmax) {
            evolvePi(lat, param, dtau / 2., (it)*dtau);
            // the last argument is the current time tau.

            evolveE(lat, param, dtau / 2., (it)*dtau);
        }

        if (it == 1 && param->getWriteOutputs() == 3) {
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
                        if (pos > 0 && pos < (N - 1) * N + N - 1) {
                            g2mu2A = lat->cells[pos]->getg2mu2A();
                        } else
                            g2mu2A = 0;

                        if (pos > 0 && pos < (N - 1) * N + N - 1) {
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
                        if (pos > 0 && pos < (N - 1) * N + N - 1) {
                            g2mu2A = lat->cells[pos]->getg2mu2A();
                        } else
                            g2mu2A = 0;

                        if (pos > 0 && pos < (N - 1) * N + N - 1) {
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
        // if (it == 1 || it == it0 || it == it1 || it == it2
        //     || it == it3 || it == itmax) {
        if (it == 1 || it == itmax) {
            eccentricity(lat, param, it, 0.0, 0);
            // eccentricity(lat, param, it, 0.1, 0);
            // eccentricity(lat, param, it, 1., 0);
            // eccentricity(lat, param, it, 10., 0);

            success = multiplicity(lat, group, param, it);
        }

        if (success == 0) break;
    }
}

void Evolution::Tmunu(Lattice *lat, Parameters *param, int it) {
    double averageTtautau = 0.;
    double averageTtaueta = 0.;
    double averageTxx = 0.;

    int N = param->getSize();
    int Nc = param->getNc();
    int pos, posX, posY, posmX, posmY, posXY, posmXpY, pospXmY, pos2X, pos2Y,
        posX2Y, pos2XY;
    double L = param->getL();
    double a = L / N;  // lattice spacing in fm
    double g = param->getg();
    double dtau = param->getdtau();
    Matrix one(Nc, 1.);

    Matrix Ux(Nc);
    Matrix Uy(Nc);
    Matrix UxmX(Nc);
    Matrix UymY(Nc);
    Matrix UDx(Nc);
    Matrix UDy(Nc);
    Matrix UDxmX(Nc);
    Matrix UDymY(Nc);
    Matrix UDxmXpY(Nc);
    Matrix UDxpXpY(Nc);
    Matrix UxpX(Nc);
    Matrix UxpY(Nc);
    Matrix UDxpY(Nc);
    Matrix UxpXpY(Nc);
    Matrix UDypXmY(Nc);
    Matrix UypY(Nc);
    Matrix UypX(Nc);
    Matrix UDypX(Nc);
    Matrix UypXpY(Nc);
    Matrix UDypXpY(Nc);
    Matrix UymX(Nc);
    Matrix UxmXpY(Nc);
    Matrix UxmY(Nc);
    Matrix UDxmY(Nc);
    Matrix UypXmY(Nc);
    Matrix UDyp2X(Nc);
    Matrix Uyp2X(Nc);
    Matrix UDxpX(Nc);
    Matrix Uxp2Y(Nc);
    Matrix UDxp2Y(Nc);
    Matrix UDypY(Nc);
    Matrix UDymX(Nc);
    Matrix Uplaq(Nc), UplaqD(Nc), Uplaq1(Nc), Uplaq1D(Nc), Uplaq2(Nc);
    Matrix E1(Nc);
    Matrix E2(Nc);
    Matrix E1p(Nc);
    Matrix E2p(Nc);
    Matrix pi(Nc);
    Matrix piX(Nc);
    Matrix piY(Nc);
    Matrix piXY(Nc);
    Matrix phi(Nc);
    Matrix phiX(Nc);
    Matrix phiY(Nc);
    Matrix phiXY(Nc);
    Matrix phimX(Nc);
    Matrix phimY(Nc);
    Matrix phi2XY(Nc);
    Matrix phiX2Y(Nc);
    Matrix phi2X(Nc);
    Matrix phi2Y(Nc);
    Matrix phimXpY(Nc);
    Matrix phipXmY(Nc);
    Matrix phiTildeX(Nc);
    Matrix phiTildeY(Nc);
    Matrix phiTildeXY1(Nc);
    Matrix phiTildeXY2(Nc);

    // set plaquette in every cell
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;

            posX = lat->pospX[pos];
            posY = lat->pospY[pos];

            UDx = lat->cells[posY]->getUx();
            UDy = lat->cells[pos]->getUy();
            UDx.conjg();
            UDy.conjg();

            Uplaq = lat->cells[pos]->getUx()
                    * (lat->cells[posX]->getUy() * (UDx * UDy));
            lat->cells[pos]->setUplaq(Uplaq);
            if (i == N - 1) lat->cells[pos]->setUplaq(one);
        }
    }

    // T^\tau\tau, Txx, Tyy, Tetaeta:
    // electric part:
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;
            posX = lat->pospX[pos];
            posY = lat->pospY[pos];

            posXY = std::min(N - 1, i + 1) * N + std::min(N - 1, j + 1);

            E1 = lat->cells[pos]->getE1();
            E2 = lat->cells[pos]->getE2();
            E1p = lat->cells[posY]->getE1();
            E2p = lat->cells[posX]->getE2();  // shift y value in x direction

            pi = lat->cells[pos]->getpi();
            piX = lat->cells[posX]->getpi();
            piY = lat->cells[posY]->getpi();
            piXY = lat->cells[posXY]->getpi();

            lat->cells[pos]->setTtautau(
                (g * g / (it * dtau) / (it * dtau)
                     * real(
                         (E1 * E1).trace() + (E1p * E1p).trace()
                         + (E2 * E2).trace() + (E2p * E2p).trace())
                     / 2.  // trans.
                 + (((pi * pi).trace()).real() + ((piX * piX).trace()).real()
                    + ((piY * piY).trace()).real()
                    + ((piXY * piXY).trace()).real())
                       / 4.));  // long.
            lat->cells[pos]->setTxx(
                (g * g / (it * dtau) / (it * dtau)
                     * real(
                         -1. * (E1 * E1).trace() - (E1p * E1p).trace()
                         + (E2 * E2).trace() + (E2p * E2p).trace())
                     / 2.  // trans.
                 + (((pi * pi).trace()).real() + ((piX * piX).trace()).real()
                    + ((piY * piY).trace()).real()
                    + ((piXY * piXY).trace()).real())
                       / 4.));  // long.
            lat->cells[pos]->setTyy(
                (g * g / (it * dtau) / (it * dtau)
                     * real(
                         (E1 * E1).trace() + (E1p * E1p).trace()
                         - (E2 * E2).trace() - (E2p * E2p).trace())
                     / 2.  // trans.
                 + (((pi * pi).trace()).real() + ((piX * piX).trace()).real()
                    + ((piY * piY).trace()).real()
                    + ((piXY * piXY).trace()).real())
                       / 4.));  // long.
            lat->cells[pos]->setTetaeta(
                1. / (it * dtau) / (it * dtau)
                * ((g * g / (it * dtau) / (it * dtau)
                        * real(
                            (E1 * E1).trace() + (E1p * E1p).trace()
                            + (E2 * E2).trace() + (E2p * E2p).trace())
                        / 2.  // trans.
                    - (((pi * pi).trace()).real() + ((piX * piX).trace()).real()
                       + ((piY * piY).trace()).real()
                       + ((piXY * piXY).trace()).real())
                          / 4.)));  // long.
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;

            posX = lat->pospX[pos];
            posY = lat->pospY[pos];

            posXY = std::min(N - 1, i + 1) * N + std::min(N - 1, j + 1);

            Uplaq = lat->cells[pos]->getUplaq();

            phi = lat->cells[pos]->getphi();
            phiX = lat->cells[posX]->getphi();
            phiY = lat->cells[posY]->getphi();
            phiXY = lat->cells[posXY]->getphi();

            Ux = lat->cells[pos]->getUx();
            Uy = lat->cells[pos]->getUy();
            UDx = Ux;
            UDx.conjg();
            UDy = Uy;
            UDy.conjg();

            phiTildeX = Ux * phiX * UDx;
            phiTildeY = Uy * phiY * UDy;

            // same at one up in the other direction
            Ux = lat->cells[posY]->getUx();
            Uy = lat->cells[posX]->getUy();
            UDx = Ux;
            UDx.conjg();
            UDy = Uy;
            UDy.conjg();

            phiTildeXY1 = Ux * phiXY * UDx;
            phiTildeXY2 = Uy * phiXY * UDy;

            lat->cells[pos]->setTtautau(
                lat->cells[pos]->getTtautau()
                + 2. / pow(g, 2.)
                      * (static_cast<double>(Nc) - (Uplaq.trace()).real())
                + 0.5 / (it * dtau) / (it * dtau)
                      * (real(((phi - phiTildeX) * (phi - phiTildeX)).trace())
                         + real(((phiY - phiTildeXY1) * (phiY - phiTildeXY1))
                                    .trace())
                         + real(((phi - phiTildeY) * (phi - phiTildeY)).trace())
                         + real(((phiX - phiTildeXY2) * (phiX - phiTildeXY2))
                                    .trace())));

            lat->cells[pos]->setTxx(
                lat->cells[pos]->getTxx()
                + 2. / pow(g, 2.)
                      * (static_cast<double>(Nc) - (Uplaq.trace()).real())
                + 0.5 / (it * dtau) / (it * dtau)
                      * (real(((phi - phiTildeX) * (phi - phiTildeX)).trace())
                         + real(((phiY - phiTildeXY1) * (phiY - phiTildeXY1))
                                    .trace())
                         - real(((phi - phiTildeY) * (phi - phiTildeY)).trace())
                         - real(((phiX - phiTildeXY2) * (phiX - phiTildeXY2))
                                    .trace())));

            lat->cells[pos]->setTyy(
                lat->cells[pos]->getTyy()
                + 2. / pow(g, 2.)
                      * (static_cast<double>(Nc) - (Uplaq.trace()).real())
                + 0.5 / (it * dtau) / (it * dtau)
                      * (-real(((phi - phiTildeX) * (phi - phiTildeX)).trace())
                         - real(((phiY - phiTildeXY1) * (phiY - phiTildeXY1))
                                    .trace())
                         + real(((phi - phiTildeY) * (phi - phiTildeY)).trace())
                         + real(((phiX - phiTildeXY2) * (phiX - phiTildeXY2))
                                    .trace())));

            lat->cells[pos]->setTetaeta(
                lat->cells[pos]->getTetaeta()
                + 1. / (it * dtau) / (it * dtau)
                      * (-2. / pow(g, 2.) * (Nc - (Uplaq.trace()).real())
                         + 0.5 / (it * dtau) / (it * dtau)
                               * (+real(((phi - phiTildeX) * (phi - phiTildeX))
                                            .trace())
                                  + real(((phiY - phiTildeXY1)
                                          * (phiY - phiTildeXY1))
                                             .trace())
                                  + real(((phi - phiTildeY) * (phi - phiTildeY))
                                             .trace())
                                  + real(((phiX - phiTildeXY2)
                                          * (phiX - phiTildeXY2))
                                             .trace()))));
        }
    }

    for (pos = 0; pos < N * N; pos++) {
        // clean up numerical noise outside the interaction region
        if (lat->cells[pos]->getg2mu2A() < 1e-12
            || lat->cells[pos]->getg2mu2B() < 1e-12) {
            lat->cells[pos]->setEpsilon(0.);
            lat->cells[pos]->setTtautau(0.);
            lat->cells[pos]->setTxx(0.);
            lat->cells[pos]->setTyy(0.);
            lat->cells[pos]->setTetaeta(0.);
        } else {
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

    // T^\tau x, T^\tau y
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;
            posX = lat->pospX[pos];
            posY = lat->pospY[pos];
            posXY = std::min(N - 1, i + 1) * N + std::min(N - 1, j + 1);

            posmX = lat->posmX[pos];
            posmY = lat->posmY[pos];

            posmXpY = std::max(0, i - 1) * N + std::min(N - 1, j + 1);
            pospXmY = std::min(N - 1, i + 1) * N + std::max(0, j - 1);

            pos2X = std::min(N - 1, i + 2) * N + j;
            pos2Y = i * N + std::min(N - 1, j + 2);

            pos2XY = std::min(N - 1, i + 2) * N + std::min(N - 1, j + 1);
            posX2Y = std::min(N - 1, i + 1) * N + std::min(N - 1, j + 2);

            E1 = lat->cells[pos]->getE1();
            E2 = lat->cells[pos]->getE2();
            E1p = lat->cells[posY]->getE1();  // shift x value in y direction
            E2p = lat->cells[posX]->getE2();  // shift y value in x direction

            Uplaq = lat->cells[pos]->getUplaq();
            Uplaq1 = lat->cells[posmX]->getUplaq();
            Uplaq2 = lat->cells[posmY]->getUplaq();
            UplaqD = Uplaq;
            UplaqD.conjg();
            Uplaq1D = Uplaq1;
            Uplaq1D.conjg();

            pi = lat->cells[pos]->getpi();
            piX = lat->cells[posX]->getpi();
            piY = lat->cells[posY]->getpi();
            piXY = lat->cells[posXY]->getpi();

            phi = lat->cells[pos]->getphi();
            phimX = lat->cells[posmX]->getphi();
            phiX = lat->cells[posX]->getphi();
            phimY = lat->cells[posmY]->getphi();
            phiY = lat->cells[posY]->getphi();
            phiXY = lat->cells[posXY]->getphi();
            phimXpY = lat->cells[posmXpY]->getphi();
            phipXmY = lat->cells[pospXmY]->getphi();
            phi2X = lat->cells[pos2X]->getphi();
            phi2XY = lat->cells[pos2XY]->getphi();
            phi2Y = lat->cells[pos2Y]->getphi();
            phiX2Y = lat->cells[posX2Y]->getphi();

            Ux = lat->cells[pos]->getUx();
            UDx = Ux;
            UDx.conjg();

            UxmX = lat->cells[posmX]->getUx();
            UDxmX = lat->cells[posmX]->getUx();
            UDxmX.conjg();
            UxmXpY = lat->cells[posmXpY]->getUx();
            UDxmXpY = lat->cells[posmXpY]->getUx();
            UDxmXpY.conjg();

            UxpX = lat->cells[posX]->getUx();
            UxpY = lat->cells[posY]->getUx();
            UDxpX = UxpX;
            UDxpX.conjg();
            UDxpY = UxpY;
            UDxpY.conjg();

            UxpXpY = lat->cells[posXY]->getUx();
            UDxpXpY = lat->cells[posXY]->getUx();
            UDxpXpY.conjg();
            UxmXpY = lat->cells[posmXpY]->getUx();
            UDxmXpY = UxmXpY;
            UDxmXpY.conjg();

            Uy = lat->cells[pos]->getUy();
            UDy = Uy;
            UDy.conjg();

            UymY = lat->cells[posmY]->getUy();
            UDymY = lat->cells[posmY]->getUy();
            UDymY.conjg();
            UypXmY = lat->cells[pospXmY]->getUy();
            UDypXmY = lat->cells[pospXmY]->getUy();
            UDypXmY.conjg();

            UypY = lat->cells[posY]->getUy();
            UypX = lat->cells[posX]->getUy();
            UDypX = UypX;
            UDypX.conjg();

            UDypY = lat->cells[posY]->getUy();
            UDypY.conjg();

            UDxpX = lat->cells[posX]->getUx();
            UDxpX.conjg();

            Uyp2X = lat->cells[pos2X]->getUy();
            UDyp2X = Uyp2X;
            UDyp2X.conjg();
            Uxp2Y = lat->cells[pos2Y]->getUx();
            UDxp2Y = Uxp2Y;
            UDxp2Y.conjg();

            UypXpY = lat->cells[posXY]->getUy();
            UDypXpY = lat->cells[posXY]->getUy();
            UDypXpY.conjg();
            UymX = lat->cells[posmX]->getUy();
            UDymX = UymX;
            UDymX.conjg();
            UxmY = lat->cells[posmY]->getUx();
            UDxmY = UxmY;
            UDxmY.conjg();
            UypXmY = lat->cells[pospXmY]->getUy();

            // note that the minus sign of the first terms in T^\taux and
            // T^\tauy comes from the direction of the plaquettes - I am using
            // +F^{yx} instead of -F^{xy} if you like.
            lat->cells[pos]->setTtaux(
                +2. / (it * dtau) / 8.
                    * (E2
                           * (Uy * UxpY * UDypX * UDx - Ux * UypX * UDxpY * UDy
                              - (Uy * UxpY * UDypX * UDx
                                 - Ux * UypX * UDxpY * UDy)
                                        .trace()
                                    / static_cast<double>(Nc) * one
                              + UDxmX * UymX * UxmXpY * UDy
                              - Uy * UDxmXpY * UDymX * UxmX
                              - (UDxmX * UymX * UxmXpY * UDy
                                 - Uy * UDxmXpY * UDymX * UxmX)
                                        .trace()
                                    / static_cast<double>(Nc) * one)
                       + E2p
                             * (UypX * UxpXpY * UDyp2X * UDxpX
                                - UxpX * Uyp2X * UDxpXpY * UDypX
                                - (UypX * UxpXpY * UDyp2X * UDxpX
                                   - UxpX * Uyp2X * UDxpXpY * UDypX)
                                          .trace()
                                      / static_cast<double>(Nc) * one
                                + UDx * Uy * UxpY * UDypX
                                - UypX * UDxpY * UDy * Ux
                                - (UDx * Uy * UxpY * UDypX
                                   - UypX * UDxpY * UDy * Ux)
                                          .trace()
                                      / static_cast<double>(Nc) * one))
                          .trace()
                          .imag()
                - 2. / 8. / (it * dtau)
                      * (pi * (Ux * phiX * UDx - UDxmX * phimX * UxmX)
                         + piY
                               * (UxpY * phiXY * UDxpY
                                  - UDxmXpY * phimXpY * UxmXpY)
                         + piX * (UxpX * phi2X * UDxpX - UDx * phi * Ux)
                         + piXY
                               * (UxpXpY * phi2XY * UDxpXpY
                                  - UDxpY * phiY * UxpY))
                            .trace()
                            .real());

            lat->cells[pos]->setTtauy(
                +2. / (it * dtau) / 8.
                    * (E1
                           * (Ux * UypX * UDxpY * UDy - Uy * UxpY * UDypX * UDx
                              - (Ux * UypX * UDxpY * UDy
                                 - Uy * UxpY * UDypX * UDx)
                                        .trace()
                                    / static_cast<double>(Nc) * one
                              + UDymY * UxmY * UypXmY * UDx
                              - Ux * UDypXmY * UDxmY * UymY
                              - (UDymY * UxmY * UypXmY * UDx
                                 - Ux * UDypXmY * UDxmY * UymY)
                                        .trace()
                                    / static_cast<double>(Nc) * one)
                       + E1p
                             * (UxpY * UypXpY * UDxp2Y * UDypY
                                - UypY * Uxp2Y * UDypXpY * UDxpY
                                - (UxpY * UypXpY * UDxp2Y * UDypY
                                   - UypY * Uxp2Y * UDypXpY * UDxpY)
                                          .trace()
                                      / static_cast<double>(Nc) * one
                                + UDy * Ux * UypX * UDxpY
                                - UxpY * UDypX * UDx * Uy
                                - (UDy * Ux * UypX * UDxpY
                                   - UxpY * UDypX * UDx * Uy)
                                          .trace()
                                      / static_cast<double>(Nc) * one))
                          .trace()
                          .imag()
                - 2. / 8. / (it * dtau)
                      * (pi * (Uy * phiY * UDy - UDymY * phimY * UymY)
                         + piX
                               * (UypX * phiXY * UDypX
                                  - UDypXmY * phipXmY * UypXmY)
                         + piY * (UypY * phi2Y * UDypY - UDy * phi * Uy)
                         + piXY
                               * (UypXpY * phiX2Y * UDypXpY
                                  - UDypX * phiX * UypX))
                            .trace()
                            .real());

            lat->cells[pos]->setTtaueta(
                g / (it * dtau) / (it * dtau) / (it * dtau)
                * (E1 * (Ux * phiX * UDx - phi)
                   + E1p * (UxpY * phiXY * UDxpY - phiY)
                   + E2 * (Uy * phiY * UDy - phi)
                   + E2p * (UypX * phiXY * UDypX - phiX))
                      .trace()
                      .real());

            // T^xy
            lat->cells[pos]->setTxy(
                2. / (it * dtau) / (it * dtau)
                * (-1. / 4. * g * g * (E1 + Uy * E1p * UDy)
                       * (E2 + Ux * E2p * UDx)
                   + 1. / 4.
                         * ((Ux * phiX * UDx - phi) * (Uy * phiY * UDy - phi)
                            + Uy * (UxpY * phiXY * UDxpY - phiY) * UDy
                                  * (Uy * phiY * UDy - phi)
                            + (Ux * phiX * UDx - phi) * Ux
                                  * (UypX * phiXY * UDypX - phiX) * UDx
                            + Uy * (UxpY * phiXY * UDxpY - phiY) * UDy * Ux
                                  * (UypX * phiXY * UDypX - phiX) * UDx))
                      .trace()
                      .real());

            lat->cells[pos]->setTxeta(
                -2. / (it * dtau) / (it * dtau)
                * (1. / 4. * g
                       * (E1 * (pi + Ux * piX * UDx)
                          + E1p * (piY + UxpY * piXY * UDxpY))
                             .trace()
                             .real()
                   + 1. / 8. / g
                         * ((Ux * UypX * UDxpY * UDy - Uy * UxpY * UDypX * UDx
                             - (Ux * UypX * UDxpY * UDy
                                - Uy * UxpY * UDypX * UDx)
                                       .trace()
                                   / static_cast<double>(Nc) * one
                             + Uy * UDxmXpY * UDymX * UxmX
                             - UDxmX * UymX * UxmXpY * UDy
                             - (Uy * UDxmXpY * UDymX * UxmX
                                - UDxmX * UymX * UxmXpY * UDy)
                                       .trace()
                                   / static_cast<double>(Nc) * one)
                                * (Uy * phiY * UDy - phi)
                            + (UypX * UDxpY * UDy * Ux - UDx * Uy * UxpY * UDypX
                               - (UypX * UDxpY * UDy * Ux
                                  - UDx * Uy * UxpY * UDypX)
                                         .trace()
                                     / static_cast<double>(Nc) * one
                               + UxpX * Uyp2X * UDxpXpY * UDypX
                               - UypX * UxpXpY * UDyp2X * UDxpX
                               - (UxpX * Uyp2X * UDxpXpY * UDypX
                                  - UypX * UxpXpY * UDyp2X * UDxpX)
                                         .trace()
                                     / static_cast<double>(Nc) * one)
                                  * (UypX * phiXY * UDypX - phiX))
                               .trace()
                               .imag()));

            lat->cells[pos]->setTyeta(
                -2. / (it * dtau) / (it * dtau)
                * (1. / 4. * g
                       * (E2 * (pi + Uy * piY * UDy)
                          + E2p * (piX + UypX * piXY * UDypX))
                             .trace()
                             .real()
                   + 1. / 8. / g
                         * ((Uy * UxpY * UDypX * UDx - Ux * UypX * UDxpY * UDy
                             - (Uy * UxpY * UDypX * UDx
                                - Ux * UypX * UDxpY * UDy)
                                       .trace()
                                   / static_cast<double>(Nc) * one
                             + Ux * UDypXmY * UDxmY * UymY
                             - UDymY * UxmY * UypXmY * UDx
                             - (Ux * UDypXmY * UDxmY * UymY
                                - UDymY * UxmY * UypXmY * UDx)
                                       .trace()
                                   / static_cast<double>(Nc) * one)
                                * (Ux * phiX * UDx - phi)
                            + (UxpY * UDypX * UDx * Uy - UDy * Ux * UypX * UDxpY
                               - (UxpY * UDypX * UDx * Uy
                                  - UDy * Ux * UypX * UDxpY)
                                         .trace()
                                     / static_cast<double>(Nc) * one
                               + UypY * Uxp2Y * UDypXpY * UDxpY
                               - UxpY * UypXpY * UDxp2Y * UDypY
                               - (UypY * Uxp2Y * UDypXpY * UDxpY
                                  - UxpY * UypXpY * UDxp2Y * UDypY)
                                         .trace()
                                     / static_cast<double>(Nc) * one)
                                  * (UxpY * phiXY * UDxpY - phiY))
                               .trace()
                               .imag()));

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

            averageTtautau +=
                lat->cells[pos]->getTtautau() * lat->cells[pos]->getTtautau();
            averageTtaueta +=
                lat->cells[pos]->getTtaueta() * lat->cells[pos]->getTtaueta();
            averageTxx += lat->cells[pos]->getTxx() * lat->cells[pos]->getTxx();
        }
    }
    averageTtautau /= double(N);
    averageTtaueta /= double(N);
    averageTxx /= double(N);
}

void Evolution::u(Lattice *lat, Parameters *param, int it, bool finalFlag) {
    MyEigen myeigen;
    myeigen.flowVelocity4D(lat, param, it, finalFlag);
}

void Evolution::anisotropy(Lattice *lat, Parameters *param, int it) {
    stringstream straniso_name;
    straniso_name << "anisotropy" << param->getEventId() << ".dat";
    string aniso_name;
    aniso_name = straniso_name.str();

    ofstream foutAniso(aniso_name.c_str(), std::ios::app);
    int N = param->getSize();
    double L = param->getL();
    double a = L / N;  // lattice spacing in fm

    double num = 0., den = 0.;
    int pos;
    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            pos = ix * N + iy;
            if (lat->cells[pos]->getTtautau() > 10.) {
                num += lat->cells[pos]->getTxx() - lat->cells[pos]->getTyy();
                den += lat->cells[pos]->getTxx() + lat->cells[pos]->getTyy();
            }
        }
    }

    foutAniso << it * a * param->getdtau() << " " << num / den << endl;
    foutAniso.close();
}

void Evolution::eccentricity(
    Lattice *lat, Parameters *param, int it, double cutoff, int doAniso) {
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

            if (param->getRunningCoupling()) {
                if (pos > 0 && pos < (N - 1) * N + N - 1) {
                    g2mu2A = lat->cells[pos]->getg2mu2A();
                } else
                    g2mu2A = 0;

                if (pos > 0 && pos < (N - 1) * N + N - 1) {
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

                    gfactor = g * g / (4. * M_PI * alphas);
                }
            } else
                gfactor = 1.;

            if (lat->cells[pos]->getEpsilon() * gfactor
                < cutoff)  // this is 1/fm^4, so Lambda_QCD^{-4} (because
                           // \Lambda_QCD is roughly 1/fm)
            {
                weight = 0.;
            } else {
                // weight = lat->cells[pos]->getEpsilon() * gfactor;
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

            if (param->getRunningCoupling()) {
                if (pos > 0 && pos < (N - 1) * N + N - 1) {
                    g2mu2A = lat->cells[pos]->getg2mu2A();
                } else
                    g2mu2A = 0;

                if (pos > 0 && pos < (N - 1) * N + N - 1) {
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

                    gfactor = g * g / (4. * M_PI * alphas);
                }
            } else
                gfactor = 1.;

            if (lat->cells[pos]->getEpsilon() * gfactor
                < cutoff)  // this is 1/fm^4, so Lambda_QCD^{-4}
            {
                weight = 0.;
            } else {
                // weight = lat->cells[pos]->getEpsilon() * gfactor;
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

            if (param->getRunningCoupling()) {
                if (pos > 0 && pos < (N - 1) * N + N - 1) {
                    g2mu2A = lat->cells[pos]->getg2mu2A();
                } else
                    g2mu2A = 0;

                if (pos > 0 && pos < (N - 1) * N + N - 1) {
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

                    gfactor = g * g / (4. * M_PI * alphas);
                }
            } else
                gfactor = 1.;

            if (lat->cells[pos]->getEpsilon() * gfactor
                < cutoff)  // this is 1/fm^4, so Lambda_QCD^{-4}
            {
                weight = 0.;
            } else {
                // weight = lat->cells[pos]->getEpsilon() * gfactor;
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
    param->setEccentricity2(eccentricity2);
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
    cout << "Reading n(k_T) from file ";
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
    cout << mult_name.c_str() << " ... ";

    // open file

    ifstream fin2;
    stringstream strmult_name2;
    strmult_name2 << "NpartdNdy" << param->getEventId() << ".dat";
    string mult_name2;
    mult_name2 = strmult_name2.str();
    fin2.open(mult_name2.c_str());
    cout << mult_name2.c_str() << " ... ";

    // read file

    if (fin) {
        for (int ikt = 0; ikt < 100; ikt++) {
            if (!fin.eof()) {
                fin >> dummy;
                fin >> kt;
                fin >> nkt;
                nIn[ikt] = atof(nkt.c_str());
                fin >> dummy >> Tpp >> b >> Npart;
                if (ikt == 0) dkt = atof(kt.c_str());
                if (ikt == 1) dkt = dkt - atof(kt.c_str());
            }
            cout << nIn[ikt] << endl;
        }
        fin.close();
        cout << " done." << endl;
    } else {
        cout << "[Evolution.cpp:readNkt]: File " << mult_name.c_str()
             << " does not exist. Exiting." << endl;
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
        cout << " done." << endl;
    } else {
        cout << "[Evolution.cpp:readNkt]: File " << mult_name2.c_str()
             << " does not exist. Exiting." << endl;
        exit(1);
    }

    double m, P;
    m = param->getJacobianm();                                // in GeV
    P = 0.13 + 0.32 * pow(param->getRoots() / 1000., 0.115);  // in GeV
    double dNdeta2;
    dNdeta2 = 0.;

    for (int ik = 0; ik < 100; ik++) {
        if (param->getUsePseudoRapidity() == 0) {
            dNdeta2 += nIn[ik] * (ik + 0.5) * dkt * dkt * 2.
                       * M_PI;  // integrate, gives a ik*dkt*2pi*dkt
        } else {
            dNdeta2 +=
                nIn[ik] * (ik + 0.5) * dkt * dkt * 2. * M_PI
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
    int N = param->getSize();
    int Nc = param->getNc();
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
    cout << "Measuring multiplicity ... " << endl;

    // fix transverse Coulomb gauge
    GaugeFix gaugefix;

    double maxtime;
    if (param->getInverseQsForMaxTime() == 1) {
        maxtime = 1. / param->getAverageQs() * hbarc;
        cout << "maximal evolution time = " << maxtime << " fm" << endl;
    } else {
        maxtime = param->getMaxtime();  // maxtime is in fm
    }

    int itmax = static_cast<int>(floor(maxtime / (a * dtau) + 1e-10));

    gaugefix.FFTChi(fft, lat, group, param, 4000);
    // gauge is fixed

    Matrix **E1;
    E1 = new Matrix *[N * N];

    for (int i = 0; i < N * N; i++) {
        E1[i] = new Matrix(Nc, 0.);
    }

    double g2mu2A, g2mu2B, gfactor, alphas = 0., Qs = 0.;
    double c = param->getc();
    double muZero = param->getMuZero();

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;

            if (param->getRunningCoupling()) {
                if (pos > 0 && pos < (N - 1) * N + N - 1) {
                    g2mu2A = lat->cells[pos]->getg2mu2A();
                } else
                    g2mu2A = 0;

                if (pos > 0 && pos < (N - 1) * N + N - 1) {
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

                    gfactor = g * g / (4. * M_PI * alphas);
                }
            } else
                gfactor = 1.;

            if (param->getRunWithkt() == 0) {
                *E1[pos] = lat->cells[pos]->getE1()
                           * sqrt(gfactor);  // replace one of the 1/g in the
                                             // lattice E^i by the running one
            } else {
                *E1[pos] = lat->cells[pos]->getE1();
            }
        }
    }

    // do Fourier transforms
    fft->fftn(E1, E1, nn, 1);

    for (int ik = 0; ik < bins; ik++) {
        n[ik] = 0.;
        E[ik] = 0.;
        n2[ik] = 0.;
        counter[ik] = 0;
    }

    const int hbins = 2000;
    // double Nh[hbins+1], Eh[hbins+1], Ehgsl[hbins+1], NhL[hbins+1],
    // NhLgsl[hbins+1], NhH[hbins+1], NhHgsl[hbins+1];
    double Nhgsl[hbins + 1];
    double Ng;
    // for (int ih=0; ih<=hbins; ih++)
    //   {
    //     Nh[ih]=0.;
    //     Eh[ih]=0.;
    //     NhL[ih]=0.;
    //     NhH[ih]=0.;
    //   }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            nkt = 0.;
            pos = i * N + j;
            npos = (N - i) * N + (N - j);

            kx = 2. * M_PI
                 * (-0.5 + static_cast<double>(i) / static_cast<double>(N));
            ky = 2. * M_PI
                 * (-0.5 + static_cast<double>(j) / static_cast<double>(N));
            kt2 = 4.
                  * (sin(kx / 2.) * sin(kx / 2.)
                     + sin(ky / 2.) * sin(ky / 2.));  //
            omega2 = 4.
                     * (sin(kx / 2.) * sin(kx / 2.)
                        + sin(ky / 2.)
                              * sin(ky / 2.));  // lattice dispersion relation
                                                // (this is omega squared)
            if (i != 0 && j != 0) {
                if (omega2 != 0) {
                    nkt = 2. / sqrt(omega2) / static_cast<double>(N * N)
                          * (g * g / ((it - 0.5) * dtau)
                             * ((((*E1[pos]) * (*E1[npos])).trace()).real()));
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
                        E[ik] += sqrt(omega2) * hbarc / a * nkt / dkt / 2 / M_PI
                                 / sqrt(kt2) * 2 * M_PI * sqrt(kt2) * dkt * N
                                 * N / M_PI / M_PI / 2. / 2.;
                        n2[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2);
                        // dividing by bin size; bin is dkt times Jacobian
                        // k(=ik*dkt) times 2Pi in phi times the correct number
                        // of counts for an infinite lattice: area in bin
                        // divided by total area
                        counter[ik] += 1;  // number of entries in n[ik]
                    }
                }
            }
        }
    }

    /// -------- 2 ---------

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;

            if (param->getRunningCoupling()) {
                if (pos > 0 && pos < (N - 1) * N + N - 1) {
                    g2mu2A = lat->cells[pos]->getg2mu2A();
                } else
                    g2mu2A = 0;

                if (pos > 0 && pos < (N - 1) * N + N - 1) {
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

                    gfactor = g * g / (4. * M_PI * alphas);
                }
            } else
                gfactor = 1.;

            if (param->getRunWithkt() == 0) {
                *E1[pos] = lat->cells[pos]->getE2() * sqrt(gfactor);  // "
            } else {
                *E1[pos] = lat->cells[pos]->getE2();
            }
        }
    }

    fft->fftn(E1, E1, nn, 1);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            nkt = 0.;
            pos = i * N + j;
            npos = (N - i) * N + (N - j);

            kx = 2. * M_PI
                 * (-0.5 + static_cast<double>(i) / static_cast<double>(N));
            ky = 2. * M_PI
                 * (-0.5 + static_cast<double>(j) / static_cast<double>(N));
            kt2 = 4.
                  * (sin(kx / 2.) * sin(kx / 2.)
                     + sin(ky / 2.) * sin(ky / 2.));  //
            omega2 = 4.
                     * (sin(kx / 2.) * sin(kx / 2.)
                        + sin(ky / 2.)
                              * sin(ky / 2.));  // lattice dispersion relation
                                                // (this is omega squared)

            // i=0 or j=0 have no negative k_T value available

            if (i != 0 && j != 0) {
                if (omega2 != 0) {
                    nkt = 2. / sqrt(omega2) / static_cast<double>(N * N)
                          * (g * g / ((it - 0.5) * dtau)
                             * (((((*E1[pos]) * (*E1[npos])).trace()).real())));
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
                        E[ik] += sqrt(omega2) * hbarc / a * nkt / dkt / 2 / M_PI
                                 / sqrt(kt2) * 2 * M_PI * sqrt(kt2) * dkt * N
                                 * N / M_PI / M_PI / 2. / 2.;
                        n2[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2);
                        // dividing by bin size; bin is dkt times Jacobian
                        // k(=ik*dkt) times 2Pi in phi times the correct number
                        // of counts for an infinite lattice: area in bin
                        // divided by total area
                    }
                }
            }
        }
    }

    /// ------3 --------

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;

            if (param->getRunningCoupling()) {
                if (pos > 0 && pos < (N - 1) * N + N - 1) {
                    g2mu2A = lat->cells[pos]->getg2mu2A();
                } else
                    g2mu2A = 0;

                if (pos > 0 && pos < (N - 1) * N + N - 1) {
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

                    gfactor = g * g / (4. * M_PI * alphas);
                }
            } else
                gfactor = 1.;

            if (param->getRunWithkt() == 0) {
                *E1[pos] = lat->cells[pos]->getpi()
                           * sqrt(gfactor);  // replace the only 1/g by the
                                             // running one (physical pi goes
                                             // like 1/g, like physical E^i)
            } else {
                *E1[pos] = lat->cells[pos]->getpi();
            }
        }
    }

    // do Fourier transforms
    fft->fftn(E1, E1, nn, 1);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            nkt = 0.;
            pos = i * N + j;
            npos = (N - i) * N + (N - j);

            kx = 2. * M_PI
                 * (-0.5 + static_cast<double>(i) / static_cast<double>(N));
            ky = 2. * M_PI
                 * (-0.5 + static_cast<double>(j) / static_cast<double>(N));
            kt2 = 4.
                  * (sin(kx / 2.) * sin(kx / 2.)
                     + sin(ky / 2.) * sin(ky / 2.));  //
            omega2 = 4.
                     * (sin(kx / 2.) * sin(kx / 2.)
                        + sin(ky / 2.)
                              * sin(ky / 2.));  // lattice dispersion relation
                                                // (this is omega squared)

            // i=0 or j=0 have no negative k_T value available

            if (i != 0 && j != 0) {
                if (omega2 != 0) {
                    nkt = 2. / sqrt(omega2) / static_cast<double>(N * N)
                          * (((it - 0.5) * dtau)
                             * ((((*E1[pos]) * (*E1[npos])).trace()).real()));
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
                        E[ik] += sqrt(omega2) * hbarc / a * nkt / dkt / 2 / M_PI
                                 / sqrt(kt2) * 2 * M_PI * sqrt(kt2) * dkt * N
                                 * N / M_PI / M_PI / 2. / 2.;
                        n2[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2);
                        // dividing by bin size; bin is dkt times Jacobian
                        // k(=ik*dkt) times 2Pi in phi times the correct number
                        // of counts for an infinite lattice: area in bin
                        // divided by total area
                    }
                }
            }
        }
    }

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

    //  double dNdetaHadrons, dNdetaHadronsCut, dNdetaHadronsCut2;
    //  double dEdetaHadrons, dEdetaHadronsCut, dEdetaHadronsCut2;

    // compute hadrons using fragmentation function
    if (it == itmax && param->getWriteOutputs() == 3) {
        cout << " Hadronizing ... " << endl;
        double z, frac;
        double mypt, kt;
        int ik;
        const int steps = 6000;
        double dz = 0.95 / static_cast<double>(steps);
        double zValues[steps + 1];
        double zintegrand[steps + 1];
        // double Ezintegrand[steps+1];
        // double Lzintegrand[steps+1];
        // double Hzintegrand[steps+1];
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
                    // Ezintegrand[iz] = mypt * 1./(z*z) * Ng * kkp(7,1,z,kt);
                    // Lzintegrand[iz] = 1./(z*z) * Ng * kkp(7,1,z,kt/2.);
                    // Hzintegrand[iz] = 1./(z*z) * Ng * kkp(7,1,z,kt*2.);
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

                    // Ezintegrand[iz] =  mypt * 1./(z*z) * Ng *
                    // 	2. *
                    // (kkp(1,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_pion*m_pion/(mypt*mypt)))
                    // 	      +kkp(2,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_kaon*m_kaon/(mypt*mypt)))
                    // 	      +kkp(4,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_proton*m_proton/(mypt*mypt))));

                    // Lzintegrand[iz] =  1./(z*z) * Ng *
                    // 	2. *
                    // (kkp(1,1,z,kt/2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_pion*m_pion/(mypt*mypt)))
                    // 	      +kkp(2,1,z,kt/2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_kaon*m_kaon/(mypt*mypt)))
                    // 	      +kkp(4,1,z,kt/2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_proton*m_proton/(mypt*mypt))));

                    // Hzintegrand[iz] =  1./(z*z) * Ng *
                    // 	2. *
                    // (kkp(1,1,z,kt*2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_pion*m_pion/(mypt*mypt)))
                    // 	      +kkp(2,1,z,kt*2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_kaon*m_kaon/(mypt*mypt)))
                    // 	      +kkp(4,1,z,kt*2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_proton*m_proton/(mypt*mypt))));
                }
            }

            zValues[steps] = 1.;  // set exactly 1

            gsl_spline_init(zspline, zValues, zintegrand, steps + 1);
            Nhgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);

            //	  gsl_spline_init (zspline, zValues, Lzintegrand, steps+1);
            // NhLgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);

            // gsl_spline_init (zspline, zValues, Hzintegrand, steps+1);
            // NhHgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);
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

        cout << " done." << endl;

        // integrate over pT using gsl
        double pt[hbins + 1];
        double integrand[hbins + 1];
        double Eintegrand[hbins + 1];
        for (int ih = 0; ih <= hbins; ih++) {
            pt[ih] = ih * 20. / static_cast<double>(hbins);
            integrand[ih] = Nhgsl[ih] * pt[ih];
            Eintegrand[ih] = Nhgsl[ih] * pt[ih] * pt[ih];
        }

        gsl_interp_accel *ptacc = gsl_interp_accel_alloc();
        gsl_spline *ptspline = gsl_spline_alloc(gsl_interp_cspline, hbins + 1);
        gsl_spline_init(ptspline, pt, integrand, hbins + 1);
        //   dNdetaHadrons = 2*M_PI*gsl_spline_eval_integ(ptspline, 0.25, 19.,
        //   ptacc);
        // dNdetaHadronsCut = 2*M_PI*gsl_spline_eval_integ(ptspline, 3., 19.,
        // ptacc); dNdetaHadronsCut2 =
        // 2*M_PI*gsl_spline_eval_integ(ptspline, 6., 19., ptacc);

        gsl_spline_init(ptspline, pt, Eintegrand, hbins + 1);
        // dEdetaHadrons = 2*M_PI*gsl_spline_eval_integ(ptspline, 0.25, 19.,
        // ptacc); dEdetaHadronsCut =
        // 2*M_PI*gsl_spline_eval_integ(ptspline, 3., 19., ptacc);
        // dEdetaHadronsCut2 = 2*M_PI*gsl_spline_eval_integ(ptspline, 6., 19.,
        // ptacc);

        gsl_spline_free(ptspline);
        gsl_interp_accel_free(ptacc);
    }

    if (param->getUsePseudoRapidity() == 0 && param->getMPIRank() == 0) {
        cout << "dN/dy 1 = " << dNdeta << ", dE/dy 1 = " << dEdeta << endl;
        cout << "dN/dy 2 = " << dNdeta2 << ", dE/dy 2 = " << dEdeta2 << endl;
        cout << "gluon <p_T> = " << dEdeta / dNdeta << endl;
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
            cout << "dN/deta 1 = " << dNdeta << ", dE/deta 1 = " << dEdeta
                 << endl;
            cout << "dN/deta 2 = " << dNdeta2 << ", dE/deta 2 = " << dEdeta2
                 << endl;
            cout << "dN/deta_cut 1 = " << dNdetaCut << endl;
            cout << "dN/deta_cut 2 = " << dNdetaCut2 << endl;
            cout << "gluon <p_T> = " << dEdeta / dNdeta << endl;
        }
    }

    if (dNdeta == 0.) {
        cout << "No collision happened on rank " << param->getMPIRank()
             << ". Restarting with new random number..." << endl;
        for (int i = 0; i < N * N; i++) {
            delete E1[i];
        }

        delete[] E1;
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
    }

    for (int i = 0; i < N * N; i++) {
        delete E1[i];
    }

    delete[] E1;

    cout << " done." << endl;
    param->setSuccess(1);
    return 1;
}

int Evolution::multiplicitynkxky(
    Lattice *lat, Group *group, Parameters *param, int it) {
    const int N = param->getSize();
    const int Nc = param->getNc();
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
    vector<double> Nkxky(N * N, 0);

    stringstream strnkxky_name;
    strnkxky_name << "nkxky-t" << it * dtau * a << "-" << param->getEventId()
                  << ".dat";
    string nkxky_name;
    nkxky_name = strnkxky_name.str();

    stringstream strNpartdNdy_name;
    strNpartdNdy_name << "NpartdNdy-t" << it * dtau * a << "-"
                      << param->getEventId() << ".dat";
    string NpartdNdy_name;
    NpartdNdy_name = strNpartdNdy_name.str();

    stringstream strNpartdNdyH_name;
    strNpartdNdyH_name << "NpartdNdyHadrons-t" << it * dtau * a << "-"
                       << param->getEventId() << ".dat";
    string NpartdNdyH_name;
    NpartdNdyH_name = strNpartdNdyH_name.str();

    stringstream strmult_name;
    strmult_name << "multiplicity-t" << it * dtau * a << "-"
                 << param->getEventId() << ".dat";
    string mult_name;
    mult_name = strmult_name.str();

    stringstream strdNdy_name;
    strdNdy_name << "dNdy-t" << it * dtau * a << "-" << param->getEventId()
                 << ".dat";
    string dNdy_name;
    dNdy_name = strdNdy_name.str();

    cout << "Measuring multiplicity ... " << endl;

    // fix transverse Coulomb gauge
    GaugeFix gaugefix;

    double maxtime;
    if (param->getInverseQsForMaxTime() == 1) {
        maxtime = 1. / param->getAverageQs() * hbarc;
        cout << "maximal evolution time = " << maxtime << " fm" << endl;
    } else {
        maxtime = param->getMaxtime();  // maxtime is in fm
    }

    int itmax = static_cast<int>(floor(maxtime / (a * dtau) + 1e-10));

    gaugefix.FFTChi(fft, lat, group, param, 4000);
    // gauge is fixed

    Matrix **E1;
    E1 = new Matrix *[N * N];

    for (int i = 0; i < N * N; i++) {
        E1[i] = new Matrix(Nc, 0.);
    }

    double g2mu2A, g2mu2B, gfactor, alphas = 0., Qs = 0.;
    double c = param->getc();
    double muZero = param->getMuZero();

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;

            if (param->getRunningCoupling()) {
                if (pos > 0 && pos < (N - 1) * N + N - 1) {
                    g2mu2A = lat->cells[pos]->getg2mu2A();
                } else
                    g2mu2A = 0;

                if (pos > 0 && pos < (N - 1) * N + N - 1) {
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

                    gfactor = g * g / (4. * M_PI * alphas);
                }
            } else
                gfactor = 1.;

            if (param->getRunWithkt() == 0) {
                *E1[pos] = lat->cells[pos]->getE1()
                           * sqrt(gfactor);  // replace one of the 1/g in the
                                             // lattice E^i by the running one
            } else {
                *E1[pos] = lat->cells[pos]->getE1();
            }
        }
    }

    // do Fourier transforms
    fft->fftn(E1, E1, nn, 1);

    for (int ik = 0; ik < bins; ik++) {
        n[ik] = 0.;
        E[ik] = 0.;
        n2[ik] = 0.;
        counter[ik] = 0;
    }

    const int hbins = 2000;

    //  double Nh[hbins+1], Eh[hbins+1], Ehgsl[hbins+1], NhL[hbins+1],
    //  NhLgsl[hbins+1], NhH[hbins+1], NhHgsl[hbins+1];
    double Nhgsl[hbins + 1], Ng;
    // for (int ih=0; ih<=hbins; ih++)
    //   {
    //     Nh[ih]=0.;
    //     Eh[ih]=0.;
    //     NhL[ih]=0.;
    //     NhH[ih]=0.;
    //   }

    ofstream foutNkxky((nkxky_name).c_str(), std::ios::out);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;
            Nkxky[pos] = 0.;
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            nkt = 0.;
            pos = i * N + j;
            npos = (N - i) * N + (N - j);

            kx = 2. * M_PI
                 * (-0.5 + static_cast<double>(i) / static_cast<double>(N));
            ky = 2. * M_PI
                 * (-0.5 + static_cast<double>(j) / static_cast<double>(N));
            kt2 = 4.
                  * (sin(kx / 2.) * sin(kx / 2.)
                     + sin(ky / 2.) * sin(ky / 2.));  //
            omega2 = 4.
                     * (sin(kx / 2.) * sin(kx / 2.)
                        + sin(ky / 2.)
                              * sin(ky / 2.));  // lattice dispersion relation
                                                // (this is omega squared)

            // i=0 or j=0 have no negative k_T value available

            if (i != 0 && j != 0) {
                if (omega2 != 0) {
                    nkt = 2. / sqrt(omega2) / static_cast<double>(N * N)
                          * (g * g / ((it - 0.5) * dtau)
                             * ((((*E1[pos]) * (*E1[npos])).trace()).real()));
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
                        E[ik] += sqrt(omega2) * hbarc / a * nkt / dkt / 2 / M_PI
                                 / sqrt(kt2) * 2 * M_PI * sqrt(kt2) * dkt * N
                                 * N / M_PI / M_PI / 2. / 2.;
                        n2[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2);
                        // dividing by bin size; bin is dkt times Jacobian
                        // k(=ik*dkt) times 2Pi in phi times the correct number
                        // of counts for an infinite lattice: area in bin
                        // divided by total area
                        counter[ik] += 1;  // number of entries in n[ik]
                    }
                }
            }
            if (i != 0 && j != 0) {
                Nkxky[pos] = nkt * N * N / M_PI / M_PI / 2. / 2.;
            }
        }
    }

    /// -------- 2 ---------

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;

            if (param->getRunningCoupling()) {
                if (pos > 0 && pos < (N - 1) * N + N - 1) {
                    g2mu2A = lat->cells[pos]->getg2mu2A();
                } else
                    g2mu2A = 0;

                if (pos > 0 && pos < (N - 1) * N + N - 1) {
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

                    gfactor = g * g / (4. * M_PI * alphas);
                }
            } else
                gfactor = 1.;

            if (param->getRunWithkt() == 0) {
                *E1[pos] = lat->cells[pos]->getE2() * sqrt(gfactor);  // "
            } else {
                *E1[pos] = lat->cells[pos]->getE2();
            }
        }
    }

    fft->fftn(E1, E1, nn, 1);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            nkt = 0.;
            pos = i * N + j;
            npos = (N - i) * N + (N - j);

            kx = 2. * M_PI
                 * (-0.5 + static_cast<double>(i) / static_cast<double>(N));
            ky = 2. * M_PI
                 * (-0.5 + static_cast<double>(j) / static_cast<double>(N));
            kt2 = 4.
                  * (sin(kx / 2.) * sin(kx / 2.)
                     + sin(ky / 2.) * sin(ky / 2.));  //
            omega2 = 4.
                     * (sin(kx / 2.) * sin(kx / 2.)
                        + sin(ky / 2.)
                              * sin(ky / 2.));  // lattice dispersion relation
                                                // (this is omega squared)

            // i=0 or j=0 have no negative k_T value available

            if (i != 0 && j != 0) {
                if (omega2 != 0) {
                    nkt = 2. / sqrt(omega2) / static_cast<double>(N * N)
                          * (g * g / ((it - 0.5) * dtau)
                             * (((((*E1[pos]) * (*E1[npos])).trace()).real())));
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
                        E[ik] += sqrt(omega2) * hbarc / a * nkt / dkt / 2 / M_PI
                                 / sqrt(kt2) * 2 * M_PI * sqrt(kt2) * dkt * N
                                 * N / M_PI / M_PI / 2. / 2.;
                        n2[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2);
                        // dividing by bin size; bin is dkt times Jacobian
                        // k(=ik*dkt) times 2Pi in phi times the correct number
                        // of counts for an infinite lattice: area in bin
                        // divided by total area
                    }
                }
            }
            if (i != 0 && j != 0) {
                Nkxky[pos] += nkt * N * N / M_PI / M_PI / 2. / 2.;
            }
        }
    }

    /// ------3 --------

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;

            if (param->getRunningCoupling()) {
                if (pos > 0 && pos < (N - 1) * N + N - 1) {
                    g2mu2A = lat->cells[pos]->getg2mu2A();
                } else
                    g2mu2A = 0;

                if (pos > 0 && pos < (N - 1) * N + N - 1) {
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

                    gfactor = g * g / (4. * M_PI * alphas);
                }
            } else
                gfactor = 1.;

            if (param->getRunWithkt() == 0) {
                *E1[pos] = lat->cells[pos]->getpi()
                           * sqrt(gfactor);  // replace the only 1/g by the
                                             // running one (physical pi goes
                                             // like 1/g, like physical E^i)
            } else {
                *E1[pos] = lat->cells[pos]->getpi();
            }
        }
    }

    // do Fourier transforms
    fft->fftn(E1, E1, nn, 1);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            nkt = 0.;
            pos = i * N + j;
            npos = (N - i) * N + (N - j);

            kx = 2. * M_PI
                 * (-0.5 + static_cast<double>(i) / static_cast<double>(N));
            ky = 2. * M_PI
                 * (-0.5 + static_cast<double>(j) / static_cast<double>(N));
            kt2 = 4.
                  * (sin(kx / 2.) * sin(kx / 2.)
                     + sin(ky / 2.) * sin(ky / 2.));  //
            omega2 = 4.
                     * (sin(kx / 2.) * sin(kx / 2.)
                        + sin(ky / 2.)
                              * sin(ky / 2.));  // lattice dispersion relation
                                                // (this is omega squared)

            // i=0 or j=0 have no negative k_T value available

            if (i != 0 && j != 0) {
                if (omega2 != 0) {
                    nkt = 2. / sqrt(omega2) / static_cast<double>(N * N)
                          * (((it - 0.5) * dtau)
                             * ((((*E1[pos]) * (*E1[npos])).trace()).real()));
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
                        E[ik] += sqrt(omega2) * hbarc / a * nkt / dkt / 2 / M_PI
                                 / sqrt(kt2) * 2 * M_PI * sqrt(kt2) * dkt * N
                                 * N / M_PI / M_PI / 2. / 2.;
                        n2[ik] += nkt / dkt / 2 / M_PI / sqrt(kt2);
                        // dividing by bin size; bin is dkt times Jacobian
                        // k(=ik*dkt) times 2Pi in phi times the correct number
                        // of counts for an infinite lattice: area in bin
                        // divided by total area
                    }
                }
            }
            if (i != 0 && j != 0) {
                Nkxky[pos] += nkt * N * N / M_PI / M_PI / 2. / 2.;
            }
            if (param->getWriteOutputs() == 2)
                foutNkxky << 2. * sin(kx / 2.) / a * hbarc << " "
                          << 2. * sin(ky / 2.) / a * hbarc << " "
                          << Nkxky[pos] * a / hbarc * a / hbarc << "\n";
        }
        if (param->getWriteOutputs() == 2) foutNkxky << endl;
    }
    foutNkxky.close();

    double m, P;
    m = param->getJacobianm();                                // in GeV
    P = 0.13 + 0.32 * pow(param->getRoots() / 1000., 0.115);  // in GeV

    ofstream foutMult(mult_name.c_str(), std::ios::out);
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
            // integrate, gives a ik*dkt*2pi*dkt, in |eta|<2.4, 0.4 GeV p_T cut,
            // charged N_track (offline, factor 0.83)
        }

        // output dN/d^2k
        if (it > 0) {
            foutMult << it * dtau * a << " " << ik * dkt / a * hbarc << " "
                     << n[ik] * a / hbarc * a / hbarc << " "
                     << n2[ik] * a / hbarc * a / hbarc << " " << param->getTpp()
                     << " " << param->getb() << " " << param->getNpart()
                     << endl;
        }
    }

    foutMult.close();

    double dNdetaHadrons = 0;
    double dNdetaHadronsCut = 0;
    double dNdetaHadronsCut2 = 0;
    double dEdetaHadrons = 0;
    double dEdetaHadronsCut = 0;
    double dEdetaHadronsCut2 = 0;

    // compute hadrons using fragmentation function
    if (it == itmax && param->getWriteOutputs() == 3) {
        cout << " Hadronizing ... " << endl;
        double z, frac;
        double mypt, kt;
        int ik;
        const int steps = 6000;
        double dz = 0.95 / static_cast<double>(steps);
        double zValues[steps + 1];
        double zintegrand[steps + 1];
        // double Ezintegrand[steps+1];
        // double Lzintegrand[steps+1];
        // double Hzintegrand[steps+1];
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
                    // Ezintegrand[iz] = mypt * 1./(z*z) * Ng * kkp(7,1,z,kt);
                    // Lzintegrand[iz] = 1./(z*z) * Ng * kkp(7,1,z,kt/2.);
                    // Hzintegrand[iz] = 1./(z*z) * Ng * kkp(7,1,z,kt*2.);
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

                    // Ezintegrand[iz] =  mypt * 1./(z*z) * Ng *
                    // 	2. *
                    // (kkp(1,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_pion*m_pion/(mypt*mypt)))
                    // 	      +kkp(2,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_kaon*m_kaon/(mypt*mypt)))
                    // 	      +kkp(4,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_proton*m_proton/(mypt*mypt))));

                    // Lzintegrand[iz] =  1./(z*z) * Ng *
                    // 	2. *
                    // (kkp(1,1,z,kt/2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_pion*m_pion/(mypt*mypt)))
                    // 	      +kkp(2,1,z,kt/2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_kaon*m_kaon/(mypt*mypt)))
                    // 	      +kkp(4,1,z,kt/2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_proton*m_proton/(mypt*mypt))));

                    // Hzintegrand[iz] =  1./(z*z) * Ng *
                    // 	2. *
                    // (kkp(1,1,z,kt*2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_pion*m_pion/(mypt*mypt)))
                    // 	      +kkp(2,1,z,kt*2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_kaon*m_kaon/(mypt*mypt)))
                    // 	      +kkp(4,1,z,kt*2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m_proton*m_proton/(mypt*mypt))));
                }
            }

            zValues[steps] = 1.;  // set exactly 1

            gsl_spline_init(zspline, zValues, zintegrand, steps + 1);
            Nhgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);

            //	  gsl_spline_init (zspline, zValues, Lzintegrand, steps+1);
            // NhLgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);

            // gsl_spline_init (zspline, zValues, Hzintegrand, steps+1);
            // NhHgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);
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

        cout << " done." << endl;

        // integrate over pT using gsl
        double pt[hbins + 1];
        double integrand[hbins + 1];
        double Eintegrand[hbins + 1];
        for (int ih = 0; ih <= hbins; ih++) {
            pt[ih] = ih * 20. / static_cast<double>(hbins);
            integrand[ih] = Nhgsl[ih] * pt[ih];
            Eintegrand[ih] = Nhgsl[ih] * pt[ih] * pt[ih];
        }

        gsl_interp_accel *ptacc = gsl_interp_accel_alloc();
        gsl_spline *ptspline = gsl_spline_alloc(gsl_interp_cspline, hbins + 1);
        gsl_spline_init(ptspline, pt, integrand, hbins + 1);
        dNdetaHadrons =
            2 * M_PI * gsl_spline_eval_integ(ptspline, 0.25, 19., ptacc);
        dNdetaHadronsCut =
            2 * M_PI * gsl_spline_eval_integ(ptspline, 3., 19., ptacc);
        dNdetaHadronsCut2 =
            2 * M_PI * gsl_spline_eval_integ(ptspline, 6., 19., ptacc);

        gsl_spline_init(ptspline, pt, Eintegrand, hbins + 1);
        dEdetaHadrons =
            2 * M_PI * gsl_spline_eval_integ(ptspline, 0.25, 19., ptacc);
        dEdetaHadronsCut =
            2 * M_PI * gsl_spline_eval_integ(ptspline, 3., 19., ptacc);
        dEdetaHadronsCut2 =
            2 * M_PI * gsl_spline_eval_integ(ptspline, 6., 19., ptacc);

        gsl_spline_free(ptspline);
        gsl_interp_accel_free(ptacc);
    }

    if (param->getUsePseudoRapidity() == 0 && param->getMPIRank() == 0) {
        cout << "dN/dy 1 = " << dNdeta << ", dE/dy 1 = " << dEdeta << endl;
        cout << "dN/dy 2 = " << dNdeta2 << ", dE/dy 2 = " << dEdeta2 << endl;
        cout << "gluon <p_T> = " << dEdeta / dNdeta << endl;
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
            cout << "dN/deta 1 = " << dNdeta << ", dE/deta 1 = " << dEdeta
                 << endl;
            cout << "dN/deta 2 = " << dNdeta2 << ", dE/deta 2 = " << dEdeta2
                 << endl;
            cout << "dN/deta_cut 1 = " << dNdetaCut << endl;
            cout << "dN/deta_cut 2 = " << dNdetaCut2 << endl;
            cout << "gluon <p_T> = " << dEdeta / dNdeta << endl;
        }
    }

    if (dNdeta == 0.) {
        cout << "No collision happened on rank " << param->getMPIRank()
             << ". Restarting with new random number..." << endl;
        for (int i = 0; i < N * N; i++) {
            delete E1[i];
        }

        delete[] E1;
        return 0;
    }

    if (it == itmax) {
        cout << "hadron <p_T> = " << dEdetaHadrons / dNdetaHadrons << endl;
        cout << "Hadrons: dN/dy(p_T>250 MeV)=" << dNdetaHadrons
             << ", dE/dy(p_T>250 MeV)=" << dEdetaHadrons << endl;

        stringstream strmeanpt_name;
        strmeanpt_name << "meanpt" << param->getEventId() << ".dat";
        string meanpt_name;
        meanpt_name = strmeanpt_name.str();

        ofstream foutNch(meanpt_name.c_str(), std::ios::out);
        foutNch << dNdeta << " " << dEdeta / dNdeta << " " << dNdetaHadrons
                << " " << dEdetaHadrons / dNdetaHadrons << endl;
        foutNch.close();

        ofstream foutNN(NpartdNdy_name.c_str(), std::ios::app);
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

        ofstream foutNNH(NpartdNdyH_name.c_str(), std::ios::app);
        foutNNH << param->getNpart() << " " << dNdetaHadrons << " "
                << param->getTpp() << " " << param->getb() << " "
                << dEdetaHadrons << " " << param->getRandomSeed() << " "
                << "N/A"
                << " "
                << "N/A"
                << " "
                << "N/A"
                << " " << dNdetaHadronsCut << " " << dEdetaHadronsCut << " "
                << dNdetaHadronsCut2 << " " << dEdetaHadronsCut2 << " "
                << g * g
                       / (4. * M_PI * 4. * M_PI
                          / (9.
                             * log(pow(
                                 pow(muZero / 0.2, 2. / c)
                                     + pow(
                                         param->getRunWithThisFactorTimesQs()
                                             * param->getAverageQs() / 0.2,
                                         2. / c),
                                 c))))
                << endl;
        foutNNH.close();
    }

    for (int i = 0; i < N * N; i++) {
        delete E1[i];
    }

    delete[] E1;

    cout << " done." << endl;
    param->setSuccess(1);
    return 1;
}

int Evolution::correlations(
    Lattice *lat, Group *group, Parameters *param, int it) {
    const int N = param->getSize();
    const int Nc = param->getNc();
    int npos, pos;
    double L = param->getL();
    double a = L / N;  // lattice spacing in fm
    double kx, ky, kt2, omega2;
    double g = param->getg();
    int nn[2];
    nn[0] = N;
    nn[1] = N;
    double dtau = param->getdtau();
    double nkt, nkt1, nkt2, nkt3, nkt4, nkt5, nkt6;
    const int bins = 40;
    const int phiBins = 16;
    double n[bins][phiBins];  // |k_T|, phi array
    //  double n2[bins][phiBins]; // |k_T|, phi array
    double nkxky[N][N];  // kx, ky array
    double nk[bins];     //|k_T| array
    // double nNoMixedTerms[bins][phiBins]; //|k_T|, phi array
    // double nkNoMixedTerms[bins]; //|k_T| array
    // int counter[bins][phiBins];
    int counterk[bins];
    double dkt = 2.83 / static_cast<double>(bins);
    double dNdeta = 0.;
    double dNdetaNoMixedTerms = 0.;
    double dNdeta1 = 0.;
    double dNdeta2 = 0.;
    double dNdeta3 = 0.;
    double dNdeta4 = 0.;
    double dNdeta5 = 0.;
    double dNdeta6 = 0.;
    double anglePhi;
    double k;
    double deltaPhi = 2. * M_PI / static_cast<double>(phiBins);

    stringstream strCorr_name;
    strCorr_name << "Corr" << param->getEventId() << ".dat";
    string Corr_name;
    Corr_name = strCorr_name.str();

    stringstream strPhiMult_name;
    strPhiMult_name << "MultPhi" << param->getEventId() << ".dat";
    string PhiMult_name;
    PhiMult_name = strPhiMult_name.str();

    stringstream strPhiMultHad_name;
    strPhiMultHad_name << "MultPhiHadrons" << param->getEventId() << ".dat";
    string PhiMultHad_name;
    PhiMultHad_name = strPhiMultHad_name.str();

    cout << "Measuring multiplicity version 2... " << endl;

    // fix transverse Coulomb gauge
    GaugeFix gaugefix;

    double maxtime;
    if (param->getInverseQsForMaxTime() == 1) {
        maxtime = 1. / param->getAverageQs() * hbarc;
        cout << "maximal evolution time = " << maxtime << " fm" << endl;
    } else {
        maxtime = param->getMaxtime();  // maxtime is in fm
    }

    //  int itmax = static_cast<int>(floor(maxtime/(a*dtau)+1e-10));
    gaugefix.FFTChi(fft, lat, group, param, 4000);
    // gauge is fixed

    Matrix U1(Nc, 1.);
    Matrix U2(Nc, 1.);
    Matrix U1dag(Nc, 1.);
    Matrix U2dag(Nc, 1.);

    Matrix **A1;
    A1 = new Matrix *[N * N];
    Matrix **A2;
    A2 = new Matrix *[N * N];
    Matrix **phi;
    phi = new Matrix *[N * N];

    Matrix **E1;
    Matrix **E2;
    Matrix **pi;
    E1 = new Matrix *[N * N];
    E2 = new Matrix *[N * N];
    pi = new Matrix *[N * N];

    for (int i = 0; i < N * N; i++) {
        A1[i] = new Matrix(Nc, 0.);
        A2[i] = new Matrix(Nc, 0.);
        E1[i] = new Matrix(Nc, 0.);
        E2[i] = new Matrix(Nc, 0.);
        pi[i] = new Matrix(Nc, 0.);
        phi[i] = new Matrix(Nc, 0.);
    }

    // version that determines the exact log of U1 and U2:
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;
            U1 = lat->cells[pos]->getUx();
            U2 = lat->cells[pos]->getUy();

            U1.logm();
            U2.logm();

            *A1[pos] = complex<double>(0., -1.) * U1;
            *A2[pos] = complex<double>(0., -1.) * U2;
        }
    }

    double g2mu2A, g2mu2B, gfactor, alphas = 0., Qs = 0.;
    double c = param->getc();
    double muZero = param->getMuZero();

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;

            if (param->getRunningCoupling()) {
                if (pos > 0 && pos < (N - 1) * N + N - 1) {
                    g2mu2A = lat->cells[pos]->getg2mu2A();
                } else
                    g2mu2A = 0;

                if (pos > 0 && pos < (N - 1) * N + N - 1) {
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

                    gfactor = g * g / (4. * M_PI * alphas);
                }
            } else
                gfactor = 1.;

            if (param->getRunWithkt() == 0) {
                *E1[pos] = lat->cells[pos]->getE1()
                           * sqrt(gfactor);  // replace one of the 1/g in the
                                             // lattice E^i by the running one
                *E2[pos] = lat->cells[pos]->getE2() * sqrt(gfactor);  // "
                *pi[pos] = lat->cells[pos]->getpi()
                           * sqrt(gfactor);  // replace the only 1/g by the
                                             // running one (physical pi goes
                                             // like 1/g, like physical E^i)
                *A1[pos] = *A1[pos] * sqrt(gfactor);  //"
                *A2[pos] = *A2[pos] * sqrt(gfactor);  // "
                *phi[pos] = lat->cells[pos]->getphi()
                            * sqrt(gfactor);  // replace the only 1/g by the
                                              // running one (physical pi goes
                                              // like 1/g, like physical E^i)
            } else {
                *E1[pos] = lat->cells[pos]->getE1();
                *E2[pos] = lat->cells[pos]->getE2();
                *pi[pos] = lat->cells[pos]->getpi();
                *phi[pos] = lat->cells[pos]->getphi();
            }
        }
    }

    // do Fourier transforms

    fft->fftn(A1, A1, nn, 1);
    fft->fftn(A2, A2, nn, 1);
    fft->fftn(phi, phi, nn, 1);

    fft->fftn(E1, E1, nn, 1);
    fft->fftn(E2, E2, nn, 1);
    fft->fftn(pi, pi, nn, 1);

    for (int ik = 0; ik < bins; ik++) {
        nk[ik] = 0.;
        // nkNoMixedTerms[ik] = 0.;
        counterk[ik] = 0;
        for (int iphi = 0; iphi < phiBins; iphi++) {
            n[ik][iphi] = 0.;
            // n2[ik][iphi] = 0.;
            // nNoMixedTerms[ik][iphi] = 0.;
            //	  counter[ik][iphi]=0;
        }
    }

    // leave out the first cell to make it symmetric
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            pos = i * N + j;
            npos = (N - i) * N + (N - j);

            kx = 2. * M_PI
                 * (-0.5 + static_cast<double>(i) / static_cast<double>(N));
            ky = 2. * M_PI
                 * (-0.5 + static_cast<double>(j) / static_cast<double>(N));
            kt2 = 4.
                  * (sin(kx / 2.) * sin(kx / 2.)
                     + sin(ky / 2.) * sin(ky / 2.));  //
            omega2 = 4.
                     * (sin(kx / 2.) * sin(kx / 2.)
                        + sin(ky / 2.)
                              * sin(ky / 2.));  // lattice dispersion relation
                                                // (this is omega squared)

            // i=0 or j=0 have no negative k_T value available

            if (i != 0 && j != 0) {
                if (omega2 != 0) {
                    nkt1 =
                        1. / sqrt(omega2) / static_cast<double>(N * N)
                        * (1. / ((it - 0.5) * dtau)
                           * ((((*E1[pos]) * (*E1[npos])).trace()).real()
                              + (((*E2[pos]) * (*E2[npos])).trace()).real()));

                    nkt2 = 1. / sqrt(omega2) / static_cast<double>(N * N)
                           * (((it - 0.5) * dtau)
                              * ((((*pi[pos]) * (*pi[npos])).trace()).real()));

                    nkt3 =
                        sqrt(omega2) / static_cast<double>(N * N)
                        * ((it)*dtau
                           * ((((*A1[pos]) * (*A1[npos])).trace()).real()
                              + (((*A2[pos]) * (*A2[npos])).trace()).real()));

                    nkt4 =
                        sqrt(omega2) / static_cast<double>(N * N)
                        * (1. / ((it)*dtau)
                           * ((((*phi[pos]) * (*phi[npos])).trace()).real()));

                    nkt5 = 1. / static_cast<double>(N * N)
                           * (complex<double>(0., 1.)
                              * ((*E1[pos]) * (*A1[npos])
                                 - (*A1[pos]) * (*E1[npos])
                                 + (*E2[pos]) * (*A2[npos])
                                 - (*A2[pos]) * (*E2[npos]))
                                    .trace())
                                 .real();

                    nkt6 = 1. / static_cast<double>(N * N)
                           * (complex<double>(0., 1.)
                              * ((*pi[pos]) * (*phi[npos])
                                 - (*phi[pos]) * (*pi[npos]))
                                    .trace())
                                 .real();

                    if (param->getRunWithkt() == 1) {
                        nkt1 *=
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
                        nkt2 *=
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
                        nkt3 *=
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
                        nkt4 *=
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
                        nkt5 *=
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
                        nkt6 *=
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
                } else {
                    nkt1 = 0.;
                    nkt2 = 0.;
                    nkt3 = 0.;
                    nkt4 = 0.;
                    nkt5 = 0.;
                    nkt6 = 0.;
                }

                dNdeta1 += nkt1;
                dNdeta2 += nkt2;
                dNdeta3 += nkt3;
                dNdeta4 += nkt4;
                dNdeta5 += nkt5;
                dNdeta6 += nkt6;

                nkt = nkt1 + nkt2 + nkt3 + nkt4 + nkt5 + nkt6;

                dNdeta += nkt;  // total multiplicity

                nkxky[i][j] = nkt;
            }
        }
    }

    double latkx, latky;
    double fracX, fracY;
    for (int ik = 0; ik < bins; ik++) {
        k = ik * dkt;
        for (int iphi = 0; iphi < phiBins; iphi++) {
            anglePhi = deltaPhi * iphi;

            kx = k * cos(anglePhi);
            ky = k * sin(anglePhi);

            int i = floor(((kx) / 2 / M_PI + 0.5) * N + 1e-10);
            int j = floor(((ky) / 2 / M_PI + 0.5) * N + 1e-10);

            latkx =
                (2. * M_PI
                 * (-0.5 + static_cast<double>(i) / static_cast<double>(N)));
            latky =
                (2. * M_PI
                 * (-0.5 + static_cast<double>(j) / static_cast<double>(N)));

            fracX = (kx - latkx) / (2 * M_PI / static_cast<double>(N));
            fracY = (ky - latky) / (2 * M_PI / static_cast<double>(N));

            if (i + 1 < N && j + 1 < N)
                n[ik][iphi] = ((1. - fracX) * (1. - fracY) * nkxky[i][j]
                               + (fracX) * (1. - fracY) * nkxky[i + 1][j]
                               + (1. - fracX) * (fracY)*nkxky[i][j + 1]
                               + (fracX) * (fracY)*nkxky[i + 1][j + 1])
                              / 2. / M_PI / 2. / M_PI * N
                              * N;  // dkt/2/Pi/sqrt(kx*kx+ky*ky);
            else
                n[ik][iphi] = 0.;

            if (k == 0) n[ik][iphi] = 0.;
        }
    }

    //  double m,P;
    // m=param->getJacobianm(); // in GeV
    // P=0.13+0.32*pow(param->getRoots()/1000.,0.115); //in GeV
    double result, fullResult, fullResult2;
    fullResult = 0.;
    fullResult2 = 0.;

    for (int ik = 1; ik < bins; ik++) {
        if (counterk[ik] > 0) {
            nk[ik] = nk[ik] / static_cast<double>(counterk[ik]);
        }
        result = 0.;
        for (int iphi = 0; iphi < phiBins; iphi++) {
            result += n[ik][iphi] * deltaPhi;
        }
        fullResult += result * (ik)*dkt * dkt;
        fullResult2 += nk[ik] * 2. * M_PI * (ik + 0.5) * dkt * dkt;
    }
    cout << "N=" << dNdeta << ", k integrated N=" << fullResult2 << endl;
    cout << "N=" << dNdeta << ", k and phi integrated N=" << fullResult << endl;

    // output dN/d^2k
    ofstream foutPhiMult(PhiMult_name.c_str(), std::ios::out);
    if (it == 1) {
        foutPhiMult
            << "3"
            << " " << bins << " " << phiBins
            << endl;  // 3 is the number of times we read out. modify if needed.
    }
    for (int ik = 1; ik < bins; ik += 4) {
        for (int iphi = 0; iphi < phiBins; iphi++) {
            foutPhiMult
                << it * dtau * a << " " << ik * dkt / a * hbarc << " "
                << iphi * deltaPhi << " " << n[ik][iphi] * a / hbarc * a / hbarc
                << endl;  //<< " " << nNoMixedTerms[ik][iphi]*a/hbarc*a/hbarc <<
                          // endl;
        }
    }
    foutPhiMult.close();

    // compute hadrons using fragmentation function

    const int hbins = 40;
    double Nh[hbins + 1][phiBins], Ng;

    for (int ih = 0; ih <= hbins; ih++) {
        for (int iphi = 0; iphi < phiBins; iphi++) {
            Nh[ih][iphi] = 0.;
        }
    }
    double z, frac;
    double mypt, kt;
    int ik;
    int steps = 200;
    double dz = 0.95 / static_cast<double>(steps);

    for (int iphi = 0; iphi < phiBins; iphi++) {
        for (int ih = 0; ih <= hbins; ih++) {
            mypt = ih * dkt / a * hbarc;  //(10./static_cast<double>(hbins)); //
                                          // the hadron's p_T

            for (int iz = 0; iz < steps; iz++) {
                z = 0.05 + iz * dz;

                kt = mypt / z;  // the gluon's k_T

                ik = static_cast<int>(
                    floor(kt * a / hbarc / dkt - 0.5 + 0.00000001));

                frac = (kt - (ik + 0.5) * dkt / a * hbarc) / (dkt / a * hbarc);

                if (ik + 1 < bins && ik >= 0) {
                    Ng = ((1. - frac) * n[ik][iphi] + frac * n[ik + 1][iphi])
                         * a / hbarc * a
                         / hbarc;  // to make dN/d^2k_T fo k_T in GeV
                    if (kt > 2) Ng *= exp(-(kt - 2) * 0.5);
                } else
                    Ng = 0.;

                if (param->getUsePseudoRapidity() == 0) {
                    if (z == 0.05 || z == 1.) {
                        Nh[ih][iphi] +=
                            1. / (z * z) * Ng * kkp(7, 1, z, kt) * dz * 0.5;
                    } else {
                        Nh[ih][iphi] +=
                            1. / (z * z) * Ng * kkp(7, 1, z, kt) * dz;
                    }
                } else {
                    if (z == 0.05 || z == 1.) {
                        Nh[ih][iphi] +=
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
                                         + m_proton * m_proton
                                               / (mypt * mypt))))
                            * dz * 0.5;
                    } else {
                        Nh[ih][iphi] +=
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
                                         + m_proton * m_proton
                                               / (mypt * mypt))))
                            * dz;
                    }
                }
            }
        }
    }
    // output dN/d^2k
    ofstream foutPhiMultHad(PhiMultHad_name.c_str(), std::ios::out);
    if (it == 1) {
        foutPhiMultHad
            << "3"
            << " " << hbins << " " << phiBins
            << endl;  // 3 is the number of times we read out. modify if needed.
    }
    for (int ih = 0; ih < hbins; ih++) {
        for (int iphi = 0; iphi < phiBins; iphi++) {
            foutPhiMultHad << it * dtau * a << " " << ih * dkt / a * hbarc
                           << " " << iphi * deltaPhi << " " << Nh[ih][iphi]
                           << endl;
        }
    }
    foutPhiMultHad.close();

    ofstream foutCorr(Corr_name.c_str(), std::ios::out);
    foutCorr << it * dtau * a << " " << dNdeta1 << " " << dNdeta2 << " "
             << dNdeta3 << " " << dNdeta4 << " " << dNdeta5 << " " << dNdeta6
             << " " << dNdeta << " " << dNdetaNoMixedTerms << endl;
    foutCorr.close();

    for (int i = 0; i < N * N; i++) {
        delete E1[i];
        delete E2[i];
        delete pi[i];
        delete A1[i];
        delete A2[i];
        delete phi[i];
    }

    delete[] E1;
    delete[] E2;
    delete[] pi;
    delete[] A1;
    delete[] A2;
    delete[] phi;

    cout << " done." << endl;
    param->setSuccess(1);
    return 1;
}

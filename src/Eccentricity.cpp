// Eccentricity.cpp is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Eccentricity.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "Instrumentation.h"
#include "PhysConst.h"
#include "RunningCoupling.h"

using PhysConst::hbarc;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;

namespace {
/// Result of computeRotatedAnisotropy(): the rotated- and
/// unrotated-frame \f$T^{xx}-T^{yy}\f$ spatial anisotropy sums.
struct AnisotropyResult {
    /// \f$\sum (T^{xx}_{\text{rot}}-T^{yy}_{\text{rot}})\f$ at the
    /// sampled angle \f$\Psi\f$.
    double num;
    /// \f$\sum (T^{xx}_{\text{rot}}+T^{yy}_{\text{rot}})\f$ at the
    /// sampled angle \f$\Psi\f$.
    double den;
    /// \f$\sum (T^{xx}-T^{yy})\f$ in the unrotated frame.
    double num2;
    /// \f$\sum (T^{xx}+T^{yy})\f$ in the unrotated frame.
    double den2;
};

/**
 * Eccentricity::compute()'s `doAniso==1` branch samples this at
 * ten values of \p Psi (previously ten copy-pasted ~25-line blocks,
 * differing only in \p Psi): sums \f$T^{xx}-T^{yy}\f$ and
 * \f$T^{xx}+T^{yy}\f$ over the whole lattice, both in the unrotated
 * frame and after rotating \f$T^{xx}\f$/\f$T^{xy}\f$/\f$T^{yy}\f$ by
 * \p Psi.
 * \param[in] lat Lattice to read `Txx`/`Txy`/`Tyy` from.
 * \param[in] N Lattice side length.
 * \param[in] Psi Rotation angle [rad].
 * \return The summed rotated- and unrotated-frame anisotropy
 * numerators/denominators.
 */
AnisotropyResult computeRotatedAnisotropy(Lattice *lat, int N, double Psi) {
    double num = 0., den = 0., num2 = 0., den2 = 0.;
    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            int pos = lat->positionFromXY(ix, iy);

            double TxxRot = cos(Psi)
                                * (cos(Psi) * lat->cells[pos]->getTxx()
                                   - sin(Psi) * lat->cells[pos]->getTxy())
                            - sin(Psi)
                                  * (cos(Psi) * lat->cells[pos]->getTxy()
                                     - sin(Psi) * lat->cells[pos]->getTyy());
            double TyyRot = sin(Psi)
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
    return {num, den, num2, den2};
}
}  // namespace

void Eccentricity::compute(
    Lattice *lat, Parameters *param, int it, double cutoff, int doAniso) {
    IPG_PROFILE_SCOPE("observables.eccentricity");
    stringstream strecc_name;
    strecc_name << "eccentricities" << param->event.eventId << ".dat";
    string ecc_name;
    ecc_name = strecc_name.str();

    // cutoff on energy density is 'cutoff' times Lambda_QCD^4
    int N = param->lattice.size;
    int pos;
    double rA, phiA, x, y;
    double L = param->lattice.L;
    double a = L / N;  // lattice spacing in fm
    double eccentricity1, eccentricity2, eccentricity3, eccentricity4,
        eccentricity5, eccentricity6;
    double avcos, avsin, avcos1, avsin1, avcos3, avsin3, avrSq, avxSq, avySq,
        avr1, avr3, avcos4, avsin4, avr4, avcos5, avsin5, avr5, avcos6, avsin6,
        avr6;
    double Rbar;
    double Psi1, Psi2, Psi3, Psi4, Psi5, Psi6;
    double maxEps = 0;
    double g = param->coupling.g;

    double g2mu2A, g2mu2B, gfactor, alphas = 0., Qs = 0.;
    double c = param->coupling.c;
    double muZero = param->coupling.mu0;

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
            pos = lat->positionFromXY(ix, iy);
            maxEps = std::max(lat->cells[pos]->getEpsilon(), maxEps);
        }
    }

    // first shift to the center
    for (int ix = 0; ix < N; ix++) {
        x = -L / 2. + a * ix;
        for (int iy = 0; iy < N; iy++) {
            y = -L / 2. + a * iy;
            pos = lat->positionFromXY(ix, iy);

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
                avgQs2AQs2B += g2mu2A * param->colorCharge.QsMuRatio
                               * param->colorCharge.QsMuRatio * g2mu2B
                               * param->colorCharge.QsMuRatio
                               * param->colorCharge.QsMuRatio / a / a / a / a;
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
    param->event.area = area;

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
            pos = lat->positionFromXY(ix, iy);
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
            pos = lat->positionFromXY(ix, iy);

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
    if (it == 1) param->event.psi = Psi2;

    if (doAniso == 0) {
        ofstream foutEcc(ecc_name.c_str(), std::ios::app);
        foutEcc << it * a * param->run.dtau << " " << eccentricity1 << " "
                << Psi1 << " " << eccentricity2 << " " << Psi2 << " "
                << eccentricity3 << " " << Psi3 << " " << eccentricity4 << " "
                << Psi4 << " " << eccentricity5 << " " << Psi5 << " "
                << eccentricity6 << " " << Psi6 << " " << cutoff << " "
                << sqrt(avrSq) << " " << maxX << " " << maxY << " "
                << param->event.b << " " << param->event.Tpp << " "
                << param->event.area << " " << Rbar << " " << avgeden << " "
                << avgQs2AQs2B * hbarc << endl;
        foutEcc.close();
    }

    if (doAniso == 1) {
        stringstream straniso_name;
        straniso_name << "anisotropy" << param->event.eventId << ".dat";
        string aniso_name;
        aniso_name = straniso_name.str();

        ofstream foutAniso(aniso_name.c_str(), std::ios::app);

        double ux, uy, PsiU;
        double unum = 0., uden = 0.;

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                pos = lat->positionFromXY(ix, iy);
                ux = lat->cells[pos]->getux();
                uy = lat->cells[pos]->getuy();
                unum += sqrt(ux * ux + uy * uy) * sin(2. * atan2(uy, ux));
                uden += sqrt(ux * ux + uy * uy) * cos(2. * atan2(uy, ux));
            }
        }

        PsiU = atan2(unum, uden) / 2.;

        foutAniso << "Psi2=" << Psi2 << ", cos(Psi2)=" << cos(Psi2)
                  << ", sin(Psi2)=" << sin(Psi2) << endl;
        foutAniso << "PsiU=" << PsiU << ", cos(PsiU)=" << cos(PsiU)
                  << ", sin(PsiU)=" << sin(PsiU) << endl;

        // Sample the rotated-tensor anisotropy at Psi = PsiU + k*Pi/8 for
        // k=0..9 (k=0: Psi = PsiU;  // param->event.psi;//-Pi/2.;).
        for (int k = 0; k < 10; ++k) {
            const double Psi = PsiU + static_cast<double>(k) * M_PI / 8.;
            const AnisotropyResult result =
                computeRotatedAnisotropy(lat, N, Psi);
            foutAniso << it * a * param->run.dtau << " "
                      << result.num / result.den << " "
                      << result.num2 / result.den2 << " angle=" << Psi << endl;
        }

        foutAniso.close();
    }
}

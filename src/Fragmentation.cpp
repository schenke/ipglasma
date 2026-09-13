// Fragmentation.cpp is part of the IP-Glasma solver.
// Copyright (C) 2013 Bjoern Schenke.
#include "Fragmentation.h"

#include <cmath>

#include "PrettyOstream.h"

namespace Fragmentation {
//**************************************************************************
// Fragmentation functions

// using KKP from http://www.desy.de/~poetter/kkp.html
// rewritten into C++ by Bjoern Schenke 2013
//**************************************************************************

//=====================================================================
//
//     ------------------------------------------------------------
//     Fragmentation functions for: Pions, Kaons, Protons, Neutrons
//            (includes mass-threshholds for c and b quarks)
//
//     Reference: B.A.Kniehl, G.Kramer, B.Potter, NPB582 (2000) 514
//     ------------------------------------------------------------
//
//     ih, iset, x, qs are input; dh is output.
//     ih   = 1 : (pi^+ + pi^-)  /2
//     ih   = 2 : (K^+ + K^-)    /2
//     ih   = 3 : (K^0 + K^0_bar)/2
//     ih   = 4 : (p + p_bar)    /2
//     ih   = 5 : (pi^0)
//     ih   = 6 : (n + n_bar)    /2
//     ih   = 7 : (h^+ + h^-)         [as sum of pions, kaons and protons]
//
//     iset = 0 : LO
//     iset = 1 : NLO
//
//     x    = longitudinal-momentum fraction
//     qs   = fragmentation scale (in GeV)
//
//     Parton label:
//     0    1    2    3    4    5    6    7    8     9    10
//     g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
//
//     Lambda_QCD (in GeV):
//     0.088 in LO
//     0.213 in NLO
//
//=====================================================================

namespace {

// One (b1,b2,b3,a1..a11) coefficient set for the KKP parametrization's
// fitted functional form (see evalKkpFit); extraA3rd is the LO
// proton-gluon fragmentation's unique extra +extraA3rd*S^3 term inside
// the last factor (0 for every other species/flavor).
struct KkpFitCoeffs {
    double b1, b2, b3;
    double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11;
    double extraA3rd = 0.;
};

// One coefficient set per parton flavor (gluon/up/down-or-analog/charm/
// bottom) for a fixed hadron species (pion/kaon/proton) and iset
// (LO or NLO).
struct KkpSpeciesFits {
    KkpFitCoeffs pionG, pionU, pionS, pionC, pionB;
    KkpFitCoeffs kaonG, kaonU, kaonD, kaonC, kaonB;
    KkpFitCoeffs protonG, protonU, protonS, protonC, protonB;
};

const KkpSpeciesFits kLoFits = {
    /*pionG*/
    {6.04510, -0.71378, 2.92133, -6.61523, -1.64978, 2.68223, 0.14705, -1.08423,
     -0.43182, 1.48429, 1.32887, -1.78696, 0.23086, -0.29182},
    /*pionU*/
    {0.54610, -1.46616, 1.01864, -0.22946, -0.22594, 0.21119, -0.45404,
     -0.12684, 0.27646, 0.95367, -1.09835, 0.74657, -0.01877, 0.02949},
    /*pionS*/
    {22.2815, 0.12732, 6.13697, -20.8125, -11.5725, 15.5372, 0.23075, -2.71424,
     1.72456, 2.18849, -5.04475, 3.29117, 0.09044, -0.07589},
    /*pionC*/
    {8.75500, -0.38611, 5.61846, -9.32277, 1.80600, 2.02179, -0.41190, -0.48496,
     0.42525, 0.74035, -0.64929, 0.66788, 0.06652, -0.05531},
    /*pionB*/
    {0.31147, -1.92993, 3.47086, -0.19319, -0.10487, 0.18824, -0.44692,
     -0.08271, 0.30441, 0.79775, -0.28091, 0.39504, -0.04887, 0.03212},
    /*kaonG*/
    {0.02862, -2.94091, 2.73474, -0.02113, 0.00389, 0.00901, 0.66881, -0.29670,
     0.20574, -0.58222, 0.04329, 0.78033, 0.03586, -0.01220},
    /*kaonU*/
    {0.25937, -0.61925, 0.85946, -0.10502, 0.00572, -0.00269, 0.09956, 0.07389,
     -0.00070, 0.57965, 0.26397, -0.12764, 0.15303, 0.14807},
    /*kaonD*/
    {5.38115, -0.00321, 3.07632, -3.05084, -1.10056, 1.31207, -0.25889,
     -0.18494, 0.13994, 1.13745, -0.90413, 0.56581, 0.05141, -0.00697},
    /*kaonC*/
    {5.18266, -0.17751, 4.30306, -3.48519, -1.00982, 1.17996, 0.02309, -0.61327,
     -0.03532, 1.00547, -0.51779, 0.20683, 0.13514, -0.17778},
    /*kaonB*/
    {1.57044, -0.84143, 6.01488, -1.78340, 0.57100, 0.15469, -0.43448, -0.05314,
     -0.36621, 0.72953, -0.64433, 0.92351, 0.01024, -0.06160},
    /*protonG*/
    {0.73953, -0.76986, 7.69079, -1.64519, 1.01189, -0.10175, -3.58787, 13.8025,
     -13.8902, -2.84470, -0.36719, -2.21825, 1.26515, -1.96117,
     /*extraA3rd=*/0.54769},
    /*protonU*/
    {0.40211, -0.85973, 2.80160, -0.21633, -0.07045, 0.07831, 0.13987, -0.82412,
     0.43114, 0.78923, -0.05344, 0.01460, 0.05198, -0.04623},
    /*protonS*/
    {4.07885, -0.09735, 4.99191, -2.97392, -0.92973, 1.23517, 0.25834, -1.52246,
     0.77060, 1.14379, -0.85320, 0.45607, 0.07174, -0.08321},
    /*protonC*/
    {0.11061, -1.54340, 2.20681, -0.07726, 0.05422, -0.03364, -0.20804, 0.29038,
     -0.23662, 0.62274, 0.29713, -0.21861, 0.00831, 0.00065},
    /*protonB*/
    {40.0971, 0.74249, 12.3729, -123.531, 128.666, -29.1808, -1.29639, -3.65003,
     3.05340, -1.04932, 0.34662, -1.34412, -0.04290, -0.30359},
};

const KkpSpeciesFits kNloFits = {
    /*pionG*/
    {3.73331, -0.74159, 2.33092, -3.16946, -0.47683, 0.70270, -0.51377,
     -0.19705, -0.17917, 2.03394, -0.50764, -0.08565, 0.09466, -0.10222},
    /*pionU*/
    {0.44809, -1.47598, 0.91338, -0.13828, -0.06951, 0.01354, -0.30498,
     -0.01863, -0.12529, 0.64145, 0.07270, -0.16989, 0.07396, -0.07757},
    /*pionS*/
    {16.5987, 0.13345, 5.89903, -18.3856, 2.44225, 2.13225, 0.22712, -0.83625,
     0.38526, -0.16911, 0.59886, -0.25630, -0.18619, 0.87362},
    /*pionC*/
    {6.17173, -0.53618, 5.60108, -4.82450, -1.30844, 1.95527, -0.27879,
     -0.51337, 0.10900, 0.83571, -1.15141, 0.77027, 0.09268, -0.11267},
    /*pionB*/
    {0.25944, -1.98713, 3.52857, -0.11449, 0.03733, -0.18028, -0.35858, 0.22277,
     -0.66413, 0.72303, 0.46260, -0.99235, -0.02701, -0.02089},
    /*kaonG*/
    {0.23140, -1.36400, 1.79761, -0.33644, 0.16204, -0.02598, 0.97182, -0.02908,
     -0.43195, 1.57116, 0.71847, -0.68331, 0.36906, 2.39060},
    /*kaonU*/
    {0.17806, -0.53733, 0.75940, -0.10988, -0.02524, 0.03142, -0.60058, 0.07863,
     0.13276, 0.61356, -0.43886, 0.23942, 0.10742, 0.12800},
    /*kaonD*/
    {4.96269, 0.05562, 2.79926, 1.54098, -9.06376, 4.94791, 1.88660, -2.94350,
     1.04227, 3.02991, -4.14807, 1.91494, 0.85450, -0.61016},
    /*kaonC*/
    {4.25954, -0.24144, 4.21265, -5.44309, 6.11031, -3.13973, -1.07757, 1.52364,
     -0.74308, 0.25590, 0.98423, -0.52839, -0.04000, 0.08695},
    /*kaonB*/
    {1.32443, -0.88351, 6.15221, -1.41156, -0.04809, 0.79066, -0.44818,
     -0.60073, 0.45526, 0.46679, -0.50792, 0.67006, -0.00477, -0.05503},
    /*protonG*/
    {1.56255, 0.01567, 3.57583, -1.48158, -0.39439, 0.51249, -2.16232, 2.47127,
     -0.93259, 3.33958, -3.05265, 1.21042, -0.84816, 1.23583},
    /*protonU*/
    {1.25946, 0.07124, 4.12795, -1.17505, 0.37550, -0.01416, -0.29533, -0.24540,
     0.16543, 0.98867, -0.46846, 0.20750, 0.18957, -0.01116},
    /*protonS*/
    {4.01135, 0.17258, 5.20766, 8.67124, -22.7888, 11.4720, 4.57608, -9.64835,
     4.61792, 7.25144, -12.6313, 6.07314, 0.16931, -0.09541},
    /*protonC*/
    {0.08250, -1.61290, 2.01255, -0.04512, -0.00565, 0.00900, -0.38012,
     -0.06840, 0.08888, 0.63782, -0.14146, 0.06083, -0.02958, 0.01130},
    /*protonB*/
    {24.2916, 0.57939, 12.1207, -88.3524, 93.1056, -17.4089, -0.80783, -5.07200,
     -2.45377, -3.27370, 1.21188, -5.50374, 0.14628, -0.78634},
};

double evalKkpFit(const KkpFitCoeffs &c, double s, double x) {
    return (c.b1 + c.a1 * s + c.a2 * s * s + c.a3 * s * s * s)
           * pow(x, c.b2 + c.a4 * s + c.a5 * s * s + c.a6 * s * s * s)
           * pow((1. - x), c.b3 + c.a7 * s + c.a8 * s * s + c.a9 * s * s * s)
           * (1. + (c.a10 * s + c.a11 * s * s + c.extraA3rd * s * s * s) / x);
}

}  // namespace

double kkp(int ih, int iset, double x, double qs) {
    // --- Mass-thresholds:
    const double rmcc = 2.9788;
    const double rmbb = 9.46037;
    // --- Q_0 (in GeV):
    const double q0 = sqrt(2.);
    if (qs < q0) qs = q0;

    double rlam;
    if (iset == 0) {
        rlam = 0.088;
    } else {
        if (iset != 1) {
            PrettyOstream messager;
            messager.warning(
                "[Fragmentation::kkp]: iset should be 0 (LO) or 1 (NLO); "
                "got a different value, proceeding with the NLO "
                "fragmentation functions.");
        }
        rlam = 0.213;
    }
    const double s =
        log(log(qs * qs / (rlam * rlam)) / log(q0 * q0 / (rlam * rlam)));
    const double sc =
        log(log(qs * qs / (rlam * rlam)) / log(rmcc * rmcc / (rlam * rlam)));
    const double sb =
        log(log(qs * qs / (rlam * rlam)) / log(rmbb * rmbb / (rlam * rlam)));

    const KkpSpeciesFits &fits = (iset == 0) ? kLoFits : kNloFits;

    double dpg = evalKkpFit(fits.pionG, s, x);
    double dpu = evalKkpFit(fits.pionU, s, x);
    double dps = evalKkpFit(fits.pionS, s, x);
    double dpc = evalKkpFit(fits.pionC, sc, x);
    double dpb = evalKkpFit(fits.pionB, sb, x);

    double dkg = evalKkpFit(fits.kaonG, s, x);
    double dku = evalKkpFit(fits.kaonU, s, x);
    double dkd = evalKkpFit(fits.kaonD, s, x);
    double dkc = evalKkpFit(fits.kaonC, sc, x);
    double dkb = evalKkpFit(fits.kaonB, sb, x);

    double dprg = evalKkpFit(fits.protonG, s, x);
    double dpru = evalKkpFit(fits.protonU, s, x);
    double dprs = evalKkpFit(fits.protonS, s, x);
    double dprc = evalKkpFit(fits.protonC, sc, x);
    double dprb = evalKkpFit(fits.protonB, sb, x);

    double dks;
    double dh[11];

    // --- Evaluate different contributions
    double dpd = dpu;
    dks = dku;
    double dprd = 0.5 * dpru;
    if (qs < rmbb) {
        dpb = 0.;
        dkb = 0.;
        dprb = 0.;
    }
    if (qs < rmcc) {
        dpc = 0.;
        dkc = 0.;
        dprc = 0.;
    }
    if (ih == 1) {
        dh[0] = dpg / 2.;
        dh[1] = dpu / 2.;
        dh[2] = dpu / 2.;
        dh[3] = dpd / 2.;
        dh[4] = dpd / 2.;
        dh[5] = dps / 2.;
        dh[6] = dps / 2.;
        dh[7] = dpc / 2.;
        dh[8] = dpc / 2.;
        dh[9] = dpb / 2.;
        dh[10] = dpb / 2.;
    } else if (ih == 2) {
        dh[0] = dkg / 2.;
        dh[1] = dku / 2.;
        dh[2] = dku / 2.;
        dh[3] = dkd / 2.;
        dh[4] = dkd / 2.;
        dh[5] = dks / 2.;
        dh[6] = dks / 2.;
        dh[7] = dkc / 2.;
        dh[8] = dkc / 2.;
        dh[9] = dkb / 2.;
        dh[10] = dkb / 2.;
    } else if (ih == 3) {
        dh[0] = dkg / 2.;
        dh[1] = dkd / 2.;
        dh[2] = dkd / 2.;
        dh[3] = dku / 2.;
        dh[4] = dku / 2.;
        dh[5] = dks / 2.;
        dh[6] = dks / 2.;
        dh[7] = dkc / 2.;
        dh[8] = dkc / 2.;
        dh[9] = dkb / 2.;
        dh[10] = dkb / 2.;
    } else if (ih == 4) {
        dh[0] = dprg / 2.;
        dh[1] = dpru / 2.;
        dh[2] = dpru / 2.;
        dh[3] = dprd / 2.;
        dh[4] = dprd / 2.;
        dh[5] = dprs / 2.;
        dh[6] = dprs / 2.;
        dh[7] = dprc / 2.;
        dh[8] = dprc / 2.;
        dh[9] = dprb / 2.;
        dh[10] = dprb / 2.;
    } else if (ih == 5) {
        dh[0] = dpg / 2.;
        dh[1] = dpu / 2.;
        dh[2] = dpu / 2.;
        dh[3] = dpd / 2.;
        dh[4] = dpd / 2.;
        dh[5] = dps / 2.;
        dh[6] = dps / 2.;
        dh[7] = dpc / 2.;
        dh[8] = dpc / 2.;
        dh[9] = dpb / 2.;
        dh[10] = dpb / 2.;
    } else if (ih == 6) {
        dh[0] = dprg / 2.;
        dh[1] = dpru / 4.;
        dh[2] = dpru / 4.;
        dh[3] = dprd;
        dh[4] = dprd;
        dh[5] = dprs / 2.;
        dh[6] = dprs / 2.;
        dh[7] = dprc / 2.;
        dh[8] = dprc / 2.;
        dh[9] = dprb / 2.;
        dh[10] = dprb / 2.;
    } else {
        dh[0] = dpg + dkg + dprg;
        dh[1] = dpu + dku + dpru;
        dh[2] = dpu + dku + dpru;
        dh[3] = dpd + dkd + dprd;
        dh[4] = dpd + dkd + dprd;
        dh[5] = dps + dks + dprs;
        dh[6] = dps + dks + dprs;
        dh[7] = dpc + dkc + dprc;
        dh[8] = dpc + dkc + dprc;
        dh[9] = dpb + dkb + dprb;
        dh[10] = dpb + dkb + dprb;
    }

    return dh[0];
}

}  // namespace Fragmentation

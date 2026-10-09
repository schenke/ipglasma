// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include "NuclearQsTable.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

#include "Instrumentation.h"

using std::ifstream;
using std::string;

void NuclearQsTable::read(const std::string &fileName) {
    IPG_PROFILE_SCOPE("initialization.read_qs_table");
    // steps in qs0 and Y in the file
    string dummy;
    string T, Qs;
    // open file

    messager_ << "[NuclearQsTable::read]: Reading Q_s(sum(T_p),y) from file ";
    messager_ << fileName << " ... ";
    messager_.flush("info");

    ifstream fin;
    fin.open((fileName).c_str());
    if (fin) {
        for (int iT = 0; iT < nT_; iT++) {
            for (int iy = 0; iy < nY_; iy++) {
                // checking the extraction (not eof() beforehand) also catches
                // a truncated last line
                if (fin >> dummy >> T >> Qs) {
                    T_[iT] = atof(T.c_str());
                    Qs2_[iT][iy] = atof(Qs.c_str());
                } else {
                    messager_ << "[NuclearQsTable::read]: Could not read "
                                 "entry "
                              << iT * nY_ + iy + 1
                              << " (the file ends early or is malformed) -- "
                                 "did the Q_s table file change? Exiting.";
                    messager_.flush("error");
                    exit(1);
                }
            }
        }
        fin.close();
    } else {
        messager_ << "[NuclearQsTable::read]: File " << fileName
                  << " does not exist. Exiting.";
        messager_.flush("error");
        exit(1);
    }
}

double NuclearQsTable::qs2(double T, double x) const {
    double value, fracy, fracT, QsYdown, QsYup;
    int posy, check = 0;
    fracy = 0.;

    // x > 0.01 only occurs with a fixed x, which the input validation
    // rejects (Parameters::validationErrors()); also catches NaN.
    if (!(x > 0. && x <= 0.01)) {
        // qs2() is called from inside an omp parallel for loop
        // (Init::setColorChargeDensity), so use a fresh, stack-local
        // instance rather than sharing messager_, which is not thread-safe.
        PrettyOstream localMessager;
        localMessager << "[NuclearQsTable::qs2]: x=" << x
                      << " is outside the tabulated range (0 < x <= 0.01). "
                         "Exiting.";
        localMessager.flush("error");
        exit(1);
    }
    double y = std::log(0.01 / x);
    // Above the largest tabulated rapidity (dilute cells with a fluctuating
    // x), use the last tabulated rapidity, as for T above the table.
    const double yMax = (nY_ - 1) * deltaY_;
    if (y > yMax) {
        if (!warnedAboveYMax_.exchange(true)) {
            // Local instance for the same omp thread-safety reason as above.
            PrettyOstream localMessager;
            localMessager << "[NuclearQsTable::qs2]: x=" << x << " (y=" << y
                          << ") is below the tabulated range (min x="
                          << 0.01 * std::exp(-yMax) << ", max y=" << yMax
                          << "); clamping to the minimal tabulated x (further "
                             "occurrences are not reported).";
            localMessager.flush("warning");
        }
        y = yMax;
    }
    posy = static_cast<int>(floor(y / deltaY_ + 0.0000001));
    // at the last tabulated y, interpolate from the bin below with fracy = 1
    // so that posy + 1 stays inside the table
    if (posy > nY_ - 2) posy = nY_ - 2;

    // at the last tabulated T there is no bin above to interpolate in
    if (T == T_[nT_ - 1]) {
        fracy = (y - static_cast<double>(posy) * deltaY_) / deltaY_;
        return fracy * Qs2_[nT_ - 1][posy + 1]
               + (1. - fracy) * Qs2_[nT_ - 1][posy];
    }

    //  if ( T > Qs2_[nT_-1][nY_-1] )
    if (T > T_[nT_ - 1]) {
        // Local instance for the same omp thread-safety reason as above.
        PrettyOstream localMessager;
        localMessager << "[NuclearQsTable::qs2]: T=" << T
                      << " exceeds the tabulated range (max T=" << T_[nT_ - 1]
                      << "); clamping to the maximal tabulated T.";
        localMessager.flush("warning");
        check = 1;
        fracy = (y - static_cast<double>(posy) * deltaY_) / deltaY_;
        QsYdown = (Qs2_[nT_ - 1][posy]);
        QsYup = (Qs2_[nT_ - 1][posy + 1]);
        value = (fracy * QsYup + (1. - fracy) * QsYdown);  //*hbarc*hbarc;
        return value;
    }

    if (T < T_[0]) {
        check = 1;
        return 0.;
    }

    for (int iT = 0; iT < nT_ - 1; iT++) {
        if (T >= T_[iT] && T < T_[iT + 1]) {
            fracT = (T - T_[iT]) / (T_[iT + 1] - T_[iT]);
            fracy = (y - static_cast<double>(posy) * deltaY_) / deltaY_;

            QsYdown = (fracT) * (Qs2_[iT + 1][posy])
                      + (1. - fracT) * (Qs2_[iT][posy]);
            QsYup = (fracT) * (Qs2_[iT + 1][posy + 1])
                    + (1. - fracT) * (Qs2_[iT][posy + 1]);
            value = (fracy * QsYup + (1. - fracy) * QsYdown);  //*hbarc*hbarc;

            check++;
            continue;
        }
    }

    if (check != 1) {
        // Local instance: same omp thread-safety reason as above.
        PrettyOstream localMessager;
        localMessager << "[NuclearQsTable::qs2]: could not uniquely "
                         "determine Qs^2 (check="
                      << check << ", T=" << T
                      << "); falling back to the maximal tabulated T_p.";
        localMessager.flush("warning");
        value =
            (fracy * Qs2_[nT_ - 1][posy + 1]
             + (1. - fracy) * Qs2_[nT_ - 1][posy]);
    }

    return value;
}

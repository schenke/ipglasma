// NuclearQsTable.cpp is part of the IP-Glasma solver.

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
                if (!fin.eof()) {
                    fin >> dummy;
                    fin >> T;
                    T_[iT] = atof(T.c_str());
                    fin >> Qs;
                    Qs2_[iT][iy] = atof(Qs.c_str());
                } else {
                    messager_ << "[NuclearQsTable::read]: End of file reached "
                                 "prematurely -- did the Q_s table file "
                                 "change? Exiting.";
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

double NuclearQsTable::qs2(double T, double y) const {
    double value, fracy, fracT, QsYdown, QsYup;
    int posy, check = 0;
    fracy = 0.;
    posy = static_cast<int>(floor(y / deltaY_ + 0.0000001));

    if (y > nY_ * deltaY_) {
        // qs2() is called from inside an omp parallel for loop
        // (Init::setColorChargeDensity), so use a fresh, stack-local
        // instance rather than sharing messager_, which is not thread-safe.
        PrettyOstream localMessager;
        localMessager << "[NuclearQsTable::qs2]: y=" << y
                      << " is above the tabulated range (max y="
                      << nY_ * deltaY_ << "). Exiting.";
        localMessager.flush("error");
        exit(1);
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

    for (int iT = 0; iT < nT_; iT++) {
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

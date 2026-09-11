#include "Lattice.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Glauber.h"
#include "Instrumentation.h"

Lattice::Lattice(Parameters *param, int length) {
    IPG_PROFILE_SCOPE("lattice.allocate");
    size_ = length * length;
    const double a = param->getL() / static_cast<double>(length);

    messager_ << "[Lattice::Lattice]: Allocating square lattice of size "
              << length << "x" << length << " with a=" << a << " fm ...";

    // Each vector is one contiguous field of fixed 3x3 matrices.  Preserve the
    // original Cell constructor semantics: all eight matrices start as I_3.
    const Matrix identity(1.0);
    U.assign(size_, identity);
    U2.assign(size_, identity);
    Ux.assign(size_, identity);
    Uy.assign(size_, identity);
    Ux1.assign(size_, identity);
    Uy1.assign(size_, identity);
    Ux2.assign(size_, identity);
    Uy2.assign(size_, identity);

    cellStorage.reserve(size_);
    cells.reserve(size_);
    for (int i = 0; i < size_; ++i) cellStorage.emplace_back();
    for (int i = 0; i < size_; ++i) cells.push_back(&cellStorage[i]);

    posmX.reserve(size_);
    pospX.reserve(size_);
    posmY.reserve(size_);
    pospY.reserve(size_);
    posmXpY.reserve(size_);
    pospXmY.reserve(size_);

    for (int i = 0; i < length; ++i) {
        const int im = std::max(0, i - 1);
        const int ip = std::min(length - 1, i + 1);
        for (int j = 0; j < length; ++j) {
            const int jm = std::max(0, j - 1);
            const int jp = std::min(length - 1, j + 1);

            pospX.push_back(ip * length + j);
            pospY.push_back(i * length + jp);
            posmX.push_back(im * length + j);
            posmY.push_back(i * length + jm);
            posmXpY.push_back(im * length + jp);
            pospXmY.push_back(ip * length + jm);
        }
    }

    messager_ << " done on rank " << param->getMPIRank() << ".";
    messager_.flush("info");
}

void Lattice::writeSU3Matrices(std::string fileprefix, Parameters *param) {
    // Logical aliases (see the field comment in Lattice.h): Ux2 <-> pi,
    // Uy2 <-> phi.
    const int N = param->getSize();

    std::stringstream strVOne_name;
    strVOne_name << fileprefix << "Phi-"
                 << param->getEventId()
                        + 2 * param->getSeed() * param->getMPISize()
                 << ".txt";

    std::stringstream strVTwo_name;
    strVTwo_name << fileprefix << "Pi-"
                 << param->getEventId()
                        + (1 + 2 * param->getSeed()) * param->getMPISize()
                 << ".txt";

    // Output in text
    std::ofstream foutU(strVOne_name.str().c_str(), std::ios::out);
    foutU.precision(15);

    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            int pos = ix * N + iy;
            foutU << ix << " " << iy << " " << Uy2[pos].MatrixToString()
                  << std::endl;
        }
        foutU << std::endl;
    }
    foutU.close();

    messager_ << "[Lattice::writeSU3Matrices]: wrote " << strVOne_name.str();
    messager_.flush("info");

    std::ofstream foutU2(strVTwo_name.str().c_str(), std::ios::out);
    foutU2.precision(15);
    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            int pos = ix * N + iy;
            foutU2 << ix << " " << iy << " " << Ux2[pos].MatrixToString()
                   << std::endl;
        }
        foutU2 << std::endl;
    }
    foutU2.close();

    messager_ << "[Lattice::writeSU3Matrices]: wrote " << strVTwo_name.str();
    messager_.flush("info");
}

void Lattice::writeWilsonLines(
    std::string fileprefix, Parameters *param, NucleusRole nucleus) {
    const int N = param->getSize();
    const double L = param->getL();
    const double a = L / static_cast<double>(N);  // lattice spacing in fm
    const bool isProjectile = (nucleus == NucleusRole::Projectile);
    // Preserves the historical iA=1 (projectile) / iA=2 (target) numbering
    // used in the output filename below.
    const int iA = isProjectile ? 1 : 2;

    std::stringstream strVOne_name;
    strVOne_name << fileprefix << "V-"
                 << param->getEventId()
                        + (iA + 2 * param->getSeed()) * param->getMPISize();
    if (param->getWriteWilsonLines() == 1) strVOne_name << ".txt";

    // Output in text
    if (param->getWriteWilsonLines() == 1) {
        std::ofstream foutU(strVOne_name.str().c_str(), std::ios::out);
        foutU.precision(15);

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                int pos = ix * N + iy;
                if (isProjectile) {
                    foutU << ix << " " << iy << " " << U[pos].MatrixToString()
                          << std::endl;
                } else {
                    foutU << ix << " " << iy << " " << U2[pos].MatrixToString()
                          << std::endl;
                }
            }
            foutU << std::endl;
        }
        foutU.close();

        messager_ << "[Lattice::writeWilsonLines]: wrote "
                  << strVOne_name.str();
        messager_.flush("info");
    } else if (param->getWriteWilsonLines() == 2) {
        std::ofstream Outfile1;
        Outfile1.open(
            strVOne_name.str().c_str(), std::ios::out | std::ios::binary);

        double temp = param->getRapidityA();
        if (!isProjectile) temp = param->getRapidityB();

        // print header ------------- //
        Outfile1.write((char *)&N, sizeof(int));
        Outfile1.write((char *)&Nc_, sizeof(int));
        Outfile1.write((char *)&L, sizeof(double));
        Outfile1.write((char *)&a, sizeof(double));
        Outfile1.write((char *)&temp, sizeof(double));

        double *val1 = new double[2];

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                for (int a1 = 0; a1 < 3; a1++) {
                    for (int b = 0; b < 3; b++) {
                        int indx = N * iy + ix;
                        int SU3indx = a1 * Nc_ + b;
                        if (isProjectile) {
                            val1[0] = U[indx].getRe(SU3indx);
                            val1[1] = U[indx].getIm(SU3indx);
                        } else {
                            val1[0] = U2[indx].getRe(SU3indx);
                            val1[1] = U2[indx].getIm(SU3indx);
                        }
                        Outfile1.write((char *)val1, 2 * sizeof(double));
                    }
                }
            }
        }

        if (Outfile1.good() == false) {
            messager_.error(
                "[Lattice::writeWilsonLines]: CRITICAL ERROR -- BINARY "
                "OUTPUT OF VECTOR CURRENTS FAILED");
            exit(1);
        }

        delete[] val1;

        Outfile1.close();
        messager_ << "[Lattice::writeWilsonLines]: wrote "
                  << strVOne_name.str();
        messager_.flush("info");
    } else {
        std::stringstream errorMsg;
        errorMsg << "[Lattice::writeWilsonLines]: Unknown option "
                    "param->getWriteWilsonLines()=="
                 << param->getWriteWilsonLines();
        messager_.error(errorMsg.str());
        exit(1);
    }
}

// constructor
BufferLattice::BufferLattice(int length) {
    size_ = length * length;

    const Matrix identity(1.0);
    buffer1.assign(size_, identity);
    buffer2.assign(size_, identity);
}

#include "Lattice.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Instrumentation.h"

namespace {

inline void requireSU3Lattice(int Nc) {
    if (Nc != 3) {
        std::cerr << "Error: lattice matrix storage is SU(3)-only; received Nc="
                  << Nc << ". Exiting." << std::endl;
        std::exit(1);
    }
}

}  // namespace

Lattice::Lattice(Parameters *param, int N, int length) {
    IPG_PROFILE_SCOPE("lattice.allocate");
    Nc = N;
    requireSU3Lattice(Nc);
    size = length * length;
    const double a = param->getL() / static_cast<double>(length);

    std::cout << "Allocating square lattice of size " << length << "x" << length
              << " with a=" << a << " fm ...";

    // Each vector is one contiguous field of fixed 3x3 matrices.  Preserve the
    // original Cell constructor semantics: all eight matrices start as I_3.
    const Matrix identity(3, 1.0);
    U.assign(size, identity);
    U2.assign(size, identity);
    Ux.assign(size, identity);
    Uy.assign(size, identity);
    Ux1.assign(size, identity);
    Uy1.assign(size, identity);
    Ux2.assign(size, identity);
    Uy2.assign(size, identity);

    cellStorage.reserve(size);
    cells.reserve(size);
    for (int i = 0; i < size; ++i) cellStorage.emplace_back(Nc);
    for (int i = 0; i < size; ++i) cells.push_back(&cellStorage[i]);

    posmX.reserve(size);
    pospX.reserve(size);
    posmY.reserve(size);
    pospY.reserve(size);
    posmXpY.reserve(size);
    pospXmY.reserve(size);

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

    std::cout << " done on rank " << param->getMPIRank() << "." << std::endl;
}

void Lattice::WriteSU3Matricies(std::string fileprefix, Parameters *param) {
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

    std::cout << "wrote " << strVOne_name.str() << std::endl;

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

    std::cout << "wrote " << strVTwo_name.str() << std::endl;
}

void Lattice::WriteWilsonLines(
    std::string fileprefix, Parameters *param, const int iA) {
    const int N = param->getSize();
    const double L = param->getL();
    const double a = L / static_cast<double>(N);  // lattice spacing in fm

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
                if (iA == 1) {
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

        std::cout << "wrote " << strVOne_name.str() << std::endl;
    } else if (param->getWriteWilsonLines() == 2) {
        std::ofstream Outfile1;
        Outfile1.open(
            strVOne_name.str().c_str(), std::ios::out | std::ios::binary);

        double temp = param->getRapidityA();
        if (iA == 2) temp = param->getRapidityB();

        // print header ------------- //
        Outfile1.write((char *)&N, sizeof(int));
        Outfile1.write((char *)&Nc, sizeof(int));
        Outfile1.write((char *)&L, sizeof(double));
        Outfile1.write((char *)&a, sizeof(double));
        Outfile1.write((char *)&temp, sizeof(double));

        double *val1 = new double[2];

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                for (int a1 = 0; a1 < 3; a1++) {
                    for (int b = 0; b < 3; b++) {
                        int indx = N * iy + ix;
                        int SU3indx = a1 * Nc + b;
                        if (iA == 1) {
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
            std::cerr << "#CRTICAL ERROR -- BINARY OUTPUT OF VECTOR "
                         "CURRENTS FAILED"
                      << std::endl;
            exit(1);
        }

        delete[] val1;

        Outfile1.close();
        std::cout << "wrote " << strVOne_name.str() << std::endl;
    } else {
        std::cerr << "# Unknwon option param->getWriteWilsonLines()=="
                  << param->getWriteWilsonLines() << std::endl;
        exit(1);
    }
}

// constructor
BufferLattice::BufferLattice(int N, int length) {
    Nc = N;
    size = length * length;

    const Matrix identity(3, 1.0);
    buffer1.assign(size, identity);
    buffer2.assign(size, identity);
}

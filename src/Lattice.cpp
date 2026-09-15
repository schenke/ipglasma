#include "Lattice.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
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

void Lattice::writeMatrixArrayText(
    const std::string &file_name, std::vector<Matrix> &field, int N) {
    std::ofstream fout(file_name.c_str(), std::ios::out);
    fout.precision(15);

    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            int pos = ix * N + iy;
            fout << ix << " " << iy << " " << field[pos].MatrixToString()
                 << std::endl;
        }
        fout << std::endl;
    }
    fout.close();

    messager_ << "[Lattice::writeSU3Matrices]: wrote " << file_name;
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

    writeMatrixArrayText(strVOne_name.str(), Uy2, N);
    writeMatrixArrayText(strVTwo_name.str(), Ux2, N);
}

void Lattice::writeWilsonLines(
    Parameters *param, NucleusRole nucleus, double x) {
    std::string wLineFile = generateWilsonLineDataFileName(param, x, nucleus);
    const int N = param->getSize();

    // Output in text
    if (param->getWriteWilsonLines() == 1) {
        std::ofstream foutU(wLineFile, std::ios::out);
        foutU.precision(15);

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                int pos = ix * N + iy;
                if (nucleus == NucleusRole::Projectile) {
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
                  << wLineFile;
        messager_.flush("info");
    } else if (param->getWriteWilsonLines() == 2) {
         
        const double L = param->getL();
        const double a = L / static_cast<double>(N);  // lattice spacing in fm
        
        std::ofstream Outfile1;
        Outfile1.open(wLineFile, std::ios::out | std::ios::binary);


        double temp = (nucleus == NucleusRole::Projectile) ? param->getRapidityA() : param->getRapidityB();

        // print header ------------- //
        Outfile1.write((char *)&N, sizeof(int));
        Outfile1.write((char *)&Nc_, sizeof(int));
        Outfile1.write((char *)&L, sizeof(double));
        Outfile1.write((char *)&a, sizeof(double));
        Outfile1.write((char *)&temp, sizeof(double));

        double val1[2];

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                for (int a1 = 0; a1 < 3; a1++) {
                    for (int b = 0; b < 3; b++) {
                        // Matches the text branch above and every other
                        // U/U2 indexing in the codebase; previously this
                        // read N * iy + ix, transposing the lattice in the
                        // binary output whenever ix != iy (invisible on a
                        // symmetric/identity lattice, which is why no test
                        // caught it -- see Init::readVFromFile's matching
                        // fix on the read side).
                        int indx = N * ix + iy;
                        int SU3indx = a1 * Nc_ + b;
                        if (nucleus == NucleusRole::Projectile) {
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
                "[Lattice::writeWilsonLines]: Failed to write the Wilson "
                "line binary output file.");
            exit(1);
        }

        Outfile1.close();
        messager_ << "[Lattice::writeWilsonLines]: wrote "
                  << wLineFile;
        messager_.flush("info");
    } else {
        std::stringstream errorMsg;
        errorMsg << "[Lattice::writeWilsonLines]: Unknown writeWilsonLines "
                    "value "
                 << param->getWriteWilsonLines() << ". Exiting.";
        messager_.error(errorMsg.str());
        exit(1);
    }
}

std::string Lattice::generateWilsonLineDataFileName(Parameters *param,
        const double x, NucleusRole nucleus)
{
    const bool isProjectile = (nucleus == NucleusRole::Projectile);
    const int iA = isProjectile ? 1 : 2;

    std::stringstream Vname;
    Vname << param->getWilsonLinePath() << "/WilsonLine";
    if (x >= 0) Vname << "_x_" << std::scientific
                 << std::setprecision(5) << x;
    Vname << "_" << param->getEventId()
                        + (iA + 2 * param->getSeed()) * param->getMPISize();

    if (param->getWriteWilsonLines() == 1) Vname << ".txt";

    return Vname.str();
}
// constructor
BufferLattice::BufferLattice(int length) {
    size_ = length * length;

    const Matrix identity(1.0);
    buffer1.assign(size_, identity);
    buffer2.assign(size_, identity);
}

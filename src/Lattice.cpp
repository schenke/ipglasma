// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include "Lattice.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "Glauber.h"
#include "Instrumentation.h"
#include "PhysConst.h"

using PhysConst::Nc;

Lattice::Lattice(Parameters *param, int length) {
    IPG_PROFILE_SCOPE("lattice.allocate");
    size_ = length * length;
    length_ = length;
    const double a = param->lattice.L / static_cast<double>(length);

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

            pospX.push_back(latticeIndex(ip, j, length));
            pospY.push_back(latticeIndex(i, jp, length));
            posmX.push_back(latticeIndex(im, j, length));
            posmY.push_back(latticeIndex(i, jm, length));
            posmXpY.push_back(latticeIndex(im, jp, length));
            pospXmY.push_back(latticeIndex(ip, jm, length));
        }
    }

    messager_ << " done on rank " << param->run.MPIRank << ".";
    messager_.flush("info");
}

void Lattice::writeMatrixArrayText(
    const std::string &file_name, std::vector<Matrix> &field, int N) {
    std::ofstream fout(file_name.c_str(), std::ios::out);
    fout.precision(15);

    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            int pos = positionFromXY(ix, iy);
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
    const int N = param->lattice.size;

    std::stringstream strVOne_name;
    strVOne_name << fileprefix << "Phi-"
                 << param->event.eventId
                        + 2 * param->random.fileSeed() * param->run.MPISize
                 << ".txt";

    std::stringstream strVTwo_name;
    strVTwo_name << fileprefix << "Pi-"
                 << param->event.eventId
                        + (1 + 2 * param->random.fileSeed())
                              * param->run.MPISize
                 << ".txt";

    writeMatrixArrayText(strVOne_name.str(), Uy2, N);
    writeMatrixArrayText(strVTwo_name.str(), Ux2, N);
}

// constructor
BufferLattice::BufferLattice(int length) {
    size_ = length * length;

    const Matrix identity(1.0);
    buffer1.assign(size_, identity);
    buffer2.assign(size_, identity);
}

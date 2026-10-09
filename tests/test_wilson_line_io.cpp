// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <set>
#include <string>
#include <vector>

#include "Glauber.h"  // for NucleusRole
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "WilsonLineIO.h"
#include "doctest.h"
#include "test_helpers.h"

TEST_CASE(
    "WilsonLineIO::isValidFormat accepts only 1 (text) and 2 "
    "(binary)") {
    CHECK(WilsonLineIO::isValidFormat(1) == true);
    CHECK(WilsonLineIO::isValidFormat(2) == true);
    CHECK(WilsonLineIO::isValidFormat(0) == false);
    CHECK(WilsonLineIO::isValidFormat(3) == false);
    CHECK(WilsonLineIO::isValidFormat(-1) == false);
}

TEST_CASE("WilsonLineIO::write (text format) writes a non-empty file") {
    const int length = 4;
    Parameters param;
    makeLatticeParam(param, length);
    param.wilsonLines.writeWilsonLines = 1;  // text
    param.colorCharge.useFluctuatingX =
        true;  // with this option, no x value in the generated filename
    Lattice lat(&param, length);

    WilsonLineIO().write(&lat, &param, NucleusRole::Projectile);

    // filename suffix = eventId + (iA + 2*seed)*MPISize; iA=1 for
    // Projectile, and eventId=seed=0, MPISize=1 here, so suffix = 1.
    const std::string path = "./WilsonLine_1.txt";
    std::ifstream in(path);
    REQUIRE(in.good());
    std::string firstLine;
    std::getline(in, firstLine);
    CHECK(!firstLine.empty());
    in.close();
    std::remove(path.c_str());
}

TEST_CASE(
    "WilsonLineIO::write (binary format) writes a matching header and "
    "data") {
    const int length = 4;
    Parameters param;
    makeLatticeParam(param, length);
    param.wilsonLines.writeWilsonLines = 2;  // binary
    param.colorCharge.useFluctuatingX =
        true;  // with this option, no x value in the generated filename
    Lattice lat(&param, length);

    // Give two off-diagonal sites (ix != iy, swapped between them) distinct
    // values, so a transposed site index (previously N*iy+ix instead of
    // N*ix+iy in the binary branch -- see Lattice.cpp) would read back the
    // wrong site's data here, instead of going unnoticed as it would on an
    // all-identity lattice.
    Matrix markerA(Matrix::noInit), markerB(Matrix::noInit);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            markerA.set(i, j, std::complex<double>(100 + 10 * i + j, 0.0));
            markerB.set(i, j, std::complex<double>(200 + 10 * i + j, 0.0));
        }
    }
    lat.U[lat.positionFromXY(1, 2)] = markerA;
    lat.U[lat.positionFromXY(2, 1)] = markerB;

    WilsonLineIO().write(&lat, &param, NucleusRole::Projectile);

    // filename suffix = eventId + (iA + 2*seed)*MPISize = 1, as above; no
    // extension is appended for the binary format.
    const std::string path = "./WilsonLine_1";
    std::ifstream in(path, std::ios::binary);
    REQUIRE(in.good());

    int n = 0, nc = 0;
    double L = 0.0, a = 0.0, x = 0.0;
    in.read(reinterpret_cast<char *>(&n), sizeof(int));
    in.read(reinterpret_cast<char *>(&nc), sizeof(int));
    in.read(reinterpret_cast<char *>(&L), sizeof(double));
    in.read(reinterpret_cast<char *>(&a), sizeof(double));
    in.read(reinterpret_cast<char *>(&x), sizeof(double));

    CHECK(n == length);
    CHECK(nc == 3);
    CHECK(L == doctest::Approx(param.lattice.L));
    CHECK(a == doctest::Approx(param.lattice.L / length));
    CHECK(x == -1.);  // no fixed x

    // The writer nests ix outer / iy inner (see Lattice.cpp), so the
    // site-th 9-(re, im)-pair block in the file must be lat.U[site], i.e.
    // lat.U[ix*length+iy] -- not lat.U[iy*length+ix]. Checking against
    // lat.U directly (rather than a hardcoded identity pattern) also
    // covers the two marker sites planted above.
    for (int site = 0; site < length * length; ++site) {
        for (int k = 0; k < 9; ++k) {
            double re = 0.0, im = 0.0;
            in.read(reinterpret_cast<char *>(&re), sizeof(double));
            in.read(reinterpret_cast<char *>(&im), sizeof(double));
            CHECK(re == doctest::Approx(lat.U[site].get(k).real()));
            CHECK(im == doctest::Approx(lat.U[site].get(k).imag()));
        }
    }
    CHECK(in.good());
    in.close();
    std::remove(path.c_str());
}

TEST_CASE(
    "WilsonLineIO::write (binary format) stores the x of the Wilson lines in "
    "the header") {
    const int length = 4;
    Parameters param;
    makeLatticeParam(param, length);
    param.wilsonLines.writeWilsonLines = 2;  // binary
    Lattice lat(&param, length);

    for (const auto &[nucleus, x] :
         {std::pair<NucleusRole, double> {NucleusRole::Projectile, 1e-3},
          {NucleusRole::Target, 2.5e-4}}) {
        const double xValue = x;
        CAPTURE(xValue);
        WilsonLineIO().write(&lat, &param, nucleus, x);
        const std::string path = WilsonLineIO::fileName(&param, x, nucleus);
        std::ifstream in(path, std::ios::binary);
        REQUIRE(in.good());
        // skip N, Nc, L and a
        in.seekg(2 * sizeof(int) + 2 * sizeof(double));
        double xInFile = 0.;
        in.read(reinterpret_cast<char *>(&xInFile), sizeof(double));
        CHECK(xInFile == x);
        in.close();
        std::remove(path.c_str());
    }
}

TEST_CASE(
    "WilsonLineIO::fileName: an explicit format overrides "
    "writeWilsonLines") {
    Parameters param;
    makeLatticeParam(param, 4);
    param.wilsonLines.wilsonLinePath = ".";

    // Reading binary (format 2) while writing text must not look for .txt
    param.wilsonLines.writeWilsonLines = 1;
    CHECK(
        WilsonLineIO::fileName(&param, -1., NucleusRole::Projectile, 2)
        == "./WilsonLine_1");
    // and vice versa
    param.wilsonLines.writeWilsonLines = 2;
    CHECK(
        WilsonLineIO::fileName(&param, -1., NucleusRole::Projectile, 1)
        == "./WilsonLine_1.txt");
    // the default (format < 0) follows writeWilsonLines
    CHECK(
        WilsonLineIO::fileName(&param, -1., NucleusRole::Projectile)
        == "./WilsonLine_1");
    param.wilsonLines.writeWilsonLines = 1;
    CHECK(
        WilsonLineIO::fileName(&param, -1., NucleusRole::Projectile)
        == "./WilsonLine_1.txt");
}

TEST_CASE(
    "WilsonLineIO::readText/readBinary put every site of a synthetic file "
    "in place, whatever the impact parameter") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    // the readers do not shift; Init::shiftFieldsWithImpactParameter() does
    param.event.b = 3.0;

    WilsonLineIO io;

    // Text format: N*N lines of "i j Re0 Im0 ... Re8 Im8" (dummies i,j are
    // discarded by the reader; loop order alone fixes position).
    const std::string textPath = "ipglasma_test_init_wilson.txt";
    {
        std::ofstream out(textPath);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                out << i << " " << j;
                for (int k = 0; k < 9; ++k) {
                    out << " " << (1000. * i + 100. * j + k) << " "
                        << -(1000. * i + 100. * j + k);
                }
                out << "\n";
            }
        }
    }

    Lattice lat(&param, N);
    io.readText(textPath, &param, lat.U);
    std::remove(textPath.c_str());

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            CHECK(
                lat.U[latticeIndex(i, j, N)].get(0).real()
                == doctest::Approx(1000. * i + 100. * j));
        }
    }

    // Binary format: header (N, Nc, L, a, dummy) then N*N*9 (re,im) pairs in
    // ix-outer/iy-inner block order.
    const std::string binPath = "ipglasma_test_init_wilson_bin";
    {
        std::ofstream out(binPath, std::ios::out | std::ios::binary);
        int Nc = 3;
        double L = param.lattice.L;
        double a = L / N;
        double dummy = 0.;
        out.write(reinterpret_cast<const char *>(&N), sizeof(int));
        out.write(reinterpret_cast<const char *>(&Nc), sizeof(int));
        out.write(reinterpret_cast<const char *>(&L), sizeof(double));
        out.write(reinterpret_cast<const char *>(&a), sizeof(double));
        out.write(reinterpret_cast<const char *>(&dummy), sizeof(double));
        for (int ix = 0; ix < N; ++ix) {
            for (int iy = 0; iy < N; ++iy) {
                for (int row = 0; row < 3; ++row) {
                    for (int col = 0; col < 3; ++col) {
                        double re = 1000. * ix + 100. * iy + 10. * row + col;
                        double im = -re;
                        out.write(
                            reinterpret_cast<const char *>(&re),
                            sizeof(double));
                        out.write(
                            reinterpret_cast<const char *>(&im),
                            sizeof(double));
                    }
                }
            }
        }
    }

    Lattice lat2(&param, N);
    io.readBinary(binPath, &param, lat2.U2);
    std::remove(binPath.c_str());

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            const Matrix &matrix = lat2.U2[latticeIndex(ix, iy, N)];
            CHECK(
                matrix.get(0).real()
                == doctest::Approx(1000. * ix + 100. * iy));
            CHECK(
                matrix.get(5).imag()
                == doctest::Approx(-(1000. * ix + 100. * iy + 12.)));
        }
    }
}

TEST_CASE(
    "WilsonLineIO::write -> WilsonLineIO::read{Text,Binary} "
    "round-trips every matrix element unchanged") {
    const int N = 4;
    for (int format : {1, 2}) {
        CAPTURE(format);
        Parameters param;
        makeInitTestParam(param, N);
        param.wilsonLines.wilsonLinePath = ".";
        param.wilsonLines.writeWilsonLines = format;

        Lattice lat(&param, N);
        for (int pos = 0; pos < N * N; ++pos) {
            lat.U[pos] = makeTestMatrix(pos);
        }
        WilsonLineIO().write(&lat, &param, NucleusRole::Projectile);
        const std::string path = WilsonLineIO::fileName(
            &param, -1., NucleusRole::Projectile, format);

        WilsonLineIO io;
        Lattice lat2(&param, N);
        if (format == 1) {
            io.readText(path, &param, lat2.U);
        } else {
            io.readBinary(path, &param, lat2.U);
        }
        std::remove(path.c_str());

        for (int pos = 0; pos < N * N; ++pos) {
            CAPTURE(pos);
            for (int k = 0; k < 9; ++k) {
                CHECK(
                    lat2.U[pos].get(k).real()
                    == doctest::Approx(lat.U[pos].get(k).real()));
                CHECK(
                    lat2.U[pos].get(k).imag()
                    == doctest::Approx(lat.U[pos].get(k).imag()));
            }
        }
    }
}

TEST_CASE(
    "WilsonLineIO: both formats put every column back in place when the "
    "lattice spacing is not a binary fraction") {
    // N = 10, L = 7: a = 0.7, where (a * i) / a truncates to i - 1 for
    // i = 3 and 6; the text reader used to truncate instead of rounding
    const int N = 10;
    for (int format : {1, 2}) {
        CAPTURE(format);
        Parameters param;
        makeInitTestParam(param, N);
        param.lattice.L = 7.;
        param.wilsonLines.wilsonLinePath = ".";
        param.wilsonLines.writeWilsonLines = format;

        Lattice lat(&param, N);
        for (int pos = 0; pos < N * N; ++pos) {
            lat.U2[pos] = makeTestMatrix(pos);  // the target's Wilson line
        }
        WilsonLineIO().write(&lat, &param, NucleusRole::Target);
        const std::string path =
            WilsonLineIO::fileName(&param, -1., NucleusRole::Target, format);

        WilsonLineIO io;
        Lattice lat2(&param, N);
        if (format == 1) {
            io.readText(path, &param, lat2.U2);
        } else {
            io.readBinary(path, &param, lat2.U2);
        }
        std::remove(path.c_str());

        for (int pos = 0; pos < N * N; ++pos) {
            CAPTURE(pos);
            CHECK(matricesClose(lat2.U2[pos], lat.U2[pos], 1e-12));
        }
    }
}

TEST_CASE(
    "WilsonLineIO::fileName: every nucleus and event of a run, and runs "
    "with another seed, get different numbers") {
    Parameters param;
    makeLatticeParam(param, 4);
    param.wilsonLines.writeWilsonLines = 2;
    // one event on one rank: 2 seed + 1 and 2 seed + 2, as before
    param.random.seed = 7;
    param.run.eventsPerRank = 1;
    CHECK(
        WilsonLineIO::fileName(&param, -1., NucleusRole::Projectile)
        == "./WilsonLine_15");
    CHECK(
        WilsonLineIO::fileName(&param, -1., NucleusRole::Target)
        == "./WilsonLine_16");

    // 3 events on each of 2 ranks, with seeds 7 and 8
    param.run.MPISize = 2;
    param.run.eventsPerRank = 3;
    std::set<std::string> names;
    for (const unsigned long long seed : {7ULL, 8ULL}) {
        param.random.seed = seed;
        for (int eventId = 0; eventId < 6; eventId++) {
            param.event.eventId = eventId;
            for (const NucleusRole nucleus :
                 {NucleusRole::Projectile, NucleusRole::Target}) {
                names.insert(WilsonLineIO::fileName(&param, 0.01, nucleus));
            }
        }
    }
    CHECK(names.size() == 2 * 6 * 2);
}

TEST_CASE(
    "WilsonLineIO::writeGeometry -> readGeometry restores the nucleons, the "
    "color-charge and thickness maps and QsMuRatio exactly") {
    const int N = 4;
    Parameters param;
    makeInitTestParam(param, N);
    param.wilsonLines.wilsonLinePath = ".";
    param.run.eventsPerRank = 1;
    param.collision.projectile = "Au";
    param.collision.target = "Pb";
    param.coupling.g = 1.;
    param.colorCharge.QsMuRatio = 0.643;

    Lattice lat(&param, N);
    for (int pos = 0; pos < N * N; ++pos) {
        lat.cells[pos]->setg2mu2A(0.1 * pos + 1. / 3.);
        lat.cells[pos]->setg2mu2B(0.2 * pos + 1. / 7.);
        lat.cells[pos]->setTpA(0.01 * pos);
        lat.cells[pos]->setTpB(0.02 * pos + 1e-9);
    }
    std::vector<ReturnValue> nucleons(3);
    for (int i = 0; i < 3; ++i) {
        nucleons[i].x = 1.1 * i - 0.3;
        nucleons[i].y = -0.7 * i + 1. / 3.;
        nucleons[i].z = 0.5 * i;
        nucleons[i].proton = (i != 1);
        nucleons[i].collided = 1;
    }

    WilsonLineIO io;
    // the same number as the Wilson-line files
    CHECK(
        WilsonLineIO::geometryFileName(&param, NucleusRole::Target)
        == "./WilsonLineGeometry_2");
    io.writeGeometry(&lat, &param, NucleusRole::Projectile, nucleons);
    io.writeGeometry(&lat, &param, NucleusRole::Target, nucleons);

    param.colorCharge.QsMuRatio = 0.5;  // the file's value wins
    Lattice lat2(&param, N);
    const NucleusGeometry projectile =
        io.readGeometry(&lat2, &param, NucleusRole::Projectile);
    const NucleusGeometry target =
        io.readGeometry(&lat2, &param, NucleusRole::Target);
    std::remove(WilsonLineIO::geometryFileName(&param, NucleusRole::Projectile)
                    .c_str());
    std::remove(
        WilsonLineIO::geometryFileName(&param, NucleusRole::Target).c_str());

    CHECK(projectile.QsMuRatio == 0.643);
    CHECK(target.QsMuRatio == 0.643);
    for (const NucleusGeometry *geometry : {&projectile, &target}) {
        REQUIRE(geometry->nucleons.size() == 3);
        for (int i = 0; i < 3; ++i) {
            CAPTURE(i);
            CHECK(geometry->nucleons[i].x == nucleons[i].x);
            CHECK(geometry->nucleons[i].y == nucleons[i].y);
            CHECK(geometry->nucleons[i].z == nucleons[i].z);
            CHECK(geometry->nucleons[i].proton == nucleons[i].proton);
            CHECK(geometry->nucleons[i].collided == 0);
        }
    }
    for (int pos = 0; pos < N * N; ++pos) {
        CAPTURE(pos);
        CHECK(lat2.cells[pos]->getg2mu2A() == lat.cells[pos]->getg2mu2A());
        CHECK(lat2.cells[pos]->getg2mu2B() == lat.cells[pos]->getg2mu2B());
        CHECK(lat2.cells[pos]->getTpA() == lat.cells[pos]->getTpA());
        CHECK(lat2.cells[pos]->getTpB() == lat.cells[pos]->getTpB());
    }
}

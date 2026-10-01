#include <complex>
#include <cstdio>
#include <fstream>
#include <string>

#include "Glauber.h"  // for NucleusRole
#include "Lattice.h"
#include "Parameters.h"
#include "doctest.h"

TEST_CASE(
    "Lattice::IsValidWilsonLineDataFormat accepts only 1 (text) and 2 "
    "(binary)") {
    CHECK(Lattice::IsValidWilsonLineDataFormat(1) == true);
    CHECK(Lattice::IsValidWilsonLineDataFormat(2) == true);
    CHECK(Lattice::IsValidWilsonLineDataFormat(0) == false);
    CHECK(Lattice::IsValidWilsonLineDataFormat(3) == false);
    CHECK(Lattice::IsValidWilsonLineDataFormat(-1) == false);
}

namespace {
// Takes an out-parameter rather than returning by value: Parameters holds
// a PrettyOstream member, which holds a non-copyable/non-movable
// std::ostringstream -- see the identical note in test_parameters.cpp.
void makeLatticeParam(Parameters &param, int size) {
    param.lattice.L = 10.0;
    param.run.MPIRank = 0;
    param.lattice.size = size;
    param.event.eventId = 0;
    param.random.seed = 0;
    param.run.MPISize = 1;
    param.colorCharge.RapidityA = 0.0;
    param.colorCharge.RapidityB = 0.0;
    param.wilsonLines.wilsonLinePath = ".";
}
}  // namespace

TEST_CASE(
    "Lattice: constructor allocates fields as identity and sizes them "
    "correctly") {
    const int length = 4;
    Parameters param;
    makeLatticeParam(param, length);
    Lattice lat(&param, length);

    CHECK(lat.getSize() == length * length);
    CHECK(lat.U.size() == static_cast<std::size_t>(length * length));
    CHECK(lat.cells.size() == static_cast<std::size_t>(length * length));

    const Matrix identity(1.0);
    for (int pos = 0; pos < length * length; ++pos) {
        for (int k = 0; k < 9; ++k) {
            CHECK(std::abs(lat.U[pos].get(k) - identity.get(k)) < 1e-14);
            CHECK(std::abs(lat.Uy2[pos].get(k) - identity.get(k)) < 1e-14);
        }
    }
}

TEST_CASE(
    "Lattice: neighbor index arrays wrap at the boundary (clamped, not "
    "periodic)") {
    const int length = 4;
    Parameters param;
    makeLatticeParam(param, length);
    Lattice lat(&param, length);

    // Site (0, 0): both minus-neighbors clamp back to the site itself.
    CHECK(lat.posmX[0] == 0);
    CHECK(lat.posmY[0] == 0);
    // Site (0, 0): plus-neighbors move one step in x or y.
    CHECK(lat.pospX[0] == length);  // (1, 0)
    CHECK(lat.pospY[0] == 1);       // (0, 1)

    // Site (length-1, length-1): both plus-neighbors clamp back to itself.
    const int last = length * length - 1;
    CHECK(lat.pospX[last] == last);
    CHECK(lat.pospY[last] == last);
}

TEST_CASE("Lattice::writeWilsonLines (text format) writes a non-empty file") {
    const int length = 4;
    Parameters param;
    makeLatticeParam(param, length);
    param.wilsonLines.writeWilsonLines = 1;  // text
    param.colorCharge.useFluctuatingx =
        true;  // with this option, no x value in the generated filename
    Lattice lat(&param, length);

    lat.writeWilsonLines(&param, NucleusRole::Projectile);

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
    "Lattice::writeWilsonLines (binary format) writes a matching header and "
    "data") {
    const int length = 4;
    Parameters param;
    makeLatticeParam(param, length);
    param.wilsonLines.writeWilsonLines = 2;  // binary
    param.colorCharge.useFluctuatingx =
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

    lat.writeWilsonLines(&param, NucleusRole::Projectile);

    // filename suffix = eventId + (iA + 2*seed)*MPISize = 1, as above; no
    // extension is appended for the binary format.
    const std::string path = "./WilsonLine_1";
    std::ifstream in(path, std::ios::binary);
    REQUIRE(in.good());

    int n = 0, nc = 0;
    double L = 0.0, a = 0.0, rapidity = 0.0;
    in.read(reinterpret_cast<char *>(&n), sizeof(int));
    in.read(reinterpret_cast<char *>(&nc), sizeof(int));
    in.read(reinterpret_cast<char *>(&L), sizeof(double));
    in.read(reinterpret_cast<char *>(&a), sizeof(double));
    in.read(reinterpret_cast<char *>(&rapidity), sizeof(double));

    CHECK(n == length);
    CHECK(nc == 3);
    CHECK(L == doctest::Approx(param.lattice.L));
    CHECK(a == doctest::Approx(param.lattice.L / length));
    CHECK(rapidity == doctest::Approx(param.colorCharge.RapidityA));

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
    "Lattice::generateWilsonLineDataFileName: an explicit format overrides "
    "writeWilsonLines") {
    Parameters param;
    makeLatticeParam(param, 4);
    param.wilsonLines.wilsonLinePath = ".";

    // Reading binary (format 2) while writing text must not look for .txt
    param.wilsonLines.writeWilsonLines = 1;
    CHECK(
        Lattice::generateWilsonLineDataFileName(
            &param, -1., NucleusRole::Projectile, 2)
        == "./WilsonLine_1");
    // and vice versa
    param.wilsonLines.writeWilsonLines = 2;
    CHECK(
        Lattice::generateWilsonLineDataFileName(
            &param, -1., NucleusRole::Projectile, 1)
        == "./WilsonLine_1.txt");
    // the default (format < 0) follows writeWilsonLines
    CHECK(
        Lattice::generateWilsonLineDataFileName(
            &param, -1., NucleusRole::Projectile)
        == "./WilsonLine_1");
    param.wilsonLines.writeWilsonLines = 1;
    CHECK(
        Lattice::generateWilsonLineDataFileName(
            &param, -1., NucleusRole::Projectile)
        == "./WilsonLine_1.txt");
}

TEST_CASE("Lattice::writeSU3Matrices writes non-empty Phi/Pi files") {
    const int length = 4;
    Parameters param;
    makeLatticeParam(param, length);
    Lattice lat(&param, length);

    const std::string prefix = "ipglasma_test_lattice_su3_";
    lat.writeSU3Matrices(prefix, &param);

    // Phi suffix = eventId + 2*seed*MPISize = 0; Pi suffix =
    // eventId + (1 + 2*seed)*MPISize = 1, with eventId=seed=0, MPISize=1.
    for (const char *name : {"Phi-0.txt", "Pi-1.txt"}) {
        const std::string path = prefix + name;
        std::ifstream in(path);
        REQUIRE(in.good());
        std::string firstLine;
        std::getline(in, firstLine);
        CHECK(!firstLine.empty());
        in.close();
        std::remove(path.c_str());
    }
}

TEST_CASE(
    "BufferLattice: allocates buffer1/buffer2 as identity, sized correctly") {
    const int length = 4;
    BufferLattice buf(length);
    CHECK(buf.buffer1.size() == static_cast<std::size_t>(length * length));
    CHECK(buf.buffer2.size() == static_cast<std::size_t>(length * length));

    const Matrix identity(1.0);
    for (int k = 0; k < 9; ++k) {
        CHECK(std::abs(buf.buffer1[0].get(k) - identity.get(k)) < 1e-14);
    }
}

TEST_CASE(
    "Lattice::positionFromXY/xFromPosition/yFromPosition are consistent "
    "with each other, with latticeIndex() and with the neighbor tables") {
    const int length = 5;
    Parameters param;
    makeLatticeParam(param, length);
    Lattice lat(&param, length);
    CHECK(lat.getLength() == length);

    for (int ix = 0; ix < length; ++ix) {
        for (int iy = 0; iy < length; ++iy) {
            const int pos = lat.positionFromXY(ix, iy);
            CHECK(pos == latticeIndex(ix, iy, length));
            CHECK(lat.xFromPosition(pos) == ix);
            CHECK(lat.yFromPosition(pos) == iy);
            CHECK(latticeX(pos, length) == ix);
            CHECK(latticeY(pos, length) == iy);
            // interior sites: the neighbor tables step one site in x / y
            if (ix > 0 && ix < length - 1 && iy > 0 && iy < length - 1) {
                CHECK(lat.pospX[pos] == lat.positionFromXY(ix + 1, iy));
                CHECK(lat.posmX[pos] == lat.positionFromXY(ix - 1, iy));
                CHECK(lat.pospY[pos] == lat.positionFromXY(ix, iy + 1));
                CHECK(lat.posmY[pos] == lat.positionFromXY(ix, iy - 1));
            }
        }
    }
}

// Checks the output-file layouts described in OUTPUT.md, so the document and
// the writers cannot drift apart unnoticed.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Cell.h"
#include "Lattice.h"
#include "MyEigen.h"
#include "Parameters.h"
#include "PhysConst.h"
#include "WilsonLineIO.h"
#include "doctest.h"
#include "test_helpers.h"

using PhysConst::hbarc;

namespace {
const int N = 8;
const int it = 5;  // tau = it * dtau * a = 0.5 fm
const int gridPoints = 6;
const double e = 1.0;  // energy density in the fluid rest frame [fm^-4]
const double p = e / 3.;
const double v = 0.2;  // flow velocity along x

// An 8x8 lattice with a = 1 fm, filled with an ideal fluid moving along x,
// and a 6x6 output grid with spacing 1 fm inside it.
void makeOutputTestParam(Parameters &param) {
    makeLatticeParam(param, N);
    param.lattice.L = N;
    param.run.dtau = 0.1;
    param.coupling.g = 1.;
    param.coupling.runningCoupling = false;  // no running-coupling factor
    param.output.sizeOutput = gridPoints;
    param.output.LOutput = gridPoints;
    param.output.etaSizeOutput = 1;
    param.output.dEtaOutput = 0.1;
    param.output.computeGluonMultiplicity = false;
}

double lorentzGamma() { return 1. / std::sqrt(1. - v * v); }

void fillIdealFluid(Lattice &lat, const Parameters &param) {
    const double a = param.lattice.L / N;
    const double tau = it * param.run.dtau * a;
    const double g2 = lorentzGamma() * lorentzGamma();
    for (int pos = 0; pos < N * N; ++pos) {
        Cell *cell = lat.cells[pos];
        cell->setTtautau((e + p) * g2 - p);
        cell->setTtaux((e + p) * g2 * v);
        cell->setTxx((e + p) * g2 * v * v + p);
        cell->setTyy(p);
        cell->setTetaeta(p / (tau * tau));
        cell->setg2mu2A(1.);
        cell->setg2mu2B(1.);
    }
}

std::vector<std::string> tokens(const std::string &line) {
    std::istringstream in(line);
    std::vector<std::string> result;
    for (std::string token; in >> token;) result.push_back(token);
    return result;
}

// Data lines (non-empty, not starting with '#') and their column counts.
struct TextTable {
    std::string header;
    std::vector<std::vector<std::string>> rows;
    int blankLines = 0;
};

TextTable readTextTable(const std::string &fileName) {
    std::ifstream in(fileName);
    REQUIRE(in.good());
    TextTable table;
    for (std::string line; std::getline(in, line);) {
        if (line.empty()) {
            table.blankLines++;
        } else if (line[0] == '#') {
            table.header = line;
        } else {
            table.rows.push_back(tokens(line));
        }
    }
    return table;
}

std::uint64_t readLittleEndian(std::ifstream &in, int bytes) {
    std::uint64_t value = 0;
    for (int i = 0; i < bytes; ++i) {
        const int byte = in.get();
        value |= static_cast<std::uint64_t>(byte & 0xff) << (8 * i);
    }
    return value;
}

// Every T^{mu nu} file component, in the order documented in OUTPUT.md.
const char *const tmunuComponents =
    "\"components\":[\"T00\",\"Txx\",\"Tyy\",\"tau2_Tetaeta\",\"neg_T0x\","
    "\"neg_T0y\",\"neg_tau_T0eta\",\"neg_Txy\",\"neg_tau_Tyeta\","
    "\"neg_tau_Txeta\"]";
}  // namespace

TEST_CASE(
    "Output format: binary T^{mu nu} (.ipgt) has the documented header and "
    "[y][x][10] float32 data in GeV/fm^3") {
    unsetenv("IPGLASMA_BINARY_TMUNU");
    Parameters param;
    makeOutputTestParam(param);
    param.output.writeOutputs = 4;
    param.output.writeTmunuBinary = true;
    Lattice lat(&param, N);
    fillIdealFluid(lat, param);
    MyEigen().writeTmunu4D(&lat, &param, it);

    const std::string fileName = "Tmunu-t0.5-0.ipgt";
    std::ifstream in(fileName, std::ios::binary);
    REQUIRE(in.good());
    char magic[8];
    in.read(magic, 8);
    CHECK(std::string(magic, 8) == "IPGTMU01");
    const std::uint64_t metadataLength = readLittleEndian(in, 4);
    std::string metadata(metadataLength, ' ');
    in.read(&metadata[0], static_cast<std::streamsize>(metadataLength));
    CHECK(metadata.find("\"format\":\"ipglasma-tmunu\"") != std::string::npos);
    CHECK(metadata.find("\"shape\":[6,6,10]") != std::string::npos);
    CHECK(
        metadata.find("\"axis_order\":[\"y\",\"x\",\"component\"]")
        != std::string::npos);
    CHECK(metadata.find(tmunuComponents) != std::string::npos);
    CHECK(metadata.find("\"tau_fm\":0.5") != std::string::npos);

    std::vector<float> data(gridPoints * gridPoints * 10);
    in.read(
        reinterpret_cast<char *>(data.data()),
        static_cast<std::streamsize>(data.size() * sizeof(float)));
    CHECK(in.gcount() == static_cast<std::streamsize>(data.size() * 4));
    CHECK(in.peek() == EOF);  // nothing after the data
    in.close();
    std::remove(fileName.c_str());

    const double g2 = lorentzGamma() * lorentzGamma();
    const float *point = &data[(2 * gridPoints + 3) * 10];  // y = 2, x = 3
    CHECK(point[0] == doctest::Approx(((e + p) * g2 - p) * hbarc));
    CHECK(point[1] == doctest::Approx(((e + p) * g2 * v * v + p) * hbarc));
    CHECK(point[2] == doctest::Approx(p * hbarc));
    CHECK(point[3] == doctest::Approx(p * hbarc));  // tau^2 T^{eta eta}
    CHECK(point[4] == doctest::Approx(-(e + p) * g2 * v * hbarc));
}

TEST_CASE(
    "Output format: text T^{mu nu} (.dat) has the documented header and "
    "12-column rows, y outer") {
    unsetenv("IPGLASMA_BINARY_TMUNU");
    Parameters param;
    makeOutputTestParam(param);
    param.output.writeOutputs = 4;
    param.output.writeTmunuBinary = false;
    Lattice lat(&param, N);
    fillIdealFluid(lat, param);
    MyEigen().writeTmunu4D(&lat, &param, it);

    const std::string fileName = "Tmunu-t0.5-0.dat";
    const TextTable table = readTextTable(fileName);
    std::remove(fileName.c_str());
    CHECK(
        tokens(table.header)
        == std::vector<std::string> {
            "#", "dummy", "1", "etamax=", "1", "xmax=", "6", "ymax=", "6",
            "deta=", "0.1", "dx=", "1", "dy=", "1"});
    REQUIRE(table.rows.size() == gridPoints * gridPoints);
    CHECK(table.blankLines == gridPoints);
    for (const auto &row : table.rows) CHECK(row.size() == 12);
    // y outer, x inner: the second row is (ix, iy) = (1, 0)
    CHECK(table.rows[1][0] == "1");
    CHECK(table.rows[1][1] == "0");
}

TEST_CASE(
    "Output format: the hydro file has the documented header and 18 columns "
    "(eps in GeV/fm^3, u^mu)") {
    Parameters param;
    makeOutputTestParam(param);
    param.output.writeOutputs = 1;
    param.output.writeEpsilonUHydro = true;
    Lattice lat(&param, N);
    fillIdealFluid(lat, param);
    MyEigen().flowVelocity4D(&lat, &param, it, /*finalFlag=*/true);

    const std::string fileName = "epsilon-u-Hydro-TauHydro-0.dat";
    const TextTable table = readTextTable(fileName);
    std::remove(fileName.c_str());
    const std::vector<std::string> header = tokens(table.header);
    REQUIRE(header.size() == 17);
    CHECK(header[15] == "tau=");
    CHECK(header[16] == "0.5");
    REQUIRE(table.rows.size() == gridPoints * gridPoints);
    CHECK(table.blankLines == 1);  // after the single eta slice
    for (const auto &row : table.rows) {
        REQUIRE(row.size() == 18);
        CHECK(std::stod(row[3]) == doctest::Approx(e * hbarc));       // eps
        CHECK(std::stod(row[4]) == doctest::Approx(lorentzGamma()));  // u^tau
        CHECK(std::stod(row[5]) == doctest::Approx(lorentzGamma() * v));  // u^x
        for (int k = 8; k < 18; ++k) {  // an ideal fluid has pi = 0
            CHECK(std::abs(std::stod(row[k])) < 1e-8);
        }
    }
    // x outer, y inner: the second row has the same x and the next y
    CHECK(table.rows[1][1] == table.rows[0][1]);
    CHECK(std::stod(table.rows[1][2]) == doctest::Approx(-2.));
}

TEST_CASE(
    "Output format: the Jazma file has 18 columns with u = (1, 0, 0, 0) and "
    "pi = 0") {
    Parameters param;
    makeOutputTestParam(param);
    param.output.writeOutputs = 2;
    param.output.writeEpsilonUHydro = true;
    Lattice lat(&param, N);
    fillIdealFluid(lat, param);
    MyEigen().flowVelocity4D(&lat, &param, it, /*finalFlag=*/true);

    const std::string fileName = "Jazma-Hydro-t0.5-0.dat";
    const TextTable table = readTextTable(fileName);
    std::remove(fileName.c_str());
    CHECK(tokens(table.header).size() == 15);  // no tau= entry
    REQUIRE(table.rows.size() == gridPoints * gridPoints);
    for (const auto &row : table.rows) {
        REQUIRE(row.size() == 18);
        CHECK(std::stod(row[3]) > 0.);
        CHECK(std::stod(row[4]) == 1.);
        for (int k = 5; k < 18; ++k) CHECK(std::stod(row[k]) == 0.);
    }
}

TEST_CASE(
    "Output format: initialWilsonLines.ipgw has the documented header and "
    "[2][2][N][N][3][3] float32 data") {
    Parameters param;
    makeOutputTestParam(param);
    Lattice lat(&param, N);
    for (int pos = 0; pos < N * N; ++pos) {
        lat.U[pos] = makeTestMatrix(pos);
        lat.U2[pos] = makeTestMatrix(pos + 1000);
    }
    WilsonLineIO().writeTrainingData(&lat, &param);

    const std::string fileName = "initialWilsonLines0.ipgw";
    std::ifstream in(fileName, std::ios::binary);
    REQUIRE(in.good());
    char magic[8];
    in.read(magic, 8);
    CHECK(std::memcmp(magic, "IPGWIL1\0", 8) == 0);
    const std::uint64_t metadataLength = readLittleEndian(in, 8);
    std::string metadata(metadataLength, ' ');
    in.read(&metadata[0], static_cast<std::streamsize>(metadataLength));
    CHECK(
        metadata.find("\"format\":\"ipglasma-initial-wilson-lines\"")
        != std::string::npos);
    CHECK(metadata.find("\"shape\":[2,2,8,8,3,3]") != std::string::npos);

    const std::size_t elements = 2 * 2 * N * N * 9;
    std::vector<float> data(elements);
    in.read(
        reinterpret_cast<char *>(data.data()),
        static_cast<std::streamsize>(elements * sizeof(float)));
    CHECK(in.gcount() == static_cast<std::streamsize>(elements * 4));
    CHECK(in.peek() == EOF);
    in.close();
    std::remove(fileName.c_str());

    // [beam][complex part][x][y][row][col]
    auto at = [&](int beam, int part, int x, int y, int row, int col) {
        return data
            [((((beam * 2 + part) * N + x) * N + y) * 3 + row) * 3 + col];
    };
    const int x = 2, y = 5, row = 1, col = 2;
    const int pos = lat.positionFromXY(x, y);
    CHECK(
        at(0, 0, x, y, row, col)
        == doctest::Approx(lat.U[pos].get(row * 3 + col).real()));
    CHECK(
        at(0, 1, x, y, row, col)
        == doctest::Approx(lat.U[pos].get(row * 3 + col).imag()));
    CHECK(
        at(1, 0, x, y, row, col)
        == doctest::Approx(lat.U2[pos].get(row * 3 + col).real()));
}

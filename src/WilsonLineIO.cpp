// WilsonLineIO.cpp is part of the IP-Glasma solver.

#include "WilsonLineIO.h"

#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Instrumentation.h"
#include "PhysConst.h"

using PhysConst::Nc;
using std::complex;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;

std::string WilsonLineIO::fileName(
    Parameters *param, const double x, NucleusRole nucleus, int format) {
    const bool isProjectile = (nucleus == NucleusRole::Projectile);
    const int iA = isProjectile ? 1 : 2;

    std::stringstream Vname;
    Vname << param->wilsonLines.wilsonLinePath << "/WilsonLine";
    if (x >= 0) Vname << "_x_" << std::scientific << std::setprecision(5) << x;
    Vname << "_"
          << param->event.eventId
                 + (iA + 2 * param->random.seed) * param->run.MPISize;

    const int fileFormat =
        (format < 0) ? param->wilsonLines.writeWilsonLines : format;
    if (fileFormat == 1) Vname << ".txt";

    return Vname.str();
}

void WilsonLineIO::write(
    Lattice *lat, Parameters *param, NucleusRole nucleus, double x) {
    std::string wLineFile = fileName(param, x, nucleus);
    const int N = param->lattice.size;

    // Output in text
    if (param->wilsonLines.writeWilsonLines == 1) {
        std::ofstream foutU(wLineFile, std::ios::out);
        foutU.precision(15);

        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                int pos = lat->positionFromXY(ix, iy);
                if (nucleus == NucleusRole::Projectile) {
                    foutU << ix << " " << iy << " "
                          << lat->U[pos].MatrixToString() << std::endl;
                } else {
                    foutU << ix << " " << iy << " "
                          << lat->U2[pos].MatrixToString() << std::endl;
                }
            }
            foutU << std::endl;
        }
        foutU.close();

        messager_ << "[WilsonLineIO::write]: wrote " << wLineFile;
        messager_.flush("info");
    } else if (param->wilsonLines.writeWilsonLines == 2) {
        const double L = param->lattice.L;
        const double a = L / static_cast<double>(N);  // lattice spacing in fm

        std::ofstream Outfile1;
        Outfile1.open(wLineFile, std::ios::out | std::ios::binary);

        double temp = (nucleus == NucleusRole::Projectile)
                          ? param->colorCharge.rapidityA
                          : param->colorCharge.rapidityB;

        // print header ------------- //
        Outfile1.write((char *)&N, sizeof(int));
        Outfile1.write((char *)&Nc, sizeof(int));
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
                        // caught it -- see WilsonLineIO::read's matching
                        // fix on the read side).
                        int indx = lat->positionFromXY(ix, iy);
                        int SU3indx = a1 * Nc + b;
                        if (nucleus == NucleusRole::Projectile) {
                            val1[0] = lat->U[indx].getRe(SU3indx);
                            val1[1] = lat->U[indx].getIm(SU3indx);
                        } else {
                            val1[0] = lat->U2[indx].getRe(SU3indx);
                            val1[1] = lat->U2[indx].getIm(SU3indx);
                        }
                        Outfile1.write((char *)val1, 2 * sizeof(double));
                    }
                }
            }
        }

        if (Outfile1.good() == false) {
            messager_.error(
                "[WilsonLineIO::write]: Failed to write the Wilson "
                "line binary output file.");
            exit(1);
        }

        Outfile1.close();
        messager_ << "[WilsonLineIO::write]: wrote " << wLineFile;
        messager_.flush("info");
    } else {
        std::stringstream errorMsg;
        errorMsg << "[WilsonLineIO::write]: Unknown writeWilsonLines "
                    "value "
                 << param->wilsonLines.writeWilsonLines << ". Exiting.";
        messager_.error(errorMsg.str());
        exit(1);
    }
}

void WilsonLineIO::read(Lattice *lat, Parameters *param, int format, double x) {
    IPG_PROFILE_SCOPE("initialization.read_wilson_lines");
    if (!isValidFormat(format)) {
        messager_ << "[WilsonLineIO::read]: Unknown format " << format
                  << " when reading the initial Wilson lines, supported "
                     "formats: 1,2";
        messager_.flush("error");
        exit(1);
    }

    string VOne_name = fileName(param, x, NucleusRole::Projectile, format);
    string VTwo_name = fileName(param, x, NucleusRole::Target, format);

    messager_ << "[WilsonLineIO::read]: Reading Wilson lines from files "
              << VOne_name << " and " << VTwo_name;
    messager_.flush("info");

    if (format == 1) {
        readText(VOne_name, param, NucleusRole::Projectile, lat->U);
        readText(VTwo_name, param, NucleusRole::Target, lat->U2);
    } else if (format == 2) {
        readBinary(VOne_name, param, NucleusRole::Projectile, lat->U);
        readBinary(VTwo_name, param, NucleusRole::Target, lat->U2);
    }

    messager_ << "[WilsonLineIO::read]: Wilson lines V_A and V_B set on rank "
              << param->run.MPIRank << ". ";
    messager_.flush("info");
}

void WilsonLineIO::readText(
    const std::string &fileName, Parameters *param, NucleusRole role,
    std::vector<Matrix> &U) {
    const bool isProjectile = (role == NucleusRole::Projectile);
    int N = param->lattice.size;

    double L = param->lattice.L;
    double a = L / static_cast<double>(N);

    Matrix temp(1.);

    double Re[9], Im[9];
    double dummy;

    ifstream fin(fileName.c_str(), std::ios::in);

    if (!fin) {
        messager_ << "[WilsonLineIO::read]: File " << fileName
                  << " not found. Exiting.";
        messager_.flush("error");
        exit(1);
    }

    messager_ << "[WilsonLineIO::read]: Reading Wilson line from file "
              << fileName << " ...";
    messager_.flush("info");

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fin >> dummy >> dummy >> Re[0] >> Im[0] >> Re[1] >> Im[1] >> Re[2]
                >> Im[2] >> Re[3] >> Im[3] >> Re[4] >> Im[4] >> Re[5] >> Im[5]
                >> Re[6] >> Im[6] >> Re[7] >> Im[7] >> Re[8] >> Im[8];

            temp.set(0, 0, complex<double>(Re[0], Im[0]));
            temp.set(0, 1, complex<double>(Re[1], Im[1]));
            temp.set(0, 2, complex<double>(Re[2], Im[2]));
            temp.set(1, 0, complex<double>(Re[3], Im[3]));
            temp.set(1, 1, complex<double>(Re[4], Im[4]));
            temp.set(1, 2, complex<double>(Re[5], Im[5]));
            temp.set(2, 0, complex<double>(Re[6], Im[6]));
            temp.set(2, 1, complex<double>(Re[7], Im[7]));
            temp.set(2, 2, complex<double>(Re[8], Im[8]));

            double bb = param->event.b;
            a = L / static_cast<double>(N);

            double xtemp = isProjectile ? (a * i - bb / 2.) : (a * i + bb / 2.);
            int ix = xtemp / a;

            if (isProjectile) {
                if (ix < 0) continue;
            } else {
                if (ix >= N) continue;
            }

            int pos = latticeIndex(ix, j, N);
            U[pos] = (temp);
        }
    }

    fin.close();
}

void WilsonLineIO::readBinary(
    const std::string &fileName, Parameters *param, NucleusRole role,
    std::vector<Matrix> &U) {
    const bool isProjectile = (role == NucleusRole::Projectile);
    std::ifstream InStream;
    InStream.precision(15);
    InStream.open(fileName.c_str(), std::ios::in | std::ios::binary);
    int N;
    int NcInFile;
    double L, a, temp;

    if (!InStream.good()) {
        messager_ << "[WilsonLineIO::read]: File " << fileName.c_str()
                  << " does not exist!";
        messager_.flush("error");
        exit(1);
    }

    if (!InStream.is_open()) return;

    // READING IN PARAMETERS
    InStream.read(reinterpret_cast<char *>(&N), sizeof(int));
    InStream.read(reinterpret_cast<char *>(&NcInFile), sizeof(int));
    InStream.read(reinterpret_cast<char *>(&L), sizeof(double));
    InStream.read(reinterpret_cast<char *>(&a), sizeof(double));
    InStream.read(reinterpret_cast<char *>(&temp), sizeof(double));

    if (N != param->lattice.size) {
        messager_ << "[WilsonLineIO::read]: wrong lattice "
                     "size, data is "
                  << N << " but you have specified " << param->lattice.size;
        messager_.flush("error");
        exit(1);
    }
    if (std::abs(L - param->lattice.L) > 1e-5) {
        messager_ << "[WilsonLineIO::read]: wrong grid length, "
                     "data has "
                  << L << " but you have specified " << param->lattice.L;
        messager_.flush("error");
        exit(1);
    }

    // READING ACTUAL DATA
    double ValueBuffer;
    int INPUT_CTR = 0;
    double re, im;
    re = 0.;
    im = 0.;

    while (
        InStream.read(reinterpret_cast<char *>(&ValueBuffer), sizeof(double))) {
        if (INPUT_CTR % 2 == 0)  // this is the real part
        {
            re = ValueBuffer;
        } else  // this is the imaginary part, write then to
                // variable //
        {
            im = ValueBuffer;
            int TEMPINDX = ((INPUT_CTR - 1) / 2);
            int PositionIndx = TEMPINDX / 9;

            // PositionIndx enumerates the sites in the writer's order
            // (see WilsonLineIO::write), i.e. latticeIndex().
            int ixRaw = latticeX(PositionIndx, N);
            int iy = latticeY(PositionIndx, N);

            double bb = param->event.b;
            a = L / static_cast<double>(N);

            // shift here by half an impact parameter
            double xtemp =
                isProjectile ? (a * ixRaw - bb / 2.) : (a * ixRaw + bb / 2.);

            int ix = round(xtemp / a);

            int MatrixIndx = TEMPINDX - PositionIndx * 9;
            int j = MatrixIndx / 3;
            int k = MatrixIndx - j * 3;

            int indx = latticeIndex(ix, iy, N);
            if (indx >= N * N || indx < 0) {
                if (bb == 0) {
                    messager_ << "[WilsonLineIO::read]: datafile " << fileName
                              << " has an element " << indx << " (iy=" << iy
                              << ", ix=" << ix << "), but the grid is N=" << N
                              << ". Element is (" << re << " + " << im
                              << "i), skipping it.";
                    messager_.flush("warning");
                }
                INPUT_CTR++;
                continue;
            }
            U[indx].set(j, k, complex<double>(re, im));
        }
        INPUT_CTR++;
    }

    InStream.close();
}

void WilsonLineIO::writeTrainingData(Lattice *lat, Parameters *param) {
    const int N = param->lattice.size;
    const double L = param->lattice.L;
    const double a = L / static_cast<double>(N);

    // Payload: [beam, real_or_imag, x, y, row, col], C-order.
    constexpr int nBeams = 2;
    const std::size_t matrixElements =
        static_cast<std::size_t>(N) * N * Nc * Nc;
    std::vector<float> payload(
        static_cast<std::size_t>(nBeams) * 2 * matrixElements);

    for (int beam = 0; beam < nBeams; ++beam) {
        const std::size_t realOffset =
            static_cast<std::size_t>(2 * beam) * matrixElements;
        const std::size_t imagOffset = realOffset + matrixElements;
        for (int x = 0; x < N; ++x) {
            for (int y = 0; y < N; ++y) {
                const int pos = lat->positionFromXY(x, y);
                const Matrix &matrix = (beam == 0) ? lat->U[pos] : lat->U2[pos];
                const std::complex<double> *elements = matrix.data();
                const std::size_t siteOffset =
                    static_cast<std::size_t>(pos) * Nc * Nc;
                for (int row = 0; row < Nc; ++row) {
                    for (int col = 0; col < Nc; ++col) {
                        const std::size_t element =
                            static_cast<std::size_t>(row) * Nc + col;
                        payload[realOffset + siteOffset + element] =
                            static_cast<float>(elements[element].real());
                        payload[imagOffset + siteOffset + element] =
                            static_cast<float>(elements[element].imag());
                    }
                }
            }
        }
    }

    const std::uint16_t endianProbe = 1;
    if (*reinterpret_cast<const unsigned char *>(&endianProbe) != 1) {
        throw std::runtime_error(
            "WilsonLineIO::writeTrainingData requires a little-endian host");
    }

    std::stringstream metadata;
    metadata << std::setprecision(17)
             << "{\"format\":\"ipglasma-initial-wilson-lines\","
             << "\"version\":1,"
             << "\"dtype\":\"<f4\","
             << "\"shape\":[2,2," << N << "," << N << "," << Nc << "," << Nc
             << "],"
             << "\"axis_order\":[\"beam\",\"complex_part\",\"x\",\"y\","
                "\"row\",\"col\"],"
             << "\"fields\":[\"VA\",\"VB\"],"
             << "\"complex_part\":[\"real\",\"imag\"],"
             << "\"native_site_index\":\"pos=x*N+y\","
             << "\"event_id\":" << param->event.eventId << ","
             << "\"N\":" << N << ","
             << "\"Nc\":" << Nc << ","
             << "\"L_fm\":" << L << ","
             << "\"a_fm\":" << a << ","
             << "\"rapidity\":" << param->colorCharge.rapidity() << "}";
    const std::string metadataString = metadata.str();

    std::stringstream filename;
    filename << "initialWilsonLines" << param->event.eventId << ".ipgw";
    std::ofstream output(
        filename.str().c_str(),
        std::ios::out | std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "could not open initial-Wilson snapshot " + filename.str());
    }
    const char magic[8] = {'I', 'P', 'G', 'W', 'I', 'L', '1', '\0'};
    const std::uint64_t metadataBytes =
        static_cast<std::uint64_t>(metadataString.size());
    output.write(magic, sizeof(magic));
    output.write(
        reinterpret_cast<const char *>(&metadataBytes), sizeof(metadataBytes));
    output.write(metadataString.data(), metadataString.size());
    output.write(
        reinterpret_cast<const char *>(payload.data()),
        static_cast<std::streamsize>(payload.size() * sizeof(float)));
    output.close();
    if (!output) {
        throw std::runtime_error(
            "failed while writing initial-Wilson snapshot " + filename.str());
    }
    messager_
        << "[WilsonLineIO::writeTrainingData]: Wrote incoming Wilson lines to "
        << filename.str();
    messager_.flush("info");
}

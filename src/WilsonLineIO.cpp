// WilsonLineIO.cpp is part of the IP-Glasma solver.

#include "WilsonLineIO.h"

#include <algorithm>
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

namespace {
/**
 * The number in the names of a nucleus' Wilson-line and geometry files:
 * \f$2(sN + i) + j\f$, see WilsonLineIO::fileName().
 * \param[in] param Simulation parameters.
 * \param[in] nucleus Which nucleus.
 * \return The number.
 */
unsigned long long fileNumber(Parameters *param, NucleusRole nucleus) {
    const int iA = (nucleus == NucleusRole::Projectile) ? 1 : 2;
    // 2 (seed N + eventId) + iA with N events in the run: different for
    // every nucleus and event of a run, and for runs with different seeds
    // and the same N
    const unsigned long long eventsPerRun =
        static_cast<unsigned long long>(param->run.eventsPerRank)
        * param->run.MPISize;
    return 2 * (param->random.seed * eventsPerRun + param->event.eventId) + iA;
}

/**
 * The value of \p key in the flat JSON object of a geometry file's header
 * (as written by WilsonLineIO::writeGeometry()), without quotes.
 * \param[in] json The header.
 * \param[in] key The key.
 * \return The value, or an empty string if \p key is missing.
 */
std::string jsonValue(const std::string &json, const std::string &key) {
    const std::string pattern = "\"" + key + "\":";
    const std::size_t start = json.find(pattern);
    if (start == std::string::npos) return "";
    const std::size_t begin = start + pattern.size();
    if (begin < json.size() && json[begin] == '"') {
        const std::size_t end = json.find('"', begin + 1);
        return json.substr(begin + 1, end - begin - 1);
    }
    const std::size_t end = json.find_first_of(",}", begin);
    return json.substr(begin, end - begin);
}

/// Magic string at the start of a geometry file.
const char kGeometryMagic[8] = {'I', 'P', 'G', 'G', 'E', 'O', '1', '\0'};

/**
 * Checks the host's byte order; the geometry files are little-endian by
 * definition.
 * \return `true` if this host is little-endian.
 */
bool littleEndianHost() {
    const std::uint16_t probe = 1;
    return *reinterpret_cast<const unsigned char *>(&probe) == 1;
}
}  // namespace

std::string WilsonLineIO::fileName(
    Parameters *param, const double x, NucleusRole nucleus, int format) {
    std::stringstream Vname;
    Vname << param->wilsonLines.wilsonLinePath << "/WilsonLine";
    if (x >= 0) Vname << "_x_" << std::scientific << std::setprecision(5) << x;
    Vname << "_" << fileNumber(param, nucleus);

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

        double temp = param->colorCharge.rapidity;

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

std::string WilsonLineIO::geometryFileName(
    Parameters *param, NucleusRole nucleus) {
    std::stringstream name;
    name << param->wilsonLines.wilsonLinePath << "/WilsonLineGeometry_"
         << fileNumber(param, nucleus);
    return name.str();
}

void WilsonLineIO::read(Lattice *lat, Parameters *param, int format) {
    IPG_PROFILE_SCOPE("initialization.read_wilson_lines");
    if (!isValidFormat(format)) {
        messager_ << "[WilsonLineIO::read]: Unknown format " << format
                  << " when reading the initial Wilson lines, supported "
                     "formats: 1,2";
        messager_.flush("error");
        exit(1);
    }

    string VOne_name = fileName(
        param, param->xBeforeJimwlk(NucleusRole::Projectile),
        NucleusRole::Projectile, format);
    string VTwo_name = fileName(
        param, param->xBeforeJimwlk(NucleusRole::Target), NucleusRole::Target,
        format);

    messager_ << "[WilsonLineIO::read]: Reading Wilson lines from files "
              << VOne_name << " and " << VTwo_name;
    messager_.flush("info");

    if (format == 1) {
        readText(VOne_name, param, lat->U);
        readText(VTwo_name, param, lat->U2);
    } else if (format == 2) {
        readBinary(VOne_name, param, lat->U);
        readBinary(VTwo_name, param, lat->U2);
    }

    messager_ << "[WilsonLineIO::read]: Wilson lines V_A and V_B set on rank "
              << param->run.MPIRank << ". ";
    messager_.flush("info");
}

void WilsonLineIO::readText(
    const std::string &fileName, Parameters *param, std::vector<Matrix> &U) {
    const int N = param->lattice.size;

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

            U[latticeIndex(i, j, N)] = temp;
        }
    }

    fin.close();
}

void WilsonLineIO::readBinary(
    const std::string &fileName, Parameters *param, std::vector<Matrix> &U) {
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
            const int ix = latticeX(PositionIndx, N);
            const int iy = latticeY(PositionIndx, N);

            int MatrixIndx = TEMPINDX - PositionIndx * 9;
            int j = MatrixIndx / 3;
            int k = MatrixIndx - j * 3;

            if (ix >= N) {
                messager_ << "[WilsonLineIO::read]: datafile " << fileName
                          << " has more than N^2 sites (N=" << N
                          << "); skipping element (" << re << " + " << im
                          << "i).";
                messager_.flush("warning");
                INPUT_CTR++;
                continue;
            }
            U[latticeIndex(ix, iy, N)].set(j, k, complex<double>(re, im));
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
             << "\"rapidity\":" << param->colorCharge.rapidity << "}";
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

void WilsonLineIO::writeGeometry(
    Lattice *lat, Parameters *param, NucleusRole nucleus,
    const std::vector<ReturnValue> &nucleons) {
    if (!littleEndianHost()) {
        throw std::runtime_error(
            "WilsonLineIO::writeGeometry requires a little-endian host");
    }
    const bool isProjectile = (nucleus == NucleusRole::Projectile);
    const int N = param->lattice.size;
    const std::size_t sites = static_cast<std::size_t>(N) * N;

    // nucleons (x, y, z, proton), then the g^2 mu^2 and T_p maps
    std::vector<double> payload;
    payload.reserve(4 * nucleons.size() + 2 * sites);
    for (const ReturnValue &nucleon : nucleons) {
        payload.push_back(nucleon.x);
        payload.push_back(nucleon.y);
        payload.push_back(nucleon.z);
        payload.push_back(nucleon.proton ? 1. : 0.);
    }
    for (std::size_t pos = 0; pos < sites; pos++) {
        payload.push_back(
            isProjectile ? lat->cells[pos]->getg2mu2A()
                         : lat->cells[pos]->getg2mu2B());
    }
    for (std::size_t pos = 0; pos < sites; pos++) {
        payload.push_back(
            isProjectile ? lat->cells[pos]->getTpA()
                         : lat->cells[pos]->getTpB());
    }

    std::stringstream metadata;
    metadata << std::setprecision(17)
             << "{\"format\":\"ipglasma-nucleus-geometry\","
             << "\"version\":1,"
             << "\"dtype\":\"<f8\","
             << "\"nucleus\":\"" << (isProjectile ? "projectile" : "target")
             << "\","
             << "\"species\":\""
             << (isProjectile ? param->collision.projectile
                              : param->collision.target)
             << "\","
             << "\"nucleons\":" << nucleons.size() << ","
             << "\"N\":" << N << ","
             << "\"L_fm\":" << param->lattice.L << ","
             << "\"g\":" << param->coupling.g << ","
             << "\"QsMuRatio\":" << param->colorCharge.QsMuRatio << ","
             << "\"blocks\":[\"nucleons[nucleons][x_fm,y_fm,z_fm,proton]\","
                "\"g2mu2[N*N]\",\"Tp_GeV2[N*N]\"],"
             << "\"native_site_index\":\"pos=x*N+y\"}";
    const std::string metadataString = metadata.str();

    const std::string name = geometryFileName(param, nucleus);
    std::ofstream output(
        name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    const std::uint64_t metadataBytes =
        static_cast<std::uint64_t>(metadataString.size());
    output.write(kGeometryMagic, sizeof(kGeometryMagic));
    output.write(
        reinterpret_cast<const char *>(&metadataBytes), sizeof(metadataBytes));
    output.write(metadataString.data(), metadataString.size());
    output.write(
        reinterpret_cast<const char *>(payload.data()),
        static_cast<std::streamsize>(payload.size() * sizeof(double)));
    output.close();
    if (!output) {
        messager_ << "[WilsonLineIO::writeGeometry]: could not write " << name
                  << ". Exiting.";
        messager_.flush("error");
        exit(1);
    }
}

NucleusGeometry WilsonLineIO::readGeometry(
    Lattice *lat, Parameters *param, NucleusRole nucleus) {
    const bool isProjectile = (nucleus == NucleusRole::Projectile);
    const std::string name = geometryFileName(param, nucleus);
    auto fail = [this, &name](const std::string &reason) {
        messager_ << "[WilsonLineIO::readGeometry]: " << name << ": " << reason
                  << ". Exiting.";
        messager_.flush("error");
        exit(1);
    };
    if (!littleEndianHost()) fail("reading needs a little-endian host");

    std::ifstream input(name.c_str(), std::ios::in | std::ios::binary);
    if (!input) {
        fail(
            "not found (reading Wilson lines needs the geometry file written "
            "with them)");
    }
    char magic[8];
    std::uint64_t metadataBytes = 0;
    input.read(magic, sizeof(magic));
    input.read(reinterpret_cast<char *>(&metadataBytes), sizeof(metadataBytes));
    if (!input || !std::equal(magic, magic + sizeof(magic), kGeometryMagic)) {
        fail("not a geometry file");
    }
    std::string metadata(static_cast<std::size_t>(metadataBytes), '\0');
    input.read(&metadata[0], static_cast<std::streamsize>(metadataBytes));

    // the geometry must belong to this lattice, coupling and species
    const int N = param->lattice.size;
    const std::string species =
        isProjectile ? param->collision.projectile : param->collision.target;
    if (jsonValue(metadata, "nucleus")
        != (isProjectile ? "projectile" : "target")) {
        fail("holds the other nucleus");
    }
    if (jsonValue(metadata, "species") != species) {
        fail(
            "holds " + jsonValue(metadata, "species") + ", but this run has "
            + species);
    }
    if (std::atoi(jsonValue(metadata, "N").c_str()) != N) {
        fail(
            "has size " + jsonValue(metadata, "N") + ", but this run has "
            + std::to_string(N));
    }
    if (std::abs(
            std::atof(jsonValue(metadata, "L_fm").c_str()) - param->lattice.L)
        > 1e-5) {
        fail("has another L than this run");
    }
    if (!PhysConst::isClose(
            std::atof(jsonValue(metadata, "g").c_str()), param->coupling.g)) {
        fail("was written with another g than this run");
    }

    NucleusGeometry geometry;
    geometry.QsMuRatio = std::atof(jsonValue(metadata, "QsMuRatio").c_str());
    const std::size_t count = static_cast<std::size_t>(
        std::atoll(jsonValue(metadata, "nucleons").c_str()));
    const std::size_t sites = static_cast<std::size_t>(N) * N;
    std::vector<double> payload(4 * count + 2 * sites);
    input.read(
        reinterpret_cast<char *>(payload.data()),
        static_cast<std::streamsize>(payload.size() * sizeof(double)));
    if (!input) fail("is truncated");

    geometry.nucleons.resize(count);
    for (std::size_t i = 0; i < count; i++) {
        ReturnValue &nucleon = geometry.nucleons[i];
        nucleon.x = payload[4 * i];
        nucleon.y = payload[4 * i + 1];
        nucleon.z = payload[4 * i + 2];
        nucleon.proton = payload[4 * i + 3] != 0.;
        nucleon.phi = 0.;
        nucleon.collided = 0;
    }
    const double *g2mu2 = payload.data() + 4 * count;
    const double *Tp = g2mu2 + sites;
    for (std::size_t pos = 0; pos < sites; pos++) {
        if (isProjectile) {
            lat->cells[pos]->setg2mu2A(g2mu2[pos]);
            lat->cells[pos]->setTpA(Tp[pos]);
        } else {
            lat->cells[pos]->setg2mu2B(g2mu2[pos]);
            lat->cells[pos]->setTpB(Tp[pos]);
        }
    }
    messager_ << "[WilsonLineIO::readGeometry]: Read " << count
              << " nucleons and the color-charge densities from " << name;
    messager_.flush("info");
    return geometry;
}

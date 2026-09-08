// MyEigen.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "MyEigen.h"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Instrumentation.h"
#include "Phys_consts.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_eigen.h"

using PhysConst::hbarc;
using PhysConst::small_eps;

using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::vector;

namespace {

const std::size_t kTextOutputBufferBytes = 4u * 1024u * 1024u;

void openBufferedTextOutput(
    ofstream &output, vector<char> &buffer, const string &filename) {
    output.rdbuf()->pubsetbuf(
        buffer.data(), static_cast<std::streamsize>(buffer.size()));
    output.open(filename.c_str(), std::ios::out);
    if (!output) {
        throw std::runtime_error("could not open output file " + filename);
    }
}

void closeBufferedTextOutput(ofstream &output, const string &filename) {
    output.close();
    if (!output) {
        throw std::runtime_error(
            "failed while writing output file " + filename);
    }
}

bool binaryTmunuEnabled(Parameters *param) {
    const bool inputDefault = param->getWriteTmunuBinary() != 0;
    const char *value = std::getenv("IPGLASMA_BINARY_TMUNU");
    if (value == NULL || value[0] == '\0') return inputDefault;

    // An explicitly set environment variable remains a convenient runtime
    // override for benchmarking and existing launch scripts.
    const string text(value);
    return !(
        text == "0" || text == "false" || text == "FALSE" || text == "off"
        || text == "OFF" || text == "no" || text == "NO");
}

bool littleEndianHost() {
    const std::uint16_t probe = 1;
    return *reinterpret_cast<const unsigned char *>(&probe) == 1;
}

void writeUint32LittleEndian(ofstream &output, std::uint32_t value) {
    const unsigned char bytes[4] = {
        static_cast<unsigned char>(value & 0xffu),
        static_cast<unsigned char>((value >> 8) & 0xffu),
        static_cast<unsigned char>((value >> 16) & 0xffu),
        static_cast<unsigned char>((value >> 24) & 0xffu)};
    output.write(reinterpret_cast<const char *>(bytes), sizeof(bytes));
}

void openTmunuBinaryOutput(
    ofstream &output, const string &filename, int hx, int hy, int heta,
    double tauFm, double deta, double dxFm, int eventId) {
    if (!littleEndianHost()) {
        throw std::runtime_error(
            "binary Tmunu output currently requires a little-endian host");
    }
    output.open(
        filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("could not open output file " + filename);
    }

    std::ostringstream metadata;
    metadata << std::setprecision(17) << "{\"format\":\"ipglasma-tmunu\","
             << "\"version\":1,"
             << "\"dtype\":\"<f4\","
             << "\"shape\":[" << hy << "," << hx << ",10],"
             << "\"axis_order\":[\"y\",\"x\",\"component\"],"
             << "\"components\":[\"T00\",\"Txx\",\"Tyy\","
                "\"tau2_Tetaeta\",\"neg_T0x\",\"neg_T0y\","
                "\"neg_tau_T0eta\",\"neg_Txy\",\"neg_tau_Tyeta\","
                "\"neg_tau_Txeta\"],"
             << "\"tau_fm\":" << tauFm << ","
             << "\"eta_points\":" << heta << ","
             << "\"deta\":" << deta << ","
             << "\"dx_fm\":" << dxFm << ","
             << "\"dy_fm\":" << dxFm << ","
             << "\"event_id\":" << eventId << "}";
    const string metadataText = metadata.str();
    if (metadataText.size() > 0xffffffffu) {
        throw std::runtime_error("binary Tmunu metadata is unexpectedly large");
    }

    static const char magic[8] = {'I', 'P', 'G', 'T', 'M', 'U', '0', '1'};
    output.write(magic, sizeof(magic));
    writeUint32LittleEndian(
        output, static_cast<std::uint32_t>(metadataText.size()));
    output.write(
        metadataText.data(), static_cast<std::streamsize>(metadataText.size()));
    if (!output) {
        throw std::runtime_error(
            "failed while writing Tmunu header " + filename);
    }
}

}  // namespace

//**************************************************************************
// MyEigen class.

void MyEigen::flowVelocity4D(
    Lattice *lat, Parameters *param, int it, bool finalFlag) {
    flowVelocity4DImpl(lat, param, it, finalFlag, false);
}

void MyEigen::writeTmunu4D(Lattice *lat, Parameters *param, int it) {
    flowVelocity4DImpl(lat, param, it, false, true);
}

void MyEigen::flowVelocity4DImpl(
    Lattice *lat, Parameters *param, int it, bool finalFlag, bool tmunuOnly) {
    int N = param->getSize();
    int count;
    double L = param->getL();
    double a = L / N;  // lattice spacing in fm
    double x, y;
    double dtau = param->getdtau();
    double averageux = 0.;
    double averageuy = 0.;
    double averageueta = 0.;
    double averageeps = 0.;
    count = 0;

    if (!tmunuOnly) {
        // The per-cell flow-velocity solve is now independent across cells (the
        // velocity carryover that previously serialized it has been removed via
        // the rest-frame reset below), so it is parallelized over the lattice.
        // Each thread allocates its own GSL workspace/eigenvalue/eigenvector
        // objects once. All per-cell writes touch only cells[pos], so every
        // DATA-FILE output is bit-for-bit identical to the serial version
        // regardless of thread count. The four diagnostic averages (printed to
        // stdout only) use an OpenMP reduction whose summation order differs
        // from serial, so those logged numbers may differ in their last bits --
        // no data file is affected.
#pragma omp parallel reduction( \
        + : averageux, averageuy, averageueta, averageeps, count)
        {
            gsl_vector_complex *eval_ws = gsl_vector_complex_alloc(4);
            gsl_matrix_complex *evec_ws = gsl_matrix_complex_alloc(4, 4);
            gsl_eigen_nonsymmv_workspace *w_ws = gsl_eigen_nonsymmv_alloc(4);

#pragma omp for
            for (int posLoop = 0; posLoop < N * N; posLoop++) {
                const int si = posLoop / N;
                const int sj = posLoop % N;
                const int pos = posLoop;
                // Flow velocity defaults to the local rest frame (0,0,0,1):
                // zero spatial flow, u^tau = 1. (See note in the serial
                // version: this also removes the old cross-cell carryover bug.)
                double ux = 0., uy = 0., ueta = 0., utau = 1.;
                double eps = 0.;
                int changeSign;
                int foundU;
                gsl_complex square;
                gsl_complex factor;
                gsl_complex euklidiansquare;
                gsl_complex z_aux;
                gsl_complex tau2;

                GSL_SET_COMPLEX(&square, 0, 0);
                // one upper, one lower index
                double data[] = {
                    lat->cells[pos]->getTtautau(),
                    -lat->cells[pos]->getTtaux(),
                    -lat->cells[pos]->getTtauy(),
                    -(it * dtau * a) * (it * dtau * a)
                        * lat->cells[pos]->getTtaueta(),
                    lat->cells[pos]->getTtaux(),
                    -lat->cells[pos]->getTxx(),
                    -lat->cells[pos]->getTxy(),
                    -(it * dtau * a) * (it * dtau * a)
                        * lat->cells[pos]->getTxeta(),
                    lat->cells[pos]->getTtauy(),
                    -lat->cells[pos]->getTxy(),
                    -lat->cells[pos]->getTyy(),
                    -(it * dtau * a) * (it * dtau * a)
                        * lat->cells[pos]->getTyeta(),
                    lat->cells[pos]->getTtaueta(),
                    -lat->cells[pos]->getTxeta(),
                    -lat->cells[pos]->getTyeta(),
                    -(it * dtau * a) * (it * dtau * a)
                        * lat->cells[pos]->getTetaeta()};

                gsl_matrix_view m =
                    gsl_matrix_view_array(data, 4, 4);  // matrix

                gsl_vector_complex *eval = eval_ws;
                gsl_matrix_complex *evec = evec_ws;

                gsl_eigen_nonsymmv(
                    &m.matrix, eval, evec,
                    w_ws);  // solve for eigenvalues and eigenvectors (without
                            // 'v' only compute eigenvalues)

                // set to 'zero'
                lat->cells[pos]->setEpsilon(lat->cells[pos]->getTtautau());
                lat->cells[pos]->setutau(1.);
                lat->cells[pos]->setux(0.);
                lat->cells[pos]->setuy(0.);
                lat->cells[pos]->setueta(0.);

                {  // output:
                    int i;
                    foundU = 0;
                    for (i = 0; i < 4; i++) {
                        gsl_complex eval_i = gsl_vector_complex_get(eval, i);
                        gsl_vector_complex_view evec_i =
                            gsl_matrix_complex_column(evec, i);

                        GSL_SET_COMPLEX(&square, 0, 0);
                        GSL_SET_COMPLEX(
                            &tau2, a * it * dtau * a * it * dtau, 0);
                        GSL_SET_COMPLEX(&euklidiansquare, 0, 0);

                        for (int j = 0; j < 4; ++j) {
                            gsl_complex z =
                                gsl_vector_complex_get(&evec_i.vector, j);
                            z_aux = gsl_complex_mul(tau2, z);
                            euklidiansquare = gsl_complex_add(
                                euklidiansquare, gsl_complex_mul(z, z));

                            if (j == 0)
                                square = gsl_complex_add(
                                    square, gsl_complex_mul(z, z));
                            else if (j < 3)
                                square = gsl_complex_sub(
                                    square, gsl_complex_mul(z, z));
                            else
                                square = gsl_complex_sub(
                                    square, gsl_complex_mul(z_aux, z));
                        }

                        GSL_SET_COMPLEX(
                            &factor,
                            sqrt(abs(
                                GSL_REAL(euklidiansquare) / GSL_REAL(square))),
                            0);
                        if (GSL_REAL(square) > 0) {
                            eps = GSL_REAL(eval_i);
                            if (abs(GSL_IMAG(eval_i)) > 0.001 && si > 0
                                && sj > 0 && si < N - 5 && sj < N - 5) {
                                eps = lat->cells[pos]->getTtautau();
                            }
                        }

                        GSL_SET_COMPLEX(&square, 0, 0);

                        for (int j = 0; j < 4; ++j) {
                            gsl_complex z =
                                gsl_vector_complex_get(&evec_i.vector, j);
                            z = gsl_complex_mul(z, factor);
                            z_aux = gsl_complex_mul(tau2, z);

                            if (j == 0)
                                square = gsl_complex_add(
                                    square, gsl_complex_mul(z, z));
                            else if (j < 3)
                                square = gsl_complex_sub(
                                    square, gsl_complex_mul(z, z));
                            else
                                square = gsl_complex_sub(
                                    square, gsl_complex_mul(z_aux, z));
                        }
                        changeSign = 0;
                        // for the time-like eigenvector do the following (this
                        // is the flow velocity)
                        if (GSL_REAL(square) > 0) {
                            foundU += 1;
                            for (int j = 0; j < 4; ++j) {
                                gsl_complex z =
                                    gsl_vector_complex_get(&evec_i.vector, j);
                                z = gsl_complex_mul(z, factor);

                                if (j == 0 && GSL_REAL(z) < 0) {
                                    changeSign = 1;
                                    GSL_SET_COMPLEX(
                                        &z, -1. * GSL_REAL(z),
                                        -1. * GSL_IMAG(z));
                                }

                                if (j > 0 && changeSign == 1) {
                                    GSL_SET_COMPLEX(
                                        &z, -1. * GSL_REAL(z),
                                        -1. * GSL_IMAG(z));
                                }

                                if (j == 0) {
                                    if (eps > 0.1)
                                        utau = GSL_REAL(z);
                                    else
                                        utau = 1.;
                                }
                                if (j == 1) {
                                    if (eps > 0.1)
                                        ux = GSL_REAL(z);
                                    else
                                        ux = 0.;
                                }
                                if (j == 2) {
                                    if (eps > 0.1)
                                        uy = GSL_REAL(z);
                                    else
                                        uy = 0.;
                                }
                                if (j == 3) {
                                    if (eps > 0.1)
                                        ueta = GSL_REAL(z);
                                    else
                                        ueta = 0.;
                                }
                                if (abs(GSL_IMAG(z)) > 0.001 && eps > 0.001
                                    && si > 10 && sj > 10 && si < N - 10
                                    && sj < N - 10) {
                                    utau = 1.;
                                    ux = 0.;
                                    uy = 0.;
                                    ueta = 0.;
                                    eps = lat->cells[pos]->getTtautau();
                                }
                            }

                            lat->cells[pos]->setutau(utau);
                            lat->cells[pos]->setux(ux);
                            lat->cells[pos]->setuy(uy);
                            lat->cells[pos]->setueta(ueta);
                            lat->cells[pos]->setEpsilon(eps);

                            // //clean up numerical noise outside the
                            // interaction region if
                            // (lat->cells[pos]->getg2mu2A() < 1e-12 ||
                            //     lat->cells[pos]->getg2mu2B() < 1e-12)
                            //   lat->cells[pos]->setEpsilon(0.);

                            averageux += ux * ux * eps;
                            averageuy += uy * uy * eps;
                            averageueta += ueta * ueta * eps * it * dtau * a
                                           * it * dtau * a;
                            averageeps += eps;

                            count++;
                        }
                    }
                }

                // write Tmunu in case no u was found
                if (foundU == 0) {
                    if (si == N / 2 && sj == N / 2) {
                        cout << si << " " << sj << endl << endl;
                        cout << lat->cells[pos]->getTtautau() << " "
                             << lat->cells[pos]->getTtaux() << " "
                             << lat->cells[pos]->getTtauy() << " "
                             << lat->cells[pos]->getTtaueta() << endl;
                        cout << lat->cells[pos]->getTtaux() << " "
                             << lat->cells[pos]->getTxx() << " "
                             << lat->cells[pos]->getTxy() << " "
                             << lat->cells[pos]->getTxeta() << endl;
                        cout << lat->cells[pos]->getTtauy() << " "
                             << lat->cells[pos]->getTxy() << " "
                             << lat->cells[pos]->getTyy() << " "
                             << lat->cells[pos]->getTyeta() << endl;
                        cout << lat->cells[pos]->getTtaueta() << " "
                             << lat->cells[pos]->getTxeta() << " "
                             << lat->cells[pos]->getTyeta() << " "
                             << lat->cells[pos]->getTetaeta() << endl;
                        cout << "ux=" << ux << endl;
                        cout << "uy=" << uy << endl;
                        cout << "ueta=" << ueta << endl;
                    }
                }

                // compute pi^{\mu\nu}
                if (utau == 1 && ux == 0 && uy == 0 && ueta == 0) {
                    lat->cells[pos]->setpitautau(0.);
                    lat->cells[pos]->setpixx(0.);
                    lat->cells[pos]->setpiyy(0.);
                    lat->cells[pos]->setpietaeta(0.);

                    lat->cells[pos]->setpitaux(0.);
                    lat->cells[pos]->setpitauy(0.);
                    lat->cells[pos]->setpitaueta(0.);

                    lat->cells[pos]->setpixeta(0.);
                    lat->cells[pos]->setpixy(0.);
                    lat->cells[pos]->setpiyeta(0.);
                } else {
                    lat->cells[pos]->setpitautau(
                        lat->cells[pos]->getTtautau()
                        - 4. / 3. * eps * utau * utau + eps / 3.);
                    lat->cells[pos]->setpixx(
                        lat->cells[pos]->getTxx() - 4. / 3. * eps * ux * ux
                        - eps / 3.);
                    lat->cells[pos]->setpiyy(
                        lat->cells[pos]->getTyy() - 4. / 3. * eps * uy * uy
                        - eps / 3.);
                    lat->cells[pos]->setpietaeta(
                        lat->cells[pos]->getTetaeta()
                        - 4. / 3. * eps * ueta * ueta
                        - eps / 3. / it / dtau / a / it / dtau / a);

                    lat->cells[pos]->setpitaux(
                        lat->cells[pos]->getTtaux()
                        - 4. / 3. * eps * utau * ux);
                    lat->cells[pos]->setpitauy(
                        lat->cells[pos]->getTtauy()
                        - 4. / 3. * eps * utau * uy);
                    lat->cells[pos]->setpitaueta(
                        lat->cells[pos]->getTtaueta()
                        - 4. / 3. * eps * utau * ueta);

                    lat->cells[pos]->setpixeta(
                        lat->cells[pos]->getTxeta()
                        - 4. / 3. * eps * ux * ueta);
                    lat->cells[pos]->setpixy(
                        lat->cells[pos]->getTxy() - 4. / 3. * eps * ux * uy);
                    lat->cells[pos]->setpiyeta(
                        lat->cells[pos]->getTyeta()
                        - 4. / 3. * eps * uy * ueta);
                }
            }  // omp for over posLoop

            gsl_eigen_nonsymmv_free(w_ws);
            gsl_vector_complex_free(eval_ws);
            gsl_matrix_complex_free(evec_ws);
        }  // omp parallel

        cout << it * dtau * a << " average u^x=" << sqrt(averageux / averageeps)
             << endl;
        cout << it * dtau * a << " average u^y=" << sqrt(averageuy / averageeps)
             << endl;
        cout << it * dtau * a
             << " average tau u^eta=" << sqrt(averageueta / averageeps) << endl;
    }

    //  double maxtime = param->getMaxtime(); // maxtime is in fm
    //  int itmax = static_cast<int>(floor(maxtime/(a*dtau)+1e-10));

    double Etot = 0.;

    // output for hydro
    if (param->getWriteOutputs() <= 0) return;

    double alphas = param->getalphas();
    double g = param->getg();
    double gfactor;
    int hx = param->getSizeOutput();
    int hy = hx;
    int heta = param->getEtaSizeOutput();
    double hL = param->getLOutput();
    double deta = param->getDetaOutput();
    double c = param->getc();
    double muZero = param->getMuZero();

    if (param->getRunningCoupling()) {
        // run with average Q_s only ! local makes no sense here (stuff has
        // moved in the mean time)
        alphas = 4. * M_PI
                 / (9.
                    * log(
                        pow(pow(muZero / 0.2, 2. / c)
                                + pow(
                                    param->getRunWithThisFactorTimesQs()
                                        * param->getAverageQsmin() / 0.2,
                                    2. / c),
                            c)));
        gfactor = g * g / (4. * M_PI * alphas);
    } else {
        gfactor = 1.;
    }

    if (hL > L) {
        cout << "WARNING: hydro grid length larger than the computed one."
             << endl;
    }

    int xpos, ypos, xposUp, yposUp, pos1, pos2, pos3;
    double fracx, fracy, x1, x2;
    double xlow, ylow;
    int pos4;
    double resultE, resultutau, resultux, resultuy, resultueta;
    double resultpi00, resultpi0x, resultpi0y, resultpi0eta;
    double resultpixy, resultpixeta, resultpiyeta, resultpixx, resultpiyy,
        resultpietaeta;
    double g2mu2A, g2mu2B;

    double tau0 = it * dtau * a;

    double ha;
    ha = hL / static_cast<double>(hx);

    stringstream streuH_name;
    if (finalFlag) {
        streuH_name << "epsilon-u-Hydro-TauHydro-" << param->getEventId()
                    << ".dat";
    } else {
        streuH_name << "epsilon-u-Hydro-t" << it * dtau * a << "-"
                    << param->getEventId() << ".dat";
    }

    // stringstream strEtot_name;
    // strEtot_name << "Etot-t" << it*dtau*a << "-" << param->getEventId() <<
    // ".dat"; string Etot_name; Etot_name = strEtot_name.str();

    if (!tmunuOnly && param->getWriteOutputs() % 2 == 1) {
        IPG_PROFILE_SCOPE("output.hydro_text");
        const string outputFilename = streuH_name.str();
        vector<char> outputBuffer(kTextOutputBufferBytes);
        ofstream foutEps2;
        openBufferedTextOutput(foutEps2, outputBuffer, outputFilename);
        //      ofstream foutEtot(Etot_name.c_str(),ios::out);

        foutEps2 << "# dummy " << 1 << " etamax= " << heta << " xmax= " << hx
                 << " ymax= " << hy << " deta= " << deta << " dx= " << ha
                 << " dy= " << ha << " tau= " << tau0 << '\n';

        for (int ieta = 0; ieta < heta; ieta++)  // loop over all positions
        {
            for (int ix = 0; ix < hx; ix++)  // loop over all positions
            {
                for (int iy = 0; iy < hy; iy++) {
                    x = -hL / 2. + ha * ix;
                    y = -hL / 2. + ha * iy;

                    if (std::abs(x) < (L / 2. - 0.5)
                        && std::abs(y) < (L / 2. - 0.5)) {
                        // if (std::abs(x) < L / 2. && std::abs(y) < L / 2.) {
                        xpos = static_cast<int>(
                            floor((x + L / 2.) / a + 0.0000000001));
                        ypos = static_cast<int>(
                            floor((y + L / 2.) / a + 0.0000000001));

                        if (xpos < N - 1)
                            xposUp = xpos + 1;
                        else
                            xposUp = xpos;

                        if (ypos < N - 1)
                            yposUp = ypos + 1;
                        else
                            yposUp = ypos;

                        xlow = -L / 2. + a * xpos;
                        ylow = -L / 2. + a * ypos;

                        // xhigh = -L/2.+a*xposUp;
                        // yhigh = -L/2.+a*yposUp;

                        fracx = (x - xlow) / a;

                        pos1 = xpos * N + ypos;
                        pos2 = xposUp * N + ypos;
                        pos3 = xpos * N + yposUp;
                        pos4 = xposUp * N + yposUp;

                        // -----------------------------epsilon---------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx)
                                     * abs(lat->cells[pos1]->getEpsilon())
                                 + fracx * abs(lat->cells[pos2]->getEpsilon());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx)
                                     * abs(lat->cells[pos3]->getEpsilon())
                                 + fracx * abs(lat->cells[pos4]->getEpsilon());
                        else
                            x2 = 0.;

                        fracy = (y - ylow) / a;

                        resultE = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------g2mu2A----------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 =
                                (1 - fracx) * abs(lat->cells[pos1]->getg2mu2A())
                                + fracx * abs(lat->cells[pos2]->getg2mu2A());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 =
                                (1 - fracx) * abs(lat->cells[pos3]->getg2mu2A())
                                + fracx * abs(lat->cells[pos4]->getg2mu2A());
                        else
                            x2 = 0.;

                        g2mu2A = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------g2mu2B----------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 =
                                (1 - fracx) * abs(lat->cells[pos1]->getg2mu2B())
                                + fracx * abs(lat->cells[pos2]->getg2mu2B());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 =
                                (1 - fracx) * abs(lat->cells[pos3]->getg2mu2B())
                                + fracx * abs(lat->cells[pos4]->getg2mu2B());
                        else
                            x2 = 0.;

                        g2mu2B = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------utau------------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getutau())
                                 + fracx * (lat->cells[pos2]->getutau());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getutau())
                                 + fracx * (lat->cells[pos4]->getutau());
                        else
                            x2 = 0.;

                        resultutau = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------ux--------------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getux())
                                 + fracx * (lat->cells[pos2]->getux());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getux())
                                 + fracx * (lat->cells[pos4]->getux());
                        else
                            x2 = 0.;

                        resultux = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------uy--------------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getuy())
                                 + fracx * (lat->cells[pos2]->getuy());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getuy())
                                 + fracx * (lat->cells[pos4]->getuy());
                        else
                            x2 = 0.;

                        resultuy = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------ueta------------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getueta())
                                 + fracx * (lat->cells[pos2]->getueta());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getueta())
                                 + fracx * (lat->cells[pos4]->getueta());
                        else
                            x2 = 0.;

                        resultueta = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------pitautau--------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getpitautau())
                                 + fracx * (lat->cells[pos2]->getpitautau());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getpitautau())
                                 + fracx * (lat->cells[pos4]->getpitautau());
                        else
                            x2 = 0.;

                        resultpi00 = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------pitaux----------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getpitaux())
                                 + fracx * (lat->cells[pos2]->getpitaux());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getpitaux())
                                 + fracx * (lat->cells[pos4]->getpitaux());
                        else
                            x2 = 0.;

                        resultpi0x = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------pitauy----------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getpitauy())
                                 + fracx * (lat->cells[pos2]->getpitauy());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getpitauy())
                                 + fracx * (lat->cells[pos4]->getpitauy());
                        else
                            x2 = 0.;

                        resultpi0y = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------pitaueta--------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getpitaueta())
                                 + fracx * (lat->cells[pos2]->getpitaueta());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getpitaueta())
                                 + fracx * (lat->cells[pos4]->getpitaueta());
                        else
                            x2 = 0.;

                        resultpi0eta = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------pixy----------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getpixy())
                                 + fracx * (lat->cells[pos2]->getpixy());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getpixy())
                                 + fracx * (lat->cells[pos4]->getpixy());
                        else
                            x2 = 0.;

                        resultpixy = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------pixeta----------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getpixeta())
                                 + fracx * (lat->cells[pos2]->getpixeta());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getpixeta())
                                 + fracx * (lat->cells[pos4]->getpixeta());
                        else
                            x2 = 0.;

                        resultpixeta = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------piyeta----------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getpiyeta())
                                 + fracx * (lat->cells[pos2]->getpiyeta());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getpiyeta())
                                 + fracx * (lat->cells[pos4]->getpiyeta());
                        else
                            x2 = 0.;

                        resultpiyeta = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------pixx------------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getpixx())
                                 + fracx * (lat->cells[pos2]->getpixx());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getpixx())
                                 + fracx * (lat->cells[pos4]->getpixx());
                        else
                            x2 = 0.;

                        resultpixx = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------piyy------------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getpiyy())
                                 + fracx * (lat->cells[pos2]->getpiyy());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getpiyy())
                                 + fracx * (lat->cells[pos4]->getpiyy());
                        else
                            x2 = 0.;

                        resultpiyy = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------pietaeta--------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 = (1 - fracx) * (lat->cells[pos1]->getpietaeta())
                                 + fracx * (lat->cells[pos2]->getpietaeta());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 = (1 - fracx) * (lat->cells[pos3]->getpietaeta())
                                 + fracx * (lat->cells[pos4]->getpietaeta());
                        else
                            x2 = 0.;

                        resultpietaeta = (1. - fracy) * x1 + fracy * x2;

                        resultutau = sqrt(
                            1. + resultux * resultux + resultuy * resultuy
                            + tau0 * tau0 * resultueta * resultueta);

                        Etot += abs(hbarc * resultE * gfactor) * ha * ha * it
                                * dtau * a;
                        if (abs(hbarc * resultE * gfactor) > 0.0000000001) {
                            foutEps2 << -(heta - 1) / 2. * deta + deta * ieta
                                     << " " << x << " " << y << " "
                                     << abs(hbarc * resultE * gfactor) << " "
                                     << resultutau << " " << resultux << " "
                                     << resultuy << " " << resultueta << " "
                                     << resultpi00 * gfactor << " "
                                     << resultpi0x * gfactor << " "
                                     << resultpi0y * gfactor << " "
                                     << resultpi0eta * gfactor << " "
                                     << resultpixx * gfactor << " "
                                     << resultpixy * gfactor << " "
                                     << resultpixeta * gfactor << " "
                                     << resultpiyy * gfactor << " "
                                     << resultpiyeta * gfactor << " "
                                     << resultpietaeta * gfactor << '\n';
                        } else {
                            foutEps2 << -(heta - 1) / 2. * deta + deta * ieta
                                     << " " << x << " " << y << " " << 0. << " "
                                     << 1. << " " << 0. << " " << 0. << " "
                                     << 0. << " " << 0. << " " << 0. << " "
                                     << 0. << " " << 0. << " " << 0. << " "
                                     << 0. << " " << 0. << " " << 0. << " "
                                     << 0. << " " << 0. << '\n';
                        }
                    } else {
                        foutEps2 << -(heta - 1) / 2. * deta + deta * ieta << " "
                                 << x << " " << y << " " << 0. << " " << 1.
                                 << " " << 0. << " " << 0. << " " << 0. << " "
                                 << 0. << " " << 0. << " " << 0. << " " << 0.
                                 << " " << 0. << " " << 0. << " " << 0. << " "
                                 << 0. << " " << 0. << " " << 0. << '\n';
                    }
                }
            }
            foutEps2 << '\n';
        }

        closeBufferedTextOutput(foutEps2, outputFilename);
        cout << "Etot = " << Etot << " GeV " << endl;
    }
    //       foutEtot <<  Etot << endl;
    // foutEtot.close();

    if (static_cast<int>(param->getWriteOutputs() / 4) == 1) {
        double resultT00, resultT0x, resultT0y, resultT0eta, resultTxx,
            resultTxy;
        double resultTxeta, resultTyy, resultTyeta, resultTetaeta;
        const bool writeBinaryTmunu = binaryTmunuEnabled(param);
        stringstream strTmunu_name;
        strTmunu_name << "Tmunu-t" << it * dtau * a << "-"
                      << param->getEventId()
                      << (writeBinaryTmunu ? ".ipgt" : ".dat");
        IPG_PROFILE_SCOPE(
            writeBinaryTmunu ? "output.tmunu_binary" : "output.tmunu_text");
        const string outputFilename = strTmunu_name.str();
        ofstream foutEps1;
        vector<char> outputBuffer;
        vector<float> binaryRow;
        if (writeBinaryTmunu) {
            openTmunuBinaryOutput(
                foutEps1, outputFilename, hx, hy, heta, tau0, deta, ha,
                param->getEventId());
            binaryRow.resize(static_cast<std::size_t>(hx) * 10u);
        } else {
            outputBuffer.resize(kTextOutputBufferBytes);
            openBufferedTextOutput(foutEps1, outputBuffer, outputFilename);
            foutEps1 << "# dummy " << 1 << " etamax= " << heta
                     << " xmax= " << hx << " ymax= " << hy << " deta= " << deta
                     << " dx= " << ha << " dy= " << ha << '\n';
        }
        // loop over all positions
        for (int iy = 0; iy < hy; iy++) {
            for (int ix = 0; ix < hx; ix++) {
                x = -hL / 2. + ha * ix;
                y = -hL / 2. + ha * iy;
                if (abs(x) < L / 2. && abs(y) < L / 2.) {
                    xpos = static_cast<int>(
                        floor((x + L / 2.) / a + 0.0000000001));
                    ypos = static_cast<int>(
                        floor((y + L / 2.) / a + 0.0000000001));

                    if (xpos < N - 1) {
                        xposUp = xpos + 1;
                    } else {
                        xposUp = xpos;
                    }
                    if (ypos < N - 1) {
                        yposUp = ypos + 1;
                    } else {
                        yposUp = ypos;
                    }

                    xlow = -L / 2. + a * xpos;
                    ylow = -L / 2. + a * ypos;

                    pos1 = xpos * N + ypos;
                    pos2 = xposUp * N + ypos;
                    pos3 = xpos * N + yposUp;
                    pos4 = xposUp * N + yposUp;

                    fracx = (x - xlow) / a;
                    fracy = (y - ylow) / a;

                    // ---------------------T^tautau----------------------- //
                    if (pos1 >= 0 && pos1 < N * N && pos2 >= 0
                        && pos2 < N * N) {
                        x1 =
                            ((1 - fracx) * (lat->cells[pos1]->getTtautau())
                             + fracx * (lat->cells[pos2]->getTtautau()));
                    } else {
                        x1 = 0.;
                    }
                    if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                        && pos4 < N * N) {
                        x2 =
                            ((1 - fracx) * (lat->cells[pos3]->getTtautau())
                             + fracx * (lat->cells[pos4]->getTtautau()));
                    } else {
                        x2 = 0.;
                    }
                    resultT00 = (1. - fracy) * x1 + fracy * x2;

                    // ---------------------T^taux----------------------- //
                    if (pos1 >= 0 && pos1 < N * N && pos2 >= 0
                        && pos2 < N * N) {
                        x1 =
                            ((1 - fracx) * (lat->cells[pos1]->getTtaux())
                             + fracx * (lat->cells[pos2]->getTtaux()));
                    } else {
                        x1 = 0.;
                    }
                    if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                        && pos4 < N * N) {
                        x2 =
                            ((1 - fracx) * (lat->cells[pos3]->getTtaux())
                             + fracx * (lat->cells[pos4]->getTtaux()));
                    } else {
                        x2 = 0.;
                    }
                    resultT0x = (1. - fracy) * x1 + fracy * x2;

                    // ---------------------T^tauy----------------------- //
                    if (pos1 >= 0 && pos1 < N * N && pos2 >= 0
                        && pos2 < N * N) {
                        x1 =
                            ((1 - fracx) * (lat->cells[pos1]->getTtauy())
                             + fracx * (lat->cells[pos2]->getTtauy()));
                    } else {
                        x1 = 0.;
                    }
                    if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                        && pos4 < N * N) {
                        x2 =
                            ((1 - fracx) * (lat->cells[pos3]->getTtauy())
                             + fracx * (lat->cells[pos4]->getTtauy()));
                    } else {
                        x2 = 0.;
                    }
                    resultT0y = (1. - fracy) * x1 + fracy * x2;

                    // ---------------------T^taueta----------------------- //
                    if (pos1 >= 0 && pos1 < N * N && pos2 >= 0
                        && pos2 < N * N) {
                        x1 =
                            ((1 - fracx) * (lat->cells[pos1]->getTtaueta())
                             + fracx * (lat->cells[pos2]->getTtaueta()));
                    } else {
                        x1 = 0.;
                    }
                    if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                        && pos4 < N * N) {
                        x2 =
                            ((1 - fracx) * (lat->cells[pos3]->getTtaueta())
                             + fracx * (lat->cells[pos4]->getTtaueta()));
                    } else {
                        x2 = 0.;
                    }
                    resultT0eta = (1. - fracy) * x1 + fracy * x2;

                    // ---------------------T^xx----------------------- //
                    if (pos1 >= 0 && pos1 < N * N && pos2 >= 0
                        && pos2 < N * N) {
                        x1 =
                            ((1 - fracx) * (lat->cells[pos1]->getTxx())
                             + fracx * (lat->cells[pos2]->getTxx()));
                    } else {
                        x1 = 0.;
                    }
                    if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                        && pos4 < N * N) {
                        x2 =
                            ((1 - fracx) * (lat->cells[pos3]->getTxx())
                             + fracx * (lat->cells[pos4]->getTxx()));
                    } else {
                        x2 = 0.;
                    }
                    resultTxx = (1. - fracy) * x1 + fracy * x2;

                    // ---------------------T^xy----------------------- //
                    if (pos1 >= 0 && pos1 < N * N && pos2 >= 0
                        && pos2 < N * N) {
                        x1 =
                            ((1 - fracx) * (lat->cells[pos1]->getTxy())
                             + fracx * (lat->cells[pos2]->getTxy()));
                    } else {
                        x1 = 0.;
                    }
                    if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                        && pos4 < N * N) {
                        x2 =
                            ((1 - fracx) * (lat->cells[pos3]->getTxy())
                             + fracx * (lat->cells[pos4]->getTxy()));
                    } else {
                        x2 = 0.;
                    }
                    resultTxy = (1. - fracy) * x1 + fracy * x2;

                    // ---------------------T^xeta----------------------- //
                    if (pos1 >= 0 && pos1 < N * N && pos2 >= 0
                        && pos2 < N * N) {
                        x1 =
                            ((1 - fracx) * (lat->cells[pos1]->getTxeta())
                             + fracx * (lat->cells[pos2]->getTxeta()));
                    } else {
                        x1 = 0.;
                    }
                    if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                        && pos4 < N * N) {
                        x2 =
                            ((1 - fracx) * (lat->cells[pos3]->getTxeta())
                             + fracx * (lat->cells[pos4]->getTxeta()));
                    } else {
                        x2 = 0.;
                    }
                    resultTxeta = (1. - fracy) * x1 + fracy * x2;

                    // ---------------------T^yy----------------------- //
                    if (pos1 >= 0 && pos1 < N * N && pos2 >= 0
                        && pos2 < N * N) {
                        x1 =
                            ((1 - fracx) * (lat->cells[pos1]->getTyy())
                             + fracx * (lat->cells[pos2]->getTyy()));
                    } else {
                        x1 = 0.;
                    }
                    if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                        && pos4 < N * N) {
                        x2 =
                            ((1 - fracx) * (lat->cells[pos3]->getTyy())
                             + fracx * (lat->cells[pos4]->getTyy()));
                    } else {
                        x2 = 0.;
                    }
                    resultTyy = (1. - fracy) * x1 + fracy * x2;

                    // ---------------------T^yeta----------------------- //
                    if (pos1 >= 0 && pos1 < N * N && pos2 >= 0
                        && pos2 < N * N) {
                        x1 =
                            ((1 - fracx) * (lat->cells[pos1]->getTyeta())
                             + fracx * (lat->cells[pos2]->getTyeta()));
                    } else {
                        x1 = 0.;
                    }
                    if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                        && pos4 < N * N) {
                        x2 =
                            ((1 - fracx) * (lat->cells[pos3]->getTyeta())
                             + fracx * (lat->cells[pos4]->getTyeta()));
                    } else {
                        x2 = 0.;
                    }
                    resultTyeta = (1. - fracy) * x1 + fracy * x2;

                    // ---------------------T^etaeta----------------------- //
                    if (pos1 >= 0 && pos1 < N * N && pos2 >= 0
                        && pos2 < N * N) {
                        x1 =
                            ((1 - fracx) * (lat->cells[pos1]->getTetaeta())
                             + fracx * (lat->cells[pos2]->getTetaeta()));
                    } else {
                        x1 = 0.;
                    }
                    if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                        && pos4 < N * N) {
                        x2 =
                            ((1 - fracx) * (lat->cells[pos3]->getTetaeta())
                             + fracx * (lat->cells[pos4]->getTetaeta()));
                    } else {
                        x2 = 0.;
                    }
                    resultTetaeta = (1. - fracy) * x1 + fracy * x2;

                    double values[10];
                    if (resultT00 * gfactor * hbarc > small_eps) {
                        values[0] = resultT00 * gfactor * hbarc;
                        values[1] = resultTxx * gfactor * hbarc;
                        values[2] = resultTyy * gfactor * hbarc;
                        values[3] =
                            tau0 * tau0 * resultTetaeta * gfactor * hbarc;
                        values[4] = -resultT0x * gfactor * hbarc;
                        values[5] = -resultT0y * gfactor * hbarc;
                        values[6] = -tau0 * resultT0eta * gfactor * hbarc;
                        values[7] = -resultTxy * gfactor * hbarc;
                        values[8] = -tau0 * resultTyeta * gfactor * hbarc;
                        values[9] = -tau0 * resultTxeta * gfactor * hbarc;
                    } else {
                        values[0] = small_eps;
                        values[1] = small_eps / 2.;
                        values[2] = small_eps / 2.;
                        for (int component = 3; component < 10; ++component) {
                            values[component] = 0.0;
                        }
                    }
                    if (writeBinaryTmunu) {
                        const std::size_t offset =
                            static_cast<std::size_t>(ix) * 10u;
                        for (int component = 0; component < 10; ++component) {
                            binaryRow[offset + component] =
                                static_cast<float>(values[component]);
                        }
                    } else {
                        foutEps1 << ix << " " << iy;
                        for (int component = 0; component < 10; ++component) {
                            foutEps1 << " " << values[component];
                        }
                        foutEps1 << '\n';
                    }
                } else {
                    double values[10] = {small_eps,
                                         small_eps / 2.,
                                         small_eps / 2.,
                                         0.0,
                                         0.0,
                                         0.0,
                                         0.0,
                                         0.0,
                                         0.0,
                                         0.0};
                    if (writeBinaryTmunu) {
                        const std::size_t offset =
                            static_cast<std::size_t>(ix) * 10u;
                        for (int component = 0; component < 10; ++component) {
                            binaryRow[offset + component] =
                                static_cast<float>(values[component]);
                        }
                    } else {
                        foutEps1 << ix << " " << iy;
                        for (int component = 0; component < 10; ++component) {
                            foutEps1 << " " << values[component];
                        }
                        foutEps1 << '\n';
                    }
                }
            }
            if (writeBinaryTmunu) {
                foutEps1.write(
                    reinterpret_cast<const char *>(binaryRow.data()),
                    static_cast<std::streamsize>(
                        binaryRow.size() * sizeof(binaryRow[0])));
            } else {
                foutEps1 << '\n';
            }
        }
        closeBufferedTextOutput(foutEps1, outputFilename);
    }

    if (tmunuOnly) return;

    if (static_cast<int>((param->getWriteOutputs() % 4) / 2) == 1) {
        double Jaztot = 0.;
        // Jazma output:
        // compute sum first for normalization
        for (int ieta = 0; ieta < heta; ieta++)  // loop over all positions
        {
            for (int ix = 0; ix < hx; ix++)  // loop over all positions
            {
                for (int iy = 0; iy < hy; iy++) {
                    x = -hL / 2. + ha * ix;
                    y = -hL / 2. + ha * iy;

                    if (abs(x) < L / 2. && abs(y) < L / 2.) {
                        xpos = static_cast<int>(
                            floor((x + L / 2.) / a + 0.0000000001));
                        ypos = static_cast<int>(
                            floor((y + L / 2.) / a + 0.0000000001));

                        if (xpos < N - 1)
                            xposUp = xpos + 1;
                        else
                            xposUp = xpos;

                        if (ypos < N - 1)
                            yposUp = ypos + 1;
                        else
                            yposUp = ypos;

                        xlow = -L / 2. + a * xpos;
                        ylow = -L / 2. + a * ypos;

                        fracx = (x - xlow) / a;

                        pos1 = xpos * N + ypos;
                        pos2 = xposUp * N + ypos;
                        pos3 = xpos * N + yposUp;
                        pos4 = xposUp * N + yposUp;

                        // -----------------------------g2mu2A----------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 =
                                (1 - fracx) * abs(lat->cells[pos1]->getg2mu2A())
                                + fracx * abs(lat->cells[pos2]->getg2mu2A());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 =
                                (1 - fracx) * abs(lat->cells[pos3]->getg2mu2A())
                                + fracx * abs(lat->cells[pos4]->getg2mu2A());
                        else
                            x2 = 0.;

                        fracy = (y - ylow) / a;

                        g2mu2A = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------g2mu2B----------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 =
                                (1 - fracx) * abs(lat->cells[pos1]->getg2mu2B())
                                + fracx * abs(lat->cells[pos2]->getg2mu2B());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 =
                                (1 - fracx) * abs(lat->cells[pos3]->getg2mu2B())
                                + fracx * abs(lat->cells[pos4]->getg2mu2B());
                        else
                            x2 = 0.;

                        g2mu2B = (1. - fracy) * x1 + fracy * x2;

                        Jaztot += g2mu2A * g2mu2B * ha * ha * it * dtau
                                  * a;  // same units as in Etot above
                    }
                }
            }
        }

        stringstream strJaz_name;
        strJaz_name << "Jazma-Hydro-t" << it * dtau * a << "-"
                    << param->getEventId() << ".dat";
        std::string Jaz_name;
        Jaz_name = strJaz_name.str();

        // stringstream strtwo_name;
        // strtwo_name << "twopointfct-t" << it*dtau*a << "-" <<
        // param->getEventId()
        // << ".dat"; string two_name; two_name = strtwo_name.str();

        ofstream foutEps3(Jaz_name.c_str(), std::ios::out);

        foutEps3 << "# dummy " << 1 << " etamax= " << heta << " xmax= " << hx
                 << " ymax= " << hy << " deta= " << deta << " dx= " << ha
                 << " dy= " << ha << endl;

        // ofstream foutEps4(two_name.c_str(),ios::out);

        // foutEps4 << "# dummy " << 1 << " etamax= " << heta
        //          << " xmax= " << hx << " ymax= " << hy << " deta= " << deta
        //          << " dx= " << ha << " dy= " << ha << endl;

        for (int ieta = 0; ieta < heta; ieta++)  // loop over all positions
        {
            for (int ix = 0; ix < hx; ix++)  // loop over all positions
            {
                for (int iy = 0; iy < hy; iy++) {
                    x = -hL / 2. + ha * ix;
                    y = -hL / 2. + ha * iy;

                    if (abs(x) < L / 2. && abs(y) < L / 2.) {
                        xpos = static_cast<int>(
                            floor((x + L / 2.) / a + 0.0000000001));
                        ypos = static_cast<int>(
                            floor((y + L / 2.) / a + 0.0000000001));

                        if (xpos < N - 1)
                            xposUp = xpos + 1;
                        else
                            xposUp = xpos;

                        if (ypos < N - 1)
                            yposUp = ypos + 1;
                        else
                            yposUp = ypos;

                        xlow = -L / 2. + a * xpos;
                        ylow = -L / 2. + a * ypos;

                        fracx = (x - xlow) / a;

                        pos1 = xpos * N + ypos;
                        pos2 = xposUp * N + ypos;
                        pos3 = xpos * N + yposUp;
                        pos4 = xposUp * N + yposUp;

                        // -----------------------------g2mu2A----------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 =
                                (1 - fracx) * abs(lat->cells[pos1]->getg2mu2A())
                                + fracx * abs(lat->cells[pos2]->getg2mu2A());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 =
                                (1 - fracx) * abs(lat->cells[pos3]->getg2mu2A())
                                + fracx * abs(lat->cells[pos4]->getg2mu2A());
                        else
                            x2 = 0.;

                        fracy = (y - ylow) / a;

                        g2mu2A = (1. - fracy) * x1 + fracy * x2;

                        // -----------------------------g2mu2B----------------------------------
                        // //
                        if (pos1 >= 0 && pos1 < (N) * (N) && pos2 >= 0
                            && pos2 < (N) * (N))
                            x1 =
                                (1 - fracx) * abs(lat->cells[pos1]->getg2mu2B())
                                + fracx * abs(lat->cells[pos2]->getg2mu2B());
                        else
                            x1 = 0.;

                        if (pos3 >= 0 && pos3 < N * N && pos4 >= 0
                            && pos4 < N * N)
                            x2 =
                                (1 - fracx) * abs(lat->cells[pos3]->getg2mu2B())
                                + fracx * abs(lat->cells[pos4]->getg2mu2B());
                        else
                            x2 = 0.;

                        g2mu2B = (1. - fracy) * x1 + fracy * x2;

                        // QsAsqr=
                        // g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg();
                        // QsBsqr=
                        // g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*hbarc*hbarc*param->getg()*param->getg();

                        foutEps3 << -(heta - 1) / 2. * deta + deta * ieta << " "
                                 << x << " " << y << " "
                                 << g2mu2A * g2mu2B / Jaztot * Etot << " " << 1.
                                 << " " << 0. << " " << 0. << " " << 0. << " "
                                 << 0. << " " << 0. << " " << 0. << " " << 0.
                                 << " " << 0. << " " << 0. << " " << 0. << " "
                                 << 0. << " " << 0. << " " << 0. << endl;

                        // // write two point and one point functions in
                        // [1/fm^6] and [1/fm^4] from
                        // https://arxiv.org/pdf/1902.07168.pdf if (QsAsqr>0
                        // && QsBsqr>0)
                        //   {
                        //     foutEps4 << -(heta-1)/2.*deta+deta*ieta << " " <<
                        //     x << " " << y << " "
                        //              << 16.*PI/9.*QsAsqr*QsBsqr/hbarc/hbarc/hbarc/hbarc*(QsAsqr/hbarc/hbarc*log(QsBsqr/pow(param->getm(),2.))
                        //                                                                  +QsBsqr/hbarc/hbarc*log(QsAsqr/pow(param->getm(),2.)))
                        //              << " "
                        //              << 4./3.*QsAsqr*QsBsqr/hbarc/hbarc/hbarc/hbarc
                        //              << " " << sqrt(QsAsqr) << " " <<
                        //              sqrt(QsBsqr) << " "
                        //              << 16.*PI/9.*QsAsqr*QsBsqr/hbarc/hbarc/hbarc/hbarc*(QsAsqr/hbarc/hbarc*log(QsBsqr/pow(param->getm(),2.)+1.)
                        //                                                                  +QsBsqr/hbarc/hbarc*log(QsAsqr/pow(param->getm(),2.)+1.)) << endl;
                        //     //        cout << QsAsqr << " " << QsBsqr << " "
                        //     << " " << param->getm()*param->getm() << endl;
                        //   }
                        // else
                        //   {
                        //     foutEps4 << -(heta-1)/2.*deta+deta*ieta << " " <<
                        //     x << " " << y << " "
                        //              << 0. << " "
                        //              << 4./3.*QsAsqr*QsBsqr/hbarc/hbarc/hbarc/hbarc
                        //              <<
                        //       " " << sqrt(QsAsqr) << " " << sqrt(QsBsqr) << "
                        //       " << 0. << endl;
                        //   }
                    } else {
                        foutEps3 << -(heta - 1) / 2. * deta + deta * ieta << " "
                                 << x << " " << y << " " << 0. << " " << 1.
                                 << " " << 0. << " " << 0. << " " << 0. << " "
                                 << 0. << " " << 0. << " " << 0. << " " << 0.
                                 << " " << 0. << " " << 0. << " " << 0. << " "
                                 << 0. << " " << 0. << " " << 0. << endl;
                        // foutEps4 << -(heta-1)/2.*deta+deta*ieta << " " << x
                        // << " " << y
                        // << " "
                        //          << 0. << " " << 0. << " "  << 0. << " " <<
                        //          0. << " " << 0. << endl;
                    }
                }
            }
        }
    }
    cout << "Wrote outputs" << endl;
    // done output for hydro
}

void MyEigen::test() {
    gsl_complex square;
    gsl_complex factor;
    gsl_complex euklidiansquare;
    GSL_SET_COMPLEX(&square, 0, 0);

    double data[] = {0.1,     -0.001, -0.01,   0.001, -0.1,
                     -0.0001, 0.01,   -0.0001, -0.1};

    gsl_matrix_view m = gsl_matrix_view_array(data, 3, 3);  // matrix

    gsl_vector_complex *eval = gsl_vector_complex_alloc(
        3);  // eigenvalues are components of this vector
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(
        3, 3);  // eigenvectors are columns of this matrix

    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(3);  // workspace

    gsl_eigen_nonsymmv(
        &m.matrix, eval, evec,
        w);  // solve for eigenvalues and eigenvectors (without 'v'
             // only compute eigenvalues)

    gsl_eigen_nonsymmv_free(w);  // free memory associated with workspace

    {  // output:
        int i;

        for (i = 0; i < 3; i++) {
            gsl_complex eval_i = gsl_vector_complex_get(eval, i);
            gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);

            printf(
                "before eigenvalue = %g + %gi\n", GSL_REAL(eval_i),
                GSL_IMAG(eval_i));
            GSL_SET_COMPLEX(&square, 0, 0);
            GSL_SET_COMPLEX(&euklidiansquare, 0, 0);

            for (int j = 0; j < 3; ++j) {
                gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
                printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
                euklidiansquare =
                    gsl_complex_add(euklidiansquare, gsl_complex_mul(z, z));

                if (j == 0)
                    square = gsl_complex_add(square, gsl_complex_mul(z, z));
                else
                    square = gsl_complex_sub(square, gsl_complex_mul(z, z));
            }

            if (GSL_REAL(square) > 0)
                printf(
                    "before u_mu u^mu = %g + %gi\n", GSL_REAL(square),
                    GSL_IMAG(square));
            else
                printf(
                    "before z_mu z^mu = %g + %gi\n", GSL_REAL(square),
                    GSL_IMAG(square));

            GSL_SET_COMPLEX(
                &factor,
                sqrt(abs(GSL_REAL(euklidiansquare) / GSL_REAL(square))), 0);
            printf(
                "eigenvalue = %g + %gi\n", GSL_REAL(eval_i), GSL_IMAG(eval_i));

            GSL_SET_COMPLEX(&square, 0, 0);

            for (int j = 0; j < 3; ++j) {
                gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
                z = gsl_complex_mul(z, factor);

                printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));

                if (j == 0)
                    square = gsl_complex_add(square, gsl_complex_mul(z, z));
                else
                    square = gsl_complex_sub(square, gsl_complex_mul(z, z));
            }

            if (GSL_REAL(square) > 0)
                printf(
                    "u_mu u^mu = %g + %gi\n", GSL_REAL(square),
                    GSL_IMAG(square));
            else
                printf(
                    "z_mu z^mu = %g + %gi\n", GSL_REAL(square),
                    GSL_IMAG(square));
        }
    }

    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);

    exit(1);
}

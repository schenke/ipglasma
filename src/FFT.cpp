// FFT.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.
// This version uses FFTW
#include "FFT.h"

#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>

#include "Instrumentation.h"
#include "Matrix.h"

//**************************************************************************
// FFT class.

//**************************************************************************

void FFT::fftnVector(
    vector<complex<double>> **data, vector<complex<double>> **outdata,
    const int nn[], const int isign) {
    IPG_PROFILE_SCOPE("fft.total");
    unsigned ntot = nn[0] * nn[1];
    int pos, newpos;
    vector<complex<double>>::iterator position;

    int mDim = data[0]->size();
    // mDim is the size of the vector (how many rows)

    // ndim is the dimension of the FFT (always 2 here)

    //    cout << "mDim=" << mDim << endl;
    // for each component of the vector fill the input array for the FFT (resort
    // as you fill in)
    for (int k = 0; k < mDim; k++) {
        //	oo   ->  xo
        //      ox       oo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j - nn[1] / 2;
                input[newpos][0] = real(data[pos]->at(k));
                input[newpos][1] = imag(data[pos]->at(k));
            }
        }
        //	xo   ->  oo
        //      oo       ox
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + nn[1] / 2 + j;
                input[newpos][0] = real(data[pos]->at(k));
                input[newpos][1] = imag(data[pos]->at(k));
            }
        }
        //	ox   ->  oo
        //      oo       xo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j + nn[1] / 2;
                input[newpos][0] = real(data[pos]->at(k));
                input[newpos][1] = imag(data[pos]->at(k));
            }
        }
        //	oo   ->  ox
        //      xo       oo
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + j - nn[1] / 2;
                input[newpos][0] = real(data[pos]->at(k));
                input[newpos][1] = imag(data[pos]->at(k));
            }
        }

        if (isign == 1)
            fftw_execute(p);
        else
            fftw_execute(pback);

        // if this is inverse transform, normalize.
        if (isign == -1) {
            for (unsigned i = 0; i < ntot; i++) {
                output[i][0] /= static_cast<double>(ntot);
                output[i][1] /= static_cast<double>(ntot);
            }
        }

        //	xo   ->  oo
        //      oo       ox
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j - nn[1] / 2;
                position = outdata[pos]->begin() + k;
                *position =
                    complex<double>(output[newpos][0], output[newpos][1]);
                //	cout << output[newpos][0] << " " << output[newpos][1] <<
                // endl;
                // cout << "data=" << data[pos]->at(k) << endl;
            }
        }
        //	oo   ->  xo
        //      ox       oo
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + nn[1] / 2 + j;
                position = outdata[pos]->begin() + k;
                *position =
                    complex<double>(output[newpos][0], output[newpos][1]);
            }
        }
        //	oo   ->  ox
        //      xo       oo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j + nn[1] / 2;
                position = outdata[pos]->begin() + k;
                *position =
                    complex<double>(output[newpos][0], output[newpos][1]);
            }
        }
        //	ox   ->  oo
        //      oo       xo

        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + j - nn[1] / 2;
                position = outdata[pos]->begin() + k;
                *position =
                    complex<double>(output[newpos][0], output[newpos][1]);
            }
        }
    }
}

void FFT::fftnArray(
    complex<double> **data, complex<double> **outdata, const int nn[],
    const int isign, const int mDim) {
    IPG_PROFILE_SCOPE("fft.total");
    unsigned ntot = nn[0] * nn[1];
    int pos, newpos;

    //    cout << "size=" << mDim << endl;
    // mDim is the size of the vector (how many rows)

    // ndim is the dimension of the FFT (always 2 here)

    //    cout << "mDim=" << mDim << endl;
    // for each component of the vector fill the input array for the FFT (resort
    // as you fill in)
    for (int k = 0; k < mDim; k++) {
        //	oo   ->  xo
        //      ox       oo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j - nn[1] / 2;
                input[newpos][0] = real(data[pos][k]);
                input[newpos][1] = imag(data[pos][k]);
            }
        }
        //	xo   ->  oo
        //      oo       ox
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + nn[1] / 2 + j;
                input[newpos][0] = real(data[pos][k]);
                input[newpos][1] = imag(data[pos][k]);
            }
        }
        //	ox   ->  oo
        //      oo       xo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j + nn[1] / 2;
                input[newpos][0] = real(data[pos][k]);
                input[newpos][1] = imag(data[pos][k]);
            }
        }
        //	oo   ->  ox
        //      xo       oo
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + j - nn[1] / 2;
                input[newpos][0] = real(data[pos][k]);
                input[newpos][1] = imag(data[pos][k]);
            }
        }

        if (isign == 1)
            fftw_execute(p);
        else
            fftw_execute(pback);

        // if this is inverse transform, normalize.
        if (isign == -1) {
            for (unsigned i = 0; i < ntot; i++) {
                output[i][0] /= static_cast<double>(ntot);
                output[i][1] /= static_cast<double>(ntot);
            }
        }

        //	xo   ->  oo
        //      oo       ox
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j - nn[1] / 2;
                outdata[pos][k] =
                    complex<double>(output[newpos][0], output[newpos][1]);
                //	cout << output[newpos][0] << " " << output[newpos][1] <<
                // endl;
                // cout << "data=" << data[pos][k] << endl;
            }
        }
        //	oo   ->  xo
        //      ox       oo
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + nn[1] / 2 + j;
                outdata[pos][k] =
                    complex<double>(output[newpos][0], output[newpos][1]);
            }
        }
        //	oo   ->  ox
        //      xo       oo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j + nn[1] / 2;
                outdata[pos][k] =
                    complex<double>(output[newpos][0], output[newpos][1]);
            }
        }
        //	ox   ->  oo
        //      oo       xo

        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + j - nn[1] / 2;
                outdata[pos][k] =
                    complex<double>(output[newpos][0], output[newpos][1]);
            }
        }
    }
}

// Performs Fast Fourier Transform of any object of class "T" (matrix or
// something else) using a wrapper for FFTW This routine takes data as a
// function of -x_max/2 to x_max/2 and returns it ordered similarly - no need to
// resort before or after!
template <class T>
void FFT::fftn(T **data, T **outdata, const int nn[], const int isign) {
    IPG_PROFILE_SCOPE("fft.total");
    const unsigned ntot = static_cast<unsigned>(nn[0] * nn[1]);
    int mDim = data[0]->getNDim();
    mDim *= mDim;

    // The historical centered FFT explicitly quadrant-shifted both the input
    // and output arrays.  For even lattice dimensions,
    //
    //     S F S = G M F M,
    //
    // where S shifts by half a lattice, M(i,j)=(-1)^(i+j), and
    // G=(-1)^(N0/2+N1/2).  Applying these checkerboard signs lets us keep the
    // FFTW buffers in natural row-major order, eliminating the four-quadrant
    // index remapping in both the pack and unpack passes.
    const double outputGlobalSign =
        (((nn[0] / 2 + nn[1] / 2) & 1) != 0) ? -1.0 : 1.0;
    const double inverseNorm =
        (isign == -1) ? 1.0 / static_cast<double>(ntot) : 1.0;

    const bool profileMatrixFFT = ipg::Profiler::instance().enabled();
    double packSeconds = 0.;
    double executeSeconds = 0.;
    double unpackSeconds = 0.;

    // All nine SU(3) matrix components are independent FFTs.  Pack each
    // component into its own plane of the existing 9x scratch buffers, then
    // execute the same FFTW plan concurrently on different planes using the
    // new-array interface.  FFTW guarantees concurrent execution of the same
    // plan is thread-safe when the input/output arrays are distinct.
    double segmentStart = 0.;
#pragma omp parallel shared( \
        segmentStart, packSeconds, executeSeconds, unpackSeconds)
    {
#pragma omp single
        {
            if (profileMatrixFFT) segmentStart = ipg::wallSeconds();
        }

#pragma omp for schedule(static)
        for (int k = 0; k < mDim; ++k) {
            fftw_complex *localInput =
                inputMany + static_cast<std::size_t>(k) * ntot;
            for (int i = 0; i < nn[0]; ++i) {
                for (int j = 0; j < nn[1]; ++j) {
                    const int pos = i * nn[1] + j;
                    const double sign = ((i + j) & 1) ? -1.0 : 1.0;
                    const complex<double> value = data[pos]->data()[k];
                    localInput[pos][0] = sign * value.real();
                    localInput[pos][1] = sign * value.imag();
                }
            }
        }

#pragma omp single
        {
            if (profileMatrixFFT) {
                packSeconds = ipg::wallSeconds() - segmentStart;
                segmentStart = ipg::wallSeconds();
            }
        }

#pragma omp for schedule(static)
        for (int k = 0; k < mDim; ++k) {
            fftw_complex *localInput =
                inputMany + static_cast<std::size_t>(k) * ntot;
            fftw_complex *localOutput =
                outputMany + static_cast<std::size_t>(k) * ntot;
            if (isign == 1)
                fftw_execute_dft(p, localInput, localOutput);
            else
                fftw_execute_dft(pback, localInput, localOutput);
        }

#pragma omp single
        {
            if (profileMatrixFFT) {
                executeSeconds = ipg::wallSeconds() - segmentStart;
                segmentStart = ipg::wallSeconds();
            }
        }

#pragma omp for schedule(static)
        for (int k = 0; k < mDim; ++k) {
            fftw_complex *localOutput =
                outputMany + static_cast<std::size_t>(k) * ntot;
            for (int i = 0; i < nn[0]; ++i) {
                for (int j = 0; j < nn[1]; ++j) {
                    const int pos = i * nn[1] + j;
                    const double sign =
                        ((i + j) & 1) ? -outputGlobalSign : outputGlobalSign;
                    outdata[pos]->data()[k] = complex<double>(
                        sign * localOutput[pos][0] * inverseNorm,
                        sign * localOutput[pos][1] * inverseNorm);
                }
            }
        }

#pragma omp single
        {
            if (profileMatrixFFT) {
                unpackSeconds = ipg::wallSeconds() - segmentStart;
            }
        }
    }

    if (profileMatrixFFT) {
        ipg::Profiler &profiler = ipg::Profiler::instance();
        profiler.add("fft.matrix.pack", packSeconds);
        profiler.add("fft.matrix.execute", executeSeconds);
        profiler.add("fft.matrix.unpack", unpackSeconds);
    }
}

template <class T>
void FFT::fftnMany(T **data, T **outdata, const int nn[], const int isign) {
    IPG_PROFILE_SCOPE("fft.total");
    unsigned ntot = nn[0] * nn[1];
    int mDim, pos, newpos;  // matrix dimension
    mDim = data[0]->getNDim();

    // mDim is the size of the matrix (how many rows)

    mDim *= mDim;
    // ndim is the dimension of the FFT (always 2 here)

    // for each component of the matrix fill the input array for the FFT (resort
    // as you fill in)

    //    cout << "mDim=" << mDim << endl;

    for (int k = 0; k < mDim; k++) {
        //	oo   ->  xo
        //      ox       oo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j - nn[1] / 2 + k * ntot;
                input[newpos][0] = data[pos]->getRe(k);
                input[newpos][1] = data[pos]->getIm(k);
            }
        }
        //	xo   ->  oo
        //      oo       ox
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j + k;
                newpos = (i + nn[0] / 2) * nn[1] + nn[1] / 2 + j + k * ntot;
                input[newpos][0] = data[pos]->getRe(k);
                input[newpos][1] = data[pos]->getIm(k);
            }
        }
        //	ox   ->  oo
        //      oo       xo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j + nn[1] / 2 + k * ntot;
                input[newpos][0] = data[pos]->getRe(k);
                input[newpos][1] = data[pos]->getIm(k);
            }
        }
        //	oo   ->  ox
        //      xo       oo
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + j - nn[1] / 2 + k * ntot;
                input[newpos][0] = data[pos]->getRe(k);
                input[newpos][1] = data[pos]->getIm(k);
            }
        }
    }

    if (isign == 1)
        fftw_execute(pmany);
    else
        fftw_execute(pmanyback);

    // if this is inverse transform, normalize.
    if (isign == -1) {
        for (unsigned i = 0; i < ntot * 9; i++) {
            output[i][0] /= static_cast<double>(ntot);
            output[i][1] /= static_cast<double>(ntot);
        }
    }

    for (int k = 0; k < mDim; k++) {
        //	oo   ->  xo
        //      ox       oo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j - nn[1] / 2 + k * ntot;
                outdata[pos]->setRe(k, output[newpos][0]);
                outdata[pos]->setIm(k, output[newpos][1]);
            }
        }
        //	xo   ->  oo
        //      oo       ox
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + nn[1] / 2 + j + k * ntot;
                outdata[pos]->setRe(k, output[newpos][0]);
                outdata[pos]->setIm(k, output[newpos][1]);
            }
        }
        //	ox   ->  oo
        //      oo       xo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j + nn[1] / 2 + k * ntot;
                outdata[pos]->setRe(k, output[newpos][0]);
                outdata[pos]->setIm(k, output[newpos][1]);
            }
        }
        //	oo   ->  ox
        //      xo       oo

        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + j - nn[1] / 2 + k * ntot;
                outdata[pos]->setRe(k, output[newpos][0]);
                outdata[pos]->setIm(k, output[newpos][1]);
            }
        }
    }

    //------
}

void FFT::fftnComplex(
    complex<double> *data, complex<double> *outdata, const int nn[],
    const int isign) {
    IPG_PROFILE_SCOPE("fft.total");
    unsigned ntot = nn[0] * nn[1];
    int mDim, pos, newpos;  // matrix dimension
    mDim = 1;

    // mDim is the size of the matrix (how many rows)

    mDim *= mDim;
    // ndim is the dimension of the FFT (always 2 here)

    // for each component of the matrix fill the input array for the FFT (resort
    // as you fill in)

    //    cout << "mDim=" << mDim << endl;
    for (int k = 0; k < mDim; k++) {
        //	oo   ->  xo
        //      ox       oo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j - nn[1] / 2;
                input[newpos][0] = data[pos].real();
                input[newpos][1] = data[pos].imag();
            }
        }
        //	xo   ->  oo
        //      oo       ox
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + nn[1] / 2 + j;
                input[newpos][0] = data[pos].real();
                input[newpos][1] = data[pos].imag();
            }
        }
        //	ox   ->  oo
        //      oo       xo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j + nn[1] / 2;
                input[newpos][0] = data[pos].real();
                input[newpos][1] = data[pos].imag();
            }
        }
        //	oo   ->  ox
        //      xo       oo
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + j - nn[1] / 2;
                input[newpos][0] = data[pos].real();
                input[newpos][1] = data[pos].imag();
            }
        }

        if (isign == 1)
            fftw_execute(p);
        else
            fftw_execute(pback);

        // if this is inverse transform, normalize.
        if (isign == -1) {
            for (unsigned i = 0; i < ntot; i++) {
                output[i][0] /= static_cast<double>(ntot);
                output[i][1] /= static_cast<double>(ntot);
            }
        }

        //	oo   ->  xo
        //      ox       oo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j - nn[1] / 2;
                outdata[pos] =
                    complex<double>(output[newpos][0], output[newpos][1]);
                //	outdata[pos]->setIm(k,output[newpos][1]);
            }
        }
        //	xo   ->  oo
        //      oo       ox
        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + nn[1] / 2 + j;
                outdata[pos] =
                    complex<double>(output[newpos][0], output[newpos][1]);
                //	outdata[pos]->setIm(k,output[newpos][1]);
            }
        }
        //	ox   ->  oo
        //      oo       xo
        for (int i = nn[0] / 2; i < nn[0]; i++) {
            for (int j = 0; j < nn[1] / 2; j++) {
                pos = i * nn[1] + j;
                newpos = (i - nn[0] / 2) * nn[1] + j + nn[1] / 2;
                outdata[pos] =
                    complex<double>(output[newpos][0], output[newpos][1]);
                //	outdata[pos]->setIm(k,output[newpos][1]);
            }
        }
        //	oo   ->  ox
        //      xo       oo

        for (int i = 0; i < nn[0] / 2; i++) {
            for (int j = nn[1] / 2; j < nn[1]; j++) {
                pos = i * nn[1] + j;
                newpos = (i + nn[0] / 2) * nn[1] + j - nn[1] / 2;
                outdata[pos] =
                    complex<double>(output[newpos][0], output[newpos][1]);
                //	outdata[pos]->setIm(k,output[newpos][1]);
            }
        }
    }
}

// Define specializations of the template:
template void FFT::fftn(
    Matrix **data, Matrix **outdata, const int nn[], const int isign);
template void FFT::fftnMany(
    Matrix **data, Matrix **outdata, const int nn[], const int isign);

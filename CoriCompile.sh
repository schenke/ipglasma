#!/bin/bash

module swap PrgEnv-intel PrgEnv-gnu
module load gsl
module load cray-fftw

rm *exe
cd src
rm *exe
CC -std=c++11 NIKLAS_MAIN.cpp -o GENERATE_IPGLASMA_DIJET.exe Fragmentation.cpp FFT.cpp Matrix.cpp Setup.cpp Init.cpp Random.cpp Group.cpp Lattice.cpp Cell.cpp Glauber.cpp Util.cpp Evolution.cpp GaugeFix.cpp Spinor.cpp MyEigen.cpp -I/usr/local/include -L/usr/local/include -L/usr/local/lib `gsl-config --libs` `gsl-config --cflags`  -lfftw3 -lfftw3_threads -lm -lgsl -lgslcblas -fopenmp
cp GENERATE_IPGLASMA_DIJET.exe ..
cd ..



#cat NiklasCoriMakefile>Makefile
#make clean
#make install
#rm GENERATE_IPGLASMA_DIJET.exe
#cd ..

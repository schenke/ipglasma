cd src
mpicxx -I/opt/local/include/openmpi-gcc5 -L/opt/local/lib/openmpi-gcc5 -lmpi_cxx -std=c++11 -O3 NIKLAS_MAIN.cpp -o GENERATE_IPGLASMA_DIJET.exe Fragmentation.cpp FFT.cpp Matrix.cpp Setup.cpp Init.cpp Random.cpp Group.cpp Lattice.cpp Cell.cpp Glauber.cpp Util.cpp Evolution.cpp GaugeFix.cpp Spinor.cpp MyEigen.cpp -I/usr/local/include -L/usr/local/include -L/usr/local/lib `gsl-config --libs` `gsl-config --cflags`  -lfftw3 -lfftw3_threads -lm -lgsl -lgslcblas -fopenmp -finline-functions -I/usr/local/include
cp GENERATE_IPGLASMA_DIJET.exe ..
rm GENERATE_IPGLASMA_DIJET.exe
cd ..

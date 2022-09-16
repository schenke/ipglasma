#/usr/bin/env bash

Flag=$1

mkdir -p build
cd build
rm -fr *
if [ "$Flag" == "KNL" ]; then
    CXX=mpiicpc cmake .. -DKNL=ON
elif [ "$Flag" == "noMPI" ]; then
    cmake .. -DdisableMPI=ON
else
    CXX=g++ cmake ..
fi

make -j4
make install
cd ..
rm -fr build

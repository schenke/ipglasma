#!/usr/bin/env bash
# Builds IP-Glasma with CMake in build/ and installs the executable ipglasma
# in the repository root.
# Usage: ./compile_IPGlasma.sh [noMPI|KNL]
set -euo pipefail

Flag=${1:-}
case "$Flag" in
    "" | noMPI | KNL) ;;
    *)
        echo "Unknown option '$Flag'. Usage: $0 [noMPI|KNL]" >&2
        exit 1
        ;;
esac

# work in the repository root, wherever the script is called from
cd "$(dirname "$0")"

mkdir -p build
cd build
rm -fr ./*
if [ "$Flag" == "KNL" ]; then
    CXX=mpiicpc cmake .. -DKNL=ON
elif [ "$Flag" == "noMPI" ]; then
    cmake .. -DdisableMPI=ON
else
    cmake ..
fi

make -j4
make install
cd ..
rm -fr build

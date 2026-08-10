#!/usr/bin/env bash
set -euo pipefail

# Build IP-Glasma with GNU GCC/G++ and OpenMP, then copy the executable
# into the source-tree root as ./ipglasma.
#
# Usage:
#   ./compile_gnu_openmp.sh
#   ./compile_gnu_openmp.sh --clean
#
# Optional environment overrides:
#   GCC=/opt/homebrew/bin/gcc-15 GXX=/opt/homebrew/bin/g++-15 ./compile_gnu_openmp.sh
#   BUILD_JOBS=8 ./compile_gnu_openmp.sh

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${ROOT}/build-gnu-openmp"
OUTPUT_EXE="${ROOT}/ipglasma"

clean=0
if [[ "${1:-}" == "--clean" ]]; then
    clean=1
elif [[ $# -gt 0 ]]; then
    echo "Usage: $0 [--clean]" >&2
    exit 2
fi

find_gnu_compilers() {
    # Respect explicit overrides first.
    if [[ -n "${GCC:-}" || -n "${GXX:-}" ]]; then
        if [[ -z "${GCC:-}" || -z "${GXX:-}" ]]; then
            echo "Error: set both GCC and GXX, or neither." >&2
            exit 1
        fi
        CC_BIN="$GCC"
        CXX_BIN="$GXX"
        return
    fi

    # Homebrew/macOS typically installs GNU compilers with version suffixes.
    for v in 16 15 14 13 12 11; do
        if command -v "gcc-${v}" >/dev/null 2>&1 &&
           command -v "g++-${v}" >/dev/null 2>&1; then
            CC_BIN="$(command -v "gcc-${v}")"
            CXX_BIN="$(command -v "g++-${v}")"
            return
        fi
    done

    # Fall back to unversioned gcc/g++ only if they really are GNU.
    if command -v gcc >/dev/null 2>&1 &&
       command -v g++ >/dev/null 2>&1 &&
       gcc --version 2>/dev/null | head -n 1 | grep -qi 'gcc'; then
        CC_BIN="$(command -v gcc)"
        CXX_BIN="$(command -v g++)"
        return
    fi

    echo "Error: GNU GCC/G++ not found." >&2
    echo "Install GCC with Homebrew, e.g.:" >&2
    echo "  brew install gcc" >&2
    exit 1
}

find_gnu_compilers

echo "Source tree : $ROOT"
echo "Build dir   : $BUILD_DIR"
echo "C compiler  : $CC_BIN"
echo "C++ compiler: $CXX_BIN"
echo
"$CC_BIN" --version | head -n 1
"$CXX_BIN" --version | head -n 1
echo

if [[ $clean -eq 1 ]]; then
    echo "Removing old build directory..."
    rm -rf "$BUILD_DIR"
fi


cmake -S "$ROOT" -B "$BUILD_DIR" \
    -DdisableMPI=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$CC_BIN" \
    -DCMAKE_CXX_COMPILER="$CXX_BIN" \
    -DCMAKE_CXX_FLAGS="-DIPGLASMA_DETERMINISTIC_FFT"

# CMake's OpenMP package/project logic should select -fopenmp with GNU GCC.
# Verify that OpenMP was detected/configured before building.
if ! grep -Eq 'OpenMP_CXX_FLAGS(:[^=]*)?=.*-fopenmp|OpenMP_CXX_FOUND(:[^=]*)?=TRUE' \
        "$BUILD_DIR/CMakeCache.txt" 2>/dev/null; then
    echo
    echo "Warning: could not confirm OpenMP from CMakeCache.txt."
    echo "The project may still configure OpenMP internally; inspect the verbose build if needed:"
    echo "  cmake --build \"$BUILD_DIR\" --verbose"
    echo
fi

JOBS="${BUILD_JOBS:-}"
if [[ -n "$JOBS" ]]; then
    cmake --build "$BUILD_DIR" -j "$JOBS"
else
    cmake --build "$BUILD_DIR" -j
fi

BUILT_EXE="$BUILD_DIR/src/ipglasma"

if [[ ! -x "$BUILT_EXE" ]]; then
    echo "Error: expected executable not found at:" >&2
    echo "  $BUILT_EXE" >&2
    exit 1
fi

cp -f "$BUILT_EXE" "$OUTPUT_EXE"
chmod +x "$OUTPUT_EXE"

echo
echo "Build complete."
echo "Executable copied to:"
echo "  $OUTPUT_EXE"
echo

# On macOS, show whether libgomp is linked.
if command -v otool >/dev/null 2>&1; then
    echo "OpenMP runtime linkage:"
    if otool -L "$OUTPUT_EXE" | grep -i 'gomp'; then
        :
    else
        echo "  Warning: libgomp was not found in otool output."
    fi
fi

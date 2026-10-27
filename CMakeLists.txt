cmake_minimum_required(VERSION 3.10)
project(IPGlasma VERSION 2.0.0 LANGUAGES CXX)

set(libname "ipglasma_lib")
set(exename "ipglasma")

# require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Tell cmake where to find modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/CMakeModules")

set(CMAKE_INSTALL_RPATH ${CMAKE_HOME_DIRECTORY})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Find GSL and FFTW for this project:
find_package(FFTW REQUIRED)
if (FFTW_FOUND)
    message("Found FFTW library ${FFTW_INCLUDE_DIRS}")
    include_directories(${FFTW_INCLUDE_DIRS})
endif()

find_package(GSL REQUIRED)
if (GSL_FOUND)
    message("Found GSL library ${GSL_INCLUDE_DIR}")
    include_directories(${GSL_INCLUDE_DIR})
endif()

# PrettyOstream uses std::mutex to serialize terminal output across
# threads; find it explicitly rather than relying on OpenMP to pull in
# pthread as a side effect, since OpenMP is unavailable/disabled on some
# platforms (e.g. AppleClang below).
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
        message(FATAL_ERROR "Insufficient gcc version")
    endif()
endif()

# OpenMP is linked as a target (OPENMP_LIBS), which adds its compile and
# link flags
set(OPENMP_LIBS "")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    message("Using Clang compiler without OpenMP parallelization... ")
else ()
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        set(OPENMP_LIBS OpenMP::OpenMP_CXX)
    endif ()
endif()

# default build type: Release (-O3 -DNDEBUG); Debug (-g) and RelWithDebInfo
# (-O2 -g -DNDEBUG) are the other common choices. The optimization level
# only comes from the build type.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_HOME_DIRECTORY}/"
        CACHE PATH "Install path prefix, prepended onto install directories."
        FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

option (KNL "Build executable on KNL" OFF)
option (unittest "Build Unit tests" OFF)
option (build_lib "Build library file" OFF)
option (disableMPI "Build executable without MPI" OFF)
if (NOT disableMPI)
    find_package(MPI)
    if (MPI_FOUND)
        message("Found MPI library ${MPI_CXX_INCLUDE_DIRS}")
        include_directories(${MPI_CXX_INCLUDE_PATH})
    else ()
        option (disableMPI "Build executable without MPI" ON)
    endif()
endif()

# Compiler-specific flags are appended, so flags given with -DCMAKE_CXX_FLAGS
# are kept. There is no -Ofast or -ffast-math: with them, std::isnan is
# always false, which disables the NaN checks of the code.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    if (KNL)
        string(APPEND CMAKE_CXX_FLAGS " -xMIC-AVX512 -march=knl -finline-functions")
    else ()
        string(APPEND CMAKE_CXX_FLAGS " -ipo -qopt-report=5 -qopt-report-phase:vec -fimf-use-svml -fma -align -finline-functions")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    string(APPEND CMAKE_CXX_FLAGS " -stdlib=libc++")
endif()

string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra -Woverloaded-virtual -Wredundant-decls -Wshadow -Wundef -Wno-unused")

# std::filesystem (used in Parameters.cpp) is part of libstdc++/libc++ proper
# starting with GCC 9 / Clang 9, but on older compilers that still ship it
# under C++17 it lives in a separate runtime library that must be linked
# explicitly.
set(FILESYSTEM_LIBS "")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    set(FILESYSTEM_LIBS stdc++fs)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    set(FILESYSTEM_LIBS c++fs)
endif()

# Pins FFTW plan selection via an on-disk wisdom cache (see FFT.h) so that
# repeated runs with the same seed are bit-reproducible, at the cost of the
# first run in a given directory measuring (rather than estimating) its FFT
# plan. Off by default. Enable with -DIPGLASMA_DETERMINISTIC_FFT=ON.
option(IPGLASMA_DETERMINISTIC_FFT
    "Pin FFTW plan selection via an on-disk wisdom cache for bit-reproducible runs"
    OFF)
if (IPGLASMA_DETERMINISTIC_FFT)
    # add_definitions rather than add_compile_definitions, which needs
    # CMake 3.12
    add_definitions(-DIPGLASMA_DETERMINISTIC_FFT)
endif()

add_subdirectory (src)

# Doctest-based unit test suite (tests/), covering small, fast, deterministic
# pieces (matrix algebra, RNG stream equivalence, input parsing, ...) rather
# than full lattice evolution. Off by default; enable with -Dunittest=ON.
# Vendors its own copy of doctest.h (tests/third_party/), so this needs no
# extra dependency beyond what the main build already requires.
if (unittest)
    enable_testing()
    add_subdirectory (tests)
endif()

# Doxygen developer-reference documentation (`make doc`); see
# doc/CMakeLists.txt. A no-op with a message if Doxygen isn't installed.
add_subdirectory (doc)

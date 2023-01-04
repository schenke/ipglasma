cmake_minimum_required(VERSION 3.4)
project(IPGlasma VERSION 1.0 LANGUAGES CXX)

set(libname "ipglasma_lib")
set(exename "ipglasma")

# require C++11
set(CMAKE_CXX_STANDARD 11)
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

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
        message(FATAL_ERROR "Insufficient gcc version")
    endif()
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    message("Using Clang compiler without OpenMP parallelization... ")
else ()
    find_package(OpenMP)
    if (OpenMP_FOUND)
        set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS})
        set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})
    endif ()
endif()

# default build type: Release
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

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    if (KNL)
       set(CMAKE_CXX_FLAGS "-g -Wall -xMIC-AVX512 -std=c++11 -O3 -fopenmp -finline-functions -march=knl  -DDOCTEST_CONFIG_DISABLE")
       #set(CMAKE_CXX_FLAGS "-g -vecabi=cmdtarget -debug inline-debug-info -fimf-use-svml -xmic-avx512 -qopt-report=5 -qopt-report-phase:vec -Ofast ${OpenMP_CXX_FLAGS} -std=c++11 -fma -align -finline-functions -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
    else ()
        set(CMAKE_CXX_FLAGS "-g -ipo -qopt-report=5 -qopt-report-phase:vec -Ofast ${OpenMP_CXX_FLAGS} -std=c++11 -fimf-use-svml -fma -align -finline-functions -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-g -Ofast ${OpenMP_CXX_FLAGS} -std=c++11 -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "-g -O3 -std=c++11  -stdlib=libc++ -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
endif()

if (unittest)
    set(CMAKE_CXX_FLAGS "-g ${OpenMP_CXX_FLAGS} -std=c++11 -DDOCTEST_CONFIG_IMPLEMENT_WITH_MAIN")
endif()

string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra -Woverloaded-virtual -Wredundant-decls -Wshadow -Wundef -Wno-unused")

add_subdirectory (src)

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

# The GPU backends of the device layer (src/device/). The layer itself and
# its host backend are part of IPGLASMA_LIB_SOURCES (src/Sources.cmake) and
# always built. This file adds, for the executable and the unit tests alike,
#   IPGLASMA_DEVICE_SOURCES      extra sources (absolute paths)
#   IPGLASMA_DEVICE_LIBS         libraries to link
#   IPGLASMA_DEVICE_DEFINITIONS  compile definitions
#   IPGLASMA_DEVICE_INCLUDES     include directories
# for
#   -DIPGLASMA_CUDA=ON    NVIDIA GPUs (CUDA toolkit with cuFFT)
#   -DIPGLASMA_OPENCL=ON  OpenCL 1.2 devices, e.g. the GPU of an Apple M
#                         processor (macOS ships OpenCL; nothing to install)
# Both can be on; IPGLASMA_DEVICE selects the backend at run time.

option(IPGLASMA_CUDA "Build the CUDA backend of the device layer" OFF)
option(IPGLASMA_OPENCL
    "Build the OpenCL backend of the device layer (Apple M GPUs and others)"
    OFF)

set(IPGLASMA_DEVICE_SOURCES "")
set(IPGLASMA_DEVICE_LIBS "")
set(IPGLASMA_DEVICE_DEFINITIONS "")
set(IPGLASMA_DEVICE_INCLUDES "")
set(ipglasmaDeviceDir "${CMAKE_CURRENT_LIST_DIR}")

if (IPGLASMA_CUDA)
    if (CMAKE_VERSION VERSION_LESS 3.17)
        message(FATAL_ERROR "IPGLASMA_CUDA needs CMake 3.17 or newer")
    endif()
    enable_language(CUDA)
    # let -DCUDAToolkit_ROOT=... (or the environment variable) point to the
    # toolkit
    if (POLICY CMP0074)
        cmake_policy(SET CMP0074 NEW)
    endif()
    find_package(CUDAToolkit REQUIRED)
    if (NOT DEFINED CMAKE_CUDA_STANDARD)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    endif()
    list(APPEND IPGLASMA_DEVICE_SOURCES "${ipglasmaDeviceDir}/CudaBackend.cu")
    list(APPEND IPGLASMA_DEVICE_LIBS CUDA::cufft CUDA::cudart)
    list(APPEND IPGLASMA_DEVICE_DEFINITIONS IPGLASMA_CUDA)
    message(STATUS "Device layer: CUDA backend (CUDA ${CUDAToolkit_VERSION})")
endif()

if (IPGLASMA_OPENCL)
    find_package(OpenCL REQUIRED)
    # the OpenCL driver compiles the kernels at start-up from their source:
    # ipg_kernels.h as a C++ raw string literal, regenerated when it changes
    set(kernelSourceInc "${CMAKE_BINARY_DIR}/device/ipg_kernels_source.inc")
    file(READ "${ipglasmaDeviceDir}/ipg_kernels.h" kernelSource)
    file(WRITE "${kernelSourceInc}.new"
        "R\"IPGKERNELS(\n${kernelSource})IPGKERNELS\"\n")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${kernelSourceInc}.new" "${kernelSourceInc}")
    file(REMOVE "${kernelSourceInc}.new")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
        "${ipglasmaDeviceDir}/ipg_kernels.h")
    list(APPEND IPGLASMA_DEVICE_SOURCES
        "${ipglasmaDeviceDir}/OpenCLBackend.cpp")
    list(APPEND IPGLASMA_DEVICE_LIBS OpenCL::OpenCL)
    list(APPEND IPGLASMA_DEVICE_DEFINITIONS IPGLASMA_OPENCL)
    list(APPEND IPGLASMA_DEVICE_INCLUDES "${CMAKE_BINARY_DIR}/device")
    message(STATUS "Device layer: OpenCL backend (${OpenCL_LIBRARIES})")
endif()

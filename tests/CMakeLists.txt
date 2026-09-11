# Unit test executable: links the same IP-Glasma library sources as the
# main "ipglasma" target (see src/Sources.cmake), minus main.cpp, plus the
# doctest test files below. None of the library sources have a hard MPI
# dependency (Init.cpp's MPI calls are commented out), so this target never
# links MPI, regardless of the disableMPI setting used for the main build.
include(${CMAKE_SOURCE_DIR}/src/Sources.cmake)

set(LIB_SOURCES ${IPGLASMA_LIB_SOURCES})
list(TRANSFORM LIB_SOURCES PREPEND "${CMAKE_SOURCE_DIR}/src/")

set(TEST_SOURCES
    test_main.cpp
    test_matrix.cpp
    test_su3.cpp
    test_group.cpp
    test_random.cpp
    test_setup.cpp
    test_glauber.cpp
    test_parameters.cpp
    test_fft.cpp
    test_lattice.cpp
    test_cell.cpp
    test_util.cpp
    test_prettyostream.cpp
    test_instrumentation.cpp
    test_fragmentation.cpp
    )

add_executable(ipglasma_tests ${TEST_SOURCES} ${LIB_SOURCES})
target_include_directories(ipglasma_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party
    )
target_link_libraries(ipglasma_tests ${GSL_LIBRARIES} ${FFTW_LIBRARIES} Threads::Threads)
if (APPLE)
    # Matches src/CMakeLists.txt: PrettyOstream::getMemoryUsage() reads
    # ru_maxrss in different units on macOS vs. Linux and picks the
    # conversion based on this define.
    target_compile_definitions(ipglasma_tests PRIVATE APPLE)
endif (APPLE)

add_test(NAME ipglasma_unit_tests COMMAND ipglasma_tests)

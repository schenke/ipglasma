include(${CMAKE_CURRENT_SOURCE_DIR}/Sources.cmake)
set (SOURCES ${IPGLASMA_LIB_SOURCES})

set(CompileFlags "")
if (APPLE)
    set(CompileFlags "${CompileFlags} -DAPPLE")
endif (APPLE)
if (disableMPI)
    set(CompileFlags "${CompileFlags} -DDISABLEMPI")
endif (disableMPI)

if (build_lib)
    add_library(${libname} SHARED ${SOURCES})
    set_target_properties (${libname} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
    target_link_libraries (${libname} ${GSL_LIBRARIES} ${MPI_CXX_LIBRARIES} ${FFTW_LIBRARIES} Threads::Threads ${FILESYSTEM_LIBS} ${OPENMP_LIBS})
    install(TARGETS ${libname} DESTINATION ${CMAKE_HOME_DIRECTORY})

    add_executable (${exename} main.cpp)
    set_target_properties (${exename} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
    target_link_libraries (${exename} ${libname})
    install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})
else ()
    add_executable (${exename} main.cpp ${SOURCES})
    set_target_properties (${exename} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
    target_link_libraries (${exename} ${GSL_LIBRARIES} ${MPI_CXX_LIBRARIES} ${FFTW_LIBRARIES} Threads::Threads ${FILESYSTEM_LIBS} ${OPENMP_LIBS})
    install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})
endif ()

# the HDF5 collection (writeOutputsToHDF5) calls this script; the version
# is printed at start-up and by --version and recorded in
# usedParameters<id>.dat
target_compile_definitions(${exename} PRIVATE
    IPGLASMA_HDF5_SCRIPT="${PROJECT_SOURCE_DIR}/utilities/combine_events_into_hdf5.py"
    IPGLASMA_VERSION="${PROJECT_VERSION}")

# the git description of the source tree, next to the version; recorded at
# every build (CMakeModules/GitDescribe.cmake), so that it follows commits
# without a new configure
set(gitDescribeHeader "${CMAKE_CURRENT_BINARY_DIR}/GitDescribe.h")
add_custom_target(ipglasma_git_describe
    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
        -DOUTPUT=${gitDescribeHeader}
        -P ${PROJECT_SOURCE_DIR}/CMakeModules/GitDescribe.cmake
    BYPRODUCTS ${gitDescribeHeader}
    COMMENT "Recording the git description of the source tree"
    VERBATIM)
add_dependencies(${exename} ipglasma_git_describe)
target_include_directories(${exename} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

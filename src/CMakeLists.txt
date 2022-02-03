set (SOURCES
    Fragmentation.cpp
    FFT.cpp
    Matrix.cpp
    Setup.cpp
    Init.cpp
    Random.cpp
    Group.cpp
    Lattice.cpp
    Cell.cpp
    Glauber.cpp
    Util.cpp
    Evolution.cpp
    GaugeFix.cpp
    Spinor.cpp
    MyEigen.cpp
    pretty_ostream.cpp
    )

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
    target_link_libraries (${libname} ${GSL_LIBRARIES} ${MPI_CXX_LIBRARIES} ${FFTW_LIBRARIES})
    install(TARGETS ${libname} DESTINATION ${CMAKE_HOME_DIRECTORY})

    add_executable (${exename} main.cpp)
    set_target_properties (${exename} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
    target_link_libraries (${exename} ${libname})
    install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})
else ()
    add_executable (${exename} main.cpp ${SOURCES})
    set_target_properties (${exename} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
    target_link_libraries (${exename} ${GSL_LIBRARIES} ${MPI_CXX_LIBRARIES} ${FFTW_LIBRARIES})
    install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})
endif ()

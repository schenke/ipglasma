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
    )

if (build_lib)
    add_library(${libname} SHARED ${SOURCES})
    target_link_libraries (${libname} ${GSL_LIBRARIES} ${MPI_CXX_LIBRARIES} ${FFTW_LIBRARIES})
    install(TARGETS ${libname} DESTINATION ${CMAKE_HOME_DIRECTORY})

    add_executable (${exename} main.cpp)
    target_link_libraries (${exename} ${libname})
    install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})
else ()
    add_executable (${exename} main.cpp ${SOURCES})
    target_link_libraries (${exename} ${GSL_LIBRARIES} ${MPI_CXX_LIBRARIES} ${FFTW_LIBRARIES})
    install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})
endif ()

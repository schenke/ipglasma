# Shared list of IP-Glasma library sources (everything except main.cpp).
# Included by both src/CMakeLists.txt (the ipglasma executable/library) and
# tests/CMakeLists.txt (the unit test executable), so the two targets can
# never drift out of sync with each other.
set (IPGLASMA_LIB_SOURCES
    Fragmentation.cpp
    FFT.cpp
    Matrix.cpp
    Setup.cpp
    Init.cpp
    JIMWLK.cpp
    Random.cpp
    Group.cpp
    Lattice.cpp
    Cell.cpp
    Glauber.cpp
    Util.cpp
    Evolution.cpp
    GaugeFix.cpp
    MyEigen.cpp
    PrettyOstream.cpp
    Parameters.cpp
    Instrumentation.cpp
    )

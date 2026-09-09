#ifndef Lattice_h
#define Lattice_h

#include <string>
#include <vector>

#include "Cell.h"
#include "Matrix.h"
#include "Parameters.h"

// Lattice matrix state is stored structure-of-arrays: every fundamental SU(3)
// field is one contiguous std::vector<Matrix>, and Matrix itself is exactly
// complex<double>[9]. Cell now contains scalar observables only; matrix hot
// paths access these field arrays directly, without Cell pointer chasing.
class Lattice {
  private:
    int size;
    int Nc;

  public:
    Lattice(Parameters *param, int N, int length);
    ~Lattice() = default;
    Lattice(const Lattice &) = delete;
    Lattice &operator=(const Lattice &) = delete;

    int getSize() const { return size; }

    // Fundamental matrix lattice fields. Logical aliases are:
    // U/E1, U2/E2, Ux1/g, Uy1/Uplaq, Ux2/pi, Uy2/phi.
    std::vector<Matrix> U;
    std::vector<Matrix> U2;
    std::vector<Matrix> Ux;
    std::vector<Matrix> Uy;
    std::vector<Matrix> Ux1;
    std::vector<Matrix> Uy1;
    std::vector<Matrix> Ux2;
    std::vector<Matrix> Uy2;

    std::vector<Cell *> cells;
    std::vector<Cell> cellStorage;

    std::vector<int> posmX;
    std::vector<int> pospX;
    std::vector<int> posmY;
    std::vector<int> pospY;

    void WriteWilsonLines(
        std::string fileprefix, Parameters *param, const int iA);
    void WriteSU3Matricies(std::string fileprefix, Parameters *param);
    std::vector<int> posmXpY;
    std::vector<int> pospXmY;
};

class BufferLattice {
  private:
    int size;
    int Nc;

  public:
    BufferLattice(int N, int length);
    ~BufferLattice() = default;
    BufferLattice(const BufferLattice &) = delete;
    BufferLattice &operator=(const BufferLattice &) = delete;

    std::vector<Matrix> buffer1;
    std::vector<Matrix> buffer2;
};

#endif

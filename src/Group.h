#ifndef Group_h
#define Group_h

#include <array>

#include "Matrix.h"

class Group {
  private:
    std::array<Matrix, 8> t;

  public:
    explicit Group(int N);
    ~Group() = default;

    Matrix &getT(int i) { return t[static_cast<std::size_t>(i)]; }
    const Matrix &getT(int i) const { return t[static_cast<std::size_t>(i)]; }
};
#endif

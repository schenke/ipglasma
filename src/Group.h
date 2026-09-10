#ifndef SRC_GROUP_H_
#define SRC_GROUP_H_

#include <array>

#include "Matrix.h"

class Group {
  private:
    std::array<Matrix, 8> t_;

  public:
    Group();
    ~Group() = default;

    Matrix &getT(int i) { return t_[static_cast<std::size_t>(i)]; }
    const Matrix &getT(int i) const { return t_[static_cast<std::size_t>(i)]; }
};
#endif  // SRC_GROUP_H_

#ifndef SRC_GROUP_H_
#define SRC_GROUP_H_

#include <array>

#include "Matrix.h"

/**
 * The eight fundamental-representation SU(3) generators
 * \f$t_a = \lambda_a/2\f$ (\f$\lambda_a\f$ the Gell-Mann matrices),
 * normalized so \f$\mathrm{Tr}(t_a t_b) = \tfrac{1}{2}\delta_{ab}\f$.
 *
 * Constructed once in main.cpp and passed by pointer everywhere a
 * projection onto (or rotation by) a color generator is needed: JIMWLK's
 * Langevin noise generation, Init.cpp's color-charge-density sampling,
 * and Evolution's evolution/multiplicity entry points (though not every
 * holder of a Group* actually uses it -- see e.g. GaugeFix::fftChi()).
 */
class Group {
  private:
    /// The eight generators \f$t_0, \ldots, t_7\f$ (physics convention
    /// \f$t_1, \ldots, t_8\f$), each a \f$3\times3\f$ Hermitian,
    /// traceless matrix.
    std::array<Matrix, 8> t_;

  public:
    /**
     * Constructs a Group with all eight generators set to their
     * standard Gell-Mann-matrix values.
     */
    Group();
    ~Group() = default;

    /**
     * Returns a mutable reference to one generator.
     * \param[in] i Generator index, `0`-`7` (physics convention
     * \f$t_{i+1}\f$).
     * \return Reference to \f$t_i\f$.
     */
    Matrix &getT(int i) { return t_[static_cast<std::size_t>(i)]; }
    /**
     * Returns a read-only reference to one generator.
     * \param[in] i Generator index, `0`-`7` (physics convention
     * \f$t_{i+1}\f$).
     * \return Const reference to \f$t_i\f$.
     */
    const Matrix &getT(int i) const { return t_[static_cast<std::size_t>(i)]; }
};
#endif  // SRC_GROUP_H_

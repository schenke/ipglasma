// NuclearQsTable.h is part of the IP-Glasma solver.

#ifndef SRC_NUCLEARQSTABLE_H_
#define SRC_NUCLEARQSTABLE_H_

#include <string>

#include "PrettyOstream.h"

/**
 * The tabulated saturation scale \f$Q_s^2(T_p, y)\f$ of a nucleus as a
 * function of the summed nucleon thickness \f$T_p\f$ and the rapidity
 * \f$y\f$ (the file `colorCharge.nucleusQsTableFileName`, from IP-Sat).
 * read() loads it, qs2() interpolates it.
 *
 * File format: one line per point, `y T_p Qs^2`, for nT_ values of
 * \f$T_p\f$ (outer loop) times nY_ values of \f$y\f$ (inner loop) in
 * steps of deltaY_.
 */
class NuclearQsTable {
  public:
    /**
     * Reads the table.
     * \param[in] fileName Path of the table file; exits with an error if
     * it doesn't exist or ends prematurely.
     */
    void read(const std::string &fileName);
    /**
     * Bilinearly interpolates the table.
     * \param[in] T Nuclear thickness \f$T_p\f$ to interpolate at.
     * \param[in] y Rapidity to interpolate at; exits with an error if
     * outside the tabulated range \f$0 \le y \le 10.75\f$.
     * \return Interpolated \f$Q_s^2\f$; `0` if \p T is below the
     * tabulated range, clamped to the maximal tabulated \f$T_p\f$ (with
     * a warning) if above it.
     */
    double qs2(double T, double y) const;

  private:
    /// Number of tabulated rapidities.
    static constexpr int nY_ = 44;
    /// Number of tabulated \f$T_p\f$ values (240 since Sep 2026, a 10x
    /// extended \f$T_p\f$ range).
    static constexpr int nT_ = 240;
    /// Rapidity step of the table.
    static constexpr double deltaY_ = 0.25;
    /// \f$Q_s^2\f$ [GeV\f$^2\f$], indexed `[iT][iy]`.
    double Qs2_[nT_][nY_];
    /// The tabulated \f$T_p\f$ values.
    double T_[nT_];
    /// Log sink for progress/error messages.
    PrettyOstream messager_;
};

#endif  // SRC_NUCLEARQSTABLE_H_

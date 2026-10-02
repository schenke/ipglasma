// NucleonModel.h is part of the IP-Glasma solver.

#ifndef SRC_NUCLEONMODEL_H_
#define SRC_NUCLEONMODEL_H_

#include <memory>
#include <string>
#include <vector>

#include "Glauber.h"
#include "Parameters.h"
#include "Random.h"

/**
 * The transverse density \f$T_p(x, y)\f$ of one sampled nucleon, as
 * produced by a NucleonModel. It knows the nucleon's center, so
 * thickness() takes positions on the lattice.
 */
class NucleonProfile {
  public:
    virtual ~NucleonProfile() = default;
    /**
     * The nucleon's thickness at a transverse position.
     * \param[in] x Transverse \f$x\f$ position [fm] on the lattice.
     * \param[in] y Transverse \f$y\f$ position [fm] on the lattice.
     * \return \f$T_p\f$ [GeV\f$^2\f$]; its integral over \f$d^2b\f$
     * [GeV\f$^{-2}\f$] is the nucleon's Qs normalization (1 without
     * `smearQs`).
     */
    virtual double thickness(double x, double y) const = 0;
};

/**
 * A single Gaussian nucleon of width \f$B_G\f$, optionally stretched
 * along the angle \c phi by `protonAnisotropy`.
 */
class GaussianProfile : public NucleonProfile {
  public:
    /**
     * \param[in] xm Nucleon center \f$x\f$ [fm].
     * \param[in] ym Nucleon center \f$y\f$ [fm].
     * \param[in] BG Gaussian width \f$B_G\f$ [GeV\f$^{-2}\f$].
     * \param[in] anisotropy Elongation \f$\xi\f$ (`protonAnisotropy`).
     * \param[in] phi Orientation of the elongation [rad].
     * \param[in] normalization Qs normalization factor.
     */
    GaussianProfile(
        double xm, double ym, double BG, double anisotropy, double phi,
        double normalization)
        : xm_(xm),
          ym_(ym),
          BG_(BG),
          xi_(anisotropy),
          phi_(phi),
          normalization_(normalization) {}
    double thickness(double x, double y) const override;

  private:
    /// Nucleon center \f$x\f$ [fm].
    double xm_;
    /// Nucleon center \f$y\f$ [fm].
    double ym_;
    /// Gaussian width \f$B_G\f$ [GeV\f$^{-2}\f$].
    double BG_;
    /// Elongation \f$\xi\f$.
    double xi_;
    /// Orientation of the elongation [rad].
    double phi_;
    /// Qs normalization factor.
    double normalization_;
};

/**
 * A nucleon made of Gaussian hot spots (constituent quarks) of width
 * \f$B_q\f$; the hot spots share the nucleon's Qs normalization equally.
 */
class HotSpotProfile : public NucleonProfile {
  public:
    /**
     * \param[in] xm Nucleon center \f$x\f$ [fm].
     * \param[in] ym Nucleon center \f$y\f$ [fm].
     * \param[in] xq Hot-spot \f$x\f$ offsets from the center [fm].
     * \param[in] yq Hot-spot \f$y\f$ offsets from the center [fm].
     * \param[in] BGq Hot-spot widths \f$B_q\f$ [GeV\f$^{-2}\f$].
     * \param[in] normalization Qs normalization factor of each hot spot.
     */
    HotSpotProfile(
        double xm, double ym, std::vector<double> xq, std::vector<double> yq,
        std::vector<double> BGq, std::vector<double> normalization)
        : xm_(xm),
          ym_(ym),
          xq_(std::move(xq)),
          yq_(std::move(yq)),
          BGq_(std::move(BGq)),
          normalization_(std::move(normalization)) {}
    double thickness(double x, double y) const override;
    /**
     * Returns the number of hot spots.
     * \return The number of hot spots of this nucleon.
     */
    int numberOfHotSpots() const { return static_cast<int>(xq_.size()); }

  private:
    /// Nucleon center \f$x\f$ [fm].
    double xm_;
    /// Nucleon center \f$y\f$ [fm].
    double ym_;
    /// Hot-spot \f$x\f$ offsets from the center [fm].
    std::vector<double> xq_;
    /// Hot-spot \f$y\f$ offsets from the center [fm].
    std::vector<double> yq_;
    /// Hot-spot widths \f$B_q\f$ [GeV\f$^{-2}\f$].
    std::vector<double> BGq_;
    /// Qs normalization factor of each hot spot.
    std::vector<double> normalization_;
};

/**
 * How a nucleon's transverse structure is sampled, selected by the
 * input parameter `nucleonModel`. A new model (e.g. a stringy proton)
 * only needs a NucleonModel and a NucleonProfile subclass and an entry
 * in create().
 */
class NucleonModel {
  public:
    virtual ~NucleonModel() = default;
    /**
     * Samples the structure of one nucleon.
     * \param[in,out] random Random-number source.
     * \param[in] nucleon The nucleon (its position and orientation
     * \c phi).
     * \return The nucleon's profile.
     */
    virtual std::unique_ptr<NucleonProfile> sample(
        Random &random, const ReturnValue &nucleon) const = 0;
    /**
     * Creates the model selected by `subnucleon.nucleonModel`.
     * \param[in] param Simulation parameters; their current values
     * (including an event's posterior parameter set) are copied.
     * \return The model.
     */
    static std::unique_ptr<NucleonModel> create(const Parameters &param);
};

/// `nucleonModel gaussian`: GaussianProfile nucleons.
class GaussianNucleon : public NucleonModel {
  public:
    /**
     * \param[in] param Simulation parameters (`BG`, `protonAnisotropy`,
     * `smearQs`, `smearingWidth`).
     */
    explicit GaussianNucleon(const Parameters &param);
    std::unique_ptr<NucleonProfile> sample(
        Random &random, const ReturnValue &nucleon) const override;

  private:
    /// Gaussian width \f$B_G\f$ [GeV\f$^{-2}\f$].
    double BG_;
    /// Elongation \f$\xi\f$ (`protonAnisotropy`).
    double anisotropy_;
    /// Whether the Qs normalization fluctuates (`smearQs`).
    bool smearQs_;
    /// Width of the log-normal Qs fluctuations (`smearingWidth`).
    double smearingWidth_;
};

/// `nucleonModel hotspots`: HotSpotProfile nucleons.
class HotSpotNucleon : public NucleonModel {
  public:
    /**
     * \param[in] param Simulation parameters (`BG`, `BGq`, `BGqVar`,
     * `dqMin`, `omega`, `NqBase`, `NqFluc`,
     * `shiftConstituentQuarkProtonOrigin`, `smearQs`, `smearingWidth`).
     */
    explicit HotSpotNucleon(const Parameters &param);
    std::unique_ptr<NucleonProfile> sample(
        Random &random, const ReturnValue &nucleon) const override;
    /**
     * Samples the number of hot spots: floor(NqBase), one more with
     * probability equal to the fractional part, plus a Poisson
     * fluctuation of mean NqFluc; at least 1.
     * \param[in,out] random Random-number source.
     * \return The number of hot spots.
     */
    int sampleNumberOfPartons(Random &random) const;

  private:
    /// Width \f$B_G\f$ of the hot-spot position distribution
    /// [GeV\f$^{-2}\f$].
    double BG_;
    /// Mean hot-spot width \f$B_q\f$ [GeV\f$^{-2}\f$].
    double BGq_;
    /// Variance of the log-normal hot-spot width distribution.
    double BGqVar_;
    /// Minimum distance between hot spots [fm].
    double dqMin_;
    /// Gamma-distribution shape of the radial hot-spot positions; `1`
    /// is a 3D Gaussian.
    double omega_;
    /// Mean number of hot spots.
    double NqBase_;
    /// Mean of the Poisson fluctuation of the number of hot spots.
    double NqFluc_;
    /// Whether the hot spots' center of mass is moved to the center.
    bool shiftOrigin_;
    /// Whether the Qs normalization fluctuates (`smearQs`).
    bool smearQs_;
    /// Width of the log-normal Qs fluctuations (`smearingWidth`).
    double smearingWidth_;
};

#endif  // SRC_NUCLEONMODEL_H_

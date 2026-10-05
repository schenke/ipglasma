// NucleonModel.h is part of the IP-Glasma solver.

#ifndef SRC_NUCLEONMODEL_H_
#define SRC_NUCLEONMODEL_H_

#include <array>
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
    /**
     * \param[in] xm Nucleon center \f$x\f$ [fm].
     * \param[in] ym Nucleon center \f$y\f$ [fm].
     */
    NucleonProfile(double xm, double ym) : xm_(xm), ym_(ym) {}
    /// Destroys the profile.
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

  protected:
    /// Nucleon center \f$x\f$ [fm].
    double xm_;
    /// Nucleon center \f$y\f$ [fm].
    double ym_;
};

/**
 * A single Gaussian nucleon of width \f$B_G\f$, optionally made
 * anisotropic by `protonAnisotropy` \f$\xi\f$: the width along the
 * direction at the angle \c phi becomes \f$B_G/(1+\xi)\f$ (narrower for
 * \f$\xi>0\f$), at the same normalization.
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
        : NucleonProfile(xm, ym),
          BG_(BG),
          xi_(anisotropy),
          phi_(phi),
          normalization_(normalization) {}
    double thickness(double x, double y) const override;

  private:
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
        : NucleonProfile(xm, ym),
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
    /**
     * Returns the hot spots' \f$x\f$ offsets from the nucleon center.
     * \return The offsets [fm], one per hot spot.
     */
    const std::vector<double> &offsetsX() const { return xq_; }
    /**
     * Returns the hot spots' \f$y\f$ offsets from the nucleon center.
     * \return The offsets [fm], one per hot spot.
     */
    const std::vector<double> &offsetsY() const { return yq_; }

  private:
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
 * input parameter `nucleonModel`. A new model only needs a NucleonModel
 * (and, unless an existing one fits, a NucleonProfile) subclass and an
 * entry in create().
 */
class NucleonModel {
  public:
    /// Destroys the model.
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
     * Creates the model selected by `param.subnucleon.nucleonModel`.
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

/**
 * The hot spots of one sampled nucleon, before they are turned into a
 * profile (see HotSpotNucleon::sampleHotSpots()).
 */
struct HotSpotConfiguration {
    /// \f$x\f$ offsets from the nucleon center [fm].
    std::vector<double> x;
    /// \f$y\f$ offsets from the nucleon center [fm].
    std::vector<double> y;
    /// Longitudinal (\f$z\f$) offsets [fm]; `0` for `omega != 1`.
    std::vector<double> z;
    /// Width \f$B_q\f$ shared by the nucleon's hot spots
    /// [GeV\f$^{-2}\f$].
    double BGq;
    /// Qs normalization factor of each hot spot.
    std::vector<double> normalization;
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
    /**
     * Samples the hot spots of one nucleon, in this order: their common
     * width (0.09 GeV\f$^{-2}\f$ plus a log-normal number with mean
     * `BGq` - 0.09 and variance `BGqVar`), their number, their positions
     * (for `omega 1` in 3D, each coordinate Gaussian with variance `BG`;
     * otherwise in the transverse plane at radius
     * \f$\sqrt{\omega x B_G}\f$ with \f$x\f$ from Random::sampleGammaInc()
     * and a uniform angle; keeping a best-effort minimum distance
     * `dqMin`), the optional shift of their center of mass to the nucleon
     * center, and their Qs normalization.
     * \param[in,out] random Random-number source.
     * \param[in] number Number of hot spots; `0` samples it with
     * sampleNumberOfPartons().
     * \return The hot spots, relative to the nucleon center.
     */
    HotSpotConfiguration sampleHotSpots(Random &random, int number = 0) const;

  private:
    /// Width \f$B_G\f$ of the hot-spot position distribution
    /// [GeV\f$^{-2}\f$].
    double BG_;
    /// Mean hot-spot width \f$B_q\f$ [GeV\f$^{-2}\f$], at least 0.09.
    double BGq_;
    /// Variance of the log-normal hot-spot width distribution.
    double BGqVar_;
    /// Minimum distance between hot spots [fm].
    double dqMin_;
    /// Gamma-distribution shape of the radial hot-spot positions; `1`
    /// is a 3D Gaussian, other values place the hot spots in the
    /// transverse plane.
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

/// A point in 3D [fm].
using Point3 = std::array<double, 3>;

/**
 * The strings of one sampled `nucleonModel strings` nucleon (see
 * StringyNucleon::sampleStrings()).
 */
struct StringConfiguration {
    /// The three hot spots the strings end on, relative to the nucleon
    /// center.
    HotSpotConfiguration hotSpots;
    /// Junction where the three strings meet: the Fermat point of the
    /// hot spots [fm].
    Point3 junction;
    /// Position of each moved hot spot along its string, from 0 at the
    /// junction to 1 at the original hot spot.
    std::vector<double> t;
    /// Transverse \f$x\f$ offsets of the moved hot spots from the
    /// nucleon center [fm].
    std::vector<double> x;
    /// Transverse \f$y\f$ offsets of the moved hot spots [fm].
    std::vector<double> y;
};

/**
 * `nucleonModel strings`: a nucleon of three hot spots connected by
 * strings that meet at a junction. The three hot spots are sampled as
 * for `nucleonModel hotspots` (in 3D for `omega 1`), and the junction is
 * their Fermat point, which minimizes the total string length. Each hot
 * spot is then moved to a uniformly random point on its string, from
 * the junction to its sampled position, and the nucleon's thickness is
 * that of hot spots at the projections of these points onto the
 * transverse plane (a HotSpotProfile). The minimum distance `dqMin`
 * applies to the sampled hot spots, the ends of the strings, not to the
 * moved ones, which can come arbitrarily close to each other.
 */
class StringyNucleon : public NucleonModel {
  public:
    /// Number of hot spots, and so of strings, of every nucleon.
    static constexpr int numberOfHotSpots = 3;
    /**
     * \param[in] param Simulation parameters (`BG`, `BGq`, `BGqVar`,
     * `dqMin`, `omega`, `shiftConstituentQuarkProtonOrigin`, `smearQs`,
     * `smearingWidth`).
     */
    explicit StringyNucleon(const Parameters &param);
    std::unique_ptr<NucleonProfile> sample(
        Random &random, const ReturnValue &nucleon) const override;
    /**
     * Samples the strings of one nucleon: the three hot spots (see
     * HotSpotNucleon::sampleHotSpots()), their junction, and the
     * position of each moved hot spot along its string.
     * \param[in,out] random Random-number source.
     * \return The strings, relative to the nucleon center.
     */
    StringConfiguration sampleStrings(Random &random) const;
    /**
     * The Fermat point of three points: the point with the smallest
     * sum of distances to them. If one angle of their triangle is at
     * least 120 degrees, it is that vertex; otherwise it lies in the
     * triangle's plane, where each pair of points is seen at 120
     * degrees, with barycentric weights \f$a/\sin(A+\pi/3)\f$ (\f$a\f$ the
     * side opposite to the vertex with angle \f$A\f$). If two points
     * coincide, it is that point.
     * \param[in] points The three points.
     * \return The Fermat point.
     */
    static Point3 fermatPoint(const std::array<Point3, 3> &points);

  private:
    /// Samples the hot spots the strings end on.
    HotSpotNucleon hotSpots_;
};

#endif  // SRC_NUCLEONMODEL_H_

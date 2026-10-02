// CollisionGeometry.h is part of the IP-Glasma solver.

#ifndef SRC_COLLISIONGEOMETRY_H_
#define SRC_COLLISIONGEOMETRY_H_

#include <vector>

#include "Glauber.h"
#include "Lattice.h"
#include "Parameters.h"
#include "PrettyOstream.h"
#include "Random.h"

/**
 * The collision geometry of one event: samples the impact parameter and
 * reaction-plane angle, determines the wounded nucleons
 * (\f$N_{\text{part}}\f$, \f$N_{\text{coll}}\f$) from the sampled nucleon
 * positions, and computes the overlap-region \f$Q_s\f$ averages,
 * \f$T_{pp}\f$ and the running \f$\alpha_s\f$ from the lattice's color-charge
 * densities, deciding whether the event is accepted.
 *
 * Works on the nucleon lists of the two nuclei passed to the constructor
 * (both centered at the origin); the nucleus A (projectile) is placed at
 * \f$+b/2\f$, the nucleus B (target) at \f$-b/2\f$ along the reaction
 * plane.
 */
class CollisionGeometry {
  public:
    /**
     * Constructs the geometry for two nucleon lists, which must outlive it.
     * \param[in,out] nucleusA Nucleons of nucleus A (projectile); their
     * `.collided` flags are set by computeQuantities().
     * \param[in,out] nucleusB Nucleons of nucleus B (target).
     */
    CollisionGeometry(
        std::vector<ReturnValue> &nucleusA, std::vector<ReturnValue> &nucleusB)
        : nucleusA_(nucleusA), nucleusB_(nucleusB) {}

    /**
     * Samples this event's impact parameter \f$b\f$ (linearly or
     * uniformly distributed between `collision.bMin`/`collision.bMax`, or
     * `0` for the constant-color-charge-density case) and reaction-plane
     * angle, and resets every nucleon's `.collided` flag to `0`.
     * \param[in,out] param Simulation parameters; `event.b`/`event.phiRP`
     * store the sampled values.
     * \param[in,out] random Random-number source.
     */
    void sampleImpactParameter(Parameters *param, Random *random);
    /**
     * Computes and logs this event's collision-geometry summary
     * (\f$N_{\text{part}}\f$, \f$N_{\text{coll}}\f$, \f$T_{pp}\f$,
     * average \f$Q_s\f$, \f$\alpha_s\f$), writes the
     * `usedParameters*.dat`/`NgluonEstimators*.dat` files, and marks
     * the event a success or failure (e.g. no overlap region, no
     * physical \f$Q_s\f$, or \f$Q_{s,\min}^2 S_T\f$ below
     * `param->colorCharge.minimumQs2ST`).
     * \param[in] lat Lattice to read color-charge densities from.
     * \param[in,out] param Simulation parameters.
     * \param[in,out] random Random-number source (Gaussian wounding).
     */
    void computeQuantities(Lattice *lat, Parameters *param, Random *random);
    /**
     * Determines \f$N_{\text{part}}\f$/\f$N_{\text{coll}}\f$ from the
     * nucleon positions, writes `NcollList*.dat`/`NpartList*.dat`, and
     * sets `param->event.Npart`.
     * \param[in,out] param Simulation parameters.
     * \param[in,out] random Random-number source (Gaussian wounding).
     * \param[out] Npart Number of participants.
     * \param[out] Ncoll Number of binary collisions.
     * \return `false` (having set `param->event.success = 0`) if
     * `useFixedNpart` is set and this event's \f$N_{\text{part}}\f$
     * doesn't match, signaling the caller to abort and resample;
     * `true` otherwise.
     */
    bool determineNpartAndNcoll(
        Parameters *param, Random *random, int &Npart, int &Ncoll);
    /**
     * determineNpartAndNcoll()'s binary-collision pair loop: writes
     * `NcollList<id>.dat` and marks each colliding nucleon pair's
     * `.collided`, using either a hard-sphere (\f$d_{ij}^2 <\f$ \p d2)
     * or Gaussian-profile wounding criterion depending on
     * `param->collision.gaussianWounding`.
     * \param[in] param Simulation parameters.
     * \param[in,out] random Random-number source (Gaussian wounding only).
     * \param[in] d2 Squared wounding distance
     * (\f$\sigma_{NN}/(10\pi)\f$) [fm\f$^2\f$].
     * \param[in] b Impact parameter [fm].
     * \param[in] phiRP Reaction-plane angle [rad].
     * \param[in,out] Ncoll Incremented for each colliding pair found.
     */
    void computeNcollList(
        Parameters *param, Random *random, double d2, double b, double phiRP,
        int &Ncoll);
    /**
     * Scans the full lattice, accumulating the \f$Q_s\f$/\f$T_{pp}\f$
     * averages computeQuantities() reports and stores (only over cells
     * within the wounding distance of at least one collided nucleon from
     * each nucleus, or with both thickness profiles nonzero for a smooth
     * nucleus).
     * \param[in] lat Lattice to read color-charge densities from.
     * \param[in] param Simulation parameters.
     * \param[in] N Lattice side length.
     * \param[in] a Lattice spacing [fm].
     * \param[in] b Impact parameter [fm].
     * \param[in] phiRP Reaction-plane angle [rad].
     * \param[out] averageQs Running sum of \f$Q_s\f$ (max of the two
     * nuclei at each cell).
     * \param[out] averageQs2 Running sum of \f$Q_s^2\f$ (max).
     * \param[out] averageQs2Avg Running sum of \f$Q_s^2\f$ (average of
     * the two nuclei).
     * \param[out] averageQs2min Running sum of \f$Q_s^2\f$ (min).
     * \param[out] averageQs2min2 Running sum of \f$Q_s^2\f$ (min),
     * accumulated over every lattice cell rather than just the
     * overlap region.
     * \param[out] Tpp Running sum of \f$T_p^A T_p^B\f$ [fm\f$^{-2}\f$].
     * \param[out] count Number of cells included in the overlap-region
     * averages.
     */
    void scanOverlap(
        Lattice *lat, Parameters *param, int N, double a, double b,
        double phiRP, double &averageQs, double &averageQs2,
        double &averageQs2Avg, double &averageQs2min, double &averageQs2min2,
        double &Tpp, int &count);
    /**
     * Sets \p param's running-coupling \f$\alpha_s\f$ from whichever
     * \f$Q_s\f$ choice `param->coupling.runWithQs` selects, or a
     * fixed value if running coupling is disabled or \f$\alpha_s\f$ runs with
     * \f$k_T\f$ instead (handled per-cell elsewhere via
     * `computeRunningCouplingGfactor`, which shares computeAlphaS()
     * with this function).
     * \param[in,out] param Simulation parameters;
     * `event.alphas` stores the result.
     */
    void computeAndSetRunningAlphaS(Parameters *param);
    /**
     * Logs computeQuantities()'s
     * \f$N_{\text{part}}\f$/\f$N_{\text{coll}}\f$/\f$T_{pp}\f$/\f$Q_s\f$
     * summary.
     * \param[in] param Simulation parameters.
     * \param[in] Npart Number of participants.
     * \param[in] Ncoll Number of binary collisions.
     * \param[in] Tpp \f$T_{pp}\f$ [fm\f$^{-2}\f$].
     * \param[in] a Lattice spacing [fm].
     * \param[in] averageQs2 Average \f$Q_s^2\f$ (max).
     * \param[in] averageQs2Avg Average \f$Q_s^2\f$ (average).
     * \param[in] averageQs2min Average \f$Q_s^2\f$ (min).
     * \param[in] averageQs2min2 Average \f$Q_s^2\f$ (min), whole
     * lattice.
     * \param[in] count Number of cells the overlap-region averages were
     * computed over.
     */
    void logQuantities(
        Parameters *param, int Npart, int Ncoll, double Tpp, double a,
        double averageQs2, double averageQs2Avg, double averageQs2min,
        double averageQs2min2, int count);
    /**
     * Appends this event's collision geometry to `usedParameters<id>.dat`
     * as comment lines (called only when computeQuantities() marks the
     * event a success).
     * \param[in] param Simulation parameters.
     * \param[in] phiRP Reaction-plane angle [rad].
     * \param[in] Npart Number of participants.
     * \param[in] Ncoll Number of binary collisions.
     */
    void writeUsedParametersFile(
        Parameters *param, double phiRP, int Npart, int Ncoll);
    /**
     * Writes this event's `NgluonEstimators<id>.dat` file (rough
     * gluon-multiplicity estimator inputs: \f$Q_s^2 S_T\f$ combinations
     * for the min/average/max \f$Q_s\f$ choices).
     * \param[in] param Simulation parameters.
     * \param[in] a Lattice spacing [fm].
     * \param[in] averageQs2 Average \f$Q_s^2\f$ (max).
     * \param[in] averageQs2Avg Average \f$Q_s^2\f$ (average).
     * \param[in] averageQs2min2 Average \f$Q_s^2\f$ (min), whole
     * lattice.
     * \param[in] count Number of cells the overlap-region averages were
     * computed over.
     */
    void writeNgluonEstimatorsFile(
        Parameters *param, double a, double averageQs2, double averageQs2Avg,
        double averageQs2min2, int count);

  private:
    /// Nucleons of nucleus A (projectile), not owned.
    std::vector<ReturnValue> &nucleusA_;
    /// Nucleons of nucleus B (target), not owned.
    std::vector<ReturnValue> &nucleusB_;
    /// Log sink for progress/warning messages.
    PrettyOstream messager_;
};

#endif  // SRC_COLLISIONGEOMETRY_H_

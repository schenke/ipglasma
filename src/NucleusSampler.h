// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#ifndef SRC_NUCLEUSSAMPLER_H_
#define SRC_NUCLEUSSAMPLER_H_

#include <string>
#include <vector>

#include "Glauber.h"
#include "Parameters.h"
#include "PrettyOstream.h"
#include "Random.h"

/// The sampled nucleons of both nuclei, see NucleusSampler::sample().
struct Nuclei {
    /// Nucleons of the projectile (nucleus A).
    std::vector<ReturnValue> projectile;
    /// Nucleons of the target (nucleus B).
    std::vector<ReturnValue> target;
};

/// A radius and polar angle sampled together, see
/// NucleusSampler::sampleRAndCosthetaFromDeformedWoodsSaxon().
struct RadiusAndCosTheta {
    /// Radius [fm].
    double r;
    /// \f$\cos\theta\f$.
    double cosTheta;
};

/**
 * The configuration file `nucleonPositionsFromFile 1` reads for a
 * nucleus other than the deuteron.
 * \param[in] nucleusA Mass number.
 * \param[in] lightNucleusOption The input `lightNucleusOption`; only
 * matters for \f$A\f$ = 3, 12, 16, 20 and 40.
 * \return The file name, without the directory; empty if there is no
 * file for this mass number or the option does not exist for it.
 */
std::string configurationFileName(int nucleusA, int lightNucleusOption);

/**
 * Checks that `lightNucleusOption` selects a configuration file of a
 * nucleus (with `nucleonPositionsFromFile 1`).
 * \param[in] nucleusA Mass number.
 * \param[in] lightNucleusOption The input `lightNucleusOption`.
 * \return "" if the option exists for \p nucleusA, or if the nucleus has
 * only one or no configuration file; otherwise a message listing the
 * allowed values.
 */
std::string lightNucleusOptionError(int nucleusA, int lightNucleusOption);

/**
 * Samples the nucleon positions of the projectile and the target: from
 * (possibly deformed) Woods-Saxon distributions, or by drawing one of the
 * pre-tabulated configurations of light nuclei and Au/Pb
 * (`nucleonPositionsFromFile`), followed by the global rotation
 * selected by the nucleus' polarization. Both nuclei are centered at the
 * origin; the impact parameter is applied later.
 */
class NucleusSampler {
  public:
    /**
     * Loads the pre-tabulated configurations of both nuclei
     * (readConfigurationFile()); call once before sample() when
     * `param->nucleus.nucleonPositionsFromFile` is set.
     * \param[in] param Simulation parameters.
     * \param[in] glauber Configured Glauber instance (mass numbers).
     * \param[in,out] random Random-number source (picks the deuteron
     * polarization file for an unpolarized deuteron).
     */
    void readConfigurations(
        Parameters *param, Glauber *glauber, Random *random);
    /**
     * Samples this event's nucleon positions of both nuclei (from the
     * loaded configurations or from Woods-Saxon distributions, depending
     * on `param->nucleus.nucleonPositionsFromFile`) and applies each
     * nucleus' global polarization rotation (applyPolarizationRotation()).
     * With `param->collision.nucleiToAverage` \f$n > 1\f$, \f$n\f$
     * nuclei of each kind are sampled and rotated independently, and
     * each list holds all their nucleons (Init divides every nucleon's
     * thickness by \f$n\f$).
     * \param[in] param Simulation parameters.
     * \param[in,out] random Random-number source.
     * \param[in] glauber Configured Glauber instance providing the
     * nuclear geometry.
     * \return The nucleons of both nuclei.
     */
    Nuclei sample(Parameters *param, Random *random, Glauber *glauber);
    /**
     * Samples the nucleon positions of one nucleus, from the loaded
     * configurations (sampleFromConfigurations()) or from its Woods-Saxon
     * distribution (sampleWoodsSaxon()), depending on
     * `param->nucleus.nucleonPositionsFromFile`.
     * \param[in] param Simulation parameters.
     * \param[in,out] random Random-number source.
     * \param[in] glauber Configured Glauber instance providing the
     * nuclear geometry.
     * \param[in] role Which nucleus to sample.
     * \return The sampled nucleons, centered at the origin.
     */
    std::vector<ReturnValue> samplePositions(
        Parameters *param, Random *random, Glauber *glauber, NucleusRole role);
    /**
     * Samples one nucleus' nucleon positions without configuration
     * files: a single proton (\f$A=1\f$) is placed at the origin, the
     * deuteron (\f$A=2\f$) uses Glauber::sampleTARejection() and its
     * neutron-proton distance convention, and everything else is
     * generate()d.
     * \param[in,out] random Random-number source.
     * \param[in] glauber Configured Glauber instance providing the
     * nuclear geometry.
     * \param[in] role Which nucleus to sample.
     * \return The sampled nucleons.
     */
    std::vector<ReturnValue> sampleWoodsSaxon(
        Random *random, Glauber *glauber, NucleusRole role);
    /**
     * Draws one nucleus' nucleon positions from a random pre-tabulated
     * configuration, assigns protons and recenters it; falls back to
     * generate() if \p configs is empty (no file for this species).
     * \param[in,out] random Random-number source.
     * \param[in] data The nucleus' species and Woods-Saxon parameters.
     * \param[in] configs The loaded configurations, three coordinates per
     * nucleon (see readConfigurationFile()).
     * \return The sampled nucleons.
     */
    std::vector<ReturnValue> sampleFromConfigurations(
        Random *random, const Nucleus &data,
        const std::vector<std::vector<float>> &configs);
    /**
     * Loads pre-tabulated binary nucleon configurations for light
     * nuclei (deuteron through Pb-208), selecting the file by \p
     * nucleusA/\p lightNucleusOption (and, for the deuteron, \p
     * polarizationFlag/\p polJz). Returns no configurations if \p
     * nucleusA isn't one of the supported species (in which case nucleon
     * positions are generated instead, see sampleFromConfigurations()).
     * \param[in] nucleusA Mass number of the nucleus to load
     * configurations for.
     * \param[in] lightNucleusOption Selects which configuration variant
     * to load for species with more than one available (Woods-Saxon,
     * variational Monte Carlo, alpha clusters, PGCM, NLEFT).
     * \param[in] polarizationFlag Deuteron only: `0` samples a random
     * polarization file per call.
     * \param[in] polJz Deuteron only: selects the \f$J_z=\pm1\f$
     * configuration file if `|polJz| == 1`.
     * \param[in,out] random Random-number source, used only for the
     * deuteron with \p polarizationFlag `0`.
     * \param[in] param Simulation parameters;
     * `param->nucleus.nuclearConfigurationsPath` gives the directory to read
     * from.
     * \return One row per configuration read from the file, three
     * coordinates per nucleon.
     */
    std::vector<std::vector<float>> readConfigurationFile(
        const int nucleusA, const int lightNucleusOption,
        const int polarizationFlag, const double polJz, Random *random,
        Parameters *param);
    /**
     * Dispatches to the appropriate nucleus-configuration generator
     * based on which deformation parameters of \p data are nonzero:
     * generateWoodsSaxon() if none are, otherwise
     * generateDeformedWoodsSaxonForceDmin() if `forceDminFlag` is set,
     * else generateTriaxialWoodsSaxon() for \f$\gamma\neq0\f$ and
     * generateDeformedWoodsSaxon() for \f$\gamma=0\f$.
     * \param[in,out] random Random-number source.
     * \param[in] data The nucleus' species and (deformed) Woods-Saxon
     * parameters.
     * \return The generated nucleons.
     */
    std::vector<ReturnValue> generate(Random *random, const Nucleus &data);
    /**
     * Generates an undeformed (spherically symmetric) Woods-Saxon
     * nucleon configuration: samples each nucleon's radius via
     * sampleRFromWoodsSaxon() (neutrons with radius and diffuseness
     * shifted by `dR_np`/`da_np`), then places them at random angles
     * subject to a best-effort (up to 100 retries) minimum-distance
     * `d_min` rejection, and recenters the result.
     * \param[in,out] random Random-number source.
     * \param[in] data The nucleus' species and Woods-Saxon parameters.
     * \return The generated nucleons.
     */
    std::vector<ReturnValue> generateWoodsSaxon(
        Random *random, const Nucleus &data);
    /**
     * Generates an axially symmetric (\f$\gamma=0\f$) deformed
     * Woods-Saxon nucleon configuration: samples each nucleon's
     * `(r, cos\theta)` jointly via
     * sampleRAndCosthetaFromDeformedWoodsSaxon(), then \f$\phi\f$
     * uniformly, subject to the same best-effort minimum-distance
     * rejection as generateWoodsSaxon().
     * \param[in,out] random Random-number source.
     * \param[in] data The nucleus' species and deformed Woods-Saxon
     * parameters.
     * \return The generated nucleons.
     */
    std::vector<ReturnValue> generateDeformedWoodsSaxon(
        Random *random, const Nucleus &data);
    /**
     * generate()'s `forceDminFlag` variant: samples each nucleon's full
     * `(r, \theta, \phi)` jointly against the (possibly triaxial,
     * \f$\gamma\f$-dependent) deformed Woods-Saxon surface, and strictly
     * enforces `d_min` by fully resampling any candidate position closer
     * than `d_min` to an already-placed nucleon (no retry cap).
     * \param[in,out] random Random-number source.
     * \param[in] data The nucleus' species and deformed Woods-Saxon
     * parameters.
     * \return The generated nucleons.
     */
    std::vector<ReturnValue> generateDeformedWoodsSaxonForceDmin(
        Random *random, const Nucleus &data);
    /**
     * generate()'s triaxial (\f$\gamma\neq0\f$), non-forced-\f$d_{\min}\f$
     * variant: like generateDeformedWoodsSaxonForceDmin()'s per-nucleon
     * `(r, \theta, \phi)` sampling against the triaxial surface, but
     * without any minimum-distance rejection. The other samplers keep
     * \f$d_{\min}\f$ on a best-effort basis by redrawing only the angles
     * the density does not depend on; a triaxial density depends on
     * \f$\phi\f$ as well, so that would bias it. For a minimum distance,
     * use `forceDMin 1`.
     * \param[in,out] random Random-number source.
     * \param[in] data The nucleus' species and deformed Woods-Saxon
     * parameters.
     * \return The generated nucleons.
     */
    std::vector<ReturnValue> generateTriaxialWoodsSaxon(
        Random *random, const Nucleus &data);
    /**
     * Samples one nucleon's radius from an undeformed Woods-Saxon
     * distribution via rejection sampling: draws \f$r\f$ with density
     * \f$\propto r^2\f$ (the correct 3D volume-element weighting, via
     * the cube root of a uniform draw) up to a generous cutoff, and
     * accepts it with probability given by fermiDistribution().
     * \param[in,out] random Random-number source.
     * \param[in] a_WS Surface diffuseness [fm].
     * \param[in] R_WS Half-density radius [fm].
     * \return The sampled radius [fm].
     */
    static double sampleRFromWoodsSaxon(
        Random *random, double a_WS, double R_WS);
    /**
     * Samples one nucleon's radius and polar angle jointly from an
     * axially symmetric (\f$\gamma=0\f$) deformed Woods-Saxon
     * distribution: draws \f$r\f$ with density \f$\propto r^2\f$ (via
     * the cube root of a uniform draw, as in sampleRFromWoodsSaxon())
     * up to a generous cutoff, \f$\cos\theta\f$ uniformly, and accepts
     * the pair with probability given by fermiDistribution() evaluated
     * at the angle-dependent surface radius \f$R(\theta) =
     * R_{WS}(1+\beta_2 Y_{20}+\beta_3 Y_{30}+\beta_4 Y_{40})\f$.
     * \param[in,out] random Random-number source.
     * \param[in] a_WS Surface diffuseness [fm].
     * \param[in] R_WS Half-density radius [fm].
     * \param[in] beta2 Quadrupole deformation.
     * \param[in] beta3 Octupole deformation.
     * \param[in] beta4 Hexadecapole deformation.
     * \return The sampled radius and \f$\cos\theta\f$.
     */
    static RadiusAndCosTheta sampleRAndCosthetaFromDeformedWoodsSaxon(
        Random *random, double a_WS, double R_WS, double beta2, double beta3,
        double beta4);
    /**
     * Evaluates the Woods-Saxon (Fermi) density profile.
     * \param[in] r Radius [fm].
     * \param[in] R_WS Half-density radius [fm].
     * \param[in] a_WS Surface diffuseness [fm].
     * \return \f$1/(1+\exp((r-R_{WS})/a_{WS}))\f$.
     */
    static double fermiDistribution(double r, double R_WS, double a_WS);
    /**
     * Evaluates the axially symmetric (\f$m=0\f$) real spherical
     * harmonic \f$Y_{l0}(\theta)\f$ for \f$l=2\f$, `3`, or `4`.
     * \param[in] l Degree; must be `2`, `3`, or `4` (returns `0`
     * otherwise).
     * \param[in] ct \f$\cos\theta\f$.
     * \return \f$Y_{l0}(\theta)\f$.
     */
    static double sphericalHarmonics(int l, double ct);
    /**
     * Evaluates the real spherical harmonic \f$Y_{22}(\theta,\phi)\f$,
     * used for triaxial (\f$\gamma\neq0\f$) nuclear deformation.
     * \param[in] ct \f$\cos\theta\f$.
     * \param[in] phi Azimuthal angle \f$\phi\f$ [rad].
     * \return \f$Y_{22}(\theta,\phi)\f$.
     */
    static double sphericalHarmonicsY22(double ct, double phi);
    /**
     * Shifts a set of coordinates so their center of mass sits at the
     * origin.
     * \param[in,out] x \f$x\f$ coordinates, shifted in place.
     * \param[in,out] y \f$y\f$ coordinates, shifted in place.
     * \param[in,out] z \f$z\f$ coordinates, shifted in place.
     */
    static void recenter(
        std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);
    /**
     * Shifts a nucleus' nucleon positions so their center of mass sits
     * at the origin.
     * \param[in,out] nucleus Nucleon positions, shifted in place.
     */
    static void recenter(std::vector<ReturnValue> &nucleus);
    /**
     * Randomly assigns \p Z of a nucleus' nucleons to be protons (the
     * rest neutrons), via a Fisher-Yates shuffle of the nucleon list
     * followed by labeling the first \p Z entries.
     * \param[in,out] random Random-number source.
     * \param[in,out] nucleus Nucleon positions; shuffled in place, with
     * `.proton` set on every entry.
     * \param[in] Z Number of protons (`std::abs(Z)` is used, in case a
     * negative value is ever passed).
     */
    static void assignProtons(
        Random *random, std::vector<ReturnValue> &nucleus, const int Z);
    /**
     * Rotates a nucleus' nucleon positions by a fixed azimuthal angle
     * \p phi_global and polar angle \p theta_global (used for
     * longitudinal and transverse polarization, where the rotation
     * axis/angle is determined rather than random).
     * \param[in] phi_global Azimuthal rotation angle [rad].
     * \param[in] theta_global Polar rotation angle [rad].
     * \param[in,out] nucleus Nucleon positions, rotated in place.
     */
    static void rotate(
        double phi_global, double theta_global,
        std::vector<ReturnValue> &nucleus);
    /**
     * Rotates a nucleus' nucleon positions by a uniformly random
     * 3D (Euler-angle) rotation, as needed for an unpolarized
     * (including triaxially deformed) nucleus.
     * \param[in,out] random Random-number source.
     * \param[in,out] nucleus Nucleon positions, rotated in place.
     */
    static void rotateRandomly(
        Random *random, std::vector<ReturnValue> &nucleus);
    /**
     * Applies sample()'s global nucleus rotation for one polarization
     * flag; called once each for the projectile and target.
     * \param[in,out] random Random-number source.
     * \param[in] polarizationFlag `0` random 3D rotation
     * (rotateRandomly()), `1` longitudinal polarization (random
     * azimuthal rotation only), `2` transverse polarization (rotates
     * \f$J\f$ to the \f$+y\f$ axis).
     * \param[in,out] nucleus Nucleon positions to rotate in place.
     */
    static void applyPolarizationRotation(
        Random *random, int polarizationFlag,
        std::vector<ReturnValue> &nucleus);

  private:
    /// Pre-tabulated nucleon configurations of the projectile, loaded by
    /// readConfigurations() (empty if there is no file for the species).
    std::vector<std::vector<float>> configsA_;
    /// Same as \c configsA_, for the target.
    std::vector<std::vector<float>> configsB_;
    /// Log sink for progress/warning/error messages.
    PrettyOstream messager_;
};

#endif  // SRC_NUCLEUSSAMPLER_H_

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

//
// The single list of every input-file parameter: its key, how it is stored
// in Parameters, its default (if it is optional) and its validity checks.
// Parameters::readInput() and Parameters::writeInputParameters() are both
// driven by this table, so adding a parameter only takes one entry here
// (plus its field in one of the groups in Parameters.h).

#include <algorithm>
#include <cctype>
#include <charconv>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "InputFile.h"
#include "Parameters.h"
#include "PhysConst.h"

namespace {

/// Predicate on the parameters read so far, deciding whether a parameter
/// is read at all (see Param::onlyIf()).
using Condition = bool (*)(const Parameters &);

/// One input parameter, with its type erased.
struct ParameterSpec {
    /// The parameter's input-file key.
    std::string name;
    /// Text of the default value; empty for a required parameter.
    std::string defaultValue;
    /// If set, the parameter is only read (and required) when this returns
    /// true for the parameters read so far; otherwise its key is accepted
    /// but ignored.
    Condition condition = nullptr;
    /// Whether a posterior file can set the parameter (see
    /// Param::posterior()).
    bool posterior = false;
    /// Parses \p text into \p param; returns an error message, or "" on
    /// success.
    std::function<std::string(Parameters &param, const std::string &text)>
        apply;
    /// Writes the current value in input-file syntax.
    std::function<void(const Parameters &param, std::ostream &out)> write;
};

/**
 * Writes \p value as the shortest text that reads back to exactly the
 * same double.
 * \param[out] out Stream to write to.
 * \param[in] value Value to write.
 */
void writeValue(std::ostream &out, double value) {
    char buffer[64];
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
    out.write(buffer, result.ptr - buffer);
}
/**
 * Writes \p values comma-separated (`none` if empty), in input-file
 * syntax.
 * \param[out] out Stream to write to.
 * \param[in] values Values to write.
 */
void writeValue(std::ostream &out, const std::vector<double> &values) {
    if (values.empty()) out << "none";
    for (std::size_t i = 0; i < values.size(); i++) {
        if (i > 0) out << ",";
        writeValue(out, values[i]);
    }
}
/**
 * Writes \p value as `0` or `1`, in input-file syntax.
 * \param[out] out Stream to write to.
 * \param[in] value Value to write.
 */
void writeValue(std::ostream &out, bool value) { out << (value ? 1 : 0); }
/**
 * Writes any other \p value (integers, strings) with `operator<<`.
 * \tparam T Type of the value.
 * \param[out] out Stream to write to.
 * \param[in] value Value to write.
 */
template <typename T>
void writeValue(std::ostream &out, const T &value) {
    out << value;
}

/**
 * Describes the values a parameter of type \p T accepts, for error
 * messages.
 * \tparam T Type of the parameter.
 * \return E.g. "an integer" or "0 or 1".
 */
template <typename T>
std::string typeName() {
    if constexpr (std::is_same_v<T, bool>) return "0 or 1";
    if constexpr (std::is_same_v<T, int>) return "an integer";
    if constexpr (std::is_same_v<T, long long>) return "an integer";
    if constexpr (std::is_same_v<T, double>) return "a number";
    if constexpr (std::is_same_v<T, std::vector<double>>)
        return "a comma-separated list of numbers, or none";
    return "a value";
}

/// A check on a parsed value: returns "" if it is valid, otherwise the
/// reason it is not (e.g. "must be positive").
template <typename T>
using Check = std::function<std::string(const T &)>;

/**
 * Typed builder for a ParameterSpec, so the table reads
 * `param("size", &P::lattice, &LatticeParameters::size).check(even())`.
 * \tparam T Type of the parameter's field.
 */
template <typename T>
class Param {
  public:
    /**
     * Describes a required parameter stored in `param.*group.*field`.
     * \tparam Group Type of the parameter group (e.g. LatticeParameters).
     * \param[in] name The parameter's input-file key.
     * \param[in] group The group member of Parameters.
     * \param[in] field The field within the group.
     */
    template <typename Group>
    Param(const char *name, Group Parameters::*group, T Group::*field)
        : name_(name),
          set_([group, field](Parameters &p, const T &v) {
              (p.*group).*field = v;
          }),
          write_([group, field](const Parameters &p, std::ostream &out) {
              writeValue(out, (p.*group).*field);
          }) {}

    /**
     * Makes the parameter optional.
     * \param[in] defaultValue Its default, as input-file text.
     * \return This builder.
     */
    Param &optional(const char *defaultValue) {
        defaultValue_ = defaultValue;
        return *this;
    }
    /**
     * Only reads (and requires) the parameter when \p condition holds;
     * otherwise its key is accepted but ignored.
     * \param[in] condition Predicate on the parameters read before it.
     * \return This builder.
     */
    Param &onlyIf(Condition condition) {
        condition_ = condition;
        return *this;
    }
    /**
     * Lets a posterior file set the parameter, as a column of its sets;
     * only for physical parameters, not for those of the lattice, the
     * time steps, the collision system or the switches between models.
     * \return This builder.
     */
    Param &posterior() {
        posterior_ = true;
        return *this;
    }
    /**
     * Adds a check that every value of the parameter must pass.
     * \param[in] check The check, e.g. positive() or even().
     * \return This builder.
     */
    Param &check(Check<T> check) {
        checks_.push_back(std::move(check));
        return *this;
    }

    /**
     * Converts the builder into the type-erased table entry.
     * \return The ParameterSpec that parses, checks, stores and writes
     * the parameter.
     */
    operator ParameterSpec() const {
        ParameterSpec spec;
        spec.name = name_;
        spec.defaultValue = defaultValue_;
        spec.condition = condition_;
        spec.posterior = posterior_;
        spec.apply = [set = set_, checks = checks_](
                         Parameters &p, const std::string &text) {
            T value {};
            if (!parseValue(text, value)) {
                return "'" + text + "' is not " + typeName<T>();
            }
            for (const auto &check : checks) {
                const std::string problem = check(value);
                if (!problem.empty()) return text + ": " + problem;
            }
            set(p, value);
            return std::string();
        };
        spec.write = write_;
        return spec;
    }

  private:
    /// The parameter's input-file key.
    std::string name_;
    /// Default as input-file text; empty for a required parameter.
    std::string defaultValue_;
    /// Condition for reading the parameter; `nullptr` to always read it.
    Condition condition_ = nullptr;
    /// Whether a posterior file can set the parameter.
    bool posterior_ = false;
    /// Stores a parsed value in its field.
    std::function<void(Parameters &, const T &)> set_;
    /// Writes the field's value in input-file syntax.
    std::function<void(const Parameters &, std::ostream &)> write_;
    /// Checks every parsed value must pass.
    std::vector<Check<T>> checks_;
};

/**
 * Starts a table entry for the parameter stored in `param.*group.*field`.
 * \tparam Group Type of the parameter group (e.g. LatticeParameters).
 * \tparam T Type of the field.
 * \param[in] name The parameter's input-file key.
 * \param[in] group The group member of Parameters.
 * \param[in] field The field within the group.
 * \return A builder for a required parameter; see Param.
 */
template <typename Group, typename T>
Param<T> param(const char *name, Group Parameters::*group, T Group::*field) {
    return Param<T>(name, group, field);
}

// ---- checks ----

/**
 * Check that a value is positive.
 * \tparam T Type of the value.
 * \return The check.
 */
template <typename T = double>
Check<T> positive() {
    return [](const T &v) { return v > 0 ? "" : "must be positive"; };
}

/**
 * Check that a value is not negative.
 * \tparam T Type of the value.
 * \return The check.
 */
template <typename T = double>
Check<T> nonNegative() {
    return [](const T &v) { return v >= 0 ? "" : "must not be negative"; };
}

/**
 * Check that an integer is even.
 * \return The check.
 */
Check<int> even() {
    return [](const int &v) { return v % 2 == 0 ? "" : "must be even"; };
}

/**
 * Check that a value is one of \p allowed.
 * \tparam T Type of the value.
 * \param[in] allowed The allowed values.
 * \param[in] note Appended to the error message, e.g. to explain the
 * values.
 * \return The check.
 */
template <typename T>
Check<T> oneOf(std::initializer_list<T> allowed, const char *note = "") {
    std::vector<T> values(allowed);
    return [values, note = std::string(note)](const T &v) {
        if (std::find(values.begin(), values.end(), v) != values.end()) {
            return std::string();
        }
        std::ostringstream message;
        message << "must be one of ";
        for (std::size_t i = 0; i < values.size(); i++) {
            message << (i == 0 ? "" : ", ") << values[i];
        }
        message << note;
        return message.str();
    };
}

/**
 * Check that a value is at least \p low.
 * \tparam T Type of the value.
 * \param[in] low Smallest allowed value.
 * \param[in] note Appended to the error message, e.g. to explain the
 * bound.
 * \return The check.
 */
template <typename T>
Check<T> atLeast(T low, const char *note = "") {
    return [low, note = std::string(note)](const T &v) {
        if (v >= low) return std::string();
        std::ostringstream message;
        message << "must be at least " << low << note;
        return message.str();
    };
}

/**
 * Check that a value lies in [\p low, \p high].
 * \tparam T Type of the value.
 * \param[in] low Smallest allowed value.
 * \param[in] high Largest allowed value.
 * \param[in] note Appended to the error message, e.g. to explain the
 * range.
 * \return The check.
 */
template <typename T>
Check<T> inRange(T low, T high, const char *note = "") {
    return [low, high, note = std::string(note)](const T &v) {
        if (v >= low && v <= high) return std::string();
        std::ostringstream message;
        message << "must be between " << low << " and " << high << note;
        return message.str();
    };
}

/**
 * Condition for the Woods-Saxon and deformation parameters.
 * \param[in] p The parameters read so far.
 * \return Whether `useInputWSParams` is set.
 */
bool wsDeformParamsSet(const Parameters &p) {
    return p.nucleus.useInputWSParams;
}
/**
 * Condition for the impact parameter and the wounded nucleons.
 * \param[in] p The parameters read so far.
 * \return Whether `useNucleus` is 1.
 */
bool nucleiCollide(const Parameters &p) { return p.collision.useNucleus; }

/**
 * Condition for sampling the nucleon positions of both nuclei, which reading
 * the Wilson lines from file skips.
 * \param[in] p The parameters read so far.
 * \return Whether `useNucleus` is 1 and `readInitialWilsonLines` is 0.
 */
bool nucleonsSampled(const Parameters &p) {
    return p.collision.useNucleus && p.wilsonLines.readInitialWilsonLines == 0;
}

/**
 * Condition for the configuration files: they are read with
 * `nucleonPositionsFromFile 1`, which a polarized nucleus also implies, when
 * the nucleon positions are sampled.
 * \param[in] p The parameters read so far.
 * \return Whether configuration files can be read.
 */
bool configurationFilesRead(const Parameters &p) {
    return nucleonsSampled(p)
           && (p.nucleus.nucleonPositionsFromFile
               || p.nucleus.polarizationProjectile != 0
               || p.nucleus.polarizationTarget != 0);
}

/**
 * Condition for the deuteron's \f$J_z\f$ of the projectile.
 * \param[in] p The parameters read so far.
 * \return Whether the projectile is polarized and its nucleons are sampled.
 */
bool projectilePolarized(const Parameters &p) {
    return nucleonsSampled(p) && p.nucleus.polarizationProjectile != 0;
}

/**
 * Condition for the deuteron's \f$J_z\f$ of the target.
 * \param[in] p The parameters read so far.
 * \return Whether the target is polarized and its nucleons are sampled.
 */
bool targetPolarized(const Parameters &p) {
    return nucleonsSampled(p) && p.nucleus.polarizationTarget != 0;
}

/**
 * Condition for the Jacobian mass of the pseudorapidity conversion.
 * \param[in] p The parameters read so far.
 * \return Whether `usePseudoRapidity` is set.
 */
bool pseudoRapidityUsed(const Parameters &p) {
    return p.colorCharge.usePseudoRapidity;
}

/**
 * Condition for \f$\sqrt{s}\f$, which enters the fluctuating \f$x\f$
 * and the pseudorapidity Jacobian.
 * \param[in] p The parameters read so far.
 * \return Whether `useFluctuatingX` or `usePseudoRapidity` is set.
 */
bool sqrtSUsed(const Parameters &p) {
    return p.colorCharge.useFluctuatingX || p.colorCharge.usePseudoRapidity;
}

/**
 * Condition for the running-coupling parameters.
 * \param[in] p The parameters read so far.
 * \return Whether `runningCoupling` is set.
 */
bool couplingRuns(const Parameters &p) { return p.coupling.runningCoupling; }

/**
 * Condition for the seed of the color charges.
 * \param[in] p The parameters read so far.
 * \return Whether the color charges are sampled (not read).
 */
bool colorChargesSampled(const Parameters &p) {
    return p.wilsonLines.readInitialWilsonLines == 0;
}

/**
 * Condition for the directory of the Wilson-line files.
 * \param[in] p The parameters read so far.
 * \return Whether Wilson lines are written or read.
 */
bool wilsonLinesWrittenOrRead(const Parameters &p) {
    return p.wilsonLines.writeWilsonLines != 0
           || p.wilsonLines.readInitialWilsonLines != 0;
}

/**
 * Condition for the JIMWLK parameters.
 * \param[in] p The parameters read so far.
 * \return Whether `useJIMWLK` is set.
 */
bool jimwlkEnabled(const Parameters &p) { return p.jimwlk.enabled; }

/**
 * Condition for the JIMWLK running-coupling prescription.
 * \param[in] p The parameters read so far.
 * \return Whether `useJIMWLK` is set and `jimwlkAlphaS` is 0.
 */
bool jimwlkCouplingRuns(const Parameters &p) {
    return p.jimwlk.enabled && PhysConst::isClose(p.jimwlk.alphaS, 0.);
}

/**
 * Condition for the constant color-charge density.
 * \param[in] p The parameters read so far.
 * \return Whether `useNucleus` is 0.
 */
bool constantColorCharge(const Parameters &p) {
    return !p.collision.useNucleus;
}

/**
 * Condition for the JIMWLK snapshot list.
 * \param[in] p The parameters read so far.
 * \return Whether `jimwlkSaveSnapshots` is set.
 */
bool saveSnapshotsSet(const Parameters &p) { return p.jimwlk.saveSnapshots; }

/**
 * Condition for the hadron spectrum, which is computed from the gluon
 * spectrum.
 * \param[in] p The parameters read so far.
 * \return Whether `computeGluonMultiplicity` is set.
 */
bool gluonMultiplicityComputed(const Parameters &p) {
    return p.output.computeGluonMultiplicity;
}

/**
 * Condition for the eccentricity cutoff.
 * \param[in] p The parameters read so far.
 * \return Whether `computeEccentricities` is set.
 */
bool eccentricitiesComputed(const Parameters &p) {
    return p.output.computeEccentricities;
}

/**
 * Condition for the output times and the output grid.
 * \param[in] p The parameters read so far.
 * \return Whether a hydro, Jazma or T^{mu nu} output is switched on.
 */
bool fieldOutputWritten(const Parameters &p) {
    return p.output.anyFieldOutput();
}

/**
 * Condition for the geometry files written with the Wilson lines.
 * \param[in] p The parameters read so far.
 * \return Whether `writeWilsonLines` is 1 or 2.
 */
bool wilsonLinesWritten(const Parameters &p) {
    return p.wilsonLines.writeWilsonLines != 0;
}

/**
 * Condition for the x of the Wilson lines to read.
 * \param[in] p The parameters read so far.
 * \return Whether `readInitialWilsonLines` is 1 or 2.
 */
bool wilsonLinesRead(const Parameters &p) {
    return p.wilsonLines.readInitialWilsonLines != 0;
}

/**
 * Condition for the fixed x of the nuclei; with a fluctuating x, x comes
 * from the local Q_s (and useJIMWLK 1 requires useFluctuatingX 0).
 * \param[in] p The parameters read so far.
 * \return Whether `useFluctuatingX` is 0.
 */
bool fixedX(const Parameters &p) { return !p.colorCharge.useFluctuatingX; }

/**
 * Condition for the factor in the fluctuating x = xQsFactor Q_s
 * e^{+-y}/sqrt(s).
 * \param[in] p The parameters read so far.
 * \return Whether `useFluctuatingX` is 1.
 */
bool fluctuatingX(const Parameters &p) { return p.colorCharge.useFluctuatingX; }

/**
 * Condition for the T^{mu nu} format.
 * \param[in] p The parameters read so far.
 * \return Whether `writeTmunu` is set.
 */
bool tmunuWritten(const Parameters &p) { return p.output.writeTmunu; }

/**
 * Condition for the parameters of the gaussian nucleon model.
 * \param[in] p The parameters read so far.
 * \return Whether `nucleonModel` is `gaussian`.
 */
bool gaussianNucleons(const Parameters &p) {
    return p.subnucleon.nucleonModel == "gaussian";
}

/**
 * Condition for the number of hot spots, which only the hot-spot nucleon
 * model samples (`strings` always has three).
 * \param[in] p The parameters read so far.
 * \return Whether `nucleonModel` is `hotspots`.
 */
bool hotSpotNucleons(const Parameters &p) {
    return p.subnucleon.nucleonModel == "hotspots";
}

/**
 * Condition for the hot-spot parameters, which the hot-spot and the
 * string nucleon models share.
 * \param[in] p The parameters read so far.
 * \return Whether `nucleonModel` is `hotspots` or `strings`.
 */
bool hotSpotSubstructure(const Parameters &p) {
    return hotSpotNucleons(p) || p.subnucleon.nucleonModel == "strings";
}

/**
 * Condition for the width of the Q_s fluctuations.
 * \param[in] p The parameters read so far.
 * \return Whether `smearQs` is set.
 */
bool qsSmeared(const Parameters &p) { return p.subnucleon.smearQs; }

/**
 * Condition for the index of the posterior set.
 * \param[in] p The parameters read so far.
 * \return Whether `posteriorFile` is not `none`.
 */
bool posteriorUsed(const Parameters &p) { return p.posterior.used(); }

/**
 * The setting a condition of a parameter that a posterior file can set
 * stands for, for the message about a column the input does not use.
 * \param[in] condition The condition.
 * \return E.g. `smearQs 1`; "" for another condition.
 */
std::string settingOf(Condition condition) {
    if (condition == qsSmeared) return "smearQs 1";
    if (condition == gaussianNucleons) return "nucleonModel gaussian";
    if (condition == hotSpotNucleons) return "nucleonModel hotspots";
    if (condition == hotSpotSubstructure) {
        return "nucleonModel hotspots or strings";
    }
    if (condition == wsDeformParamsSet) return "useInputWSParams 1";
    if (condition == constantColorCharge) return "useNucleus 0";
    if (condition == fluctuatingX) return "useFluctuatingX 1";
    if (condition == couplingRuns) return "runningCoupling 1";
    if (condition == jimwlkEnabled) return "useJIMWLK 1";
    return "";
}

/// Short name for the table entries.
using P = Parameters;

/**
 * The table of all input parameters. Its order is the order parameters
 * are read in (a condition can only depend on parameters listed before
 * it) and written to usedParameters*.dat.
 * \return The table, built on first use.
 */
const std::vector<ParameterSpec> &parameterTable() {
    static const std::vector<ParameterSpec> table = {
        // general setup
        param(
            "runEvolution", &P::evolution, &EvolutionParameters::runEvolution),
        param("size", &P::lattice, &LatticeParameters::size)
            .check(positive<int>())
            .check(even()),  // the FFTs assume even lattice dimensions
        param("L", &P::lattice, &LatticeParameters::L).check(positive()),
        // divisor of the color-charge density of each layer
        param("Ny", &P::colorCharge, &ColorChargeParameters::Ny)
            .check(positive<int>()),
        // the color-charge densities are divided by g^2
        param("g", &P::coupling, &CouplingParameters::g).check(positive()),
        param("maxTime", &P::evolution, &EvolutionParameters::maxTime)
            .check(nonNegative()),
        // the time step is about a over this number (issue #109)
        param(
            "timeStepsPerLatticeSpacing", &P::evolution,
            &EvolutionParameters::timeStepsPerLatticeSpacing)
            .optional("10")
            .check(positive<int>()),
        param(
            "inverseQsForMaxTime", &P::evolution,
            &EvolutionParameters::inverseQsForMaxTime),

        // posterior parameter sets; before every parameter marked
        // posterior(), since its columns need not be in the input
        param("posteriorFile", &P::posterior, &PosteriorParameters::file)
            .optional("none"),
        param("posteriorSet", &P::posterior, &PosteriorParameters::set)
            .onlyIf(posteriorUsed)
            .check([](const int &v) {
                return v >= -1 ? "" : "must be -1 (random) or a set index >= 0";
            }),

        // random seed
        param("seed", &P::random, &RandomParameters::seed)
            .check(atLeast(-1LL, " (-1 draws a random seed)")),
        param("useSeedList", &P::random, &RandomParameters::useSeedList),

        // collision system and geometry
        param("projectile", &P::collision, &CollisionParameters::projectile),
        param("target", &P::collision, &CollisionParameters::target),
        // without a cross section no nucleons collide, and the impact
        // parameter would be resampled forever
        param("sigmaNN", &P::collision, &CollisionParameters::sigmaNN)
            .posterior()
            .check(positive()),
        param("useNucleus", &P::collision, &CollisionParameters::useNucleus),
        // before the nucleon positions, which are only sampled (and their
        // configuration files read) when the Wilson lines are not read
        param(
            "readInitialWilsonLines", &P::wilsonLines,
            &WilsonLineParameters::readInitialWilsonLines)
            .check(oneOf({0, 1, 2}, " (0: sample, 1: text, 2: binary)")),
        // a separate stream for the color charges (issue #127)
        param("colorChargeSeed", &P::random, &RandomParameters::colorChargeSeed)
            .onlyIf(colorChargesSampled)
            .optional("-1")
            .check(atLeast(-1LL, " (-1: the stream of the run)")),
        param("bMin", &P::collision, &CollisionParameters::bMin)
            .onlyIf(nucleiCollide)
            .check(nonNegative()),
        param("bMax", &P::collision, &CollisionParameters::bMax)
            .onlyIf(nucleiCollide),
        param(
            "sampleBFromLinearDistribution", &P::collision,
            &CollisionParameters::sampleBFromLinearDistribution)
            .onlyIf(nucleiCollide),
        param(
            "rotateReactionPlane", &P::collision,
            &CollisionParameters::rotateReactionPlane)
            .onlyIf(nucleiCollide),
        param("g2muGeV", &P::collision, &CollisionParameters::g2muGeV)
            .onlyIf(constantColorCharge)
            .posterior()
            .check(positive()),
        param(
            "useSmoothNucleus", &P::nucleus,
            &NucleusParameters::useSmoothNucleus),
        param(
            "useFixedNpart", &P::collision, &CollisionParameters::useFixedNpart)
            .onlyIf(nucleiCollide)
            .check(nonNegative<int>()),
        param(
            "nucleiToAverage", &P::collision,
            &CollisionParameters::nucleiToAverage)
            .check(positive<int>()),
        param(
            "gaussianWounding", &P::collision,
            &CollisionParameters::gaussianWounding)
            .onlyIf(nucleiCollide),

        // nucleon positions
        param(
            "nucleonPositionsFromFile", &P::nucleus,
            &NucleusParameters::nucleonPositionsFromFile),
        param(
            "polarizationProjectile", &P::nucleus,
            &NucleusParameters::polarizationProjectile)
            .check(oneOf(
                {0, 1, 2},
                " (0: random orientation, 1: longitudinal, 2: transverse)")),
        param(
            "polarizationTarget", &P::nucleus,
            &NucleusParameters::polarizationTarget)
            .check(oneOf(
                {0, 1, 2},
                " (0: random orientation, 1: longitudinal, 2: transverse)")),
        param(
            "polarizationProjectileJz", &P::nucleus,
            &NucleusParameters::polarizationProjectileJz)
            .onlyIf(projectilePolarized),
        param(
            "polarizationTargetJz", &P::nucleus,
            &NucleusParameters::polarizationTargetJz)
            .onlyIf(targetPolarized),
        param(
            "nuclearConfigurationsPath", &P::nucleus,
            &NucleusParameters::nuclearConfigurationsPath)
            .onlyIf(configurationFilesRead)
            .optional("./nucleusConfigurations"),
        param(
            "lightNucleusOption", &P::nucleus,
            &NucleusParameters::lightNucleusOption)
            .onlyIf(configurationFilesRead),

        // Woods-Saxon deformation
        param(
            "useInputWSParams", &P::nucleus,
            &NucleusParameters::useInputWSParams),
        param("radiusWS", &P::nucleus, &NucleusParameters::radiusWS)
            .onlyIf(wsDeformParamsSet)
            .posterior()
            .check(positive()),
        param("diffusenessWS", &P::nucleus, &NucleusParameters::diffusenessWS)
            .onlyIf(wsDeformParamsSet)
            .posterior()
            .check(positive()),
        param("beta2", &P::nucleus, &NucleusParameters::beta2)
            .onlyIf(wsDeformParamsSet)
            .posterior(),
        param("beta3", &P::nucleus, &NucleusParameters::beta3)
            .onlyIf(wsDeformParamsSet)
            .posterior(),
        param("beta4", &P::nucleus, &NucleusParameters::beta4)
            .onlyIf(wsDeformParamsSet)
            .posterior(),
        param("gamma", &P::nucleus, &NucleusParameters::gamma)
            .onlyIf(wsDeformParamsSet)
            .posterior(),
        param("deltaRnp", &P::nucleus, &NucleusParameters::deltaRnp)
            .onlyIf(wsDeformParamsSet)
            .posterior(),
        param("deltaAnp", &P::nucleus, &NucleusParameters::deltaAnp)
            .onlyIf(wsDeformParamsSet)
            .posterior(),
        // Glauber::findNucleusData applies these regardless of
        // setWSDeformParams
        param("forceDMin", &P::nucleus, &NucleusParameters::forceDMin),
        param("dMin", &P::nucleus, &NucleusParameters::dMin).posterior(),

        // nucleon substructure and color charges
        param(
            "nucleonModel", &P::subnucleon, &SubnucleonParameters::nucleonModel)
            .check(oneOf<std::string>({"gaussian", "hotspots", "strings"})),
        param("m", &P::subnucleon, &SubnucleonParameters::m).posterior(),
        param("BG", &P::subnucleon, &SubnucleonParameters::BG)
            .posterior()
            .check(positive()),
        // gaussian nucleons
        param(
            "protonAnisotropy", &P::subnucleon,
            &SubnucleonParameters::protonAnisotropy)
            .onlyIf(gaussianNucleons)
            .posterior()
            // the thickness is normalized with sqrt(1 + protonAnisotropy)
            .check([](const double &v) {
                return v > -1. ? "" : "must be larger than -1";
            }),
        // hot spots (also the ends of the strings of `strings`)
        // the hot-spot width is 0.09 plus a log-normal number with mean
        // BGq - 0.09, which needs a positive mean
        param("BGq", &P::subnucleon, &SubnucleonParameters::BGq)
            .onlyIf(hotSpotSubstructure)
            .posterior()
            .check([](const double &v) {
                return v > 0.09 ? "" : "must be larger than 0.09 (GeV^-2)";
            }),
        param("BGqVar", &P::subnucleon, &SubnucleonParameters::BGqVar)
            .onlyIf(hotSpotSubstructure)
            .posterior()
            .check(nonNegative()),
        param("dqMin", &P::subnucleon, &SubnucleonParameters::dqMin)
            .onlyIf(hotSpotSubstructure)
            .posterior()
            .check(nonNegative()),
        param("omega", &P::subnucleon, &SubnucleonParameters::omega)
            .onlyIf(hotSpotSubstructure)
            .posterior()
            .check(positive()),
        // the number of hot spots (`strings` always has three)
        param("Nq", &P::subnucleon, &SubnucleonParameters::Nq)
            .onlyIf(hotSpotNucleons)
            .posterior()
            .check([](const double &v) {
                return v >= 1. ? "" : "must be at least 1";
            }),
        param("NqFluc", &P::subnucleon, &SubnucleonParameters::NqFluc)
            .onlyIf(hotSpotNucleons)
            .posterior(),
        param(
            "shiftConstituentQuarkProtonOrigin", &P::subnucleon,
            &SubnucleonParameters::shiftConstituentQuarkProtonOrigin)
            .onlyIf(hotSpotSubstructure),
        // g^2 mu = Q_s / QsMuRatio
        param("QsMuRatio", &P::colorCharge, &ColorChargeParameters::QsMuRatio)
            .posterior()
            .check(positive()),
        param("smearQs", &P::subnucleon, &SubnucleonParameters::smearQs),
        param(
            "smearingWidth", &P::subnucleon,
            &SubnucleonParameters::smearingWidth)
            .onlyIf(qsSmeared)
            .posterior()
            .check(nonNegative()),
        // a negative damping length would amplify the UV
        param("UVDamp", &P::subnucleon, &SubnucleonParameters::UVDamp)
            .posterior()
            .check(nonNegative()),
        param(
            "removeZeroMode", &P::subnucleon,
            &SubnucleonParameters::removeZeroMode)
            .optional("0"),
        param(
            "minimumQs2ST", &P::colorCharge,
            &ColorChargeParameters::minimumQs2ST)
            .posterior()
            .check(nonNegative()),
        param(
            "nucleusQsTableFileName", &P::colorCharge,
            &ColorChargeParameters::nucleusQsTableFileName),

        // rapidity and x
        param("rapidity", &P::colorCharge, &ColorChargeParameters::rapidity),
        param(
            "usePseudoRapidity", &P::colorCharge,
            &ColorChargeParameters::usePseudoRapidity),
        param(
            "jacobianMass", &P::colorCharge,
            &ColorChargeParameters::jacobianMass)
            .onlyIf(pseudoRapidityUsed),
        param(
            "useFluctuatingX", &P::colorCharge,
            &ColorChargeParameters::useFluctuatingX),
        param("sqrtS", &P::collision, &CollisionParameters::sqrtS)
            .onlyIf(sqrtSUsed)
            .check(positive()),
        param(
            "projectileX", &P::colorCharge, &ColorChargeParameters::projectileX)
            .onlyIf(fixedX)
            .check(positive()),
        param("targetX", &P::colorCharge, &ColorChargeParameters::targetX)
            .onlyIf(fixedX)
            .check(positive()),
        param("xQsFactor", &P::colorCharge, &ColorChargeParameters::xQsFactor)
            .onlyIf(fluctuatingX)
            .posterior()
            .check(positive()),

        // running coupling
        param(
            "runningCoupling", &P::coupling,
            &CouplingParameters::runningCoupling),
        param("mu0", &P::coupling, &CouplingParameters::mu0)
            .onlyIf(couplingRuns)
            .posterior(),
        param("c", &P::coupling, &CouplingParameters::c)
            .onlyIf(couplingRuns)
            .posterior()
            .check(positive()),
        param("nFlavors", &P::coupling, &CouplingParameters::nFlavors)
            .optional("3")
            .check(inRange(
                0, 16,
                " (the one-loop beta-function coefficient 11*Nc - "
                "2*nFlavors must be positive)")),
        param("LambdaQCD", &P::coupling, &CouplingParameters::LambdaQCD)
            .onlyIf(couplingRuns)
            .posterior()
            .optional("0.2")
            .check(positive()),
        param("runWithQs", &P::coupling, &CouplingParameters::runWithQs)
            .onlyIf(couplingRuns)
            .check(oneOf({0, 1, 2}, " (0: min, 1: average, 2: max Qs)")),
        param(
            "runningCouplingQsFactor", &P::coupling,
            &CouplingParameters::runningCouplingQsFactor)
            .onlyIf(couplingRuns)
            .posterior(),
        param(
            "runWithLocalQs", &P::coupling, &CouplingParameters::runWithLocalQs)
            .onlyIf(couplingRuns),
        param("runWithKt", &P::coupling, &CouplingParameters::runWithKt)
            .onlyIf(couplingRuns),

        // observables
        param(
            "computeGluonMultiplicity", &P::output,
            &OutputParameters::computeGluonMultiplicity),
        param(
            "computeEccentricities", &P::output,
            &OutputParameters::computeEccentricities),
        param(
            "eccentricityCutoff", &P::output,
            &OutputParameters::eccentricityCutoff)
            .onlyIf(eccentricitiesComputed)
            .optional("0")
            .check(nonNegative()),

        // output
        param("writeHydro", &P::output, &OutputParameters::writeHydro),
        param("writeJazma", &P::output, &OutputParameters::writeJazma),
        param("writeTmunu", &P::output, &OutputParameters::writeTmunu),
        param(
            "writeTmunuBinary", &P::output, &OutputParameters::writeTmunuBinary)
            .onlyIf(tmunuWritten)
            .optional("1"),
        param("outputTimes", &P::output, &OutputParameters::outputTimes)
            .onlyIf(fieldOutputWritten)
            .optional("none")
            .check([](const std::vector<double> &times) {
                for (const double t : times) {
                    if (t <= 0.) return "every time must be positive";
                }
                return "";
            }),
        param(
            "writeHadronSpectrum", &P::output,
            &OutputParameters::writeHadronSpectrum)
            .onlyIf(gluonMultiplicityComputed)
            .optional("0"),
        param(
            "gluonSpectrumEstimator", &P::output,
            &OutputParameters::gluonSpectrumEstimator)
            .onlyIf(gluonMultiplicityComputed)
            .optional("1")
            .check(oneOf(
                {0, 1},
                " (0: 2|E|^2/omega, 1: projection onto the free "
                "Bjorken modes)")),
        param(
            "writeWilsonLineSnapshot", &P::output,
            &OutputParameters::writeWilsonLineSnapshot)
            .optional("0"),
        param("writeNpartList", &P::output, &OutputParameters::writeNpartList)
            .optional("1"),
        param("writeNcollList", &P::output, &OutputParameters::writeNcollList)
            .optional("1"),
        param(
            "writeNgluonEstimators", &P::output,
            &OutputParameters::writeNgluonEstimators)
            .optional("1"),
        param(
            "writeOutputsToHDF5", &P::output,
            &OutputParameters::writeOutputsToHDF5),
        param("LOutput", &P::output, &OutputParameters::LOutput)
            .onlyIf(fieldOutputWritten)
            .check(positive()),
        param("sizeOutput", &P::output, &OutputParameters::sizeOutput)
            .onlyIf(fieldOutputWritten)
            .check(positive<int>()),

        // Wilson lines
        param(
            "writeWilsonLines", &P::wilsonLines,
            &WilsonLineParameters::writeWilsonLines)
            .check(oneOf({0, 1, 2}, " (0: none, 1: text, 2: binary)")),
        param(
            "writeWilsonLineGeometry", &P::wilsonLines,
            &WilsonLineParameters::writeGeometry)
            .onlyIf(wilsonLinesWritten)
            .optional("1"),
        param(
            "wilsonLinePath", &P::wilsonLines,
            &WilsonLineParameters::wilsonLinePath)
            .onlyIf(wilsonLinesWrittenOrRead)
            .optional("./"),
        param("readWilsonLinesX", &P::wilsonLines, &WilsonLineParameters::readX)
            .onlyIf(wilsonLinesRead)
            .optional("0")
            .check(nonNegative()),

        // JIMWLK
        param("useJIMWLK", &P::jimwlk, &JimwlkParameters::enabled),
        param("jimwlkMu0", &P::jimwlk, &JimwlkParameters::mu0)
            .onlyIf(jimwlkEnabled)
            .posterior(),
        param(
            "jimwlkLambdaQCD", &P::jimwlk,
            &JimwlkParameters::LambdaQCD)
            .onlyIf(jimwlkEnabled)
            .posterior()
            .check(positive()),  // in GeV
        param("jimwlkC", &P::jimwlk, &JimwlkParameters::c)
            .onlyIf(jimwlkEnabled)
            .posterior()
            .optional("0.2")
            .check(positive()),
        param("jimwlkMass", &P::jimwlk, &JimwlkParameters::mass)
            .onlyIf(jimwlkEnabled)
            .posterior(),
        // 0 selects the running coupling, > 0 a fixed coupling; a posterior
        // set cannot switch between the two (readInput())
        param("jimwlkAlphaS", &P::jimwlk, &JimwlkParameters::alphaS)
            .onlyIf(jimwlkEnabled)
            .posterior()
            .check(nonNegative()),
        param(
            "jimwlkRunningCoupling", &P::jimwlk,
            &JimwlkParameters::runningCoupling)
            .onlyIf(jimwlkCouplingRuns)
            .optional("1")
            .check(oneOf({1, 2}, " (1: square root, 2: noise)")),
        // divisor of the step count and under a square root in the step
        param("jimwlkDs", &P::jimwlk, &JimwlkParameters::Ds)
            .onlyIf(jimwlkEnabled)
            .check(positive()),
        param("jimwlkInitialX", &P::jimwlk, &JimwlkParameters::initialX)
            .onlyIf(jimwlkEnabled)
            .posterior()
            .check(positive()),
        param(
            "jimwlkSaveSnapshots", &P::jimwlk, &JimwlkParameters::saveSnapshots)
            .onlyIf(jimwlkEnabled),
        param(
            "jimwlkXSnapshotList", &P::jimwlk, &JimwlkParameters::xSnapshotList)
            .onlyIf(saveSnapshotsSet),
    };
    return table;
}

/**
 * Input keys renamed in IP-Glasma 2.0 (issue #32), so an old input file
 * gets told the new name instead of only "unknown parameter".
 * \return Map from each old key to its new key.
 */
const std::map<std::string, std::string> &renamedKeys() {
    static const std::map<std::string, std::string> renamed = {
        {"maxtime", "maxTime"},
        {"Projectile", "projectile"},
        {"Target", "target"},
        {"roots", "sqrtS"},
        {"SigmaNN", "sigmaNN"},
        {"bmin", "bMin"},
        {"bmax", "bMax"},
        {"samplebFromLinearDistribution", "sampleBFromLinearDistribution"},
        {"averageOverThisManyNuclei", "nucleiToAverage"},
        {"polariztionProjectile", "polarizationProjectile"},
        {"polariztionTarget", "polarizationTarget"},
        {"setWSDeformParams", "useInputWSParams"},
        {"R_WS", "radiusWS"},
        {"a_WS", "diffusenessWS"},
        {"dR_np", "deltaRnp"},
        {"da_np", "deltaAnp"},
        {"force_dmin_flag", "forceDMin"},
        {"d_min", "dMin"},
        {"SubNucleonParamSet", "posteriorSet"},
        {"subNucleonParamSet", "posteriorSet"},
        {"UVdamp", "UVDamp"},
        {"QsmuRatio", "QsMuRatio"},
        {"NucleusQsTableFileName", "nucleusQsTableFileName"},
        {"Jacobianm", "jacobianMass"},
        {"useFluctuatingx", "useFluctuatingX"},
        {"xFromThisFactorTimesQs", "xQsFactor"},
        {"muZero", "mu0"},
        {"runWith0Min1Avg2MaxQs", "runWithQs"},
        {"runWithThisFactorTimesQs", "runningCouplingQsFactor"},
        {"runWithkt", "runWithKt"},
        {"writeInitialWilsonLines", "writeWilsonLines"},
    };
    return renamed;
}

/**
 * Input keys of IP-Glasma before 2.0 that have no single new name or were
 * removed, with a hint at what replaced them.
 * \return Map from each old key to the hint.
 */
const std::map<std::string, std::string> &replacedKeys() {
    static const std::map<std::string, std::string> replaced = {
        {"subNucleonParamType",
         "replaced by posteriorFile: none for 0, tables/posterior.csv for 1, "
         "tables/posterior_Nq3.csv for 2 and tables/posterior5020_Nq3.csv "
         "for 4"},
        {"SubNucleonParamType",
         "replaced by posteriorFile: none for 0, tables/posterior.csv for 1, "
         "tables/posterior_Nq3.csv for 2 and tables/posterior5020_Nq3.csv "
         "for 4"},
        {"useConstituentQuarkProton",
         "replaced by nucleonModel: gaussian, or hotspots with Nq hot spots"},
        {"Rapidity",
         "replaced by projectileX and targetX (the x of the nuclei; "
         "Rapidity y set x = 0.01 exp(-y)) and by rapidity, which is "
         "the rapidity of "
         "the spectra"},
        {"rapidityA",
         "replaced by projectileX (the x of the projectile) and rapidity "
         "(the rapidity of the spectra)"},
        {"rapidityB",
         "replaced by targetX (the x of the target) and rapidity (the "
         "rapidity of the spectra)"},
        {"RapidityA", "replaced by projectileX and rapidity"},
        {"RapidityB", "replaced by targetX and rapidity"},
        {"writeOutputs",
         "replaced by writeHydro, writeJazma, writeTmunu, outputTimes, "
         "writeHadronSpectrum and writeWilsonLineSnapshot"},
        {"Nc", "removed: the code is SU(3) only"},
        {"rmax", "removed: color charges cover the whole lattice"},
        {"dtau", "removed: the time step follows from maxTime, L and size"},
        {"tDistNu", "removed: it had no effect"},
        {"useFatTails", "removed: it had no effect"},
        {"writeEvolution", "removed: it had no effect"},
        {"readMultFromFile",
         "removed: no version of the code writes the files it read"},
        {"useRandomSeed", "replaced by seed -1, which draws a random seed"},
        {"useTimeForSeed", "replaced by seed -1, which draws a random seed"},
        {"mode",
         "replaced by runEvolution: 1 for mode 1, 0 for any other mode"},
        {"useGaussian", "removed: useNucleus 0 always uses a constant g2muGeV"},
        {"g2mu",
         "replaced by g2muGeV, in GeV instead of lattice units: g2muGeV = "
         "g2mu * 0.19733 / (L/size)"},
        {"etaSizeOutput",
         "removed: the fields are boost invariant, so the outputs hold one "
         "eta slice"},
        {"dEtaOutput",
         "removed: the fields are boost invariant, so the outputs hold one "
         "eta slice"},
        {"detaOutput",
         "removed: the fields are boost invariant, so the outputs hold one "
         "eta slice"},
    };
    return replaced;
}

/**
 * Looks up the old names of a parameter.
 * \param[in] name Current input-file key.
 * \return Its old keys, empty if it was not renamed.
 */
std::vector<std::string> oldNamesOf(const std::string &name) {
    std::vector<std::string> oldNames;
    for (const auto &[oldName, newName] : renamedKeys()) {
        if (newName == name) oldNames.push_back(oldName);
    }
    return oldNames;
}

/**
 * Case-insensitive Levenshtein distance between two keys.
 * \param[in] a First key.
 * \param[in] b Second key.
 * \return Number of single-character insertions, deletions and
 * substitutions turning \p a into \p b.
 */
std::size_t editDistance(const std::string &a, const std::string &b) {
    std::vector<std::size_t> row(b.size() + 1);
    for (std::size_t j = 0; j <= b.size(); j++) row[j] = j;
    for (std::size_t i = 1; i <= a.size(); i++) {
        std::size_t diagonal = row[0];
        row[0] = i;
        for (std::size_t j = 1; j <= b.size(); j++) {
            const std::size_t above = row[j];
            const bool same =
                std::tolower(static_cast<unsigned char>(a[i - 1]))
                == std::tolower(static_cast<unsigned char>(b[j - 1]));
            row[j] = std::min({row[j] + 1, row[j - 1] + 1, diagonal + !same});
            diagonal = above;
        }
    }
    return row[b.size()];
}

/**
 * Finds a suggestion for an unknown key.
 * \param[in] key The unknown key.
 * \return The known parameter name closest to \p key, or "" if none is
 * close enough to be a likely typo.
 */
std::string closestName(const std::string &key) {
    std::string best;
    std::size_t bestDistance = std::numeric_limits<std::size_t>::max();
    for (const ParameterSpec &spec : parameterTable()) {
        const std::size_t distance = editDistance(key, spec.name);
        // a suggestion needs some characters in common with both names, so
        // short keys don't match the one-letter parameters (L, g, m, c)
        if (distance >= std::min(key.size(), spec.name.size())) continue;
        if (distance < bestDistance) {
            bestDistance = distance;
            best = spec.name;
        }
    }
    return (bestDistance <= std::max<std::size_t>(2, key.size() / 4)) ? best
                                                                      : "";
}

/**
 * Removes the spaces, tabs and carriage returns around a text.
 * \param[in] text The text.
 * \return \p text without them.
 */
std::string trimmed(const std::string &text) {
    const std::size_t first = text.find_first_not_of(" \t\r");
    if (first == std::string::npos) return "";
    const std::size_t last = text.find_last_not_of(" \t\r");
    return text.substr(first, last - first + 1);
}

/**
 * Splits a line of a posterior file into its comma-separated values.
 * \param[in] line The line.
 * \return Its values, without the spaces around them.
 */
std::vector<std::string> commaSeparated(const std::string &line) {
    std::vector<std::string> values;
    std::size_t start = 0;
    while (true) {
        const std::size_t comma = line.find(',', start);
        values.push_back(trimmed(line.substr(start, comma - start)));
        if (comma == std::string::npos) return values;
        start = comma + 1;
    }
}

/// Number of sets of a posterior file whose problems are reported; the
/// others are only counted.
constexpr int reportedPosteriorSets = 10;

}  // namespace

std::vector<std::string> Parameters::loadPosteriorTable() {
    posteriorTable_.reset();
    const std::string &fileName = posterior.file;
    std::ifstream file(fileName);
    if (!file.is_open()) {
        return {"cannot open the posterior file " + fileName};
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        // skip a UTF-8 byte order mark written by some editors
        if (lines.empty() && line.compare(0, 3, "\xEF\xBB\xBF") == 0) {
            line.erase(0, 3);
        }
        lines.push_back(trimmed(line));
    }
    // The columns are named on the first line that does not start with #,
    // unless that line is already a set (all parameters a posterior can set
    // are numbers); then on the last # line before it. The # lines before
    // the names are comments.
    std::size_t first = 0;
    while (first < lines.size()
           && (lines[first].empty() || lines[first][0] == '#')) {
        first++;
    }
    std::size_t header = lines.size();  // none
    double number = 0.;
    if (first < lines.size()
        && !parseValue(commaSeparated(lines[first])[0], number)) {
        header = first;
    } else {
        for (std::size_t i = first; i-- > 0;) {
            if (!lines[i].empty()) {
                header = i;
                break;
            }
        }
    }
    if (header == lines.size()) {
        return {fileName + " has no line naming the columns"};
    }

    const std::vector<ParameterSpec> &specs = parameterTable();
    auto table = std::make_shared<PosteriorTable>();
    std::vector<std::string> errors;
    // units in brackets after the names
    std::string content = lines[header];
    if (content[0] == '#') content.erase(0, 1);
    const std::vector<std::string> names = commaSeparated(content);
    std::string where = fileName + ":" + std::to_string(header + 1) + ": ";
    for (const std::string &column : names) {
        const std::string name = trimmed(column.substr(0, column.find('[')));
        const auto spec = std::find_if(
            specs.begin(), specs.end(),
            [&name](const ParameterSpec &s) { return s.name == name; });
        if (name.empty()) {
            errors.push_back(where + "a column has no name");
        } else if (spec == specs.end()) {
            std::string message = where + "unknown parameter " + name;
            const std::string suggestion = closestName(name);
            if (!suggestion.empty())
                message += " (did you mean " + suggestion + "?)";
            errors.push_back(message);
        } else if (!spec->posterior) {
            errors.push_back(
                where + name
                + " cannot be set by a posterior file, only the physical "
                  "parameters can");
        } else if (
            std::find(table->columns.begin(), table->columns.end(), name)
            != table->columns.end()) {
            errors.push_back(where + "column " + name + " appears twice");
        } else {
            table->columns.push_back(name);
            table->specs.push_back(
                static_cast<std::size_t>(spec - specs.begin()));
        }
    }

    int malformedSets = 0;
    for (std::size_t i = header + 1; i < lines.size(); i++) {
        if (lines[i].empty() || lines[i][0] == '#') continue;
        std::vector<std::string> values = commaSeparated(lines[i]);
        if (values.size() != names.size()) {
            if (++malformedSets <= reportedPosteriorSets) {
                errors.push_back(
                    fileName + ":" + std::to_string(i + 1) + ": "
                    + std::to_string(values.size()) + " values, but line "
                    + std::to_string(header + 1) + " names "
                    + std::to_string(names.size()) + " columns");
            }
            continue;
        }
        table->sets.push_back(std::move(values));
        table->lines.push_back(static_cast<int>(i + 1));
    }
    if (malformedSets > reportedPosteriorSets) {
        errors.push_back(
            fileName + ": "
            + std::to_string(malformedSets - reportedPosteriorSets)
            + " more sets with the wrong number of values");
    }
    if (errors.empty() && table->sets.empty()) {
        errors.push_back(fileName + " contains no parameter sets");
    }
    if (errors.empty()) posteriorTable_ = std::move(table);
    return errors;
}

std::string Parameters::applyPosteriorValues(std::size_t set) {
    const PosteriorTable &table = *posteriorTable_;
    for (std::size_t column = 0; column < table.columns.size(); column++) {
        const std::string problem = parameterTable()[table.specs[column]].apply(
            *this, table.sets[set][column]);
        if (!problem.empty()) return table.columns[column] + " " + problem;
    }
    return "";
}

void Parameters::applyPosteriorSet(int set) {
    if (!posteriorTable_) return;
    set %= static_cast<int>(posteriorTable_->sets.size());
    event.posteriorSet = set;
    messager_ << "[Parameters::applyPosteriorSet]: Using parameter set " << set
              << " of the posterior file " << posterior.file << ".";
    messager_.flush("info");
    // readInput() checked the values of every set
    applyPosteriorValues(static_cast<std::size_t>(set));
}

std::vector<std::string> Parameters::readInput(const InputFile &input) {
    std::vector<std::string> errors = input.errors();
    if (!input.isOpen()) return errors;
    const std::string &source = input.sourceName();
    inputWarnings_.clear();
    posteriorTable_.reset();

    // old keys already reported as "renamed", so they are not reported
    // again as unknown below
    std::vector<std::string> reportedOldKeys;
    // the columns of the posterior file whose parameter this input reads
    std::vector<bool> columnRead;
    const std::vector<ParameterSpec> &table = parameterTable();
    for (std::size_t index = 0; index < table.size(); index++) {
        const ParameterSpec &spec = table[index];
        if (spec.condition && !spec.condition(*this)) continue;
        // the column of the posterior file that sets this parameter, if any
        int column = -1;
        if (posteriorTable_) {
            const std::vector<std::size_t> &specs = posteriorTable_->specs;
            const auto found = std::find(specs.begin(), specs.end(), index);
            if (found != specs.end())
                column = static_cast<int>(found - specs.begin());
        }
        const InputFile::Entry *entry = input.find(spec.name);
        std::string oldName;
        const InputFile::Entry *oldEntry = nullptr;
        for (const std::string &name : oldNamesOf(spec.name)) {
            oldEntry = input.find(name);
            if (oldEntry) {
                oldName = name;
                break;
            }
        }
        std::string problem;
        std::string where = source + ": ";
        if (entry) {
            where = source + ":" + std::to_string(entry->line) + ": ";
            problem = spec.apply(*this, entry->value);
            if (problem.empty() && column >= 0) {
                inputWarnings_.push_back(
                    where + spec.name + " " + entry->value
                    + " is replaced by the posterior file " + posterior.file);
            }
        } else if (oldEntry) {
            errors.push_back(
                source + ":" + std::to_string(oldEntry->line) + ": " + oldName
                + " was renamed to " + spec.name);
            reportedOldKeys.push_back(oldName);
            // read, under its old name
            if (column >= 0) columnRead[column] = true;
            continue;
        } else if (column >= 0) {
            // set from the posterior file below
        } else if (!spec.defaultValue.empty()) {
            problem = spec.apply(*this, spec.defaultValue);
        } else {
            problem = "is required but not given";
        }
        if (!problem.empty())
            errors.push_back(where + spec.name + " " + problem);
        if (column >= 0) {
            // the first set, so the parameters read after this one see its
            // value, as they will see that of every set (checked below)
            columnRead[column] = true;
            const std::string valueProblem = spec.apply(
                *this, posteriorTable_->sets[0][column]);
            if (!valueProblem.empty()) {
                errors.push_back(
                    posterior.file + ":"
                    + std::to_string(posteriorTable_->lines[0]) + ": "
                    + spec.name + " " + valueProblem);
            }
        }
        // the parameters of the columns are read after the posterior file
        if (spec.name == "posteriorFile" && posterior.used()) {
            for (const std::string &error : loadPosteriorTable())
                errors.push_back(error);
            if (posteriorTable_)
                columnRead.assign(posteriorTable_->columns.size(), false);
        }
    }
    if (posteriorTable_) {
        for (std::size_t column = 0; column < columnRead.size(); column++) {
            if (columnRead[column]) continue;
            const std::string &name = posteriorTable_->columns[column];
            const std::string setting =
                settingOf(table[posteriorTable_->specs[column]].condition);
            errors.push_back(
                posterior.file + ": column " + name
                + (setting.empty()
                       ? " sets a parameter that this input does not use"
                       : " is only used with " + setting
                             + ", which this input does not set"));
        }
    }

    for (const auto &[key, entry] : input.entries()) {
        const bool known = std::any_of(
            table.begin(), table.end(),
            [&key](const ParameterSpec &spec) { return spec.name == key; });
        if (known) continue;
        if (std::find(reportedOldKeys.begin(), reportedOldKeys.end(), key)
            != reportedOldKeys.end()) {
            continue;
        }
        std::string message = source + ":" + std::to_string(entry.line)
                              + ": unknown parameter " + key;
        const auto renamed = renamedKeys().find(key);
        const auto replaced = replacedKeys().find(key);
        if (renamed != renamedKeys().end()) {
            message += " (renamed to " + renamed->second + ")";
        } else if (replaced != replacedKeys().end()) {
            message += " (" + replaced->second + ")";
        } else {
            const std::string suggestion = closestName(key);
            if (!suggestion.empty())
                message += " (did you mean " + suggestion + "?)";
        }
        errors.push_back(message);
    }
    if (!errors.empty()) return errors;

    // Every set must pass the checks of its values and use the same
    // parameters as the first: a set cannot switch between models (e.g.
    // jimwlkAlphaS between a fixed and the running coupling).
    if (posteriorTable_) {
        Parameters withSet(*this);
        int failingSets = 0;
        for (std::size_t set = 1; set < posteriorTable_->sets.size(); set++) {
            std::string problem = withSet.applyPosteriorValues(set);
            for (std::size_t index = 0; index < table.size() && problem.empty();
                 index++) {
                const Condition condition = table[index].condition;
                if (condition && condition(withSet) != condition(*this)) {
                    problem = "the set changes whether " + table[index].name
                              + " is used";
                }
            }
            if (problem.empty()) continue;
            if (++failingSets <= reportedPosteriorSets) {
                errors.push_back(
                    posterior.file + ":"
                    + std::to_string(posteriorTable_->lines[set]) + ": "
                    + problem);
            }
        }
        if (failingSets > reportedPosteriorSets) {
            errors.push_back(
                posterior.file + ": "
                + std::to_string(failingSets - reportedPosteriorSets)
                + " more sets with problems");
        }
        if (!errors.empty()) return errors;
    }

    // Derived settings.
    // dtau (in lattice units) is ~1/timeStepsPerLatticeSpacing, adjusted so
    // maxTime is a whole number of steps. With fewer than one step (e.g.
    // maxTime 0 to only produce Wilson lines) use 1/timeStepsPerLatticeSpacing
    // instead of dividing by zero.
    const double latticeSpacing = lattice.L / static_cast<double>(lattice.size);
    const double stepsPerSpacing = evolution.timeStepsPerLatticeSpacing;
    const double steps = stepsPerSpacing * evolution.maxTime / latticeSpacing;
    if (steps > 1e8) {
        errors.push_back(
            source + ": maxTime " + std::to_string(evolution.maxTime)
            + " needs more than 1e8 time steps on this lattice");
        return errors;
    }
    const int timeSteps = static_cast<int>(steps);
    run.dtau = (timeSteps > 0)
                   ? evolution.maxTime / (timeSteps * latticeSpacing)
                   : 1. / stepsPerSpacing;
    // polarized nuclei are only available as configuration files
    if (nucleus.polarizationProjectile != 0
        || nucleus.polarizationTarget != 0) {
        nucleus.nucleonPositionsFromFile = true;
    }
    return errors;
}

void Parameters::writeInputParameters(std::ostream &out) const {
    for (const ParameterSpec &spec : parameterTable()) {
        if (spec.condition && !spec.condition(*this)) continue;
        out << spec.name << " ";
        // the set the event used, whose values are those written
        if (spec.name == "posteriorSet" && event.posteriorSet >= 0) {
            out << event.posteriorSet;
        } else {
            spec.write(*this, out);
        }
        out << "\n";
    }
}

// ParameterTable.cpp is part of the IP-Glasma solver.
//
// The single list of every input-file parameter: its key, how it is stored
// in Parameters, its default (if it is optional) and its validity checks.
// Parameters::readInput() and Parameters::writeInputParameters() are both
// driven by this table, so adding a parameter only takes one entry here
// (plus its field in one of the groups in Parameters.h).

#include <algorithm>
#include <cctype>
#include <charconv>
#include <functional>
#include <initializer_list>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "InputFile.h"
#include "Parameters.h"

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
 * Writes \p values comma-separated, in input-file syntax.
 * \param[out] out Stream to write to.
 * \param[in] values Values to write.
 */
void writeValue(std::ostream &out, const std::vector<double> &values) {
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
    if constexpr (std::is_same_v<T, unsigned long long>)
        return "a non-negative integer";
    if constexpr (std::is_same_v<T, double>) return "a number";
    if constexpr (std::is_same_v<T, std::vector<double>>)
        return "a comma-separated list of numbers";
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
 * Condition for the JIMWLK snapshot list.
 * \param[in] p The parameters read so far.
 * \return Whether `jimwlkSaveSnapshots` is set.
 */
bool saveSnapshotsSet(const Parameters &p) { return p.jimwlk.saveSnapshots; }

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
        param("mode", &P::evolution, &EvolutionParameters::mode),
        param("size", &P::lattice, &LatticeParameters::size)
            .check(positive<int>())
            .check(even()),  // the FFTs assume even lattice dimensions
        param("L", &P::lattice, &LatticeParameters::L).check(positive()),
        param("Ny", &P::colorCharge, &ColorChargeParameters::Ny),
        param("sqrtS", &P::collision, &CollisionParameters::sqrtS),
        param("g", &P::coupling, &CouplingParameters::g),
        param("g2mu", &P::collision, &CollisionParameters::g2mu),
        param("maxTime", &P::evolution, &EvolutionParameters::maxTime)
            .check(nonNegative()),
        param(
            "inverseQsForMaxTime", &P::evolution,
            &EvolutionParameters::inverseQsForMaxTime),

        // random seed
        param("seed", &P::random, &RandomParameters::seed),
        param("useSeedList", &P::random, &RandomParameters::useSeedList),
        param("useTimeForSeed", &P::random, &RandomParameters::useTimeForSeed),

        // collision system and geometry
        param("projectile", &P::collision, &CollisionParameters::projectile),
        param("target", &P::collision, &CollisionParameters::target),
        param("sigmaNN", &P::collision, &CollisionParameters::sigmaNN),
        param("bMin", &P::collision, &CollisionParameters::bMin),
        param("bMax", &P::collision, &CollisionParameters::bMax),
        param(
            "sampleBFromLinearDistribution", &P::collision,
            &CollisionParameters::sampleBFromLinearDistribution),
        param(
            "rotateReactionPlane", &P::collision,
            &CollisionParameters::rotateReactionPlane),
        param("useNucleus", &P::collision, &CollisionParameters::useNucleus),
        param("useGaussian", &P::collision, &CollisionParameters::useGaussian),
        param(
            "useSmoothNucleus", &P::nucleus,
            &NucleusParameters::useSmoothNucleus),
        param(
            "useFixedNpart", &P::collision,
            &CollisionParameters::useFixedNpart),
        param(
            "nucleiToAverage", &P::collision,
            &CollisionParameters::nucleiToAverage)
            .check(positive<int>()),
        param(
            "gaussianWounding", &P::collision,
            &CollisionParameters::gaussianWounding),

        // nucleon positions
        param(
            "nucleonPositionsFromFile", &P::nucleus,
            &NucleusParameters::nucleonPositionsFromFile),
        param(
            "nuclearConfigurationsPath", &P::nucleus,
            &NucleusParameters::nuclearConfigurationsPath)
            .optional("./nucleusConfigurations"),
        param(
            "lightNucleusOption", &P::nucleus,
            &NucleusParameters::lightNucleusOption),
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
            &NucleusParameters::polarizationProjectileJz),
        param(
            "polarizationTargetJz", &P::nucleus,
            &NucleusParameters::polarizationTargetJz),

        // Woods-Saxon deformation
        param(
            "useInputWSParams", &P::nucleus,
            &NucleusParameters::useInputWSParams),
        param("radiusWS", &P::nucleus, &NucleusParameters::radiusWS)
            .onlyIf(wsDeformParamsSet),
        param("diffusenessWS", &P::nucleus, &NucleusParameters::diffusenessWS)
            .onlyIf(wsDeformParamsSet),
        param("beta2", &P::nucleus, &NucleusParameters::beta2)
            .onlyIf(wsDeformParamsSet),
        param("beta3", &P::nucleus, &NucleusParameters::beta3)
            .onlyIf(wsDeformParamsSet),
        param("beta4", &P::nucleus, &NucleusParameters::beta4)
            .onlyIf(wsDeformParamsSet),
        param("gamma", &P::nucleus, &NucleusParameters::gamma)
            .onlyIf(wsDeformParamsSet),
        param("deltaRnp", &P::nucleus, &NucleusParameters::deltaRnp)
            .onlyIf(wsDeformParamsSet),
        param("deltaAnp", &P::nucleus, &NucleusParameters::deltaAnp)
            .onlyIf(wsDeformParamsSet),
        // Glauber::findNucleusData applies these regardless of
        // setWSDeformParams
        param("forceDMin", &P::nucleus, &NucleusParameters::forceDMin),
        param("dMin", &P::nucleus, &NucleusParameters::dMin),

        // nucleon substructure and color charges
        param("m", &P::subnucleon, &SubnucleonParameters::m),
        param("BG", &P::subnucleon, &SubnucleonParameters::BG),
        param("BGq", &P::subnucleon, &SubnucleonParameters::BGq),
        param("BGqVar", &P::subnucleon, &SubnucleonParameters::BGqVar),
        param("dqMin", &P::subnucleon, &SubnucleonParameters::dqMin),
        param("omega", &P::subnucleon, &SubnucleonParameters::omega)
            .check(positive()),
        param("Nq", &P::subnucleon, &SubnucleonParameters::Nq)
            .check(nonNegative()),
        param("NqFluc", &P::subnucleon, &SubnucleonParameters::NqFluc),
        param(
            "shiftConstituentQuarkProtonOrigin", &P::subnucleon,
            &SubnucleonParameters::shiftConstituentQuarkProtonOrigin),
        param(
            "protonAnisotropy", &P::subnucleon,
            &SubnucleonParameters::protonAnisotropy),
        param(
            "subNucleonParamType", &P::subnucleon,
            &SubnucleonParameters::subNucleonParamType)
            .check(oneOf({0, 1, 2, 4}, " (0: use the input values)")),
        param(
            "subNucleonParamSet", &P::subnucleon,
            &SubnucleonParameters::subNucleonParamSet),
        param("QsMuRatio", &P::colorCharge, &ColorChargeParameters::QsMuRatio),
        param("smearQs", &P::subnucleon, &SubnucleonParameters::smearQs),
        param(
            "smearingWidth", &P::subnucleon,
            &SubnucleonParameters::smearingWidth),
        param("UVDamp", &P::subnucleon, &SubnucleonParameters::UVDamp),
        param(
            "minimumQs2ST", &P::colorCharge,
            &ColorChargeParameters::minimumQs2ST),
        param(
            "nucleusQsTableFileName", &P::colorCharge,
            &ColorChargeParameters::nucleusQsTableFileName),

        // rapidity and x
        param("rapidityA", &P::colorCharge, &ColorChargeParameters::rapidityA),
        param("rapidityB", &P::colorCharge, &ColorChargeParameters::rapidityB),
        param(
            "usePseudoRapidity", &P::colorCharge,
            &ColorChargeParameters::usePseudoRapidity),
        param(
            "jacobianMass", &P::colorCharge,
            &ColorChargeParameters::jacobianMass),
        param(
            "useFluctuatingX", &P::colorCharge,
            &ColorChargeParameters::useFluctuatingX),
        param("xQsFactor", &P::colorCharge, &ColorChargeParameters::xQsFactor),

        // running coupling
        param(
            "runningCoupling", &P::coupling,
            &CouplingParameters::runningCoupling),
        param("mu0", &P::coupling, &CouplingParameters::mu0),
        param("c", &P::coupling, &CouplingParameters::c).check(positive()),
        param("nFlavors", &P::coupling, &CouplingParameters::nFlavors)
            .optional("3")
            .check(inRange(
                0, 16,
                " (the one-loop beta-function coefficient 11*Nc - "
                "2*nFlavors must be positive)")),
        param("LambdaQCD", &P::coupling, &CouplingParameters::LambdaQCD)
            .optional("0.2")
            .check(positive()),
        param("runWithQs", &P::coupling, &CouplingParameters::runWithQs)
            .check(oneOf({0, 1, 2}, " (0: min, 1: average, 2: max Qs)")),
        param(
            "runningCouplingQsFactor", &P::coupling,
            &CouplingParameters::runningCouplingQsFactor),
        param(
            "runWithLocalQs", &P::coupling,
            &CouplingParameters::runWithLocalQs),
        param("runWithKt", &P::coupling, &CouplingParameters::runWithKt),

        // observables
        param(
            "computeGluonMultiplicity", &P::output,
            &OutputParameters::computeGluonMultiplicity),
        param(
            "readMultFromFile", &P::output,
            &OutputParameters::readMultFromFile),

        // output
        param("writeOutputs", &P::output, &OutputParameters::writeOutputs),
        param(
            "writeEpsilonUHydro", &P::output,
            &OutputParameters::writeEpsilonUHydro)
            .optional("1"),
        param(
            "writeTmunuBinary", &P::output, &OutputParameters::writeTmunuBinary)
            .optional("1"),
        param(
            "writeOutputsToHDF5", &P::output,
            &OutputParameters::writeOutputsToHDF5),
        param("LOutput", &P::output, &OutputParameters::LOutput),
        param("sizeOutput", &P::output, &OutputParameters::sizeOutput),
        param("etaSizeOutput", &P::output, &OutputParameters::etaSizeOutput),
        param("dEtaOutput", &P::output, &OutputParameters::dEtaOutput),

        // Wilson lines
        param(
            "writeWilsonLines", &P::wilsonLines,
            &WilsonLineParameters::writeWilsonLines)
            .check(oneOf({0, 1, 2}, " (0: none, 1: text, 2: binary)")),
        param(
            "wilsonLinePath", &P::wilsonLines,
            &WilsonLineParameters::wilsonLinePath)
            .optional("./"),
        param(
            "readInitialWilsonLines", &P::wilsonLines,
            &WilsonLineParameters::readInitialWilsonLines)
            .check(oneOf({0, 1, 2}, " (0: sample, 1: text, 2: binary)")),

        // JIMWLK
        param("useJIMWLK", &P::jimwlk, &JimwlkParameters::enabled),
        param("jimwlkMu0", &P::jimwlk, &JimwlkParameters::mu0),
        param(
            "jimwlkLambdaQCD", &P::jimwlk,
            &JimwlkParameters::LambdaQCD)
            .check(positive()),  // in GeV
        param("jimwlkC", &P::jimwlk, &JimwlkParameters::c)
            .optional("0.2")
            .check(positive()),
        param("jimwlkMass", &P::jimwlk, &JimwlkParameters::mass),
        // 0 selects the running coupling, > 0 a fixed coupling
        param("jimwlkAlphaS", &P::jimwlk, &JimwlkParameters::alphaS)
            .check(nonNegative()),
        param("jimwlkDs", &P::jimwlk, &JimwlkParameters::Ds),
        param("jimwlkInitialX", &P::jimwlk, &JimwlkParameters::initialX),
        param("jimwlkXProjectile", &P::jimwlk, &JimwlkParameters::xProjectile),
        param("jimwlkXTarget", &P::jimwlk, &JimwlkParameters::xTarget),
        param(
            "jimwlkSaveSnapshots", &P::jimwlk,
            &JimwlkParameters::saveSnapshots),
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
        {"useConstituentQuarkProton", "Nq"},
        {"SubNucleonParamType", "subNucleonParamType"},
        {"SubNucleonParamSet", "subNucleonParamSet"},
        {"UVdamp", "UVDamp"},
        {"QsmuRatio", "QsMuRatio"},
        {"NucleusQsTableFileName", "nucleusQsTableFileName"},
        {"RapidityA", "rapidityA"},
        {"RapidityB", "rapidityB"},
        {"Jacobianm", "jacobianMass"},
        {"useFluctuatingx", "useFluctuatingX"},
        {"xFromThisFactorTimesQs", "xQsFactor"},
        {"muZero", "mu0"},
        {"runWith0Min1Avg2MaxQs", "runWithQs"},
        {"runWithThisFactorTimesQs", "runningCouplingQsFactor"},
        {"runWithkt", "runWithKt"},
        {"detaOutput", "dEtaOutput"},
        {"mu0_jimwlk", "jimwlkMu0"},
        {"Lambda_QCD_jimwlk", "jimwlkLambdaQCD"},
        {"c_jimwlk", "jimwlkC"},
        {"m_jimwlk", "jimwlkMass"},
        {"alphas_jimwlk", "jimwlkAlphaS"},
        {"Ds_jimwlk", "jimwlkDs"},
        {"jimwlk_ic_x", "jimwlkInitialX"},
        {"x_projectile_jimwlk", "jimwlkXProjectile"},
        {"x_target_jimwlk", "jimwlkXTarget"},
        {"saveSnapshots", "jimwlkSaveSnapshots"},
        {"xSnapshotList", "jimwlkXSnapshotList"},
    };
    return renamed;
}

/**
 * Looks up the pre-2.0 name of a parameter.
 * \param[in] name Current input-file key.
 * \return Its old key, or "" if it was not renamed.
 */
std::string oldNameOf(const std::string &name) {
    for (const auto &[oldName, newName] : renamedKeys()) {
        if (newName == name) return oldName;
    }
    return "";
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

}  // namespace

std::vector<std::string> Parameters::readInput(const InputFile &input) {
    std::vector<std::string> errors = input.errors();
    if (!input.isOpen()) return errors;
    const std::string &source = input.sourceName();

    // old keys already reported as "renamed", so they are not reported
    // again as unknown below
    std::vector<std::string> reportedOldKeys;
    for (const ParameterSpec &spec : parameterTable()) {
        if (spec.condition && !spec.condition(*this)) continue;
        const InputFile::Entry *entry = input.find(spec.name);
        const std::string oldName = oldNameOf(spec.name);
        const InputFile::Entry *oldEntry =
            oldName.empty() ? nullptr : input.find(oldName);
        std::string problem;
        std::string where = source + ": ";
        if (entry) {
            where = source + ":" + std::to_string(entry->line) + ": ";
            problem = spec.apply(*this, entry->value);
        } else if (oldEntry) {
            errors.push_back(
                source + ":" + std::to_string(oldEntry->line) + ": " + oldName
                + " was renamed to " + spec.name);
            reportedOldKeys.push_back(oldName);
            continue;
        } else if (!spec.defaultValue.empty()) {
            problem = spec.apply(*this, spec.defaultValue);
        } else {
            problem = "is required but not given";
        }
        if (!problem.empty())
            errors.push_back(where + spec.name + " " + problem);
    }

    for (const auto &[key, entry] : input.entries()) {
        const auto &table = parameterTable();
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
        if (renamed != renamedKeys().end()) {
            message += " (renamed to " + renamed->second + ")";
        } else {
            const std::string suggestion = closestName(key);
            if (!suggestion.empty())
                message += " (did you mean " + suggestion + "?)";
        }
        errors.push_back(message);
    }
    if (!errors.empty()) return errors;

    // Derived settings.
    // dtau (in lattice units) is ~0.1, adjusted so maxtime is a whole number
    // of steps. With fewer than one step (e.g. maxTime 0 to only produce
    // Wilson lines) use 0.1 instead of dividing by zero.
    const double latticeSpacing = lattice.L / static_cast<double>(lattice.size);
    const double steps = 10 * evolution.maxTime / latticeSpacing;
    if (steps > 1e8) {
        errors.push_back(
            source + ": maxTime " + std::to_string(evolution.maxTime)
            + " needs more than 1e8 time steps on this lattice");
        return errors;
    }
    const int timeSteps = static_cast<int>(steps);
    run.dtau = (timeSteps > 0)
                   ? evolution.maxTime / (timeSteps * latticeSpacing)
                   : 0.1;
    // polarized nuclei are only available as configuration files
    if (nucleus.polarizationProjectile != 0
        || nucleus.polarizationTarget != 0) {
        nucleus.nucleonPositionsFromFile = true;
    }
    subnucleon.NqBase = subnucleon.Nq;
    if (subnucleon.subNucleonParamType > 0) {
        const std::string problem =
            loadPosteriorParameterSets(subnucleon.subNucleonParamType);
        if (!problem.empty()) errors.push_back(problem);
    }
    return errors;
}

void Parameters::writeInputParameters(std::ostream &out) const {
    for (const ParameterSpec &spec : parameterTable()) {
        if (spec.condition && !spec.condition(*this)) continue;
        out << spec.name << " ";
        spec.write(*this, out);
        out << "\n";
    }
}

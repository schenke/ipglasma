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
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "InputFile.h"
#include "Parameters.h"

namespace {

using Condition = bool (*)(const Parameters &);

/// One input parameter, with its type erased.
struct ParameterSpec {
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

void writeValue(std::ostream &out, double value) {
    // shortest text that reads back to exactly the same double
    char buffer[64];
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
    out.write(buffer, result.ptr - buffer);
}
void writeValue(std::ostream &out, const std::vector<double> &values) {
    for (std::size_t i = 0; i < values.size(); i++) {
        if (i > 0) out << ",";
        writeValue(out, values[i]);
    }
}
void writeValue(std::ostream &out, bool value) { out << (value ? 1 : 0); }
template <typename T>
void writeValue(std::ostream &out, const T &value) {
    out << value;
}

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

/// Typed builder for a ParameterSpec, so the table reads
/// `param("size", &P::lattice, &LatticeParameters::size).check(even())`.
template <typename T>
class Param {
  public:
    template <typename Group>
    Param(const char *name, Group Parameters::*group, T Group::*field)
        : name_(name),
          set_([group, field](Parameters &p, const T &v) {
              (p.*group).*field = v;
          }),
          write_([group, field](const Parameters &p, std::ostream &out) {
              writeValue(out, (p.*group).*field);
          }) {}

    /// Makes the parameter optional, with this default (as input text).
    Param &optional(const char *defaultValue) {
        defaultValue_ = defaultValue;
        return *this;
    }
    /// Only reads the parameter when \p condition holds.
    Param &onlyIf(Condition condition) {
        condition_ = condition;
        return *this;
    }
    Param &check(Check<T> check) {
        checks_.push_back(std::move(check));
        return *this;
    }

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
    std::string name_;
    std::string defaultValue_;
    Condition condition_ = nullptr;
    std::function<void(Parameters &, const T &)> set_;
    std::function<void(const Parameters &, std::ostream &)> write_;
    std::vector<Check<T>> checks_;
};

/// One input parameter, stored in the field `param.*group.*field`.
template <typename Group, typename T>
Param<T> param(const char *name, Group Parameters::*group, T Group::*field) {
    return Param<T>(name, group, field);
}

// ---- checks ----

template <typename T = double>
Check<T> positive() {
    return [](const T &v) { return v > 0 ? "" : "must be positive"; };
}

Check<int> even() {
    return [](const int &v) { return v % 2 == 0 ? "" : "must be even"; };
}

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

template <typename T>
Check<T> inRange(T low, T high, const char *note = "") {
    return [low, high, note = std::string(note)](const T &v) {
        if (v >= low && v <= high) return std::string();
        std::ostringstream message;
        message << "must be between " << low << " and " << high << note;
        return message.str();
    };
}

bool wsDeformParamsSet(const Parameters &p) {
    return p.nucleus.setWSDeformParams;
}
bool saveSnapshotsSet(const Parameters &p) { return p.jimwlk.saveSnapshots; }

using P = Parameters;

// Table order is the order parameters are read in (a condition can only
// depend on parameters listed before it) and written to usedParameters*.dat.
const std::vector<ParameterSpec> &parameterTable() {
    static const std::vector<ParameterSpec> table = {
        // general setup
        param("mode", &P::evolution, &EvolutionParameters::mode),
        param("size", &P::lattice, &LatticeParameters::size)
            .check(positive<int>())
            .check(even()),  // the FFTs assume even lattice dimensions
        param("L", &P::lattice, &LatticeParameters::L),
        param("Ny", &P::colorCharge, &ColorChargeParameters::Ny),
        param("roots", &P::collision, &CollisionParameters::roots),
        param("g", &P::coupling, &CouplingParameters::g),
        param("g2mu", &P::collision, &CollisionParameters::g2mu),
        param("maxtime", &P::evolution, &EvolutionParameters::maxtime),
        param(
            "inverseQsForMaxTime", &P::evolution,
            &EvolutionParameters::inverseQsForMaxTime),

        // random seed
        param("seed", &P::random, &RandomParameters::seed),
        param("useSeedList", &P::random, &RandomParameters::useSeedList),
        param("useTimeForSeed", &P::random, &RandomParameters::useTimeForSeed),

        // collision system and geometry
        param("Projectile", &P::collision, &CollisionParameters::Projectile),
        param("Target", &P::collision, &CollisionParameters::Target),
        param("SigmaNN", &P::collision, &CollisionParameters::SigmaNN),
        param("bmin", &P::collision, &CollisionParameters::bmin),
        param("bmax", &P::collision, &CollisionParameters::bmax),
        param(
            "samplebFromLinearDistribution", &P::collision,
            &CollisionParameters::samplebFromLinearDistribution),
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
            "averageOverThisManyNuclei", &P::collision,
            &CollisionParameters::averageOverThisManyNuclei),
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
            "polariztionProjectile", &P::nucleus,
            &NucleusParameters::polariztionProjectile),
        param(
            "polariztionTarget", &P::nucleus,
            &NucleusParameters::polariztionTarget),
        param(
            "polarizationProjectileJz", &P::nucleus,
            &NucleusParameters::polarizationProjectileJz),
        param(
            "polarizationTargetJz", &P::nucleus,
            &NucleusParameters::polarizationTargetJz),

        // Woods-Saxon deformation
        param(
            "setWSDeformParams", &P::nucleus,
            &NucleusParameters::setWSDeformParams),
        param("R_WS", &P::nucleus, &NucleusParameters::R_WS)
            .onlyIf(wsDeformParamsSet),
        param("a_WS", &P::nucleus, &NucleusParameters::a_WS)
            .onlyIf(wsDeformParamsSet),
        param("beta2", &P::nucleus, &NucleusParameters::beta2)
            .onlyIf(wsDeformParamsSet),
        param("beta3", &P::nucleus, &NucleusParameters::beta3)
            .onlyIf(wsDeformParamsSet),
        param("beta4", &P::nucleus, &NucleusParameters::beta4)
            .onlyIf(wsDeformParamsSet),
        param("gamma", &P::nucleus, &NucleusParameters::gamma)
            .onlyIf(wsDeformParamsSet),
        param("dR_np", &P::nucleus, &NucleusParameters::dR_np)
            .onlyIf(wsDeformParamsSet),
        param("da_np", &P::nucleus, &NucleusParameters::da_np)
            .onlyIf(wsDeformParamsSet),
        // Glauber::findNucleusData applies these regardless of
        // setWSDeformParams
        param(
            "force_dmin_flag", &P::nucleus,
            &NucleusParameters::force_dmin_flag),
        param("d_min", &P::nucleus, &NucleusParameters::d_min),

        // nucleon substructure and color charges
        param("m", &P::subnucleon, &SubnucleonParameters::m),
        param("BG", &P::subnucleon, &SubnucleonParameters::BG),
        param("BGq", &P::subnucleon, &SubnucleonParameters::BGq),
        param("BGqVar", &P::subnucleon, &SubnucleonParameters::BGqVar),
        param("dqMin", &P::subnucleon, &SubnucleonParameters::dqMin),
        param("omega", &P::subnucleon, &SubnucleonParameters::omega)
            .check(positive()),
        param(
            "useConstituentQuarkProton", &P::subnucleon,
            &SubnucleonParameters::useConstituentQuarkProton),
        param("NqFluc", &P::subnucleon, &SubnucleonParameters::NqFluc),
        param(
            "shiftConstituentQuarkProtonOrigin", &P::subnucleon,
            &SubnucleonParameters::shiftConstituentQuarkProtonOrigin),
        param(
            "protonAnisotropy", &P::subnucleon,
            &SubnucleonParameters::protonAnisotropy),
        param(
            "SubNucleonParamType", &P::subnucleon,
            &SubnucleonParameters::SubNucleonParamType)
            .check(oneOf({0, 1, 2, 4}, " (0: use the input values)")),
        param(
            "SubNucleonParamSet", &P::subnucleon,
            &SubnucleonParameters::SubNucleonParamSet),
        param("QsmuRatio", &P::colorCharge, &ColorChargeParameters::QsmuRatio),
        param("smearQs", &P::subnucleon, &SubnucleonParameters::smearQs),
        param(
            "smearingWidth", &P::subnucleon,
            &SubnucleonParameters::smearingWidth),
        param("UVdamp", &P::subnucleon, &SubnucleonParameters::UVdamp),
        param(
            "minimumQs2ST", &P::colorCharge,
            &ColorChargeParameters::minimumQs2ST),
        param(
            "NucleusQsTableFileName", &P::colorCharge,
            &ColorChargeParameters::NucleusQsTableFileName),

        // rapidity and x
        param("RapidityA", &P::colorCharge, &ColorChargeParameters::RapidityA),
        param("RapidityB", &P::colorCharge, &ColorChargeParameters::RapidityB),
        param(
            "usePseudoRapidity", &P::colorCharge,
            &ColorChargeParameters::usePseudoRapidity),
        param("Jacobianm", &P::colorCharge, &ColorChargeParameters::Jacobianm),
        param(
            "useFluctuatingx", &P::colorCharge,
            &ColorChargeParameters::useFluctuatingx),
        param(
            "xFromThisFactorTimesQs", &P::colorCharge,
            &ColorChargeParameters::xFromThisFactorTimesQs),

        // running coupling
        param(
            "runningCoupling", &P::coupling,
            &CouplingParameters::runningCoupling),
        param("muZero", &P::coupling, &CouplingParameters::muZero),
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
        param(
            "runWith0Min1Avg2MaxQs", &P::coupling,
            &CouplingParameters::runWith0Min1Avg2MaxQs)
            .check(oneOf({0, 1, 2}, " (0: min, 1: average, 2: max Qs)")),
        param(
            "runWithThisFactorTimesQs", &P::coupling,
            &CouplingParameters::runWithThisFactorTimesQs),
        param(
            "runWithLocalQs", &P::coupling,
            &CouplingParameters::runWithLocalQs),
        // read by the gluon spectrum even without running coupling
        param("runWithkt", &P::coupling, &CouplingParameters::runWithkt)
            .check(oneOf({0, 1})),

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
        param("detaOutput", &P::output, &OutputParameters::detaOutput),

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
        param("useJIMWLK", &P::jimwlk, &JimwlkParameters::useJIMWLK),
        param("mu0_jimwlk", &P::jimwlk, &JimwlkParameters::mu0_jimwlk),
        param(
            "Lambda_QCD_jimwlk", &P::jimwlk,
            &JimwlkParameters::Lambda_QCD_jimwlk)
            .check(positive()),  // in GeV
        param("c_jimwlk", &P::jimwlk, &JimwlkParameters::c_jimwlk)
            .optional("0.2")
            .check(positive()),
        param("m_jimwlk", &P::jimwlk, &JimwlkParameters::m_jimwlk),
        param("alphas_jimwlk", &P::jimwlk, &JimwlkParameters::alphas_jimwlk)
            .check([](const double &v) {
                // 0 selects the running coupling, > 0 a fixed coupling
                return v >= 0 ? "" : "must not be negative";
            }),
        param("Ds_jimwlk", &P::jimwlk, &JimwlkParameters::Ds_jimwlk),
        param("jimwlk_ic_x", &P::jimwlk, &JimwlkParameters::jimwlk_ic_x),
        param(
            "x_projectile_jimwlk", &P::jimwlk,
            &JimwlkParameters::x_projectile_jimwlk),
        param(
            "x_target_jimwlk", &P::jimwlk, &JimwlkParameters::x_target_jimwlk),
        param("saveSnapshots", &P::jimwlk, &JimwlkParameters::saveSnapshots),
        param("xSnapshotList", &P::jimwlk, &JimwlkParameters::xSnapshotList)
            .onlyIf(saveSnapshotsSet),
    };
    return table;
}

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

/// The known parameter name closest to \p key, or "" if none is close.
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

    for (const ParameterSpec &spec : parameterTable()) {
        if (spec.condition && !spec.condition(*this)) continue;
        const InputFile::Entry *entry = input.find(spec.name);
        std::string problem;
        std::string where = source + ": ";
        if (entry) {
            where = source + ":" + std::to_string(entry->line) + ": ";
            problem = spec.apply(*this, entry->value);
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
        std::string message = source + ":" + std::to_string(entry.line)
                              + ": unknown parameter " + key;
        const std::string suggestion = closestName(key);
        if (!suggestion.empty())
            message += " (did you mean " + suggestion + "?)";
        errors.push_back(message);
    }
    if (!errors.empty()) return errors;

    // Derived settings.
    // dtau (in lattice units) is ~0.1, adjusted so maxtime is a whole number
    // of steps. With fewer than one step (e.g. maxtime 0 to only produce
    // Wilson lines) use 0.1 instead of dividing by zero.
    const double latticeSpacing = lattice.L / static_cast<double>(lattice.size);
    const int timeSteps =
        static_cast<int>(10 * evolution.maxtime / latticeSpacing);
    run.dtau = (timeSteps > 0)
                   ? evolution.maxtime / (timeSteps * latticeSpacing)
                   : 0.1;
    // polarized nuclei are only available as configuration files
    if (nucleus.polariztionProjectile != 0 || nucleus.polariztionTarget != 0) {
        nucleus.nucleonPositionsFromFile = 1;
    }
    subnucleon.NqBase = subnucleon.useConstituentQuarkProton;
    if (subnucleon.SubNucleonParamType > 0) {
        loadPosteriorParameterSets(subnucleon.SubNucleonParamType);
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

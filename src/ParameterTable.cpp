// ParameterTable.cpp is part of the IP-Glasma solver.
//
// The single list of every input-file parameter: its key, how it is stored
// in Parameters, its default (if it is optional) and its validity checks.
// Parameters::readInput() and Parameters::writeInputParameters() are both
// driven by this table, so adding a parameter only takes one entry here
// (plus its member/getter/setter in Parameters.h).

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
/// `param("size", &Parameters::setSize, &Parameters::getSize).check(even())`.
template <typename T>
class Param {
  public:
    template <typename Arg, typename Ret>
    Param(
        const char *name, void (Parameters::*set)(Arg),
        Ret (Parameters::*get)() const)
        : name_(name),
          set_([set](Parameters &p, const T &v) { (p.*set)(v); }),
          write_([get](const Parameters &p, std::ostream &out) {
              writeValue(out, (p.*get)());
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

template <typename Arg, typename Ret>
Param<std::decay_t<Arg>> param(
    const char *name, void (Parameters::*set)(Arg),
    Ret (Parameters::*get)() const) {
    return Param<std::decay_t<Arg>>(name, set, get);
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

bool wsDeformParamsSet(const Parameters &p) { return p.getSetWSDeformParams(); }
bool saveSnapshotsSet(const Parameters &p) { return p.getSaveSnapshots(); }

using P = Parameters;

// Table order is the order parameters are read in (a condition can only
// depend on parameters listed before it) and written to usedParameters*.dat.
const std::vector<ParameterSpec> &parameterTable() {
    static const std::vector<ParameterSpec> table = {
        // general setup
        param("mode", &P::setMode, &P::getMode),
        param("size", &P::setSize, &P::getSize)
            .check(positive<int>())
            .check(even()),  // the FFTs assume even lattice dimensions
        param("L", &P::setL, &P::getL),
        param("Ny", &P::setNy, &P::getNy),
        param("roots", &P::setRoots, &P::getRoots),
        param("g", &P::setg, &P::getg),
        param("g2mu", &P::setg2mu, &P::getg2mu),
        param("maxtime", &P::setMaxtime, &P::getMaxtime),
        param(
            "inverseQsForMaxTime", &P::setInverseQsForMaxTime,
            &P::getInverseQsForMaxTime),

        // random seed
        param("seed", &P::setSeed, &P::getSeed),
        param("useSeedList", &P::setUseSeedList, &P::getUseSeedList),
        param("useTimeForSeed", &P::setUseTimeForSeed, &P::getUseTimeForSeed),

        // collision system and geometry
        param("Projectile", &P::setProjectile, &P::getProjectile),
        param("Target", &P::setTarget, &P::getTarget),
        param("SigmaNN", &P::setSigmaNN, &P::getSigmaNN),
        param("bmin", &P::setbmin, &P::getbmin),
        param("bmax", &P::setbmax, &P::getbmax),
        param("samplebFromLinearDistribution", &P::setLinearb, &P::getLinearb),
        param(
            "rotateReactionPlane", &P::setRotateReactionPlane,
            &P::getRotateReactionPlane),
        param("useNucleus", &P::setUseNucleus, &P::getUseNucleus),
        param("useGaussian", &P::setUseGaussian, &P::getUseGaussian),
        param(
            "useSmoothNucleus", &P::setUseSmoothNucleus,
            &P::getUseSmoothNucleus),
        param("useFixedNpart", &P::setUseFixedNpart, &P::getUseFixedNpart),
        param(
            "averageOverThisManyNuclei", &P::setAverageOverNuclei,
            &P::getAverageOverNuclei),
        param(
            "gaussianWounding", &P::setGaussianWounding,
            &P::getGaussianWounding),

        // nucleon positions
        param(
            "nucleonPositionsFromFile", &P::setNucleonPositionsFromFile,
            &P::getNucleonPositionsFromFile),
        param(
            "nuclearConfigurationsPath", &P::setNuclearConfigurationsPath,
            &P::getNuclearConfigurationsPath)
            .optional("./nucleusConfigurations"),
        param(
            "lightNucleusOption", &P::setlightNucleusOption,
            &P::getlightNucleusOption),
        param(
            "polariztionProjectile", &P::setPolarizationProjectile,
            &P::getPolarizationProjectile),
        param(
            "polariztionTarget", &P::setPolarizationTarget,
            &P::getPolarizationTarget),
        param(
            "polarizationProjectileJz", &P::setPolarizationProjectileJz,
            &P::getPolarizationProjectileJz),
        param(
            "polarizationTargetJz", &P::setPolarizationTargetJz,
            &P::getPolarizationTargetJz),

        // Woods-Saxon deformation
        param(
            "setWSDeformParams", &P::setSetWSDeformParams,
            &P::getSetWSDeformParams),
        param("R_WS", &P::setR_WS, &P::getR_WS).onlyIf(wsDeformParamsSet),
        param("a_WS", &P::setA_WS, &P::getA_WS).onlyIf(wsDeformParamsSet),
        param("beta2", &P::setBeta2, &P::getBeta2).onlyIf(wsDeformParamsSet),
        param("beta3", &P::setBeta3, &P::getBeta3).onlyIf(wsDeformParamsSet),
        param("beta4", &P::setBeta4, &P::getBeta4).onlyIf(wsDeformParamsSet),
        param("gamma", &P::setGamma, &P::getGamma).onlyIf(wsDeformParamsSet),
        param("dR_np", &P::setWSdR_np, &P::getWSdR_np)
            .onlyIf(wsDeformParamsSet),
        param("da_np", &P::setWSda_np, &P::getWSda_np)
            .onlyIf(wsDeformParamsSet),
        // Glauber::findNucleusData applies these regardless of
        // setWSDeformParams
        param("force_dmin_flag", &P::setForceDmin, &P::getForceDmin),
        param("d_min", &P::setDmin, &P::getDmin),

        // nucleon substructure and color charges
        param("m", &P::setm, &P::getm),
        param("BG", &P::setBG, &P::getBG),
        param("BGq", &P::setBGq, &P::getBGq),
        param("BGqVar", &P::setBGqVar, &P::getBGqVar),
        param("dqMin", &P::setDqmin, &P::getDqmin),
        param("omega", &P::setOmega, &P::getOmega).check(positive()),
        param(
            "useConstituentQuarkProton", &P::setUseConstituentQuarkProton,
            &P::getUseConstituentQuarkProton),
        param("NqFluc", &P::setNqFluc, &P::getNqFluc),
        param(
            "shiftConstituentQuarkProtonOrigin",
            &P::setShiftConstituentQuarkProtonOrigin,
            &P::getShiftConstituentQuarkProtonOrigin),
        param(
            "protonAnisotropy", &P::setProtonAnisotropy,
            &P::getProtonAnisotropy),
        param(
            "SubNucleonParamType", &P::setSubNucleonParamType,
            &P::getSubNucleonParamType)
            .check(oneOf({0, 1, 2, 4}, " (0: use the input values)")),
        param(
            "SubNucleonParamSet", &P::setSubNucleonParamSet,
            &P::getSubNucleonParamSet),
        param("QsmuRatio", &P::setQsmuRatio, &P::getQsmuRatio),
        param("smearQs", &P::setSmearQs, &P::getSmearQs),
        param("smearingWidth", &P::setSmearingWidth, &P::getSmearingWidth),
        param("UVdamp", &P::setUVdamp, &P::getUVdamp),
        param("minimumQs2ST", &P::setMinimumQs2ST, &P::getMinimumQs2ST),
        param(
            "NucleusQsTableFileName", &P::setNucleusQsTableFileName,
            &P::getNucleusQsTableFileName),

        // rapidity and x
        param("RapidityA", &P::setRapidityA, &P::getRapidityA),
        param("RapidityB", &P::setRapidityB, &P::getRapidityB),
        param(
            "usePseudoRapidity", &P::setUsePseudoRapidity,
            &P::getUsePseudoRapidity),
        param("Jacobianm", &P::setJacobianm, &P::getJacobianm),
        param(
            "useFluctuatingx", &P::setUseFluctuatingx, &P::getUseFluctuatingx),
        param(
            "xFromThisFactorTimesQs", &P::setxFromThisFactorTimesQs,
            &P::getxFromThisFactorTimesQs),

        // running coupling
        param(
            "runningCoupling", &P::setRunningCoupling, &P::getRunningCoupling),
        param("muZero", &P::setMuZero, &P::getMuZero),
        param("c", &P::setc, &P::getc).check(positive()),
        param("nFlavors", &P::setNFlavors, &P::getNFlavors)
            .optional("3")
            .check(inRange(
                0, 16,
                " (the one-loop beta-function coefficient 11*Nc - "
                "2*nFlavors must be positive)")),
        param("LambdaQCD", &P::setLambdaQCD, &P::getLambdaQCD)
            .optional("0.2")
            .check(positive()),
        param("runWith0Min1Avg2MaxQs", &P::setRunWithQs, &P::getRunWithQs)
            .check(oneOf({0, 1, 2}, " (0: min, 1: average, 2: max Qs)")),
        param(
            "runWithThisFactorTimesQs", &P::setRunWithThisFactorTimesQs,
            &P::getRunWithThisFactorTimesQs),
        param("runWithLocalQs", &P::setRunWithLocalQs, &P::getRunWithLocalQs),
        // read by the gluon spectrum even without running coupling
        param("runWithkt", &P::setRunWithkt, &P::getRunWithkt)
            .check(oneOf({0, 1})),

        // observables
        param(
            "computeGluonMultiplicity", &P::setComputeGluonMultiplicity,
            &P::getComputeGluonMultiplicity),
        param(
            "readMultFromFile", &P::setReadMultFromFile,
            &P::getReadMultFromFile),

        // output
        param("writeOutputs", &P::setWriteOutputs, &P::getWriteOutputs),
        param(
            "writeEpsilonUHydro", &P::setWriteEpsilonUHydro,
            &P::getWriteEpsilonUHydro)
            .optional("1"),
        param(
            "writeTmunuBinary", &P::setWriteTmunuBinary,
            &P::getWriteTmunuBinary)
            .optional("1"),
        param(
            "writeOutputsToHDF5", &P::setWriteOutputsToHDF5,
            &P::getWriteOutputsToHDF5),
        param("LOutput", &P::setLOutput, &P::getLOutput),
        param("sizeOutput", &P::setSizeOutput, &P::getSizeOutput),
        param("etaSizeOutput", &P::setEtaSizeOutput, &P::getEtaSizeOutput),
        param("detaOutput", &P::setDetaOutput, &P::getDetaOutput),

        // Wilson lines
        param(
            "writeWilsonLines", &P::setWriteWilsonLines,
            &P::getWriteWilsonLines)
            .check(oneOf({0, 1, 2}, " (0: none, 1: text, 2: binary)")),
        param("wilsonLinePath", &P::setWilsonLinePath, &P::getWilsonLinePath)
            .optional("./"),
        param(
            "readInitialWilsonLines", &P::setReadInitialWilsonLines,
            &P::getReadInitialWilsonLines)
            .check(oneOf({0, 1, 2}, " (0: sample, 1: text, 2: binary)")),

        // JIMWLK
        param("useJIMWLK", &P::setUseJIMWLK, &P::getUseJIMWLK),
        param("mu0_jimwlk", &P::setMu0_jimwlk, &P::getMu0_jimwlk),
        param(
            "Lambda_QCD_jimwlk", &P::setLambdaQCD_jimwlk,
            &P::getLambdaQCD_jimwlk)
            .check(positive()),  // in GeV
        param("c_jimwlk", &P::setc_jimwlk, &P::getc_jimwlk)
            .optional("0.2")
            .check(positive()),
        param("m_jimwlk", &P::setm_jimwlk, &P::getm_jimwlk),
        param("alphas_jimwlk", &P::setJimwlk_alphas, &P::getJimwlk_alphas)
            .check([](const double &v) {
                // 0 selects the running coupling, > 0 a fixed coupling
                return v >= 0 ? "" : "must not be negative";
            }),
        param("Ds_jimwlk", &P::setDs_jimwlk, &P::getDs_jimwlk),
        param("jimwlk_ic_x", &P::setJimwlk_x0, &P::getJimwlk_x0),
        param(
            "x_projectile_jimwlk", &P::setJimwlk_x_projectile,
            &P::getJimwlk_x_projectile),
        param(
            "x_target_jimwlk", &P::setJimwlk_x_target, &P::getJimwlk_x_target),
        param("saveSnapshots", &P::setSaveSnapshots, &P::getSaveSnapshots),
        param("xSnapshotList", &P::setxSnapshotList, &P::getxSnapshotList)
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
    const double latticeSpacing = getL() / static_cast<double>(getSize());
    const int timeSteps = static_cast<int>(10 * getMaxtime() / latticeSpacing);
    setdtau(
        (timeSteps > 0) ? getMaxtime() / (timeSteps * latticeSpacing) : 0.1);
    // polarized nuclei are only available as configuration files
    if (getPolarizationProjectile() != 0 || getPolarizationTarget() != 0) {
        setNucleonPositionsFromFile(1);
    }
    setNqBase(getUseConstituentQuarkProton());
    if (getSubNucleonParamType() > 0) {
        loadPosteriorParameterSets(getSubNucleonParamType());
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

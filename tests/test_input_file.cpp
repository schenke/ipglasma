#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include "InputFile.h"
#include "Parameters.h"
#include "doctest.h"

namespace {
InputFile inputFromText(const std::string &text) {
    std::istringstream in(text);
    return InputFile(in, "test");
}

std::string readSourceFile(const std::string &relativePath) {
    std::ifstream in(std::string(IPGLASMA_SOURCE_DIR) + "/" + relativePath);
    REQUIRE(in.good());
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

// `text` with `key`'s line replaced by `key value` (or removed if value
// is empty).
std::string inputWith(
    const std::string &text, const std::string &key, const std::string &value) {
    std::istringstream in(text);
    std::string line, result;
    bool found = false;
    while (std::getline(in, line)) {
        std::istringstream tokens(line);
        std::string first;
        tokens >> first;
        if (first == key) {
            found = true;
            if (!value.empty()) result += key + " " + value + "\n";
        } else {
            result += line + "\n";
        }
    }
    REQUIRE(found);
    return result;
}

// The shipped example input, with `key`'s line replaced by `key value`
// (or removed if value is empty).
std::string exampleInputWith(const std::string &key, const std::string &value) {
    return inputWith(readSourceFile("input"), key, value);
}

// Inserts extra lines before the EndOfFile line of an input text.
std::string insertBeforeEndOfFile(
    const std::string &text, const std::string &extra) {
    const std::size_t end = text.find("\nEndOfFile");
    REQUIRE(end != std::string::npos);
    return text.substr(0, end + 1) + extra + text.substr(end + 1);
}

std::vector<std::string> readErrors(const std::string &text) {
    Parameters param;
    return param.readInput(inputFromText(text));
}

bool anyContains(
    const std::vector<std::string> &errors, const std::string &part) {
    for (const std::string &error : errors) {
        if (error.find(part) != std::string::npos) return true;
    }
    return false;
}
}  // namespace

// ---- InputFile ----

TEST_CASE("InputFile: reads key/value pairs, comments and EndOfFile") {
    const InputFile input = inputFromText(
        "# full-line comment\n"
        "size 256\n"
        "\n"
        "L 30.5   # trailing comment\n"
        "EndOfFile\n"
        "ignored after EndOfFile\n");
    CHECK(input.errors().empty());
    CHECK(input.entries().size() == 2);
    REQUIRE(input.find("size") != nullptr);
    CHECK(input.find("size")->value == "256");
    CHECK(input.find("size")->line == 2);
    CHECK(input.find("L")->value == "30.5");
    CHECK(input.find("ignored") == nullptr);
    CHECK(input.find("missing") == nullptr);
}

TEST_CASE("InputFile: EndOfFile must be alone on its line") {
    const InputFile input = inputFromText("size 256\nEndOfFile typo\nL 30\n");
    REQUIRE(input.errors().size() == 1);
    CHECK(
        input.errors()[0]
        == "test:2: unexpected text after EndOfFile ('typo' ...)");
    CHECK(input.find("L") == nullptr);  // still the end of the input
}

TEST_CASE("InputFile: EndOfFile is optional") {
    const InputFile input = inputFromText("size 256\nL 30");
    CHECK(input.errors().empty());
    CHECK(input.entries().size() == 2);
}

TEST_CASE(
    "InputFile: reports missing values, extra values and duplicate keys "
    "with their line numbers") {
    const InputFile input = inputFromText(
        "size\n"
        "L 30 40\n"
        "m 0.4\n"
        "m 0.5\n");
    REQUIRE(input.errors().size() == 3);
    CHECK(input.errors()[0] == "test:1: no value given for size");
    CHECK(
        input.errors()[1]
        == "test:2: more than one value given for L ('40' ...)");
    CHECK(input.errors()[2] == "test:4: m is already set on line 3");
    // the first occurrence is kept
    CHECK(input.find("m")->value == "0.4");
}

TEST_CASE("InputFile: reports a missing file") {
    const InputFile input("this_input_file_does_not_exist");
    REQUIRE(input.errors().size() == 1);
    CHECK(anyContains(input.errors(), "cannot open"));
}

// ---- parseValue ----

TEST_CASE("parseValue: integers must be written as integers") {
    int i = -7;
    CHECK(parseValue("42", i));
    CHECK(i == 42);
    CHECK(parseValue("-3", i));
    CHECK(i == -3);
    CHECK(parseValue("+5", i));
    CHECK(i == 5);
    for (const char *bad :
         {"3.0", "3.", "1e3", "abc", "", "4x", "99999999999", "+-5", "+"}) {
        CAPTURE(bad);
        CHECK_FALSE(parseValue(bad, i));
    }
    CHECK(i == 5);  // unchanged by the failed parses

    long long l = 0;
    CHECK(parseValue("9223372036854775807", l));
    CHECK(l == 9223372036854775807LL);
    CHECK(parseValue("-1", l));
    CHECK(l == -1);
    CHECK_FALSE(parseValue("9223372036854775808", l));
}

TEST_CASE("parseValue: flags must be 0 or 1") {
    bool flag = false;
    CHECK(parseValue("1", flag));
    CHECK(flag == true);
    CHECK(parseValue("0", flag));
    CHECK(flag == false);
    for (const char *bad : {"2", "-1", "true", "1.0", ""}) {
        CAPTURE(bad);
        CHECK_FALSE(parseValue(bad, flag));
    }
}

TEST_CASE("parseValue: doubles must be finite numbers") {
    double d = 0.;
    CHECK(parseValue("0.5", d));
    CHECK(d == 0.5);
    CHECK(parseValue("1.", d));
    CHECK(d == 1.);
    CHECK(parseValue("-2e-3", d));
    CHECK(d == -2e-3);
    for (const char *bad : {"nan", "inf", "1x", "", "0.5.1"}) {
        CAPTURE(bad);
        CHECK_FALSE(parseValue(bad, d));
    }
}

TEST_CASE("parseValue: lists are comma-separated doubles") {
    std::vector<double> list;
    CHECK(parseValue("5e-3,2e-3,0.0001", list));
    CHECK(list == std::vector<double> {5e-3, 2e-3, 0.0001});
    CHECK(parseValue("7", list));
    CHECK(list == std::vector<double> {7.});
    CHECK(parseValue("none", list));
    CHECK(list.empty());
    for (const char *bad : {"1,,2", "1,", ",1", "1;2", "", "none,1"}) {
        CAPTURE(bad);
        CHECK_FALSE(parseValue(bad, list));
    }
}

// ---- Parameters::readInput (the parameter table) ----

TEST_CASE("Parameters::readInput: the shipped input files are valid") {
    for (const char *file :
         {"input", "validations/input_eccentricity",
          "validations/input_vm_proton"}) {
        CAPTURE(file);
        Parameters param;
        const std::vector<std::string> errors =
            param.readInput(inputFromText(readSourceFile(file)));
        for (const std::string &error : errors) CAPTURE(error);
        CHECK(errors.empty());
        CHECK(param.validationErrors().empty());
    }
}

TEST_CASE("Parameters::readInput: reads values and sets derived ones") {
    Parameters param;
    REQUIRE(param.readInput(inputFromText(readSourceFile("input"))).empty());
    CHECK(param.lattice.size == 256);
    CHECK(param.lattice.L == 30.);
    CHECK(param.collision.projectile == "Pb");
    CHECK(param.jimwlk.xSnapshotList.size() == 5);
    CHECK(param.jimwlk.enabled == true);
    // derived: NqBase from Nq, dtau ~0.1 with
    // maxTime a whole number of steps
    CHECK(param.subnucleon.NqBase == param.subnucleon.Nq);
    const double a = param.lattice.L / param.lattice.size;
    const double steps = param.evolution.maxTime / (a * param.run.dtau);
    CHECK(steps == doctest::Approx(static_cast<int>(steps + 0.5)));
    CHECK(param.run.dtau == doctest::Approx(0.1).epsilon(0.05));
}

TEST_CASE("Parameters::readInput: maxTime 0 gives dtau 0.1, not NaN") {
    Parameters param;
    REQUIRE(param.readInput(inputFromText(exampleInputWith("maxTime", "0")))
                .empty());
    CHECK(param.run.dtau == 0.1);
}

TEST_CASE("Parameters::readInput: unknown keys are errors, with a suggestion") {
    const std::vector<std::string> errors = readErrors(insertBeforeEndOfFile(
        readSourceFile("input"), "writeTmunuBinry 0\nfoo 1\n"));
    REQUIRE(errors.size() == 2);
    CHECK(anyContains(
        errors,
        "unknown parameter writeTmunuBinry (did you mean writeTmunuBinary?)"));
    CHECK(anyContains(errors, "unknown parameter foo"));
    CHECK_FALSE(anyContains(errors, "foo (did you mean"));
}

TEST_CASE("Parameters::readInput: required keys must be given") {
    const std::vector<std::string> errors =
        readErrors(exampleInputWith("sigmaNN", ""));
    REQUIRE(errors.size() == 1);
    CHECK(errors[0] == "test: sigmaNN is required but not given");
}

TEST_CASE(
    "Parameters::readInput: a conditional key is only required, and only "
    "read, when its setting is on") {
    // the example input has runningCoupling 0
    CHECK(readErrors(exampleInputWith("mu0", "")).empty());
    const std::vector<std::string> errors = readErrors(
        inputWith(exampleInputWith("mu0", ""), "runningCoupling", "1"));
    REQUIRE(errors.size() == 1);
    CHECK(errors[0] == "test: mu0 is required but not given");

    // given while unused, it is accepted, ignored and not written out
    Parameters param;
    REQUIRE(
        param.readInput(inputFromText(exampleInputWith("mu0", "0.7"))).empty());
    CHECK(param.coupling.mu0 == CouplingParameters {}.mu0);
    std::ostringstream used;
    param.writeInputParameters(used);
    CHECK(used.str().find("\nmu0 ") == std::string::npos);

    // the JIMWLK parameters are not needed without JIMWLK
    std::string text = inputWith(
        exampleInputWith("useJIMWLK", "0"), "jimwlkSaveSnapshots", "");
    for (const char *key :
         {"jimwlkMu0", "jimwlkLambdaQCD", "jimwlkMass", "jimwlkAlphaS",
          "jimwlkDs", "jimwlkInitialX", "jimwlkXSnapshotList"}) {
        text = inputWith(text, key, "");
    }
    CHECK(readErrors(text).empty());
}

TEST_CASE("Parameters::readInput: optional keys fall back to their default") {
    Parameters param;
    REQUIRE(param.readInput(inputFromText(exampleInputWith("nFlavors", "")))
                .empty());
    CHECK(param.coupling.nFlavors == 3);
    Parameters param2;
    REQUIRE(
        param2.readInput(inputFromText(exampleInputWith("wilsonLinePath", "")))
            .empty());
    CHECK(param2.wilsonLines.wilsonLinePath == "./");
    Parameters param3;
    REQUIRE(param3
                .readInput(inputFromText(
                    exampleInputWith("writeWilsonLineGeometry", "")))
                .empty());
    CHECK(param3.wilsonLines.writeGeometry);
}

TEST_CASE("Parameters::readInput: values must have the parameter's type") {
    CHECK(anyContains(
        readErrors(exampleInputWith("useJIMWLK", "2")),
        "useJIMWLK '2' is not 0 or 1"));
    CHECK(anyContains(
        readErrors(inputWith(
            exampleInputWith("runWithKt", "2"), "runningCoupling", "1")),
        "runWithKt '2' is not 0 or 1"));
    CHECK(anyContains(
        readErrors(exampleInputWith("size", "256.0")),
        "size '256.0' is not an integer"));
    CHECK(anyContains(
        readErrors(exampleInputWith("m", "0.4GeV")),
        "m '0.4GeV' is not a number"));
    CHECK(anyContains(
        readErrors(exampleInputWith("jimwlkXSnapshotList", "1e-3;1e-4")),
        "jimwlkXSnapshotList '1e-3;1e-4' is not a comma-separated list"));
}

TEST_CASE("Parameters::readInput: per-value checks") {
    struct Case {
        const char *key;
        const char *value;
        const char *message;
    };
    for (const Case &c : std::vector<Case> {
             {"size", "0", "must be positive"},
             {"jimwlkAlphaS", "-0.3", "must not be negative"},
             {"jimwlkDs", "0", "must be positive"},
             {"nFlavors", "-1", "must be between 0 and 16"},
             {"L", "0", "must be positive"},
             {"maxTime", "-1", "must not be negative"},
             {"nucleiToAverage", "0", "must be positive"},
             {"polarizationTarget", "3", "must be one of 0, 1, 2"},
             {"Nq", "0", "must be at least 1"},
             {"Nq", "0.5", "must be at least 1"},
             {"nucleonModel", "stringy", "must be one of gaussian, hotspots"},
             {"size", "255", "must be even"},
             {"omega", "0", "must be positive"},
             {"subNucleonParamType", "3", "must be one of 0, 1, 2, 4"},
             {"runWithQs", "3", "must be one of 0, 1, 2"},
             {"nFlavors", "17", "must be between 0 and 16"},
             {"LambdaQCD", "0", "must be positive"},
             {"c", "-0.2", "must be positive"},
             {"jimwlkC", "0", "must be positive"},
             {"jimwlkLambdaQCD", "0", "must be positive"},
             {"writeWilsonLines", "3", "must be one of 0, 1, 2"},
             {"readInitialWilsonLines", "3", "must be one of 0, 1, 2"},
             {"sqrtS", "0", "must be positive"},
             {"projectileX", "0", "must be positive"},
             {"targetX", "-1e-3", "must be positive"},
             {"jimwlkInitialX", "0", "must be positive"},
             {"seed", "-2", "must be at least -1"},
         }) {
        CAPTURE(c.key);
        CAPTURE(c.value);
        // with running coupling, so its parameters are read
        const std::vector<std::string> errors = readErrors(inputWith(
            exampleInputWith(c.key, c.value), "runningCoupling", "1"));
        REQUIRE(errors.size() == 1);
        CHECK(anyContains(errors, std::string(c.key) + " " + c.value));
        CHECK(anyContains(errors, c.message));
    }
}

TEST_CASE("Parameters::readInput: reports every problem at once") {
    std::string text = exampleInputWith("size", "255");
    text = insertBeforeEndOfFile(text, "unknownKey 1\n");
    std::istringstream in(text);
    std::string withoutM;
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("m ", 0) != 0) withoutM += line + "\n";
    }
    const std::vector<std::string> errors = readErrors(withoutM);
    CHECK(errors.size() == 3);
    CHECK(anyContains(errors, "must be even"));
    CHECK(anyContains(errors, "m is required"));
    CHECK(anyContains(errors, "unknown parameter unknownKey"));
}

TEST_CASE(
    "Parameters::readInput: Woods-Saxon deformation parameters are only "
    "read with useInputWSParams 1") {
    // without them, useInputWSParams 0 is fine and 1 is not
    std::string text;
    {
        std::istringstream in(exampleInputWith("useInputWSParams", "0"));
        std::string line;
        while (std::getline(in, line)) {
            if (line.rfind("radiusWS ", 0) != 0) text += line + "\n";
        }
    }
    CHECK(readErrors(text).empty());

    std::string text1;
    {
        std::istringstream in(text);
        std::string line;
        while (std::getline(in, line)) {
            text1 += (line.rfind("useInputWSParams ", 0) == 0)
                         ? "useInputWSParams 1\n"
                         : line + "\n";
        }
    }
    const std::vector<std::string> errors = readErrors(text1);
    REQUIRE(errors.size() == 1);
    CHECK(errors[0] == "test: radiusWS is required but not given");
}

TEST_CASE(
    "Parameters::writeInputParameters output reads back to the same "
    "parameters") {
    Parameters param;
    REQUIRE(param.readInput(inputFromText(readSourceFile("input"))).empty());
    std::ostringstream written;
    param.writeInputParameters(written);

    Parameters reread;
    const std::vector<std::string> errors =
        reread.readInput(inputFromText(written.str()));
    for (const std::string &error : errors) CAPTURE(error);
    REQUIRE(errors.empty());
    std::ostringstream rewritten;
    reread.writeInputParameters(rewritten);
    CHECK(rewritten.str() == written.str());
    // exact doubles, not rounded ones
    CHECK(reread.colorCharge.projectileX == param.colorCharge.projectileX);
    CHECK(reread.run.dtau == param.run.dtau);
}

TEST_CASE("InputFile: skips a UTF-8 byte order mark") {
    const InputFile input = inputFromText("\xEF\xBB\xBFsize 256\n");
    CHECK(input.errors().empty());
    CHECK(input.find("size") != nullptr);
}

TEST_CASE("Parameters::readInput: a missing file gives a single error") {
    Parameters param;
    const std::vector<std::string> errors =
        param.readInput(InputFile("this_input_file_does_not_exist"));
    REQUIRE(errors.size() == 1);
    CHECK(anyContains(errors, "cannot open"));
}

TEST_CASE(
    "Parameters::readInput: no arbitrary suggestions for short unknown keys") {
    for (const char *key : {"xx", "R", "ab"}) {
        CAPTURE(key);
        const std::vector<std::string> errors =
            readErrors(insertBeforeEndOfFile(
                readSourceFile("input"), std::string(key) + " 1\n"));
        REQUIRE(errors.size() == 1);
        CHECK_FALSE(anyContains(errors, "did you mean"));
    }
}

TEST_CASE(
    "Parameters::readInput: outputTimes is only read with a field output, "
    "defaults to none and must be positive") {
    // the example input writes no field output and has outputTimes none
    auto withLine = [](const std::string &text, const std::string &key,
                       const std::string &value) {
        std::istringstream in(text);
        std::string line, result;
        while (std::getline(in, line)) {
            result += (line.rfind(key + " ", 0) == 0) ? key + " " + value + "\n"
                                                      : line + "\n";
        }
        return result;
    };
    const std::string withTmunu = exampleInputWith("writeTmunu", "1");
    {
        Parameters param;
        REQUIRE(param
                    .readInput(inputFromText(
                        withLine(withTmunu, "outputTimes", "0.1,0.3")))
                    .empty());
        CHECK(param.output.outputTimes == std::vector<double> {0.1, 0.3});
        std::ostringstream written;
        param.writeInputParameters(written);
        CHECK(
            written.str().find("\noutputTimes 0.1,0.3\n") != std::string::npos);
    }
    {
        // not read without a field output
        Parameters param;
        REQUIRE(param
                    .readInput(
                        inputFromText(exampleInputWith("outputTimes", "0.1")))
                    .empty());
        CHECK(param.output.outputTimes.empty());
    }
    {
        // optional, and written back as none
        Parameters param;
        REQUIRE(
            param
                .readInput(inputFromText(withLine(
                    exampleInputWith("outputTimes", ""), "writeTmunu", "1")))
                .empty());
        CHECK(param.output.outputTimes.empty());
        std::ostringstream written;
        param.writeInputParameters(written);
        CHECK(written.str().find("\noutputTimes none\n") != std::string::npos);
    }
    const std::vector<std::string> errors =
        readErrors(withLine(withTmunu, "outputTimes", "0.1,0"));
    REQUIRE(errors.size() == 1);
    CHECK(anyContains(errors, "every time must be positive"));
}

TEST_CASE(
    "Parameters::readInput: projectileX and targetX are only read with "
    "useFluctuatingX 0") {
    // fluctuating x (which excludes JIMWLK): x comes from the local Q_s
    std::string text;
    {
        std::istringstream in(
            exampleInputWith("projectileX", ""));  // drop projectileX
        std::string line;
        while (std::getline(in, line)) {
            if (line.rfind("targetX ", 0) == 0) continue;
            if (line.rfind("useFluctuatingX ", 0) == 0) {
                line = "useFluctuatingX 1";
            } else if (line.rfind("useJIMWLK ", 0) == 0) {
                line = "useJIMWLK 0";
            }
            text += line + "\n";
        }
    }
    Parameters param;
    CHECK(param.readInput(inputFromText(text)).empty());
    CHECK(param.validationErrors().empty());

    // a fixed x (the shipped input): both are required
    const std::vector<std::string> errors =
        readErrors(exampleInputWith("projectileX", ""));
    REQUIRE(errors.size() == 1);
    CHECK(errors[0] == "test: projectileX is required but not given");
}

TEST_CASE(
    "Parameters::readInput: xQsFactor is only read, and must be positive, "
    "with useFluctuatingX 1") {
    // a fixed x (the shipped input): not needed
    CHECK(readErrors(exampleInputWith("xQsFactor", "")).empty());
    CHECK(readErrors(exampleInputWith("xQsFactor", "0")).empty());

    // a fluctuating x (without JIMWLK)
    auto fluctuatingInputWith = [](const std::string &value) {
        std::istringstream in(exampleInputWith("xQsFactor", value));
        std::string line, text;
        while (std::getline(in, line)) {
            if (line.rfind("useFluctuatingX ", 0) == 0) {
                line = "useFluctuatingX 1";
            } else if (line.rfind("useJIMWLK ", 0) == 0) {
                line = "useJIMWLK 0";
            }
            text += line + "\n";
        }
        return text;
    };
    CHECK(readErrors(fluctuatingInputWith("1")).empty());
    std::vector<std::string> errors = readErrors(fluctuatingInputWith(""));
    REQUIRE(errors.size() == 1);
    CHECK(errors[0] == "test: xQsFactor is required but not given");
    for (const std::string value : {"0", "-1"}) {
        CAPTURE(value);
        errors = readErrors(fluctuatingInputWith(value));
        REQUIRE(errors.size() == 1);
        CHECK(anyContains(errors, "xQsFactor " + value));
        CHECK(anyContains(errors, "must be positive"));
    }
}

TEST_CASE(
    "Parameters::readInput: jimwlkXSnapshotList is only read with "
    "jimwlkSaveSnapshots 1") {
    std::string text;
    {
        std::istringstream in(exampleInputWith("jimwlkXSnapshotList", ""));
        std::string line;
        while (std::getline(in, line)) {
            text += (line.rfind("jimwlkSaveSnapshots ", 0) == 0)
                        ? "jimwlkSaveSnapshots 0\n"
                        : line + "\n";
        }
    }
    CHECK(readErrors(text).empty());

    const std::vector<std::string> errors = readErrors(
        exampleInputWith("jimwlkXSnapshotList", ""));  // jimwlkSaveSnapshots 1
    REQUIRE(errors.size() == 1);
    CHECK(errors[0] == "test: jimwlkXSnapshotList is required but not given");
}

TEST_CASE(
    "Parameters::readInput: a malformed condition key gives one error, not "
    "an indeterminate set of follow-up errors") {
    const std::vector<std::string> errors =
        readErrors(exampleInputWith("useInputWSParams", "1.0"));
    REQUIRE(errors.size() == 1);
    CHECK(anyContains(errors, "useInputWSParams '1.0' is not 0 or 1"));
}

TEST_CASE(
    "Parameters::readInput: keys renamed in IP-Glasma 2.0 name their new key") {
    const std::vector<std::string> errors = readErrors(
        insertBeforeEndOfFile(exampleInputWith("dMin", ""), "d_min 0.9\n"));
    REQUIRE(errors.size() == 1);
    CHECK(anyContains(errors, "d_min was renamed to dMin"));

    // given together with the new key, the old one is just unknown
    const std::vector<std::string> both = readErrors(
        insertBeforeEndOfFile(readSourceFile("input"), "d_min 0.9\n"));
    REQUIRE(both.size() == 1);
    CHECK(anyContains(both, "unknown parameter d_min (renamed to dMin)"));
}

TEST_CASE(
    "Parameters::readInput: a pre-2.0 input file reports each rename once") {
    std::istringstream in(readSourceFile("input"));
    std::string line, oldInput;
    const std::vector<std::pair<std::string, std::string>> renames = {
        {"maxTime", "maxtime"},
        {"dMin", "d_min"},
        {"writeWilsonLines", "writeInitialWilsonLines"}};
    while (std::getline(in, line)) {
        for (const auto &[newKey, oldKey] : renames) {
            if (line.rfind(newKey + " ", 0) == 0) {
                line = oldKey + line.substr(newKey.size());
            }
        }
        oldInput += line + "\n";
    }
    const std::vector<std::string> errors = readErrors(oldInput);
    CHECK(errors.size() == renames.size());
    for (const auto &[newKey, oldKey] : renames) {
        CHECK(anyContains(errors, oldKey + " was renamed to " + newKey));
    }
}

TEST_CASE(
    "Parameters::readInput: a removed or replaced pre-2.0 key gets a hint") {
    const std::vector<std::string> errors = readErrors(insertBeforeEndOfFile(
        readSourceFile("input"),
        "\nNc 3\nwriteOutputs 2\nuseGaussian 0\ng2mu 0.1\nmode "
        "1\nuseRandomSeed 0\n"));
    for (const std::string &error : errors) CAPTURE(error);
    CHECK(errors.size() == 6);
    CHECK(anyContains(errors, "unknown parameter Nc (removed: "));
    CHECK(anyContains(
        errors, "unknown parameter writeOutputs (replaced by writeHydro"));
    CHECK(anyContains(
        errors,
        "unknown parameter useGaussian (removed: useNucleus 0 always uses a "
        "constant g2muGeV)"));
    CHECK(anyContains(
        errors, "unknown parameter g2mu (replaced by g2muGeV, in GeV"));
    CHECK(anyContains(
        errors,
        "unknown parameter mode (replaced by runEvolution: 1 for mode "
        "1, 0 for any other mode)"));
    CHECK(anyContains(
        errors,
        "unknown parameter useRandomSeed (replaced by seed -1, which draws a "
        "random seed)"));
}

TEST_CASE(
    "Parameters::readInput: the removed rapidityA/rapidityB point to their "
    "replacements") {
    const std::vector<std::string> errors = readErrors(insertBeforeEndOfFile(
        readSourceFile("input"), "rapidityA 0\nrapidityB 0\n"));
    for (const std::string &error : errors) CAPTURE(error);
    CHECK(errors.size() == 2);
    CHECK(anyContains(
        errors, "unknown parameter rapidityA (replaced by projectileX"));
    CHECK(anyContains(
        errors, "unknown parameter rapidityB (replaced by targetX"));
}

TEST_CASE("Parameters::readInput: minimumQs2ST is a non-negative real number") {
    Parameters param;
    REQUIRE(
        param.readInput(inputFromText(exampleInputWith("minimumQs2ST", "12.5")))
            .empty());
    CHECK(param.colorCharge.minimumQs2ST == 12.5);
    const std::vector<std::string> errors =
        readErrors(exampleInputWith("minimumQs2ST", "-1"));
    REQUIRE(errors.size() == 1);
    CHECK(anyContains(errors, "minimumQs2ST -1: must not be negative"));
}

TEST_CASE("Parameters::readInput: a fractional Nq sets a fractional NqBase") {
    Parameters param;
    REQUIRE(
        param.readInput(inputFromText(exampleInputWith("Nq", "2.5"))).empty());
    CHECK(param.subnucleon.Nq == 2.5);
    CHECK(param.subnucleon.NqBase == 2.5);
}

TEST_CASE("Parameters::readInput: rejects a maxTime needing too many steps") {
    const std::vector<std::string> errors =
        readErrors(exampleInputWith("maxTime", "1e12"));
    REQUIRE(errors.size() == 1);
    CHECK(anyContains(errors, "more than 1e8 time steps"));
}

TEST_CASE(
    "Parameters::readInput: each nucleon model only reads its own "
    "parameters") {
    // hotspots (the shipped input): the hot-spot keys are required
    CHECK(anyContains(
        readErrors(exampleInputWith("BGq", "")), "BGq is required"));

    // gaussian: the hot-spot keys are accepted and ignored...
    {
        Parameters param;
        const std::vector<std::string> errors = param.readInput(
            inputFromText(exampleInputWith("nucleonModel", "gaussian")));
        for (const std::string &error : errors) CAPTURE(error);
        CHECK(errors.empty());
    }
    // ...and not needed
    std::string gaussian;
    {
        std::istringstream in(exampleInputWith("nucleonModel", "gaussian"));
        std::string line;
        while (std::getline(in, line)) {
            bool hotSpotKey = false;
            for (const char *key :
                 {"BGq ", "BGqVar ", "dqMin ", "omega ", "Nq ", "NqFluc ",
                  "shiftConstituentQuarkProtonOrigin "}) {
                if (line.rfind(key, 0) == 0) hotSpotKey = true;
            }
            if (!hotSpotKey) gaussian += line + "\n";
        }
    }
    Parameters param;
    const std::vector<std::string> errors =
        param.readInput(inputFromText(gaussian));
    for (const std::string &error : errors) CAPTURE(error);
    CHECK(errors.empty());
    CHECK(param.subnucleon.nucleonModel == "gaussian");
}

TEST_CASE(
    "Parameters::readInput: nucleonModel strings reads the hot-spot "
    "parameters except the number of hot spots, and allows no posterior "
    "set") {
    std::string strings;
    {
        std::istringstream in(exampleInputWith("nucleonModel", "strings"));
        std::string line;
        while (std::getline(in, line)) {
            if (line.rfind("Nq ", 0) == 0 || line.rfind("NqFluc ", 0) == 0) {
                continue;  // always three hot spots
            }
            strings += line + "\n";
        }
    }
    Parameters param;
    const std::vector<std::string> errors =
        param.readInput(inputFromText(strings));
    for (const std::string &error : errors) CAPTURE(error);
    CHECK(errors.empty());
    CHECK(param.subnucleon.nucleonModel == "strings");
    CHECK(param.validationErrors().empty());

    // the hot-spot parameters are still required
    std::string withoutBGq;
    {
        std::istringstream in(strings);
        std::string line;
        while (std::getline(in, line)) {
            if (line.rfind("BGq ", 0) != 0) withoutBGq += line + "\n";
        }
    }
    CHECK(anyContains(readErrors(withoutBGq), "BGq is required"));

    // the posterior sets are fits of hot-spot nucleons
    std::string posterior;
    {
        std::istringstream in(strings);
        std::string line;
        while (std::getline(in, line)) {
            if (line.rfind("subNucleonParamType ", 0) == 0) {
                line = "subNucleonParamType 2";
            }
            posterior += line + "\n";
        }
    }
    Parameters withPosterior;
    REQUIRE(withPosterior.readInput(inputFromText(posterior)).empty());
    CHECK(anyContains(
        withPosterior.validationErrors(),
        "requires nucleonModel hotspots, not strings"));
}

TEST_CASE(
    "Parameters::readInput: a posterior parameter set replaces m, BG, BGq, "
    "smearingWidth, QsMuRatio, dqMin and Nq, which are then not read") {
    // readInput() loads tables/posterior_Nq3.csv for type 2; provide a
    // one-row table if none exists in the working directory
    const std::string table = "tables/posterior_Nq3.csv";
    const bool ownTable = !std::filesystem::exists(table);
    const bool ownDirectory = !std::filesystem::exists("tables");
    if (ownTable) {
        std::filesystem::create_directories("tables");
        std::ofstream out(table);
        out << "m,BG,BGq,smearingWidth,QsmuRatio,dqmin\n"
            << "0.3,4.0,0.3,0.5,0.6,0.2\n";
    }

    std::string text;
    {
        std::istringstream in(exampleInputWith("subNucleonParamType", "2"));
        std::string line;
        while (std::getline(in, line)) {
            bool replaced = false;
            for (const char *key :
                 {"m ", "BG ", "BGq ", "smearingWidth ", "QsMuRatio ", "dqMin ",
                  "Nq "}) {
                if (line.rfind(key, 0) == 0) replaced = true;
            }
            if (!replaced) text += line + "\n";
        }
    }
    Parameters param;
    const std::vector<std::string> errors =
        param.readInput(inputFromText(text));
    for (const std::string &error : errors) CAPTURE(error);
    CHECK(errors.empty());
    CHECK(param.validationErrors().empty());

    // the replaced keys are accepted and ignored when present
    Parameters withReplacedKeys;
    const std::vector<std::string> replacedErrors = withReplacedKeys.readInput(
        inputFromText(exampleInputWith("subNucleonParamType", "2")));
    for (const std::string &error : replacedErrors) CAPTURE(error);
    CHECK(replacedErrors.empty());

    if (ownTable) std::remove(table.c_str());
    if (ownDirectory) {
        std::error_code ignored;  // keep a directory that is not empty
        std::filesystem::remove("tables", ignored);
    }
}

TEST_CASE(
    "Parameters::readInput: subNucleonParamSet must be -1 (random) or a set "
    "index") {
    std::string text;
    {
        std::istringstream in(exampleInputWith("subNucleonParamType", "2"));
        std::string line;
        while (std::getline(in, line)) {
            if (line.rfind("subNucleonParamSet ", 0) == 0) {
                line = "subNucleonParamSet -2";
            }
            text += line + "\n";
        }
    }
    CHECK(anyContains(
        readErrors(text),
        "subNucleonParamSet -2: must be -1 (random) or a set index >= 0"));
}

TEST_CASE(
    "Parameters::readInput: checks the ranges of the nucleon parameters") {
    struct Case {
        const char *key;
        const char *value;
        const char *error;
    };
    for (const Case &c :
         {Case {"BG", "0", "BG 0: must be positive"},
          Case {"BGq", "0.09", "BGq 0.09: must be larger than 0.09"},
          Case {"BGq", "0.05", "BGq 0.05: must be larger than 0.09"},
          Case {"BGqVar", "-0.1", "BGqVar -0.1: must not be negative"},
          Case {"dqMin", "-0.2", "dqMin -0.2: must not be negative"},
          Case {
              "smearingWidth", "-0.5",
              "smearingWidth -0.5: must not be negative"}}) {
        CAPTURE(c.key);
        CAPTURE(c.value);
        CHECK(
            anyContains(readErrors(exampleInputWith(c.key, c.value)), c.error));
    }
    CHECK(readErrors(exampleInputWith("BGq", "0.091")).empty());
}

TEST_CASE("Parameters::readInput: protonAnisotropy must be larger than -1") {
    const std::string gaussian = exampleInputWith("nucleonModel", "gaussian");
    auto withAnisotropy = [&gaussian](const std::string &value) {
        std::istringstream in(gaussian);
        std::string line, text;
        while (std::getline(in, line)) {
            if (line.rfind("protonAnisotropy ", 0) == 0) {
                line = "protonAnisotropy " + value;
            }
            text += line + "\n";
        }
        return text;
    };
    CHECK(anyContains(
        readErrors(withAnisotropy("-1")),
        "protonAnisotropy -1: must be larger than -1"));
    CHECK(readErrors(withAnisotropy("-0.5")).empty());
}

TEST_CASE(
    "Parameters::readInput: useConstituentQuarkProton points to "
    "nucleonModel") {
    std::istringstream in(readSourceFile("input"));
    std::string line, oldInput;
    while (std::getline(in, line)) {
        if (line.rfind("nucleonModel ", 0) == 0) continue;
        if (line.rfind("Nq ", 0) == 0) line = "useConstituentQuarkProton 3";
        oldInput += line + "\n";
    }
    const std::vector<std::string> errors = readErrors(oldInput);
    for (const std::string &error : errors) CAPTURE(error);
    CHECK(anyContains(
        errors,
        "unknown parameter useConstituentQuarkProton (replaced by "
        "nucleonModel"));
    CHECK(anyContains(errors, "nucleonModel is required"));
    CHECK_FALSE(anyContains(errors, "renamed to Nq"));
}

TEST_CASE(
    "Parameters::readInput: a posterior set with gaussian nucleons is "
    "reported as such, without needing the posterior table") {
    Parameters param;
    std::string text;
    {
        std::istringstream in(exampleInputWith("nucleonModel", "gaussian"));
        std::string line;
        while (std::getline(in, line)) {
            if (line.rfind("subNucleonParamType ", 0) == 0) {
                line = "subNucleonParamType 4";
            }
            text += line + "\n";
        }
    }
    const std::vector<std::string> errors =
        param.readInput(inputFromText(text));
    for (const std::string &error : errors) CAPTURE(error);
    CHECK(errors.empty());
    CHECK(anyContains(
        param.validationErrors(),
        "subNucleonParamType = 4 (a posterior parameter set) requires "
        "nucleonModel hotspots"));
}

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

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

// A file in the working directory, removed again at the end of the test.
class TemporaryFile {
  public:
    TemporaryFile(const std::string &name, const std::string &contents)
        : name_(name) {
        std::ofstream out(name_);
        out << contents;
    }
    ~TemporaryFile() { std::remove(name_.c_str()); }
    const std::string &name() const { return name_; }

  private:
    std::string name_;
};

// The shipped example input, with the first parameter set of the posterior
// file `file`.
std::string exampleInputWithPosterior(const std::string &file) {
    return inputWith(
        exampleInputWith("posteriorFile", file), "posteriorSet", "0");
}

// A posterior file shipped in tables/.
std::string shippedPosterior(const std::string &name) {
    return std::string(IPGLASMA_SOURCE_DIR) + "/tables/" + name;
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
         {"input", "validations/eccentricity_regression/input_eccentricity",
          "validations/vector_meson_production/input_vm_proton"}) {
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
    // with the JIMWLK snapshots, which the example input does not write
    Parameters param;
    REQUIRE(param
                .readInput(
                    inputFromText(exampleInputWith("jimwlkSaveSnapshots", "1")))
                .empty());
    CHECK(param.lattice.size == 800);
    CHECK(param.lattice.L == 30.);
    CHECK(param.collision.projectile == "Pb");
    CHECK(param.jimwlk.xSnapshotList.size() == 4);
    CHECK(param.jimwlk.enabled == true);
    // derived: dtau ~0.1 with maxTime a whole number of steps
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

TEST_CASE(
    "Parameters::readInput: timeStepsPerLatticeSpacing sets the time step, "
    "10 by default") {
    // the time step before the parameter existed (issue #109)
    auto dtauWith = [](const std::string &text) {
        Parameters param;
        REQUIRE(param.readInput(inputFromText(text)).empty());
        return param.run.dtau;
    };
    const std::string input = readSourceFile("input");
    Parameters example;
    REQUIRE(example.readInput(inputFromText(input)).empty());
    const double a = example.lattice.L / example.lattice.size;
    const double maxTime = example.evolution.maxTime;
    const double oldDtau = maxTime / (static_cast<int>(10 * maxTime / a) * a);
    CHECK(example.evolution.timeStepsPerLatticeSpacing == 10);
    CHECK(example.run.dtau == oldDtau);
    // optional: without the key it is 10 as well
    CHECK(
        dtauWith(inputWith(input, "timeStepsPerLatticeSpacing", ""))
        == oldDtau);
    for (const int n : {1, 20, 40, 80}) {
        CAPTURE(n);
        const double dtau = dtauWith(
            inputWith(input, "timeStepsPerLatticeSpacing", std::to_string(n)));
        // maxTime is a whole number of steps, the largest number up to
        // n maxTime / a
        const double steps = maxTime / (a * dtau);
        CHECK(steps == doctest::Approx(static_cast<int>(n * maxTime / a)));
        CHECK(steps == doctest::Approx(static_cast<int>(steps + 0.5)));
        // fewer than one step: 1/n
        CHECK(
            dtauWith(inputWith(
                exampleInputWith("maxTime", "0"), "timeStepsPerLatticeSpacing",
                std::to_string(n)))
            == 1. / n);
    }
    // the step limit counts the steps with the chosen time step: 1e8 a is
    // 3.75e6 fm/c on the lattice of the example input
    const std::vector<std::string> errors = readErrors(inputWith(
        exampleInputWith("maxTime", "3e6"), "timeStepsPerLatticeSpacing", "1"));
    CHECK(errors.empty());
    CHECK(anyContains(
        readErrors(inputWith(
            exampleInputWith("maxTime", "3e6"), "timeStepsPerLatticeSpacing",
            "1000")),
        "more than 1e8 time steps"));
    // an integer
    CHECK(anyContains(
        readErrors(exampleInputWith("timeStepsPerLatticeSpacing", "12.5")),
        "timeStepsPerLatticeSpacing"));
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

    // sqrtS is only needed with fluctuating x or pseudorapidity (both off
    // in the example input), and smearingWidth only with smearQs 1
    CHECK(readErrors(exampleInputWith("sqrtS", "")).empty());
    CHECK(readErrors(
              inputWith(exampleInputWith("smearingWidth", ""), "smearQs", "0"))
              .empty());
    CHECK(anyContains(
        readErrors(exampleInputWith("smearingWidth", "")),
        "smearingWidth is required but not given"));

    // the configuration files are only read when the nucleons are sampled,
    // not when the Wilson lines are read from file
    const std::string fromFiles = inputWith(
        inputWith(
            exampleInputWith("nucleonPositionsFromFile", "1"),
            "lightNucleusOption", ""),
        "nuclearConfigurationsPath", "");
    CHECK(anyContains(
        readErrors(fromFiles), "lightNucleusOption is required but not given"));
    Parameters reading;
    REQUIRE(reading
                .readInput(inputFromText(
                    inputWith(fromFiles, "readInitialWilsonLines", "2")))
                .empty());
    std::ostringstream readingUsed;
    reading.writeInputParameters(readingUsed);
    CHECK(readingUsed.str().find("lightNucleusOption") == std::string::npos);
}

TEST_CASE("Parameters::readInput: optional keys fall back to their default") {
    Parameters param;
    REQUIRE(param.readInput(inputFromText(exampleInputWith("nFlavors", "")))
                .empty());
    CHECK(param.coupling.nFlavors == 3);
    // wilsonLinePath is only read when Wilson lines are written or read,
    // and writeWilsonLineGeometry when they are written
    Parameters param2;
    REQUIRE(param2
                .readInput(inputFromText(inputWith(
                    exampleInputWith("wilsonLinePath", ""), "writeWilsonLines",
                    "2")))
                .empty());
    CHECK(param2.wilsonLines.wilsonLinePath == "./");
    Parameters param3;
    REQUIRE(param3
                .readInput(inputFromText(inputWith(
                    exampleInputWith("writeWilsonLineGeometry", ""),
                    "writeWilsonLines", "2")))
                .empty());
    CHECK(param3.wilsonLines.writeGeometry);
    // read with its default, not left at the initial value
    std::ostringstream used;
    param3.writeInputParameters(used);
    CHECK(
        used.str().find("\nwriteWilsonLineGeometry 1\n") != std::string::npos);
    // the example input has the JIMWLK running coupling (jimwlkAlphaS 0)
    Parameters param4;
    REQUIRE(param4
                .readInput(inputFromText(
                    exampleInputWith("jimwlkRunningCoupling", "")))
                .empty());
    CHECK(param4.jimwlk.runningCoupling == 1);
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
        readErrors(inputWith(
            exampleInputWith("jimwlkXSnapshotList", "1e-3;1e-4"),
            "jimwlkSaveSnapshots", "1")),
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
             {"timeStepsPerLatticeSpacing", "0", "must be positive"},
             {"timeStepsPerLatticeSpacing", "-10", "must be positive"},
             {"nucleiToAverage", "0", "must be positive"},
             {"polarizationTarget", "3", "must be one of 0, 1, 2"},
             {"Nq", "0", "must be at least 1"},
             {"Nq", "0.5", "must be at least 1"},
             {"nucleonModel", "stringy", "must be one of gaussian, hotspots"},
             {"size", "255", "must be even"},
             {"omega", "0", "must be positive"},
             {"runWithQs", "3", "must be one of 0, 1, 2"},
             {"nFlavors", "17", "must be between 0 and 16"},
             {"LambdaQCD", "0", "must be positive"},
             {"c", "-0.2", "must be positive"},
             {"jimwlkC", "0", "must be positive"},
             {"jimwlkLambdaQCD", "0", "must be positive"},
             {"jimwlkRunningCoupling", "3", "must be one of 1, 2"},
             {"writeWilsonLines", "3", "must be one of 0, 1, 2"},
             {"readInitialWilsonLines", "3", "must be one of 0, 1, 2"},
             {"sqrtS", "0", "must be positive"},
             {"projectileX", "0", "must be positive"},
             {"targetX", "-1e-3", "must be positive"},
             {"jimwlkInitialX", "0", "must be positive"},
             {"seed", "-2", "must be at least -1"},
             {"colorChargeSeed", "-2", "must be at least -1"},
         }) {
        CAPTURE(c.key);
        CAPTURE(c.value);
        // with running coupling and pseudorapidity, so their parameters
        // (and sqrtS) are read
        const std::vector<std::string> errors = readErrors(inputWith(
            inputWith(exampleInputWith(c.key, c.value), "runningCoupling", "1"),
            "usePseudoRapidity", "1"));
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
    // the example input writes the hydro fields (writeHydro 1) and has
    // outputTimes none
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
        REQUIRE(
            param
                .readInput(inputFromText(withLine(
                    exampleInputWith("outputTimes", "0.1"), "writeHydro", "0")))
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
    "Parameters::validationErrors: jimwlkLambdaQCD must be smaller than "
    "jimwlkMu0 with the running coupling, with either prescription") {
    for (const std::string prescription : {"1", "2"}) {
        CAPTURE(prescription);
        const std::string text =
            exampleInputWith("jimwlkRunningCoupling", prescription);
        for (const std::string mu0 : {"0.040", "0.03"}) {  // = and < Lambda
            CAPTURE(mu0);
            Parameters param;
            REQUIRE(
                param
                    .readInput(inputFromText(inputWith(text, "jimwlkMu0", mu0)))
                    .empty());
            CHECK(anyContains(
                param.validationErrors(), "must be smaller than jimwlkMu0"));
        }
        // with a fixed coupling neither is used
        Parameters fixed;
        REQUIRE(fixed
                    .readInput(inputFromText(inputWith(
                        inputWith(text, "jimwlkMu0", "0.03"), "jimwlkAlphaS",
                        "0.2")))
                    .empty());
        CHECK_FALSE(anyContains(
            fixed.validationErrors(), "must be smaller than jimwlkMu0"));
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

    const std::vector<std::string> errors = readErrors(inputWith(
        exampleInputWith("jimwlkXSnapshotList", ""), "jimwlkSaveSnapshots",
        "1"));
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

TEST_CASE("Parameters::readInput: the removed eta grid gets a hint") {
    const std::vector<std::string> errors = readErrors(insertBeforeEndOfFile(
        readSourceFile("input"), "etaSizeOutput 1\ndEtaOutput 0\n"));
    REQUIRE(errors.size() == 2);
    CHECK(anyContains(
        errors,
        "unknown parameter etaSizeOutput (removed: the fields are boost "
        "invariant"));
    CHECK(anyContains(errors, "unknown parameter dEtaOutput (removed:"));
}

TEST_CASE("Parameters::readInput: Nq can be fractional") {
    Parameters param;
    REQUIRE(
        param.readInput(inputFromText(exampleInputWith("Nq", "2.5"))).empty());
    CHECK(param.subnucleon.Nq == 2.5);
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
    "parameters except the number of hot spots") {
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

    // so a posterior file cannot set the number of hot spots
    CHECK(anyContains(
        readErrors(inputWith(
            inputWith(
                strings, "posteriorFile", shippedPosterior("posterior_Nq3.csv")),
            "posteriorSet", "0")),
        "column Nq is only used with nucleonModel hotspots, which this input "
        "does not set"));
}

TEST_CASE(
    "Parameters::readInput: a posterior file sets the parameters of its "
    "columns, which the input need not give") {
    TemporaryFile file(
        "ipglasma_test_posterior.csv",
        "# m [GeV], BG [GeV^-2], QsMuRatio\r\n"
        "0.3, 4.0, 0.6\r\n"
        "\n"
        "# a comment\n"
        "0.25,4.5,0.7\n");
    // m only in the posterior file; BG and QsMuRatio also in the input
    Parameters param;
    const std::vector<std::string> errors = param.readInput(inputFromText(
        inputWith(exampleInputWithPosterior(file.name()), "m", "")));
    for (const std::string &error : errors) CAPTURE(error);
    REQUIRE(errors.empty());
    CHECK(param.validationErrors().empty());
    // the first set, exactly as written (not rounded to float)
    CHECK(param.subnucleon.m == 0.3);
    CHECK(param.subnucleon.BG == 4.0);
    CHECK(param.colorCharge.QsMuRatio == 0.6);
    // the values of the input are replaced, with a warning
    const std::vector<std::string> &warnings = param.inputWarnings();
    REQUIRE(warnings.size() == 2);
    CHECK(anyContains(
        warnings,
        "BG 3.3 is replaced by the posterior file ipglasma_test_posterior.csv"));
    CHECK(anyContains(warnings, "QsMuRatio 0.643 is replaced"));

    param.applyPosteriorSet(1);
    CHECK(param.event.posteriorSet == 1);
    CHECK(param.subnucleon.m == 0.25);
    CHECK(param.subnucleon.BG == 4.5);
    CHECK(param.colorCharge.QsMuRatio == 0.7);
    // modulo the number of sets
    param.applyPosteriorSet(4);
    CHECK(param.event.posteriorSet == 0);
    CHECK(param.subnucleon.m == 0.3);

    // the parameters written are those of the set
    param.applyPosteriorSet(1);
    std::ostringstream written;
    param.writeInputParameters(written);
    CHECK(written.str().find("\nm 0.25\n") != std::string::npos);
    CHECK(written.str().find("\nBG 4.5\n") != std::string::npos);
    CHECK(written.str().find("\nposteriorSet 1\n") != std::string::npos);

    // without a posterior file, no warnings
    Parameters withoutPosterior;
    REQUIRE(withoutPosterior.readInput(inputFromText(readSourceFile("input")))
                .empty());
    CHECK(withoutPosterior.inputWarnings().empty());
}

TEST_CASE(
    "Parameters::writeInputParameters: with a posterior file, the output "
    "reads back to the parameters of the event, also with posteriorSet -1") {
    TemporaryFile file(
        "ipglasma_test_posterior.csv",
        "m, BG\n0.3, 4.0\n0.25, 4.5\n0.35, 3.9\n");
    Parameters param;
    REQUIRE(param
                .readInput(inputFromText(inputWith(
                    exampleInputWithPosterior(file.name()), "posteriorSet",
                    "-1")))
                .empty());
    // as main() does with the set it drew
    param.applyPosteriorSet(2);
    std::ostringstream written;
    param.writeInputParameters(written);

    Parameters reread;
    const std::vector<std::string> errors =
        reread.readInput(inputFromText(written.str()));
    for (const std::string &error : errors) CAPTURE(error);
    REQUIRE(errors.empty());
    CHECK(reread.posterior.set == 2);
    reread.applyPosteriorSet(reread.posterior.set);
    std::ostringstream rewritten;
    reread.writeInputParameters(rewritten);
    CHECK(rewritten.str() == written.str());
    CHECK(reread.subnucleon.m == 0.35);
    CHECK(reread.subnucleon.BG == 3.9);
}

TEST_CASE(
    "Parameters::readInput: the shipped posterior files set every parameter "
    "of their sets") {
    for (const char *name :
         {"posterior.csv", "posterior_Nq3.csv", "posterior5020_Nq3.csv"}) {
        CAPTURE(name);
        Parameters param;
        const std::vector<std::string> errors = param.readInput(
            inputFromText(exampleInputWithPosterior(shippedPosterior(name))));
        for (const std::string &error : errors) CAPTURE(error);
        REQUIRE(errors.empty());
        CHECK(param.validationErrors().empty());
        // m, BG, BGq, smearingWidth, Nq, QsMuRatio, dqMin, BGqVar, NqFluc
        CHECK(param.inputWarnings().size() == 9);
        CHECK(param.subnucleon.BGqVar == 0.);
        CHECK(param.subnucleon.NqFluc == 0.);
        if (std::string(name) == "posterior.csv") {
            CHECK(param.subnucleon.m == 0.506);
            CHECK(param.subnucleon.Nq == 6.79);
            CHECK(param.colorCharge.QsMuRatio == 0.598);
        } else {
            CHECK(param.subnucleon.Nq == 3.);
        }
    }
}

TEST_CASE(
    "Parameters::readInput: the columns of a posterior file are physical "
    "input keys, each once") {
    struct Case {
        const char *contents;
        const char *error;
    };
    for (const Case &c :
         {Case {
              "Qsg2muRatio\n0.6\n",
              ":1: unknown parameter Qsg2muRatio (did you mean QsMuRatio?)"},
          Case {
              "size\n128\n",
              ":1: size cannot be set by a posterior file, only the physical "
              "parameters can"},
          Case {"jimwlkDs\n0.003\n", ":1: jimwlkDs cannot be set"},
          Case {"g\n1\n", ":1: g cannot be set"},
          Case {"jacobianMass\n0.35\n", ":1: jacobianMass cannot be set"},
          Case {"posteriorSet\n1\n", ":1: posteriorSet cannot be set"},
          Case {"BG, BG\n4, 4\n", ":1: column BG appears twice"},
          Case {"m [GeV], \n0.3, 4\n", ":1: a column has no name"},
          Case {
              "m, BG\n0.3\n",
              ":2: 1 values, but line 1 names 2 columns"},
          Case {"m\n", " contains no parameter sets"},
          Case {"\n", " has no line naming the columns"},
          Case {"0.3, 4\n", " has no line naming the columns"},
          // a comment before the names is not the names
          Case {
              "# from arXiv:2502.05138\n0.3\n",
              ":1: unknown parameter from arXiv:2502.05138"}}) {
        CAPTURE(c.contents);
        TemporaryFile file("ipglasma_test_posterior.csv", c.contents);
        const std::vector<std::string> errors =
            readErrors(exampleInputWithPosterior(file.name()));
        for (const std::string &error : errors) CAPTURE(error);
        CHECK(anyContains(errors, c.error));
    }
    CHECK(anyContains(
        readErrors(exampleInputWithPosterior("no_such_posterior.csv")),
        "cannot open the posterior file no_such_posterior.csv"));
}

TEST_CASE(
    "Parameters::readInput: comment lines and a byte order mark before the "
    "names of the columns of a posterior file") {
    for (const char *contents :
         {"# sets of a posterior\n#\n# m [GeV], BG [GeV^-2]\n0.3, 4\n",
          "# sets of a posterior\n\nm, BG\n# the first set\n0.3, 4\n",
          "\xEF\xBB\xBF# m [GeV], BG [GeV^-2]\n0.3, 4\n",
          "\xEF\xBB\xBFm, BG\n0.3, 4\n"}) {
        CAPTURE(contents);
        TemporaryFile file("ipglasma_test_posterior.csv", contents);
        Parameters param;
        const std::vector<std::string> errors = param.readInput(
            inputFromText(exampleInputWithPosterior(file.name())));
        for (const std::string &error : errors) CAPTURE(error);
        REQUIRE(errors.empty());
        CHECK(param.subnucleon.m == 0.3);
        CHECK(param.subnucleon.BG == 4.);
    }
}

TEST_CASE(
    "Parameters::readInput: an old name of a parameter of a posterior column "
    "is only reported as renamed") {
    TemporaryFile file("ipglasma_test_posterior.csv", "UVDamp\n0.1\n");
    std::string text = exampleInputWithPosterior(file.name());
    text.replace(text.find("\nUVDamp "), 8, "\nUVdamp ");
    const std::vector<std::string> errors = readErrors(text);
    for (const std::string &error : errors) CAPTURE(error);
    CHECK(anyContains(errors, "UVdamp was renamed to UVDamp"));
    CHECK_FALSE(anyContains(errors, "column UVDamp"));
}

TEST_CASE(
    "Parameters::readInput: every set of a posterior file passes the checks "
    "of its values") {
    TemporaryFile file(
        "ipglasma_test_posterior.csv", "m, BG\n0.3, 4\n0.3, 0\n0.3, x\n");
    const std::vector<std::string> errors =
        readErrors(exampleInputWithPosterior(file.name()));
    for (const std::string &error : errors) CAPTURE(error);
    CHECK(errors.size() == 2);
    CHECK(anyContains(
        errors, "ipglasma_test_posterior.csv:3: BG 0: must be positive"));
    CHECK(anyContains(
        errors, "ipglasma_test_posterior.csv:4: BG 'x' is not a number"));

    // also the first set, which readInput() applies
    TemporaryFile first("ipglasma_test_posterior_first.csv", "BG\n0\n4\n");
    CHECK(anyContains(
        readErrors(exampleInputWithPosterior(first.name())),
        "ipglasma_test_posterior_first.csv:2: BG 0: must be positive"));
}

TEST_CASE(
    "Parameters::readInput: a posterior file sets only parameters the input "
    "uses, and its sets do not switch between models") {
    // gaussian nucleons have no hot spots
    TemporaryFile hotSpots("ipglasma_test_posterior.csv", "BGq\n0.3\n");
    CHECK(anyContains(
        readErrors(inputWith(
            exampleInputWithPosterior(hotSpots.name()), "nucleonModel",
            "gaussian")),
        "ipglasma_test_posterior.csv: column BGq is only used with nucleonModel "
        "hotspots or strings, which this input does not set"));

    // the Woods-Saxon parameters need useInputWSParams 1 (the example input
    // has 0)
    for (const char *name :
         {"radiusWS", "diffusenessWS", "beta2", "beta3", "beta4", "gamma",
          "deltaRnp", "deltaAnp"}) {
        CAPTURE(name);
        TemporaryFile woodsSaxon(
            "ipglasma_test_posterior_ws.csv",
            std::string(name) + "\n0.5\n");
        const std::string text = exampleInputWithPosterior(woodsSaxon.name());
        CHECK(anyContains(
            readErrors(text),
            "column " + std::string(name)
                + " is only used with useInputWSParams 1, which this input "
                  "does not set"));
        const std::vector<std::string> errors =
            readErrors(inputWith(text, "useInputWSParams", "1"));
        for (const std::string &error : errors) CAPTURE(error);
        CHECK(errors.empty());
    }

    // the shipped tables set smearingWidth, which needs smearQs 1
    CHECK(anyContains(
        readErrors(inputWith(
            exampleInputWithPosterior(shippedPosterior("posterior_Nq3.csv")),
            "smearQs", "0")),
        "column smearingWidth is only used with smearQs 1, which this input "
        "does not set"));

    // jimwlkAlphaS 0 is the running coupling of JIMWLK, > 0 a fixed one
    TemporaryFile coupling(
        "ipglasma_test_posterior_coupling.csv", "jimwlkAlphaS\n0\n0.2\n");
    CHECK(anyContains(
        readErrors(exampleInputWithPosterior(coupling.name())),
        "ipglasma_test_posterior_coupling.csv:3: the set changes whether "
        "jimwlkRunningCoupling is used"));
}

TEST_CASE(
    "Parameters::validationErrors: the checks combining parameters run with "
    "every set of a posterior file") {
    TemporaryFile file(
        "ipglasma_test_posterior.csv", "LambdaQCD\n0.2\n0.25\n0.35\n");
    Parameters param;
    const std::vector<std::string> errors =
        param.readInput(inputFromText(inputWith(
            exampleInputWithPosterior(file.name()), "runningCoupling", "1")));
    for (const std::string &error : errors) CAPTURE(error);
    REQUIRE(errors.empty());
    const std::vector<std::string> validation = param.validationErrors();
    REQUIRE(validation.size() == 1);
    CHECK(
        validation[0].rfind(
            "ipglasma_test_posterior.csv:4: with this set, LambdaQCD (0.35) "
            "must be smaller than mu0 (0.3)",
            0)
        == 0);
}

TEST_CASE(
    "Parameters::readInput: subNucleonParamType and subNucleonParamSet point "
    "to posteriorFile and posteriorSet") {
    CHECK(anyContains(
        readErrors(insertBeforeEndOfFile(
            readSourceFile("input"), "subNucleonParamType 2\n")),
        "unknown parameter subNucleonParamType (replaced by posteriorFile: "
        "none for 0, tables/posterior.csv for 1, tables/posterior_Nq3.csv "
        "for 2"));
    CHECK(anyContains(
        readErrors(insertBeforeEndOfFile(
            inputWith(
                exampleInputWith(
                    "posteriorFile", shippedPosterior("posterior_Nq3.csv")),
                "posteriorSet", ""),
            "subNucleonParamSet -1\n")),
        "subNucleonParamSet was renamed to posteriorSet"));
}

TEST_CASE(
    "Parameters::readInput: posteriorSet must be -1 (random) or a set "
    "index") {
    CHECK(anyContains(
        readErrors(inputWith(
            exampleInputWithPosterior(shippedPosterior("posterior_Nq3.csv")),
            "posteriorSet", "-2")),
        "posteriorSet -2: must be -1 (random) or a set index >= 0"));
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
        // with smearQs 1, so smearingWidth is read
        CHECK(anyContains(
            readErrors(
                inputWith(exampleInputWith(c.key, c.value), "smearQs", "1")),
            c.error));
    }
    CHECK(readErrors(exampleInputWith("BGq", "0.091")).empty());
}

TEST_CASE(
    "Parameters::readInput: rejects values that would hang the run or "
    "divide by zero") {
    struct Case {
        const char *key;
        const char *value;
        const char *error;
    };
    for (const Case &c :
         {Case {"sigmaNN", "0", "sigmaNN 0: must be positive"},
          Case {
              "useFixedNpart", "-1", "useFixedNpart -1: must not be negative"},
          Case {"Ny", "0", "Ny 0: must be positive"},
          Case {"g", "0", "g 0: must be positive"},
          Case {"QsMuRatio", "0", "QsMuRatio 0: must be positive"},
          Case {"UVDamp", "-0.1", "UVDamp -0.1: must not be negative"},
          Case {"bMin", "-1", "bMin -1: must not be negative"},
          Case {"radiusWS", "0", "radiusWS 0: must be positive"},
          Case {"diffusenessWS", "0", "diffusenessWS 0: must be positive"},
          Case {"LOutput", "0", "LOutput 0: must be positive"},
          Case {"sizeOutput", "0", "sizeOutput 0: must be positive"}}) {
        CAPTURE(c.key);
        CAPTURE(c.value);
        // with input Woods-Saxon parameters and a field output, so those
        // parameters are read
        const std::vector<std::string> errors = readErrors(inputWith(
            inputWith(
                exampleInputWith(c.key, c.value), "useInputWSParams", "1"),
            "writeHydro", "1"));
        REQUIRE(errors.size() == 1);
        CHECK(anyContains(errors, c.error));
    }
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

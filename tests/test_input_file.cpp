#include <fstream>
#include <sstream>
#include <string>
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

// The shipped example input, with `key`'s line replaced by `key value`
// (or removed if value is empty).
std::string exampleInputWith(const std::string &key, const std::string &value) {
    std::istringstream in(readSourceFile("input"));
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

    unsigned long long u = 0;
    CHECK(parseValue("18446744073709551615", u));
    CHECK(u == 18446744073709551615ULL);
    CHECK_FALSE(parseValue("-1", u));
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
    for (const char *bad : {"1,,2", "1,", ",1", "1;2", ""}) {
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
    CHECK(param.collision.Projectile == "Pb");
    CHECK(param.jimwlk.xSnapshotList.size() == 5);
    CHECK(param.jimwlk.useJIMWLK == true);
    // derived: NqBase from useConstituentQuarkProton, dtau ~0.1 with
    // maxtime a whole number of steps
    CHECK(
        param.subnucleon.NqBase == param.subnucleon.useConstituentQuarkProton);
    const double a = param.lattice.L / param.lattice.size;
    const double steps = param.evolution.maxtime / (a * param.run.dtau);
    CHECK(steps == doctest::Approx(static_cast<int>(steps + 0.5)));
    CHECK(param.run.dtau == doctest::Approx(0.1).epsilon(0.05));
}

TEST_CASE("Parameters::readInput: maxtime 0 gives dtau 0.1, not NaN") {
    Parameters param;
    REQUIRE(param.readInput(inputFromText(exampleInputWith("maxtime", "0")))
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
        readErrors(exampleInputWith("muZero", ""));
    REQUIRE(errors.size() == 1);
    CHECK(errors[0] == "test: muZero is required but not given");
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
}

TEST_CASE("Parameters::readInput: values must have the parameter's type") {
    CHECK(anyContains(
        readErrors(exampleInputWith("useJIMWLK", "2")),
        "useJIMWLK '2' is not 0 or 1"));
    CHECK(anyContains(
        readErrors(exampleInputWith("runWithkt", "2")),
        "runWithkt '2' is not 0 or 1"));
    CHECK(anyContains(
        readErrors(exampleInputWith("size", "256.0")),
        "size '256.0' is not an integer"));
    CHECK(anyContains(
        readErrors(exampleInputWith("m", "0.4GeV")),
        "m '0.4GeV' is not a number"));
    CHECK(anyContains(
        readErrors(exampleInputWith("xSnapshotList", "1e-3;1e-4")),
        "xSnapshotList '1e-3;1e-4' is not a comma-separated list"));
}

TEST_CASE("Parameters::readInput: per-value checks") {
    struct Case {
        const char *key;
        const char *value;
        const char *message;
    };
    for (const Case &c : std::vector<Case> {
             {"size", "0", "must be positive"},
             {"alphas_jimwlk", "-0.3", "must not be negative"},
             {"nFlavors", "-1", "must be between 0 and 16"},
             {"size", "255", "must be even"},
             {"omega", "0", "must be positive"},
             {"SubNucleonParamType", "3", "must be one of 0, 1, 2, 4"},
             {"runWith0Min1Avg2MaxQs", "3", "must be one of 0, 1, 2"},
             {"nFlavors", "17", "must be between 0 and 16"},
             {"LambdaQCD", "0", "must be positive"},
             {"c", "-0.2", "must be positive"},
             {"c_jimwlk", "0", "must be positive"},
             {"Lambda_QCD_jimwlk", "0", "must be positive"},
             {"writeWilsonLines", "3", "must be one of 0, 1, 2"},
             {"readInitialWilsonLines", "3", "must be one of 0, 1, 2"},
         }) {
        CAPTURE(c.key);
        CAPTURE(c.value);
        const std::vector<std::string> errors =
            readErrors(exampleInputWith(c.key, c.value));
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
    "read with setWSDeformParams 1") {
    // without them, setWSDeformParams 0 is fine and 1 is not
    std::string text;
    {
        std::istringstream in(exampleInputWith("setWSDeformParams", "0"));
        std::string line;
        while (std::getline(in, line)) {
            if (line.rfind("R_WS ", 0) != 0) text += line + "\n";
        }
    }
    CHECK(readErrors(text).empty());

    std::string text1;
    {
        std::istringstream in(text);
        std::string line;
        while (std::getline(in, line)) {
            text1 += (line.rfind("setWSDeformParams ", 0) == 0)
                         ? "setWSDeformParams 1\n"
                         : line + "\n";
        }
    }
    const std::vector<std::string> errors = readErrors(text1);
    REQUIRE(errors.size() == 1);
    CHECK(errors[0] == "test: R_WS is required but not given");
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
    CHECK(
        reread.jimwlk.x_projectile_jimwlk == param.jimwlk.x_projectile_jimwlk);
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
    "Parameters::readInput: xSnapshotList is only read with saveSnapshots 1") {
    std::string text;
    {
        std::istringstream in(exampleInputWith("xSnapshotList", ""));
        std::string line;
        while (std::getline(in, line)) {
            text += (line.rfind("saveSnapshots ", 0) == 0) ? "saveSnapshots 0\n"
                                                           : line + "\n";
        }
    }
    CHECK(readErrors(text).empty());

    const std::vector<std::string> errors =
        readErrors(exampleInputWith("xSnapshotList", ""));  // saveSnapshots 1
    REQUIRE(errors.size() == 1);
    CHECK(errors[0] == "test: xSnapshotList is required but not given");
}

TEST_CASE(
    "Parameters::readInput: a malformed condition key gives one error, not "
    "an indeterminate set of follow-up errors") {
    const std::vector<std::string> errors =
        readErrors(exampleInputWith("setWSDeformParams", "1.0"));
    REQUIRE(errors.size() == 1);
    CHECK(anyContains(errors, "setWSDeformParams '1.0' is not 0 or 1"));
}

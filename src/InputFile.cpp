// InputFile.cpp is part of the IP-Glasma solver.

#include "InputFile.h"

#include <charconv>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

InputFile::InputFile(const std::string &fileName) : sourceName_(fileName) {
    std::ifstream in(fileName);
    if (!in) {
        errors_.push_back("cannot open input file " + fileName);
        isOpen_ = false;
        return;
    }
    read(in);
}

InputFile::InputFile(std::istream &in, const std::string &sourceName)
    : sourceName_(sourceName) {
    read(in);
}

void InputFile::read(std::istream &in) {
    std::string line;
    int lineNumber = 0;
    while (std::getline(in, line)) {
        lineNumber++;
        // skip a UTF-8 byte order mark written by some editors
        if (lineNumber == 1 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) {
            line.erase(0, 3);
        }
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);

        std::istringstream tokens(line);
        std::string key, value, extra;
        if (!(tokens >> key)) continue;  // blank or comment-only line

        const std::string where =
            sourceName_ + ":" + std::to_string(lineNumber) + ": ";
        if (key == "EndOfFile") {
            if (tokens >> extra) {
                errors_.push_back(
                    where + "unexpected text after EndOfFile ('" + extra
                    + "' ...)");
            }
            break;
        }
        if (!(tokens >> value)) {
            errors_.push_back(where + "no value given for " + key);
            continue;
        }
        if (tokens >> extra) {
            errors_.push_back(
                where + "more than one value given for " + key + " ('" + extra
                + "' ...)");
            continue;
        }
        const auto [it, inserted] =
            entries_.emplace(key, Entry {value, lineNumber});
        if (!inserted) {
            errors_.push_back(
                where + key + " is already set on line "
                + std::to_string(it->second.line));
        }
    }
}

const InputFile::Entry *InputFile::find(const std::string &key) const {
    const auto it = entries_.find(key);
    return (it == entries_.end()) ? nullptr : &it->second;
}

namespace {
/**
 * Converts all of \p text to an integer of type \p T, accepting one
 * leading `+`.
 * \tparam T Integer type to convert to.
 * \param[in] text The value as written in the input.
 * \param[out] out Set to the value on success.
 * \return Whether all of \p text is an integer in the range of \p T.
 */
template <typename T>
bool parseInteger(const std::string &text, T &out) {
    T value {};
    const char *begin = text.data();
    const char *end = begin + text.size();
    // from_chars rejects a leading '+'; accept one, but not "+-5"
    if (begin + 1 < end && *begin == '+' && begin[1] != '-') begin++;
    const auto [ptr, ec] = std::from_chars(begin, end, value);
    if (ec != std::errc() || ptr != end) return false;
    out = value;
    return true;
}
}  // namespace

bool parseValue(const std::string &text, bool &out) {
    if (text != "0" && text != "1") return false;
    out = (text == "1");
    return true;
}

bool parseValue(const std::string &text, int &out) {
    return parseInteger(text, out);
}

bool parseValue(const std::string &text, long long &out) {
    return parseInteger(text, out);
}

bool parseValue(const std::string &text, double &out) {
    if (text.empty()) return false;
    char *end = nullptr;
    const double value = std::strtod(text.c_str(), &end);
    if (end != text.c_str() + text.size() || !std::isfinite(value)) {
        return false;
    }
    out = value;
    return true;
}

bool parseValue(const std::string &text, std::string &out) {
    if (text.empty()) return false;
    out = text;
    return true;
}

bool parseValue(const std::string &text, std::vector<double> &out) {
    if (text == "none") {
        out.clear();
        return true;
    }
    std::vector<double> values;
    std::size_t start = 0;
    while (true) {
        const std::size_t comma = text.find(',', start);
        double value = 0.;
        if (!parseValue(text.substr(start, comma - start), value)) {
            return false;
        }
        values.push_back(value);
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    out = values;
    return true;
}

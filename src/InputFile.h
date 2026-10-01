// InputFile.h is part of the IP-Glasma solver.

#ifndef SRC_INPUTFILE_H_
#define SRC_INPUTFILE_H_

#include <istream>
#include <map>
#include <string>
#include <vector>

/**
 * The `key value` pairs of an IP-Glasma input file, read once.
 *
 * Format: one `key value` pair per line. Blank lines are ignored, and
 * `#` starts a comment that runs to the end of the line (so values,
 * e.g. file paths, cannot contain `#`). A line holding
 * only `EndOfFile` ends the input; everything after it is ignored.
 * A key without a value, a line with more than one value, and a key
 * given twice are errors, collected in errors().
 *
 * Values are stored as text; parseValue() converts them strictly (the
 * whole value must be a valid number of the requested type).
 */
class InputFile {
  public:
    /// One `key value` entry and the line it came from.
    struct Entry {
        std::string value;
        int line;
    };

    /**
     * Reads \p fileName. A missing file is reported in errors().
     * \param[in] fileName Path of the input file.
     */
    explicit InputFile(const std::string &fileName);
    /**
     * Reads the input from \p in (for tests and in-memory input).
     * \param[in] in Stream to read from.
     * \param[in] sourceName Name used in error messages.
     */
    InputFile(std::istream &in, const std::string &sourceName);

    /// Problems found while reading (missing file, malformed lines,
    /// duplicate keys); empty if the file was read cleanly.
    const std::vector<std::string> &errors() const { return errors_; }
    /// Whether the input could be opened at all.
    bool isOpen() const { return isOpen_; }
    /// Name of the file (or stream) this was read from.
    const std::string &sourceName() const { return sourceName_; }
    /// All entries, keyed by parameter name.
    const std::map<std::string, Entry> &entries() const { return entries_; }
    /**
     * Looks up a key.
     * \param[in] key Parameter name.
     * \return The entry, or `nullptr` if the key is not in the input.
     */
    const Entry *find(const std::string &key) const;

  private:
    void read(std::istream &in);

    std::string sourceName_;
    bool isOpen_ = true;
    std::map<std::string, Entry> entries_;
    std::vector<std::string> errors_;
};

/**
 * Strict conversion of an input value. Each overload returns `false`
 * (leaving \p out unchanged) unless all of \p text is a valid value of
 * that type: an `int` or `unsigned long long` must be an integer
 * (no decimal point or exponent), a `bool` must be `0` or `1`, a
 * `double` must be finite, and a
 * list is a comma-separated sequence of doubles.
 */
bool parseValue(const std::string &text, bool &out);
bool parseValue(const std::string &text, int &out);
bool parseValue(const std::string &text, unsigned long long &out);
bool parseValue(const std::string &text, double &out);
bool parseValue(const std::string &text, std::string &out);
bool parseValue(const std::string &text, std::vector<double> &out);

#endif  // SRC_INPUTFILE_H_

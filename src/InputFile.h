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
 * A key without a value, a line with more than one value, a key given
 * twice and text after `EndOfFile` on its line are errors, collected in
 * errors().
 *
 * Values are stored as text; parseValue() converts them strictly (the
 * whole value must be a valid number of the requested type).
 */
class InputFile {
  public:
    /// One `key value` entry and the line it came from.
    struct Entry {
        /// The value, as written in the input (not yet converted).
        std::string value;
        /// Line number of the entry in the input (starting at 1).
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

    /**
     * Returns the problems found while reading.
     * \return One message per missing file, malformed line or duplicate
     * key; empty if the input was read cleanly.
     */
    const std::vector<std::string> &errors() const { return errors_; }
    /**
     * Returns whether the input could be opened at all.
     * \return `false` only if the input file could not be opened.
     */
    bool isOpen() const { return isOpen_; }
    /**
     * Returns the name used for the input in messages.
     * \return The file name (or the stream's name).
     */
    const std::string &sourceName() const { return sourceName_; }
    /**
     * Returns all entries.
     * \return The entries, keyed by parameter name.
     */
    const std::map<std::string, Entry> &entries() const { return entries_; }
    /**
     * Looks up a key.
     * \param[in] key Parameter name.
     * \return The entry, or `nullptr` if the key is not in the input.
     */
    const Entry *find(const std::string &key) const;

  private:
    /**
     * Parses every line of \p in into \c entries_, recording problems
     * in \c errors_; stops at an `EndOfFile` line.
     * \param[in] in Stream to read from.
     */
    void read(std::istream &in);

    /// Name of the file (or stream), used in messages.
    std::string sourceName_;
    /// Whether the input could be opened.
    bool isOpen_ = true;
    /// The parsed entries, keyed by parameter name.
    std::map<std::string, Entry> entries_;
    /// Problems found while reading.
    std::vector<std::string> errors_;
};

// Strict conversion of an input value: each overload returns false (leaving
// out unchanged) unless all of text is a valid value of that type.

/**
 * Converts an on/off value.
 * \param[in] text The value as written in the input.
 * \param[out] out Set to the value if \p text is `0` or `1`.
 * \return Whether \p text is `0` or `1`.
 */
bool parseValue(const std::string &text, bool &out);
/**
 * Converts an integer value (no decimal point or exponent).
 * \param[in] text The value as written in the input.
 * \param[out] out Set to the value on success.
 * \return Whether all of \p text is an integer in the range of `int`.
 */
bool parseValue(const std::string &text, int &out);
/**
 * Converts a non-negative integer value (no decimal point or exponent).
 * \param[in] text The value as written in the input.
 * \param[out] out Set to the value on success.
 * \return Whether all of \p text is an integer in the range of
 * `unsigned long long`.
 */
bool parseValue(const std::string &text, unsigned long long &out);
/**
 * Converts a floating-point value.
 * \param[in] text The value as written in the input.
 * \param[out] out Set to the value on success.
 * \return Whether all of \p text is a finite number.
 */
bool parseValue(const std::string &text, double &out);
/**
 * Takes a string value as it is.
 * \param[in] text The value as written in the input.
 * \param[out] out Set to \p text on success.
 * \return Whether \p text is non-empty.
 */
bool parseValue(const std::string &text, std::string &out);
/**
 * Converts a comma-separated list of floating-point values; `none` is
 * the empty list.
 * \param[in] text The value as written in the input, e.g. `1e-3,1e-4`.
 * \param[out] out Set to the list on success.
 * \return Whether \p text is `none` or one or more finite numbers
 * separated by single commas.
 */
bool parseValue(const std::string &text, std::vector<double> &out);

#endif  // SRC_INPUTFILE_H_

// Setup.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_SETUP_H_
#define SRC_SETUP_H_

#include <string>
#include <vector>

#include "PrettyOstream.h"

/**
 * Reads `key value` pairs from a whitespace-delimited input file (the
 * format main.cpp's parameter file uses, terminated by a line reading
 * `EndOfFile`).
 *
 * `*Find()` methods exit with an error if the file is missing or the
 * key isn't found; the `*FindOptional()` methods still exit if the
 * file is missing, but return a caller-supplied default if the key
 * isn't found instead of exiting.
 */
class Setup {
  private:
    /// Log sink for error messages.
    PrettyOstream messager_;

  public:
    /**
     * Constructs a Setup with no state of its own.
     */
    Setup() {}

    /**
     * Finds the first occurrence of \p st as a token in \p file_name
     * and returns the token immediately following it.
     * \param[in] file_name Path to the input file; exits with an error
     * if it doesn't exist.
     * \param[in] st Key to search for; exits with an error if not
     * found before an `EndOfFile` token.
     * \return The value token following the first match of \p st.
     */
    std::string stringFind(std::string file_name, std::string st);
    /**
     * Same lookup as stringFind(), but returns \p defaultValue instead
     * of exiting if \p st isn't found before an `EndOfFile` token (or
     * before the file itself runs out, if there's no `EndOfFile`
     * token). A key that appears only after `EndOfFile` is therefore
     * not found either, and yields \p defaultValue.
     * \param[in] file_name Path to the input file; exits with an error
     * if it doesn't exist.
     * \param[in] st Key to search for.
     * \param[in] defaultValue Value to return if \p st isn't found.
     * \return The value token following \p st, or \p defaultValue.
     */
    std::string stringFindOptional(
        std::string file_name, std::string st, std::string defaultValue);
    /**
     * Integer form of stringFind(), via dFind().
     * \param[in] file_name Path to the input file; exits with an error
     * if it doesn't exist.
     * \param[in] st Key to search for; exits with an error if not
     * found.
     * \return The value following \p st, truncated towards zero.
     */
    int iFind(std::string file_name, std::string st);
    /**
     * Integer form of stringFindOptional().
     * \param[in] file_name Path to the input file; exits with an error
     * if it doesn't exist.
     * \param[in] st Key to search for.
     * \param[in] defaultValue Value to return if \p st isn't found.
     * \return The value following \p st (parsed as a `double` then
     * truncated towards zero), or \p defaultValue.
     */
    int iFindOptional(std::string file_name, std::string st, int defaultValue);
    /**
     * Unsigned 64-bit integer form of stringFind(), via dFind(),
     * rounded to the nearest integer rather than truncated (unlike
     * iFind()).
     * \param[in] file_name Path to the input file; exits with an error
     * if it doesn't exist.
     * \param[in] st Key to search for; exits with an error if not
     * found.
     * \return The value following \p st, rounded to the nearest
     * integer.
     */
    unsigned long long int uLLIFind(std::string file_name, std::string st);
    /**
     * Floating-point form of stringFind().
     * \param[in] file_name Path to the input file; exits with an error
     * if it doesn't exist.
     * \param[in] st Key to search for; exits with an error if not
     * found.
     * \return The value following \p st, parsed as a `double` (`0.0`
     * if it doesn't parse as a number).
     */
    double dFind(std::string file_name, std::string st);
    /**
     * Same lookup as dFind(), but returns \p defaultValue instead of
     * exiting if \p st isn't found before an `EndOfFile` token (or
     * before the file itself runs out, if there's no `EndOfFile`
     * token). A key that appears only after `EndOfFile` is therefore
     * not found either, and yields \p defaultValue.
     * \param[in] file_name Path to the input file; exits with an error
     * if it doesn't exist.
     * \param[in] st Key to search for.
     * \param[in] defaultValue Value to return if \p st isn't found.
     * \return The value token following \p st, parsed as a `double`
     * (`0.0` if it doesn't parse as a number), or \p defaultValue.
     */
    double dFindOptional(
        std::string file_name, std::string st, double defaultValue);
    /**
     * Checks whether a file can be opened for reading.
     * \param[in] file_name Path to check.
     * \return `1` if \p file_name can be opened for reading, `0`
     * otherwise.
     */
    int isFile(std::string file_name);
    /**
     * Finds the first line containing \p st and parses it as a
     * comma-separated list of numbers (skipping the line's first
     * whitespace-delimited token, e.g. a leading key name).
     * \param[in] file_name Path to the input file; exits with an error
     * if it doesn't exist.
     * \param[in] st Substring to search for within each line.
     * \return The parsed numbers, in order; empty if no line contains
     * \p st.
     */
    std::vector<double> listFind(std::string file_name, std::string st);
};

#endif  // SRC_SETUP_H_

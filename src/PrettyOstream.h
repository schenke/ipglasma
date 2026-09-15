// Copyright Chun Shen @ 2017
// This class is inspired by the JetScapeLogger class written by Joern Putschke

#define BOLD "\033[1m"      // Bold
#define BLACK "\033[30m"    // Black
#define RED "\033[31m"      // Red
#define GREEN "\033[32m"    // Green
#define YELLOW "\033[33m"   // Yellow
#define BLUE "\033[34m"     // Blue
#define MAGENTA "\033[35m"  // Magenta
#define CYAN "\033[36m"     // Cyan
#define WHITE "\033[37m"    // White
#define ORANGE \
    "\033[38;5;208m"     // Orange (xterm 256-color; no plain ANSI
                         // code exists for orange)
#define RESET "\033[0m"  // reset

#ifndef SRC_PRETTYOSTREAM_H_
#define SRC_PRETTYOSTREAM_H_

#include <sstream>
#include <string>

/**
 * Buffered, colorized console logger: accumulates a message via
 * `operator<<` like an `ostream`, then flush() dispatches it to
 * info()/debug()/warning()/error() by category.
 *
 * Every category-specific method (info()/debug()/warning()/error())
 * locks a shared mutex around its actual write to `cout`, so instances
 * can be constructed one-per-thread (e.g. inside an
 * `#pragma omp parallel` region, since this class itself has no
 * internal state to race on beyond that shared terminal) without
 * interleaving another thread's message.
 */
class PrettyOstream {
  private:
    /// Accumulates the message being built via `operator<<`, until the
    /// next flush().
    std::ostringstream messageStream_;

  public:
    /**
     * Constructs a PrettyOstream with an empty message buffer.
     */
    PrettyOstream();
    /**
     * Destroys this PrettyOstream (nothing to release; any
     * un-flushed, buffered message is silently discarded).
     */
    ~PrettyOstream();

    /**
     * Dispatches the buffered message to info()/debug()/warning()/
     * error() based on \p type, then clears the buffer. An
     * unrecognized \p type silently discards the buffered message
     * without printing anything.
     * \param[in] type Category name, case-insensitively one of
     * `"info"`, `"debug"`, `"warning"`, `"error"`.
     */
    void flush(std::string type);

    /**
     * Prints an info-level message (uncolored) with a memory-usage
     * prefix.
     * \param[in] message Message to print.
     */
    void info(std::string message);

    /**
     * Prints a debug-level message (cyan) with a memory-usage prefix.
     * \param[in] message Message to print.
     */
    void debug(std::string message);

    /**
     * Prints a warning-level message (bold orange).
     * \param[in] message Message to print.
     */
    void warning(std::string message);

    /**
     * Prints an error-level message (bold red).
     * \param[in] message Message to print.
     */
    void error(std::string message);

    /**
     * Reads this process' peak resident memory usage via `getrusage()`.
     * \return `"<value> MB"` (4 significant digits), or an empty
     * string if `getrusage()` fails.
     */
    std::string getMemoryUsage();

    /**
     * Appends a value to the buffered message, like `ostream`'s
     * `operator<<`.
     * \param[in] value Value to append; anything `ostringstream`
     * accepts.
     * \return `*this`, for chaining.
     */
    template <typename T>
    PrettyOstream &operator<<(T const &value) {
        messageStream_ << value;
        return (*this);
    }
};

#endif  // SRC_PRETTYOSTREAM_H_

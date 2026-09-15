#ifndef SRC_INSTRUMENTATION_H_
#define SRC_INSTRUMENTATION_H_

#include <chrono>
#include <cstdint>
#include <map>
#include <mutex>
#include <string>

class Lattice;

/**
 * Opt-in performance profiling (per-phase wall-clock timing, enabled
 * via `IPGLASMA_PROFILE=1`) and a bit-reproducibility fingerprint
 * (enabled via `IPGLASMA_FINGERPRINT=1`, used to verify a refactor
 * doesn't change numeric output). Both are zero-cost when disabled
 * (the default).
 */
namespace ipg {

/**
 * Accumulated wall-clock time and call count for one named profiling
 * phase within an event.
 */
struct PhaseStat {
    /**
     * Constructs a PhaseStat with zero elapsed time and zero calls.
     */
    PhaseStat() : seconds(0.0), calls(0) {}
    /// Total accumulated wall-clock time for this phase [s].
    double seconds;
    /// Number of times this phase's timer was added to.
    std::uint64_t calls;
};

/**
 * Process-wide singleton collecting per-phase timing for the current
 * event and appending it to a per-rank TSV file when the event ends.
 * Enabled by the `IPGLASMA_PROFILE` environment variable (any value
 * other than empty/`0`/`false`/`off`/`no`, case-insensitively);
 * disabled by default, in which case every method is a cheap no-op.
 */
class Profiler {
  public:
    /**
     * Returns the process-wide Profiler instance.
     * \return Reference to the singleton.
     */
    static Profiler &instance();

    /**
     * Records this rank's MPI rank and re-reads the enabled/output-
     * directory environment variables (call once at startup, after MPI
     * is initialized).
     * \param[in] rank This process' MPI rank, used in output filenames.
     */
    void initialize(int rank);
    /**
     * Reports whether profiling is currently enabled.
     * \return `true` if `IPGLASMA_PROFILE` was set to an enabling value
     * at construction (or the last initialize() call).
     */
    bool enabled() const;
    /**
     * Starts timing a new event, clearing any phase statistics left
     * over from a previous event. A no-op if profiling is disabled.
     * \param[in] event_id Identifier recorded alongside every phase row
     * for this event (e.g. the event-loop index).
     */
    void beginEvent(int event_id);
    /**
     * Adds elapsed time to a named phase's running total for the
     * current event. A no-op if profiling is disabled or no event is
     * active.
     * \param[in] phase Phase name (dot-separated, e.g.
     * `"fft.total"`); merged with any prior calls under the same name
     * this event.
     * \param[in] seconds Elapsed time to add [s].
     */
    void add(const std::string &phase, double seconds);
    /**
     * Ends the current event, appending one row per accumulated phase
     * (plus a synthetic `"event.total"` row) to
     * `ipglasma_profile_rank<rank>.tsv` in `IPGLASMA_PROFILE_DIR`
     * (default: the current directory), writing a header line first if
     * the file is new/empty. A no-op if profiling is disabled or no
     * event is active.
     */
    void endEvent();

  private:
    Profiler();
    Profiler(const Profiler &);
    Profiler &operator=(const Profiler &);

    /// Whether `IPGLASMA_PROFILE` was set to an enabling value.
    bool enabled_;
    /// Whether beginEvent() has been called without a matching
    /// endEvent() yet.
    bool event_active_;
    /// This process' MPI rank, set by initialize().
    int rank_;
    /// Identifier of the currently (or most recently) active event.
    int event_id_;
    /// Directory the per-rank TSV file is written into.
    std::string output_dir_;
    /// Wall-clock time beginEvent() was called, used to compute
    /// `"event.total"` in endEvent().
    std::chrono::steady_clock::time_point event_start_;
    /// Running per-phase totals for the current event, keyed by phase
    /// name.
    std::map<std::string, PhaseStat> stats_;
    /// Serializes access to every member above, since add() is called
    /// concurrently from OpenMP-parallel regions.
    mutable std::mutex mutex_;
};

/**
 * RAII wall-clock timer: times its own scope and adds the elapsed time
 * to Profiler::instance() under a given phase name on destruction.
 * Reads Profiler::enabled() once at construction, so it costs one
 * branch (and nothing else) when profiling is disabled. Normally used
 * via the IPG_PROFILE_SCOPE() macro rather than named directly.
 */
class ScopedTimer {
  public:
    /**
     * Starts timing, if profiling is enabled.
     * \param[in] phase Phase name to add elapsed time to on
     * destruction; must outlive nothing (copied into \c phase_).
     */
    explicit ScopedTimer(const char *phase);
    /**
     * Starts timing, if profiling is enabled.
     * \param[in] phase Phase name to add elapsed time to on
     * destruction.
     */
    explicit ScopedTimer(const std::string &phase);
    /**
     * Adds this scope's elapsed time to Profiler::instance() under \c
     * phase_, if profiling was enabled at construction.
     */
    ~ScopedTimer();

  private:
    ScopedTimer(const ScopedTimer &);
    ScopedTimer &operator=(const ScopedTimer &);

    /// Whether profiling was enabled when this timer was constructed.
    bool active_;
    /// Phase name to report elapsed time under.
    std::string phase_;
    /// Wall-clock time this timer was constructed.
    std::chrono::steady_clock::time_point start_;
};

/**
 * Current wall-clock time, for manual timing spans that don't fit
 * ScopedTimer's scope-based lifetime (e.g. spans crossing an OpenMP
 * region boundary).
 * \return Seconds since an unspecified epoch; only differences between
 * two calls are meaningful.
 */
double wallSeconds();
/**
 * Reports whether lattice-fingerprint output is enabled.
 * \return `true` if the `IPGLASMA_FINGERPRINT` environment variable is
 * set to an enabling value (see Profiler's `IPGLASMA_PROFILE` for the
 * exact rule).
 */
bool fingerprintEnabled();
/**
 * Writes a compact, order-independent-within-value but
 * position-dependent digest (count, non-finite count, an FNV-1a hash
 * of every value's raw bits, mean, RMS, min, max) of every physically
 * meaningful lattice field to
 * `ipglasma_fingerprint_rank<rank>.tsv` in `IPGLASMA_PROFILE_DIR`,
 * plus one combined hash of all fields together -- cheap enough to run
 * every event, and precise enough that two runs' fingerprints matching
 * is strong evidence their full lattice state matches bit-for-bit
 * (used to verify a refactor changed no numeric output). A no-op if
 * fingerprinting is disabled or \p lat is `NULL`.
 * \param[in] lat Lattice to digest.
 * \param[in] rank This process' MPI rank, used in the output filename.
 * \param[in] event_id Identifier recorded alongside every row.
 */
void writeLatticeFingerprint(Lattice *lat, int rank, int event_id);

}  // namespace ipg

#define IPG_JOIN_IMPL(a, b) a##b
#define IPG_JOIN(a, b) IPG_JOIN_IMPL(a, b)
/**
 * \def IPG_PROFILE_SCOPE(phase_name)
 * Times the remainder of the enclosing scope and adds it to
 * ipg::Profiler::instance() under \p phase_name on scope exit, via a
 * uniquely-named ipg::ScopedTimer local.
 */
#define IPG_PROFILE_SCOPE(phase_name) \
    ipg::ScopedTimer IPG_JOIN(ipg_scoped_timer_, __LINE__)(phase_name)

#endif  // SRC_INSTRUMENTATION_H_

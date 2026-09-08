#ifndef IPGLASMA_INSTRUMENTATION_H
#define IPGLASMA_INSTRUMENTATION_H

#include <chrono>
#include <cstdint>
#include <map>
#include <mutex>
#include <string>

class Lattice;

namespace ipg {

struct PhaseStat {
    PhaseStat() : seconds(0.0), calls(0) {}
    double seconds;
    std::uint64_t calls;
};

class Profiler {
  public:
    static Profiler &instance();

    void initialize(int rank);
    bool enabled() const;
    void beginEvent(int event_id);
    void add(const std::string &phase, double seconds);
    void endEvent();

  private:
    Profiler();
    Profiler(const Profiler &);
    Profiler &operator=(const Profiler &);

    bool enabled_;
    bool event_active_;
    int rank_;
    int event_id_;
    std::string output_dir_;
    std::chrono::steady_clock::time_point event_start_;
    std::map<std::string, PhaseStat> stats_;
    mutable std::mutex mutex_;
};

class ScopedTimer {
  public:
    explicit ScopedTimer(const char *phase);
    explicit ScopedTimer(const std::string &phase);
    ~ScopedTimer();

  private:
    ScopedTimer(const ScopedTimer &);
    ScopedTimer &operator=(const ScopedTimer &);

    bool active_;
    std::string phase_;
    std::chrono::steady_clock::time_point start_;
};

double wallSeconds();
bool fingerprintEnabled();
void writeLatticeFingerprint(Lattice *lat, int rank, int event_id);

}  // namespace ipg

#define IPG_JOIN_IMPL(a, b) a##b
#define IPG_JOIN(a, b) IPG_JOIN_IMPL(a, b)
#define IPG_PROFILE_SCOPE(phase_name) \
    ipg::ScopedTimer IPG_JOIN(ipg_scoped_timer_, __LINE__)(phase_name)

#endif

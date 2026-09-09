#include "Instrumentation.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "Cell.h"
#include "Lattice.h"
#include "Matrix.h"

namespace {

bool envEnabled(const char *name) {
    const char *value = std::getenv(name);
    if (value == NULL || value[0] == '\0') return false;
    const std::string text(value);
    return text != "0" && text != "false" && text != "FALSE" && text != "off"
           && text != "OFF" && text != "no" && text != "NO";
}

std::string envOrDefault(const char *name, const char *fallback) {
    const char *value = std::getenv(name);
    if (value == NULL || value[0] == '\0') return std::string(fallback);
    return std::string(value);
}

std::string joinPath(const std::string &directory, const std::string &name) {
    if (directory.empty() || directory == ".") return name;
    if (directory[directory.size() - 1] == '/') return directory + name;
    return directory + "/" + name;
}

bool fileNeedsHeader(const std::string &path) {
    std::ifstream input(path.c_str(), std::ios::binary | std::ios::ate);
    return !input || input.tellg() == std::streampos(0);
}

class FieldDigest {
  public:
    FieldDigest()
        : hash_(14695981039346656037ULL),
          count_(0),
          finite_count_(0),
          nonfinite_(0),
          sum_(0.0L),
          sumsq_(0.0L),
          min_(std::numeric_limits<double>::infinity()),
          max_(-std::numeric_limits<double>::infinity()) {}

    void add(double value) {
        std::uint64_t bits = 0;
        static_assert(sizeof(bits) == sizeof(value), "unexpected double size");
        std::memcpy(&bits, &value, sizeof(bits));
        for (int shift = 0; shift < 64; shift += 8) {
            const unsigned char byte =
                static_cast<unsigned char>((bits >> shift) & 0xffU);
            hash_ ^= static_cast<std::uint64_t>(byte);
            hash_ *= 1099511628211ULL;
        }

        ++count_;
        const std::uint64_t exponent_mask = 0x7ff0000000000000ULL;
        if ((bits & exponent_mask) == exponent_mask) {
            ++nonfinite_;
            return;
        }
        ++finite_count_;
        const long double x = static_cast<long double>(value);
        sum_ += x;
        sumsq_ += x * x;
        if (value < min_) min_ = value;
        if (value > max_) max_ = value;
    }

    std::uint64_t hash() const { return hash_; }
    std::uint64_t count() const { return count_; }
    std::uint64_t nonfinite() const { return nonfinite_; }
    double mean() const {
        return finite_count_ == 0 ? std::numeric_limits<double>::quiet_NaN()
                                  : static_cast<double>(sum_ / finite_count_);
    }
    double rms() const {
        return finite_count_ == 0
                   ? std::numeric_limits<double>::quiet_NaN()
                   : std::sqrt(static_cast<double>(sumsq_ / finite_count_));
    }
    double minimum() const {
        return finite_count_ == 0 ? std::numeric_limits<double>::quiet_NaN()
                                  : min_;
    }
    double maximum() const {
        return finite_count_ == 0 ? std::numeric_limits<double>::quiet_NaN()
                                  : max_;
    }

  private:
    std::uint64_t hash_;
    std::uint64_t count_;
    std::uint64_t finite_count_;
    std::uint64_t nonfinite_;
    long double sum_;
    long double sumsq_;
    double min_;
    double max_;
};

void mixCombined(
    std::uint64_t &combined, const std::string &label,
    const std::uint64_t field_hash) {
    for (std::size_t i = 0; i < label.size(); ++i) {
        combined ^= static_cast<unsigned char>(label[i]);
        combined *= 1099511628211ULL;
    }
    for (int shift = 0; shift < 64; shift += 8) {
        combined ^= static_cast<unsigned char>((field_hash >> shift) & 0xffU);
        combined *= 1099511628211ULL;
    }
}

void writeDigestRow(
    std::ofstream &output, int rank, int event_id, const std::string &label,
    const FieldDigest &digest, std::uint64_t &combined) {
    output << rank << '\t' << event_id << '\t' << label << '\t'
           << digest.count() << '\t' << digest.nonfinite() << '\t' << "0x"
           << std::hex << std::setw(16) << std::setfill('0') << digest.hash()
           << std::dec << std::setfill(' ') << '\t' << std::setprecision(17)
           << digest.mean() << '\t' << digest.rms() << '\t' << digest.minimum()
           << '\t' << digest.maximum() << '\n';
    mixCombined(combined, label, digest.hash());
}

void digestScalarField(
    std::ofstream &output, Lattice *lat, int rank, int event_id,
    const std::string &label, double (Cell::*getter)() const,
    std::uint64_t &combined) {
    FieldDigest digest;
    for (int pos = 0; pos < lat->getSize(); ++pos) {
        digest.add((lat->cells[pos]->*getter)());
    }
    writeDigestRow(output, rank, event_id, label, digest, combined);
}

void digestMatrixField(
    std::ofstream &output, const std::vector<Matrix> &field, int rank,
    int event_id, const std::string &label, std::uint64_t &combined) {
    FieldDigest digest;
    for (std::size_t pos = 0; pos < field.size(); ++pos) {
        const std::complex<double> *values = field[pos].data();
        for (int i = 0; i < 9; ++i) {
            digest.add(values[i].real());
            digest.add(values[i].imag());
        }
    }
    writeDigestRow(output, rank, event_id, label, digest, combined);
}

}  // namespace

namespace ipg {

Profiler &Profiler::instance() {
    static Profiler profiler;
    return profiler;
}

Profiler::Profiler()
    : enabled_(envEnabled("IPGLASMA_PROFILE")),
      event_active_(false),
      rank_(0),
      event_id_(-1),
      output_dir_(envOrDefault("IPGLASMA_PROFILE_DIR", ".")) {}

void Profiler::initialize(int rank) {
    std::lock_guard<std::mutex> lock(mutex_);
    rank_ = rank;
    enabled_ = envEnabled("IPGLASMA_PROFILE");
    output_dir_ = envOrDefault("IPGLASMA_PROFILE_DIR", ".");
}

bool Profiler::enabled() const { return enabled_; }

void Profiler::beginEvent(int event_id) {
    if (!enabled_) return;
    std::lock_guard<std::mutex> lock(mutex_);
    stats_.clear();
    event_id_ = event_id;
    event_start_ = std::chrono::steady_clock::now();
    event_active_ = true;
}

void Profiler::add(const std::string &phase, double seconds) {
    if (!enabled_) return;
    std::lock_guard<std::mutex> lock(mutex_);
    if (!event_active_) return;
    PhaseStat &stat = stats_[phase];
    stat.seconds += seconds;
    ++stat.calls;
}

void Profiler::endEvent() {
    if (!enabled_) return;

    std::lock_guard<std::mutex> lock(mutex_);
    if (!event_active_) return;

    const std::chrono::duration<double> elapsed =
        std::chrono::steady_clock::now() - event_start_;
    PhaseStat &total = stats_["event.total"];
    total.seconds += elapsed.count();
    ++total.calls;

    std::ostringstream filename;
    filename << "ipglasma_profile_rank" << rank_ << ".tsv";
    const std::string path = joinPath(output_dir_, filename.str());
    const bool header = fileNeedsHeader(path);
    std::ofstream output(path.c_str(), std::ios::app);
    if (output) {
        if (header) {
            output << "rank\tevent\tphase\tseconds\tcalls\tpercent_event\n";
        }
        const double event_seconds = total.seconds;
        for (std::map<std::string, PhaseStat>::const_iterator it =
                 stats_.begin();
             it != stats_.end(); ++it) {
            const double percent =
                event_seconds > 0.0 ? 100.0 * it->second.seconds / event_seconds
                                    : 0.0;
            output << rank_ << '\t' << event_id_ << '\t' << it->first << '\t'
                   << std::setprecision(12) << it->second.seconds << '\t'
                   << it->second.calls << '\t' << percent << '\n';
        }
    }

    stats_.clear();
    event_active_ = false;
}

ScopedTimer::ScopedTimer(const char *phase)
    : active_(Profiler::instance().enabled()) {
    if (active_) {
        phase_ = phase;
        start_ = std::chrono::steady_clock::now();
    }
}

ScopedTimer::ScopedTimer(const std::string &phase)
    : active_(Profiler::instance().enabled()) {
    if (active_) {
        phase_ = phase;
        start_ = std::chrono::steady_clock::now();
    }
}

ScopedTimer::~ScopedTimer() {
    if (!active_) return;
    const std::chrono::duration<double> elapsed =
        std::chrono::steady_clock::now() - start_;
    Profiler::instance().add(phase_, elapsed.count());
}

double wallSeconds() {
    const std::chrono::duration<double> elapsed =
        std::chrono::steady_clock::now().time_since_epoch();
    return elapsed.count();
}

bool fingerprintEnabled() { return envEnabled("IPGLASMA_FINGERPRINT"); }

void writeLatticeFingerprint(Lattice *lat, int rank, int event_id) {
    if (!fingerprintEnabled() || lat == NULL) return;

    const std::string output_dir = envOrDefault("IPGLASMA_PROFILE_DIR", ".");
    std::ostringstream filename;
    filename << "ipglasma_fingerprint_rank" << rank << ".tsv";
    const std::string path = joinPath(output_dir, filename.str());
    const bool header = fileNeedsHeader(path);
    std::ofstream output(path.c_str(), std::ios::app);
    if (!output) return;

    if (header) {
        output << "rank\tevent\tfield\tcount\tnonfinite\thash_fnv1a64"
                  "\tmean\trms\tmin\tmax\n";
    }

    std::uint64_t combined = 14695981039346656037ULL;

    digestScalarField(
        output, lat, rank, event_id, "epsilon", &Cell::getEpsilon, combined);
    digestScalarField(
        output, lat, rank, event_id, "Ttautau", &Cell::getTtautau, combined);
    digestScalarField(
        output, lat, rank, event_id, "Txx", &Cell::getTxx, combined);
    digestScalarField(
        output, lat, rank, event_id, "Tyy", &Cell::getTyy, combined);
    digestScalarField(
        output, lat, rank, event_id, "Txy", &Cell::getTxy, combined);
    digestScalarField(
        output, lat, rank, event_id, "Tetaeta", &Cell::getTetaeta, combined);
    digestScalarField(
        output, lat, rank, event_id, "Ttaux", &Cell::getTtaux, combined);
    digestScalarField(
        output, lat, rank, event_id, "Ttauy", &Cell::getTtauy, combined);
    digestScalarField(
        output, lat, rank, event_id, "Ttaueta", &Cell::getTtaueta, combined);
    digestScalarField(
        output, lat, rank, event_id, "Txeta", &Cell::getTxeta, combined);
    digestScalarField(
        output, lat, rank, event_id, "Tyeta", &Cell::getTyeta, combined);
    digestScalarField(
        output, lat, rank, event_id, "utau", &Cell::getutau, combined);
    digestScalarField(
        output, lat, rank, event_id, "ux", &Cell::getux, combined);
    digestScalarField(
        output, lat, rank, event_id, "uy", &Cell::getuy, combined);
    digestScalarField(
        output, lat, rank, event_id, "ueta", &Cell::getueta, combined);
    digestScalarField(
        output, lat, rank, event_id, "g2mu2A", &Cell::getg2mu2A, combined);
    digestScalarField(
        output, lat, rank, event_id, "g2mu2B", &Cell::getg2mu2B, combined);

    digestMatrixField(output, lat->Ux, rank, event_id, "Ux.reim", combined);
    digestMatrixField(output, lat->Uy, rank, event_id, "Uy.reim", combined);
    digestMatrixField(output, lat->U, rank, event_id, "E1.reim", combined);
    digestMatrixField(output, lat->U2, rank, event_id, "E2.reim", combined);
    digestMatrixField(output, lat->Uy2, rank, event_id, "phi.reim", combined);
    digestMatrixField(output, lat->Ux2, rank, event_id, "pi.reim", combined);

    output << rank << '\t' << event_id << "\tALL_FIELDS\t0\t0\t0x" << std::hex
           << std::setw(16) << std::setfill('0') << combined << std::dec
           << std::setfill(' ') << "\tnan\tnan\tnan\tnan\n";
}

}  // namespace ipg

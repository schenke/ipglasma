#ifndef SRC_LATTICE_H_
#define SRC_LATTICE_H_

#include <string>
#include <vector>

#include "Cell.h"
#include "LatticeIndex.h"
#include "Matrix.h"
#include "Parameters.h"
#include "PrettyOstream.h"

enum class NucleusRole;

/**
 * The transverse lattice: every fundamental SU(3) matrix field (stored
 * structure-of-arrays, one contiguous `std::vector<Matrix>` per field)
 * plus the scalar Cell observables and precomputed nearest-neighbor
 * index tables.
 *
 * \note posmX/pospX/posmY/pospY/posmXpY/pospXmY use *clamped*
 * (Neumann-like) boundaries -- at an edge, the neighbor index is the
 * site itself, not a periodic wraparound. Code elsewhere that needs a
 * genuinely periodic neighbor (e.g. GaugeFix::fftChi()) computes its
 * own index instead of using these arrays.
 */
class Lattice {
  private:
    /// Total number of sites, `length*length` (despite the name, *not*
    /// the lattice side length -- see getSize()).
    int size_;
    /// Lattice side length (sites per dimension).
    int length_;
    /// Log sink for progress/error messages.
    PrettyOstream messager_;

    /**
     * Writes one \f$N\times N\f$ Matrix field to a text file, one line
     * per site (`ix iy MatrixToString()`), used by writeSU3Matrices()
     * for both its Phi and Pi passes.
     * \param[in] file_name Output file path.
     * \param[in] field Per-site matrix field to write (`Uy2` or
     * `Ux2`).
     * \param[in] N Lattice side length.
     */
    void writeMatrixArrayText(
        const std::string &file_name, std::vector<Matrix> &field, int N);

  public:
    /**
     * Allocates every field at identity/zero and builds the
     * nearest-neighbor index tables for a square lattice.
     * \param[in] param Simulation parameters; only `lattice.L` and
     * `run.MPIRank` are used (the side length comes from \p length).
     * \param[in] length Lattice side length; total site count is
     * \p length squared.
     */
    Lattice(Parameters *param, int length);
    ~Lattice() = default;
    // Copying would either deep-copy several large per-site vectors for
    // no current caller, or (if left shallow) leave `cells` dangling
    // into the wrong `cellStorage` -- disable copying (nothing needs
    // it).
    Lattice(const Lattice &) = delete;
    Lattice &operator=(const Lattice &) = delete;

    /**
     * Returns the total number of lattice sites.
     * \note Despite the name, this is `length*length` from the
     * constructor, i.e. the *total* site count, not the lattice side
     * length.
     * \return Total number of sites.
     */
    int getSize() const { return size_; }
    /**
     * Returns the lattice side length.
     * \return Number of sites per dimension.
     */
    int getLength() const { return length_; }
    /**
     * Site index of the transverse position (\p ix, \p iy); see
     * latticeIndex().
     * \param[in] ix Site index in \f$x\f$.
     * \param[in] iy Site index in \f$y\f$.
     * \return The index into the per-site vectors.
     */
    int positionFromXY(int ix, int iy) const {
        return latticeIndex(ix, iy, length_);
    }
    /**
     * The \f$x\f$ index of a site; inverse of positionFromXY().
     * \param[in] pos Index into the per-site vectors.
     * \return The site's \f$x\f$ index.
     */
    int xFromPosition(int pos) const { return latticeX(pos, length_); }
    /**
     * The \f$y\f$ index of a site; inverse of positionFromXY().
     * \param[in] pos Index into the per-site vectors.
     * \return The site's \f$y\f$ index.
     */
    int yFromPosition(int pos) const { return latticeY(pos, length_); }

    // Fundamental matrix lattice fields. Logical aliases are:
    // U/E1, U2/E2, Ux1/g, Uy1/Uplaq, Ux2/pi, Uy2/phi.
    /// Projectile Wilson line / electric-field variable \f$E_1\f$
    /// (reused for both roles across the initialization/evolution
    /// stages).
    std::vector<Matrix> U;
    /// Target Wilson line / electric-field variable \f$E_2\f$.
    std::vector<Matrix> U2;
    /// Transverse gauge link in the \f$x\f$ direction, \f$U_x\f$.
    std::vector<Matrix> Ux;
    /// Transverse gauge link in the \f$y\f$ direction, \f$U_y\f$.
    std::vector<Matrix> Uy;
    /// Scratch/gauge-transformation field \f$g(x)\f$ (see
    /// GaugeFix::fftChi()), reusing the \c Ux1 storage.
    std::vector<Matrix> Ux1;
    /// Scratch plaquette field \f$U_{\text{plaq}}\f$ (see
    /// Evolution.cpp's `tmunuPlaquetteTeam`), reusing the \c Uy1
    /// storage.
    std::vector<Matrix> Uy1;
    /// Momentum conjugate to the classical field, \f$\pi\f$ (also
    /// written `Pi` by writeSU3Matrices()), reusing the \c Ux2 storage.
    std::vector<Matrix> Ux2;
    /// Classical field variable \f$\phi\f$ (also written `Phi` by
    /// writeSU3Matrices()), reusing the \c Uy2 storage.
    std::vector<Matrix> Uy2;

    /// Pointer-per-cell view over \c cellStorage (the actual owner);
    /// exists so callers can hold/reseat per-cell pointers without
    /// invalidating \c cellStorage itself.
    std::vector<Cell *> cells;
    /// Owns the per-site scalar observables \c cells points into.
    std::vector<Cell> cellStorage;

    /// Site index one step in \f$-x\f$, clamped at the boundary (see
    /// the class-level \note).
    std::vector<int> posmX;
    /// Site index one step in \f$+x\f$, clamped at the boundary.
    std::vector<int> pospX;
    /// Site index one step in \f$-y\f$, clamped at the boundary.
    std::vector<int> posmY;
    /// Site index one step in \f$+y\f$, clamped at the boundary.
    std::vector<int> pospY;

    /**
     * Writes this instance's projectile or target Wilson line
     * (\c U or \c U2) to disk, in text or binary format depending on
     * \p param's `wilsonLines.writeWilsonLines`. The output file name is
     * generated by generateWilsonLineDataFileName().
     * \param[in] param Simulation parameters; `lattice.size`, `lattice.L`,
     * `event.eventId`, `random.seed`, `run.MPISize`,
     * `wilsonLines.wilsonLinePath`, `wilsonLines.writeWilsonLines` and (binary
     * format only) `colorCharge.RapidityA`/`colorCharge.RapidityB` are used.
     * \param[in] nucleus Which of \c U (Projectile) / \c U2 (Target) to
     * write.
     * \param[in] x If non-negative, included in the generated file name
     * (see generateWilsonLineDataFileName()); the default `-1` omits it.
     *
     * Format `1` (text): one line per site, `ix iy
     * MatrixToString()`. Format `2` (binary): a
     * `{N, Nc, L, a, rapidity}` header followed by every site's 9
     * complex components as consecutive `{real, imag}` `double` pairs,
     * indexed `N*ix+iy` (matching every other `U`/`U2` indexing in
     * the codebase; see Init::readVFromFile for the matching reader).
     * Exits with an error for any other `wilsonLines.writeWilsonLines` value,
     * or if the binary write fails.
     */
    void writeWilsonLines(
        Parameters *param, NucleusRole nucleus, double x = -1);
    /**
     * Writes this instance's \c Uy2 (`Phi`) and \c Ux2 (`Pi`) fields to
     * two separate text files via writeMatrixArrayText().
     * \param[in] fileprefix Prefix prepended to both output file names.
     * \param[in] param Simulation parameters; `lattice.size`,
     * `event.eventId`, `random.seed` and `run.MPISize` are used.
     */
    void writeSU3Matrices(std::string fileprefix, Parameters *param);
    /// Site index one step in \f$-x\f$ and one step in \f$+y\f$,
    /// clamped at the boundary.
    std::vector<int> posmXpY;
    /// Site index one step in \f$+x\f$ and one step in \f$-y\f$,
    /// clamped at the boundary.
    std::vector<int> pospXmY;

    /**
     * Checks whether a `writeWilsonLines`/`readInitialWilsonLines`
     * input value is one of the formats writeWilsonLines() actually
     * supports.
     * \param[in] format Value to check.
     * \return `true` if \p format is `1` (text) or `2` (binary).
     */
    static bool IsValidWilsonLineDataFormat(const int format) {
        return format == 1 || format == 2;
    }

    /**
     * Generates the deterministic file name writeWilsonLines() writes to
     * (and Init::readVFromFile() reads from): `<wilsonLines.wilsonLinePath>/
     * WilsonLine[_x_<x>]_<event/nucleus/seed/rank-derived index>[.txt]`.
     * \param[in] param Simulation parameters; `wilsonLines.wilsonLinePath`,
     * `event.eventId`, `random.seed`, `run.MPISize` and (when \p format is
     * negative) `wilsonLines.writeWilsonLines` are used.
     * \param[in] x If non-negative, embedded in the file name
     * (`_x_<x>`, in scientific notation); a negative value omits it.
     * \param[in] nucleus Which nucleus this file belongs to (offsets
     * the numeric index so the projectile and target never collide).
     * \param[in] format Overrides `param->wilsonLines.writeWilsonLines` for
     * deciding whether to append the `.txt` extension (`1`: text,
     * anything else: binary, no extension); a negative value (the
     * default) uses `param->wilsonLines.writeWilsonLines` instead.
     * \return The generated file path.
     */
    static std::string generateWilsonLineDataFileName(
        Parameters *param, const double x, NucleusRole nucleus,
        int format = -1);
};

/**
 * A pair of scratch \f$N\times N\f$ Matrix fields, allocated at
 * identity, used as ping-pong buffers by algorithms that need
 * temporary per-site matrix storage without touching a Lattice's own
 * fields.
 */
class BufferLattice {
  private:
    /// Total number of sites, `length*length` (see
    /// Lattice::size_/getSize() for the same naming caveat).
    int size_;

  public:
    /**
     * Allocates both buffers at identity.
     * \param[in] length Lattice side length; total site count is
     * \p length squared.
     */
    explicit BufferLattice(int length);
    ~BufferLattice() = default;
    // Deep-copying two full per-site matrix fields for no current
    // caller isn't worth supporting -- disable copying (nothing needs
    // it).
    BufferLattice(const BufferLattice &) = delete;
    BufferLattice &operator=(const BufferLattice &) = delete;

    /// First scratch matrix field.
    std::vector<Matrix> buffer1;
    /// Second scratch matrix field.
    std::vector<Matrix> buffer2;
};

#endif  // SRC_LATTICE_H_

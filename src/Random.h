#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

#include <cstddef>
#include <random>
#include <vector>

#include "gsl/gsl_rng.h"

/**
 * Random-number source for the whole simulation: a 64-bit Mersenne
 * Twister (MT19937-64, Nishimura & Matsumoto) for the reproducible
 * uniform/Gaussian streams the physics depends on bit-for-bit, plus a
 * secondary GSL generator for Poisson sampling, plus a tabulated
 * inverse-CDF sampler for a gamma-distribution-based radial profile
 * used in constituent-quark position sampling (see
 * Init::sampleConstituentQuarkGeometry()).
 */
class Random {
  private:
    // MT19937-64 constants (Nishimura & Matsumoto).
    /// State array length.
    static constexpr int NN = 312;
    /// Twist recurrence offset.
    static constexpr int MM = 156;
    /// Twist recurrence's alternating XOR mask.
    static constexpr unsigned long long MATRIX_A = 0xB5026F5AA96619E9ULL;
    static constexpr unsigned long long UM =
        0xFFFFFFFF80000000ULL;  // most significant 33 bits
    static constexpr unsigned long long LM =
        0x7FFFFFFFULL;  // least significant 31 bits

    /// `0` if the next gauss() call must generate a fresh polar pair,
    /// `1` if \c gset_ holds an already-generated, unused partner.
    int iset_;
    /// The unused partner from the last generated Box-Muller polar
    /// pair, valid only when \c iset_ is `1`.
    double gset_;

    /// MT19937-64 state array.
    unsigned long long mt_[NN];
    /// Reusable scratch buffer for gaussBulk()'s untempered raw words,
    /// owned here so repeated calls don't reallocate.
    std::vector<unsigned long long> bulkRawScratch_;
    /// Index of the next unused word in \c mt_; `NN+1` means \c mt_ is
    /// not yet initialized (twist() will seed it with 5489 on first
    /// use).
    int mti_;

    /// Secondary GSL generator, used only by poisson().
    gsl_rng *gslRandom_;

    /// Tabulated regularized-upper-incomplete-gamma CDF values built by
    /// setGammaIncCDF(), indexed in step with \c gammaIncCDFx_.
    std::vector<double> gammaIncCDF_;
    /// Abscissas corresponding to \c gammaIncCDF_'s tabulated CDF
    /// values.
    std::vector<double> gammaIncCDFx_;

    /**
     * Draws \p count raw (untempered) MT19937-64 words, advancing the
     * generator's state exactly as \p count calls to genrand64_int64()
     * would, but without tempering them (tempering is applied
     * separately and can be parallelized, since it has no
     * inter-word dependency) -- used by gaussBulk() to draw many words
     * at once while still reproducing the scalar stream bit-for-bit.
     * \param[out] out Receives \p count raw words, in generation order.
     * \param[in] count Number of words to draw.
     */
    void genrand64RawBulk(unsigned long long *out, std::size_t count);
    /**
     * Refills \c mt_[0..NN-1] with the next \c NN raw words via the
     * MT19937-64 twist recurrence and resets \c mti_ to 0.
     */
    void twist();

  public:
    /**
     * Constructs a Random with its MT19937-64 state uninitialized (see
     * \c mti_) and a freshly allocated GSL `taus` generator; call
     * init_genrand64() and gslRandomInit() before drawing any numbers.
     */
    Random() {
        iset_ = 0;
        gset_ = 0.;
        mti_ = NN + 1;
        gslRandom_ = gsl_rng_alloc(gsl_rng_taus);

    };  // constructor

    /**
     * Frees the GSL generator this instance owns.
     */
    ~Random() { gsl_rng_free(gslRandom_); };  // destructor

    // gslRandom_ is a raw handle freed in the destructor; default copies
    // would double-free it, so disable copying (nothing needs it).
    Random(const Random &) = delete;
    Random &operator=(const Random &) = delete;

    /**
     * Seeds the MT19937-64 state.
     * \param[in] seed Seed value.
     */
    void init_genrand64(unsigned long long seed);
    /**
     * Draws one raw MT19937-64 word.
     * \return A pseudorandom integer, uniform on \f$[0, 2^{64}-1]\f$.
     */
    unsigned long long genrand64_int64(void);
    /**
     * Draws one raw word and discards its lowest bit.
     * \return A pseudorandom integer, uniform on \f$[0, 2^{63}-1]\f$.
     */
    long long genrand64_int63(void);
    /**
     * Draws a uniform double using the top 53 bits of one raw word.
     * \return A pseudorandom real, uniform on \f$[0, 1]\f$ (both
     * endpoints reachable).
     */
    double genrand64_real1(void);
    /**
     * Draws a uniform double using the top 53 bits of one raw word,
     * normalized to exclude 1.
     * \return A pseudorandom real, uniform on \f$[0, 1)\f$.
     */
    double genrand64_real2(void);
    /**
     * Draws a uniform double using the top 52 bits of one raw word,
     * shifted to exclude both endpoints.
     * \return A pseudorandom real, uniform on \f$(0, 1)\f$.
     */
    double genrand64_real3(void);

    /**
     * Seeds the secondary GSL generator poisson() uses.
     * \param[in] seed Seed value.
     */
    void gslRandomInit(unsigned long long seed);
    /**
     * Draws a Poisson-distributed integer via the secondary GSL
     * generator.
     * \param[in] mean Distribution mean \f$\lambda\f$.
     * \return A pseudorandom integer, Poisson-distributed with mean
     * \p mean.
     */
    int poisson(const double mean);
    /**
     * Draws a Gaussian-distributed double via the Box-Muller polar
     * method, caching the unused partner of each generated pair in \c
     * gset_ so every other call is free.
     * \param[in] mean Distribution mean \f$\mu\f$.
     * \param[in] width Distribution standard deviation \f$\sigma\f$.
     * \return A pseudorandom real, drawn from
     * \f$\mathcal{N}(\mu, \sigma^2)\f$.
     */
    double gauss(double mean = 0., double width = 1.);
    /**
     * Draws \p count Gaussian-distributed doubles (standard normal;
     * scale/shift the result yourself), reproducing exactly the
     * sequence repeated gauss() calls would produce -- including
     * consuming (or contributing) a cached Box-Muller partner via \c
     * iset_/\c gset_ -- but with the polar-method's uniform draws and
     * rejection sampling batched and parallelized over OpenMP.
     * \param[out] out Receives \p count standard-normal values.
     * \param[in] count Number of values to draw.
     * \param[in,out] scratch Reusable scratch buffer, resized as
     * needed; pass the same buffer across calls to avoid reallocating.
     */
    void gaussBulk(
        double *out, std::size_t count, std::vector<double> &scratch);

    /**
     * Builds the tabulated inverse-CDF sampler sampleGammaInc() draws
     * from: 1000 abscissas \f$x_i = i \cdot x_{\max}/1000\f$ for
     * \f$i=0,\ldots,999\f$ (so evenly spaced over \f$[0, x_{\max})\f$,
     * \f$x_{\max} = \max(5, 5/\omega)\f$, never reaching \f$x_{\max}\f$
     * itself), paired with the regularized upper incomplete gamma
     * function \f$Q(1/\omega, x)\f$'s cumulative sum evaluated at each
     * bin's *left* edge (i.e. `gammaIncCDF_[i]` excludes bin `i`'s own
     * weight), normalized by the sum over all 1000 bins. Consequently
     * the table's last entry is strictly below `1`, not `1` itself;
     * sampleGammaInc() maps the remaining upper-tail probability (any
     * draw above that last entry) onto the last abscissa rather than
     * sampling it.
     * \param[in] omega Shape parameter; must match what
     * sampleGammaInc()'s caller subsequently scales its result by
     * (see Init::sampleConstituentQuarkGeometry()).
     */
    void setGammaIncCDF(const double omega);
    /**
     * Draws one sample from the distribution tabulated by
     * setGammaIncCDF(), via binary search on its CDF table (with the
     * upper-tail draws described there mapped to the last abscissa).
     * \return A pseudorandom real in `[0, ` \c gammaIncCDFx_'s last
     * entry `]`.
     */
    double sampleGammaInc();
};

#endif  // SRC_RANDOM_H_

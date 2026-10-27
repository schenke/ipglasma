// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

#include <cstddef>
#include <random>
#include <vector>

#include "gsl/gsl_rng.h"

/**
 * Random-number source for the whole simulation: a 64-bit Mersenne
 * Twister (MT19937-64, Nishimura & Matsumoto) for the reproducible
 * uniform, Gaussian and Gamma streams the physics depends on
 * bit-for-bit (including the hot-spot radii of HotSpotNucleon::sample()
 * via sampleGammaInc()), plus a secondary GSL generator for Poisson
 * sampling.
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
    /// Mask selecting a word's most significant 33 bits.
    static constexpr unsigned long long UM = 0xFFFFFFFF80000000ULL;
    /// Mask selecting a word's least significant 31 bits.
    static constexpr unsigned long long LM = 0x7FFFFFFFULL;

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
    /// Reusable scratch buffer for gaussBulk()'s per-thread counts of
    /// accepted pairs and their prefix sum.
    std::vector<std::size_t> bulkCountScratch_;
    /// Reusable scratch buffer for gaussBulk()'s accept decision of each
    /// candidate pair of a rejection round, made once in the counting pass
    /// and reused by the writing pass.
    std::vector<unsigned char> bulkAcceptScratch_;
    /// Index of the next unused word in \c mt_; `NN+1` means \c mt_ is
    /// not yet initialized (twist() will seed it with 5489 on first
    /// use).
    int mti_;

    /// Secondary GSL generator, used only by poisson().
    gsl_rng *gslRandom_;

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
     * rejection sampling batched: in each rejection round the raw MT words
     * are drawn serially, and their tempering, the acceptance test and the
     * compaction of the accepted pairs (with a prefix sum over the threads,
     * in the original order) run in parallel over OpenMP.
     * \param[out] out Receives \p count standard-normal values.
     * \param[in] count Number of values to draw.
     * \param[in,out] scratch Reusable scratch buffer, resized as
     * needed; pass the same buffer across calls to avoid reallocating.
     */
    void gaussBulk(
        double *out, std::size_t count, std::vector<double> &scratch);

    /**
     * Draws a Gamma-distributed double with unit scale, by the method of
     * Marsaglia and Tsang (ACM Trans. Math. Softw. 26, 363 (2000)) on the
     * Mersenne Twister: gauss() and genrand64_real3() draws, repeated
     * until one is accepted; for \p shape < 1 a draw for \p shape + 1
     * times \f$U^{1/\mathrm{shape}}\f$.
     * \param[in] shape Shape parameter \f$k > 0\f$.
     * \return A pseudorandom real with the density
     * \f$x^{k-1}e^{-x}/\Gamma(k)\f$ on \f$(0, \infty)\f$, so mean and
     * variance \f$k\f$.
     */
    double gammaVariate(double shape);
    /**
     * Draws \f$x\f$ with the density \f$\omega\,Q(1/\omega, x)\f$ on
     * \f$(0, \infty)\f$, \f$Q\f$ the regularized upper incomplete gamma
     * function, so \f$\langle x\rangle = (1+\omega)/(2\omega)\f$; the
     * hot-spot radii of HotSpotNucleon::sample() for \f$\omega \ne 1\f$.
     * Exact: with \f$s = 1/\omega\f$, \f$x = U G\f$ for \f$U\f$ uniform
     * on (0, 1) and \f$G\f$ Gamma-distributed with shape \f$s+1\f$
     * (gammaVariate()) has the density
     * \f$\int_x^\infty t^{s-1}e^{-t}\,dt/\Gamma(s+1) = Q(s, x)/s\f$.
     * \param[in] omega \f$\omega > 0\f$.
     * \return A pseudorandom real in \f$(0, \infty)\f$.
     */
    double sampleGammaInc(double omega);
};

#endif  // SRC_RANDOM_H_

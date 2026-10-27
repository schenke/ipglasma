// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#include "Random.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "Instrumentation.h"
#include "gsl/gsl_randist.h"

// This file contains the random generator, which is

/*
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/

/* initializes mt_[NN] with a seed */
void Random::init_genrand64(unsigned long long seed) {
    mt_[0] = seed;
    for (mti_ = 1; mti_ < NN; mti_++)
        mt_[mti_] =
            (6364136223846793005ULL * (mt_[mti_ - 1] ^ (mt_[mti_ - 1] >> 62))
             + mti_);
}

/* refills mt_[0..NN-1] with the next NN raw words and resets mti_ to 0 */
void Random::twist() {
    static const unsigned long long mag01[2] = {0ULL, MATRIX_A};

    if (mti_ == NN + 1) init_genrand64(5489ULL);

    int i = 0;
    unsigned long long x;
    for (; i < NN - MM; ++i) {
        x = (mt_[i] & UM) | (mt_[i + 1] & LM);
        mt_[i] = mt_[i + MM] ^ (x >> 1) ^ mag01[static_cast<int>(x & 1ULL)];
    }
    for (; i < NN - 1; ++i) {
        x = (mt_[i] & UM) | (mt_[i + 1] & LM);
        mt_[i] =
            mt_[i + (MM - NN)] ^ (x >> 1) ^ mag01[static_cast<int>(x & 1ULL)];
    }
    x = (mt_[NN - 1] & UM) | (mt_[0] & LM);
    mt_[NN - 1] = mt_[MM - 1] ^ (x >> 1) ^ mag01[static_cast<int>(x & 1ULL)];
    mti_ = 0;
}

void Random::genrand64RawBulk(unsigned long long *out, std::size_t count) {
    std::size_t produced = 0;

    while (produced < count) {
        if (mti_ >= NN) twist();

        const std::size_t available = static_cast<std::size_t>(NN - mti_);
        const std::size_t remaining = count - produced;
        const std::size_t take =
            (remaining < available) ? remaining : available;
        std::memcpy(
            out + produced, mt_ + mti_, take * sizeof(unsigned long long));
        mti_ += static_cast<int>(take);
        produced += take;
    }
}

namespace {
/**
 * MT19937-64 tempering step, applied to one raw generator word to
 * improve its equidistribution properties. Stateless and independent
 * per word, so it can run in parallel across many words (see
 * Random::gaussBulk()).
 * \param[in] x Raw (untempered) generator word.
 * \return The tempered word.
 */
unsigned long long temperMT64(unsigned long long x) {
    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);
    return x;
}
}  // namespace

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long Random::genrand64_int64(void) {
    if (mti_ >= NN) /* generate NN words at one time */
        twist();

    return temperMT64(mt_[mti_++]);
}

/* generates a random number on [0, 2^63-1]-interval */
long long Random::genrand64_int63(void) {
    return (long long)(genrand64_int64() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double Random::genrand64_real1(void) {
    return (genrand64_int64() >> 11) * (1.0 / 9007199254740991.0);
}

/* generates a random number on [0,1)-real-interval */
double Random::genrand64_real2(void) {
    return (genrand64_int64() >> 11) * (1.0 / 9007199254740992.0);
}

/* generates a random number on (0,1)-real-interval */
double Random::genrand64_real3(void) {
    return ((genrand64_int64() >> 12) + 0.5) * (1.0 / 4503599627370496.0);
}

double Random::gauss(double mean, double width) {
    // produces random numbers with distribution
    // p(y)dy  = 1/\sqrt{2\pi \sigma^2} exp(-(y-\mu)^2/(2\sigma^2)) dy

    // default is \mu=0, \sigma = 1

    // if iset_=0, generate both random numbers new, else use the unused one
    // from the previous step. this saves lots of time.
    double fac, rsq, v1, v2;

    if (iset_ == 0) {
        do {
            v1 = 2.0 * genrand64_real3() - 1.0;
            v2 = 2.0 * genrand64_real3() - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq > 1. || rsq == 0.);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset_ = v1 * fac;
        iset_ = 1;
        return mean + width * v2 * fac;
    } else {
        iset_ = 0;
        return mean + width * gset_;
    }
}

void Random::gaussBulk(
    double *out, std::size_t count, std::vector<double> &scratch) {
    if (count == 0) return;

    // Preserve the exact scalar gauss() stream, including the cached partner
    // in gset_/iset_. Generate candidate polar pairs in ordered rejection
    // rounds: if R accepted pairs are still required, draw exactly R
    // candidates.  If some are rejected, the next round draws exactly the
    // number still needed. Consequently no MT word is consumed beyond the
    // scalar stopping point.
    std::size_t outOffset = 0;
    if (iset_ != 0) {
        out[outOffset++] = gset_;
        iset_ = 0;
        if (outOffset == count) return;
    }

    const std::size_t remaining = count - outOffset;
    const std::size_t pairCount = (remaining + 1) / 2;

    // First 3*pairCount doubles hold accepted (v1,v2,rsq) triples.  The final
    // 3*pairCount doubles hold the (v1,v2,rsq) candidates of one rejection
    // round, with v = 2u - 1.  Raw MT words and the accept decisions live in
    // reusable Random-owned buffers.
    scratch.resize(6 * pairCount);
    double *accepted = scratch.data();
    double *candidates = scratch.data() + 3 * pairCount;
    bulkRawScratch_.resize(2 * pairCount);
    bulkAcceptScratch_.resize(pairCount);

    const bool profile = ipg::Profiler::instance().enabled();
    double phaseStart = profile ? ipg::wallSeconds() : 0.0;

    std::size_t acceptedPairs = 0;
    while (acceptedPairs < pairCount) {
        const std::size_t needed = pairCount - acceptedPairs;
        const std::size_t uniformCount = 2 * needed;

        // Advance the MT state exactly as the scalar generator would, but copy
        // the untempered state words out in bulk (the recurrence is serial).
        genrand64RawBulk(bulkRawScratch_.data(), uniformCount);

        // Temper the words, test the candidate pairs and keep the accepted
        // ones in one parallel region (issue #121). Each thread handles a
        // contiguous range of pairs: it computes and tests its candidates and
        // counts the accepted ones, a prefix sum over the threads gives each
        // one its offset, and the accepted pairs are copied in the original
        // RNG order, so the rejection decisions and the order are exactly
        // those of gauss(). The copy uses the stored decisions and values
        // rather than recomputing them, so the offsets of the prefix sum
        // cannot disagree with it, whatever the compiler does with the
        // floating-point arithmetic of the two loops.
#pragma omp parallel
        {
#ifdef _OPENMP
            const std::size_t thread = omp_get_thread_num();
            const std::size_t threads = omp_get_num_threads();
#else
            const std::size_t thread = 0;
            const std::size_t threads = 1;
#endif
#pragma omp single
            bulkCountScratch_.assign(threads + 1, 0);
            const std::size_t begin = needed * thread / threads;
            const std::size_t end = needed * (thread + 1) / threads;
            std::size_t acceptedHere = 0;
            for (std::size_t pair = begin; pair < end; ++pair) {
                double v[2];
                for (int k = 0; k < 2; ++k) {
                    const unsigned long long x =
                        temperMT64(bulkRawScratch_[2 * pair + k]);
                    v[k] =
                        2.0 * (((x >> 12) + 0.5) * (1.0 / 4503599627370496.0))
                        - 1.0;
                }
                const double rsq = v[0] * v[0] + v[1] * v[1];
                const bool accept = rsq <= 1.0 && rsq != 0.0;
                double *candidate = candidates + 3 * pair;
                candidate[0] = v[0];
                candidate[1] = v[1];
                candidate[2] = rsq;
                bulkAcceptScratch_[pair] = accept;
                acceptedHere += accept;
            }
            bulkCountScratch_[thread + 1] = acceptedHere;
#pragma omp barrier
#pragma omp single
            for (std::size_t t = 0; t < threads; ++t) {
                bulkCountScratch_[t + 1] += bulkCountScratch_[t];
            }
            std::size_t next = acceptedPairs + bulkCountScratch_[thread];
            for (std::size_t pair = begin; pair < end; ++pair) {
                if (bulkAcceptScratch_[pair]) {
                    std::memcpy(
                        accepted + 3 * next++, candidates + 3 * pair,
                        3 * sizeof(double));
                }
            }
        }
        acceptedPairs += bulkCountScratch_.back();
    }

    if (profile) {
        ipg::Profiler::instance().add(
            "random.gauss_bulk.generate", ipg::wallSeconds() - phaseStart);
        phaseStart = ipg::wallSeconds();
    }

    const std::size_t fullPairs = remaining / 2;
#pragma omp parallel for
    for (std::size_t pair = 0; pair < fullPairs; ++pair) {
        const std::size_t base = 3 * pair;
        const double v1 = accepted[base];
        const double v2 = accepted[base + 1];
        const double rsq = accepted[base + 2];
        const double fac = sqrt(-2.0 * log(rsq) / rsq);
        const std::size_t outBase = outOffset + 2 * pair;
        out[outBase] = v2 * fac;
        out[outBase + 1] = v1 * fac;
    }

    if ((remaining & 1U) != 0U) {
        const std::size_t pair = pairCount - 1;
        const std::size_t base = 3 * pair;
        const double v1 = accepted[base];
        const double v2 = accepted[base + 1];
        const double rsq = accepted[base + 2];
        const double fac = sqrt(-2.0 * log(rsq) / rsq);
        out[count - 1] = v2 * fac;
        gset_ = v1 * fac;
        iset_ = 1;
    }

    if (profile) {
        ipg::Profiler::instance().add(
            "random.gauss_bulk.transform", ipg::wallSeconds() - phaseStart);
    }
}

double Random::gammaVariate(const double shape) {
    if (shape < 1.) {
        // Gamma(k) = Gamma(k + 1) U^(1/k)
        const double u = genrand64_real3();
        return gammaVariate(shape + 1.) * std::pow(u, 1. / shape);
    }
    // Marsaglia & Tsang: d v with v = (1 + c z)^3 for a Gaussian z,
    // accepted with the squeeze first and the exact test second
    const double d = shape - 1. / 3.;
    const double c = 1. / std::sqrt(9. * d);
    while (true) {
        double z, v;
        do {
            z = gauss();
            v = 1. + c * z;
        } while (v <= 0.);
        v = v * v * v;
        const double u = genrand64_real3();
        const double z2 = z * z;
        if (u < 1. - 0.0331 * z2 * z2) return d * v;
        if (std::log(u) < 0.5 * z2 + d * (1. - v + std::log(v))) return d * v;
    }
}

double Random::sampleGammaInc(const double omega) {
    // x = U G, G ~ Gamma(1/omega + 1), has the density omega Q(1/omega, x);
    // two statements, so the order of the draws is fixed
    const double g = gammaVariate(1. / omega + 1.);
    return genrand64_real3() * g;
}

void Random::gslRandomInit(unsigned long long seed) {
    gsl_rng_set(gslRandom_, seed);
}

int Random::poisson(const double mean) {
    return (gsl_ran_poisson(gslRandom_, mean));
}


#include "Random.h"

#include <gsl/gsl_sf_gamma.h>

#include <cmath>
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

/* initializes mt[NN] with a seed */
void Random::init_genrand64(unsigned long long seed) {
    mt[0] = seed;
    for (mti = 1; mti < NN; mti++)
        mt[mti] =
            (6364136223846793005ULL * (mt[mti - 1] ^ (mt[mti - 1] >> 62))
             + mti);
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void Random::init_by_array64(
    unsigned long long init_key[], unsigned long long key_length) {
    unsigned long long i, j, k;
    init_genrand64(19650218ULL);
    i = 1;
    j = 0;
    k = (NN > key_length ? NN : key_length);
    for (; k; k--) {
        mt[i] =
            (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 3935559000370003845ULL))
            + init_key[j] + j; /* non linear */
        i++;
        j++;
        if (i >= NN) {
            mt[0] = mt[NN - 1];
            i = 1;
        }
        if (j >= key_length) j = 0;
    }
    for (k = NN - 1; k; k--) {
        mt[i] =
            (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 2862933555777941757ULL))
            - i; /* non linear */
        i++;
        if (i >= NN) {
            mt[0] = mt[NN - 1];
            i = 1;
        }
    }

    mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */
}

void Random::genrand64RawBulk(unsigned long long *out, std::size_t count) {
    static const unsigned long long mag01[2] = {0ULL, MATRIX_A};
    std::size_t produced = 0;

    while (produced < count) {
        if (mti >= NN) {
            if (mti == NN + 1) init_genrand64(5489ULL);

            int i = 0;
            unsigned long long x;
            for (; i < NN - MM; ++i) {
                x = (mt[i] & UM) | (mt[i + 1] & LM);
                mt[i] =
                    mt[i + MM] ^ (x >> 1) ^ mag01[static_cast<int>(x & 1ULL)];
            }
            for (; i < NN - 1; ++i) {
                x = (mt[i] & UM) | (mt[i + 1] & LM);
                mt[i] = mt[i + (MM - NN)] ^ (x >> 1)
                        ^ mag01[static_cast<int>(x & 1ULL)];
            }
            x = (mt[NN - 1] & UM) | (mt[0] & LM);
            mt[NN - 1] =
                mt[MM - 1] ^ (x >> 1) ^ mag01[static_cast<int>(x & 1ULL)];
            mti = 0;
        }

        const std::size_t available = static_cast<std::size_t>(NN - mti);
        const std::size_t remaining = count - produced;
        const std::size_t take =
            (remaining < available) ? remaining : available;
        std::memcpy(
            out + produced, mt + mti, take * sizeof(unsigned long long));
        mti += static_cast<int>(take);
        produced += take;
    }
}

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long Random::genrand64_int64(void) {
    int i;
    unsigned long long x;
    static unsigned long long mag01[2] = {0ULL, MATRIX_A};

    if (mti >= NN) { /* generate NN words at one time */

        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        if (mti == NN + 1) init_genrand64(5489ULL);

        for (i = 0; i < NN - MM; i++) {
            x = (mt[i] & UM) | (mt[i + 1] & LM);
            mt[i] = mt[i + MM] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
        }
        for (; i < NN - 1; i++) {
            x = (mt[i] & UM) | (mt[i + 1] & LM);
            mt[i] = mt[i + (MM - NN)] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
        }
        x = (mt[NN - 1] & UM) | (mt[0] & LM);
        mt[NN - 1] = mt[MM - 1] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];

        mti = 0;
    }

    x = mt[mti++];

    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);

    return x;
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

double Random::Gauss2(double mean, double sigma) {
    double x, y, z, result;
    do {
        y = genrand64_real3();
    } while (!y);
    z = genrand64_real3();
    x = z * 6.283185;
    result = mean + sigma * sin(x) * sqrt(-2. * log(y));
    return result;
}

double Random::Gauss(double mean, double width) {
    // produces random numbers with distribution
    // p(y)dy  = 1/\sqrt{2\pi \sigma^2} exp(-(y-\mu)^2/(2\sigma^2)) dy

    // default is \mu=0, \sigma = 1

    // if iset=0, generate both random numbers new, else use the unused one from
    // the previous step. this saves lots of time.
    double fac, rsq, v1, v2;

    // if(idum<0) iset=0;

    if (iset == 0) {
        do {
            v1 = 2.0 * genrand64_real3() - 1.0;
            v2 = 2.0 * genrand64_real3() - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq > 1. || rsq == 0.);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return mean + width * v2 * fac;
    } else {
        iset = 0;
        return mean + width * gset;
    }
}

void Random::GaussBulk(
    double *out, std::size_t count, std::vector<double> &scratch) {
    if (count == 0) return;

    // Preserve the exact scalar Gauss() stream, including the cached partner
    // in gset/iset. Generate candidate polar pairs in ordered rejection rounds:
    // if R accepted pairs are still required, draw exactly R candidates.  If
    // some are rejected, the next round draws exactly the number still needed.
    // Consequently no MT word is consumed beyond the scalar stopping point.
    std::size_t outOffset = 0;
    if (iset != 0) {
        out[outOffset++] = gset;
        iset = 0;
        if (outOffset == count) return;
    }

    const std::size_t remaining = count - outOffset;
    const std::size_t pairCount = (remaining + 1) / 2;

    // First 3*pairCount doubles hold accepted (v1,v2,rsq) triples.  The final
    // 2*pairCount doubles are temporary converted uniforms for one rejection
    // round.  Raw MT words live in a reusable Random-owned buffer.
    scratch.resize(5 * pairCount);
    double *accepted = scratch.data();
    double *uniforms = scratch.data() + 3 * pairCount;
    bulkRawScratch_.resize(2 * pairCount);

    const bool profile = ipg::Profiler::instance().enabled();
    double phaseStart = profile ? ipg::wallSeconds() : 0.0;

    std::size_t acceptedPairs = 0;
    while (acceptedPairs < pairCount) {
        const std::size_t needed = pairCount - acceptedPairs;
        const std::size_t uniformCount = 2 * needed;

        // Advance the MT state exactly as the scalar generator would, but copy
        // the untempered state words out in bulk. Tempering is independent for
        // every word and can therefore use the existing OpenMP team safely.
        genrand64RawBulk(bulkRawScratch_.data(), uniformCount);

#pragma omp parallel for
        for (std::size_t q = 0; q < uniformCount; ++q) {
            unsigned long long x = bulkRawScratch_[q];
            x ^= (x >> 29) & 0x5555555555555555ULL;
            x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
            x ^= (x << 37) & 0xFFF7EEE000000000ULL;
            x ^= (x >> 43);
            uniforms[q] = ((x >> 12) + 0.5) * (1.0 / 4503599627370496.0);
        }

        // Scan candidates serially in original RNG order.  This preserves the
        // precise rejection decisions and accepted-pair ordering of Gauss().
        for (std::size_t pair = 0; pair < needed; ++pair) {
            const double v1 = 2.0 * uniforms[2 * pair] - 1.0;
            const double v2 = 2.0 * uniforms[2 * pair + 1] - 1.0;
            const double rsq = v1 * v1 + v2 * v2;
            if (rsq <= 1.0 && rsq != 0.0) {
                const std::size_t base = 3 * acceptedPairs++;
                accepted[base] = v1;
                accepted[base + 1] = v2;
                accepted[base + 2] = rsq;
            }
        }
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
        gset = v1 * fac;
        iset = 1;
    }

    if (profile) {
        ipg::Profiler::instance().add(
            "random.gauss_bulk.transform", ipg::wallSeconds() - phaseStart);
    }
}

void Random::setGammaIncCDF(const double omega) {
    gammaIncCDF_.clear();
    gammaIncCDFx_.clear();
    double xmax = std::max(5., 5. / omega);
    int nX = 1000;
    gammaIncCDF_.resize(nX, 0.);
    gammaIncCDFx_.resize(nX, 0.);
    double CDF = 0.;
    for (int i = 0; i < nX; i++) {
        double x = i * xmax / nX;
        gammaIncCDFx_[i] = x;
        gammaIncCDF_[i] = CDF;
        CDF += gsl_sf_gamma_inc_Q(1. / omega, x);
    }
    for (int i = 0; i < nX; i++) {
        gammaIncCDF_[i] /= CDF;
    }
}

double Random::sampleGammaInc() {
    double u = genrand64_real1();
    int idx_l = 0;
    int idx_h = gammaIncCDF_.size() - 1;
    if (u > gammaIncCDF_[idx_h]) {
        return gammaIncCDFx_[idx_h];
    }
    int idx_m = static_cast<int>((idx_l + idx_h) / 2);
    while (idx_h - idx_l > 1) {
        if (u < gammaIncCDF_[idx_m]) {
            idx_h = idx_m;
        } else {
            idx_l = idx_m;
        }
        idx_m = static_cast<int>((idx_l + idx_h) / 2);
    }
    return gammaIncCDFx_[idx_m];
}

void Random::gslRandomInit(unsigned long long seed) {
    gsl_rng_set(gslRandom, seed);
}

double Random::NBD(double nbar, double k) {
    double p = k / (nbar + k);
    double n = k;

    return gsl_ran_negative_binomial(gslRandom, p, n);
}

int Random::Poisson(const double mean) {
    return (gsl_ran_poisson(gslRandom, mean));
}

double Random::tdist(double nu) {
    // produces random numbers with distribution
    // p(x) dx = {\Gamma((\nu + 1)/2) \over \sqrt{\pi \nu} \Gamma(\nu/2)}
    //           (1 + x^2/\nu)^{-(\nu + 1)/2} dx
    // with mean 0 and variance nu/(nu-2)
    // however, I take care of that variance when returning the value, so that
    // variance is always 1

    if (nu <= 2) {
        std::cerr << "nu has to be > 2. Exiting." << std::endl;
        exit(1);
    }

    double f;

    f = sqrt((nu - 2.) / nu) * gsl_ran_tdist(gslRandom, nu);

    return f;
}

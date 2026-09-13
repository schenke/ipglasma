#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

#include <cstddef>
#include <random>
#include <vector>

#include "gsl/gsl_rng.h"

class Random {
  private:
    // MT19937-64 constants (Nishimura & Matsumoto).
    static constexpr int NN = 312;
    static constexpr int MM = 156;
    static constexpr unsigned long long MATRIX_A = 0xB5026F5AA96619E9ULL;
    static constexpr unsigned long long UM =
        0xFFFFFFFF80000000ULL;  // most significant 33 bits
    static constexpr unsigned long long LM =
        0x7FFFFFFFULL;  // least significant 31 bits

    int iset_;
    double gset_;

    unsigned long long mt_[NN];
    std::vector<unsigned long long> bulkRawScratch_;
    /* mti_==NN+1 means mt_[NN] is not initialized */
    int mti_;

    gsl_rng *gslRandom_;

    std::vector<double> gammaIncCDF_;
    std::vector<double> gammaIncCDFx_;

    void genrand64RawBulk(unsigned long long *out, std::size_t count);
    /* refills mt_[0..NN-1] with the next NN raw words and resets mti_ */
    void twist();

  public:
    Random() {
        iset_ = 0;
        gset_ = 0.;
        mti_ = NN + 1;
        gslRandom_ = gsl_rng_alloc(gsl_rng_taus);

    };  // constructor

    ~Random() { gsl_rng_free(gslRandom_); };  // destructor

    // gslRandom_ is a raw handle freed in the destructor; default copies
    // would double-free it, so disable copying (nothing needs it).
    Random(const Random &) = delete;
    Random &operator=(const Random &) = delete;

    void init_genrand64(unsigned long long seed);
    unsigned long long genrand64_int64(void);
    long long genrand64_int63(void);
    double genrand64_real1(void);
    double genrand64_real2(void);
    double genrand64_real3(void);

    void gslRandomInit(unsigned long long seed);
    int poisson(const double mean);
    double gauss(double mean = 0., double width = 1.);
    void gaussBulk(
        double *out, std::size_t count, std::vector<double> &scratch);

    void setGammaIncCDF(const double omega);
    double sampleGammaInc();
};

#endif  // SRC_RANDOM_H_

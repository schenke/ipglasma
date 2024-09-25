#ifndef Random_h
#define Random_h

#include <complex>
#include <cstdlib>
#include <iostream>

#include "Matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"

#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL         /* Least significant 31 bits */

using namespace std;

class Random {
  private:
    int iset;
    double gset;

    unsigned long long mt[NN];
    /* mti==NN+1 means mt[NN] is not initialized */
    int mti;

    gsl_rng *gslRandom;

  public:
    Random() {
        iset = 0;
        gslRandom = gsl_rng_alloc(gsl_rng_taus);
    };  // constructor

    ~Random() { gsl_rng_free(gslRandom); };  // destructor
    void init_genrand64(unsigned long long seed);
    void init_by_array64(
        unsigned long long init_key[], unsigned long long key_length);
    unsigned long long genrand64_int64(void);
    long long genrand64_int63(void);
    double genrand64_real1(void);
    double genrand64_real2(void);
    double genrand64_real3(void);

    void gslRandomInit(unsigned long long seed);
    double tdist(double nu);
    double NBD(double nbar, double k);
    int Poisson(const double mean);
    double Gauss(double mean = 0., double width = 1.);
    double Gauss2(double mean, double sigma);
};

#endif

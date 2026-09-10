// Init.h is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#ifndef SRC_INIT_H_
#define SRC_INIT_H_

#include <cstdint>

#include "FFT.h"
#include "Glauber.h"
#include "Group.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "PrettyOstream.h"
#include "Random.h"

enum class InitializationMethod {
    SampleColorCharges,
    ReadWlineText,
    ReadWlineBinary,
    InitializeAfterJimwlk
};

class Init {
  private:
    int const static iymaxNuc_ = 44;  // for the Tp-y table

    int const static iTpmax_ =
        200;  // updated in March 2019 to a larger T_A range

    double const deltaYNuc_ = 0.25;  // for the new table
    FFT fft_;
    //  Matrix** A;
    //  Glauber *glauber;
    double Qs2Nuclear_[iTpmax_][iymaxNuc_];
    double Tlist_[iTpmax_];

    std::vector<vector<float>> nucleonPosArrA_;
    std::vector<vector<float>> nucleonPosArrB_;

    // list of x and y coordinates of nucleons in nucleus A
    std::vector<ReturnValue> nucleusA_;
    // list of x and y coordinates of nucleons in nucleus B
    std::vector<ReturnValue> nucleusB_;

    PrettyOstream messager_;

    static constexpr int Nc_ = 3;
    static constexpr int Nc2m1_ = Nc_ * Nc_ - 1;
    Group *group_ptr_;
    Random *random_ptr_;

    Matrix one_;
    vector<vector<double>> xq1_, xq2_, yq1_, yq2_, BGq1_, BGq2_, gauss1_,
        gauss2_;

  public:
    // Constructor.
    explicit Init(const int nn[]) : fft_(nn), one_(1.) {};

    ~Init() {};

    void init(
        Lattice *lat, Group *group, Parameters *param, Random *random,
        Glauber *glauber, InitializationMethod init_method);
    void shiftFieldsWithImpactParameter(Lattice *lat, Parameters *param);
    void initializeForwardLightCone(Lattice *lat, Parameters *param);
    void sampleImpactParameter(Parameters *param);
    void sampleTA(Parameters *param, Random *random, Glauber *glauber);
    void readNuclearQs(Parameters *param);
    void solveAxb(double *Jab, double *Fa, std::vector<double> &xvec);

    double getNuclearQs2(double T, double y);
    void setColorChargeDensity(
        Lattice *lat, Parameters *param, Random *random, Glauber *glauber);
    void computeCollisionGeometryQuantities(Lattice *lat, Parameters *param);
    void setV(Lattice *lat, Parameters *param, Random *random);
    void readVFromFile(Lattice *lat, Parameters *param, int format);

    // void eccentricity(Lattice *lat, Group *group, Parameters *param, Random
    // *random, Glauber *glauber);

    Matrix getUfromExponent(std::vector<double> &Q);
    bool findUInForwardLightconeChun(
        Matrix &U1, Matrix &U2, Matrix &Usol, std::uint64_t retrySeed);

    void readInNucleusConfigs(
        const int nucleusA, const int lightNucleusOption,
        const int polarizationFlag, const double polJz,
        vector<vector<float>> &nucleonPosArr);
    void generateNucleusConfiguration(
        Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
        double beta3, double beta4, double gamma, bool forceDminFlag,
        double d_min, double dR_np, double da_np,
        std::vector<ReturnValue> &nucleus);
    void generateNucleusConfigurationWithWoodsSaxon(
        Random *random, int A, int Z, double a_WS, double R_WS, double d_min,
        double dR_np, double da_np, std::vector<ReturnValue> &nucleus);
    void generateNucleusConfigurationWithDeformedWoodsSaxon(
        Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
        double beta3, double beta4, double d_min, double dR_np, double da_np,
        std::vector<ReturnValue> &nucleus);
    void generateNucleusConfigurationWithDeformedWoodsSaxon2(
        Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
        double beta3, double beta4, double gamma, double dR_np, double da_np,
        std::vector<ReturnValue> &nucleus);
    void generateNucleusConfigurationWithDeformedWoodsSaxonForceDmin(
        Random *random, int A, int Z, double a_WS, double R_WS, double beta2,
        double beta3, double beta4, double gamma, double d_min, double dR_np,
        double da_np, std::vector<ReturnValue> &nucleus);
    double sampleRFromWoodsSaxon(
        Random *random, double a_WS, double R_WS) const;
    void sampleRAndCosthetaFromDeformedWoodsSaxon(
        Random *random, double a_WS, double R_WS, double beta2, double beta3,
        double beta4, double &r, double &costheta) const;
    double fermiDistribution(double r, double R_WS, double a_WS) const;
    double sphericalHarmonics(int l, double ct) const;
    double sphericalHarmonicsY22(double ct, double phi) const;
    void recenterNucleus(
        std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);
    void recenterNucleus(std::vector<ReturnValue> &nucleus);
    void assignProtons(
        Random *random, std::vector<ReturnValue> &nucleus, const int Z);
    void rotateNucleus(
        double phi_global, double theta_global,
        std::vector<ReturnValue> &nucleus);
    void rotateNucleus3D(Random *random, std::vector<ReturnValue> &nucleus);

    void samplePartonPositions(
        Parameters *param, Random *random, std::vector<double> &x_array,
        std::vector<double> &y_array, std::vector<double> &z_array,
        std::vector<double> &BGq_array);

    double sampleLogNormalDistribution(
        Random *random, const double mean, const double variance);

    void sampleQsNormalization(
        Random *random, Parameters *param, const int Nq,
        std::vector<double> &gauss_array);
    int sampleNumberOfPartons(Random *random, Parameters *param);
};

#endif  // SRC_INIT_H_

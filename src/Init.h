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
    ReadWlineBinary
};

class Init {
  private:
    int const static iymaxNuc_ = 44;  // for the Tp-y table

    int const static iTpmax_ =
        240;  // updated in Sep 2026 to a 10x extended T_A range

    double const deltaYNuc_ = 0.25;  // for the new table
    FFT fft_;
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
    // sampleTA's four mutually-exclusive nucleonPositionsFromFile branches,
    // each filling nucleusA_/nucleusB_ by one particular method.
    void sampleTAWoodsSaxon(
        Parameters *param, Random *random, Glauber *glauber);
    void sampleTAFromConfigFiles(Random *random, Glauber *glauber);
    void sampleTAFromAlvioliFiles(Random *random, Glauber *glauber);
    // Reads one nucleus's worth of nucleon positions from a randomly
    // selected Alvioli correlated-Pb-208 configuration file, appending them
    // to nucleus. label ("A" or "B") is used only for log messages.
    void readOneAlvioliNucleus(
        Random *random, int nucleonCount, const std::string &label,
        std::vector<ReturnValue> &nucleus);
    // Applies sampleTA's global nucleus rotation for one polarization flag
    // (0: random 3D, 1: longitudinal, 2: transverse); called once each for
    // the projectile and target.
    void applyPolarizationRotation(
        Random *random, int polarizationFlag,
        std::vector<ReturnValue> &nucleus);
    void readNuclearQs(Parameters *param);
    void solveAxb(double *Jab, double *Fa, std::vector<double> &xvec);

    double getNuclearQs2(double T, double y);
    // Sets g2mu2A/g2mu2B at one cell from its already-accumulated TpA/TpB,
    // via getNuclearQs2 and (if enabled) the fluctuating-x iterative solve.
    // Called from setColorChargeDensity's per-cell loop.
    void computeCellColorCharge(
        Lattice *lat, Parameters *param, int ipos, double a,
        double rapidityA, double rapidityB);
    void setColorChargeDensity(
        Lattice *lat, Parameters *param, Random *random, Glauber *glauber);
    // Converts param's input rapidity to true rapidity when it's flagged
    // as pseudorapidity, otherwise passes it through unchanged.
    void computeEffectiveRapidities(
        Parameters *param, double &rapidityA, double &rapidityB);
    // setColorChargeDensity's useNucleus==0 (constant g^2mu background)
    // branch; returns having already called param->setSuccess(1).
    void setConstantColorChargeDensity(Lattice *lat, Parameters *param);
    // Samples each nucleon's proton-anisotropy angle phi (or sets it to 0
    // when protonAnisotropy is off).
    void sampleNucleonAnisotropyAngles(Parameters *param, Random *random);
    // Samples constituent-quark positions/widths (xq1_/yq1_/BGq1_ etc.) and
    // each nucleon's Qs-normalization gauss factor, for both nuclei.
    void sampleConstituentQuarkGeometry(Parameters *param, Random *random);
    // setColorChargeDensity's useSmoothNucleus==1 branch: sets TpA/TpB from
    // the smooth (undeformed) Woods-Saxon thickness functions.
    void computeSmoothNucleusThickness(
        Lattice *lat, Parameters *param, Glauber *glauber);
    // setColorChargeDensity's default branch: sets TpA/TpB by summing each
    // sampled nucleon's (constituent-quark or single-Gaussian) thickness.
    void computeThicknessFromNucleons(
        Lattice *lat, Parameters *param, double nucleiInAverage);
    // Determines Npart/Ncoll from the (already-sampled) nucleon positions,
    // writes NcollList*.dat/NpartList*.dat, and sets param->setNpart.
    // Returns false if useFixedNpart is set and this event's Npart doesn't
    // match, signaling the caller to abort and resample.
    bool determineNpartAndNcoll(Parameters *param, int &Npart, int &Ncoll);
    // Sets param's running-coupling alpha_s from whichever Qs choice
    // param->getRunWithQs() selects, or a fixed value if running coupling
    // is disabled or alpha_s runs with k_T instead.
    void computeAndSetRunningAlphaS(Parameters *param);
    void computeCollisionGeometryQuantities(Lattice *lat, Parameters *param);
    // Scans the full lattice, accumulating the Qs/T_pp collision-geometry
    // averages computeCollisionGeometryQuantities reports and stores.
    void scanCollisionGeometry(
        Lattice *lat, Parameters *param, int N, double a, double b,
        double phiRP, double &averageQs, double &averageQs2,
        double &averageQs2Avg, double &averageQs2min, double &averageQs2min2,
        double &Tpp, int &count);
    // Logs computeCollisionGeometryQuantities' N_part/N_coll/T_pp/Qs/alpha_s
    // summary.
    void logCollisionGeometryQuantities(
        Parameters *param, int Npart, int Ncoll, double Tpp, double a,
        double averageQs2, double averageQs2Avg, double averageQs2min,
        double averageQs2min2, int count);
    // Appends this event's usedParameters<id>.dat entry (called only when
    // computeCollisionGeometryQuantities marks the event a success).
    void writeUsedParametersFile(
        Parameters *param, double phiRP, int Npart, int Ncoll);
    // Writes this event's NgluonEstimators<id>.dat file.
    void writeNgluonEstimatorsFile(
        Parameters *param, double a, double averageQs2, double averageQs2Avg,
        double averageQs2min2, int count);
    void setV(Lattice *lat, Parameters *param, Random *random);
    void readVFromFile(Lattice *lat, Parameters *param, int format);
    // readVFromFile's format==1/format==2 branches, each called once for
    // the projectile's Wilson line file and once for the target's; role
    // selects the sign of the b/2 shift and (format 1 only) which side of
    // the lattice out-of-bounds indices are dropped on.
    void readWilsonLineText(
        const std::string &fileName, Parameters *param, NucleusRole role,
        std::vector<Matrix> &U);
    void readWilsonLineBinary(
        const std::string &fileName, Parameters *param, NucleusRole role,
        std::vector<Matrix> &U);

    Matrix getUfromExponent(std::vector<double> &Q);
    bool findUInForwardLightcone(
        Matrix &U1, Matrix &U2, Matrix &Usol, std::uint64_t retrySeed);

    void readInNucleusConfigs(
        const int nucleusA, const int lightNucleusOption,
        const int polarizationFlag, const double polJz,
        vector<vector<float>> &nucleonPosArr, Parameters *param);
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

// Parameters.h is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.

#ifndef SRC_PARAMETERS_H_
#define SRC_PARAMETERS_H_

#include <string>
#include <vector>

class Parameters {
  private:
    int subNucleonParamType_;
    int subNucleonParamSet_;
    std::vector<std::vector<float>> posteriorParamSets_;
    std::vector<std::vector<float>> posteriorParamSetsNq3_;

    double myPI_;
    double myhbarc_;

    int size_;        // the length of the lattice (make it 2^n, with n integer)
    int sizeOutput_;  // the length of the lattice for the output data
                     // (sizeOutput
                     // <= size!)
    int etaSizeOutput_;  // the length of the lattice in rapidity for the output
                        // data
    double detaOutput_;  // step size in rapidity for the output data
    int runningCoupling_;  // switch to decide if alpha_s should run (0 constant
                          // alpha_s, 1 running coupling)
    int useTimeForSeed_;  // decide if the system time should be used to generate
                         // a seed (1) or not (0)
    int useSeedList_;     // read random seeds from a file if set to (1) - this
                         // overwrites the 'use time for seed' setting
    unsigned long long int
        seed_;  // random seed that's added to the current time
                // to generate the full seed (or the full seed,
                // depending on the value of getUseTimeforSeed())
    int Ny_;     // longitudinal 'resolution' (see Lappi, Eur. Phys. J. C55,285)
    double g2mu_;       // g^2 mu [in lattice units]
    int mode_;  // mode: (1) run the evolution, (2) analysis with files from disk
    int runWithQs_;  // set whether alpha_s should run with the maximum(2),
                    // average (1) or minimum(0) of Q_s from nucleus A and B
    int runWithkt_;  // set whether alpha_s should run kt (1) or not (0) - if
                    // this is set it overwrites any running with Q_s
    int runWithLocalQs_;  // set whether alpha_s should run with the local Q_s
                         // from nucleus A and B (1) or the average (0), both
                         // use settings from runWithQs
    double runWithThisFactorTimesQs_;  // set the factor in front of Q_s under
                                      // the log in alpha_s
    double g_;  // coupling g needed in the initU3 where g^2mu does not scale out
    double m_;  // mass term in GeV to cut off the Coulomb tail - should be of
               // the order of \Lambda_QCD = 0.2 GeV
    double Jacobianm_;   // mass term in GeV in the Jacobian going from y to eta
    double QsmuRatio_;   // ratio between Qs and mu: Q_s = QsmuRatio * g^2 mu for
                        // nucleus A
    double QsmuRatioB_;  // ratio between Qs and mu: Q_s = QsmuRatio * g^2 mu for
                        // nucleus B
    double rapidityA_;  // rapidity to use when getting Q_s from IPSat.
                        // Basically to pick x for now
    double rapidityB_;  // rapidity to use when getting Q_s from IPSat.
                        // Basically to pick x for now
    int usePseudoRapidity_;  // if selected (1) the variable 'rapidity' will
                            // contain the pseudorapidity and the right
                            // conversion will be done (incl. Jacobian)
    double averageQs_;  // the average Q_s (maximum of nucleus A and B) used as
                       // scale for running coupling
    double averageQsAvg_;  // the average Q_s (average of nucleus A and B) used
                          // as scale for running coupling
    double averageQsmin_;  // the average Q_s (minimum of nucleus A and B) used
                          // as scale for running coupling
    double
        alphas_;  // the alpha_s computed at the scale given by the average Q_s
    int writeOutputs_;  // decide whether to write (1) or not write (0) large
                       // output files (like hydro input data)
    int writeEpsilonUHydro_;  // run the flow-velocity/hydro-output calculation
                             // (1) or write only Tmunu at measurement times (0)
    int writeTmunuBinary_;    // write Tmunu as compact binary .ipgt (1) or
                             // formatted text .dat (0)
    int writeOutputsToHDF5_;  // decide whether to write (1) or not write (0)
                             // output files to one hdf5 file
    int writeWilsonLines_;  // decide whether to write (1) in text or (2)
                           // in binary format or not write (0) generated
                           // Wilson lines (before any evolution)
    int readInitialWilsonLines_;  // decide wheter to generate initial Wilson
                                 // lines (0), or read these in plain text (1)
                                 // or in binary format (2)
    unsigned long long int randomSeed_;  // stores the random seed used (so the
                                        // event can be reproduced)
    std::string
        NucleusQsTableFileName_;  // the file name for the table containing Qs^2
                                 // as a function of Y and Qs^2(Y=0)
    double BG_;  // the width of the Gaussian describing the shape of the proton
                // in GeV^(-2)
    double BGq_;     // the mean width of the Gaussian describing the shape of
                     // a constituent quark in GeV^(-2)
    double BGqVar_;  // the variance of the Gaussian width describing the shape
    double omega_;
    // of a constituent quark in GeV^(-4)
    double dq_min_;  // the minimum distance between valence quarks [fm]
    double muZero_;   // mu_0 in the running coupling (makes it infrared finite)
    double c_;  // determines how smooth the cutoff in the running coupling is
    double roots_;  // square root of s: center of mass energy of the collision
                   // in GeV
    int useFluctuatingx_;  // switch to determine if the rapidity value in the
                          // input file should always be used (0) or if x should
                          // fluctuate as the local Q_s
    // x = Q_s*beta/roots (1). the value of roots is only used when this is set
    // to 1.
    double xFromThisFactorTimesQs_;  // set the factor beta in x = Q_s*beta/roots
    double Tpp_;  // This is the convolution of two T_p's to be used in the
                 // weight for different impact parameters
    // T_pp (b_T) = \sum \delta^2 x_T T_p(x_T) T_p(x_T-b_T)
    int inverseQsForMaxTime_;  // use 1/Q_s as the maximal evolution time (1) or
                              // use the manually entered maximal evolution time
                              // (0)
    double area_;  // area of the initial interaction region
    double Psi_;  // the initial angle Psi_2 that determines the event-plane
                 // (geometric/spatial one)
    // Glauber parameters:
    double SigmaNN_;      // nucleon-nucleon cross section
    double b_;            // impact parameter
    double bmin_;         // minimum impact parameter to sample from
    double bmax_;         // maximum impact parameter to sample to
    double phiRP_;       // the reaction plane angle
    int linearb_;         // sample b from a linear distribution if 1, uniform
                         // distribution otherwise
    std::string Target_;  // target nucleus' name
    std::string Projectile_;  // projectile nucleus' name
    double L_;                // lattice size in fm
    double LOutput_;          // lattice size for the output in fm
    int useNucleus_;  // use nuclei (1) or a constant g^2mu distribution over the
                     // lattice
    int lightNucleusOption_;  // for light nuclei (carbon, oxygen): 1:
                             // Woods-Saxon; 2: variational MC; 3: alpha
                             // clusters

    int polarizationFlagProjectile_;  // 0: unpolarized; 1: longitudinal
                                      // polarized; 2: transverse
    int polarizationFlagTarget_;  // 0: unpolarized; 1: longitudinal polarized;
                                  // 2: transverse
    double polJzProjectile_;      // The Jz polarization of the projectile
    double polJzTarget_;          // The Jz polarization of the target

    int useGaussian_;        // use a Gaussian profile on top of the constant
                            // background
    double dtau_;            // time step in lattice units
    double maxtime_;         // maximal evolution time in fm/c
    int Npart_;              // Number of participants
    int averageOverNuclei_;  // average over this many nuclei to get a smooth(er)
                            // distribution
    int nucleonPositionsFromFile_;  // switch to determine whether to sample
                                   // nucleon positions (0) or read them from a
                                   // file (1)
    int useFixedNpart_;  // if 0 do not demand a given N_part, if >1 sample the
                        // initial configuration until the given N_part is
                        // reached
    int smearQs_;  // decide whether to smear Q_s using a Poisson distribution
                  // around its mean at every x_T (1) or not (0)
    double smearingWidth_;  // width of the Gaussian smearing around the mean
                           // g^2mu^2
    int gaussianWounding_;  // use hard sphere profile (0) or Gaussian cross
                           // section (1) to determine whether a nucleon is
                           // wounded
    int MPIrank_;           // MPI rank
    int MPIsize_;           // MPI number of cores
    int event_id_;
    int success_;  // no collision happened (0) or collision happened (1) - used
                  // to restart if there was no collision
    int readMultFromFile_;  // if set, the gluon distribution as a function of
                           // k_T is read from file and the integrated rate
                           // computed
    double
        rmax_;  // radius at which we cut distribution for each nucleon (in fm)
    double protonAnisotropy_;  // anisotropy of the proton thickness function: xi
                              // in Exp[-(x^2 + xi y^2)/2/B]/2/Pi/B Sqrt[xi] -
                              // as a first test
    int useConstituentQuarkProton_;  // if >0, use proton made up of
                                    // useConstituentQuarkProton constituent
                                    // quarks.
    double NqBase_;
    double NqFluc_;
    int useSmoothNucleus_;  // if 1, use a smooth Woods-Saxon distribution for a
                           // heavy nucleus
    int shiftConstituentQuarkProtonOrigin_;  // if 1, move constituent quark
                                            // center of mass to origin
    double UVdamp_;                          // UV damping parameter
    int minimumQs2ST_;  // if >0 this will excludes events with Qs_min^2 S_T <
                       // minimumQs2ST. Can be used to trigger on high
                       // multiplicity events.
    double R_WS_, a_WS_;
    double beta2_;  // value of deformation parameter beta2 to test sensitivity
                   // in Uranium
    double beta3_, beta4_, gamma_;
    double d_min_;
    bool setWSDeformParams_, force_dmin_flag_;
    double WSdR_np_, WSda_np_;

    bool rotateReactionPlane_;  // flag to randomly rotate the event reaction
                                // plane

    bool computeGluonMultiplicity_;  // flag to compute gluonMultiplicity

    bool useJIMWLK_;  // flag to use JIMWLK evolution
    bool simpleLangevin_;

    double jimwlk_alphas_;  // 0 = running coupling, positive value = fixed
                           // coupling
    double m_jimwlk_;
    double mu0_jimwlk_;
    double LambdaQCD_jimwlk_;
    // int steps_jimwlk;
    // int measureSteps_jimwlk;
    double ds_jimwlk_;
    double x0_jimwlk_;  // Bjorken-x at the initial condition of the JIMLWK
                       // evolution

    double jimwlk_x1_;  // Bjorken x for the nucleus A (projectile)
    double jimwlk_x2_;  // Bjorken x for the nucleus B (target)
    bool saveSnapshots_;
    std::vector<double> xSnapshotList_;

  public:
    // constructor:
    Parameters() {};

    // functions to access the private variables:
    void setSubNucleonParamType(int paramType) {
        subNucleonParamType_ = paramType;
    }
    int getSubNucleonParamType() const { return (subNucleonParamType_); }
    void setSubNucleonParamSet(int paramSet) { subNucleonParamSet_ = paramSet; }
    int getSubNucleonParamSet() const { return (subNucleonParamSet_); }
    void setSeed(unsigned long long int x) { seed_ = x; }
    unsigned long long int getSeed() { return seed_; }
    void setNy(int x) { Ny_ = x; }
    int getNy() { return Ny_; }
    void setSize(int x) { size_ = x; }
    int getSize() { return size_; }
    void setProtonAnisotropy(double x) { protonAnisotropy_ = x; }
    double getProtonAnisotropy() { return protonAnisotropy_; }

    void setSizeOutput(int x) { sizeOutput_ = x; }
    int getSizeOutput() { return sizeOutput_; }
    void setEtaSizeOutput(int x) { etaSizeOutput_ = x; }
    int getEtaSizeOutput() { return etaSizeOutput_; }
    void setDetaOutput(double x) { detaOutput_ = x; }
    double getDetaOutput() { return detaOutput_; }

    void setAverageOverNuclei(int x) { averageOverNuclei_ = x; }
    int getAverageOverNuclei() { return averageOverNuclei_; }
    void setg2mu(double x) { g2mu_ = x; }
    double getg2mu() { return g2mu_; }
    void setMode(int x) { mode_ = x; };
    int getMode() { return mode_; }
    void setRunningCoupling(int x) { runningCoupling_ = x; };
    int getRunningCoupling() { return runningCoupling_; }
    void setg(double x) { g_ = x; }
    double getg() { return g_; }
    void setSigmaNN(double x) { SigmaNN_ = x; }
    double getSigmaNN() { return SigmaNN_; }
    void setb(double x) { b_ = x; }
    double getb() const { return b_; }
    void setPhiRP(double x) { phiRP_ = x; }
    double getPhiRP() const { return phiRP_; }
    void setbmin(double x) { bmin_ = x; }
    double getbmin() { return bmin_; }
    void setbmax(double x) { bmax_ = x; }
    double getbmax() { return bmax_; }
    void setTarget(std::string x) { Target_ = x; }
    std::string getTarget() { return Target_; }
    void setProjectile(std::string x) { Projectile_ = x; }
    std::string getProjectile() { return Projectile_; }
    void setL(double x) { L_ = x; }
    double getL() { return L_; }
    void setLOutput(double x) { LOutput_ = x; }
    double getLOutput() { return LOutput_; }
    void setm(double x) { m_ = x; }
    double getm() { return m_; }
    void setJacobianm(double x) { Jacobianm_ = x; }
    double getJacobianm() { return Jacobianm_; }
    void setQsmuRatio(double x) { QsmuRatio_ = x; }
    double getQsmuRatio() { return QsmuRatio_; }
    void setQsmuRatioB(double x) { QsmuRatioB_ = x; }
    double getQsmuRatioB() { return QsmuRatioB_; }
    void setRapidityA(double x) { rapidityA_ = x; }
    double getRapidityA() const { return rapidityA_; }
    void setRapidityB(double x) { rapidityB_ = x; }
    double getRapidityB() const { return rapidityB_; }
    double getRapidity() const { return (rapidityA_ + rapidityB_) / 2.; }
    void setMaxtime(double x) { maxtime_ = x; }
    double getMaxtime() { return maxtime_; }
    void setdtau(double x) { dtau_ = x; }
    double getdtau() { return dtau_; }
    void setNpart(int x) { Npart_ = x; };
    int getNpart() { return Npart_; }
    void setAverageQs(double x) { averageQs_ = x; }
    double getAverageQs() { return averageQs_; }
    void setAverageQsAvg(double x) { averageQsAvg_ = x; }
    double getAverageQsAvg() { return averageQsAvg_; }
    void setAverageQsmin(double x) { averageQsmin_ = x; }
    double getAverageQsmin() { return averageQsmin_; }
    void setalphas(double x) { alphas_ = x; }
    double getalphas() { return alphas_; }
    void setRandomSeed(unsigned long long int x) { randomSeed_ = x; };
    unsigned long long int getRandomSeed() { return randomSeed_; }
    void setUseTimeForSeed(int x) { useTimeForSeed_ = x; };
    int getUseTimeForSeed() { return useTimeForSeed_; }
    void setUseSeedList(int x) { useSeedList_ = x; };
    int getUseSeedList() { return useSeedList_; }
    void setNucleusQsTableFileName(std::string x) {
        NucleusQsTableFileName_ = x;
    }
    std::string getNucleusQsTableFileName() { return NucleusQsTableFileName_; }
    void setBG(double x) { BG_ = x; }
    double getBG() { return BG_; }
    void setBGq(double x) { BGq_ = x; }
    double getBGq() { return BGq_; }
    void setBGqVar(double BGqVar) { BGqVar_ = BGqVar; }
    double getBGqVar() { return BGqVar_; }
    void setOmega(double x) { omega_ = x; }
    double getOmega() const { return omega_; }
    void setDqmin(double dq_min) { dq_min_ = dq_min; }
    double getDqmin() { return dq_min_; }
    void setMuZero(double x) { muZero_ = x; }
    double getMuZero() { return muZero_; }
    void setc(double x) { c_ = x; }
    double getc() { return c_; }
    void setRoots(double x) { roots_ = x; }
    double getRoots() { return roots_; }
    void setUseFluctuatingx(int x) { useFluctuatingx_ = x; }
    int getUseFluctuatingx() { return useFluctuatingx_; }
    void setRunWithThisFactorTimesQs(double x) {
        runWithThisFactorTimesQs_ = x;
    };
    double getRunWithThisFactorTimesQs() { return runWithThisFactorTimesQs_; }
    void setxFromThisFactorTimesQs(double x) { xFromThisFactorTimesQs_ = x; };
    double getxFromThisFactorTimesQs() { return xFromThisFactorTimesQs_; }
    void setTpp(double x) { Tpp_ = x; }
    double getTpp() { return Tpp_; }
    void setUseFixedNpart(int x) { useFixedNpart_ = x; }
    int getUseFixedNpart() { return useFixedNpart_; }
    void setArea(double x) { area_ = x; }
    double getArea() { return area_; }
    void setPsi(double x) { Psi_ = x; }
    double getPsi() { return Psi_; }
    void setSmearingWidth(double x) { smearingWidth_ = x; }
    double getSmearingWidth() { return smearingWidth_; }
    void setMPIRank(int x) { MPIrank_ = x; }
    int getMPIRank() { return MPIrank_; }
    void setEventId(int x) { event_id_ = x; }
    int getEventId() { return event_id_; }
    void setMPISize(int x) { MPIsize_ = x; }
    int getMPISize() { return MPIsize_; }
    void setSuccess(int x) { success_ = x; }
    int getSuccess() { return success_; }
    void setRmax(double x) { rmax_ = x; }
    double getRmax() { return rmax_; }
    void setUVdamp(double x) { UVdamp_ = x; }
    double getUVdamp() { return UVdamp_; }
    void setSetWSDeformParams(int x) {
        if (x == 0)
            setWSDeformParams_ = false;
        else
            setWSDeformParams_ = true;
    }
    bool getSetWSDeformParams() const { return setWSDeformParams_; }
    void setR_WS(double x) { R_WS_ = x; }
    double getR_WS() const { return R_WS_; }
    void setA_WS(double x) { a_WS_ = x; }
    double getA_WS() const { return a_WS_; }
    void setBeta2(double x) { beta2_ = x; }
    double getBeta2() const { return beta2_; }
    void setBeta3(double x) { beta3_ = x; }
    double getBeta3() const { return beta3_; }
    void setBeta4(double x) { beta4_ = x; }
    double getBeta4() const { return beta4_; }
    void setGamma(double x) { gamma_ = x; }
    double getGamma() const { return gamma_; }
    void setDmin(double x) { d_min_ = x; }
    double getDmin() const { return d_min_; }
    void setForceDmin(int x) {
        if (x == 0)
            force_dmin_flag_ = false;
        else
            force_dmin_flag_ = true;
    }
    bool getForceDmin() const { return (force_dmin_flag_); }
    void setWSdR_np(double dR_np) { WSdR_np_ = dR_np; }
    double getWSdR_np() const { return WSdR_np_; }
    void setWSda_np(double da_np) { WSda_np_ = da_np; }
    double getWSda_np() const { return WSda_np_; }

    void setRotateReactionPlane(int iflag) {
        if (iflag == 0) {
            rotateReactionPlane_ = false;
        } else {
            rotateReactionPlane_ = true;
        }
    }
    bool getRotateReactionPlane() const { return rotateReactionPlane_; }

    // switches:
    void setUseNucleus(int x) { useNucleus_ = x; };
    int getUseNucleus() { return useNucleus_; }
    void setUseGaussian(int x) { useGaussian_ = x; };
    int getUseGaussian() { return useGaussian_; }
    void setlightNucleusOption(int x) { lightNucleusOption_ = x; };
    int getlightNucleusOption() { return lightNucleusOption_; }
    void setPolarizationProjectile(int x) { polarizationFlagProjectile_ = x; };
    int getPolarizationProjectile() { return polarizationFlagProjectile_; }
    void setPolarizationTarget(int x) { polarizationFlagTarget_ = x; };
    int getPolarizationTarget() { return polarizationFlagTarget_; }
    void setPolarizationProjectileJz(double x) { polJzProjectile_ = x; };
    double getPolarizationProjectileJz() { return polJzProjectile_; }
    void setPolarizationTargetJz(double x) { polJzTarget_ = x; };
    double getPolarizationTargetJz() { return polJzTarget_; }
    void setRunWithQs(int x) { runWithQs_ = x; };
    int getRunWithQs() { return runWithQs_; }
    void setRunWithkt(int x) { runWithkt_ = x; };
    int getRunWithkt() { return runWithkt_; }
    void setRunWithLocalQs(int x) { runWithLocalQs_ = x; };
    int getRunWithLocalQs() { return runWithLocalQs_; }
    void setLinearb(int x) { linearb_ = x; };
    int getLinearb() { return linearb_; }
    void setWriteOutputs(int x) { writeOutputs_ = x; };
    int getWriteOutputs() { return writeOutputs_; }
    void setWriteEpsilonUHydro(int x) { writeEpsilonUHydro_ = x; };
    int getWriteEpsilonUHydro() { return writeEpsilonUHydro_; }
    void setWriteTmunuBinary(int x) { writeTmunuBinary_ = x; };
    int getWriteTmunuBinary() { return writeTmunuBinary_; }
    void setWriteOutputsToHDF5(int x) { writeOutputsToHDF5_ = x; };
    int getWriteOutputsToHDF5() { return writeOutputsToHDF5_; }
    void setWriteWilsonLines(int x) { writeWilsonLines_ = x; }
    int getWriteWilsonLines() { return writeWilsonLines_; }
    void setReadInitialWilsonLines(int x) { readInitialWilsonLines_ = x; }
    int getReadInitialWilsonLines() { return readInitialWilsonLines_; }
    void setNucleonPositionsFromFile(int x) { nucleonPositionsFromFile_ = x; }
    int getNucleonPositionsFromFile() { return nucleonPositionsFromFile_; }
    void setInverseQsForMaxTime(int x) { inverseQsForMaxTime_ = x; };
    int getInverseQsForMaxTime() { return inverseQsForMaxTime_; }
    void setSmearQs(int x) { smearQs_ = x; }
    int getSmearQs() { return smearQs_; }
    void setReadMultFromFile(int x) { readMultFromFile_ = x; }
    int getReadMultFromFile() { return readMultFromFile_; }
    void setGaussianWounding(int x) { gaussianWounding_ = x; }
    int getGaussianWounding() { return gaussianWounding_; }
    void setUsePseudoRapidity(int x) { usePseudoRapidity_ = x; }
    int getUsePseudoRapidity() { return usePseudoRapidity_; }
    void setUseConstituentQuarkProton(int x) { useConstituentQuarkProton_ = x; }
    int getUseConstituentQuarkProton() { return useConstituentQuarkProton_; }
    void setNqBase(double NqBase) { NqBase_ = NqBase; }
    double getNqBase() { return NqBase_; }
    void setNqFluc(double NqFluc) { NqFluc_ = NqFluc; }
    double getNqFluc() { return NqFluc_; }
    void setUseSmoothNucleus(int x) { useSmoothNucleus_ = x; }
    int getUseSmoothNucleus() { return useSmoothNucleus_; }
    void setShiftConstituentQuarkProtonOrigin(int x) {
        shiftConstituentQuarkProtonOrigin_ = x;
    }
    int getShiftConstituentQuarkProtonOrigin() {
        return shiftConstituentQuarkProtonOrigin_;
    }
    void setMinimumQs2ST(int x) { minimumQs2ST_ = x; }
    int getMinimumQs2ST() { return minimumQs2ST_; }

    void setComputeGluonMultiplicity(int x) {
        if (x == 0) {
            computeGluonMultiplicity_ = false;
        } else {
            computeGluonMultiplicity_ = true;
        }
    }
    bool getComputeGluonMultiplicity() const {
        return computeGluonMultiplicity_;
    }

    void loadPosteriorParameterSetsFromFile(
        std::string posteriorFileName,
        std::vector<std::vector<float>> &ParamSet);
    void loadPosteriorParameterSets(const int itype);
    void setParamsWithPosteriorParameterSet(const int itype, int iset);

    // JIMWLK functions
    void setm_jimwlk(double x) { m_jimwlk_ = x; };
    double getm_jimwlk() { return m_jimwlk_; }
    void setMu0_jimwlk(double x) { mu0_jimwlk_ = x; }
    double getMu0_jimwlk() { return mu0_jimwlk_; }
    void setSimpleLangevin(int x) {
        if (x == 0) {
            simpleLangevin_ = false;
        } else {
            simpleLangevin_ = x;
        }
    }
    bool getSimpleLangevin() const { return simpleLangevin_; }
    void setLambdaQCD_jimwlk(double x) { LambdaQCD_jimwlk_ = x; }
    double getLambdaQCD_jimwlk() { return LambdaQCD_jimwlk_; }
    void setJimwlk_x_projectile(double x) { jimwlk_x1_ = x; }
    double getJimwlk_x_projectile() { return jimwlk_x1_; }
    void setJimwlk_x_target(double x) { jimwlk_x2_ = x; }
    double getJimwlk_x_target() { return jimwlk_x2_; }
    // void setMeasureSteps_jimwlk(int x) { measureSteps_jimwlk = x; };
    // int getMeasureSteps_jimwlk() { return measureSteps_jimwlk; }
    void setDs_jimwlk(double x) { ds_jimwlk_ = x; }
    double getDs_jimwlk() { return ds_jimwlk_; }
    void setJimwlk_alphas(double as) { jimwlk_alphas_ = as; }
    double getJimwlk_alphas() { return jimwlk_alphas_; }
    void setJimwlk_x0(double x) { x0_jimwlk_ = x; }
    double getJimwlk_x0() { return x0_jimwlk_; }
    bool getUseJIMWLK() const { return useJIMWLK_; }
    void setUseJIMWLK(int x) {
        if (x == 0) {
            useJIMWLK_ = false;
        } else {
            useJIMWLK_ = true;
        }
    }
    void setxSnapshotList(std::vector<double> xList) { xSnapshotList_ = xList; }
    std::vector<double> getxSnapshotList() { return xSnapshotList_; }
    void setSaveSnapshots(int x) {
        if (x == 0)
            saveSnapshots_ = false;
        else
            saveSnapshots_ = true;
    }
    bool getSaveSnapshots() { return saveSnapshots_; }
};
#endif  // SRC_PARAMETERS_H_

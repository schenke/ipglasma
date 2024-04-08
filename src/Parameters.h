// Parameters.h is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.

#ifndef Parameters_H
#define Parameters_H

#include <string>
#include <vector>

class Parameters {
private:
  int subNucleonParamType_;
  int subNucleonParamSet_;
  std::vector<std::vector<float>> posteriorParamSets_;
  std::vector<std::vector<float>> posteriorParamSetsNq3_;
  // switches:
  int initMethod;

  double myPI;
  double myhbarc;

  int A;          // number of nucleons;
  int Nc;         // number of colors (SU(Nc))
  int size;       // the length of the lattice (make it 2^n, with n integer)
  int sizeOutput; // the length of the lattice for the output data (sizeOutput
                  // <= size!)
  int etaSizeOutput;   // the length of the lattice in rapidity for the output
                       // data
  double detaOutput;   // step size in rapidity for the output data
  int runningCoupling; // switch to decide if alpha_s should run (0 constant
                       // alpha_s, 1 running coupling)
  double R;
  int useTimeForSeed; // decide if the system time should be used to generate a
                      // seed (1) or not (0)
  int useSeedList;    // read random seeds from a file if set to (1) - this
                      // overwrites the 'use time for seed' setting
  unsigned long long int seed; // random seed that's added to the current time
                               // to generate the full seed (or the full seed,
                               // depending on the value of getUseTimeforSeed())
  double ds;                   // 'time' step
  int Ny;      // longitudinal 'resolution' (see Lappi, Eur. Phys. J. C55,285)
  double g2mu; // g^2 mu [in lattice units]
  double Qs;   // Q_s, to be dynamically determined
  int steps;   // number of rapidity steps
  int measureSteps; // number of steps in interval between measurements
  int mode; // mode: (1) run the evolution, (2) analysis with files from disk
  int runWithQs; // set whether alpha_s should run with the maximum(2), average
                 // (1) or minimum(0) of Q_s from nucleus A and B
  int runWithkt; // set whether alpha_s should run kt (1) or not (0) - if this
                 // is set it overwrites any running with Q_s
  int runWithLocalQs; // set whether alpha_s should run with the local Q_s from
                      // nucleus A and B (1) or the average (0), both use
                      // settings from runWithQs
  double runWithThisFactorTimesQs; // set the factor in front of Q_s under the
                                   // log in alpha_s
  double mu0; // cutoff to avoid the Landau pole in the 1-loop running coupling
              // expression
  double LambdaQCD; // LambdaQCD in units of g^2mu
  double g; // coupling g needed in the initU3 where g^2mu does not scale out
  double
      kappa4Factor; // factor that multiplies the ratio of kappa4/(g^2 mu^2)^3
  double m; // mass term in GeV to cut off the Coulomb tail - should be of the
            // order of \Lambda_QCD = 0.2 GeV
  double Jacobianm;  // mass term in GeV in the Jacobian going from y to eta
  double QsmuRatio;  // ratio between Qs and mu: Q_s = QsmuRatio * g^2 mu for
                     // nucleus A
  double QsmuRatioB; // ratio between Qs and mu: Q_s = QsmuRatio * g^2 mu for
                     // nucleus B
  double rapidity; // rapidity to use when getting Q_s from IPSat. Basically to
                   // pick x for now
  int usePseudoRapidity; // if selected (1) the variable 'rapidity' will contain
                         // the pseudorapidity and the right conversion will be
                         // done (incl. Jacobian)
  double averageQs;      // the average Q_s (maximum of nucleus A and B) used as
                         // scale for running coupling
  double averageQsAvg;   // the average Q_s (average of nucleus A and B) used as
                         // scale for running coupling
  double averageQsmin;   // the average Q_s (minimum of nucleus A and B) used as
                         // scale for running coupling
  double alphas; // the alpha_s computed at the scale given by the average Q_s
  double xExponent; // - exponent with which Q_s grows with x (usually 0.31 in
                    // IP-Sat for nuclei)
  int writeOutputs; // decide whether to write (1) or not write (0) large output
                    // files (like hydro input data)
  int writeOutputsToHDF5; // decide whether to write (1) or not write (0) output
                          // files to one hdf5 file
  int writeEvolution;     // decide whether to write (1) or not write (0) time
                          // dependent quantities like the anisotropy
  int writeInitialWilsonLines; // decide whether to write (1) in text or (2) in
                               // binary format or not write (0) generated
                               // Wilson lines (before any evolution)
  int readInitialWilsonLines; // decide wheter to generate initial Wilson lines (0),
                              // or read these in plain text (1) or in binary format (2)
  unsigned long long int randomSeed; // stores the random seed used (so the
                                     // event can be reproduced)
  std::string NucleusQsTableFileName; // the file name for the table containing Qs^2
                                 // as a function of Y and Qs^2(Y=0)
  double BG;  // the width of the Gaussian describing the shape of the proton in
              // GeV^(-2)
  double BGq_; // the mean width of the Gaussian describing the shape of
               // a constituent quark in GeV^(-2)
  double BGqVar_;  // the variance of the Gaussian width describing the shape
                   // of a constituent quark in GeV^(-4)
  double dq_min_;  // the minimum distance between valence quarks [fm]
  double muZero;   // mu_0 in the running coupling (makes it infrared finite)
  double c;        // determines how smooth the cutoff in the running coupling is
  double
      roots; // square root of s: center of mass energy of the collision in GeV
  int useFluctuatingx; // switch to determine if the rapidity value in the input
                       // file should always be used (0) or if x should
                       // fluctuate as the local Q_s
  // x = Q_s*beta/roots (1). the value of roots is only used when this is set
  // to 1.
  double xFromThisFactorTimesQs; // set the factor beta in x = Q_s*beta/roots
  double Tpp; // This is the convolution of two T_p's to be used in the weight
              // for different impact parameters
  // T_pp (b_T) = \sum \delta^2 x_T T_p(x_T) T_p(x_T-b_T)
  int inverseQsForMaxTime; // use 1/Q_s as the maximal evolution time (1) or use
                           // the manually entered maximal evolution time (0)
  int useFatTails;         // if 1 use the student's t distribution instead of a
                   // Gaussian (0) to sample the rho's (standard deviation is
                   // still g2mu)
  double nu; // nu in the student's t distribution (used to produce fatter tails
             // than the Gaussian)
  double area;          // area of the initial interaction region
  double eccentricity2; // save the computed ellipticity to output together with
                        // S_T and dN/dy in the end
  double Psi; // the initial angle Psi_2 that determines the event-plane
              // (geometric/spatial one)
  // Glauber parameters:
  double SigmaNN; // nucleon-nucleon cross section
  double b;       // impact parameter
  double bmin;    // minimum impact parameter to sample from
  double bmax;    // maximum impact parameter to sample to
  int linearb; // sample b from a linear distribution if 1, uniform distribution
               // otherwise
  std::string Target;     // target nucleus' name
  std::string Projectile; // projectile nucleus' name
  double L;          // lattice size in fm
  double LOutput;    // lattice size for the output in fm
  int useNucleus;    // use nuclei (1) or a constant g^2mu distribution over the
                     // lattice
  int lightNucleusOption; // for light nuclei (carbon, oxygen): 1: Woods-Saxon;
                          // 2: variational MC; 3: alpha clusters
  int useGaussian; // use a Gaussian profile on top of the constant background
  double dtau;     // time step in lattice units
  double maxtime;  // maximal evolution time in fm/c
  int Npart;       // Number of participants
  int averageOverNuclei; // average over this many nuclei to get a smooth(er)
                         // distribution
  int nucleonPositionsFromFile; // switch to determine whether to sample nucleon
                                // positions (0) or read them from a file (1)
  int A1FromFile;    // if nuclei are read from file, store A value here
  int A2FromFile;    // if nuclei are read from file, store A value here
  int useFixedNpart; // if 0 do not demand a given N_part, if >1 sample the
                     // initial configuration until the given N_part is reached
  double rnp;        // distance between proton and neutron in the transverse
                     // projection of the deuteron
  int smearQs;       // decide whether to smear Q_s using a Poisson distribution
                     // around its mean at every x_T (1) or not (0)
  double
      smearingWidth; // width of the Gaussian smearing around the mean g^2mu^2
  int gaussianWounding; // use hard sphere profile (0) or Gaussian cross section
                        // (1) to determine whether a nucleon is wounded
  int MPIrank;          // MPI rank
  int MPIsize;          // MPI number of cores
  int event_id;
  int success; // no collision happened (0) or collision happened (1) - used to
               // restart if there was no collision
  int readMultFromFile; // if set, the gluon distribution as a function of k_T
                        // is read from file and the integrated rate computed
  double rmax; // radius at which we cut distribution for each nucleon (in fm)
  double protonAnisotropy; // anisotropy of the proton thickness function: xi in
                           // Exp[-(x^2 + xi y^2)/2/B]/2/Pi/B Sqrt[xi] - as a
                           // first test
  int useConstituentQuarkProton; // if >0, use proton made up of
                                 // useConstituentQuarkProton constituent
                                 // quarks.
  double NqBase_;
  double NqFluc_;
  int useSmoothNucleus; // if 1, use a smooth Woods-Saxon distribution for a
                        // heavy nucleus
  int shiftConstituentQuarkProtonOrigin; // if 1, move constituent quark center
                                         // of mass to origin
  double UVdamp;                         // UV damping parameter
  int minimumQs2ST;        // if >0 this will excludes events with Qs_min^2 S_T < minimumQs2ST. Can be used to trigger on high multiplicity events.
  double R_WS_, a_WS_;
  double beta2;        // value of deformation parameter beta2 to test sensitivity in Uranium
  double beta3, beta4, gamma_;
  double d_min_;
  bool setWSDeformParams_, force_dmin_flag_;
  double WSdR_np_, WSda_np_;

public:
  // constructor:
  Parameters(){};

  // functions to access the private variables:
  void setSubNucleonParamType(int paramType) {
      subNucleonParamType_ = paramType;
  }
  int getSubNucleonParamType() const { return(subNucleonParamType_); }
  void setSubNucleonParamSet(int paramSet) {
      subNucleonParamSet_ = paramSet;
  }
  int getSubNucleonParamSet() const { return(subNucleonParamSet_); }
  void setSeed(unsigned long long int x) { seed = x; }
  unsigned long long int getSeed() { return seed; }
  void setA(int x) { A = x; }
  int getA() { return A; }
  void setNc(int x) { Nc = x; }
  int getNc() { return Nc; }
  void setNy(int x) { Ny = x; }
  int getNy() { return Ny; }
  void setSize(int x) { size = x; }
  int getSize() { return size; }
  void setProtonAnisotropy(double x) { protonAnisotropy = x; }
  double getProtonAnisotropy() { return protonAnisotropy; }

  void setSizeOutput(int x) { sizeOutput = x; }
  int getSizeOutput() { return sizeOutput; }
  void setEtaSizeOutput(int x) { etaSizeOutput = x; }
  int getEtaSizeOutput() { return etaSizeOutput; }
  void setDetaOutput(double x) { detaOutput = x; }
  double getDetaOutput() { return detaOutput; }

  void setAverageOverNuclei(int x) { averageOverNuclei = x; }
  int getAverageOverNuclei() { return averageOverNuclei; }
  void setR(double x) { R = x; }
  double getR() { return R; }
  void setDs(double x) { ds = x; }
  double getDs() { return ds; }
  void setg2mu(double x) { g2mu = x; }
  double getg2mu() { return g2mu; }
  void setQs(double x) { Qs = x; }
  double getQs() { return Qs; }
  void setSteps(int x) { steps = x; };
  int getSteps() { return steps; }
  void setMeasureSteps(int x) { measureSteps = x; };
  int getMeasureSteps() { return measureSteps; }
  void setMode(int x) { mode = x; };
  int getMode() { return mode; }
  void setRunningCoupling(int x) { runningCoupling = x; };
  int getRunningCoupling() { return runningCoupling; }
  void setMu0(double x) { mu0 = x; }
  double getMu0() { return mu0; }
  void setg(double x) { g = x; }
  double getg() { return g; }
  void setLambdaQCD(double x) { LambdaQCD = x; }
  double getLambdaQCD() { return LambdaQCD; }
  void setkappa4Factor(double x) { kappa4Factor = x; }
  double getkappa4Factor() { return kappa4Factor; }
  void setSigmaNN(double x) { SigmaNN = x; }
  double getSigmaNN() { return SigmaNN; }
  void setb(double x) { b = x; }
  double getb() { return b; }
  void setbmin(double x) { bmin = x; }
  double getbmin() { return bmin; }
  void setbmax(double x) { bmax = x; }
  double getbmax() { return bmax; }
  void setTarget(std::string x) { Target = x; }
  std::string getTarget() { return Target; }
  void setProjectile(std::string x) { Projectile = x; }
  std::string getProjectile() { return Projectile; }
  void setL(double x) { L = x; }
  double getL() { return L; }
  void setLOutput(double x) { LOutput = x; }
  double getLOutput() { return LOutput; }
  void setm(double x) { m = x; }
  double getm() { return m; }
  void setJacobianm(double x) { Jacobianm = x; }
  double getJacobianm() { return Jacobianm; }
  void setQsmuRatio(double x) { QsmuRatio = x; }
  double getQsmuRatio() { return QsmuRatio; }
  void setQsmuRatioB(double x) { QsmuRatioB = x; }
  double getQsmuRatioB() { return QsmuRatioB; }
  void setRapidity(double x) { rapidity = x; }
  double getRapidity() { return rapidity; }
  void setMaxtime(double x) { maxtime = x; }
  double getMaxtime() { return maxtime; }
  void setdtau(double x) { dtau = x; }
  double getdtau() { return dtau; }
  void setNpart(int x) { Npart = x; };
  int getNpart() { return Npart; }
  void setAverageQs(double x) { averageQs = x; }
  double getAverageQs() { return averageQs; }
  void setAverageQsAvg(double x) { averageQsAvg = x; }
  double getAverageQsAvg() { return averageQsAvg; }
  void setAverageQsmin(double x) { averageQsmin = x; }
  double getAverageQsmin() { return averageQsmin; }
  void setalphas(double x) { alphas = x; }
  double getalphas() { return alphas; }
  void setxExponent(double x) { xExponent = x; }
  double getxExponent() { return xExponent; }
  void setRandomSeed(unsigned long long int x) { randomSeed = x; };
  unsigned long long int getRandomSeed() { return randomSeed; }
  void setUseTimeForSeed(int x) { useTimeForSeed = x; };
  int getUseTimeForSeed() { return useTimeForSeed; }
  void setUseSeedList(int x) { useSeedList = x; };
  int getUseSeedList() { return useSeedList; }
  void setNucleusQsTableFileName(std::string x) { NucleusQsTableFileName = x; }
  std::string getNucleusQsTableFileName() { return NucleusQsTableFileName; }
  void setA1FromFile(int x) { A1FromFile = x; }
  int getA1FromFile() { return A1FromFile; }
  void setA2FromFile(int x) { A2FromFile = x; }
  int getA2FromFile() { return A2FromFile; }
  void setBG(double x) { BG = x; }
  double getBG() { return BG; }
  void setBGq(double x) { BGq_ = x; }
  double getBGq() { return BGq_; }
  void setBGqVar(double BGqVar) { BGqVar_ = BGqVar; }
  double getBGqVar() { return BGqVar_; }
  void setDqmin(double dq_min) { dq_min_ = dq_min; }
  double getDqmin() { return dq_min_; }
  void setMuZero(double x) { muZero = x; }
  double getMuZero() { return muZero; }
  void setc(double x) { c = x; }
  double getc() { return c; }
  void setRoots(double x) { roots = x; }
  double getRoots() { return roots; }
  void setUseFluctuatingx(int x) { useFluctuatingx = x; }
  int getUseFluctuatingx() { return useFluctuatingx; }
  void setRunWithThisFactorTimesQs(double x) { runWithThisFactorTimesQs = x; };
  double getRunWithThisFactorTimesQs() { return runWithThisFactorTimesQs; }
  void setxFromThisFactorTimesQs(double x) { xFromThisFactorTimesQs = x; };
  double getxFromThisFactorTimesQs() { return xFromThisFactorTimesQs; }
  void setTpp(double x) { Tpp = x; }
  double getTpp() { return Tpp; }
  void setNu(double x) { nu = x; }
  double getNu() { return nu; }
  void setUseFixedNpart(int x) { useFixedNpart = x; }
  int getUseFixedNpart() { return useFixedNpart; }
  void setArea(double x) { area = x; }
  double getArea() { return area; }
  void setEccentricity2(double x) { eccentricity2 = x; }
  double getEccentricity2() { return eccentricity2; }
  void setRnp(double x) { rnp = x; }
  double getRnp() { return rnp; }
  void setPsi(double x) { Psi = x; }
  double getPsi() { return Psi; }
  void setSmearingWidth(double x) { smearingWidth = x; }
  double getSmearingWidth() { return smearingWidth; }
  void setMPIRank(int x) { MPIrank = x; }
  int getMPIRank() { return MPIrank; }
  void setEventId(int x) { event_id = x; }
  int getEventId() { return event_id; }
  void setMPISize(int x) { MPIsize = x; }
  int getMPISize() { return MPIsize; }
  void setSuccess(int x) { success = x; }
  int getSuccess() { return success; }
  void setRmax(double x) { rmax = x; }
  double getRmax() { return rmax; }
  void setUVdamp(double x) { UVdamp = x; }
  double getUVdamp() { return UVdamp; }
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
  void setBeta2(double x) { beta2 = x; }
  double getBeta2() const { return beta2; }
  void setBeta3(double x) { beta3 = x; }
  double getBeta3() const { return beta3; }
  void setBeta4(double x) { beta4 = x; }
  double getBeta4() const { return beta4; }
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

  // switches:
  void setInitMethod(int x) { initMethod = x; }
  int getInitMethod() { return initMethod; }
  void setUseNucleus(int x) { useNucleus = x; };
  int getUseNucleus() { return useNucleus; }
  void setUseGaussian(int x) { useGaussian = x; };
  int getUseGaussian() { return useGaussian; }
  void setlightNucleusOption(int x) { lightNucleusOption = x; };
  int getlightNucleusOption() { return lightNucleusOption; }
  void setRunWithQs(int x) { runWithQs = x; };
  int getRunWithQs() { return runWithQs; }
  void setRunWithkt(int x) { runWithkt = x; };
  int getRunWithkt() { return runWithkt; }
  void setRunWithLocalQs(int x) { runWithLocalQs = x; };
  int getRunWithLocalQs() { return runWithLocalQs; }
  void setLinearb(int x) { linearb = x; };
  int getLinearb() { return linearb; }
  void setWriteOutputs(int x) { writeOutputs = x; };
  int getWriteOutputs() { return writeOutputs; }
  void setWriteOutputsToHDF5(int x) { writeOutputsToHDF5 = x; };
  int getWriteOutputsToHDF5() { return writeOutputsToHDF5; }
  void setWriteEvolution(int x) { writeEvolution = x; };
  int getWriteEvolution() { return writeEvolution; }
  void setWriteInitialWilsonLines(int x) { writeInitialWilsonLines = x; }
  int getWriteInitialWilsonLines() { return writeInitialWilsonLines; }
  void setReadInitialWilsonLines(int x) { readInitialWilsonLines = x; }
  int getReadInitialWilsonLines() { return readInitialWilsonLines; }
  void setNucleonPositionsFromFile(int x) { nucleonPositionsFromFile = x; }
  int getNucleonPositionsFromFile() { return nucleonPositionsFromFile; }
  void setInverseQsForMaxTime(int x) { inverseQsForMaxTime = x; };
  int getInverseQsForMaxTime() { return inverseQsForMaxTime; }
  void setUseFatTails(int x) { useFatTails = x; }
  int getUseFatTails() { return useFatTails; }
  void setSmearQs(int x) { smearQs = x; }
  int getSmearQs() { return smearQs; }
  void setReadMultFromFile(int x) { readMultFromFile = x; }
  int getReadMultFromFile() { return readMultFromFile; }
  void setGaussianWounding(int x) { gaussianWounding = x; }
  int getGaussianWounding() { return gaussianWounding; }
  void setUsePseudoRapidity(int x) { usePseudoRapidity = x; }
  int getUsePseudoRapidity() { return usePseudoRapidity; }
  void setUseConstituentQuarkProton(int x) { useConstituentQuarkProton = x; }
  int getUseConstituentQuarkProton() { return useConstituentQuarkProton; }
  void setNqBase(double NqBase) { NqBase_ = NqBase; }
  double getNqBase() { return NqBase_; }
  void setNqFluc(double NqFluc) { NqFluc_ = NqFluc; }
  double getNqFluc() { return NqFluc_; }
  void setUseSmoothNucleus(int x) { useSmoothNucleus = x; }
  int getUseSmoothNucleus() { return useSmoothNucleus; }
  void setShiftConstituentQuarkProtonOrigin(int x) {
    shiftConstituentQuarkProtonOrigin = x;
  }
  int getShiftConstituentQuarkProtonOrigin() {
    return shiftConstituentQuarkProtonOrigin;
  }
  void setMinimumQs2ST(int x) { minimumQs2ST = x; }
  int getMinimumQs2ST() { return minimumQs2ST; }

  void loadPosteriorParameterSetsFromFile(std::string posteriorFileName,
                                          std::vector<std::vector<float>> &ParamSet);
  void loadPosteriorParameterSets(const int itype);
  void setParamsWithPosteriorParameterSet(const int itype, int iset);
};
#endif // Parameters_H

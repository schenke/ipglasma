#ifndef SRC_GLAUBER_H_
#define SRC_GLAUBER_H_

#include <string>

#include "Random.h"

#define TOL (1.0e-6)
#define TINY (1.0e-10)
#define LIMIT 10000

enum class NucleusRole { Projectile, Target };

struct ReturnValue {
    double x;
    double y;
    double z;
    double phi;
    int collided;
    bool proton;
    // int acceptances;
};

struct Nucleus {
    std::string name;
    int A;
    int Z;
    int anumFunc;
    int anumFuncIntegrand;
    int densityFunc;
    double w_WS;
    double a_WS;
    double R_WS;
    double rho_WS;
    double beta2;
    double beta3;
    double beta4;
    double gamma;
    bool forceDminFlag;
    double d_min;
    double dR_np;
    double da_np;
};

struct Data {
    double sigmaNN;
    Nucleus target;
    Nucleus projectile;
    double sCutoff;
    int interMax;
    /* trap door */
};

class Glauber {
  private:
    double AnumR_, NuInS_S_;
    Nucleus *Nuc_WS_;
    Data glauberData_;
    double b_;  // impact parameter
    int currentA1_;
    int currentA2_;
    int currentZ1_;
    int currentZ2_;

  public:
    Glauber() {};
    ~Glauber() { remove("tmp.dat"); }

    int nucleusA1() const { return currentA1_; }
    int nucleusA2() const { return currentA2_; }
    int nucleusZ1() const { return currentZ1_; }
    int nucleusZ2() const { return currentZ2_; }
    const Data &getGlauberData() const { return glauberData_; }
    int isFile(char *file_name);
    void findNucleusData(
        Nucleus *nucleus, std::string name, bool setWSDeformParams, double R_WS,
        double a_WS, double beta2, double beta3, double beta4, double gamma,
        bool forceDminFlag, double d_min, double dR_np, double da_np);
    void printGlauberData();
    void printNucleusData(Nucleus *nucleus);
    int linearFindXorg(double x, double *Vx, int ymax);
    double fourPtInterpolate(
        double x, double *Vx, double *Vy, double h, int x_org);
    void makeCoeff(
        double *a, double *b, double *c, double *d, double *Vy, double h,
        int x_org);
    double vInterpolate(double x, double *Vx, double *Vy, int ymax);
    double *makeVx(double down, double up, int maxi_num);
    double *makeVy(double *vx, int maxi_num);
    double *readInVx(char *file_name, int maxi_num, int quiet);
    double *readInVy(char *file_name, int maxi_num, int quiet);

    double interNuPInSP(double s);
    double interNuTInST(double s);
    void calcRho(Nucleus *nucleus);
    double nuInS(double s);

    double anum3Fermi(double R_WS);
    double anum3FermiInt(double xi);
    double nuInt3Fermi(double xi);
    double anum3Gauss(double R_WS);
    double anum3GaussInt(double xi);
    double nuInt3Gauss(double xi);
    double anum2HO();
    double anum2HOInt(double xi);
    double nuInt2HO(double xi);
    double anumHulthen();
    double nuIntHulthen(double xi);

    double integral(int id, double down, double up, double tol, int *count);
    double qnc7(
        int id, double tol, double down, double dx, double *f_of,
        double pre_sum, double area, int *count);
    double oLSIntegrand(double s);
    double tAB();
    void initGlauber(
        double sigmaNN, std::string target, std::string projectile, double inb,
        bool setWSDeformParams, double R_WS, double a_WS, double beta2,
        double beta3, double beta4, double gamma, bool forceDminFlag,
        double d_min, double dR_np, double da_np, int imax);
    double areaTA(double x, double A);
    ReturnValue sampleTARejection(Random *random, NucleusRole nucleus);
};
#endif  // SRC_GLAUBER_H_

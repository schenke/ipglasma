#ifndef SRC_GLAUBER_H_
#define SRC_GLAUBER_H_

#include <string>

#include "Random.h"

#define TOL (1.0e-6)
#define tiny (1.0e-10)
#define limit 10000

struct ReturnValue {
    double x;
    double y;
    double z;
    double phi;
    int collided;
    bool proton;
    // int acceptances;
};

typedef struct nucleus {
    std::string name;
    double A;
    double Z;
    int AnumFunc;
    int AnumFuncIntegrand;
    int DensityFunc;
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
} Nucleus;

typedef struct data {
    double SigmaNN;
    Nucleus Target;
    Nucleus Projectile;
    double SCutOff;
    int InterMax;
    /* trap door */
} Data;

class Glauber {
  private:
  public:
    typedef double (*ptr_func)(double);

    double AnumR, NuInS_S;
    Nucleus *Nuc_WS;
    Data GlauberData;
    ptr_func tempFunc;
    double b;  // impact parameter
    double currentTAB;
    double currentA1;
    double currentA2;
    double currentZ1;
    double currentZ2;

    Glauber() {};
    ~Glauber() { remove("tmp.dat"); }

    int nucleusA1() const { return static_cast<int>(currentA1); }
    int nucleusA2() const { return static_cast<int>(currentA2); }
    int nucleusZ1() const { return static_cast<int>(currentZ1); }
    int nucleusZ2() const { return static_cast<int>(currentZ2); }
    int isFile(char *file_name);
    void findNucleusData(
        Nucleus *nucleus, std::string target, std::string file_name, int rank);
    void findNucleusData2(
        Nucleus *nucleus, std::string name, bool setWSDeformParams, double R_WS,
        double a_WS, double beta2, double beta3, double beta4, double gamma,
        bool force_dmin, double d_min, double dR_np, double da_np);
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
    double *readInVx(char *, int maxi_num, int quiet);
    double *readInVy(char *, int maxi_num, int quiet);

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
        double SigmaNN, std::string Target, std::string Projectile, double inb,
        bool setWSDeformParams, double R_WS, double a_WS, double beta2,
        double beta3, double beta4, double gamma, bool force_dmin, double d_min,
        double dR_np, double da_np, int imax);
    double areaTA(double x, double A);
    ReturnValue sampleTARejection(Random *random, int PorT);
};
#endif  // SRC_GLAUBER_H_

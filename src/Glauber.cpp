#include "Glauber.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "PhysConst.h"

using PhysConst::mbToFm2;
using std::string;

namespace {

/**
 * One row of Glauber::findNucleusData()'s built-in per-species table.
 */
struct NucleusTemplate {
    /// Species name, matched case-sensitively against
    /// Glauber::findNucleusData()'s \p name argument.
    const char *name;
    /// Mass number.
    int A;
    /// Atomic number.
    int Z;
    /// Woods-Saxon half-density radius [fm].
    double R_WS;
    /// Profile-shape parameter [dimensionless or 1/fm].
    double w_WS;
    /// Surface diffuseness [fm].
    double a_WS;
    /// Quadrupole deformation [dimensionless].
    double beta2;
    /// Octupole deformation [dimensionless].
    double beta3;
    /// Hexadecapole deformation [dimensionless].
    double beta4;
    /// Triaxiality angle [rad].
    double gamma;
    /// Density-function selector, assigned identically to
    /// `Nucleus::anumFunc`/`anumFuncIntegrand`/`densityFunc` (1=2HO or
    /// readFromFile, 2=3Gauss, 3=3Fermi, 8=Hulthen).
    int funcCode;
};

/// Built-in Woods-Saxon/density-profile parameters, one row per
/// supported species name; see NucleusTemplate for the field meanings.
/// The deformations of O and Xe are from FRDM(2012) \cite Moller:2015fba;
/// U's \f$\beta_2 = 0.28\f$ and the other parameters have no recorded
/// source.
const NucleusTemplate kNucleusTemplates[] = {
    {"Au", 197, 79, 6.37, 0, 0.535, -0.13, 0., -0.03, 0., 3},
    {"Pb", 208, 82, 6.62, 0., 0.546, 0.0, 0., 0.0, 0., 3},
    {"p", 1, 1, 1., 0., 1., 0.0, 0., 0.0, 0., 3},
    {"He3", 3, 2, 0, 0, 0, 0.0, 0., 0.0, 0., 1},
    {"He4", 4, 2, 0, 0, 0, 0.0, 0., 0.0, 0., 1},
    {"d", 2, 1, 1.0, 1.18, 0.228, 0.0, 0., 0.0, 0., 8},
    {"C", 12, 6, 2.44, 1.403, 1.635, 0.0, 0., 0.0, 0., 1},
    // beta2 and beta4 from FRDM(2012), arXiv:1508.06294
    {"O", 16, 8, 2.608, -0.051, 0.513, -0.01, 0., -0.122, 0., 3},
    {"Ne", 20, 10, 2.8, 0.0, 0.57, 0.0, 0., 0.0, 0., 3},
    {"Ne22", 22, 10, 2.782, 0.0, 0.549, 0.0, 0., 0.0, 0., 3},
    {"S", 32, 16, 2.54, 0.16, 2.191, 0.0, 0., 0.0, 0., 2},
    {"Ar", 40, 18, 3.61, 0.0, 0.516, 0.1668, 0., 0.00695, 0., 3},
    {"W", 184, 74, 6.51, 0, 0.535, 0.0, 0., 0.0, 0., 3},
    {"Al", 27, 13, 3.07, 0, 0.519, 0.0, 0., 0.0, 0., 3},
    {"Ca", 40, 20, 3.766, -0.161, 0.586, 0.0, 0., 0.0, 0., 3},
    {"Cu", 63, 29, 4.163, 0, 0.606, 0.162, 0., 0.006, 0., 3},
    {"Fe", 56, 26, 4.106, 0, 0.519, 0.0, 0., 0.0, 0., 3},
    {"Pt", 195, 78, 6.78, 0, 0.54, 0.0, 0., 0.0, 0., 3},
    // beta2 = 0.28: no recorded source
    {"U", 238, 92, 6.81, 0, 0.55, 0.28, 0., 0.093, 0., 3},
    {"Ru", 96, 44, 5.085, 0, 0.46, 0.158, 0., 0.0, 0., 3},
    {"Zr", 96, 40, 5.02, 0, 0.46, 0.0, 0., 0.0, 0., 3},
    // beta2 and beta4 from FRDM(2012), arXiv:1508.06294
    {"Xe", 129, 54, 5.42, 0, 0.57, 0.162, 0., -0.003, 0., 3},
};

}  // namespace

int speciesMassNumber(const std::string &name) {
    for (const auto &candidate : kNucleusTemplates) {
        if (name == candidate.name) return candidate.A;
    }
    return 0;
}

bool speciesHasDensityProfile(const std::string &name) {
    for (const auto &candidate : kNucleusTemplates) {
        if (name == candidate.name) {
            return candidate.A <= 2 || candidate.a_WS > 0.;
        }
    }
    return true;
}

void Glauber::findNucleusData(
    Nucleus *nucleus, string name, bool setWSDeformParams, double R_WS,
    double a_WS, double beta2, double beta3, double beta4, double gamma,
    bool forceDminFlag, double d_min, double dR_np, double da_np) {
    const NucleusTemplate *tmpl = nullptr;
    for (const auto &candidate : kNucleusTemplates) {
        if (name == candidate.name) {
            tmpl = &candidate;
            break;
        }
    }
    if (tmpl == nullptr) {
        messager_ << "[Glauber::findNucleusData]: Unknown nucleus name \""
                  << name << "\". Exiting.";
        messager_.flush("error");
        exit(1);
    }

    // Not set anywhere pre-table-driven-rewrite either (findNucleusData
    // never wrote it, and no caller set it separately), leaving it
    // permanently empty for printNucleusData -- currently uncalled -- to
    // print; set it here now that we have the matched name at hand anyway.
    nucleus->name = name;
    nucleus->A = tmpl->A;
    nucleus->Z = tmpl->Z;
    nucleus->R_WS = tmpl->R_WS;
    nucleus->w_WS = tmpl->w_WS;
    nucleus->a_WS = tmpl->a_WS;
    nucleus->beta2 = tmpl->beta2;
    nucleus->beta3 = tmpl->beta3;
    nucleus->beta4 = tmpl->beta4;
    nucleus->gamma = tmpl->gamma;
    nucleus->anumFunc = tmpl->funcCode;
    nucleus->anumFuncIntegrand = tmpl->funcCode;
    nucleus->densityFunc = tmpl->funcCode;

    nucleus->rho_WS = nucleus->R_WS;

    // Setting neutron skin to 0 by default
    nucleus->dR_np = 0.;
    nucleus->da_np = 0.;

    // the input describes a Woods-Saxon nucleus, whatever the species'
    // built-in profile is; the proton and the deuteron keep their own
    // description (a nucleon at the origin, the Hulthen wave function), so
    // that e.g. d+Au with the input parameters of Au keeps its deuteron
    if (setWSDeformParams && tmpl->A > 2) {
        nucleus->anumFunc = 3;
        nucleus->anumFuncIntegrand = 3;
        nucleus->densityFunc = 3;
        nucleus->w_WS = 0.;
        nucleus->R_WS = R_WS;
        nucleus->a_WS = a_WS;
        nucleus->beta2 = beta2;
        nucleus->beta3 = beta3;
        nucleus->beta4 = beta4;
        nucleus->gamma = gamma;
        nucleus->dR_np = dR_np;
        nucleus->da_np = da_np;
    }
    nucleus->forceDminFlag = forceDminFlag;
    nucleus->d_min = d_min;
}

void Glauber::printGlauberData() {
    messager_ << "[Glauber::printGlauberData]: sigmaNN = "
              << glauberData_.sigmaNN
              << ", interMax = " << glauberData_.interMax
              << ", sCutoff = " << glauberData_.sCutoff;
    messager_.flush("debug");
}

void Glauber::printNucleusData(Nucleus *nucleus) {
    messager_ << "[Glauber::printNucleusData]: Nucleus Name: " << nucleus->name
              << "\n"
              << " Nucleus.A = " << nucleus->A << "\n"
              << " Nucleus.Z = " << nucleus->Z << "\n"
              << " Nucleus.w_WS = " << nucleus->w_WS << "\n"
              << " Nucleus.a_WS = " << nucleus->a_WS << "\n"
              << " Nucleus.R_WS = " << nucleus->R_WS;
    messager_.flush("info");
}

int Glauber::linearFindXorg(double x, double *Vx, int ymax) {
    /* finds the first of the 4 points, x is between the second and the third */

    int x_org;
    double nx;

    nx = ymax * (x - Vx[0]) / (Vx[ymax] - Vx[0]);

    x_org = (int)nx;
    x_org -= 1;

    if (x_org <= 0)
        return 0;
    else if (x_org >= ymax - 3)
        return ymax - 3;
    else
        return x_org;

} /* Linear Find Xorg */

double Glauber::fourPtInterpolate(
    double x, double *Vx, double *Vy, double h, int x_org) {
    /* interpolating points are x_org, x_org+1, x_org+2, x_org+3 */
    /* cubic polynomial approximation */

    double a, bb, c, d, f;

    makeCoeff(&a, &bb, &c, &d, Vy, h, x_org);

    f = a * pow(x - Vx[x_org], 3.);
    f += bb * pow(x - Vx[x_org], 2.);
    f += c * (x - Vx[x_org]);
    f += d;

    return f;
}

void Glauber::makeCoeff(
    double *a, double *b, double *c, double *d, double *Vy, double h,
    int x_org) {
    double f0, f1, f2, f3;

    f0 = Vy[x_org];
    f1 = Vy[x_org + 1];
    f2 = Vy[x_org + 2];
    f3 = Vy[x_org + 3];

    *a = (-f0 + 3.0 * f1 - 3.0 * f2 + f3) / (6.0 * h * h * h);

    *b = (2.0 * f0 - 5.0 * f1 + 4.0 * f2 - f3) / (2.0 * h * h);

    *c = (-11.0 * f0 + 18.0 * f1 - 9.0 * f2 + 2.0 * f3) / (6.0 * h);

    *d = f0;
}

double Glauber::vInterpolate(double x, double *Vx, double *Vy, int ymax) {
    int x_org;
    double h;

    if ((x < Vx[0]) || (x > Vx[ymax])) {
        messager_ << "[Glauber::vInterpolate]: x = " << x
                  << " is outside the tabulated range [" << Vx[0] << ", "
                  << Vx[ymax] << "]. This should never happen -- exiting.";
        messager_.flush("error");
        exit(1);
    }

    /* we only deal with evenly spaced Vx */
    /* x_org is the first of the 4 points */

    x_org = linearFindXorg(x, Vx, ymax);

    h = (Vx[ymax] - Vx[0]) / ymax;

    return fourPtInterpolate(x, Vx, Vy, h, x_org);

} /* vInterpolate */

std::vector<double> Glauber::makeVx(double down, double up, int maxi_num) {
    std::vector<double> vx(maxi_num + 1);
    double dx = (up - down) / maxi_num;
    for (int i = 0; i <= maxi_num; i++) {
        vx[i] = dx * i;
    }
    return vx;

} /* makeVx */

std::vector<double> Glauber::makeVy(const double *vx, int maxi_num) {
    std::vector<double> vy(maxi_num + 1);

    for (int i = 0; i <= maxi_num; i++) {
        vy[i] = nuInS(vx[i]);
    }

    return vy;
} /* makeVy */

std::vector<double> Glauber::readInVx(
    char *file_name, int maxi_num, int quiet) {
    std::vector<double> vx(maxi_num + 1);
    double x;
    FILE *input;
    char s[120], sx[120];

    if (quiet == 1) {
        messager_ << "[Glauber::readInVx]: Reading in Vx from " << file_name
                  << " ...";
        messager_.flush("info");
    }

    input = fopen(file_name, "r");
    if (input == nullptr) {
        messager_ << "[Glauber::readInVx]: File " << file_name
                  << " not found. Exiting.";
        messager_.flush("error");
        exit(1);
    }
    if (fscanf(input, "%s", s) != 1) {
        messager_ << "[Glauber::readInVx]: File " << file_name
                  << " is empty or malformed. Exiting.";
        messager_.flush("error");
        exit(1);
    }
    while (strcmp(s, "EndOfData") != 0) {
        if (fscanf(input, "%s", sx) != 1 || fscanf(input, "%s", s) != 1) {
            messager_ << "[Glauber::readInVx]: File " << file_name
                      << " is missing its \"EndOfData\" marker. Exiting.";
            messager_.flush("error");
            exit(1);
        }
    }

    for (int i = 0; i <= maxi_num; i++) {
        if (fscanf(input, "%lf", &x) != 1) {
            messager_ << "[Glauber::readInVx]: File " << file_name
                      << " has fewer than " << maxi_num + 1
                      << " entries. Exiting.";
            messager_.flush("error");
            exit(1);
        }
        vx[i] = x;
        if (fscanf(input, "%lf", &x) != 1) {
            messager_ << "[Glauber::readInVx]: File " << file_name
                      << " has fewer than " << maxi_num + 1
                      << " entries. Exiting.";
            messager_.flush("error");
            exit(1);
        }
    }
    fclose(input);

    return vx;

} /* readInVx */

std::vector<double> Glauber::readInVy(
    char *file_name, int maxi_num, int quiet) {
    std::vector<double> vy(maxi_num + 1);
    double y;
    FILE *input;
    char s[120], sy[120];

    if (quiet == 1) {
        messager_ << "[Glauber::readInVy]: Reading in Vy from " << file_name
                  << " ...";
        messager_.flush("info");
    }

    input = fopen(file_name, "r");
    if (input == nullptr) {
        messager_ << "[Glauber::readInVy]: File " << file_name
                  << " not found. Exiting.";
        messager_.flush("error");
        exit(1);
    }
    if (fscanf(input, "%s", s) != 1) {
        messager_ << "[Glauber::readInVy]: File " << file_name
                  << " is empty or malformed. Exiting.";
        messager_.flush("error");
        exit(1);
    }
    while (strcmp(s, "EndOfData") != 0) {
        if (fscanf(input, "%s", sy) != 1 || fscanf(input, "%s", s) != 1) {
            messager_ << "[Glauber::readInVy]: File " << file_name
                      << " is missing its \"EndOfData\" marker. Exiting.";
            messager_.flush("error");
            exit(1);
        }
    }

    for (int i = 0; i <= maxi_num; i++) {
        if (fscanf(input, "%lf", &y) != 1) {
            messager_ << "[Glauber::readInVy]: File " << file_name
                      << " has fewer than " << maxi_num + 1
                      << " entries. Exiting.";
            messager_.flush("error");
            exit(1);
        }
        if (fscanf(input, "%lf", &y) != 1) {
            messager_ << "[Glauber::readInVy]: File " << file_name
                      << " has fewer than " << maxi_num + 1
                      << " entries. Exiting.";
            messager_.flush("error");
            exit(1);
        }
        vy[i] = y;
    }
    fclose(input);

    return vy;
} /* readInVy */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double Glauber::interNuPInSP(double s) {
    double y;
    static int ind = 0;
    static double up, down;
    static int maxi_num;
    static std::vector<double> vx, vy;
    ind++;

    if (glauberData_.projectile.A == 1) return 0.0;

    if (ind == 1) {
        calcRho(&(glauberData_.projectile));
        up = 2.0 * glauberData_.sCutoff;
        down = 0.0;
        maxi_num = glauberData_.interMax;
        vx = makeVx(down, up, maxi_num);
        vy = makeVy(vx.data(), maxi_num);
    } /* if ind */

    if (s > up)
        return 0.0;
    else {
        y = vInterpolate(s, vx.data(), vy.data(), maxi_num);
        if (y < 0.0)
            return 0.0;
        else {
            return y;
        }
    }
} /* interNuPInSP */

double Glauber::interNuTInST(double s) {
    double y;
    static int ind = 0;
    static double up, down;
    static int maxi_num;
    static std::vector<double> vx, vy;

    ind++;
    if (glauberData_.target.A == 1) return 0.0;

    if (ind == 1) {
        calcRho(&(glauberData_.target));

        up = 2.0 * glauberData_.sCutoff;
        down = 0.0;
        maxi_num = glauberData_.interMax;

        vx = makeVx(down, up, maxi_num);
        vy = makeVy(vx.data(), maxi_num);
    } /* if ind */

    if (s > up)
        return 0.0;
    else {
        y = vInterpolate(s, vx.data(), vy.data(), maxi_num);
        if (y < 0.0)
            return 0.0;
        else
            return y;
    }
} /* interNuTInST */

void Glauber::calcRho(Nucleus *nucleus) {
    double f, R_WS;
    /* to pass to AnumIntegrand */

    Nuc_WS_ = nucleus;

    R_WS = nucleus->R_WS;

    if (nucleus->anumFunc == 1)
        f = anum2HO() / (nucleus->rho_WS);
    else if (nucleus->anumFunc == 2)
        f = anum3Gauss(R_WS) / (nucleus->rho_WS);
    else if (nucleus->anumFunc == 3)
        f = anum3Fermi(R_WS) / (nucleus->rho_WS);
    else if (nucleus->anumFunc == 8)
        f = anumHulthen() / (nucleus->rho_WS);
    else
        f = anum3Fermi(R_WS) / (nucleus->rho_WS);

    nucleus->rho_WS = (nucleus->A) / f;
} /* calcRho */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double Glauber::nuInS(double s) {
    double y;

    /* to pass to the densityFunc's */
    NuInS_S_ = s;

    // Nucleus::densityFunc only ever holds 1 (2HO), 2 (3Gauss), 3 (3Fermi) or
    // 8 (Hulthen) -- the same values IntegrandId uses for these cases.
    const IntegrandId id = static_cast<IntegrandId>(Nuc_WS_->densityFunc);

    y = integral(id, 0.0, 1.0);

    return y;
}

double Glauber::anum3Fermi(double R_WS) {
    double up, down, a_WS, rho, f;

    a_WS = Nuc_WS_->a_WS;
    rho = Nuc_WS_->rho_WS;

    /* to pass to Anumintegrand */
    AnumR_ = R_WS / a_WS;

    down = 0.0;
    up = 1.0;

    f = integral(IntegrandId::Anum3FermiInt, down, up);
    f *= 4.0 * M_PI * rho * pow(a_WS, 3.);

    return f;
} /* anum3Fermi */

double Glauber::anum3FermiInt(double xi) {
    double f;
    double r;
    double R_WS, w_WS;

    if (xi == 0.0) xi = TINY;
    if (xi == 1.0) xi = 1.0 - TINY;

    w_WS = Nuc_WS_->w_WS;

    /* already divided by a */
    R_WS = AnumR_;
    r = -log(xi);

    f = r * r;
    f *= (1.0 + w_WS * pow(r / R_WS, 2.));
    f /= (xi + exp(-R_WS));

    return f;
} /* anum3FermiInt */

double Glauber::nuInt3Fermi(double xi) {
    double f;
    double c;
    double z, r, s;
    double w_WS, R_WS, a_WS, rho;

    a_WS = Nuc_WS_->a_WS;
    w_WS = Nuc_WS_->w_WS;
    R_WS = Nuc_WS_->R_WS;
    rho = Nuc_WS_->rho_WS;

    /* xi = exp(-z/a), r = sqrt(z^2 + s^2) */

    if (xi == 0.0) xi = TINY;
    if (xi == 1.0) xi = 1.0 - TINY;

    /* devide by a_WS, make life simpler */

    s = NuInS_S_ / a_WS;
    z = -log(xi);
    r = sqrt(s * s + z * z);
    R_WS /= a_WS;

    c = exp(-R_WS);

    f = 2.0 * a_WS * rho * (glauberData_.sigmaNN);
    f *= 1.0 + w_WS * pow(r / R_WS, 2.);
    f /= xi + c * exp(s * s / (r + z));

    return f;
} /* nuInt3Fermi */

/* %%%%%%% 3 Parameter Gauss %%%%%%%%%%%% */

double Glauber::anum3Gauss(double R_WS) {
    double up, down, a_WS, rho, f;

    a_WS = Nuc_WS_->a_WS;
    rho = Nuc_WS_->rho_WS;

    /* to pass to Anumintegrand */
    AnumR_ = R_WS / a_WS;

    down = 0.0;
    up = 1.0;

    f = integral(IntegrandId::Anum3GaussInt, down, up);
    f *= 4.0 * M_PI * rho * pow(a_WS, 3.);

    return f;
} /* anum3Gauss */

double Glauber::anum3GaussInt(double xi) {
    double y;
    double r_sqr;
    double R_WS, w_WS;

    if (xi == 0.0) xi = TINY;
    if (xi == 1.0) xi = 1.0 - TINY;

    w_WS = Nuc_WS_->w_WS;

    /* already divided by a */

    R_WS = AnumR_;

    r_sqr = -log(xi);

    y = sqrt(r_sqr);

    y *= 1.0 + w_WS * r_sqr / pow(R_WS, 2.);

    /* 2 comes from dr^2 = 2 rdr */

    y /= 2.0 * (xi + exp(-R_WS * R_WS));

    return y;
} /* anum3GaussInt */

double Glauber::nuInt3Gauss(double xi) {
    double f;
    double c;
    double z_sqr, r_sqr, s;
    double w_WS, R_WS, a_WS, rho;

    a_WS = Nuc_WS_->a_WS;
    w_WS = Nuc_WS_->w_WS;
    R_WS = Nuc_WS_->R_WS;
    rho = Nuc_WS_->rho_WS;

    /* xi = exp(-z*z/a/a), r = sqrt(z^2 + s^2) */

    if (xi == 0.0) xi = TINY;
    if (xi == 1.0) xi = 1.0 - TINY;

    /* devide by a_WS, make life simpler */

    s = NuInS_S_ / a_WS;
    z_sqr = -log(xi);
    r_sqr = s * s + z_sqr;
    R_WS /= a_WS;

    c = exp(-R_WS * R_WS);

    f = a_WS * rho * (glauberData_.sigmaNN);
    f *= 1.0 + w_WS * r_sqr / pow(R_WS, 2.);
    f /= sqrt(z_sqr) * (xi + c * exp(s * s));

    return f;
} /* nuInt3Gauss */

/* %%%%%%% 2 Parameter HO %%%%%%%%%%%% */
double Glauber::anum2HO() {
    double up, down, a_WS, rho, f;

    a_WS = Nuc_WS_->a_WS;
    rho = Nuc_WS_->rho_WS;

    down = 0.0;
    up = 1.0;

    f = integral(IntegrandId::Anum2HOInt, down, up);
    f *= 4.0 * M_PI * rho * pow(a_WS, 3.);

    return f;
} /* anum2HO */

double Glauber::anum2HOInt(double xi) {
    double y;
    double r_sqr, r;
    double w_WS;

    if (xi == 0.0) xi = TINY;
    if (xi == 1.0) xi = 1.0 - TINY;

    w_WS = Nuc_WS_->w_WS;

    /* already divided by a */

    r_sqr = -log(xi);

    r = sqrt(r_sqr);

    /* 2 comes from dr^2 = 2 rdr */

    y = r + w_WS * r * r_sqr;
    y /= 2.0;

    return y;
} /* anum2HOInt */

double Glauber::nuInt2HO(double xi) {
    double f;
    double z_sqr, r_sqr, s;
    double w_WS, a_WS, rho;

    a_WS = Nuc_WS_->a_WS;
    w_WS = Nuc_WS_->w_WS;
    rho = Nuc_WS_->rho_WS;

    /* xi = exp(-z*z/a/a), r = sqrt(z^2 + s^2) */

    if (xi == 0.0) xi = TINY;
    if (xi == 1.0) xi = 1.0 - TINY;

    /* devide by a_WS, make life simpler */

    s = NuInS_S_ / a_WS;
    z_sqr = -log(xi);
    r_sqr = s * s + z_sqr;

    /* no need to divide by 2 here because -infty < z < infty and
       we integrate only over positive z */

    if (z_sqr < 0.0) z_sqr = TINY;
    f = a_WS * rho * (glauberData_.sigmaNN);
    f *= (1.0 + w_WS * r_sqr) * exp(-s * s) / sqrt(z_sqr);

    return f;
} /* nuInt2HO */

/* %%%%%%% Hulthen %%%%%%%%%%%% */

double Glauber::anumHulthen() {
    double a_WS, b_WS, rho, f;

    a_WS = Nuc_WS_->a_WS;
    b_WS = Nuc_WS_->w_WS; /* take this to be b */
    rho = Nuc_WS_->rho_WS;

    /* density function has the form */
    /* rho = (1/r^2) (exp(-a r) - exp(-b r))^2 */
    /* int d^3r rho is given by */

    f = 2.0 * (a_WS - b_WS) * (a_WS - b_WS) * M_PI / a_WS / b_WS
        / (a_WS + b_WS);

    f *= rho;

    /* so we return f*rho  calcRho will do f/rho and rho = A/f */
    return f;
} /* anumHulthen */

double Glauber::nuIntHulthen(double xi) {
    double f, g;
    double z, r, s;
    double b_WS, a_WS, rho;

    a_WS = Nuc_WS_->a_WS;
    b_WS = Nuc_WS_->w_WS;
    rho = Nuc_WS_->rho_WS;

    /* xi = exp(-z a), r = sqrt(z^2 + s^2) */

    if (xi == 0.0) xi = TINY;
    if (xi == 1.0) xi = 1.0 - TINY;

    /* multiply by a_WS (in fm^-1), make life simpler */

    s = NuInS_S_ * a_WS;
    z = -log(xi);
    r = sqrt(s * s + z * z);

    /* mult by 2 because the integral is originally over -infty to infty */

    f = 2.0 * a_WS * rho * (glauberData_.sigmaNN);
    g = (1.0 / r) * (exp(-r) - exp(-(b_WS / a_WS) * r));
    // dz = -dxi / (a xi)
    f *= g * g / xi;

    return f;
} /* nuIntHulthen */

double Glauber::evaluateIntegrand(IntegrandId id, double xi) {
    switch (id) {
        case IntegrandId::NuInt2HO:
            return nuInt2HO(xi);
        case IntegrandId::NuInt3Gauss:
            return nuInt3Gauss(xi);
        case IntegrandId::NuInt3Fermi:
            return nuInt3Fermi(xi);
        case IntegrandId::Anum3FermiInt:
            return anum3FermiInt(xi);
        case IntegrandId::Anum3GaussInt:
            return anum3GaussInt(xi);
        case IntegrandId::Anum2HOInt:
            return anum2HOInt(xi);
        case IntegrandId::OLSIntegrand:
            return oLSIntegrand(xi);
        case IntegrandId::NuIntHulthen:
            return nuIntHulthen(xi);
    }
    return 0.0;
}

namespace {
/// The integrand integral() passes to GSL through Glauber::integrandCallback().
struct IntegrandCall {
    /// Glauber instance whose integrand (and its current state, e.g. the
    /// transverse offset of nuInS()) is evaluated.
    Glauber *glauber;
    /// Which integrand to evaluate.
    IntegrandId id;
};
}  // namespace

double Glauber::integrandCallback(double x, void *params) {
    const IntegrandCall *call = static_cast<const IntegrandCall *>(params);
    return call->glauber->evaluateIntegrand(call->id, x);
}

double Glauber::integral(IntegrandId id, double down, double up) {
    if (down == up) return 0.0;
    gsl_integration_cquad_workspace *workspace =
        gsl_integration_cquad_workspace_alloc(QUADRATURE_INTERVALS);
    IntegrandCall call {this, id};
    gsl_function function;
    function.function = &Glauber::integrandCallback;
    function.params = &call;
    double result = 0.;
    double error = 0.;
    // report a failed integral here instead of GSL's default handler,
    // which would abort the whole run
    gsl_error_handler_t *previousHandler = gsl_set_error_handler_off();
    const int status = gsl_integration_cquad(
        &function, down, up, 0., QUADRATURE_TOLERANCE, workspace, &result,
        &error, nullptr);
    gsl_set_error_handler(previousHandler);
    gsl_integration_cquad_workspace_free(workspace);
    if (status != GSL_SUCCESS) {
        messager_ << "[Glauber::integral]: integral over [" << down << ", "
                  << up << "] did not converge (" << gsl_strerror(status)
                  << "); using " << result << " +- " << error << ".";
        messager_.flush("warning");
    }
    return result;
}

double Glauber::oLSIntegrand(double s) {
    double sum, arg, x, r;
    int k, m;
    m = 20;
    sum = 0.0;
    for (k = 1; k <= m; k++) {
        arg = M_PI * (2.0 * k - 1.0) / (2.0 * m);
        x = cos(arg);
        r = sqrt(s * s + b_ * b_ + 2.0 * s * b_ * x);
        sum += interNuTInST(r);
    } /* k */

    return s * sum * M_PI / (1.0 * m) * interNuPInSP(s);

} /* oLSIntegrand */

double Glauber::tAB() {
    double f;
    f = integral(
        IntegrandId::OLSIntegrand, 0.0,
        glauberData_.sCutoff);          // integrate oLSIntegrand(s)
    f *= 2.0 / (glauberData_.sigmaNN);  // here tAB is the number of binary
                                        // collisions, dimensionless (1/fm^4
                                        // integrated over dr_T^2 (gets rid
                                        // of 1/fm^2), divided by sigma (gets
                                        // rid of the other))
    return f;
} /* tAB */

void Glauber::initGlauber(
    double sigmaNN, string target, string projectile, double inb,
    bool setWSDeformParams, double R_WS, double a_WS, double beta2,
    double beta3, double beta4, double gamma, bool forceDminFlag, double d_min,
    double dR_np, double da_np, int imax) {
    string targetName;
    targetName = target;

    string projectileName;
    projectileName = projectile;

    string p_name;
    std::stringstream sp_name;

    const char *EOSPATH = "HYDROPROGRAMPATH";
    char *envPath = getenv(EOSPATH);
    if (envPath != 0 && *envPath != '\0') {
        sp_name << envPath;
        sp_name << "/known_nuclei.in";
        p_name = sp_name.str();
    } else {
        sp_name << "./known_nuclei.in";
        p_name = sp_name.str();
    }

    string paf;
    paf = p_name;

    findNucleusData(
        &(glauberData_.target), targetName, setWSDeformParams, R_WS, a_WS,
        beta2, beta3, beta4, gamma, forceDminFlag, d_min, dR_np, da_np);
    findNucleusData(
        &(glauberData_.projectile), projectileName, setWSDeformParams, R_WS,
        a_WS, beta2, beta3, beta4, gamma, forceDminFlag, d_min, dR_np, da_np);

    glauberData_.sigmaNN = mbToFm2 * sigmaNN;  // sigma in fm^2
    currentA1_ = glauberData_.projectile.A;
    currentA2_ = glauberData_.target.A;
    currentZ1_ = glauberData_.projectile.Z;
    currentZ2_ = glauberData_.target.Z;

    glauberData_.interMax = imax;
    glauberData_.sCutoff = 12.;

    b_ = inb;
}

double Glauber::areaTA(double x, double A) {
    double f;
    f = A * 220. * (1. - exp(-0.025 * x * x));
    return f;
}

ReturnValue Glauber::sampleTARejection(Random *random, NucleusRole nucleus) {
    ReturnValue returnVec;

    double r, x, y, tmp;
    double phi;
    double A = 1.2 * glauberData_.sigmaNN
               / 4.21325504715;  // increase the envelope for larger sigma_inel
                                 // (larger root_s) (was originally written
    // for root(s)=200 GeV, hence the cross section of 4.21325504715 fm^2
    // (=42.13 mb)
    if (nucleus == NucleusRole::Projectile) {
        do {
            phi = 2. * M_PI * random->genrand64_real1();
            r = 6.32456
                * sqrt(-log(
                    (-0.00454545
                     * (-220. * A + areaTA(15., A) * random->genrand64_real1()))
                    / A));
            // here random->genrand64_real1()*areaTA(glauberData_.sCutoff)
            // is a uniform random number on [0, area under f(x)]
            tmp = random->genrand64_real1();

            // x is uniform on [0,1]
            if (r * interNuPInSP(r) > A * r * 11. * exp(-r * r / 40.)) {
                messager_ << std::setprecision(10)
                          << "[Glauber::sampleTARejection]: rejection-"
                             "sampling envelope exceeded: TA="
                          << r * interNuPInSP(r)
                          << ", envelope=" << A * r * 11. * exp(-r * r / 40.);
                messager_.flush("warning");
            }
        } while (tmp > r * interNuPInSP(r) / (A * r * 11. * exp(-r * r / 40.)));
    } else {
        do {
            phi = 2. * M_PI * random->genrand64_real1();
            r = 6.32456
                * sqrt(-log(
                    (-0.00454545
                     * (-220. * A + areaTA(15., A) * random->genrand64_real1()))
                    / A));
            // here random->genrand64_real1()*areaTA(glauberData_.sCutoff)
            // is a uniform random number on [0, area under f(x)]
            tmp = random->genrand64_real1();
            // x is uniform on [0,1]
            if (r * interNuTInST(r) > A * r * 11. * exp(-r * r / 40.)) {
                messager_ << std::setprecision(10)
                          << "[Glauber::sampleTARejection]: rejection-"
                             "sampling envelope exceeded: TA="
                          << r * interNuTInST(r)
                          << ", envelope=" << A * r * 11. * exp(-r * r / 40.);
                messager_.flush("warning");
            }
        } while (tmp > r * interNuTInST(r) / (A * r * 11. * exp(-r * r / 40.)));
    }
    // reject if tmp is larger than the ratio p(y)/f(y),
    // f(y)=A*r*11.*exp(-r*r/40.))
    x = r * cos(phi);
    y = r * sin(phi);
    returnVec.x = x;
    returnVec.y = y;
    returnVec.collided = 0;
    return returnVec;
}

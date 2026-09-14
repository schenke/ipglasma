#include "Glauber.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

#include "Util.h"

using std::string;

namespace {

// U and Xe's beta2 historically defaults to whatever the caller passed in
// (see findNucleusData's beta2 parameter) rather than a fixed literal, so
// their table row uses this sentinel instead of a numeric default.
constexpr double kUseCallerBeta2 = std::numeric_limits<double>::quiet_NaN();

// One row per supported nucleus name. funcCode selects the density
// function (anum2HO/anum3Gauss/anum3Fermi/anumHulthen family) and is
// assigned identically to all three of nucleus->anumFunc/
// anumFuncIntegrand/densityFunc (1=2HO or readFromFile, 2=3Gauss,
// 3=3Fermi, 8=Hulthen).
struct NucleusTemplate {
    const char *name;
    int A, Z;
    double R_WS, w_WS, a_WS;
    double beta2, beta3, beta4, gamma;
    int funcCode;
};

const NucleusTemplate kNucleusTemplates[] = {
    {"Au", 197, 79, 6.37, 0, 0.535, -0.13, 0., -0.03, 0., 3},
    {"Pb", 208, 82, 6.62, 0., 0.546, 0.0, 0., 0.0, 0., 3},
    {"p", 1, 1, 1., 0., 1., 0.0, 0., 0.0, 0., 3},
    {"He3", 3, 2, 0, 0, 0, 0.0, 0., 0.0, 0., 1},
    {"He4", 4, 2, 0, 0, 0, 0.0, 0., 0.0, 0., 1},
    {"d", 2, 1, 1.0, 1.18, 0.228, 0.0, 0., 0.0, 0., 8},
    {"C", 12, 6, 2.44, 1.403, 1.635, 0.0, 0., 0.0, 0., 1},
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
    {"U", 238, 92, 6.81, 0, 0.55, kUseCallerBeta2, 0., 0.093, 0., 3},
    {"Ru", 96, 44, 5.085, 0, 0.46, 0.158, 0., 0.0, 0., 3},
    {"Zr", 96, 40, 5.02, 0, 0.46, 0.0, 0., 0.0, 0., 3},
    {"Xe", 129, 54, 5.42, 0, 0.57, kUseCallerBeta2, 0., -0.003, 0., 3},
};

}  // namespace

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
    nucleus->beta2 = std::isnan(tmpl->beta2) ? beta2 : tmpl->beta2;
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

    if (setWSDeformParams) {
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
    int count;
    int id;

    /* to pass to the densityFunc's */
    NuInS_S_ = s;

    id = Nuc_WS_->densityFunc;

    count = 0;
    y = integral(id, 0.0, 1.0, TOL, &count);

    return y;
}

double Glauber::anum3Fermi(double R_WS) {
    int count = 0;
    double up, down, a_WS, rho, f;

    a_WS = Nuc_WS_->a_WS;
    rho = Nuc_WS_->rho_WS;

    /* to pass to Anumintegrand */
    AnumR_ = R_WS / a_WS;

    down = 0.0;
    up = 1.0;

    f = integral(4, down, up, TOL, &count);
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
    int count = 0;
    double up, down, a_WS, rho, f;

    a_WS = Nuc_WS_->a_WS;
    rho = Nuc_WS_->rho_WS;

    /* to pass to Anumintegrand */
    AnumR_ = R_WS / a_WS;

    down = 0.0;
    up = 1.0;

    f = integral(5, down, up, TOL, &count);
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
    int count = 0;
    double up, down, a_WS, rho, f;

    a_WS = Nuc_WS_->a_WS;
    rho = Nuc_WS_->rho_WS;

    down = 0.0;
    up = 1.0;

    f = integral(6, down, up, TOL, &count);
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
    f *= g * g;

    return f;
} /* nuIntHulthen */

double Glauber::integral(
    int id, double down, double up, double tol, int *count) {
    double dx, y, g1[7];
    int i;

    if (down == up)
        y = 0.0;
    else {
        dx = (up - down) / 6.0;
        for (i = 0; i < 7; i++) {
            if (id == 1)
                g1[i] = nuInt2HO(down + i * dx);
            else if (id == 2)
                g1[i] = nuInt3Gauss(down + i * dx);
            else if (id == 3)
                g1[i] = nuInt3Fermi(down + i * dx);
            else if (id == 4)
                g1[i] = anum3FermiInt(down + i * dx);
            else if (id == 5)
                g1[i] = anum3GaussInt(down + i * dx);
            else if (id == 6)
                g1[i] = anum2HOInt(down + i * dx);
            else if (id == 7)
                g1[i] = oLSIntegrand(down + i * dx);
            else if (id == 8)
                g1[i] = nuIntHulthen(down + i * dx);
        }
        *count = 7;
        y = qnc7(id, tol, down, dx, g1, 0.0, 0.0, count);
    }
    return y;
} /* end of integral */

double Glauber::qnc7(
    int id, double tol, double down, double dx, double *f_of, double pre_sum,
    double area, int *count) {
    int i;
    double left_sum, right_sum, ans;
    static double w[] = {41.0 / 140.0, 54.0 / 35.0, 27.0 / 140.0, 68.0 / 35.0,
                         27.0 / 140,   54.0 / 35.0, 41.0 / 140.0};
    double fl[7];
    double fr[7];
    /*
      qnc7 calculates integral over left and right half of the given interval
      and branches
      to do so, first halve dx
      */

    dx /= 2.0;

    /*
      first calculate the left estimate
      f_of[] contains the evaluated values at down+i, 0< i <7
      store half distanced values for the left sum in fl[]
      */

    if (id == 1) {
        fl[1] = nuInt2HO(down + dx);
        fl[3] = nuInt2HO(down + 3.0 * dx);
        fl[5] = nuInt2HO(down + 5.0 * dx);
    } else if (id == 2) {
        fl[1] = nuInt3Gauss(down + dx);
        fl[3] = nuInt3Gauss(down + 3.0 * dx);
        fl[5] = nuInt3Gauss(down + 5.0 * dx);
    } else if (id == 3) {
        fl[1] = nuInt3Fermi(down + dx);
        fl[3] = nuInt3Fermi(down + 3.0 * dx);
        fl[5] = nuInt3Fermi(down + 5.0 * dx);
    } else if (id == 4) {
        fl[1] = anum3FermiInt(down + dx);
        fl[3] = anum3FermiInt(down + 3.0 * dx);
        fl[5] = anum3FermiInt(down + 5.0 * dx);
    } else if (id == 5) {
        fl[1] = anum3GaussInt(down + dx);
        fl[3] = anum3GaussInt(down + 3.0 * dx);
        fl[5] = anum3GaussInt(down + 5.0 * dx);
    } else if (id == 6) {
        fl[1] = anum2HOInt(down + dx);
        fl[3] = anum2HOInt(down + 3.0 * dx);
        fl[5] = anum2HOInt(down + 5.0 * dx);
    } else if (id == 7) {
        fl[1] = oLSIntegrand(down + dx);
        fl[3] = oLSIntegrand(down + 3.0 * dx);
        fl[5] = oLSIntegrand(down + 5.0 * dx);
    } else if (id == 8) {
        fl[1] = nuIntHulthen(down + dx);
        fl[3] = nuIntHulthen(down + 3.0 * dx);
        fl[5] = nuIntHulthen(down + 5.0 * dx);
    }

    fl[0] = f_of[0];
    fl[2] = f_of[1];
    fl[4] = f_of[2];
    fl[6] = f_of[3];

    *count += 3;

    left_sum = 0.0;
    for (i = 0; i < 7; i++) left_sum += w[i] * fl[i];
    left_sum *= dx;

    /*
      like wise, the right sum is in fr[]
      */

    if (id == 1) {
        fr[1] = nuInt2HO(down + 7.0 * dx);
        fr[3] = nuInt2HO(down + 9.0 * dx);
        fr[5] = nuInt2HO(down + 11.0 * dx);
    } else if (id == 2) {
        fr[1] = nuInt3Gauss(down + 7.0 * dx);
        fr[3] = nuInt3Gauss(down + 9.0 * dx);
        fr[5] = nuInt3Gauss(down + 11.0 * dx);
    } else if (id == 3) {
        fr[1] = nuInt3Fermi(down + 7.0 * dx);
        fr[3] = nuInt3Fermi(down + 9.0 * dx);
        fr[5] = nuInt3Fermi(down + 11.0 * dx);
    } else if (id == 4) {
        fr[1] = anum3FermiInt(down + 7.0 * dx);
        fr[3] = anum3FermiInt(down + 9.0 * dx);
        fr[5] = anum3FermiInt(down + 11.0 * dx);
    } else if (id == 5) {
        fr[1] = anum3GaussInt(down + 7.0 * dx);
        fr[3] = anum3GaussInt(down + 9.0 * dx);
        fr[5] = anum3GaussInt(down + 11.0 * dx);
    } else if (id == 6) {
        fr[1] = anum2HOInt(down + 7.0 * dx);
        fr[3] = anum2HOInt(down + 9.0 * dx);
        fr[5] = anum2HOInt(down + 11.0 * dx);
    } else if (id == 7) {
        fr[1] = oLSIntegrand(down + 7.0 * dx);
        fr[3] = oLSIntegrand(down + 9.0 * dx);
        fr[5] = oLSIntegrand(down + 11.0 * dx);
    } else if (id == 8) {
        fr[1] = nuIntHulthen(down + 7.0 * dx);
        fr[3] = nuIntHulthen(down + 9.0 * dx);
        fr[5] = nuIntHulthen(down + 11.0 * dx);
    }

    fr[0] = f_of[3];
    fr[2] = f_of[4];
    fr[4] = f_of[5];
    fr[6] = f_of[6];

    *count += 3;

    right_sum = 0.0;

    for (i = 0; i < 7; i++) {
        right_sum += w[i] * fr[i];
    }
    right_sum *= dx;

    ans = left_sum + right_sum;

    area += -fabs(pre_sum) + fabs(left_sum) + fabs(right_sum);

    if (fabs(ans - pre_sum) > tol * fabs(area) && (*count < LIMIT)) {
        /*
          branch by calling the function itself
          by calling the qnc7 twice, we are branching
          since left hand side is being calculated first, until the condition
          is satisfied, the left branch keeps branching
          when finally the condition is met by one left-most interval,
          qnc7 returns the right hand side of one up level,
          and the same process resumes
          until the criterion is met by all the branched
          intervals,
          then qnc7 returns to the original right branch and resumes halving
          until the condition is met by all intervals
          (funk, down, dx, f_of[7], pre_ans, ans)
          */

        tol /= 1.414213562;
        left_sum = qnc7(id, tol, down, dx, fl, left_sum, area, count);
        right_sum =
            qnc7(id, tol, down + dx * 6., dx, fr, right_sum, area, count);

        ans = left_sum + right_sum;

    } /* belongs to if*/

    return ans;
} /* end of qnc */

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
    int count = 0;
    f = integral(
        7, 0.0, glauberData_.sCutoff, TOL,
        &count);                        // integrate oLSIntegrand(s)
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

    glauberData_.sigmaNN = 0.1 * sigmaNN;  // sigma in fm^2
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

// Init.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.

#include "Init.h"

#include <algorithm>
#include <utility>
#include <vector>

#include "Instrumentation.h"
#include "NucleonModel.h"
#include "PhysConst.h"
#include "WilsonLineIO.h"

using PhysConst::hbarc;
using PhysConst::mbToFm2;
using PhysConst::Nc2m1;

//**************************************************************************
// Init class.

void Init::sampleTA(Parameters *param, Random *random, Glauber *glauber) {
    Nuclei nuclei = nucleusSampler_.sample(param, random, glauber);
    // move-assign: collisionGeometry_ keeps referring to these vectors
    nucleusA_ = std::move(nuclei.projectile);
    nucleusB_ = std::move(nuclei.target);
}

double Init::computeFluctuatingXG2mu2(
    Parameters *param, double a, double rapidity, double Tp, double qsmuRatio,
    double ySign) {
    // Q_s as a function of \sum T_p and y
    const double exponent = 5.6;  // see 1212.2974 Eq. (17)
    double Qs = 1.;
    double xVal = 0.;
    double g2mu2 = 0.;
    double localrapidity = rapidity;
    double yIn = rapidity;
    double Ydeviation = 10000;
    // iterative loops here to determine the fluctuating Y
    while (std::abs(Ydeviation) > 0.001) {
        if (localrapidity >= 0) {
            Qs = sqrt(qsTable_.qs2(Tp, 0.01 * exp(-localrapidity)));
        } else {
            xVal = Qs * param->colorCharge.xQsFactor / param->collision.sqrtS
                   * exp(ySign * yIn);
            if (xVal == 0)
                Qs = 0.;
            else
                Qs = sqrt(qsTable_.qs2(Tp, 0.01))
                     * sqrt(
                         pow((1 - xVal) / (1 - 0.01), exponent)
                         * pow((0.01 / xVal), 0.2));
        }
        if (Qs == 0) {
            Ydeviation = 0;
            g2mu2 = 0.;
        } else {
            g2mu2 = Qs * Qs / qsmuRatio / qsmuRatio * a * a / hbarc / hbarc
                    / param->coupling.g
                    / param->coupling.g;  // lattice units? check

            Ydeviation = localrapidity
                         - log(
                             0.01
                             / (Qs * param->colorCharge.xQsFactor
                                / param->collision.sqrtS * exp(ySign * yIn)));
            localrapidity =
                log(0.01
                    / (Qs * param->colorCharge.xQsFactor
                       / param->collision.sqrtS * exp(ySign * yIn)));
        }
    }
    if (g2mu2 != g2mu2) {
        g2mu2 = 0.;
    }
    return g2mu2;
}

void Init::computeCellColorCharge(
    // set g^2\mu^2 as the sum of the individual nucleons' g^2\mu^2, using
    // Q_s(b,y) prop to g^mu(b,y) Also compute N_part using Glauber
    Lattice *lat, Parameters *param, int ipos, double a, double rapidity) {
    if (param->colorCharge.useFluctuatingX) {  // Local Qs dependent x
        lat->cells[ipos]->setg2mu2A(computeFluctuatingXG2mu2(
            param, a, rapidity, lat->cells[ipos]->getTpA(),
            param->colorCharge.QsMuRatio, 1.));
        lat->cells[ipos]->setg2mu2B(computeFluctuatingXG2mu2(
            param, a, rapidity, lat->cells[ipos]->getTpB(),
            param->colorCharge.QsMuRatio, -1.));
    } else {  // Fixed x: jimwlkInitialX with JIMWLK, else projectileX/targetX
        // nucleus A
        lat->cells[ipos]->setg2mu2A(
            qsTable_.qs2(
                lat->cells[ipos]->getTpA(),
                param->initialX(NucleusRole::Projectile))
            / param->colorCharge.QsMuRatio / param->colorCharge.QsMuRatio * a
            * a / hbarc / hbarc / param->coupling.g
            / param->coupling.g);  // lattice units? check

        // nucleus B
        lat->cells[ipos]->setg2mu2B(
            qsTable_.qs2(
                lat->cells[ipos]->getTpB(),
                param->initialX(NucleusRole::Target))
            / param->colorCharge.QsMuRatio / param->colorCharge.QsMuRatio * a
            * a / hbarc / hbarc / param->coupling.g / param->coupling.g);
    }
}

void Init::readQsTable(const std::string &fileName) { qsTable_.read(fileName); }

void Init::setColorChargeDensity(
    Lattice *lat, Parameters *param, Random *random, Glauber *glauber) {
    IPG_PROFILE_SCOPE("initialization.color_charge_density");
    messager_.info(
        "[Init::setColorChargeDensity]: set color charge density ...");

    const int N = param->lattice.size;
    const double a = param->lattice.L / N;  // lattice spacing in fm

    double nucleiInAverage =
        static_cast<double>(param->collision.nucleiToAverage);

    if (!param->collision.useNucleus) {
        setConstantColorChargeDensity(lat, param);
        return;
    }

#pragma omp parallel for
    for (int ipos = 0; ipos < N * N; ipos++) {
        lat->cells[ipos]->setg2mu2A(0.);
        lat->cells[ipos]->setg2mu2B(0.);
    }

    // the rapidity the fluctuating-x iteration starts from
    double rapidity = 0.;
    if (param->colorCharge.useFluctuatingX) {
        rapidity = computeEffectiveRapidity(param);
    } else {
        messager_ << "[Init::setColorChargeDensity]: Q_s^2 of the projectile "
                     "at x = "
                  << param->initialX(NucleusRole::Projectile)
                  << ", of the target at x = "
                  << param->initialX(NucleusRole::Target);
        messager_.flush("info");
    }

    sampleNucleonAnisotropyAngles(param, random);
    sampleNucleonProfiles(param, random);

    // test what a smooth Woods-Saxon would give
    if (param->nucleus.useSmoothNucleus) {
        computeSmoothNucleusThickness(lat, param, glauber);
    } else {
        // Non-smooth nucleus add all T_p's (new in version 1.2)
        computeThicknessFromNucleons(lat, param, nucleiInAverage);
    }

// get Q_s^2 (and from that g^2mu^2) for a given \sum T_p and Y
#pragma omp parallel for
    for (int ipos = 0; ipos < N * N; ipos++) {
        computeCellColorCharge(lat, param, ipos, a, rapidity);
    }
    messager_.info(
        "[Init::setColorChargeDensity]: Color charge densities for nucleus A "
        "and B set. ");
}

double Init::computeEffectiveRapidity(Parameters *param) {
    const double input = param->colorCharge.rapidity;
    if (!param->colorCharge.usePseudoRapidity) return input;
    // when using pseudorapidity as input convert to rapidity here.
    // later include Jacobian in multiplicity and energy
    messager_ << "[Init::setColorChargeDensity]: Using pseudorapidity "
              << input;
    messager_.flush("info");
    double m = param->colorCharge.jacobianMass;  // in GeV
    double P =
        0.13 + 0.32 * pow(param->collision.sqrtS / 1000., 0.115);  // in GeV
    const double rapidity =
        0.5
        * log(
            (sqrt(pow(cosh(input), 2.) + m * m / (P * P)) + sinh(input))
            / (sqrt(pow(cosh(input), 2.) + m * m / (P * P)) - sinh(input)));
    messager_ << "[Init::setColorChargeDensity]: Corresponds to rapidity "
              << rapidity;
    messager_.flush("info");
    return rapidity;
}

void Init::setConstantColorChargeDensity(Lattice *lat, Parameters *param) {
    const int N = param->lattice.size;
    const double g2mu2 = param->collision.g2mu * param->collision.g2mu
                         / param->coupling.g / param->coupling.g;
    for (int pos = 0; pos < N * N; pos++) {
        lat->cells[pos]->setg2mu2A(g2mu2);
        lat->cells[pos]->setg2mu2B(g2mu2);
    }
    param->event.success = 1;
    messager_.info(
        "[Init::setColorChargeDensity]: constant color charge density set");
}

void Init::sampleNucleonAnisotropyAngles(Parameters *param, Random *random) {
    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();
    double xi = param->subnucleon.protonAnisotropy;
    if (xi != 0.) {
        for (int i = 0; i < A1; i++) {
            nucleusA_.at(i).phi = 2 * M_PI * random->genrand64_real2();
        }

        for (int i = 0; i < A2; i++) {
            nucleusB_.at(i).phi = 2 * M_PI * random->genrand64_real2();
        }
    } else {
        for (int i = 0; i < A1; i++) {
            nucleusA_.at(i).phi = 0.;
        }

        for (int i = 0; i < A2; i++) {
            nucleusB_.at(i).phi = 0.;
        }
    }
}

void Init::computeSmoothNucleusThickness(
    Lattice *lat, Parameters *param, Glauber *glauber) {
    messager_ << "[Init::setColorChargeDensity]: Using smooth nucleus for "
                 "test purposes. "
                 "Does not "
                 "include "
                 "deformation.";
    messager_.flush("info");
    const int N = param->lattice.size;
    const double L = param->lattice.L;
    const double a = L / N;
    // Both profiles are centered at the origin: the impact parameter is
    // sampled later and applied by shiftFieldsWithImpactParameter(), the
    // same as for nucleon-based nuclei. Nucleus A is the projectile.
    double normA = 0.;
    double normB = 0.;
    for (int ix = 0; ix < N; ix++) {
        const double x = -L / 2. + a * ix;
        for (int iy = 0; iy < N; iy++) {
            const double y = -L / 2. + a * iy;
            const int localpos = lat->positionFromXY(ix, iy);
            const double r = sqrt(x * x + y * y);

            double TA = glauber->interNuPInSP(r);
            double TB = glauber->interNuTInST(r);
            normA += TA * a * a;
            normB += TB * a * a;

            // remove the far tails of each nucleus separately (the overlap
            // region is not known before b is sampled)
            if (TA < 0.001) TA = 0.;
            if (TB < 0.001) TB = 0.;
            lat->cells[localpos]->setTpA(TA);
            lat->cells[localpos]->setTpB(TB);
        }
    }

    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            const int localpos = lat->positionFromXY(ix, iy);
            lat->cells[localpos]->setTpA(
                lat->cells[localpos]->getTpA() / normA * glauber->nucleusA1()
                * hbarc * hbarc);
            lat->cells[localpos]->setTpB(
                lat->cells[localpos]->getTpB() / normB * glauber->nucleusA2()
                * hbarc * hbarc);
        }
    }
}

void Init::sampleNucleonProfiles(Parameters *param, Random *random) {
    const std::unique_ptr<NucleonModel> model = NucleonModel::create(*param);
    profilesA_.clear();
    profilesB_.clear();
    for (const ReturnValue &nucleon : nucleusA_) {
        profilesA_.push_back(model->sample(*random, nucleon));
    }
    for (const ReturnValue &nucleon : nucleusB_) {
        profilesB_.push_back(model->sample(*random, nucleon));
    }
}

double Init::computeNucleonThicknessAtCell(
    const std::vector<std::unique_ptr<NucleonProfile>> &profiles, double x,
    double y, double nucleiInAverage) const {
    double Tp = 0.;
    for (const auto &profile : profiles) {
        Tp += profile->thickness(x, y) / nucleiInAverage;  // add up all T_p
    }
    return Tp;
}

void Init::computeThicknessFromNucleons(
    Lattice *lat, Parameters *param, double nucleiInAverage) {
    const int N = param->lattice.size;
    const double L = param->lattice.L;
    const double a = L / N;

#pragma omp parallel for
    for (int ipos = 0; ipos < N * N; ipos++) {
        // loop over all positions
        int iy = lat->yFromPosition(ipos);
        int ix = lat->xFromPosition(ipos);
        double x = -L / 2. + a * ix;
        double y = -L / 2. + a * iy;

        lat->cells[ipos]->setTpA(
            computeNucleonThicknessAtCell(profilesA_, x, y, nucleiInAverage));
        lat->cells[ipos]->setTpB(
            computeNucleonThicknessAtCell(profilesB_, x, y, nucleiInAverage));
    }
}

std::vector<double> Init::computeWilsonLineMomentumKernel(
    int N, int sites, double m, double UVdamp) {
    std::vector<double> momentumKernel(static_cast<std::size_t>(sites));
#pragma omp parallel for
    for (int pos = 0; pos < sites; ++pos) {
        const int i = latticeX(pos, N);
        const int j = latticeY(pos, N);
        const double kx =
            2. * M_PI
            * (-0.5 + static_cast<double>(i) / static_cast<double>(N));
        const double ky =
            2. * M_PI
            * (-0.5 + static_cast<double>(j) / static_cast<double>(N));
        const double sx = sin(kx / 2.);
        const double sy = sin(ky / 2.);
        const double kt2 = 4. * (sx * sx + sy * sy);

        if (m == 0.) {
            momentumKernel[static_cast<std::size_t>(pos)] =
                (kt2 != 0.) ? 1. / kt2 : 0.;
        } else {
            momentumKernel[static_cast<std::size_t>(pos)] =
                (1. / (kt2 + m * m)) * exp(-sqrt(kt2) * UVdamp);
        }
    }
    return momentumKernel;
}

ColorChargeScales Init::computeWilsonLineColorChargeScales(
    Lattice *lat, int sites, double g, double invNy) {
    ColorChargeScales scales;
    std::vector<double> &colorChargeScaleA = scales.projectile;
    std::vector<double> &colorChargeScaleB = scales.target;
    colorChargeScaleA.resize(static_cast<std::size_t>(sites));
    colorChargeScaleB.resize(static_cast<std::size_t>(sites));
#pragma omp parallel for
    for (int pos = 0; pos < sites; ++pos) {
        colorChargeScaleA[static_cast<std::size_t>(pos)] =
            g * sqrt(lat->cells[pos]->getg2mu2A() * invNy);
        colorChargeScaleB[static_cast<std::size_t>(pos)] =
            g * sqrt(lat->cells[pos]->getg2mu2B() * invNy);
    }
    return scales;
}

void Init::setV(Lattice *lat, Parameters *param, Random *random) {
    IPG_PROFILE_SCOPE("initialization.wilson_lines");
    messager_.info("[Init::setV]: Setting Wilson lines ...");
    const int A1 = nucleusA_.size();
    const int A2 = nucleusB_.size();

    const double d2 = param->collision.sigmaNN * mbToFm2 / M_PI;  // in fm^2
    const int N = param->lattice.size;
    const int Ny = param->colorCharge.Ny;
    const int sites = N * N;
    const int nn[2] = {N, N};
    const double L = param->lattice.L;
    const double a = L / N;  // lattice spacing in fm
    const double m = param->subnucleon.m * a / hbarc;
    const double g = param->coupling.g;
    const double invNy = 1. / static_cast<double>(Ny);
    double UVdamp = param->subnucleon.UVDamp;  // GeV^-1
    UVdamp = UVdamp / a * hbarc;

    // The lattice Poisson/UV kernel depends only on transverse momentum and
    // run parameters.  The historical implementation recomputed the same
    // sin/sqrt/exp expressions for every longitudinal sheet of both nuclei.
    std::vector<double> momentumKernel =
        computeWilsonLineMomentumKernel(N, sites, m, UVdamp);

    // rhoACoeffData owns the Nc2m1*sites backing storage; rhoACoeff is a
    // pointer-per-component view over it for FFT::fftnComplexArray's T**
    // interface.
    std::vector<complex<double>> rhoACoeffData(
        static_cast<std::size_t>(Nc2m1) * sites);
    std::vector<complex<double> *> rhoACoeff(Nc2m1);
    for (int i = 0; i < Nc2m1; i++) {
        rhoACoeff[i] = rhoACoeffData.data() + i * sites;
    }

    auto applyMomentumKernel = [&]() {
#pragma omp parallel for
        for (int n = 0; n < Nc2m1; ++n) {
            complex<double> *rho = rhoACoeff[n];
            for (int pos = 0; pos < sites; ++pos) {
                rho[pos] *= momentumKernel[static_cast<std::size_t>(pos)];
            }
        }
    };

    // Reuse the bulk Gaussian buffers for every longitudinal sheet.  The
    // linear ordering matches the historical pos-major/color-minor gauss()
    // call sequence exactly.
    std::vector<double> gaussianField(
        static_cast<std::size_t>(sites) * static_cast<std::size_t>(Nc2m1));
    std::vector<double> gaussianScratch;
    gaussianScratch.reserve(5 * ((gaussianField.size() + 1) / 2));

    // g2mu2 and Ny are fixed throughout Wilson-line construction.  Cache the
    // color-independent site scale once for each nucleus instead of repeating
    // the same sqrt in every longitudinal sheet.
    const ColorChargeScales scales =
        computeWilsonLineColorChargeScales(lat, sites, g, invNy);
    const std::vector<double> &colorChargeScaleA = scales.projectile;
    const std::vector<double> &colorChargeScaleB = scales.target;

    auto fillColorCharge = [&](const std::vector<double> &scale) {
        {
            IPG_PROFILE_SCOPE("initialization.wilson_random.gauss");
            random->gaussBulk(
                gaussianField.data(), gaussianField.size(), gaussianScratch);
        }
        {
            IPG_PROFILE_SCOPE("initialization.wilson_random.scale");
#pragma omp parallel for
            for (int pos = 0; pos < sites; ++pos) {
                const double localScale = scale[static_cast<std::size_t>(pos)];
                const std::size_t base = static_cast<std::size_t>(pos)
                                         * static_cast<std::size_t>(Nc2m1);
                for (int n = 0; n < Nc2m1; ++n) {
                    rhoACoeff[n][pos] =
                        localScale
                        * gaussianField[base + static_cast<std::size_t>(n)];
                }
            }
        }
    };

    // loop over longitudinal direction, once for each nucleus
    auto evolveNucleusWilsonLine =
        [&](const std::vector<double> &colorChargeScale,
            std::vector<Matrix> &U) {
            for (int k = 0; k < Ny; k++) {
                {
                    IPG_PROFILE_SCOPE("initialization.wilson_random");
                    fillColorCharge(colorChargeScale);
                }

                fft_.fftnComplexArray(
                    rhoACoeff.data(), rhoACoeff.data(), nn, 1, Nc2m1);

                {
                    IPG_PROFILE_SCOPE("initialization.wilson_Poisson");
                    applyMomentumKernel();
                }

                fft_.fftnComplexArray(
                    rhoACoeff.data(), rhoACoeff.data(), nn, -1, Nc2m1);

                {
                    IPG_PROFILE_SCOPE("initialization.wilson_exponent");
#pragma omp parallel
                    {
                        std::vector<double> in(Nc2m1, 0.);
                        Matrix temp(1.);
                        Matrix tempNew(0.);

#pragma omp for
                        for (int pos = 0; pos < sites; pos++) {
                            for (int aa = 0; aa < Nc2m1; aa++) {
                                // expmCoeff calculates exp(i in[a] t[a]), so
                                // multiply by -1 (not -i).
                                in[aa] = -(rhoACoeff[aa][pos]).real();
                            }
                            tempNew = Matrix::fromAlgebraExponent(in);
                            temp = tempNew * U[pos];
                            U[pos] = temp;
                        }
                    }
                }
            }
        };

    evolveNucleusWilsonLine(colorChargeScaleA, lat->U);
    evolveNucleusWilsonLine(colorChargeScaleB, lat->U2);

    if (param->output.writeWilsonLineSnapshot) {
        WilsonLineIO().writeTrainingData(lat, param);
    }

    // output U
    if (param->wilsonLines.writeWilsonLines > 0
        && (param->jimwlk.saveSnapshots || !param->jimwlk.enabled)) {
        WilsonLineIO io;
        io.write(
            lat, param, NucleusRole::Projectile,
            param->initialX(NucleusRole::Projectile));
        io.write(
            lat, param, NucleusRole::Target,
            param->initialX(NucleusRole::Target));
    }
    writeGeometry(lat, param);

    messager_ << "[Init::setV]: Wilson lines V_A and V_B set on rank "
              << param->run.MPIRank << ". ";
    messager_.flush("info");
}

void Init::sampleImpactParameter(Parameters *param) {
    collisionGeometry_.sampleImpactParameter(param, random_ptr_);
}

void Init::computeCollisionGeometryQuantities(Lattice *lat, Parameters *param) {
    collisionGeometry_.computeQuantities(lat, param, random_ptr_);
}

void Init::init(
    Lattice *lat, Parameters *param, Random *random, Glauber *glauber,
    InitializationMethod init_method) {
    random_ptr_ = random;

    messager_.info("[Init::init]: Initializing fields ... ");

    if (!param->collision.useNucleus) {
        // No real collision geometry in the constant-g^2mu case: skip
        // sampleImpactParameter() entirely, but still give b/phi_RP the
        // same defaults it would have produced.
        param->event.b = 0.;
        param->event.phiRP = 0.;
        param->event.success = 1;
    } else {
        readQsTable(param->colorCharge.nucleusQsTableFileName);
    }

    // The configuration files are only used with nucleonPositionsFromFile;
    // don't require them to exist for Woods-Saxon sampling.
    if (param->collision.useNucleus && param->nucleus.nucleonPositionsFromFile
        && init_method == InitializationMethod::SampleColorCharges) {
        nucleusSampler_.readConfigurations(param, glauber, random);
    }

    if (init_method == InitializationMethod::ReadWlineBinary
        or init_method == InitializationMethod::ReadWlineText) {
        // to read Wilson lines from file, with their nuclei's geometry so
        // that main's collision-geometry loop treats them like sampled ones
        WilsonLineIO().read(
            lat, param,
            (init_method == InitializationMethod::ReadWlineBinary) ? 2 : 1);
        // (the geometry files are not written again: a run that writes the
        // Wilson lines it read, e.g. after JIMWLK, uses the same names)
        if (param->collision.useNucleus) readGeometry(lat, param);
    } else {
        // to generate your own Wilson lines
        if (param->collision.useNucleus) {
            // populate the lists nucleusA_ and nucleusB_ with position data
            sampleTA(param, random, glauber);
        }
        setColorChargeDensity(lat, param, random, glauber);
        // sample color charges and find Wilson lines V_A and V_B
        setV(lat, param, random);
    }
}

void Init::writeGeometry(Lattice *lat, Parameters *param) {
    // with every event whose Wilson lines are written, so that they can be
    // read back and collided like sampled ones
    if (param->wilsonLines.writeWilsonLines == 0
        || !param->wilsonLines.writeGeometry || !param->collision.useNucleus) {
        return;
    }
    WilsonLineIO io;
    io.writeGeometry(lat, param, NucleusRole::Projectile, nucleusA_);
    io.writeGeometry(lat, param, NucleusRole::Target, nucleusB_);
}

void Init::readGeometry(Lattice *lat, Parameters *param) {
    WilsonLineIO io;
    NucleusGeometry projectile =
        io.readGeometry(lat, param, NucleusRole::Projectile);
    NucleusGeometry target = io.readGeometry(lat, param, NucleusRole::Target);
    // the color-charge maps were built with this ratio, which a posterior
    // parameter set may have chosen for the event
    if (!PhysConst::isClose(projectile.QsMuRatio, target.QsMuRatio)) {
        messager_ << "[Init::readGeometry]: the geometry files of the "
                     "projectile and the target were written with different "
                     "QsMuRatio ("
                  << projectile.QsMuRatio << ", " << target.QsMuRatio
                  << "). Exiting.";
        messager_.flush("error");
        exit(1);
    }
    param->colorCharge.QsMuRatio = projectile.QsMuRatio;
    // move-assign: collisionGeometry_ keeps referring to these vectors
    nucleusA_ = std::move(projectile.nucleons);
    nucleusB_ = std::move(target.nucleons);
}

void Init::shiftFieldsWithImpactParameter(Lattice *lat, Parameters *param) {
    messager_.info(
        "[Init::shiftFieldsWithImpactParameter]: Shifting fields with impact "
        "parameter...");
    const double b = param->event.b;
    const double phiRP = param->event.phiRP;
    messager_ << "[Init::shiftFieldsWithImpactParameter]: b = " << b
              << " fm, phi_RP = " << phiRP;
    messager_.flush("info");

    const int N = param->lattice.size;
    BufferLattice lat_tmp(param->lattice.size);
    for (int ipos = 0; ipos < N * N; ipos++) {
        lat_tmp.buffer1[ipos] = lat->U[ipos];
        lat_tmp.buffer2[ipos] = lat->U2[ipos];
    }

    const double L = param->lattice.L;
    const double a = L / N;  // lattice spacing in fm
    for (int ipos = 0; ipos < N * N; ipos++) {
        int ix = lat->xFromPosition(ipos);
        int iy = lat->yFromPosition(ipos);
        double x = -L / 2. + a * ix;
        double y = -L / 2. + a * iy;

        double xA = x - b / 2. * cos(phiRP);
        double yA = y - b / 2. * sin(phiRP);
        double xB = x + b / 2. * cos(phiRP);
        double yB = y + b / 2. * sin(phiRP);

        int ixA = static_cast<int>((xA + L / 2.) / a);
        int iyA = static_cast<int>((yA + L / 2.) / a);
        int ixB = static_cast<int>((xB + L / 2.) / a);
        int iyB = static_cast<int>((yB + L / 2.) / a);

        if (ixA < 0 || ixA >= N || iyA < 0 || iyA >= N) {
            lat->U[ipos] = one_;
        } else {
            int posA = lat->positionFromXY(ixA, iyA);
            lat->U[ipos] = lat_tmp.buffer1[posA];
        }
        if (ixB < 0 || ixB >= N || iyB < 0 || iyB >= N) {
            lat->U2[ipos] = one_;
        } else {
            int posB = lat->positionFromXY(ixB, iyB);
            lat->U2[ipos] = lat_tmp.buffer2[posB];
        }
    }
}

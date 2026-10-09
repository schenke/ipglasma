// Copyright (C) 2023 Chun Shen

#include "Parameters.h"

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "Glauber.h"
#include "NucleusSampler.h"
#include "PhysConst.h"

std::string Parameters::loadPosteriorParameterSetsFromFile(
    std::string posteriorFileName, std::vector<std::vector<float>> &ParamSet) {
    std::ifstream posteriorFile(posteriorFileName.c_str());
    if (!posteriorFile.is_open()) {
        return "cannot open the posterior parameter file " + posteriorFileName;
    }
    std::string tempLine;
    std::getline(posteriorFile, tempLine);
    int lineNumber = 1;
    while (std::getline(posteriorFile, tempLine)) {
        lineNumber++;
        std::stringstream lineStream(tempLine);
        std::string cell;
        std::vector<float> parsedRow;
        while (std::getline(lineStream, cell, ',')) {
            try {
                parsedRow.push_back(std::stof(cell));
            } catch (const std::exception &) {
                return posteriorFileName + ":" + std::to_string(lineNumber)
                       + ": " + cell + " is not a number";
            }
        }
        ParamSet.push_back(parsedRow);
    }
    posteriorFile.close();
    return "";
}

std::string Parameters::loadPosteriorParameterSets(const int itype) {
    std::string fileName;
    std::vector<std::vector<float>> *table = nullptr;
    std::size_t columns = 0;
    if (itype == 1) {
        fileName = "tables/posterior.csv";
        table = &posteriorParamSets_;
        columns = 7;  // m, BG, BGq, smearingWidth, Nq, QsMuRatio, dqMin
    } else if (itype == 2 || itype == 4) {
        fileName = (itype == 2) ? "tables/posterior_Nq3.csv"
                                : "tables/posterior5020_Nq3.csv";
        table = &posteriorParamSetsNq3_;
        columns = 6;  // m, BG, BGq, smearingWidth, QsMuRatio, dqMin
    } else {
        return "";
    }
    const std::string problem =
        loadPosteriorParameterSetsFromFile(fileName, *table);
    if (!problem.empty()) return problem;
    if (table->empty()) {
        return fileName + " contains no parameter sets";
    }
    for (std::size_t row = 0; row < table->size(); row++) {
        if ((*table)[row].size() < columns) {
            return fileName + ": parameter set " + std::to_string(row) + " has "
                   + std::to_string((*table)[row].size()) + " values, expected "
                   + std::to_string(columns);
        }
    }
    return "";
}

void Parameters::setParamsWithPosteriorParameterSet(const int itype, int iset) {
    if (itype == 1) {
        // variant Nq
        iset = (iset % posteriorParamSets_.size());
        event.subNucleonParamSet = iset;
        messager_ << "[Parameters::setParamsWithPosteriorParameterSet]: "
                     "Using subnucleon parameter set "
                  << iset << " (variable Nq).";
        messager_.flush("info");
        subnucleon.m = posteriorParamSets_[iset][0];
        subnucleon.BG = posteriorParamSets_[iset][1];
        subnucleon.BGq = posteriorParamSets_[iset][2];
        subnucleon.smearingWidth = posteriorParamSets_[iset][3];
        subnucleon.Nq = posteriorParamSets_[iset][4];
        colorCharge.QsMuRatio = posteriorParamSets_[iset][5];
        subnucleon.dqMin = posteriorParamSets_[iset][6];
    } else if (itype == 2 || itype == 4) {
        // fixed Nq = 3
        iset = (iset % posteriorParamSetsNq3_.size());
        event.subNucleonParamSet = iset;
        messager_ << "[Parameters::setParamsWithPosteriorParameterSet]: "
                     "Using subnucleon parameter set "
                  << iset << " (Nq = 3).";
        messager_.flush("info");
        subnucleon.m = posteriorParamSetsNq3_[iset][0];
        subnucleon.BG = posteriorParamSetsNq3_[iset][1];
        subnucleon.BGq = posteriorParamSetsNq3_[iset][2];
        subnucleon.smearingWidth = posteriorParamSetsNq3_[iset][3];
        subnucleon.Nq = 3.;
        colorCharge.QsMuRatio = posteriorParamSetsNq3_[iset][4];
        subnucleon.dqMin = posteriorParamSetsNq3_[iset][5];
    }
}

double Parameters::initialX(NucleusRole role) const {
    if (jimwlk.enabled) return jimwlk.initialX;
    // the initial condition does not correspond to a fixed x
    if (colorCharge.useFluctuatingX) return -1.;
    return (role == NucleusRole::Projectile) ? colorCharge.projectileX
                                             : colorCharge.targetX;
}

double Parameters::xBeforeJimwlk(NucleusRole role) const {
    if (wilsonLines.readInitialWilsonLines != 0 && wilsonLines.readX > 0.) {
        return wilsonLines.readX;
    }
    return initialX(role);
}

std::vector<std::string> Parameters::validationErrors() const {
    // Checks of a single value are part of the input parameter table (see
    // ParameterTable.cpp) and run while reading; these combine several.
    std::vector<std::string> errors;
    auto fail = [&errors](const std::ostringstream &message) {
        errors.push_back(message.str());
    };

    if (jimwlk.saveSnapshots && wilsonLines.writeWilsonLines == 0) {
        std::ostringstream message;
        message << "jimwlkSaveSnapshots = 1 requires writing Wilson lines "
                   "(writeWilsonLines = 1 or 2)";
        fail(message);
    }

    if (wilsonLines.writeWilsonLines != 0) {
        const std::filesystem::path outputPath(wilsonLines.wilsonLinePath);
        if (!std::filesystem::is_directory(outputPath)) {
            std::ostringstream message;
            message << "wilsonLinePath " << outputPath.string()
                    << " is not an existing directory";
            fail(message);
        }
    }

    if (coupling.runningCoupling && coupling.LambdaQCD >= coupling.mu0) {
        std::ostringstream message;
        message << "LambdaQCD (" << coupling.LambdaQCD
                << ") must be smaller than mu0 (" << coupling.mu0
                << ") with running coupling; otherwise alpha_s is singular "
                   "or negative where the local scale is zero";
        fail(message);
    }

    // JIMWLK running coupling (jimwlkAlphaS 0), independent of
    // runningCoupling
    if (jimwlk.enabled && PhysConst::isClose(jimwlk.alphaS, 0.)
        && jimwlk.LambdaQCD >= jimwlk.mu0) {
        std::ostringstream message;
        message << "jimwlkLambdaQCD (" << jimwlk.LambdaQCD
                << ") must be smaller than jimwlkMu0 (" << jimwlk.mu0
                << ") with the JIMWLK running coupling; otherwise alpha_s "
                   "is singular or negative at large dipole sizes";
        fail(message);
    }

    // <Qs> is only computed from the collision geometry of nuclei (sampled
    // or read with their Wilson lines); without it, 1/<Qs> and the running
    // coupling at the event-averaged Qs (used at least for the hydro
    // output) are undefined
    if (!collision.useNucleus && evolution.inverseQsForMaxTime) {
        std::ostringstream message;
        message << "inverseQsForMaxTime = 1 needs the event-averaged Qs, which "
                   "is not computed with useNucleus = 0";
        fail(message);
    }
    if (!collision.useNucleus && coupling.runningCoupling) {
        std::ostringstream message;
        message << "runningCoupling = 1 needs the event-averaged Qs, which is "
                   "not computed with useNucleus = 0";
        fail(message);
    }

    // JIMWLK evolves to smaller x: from readWilsonLinesX for read Wilson
    // lines, otherwise from jimwlkInitialX
    if (jimwlk.enabled
        && (colorCharge.projectileX > xBeforeJimwlk(NucleusRole::Projectile)
            || colorCharge.targetX > xBeforeJimwlk(NucleusRole::Target))) {
        const bool readLines =
            wilsonLines.readInitialWilsonLines != 0 && wilsonLines.readX > 0.;
        std::ostringstream message;
        message << "projectileX (" << colorCharge.projectileX
                << ") and targetX (" << colorCharge.targetX
                << ") must not be larger than "
                << (readLines ? "readWilsonLinesX (" : "jimwlkInitialX (")
                << xBeforeJimwlk(NucleusRole::Projectile)
                << "): JIMWLK evolves "
                << (readLines ? "the read Wilson lines"
                              : "the initial condition")
                << " from there to smaller x";
        fail(message);
    }

    // with a fixed x, Q_s^2 is looked up at y = ln(0.01/x), and the nuclear
    // Q_s table starts at y = 0
    const bool samplesColorCharges =
        collision.useNucleus && wilsonLines.readInitialWilsonLines == 0;
    // nucleons sampled from the built-in density profile need one (He3 and
    // He4 have none); a polarization already set nucleonPositionsFromFile.
    // A smooth nucleus takes its thickness from the profile even when the
    // nucleons are read from configuration files.
    if (samplesColorCharges
        && (!nucleus.nucleonPositionsFromFile || nucleus.useSmoothNucleus)
        && !nucleus.useInputWSParams) {
        const bool projectileOk =
            speciesHasDensityProfile(collision.projectile);
        const bool targetOk = speciesHasDensityProfile(collision.target);
        if (!projectileOk || !targetOk) {
            std::ostringstream message;
            if (!projectileOk && !targetOk
                && collision.projectile == collision.target) {
                message << "projectile and target " << collision.projectile;
            } else if (!projectileOk && !targetOk) {
                message << "projectile " << collision.projectile
                        << " and target " << collision.target;
            } else {
                message << (projectileOk ? "target " : "projectile ")
                        << (projectileOk ? collision.target
                                         : collision.projectile);
            }
            message << " ha" << (projectileOk || targetOk ? "s" : "ve");
            if (nucleus.useSmoothNucleus) {
                message << " no built-in density profile for the smooth "
                           "nucleus (useSmoothNucleus 1): use "
                           "useInputWSParams 1";
            } else {
                message << " no built-in density profile to sample the "
                           "nucleons from: use nucleonPositionsFromFile 1 "
                           "(configuration files) or useInputWSParams 1";
            }
            fail(message);
        }
    }

    if (collision.useNucleus && collision.bMin > collision.bMax) {
        std::ostringstream message;
        message << "bMin (" << collision.bMin
                << ") must not be larger than bMax (" << collision.bMax << ")";
        fail(message);
    }

    // lightNucleusOption must select a configuration file of both nuclei
    // when they are read (a polarization already set
    // nucleonPositionsFromFile)
    if (samplesColorCharges && nucleus.nucleonPositionsFromFile) {
        const std::string projectileError = lightNucleusOptionError(
            speciesMassNumber(collision.projectile),
            nucleus.lightNucleusOption);
        const std::string targetError = lightNucleusOptionError(
            speciesMassNumber(collision.target), nucleus.lightNucleusOption);
        std::ostringstream message;
        if (!projectileError.empty() && projectileError == targetError) {
            message << "projectile and target: " << projectileError;
            fail(message);
        } else {
            if (!projectileError.empty()) {
                message << "projectile: " << projectileError;
                fail(message);
                message.str("");
            }
            if (!targetError.empty()) {
                message << "target: " << targetError;
                fail(message);
            }
        }
    }

    if (samplesColorCharges && jimwlk.enabled && jimwlk.initialX > 0.01) {
        std::ostringstream message;
        message << "jimwlkInitialX (" << jimwlk.initialX
                << ") must not be larger than 0.01; the nuclear Q_s table "
                   "starts at x = 0.01";
        fail(message);
    }
    if (samplesColorCharges && !jimwlk.enabled && !colorCharge.useFluctuatingX
        && (colorCharge.projectileX > 0.01 || colorCharge.targetX > 0.01)) {
        std::ostringstream message;
        message << "projectileX (" << colorCharge.projectileX
                << ") and targetX (" << colorCharge.targetX
                << ") must not be larger than 0.01 with useFluctuatingX = 0; "
                   "the nuclear Q_s table starts at x = 0.01";
        fail(message);
    }

    // JIMWLK and fluctuating x are mutually exclusive:
    // JIMWLK evolves the initial condition to the fixed x of the nuclei,
    // while fluctuating x computes a local x from the local Q_s
    if (jimwlk.enabled && colorCharge.useFluctuatingX) {
        std::ostringstream message;
        message << "useJIMWLK = 1 and useFluctuatingX = 1 are mutually "
                   "exclusive: JIMWLK evolves the initial condition to the "
                   "fixed x of the nuclei, while fluctuating x computes a "
                   "local x from the local Q_s";
        fail(message);
    }

    // the final time is always written; an output time at or after it
    // would not be reached
    if (!evolution.inverseQsForMaxTime) {
        for (const double t : output.outputTimes) {
            if (t >= evolution.maxTime) {
                std::ostringstream message;
                message << "outputTimes value " << t
                        << " must be smaller than maxTime ("
                        << evolution.maxTime
                        << "); the final time is always written";
                fail(message);
            }
        }
    }

    // averaging over nuclei is not supported for protons
    if (collision.nucleiToAverage > 1
        && (collision.projectile == "p" || collision.target == "p")) {
        std::ostringstream message;
        message << "nucleiToAverage = " << collision.nucleiToAverage
                << " (averaging over nuclei) is not supported for collisions "
                   "with a proton";
        fail(message);
    }

    // the posterior parameter sets are fits of hot-spot nucleons
    if (subnucleon.subNucleonParamType != 0
        && subnucleon.nucleonModel != "hotspots") {
        std::ostringstream message;
        message << "subNucleonParamType = " << subnucleon.subNucleonParamType
                << " (a posterior parameter set) requires nucleonModel "
                   "hotspots, not "
                << subnucleon.nucleonModel;
        fail(message);
    }

    return errors;
}

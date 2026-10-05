// Copyright (C) 2023 Chun Shen

#include "Parameters.h"

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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
        columns = 7;  // m, BG, BGq, smearingWidth, NqBase, QsMuRatio, dqMin
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
        subnucleon.NqBase = posteriorParamSets_[iset][4];
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
        subnucleon.NqBase = 3.;
        colorCharge.QsMuRatio = posteriorParamSetsNq3_[iset][4];
        subnucleon.dqMin = posteriorParamSetsNq3_[iset][5];
    }
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

    // <Qs> is only computed from the collision geometry of a sampled
    // nucleus; without it, 1/<Qs> and the running coupling at the
    // event-averaged Qs (used at least for the hydro output) are undefined
    const bool noAverageQs =
        !collision.useNucleus || wilsonLines.readInitialWilsonLines != 0;
    if (noAverageQs && evolution.inverseQsForMaxTime) {
        std::ostringstream message;
        message << "inverseQsForMaxTime = 1 needs the event-averaged Qs, which "
                   "is not computed with useNucleus = 0 or "
                   "readInitialWilsonLines = 1 or 2";
        fail(message);
    }
    if (noAverageQs && coupling.runningCoupling) {
        std::ostringstream message;
        message << "runningCoupling = 1 needs the event-averaged Qs, which is "
                   "not computed with useNucleus = 0 or "
                   "readInitialWilsonLines = 1 or 2";
        fail(message);
    }

    // with a fixed x, Q_s^2 is looked up at the input rapidities, and the
    // nuclear Q_s table starts at y = 0 (a pseudorapidity keeps its sign)
    const bool samplesColorCharges =
        collision.useNucleus && wilsonLines.readInitialWilsonLines == 0;
    if (samplesColorCharges && !colorCharge.useFluctuatingX
        && (colorCharge.rapidityA < 0. || colorCharge.rapidityB < 0.)) {
        std::ostringstream message;
        message << "rapidityA (" << colorCharge.rapidityA << ") and rapidityB ("
                << colorCharge.rapidityB
                << ") must not be negative with useFluctuatingX = 0; the "
                   "nuclear Q_s table starts at y = 0";
        fail(message);
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

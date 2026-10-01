// Copyright (C) 2023 Chun Shen

#include "Parameters.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void Parameters::loadPosteriorParameterSetsFromFile(
    std::string posteriorFileName, std::vector<std::vector<float>> &ParamSet) {
    std::ifstream posteriorFile(posteriorFileName.c_str());
    if (!posteriorFile.is_open()) {
        messager_ << "[Parameters::loadPosteriorParameterSetsFromFile]: "
                     "Cannot open posterior file: "
                  << posteriorFileName;
        messager_.flush("error");
        exit(1);
    }
    std::string tempLine;
    std::getline(posteriorFile, tempLine);
    while (std::getline(posteriorFile, tempLine)) {
        std::stringstream lineStream(tempLine);
        std::string cell;
        std::vector<float> parsedRow;
        while (std::getline(lineStream, cell, ',')) {
            parsedRow.push_back(std::stof(cell));
        }
        ParamSet.push_back(parsedRow);
    }
    posteriorFile.close();
}

void Parameters::loadPosteriorParameterSets(const int itype) {
    if (itype == 1) {
        loadPosteriorParameterSetsFromFile(
            "tables/posterior.csv", posteriorParamSets_);
    } else if (itype == 2) {
        loadPosteriorParameterSetsFromFile(
            "tables/posterior_Nq3.csv", posteriorParamSetsNq3_);
    } else if (itype == 4) {
        loadPosteriorParameterSetsFromFile(
            "tables/posterior5020_Nq3.csv", posteriorParamSetsNq3_);
    }
}

void Parameters::setParamsWithPosteriorParameterSet(const int itype, int iset) {
    if (itype == 1) {
        // variant Nq
        iset = (iset % posteriorParamSets_.size());
        messager_ << "[Parameters::setParamsWithPosteriorParameterSet]: "
                     "Using subnucleon parameter set "
                  << iset << " (variable Nq).";
        messager_.flush("info");
        setm(posteriorParamSets_[iset][0]);
        setBG(posteriorParamSets_[iset][1]);
        setBGq(posteriorParamSets_[iset][2]);
        setSmearingWidth(posteriorParamSets_[iset][3]);
        setNqBase(posteriorParamSets_[iset][4]);
        setQsmuRatio(posteriorParamSets_[iset][5]);
        setDqmin(posteriorParamSets_[iset][6]);
    } else if (itype == 2 || itype == 4) {
        // fixed Nq = 3
        iset = (iset % posteriorParamSetsNq3_.size());
        messager_ << "[Parameters::setParamsWithPosteriorParameterSet]: "
                     "Using subnucleon parameter set "
                  << iset << " (Nq = 3).";
        messager_.flush("info");
        setm(posteriorParamSetsNq3_[iset][0]);
        setBG(posteriorParamSetsNq3_[iset][1]);
        setBGq(posteriorParamSetsNq3_[iset][2]);
        setSmearingWidth(posteriorParamSetsNq3_[iset][3]);
        setNqBase(3.);
        setQsmuRatio(posteriorParamSetsNq3_[iset][4]);
        setDqmin(posteriorParamSetsNq3_[iset][5]);
    }
}

std::vector<std::string> Parameters::validationErrors() const {
    // Checks of a single value are part of the input parameter table (see
    // ParameterTable.cpp) and run while reading; these combine several.
    std::vector<std::string> errors;
    auto fail = [&errors](const std::ostringstream &message) {
        errors.push_back(message.str());
    };

    if (getSaveSnapshots() && getWriteWilsonLines() == 0) {
        std::ostringstream message;
        message << "saveSnapshots = 1 requires writing Wilson lines "
                   "(writeWilsonLines = 1 or 2)";
        fail(message);
    }

    if (getWriteWilsonLines() != 0) {
        const std::filesystem::path outputPath(getWilsonLinePath());
        if (!std::filesystem::is_directory(outputPath)) {
            std::ostringstream message;
            message << "wilsonLinePath " << outputPath.string()
                    << " is not an existing directory";
            fail(message);
        }
    }

    if (getRunningCoupling() && getLambdaQCD() >= getMuZero()) {
        std::ostringstream message;
        message << "LambdaQCD (" << getLambdaQCD()
                << ") must be smaller than muZero (" << getMuZero()
                << ") with running coupling; otherwise alpha_s is singular "
                   "or negative where the local scale is zero";
        fail(message);
    }

    // JIMWLK running coupling (alphas_jimwlk 0), independent of
    // runningCoupling
    if (getUseJIMWLK() && getJimwlk_alphas() <= 1e-10
        && getLambdaQCD_jimwlk() >= getMu0_jimwlk()) {
        std::ostringstream message;
        message << "Lambda_QCD_jimwlk (" << getLambdaQCD_jimwlk()
                << ") must be smaller than mu0_jimwlk (" << getMu0_jimwlk()
                << ") with the JIMWLK running coupling; otherwise alpha_s "
                   "is singular or negative at large dipole sizes";
        fail(message);
    }

    return errors;
}

bool Parameters::ValidParameters() {
    const std::vector<std::string> errors = validationErrors();
    for (const std::string &error : errors) {
        messager_ << "[Parameters::ValidParameters]: " << error << ".";
        messager_.flush("error");
    }
    return errors.empty();
}

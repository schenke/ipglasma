// Copyright (C) 2023 Chun Shen

#include "Parameters.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Lattice.h"
#include "PhysConst.h"

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

bool Parameters::ValidParameters() {
    // Check if the parameters are valid. Return true if they are, false
    // otherwise. This function can be used to validate the parameters before
    // running the simulation.
    if (size_ <= 0) {
        messager_ << "[Parameters::ValidParameters]: Invalid lattice size "
                  << size_ << ".";
        messager_.flush("error");
        return false;
    }
    if (getWriteWilsonLines() != 0
        and !Lattice::IsValidWilsonLineDataFormat(getWriteWilsonLines())) {
        messager_ << "[Parameters::ValidParameters]: Invalid Wilson line "
                     "data format "
                  << getWriteWilsonLines();
        messager_.flush("error");
        return false;
    }
    if (getSaveSnapshots() and getWriteWilsonLines() == 0) {
        messager_ << "[Parameters::ValidParameters]: Cannot save snapshots "
                     "(saveSnapshots = "
                  << getSaveSnapshots() << ") "
                  << "without writing Wilson lines (writeWilsonLines = "
                  << getWriteWilsonLines() << ").";
        messager_.flush("error");
        return false;
    }

    if (getWriteWilsonLines() != 0) {
        const std::filesystem::path outputPath(getWilsonLinePath());
        if (!std::filesystem::exists(outputPath)
            || !std::filesystem::is_directory(outputPath)) {
            messager_ << "[Parameters::ValidParameters]: Wilson line output "
                         "directory does not exist: "
                      << outputPath.string() << ".";
            messager_.flush("error");
            return false;
        }
    }

    if (getRunningCoupling()) {
        if (getLambdaQCD() >= getMuZero()) {
            messager_ << "[Parameters::ValidParameters]: LambdaQCD ("
                      << getLambdaQCD() << ") must be smaller than muZero ("
                      << getMuZero()
                      << "); otherwise the running-coupling formula's log "
                         "argument is non-positive at the lattice edges "
                         "(where the local scale is zero), making alpha_s "
                         "singular or negative.";
            messager_.flush("error");
            return false;
        }
        if (11. * PhysConst::Nc - 2. * getNFlavors() <= 0.) {
            messager_ << "[Parameters::ValidParameters]: nFlavors ("
                      << getNFlavors()
                      << ") is too large -- the one-loop beta-function "
                         "coefficient (11*Nc - 2*nFlavors) must be "
                         "positive.";
            messager_.flush("error");
            return false;
        }
    }

    if (getUseJIMWLK() && getc_jimwlk() <= 0.) {
        messager_ << "[Parameters::ValidParameters]: c_jimwlk ("
                  << getc_jimwlk() << ") must be positive.";
        messager_.flush("error");
        return false;
    }

    return true;
}
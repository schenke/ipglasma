// Copyright (C) 2023 Chun Shen

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include "Parameters.h"

void Parameters::loadPosteriorParameterSetsFromFile(
    std::string posteriorFileName, std::vector<std::vector<float>> &ParamSet) {
    std::ifstream posteriorFile(posteriorFileName.c_str());
    if (!posteriorFile.is_open()) {
        std::cout << "[Paramters] Can not open posterior file: "
                  << posteriorFileName << std::endl;
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
        loadPosteriorParameterSetsFromFile("tables/posterior.csv",
                                           posteriorParamSets_);
    } else if (itype == 2) {
        loadPosteriorParameterSetsFromFile("tables/posterior_Nq3.csv",
                                           posteriorParamSetsNq3_);
    } else if (itype == 4) {
        loadPosteriorParameterSetsFromFile("tables/posterior5020_Nq3.csv",
                                           posteriorParamSetsNq3_);
    }
}


void Parameters::setParamsWithPosteriorParameterSet(const int itype,
                                                    int iset) {
    if (itype == 1) {
        // variant Nq
        iset = (iset % posteriorParamSets_.size());
        std::cout << "set subnucleon param set:" << iset << std::endl;
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
        std::cout << "set subnucleon param set (Nq = 3):" << iset << std::endl;
        setm(posteriorParamSetsNq3_[iset][0]);
        setBG(posteriorParamSetsNq3_[iset][1]);
        setBGq(posteriorParamSetsNq3_[iset][2]);
        setSmearingWidth(posteriorParamSetsNq3_[iset][3]);
        setNqBase(3.);
        setQsmuRatio(posteriorParamSetsNq3_[iset][4]);
        setDqmin(posteriorParamSetsNq3_[iset][5]);
    }
}

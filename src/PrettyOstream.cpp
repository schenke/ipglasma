// Copyright Chun Shen @ 2017
// This class is inspired by the JetScapeLogger class written by Joern Putschke

#include "PrettyOstream.h"

#include <sys/resource.h>
#include <sys/time.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::cout;
using std::endl;
using std::string;

PrettyOstream::PrettyOstream() {}

PrettyOstream::~PrettyOstream() {}

//! This function flushes out message to the screen
void PrettyOstream::flush(string type) {
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    if (type == "info") {
        info(message_stream.str());
    } else if (type == "warning") {
        warning(message_stream.str());
    } else if (type == "error") {
        error(message_stream.str());
    } else if (type == "debug") {
        debug(message_stream.str());
    }
    message_stream.str("");
    message_stream.clear();
}

//! This function output information message
void PrettyOstream::info(string message) {
    cout << "[Info] " << getMemoryUsage() << " " << message << endl;
}

//! This function output debug message
void PrettyOstream::debug(string message) {
    cout << CYAN << "[Debug] " << getMemoryUsage() << " " << message << RESET
         << endl;
}

//! This function output warning message
void PrettyOstream::warning(string message) {
    cout << BOLD << YELLOW << "[Warning] " << message << RESET << endl;
}

//! This function output error message
void PrettyOstream::error(string message) {
    cout << BOLD << RED << "[Error] " << message << RESET << endl;
}

//! This function returns a string for the memory usage
//! of the current program in MB
string PrettyOstream::getMemoryUsage() {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        double memory_usage_in_MB = 0.0;
#ifdef APPLE
        memory_usage_in_MB = usage.ru_maxrss / 1024. / 1024.;  // MB in Apple
#else
        memory_usage_in_MB = usage.ru_maxrss / 1024.;  // MB in linux
#endif
        std::ostringstream memory_usage;
        memory_usage << std::setprecision(4) << memory_usage_in_MB << " MB";
        return (memory_usage.str());
    } else {
        return (0);
    }
}

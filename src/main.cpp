#include <stdio.h>

#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifndef DISABLEMPI
#include "mpi.h"
#endif

#include "Evolution.h"
#include "ForwardLightCone.h"
#include "GluonMultiplicity.h"
#include "Init.h"
#include "InputFile.h"
#include "Instrumentation.h"
#include "JIMWLK.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "PrettyOstream.h"
#include "Random.h"
#include "WilsonLineIO.h"

#define _SECURE_SCL 0
#define _HAS_ITERATOR_DEBUGGING 0

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;

bool readInput(Parameters *param, int argc, char *argv[], int rank);
void display_logo();
void writeparams(Parameters *param);

// main program 1
int main(int argc, char *argv[]) {
    int rank;
    int size;

    int nev = 1;
    if (argc == 3) {
        nev = atoi(argv[2]);
    }

#ifndef DISABLEMPI
    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // get number of processes
#else
    rank = 0;
    size = 1;
#endif

    ipg::Profiler::instance().initialize(rank);

    int h5Flag = 0;
    PrettyOstream messager;

    Parameters paramStorage;
    Parameters *param = &paramStorage;
    param->run.MPIRank = rank;
    param->run.MPISize = size;

    // read and validate the parameters from the input file
    if (!readInput(param, argc, argv, rank)) {
#ifndef DISABLEMPI
        MPI_Finalize();
#endif
        return 1;
    }

    // initialize random generator using time and seed from input file
    Random randomStorage;
    Random *random = &randomStorage;
    unsigned long long int rnum;
    if (!param->random.useSeedList) {
        if (param->random.useTimeForSeed) {
            std::random_device ran_dev;
            rnum = ran_dev();
        } else {
            rnum = param->random.seed;
            messager << "[main::main]: Random seed = " << rnum + (rank * 1000)
                     << " - entered directly +rank*1000.";
            messager.flush("info");
        }
        param->run.randomSeed = rnum + rank * 1000;
        if (param->random.useTimeForSeed) {
            messager << "[main::main]: Random seed = " << param->run.randomSeed;
            messager.flush("info");
        }
        random->init_genrand64(rnum + rank * 1000);
        random->gslRandomInit(rnum + rank * 1000);
    } else {
        ifstream fin;
        fin.open("seedList");
        std::vector<unsigned long long int> seedList(size, 0);
        if (fin) {
            for (int i = 0; i < size; i++) {
                if (!(fin >> seedList[i])) {
                    messager.error(
                        "[main::main]: Not enough random seeds (or a "
                        "malformed entry) for the number of processors "
                        "selected. Exiting.");
                    exit(1);
                }
            }
        } else {
            messager.error(
                "[main::main]: Random seed file 'seedList' not found. "
                "Exiting.");
            exit(1);
        }
        fin.close();
        param->run.randomSeed = seedList[rank];
        random->init_genrand64(seedList[rank]);
        random->gslRandomInit(seedList[rank]);
        messager << "[main::main]: Random seed on rank " << rank << " = "
                 << seedList[rank] << " read from list.";
        messager.flush("info");
    }
    // only the hot-spot positions are sampled from the gamma distribution
    if (param->subnucleon.nucleonModel == "hotspots") {
        random->setGammaIncCDF(param->subnucleon.omega);
    }

    // event loop starts ...
    for (int iev = 0; iev < nev; iev++) {
        const int profiler_event_id = rank + iev * size;
        ipg::Profiler::instance().beginEvent(profiler_event_id);

        messager << "[main::main]: Generating event " << iev + 1 << " out of "
                 << nev << " ...";
        messager.flush("info");
        // welcome
        if (rank == 0) display_logo();

        if (param->subnucleon.subNucleonParamType > 0) {
            IPG_PROFILE_SCOPE("initialization.subnucleon_parameters");
            // sample the sub-nucleon parameters from the posterior distribution
            int iSubNucleonParamSet = param->subnucleon.subNucleonParamSet;
            if (iSubNucleonParamSet == -1) {
                iSubNucleonParamSet =
                    static_cast<int>(random->genrand64_int63() % 2147483647ULL);
            }
            param->setParamsWithPosteriorParameterSet(
                param->subnucleon.subNucleonParamType, iSubNucleonParamSet);
        }

        // initialize helper class objects

        param->event.eventId = rank + iev * size;
        param->event.success = 0;

        {
            IPG_PROFILE_SCOPE("parameters.write");
            writeparams(param);
        }

        int nn[2];
        nn[0] = param->lattice.size;
        nn[1] = param->lattice.size;

        stringstream strup_name;
        strup_name << "usedParameters" << param->event.eventId << ".dat";
        string up_name;
        up_name = strup_name.str();
        ofstream fout1(up_name.c_str(), std::ios::app);
        fout1 << "# Random seed used on rank " << rank << ": "
              << param->run.randomSeed << endl;
        fout1.close();

        // initialize init object
        Init init(nn);

        // initialize group
        Group group;

        // initialize Glauber class
        messager << "[main::main]: Init Glauber on rank " << param->run.MPIRank
                 << " ... ";
        messager.flush("info");
        Glauber glauber;
        {
            IPG_PROFILE_SCOPE("glauber.initialize");
            glauber.initGlauber(
                param->collision.sigmaNN, param->collision.target,
                param->collision.projectile, param->event.b,
                param->nucleus.useInputWSParams, param->nucleus.radiusWS,
                param->nucleus.diffusenessWS, param->nucleus.beta2,
                param->nucleus.beta3, param->nucleus.beta4,
                param->nucleus.gamma, param->nucleus.forceDMin,
                param->nucleus.dMin, param->nucleus.deltaRnp,
                param->nucleus.deltaAnp, 100);
        }

        // initialize evolution object
        Evolution evolution(nn);

        // either read k_T spectrum from file or do a fresh start
        if (param->output.readMultFromFile) {
            GluonMultiplicity::readNkt(param);
        }

        // Keep the lattice lifetime inside this block so destruction is timed
        // before the per-event profile is written.
        {
            // allocate lattice
            Lattice lat(param, param->lattice.size);
            messager.info("[main::main]: Lattice generated.");

            param->event.success = 0;

            // initialize U-fields on the lattice
            InitializationMethod init_method;
            if (param->wilsonLines.readInitialWilsonLines == 0) {
                init_method = InitializationMethod::SampleColorCharges;
            } else {
                init_method = (param->wilsonLines.readInitialWilsonLines == 1)
                                  ? InitializationMethod::ReadWlineText
                                  : InitializationMethod::ReadWlineBinary;
            }
            // First generate the V
            init.init(&lat, param, random, &glauber, init_method);

            if (param->jimwlk.enabled) {
                messager.info("[main::main]: Start JIMWLK");
                JIMWLK jimwlkSolver(*param, &group, &lat, random);
                jimwlkSolver.evolution();
                messager.info("[main::main]: Finish JIMWLK");

                // Store final Wilson lines after JIMWLK evolution
                if (param->wilsonLines.writeWilsonLines > 0) {
                    WilsonLineIO io;
                    io.write(
                        &lat, param, NucleusRole::Projectile,
                        param->jimwlk.xProjectile);
                    io.write(
                        &lat, param, NucleusRole::Target,
                        param->jimwlk.xTarget);
                }
            }

            if (param->evolution.mode == 1) {
                while (param->event.success == 0) {
                    // sample collision impact parameter
                    // and compute Npart, Ncoll,etc, and check if there was a
                    // collision
                    init.sampleImpactParameter(param);
                    init.computeCollisionGeometryQuantities(&lat, param);
                }
                init.shiftFieldsWithImpactParameter(&lat, param);
                ForwardLightCone(&group).initialize(&lat, param);
                messager.info("[main::main]: Start CYM evolution");
                // do the CYM evolution of the initialized fields using
                // parmeters in param
                evolution.run(&lat, &group, param);
            }

#ifndef DISABLEMPI
            {
                IPG_PROFILE_SCOPE("mpi.barrier");
                MPI_Barrier(MPI_COMM_WORLD);
            }
#endif

            messager.info("[main::main]: One event finished");
            if (param->output.writeOutputsToHDF5) {
                IPG_PROFILE_SCOPE("output.hdf5_collect_event");
                int status = 0;
                stringstream h5output_filename;
                h5output_filename << "RESULTS_rank" << rank;
                stringstream collect_command;
                collect_command
                    << "python3 utilities/combine_events_into_hdf5.py ."
                    << " --output_filename " << h5output_filename.str()
                    << " --event_id " << param->event.eventId;
                status = system(collect_command.str().c_str());
                if (status == 0) {
                    messager << "[main::main]: Collected this event's "
                                "output into an HDF5 file.";
                    messager.flush("info");
                } else {
                    messager << "[main::main]: combine_events_into_hdf5.py "
                                "exited with status "
                             << status
                             << " while collecting this event's "
                                "output.";
                    messager.flush("warning");
                }
                h5Flag = 1;
            }

            {
                IPG_PROFILE_SCOPE("correctness.fingerprint");
                ipg::writeLatticeFingerprint(&lat, rank, param->event.eventId);
            }
        }  // lattice lifetime

        ipg::Profiler::instance().endEvent();
    }

#ifndef DISABLEMPI
    // Every rank must have finished appending its last event to its
    // RESULTS_rank*.h5 before rank 0 merges (and deletes) those files.
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (h5Flag == 1 && rank == 0) {
        int status = 0;
        stringstream collect_command;
        collect_command << "python3 utilities/combine_events_into_hdf5.py ."
                        << " --output_filename RESULTS"
                        << " --combine_hdf5_files_only";
        status = system(collect_command.str().c_str());
        if (status == 0) {
            messager << "[main::main]: Combined all per-rank HDF5 files "
                        "into RESULTS.h5.";
            messager.flush("info");
        } else {
            messager << "[main::main]: combine_events_into_hdf5.py exited "
                        "with status "
                     << status
                     << " while combining the per-rank HDF5 "
                        "files.";
            messager.flush("warning");
        }
    }

#ifndef DISABLEMPI
    MPI_Finalize();
#endif

    return 0;
}

void display_logo() {
    cout << endl;
    cout << "--------------------------------------------------------------"
            "----"
            "--"
            "---------"
         << endl;
    cout << "| Classical Yang-Mills evolution with IP-Glasma initial "
            "configurations      |"
         << endl;
    cout << "--------------------------------------------------------------"
            "----"
            "--"
            "---------"
         << endl;
    cout << "| References:                                                 "
            "    "
            "  "
            "        |"
         << endl;
    cout << "| B. Schenke, P. Tribedy, R. Venugopalan                      "
            "    "
            "  "
            "        |"
         << endl;
    cout << "| Phys. Rev. Lett. 108, 252301 (2012) and Phys. Rev. C86, "
            "034908 "
            "(2012)     |"
         << endl;
    cout << "| H. Mäntysaari, B. Schenke, C. Shen and W. Zhao "
            "           "
            "        "
            "        |"
         << endl;
    cout << "| Phys. Rev. Lett. 135, 022302 (2025)                             "
            "          |"
         << endl;
    cout << "--------------------------------------------------------------"
            "----"
            "--"
            "---------"
         << endl;

    cout << "This version uses Qs as obtained from IP-Sat using the sum "
            "over "
            "proton T_p(b)"
         << endl;
    cout << "This is a simple MPI version that runs many events in one "
            "job. No "
            "communication."
         << endl;

    cout << "Run using large lattices to improve convergence of the root "
            "finder "
            "in initial condition. "
         << "Recommended: 600x600 using L=30fm" << endl;
    cout << endl;
}

bool readInput(Parameters *param, int argc, char *argv[], int rank) {
    // the first given argument is taken to be the input file name
    // if none is given, that file name is "input"
    PrettyOstream messager;
    const string file_name = (argc > 1) ? argv[1] : "input";
    if (rank == 0) {
        messager << "[main::readInput]: Reading parameters from \"" << file_name
                 << "\".";
        messager.flush("info");
    }

    const InputFile input(file_name);
    std::vector<string> errors = param->readInput(input);
    // checks combining several parameters need all values read
    if (errors.empty()) errors = param->validationErrors();
    if (!errors.empty()) {
        if (rank == 0) {
            for (const string &error : errors) {
                messager << "[main::readInput]: " << error;
                messager.flush("error");
            }
            messager << "[main::readInput]: Invalid input parameters. "
                        "Exiting.";
            messager.flush("error");
        }
        return false;
    }
    return true;
}

void writeparams(Parameters *param) {
    // write the values of all input parameters this event used to
    // "usedParameters<id>.dat", in input-file syntax
    stringstream strup_name;
    strup_name << "usedParameters" << param->event.eventId << ".dat";
    ofstream fout1(strup_name.str());
    time_t rawtime = time(0);
    fout1 << "# Input parameters used by IP-Glasma for event "
          << param->event.eventId << ", written " << ctime(&rawtime);
    fout1 << "# This file is a valid input file. Running it does not "
             "reproduce this event:\n"
             "# the random numbers also depend on the MPI rank and the "
             "event's position in\n"
             "# the run (and on the time with useTimeForSeed 1), and "
             "subNucleonParamSet -1\n"
             "# draws a new posterior parameter set.\n";
    param->writeInputParameters(fout1);
    if (param->subnucleon.subNucleonParamType > 0) {
        // these values are not input parameters in this case, so they are
        // not in the list above
        fout1 << std::setprecision(9) << "# Posterior parameter set used: "
              << param->event.subNucleonParamSet << "\n#   m "
              << param->subnucleon.m << ", BG " << param->subnucleon.BG
              << ", BGq " << param->subnucleon.BGq << ", smearingWidth "
              << param->subnucleon.smearingWidth << ", NqBase "
              << param->subnucleon.NqBase << ", QsMuRatio "
              << param->colorCharge.QsMuRatio << ", dqMin "
              << param->subnucleon.dqMin << "\n";
    }
}

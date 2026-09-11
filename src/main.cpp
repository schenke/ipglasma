#include <stdio.h>

#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifndef DISABLEMPI
#include "mpi.h"
#endif

#include "Evolution.h"
#include "FFT.h"
#include "Init.h"
#include "Instrumentation.h"
#include "JIMWLK.h"
#include "Lattice.h"
#include "Matrix.h"
#include "Parameters.h"
#include "PrettyOstream.h"
#include "Random.h"
#include "Setup.h"

#define _SECURE_SCL 0
#define _HAS_ITERATOR_DEBUGGING 0

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;

int readInput(
    Setup *setup, Parameters *param, int argc, char *argv[], int rank);
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

    Parameters *param = new Parameters();
    param->setMPIRank(rank);
    param->setMPISize(size);
    Setup setup;

    // read parameters from file
    readInput(&setup, param, argc, argv, rank);

    // Validate parameters before proceeding
    if (!param->ValidParameters()) {
        messager << "[main::main]: Invalid parameters detected. Exiting.";
        messager.flush("error");
        return 1;
    }

    // initialize random generator using time and seed from input file
    Random *random = new Random();
    unsigned long long int rnum;
    if (param->getUseSeedList() == 0) {
        if (param->getUseTimeForSeed() == 1) {
            std::random_device ran_dev;
            rnum = ran_dev();
            // rnum = time(0) + param->getSeed() * 10000;
        } else {
            rnum = param->getSeed();
            messager << "[main::main]: Random seed = " << rnum + (rank * 1000)
                     << " - entered directly +rank*1000.";
            messager.flush("info");
        }
        param->setRandomSeed(rnum + rank * 1000);
        if (param->getUseTimeForSeed() == 1) {
            messager << "[main::main]: Random seed = "
                     << param->getRandomSeed();
            //<< " made from time " << rnum - param->getSeed() - (rank * 1000)
            //<< " and argument (+1000*rank) "
            //<< param->getSeed() + (rank * 1000);
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
                if (!fin.eof()) {
                    fin >> seedList[i];
                } else {
                    messager.error(
                        "[main::main]: Not enough random seeds for the number "
                        "of processors selected. Exiting.");
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
        param->setRandomSeed(seedList[rank]);
        random->init_genrand64(seedList[rank]);
        random->gslRandomInit(seedList[rank]);
        messager << "[main::main]: Random seed on rank " << rank << " = "
                 << seedList[rank] << " read from list.";
        messager.flush("info");
    }
    random->setGammaIncCDF(param->getOmega());

    // event loop starts ...
    for (int iev = 0; iev < nev; iev++) {
        const int profiler_event_id = rank + iev * size;
        ipg::Profiler::instance().beginEvent(profiler_event_id);

        messager << "[main::main]: Generating event " << iev + 1 << " out of "
                 << nev << " ...";
        messager.flush("info");
        // welcome
        if (rank == 0) display_logo();

        if (param->getSubNucleonParamType() > 0) {
            IPG_PROFILE_SCOPE("initialization.subnucleon_parameters");
            // sample the sub-nucleon parameters from the posterior distribution
            int iSubNucleonParamSet = param->getSubNucleonParamSet();
            if (iSubNucleonParamSet == -1) {
                iSubNucleonParamSet = random->genrand64_int63();
            }
            param->setParamsWithPosteriorParameterSet(
                param->getSubNucleonParamType(), iSubNucleonParamSet);
        }

        // initialize helper class objects

        param->setEventId(rank + iev * size);
        param->setSuccess(0);

        {
            IPG_PROFILE_SCOPE("parameters.write");
            writeparams(param);
        }

        int nn[2];
        nn[0] = param->getSize();
        nn[1] = param->getSize();

        stringstream strup_name;
        strup_name << "usedParameters" << param->getEventId() << ".dat";
        string up_name;
        up_name = strup_name.str();
        ofstream fout1(up_name.c_str(), std::ios::app);
        fout1 << "Random seed used on rank " << rank << ": "
              << param->getRandomSeed() << endl;
        fout1.close();

        // initialize init object
        Init init(nn);

        // initialize group
        Group group;

        // initialize Glauber class
        messager << "[main::main]: Init Glauber on rank " << param->getMPIRank()
                 << " ... ";
        messager.flush("info");
        Glauber glauber;
        {
            IPG_PROFILE_SCOPE("glauber.initialize");
            glauber.initGlauber(
                param->getSigmaNN(), param->getTarget(), param->getProjectile(),
                param->getb(), param->getSetWSDeformParams(), param->getR_WS(),
                param->getA_WS(), param->getBeta2(), param->getBeta3(),
                param->getBeta4(), param->getGamma(), param->getForceDmin(),
                param->getDmin(), param->getWSdR_np(), param->getWSda_np(),
                100);
        }

        // measure and output eccentricity, triangularity
        // init.eccentricity(lat, &group, param, random, glauber);

        // initialize evolution object
        Evolution evolution(nn);

        // either read k_T spectrum from file or do a fresh start
        if (param->getReadMultFromFile() == 1) {
            evolution.readNkt(param);
        } else {
            // clean files
            // stringstream strNpartdNdy_name;
            // strNpartdNdy_name << "NpartdNdy" << rank << ".dat";
            // string NpartdNdy_name;
            // NpartdNdy_name = strNpartdNdy_name.str();

            // ofstream foutNN(NpartdNdy_name.c_str(),ios::out);
            // foutNN.close();

            // stringstream strNpartdNdyH_name;
            // strNpartdNdyH_name << "NpartdNdyHadrons" << rank << ".dat";
            // string NpartdNdyH_name;
            // NpartdNdyH_name = strNpartdNdyH_name.str();

            // ofstream foutNNH(NpartdNdyH_name.c_str(),ios::out);
            // foutNNH.close();

            // stringstream strNpartdEdy_name;
            // strNpartdEdy_name << "NpartdEdy" << param->getEventId() <<
            // ".dat"; string NpartdEdy_name; NpartdEdy_name =
            // strNpartdEdy_name.str();

            // ofstream foutE(NpartdEdy_name.c_str(),ios::out);
            // foutE.close();

            // stringstream strdNdy_name;
            // strdNdy_name << "dNdy" << param->getEventId() << ".dat";
            // string dNdy_name;
            // dNdy_name = strdNdy_name.str();

            // ofstream foutN(dNdy_name.c_str(),ios::out);
            // foutN.close();

            // stringstream strCorr_name;
            // strCorr_name << "Corr" << param->getEventId() << ".dat";
            // string Corr_name;
            // Corr_name = strCorr_name.str();

            // ofstream foutCorr(Corr_name.c_str(),ios::out);
            // foutCorr.close();

            // stringstream strPhiMult_name;
            // strPhiMult_name << "MultPhi" << param->getEventId() << ".dat";
            // string PhiMult_name;
            // PhiMult_name = strPhiMult_name.str();

            // ofstream foutPhiMult(PhiMult_name.c_str(),ios::out);
            // foutPhiMult.close();

            // stringstream strPhi2ParticleMult_name;
            // strPhi2ParticleMult_name << "MultPhi2Particle" <<
            // param->getEventId()
            // << ".dat"; string Phi2ParticleMult_name; Phi2ParticleMult_name =
            // strPhi2ParticleMult_name.str();

            // ofstream
            // foutPhi2ParticleMult(Phi2ParticleMult_name.c_str(),ios::out);
            // foutPhi2ParticleMult.close();

            // stringstream strPhiMultHad_name;
            // strPhiMultHad_name << "MultPhiHadrons" << param->getEventId() <<
            // ".dat"; string PhiMultHad_name; PhiMultHad_name =
            // strPhiMultHad_name.str();

            // ofstream foutPhiMultHad(PhiMultHad_name.c_str(),ios::out);
            // foutPhiMultHad.close();

            // stringstream strPhi2ParticleMultHad_name;
            // strPhi2ParticleMultHad_name << "MultPhiHadrons2Particle" <<
            // param->getEventId() << ".dat"; string Phi2ParticleMultHad_name;
            // Phi2ParticleMultHad_name = strPhi2ParticleMultHad_name.str();

            // ofstream
            // foutPhi2ParticleMultHad(Phi2ParticleMultHad_name.c_str(),ios::out);
            // foutPhi2ParticleMultHad.close();

            // stringstream strame_name;
            // strame_name << "AverageMaximalEpsilon" << param->getEventId() <<
            // ".dat"; string ame_name; ame_name = strame_name.str();

            // ofstream foutEpsA(ame_name.c_str(),ios::out);
            // foutEpsA.close();

            // stringstream strepsx_name;
            // strepsx_name << "eps-x" << param->getEventId() << ".dat";
            // string epsx_name;
            // epsx_name = strepsx_name.str();

            // ofstream foutEpsX(epsx_name.c_str(),ios::out);
            // foutEpsX.close();

            // stringstream strdEdy_name;
            // strdEdy_name << "dEdy" << param->getEventId() << ".dat";
            // string dEdy_name;
            // dEdy_name = strdEdy_name.str();

            // ofstream foutdE(dEdy_name.c_str(),ios::out);
            // foutdE.close();

            // stringstream straniso_name;
            // straniso_name << "anisotropy" << param->getEventId() << ".dat";
            // string aniso_name;
            // aniso_name = straniso_name.str();

            // ofstream foutAni(aniso_name.c_str(),ios::out);
            // foutAni.close();

            // stringstream strecc_name;
            // strecc_name << "eccentricities" << param->getEventId() << ".dat";
            // string ecc_name;
            // ecc_name = strecc_name.str();

            // ofstream foutEcc(ecc_name.c_str(),ios::out);
            // foutEcc.close();

            // stringstream strmult_name;
            // strmult_name << "multiplicity" << param->getEventId() << ".dat";
            // string mult_name;
            // mult_name = strmult_name.str();
            // ofstream foutmult(mult_name.c_str(),ios::out);
            // foutmult.close();

            // stringstream strmult2_name;
            // strmult2_name << "multiplicityCorr" << param->getEventId() <<
            // ".dat"; string mult2_name; mult2_name = strmult2_name.str();
            // ofstream foutmult2(mult2_name.c_str(),ios::out);
            // foutmult2.close();

            // stringstream strmult3_name;
            // strmult3_name << "multiplicityCorrFromPhi" << param->getEventId()
            // <<
            // ".dat"; string mult3_name; mult3_name = strmult3_name.str();
            // ofstream foutmult3(mult3_name.c_str(),ios::out);
            // foutmult3.close();

            // stringstream strmult4_name;
            // strmult4_name << "multiplicityCorrFromPhiHadrons" <<
            // param->getEventId() << ".dat"; string mult4_name; mult4_name =
            // strmult4_name.str(); ofstream
            // foutmult4(mult4_name.c_str(),ios::out); foutmult4.close();
        }

        // Keep the lattice lifetime inside this block so destruction is timed
        // before the per-event profile is written.
        {
            // allocate lattice
            Lattice lat(param, param->getSize());
            messager.info("[main::main]: Lattice generated.");

            param->setSuccess(0);

            // initialize U-fields on the lattice
            InitializationMethod init_method;
            if (param->getReadInitialWilsonLines() == 0) {
                init_method = InitializationMethod::SampleColorCharges;
            } else {
                init_method = (param->getReadInitialWilsonLines() == 1)
                                  ? InitializationMethod::ReadWlineText
                                  : InitializationMethod::ReadWlineBinary;
            }
            // First generate the V
            init.init(&lat, &group, param, random, &glauber, init_method);

            if (param->getUseJIMWLK()) {
                messager.info("[main::main]: Start JIMWLK");
                JIMWLK jimwlkSolver(*param, &group, &lat, random);
                jimwlkSolver.evolution();
                messager.info("[main::main]: Finish JIMWLK");

                if (param->getWriteWilsonLines() > 0) {
                    std::stringstream s1;
                    s1 << "Final_x_"
                       << std::to_string(param->getJimwlk_x_projectile())
                       << "_";
                    lat.writeWilsonLines(
                        s1.str(), param, NucleusRole::Projectile);
                    std::stringstream s2;
                    s2 << "Final_x_"
                       << std::to_string(param->getJimwlk_x_target()) << "_";
                    lat.writeWilsonLines(s2.str(), param, NucleusRole::Target);
                }
            }

            if (param->getMode() == 1) {
                while (param->getSuccess() == 0) {
                    // sample collision impact parameter
                    // and compute Npart, Ncoll,etc, and check if there was a
                    // collision
                    init.sampleImpactParameter(param);
                    init.computeCollisionGeometryQuantities(&lat, param);
                }
                init.shiftFieldsWithImpactParameter(&lat, param);
                init.initializeForwardLightCone(&lat, param);
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
            if (param->getWriteOutputsToHDF5() == 1) {
                IPG_PROFILE_SCOPE("output.hdf5_collect_event");
                int status = 0;
                stringstream h5output_filename;
                h5output_filename << "RESULTS_rank" << rank;
                stringstream collect_command;
                collect_command
                    << "python3 utilities/combine_events_into_hdf5.py ."
                    << " --output_filename " << h5output_filename.str()
                    << " --event_id " << param->getEventId();
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
                ipg::writeLatticeFingerprint(&lat, rank, param->getEventId());
            }
        }  // lattice lifetime

        ipg::Profiler::instance().endEvent();
    }

    delete random;
    delete param;

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

int readInput(
    Setup *setup, Parameters *param, int argc, char *argv[], int rank) {
    // the first given argument is taken to be the input file name
    // if none is given, that file name is "input"
    PrettyOstream messager;
    string file_name;
    if (argc > 1) {
        file_name = argv[1];
        if (rank == 0) {
            messager << "[main::readInput]: Using file name \"" << file_name
                     << "\".";
            messager.flush("info");
        }
    } else {
        file_name = "input";
        if (rank == 0) {
            messager << "[main::readInput]: No input file name given. Using "
                        "default \""
                     << file_name << "\".";
            messager.flush("info");
        }
    }

    // read and set all the parameters in the "param" object of class
    // "Parameters"
    if (rank == 0) {
        messager << "[main::readInput]: Reading parameters from file ... ";
        // Flush immediately rather than deferring to the "done." message
        // far below: if any of the reads that follow hits a missing
        // key/file and calls exit(1), this is the only indication that
        // parameter parsing had even started.
        messager.flush("info");
    }
    param->setNucleusQsTableFileName(
        setup->stringFind(file_name, "NucleusQsTableFileName"));
    param->setNucleonPositionsFromFile(
        setup->iFind(file_name, "nucleonPositionsFromFile"));
    param->setNuclearConfigurationsPath(setup->stringFindOptional(
        file_name, "nuclearConfigurationsPath", "./nucleusConfigurations"));
    param->setTarget(setup->stringFind(file_name, "Target"));
    param->setProjectile(setup->stringFind(file_name, "Projectile"));
    param->setMode(setup->iFind(file_name, "mode"));
    param->setRunningCoupling(setup->iFind(file_name, "runningCoupling"));
    param->setL(setup->dFind(file_name, "L"));
    param->setLOutput(setup->dFind(file_name, "LOutput"));
    param->setBG(setup->dFind(file_name, "BG"));
    param->setBGq(setup->dFind(file_name, "BGq"));
    param->setBGqVar(setup->dFind(file_name, "BGqVar"));
    param->setDqmin(setup->dFind(file_name, "dqMin"));
    param->setOmega(setup->dFind(file_name, "omega"));
    param->setMuZero(setup->dFind(file_name, "muZero"));
    param->setc(setup->dFind(file_name, "c"));
    param->setSize(setup->iFind(file_name, "size"));
    param->setSizeOutput(setup->iFind(file_name, "sizeOutput"));
    param->setEtaSizeOutput(setup->iFind(file_name, "etaSizeOutput"));
    param->setDetaOutput(setup->dFind(file_name, "detaOutput"));
    param->setUseFluctuatingx(setup->iFind(file_name, "useFluctuatingx"));
    param->setInverseQsForMaxTime(
        setup->iFind(file_name, "inverseQsForMaxTime"));
    param->setSeed(setup->uLLIFind(file_name, "seed"));
    param->setUseSeedList(setup->iFind(file_name, "useSeedList"));
    param->setNy(setup->iFind(file_name, "Ny"));
    param->setRoots(setup->dFind(file_name, "roots"));
    param->setg(setup->dFind(file_name, "g"));
    param->setm(setup->dFind(file_name, "m"));
    param->setJacobianm(setup->dFind(file_name, "Jacobianm"));
    param->setSigmaNN(setup->dFind(file_name, "SigmaNN"));
    param->setRmax(setup->dFind(file_name, "rmax"));
    param->setUVdamp(setup->dFind(file_name, "UVdamp"));
    param->setSetWSDeformParams(setup->iFind(file_name, "setWSDeformParams"));
    if (param->getSetWSDeformParams()) {
        param->setR_WS(setup->dFind(file_name, "R_WS"));
        param->setA_WS(setup->dFind(file_name, "a_WS"));
        param->setBeta2(setup->dFind(file_name, "beta2"));
        param->setBeta3(setup->dFind(file_name, "beta3"));
        param->setBeta4(setup->dFind(file_name, "beta4"));
        param->setGamma(setup->dFind(file_name, "gamma"));
        param->setWSdR_np(setup->dFind(file_name, "dR_np"));
        param->setWSda_np(setup->dFind(file_name, "da_np"));
    }
    // Glauber::findNucleusData applies forceDminFlag/d_min unconditionally
    // (unlike the other deform params above, which it only applies when
    // setWSDeformParams is set), so these must always be read.
    param->setForceDmin(setup->dFind(file_name, "force_dmin_flag"));
    param->setDmin(setup->dFind(file_name, "d_min"));
    param->setbmin(setup->dFind(file_name, "bmin"));
    param->setbmax(setup->dFind(file_name, "bmax"));
    param->setRotateReactionPlane(
        setup->iFind(file_name, "rotateReactionPlane"));
    param->setComputeGluonMultiplicity(
        setup->iFind(file_name, "computeGluonMultiplicity"));
    param->setQsmuRatio(setup->dFind(file_name, "QsmuRatio"));
    param->setUsePseudoRapidity(setup->dFind(file_name, "usePseudoRapidity"));
    param->setRapidityA(setup->dFind(file_name, "RapidityA"));
    param->setRapidityB(setup->dFind(file_name, "RapidityB"));
    param->setUseNucleus(setup->iFind(file_name, "useNucleus"));
    param->setUseGaussian(setup->iFind(file_name, "useGaussian"));
    param->setlightNucleusOption(setup->iFind(file_name, "lightNucleusOption"));
    param->setPolarizationProjectile(
        setup->iFind(file_name, "polariztionProjectile"));
    param->setPolarizationTarget(setup->iFind(file_name, "polariztionTarget"));
    param->setPolarizationProjectileJz(
        setup->iFind(file_name, "polarizationProjectileJz"));
    param->setPolarizationTargetJz(
        setup->iFind(file_name, "polarizationTargetJz"));
    if (param->getPolarizationProjectile() != 0
        || param->getPolarizationTarget() != 0) {
        param->setNucleonPositionsFromFile(1);
    }
    param->setg2mu(setup->dFind(file_name, "g2mu"));
    param->setMaxtime(setup->dFind(file_name, "maxtime"));
    double lattice_a = param->getL() / static_cast<double>(param->getSize());
    // param->setdtau(setup->dFind(file_name, "dtau"));
    //   int iTimeSteps = static_cast<int>(10 * param->getMaxtime() /
    //   lattice_a) + 1;
    int iTimeSteps = static_cast<int>(10 * param->getMaxtime() / lattice_a);
    param->setdtau(param->getMaxtime() / (iTimeSteps * lattice_a));
    param->setRunWithQs(setup->iFind(file_name, "runWith0Min1Avg2MaxQs"));
    param->setRunWithkt(setup->iFind(file_name, "runWithkt"));
    param->setRunWithLocalQs(setup->iFind(file_name, "runWithLocalQs"));
    param->setRunWithThisFactorTimesQs(
        setup->dFind(file_name, "runWithThisFactorTimesQs"));
    param->setxFromThisFactorTimesQs(
        setup->dFind(file_name, "xFromThisFactorTimesQs"));
    param->setLinearb(setup->iFind(file_name, "samplebFromLinearDistribution"));
    param->setWriteOutputs(setup->iFind(file_name, "writeOutputs"));
    param->setWriteEpsilonUHydro(
        setup->iFindOptional(file_name, "writeEpsilonUHydro", 1));
    param->setWriteTmunuBinary(
        setup->iFindOptional(file_name, "writeTmunuBinary", 1));
    param->setWriteOutputsToHDF5(setup->iFind(file_name, "writeOutputsToHDF5"));
    param->setWriteWilsonLines(setup->iFind(file_name, "writeWilsonLines"));
    param->setReadInitialWilsonLines(
        setup->iFind(file_name, "readInitialWilsonLines"));
    param->setAverageOverNuclei(
        setup->iFind(file_name, "averageOverThisManyNuclei"));
    param->setUseTimeForSeed(setup->iFind(file_name, "useTimeForSeed"));
    param->setUseFixedNpart(setup->iFind(file_name, "useFixedNpart"));
    param->setSmearQs(setup->iFind(file_name, "smearQs"));
    param->setSmearingWidth(setup->dFind(file_name, "smearingWidth"));
    param->setGaussianWounding(setup->iFind(file_name, "gaussianWounding"));
    param->setReadMultFromFile(setup->iFind(file_name, "readMultFromFile"));
    param->setProtonAnisotropy(setup->dFind(file_name, "protonAnisotropy"));
    param->setUseConstituentQuarkProton(
        setup->dFind(file_name, "useConstituentQuarkProton"));
    param->setNqBase(setup->dFind(file_name, "useConstituentQuarkProton"));
    param->setNqFluc(setup->dFind(file_name, "NqFluc"));
    param->setUseSmoothNucleus(setup->iFind(file_name, "useSmoothNucleus"));
    param->setShiftConstituentQuarkProtonOrigin(
        setup->dFind(file_name, "shiftConstituentQuarkProtonOrigin"));
    param->setMinimumQs2ST(setup->iFind(file_name, "minimumQs2ST"));
    param->setSubNucleonParamType(
        setup->iFind(file_name, "SubNucleonParamType"));
    param->setSubNucleonParamSet(setup->iFind(file_name, "SubNucleonParamSet"));
    if (param->getSubNucleonParamType() > 0) {
        param->loadPosteriorParameterSets(param->getSubNucleonParamType());
    }

    // JIMWLK parameters
    param->setUseJIMWLK(setup->iFind(file_name, "useJIMWLK"));
    param->setSimpleLangevin(setup->iFind(file_name, "simpleLangevin"));
    param->setMu0_jimwlk(setup->dFind(file_name, "mu0_jimwlk"));
    param->setLambdaQCD_jimwlk(
        setup->dFind(file_name, "Lambda_QCD_jimwlk"));  // in units of g^2mu
    param->setm_jimwlk(setup->dFind(file_name, "m_jimwlk"));
    param->setJimwlk_alphas(setup->dFind(file_name, "alphas_jimwlk"));
    param->setDs_jimwlk(setup->dFind(file_name, "Ds_jimwlk"));
    param->setJimwlk_x_projectile(
        setup->dFind(file_name, "x_projectile_jimwlk"));
    param->setJimwlk_x_target(setup->dFind(file_name, "x_target_jimwlk"));
    param->setJimwlk_x0(setup->dFind(file_name, "jimwlk_ic_x"));
    param->setSaveSnapshots(setup->iFind(file_name, "saveSnapshots"));
    param->setxSnapshotList(setup->listFind(file_name, "xSnapshotList"));

    if (rank == 0) {
        messager << "[main::readInput]: Finished reading parameters.";
        messager.flush("info");
    }

    return 0;
}

void writeparams(Parameters *param) {
    // write the used parameters into file "usedParameters.dat" as a double
    // check for later
    time_t rawtime = time(0);
    stringstream strup_name;
    strup_name << "usedParameters" << param->getEventId() << ".dat";
    string up_name;
    up_name = strup_name.str();

    ofstream fout1(up_name.c_str(), std::ios::out);
    char *timestring = ctime(&rawtime);
    fout1 << "File created on " << timestring << endl;
    fout1 << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
    fout1 << "Used parameters by IP-Glasma v1.3" << endl;
    fout1 << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
    fout1 << " " << endl;
    fout1 << " Output by readInput in main.cpp: " << endl;
    fout1 << " " << endl;
    fout1 << "Program run in mode " << param->getMode() << endl;
    fout1 << "Nc 3" << endl;
    fout1 << "size " << param->getSize() << endl;
    fout1 << "lattice spacing a "
          << param->getL() / static_cast<double>(param->getSize()) << " fm "
          << endl;
    fout1 << "Ny " << param->getNy() << endl;
    fout1 << "Projectile " << param->getProjectile() << endl;
    fout1 << "Target " << param->getTarget() << endl;
    if (param->getUseConstituentQuarkProton() > 0) {
        fout1 << "Nucleons consists of "
              << param->getUseConstituentQuarkProton() << " constituent quarks"
              << endl;
        if (param->getShiftConstituentQuarkProtonOrigin())
            fout1 << "... constituent quark center of mass moved to origin"
                  << endl;
    }
    fout1 << "Smooth nucleus " << param->getUseSmoothNucleus() << endl;
    fout1 << "Gaussian wounding " << param->getGaussianWounding() << endl;
    fout1 << "Using fluctuating x=Qs/root(s) " << param->getUseFluctuatingx()
          << endl;
    if (param->getRunWithkt() == 0)
        fout1 << "Using local Qs to run " << param->getRunWithLocalQs() << endl;
    else
        fout1 << "running alpha_s with k_T" << endl;
    fout1 << "QsmuRatio " << param->getQsmuRatio() << endl;
    fout1 << "smeared mu " << param->getSmearQs() << endl;
    fout1 << "m " << param->getm() << endl;
    fout1 << "rmax " << param->getRmax() << endl;
    fout1 << "UVdamp " << param->getUVdamp() << endl;
    fout1 << "writeTmunuBinary " << param->getWriteTmunuBinary() << endl;
    if (param->getSetWSDeformParams()) {
        fout1 << "setWSDeformParams " << param->getSetWSDeformParams() << endl;
        fout1 << "R_WS " << param->getR_WS() << endl;
        fout1 << "a_WS " << param->getA_WS() << endl;
        fout1 << "beta2 " << param->getBeta2() << endl;
        fout1 << "beta3 " << param->getBeta3() << endl;
        fout1 << "beta4 " << param->getBeta4() << endl;
        fout1 << "gamma " << param->getGamma() << endl;
    }
    if (param->getSmearQs() == 1) {
        fout1 << "smearing width " << param->getSmearingWidth() << endl;
    }
    fout1.close();
}

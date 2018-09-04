#include "mpi.h"
#include <stdio.h>
#include <string>
#include <cmath>
#include <ctime>
#include <iostream>
#include <complex>
#include <fstream>
#include <vector>

#include "Setup.h"
#include "Init.h"
#include "Random.h"
#include "FFT.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Lattice.h"
#include "Evolution.h"
#include "Spinor.h"
#include "MyEigen.h"
#include "Fragmentation.h"
#include "NIKLAS_grids.h"
#include "NIKLAS_kinematics.h"
#include "NIKLAS_interpolation.h"
#include "NIKLAS_random_seed.h"
#include "NIKLAS_wilsonlines_IO.h"
#include "NIKLAS_compute_dipole_amplitude.h"


#define _SECURE_SCL 0
#define _HAS_ITERATOR_DEBUGGING 0
using namespace std;

int readInput(Setup *setup, Parameters *param, int argc, char *argv[], int MPI_ID);

// main program
int main(int argc, char *argv[])
{

    MPI::Init(argc, argv);
    
    int MPI_ID = MPI::COMM_WORLD.Get_rank(); //number of current processor
    int MPI_SIZE = MPI::COMM_WORLD.Get_size(); //total number of processors
    
    if(MPI_ID==0)
    {
        std::cout << "# NIKLAS' BASTARD VERSION OF IPGLASMA" << std::endl;
    }
    
    // initialize helper class objects
    MyEigen *myeigen;
    myeigen = new MyEigen();
    //myeigen->test();
    
    Setup *setup;
    setup = new Setup();
    Random *random;
    random = new Random();
    Parameters *param;
    param = new Parameters();
    
    param->setMPIRank(MPI_ID);
    param->setMPISize(MPI_SIZE);
    param->setSuccess(0);
    
    // read parameters from file
    readInput(setup, param, argc, argv, MPI_ID);
    
    
    
    int cells = param->getSize()*param->getSize();
    int Nc2m1 = param->getNc()*param->getNc()-1; // N_c^2-1
    int nn[2];
    int pos;
    double x,y;
    double ds = param->getDs();
    nn[0]=param->getSize();
    nn[1]=param->getSize();
    
    stringstream strup_name;
    strup_name << "usedParameters" << param->getMPIRank() << ".dat";
    string up_name;
    up_name = strup_name.str();
    
    //initialize init object
    Init *init;
    init = new Init(nn);
    
    // initialize group
    Group *group;
    group = new Group(param->getNc());
    
    // initialize Glauber class
    cout << "Init Glauber on MPI_ID " << param->getMPIRank() << " ... ";
    Glauber *glauber;
    glauber = new Glauber;
    
    // measure and output eccentricity, triangularity
    // init->eccentricity(lat, group, param, random, glauber);
    
    // initialize evolution object
    Evolution *evolution;
    evolution = new Evolution(nn);
    
    // initialize Glauber
    glauber->initGlauber(param->getSigmaNN(), param->getTarget(), param->getProjectile(), param->getb(), 100, MPI_ID);
    //cout << "done." << endl;
    
    
    // either read k_T spectrum from file or do a fresh start
    if (param->getReadMultFromFile()==1)
    evolution->readNkt(param);
    else
    {
        // clean files
        stringstream strNpartdNdy_name;
        strNpartdNdy_name << "NpartdNdy" << MPI_ID << ".dat";
        string NpartdNdy_name;
        NpartdNdy_name = strNpartdNdy_name.str();
        
        ofstream foutNN(NpartdNdy_name.c_str(),ios::out);
        foutNN.close();
        
        stringstream strNpartdNdyH_name;
        strNpartdNdyH_name << "NpartdNdyHadrons" << MPI_ID << ".dat";
        string NpartdNdyH_name;
        NpartdNdyH_name = strNpartdNdyH_name.str();
        
        ofstream foutNNH(NpartdNdyH_name.c_str(),ios::out);
        foutNNH.close();
        
        stringstream strmult_name;
        strmult_name << "multiplicity" << param->getMPIRank() << ".dat";
        string mult_name;
        mult_name = strmult_name.str();
        ofstream foutmult(mult_name.c_str(),ios::out);
        foutmult.close();
    }
    
    while(param->getSuccess()==0)
    {
        param->setSuccess(0);
        // allocate lattice
        Lattice *lat;
        lat = new Lattice(param, param->getNc(), param->getSize());
        
        //initialize random generator using time and seed from input file
        unsigned long long int rnum;
        
        if(param->getUseSeedList()==0)
        {
            if(param->getUseTimeForSeed()==1)
            {
                rnum=time(0)+param->getSeed()*10000;
            }
            else
            {
                rnum = param->getSeed();
                cout << "Random seed = " << rnum+(MPI_ID*1000) << " - entered directly +MPI_ID*1000."  << endl;
            }
            
            param->setRandomSeed(rnum+MPI_ID*1000);
            if(param->getUseTimeForSeed()==1)
            {
                cout << "Random seed = " << param->getRandomSeed() << " made from time "
                << rnum-param->getSeed()-(MPI_ID*1000) << " and argument (+1000*MPI_ID) "
                << param->getSeed()+(MPI_ID*1000) << endl;
            }
            
            random->init_genrand64(rnum+MPI_ID*1000);
        }
        else
        {
            ifstream fin;
            fin.open("seedList");
            unsigned long long int seedList[MPI_SIZE];
            if(fin)
            {
                for (int i=0; i<MPI_SIZE; i++)
                {
                    if (!fin.eof())
                    {
                        fin >> seedList[i];
                    }
                    else
                    {
                        cerr << "Error: Not enough random seeds for the number of processors selected. Exiting." << endl;
                        exit(1);
                    }
                }
            }
            else
            {
                cerr << "Random seed file 'seedList' not found. Exiting." << endl;
                exit(1);
            }
            fin.close();
            
            param->setRandomSeed(seedList[MPI_ID]);
            random->init_genrand64(seedList[MPI_ID]);
            cout << "Random seed on MPI_ID " << MPI_ID << " = " << seedList[MPI_ID] << " read from list."  << endl;
        }
        
        ofstream fout1(up_name.c_str(),ios::app);
        fout1 << "Random seed used on MPI_ID " << MPI_ID << ": " << param->getRandomSeed() << endl;
        fout1.close();
        // initialize gsl random number generator (used for non-Gaussian distributions)
        random->gslRandomInit(rnum);
        // initialize U-fields on the lattice
        //      init->init(lat, group, param, random, glauber);
        int READFROMFILE = 0;
        
        // ------------------------------------------------------------------------------------
        // THIS IS WHAT I WANT ----------------------------------------------------------------
        // -----------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------
        // ----------------------------------------------------------------------------------- 
        // COMMPUTING WILSON LINES //
        std::cout << "# COMPUTING WILSON LINES" << std::endl;
        init->init(lat, group, param, random, glauber, READFROMFILE);
        
        // GENERATE SOME OUTPUT HERE // 
        WILSONLINES_IO::PrintWilsonLines(lat, param, "WilsonLine.txt");
        
        INT N, Nc;
        DOUBLE L, a;
        
        Parameters *Myparam;
        Lattice *Mylat=WILSONLINES_IO::ReadWilsonLines(Myparam,"WilsonLine.txt", N, Nc, L, a, MPI_ID);
        
        WILSONLINES_IO::CompareLattice(lat, Mylat, param);

        std::cout << "# DONE" << std::endl;
        
        // -----------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------
        
        if(param->getSuccess()==0)
        {
            delete lat;
            continue;
        }

        delete init;
        delete random;
        delete glauber;
        delete lat;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    delete group;
    delete evolution;
    delete param;
    delete setup;
    delete myeigen;
    MPI::Finalize();
    return 0;
}





























int readInput(Setup *setup, Parameters *param, int argc, char *argv[], int MPI_ID)
{
    // the first given argument is taken to be the input file name
    // if none is given, that file name is "input"
    //cout << "Opening input file ... " << endl;
    string file_name;
    if (argc>1)
    {
        file_name = argv[1];
        if(MPI_ID==0)
        cout << "Using file name \"" << file_name << "\"." << endl;
    }
    else
    {
        file_name = "input";
        if(MPI_ID==0)
        cout << "No input file name given. Using default \"" << file_name << "\"." << endl;
    }
    
    
    // read and set all the parameters in the "param" object of class "Parameters"
    if(MPI_ID==0)
    cout << "Reading parameters from file ... ";
    param->setNucleusQsTableFileName(setup->StringFind(file_name,"NucleusQsTableFileName"));
    param->setNucleonPositionsFromFile(setup->IFind(file_name,"nucleonPositionsFromFile"));
    param->setTarget(setup->StringFind(file_name,"Target"));
    param->setProjectile(setup->StringFind(file_name,"Projectile"));
    param->setMode(setup->IFind(file_name,"mode"));
    param->setRunningCoupling(setup->IFind(file_name,"runningCoupling"));
    param->setL(setup->DFind(file_name,"L"));
    param->setLOutput(setup->DFind(file_name,"LOutput"));
    param->setBG(setup->DFind(file_name,"BG"));
    param->setBGq(setup->DFind(file_name,"BGq"));
    param->setMuZero(setup->DFind(file_name,"muZero"));
    param->setc(setup->DFind(file_name,"c"));
    param->setSize(setup->IFind(file_name,"size"));
    param->setSizeOutput(setup->IFind(file_name,"sizeOutput"));
    param->setEtaSizeOutput(setup->IFind(file_name,"etaSizeOutput"));
    param->setDetaOutput(setup->DFind(file_name,"detaOutput"));
    param->setUseFluctuatingx(setup->IFind(file_name,"useFluctuatingx"));
    param->setNc(setup->IFind(file_name,"Nc"));
    param->setInverseQsForMaxTime(setup->IFind(file_name,"inverseQsForMaxTime"));
    param->setSeed(setup->ULLIFind(file_name,"seed"));
    param->setUseSeedList(setup->IFind(file_name,"useSeedList"));
    param->setNy(setup->IFind(file_name,"Ny"));
    param->setRoots(setup->DFind(file_name,"roots"));
    param->setNu(setup->DFind(file_name,"tDistNu"));
    param->setUseFatTails(setup->IFind(file_name,"useFatTails"));
    param->setg(setup->DFind(file_name,"g"));
    param->setm(setup->DFind(file_name,"m"));
    param->setJacobianm(setup->DFind(file_name,"Jacobianm"));
    param->setSigmaNN(setup->DFind(file_name,"SigmaNN"));
    param->setRmax(setup->DFind(file_name,"rmax"));
    param->setUVdamp(setup->DFind(file_name,"UVdamp"));
    param->setbmin(setup->DFind(file_name,"bmin"));
    param->setbmax(setup->DFind(file_name,"bmax"));
    param->setQsmuRatio(setup->DFind(file_name,"QsmuRatio"));
    param->setUsePseudoRapidity(setup->DFind(file_name,"usePseudoRapidity"));
    param->setRapidity(setup->DFind(file_name,"Rapidity"));
    param->setUseNucleus(setup->IFind(file_name,"useNucleus"));
    param->setg2mu(setup->DFind(file_name,"g2mu"));
    param->setMaxtime(setup->DFind(file_name,"maxtime"));
    param->setdtau(setup->DFind(file_name,"dtau"));
    // param->setxExponent(setup->DFind(file_name,"xExponent")); //  is now obsolete
    param->setRunWithQs(setup->IFind(file_name,"runWith0Min1Avg2MaxQs"));
    param->setRunWithkt(setup->IFind(file_name,"runWithkt"));
    param->setRunWithLocalQs(setup->IFind(file_name,"runWithLocalQs"));
    param->setRunWithThisFactorTimesQs(setup->DFind(file_name,"runWithThisFactorTimesQs"));
    param->setxFromThisFactorTimesQs(setup->DFind(file_name,"xFromThisFactorTimesQs"));
    param->setLinearb(setup->IFind(file_name,"samplebFromLinearDistribution"));
    param->setWriteOutputs(setup->IFind(file_name,"writeOutputs"));
    param->setWriteEvolution(setup->IFind(file_name,"writeEvolution"));
    param->setAverageOverNuclei(setup->IFind(file_name,"averageOverThisManyNuclei"));
    param->setUseTimeForSeed(setup->IFind(file_name,"useTimeForSeed"));
    param->setUseFixedNpart(setup->IFind(file_name,"useFixedNpart"));
    param->setSmearQs(setup->IFind(file_name,"smearQs"));
    param->setSmearingWidth(setup->DFind(file_name,"smearingWidth"));
    param->setGaussianWounding(setup->IFind(file_name,"gaussianWounding"));
    param->setReadMultFromFile(setup->IFind(file_name,"readMultFromFile"));
    param->setProtonAnisotropy(setup->DFind(file_name,"protonAnisotropy"));
    param->setUseConstituentQuarkProton(setup->DFind(file_name,"useConstituentQuarkProton"));
    if(MPI_ID==0)
    cout << "done." << endl;
    
    // write the used parameters into file "usedParameters.dat" as a double check for later
    time_t rawtime = time(0);
    stringstream strup_name;
    strup_name << "usedParameters" << param->getMPIRank() << ".dat";
    string up_name;
    up_name = strup_name.str();
    
    fstream fout1(up_name.c_str(),ios::out);
    char * timestring = ctime(&rawtime);
    fout1 << "File created on " << timestring << endl;
    fout1 << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
    fout1 << "Used parameters by IP-Glasma v1.3" << endl;
    fout1 << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
    fout1 << " " << endl;
    fout1 << " Output by readInput in main.cpp: " << endl;
    fout1 << " " << endl;
    fout1 << "Program run in mode " << param->getMode() << endl;
    fout1 << "Nc " << param->getNc() << endl;
    fout1 << "size " << param->getSize() << endl;
    fout1 << "lattice spacing a " << param->getL()/static_cast<double>(param->getSize()) << " fm " << endl;
    fout1 << "Ny " << param->getNy() << endl;
    fout1 << "Projectile " << param->getProjectile() << endl;
    fout1 << "Target " << param->getTarget() << endl;
    fout1 << "Gaussian wounding " << param->getGaussianWounding() << endl;
    fout1 << "Using fluctuating x=Qs/root(s) " << param->getUseFluctuatingx() << endl;
    if( param->getRunWithkt()==0)
    fout1 << "Using local Qs to run " << param->getRunWithLocalQs() << endl;
    else
    fout1 << "running alpha_s with k_T" << endl;
    fout1 << "QsmuRatio " << param->getQsmuRatio() << endl;
    fout1 << "smeared mu " << param->getSmearQs() << endl;
    fout1 << "m " << param->getm() << endl;
    fout1 << "rmax " << param->getRmax() << endl;
    fout1 << "UVdamp " << param->getUVdamp() << endl;
    if (param->getSmearQs()==1)
    {
        fout1 << "smearing width " << param->getSmearingWidth() << endl;
    }
    fout1 << "Using fat tailed distribution " << param->getUseFatTails() << endl;
    fout1.close();
    
    return 0;
}

//
//  commandlineargs.h
//  JETS
//
//  Created by Niklas Mueller on 2/28/18.
//  Copyright © 2018 Niklas Mueller. All rights reserved.
//

#ifndef commandlineargs_h
#define commandlineargs_h

//#include "../src_niklas/IO/StringManipulation.h"
//#include "../src_niklas/IO/OutputManagment.h"

#include "StringManipulation.h"
#include "OutputManagment.h"


void ProcessCommandlineArguments(int argc,char **argv){
    
    Konfig arguments(argc,argv);
    
    char OutDir[768]="OUTPUT";
    arguments.Getval("o",OutDir);
    IO::SetOutputDirectory(OutDir);
    char OutFile[768]="test.txt";
    arguments.Getval("f",OutFile);
    IO::SetOutputFile(OutFile);
    
    // GLOBAL KINEMATICS //
    arguments.Getval("W",PAR::W);
    arguments.Getval("Q2",PAR::PHOTON::Q2);
    arguments.Getval("Ep",PAR::TARGET::Ep);
    arguments.Getval("mq",PAR::JET::mquark);
    arguments.Getval("dipole",PAR::dipole_model);
    arguments.Getval("mN",PAR::TARGET::mNucleon);
    
    // MONTE CARLO POINTS //
    arguments.Getval("MC",SIMPAR::NMCSAMPLE);
    arguments.Getval("SEED",SIMPAR::SEED);
    
    // -------------------------------------------//
    // GRIDS IN TRANSVERSE AND RAPIDITY DIRECTION //
    // -------------------------------------------//
    arguments.Getval("Ny0",SIMPAR::Ny0);
    arguments.Getval("Ny1",SIMPAR::Ny1);
    arguments.Getval("Np0",SIMPAR::Np0);
    arguments.Getval("Np1",SIMPAR::Np1);
    arguments.Getval("Ntheta",SIMPAR::Ntheta);
    
    // IN TRANSFORMED COORDINATES //
    arguments.Getval("N_abs_Delta",SIMPAR::N_abs_Delta);
    arguments.Getval("N_abs_k",SIMPAR::N_abs_k);
    
    
    // ACTUAL LIMITS //
    arguments.Getval("y0min",SIMPAR::y0min);
    arguments.Getval("y0max",SIMPAR::y0max);
    arguments.Getval("z0min",SIMPAR::z0min);
    arguments.Getval("z0max",SIMPAR::z0max);
    arguments.Getval("y1min",SIMPAR::y1min);
    arguments.Getval("y1max",SIMPAR::y1max);
    arguments.Getval("z1min",SIMPAR::z1min);
    arguments.Getval("z1max",SIMPAR::z1max);
    arguments.Getval("p0min",SIMPAR::p0min);
    arguments.Getval("p0max",SIMPAR::p0max);
    arguments.Getval("p1min",SIMPAR::p1min);
    arguments.Getval("p1max",SIMPAR::p1max);
    arguments.Getval("thetamin",SIMPAR::thetamin);
    arguments.Getval("thetamax",SIMPAR::thetamax);
    
    // ARTIFICIAL Qs DEPENDENCE //
    arguments.Getval("Qs",SIMPAR::Qsfactor);
    
    // HISTORIC NOT USED ANY MORE //
    arguments.Getval("p0",PAR::JET::abs_p0);
    arguments.Getval("p1",PAR::JET::abs_p1);
    arguments.Getval("NTrapez",SIMPAR::NTrapezodialPoints);
    
    // KINEMATIC CUTS //
    arguments.Getval("Res",CUTS::ResolutionFactor);
    arguments.Getval("Cp0min",CUTS::p0min);
    arguments.Getval("Cp0max",CUTS::p0max);
    arguments.Getval("Cp1min",CUTS::p1min);
    arguments.Getval("Cp1max",CUTS::p1max);
    arguments.Getval("Cy0min",CUTS::y0min);
    arguments.Getval("Cy0max",CUTS::y0max);
    arguments.Getval("Cy1min",CUTS::y1min);
    arguments.Getval("Cy1max",CUTS::y1max);
    arguments.Getval("CDeltamin",CUTS::Delta_min);
    arguments.Getval("CDeltamax",CUTS::Delta_max);
    arguments.Getval("Ckmin",CUTS::k_min);
    arguments.Getval("Ckmax",CUTS::k_max);
    
}






#endif /* commandlineargs_h */

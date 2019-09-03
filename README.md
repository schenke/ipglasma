# README #

IP-Glasma code with improved matrix exponential - openmp development


### openmp IP-Glasma ###

 * this is version 0.1
 * work on openmp fftw (http://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html)
 
 
## Input parameters

 - **writeOutputs**: this parameter controls output files
 	- 0: no output
 	- 1: output initial conditions e, u^\mu, and pi^{\mu\nu} for hydrodynamic simulations
 	- 2: output the initial condition for energy density according to Jazma
 	- 3: output 1 & 2
 	- 4: output initial T^{\mu\nu} for effective kinetic theory (KomPoST) simulations
 	- 5: output 1 & 4
 	- 6: output 2 & 4
 	- 7: output 1 & 2 & 4
 
 
 
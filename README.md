# README 

IP-Glasma initial condition with JIMWLK evolution

References
* Original IP-Glasma: [Schenke, Tribedy, Venugopalan, PRL 108 (2012) 252301](https://doi.org/10.1103/PhysRevLett.108.252301), [arXiv:1202.6646](https://arxiv.org/abs/1202.6646) and [Schenke, Tribedy, Venugopalan, PRC 86 (2012) 034908](https://doi.org/10.1103/PhysRevC.86.034908), [arXiv:1206.6805](https://arxiv.org/abs/1206.6805)
* JIMWLK evolution implementation: Mäntysaari, Schenke, Shen, Zhao, in preparation
* Standalone version of the JIMWLK code: https://github.com/hejajama/jimwlk


### openmp IP-Glasma 
 * this is version 0.1
 * work on openmp fftw (http://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html)
 
## Compile
To compile IP-Glasma, run `./compile_IPGlasma.sh`
Dependencies
* CMake
* FFTW

 
## Input parameters
See `src/Parameters.h` for a more detailed description of parameters that are specified in the file `input`. If a command line argument is provided, that refers to the input file that will be used.

### Lattice
- **size**: controls the size of the lattice that is `size`$^2$. 
  - Recommended to be of the form $2^n$
- **L**: the total physical extend of hte lattice (in fm)

### Initial state
- **useNucleus**
  - 1: nucleus with finite geometry
  - 0: infite target with constant color charge density controlled by `g2mu` (in lattice units)
- **Projectile** and **Target**: specify nuclei
  - Typical values: `p`, `Pb`, `Au`
  - See `src/Glauber.cpp` for all supported nuclei and details
- **m**: infrared regulator in GeV
- **BG**: controls the nucleon width, density profile is $T \sim e^{-b^2/(2B)}$
- **BGq**: controls the hot spot width (if nucleon substructure is enabled), hot spot density profile is $T_q \sim e^{-b^2/(2B_{Gq})}$
- **useConstituentQuarkProton**: control nucleon substructure
  - 0: no substructure
  - Positive value: number of hot spots 
- **shiftConstituentQuarkProtonOrigin** whether to shift the center-of-mass to origin (1) or not (0) after sampling the hot spot positions
- **smearQs**: enable (1) or disable (0) saturation scale fluctuations
- **smearingWidth**: width of the saturation scale fluctuations, parameter $\sigma$ in Eq. (23) of [arXiv:1607.01711](https://arxiv.org/pdf/1607.01711)
- **useFluctuatingx**: controls how to determine Bjorken-$x$ when generating the initial condition
  - 1: Dynanically determined $b_\perp$ dependent $x$ 
  - 0: Fixed $x$
- **Rapidity**:
  - If `useFluctuatingx 0`, then $x = 0.01 e^{-\mathrm{Rapidity}}$ (for both projectile and target)
  - If `useFluctuatingx 1`, consider particle production at rapidity $y$


### Output
 - **writeOutputs**: this parameter controls output files
 	- 0: no output
 	- 1: output initial conditions e, u^\mu, and pi^{\mu\nu} for hydrodynamic simulations
 	- 2: output the initial condition for energy density according to the Jazma model
 	- 3: output 1 & 2
 	- 4: output initial T^{\mu\nu} for the effective kinetic theory (KoMPoST) simulations
 	- 5: output 1 & 4
 	- 6: output 2 & 4
 	- 7: output 1 & 2 & 4
 
 - **writeOutputsToHDF5**: this parameter decides whether to collect all the IPGlasma output files into a hdf5 data file
 	- 0: no
 	- 1: yes	
 - **writeInitialWilsonLines**: controls if the generated Wilson lines for the saved on disc. File names depend on random seed (parameter `seed`). Wilson lines at the initial condition and after the evolution are saved.
    - 0: do not save Wilson lines
	- 1: save in text format
	- 2: save in binary format (faster I/O, smaller file size)

### JIMWLK evolution
Note that when using the JIMWLK evolution, one should use `useFluctuatingx 0` which corresponds to having a fixed $x$ at the initial state of the evoluiton.

- **useJIMWLK**: with JIMWLK (1), or no JIMWLK (0)
- **jimwlk_ic_x**: Bjorken-x at the initial condition
- **x_projectile_jimwlk**: Bjorken-$x$ to which the projectile is evolved
- **x_target_jimwlk**: Bjorken-$x$ to which the target is evolved
- **m_jimwlk**: Infrared regulator in GeV in the JIMWLK kernel, see (21) in [https://arxiv.org/pdf/2207.03712](arXiv:2207.03712)
- **alphas_jimwlk**: Coupling constant in the JIMWLK evolution
  - 0: Use running coupling 
- **Lambda_QCD_jimwlk** $\Lambda_\mathrm{QCD}$ in $\alpha_s(r)$ in GeV as in Eq. (22) of [https://arxiv.org/pdf/2207.03712](arXiv:2207.03712)
- **mu0_jimwlk**: Regulator in $\alpha_s(r)$ as in Eq. (22) of [https://arxiv.org/pdf/2207.03712](arXiv:2207.03712)
- **Ds_jimwlk**: step size in JIMWLK evolution. Recommended values
  - 0.005 with running coupling
  - 0.0005 with fixed coupling

- **simpleLangevin**: JIMWLK discretization method
  - 1: Use the simple Langevin step developed in [arXiv:1212.4825](https://arxiv.org/abs/1212.4825)

  Default parameters for the JIMWLK evolution with fluctuating proton at initial $x=0.01$ fitted to HERA vector meson production data are reported in [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712)
# README

IP-Glasma initial condition with JIMWLK evolution

References
* Original IP-Glasma: [Schenke, Tribedy, Venugopalan, PRL 108 (2012) 252301](https://doi.org/10.1103/PhysRevLett.108.252301), [arXiv:1202.6646](https://arxiv.org/abs/1202.6646) and [Schenke, Tribedy, Venugopalan, PRC 86 (2012) 034908](https://doi.org/10.1103/PhysRevC.86.034908), [arXiv:1206.6805](https://arxiv.org/abs/1206.6805)
* JIMWLK evolution implementation: Mäntysaari, Schenke, Shen, Zhao, [Phys.Rev.Lett. 135 (2025) 2, 022302](https://doi.org/10.1103/gf4y-p5j7), [arXiv:2502.05138](https://arxiv.org/abs/2502.05138) and [Phys.Rev.C 113 (2026) 3, 034914](https://doi.org/10.1103/qsd7-zwmw), [arXiv:2511.03588](https://arxiv.org/abs/2511.03588)
* Standalone version of the JIMWLK code: https://github.com/hejajama/jimwlk



## Compile
To compile IP-Glasma, run `./compile_IPGlasma.sh`. It builds the code with CMake in the directory `build/` of the repository, which it empties first and deletes afterwards, and installs the executable `ipglasma` in the repository root. If the build fails, it stops and keeps `build/`. `./compile_IPGlasma.sh noMPI` builds the code without MPI.

Dependencies
* a C++17 compiler and CMake
* FFTW 3
* GSL
* optional: MPI (used if CMake finds it) and OpenMP (not with Apple Clang)

The code can also be built with CMake directly, e.g. `cmake -B build` and `cmake --build build`; the executable is then `build/src/ipglasma`.

### Reproducible runs
By default, FFTW picks its FFT algorithms by timing them when the program starts, so two runs with the same seed can differ by floating-point rounding (around 1e-14). For bit-identical output, build with
```
cmake -B build -DIPGLASMA_DETERMINISTIC_FFT=ON
cmake --build build
```
The first run then stores the chosen FFT plans in `ipglasma_fftw_wisdom.dat` in its working directory, and later runs that find this file reuse them. Runs are only bit-identical if they use the same wisdom file, so run them in the same directory or copy the file there. Different machines or FFTW versions can still choose different plans.

## Unit tests
A [doctest](https://github.com/doctest/doctest)-based unit test suite lives under `tests/` and covers small, fast, deterministic pieces of the code (matrix algebra, RNG stream properties, input parsing, ...) rather than full lattice evolution, so it runs in a fraction of a second. It is off by default; build and run it with:
```
cmake -B build -Dunittest=ON
cmake --build build
ctest --test-dir build --output-on-failure
```
It also runs automatically on every push/PR to `devel`/`main` via GitHub Actions (see `.github/workflows/tests.yml`).

## Running
```
./ipglasma [input file] [events per rank]
mpirun -np <ranks> ./ipglasma [input file] [events per rank]
```
The input file defaults to `input`, and the number of events per MPI rank to 1; with MPI, the run generates ranks × events per rank events. The number of events per rank also enters the names of the Wilson-line files. File names in the input file (e.g. `nucleusQsTableFileName`) are relative to the working directory, and the output files are written there, apart from the Wilson-line files (see `wilsonLinePath` and [OUTPUT.md](OUTPUT.md)).


## Input parameters
Every parameter is listed once in `src/ParameterTable.cpp`, together with its default value (if it is optional) and its validity checks; see `src/Parameters.h` for a more detailed description of each parameter.

The input file has one `key value` pair per line:
- `#` starts a comment that runs to the end of the line, and blank lines are ignored.
- A line with only `EndOfFile` ends the input; everything after it is ignored. It is optional.
- Unknown keys, keys given twice, missing required keys and malformed values (e.g. `256.0` for an integer) are errors. All problems are reported at once before the run stops.
- A parameter marked "only used with …" below is only read, and only required, when that setting is on. Otherwise it can be left out, and if it is given, it is ignored.

Each event writes the values of all input parameters it used to `usedParameters<event>.dat`, followed by its random seed and collision geometry as comments. The file is itself a valid input file. Running it does not reproduce the same event, though: the random numbers also depend on the MPI rank and on the event's position in the run, `seed -1` draws a new seed, and with `subNucleonParamSet -1` a new posterior parameter set is drawn.

### Run control
- **runEvolution**: `1` samples the collision, runs the classical Yang-Mills evolution and writes the outputs. `0` stops after the Wilson lines of the two nuclei are built (and evolved with JIMWLK, see `useJIMWLK`); with `writeWilsonLines` they are written to disk
- **maxTime**: proper time $\tau$ in fm/c at which the classical Yang-Mills evolution stops
- **inverseQsForMaxTime**: `1` stops the evolution at $\tau = 1/\langle Q_s \rangle$ instead of `maxTime`, with $\langle Q_s \rangle$ the larger of the two nuclei's $Q_s$ averaged over the overlap region (see `runWithQs`); with `useJIMWLK 1` that of the initial condition, before the JIMWLK evolution. Not possible with `useNucleus 0`, which does not compute $\langle Q_s \rangle$

### Random seed
- **seed**: random seed; MPI rank $r$ uses `seed` $+ 1000 r$. `-1` draws the seed from `std::random_device` (rank $r$ again adds $1000 r$). It also enters the names of the Wilson-line files (see `writeWilsonLines`), as `0` for `seed -1`, and also with `useSeedList 1`
- **useSeedList**: `1` reads one seed per MPI rank from the file `seedList` in the working directory (rank $r$ uses the $(r+1)$-th number); it overrides `seed`

### Lattice
- **size**: controls the size of the lattice that is `size`$^2$.
  - Recommended to be of the form $2^n$
- **L**: the total physical extent of the lattice (in fm)
- **Ny**: number of longitudinal layers of color charges per nucleus, at least 1: each Wilson line is a product of `Ny` factors, each sampled with the color-charge density $g^2\mu^2/N_y$ (T. Lappi, Eur. Phys. J. C 55 (2008) 285)

### Initial state
- **useNucleus**
  - 1: nucleus with finite geometry
  - 0: infinite target with constant color charge density controlled by `g2muGeV`
- **g2muGeV** (only used with `useNucleus 0`): the constant $g^2\mu$ in GeV, positive; it is converted to lattice units with the lattice spacing $a$ = `L`/`size`, so the physical system does not depend on the lattice
- **projectile** and **target**: specify nuclei
  - Typical values: `p`, `Pb`, `Au`
  - See `src/Glauber.cpp` for all supported nuclei and details
- **sqrtS** (only used with `useFluctuatingX 1` or `usePseudoRapidity 1`): center-of-mass energy per nucleon pair $\sqrt{s}$ in GeV. It enters Bjorken $x$ with `useFluctuatingX 1` (see `xQsFactor`) and the pseudorapidity Jacobian (see `usePseudoRapidity`); it does not set `sigmaNN`
- **m**: infrared regulator in GeV
- **UVDamp**: UV damping length in GeV$^{-1}$, not negative: the propagator $1/(k^2+m^2)$ that gives the gauge field of the color charges is multiplied by $e^{-|k| \cdot \mathrm{UVDamp}}$. Only applied with `m` $\neq 0$; `0` for no damping
- **BG**: nucleon width in GeV$^{-2}$. With `nucleonModel gaussian` the nucleon's thickness is $T \sim e^{-b^2/(2B_G)}$; with `hotspots` and `strings` it is the width of the hot-spot position distribution (see `omega`)
- **nucleonModel**: transverse structure of a nucleon
  - `gaussian`: a single Gaussian of width `BG` (optionally elongated by `protonAnisotropy`)
  - `hotspots`: `Nq` Gaussian hot spots of width `BGq`, whose positions are distributed with width `BG`
  - `strings`: three hot spots, sampled as for `hotspots`, connected by strings that meet at a junction, the point with the shortest total string length (Fermat point, computed in 3D). Each hot spot is moved to a uniformly random point on its string, between the junction and its sampled position, so the nucleon is smaller than with `hotspots` for the same `BG`. The minimum distance `dqMin` applies to the sampled hot spots (the ends of the strings), not to the moved ones. Uses the hot-spot parameters except `Nq` and `NqFluc`.

  The parameters of the other models are not read.
- **subNucleonParamType**: `0` uses the input values for the nucleon substructure; `1`, `2` or `4` draw them every event from a Bayesian posterior parameter set (`tables/posterior.csv` with variable Nq, `tables/posterior_Nq3.csv` or `tables/posterior5020_Nq3.csv` with Nq = 3), selected by **subNucleonParamSet** (`-1`: random). The posterior sets are fits of hot-spot nucleons, so they require `nucleonModel hotspots`. They replace `m`, `BG`, `BGq`, `smearingWidth`, `QsMuRatio`, `dqMin` and the number of hot spots, which are then not read from the input; `usedParameters<event>.dat` lists the values an event used. `BGqVar` and `NqFluc`, which the fits did not vary, are then 0 and not read either, so an event uses exactly the fitted model.
- **protonAnisotropy** (`gaussian`): anisotropy $\xi > -1$ of the nucleon, $T \propto \sqrt{1+\xi} e^{-(b^2 + \xi (\vec b \cdot \hat n)^2)/(2B_G)}$ with a random direction $\hat n$ per nucleon, so the nucleon is narrower along $\hat n$ for $\xi > 0$; `0` gives a round nucleon
- **BGq** (`hotspots`, `strings`): mean width $B_q$ in GeV$^{-2}$ of the hot spots, whose density profile is $T_q \sim e^{-b^2/(2B_q)}$; the posterior parameter sets fit it (see `subNucleonParamType`). All hot spots of a nucleon share one width. With `BGqVar 0` it is exactly `BGq`; otherwise it fluctuates from nucleon to nucleon around the mean `BGq` but never falls below 0.09 GeV$^{-2}$: it is 0.09 GeV$^{-2}$ plus a log-normal number with mean `BGq` $-$ 0.09 and variance `BGqVar`. So `BGq` must be larger than 0.09 GeV$^{-2}$
- **BGqVar** (`hotspots`, `strings`): variance in GeV$^{-4}$ of the hot-spot width from nucleon to nucleon (see `BGq`); `0` gives every nucleon the width `BGq`. With a posterior parameter set it is 0 (see `subNucleonParamType`)
- **dqMin** (`hotspots`, `strings`): minimum distance in fm between the hot spots of a nucleon, in 3D for `omega 1` and in the transverse plane otherwise. It is kept on a best-effort basis: a hot spot keeps its sampled radius and only its direction is redrawn, up to 100 times. With `strings` it applies to the ends of the strings
- **omega** (`hotspots`, `strings`): shape and width of the distribution of the hot-spot positions; it changes both, not only 2D against 3D
  - `1` (more precisely $|\omega - 1| < 10^{-8}$): positions in 3D, each coordinate Gaussian with variance `BG`, so the transverse distribution is a Gaussian with $\langle b^2 \rangle = 2 B_G$
  - any other value: positions in the transverse plane, at a uniformly random angle and radius $b = \sqrt{2\omega x B_G}$, where $x$ follows the density $\propto Q(1/\omega, x)$ (regularized upper incomplete gamma function), so $\langle b^2 \rangle = (1+\omega) B_G$. The shape is not a Gaussian: for $\omega < 1$ the distribution is flatter and more compact, approaching a uniform disk of radius $\sqrt{2 B_G}$ for $\omega \to 0$; for $\omega > 1$ it is more peaked at the center with a longer tail. For $\omega \to 1$ it approaches the Gaussian of `omega 1`
- **Nq** (`hotspots`): mean number of hot spots per nucleon, at least 1; a fractional value such as 2.5 gives 2 or 3 hot spots with the corresponding probabilities (plus the fluctuation set by `NqFluc`)
- **NqFluc** (`hotspots`): mean of a Poisson-distributed number of hot spots added to the `Nq` ones: a nucleon gets `Nq` hot spots (rounded as described for `Nq`) plus a Poisson number with mean `NqFluc`, so the mean is `Nq` + `NqFluc`; e.g. `Nq 2` and `NqFluc 1` give 3 on average. `0` for no fluctuation. Every nucleon has at least one hot spot
- **shiftConstituentQuarkProtonOrigin** (`hotspots`, `strings`): whether to shift the center-of-mass to origin (1) or not (0) after sampling the hot spot positions
- **smearQs**: enable (1) or disable (0) saturation scale fluctuations: the thickness of each nucleon (`gaussian`) or each hot spot (`hotspots`, `strings`) is multiplied by a log-normal factor $e^{X}/e^{\sigma^2/2}$ with mean 1, where $X$ is Gaussian with width $\sigma$ = `smearingWidth`
- **smearingWidth** (only used with `smearQs 1`): width $\sigma$ of the saturation scale fluctuations (see `smearQs`), parameter $\sigma$ in Eq. (23) of [arXiv:1607.01711](https://arxiv.org/pdf/1607.01711)
- **QsMuRatio**: ratio $Q_s/(g^2\mu)$, positive, the same for both nuclei: a cell's color-charge density is $g^2\mu = Q_s/$`QsMuRatio`, with $Q_s$ from `nucleusQsTableFileName`
- **nucleusQsTableFileName**: file with the table of $Q_s^2$ as a function of the summed nucleon thickness $T_p$ and rapidity $y$ (e.g. `qs2Adj_vs_Tp_vs_Y_240.in` in the repository root), relative to the working directory; the file is only read with `useNucleus 1`
- **projectileX** and **targetX** (only used with `useFluctuatingX 0`): Bjorken $x$ of the projectile and of the target. Without JIMWLK and with `useFluctuatingX 0`, $Q_s^2$ of each nucleus is read from the nuclear $Q_s$ table at $y = \ln(0.01/x)$, so both must not be larger than 0.01 (the table covers $0 \le y \le 10.75$, larger values use $Q_s$ at $y = 10.75$). With `useJIMWLK 1`, the nuclei are evolved from `jimwlkInitialX` (or from `readWilsonLinesX` for read Wilson lines) to these $x$, which must then not be larger than it.
- **useFluctuatingX**: controls how to determine Bjorken-$x$ when generating the initial condition; must be `0` with `useJIMWLK 1`, where $Q_s^2$ is evaluated at the fixed $x$ = `jimwlkInitialX`
  - 1: Dynamically determined $b_\perp$ dependent $x$, see `xQsFactor`
  - 0: Fixed $x$, `projectileX` and `targetX`
- **rapidity**: rapidity $y$ at which the gluon and hadron spectra are computed (`computeGluonMultiplicity`, `writeHadronSpectrum`). It does not set the $x$ of the nuclei, except with `useFluctuatingX 1` (see `xQsFactor`)
- **xQsFactor** (only used with `useFluctuatingX 1`): the factor $\beta$ in $x = \beta Q_s e^{\pm y}/\sqrt{s}$ ($+$ for the projectile, $-$ for the target), with $y$ = `rapidity`. $x$ and the local $Q_s$ are solved for self-consistently in every cell; for $x > 0.01$ the table's $Q_s^2$ at $y = 0$ is extrapolated with the factor $[(1-x)/0.99]^{5.6} (0.01/x)^{0.2}$
- **usePseudoRapidity**: `1`: `rapidity` is a pseudorapidity $\eta$, and the gluon multiplicity and hadron spectrum are given per unit $\eta$ instead of $y$, with the Jacobian for a particle of mass `jacobianMass` and transverse momentum $P = 0.13 + 0.32 (\sqrt{s}/\mathrm{TeV})^{0.115}$ GeV. With `useFluctuatingX 1`, $\eta$ is converted to $y$ with the same mass and momentum for the fluctuating $x$
- **jacobianMass** (only used with `usePseudoRapidity 1`): mass in GeV in the rapidity-pseudorapidity Jacobian

### Collision geometry
- **sigmaNN**: inelastic nucleon-nucleon cross section in mb, positive, used to decide which nucleons collide
- **gaussianWounding** (only used with `useNucleus 1`): how it is decided whether two nucleons collide
  - 0: hard sphere, they collide if their transverse distance $d$ is below $\sqrt{\sigma_{NN}/\pi}$
  - 1: they collide with the probability $p(d) = G e^{-G \pi d^2/\sigma_{NN}}$, $G = 0.92$, the Gaussian wounding profile of GLISSANDO (Eq. (13) of [arXiv:0710.5731](https://arxiv.org/abs/0710.5731)); its value of $G$ is taken from analyses of pp scattering at ISR energies (U. Amaldi and K. R. Schubert, Nucl. Phys. B 166 (1980) 301) and from A. Białas and A. Bzdak, Acta Phys. Polon. B 38 (2007) 159. Both profiles integrate to $\sigma_{NN}$
- **bMin** and **bMax** (only used with `useNucleus 1`): range of the impact parameter $b$ in fm, $0 \le$ `bMin` $\le$ `bMax`
- **sampleBFromLinearDistribution** (only used with `useNucleus 1`): `1` samples $b$ with a probability density $\propto b$ (uniform in the transverse plane, as for minimum-bias events), `0` uniformly in $b$
- **rotateReactionPlane** (only used with `useNucleus 1`): `1` points the impact parameter in a uniformly random direction (reaction-plane angle $\phi_{RP}$ in $[0, 2\pi)$), `0` along $x$
- **useFixedNpart** (only used with `useNucleus 1`): not negative; if not `0`, the impact parameter and reaction-plane angle are resampled, keeping the nucleon positions, until the event has exactly this number of participants; the value must be reachable for the sampled nuclei
- **minimumQs2ST**: trigger on high-multiplicity events: the impact parameter is resampled, keeping the nuclei, until $Q_{s,\min}^2 S_T$ exceeds this value (a non-negative number; `0` for no trigger). $Q_{s,\min}^2 S_T$ is the sum over all lattice cells of the smaller of the two nuclei's $Q_s^2$ times the cell area

### Nucleon positions
- **nucleonPositionsFromFile**: `1` takes each nucleus' nucleon positions from a randomly chosen configuration in the files in `nuclearConfigurationsPath` (see `lightNucleusOption`) and assigns the protons randomly; species without a file are sampled as with `0`. `0` samples them: a proton is a single nucleon, a deuteron is sampled from the Hulthén wave function, and heavier nuclei from the density profile of their species (see [Nuclear density profiles](#nuclear-density-profiles))
- **nuclearConfigurationsPath** (optional, default `./nucleusConfigurations`, only used with `useNucleus 1`, `readInitialWilsonLines 0`, and `nucleonPositionsFromFile 1` or a nonzero polarization): directory of the configuration files; running `download_nucleusTables.sh` inside `nucleusConfigurations/` downloads them there
- **lightNucleusOption** (only used with `useNucleus 1`, `readInitialWilsonLines 0`, and `nucleonPositionsFromFile 1` or a nonzero polarization): which configuration file `nucleonPositionsFromFile 1` uses. The file is chosen by the mass number $A$:
  - $A = 3$: `0` ³He, `1` triton
  - $A = 12$ (C): `0` variational Monte Carlo (VMC), `1` alpha clusters
  - $A = 16$ (O): `0` VMC, `1` alpha clusters, `2`/`3` clustered/uniform PGCM, `4`/`5` NLEFT with positive/negative weights
  - $A = 20$ (Ne): `0` or `2` clustered PGCM, `3` uniform PGCM, `4`/`5` NLEFT with positive/negative weights
  - $A = 40$: `0` VMC, `4` NLEFT (Ar files, also used for Ca)
  - d, ⁴He, Ne22, Au and Pb have one file each (for the deuteron see `polarizationProjectileJz`) and ignore this parameter. For the species above, a value that is not listed is rejected when reading the input
- **polarizationProjectile** and **polarizationTarget**: orientation of the nucleus, applied after its nucleon positions are sampled or read
  - 0: random orientation
  - 1: longitudinal: the nucleus' $z$ axis (the symmetry axis of a deformed Woods-Saxon nucleus) stays along the beam, with a random rotation about it
  - 2: transverse: the $z$ axis is turned to the $+y$ direction
- **polarizationProjectileJz** and **polarizationTargetJz** (only used with `useNucleus 1`, `readInitialWilsonLines 0`, and a nonzero `polarizationProjectile` or `polarizationTarget`, respectively): for a deuteron read from file (`nucleonPositionsFromFile 1`), $|J_z| = 1$ (any value within $10^{-8}$ of 1) uses the $J_z = \pm 1$ configurations and any other value the $J_z = 0$ ones. With polarization `0` they are not used: each event takes $J_z = 0$ with probability 1/3 and $J_z = \pm 1$ otherwise
- **useSmoothNucleus**: test option: `1` replaces the nucleons by the smooth thickness of each nucleus' density profile (without deformation, nucleons or hot spots); $N_\text{part}$ and $N_\text{coll}$ are set to 2, and the `NpartList`/`NcollList` files are not written
- **nucleiToAverage**: number $n$ of nuclei to average over, at least 1, `1` for normal runs. With $n > 1$, $n$ projectile and $n$ target nuclei are sampled (or read from file) and oriented independently, and every nucleon's thickness is divided by $n$, which gives a smoother thickness. $N_\text{part}$, $N_\text{coll}$ and the `NpartList`/`NcollList` files then count the nucleons of all $n$ nuclei, colliding every projectile with every target nucleon, so they are roughly $n N_\text{part}$ and $n^2 N_\text{coll}$; `useFixedNpart` refers to these totals. Not allowed for collisions with a proton

### Nuclear density profiles
Without configuration files, the nucleons of a nucleus with $A > 2$ are sampled from the radial density profile of its species, the profile the nuclear thickness integrals use too. `src/Glauber.cpp` has the built-in profile and parameters of each species:

| Profile | $\rho(r) \propto$ | Species |
|---|---|---|
| 3-parameter Fermi, a Woods-Saxon for $w = 0$ | $(1 + w r^2/R^2) / (1 + e^{(r-R)/a})$ | O and Ca ($w \neq 0$); Ne, Ne22, Al, Ar, Fe, Cu, Ru, Zr, Xe, W, Pt, Au, Pb, U ($w = 0$) |
| 3-parameter Gauss | $(1 + w r^2/R^2) / (1 + e^{(r^2-R^2)/a^2})$ | S |
| harmonic oscillator | $(1 + w r^2/a^2) e^{-r^2/a^2}$ | C |

He3 and He4 have no built-in profile, so they need configuration files (`nucleonPositionsFromFile 1`) or Woods-Saxon parameters from the input (`useInputWSParams 1`). A deformed nucleus is sampled from the 3-parameter Fermi profile with the angle-dependent radius $R(\theta, \phi) = R [1 + \beta_2 (\cos(\gamma) Y_{20} + \sin(\gamma) Y_{22}) + \beta_3 Y_{30} + \beta_4 Y_{40}]$. These species are deformed by default (all with $\beta_3 = \gamma = 0$):

| Species | $\beta_2$ | $\beta_4$ |
|---|---|---|
| O | −0.01 | −0.122 |
| Ar | 0.1668 | 0.00695 |
| Cu | 0.162 | 0.006 |
| Ru | 0.158 | 0 |
| Xe | 0.162 | −0.003 |
| Au | −0.13 | −0.03 |
| U | 0.28 | 0.093 |

The deformations of O and Xe are from FRDM(2012) ([arXiv:1508.06294](https://arxiv.org/abs/1508.06294)); the other parameters have no recorded source.

- **useInputWSParams**: `1` makes both nuclei Woods-Saxon nuclei (3-parameter Fermi with $w = 0$, whatever their built-in profile) with the input values of `radiusWS`, `diffusenessWS`, `beta2`, `beta3`, `beta4`, `gamma`, `deltaRnp` and `deltaAnp`, which are only used with `1`, instead of the built-in ones. A nucleus whose deformation parameters are all 0 (within $10^{-8}$) is sampled as spherical
- **radiusWS**: radius $R$ in fm, positive
- **diffusenessWS**: diffuseness $a$ in fm, positive
- **beta2**, **beta3** and **beta4**: deformation parameters $\beta_2$, $\beta_3$ and $\beta_4$
- **gamma**: triaxiality angle $\gamma$ in radians
- **deltaRnp** and **deltaAnp**: neutron skin: neutrons are sampled with radius $R$ + `deltaRnp` and diffuseness $a$ + `deltaAnp` (in fm); both are 0 without `useInputWSParams 1`
- **forceDMin**: `1` enforces `dMin` strictly for deformed nuclei: a nucleon's position is redrawn until it is at least `dMin` away from all others. `0` keeps `dMin` on a best-effort basis (see `dMin`)
- **dMin**: minimum distance in fm between two nucleons of a nucleus sampled from its density profile, also used with `useInputWSParams 0`. Without `forceDMin` it is kept on a best-effort basis: a nucleon keeps its sampled radius and only its direction is redrawn, up to 100 times; triaxial nuclei ($\gamma \neq 0$) then ignore it, since their density depends on every angle (use `forceDMin 1` for them). `0` for no minimum distance

### Coupling
- **g**: coupling constant $g$ of the classical Yang-Mills fields, positive; with fixed coupling, $\alpha_s = g^2/(4\pi)$
- **runningCoupling**: `0` for the fixed coupling. `1` rescales the energy density and $T^{\mu\nu}$ output, the gluon multiplicity and the eccentricities by $g^2/(4\pi\alpha_s(Q))$, with $\alpha_s(Q) = \frac{4\pi}{\beta_0 c \ln[(\mu_0/\Lambda_\mathrm{QCD})^{2/c} + (Q/\Lambda_\mathrm{QCD})^{2/c}]}$, $\beta_0 = (11 N_c - 2 N_f)/3$ and $Q$ = `runningCouplingQsFactor` times the $Q_s$ chosen by `runWithQs` and `runWithLocalQs`. All of them use the same factor in each cell, except that the gluon multiplicity uses $k_T$ with `runWithKt 1`. With `useJIMWLK 1`, $Q_s$ (averaged or local) is that of the initial condition, computed from the color-charge densities before the JIMWLK evolution, not of the evolved Wilson lines. Needs `useNucleus 1` and `LambdaQCD` < `mu0`
- **mu0** (only used with `runningCoupling 1`): $\mu_0$ in GeV, keeps $\alpha_s$ finite for $Q \to 0$
- **c** (only used with `runningCoupling 1`): how sharply $\alpha_s$ changes over from the $\mu_0$ to the $Q$ regime; must be positive
- **LambdaQCD** (optional, default `0.2`, only used with `runningCoupling 1`): $\Lambda_\mathrm{QCD}$ in GeV
- **nFlavors** (optional, default `3`): number of quark flavors $N_f$ in $\beta_0$, between 0 and 16; also used by the JIMWLK running coupling
- **runWithQs** (only used with `runningCoupling 1`): which of the two nuclei's $Q_s$ sets the scale: `0` the smaller, `1` the average, `2` the larger, averaged ($\sqrt{\langle Q_s^2 \rangle}$) over the overlap region, or per cell with `runWithLocalQs 1`
- **runningCouplingQsFactor** (only used with `runningCoupling 1`): factor multiplying $Q_s$ (or $k_T$ with `runWithKt 1`) to give the scale $Q$
- **runWithLocalQs** (only used with `runningCoupling 1`): `1` uses the local $Q_s$ of each cell instead of the overlap-region average (for the energy density and $T^{\mu\nu}$ output, the gluon multiplicity and the eccentricities)
- **runWithKt** (only used with `runningCoupling 1`): `1` computes the gluon multiplicity with $\alpha_s(Q)$, $Q$ = `runningCouplingQsFactor` $k_T$, in each $k_T$ bin instead of with the $Q_s$-based coupling; the energy density and $T^{\mu\nu}$ output and the eccentricities keep the $Q_s$-based coupling


### Output
The files themselves (names, order, layout, columns and units) are described in [OUTPUT.md](OUTPUT.md). Every output has its own switch, except `usedParameters<event>.dat`, which is always written (`NpartdNdy-t*` and `gluonMultiplicity*.json` share `computeGluonMultiplicity`). The hydro, Jazma and $T^{\mu\nu}$ files, the multiplicity, the eccentricities and the collision-geometry files (`NpartList`, `NcollList`, `NgluonEstimators`) are only written with `runEvolution 1`, the collision-geometry files and the Wilson-line geometry files only with `useNucleus 1`.

 - **writeHydro**: `1` writes the initial condition for hydrodynamic simulations, $\epsilon$, $u^\mu$ and $\pi^{\mu\nu}$ (`epsilon-u-Hydro-*.dat`), at the final time and at `outputTimes`
 - **writeJazma**: `1` writes the energy density of the Jazma model (`Jazma-Hydro-*.dat`), at the final time and at `outputTimes`
 - **writeTmunu**: `1` writes $T^{\mu\nu}$, e.g. for effective kinetic theory (KoMPoST) simulations (`Tmunu-*`), at the final time and at `outputTimes`
 - **writeTmunuBinary** (optional, default `1`, only used with `writeTmunu 1`): $T^{\mu\nu}$ in binary (`.ipgt`, `1`) or text (`.dat`, `0`) format
 - **outputTimes** (optional, default `none`, only used with `writeHydro`, `writeJazma` or `writeTmunu` 1): comma-separated proper times in fm/c (no spaces), e.g. `0.1,0.2,0.3,0.4`, at which these files are also written before the final time; `none` for only the final time. Each time is rounded down to a time step and must be smaller than `maxTime`; with `inverseQsForMaxTime 1`, times that are not before the final time are skipped
 - **sizeOutput**, **LOutput** (only used with `writeHydro`, `writeJazma` or `writeTmunu` 1): number of grid points per direction and side length [fm] of the transverse output grid the fields are interpolated to, both positive
 - **computeGluonMultiplicity**: at the final time, measure the gluon spectrum and multiplicity (files `NpartdNdy-t*` and `gluonMultiplicity*.json`)
 - **writeHadronSpectrum** (optional, default `0`, only used with `computeGluonMultiplicity 1`): `1` also writes the hadron spectrum from fragmenting the gluon spectrum (`multiplicityHadrons<id>.dat`)
 - **computeEccentricities**: at the final time, compute the eccentricities of the energy density (`eccentricities<id>.dat`)
 - **eccentricityCutoff** (optional, default `0`, only used with `computeEccentricities 1`): energy density in GeV/fm$^3$ below which a cell is left out of the eccentricities
 - **writeNpartList**, **writeNcollList** and **writeNgluonEstimators** (optional, default `1`): write the positions of the nucleons and whether they collided (`NpartList<id>.dat`), the binary collisions (`NcollList<id>.dat`) and the gluon-number estimators (`NgluonEstimators<id>.dat`)
 - **writeWilsonLineSnapshot** (optional, default `0`): `1` writes the initial Wilson lines of both nuclei as one binary file with a JSON header (`initialWilsonLines<id>.ipgw`) when they are built, with any `runEvolution` (not with `readInitialWilsonLines` 1 or 2)
 - **writeOutputsToHDF5**: this parameter decides whether to collect output files into an HDF5 file
   - 0: no
   - 1: yes; after each event, its `usedParameters`, `NpartList`, `NcollList` and `NpartdNdy-t*` files, the hydro files at `outputTimes` and the text $T^{\mu\nu}$ files are moved into `RESULTS_rank<rank>.h5` (the originals are deleted), and at the end of the run these are merged into `RESULTS.h5`. The other files stay on disk, and so does a per-rank file that cannot be merged. Needs `python3` with `h5py` and `numpy` (it calls `utilities/combine_events_into_hdf5.py` of the source tree, so the run can start in any directory)
 - **writeWilsonLines**: controls whether the generated Wilson lines are saved to disk. The file names (see OUTPUT.md) depend on the seed, the event, the number of events of the run and the nucleus. The initial Wilson lines are saved (with `useJIMWLK 1` only with `jimwlkSaveSnapshots 1`), and those after the JIMWLK evolution. With `useNucleus 1`, each nucleus' geometry (nucleon positions, color-charge densities and thicknesses) goes to a file `WilsonLineGeometry_<n>` next to its Wilson lines, which `readInitialWilsonLines` needs (see `writeWilsonLineGeometry`).
   - 0: do not save Wilson lines
   - 1: save in text format
   - 2: save in binary format (faster I/O, smaller file size)
 - **wilsonLinePath** (optional, default `./`, only used with `writeWilsonLines` or `readInitialWilsonLines` 1 or 2): directory used both when writing Wilson lines (`writeWilsonLines` is 1 or 2) and when reading them back in (`readInitialWilsonLines` is 1 or 2). When Wilson lines are written, the directory must already exist, otherwise the run fails at startup.
 - **writeWilsonLineGeometry** (optional, default `1`, only used with `writeWilsonLines` 1 or 2): `1` writes, next to each nucleus' Wilson-line files, a geometry file `WilsonLineGeometry_<n>` with its nucleon positions, color-charge densities and thicknesses; it is about 12% of the size of a binary Wilson-line file of the same lattice. Reading the Wilson lines back with `readInitialWilsonLines` and `useNucleus 1` needs them; runs that only produce Wilson lines for other codes, e.g. for vector-meson production, can set `0`
 - **readInitialWilsonLines**: `0` samples the color charges and builds the Wilson lines; `1` (text) or `2` (binary) instead reads the Wilson lines of both nuclei from files in `wilsonLinePath`, which an earlier run wrote with `writeWilsonLines`. The file names contain a number made from the `seed`, the number of events of the run and the event (see OUTPUT.md), so the reading run must use the same `seed`, events per rank and MPI ranks as the writing run to find the files of its events; their x is the one of `readWilsonLinesX`. With `useNucleus 1` it also reads their geometry files, so the collision is treated like one of sampled nuclei: the impact parameter is sampled with the same collision criterion, and $N_\text{part}$, $N_\text{coll}$, $\langle Q_s \rangle$ and all outputs are computed as usual. The geometry files must match the run's `size`, `L`, `g`, `projectile` and `target`; `QsMuRatio` is taken from them, since the color-charge densities were built with it. A run that reads the Wilson lines with the parameters of the run that wrote them, and with a fixed impact parameter and reaction plane and `gaussianWounding 0`, reproduces that run's outputs
 - **readWilsonLinesX** (optional, default `0`, only used with `readInitialWilsonLines` 1 or 2): Bjorken $x$ in the names of the Wilson-line files to read, e.g. the final $x$ of a JIMWLK evolution or a snapshot's $x$; both nuclei are read at this $x$. `0` reads the initial Wilson lines, at x = `jimwlkInitialX` with `useJIMWLK 1`, no x with `useFluctuatingX 1`, and `projectileX`/`targetX` otherwise. With `useJIMWLK 1`, the read Wilson lines are evolved from this $x$ (which must then not be smaller than `projectileX` and `targetX`)

### JIMWLK evolution
The JIMWLK evolution starts from a fixed $x$, `jimwlkInitialX`, and evolves the projectile to `projectileX` and the target to `targetX`; it requires `useFluctuatingX 0`.

- **useJIMWLK**: with JIMWLK (1), or no JIMWLK (0)
- **jimwlkInitialX** (only used with `useJIMWLK 1`): Bjorken-$x$ at the initial condition, at which $Q_s^2$ of both nuclei is read from the nuclear $Q_s$ table ($y = \ln(0.01/x)$, so it must not be larger than 0.01)
- **jimwlkSaveSnapshots** (only used with `useJIMWLK 1`): `1` also writes the Wilson lines during the evolution, at the $x$ values in `jimwlkXSnapshotList`, in the format of `writeWilsonLines` (which must be 1 or 2); the initial Wilson lines are then written too (unless they were read from file)
- **jimwlkXSnapshotList**: comma-separated $x$ values (no spaces), in any order, to save snapshots at, only used with `jimwlkSaveSnapshots 1`. Each is saved at the evolution step whose $x$ is closest to it in $\ln x$ (at most half a step away), in a file named with the requested value. A value more than half a step outside a nucleus' evolution (above its starting $x$, `jimwlkInitialX` or the `readWilsonLinesX` of read Wilson lines, or below its final $x$) is skipped with a warning
- **jimwlkMass** (only used with `useJIMWLK 1`): Infrared regulator in GeV in the JIMWLK kernel, see (21) in [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712)
- **jimwlkAlphaS** (only used with `useJIMWLK 1`): Coupling constant in the JIMWLK evolution
  - 0 (any value within $10^{-8}$ of 0): Use running coupling
- **jimwlkLambdaQCD** (only used with `useJIMWLK 1`): $\Lambda_\mathrm{QCD}$ in $\alpha_s(r)$ in GeV as in Eq. (22) of [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712)
- **jimwlkMu0** (only used with `useJIMWLK 1`): Regulator in $\alpha_s(r)$ as in Eq. (22) of [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712)
- **jimwlkC** (optional, default `0.2`, only used with `useJIMWLK 1`): parameter $c$ in $\alpha_s(r)$ as in Eq. (22) of [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712); $N_f$ there is `nFlavors`
- **jimwlkDs** (only used with `useJIMWLK 1`): step size in JIMWLK evolution. Recommended values
  - 0.005 with running coupling
  - 0.0005 with fixed coupling

  Default parameters for the JIMWLK evolution with fluctuating proton at initial $x=0.01$ fitted to HERA vector meson production data are reported in [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712)

## Utilities
Python scripts in `utilities/` for working with the output files (see [OUTPUT.md](OUTPUT.md)):
- `read_tmunu.py`: reads the $T^{\mu\nu}$ files, binary (`.ipgt`) and text (`.dat`), into numpy arrays, e.g. the energy density with `get_energy_density()`.
- `eccentricity.py`: the eccentricities $\epsilon_n$ of an energy-density grid, e.g. from `read_tmunu.py` (`eccentricity_from_tmunu_file()`).
- `combine_events_into_hdf5.py`: collects the text outputs of events into an HDF5 file. The program calls it after each event with `writeOutputsToHDF5 1`; it can also be run by hand on a results folder (`--help` lists the options).
- `fetch_IPGlasma_event_from_hdf5_database.py`: writes the hydro file (the MUSIC input) or the text $T^{\mu\nu}$ file of one event at a given proper time back out of such an HDF5 file: `fetch_IPGlasma_event_from_hdf5_database.py RESULTS.h5 <event id> <0: hydro, 1: Tmunu> <tau>`. Without the proper time it lists the available ones.
- `saveToBinaryFile.py`: an example of converting a text table of nucleon configurations into the binary format of the configuration files.

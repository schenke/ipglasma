# README

IP-Glasma initial condition with JIMWLK evolution

References
* Original IP-Glasma: [Schenke, Tribedy, Venugopalan, PRL 108 (2012) 252301](https://doi.org/10.1103/PhysRevLett.108.252301), [arXiv:1202.6646](https://arxiv.org/abs/1202.6646) and [Schenke, Tribedy, Venugopalan, PRC 86 (2012) 034908](https://doi.org/10.1103/PhysRevC.86.034908), [arXiv:1206.6805](https://arxiv.org/abs/1206.6805)
* JIMWLK evolution implementation: Mäntysaari, Schenke, Shen, Zhao, [Phys.Rev.Lett. 135 (2025) 2, 022302](https://doi.org/10.1103/gf4y-p5j7), [arXiv:2502.05138](https://arxiv.org/abs/2502.05138) and [Phys.Rev.C 113 (2026) 3, 034914](https://doi.org/10.1103/qsd7-zwmw), [arXiv:2511.03588](https://arxiv.org/abs/2511.03588)
* Standalone version of the JIMWLK code: https://github.com/hejajama/jimwlk



## Compile
To compile IP-Glasma, run `./compile_IPGlasma.sh`
Dependencies
* CMake
* FFTW
* GSL

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


## Input parameters
The input file is given as the first command line argument (default: `input`). Every parameter is listed once in `src/ParameterTable.cpp`, together with its default value (if it is optional) and its validity checks; see `src/Parameters.h` for a more detailed description of each parameter.

The input file has one `key value` pair per line:
- `#` starts a comment that runs to the end of the line, and blank lines are ignored.
- A line with only `EndOfFile` ends the input; everything after it is ignored. It is optional.
- Unknown keys, keys given twice, missing required keys and malformed values (e.g. `256.0` for an integer) are errors. All problems are reported at once before the run stops.

Each event writes the values of all input parameters it used to `usedParameters<event>.dat`, followed by its random seed and collision geometry as comments. The file is itself a valid input file. Running it does not reproduce the same event, though: the random numbers also depend on the MPI rank and on the event's position in the run, and with `subNucleonParamSet -1` a new posterior parameter set is drawn.

### Lattice
- **size**: controls the size of the lattice that is `size`$^2$.
  - Recommended to be of the form $2^n$
- **L**: the total physical extent of the lattice (in fm)

### Initial state
- **useNucleus**
  - 1: nucleus with finite geometry
  - 0: infinite target with constant color charge density controlled by `g2mu` (in lattice units)
- **projectile** and **target**: specify nuclei
  - Typical values: `p`, `Pb`, `Au`
  - See `src/Glauber.cpp` for all supported nuclei and details
- **m**: infrared regulator in GeV
- **BG**: nucleon width in GeV$^{-2}$. With `nucleonModel gaussian` the nucleon's thickness is $T \sim e^{-b^2/(2B_G)}$; with `hotspots` and `strings` it is the width of the hot-spot position distribution (see `omega`)
- **nucleonModel**: transverse structure of a nucleon
  - `gaussian`: a single Gaussian of width `BG` (optionally elongated by `protonAnisotropy`)
  - `hotspots`: `Nq` Gaussian hot spots of width `BGq`, whose positions are distributed with width `BG`
  - `strings`: three hot spots, sampled as for `hotspots`, connected by strings that meet at a junction, the point with the shortest total string length (Fermat point, computed in 3D). Each hot spot is moved to a uniformly random point on its string, between the junction and its sampled position, so the nucleon is smaller than with `hotspots` for the same `BG`. The minimum distance `dqMin` applies to the sampled hot spots (the ends of the strings), not to the moved ones. Uses the hot-spot parameters except `Nq` and `NqFluc`.

  The parameters of the other models are not read.
- **subNucleonParamType**: `0` uses the input values for the nucleon substructure; `1`, `2` or `4` draw them every event from a Bayesian posterior parameter set (`tables/posterior.csv` with variable Nq, `tables/posterior_Nq3.csv` or `tables/posterior5020_Nq3.csv` with Nq = 3), selected by **subNucleonParamSet** (`-1`: random). The posterior sets are fits of hot-spot nucleons, so they require `nucleonModel hotspots`. They replace `m`, `BG`, `BGq`, `smearingWidth`, `QsMuRatio`, `dqMin` and the number of hot spots, which are then not read from the input; `usedParameters<event>.dat` lists the values an event used.
- **protonAnisotropy** (`gaussian`): anisotropy $\xi > -1$ of the nucleon, $T \propto \sqrt{1+\xi}\, e^{-(b^2 + \xi (\vec b \cdot \hat n)^2)/(2B_G)}$ with a random direction $\hat n$ per nucleon, so the nucleon is narrower along $\hat n$ for $\xi > 0$; `0` gives a round nucleon
- **BGq** (`hotspots`, `strings`): mean hot-spot width $B_q$ in GeV$^{-2}$, hot spot density profile is $T_q \sim e^{-b^2/(2B_q)}$. The hot spots of a nucleon share one width, 0.09 GeV$^{-2}$ plus a log-normal number with mean `BGq` $-$ 0.09 and variance `BGqVar`, so `BGq` must be larger than 0.09 GeV$^{-2}$
- **BGqVar** (`hotspots`, `strings`): variance of the hot-spot width in GeV$^{-4}$ (see `BGq`); `0` gives every nucleon the width `BGq`
- **dqMin** (`hotspots`, `strings`): minimum distance in fm between the hot spots of a nucleon, in 3D for `omega 1` and in the transverse plane otherwise. It is kept on a best-effort basis: a hot spot keeps its sampled radius and only its direction is redrawn, up to 100 times. With `strings` it applies to the ends of the strings
- **omega** (`hotspots`, `strings`): radial distribution of the hot spots
  - `1` (any value within $10^{-8}$ of 1): positions in 3D, each coordinate Gaussian with variance `BG` (transverse $\langle b^2 \rangle = 2 B_G$)
  - otherwise: positions in the transverse plane, at a uniformly random angle and radius $b = \sqrt{\omega x B_G}$, where $x$ follows the density $\propto Q(1/\omega, x)$ (regularized upper incomplete gamma function), so $\langle b^2 \rangle \approx (1+\omega) B_G/2$. Note that this does not approach the `omega 1` case for $\omega \to 1$
- **Nq** (`hotspots`): mean number of hot spots per nucleon, at least 1; a fractional value such as 2.5 gives 2 or 3 hot spots with the corresponding probabilities (plus the fluctuation set by `NqFluc`)
- **NqFluc** (`hotspots`): mean of a Poisson-distributed number of additional hot spots per nucleon; `0` for no fluctuation. Every nucleon has at least one hot spot
- **shiftConstituentQuarkProtonOrigin** (`hotspots`, `strings`): whether to shift the center-of-mass to origin (1) or not (0) after sampling the hot spot positions
- **smearQs**: enable (1) or disable (0) saturation scale fluctuations: the thickness of each nucleon (`gaussian`) or each hot spot (`hotspots`, `strings`) is multiplied by a log-normal factor $e^{X}/e^{\sigma^2/2}$ with mean 1, where $X$ is Gaussian with width $\sigma$ = `smearingWidth`
- **smearingWidth**: width $\sigma$ of the saturation scale fluctuations (see `smearQs`), parameter $\sigma$ in Eq. (23) of [arXiv:1607.01711](https://arxiv.org/pdf/1607.01711)
- **useFluctuatingX**: controls how to determine Bjorken-$x$ when generating the initial condition
  - 1: Dynamically determined $b_\perp$ dependent $x$
  - 0: Fixed $x$
- **rapidityA** and **rapidityB**:
  - If `useFluctuatingX 0`, then $x = 0.01 e^{-\mathrm{rapidityA}}$ for the projectile and $x = 0.01 e^{-\mathrm{rapidityB}}$ for the target; both must not be negative (the nuclear $Q_s$ table covers $0 \le y \le 10.75$, larger values use $Q_s$ at $y = 10.75$)
  - If `useFluctuatingX 1`, consider particle production at rapidity $y$


### Output
The files themselves (names, order, layout, columns and units) are described in [OUTPUT.md](OUTPUT.md). The parameters that switch them on:

 - **writeOutputs**: this parameter controls output files (written in `mode 1`; the values add up as bits)
   - 0: no output
   - 1: output initial conditions $\epsilon$, $u^\mu$, and $\pi^{\mu\nu}$ for hydrodynamic simulations (needs `writeEpsilonUHydro 1`)
   - 2: output the initial condition for energy density according to the Jazma model (needs `writeEpsilonUHydro 1`)
   - 3: output 1 & 2; with `computeGluonMultiplicity 1` also the hadron spectrum
   - 4: output initial $T^{\mu\nu}$ for the effective kinetic theory (KoMPoST) simulations
   - 5: output 1 & 4; in addition the same outputs at $\tau \approx$ 0.1, 0.2, 0.3 and 0.4 fm/c, and a snapshot of the initial Wilson lines (`initialWilsonLines<id>.ipgw`)
   - 6: output 2 & 4
   - 7: output 1 & 2 & 4
 - **writeEpsilonUHydro** (optional, default `1`): write the hydro and Jazma files, which needs the Landau matching (flow velocity); with `0`, only $T^{\mu\nu}$ is written
 - **writeTmunuBinary** (optional, default `1`): $T^{\mu\nu}$ in binary (`.ipgt`, `1`) or text (`.dat`, `0`) format; the environment variable `IPGLASMA_BINARY_TMUNU` overrides it
 - **sizeOutput**, **LOutput**: number of grid points per direction and side length [fm] of the transverse output grid the fields are interpolated to
 - **etaSizeOutput**, **dEtaOutput**: number of points and spacing of the (boost-invariant) $\eta$ grid in the hydro and Jazma files
 - **computeGluonMultiplicity**: at the final time, measure the gluon spectrum and multiplicity and the eccentricities (files `NpartdNdy-t*`, `gluonMultiplicity*.json`, `eccentricities*.dat`)
 - **eccentricityCutoff** (optional, default `0`): energy density in GeV/fm$^3$ below which a cell is left out of the eccentricities
 - **readMultFromFile**: post-processing mode that rescales the multiplicity of an earlier run (see OUTPUT.md) and then stops; `0` for normal runs
 - **writeOutputsToHDF5**: this parameter decides whether to collect all the IPGlasma output files into a hdf5 data file
   - 0: no
   - 1: yes; after each event its text files are moved into `RESULTS_rank<rank>.h5` (the originals are deleted)
 - **writeWilsonLines**: controls whether the generated Wilson lines are saved to disk. File names depend on random seed (parameter `seed`), see `WilsonLineIO::fileName()`. Wilson lines at the initial condition and after the JIMWLK evolution are saved.
   - 0: do not save Wilson lines
   - 1: save in text format
   - 2: save in binary format (faster I/O, smaller file size)
 - **wilsonLinePath** (optional): directory used both when writing Wilson lines (`writeWilsonLines` is 1 or 2) and when reading them back in (`readInitialWilsonLines` is 1 or 2). Defaults to `./`. The directory must already exist, otherwise the run fails at startup.

### JIMWLK evolution
Note that when using the JIMWLK evolution, one should use `useFluctuatingX 0` which corresponds to having a fixed $x$ at the initial state of the evolution.

- **useJIMWLK**: with JIMWLK (1), or no JIMWLK (0)
- **jimwlkInitialX**: Bjorken-x at the initial condition
- **jimwlkXProjectile**: Bjorken-$x$ to which the projectile is evolved
- **jimwlkXTarget**: Bjorken-$x$ to which the target is evolved
- **jimwlkMass**: Infrared regulator in GeV in the JIMWLK kernel, see (21) in [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712)
- **jimwlkAlphaS**: Coupling constant in the JIMWLK evolution
  - 0 (any value within $10^{-8}$ of 0): Use running coupling
- **jimwlkLambdaQCD** $\Lambda_\mathrm{QCD}$ in $\alpha_s(r)$ in GeV as in Eq. (22) of [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712)
- **jimwlkMu0**: Regulator in $\alpha_s(r)$ as in Eq. (22) of [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712)
- **jimwlkDs**: step size in JIMWLK evolution. Recommended values
  - 0.005 with running coupling
  - 0.0005 with fixed coupling

  Default parameters for the JIMWLK evolution with fluctuating proton at initial $x=0.01$ fitted to HERA vector meson production data are reported in [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712)
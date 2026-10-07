# Changelog

All notable changes to this project will be documented in this file.

The versioning of the codebase is inspired by [Semantic Versioning](https://semver.org/spec/v2.0.0.html) with a version number `X.Y.Z`, where

* `X` is incremented for major changes (e.g. large backwards incompatible updates),
* `Y` is incremented for minor changes (e.g. a new feature) and
* `Z` is incremented for the indication of a bug fix or other very small changes that are not backwards incompatible.

The main categories for changes in this file are:

* `Input / Output` for all, in particular breaking, changes, fixes and additions to the in- and output files;
* `Added` for new features;
* `Changed` for changes in existing functionality;
* `Fixed` for any bug fixes;
* `Removed` for now removed features.

## 2.0.0

The changes since b36b2a9, the last commit on the default branch (`main`, formerly `master`) before this release; see the 1.x section below. Input parameters are given by their new names; the first entries list the renames.

### Input / Output
* **Not backwards compatible:** the nucleon substructure is now selected by the new required parameter `nucleonModel` (`gaussian`, `hotspots` or `strings`) instead of by `useConstituentQuarkProton > 0`, and the number of hot spots is `Nq`. `Nq` (at least 1; a fractional value is the mean number of hot spots) and `NqFluc` are only read with `nucleonModel hotspots`, `BGq`, `BGqVar`, `dqMin`, `omega` and `shiftConstituentQuarkProtonOrigin` with `hotspots` or `strings`, and `protonAnisotropy` only with `gaussian`. The posterior parameter sets (`subNucleonParamType` 1, 2, 4) require `nucleonModel hotspots`; they replace `m`, `BG`, `BGq`, `smearingWidth`, `QsMuRatio`, `dqMin` and the number of hot spots every event, so those keys (and `Nq`) are then not read from the input, and `usedParameters<event>.dat` lists the values used as a comment.
* **Not backwards compatible:** rename input parameters to one naming convention (issue #32): lowerCamelCase, with physics symbols keeping their conventional spelling (`Qs`, `L`, `BG`, `sigmaNN`, `LambdaQCD`), no underscores, and a `jimwlk` prefix for the parameters of the new JIMWLK evolution. Renamed: `maxtime` → `maxTime`, `Projectile` → `projectile`, `Target` → `target`, `roots` → `sqrtS`, `SigmaNN` → `sigmaNN`, `bmin` → `bMin`, `bmax` → `bMax`, `samplebFromLinearDistribution` → `sampleBFromLinearDistribution`, `averageOverThisManyNuclei` → `nucleiToAverage`, `polariztionProjectile` → `polarizationProjectile`, `polariztionTarget` → `polarizationTarget`, `setWSDeformParams` → `useInputWSParams`, `R_WS` → `radiusWS`, `a_WS` → `diffusenessWS`, `dR_np` → `deltaRnp`, `da_np` → `deltaAnp`, `force_dmin_flag` → `forceDMin`, `d_min` → `dMin`, `SubNucleonParamType` → `subNucleonParamType`, `SubNucleonParamSet` → `subNucleonParamSet`, `UVdamp` → `UVDamp`, `QsmuRatio` → `QsMuRatio`, `NucleusQsTableFileName` → `nucleusQsTableFileName`, `Jacobianm` → `jacobianMass`, `useFluctuatingx` → `useFluctuatingX`, `xFromThisFactorTimesQs` → `xQsFactor`, `muZero` → `mu0`, `runWith0Min1Avg2MaxQs` → `runWithQs`, `runWithThisFactorTimesQs` → `runningCouplingQsFactor`, `runWithkt` → `runWithKt`, `detaOutput` → `dEtaOutput`, `writeInitialWilsonLines` → `writeWilsonLines`, `useTimeForSeed` → `useRandomSeed` (it draws the seed from `std::random_device`, not from the time). An input file that still uses an old name is rejected with a message naming the new one (e.g. `maxtime was renamed to maxTime`), and a removed or replaced key with a hint at what replaced it.
* **Not backwards compatible:** replace `Rapidity` by `projectileX` and `targetX`, the Bjorken x of the projectile and of the target, and `rapidity`, the rapidity of the gluon and hadron spectra. Without JIMWLK and with `useFluctuatingX 0`, Q_s² is evaluated at `projectileX` and `targetX` (the old `Rapidity` y corresponds to x = 0.01 e^{-y}; x must not be above 0.01); with `useJIMWLK 1`, which requires `useFluctuatingX 0` (both set to 1 is rejected), it is evaluated at `jimwlkInitialX`, and the nuclei are evolved to `projectileX` and `targetX` (`jimwlkXProjectile` and `jimwlkXTarget` in earlier development versions). With `useFluctuatingX 1` and no JIMWLK, `rapidity` is the y in x = `xQsFactor` Q_s e^{±y}/√s, as `Rapidity` was. The Wilson-line binary header and `initialWilsonLines<id>.ipgw` store `rapidity`.
* **Not backwards compatible:** remove the input parameters `Nc` (the code is SU(3) only), `rmax` (every shipped input file disabled its cutoff with `10000000.`; color charges are now always assigned over the whole lattice), `dtau` (it was ignored: the time step follows from `maxTime`, `L` and `size`), and `tDistNu`, `useFatTails` and `writeEvolution`, which had no effect.
* **Not backwards compatible:** replace `writeOutputs`, whose values 0-7 selected combinations of output files, by one switch per output file: `writeHydro`, `writeJazma` and `writeTmunu` for the field outputs, written at the final time and at the proper times of the new list `outputTimes` (optional, default `none`), which replaces the fixed 0.1, 0.2, 0.3 and 0.4 fm/c of `writeOutputs 5`; `writeHadronSpectrum` (was `writeOutputs 3`); `writeWilsonLineSnapshot` for `initialWilsonLines<id>.ipgw` (was `writeOutputs 5`); and `writeNpartList`, `writeNcollList` and `writeNgluonEstimators` (optional, default 1) for the files that used to be written always. Any combination is possible now, e.g. the hadron spectrum together with T^μν or the Jazma file at intermediate times. A list value can be `none`. T^μν is written in the binary `.ipgt` format by default; the new `writeTmunuBinary 0` (optional, default 1) writes the text `.dat` file as before. The eccentricities have their own switch, `computeEccentricities`, instead of coming with `computeGluonMultiplicity`, and `eccentricities<id>.dat` holds only the final time (it also had a row for the first time step).
* **Not backwards compatible:** the Wilson-line files are named `WilsonLine_x_<x>_<n>`, with the Bjorken x of the Wilson lines (`WilsonLine_<n>` with `useFluctuatingX 1`, where x varies across the nucleus) and `.txt` for text, instead of `V-<n>`, and are written to and read from the directory `wilsonLinePath` (optional, default `./`). They are numbered `2 (seed N + eventId) + iA`, with N the events per rank times the MPI ranks and iA 1 for the projectile and 2 for the target. The old numbers, `eventId + 2 seed nRanks` for the projectile and `eventId + (1 + 2 seed) nRanks` for the target, gave the target of one event the number of the projectile of another when a rank ran several events, so their files overwrote each other.
* **Not backwards compatible:** the Wilson-line files are now stored in the same site and element order as the in-memory lattice (site index `N*ix+iy`, ix outer; row-major 3x3 elements), for both binary (`writeWilsonLines 2`) and text (`writeWilsonLines 1`) output. Binary files used to be stored with x and y swapped, and text files with each matrix transposed (column-major). Files written by older versions are read transposed by this version, and external readers that compensated for the old layout (e.g. subnucleondiffraction, which transposes text Wilson lines and swaps x/y in binary ones) need to be updated. The text `Phi`/`Pi` matrices written by `Lattice::writeSU3Matrices` are now row-major as well.
* **Not backwards compatible:** reading Wilson lines (`readInitialWilsonLines 1/2`) with `useNucleus 1` needs the new geometry files `WilsonLineGeometry_<n>` (nucleon positions, color-charge densities and thicknesses of each nucleus), which are written next to the Wilson lines (about 12% of their size; `writeWilsonLineGeometry 0` switches them off for runs that never read the Wilson lines back). With them, read Wilson lines are collided like sampled nuclei (issue #52): the impact parameter and reaction plane are sampled with the same collision criterion, N_part, N_coll and ⟨Q_s⟩ are computed, and `runningCoupling 1` and `inverseQsForMaxTime 1` are possible. Before, the impact parameter was sampled without the collision criterion, the reader shifted the read fields by it, and N_part, N_coll and ⟨Q_s⟩ were 0. The new optional `readWilsonLinesX` reads the Wilson lines at a given x, e.g. after a JIMWLK evolution, which JIMWLK then continues from.
* **Not backwards compatible:** the input file is now read once and checked strictly. Unknown keys (e.g. typos, or parameters removed in this release), keys given twice, missing required keys and malformed values (e.g. `256.0` for an integer parameter, `0.4GeV` for a number) are errors, and all of them are reported together before the run stops. `#` comments are allowed (so values cannot contain `#`) and the `EndOfFile` line is optional.
* `usedParameters<event>.dat` now lists every input parameter used by the event (not a hand-picked subset) in input-file syntax, with the per-event information (random seed, `b`, `Npart`, ...) as `#` comments, so it is a valid input file with the same parameter values (it does not reproduce the same event, since the random numbers also depend on the MPI rank and the event's position in the run).
* Add `nuclearConfigurationsPath` input parameter (optional, default `./nucleusConfigurations`) to configure the directory nucleon configuration tables are read from.
* Add the Ne22 (NLEFT) nucleon configurations and the Ne22 species.
* Add `nFlavors` (optional, default `3`) and `LambdaQCD` (optional, default `0.2` GeV), which set the running coupling of the classical evolution and of the hydro output; they were hard-coded. `nFlavors` also sets the JIMWLK running coupling.
* Replace the default nuclear Q_s² table `qs2Adj_vs_Tp_vs_Y_200.in` by `qs2Adj_vs_Tp_vs_Y_240.in`, which covers a 10 times larger T_p range, so the "T out of range, using maximal T in table" clamping at a high local thickness is avoided.
* Remove the `nucleonPositionsFromFile 2` option (correlated Pb-208 configurations of Alvioli et al.), which read from a hard-coded, machine-specific path and was never packaged with the code; `nucleonPositionsFromFile` now only accepts `0` (sample) or `1` (read the configuration files).
* Remove the `epsilonInitialPlot<id>.dat` and `epsilonIntermediatePlot<id>.dat` energy-density maps of `writeOutputs 3`, written at the initial and the halfway time step.

### Added
* Add a JIMWLK small-x evolution of the projectile and target Wilson lines before the classical Yang-Mills evolution (`useJIMWLK 1`), with the parameters `jimwlkMu0`, `jimwlkLambdaQCD`, `jimwlkC`, `jimwlkMass`, `jimwlkAlphaS` (`0` for running coupling), `jimwlkDs` and `jimwlkInitialX` (evolving to `projectileX` and `targetX`); the running coupling uses `nFlavors`. With `jimwlkSaveSnapshots 1`, the Wilson lines are also written at the evolution step closest (in ln x) to each value of `jimwlkXSnapshotList`, named with that value; a value outside the evolution is skipped with a warning. Inputs that would give a NaN or negative coupling (e.g. a `jimwlkLambdaQCD` that is not positive or not below `jimwlkMu0`, a non-positive `jimwlkC`, a negative `jimwlkAlphaS`) are rejected.
* Add the string nucleon model, `nucleonModel strings`: three hot spots connected by strings that meet at their Fermat point (the junction, computed in 3D), each hot spot moved to a uniformly random point on its string. It corresponds to `Use_stringy_proton 2` of the stringy-proton branch, but samples the junction and the positions on the strings once per nucleon instead of anew for every lattice cell, and computes the junction in closed form instead of iteratively, so no event is rejected. The junction shift (`GeoM_shift`) is not included.
* Add the optional input parameter `eccentricityCutoff` (GeV/fm³, default `0`): cells with a lower energy density are left out of the eccentricities. It replaces the hard-coded cutoff of 0, which was described in units of Λ_QCD⁴.
* Add the gluon spectrum output `gluonMultiplicity<id>.json` (with `mode 1` and `computeGluonMultiplicity 1`): the gluon number and energy, ⟨k_T⟩, the yields above 3 and 6 GeV, and the event's N_part, T_pp, impact parameter and seed.
* Add `OUTPUT.md`, which describes every output file (when it is written, in which order, its name, layout, columns and units); it is also part of the Doxygen documentation, and the functions writing the files point to it. Document the output parameters that were missing from the README (`writeTmunuBinary`, `computeGluonMultiplicity`, the output grid), the hadron spectrum and the outputs at intermediate times. New tests check the documented layouts.
* Add the Doxygen documentation of all classes and functions (CMake target `doc`, with the README and `OUTPUT.md` as pages and the formulas rendered by MathJax); the `undocumented_test` target fails if anything is undocumented.
* Add the CMake option `-DIPGLASMA_DETERMINISTIC_FFT=ON` (off by default), which keeps the FFTW plans in a wisdom file, so repeated runs with the same seed give identical results; `FFTW_MEASURE` can otherwise pick a different plan, and so a different rounding, in every run. Flags passed with `-DCMAKE_CXX_FLAGS` are no longer overwritten by `CMakeLists.txt`.
* Add profiling and fingerprint output of the initialization and evolution, switched on with the environment variables `IPGLASMA_PROFILE` and `IPGLASMA_FINGERPRINT` (see `OUTPUT.md`).
* Add the validation scripts in `validations/`: the coherent and incoherent J/ψ cross section from IP-Glasma, JIMWLK and subnucleondiffraction, compared to arXiv:2207.03712, and a regression check of the ε₂ distribution.
* Add a `doctest` unit test suite under `tests/` (built with `-Dunittest=ON`), which covers the matrix algebra, random numbers, input reading, nucleus and nucleon sampling, collision geometry, Wilson-line files, forward light cone, energy-momentum tensor, multiplicity and eccentricities.
* Add GitHub Actions workflows that build the code and run the unit tests on Ubuntu and macOS, with and without MPI (`build.yml`, `tests.yml`), check the Doxygen documentation (`doxygen.yml`), and apply `clang-format` to pull requests from branches of this repository (`clang-format.yml`).

### Changed
* Rewrite `Matrix`, `Cell` and `Lattice` as a structure-of-arrays layout with a fixed 3x3 `Matrix`, generate the random numbers in bulk, and batch and parallelize the FFTs, which makes the code substantially faster.
* Split `Init` and `Evolution` into classes with one task each: `ForwardLightCone`, `WilsonLineIO`, `NucleonModel` (a new nucleon model only needs a new model and profile class), `NuclearQsTable`, `NucleusSampler`, `CollisionGeometry`, `EnergyMomentumTensor`, `GluonMultiplicity`, `RunningCoupling` and `Eccentricity`; define the lattice site layout once (`LatticeIndex.h`); and merge duplicated code, e.g. the projectile and target branches, the running-coupling formula, the KKP fragmentation fits and the nucleus parameters (now tables).
* Read the input file once (`InputFile`) and define each parameter's key, default and checks in one table (`src/ParameterTable.cpp`), which drives reading, validation and `usedParameters<event>.dat`. `Parameters` holds the values as public fields grouped by topic instead of ~120 getter/setter pairs. Checks of a single value (e.g. an even `size`, `omega > 0`) apply even when the feature using it is off, and on/off parameters only accept `0` or `1`.
* Use one naming convention for files, classes, functions and variables, one include-guard style, and format all code with `clang-format` (`formatCode.sh`). Renamed files include `pretty_ostream` → `PrettyOstream` and `Phys_consts.h` → `PhysConst.h`.
* Print all messages through `PrettyOstream`, with a `[Class::function]` tag and colors for warnings (orange) and errors (red).
* Compute the nuclear thickness and normalization integrals with GSL's adaptive CQUAD quadrature (relative tolerance 1e-6) instead of a hand-written recursive Newton-Cotes rule (`Glauber::qnc7`). This corrects the thickness of the 3-parameter Gauss profile (e.g. S) by up to 7.5e-4 relative, an error of the old rule's convergence criterion. A failed integral gives a warning instead of being able to abort the run.
* Replace raw `new[]`/`malloc` arrays by `std::vector`.
* Build the nuclei and their Wilson lines centered at the origin and shift them by the impact parameter only afterwards, so the Wilson lines can be evolved with JIMWLK and written or read independently of b. The random numbers are drawn in a different order, so a given seed gives different events than 1.x.
* Average over nuclei (`nucleiToAverage` n > 1) with n independently sampled and oriented projectile and target nuclei, which also works with `nucleonPositionsFromFile 1` (n configurations). It used to sample one nucleus with n·A nucleons from the same Woods-Saxon distribution, and with configuration files only divided every nucleon's thickness by n. A collision with a proton is rejected when reading the input.
* Compare input values that select a code path (`omega 1`, `jimwlkAlphaS 0`, a deuteron's `|Jz| = 1`, vanishing deformation parameters) with one tolerance, 1e-8, instead of separate hard-coded values (1e-8, 1e-10, 1e-15).
* Deform U and Xe by default, with β2 = 0.28 and 0.162 (their β4 = 0.093 and −0.003 were built in already); their β2 used to come from the input `beta2`, which is only read with `useInputWSParams 1`, so it was uninitialized otherwise.
* Change `QsMuRatio` in the example input file from 0.8 to 0.643, following arXiv:2207.03712.
* Rename the GitHub repository's default branch from `master` to `main`.

### Fixed
* With `omega` different from 1, the hot-spot radius is now $b=\sqrt{2\omega x B_G}$ instead of $\sqrt{\omega x B_G}$ (issue #59), so $\langle b^2\rangle = (1+\omega)B_G$ and the transverse hot-spot distribution approaches that of `omega 1` for $\omega\to1$; the old radii were too small by $\sqrt2$ and $\langle b^2\rangle$ halved when $\omega$ moved away from 1. The table $x$ is drawn from now extends to $x=\max(20, 20/\omega)$ instead of $\max(5, 5/\omega)$, whose cut lowered $\langle b^2\rangle$ by another ~4% near $\omega=1$.
* Fix a polarized projectile with an unpolarized target having no effect with `nucleonPositionsFromFile 0`: polarized nuclei only exist as configuration files, but the check that switches to them looked at the target polarization twice.
* Fix a missing nucleon-configuration file making the run loop forever instead of stopping with an error. The files are now only read when they are used, so Woods-Saxon runs (`nucleonPositionsFromFile 0`) of species that have configuration files (e.g. Pb, Au, O) no longer need `nucleusConfigurations/`.
* Fix `forceDMin` and `dMin` being read only with `useInputWSParams 1` although they are always used: otherwise they were uninitialized, so the nucleon sampling of deformed species (e.g. Au, Cu) depended on garbage memory and was not reproducible for a fixed seed.
* Fix the Au197 and Pb208 nucleon configurations (`nucleonPositionsFromFile 1`): the files store 4 numbers per nucleon, but were read as if they stored 3, so all nucleons after the first were wrong.
* Fix asymmetric collisions of smooth nuclei (`useSmoothNucleus 1`), whose projectile and target thicknesses were swapped.
* Fix polarized deuterons with `Jz = 0` loading the `Jz = +-1` configurations.
* Fix the deuteron (Hulthen) thickness function missing the Jacobian of its integration variable, which made `T(s)` 17-44% too low and too compact (it integrated to 1.45 instead of A = 2).
* Fix `useConstituentQuarkProton` values between 0 and 1 switching the hot spots off, since the on/off flag was the truncated value; `Nq` must now be at least 1, and a fractional value sets the mean number of hot spots.
* Fix the event-averaged running coupling (`runningCoupling 1`), which decides whether an event is accepted and is written to `usedParameters<event>.dat`, using a different formula than the evolution (3 flavors, Λ_QCD = 0.2 GeV, no `mu0` regulator): it could be negative or singular at a small ⟨Q_s⟩ and reject valid events. It now uses the formula of the evolution, so with running coupling other events can be accepted than in 1.x (issue #53).
* Fix `runningCoupling 1` with `runWithKt 1` rejecting every event (alpha_s is evaluated per k_T bin there, but the event-acceptance cut required the unset global alpha_s to be positive), so the run never finished.
* Fix the fixed coupling alpha_s in `usedParameters<event>.dat`, which was written before it was set.
* Fix the eccentricity computation using an unset Qs/mu ratio for nucleus B when Wilson lines are read from file.
* Fix out-of-bounds reads of the nuclear Qs table: a rapidity between 10.75 (the largest tabulated value) and 11 read past the table, and so did a T_p equal to the largest tabulated one and a negative rapidity with a fixed x. A rapidity above the table now uses Q_s at y = 10.75 (with a warning), like T_p above the table, so fluctuating-x events with dilute cells at forward rapidity, which used to stop above y = 11, now run; an x above 0.01 (a negative rapidity) with `useFluctuatingX 0` is rejected when reading the input.
* Fix the overlap-region scan of the collision geometry skipping the corner cell (0, 0) and, for impact parameters reaching the lattice edge, reading cells on the opposite edge (an out-of-range coordinate aliased a valid site).
* Fix reading Wilson lines in text format (`readInitialWilsonLines 1`) when the lattice spacing `L/size` is not a binary fraction (e.g. `L 21`, `size 60`): the column index was truncated instead of rounded, so a few percent of the columns were stored one column too far left. Both formats now place the columns the same way.
* Fix the matrix exponential's test for its NaN fallback, a vanishing real part of the identity coefficient, which a valid SU(3) matrix can also have; it now tests that all coefficients vanish.
* Fix a NaN time step when `maxTime` is shorter than one step (e.g. `maxTime 0` to only produce Wilson lines); `dtau` now falls back to 0.1 and no evolution step is run.
* Reject inputs that used to fail silently, give NaN or undefined behaviour: an odd `size` (the FFTs assume an even size), `L <= 0`, `omega <= 0`, a negative `maxTime` or one needing more than 1e8 time steps, `nucleiToAverage < 1`, polarization values other than 0/1/2, a negative `Nq`, `subNucleonParamSet < -1`, an unknown `subNucleonParamType`, `runWithQs` outside 0-2, `runWithKt` other than 0/1, `protonAnisotropy <= -1`, `BG <= 0`, `BGq` up to 0.09 GeV⁻² (the hot-spot width is 0.09 plus a log-normal number with mean `BGq` − 0.09, so it was wrong or NaN), negative `BGqVar`, `dqMin` and `smearingWidth`, and, with running coupling, `LambdaQCD <= 0`, `c <= 0`, `LambdaQCD >= mu0` and an `nFlavors` that makes the one-loop β-function coefficient non-positive. `inverseQsForMaxTime 1` and `runningCoupling 1` are rejected with `useNucleus 0`, where the event-averaged Q_s is not computed. Empty, truncated or malformed posterior and nuclear Q_s tables are reported instead of being read silently.
* Fix a malformed entry in the random-seed list (`useSeedList 1`) silently setting that and all later seeds to 0; it is now an error.
* Fix the input parser hanging forever instead of reporting a missing required key when the input file has no `EndOfFile` line.
* Fix error messages that were never printed before the program stopped (e.g. a wrong lattice size or grid length of a read Wilson-line file), or printed as info; errors are now printed as errors, and problems the run continues after as warnings.
* Fix messages of parallel OpenMP threads interleaving in the output.
* Fix `Glauber`'s destructor deleting any `tmp.dat` file in the working directory.
* Fix a memory leak in the matrix logarithm (`Matrix::logmPade`), which leaked a GSL integration table on every call.
* Fix a race on MPI runs with `writeOutputsToHDF5 1` where rank 0 merged and deleted the per-rank HDF5 files while other ranks were still writing their last event.
* Fix the HDF5 collection (`writeOutputsToHDF5 1`): the per-rank files were merged into `RESULTS.h5` with `h5copy` and deleted even when that failed (e.g. without `h5copy` installed), losing every collected event; they are now copied with `h5py`, and a file that cannot be merged is kept. The script also failed with NumPy 2 (`np.string_`), and was only found when the run started in the repository root.

### Removed
* Remove the input parameter `readMultFromFile`: its post-processing mode read `multiplicity<id>.dat` and `NpartdNdy<id>.dat`, which no version of the code writes under these names. `GluonMultiplicity::readNkt()` is kept for such files.
* Remove code that was never called or had no effect: unused functions (e.g. `Evolution::evolveUfast`, `GaugeFix::gaugeTransform`, `Init::findUInForwardLightconeBjoern`, `FFT::fftnMany`), the `Spinor` class, `Util.cpp`/`Util.h` (duplicates of `Setup`'s functions), unused member variables, `Matrix` operators that were never defined or used, and commented-out debug output.
* Remove `src/INTELmakefile`, an unmaintained alternate build file (last touched in 2018) that referenced a `Spinor.cpp`/`Spinor.h` pair no longer in the repo and was not wired into the CMake build.
* Remove the make-based build (`GNUmakefile`, `src/GNUmakefile`); build with CMake as described in the README.

[Link to diff from the last 1.x commit](https://github.com/schenke/ipglasma/compare/b36b2a9...2.0.0)

## 1.x (after 1.0, not released)

The changes on the default branch (`master`, now `main`) after the 1.0 tag, up to b36b2a9. Input parameters are given by their names at that time.

### Input / Output
* Add binary output of the initial Wilson lines (`writeInitialWilsonLines 2`; `1` writes text). They are only written when `writeInitialWilsonLines` is set.
* Add `readInitialWilsonLines` (`1` text, `2` binary) to read the initial Wilson lines `V-<n>` from file instead of computing them; they are shifted by the sampled impact parameter. Writing Wilson lines at b ≠ 0 gives a warning, and a missing file an error.
* Add the output `NgluonEstimators<id>.dat` with estimators of the gluon number, e.g. to trigger on high multiplicities.
* Add the output `eccentricities<id>.dat`, with the eccentricities weighted by ε u^τ and without an energy-density cutoff.
* Add `minimumQs2ST`: events with a smaller Q_s,min² S_T are rejected, to trigger on high-multiplicity events.
* Name the hydro file at the final time `epsilon-u-Hydro-TauHydro-<id>.dat` instead of `epsilon-u-Hydro-t<tau>-<id>.dat`, and no longer write it at the first time step and at 0.1, 0.2, 0.4 and 0.6 fm/c.

### Added
* Add deformed Woods-Saxon nuclei with β2, β3, β4 and the triaxiality γ (`beta2`, `beta3`, `beta4`, `gamma`), and the Woods-Saxon radius and diffuseness (`R_WS`, `a_WS`), read with `setWSDeformParams 1`.
* Add `force_dmin_flag` and `d_min` to impose a minimum distance between the nucleons of deformed nuclei.
* Add a neutron skin: protons and neutrons are sampled from Woods-Saxon distributions that differ by `dR_np` and `da_np`.
* Add a random 3D rotation of each nucleus (all three Euler angles) before the collision, and `rotateReactionPlane` for a random reaction-plane angle.
* Add polarized deuterons (`polariztionProjectile`, `polariztionTarget`, `polarizationProjectileJz`, `polarizationTargetJz`).
* Add nucleon configuration files (`nucleonPositionsFromFile 1`) for He3, He4, C12, O16 (including PGCM and NLEFT configurations), Ne20, Ar40, Au197 and Pb208, in a compact binary format, with a script that converts text tables to it.
* Add fluctuating hot spots: the number of hot spots per nucleon (`NqFluc`, with mean `useConstituentQuarkProton`) and the width of each hot spot (`BGqVar`) fluctuate, with a minimum distance between hot spots (`dqMin`).
* Add posterior parameter sets of the nucleon substructure (`SubNucleonParamType`, `SubNucleonParamSet`), which set the substructure parameters of each event from a Bayesian posterior.
* Add `computeGluonMultiplicity` to switch the gluon multiplicity computation off.
* Make MPI optional at compile time.
* Add a Dockerfile and Singularity files (`docker/`).

### Changed
* Draw the random seed of `useTimeForSeed 1` from `std::random_device` instead of the time.
* Set the time step from the lattice spacing; the input `dtau` is ignored.
* Solve for the gauge links in the forward light cone with a new root finder, which converges robustly also at the lattice edges and in the low-density region, with a higher precision, so there is no noise outside the reaction region.
* No longer use periodic boundary conditions at the lattice edges in the classical evolution (the gauge fixing still does).

### Fixed
* Fix a NaN in the matrix exponential in the very-low-density region.
* Fix the random-number generator being re-initialized every event; it is now initialized once so a fixed seed reproduces a multi-event run.
* Fix a sporadic segfault in the forward-light-cone solver.
* Fix the T^μν output interpolation to use the simulation-lattice spacing and correct cell-boundary handling.
* Fix the hot-spot widths and normalizations of a nucleon keeping the values of the previous nucleon when the number of hot spots fluctuates, and the Frobenius and 1-norms of a matrix being computed from single elements.
* Fix the target's hot-spot normalization using the projectile's number of hot spots in asymmetric collisions.
* Fix the target's nucleon configurations being read with the projectile's mass number, and assign the protons of configuration-file nuclei with the seeded random-number generator.
* Fix the fluctuating-x solve (`useFluctuatingx 1`) using |y| instead of y, and use the floating-point `std::abs` in its convergence test. Events without an overlap region are rejected.
* Fix the `rmax` cutoff when a nucleon's thickness vanishes (the logarithm of 0).
* Replace Wilson lines that are 0 at the lattice edge by the identity.
* Fix the parameters being written before the event ID was set.
* Fix an uninitialized variable when reading Wilson lines in binary format.
* Fix memory leaks from raw-pointer use and from re-allocating the lattice in every loop iteration.

[Link to diff from 1.0](https://github.com/schenke/ipglasma/compare/1.0...b36b2a9)

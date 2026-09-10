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

### Input / Output
* Add input parameters to configure a JIMWLK small-x evolution stage (`useJIMWLK`, `alphas_jimwlk`, `mu0_jimwlk`, `Lambda_QCD_jimwlk`, `m_jimwlk`, `Ds_jimwlk`, `x_projectile_jimwlk`, `x_target_jimwlk`, `simpleLangevin`, `saveSnapshots`).
* Add input parameters for deformed Woods-Saxon nuclei (`beta2`, `beta3`, `beta4`, `gamma`, `setWSDeformParams`, `force_dmin_flag`, `d_min`, `dR_np`, `da_np`).
* Add input parameters for fluctuating sub-nucleon (constituent-quark) structure (`dqMin`, `BGqVar`, `NqFluc`, `useSmoothNucleus`, `shiftConstituentQuarkProtonOrigin`).
* Add input parameters for random reaction-plane/nucleus rotation (`rotateReactionPlane`) and deuteron/light-nucleus polarization (`polariztionProjectile`, `polariztionTarget`, `polarizationProjectileJz`, `polarizationTargetJz`).
* Add `minimumQs2ST` input parameter to trigger only on high-multiplicity events.
* Add binary format support for reading and writing initial Wilson lines (`writeWilsonLines`/`readInitialWilsonLines` = 2) and for JIMWLK snapshots.
* Add `evolvedFields*.ipgf` binary snapshot output and `NgluonEstimators.dat`/eccentricity output files.
* Add nucleon configuration tables for He3, He4, C12, O16, Ne, Ne22 and Ar, including deformed and ab initio (PGCM/NLEFT) variants.
* Replace the default nuclear-`Qs^2` table `qs2Adj_vs_Tp_vs_Y_200.in` with `qs2Adj_vs_Tp_vs_Y_240.in`, which extends the covered `T_p` range 10x to avoid the "T out of range, using maximal T in table" clamping warning at high local thickness; `Init::iTpmax_` is updated accordingly (200 → 240).
* Remove the `Nc` input parameter; the code has always been SU(3)-only and now hardcodes it internally.
* Remove the dead `tDistNu`, `useFatTails` and `writeEvolution` input parameters, which never had any effect on the simulation.

### Added
* Add a JIMWLK small-x evolution stage, run on the projectile and target Wilson lines before the classical Yang-Mills evolution.
* Add a new, more robust forward-lightcone Wilson-line solver (`Init::findUInForwardLightconeChun`).
* Add support for light nuclei with realistic nucleon configurations sampled from file, plus a helper script to convert the configuration tables to a compact binary format.
* Add deformed Woods-Saxon nucleus sampling and full 3D nucleus rotation before collisions.
* Add fluctuating constituent-quark sub-structure: number and width of constituent quarks sampled per nucleon from posterior parameter sets, with a minimum inter-quark distance.
* Add eccentricity/anisotropy computation and output as a function of time.
* Add a flag to enable or disable the gluon-multiplicity calculation.
* Add profiling instrumentation (`Instrumentation.h`/`.cpp`) used throughout initialization and evolution.
* Add validation tools for J/Psi and vector-meson production spectra under `utilities/`.

### Changed
* Rewrite `Matrix`, `Cell` and `Lattice` as a structure-of-arrays layout with a fixed 3x3 `Matrix`, and switch random-number sampling to bulk generation, substantially speeding up both the classical Yang-Mills and JIMWLK evolution.
* Optimize the JIMWLK evolution kernel and noise generation.
* Parallelize and batch `FFT::fftnArray` (used by the per-step JIMWLK noise/kernel transforms) the same way `FFT::fftn` already was: pack all planes into the shared scratch buffer and execute them concurrently via FFTW's thread-safe new-array interface, instead of looping over them serially. Profiling on a 256x256 lattice showed this cut the FFT phase from 50% to 34% of total event time (~24% faster overall).
* Replace `JIMWLK::evolutionStep`'s per-step noise generation, a serial loop of scalar `Random::gauss()` calls, with `Random::gaussBulk` (already used the same way in `Init.cpp`), which reproduces the identical value stream but lets the scatter into `xi2_` run in parallel. Profiling showed the scalar loop was 18% of total event time; this step alone gave a further ~8% speedup on top of the `fftnArray` change above, with bit-identical output.
* Change the nucleus/impact-parameter sampling order and initialization to a two-step field-shifting procedure.
* Update the default `QsmuRatio` to 0.643, following arXiv:2207.03712.
* Rename several input parameters and internal variables for clarity (e.g. add an enum class for `NucleusRole`).
* `setWriteEpsilonUHydro` now uses `iFindOptional`, so older input files without this key still work.
* Make the MPI dependency optional at compile time.
* Reformat the whole codebase to a consistent style.
* Standardize header include guards to a single `SRC_<FILE>_H_` style across all headers.
* Rename `jimwlk.cpp`/`.h` to `JIMWLK.cpp`/`.h` and the `pretty_ostream` class/files to `PrettyOstream`, to match their class names.
* Rename `utilities/read_test.cc`/`save_to_binary.cc` to `.cpp`, matching the rest of the codebase.
* Rename `Setup`/`Util`/`Glauber`'s `PascalCase` functions (`IFind`, `PrintGlauberData`, `ReadInVx`, ...) and `Init`'s `snake_case` nucleus-generation functions to `camelCase`, for a single consistent naming convention.
* Standardize private member-variable names to a trailing-underscore convention (`size_`, `mode_`, ...) across `Parameters`, `Init`, `Random`, `Cell`, `Evolution`, `Lattice`, and `Matrix`.
* Encapsulate `Glauber`: move its data members from `public` to `private` (with trailing underscores) and add a `getGlauberData()` accessor for the one member that was read from outside the class.
* Rename the remaining `PascalCase`/`snake_case` methods to `camelCase`: `Matrix::FrobeniusNorm`/`OneNorm`/`logm_pade`, `GaugeFix::FFTChi`, `Random::Gauss`/`GaussBulk`/`Poisson`, `Evolution::Tmunu`, `Lattice::WriteWilsonLines`/`WriteSU3Matricies` (also fixing a "Matricies" typo to "Matrices"), and `PrettyOstream::get_memory_usage`.
* Rename `Glauber`'s lowercase `tiny`/`limit` macros to `TINY`/`LIMIT`, matching the `ALL_CAPS` convention used by every other macro in the codebase.
* Modernize `Init`'s `Initialization_method` to a scoped `enum class InitializationMethod` with `PascalCase` enumerators, matching `NucleusRole`.
* Simplify `Glauber`'s `Nucleus`/`Data` from C-style `typedef struct` to plain `struct` declarations.
* Standardize function-parameter naming: `Cell`'s setters now use `x` like every other class' setters (instead of `in`), and the "force minimum inter-nucleon distance" flag is now consistently named `forceDminFlag` everywhere it appears (`Glauber`, `Init`), instead of drifting between `force_dmin`, `forceDminFlag`, and `force_dmin_flag` depending on the file.
* Replace `Glauber::sampleTARejection`'s magic-number `int PorT` parameter with the existing `NucleusRole` enum (moved from `JIMWLK.h` to `Glauber.h`, where it more naturally belongs), removing the last raw `1`/`2` projectile/target literals from `Init`.
* Fix several declaration/definition parameter-name mismatches found by cross-checking every header against its `.cpp`: `Setup::listFind` (`fileName`/`paramName` → `file_name`/`st`, matching its own declaration and sibling functions), `Glauber::makeCoeff` (`bb` → `b`, no longer needed now that the colliding public member is `b_`), `Glauber::readInVx`/`readInVy` (name the previously-unnamed `char *` parameter `file_name`), `Init::getNuclearQs2` (declaration said `Qs2atZeroY`, but the parameter is actually a temperature-like table lookup value, named `T` in the implementation and now in the declaration too), and `Matrix::expmCoeff` (`result` → `out`, matching its declaration).
* Rename `Matrix::traceOfProdcutOfMatrix`'s parameters `M1`/`M2` to `a`/`b`, matching every other two-matrix function in `Matrix`/`SU3.h`, and rename `Init::getUfromExponent`'s parameter `in` to `Q`, matching the identically-meaning parameter of `Matrix::expmCoeff` that it forwards to.
* Fix three more private member variables that were missed by the earlier trailing-underscore pass: `Group::t` → `t_`, `FFT::p`/`pback` → `p_`/`pback_`, and `PrettyOstream::message_stream` → `messageStream_` (also fixing its casing). Also rename `Parameters`' backing member for the "force minimum inter-nucleon distance" flag from `force_dmin_flag_` to `forceDminFlag_`, matching the parameter name it was already unified to everywhere else.
* Fix a batch of `PascalCase` quantity names left over from before the member-variable and struct-field naming passes, none of which were physics-notation exceptions: `Parameters`' `Target_`/`Projectile_`/`SigmaNN_`/`Psi_`/`NucleusQsTableFileName_` members, `Glauber`'s `GlauberData_` member, its `Data` struct's `SigmaNN`/`Target`/`Projectile`/`SCutOff`/`InterMax` fields and `Nucleus`'s `AnumFunc`/`AnumFuncIntegrand`/`DensityFunc` fields, and `Glauber::initGlauber`'s `SigmaNN`/`Target`/`Projectile` parameters (and its local `Target_Name`/`Projectile_Name` copies) are now `sigmaNN`/`target`/`projectile`/etc., matching the lowercase-leading convention already used for every other member, field and parameter (e.g. `beta2_`/`gamma_` were already lowercase Greek letters, unlike the capitalized `SigmaNN_`).
* Rename `Phys_consts.h` to `PhysConst.h`, matching its `PhysConst` namespace, and rename its `small_eps` constant to `smallEps` (its sibling constants, `hbarc`/`m_pion`/`m_kaon`/`m_proton`, already followed the codebase's symbol/symbol_subscript physics-notation convention).

### Fixed
* Fix a NaN in the matrix exponential in the very-low-density region.
* Fix the random-number generator being re-initialized every event; it is now initialized once so a fixed seed reproduces a multi-event run.
* Fix a sporadic segfault in the forward-lightcone Wilson-line solver.
* Fix the Tmunu output interpolation to use the simulation-lattice spacing and correct cell-boundary handling.
* Fix several bugs in fluctuating-Nq and fluctuating-Qs sampling.
* Fix the classical Yang-Mills evolution when initial Wilson lines are read back in from disk.
* Fix an uninitialized-variable bug when reading Wilson lines in binary format.
* Fix memory leaks from raw-pointer use and from re-allocating the lattice on every loop iteration.
* Fix a duplicated-getter bug where setting `polariztionProjectile` alone had no effect, because the target polarization flag was checked twice instead of the projectile flag.
* Fix `alphas_jimwlk`: it was parsed as an integer (silently truncating fractional fixed-coupling values to 0) and was ignored by the JIMWLK kernel even when set to a valid value.
* Fix `b`/`phi_RP` being left uninitialized when `useNucleus=0`.
* Fix `Init::readInNucleusConfigs` hanging forever in an infinite loop instead of exiting with an error message when the requested nucleon-configuration file doesn't exist.
* Fix `Glauber::readInVx`/`readInVy` ignoring every `fscanf` return value: they now exit with an error message on a missing file, a missing `"EndOfData"` marker (which previously could loop forever on truncated input, the same bug as `readInNucleusConfigs`), or fewer entries than requested (which previously read silently stale/garbage values past EOF).
* Fix `Matrix::logmPade` leaking its GSL Gauss-Legendre integration table on every call.
* Fix `force_dmin_flag`/`d_min` only being read from the input file when `setWSDeformParams=1`, even though `Glauber::findNucleusData` applies them unconditionally; with `setWSDeformParams=0` (the common case) and a species whose built-in deformation is nonzero (e.g. Au, Cu), this left `Parameters`' `forceDminFlag_`/`d_min_` uninitialized, making nucleon sampling depend on garbage memory and non-reproducible even for a fixed seed.
* Fix `Init::computeCollisionGeometryQuantities` printing `param->getalphas()` to `usedParameters*.dat` before `param->setalphas()` was ever called, so the fixed-coupling `alpha_s` value it recorded was always stale/uninitialized.
* Fix run-to-run non-reproducibility at fixed seed: `FFTW_MEASURE` times candidate FFT algorithms against the current machine state and can pick a different plan (hence different floating-point rounding) between separate invocations of the same binary on the same problem size. Building with the new `-DIPGLASMA_DETERMINISTIC_FFT=ON` CMake option pins the plan via an on-disk FFTW wisdom cache (`FFT.h`), so repeated runs converge onto the same, previously-measured plan instead of re-deciding (off by default, since the first run in a directory pays the same one-time measuring cost `FFTW_MEASURE` always has). Also fixes `-DIPGLASMA_DETERMINISTIC_FFT` never actually reaching the compiler when passed the previously-documented way: every branch in `CMakeLists.txt` unconditionally overwrote `CMAKE_CXX_FLAGS`, silently discarding anything passed via `-DCMAKE_CXX_FLAGS` on the command line.

### Removed
* Remove functions that were declared or defined but never called, including `Evolution::evolveUfast`/`multiplicitynkxky`/`correlations`/`anisotropy`, `GaugeFix::gaugeTransform`, the `Spinor` class and `Matrix::reu`/`reu2`/`imag`, `Init::solveAxbComplex`/`multiplicity`/the 2-argument `rotate_nucleus` overload/`findUInForwardLightconeBjoern`, `MyEigen::test`, `FFT::fftnMany`, `Glauber::FindXorg`/`PAB`/`AnumHulthenInt`, and about a dozen unused `Parameters` getter/setter pairs.
* Remove `Glauber`'s dead `findNucleusData` overload (declared but never defined or called) and rename the surviving `findNucleusData2` to `findNucleusData`.
* Remove the old serial validation script and other now-unused files.
* Remove unused private member variables found by checking every class for members with zero references outside their own declaration: `Parameters::myPI_`/`myhbarc_` (orphaned duplicates of `PhysConst::hbarc` and a pi constant, both superseded), `Glauber::tempFunc_` (and its now-unused `ptr_func` typedef) and `currentTAB_`, and `Init::As_`.

[Link to diff from previous version](https://github.com/schenke/ipglasma/compare/1.0...2.0.0)

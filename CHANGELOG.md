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
Date: 2026-XX-XX

### Input / Output
* Add input parameters to configure a JIMWLK small-x evolution stage (`useJIMWLK`, `alphas_jimwlk`, `mu0_jimwlk`, `Lambda_QCD_jimwlk`, `m_jimwlk`, `Ds_jimwlk`, `x_projectile_jimwlk`, `x_target_jimwlk`, `simpleLangevin`, `saveSnapshots`).
* Add input parameters for deformed Woods-Saxon nuclei (`beta2`, `beta3`, `beta4`, `gamma`, `setWSDeformParams`, `force_dmin_flag`, `d_min`, `dR_np`, `da_np`).
* Add input parameters for fluctuating sub-nucleon (constituent-quark) structure (`dqMin`, `BGqVar`, `NqFluc`, `useSmoothNucleus`, `shiftConstituentQuarkProtonOrigin`).
* Add input parameters for random reaction-plane/nucleus rotation (`rotateReactionPlane`) and deuteron/light-nucleus polarization (`polariztionProjectile`, `polariztionTarget`, `polarizationProjectileJz`, `polarizationTargetJz`).
* Add `minimumQs2ST` input parameter to trigger only on high-multiplicity events.
* Add binary format support for reading and writing initial Wilson lines (`writeWilsonLines`/`readInitialWilsonLines` = 2) and for JIMWLK snapshots.
* Add `evolvedFields*.ipgf` binary snapshot output and `NgluonEstimators.dat`/eccentricity output files.
* Add nucleon configuration tables for He3, He4, C12, O16, Ne, Ne22 and Ar, including deformed and ab initio (PGCM/NLEFT) variants.
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

### Removed
* Remove functions that were declared or defined but never called, including `Evolution::evolveUfast`/`multiplicitynkxky`/`correlations`/`anisotropy`, `GaugeFix::gaugeTransform`, the `Spinor` class and `Matrix::reu`/`reu2`/`imag`, `Init::solveAxbComplex`/`multiplicity`/the 2-argument `rotate_nucleus` overload/`findUInForwardLightconeBjoern`, `MyEigen::test`, `FFT::fftnMany`, `Glauber::FindXorg`/`PAB`/`AnumHulthenInt`, and about a dozen unused `Parameters` getter/setter pairs.
* Remove the old serial validation script and other now-unused files.

[Link to diff from previous version](https://github.com/schenke/ipglasma/compare/1.0...2.0.0)

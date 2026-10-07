# Output files

This page lists every file IP-Glasma writes: when it is written, in which
order, and what it contains. The input parameters that switch the outputs on
are described in the [README](README.md#output); this page describes the
files themselves.

## Conventions

- **Location.** All files are written to the working directory, except the
  Wilson-line and geometry files, which go to `wilsonLinePath`, and the
  profile and fingerprint files, which go to `IPGLASMA_PROFILE_DIR`.
- **Event id.** `<id>` in a file name is the event id
  `rank + iev * nRanks`, where `iev` = 0, 1, … counts the events of an MPI
  rank and `nRanks` is the number of MPI ranks.
- **Times.** `<tau>` in a file name is the proper time in fm/c, written with
  the default stream precision (6 significant digits, e.g. `0.4` or
  `0.0996094`).
- **Text precision.** Text files use 6 significant digits unless stated
  otherwise. Exceptions: the parameter values in `usedParameters` are
  written exactly (shortest round-trip form), the posterior-set values there
  with 9 digits, the Wilson-line text files with 15, the gluon-spectrum JSON
  and the fingerprint file with 17, and the profile file with 12.
- **Lattice coordinates.** The lattice has `size` × `size` sites with spacing
  `a = L/size`. Site `(ix, iy)` sits at `x = -L/2 + a*ix`, `y = -L/2 + a*iy`
  [fm]. In memory and in the Wilson-line files the site index is
  `ix*size + iy` (`ix` outer, `iy` inner).
- **Output grid.** The hydro, Jazma and T^μν files are interpolated (bilinear)
  onto a separate grid:
  - `sizeOutput` × `sizeOutput` points with spacing
    `dx = LOutput/sizeOutput`;
  - point `(ix, iy)` sits at `x = -LOutput/2 + dx*ix` (same for y);
  - `etaSizeOutput` copies in η with spacing `dEtaOutput`, centred on η = 0.
    The fields are boost invariant, so all η slices are identical.
- **Running coupling.** With `runningCoupling 1`, the hydro and T^μν outputs
  are multiplied by the factor g²/(4π α_s) of the gluon spectrum and the
  eccentricities, with α_s evaluated at `runningCouplingQsFactor` × the Q_s
  that `runWithQs` selects: averaged over the overlap region, or with
  `runWithLocalQs 1` the one of each cell, interpolated to the output grid
  like the fields.
- **Byte order.** The binary formats with a JSON header (`.ipgt`, `.ipgw`,
  `.ipgf` and the Wilson-line geometry files) are little-endian by
  definition, and the code refuses to write (or, for the geometry files, to
  read) them on a big-endian host. The binary Wilson-line files are written
  in the host's byte order.
- **Short runs.** With a `maxTime` below one time step (0.1 `a`), there is no
  time step, so none of the files written at the final time (hydro, Jazma,
  T^μν, eccentricities, multiplicity) is written.

## Overview

| File | Section | Written when | Format |
|---|---|---|---|
| `usedParameters<id>.dat` | usedParameters | always | text, valid input file |
| `WilsonLine[_x_<x>]_<n>[.txt]` | Wilson lines | `writeWilsonLines 1` or `2` | text or binary |
| `WilsonLineGeometry_<n>` | Wilson-line geometry | `writeWilsonLines 1` or `2`, `writeWilsonLineGeometry 1` (default), `useNucleus 1` | binary with JSON header |
| `initialWilsonLines<id>.ipgw` | Initial Wilson lines snapshot | `writeWilsonLineSnapshot 1`, color charges sampled | binary with JSON header |
| `NpartList<id>.dat`, `NcollList<id>.dat` | Participants and binary collisions | `mode 1`, `useNucleus 1`, `useSmoothNucleus 0`, `writeNpartList 1`, `writeNcollList 1` | text |
| `NgluonEstimators<id>.dat` | Gluon number estimators | `mode 1`, `useNucleus 1`, `writeNgluonEstimators 1` | text |
| `Tmunu-t<tau>-<id>.ipgt` / `.dat` | Energy-momentum tensor | `mode 1`, `writeTmunu 1` | binary or text |
| `epsilon-u-Hydro-t<tau>-<id>.dat`, `epsilon-u-Hydro-TauHydro-<id>.dat` | Hydro initial conditions | `mode 1`, `writeHydro 1` | text |
| `Jazma-Hydro-t<tau>-<id>.dat` | Jazma energy density | `mode 1`, `writeJazma 1` | text |
| `eccentricities<id>.dat` | Eccentricities | `mode 1`, `computeEccentricities 1` | text, appended |
| `NpartdNdy-t<tau>-<id>.dat` | Multiplicity summary | `mode 1`, `computeGluonMultiplicity 1` | text |
| `gluonMultiplicity<id>.json` | Gluon spectrum | `mode 1`, `computeGluonMultiplicity 1` | JSON |
| `multiplicityHadrons<id>.dat` | Hadron spectrum | as above and `writeHadronSpectrum 1` | text |
| `RESULTS_rank<rank>.h5`, `RESULTS.h5` | HDF5 collection | `writeOutputsToHDF5 1` | HDF5 |
| `ipglasma_fftw_wisdom.dat` | Diagnostic files | built with `-DIPGLASMA_DETERMINISTIC_FFT=ON` | FFTW wisdom |
| `ipglasma_profile_rank<rank>.tsv`, `ipglasma_fingerprint_rank<rank>.tsv` | Diagnostic files | environment variables | tab-separated text |

## Order within an event

1. `usedParameters<id>.dat`: the input parameters and the random seed.
2. With sampled color charges, after the Wilson lines are built:
   1. `initialWilsonLines<id>.ipgw` (`writeWilsonLineSnapshot 1`);
   2. the initial Wilson lines (`writeWilsonLines` > 0, and without JIMWLK or
      with `jimwlkSaveSnapshots 1`);
   3. the geometry files of both nuclei (`writeWilsonLines` > 0,
      `writeWilsonLineGeometry 1`, `useNucleus 1`). A run that reads Wilson
      lines does not write them again: they have the same names as the files
      it read.
3. With JIMWLK:
   1. the Wilson-line snapshots at the x values of `jimwlkXSnapshotList`
      (`jimwlkSaveSnapshots 1`), each at the evolution step closest to it;
   2. the final Wilson lines at `projectileX`/`targetX`
      (`writeWilsonLines` > 0).
4. `mode 1` with `useNucleus 1` (sampled nuclei, or nuclei read with their
   Wilson lines), for each impact parameter tried, as switched on:
   `NcollList<id>.dat`, then `NpartList<id>.dat`, then, once an impact
   parameter is accepted, its collision geometry is appended to
   `usedParameters<id>.dat`, then `NgluonEstimators<id>.dat` (not for a try
   rejected by `useFixedNpart` or for lack of overlap). Each try overwrites
   the files of the previous one.
5. `mode 1`, during the evolution: the hydro, T^μν and Jazma files that are
   switched on, in this order, at each of the `outputTimes` and at the final
   time. Each output time is rounded down to a time step; times below one
   time step or not before the final step are skipped, and times that round
   to the same step are written once.
6. `mode 1`, after the final time: `eccentricities<id>.dat`
   (`computeEccentricities 1`), then, with `computeGluonMultiplicity 1`,
   `multiplicityHadrons<id>.dat` (`writeHadronSpectrum 1`),
   `NpartdNdy-t<tau>-<id>.dat` and `gluonMultiplicity<id>.json`.
7. With `writeOutputsToHDF5 1`, some of the event's text files are moved into
   `RESULTS_rank<rank>.h5` (see "HDF5 collection"); at the end of the run
   these are merged into `RESULTS.h5`.
8. The fingerprint row (`IPGLASMA_FINGERPRINT`) and, at the end of the
   event, the profile rows (`IPGLASMA_PROFILE`).

## usedParameters

`usedParameters<id>.dat` (`main.cpp`, `CollisionGeometry::writeUsedParametersFile()`).

A valid input file with the values of every input parameter the event used,
in input-file syntax (`key value`, one per line), in the order of the
parameter table (`src/ParameterTable.cpp`). Parameters that are not read
under the chosen settings are left out. Comment lines (`#`) carry the event
information:

- a header with the event id and the date;
- with a posterior parameter set (`subNucleonParamType` > 0), the set that
  was used and its values of `m`, `BG`, `BGq`, `smearingWidth`, `NqBase`,
  `QsMuRatio` and `dqMin`;
- `# Random seed used on rank <rank>: <seed>`;
- with `mode 1` and `useNucleus 1`, once an impact parameter is accepted, a
  block starting with `# Collision geometry of this event:`:
  - `b` [fm], `phiRP`, `Npart` and `Ncoll`;
  - with running coupling, the event-averaged Q_s [GeV] selected by
    `runWithQs` and α_s (0 with `runWithKt 1`); otherwise the fixed α_s;
  - with `readInitialWilsonLines` 1 or 2, `# QsMuRatio = <value> (from the
    Wilson-line geometry files)`, the value the event used instead of the
    input one.

Running the file again gives the same parameters, but not the same event:
the random numbers also depend on the MPI rank and on the event's position
in the run, `useRandomSeed 1` draws a new seed, and `subNucleonParamSet -1`
draws a new posterior set.

## Wilson lines

`<wilsonLinePath>/WilsonLine[_x_<x>]_<n>[.txt]` (`WilsonLineIO::write()`).

- **Name.**
  - `_x_<x>` gives Bjorken x in scientific notation with 5 decimals. Only
    the initial Wilson lines with `useFluctuatingX 1` (no JIMWLK) have
    no fixed x and leave it out; JIMWLK snapshots and final lines always
    have it.
  - `<n> = 2 (seed · N + <id>) + iA`, where seed is the input parameter
    `seed` (also with `useRandomSeed 1` or `useSeedList 1`), N is the number
    of events of the run (events per rank times MPI ranks) and `iA` is 1 for
    the projectile and 2 for the target. So no two files of a run, nor of runs with different
    seeds and the same N, have the same number; one event on one rank gives
    2·seed + 1 and 2·seed + 2.
  - Text files end in `.txt`; binary files have no extension.
- **When.**
  - The initial Wilson lines: without JIMWLK at x = `projectileX`/`targetX`
    (no x with `useFluctuatingX 1`), with JIMWLK at `jimwlkInitialX` (and
    only with `jimwlkSaveSnapshots 1`).
  - JIMWLK snapshots, at the step closest to each value of
    `jimwlkXSnapshotList` and named with that value, and the final Wilson
    lines, see "Order within an event" above.
  - `readInitialWilsonLines 1`/`2` reads the files back, under the names a
    run with the same parameters writes, at the x of `readWilsonLinesX`
    (default: the initial Wilson lines), with `useNucleus 1` together with
    the geometry files.
    The fields are stored as built, centered at the origin; the impact
    parameter is applied after reading, as for sampled nuclei.

**Text format** (`writeWilsonLines 1`): one line per site, `ix` outer and
`iy` inner, with a blank line after each `ix`. Each line holds 20 columns:

```
ix iy Re(V11) Im(V11) Re(V12) Im(V12) Re(V13) Im(V13) Re(V21) ... Re(V33) Im(V33)
```

The matrix elements are row-major and written with 15 significant digits.

**Binary format** (`writeWilsonLines 2`), in the host's byte order:

| Bytes | Type | Content |
|---|---|---|
| 4 | `int32` | `size` (N) |
| 4 | `int32` | N_c = 3 |
| 8 | `double` | `L` [fm] |
| 8 | `double` | `a` [fm] |
| 8 | `double` | Bjorken x of the Wilson lines, the x in the file name; −1 with `useFluctuatingX 1` (no fixed x) |
| N² × 9 × 16 | `double` pairs | (Re, Im) of each matrix element, sites `ix` outer and `iy` inner, elements row-major |

## Wilson-line geometry

`<wilsonLinePath>/WilsonLineGeometry_<n>` (`WilsonLineIO::writeGeometry()`),
one per nucleus, with the same number `<n>` as its Wilson-line files and no
x, since one geometry serves the Wilson lines at every x. Written right
after the Wilson lines are built, with `writeWilsonLines` 1 or 2 and
`useNucleus 1`, unless `writeWilsonLineGeometry 0`, also with JIMWLK when
only the final Wilson lines are written. A run that reads Wilson lines does
not write it again. `readInitialWilsonLines` reads it to sample the
collision geometry of read Wilson lines; it must match the run's `size`,
`L`, `g`, `projectile` and `target`.

Little-endian by definition:

| Bytes | Content |
|---|---|
| 8 | magic `IPGGEO1\0` |
| 8 | length M of the JSON metadata, `uint64` |
| M | JSON metadata: `format`, `version`, `dtype` (`<f8`), `nucleus` (`projectile` or `target`), `species`, `nucleons`, `N`, `L_fm`, `g`, `QsMuRatio`, `blocks`, `native_site_index` |
| nucleons × 4 × 8 | `float64` x, y, z [fm] and proton (1 or 0) of each nucleon, centered at the origin |
| N² × 8 | `float64` g²μ² of each site, as stored on the lattice, sites `ix` outer and `iy` inner |
| N² × 8 | `float64` T_p [GeV²] (the summed nucleon thickness, as stored on the lattice) of each site, in the same order |

## Initial Wilson lines snapshot

`initialWilsonLines<id>.ipgw` (`WilsonLineIO::writeTrainingData()`), with
`writeWilsonLineSnapshot 1`, in any `mode`, but not when the Wilson lines are
read from file (`readInitialWilsonLines` 1 or 2). It holds both nuclei's
Wilson lines right after they are built, before JIMWLK and the
impact-parameter shift.

| Bytes | Content |
|---|---|
| 8 | magic `IPGWIL1\0` |
| 8 | metadata length `n` as little-endian `uint64` |
| `n` | JSON metadata |
| rest | `float32` (little-endian) data of shape `[2, 2, N, N, 3, 3]` |

- **Data axes:** `[beam, complex_part, x, y, row, col]`, with beam 0 = V_A
  (projectile), beam 1 = V_B (target), and complex part 0 = real,
  1 = imaginary.
- **Metadata:** repeats the format name (`ipglasma-initial-wilson-lines`),
  version, dtype, shape and axis order, and adds `fields` (`VA`, `VB`),
  `complex_part`, `native_site_index`, `event_id`, `N`, `Nc`, `L_fm`, `a_fm`
  and `x_A`/`x_B`, the Bjorken x of the projectile's and the target's
  Wilson lines (`jimwlkInitialX` with JIMWLK, otherwise `projectileX`/
  `targetX`; −1 with `useFluctuatingX 1`, where x is not fixed).

## Participants and binary collisions

`NpartList<id>.dat` and `NcollList<id>.dat`
(`CollisionGeometry::determineNpartAndNcoll()`, `computeNcollList()`), with
`writeNpartList 1` and `writeNcollList 1` (the default), in `mode 1` with
`useNucleus 1` (also for nuclei read with their Wilson lines). They are not
written with `useSmoothNucleus 1`. Positions are in fm, in the frame
of the collision: the projectile is centred at +b/2 and the target at −b/2
along the reaction plane. Two nucleons collide if their transverse distance
is below √(σ_NN/π) (`gaussianWounding 0`) or with the Gaussian probability
of GLISSANDO (`gaussianWounding 1`, see the README). With `nucleiToAverage`
n > 1 the lists hold the nucleons of all n nuclei of each kind, first those
of the first nucleus.

**`NpartList<id>.dat`:** one line per nucleon, first all projectile nucleons,
then a blank line, then all target nucleons:

| Column | Content |
|---|---|
| 1, 2 | x, y [fm] |
| 3 | 1 for a proton, 0 for a neutron |
| 4 | 1 if the nucleon collided, 0 otherwise |

In p+p both protons count as participants even if column 4 is 0.

**`NcollList<id>.dat`:** one line per binary collision, with the midpoint
x, y [fm] of the two colliding nucleons.

## Gluon number estimators

`NgluonEstimators<id>.dat` (`CollisionGeometry::writeNgluonEstimatorsFile()`),
with `writeNgluonEstimators 1` (the default), in `mode 1` with
`useNucleus 1`, also with `useSmoothNucleus 1`. One header line (`#`) and one
line of four numbers:

| Column | Content |
|---|---|
| 1 | Q_s²(min) S_T, summed over the whole lattice |
| 2 | Q_s²(avg) S_T over the overlap region |
| 3 | Q_s²(max) S_T over the overlap region |
| 4 | column 1 × ln²(column 3 / column 1) |

Q_s(min/avg/max) is the smaller, average or larger of the two nuclei's Q_s
in a cell, and S_T is the area. All four values are dimensionless.

## Energy-momentum tensor

`Tmunu-t<tau>-<id>.ipgt` (binary, `writeTmunuBinary 1`, the default) or
`Tmunu-t<tau>-<id>.dat` (text) (`MyEigen::writeRawTmunu()`), with
`writeTmunu 1`, at the final time and at `outputTimes`.

Each output-grid point holds ten components, in GeV/fm³:

| # | Name in the metadata | Component |
|---|---|---|
| 0 | `T00` | T^ττ |
| 1 | `Txx` | T^xx |
| 2 | `Tyy` | T^yy |
| 3 | `tau2_Tetaeta` | τ² T^ηη |
| 4 | `neg_T0x` | −T^τx |
| 5 | `neg_T0y` | −T^τy |
| 6 | `neg_tau_T0eta` | −τ T^τη |
| 7 | `neg_Txy` | −T^xy |
| 8 | `neg_tau_Tyeta` | −τ T^yη |
| 9 | `neg_tau_Txeta` | −τ T^xη |

Points outside the lattice, or with T^ττ below 10⁻¹⁶ GeV/fm³, get
(10⁻¹⁶, 0.5·10⁻¹⁶, 0.5·10⁻¹⁶, 0, …, 0).

**Binary format** (`.ipgt`):

| Bytes | Content |
|---|---|
| 8 | magic `IPGTMU01` |
| 4 | metadata length `n` as little-endian `uint32` |
| `n` | JSON metadata |
| rest | `float32` (little-endian) data of shape `[sizeOutput (y), sizeOutput (x), 10]` |

- **Data order:** y is the outer axis, then x, then the component.
- **Metadata:** `format` (`ipglasma-tmunu`), `version`, `dtype`, `shape`,
  `axis_order`, `components`, `tau_fm`, `eta_points`, `deta`, `dx_fm`,
  `dy_fm` and `event_id`.
- The file holds one η slice; `eta_points` and `deta` describe the intended
  η grid.
- `utilities/read_tmunu.py` reads both formats.

**Text format** (`.dat`):
- One header line: `# dummy 1 etamax= <etaSizeOutput> xmax= <sizeOutput> ymax= <sizeOutput> deta= <dEtaOutput> dx= <dx> dy= <dx>`.
- Then one line per grid point, `iy` outer and `ix` inner, with a blank line
  after each `iy`. Each line has 12 columns: `ix iy` and the ten components.

## Hydro initial conditions

`epsilon-u-Hydro-t<tau>-<id>.dat` (at `outputTimes`) or
`epsilon-u-Hydro-TauHydro-<id>.dat` (final time) (`MyEigen::writeHydroText()`),
with `writeHydro 1`. Energy density, flow
velocity and shear-stress tensor after Landau matching, in the format read by
MUSIC.

- **Header line:** `# dummy 1 etamax= <etaSizeOutput> xmax= <sizeOutput> ymax= <sizeOutput> deta= <dEtaOutput> dx= <dx> dy= <dx> tau= <tau>`.
- **Order:** one line per grid point, η outer, then x, then y, with a blank
  line after each η slice.
- **Columns (18):**

| Column | Content | Unit |
|---|---|---|
| 1 | η | |
| 2, 3 | x, y | fm |
| 4 | ε | GeV/fm³ |
| 5–8 | u^τ, u^x, u^y, u^η | u^η in fm⁻¹ |
| 9–18 | π^ττ, π^τx, π^τy, π^τη, π^xx, π^xy, π^xη, π^yy, π^yη, π^ηη | fm⁻⁴; components with an η index carry an extra fm⁻¹ each |

Note that ε is in GeV/fm³, but π^μν is in fm⁻⁴ (multiply by ħc for GeV/fm³).

Points with ε ≤ 10⁻¹⁰ GeV/fm³, and points within 0.5 fm of the lattice edge,
are written as vacuum: ε = 0, u = (1, 0, 0, 0), π = 0.

## Jazma energy density

`Jazma-Hydro-t<tau>-<id>.dat` (`MyEigen::writeJazma()`), with
`writeJazma 1`, at the final time and at `outputTimes`; `<tau>` is the time
also for the final one. The same 18 columns and point order (η outer, then
x, then y) as the hydro file, but without the `tau=` entry in the header and
without a blank line between η slices. Here ε ∝ g²μ²_A · g²μ²_B, normalized
to the same total energy as the hydro output, with u = (1, 0, 0, 0) and
π = 0. Only points outside the lattice get ε = 0; the hydro file's 0.5 fm
edge margin and ε threshold do not apply.

## Eccentricities

`eccentricities<id>.dat` (`Eccentricity::compute()`), at the final time with
`computeEccentricities 1`. The line is *appended*, so rerunning an event
id in the same directory adds lines. The weights are ε u^τ (times the
running-coupling factor), over cells with ε at least `eccentricityCutoff`
(default 0), about the energy-weighted centroid.

| Column | Content |
|---|---|
| 1 | τ [fm/c] |
| 2–13 | ε_n, Ψ_n for n = 1…6 (ε_1, Ψ_1, ε_2, Ψ_2, …); ε_n = \|⟨r^m e^{inφ}⟩\| / ⟨r^m⟩ with m = 3 for n = 1 and m = n otherwise |
| 14 | energy-density cutoff `eccentricityCutoff` [GeV/fm³] |
| 15 | √⟨r²⟩ [fm] |
| 16, 17 | largest x and y offset from the centroid, along the row and column through it, where the weight exceeds the cutoff [fm] |
| 18 | impact parameter b [fm] |
| 19 | T_pp [fm⁻²] |
| 20 | area of the cells above the cutoff [fm²] |
| 21 | R̄ = 1/√(1/⟨x²⟩ + 1/⟨y²⟩) [fm] |
| 22 | average energy density of the cells above the cutoff [GeV/fm³] |
| 23 | average of Q_{s,A}² Q_{s,B}² / (g⁴ (ħc)³) over the same cells |

## Multiplicity summary

`NpartdNdy-t<tau>-<id>.dat` (`GluonMultiplicity::compute()`), at the final
time with `computeGluonMultiplicity 1`. It is not written if the event has
no gluons (dN/dy = 0), which ends the event. One line:

| Column | Content |
|---|---|
| 1 | N_part |
| 2 | dN/dy (dN/dη with `usePseudoRapidity 1`) of gluons |
| 3 | T_pp [fm⁻²] |
| 4 | b [fm] |
| 5 | dE/dy [GeV] |
| 6 | random seed |
| 7–9 | `N/A` (placeholders) |
| 10, 11 | dN/dy, dE/dy [GeV] for k_T > 3 GeV |
| 12, 13 | dN/dy, dE/dy [GeV] for k_T > 6 GeV |
| 14 | g²/(4π α_s) at `runningCouplingQsFactor` × the ⟨Q_s⟩ that `runWithQs` selects; 1 with `runningCoupling 0`. With `runWithLocalQs 1` or `runWithKt 1`, the factor applied to the spectrum varies per cell or k_T bin, and this is the one at the event-averaged ⟨Q_s⟩ |

## Gluon spectrum

`gluonMultiplicity<id>.json` (`GluonMultiplicity::writeTarget()`), written
together with the multiplicity summary. A JSON object with:

- `format` (`ipglasma-gluon-target`) and `version`;
- `event_id`, `step`, `tau_fm`;
- `rapidity_variable` (`y`, or `eta` with `usePseudoRapidity 1`);
- `dN`, `dE_GeV` and their binned cross-checks `dN_binned_check`,
  `dE_binned_check_GeV`;
- `mean_kT_GeV`;
- `dN_kT_gt_3_GeV`, `dE_kT_gt_3_GeV`, `dN_kT_gt_6_GeV`, `dE_kT_gt_6_GeV`;
- `Npart`, `Tpp`, `impact_parameter_fm`, `random_seed`;
- `spectrum_definition` (a description);
- four arrays of 100 k_T bins:
  - `kt_GeV`: the bin centres;
  - `dN_d2k_GeV_minus2`: dN/(dy d²k_T) [GeV⁻²];
  - `dE_d2k_GeV_minus1`: dE/(dy d²k_T) [GeV⁻¹];
  - `lattice_bin_counts`: the number of lattice momentum modes in each bin.

The gluon spectrum is measured in transverse Coulomb gauge.

## Hadron spectrum

`multiplicityHadrons<id>.dat` (`GluonMultiplicity::hadronizeAndWrite()`),
with `computeGluonMultiplicity 1` and `writeHadronSpectrum 1`. It is written
before the multiplicity summary, also for an event without gluons. The gluon
spectrum
convolved with KKP fragmentation functions. One line per p_T from 0 to 20 GeV
in steps of 0.1 GeV:

| Column | Content |
|---|---|
| 1 | p_T [GeV] |
| 2 | charged-hadron spectrum |
| 3, 4 | 0 (placeholders) |
| 5 | T_pp [fm⁻²] |
| 6 | b [fm] |

## HDF5 collection

`RESULTS_rank<rank>.h5` and `RESULTS.h5`, with `writeOutputsToHDF5 1`.
After each event, `utilities/combine_events_into_hdf5.py` collects some of the
event's text files into the group `event-<id>` of `RESULTS_rank<rank>.h5` and
then **deletes them**:

- `usedParameters<id>.dat`, each line as a group attribute named by its
  line number;
- as datasets, named like the files:
  - `NcollList<id>.dat` (N_coll × 2) and `NpartList<id>.dat` (one row of 4
    per nucleon, without the blank line between the nuclei);
  - `NpartdNdy-t*-<id>.dat`, with the `N/A` placeholders stored as 0;
  - `epsilon-u-Hydro-t*-<id>.dat` (the hydro files at `outputTimes`, not the
    final `epsilon-u-Hydro-TauHydro-<id>.dat`), without the η, x, y columns
    (15 columns), with the attributes `header`, `x_size`, `y_size`, `dx`,
    `dy`, `nx` and `ny`;
  - the text `Tmunu-t*-<id>.dat`, without the `ix`, `iy` columns (10
    columns), with the same attributes.

The other files, among them the binary `.ipgt` files, stay on disk. At the
end of the run, rank 0 copies the groups of every `RESULTS_rank<rank>.h5`
into `RESULTS.h5`, which is opened for appending, and deletes the per-rank
file; a file that cannot be read
or copied, or one with a group that `RESULTS.h5` already holds, is kept (and
the program warns), and the groups already copied from it are removed again,
so a later run can merge it. The script needs `python3` with `h5py` and `numpy`; the program calls
it from the source tree (`utilities/combine_events_into_hdf5.py`), so the run
can start in any directory.

## Diagnostic files

- **`ipglasma_fftw_wisdom.dat`:** only when built with
  `-DIPGLASMA_DETERMINISTIC_FFT=ON`. FFTW's plan cache, read and written
  again whenever an FFT is set up (several times per event), via a temporary
  `ipglasma_fftw_wisdom.dat.tmp.<pid>`, so later runs choose the same FFT
  algorithms (see "Reproducible runs" in the README).
- **`ipglasma_profile_rank<rank>.tsv`:** with the environment variable
  `IPGLASMA_PROFILE=1`, in `IPGLASMA_PROFILE_DIR` (default `.`).
  - Appended per event, tab-separated, with the columns `rank`, `event`,
    `phase`, `seconds`, `calls`, `percent_event`, including an `event.total`
    row.
- **`ipglasma_fingerprint_rank<rank>.tsv`:** with
  `IPGLASMA_FINGERPRINT=1`, also in `IPGLASMA_PROFILE_DIR`.
  - Per event and field (ε, the T^μν components, …), the columns `rank`,
    `event`, `field`, `count`, `nonfinite`, `hash_fnv1a64`, `mean`, `rms`,
    `min`, `max`, and an `ALL_FIELDS` row with the hash of all fields.
  - Used to compare runs for bit-identical results.

## Writers not used in a normal run

These writers exist for debugging and tests. No input parameter switches them
on.

- **`anisotropy<id>.dat`** (`Eccentricity::compute()` with `doAniso = 1`):
  - two header lines with Ψ_2 and the flow-velocity angle Ψ_U;
  - then ten lines `tau ratio ratio2 angle=<Ψ>` for Ψ = Ψ_U + kπ/8. `ratio`
    is ⟨T^xx − T^yy⟩/⟨T^xx + T^yy⟩ in the frame rotated by Ψ; `ratio2` is the
    same unrotated.
- **`evolvedFields<id>_it<step>.ipgf`** (`Evolution::writeEvolvedFields()`),
  with the step zero-padded to 8 digits (e.g. `evolvedFields0_it00000123.ipgf`):
  - magic `IPGFLD1\0`, then the metadata length as a little-endian `uint64`,
    then JSON metadata;
  - then `float32` data of shape `[6, 2, N, N, 3, 3]`, with axes `[field,
    complex_part, x, y, row, col]` and fields φ, π, E1, E2, U_x, U_y;
  - U_x, U_y and φ are at τ, E1, E2 and π at τ − dτ/2.
- **`<prefix>Phi-<n>.txt`, `<prefix>Pi-<n>.txt`**
  (`Lattice::writeSU3Matrices()`): φ and π in the Wilson-line text format,
  with `<n>` = `<id>` + 2·seed·nRanks for φ and `<id>` + (2·seed + 1)·nRanks
  for π (not the Wilson-line numbering).
- **`NpartdNdy-mod.dat`** (`GluonMultiplicity::readNkt()`): reads
  `multiplicity<id>.dat` and `NpartdNdy<id>.dat` (names no version of the
  program writes), converts dN/dy to dN/dη with the `jacobianMass`/`sqrtS`
  Jacobian, writes one line (N_part, the two dN/dη, T_pp, b) and stops the
  program with exit status 1.

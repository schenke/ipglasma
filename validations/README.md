# Running the validation scripts

This folder contains, each in its own subfolder, two regression tests for
IP-Glasma and the small validation setups of the code paper, which test
parts of the code against exact or published results:

- **`vector_meson_production/`** (`parallel_test_vector_meson_production.py`)
  — runs many independent
  IP-Glasma + JIMWLK events, feeds the resulting Wilson lines to the external
  [`subnucleondiffraction`](https://github.com/hejajama/subnucleondiffraction)
  code, and compares the resulting coherent/incoherent J/ψ cross section to a
  stored reference (`validation_spectra`) and to Fig. 1 ("CGC+shape fluct") of
  [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712).
- **`eccentricity_regression/`** (`parallel_test_eccentricity.py`) — runs
  many independent IP-Glasma
  events (no external dependency), computes the energy-density-weighted
  eccentricity ε₂ of the final-state `Tmunu` snapshot for each event, and
  checks the resulting ε₂ distribution against a saved baseline. This is the
  lighter-weight "did the physics silently change" check and is a good first
  thing to run after touching the code.
- **`dipole_correlator/`** — the dipole correlator of the MV model (constant
  color-charge density) against its exact lattice prediction, for a scan of
  `g2muGeV`, `m` and `Ny`, optionally also of `UVDamp` (see below).
- **`lattice_symmetries/`** — two exact symmetries of the lattice
  formulation, the transposition x ↔ y and the exchange of projectile and
  target, applied to the Wilson lines of a Pb+Pb event (see below).
- **`dilute_limit/`** — the dilute limit of the glasma: the (g²μ)⁴/g²
  scaling of the energy density and the gluon number (see below).
- **`numerical_parameters/`** — how the observables of a Pb+Pb event depend
  on the lattice spacing, `Ny` and `m` at a fixed geometry (see below).
- **`constant_glasma/`** — the glasma of a constant color-charge density
  against the published classical Yang-Mills results of Lappi (see below).
- **`performance/`** — the time per stage of an event and the peak memory
  against the lattice size and the number of threads (see below).
- **`jimwlk_step_size/`** — the convergence of the JIMWLK evolution in its
  step size `jimwlkDs`, and its error at the recommended step sizes, in the
  setup of arXiv:1212.4825 (with BK and the curves of that paper) and in the
  example input (see below).

The two regression tests parallelize over random seeds using
`concurrent.futures.ProcessPoolExecutor`, so each event runs in its own
working directory (`<datadir>/seed_<n>/`) and multiple events run
simultaneously.

## Prerequisites

1. Build the `ipglasma` binary from the repository root:
   ```bash
   ./compile_IPGlasma.sh
   ```
   (dependencies: CMake, FFTW and GSL, see the main README). The scripts
   locate the binary via `--ipglasma-path` (default `../..`, the repository
   root seen from the folder of the script), expecting it at
   `<ipglasma-path>/ipglasma`.
2. For `parallel_test_vector_meson_production.py` only: build
   [`subnucleondiffraction`](https://github.com/hejajama/subnucleondiffraction)
   separately and note its path (needed for `--subnucleondiffraction-path`).
3. Python 3 with `numpy` and `matplotlib` (and optionally `scipy`, used for
   the two-sample KS test in the eccentricity script if available; the
   JIMWLK step-size script needs it).

## Running interactively

From the folder of each script (paths below assume the scripts are run
from there; both scripts also accept absolute paths for `--ipglasma-path`
etc. if you invoke them from elsewhere; their input templates and
references are found next to the script):

### Eccentricity regression test (quick sanity check)

```bash
cd eccentricity_regression
# One-time: record a baseline from a known-good build
python3 parallel_test_eccentricity.py --input-template input_eccentricity \
    --ipglasma-path ../.. --maxevents 40 --max-workers 8 --save-reference

# Later, e.g. after changing the code: check nothing drifted
python3 parallel_test_eccentricity.py --input-template input_eccentricity \
    --ipglasma-path ../.. --maxevents 40 --max-workers 8
```

Useful flags:
- `--datadir DIR` — where per-event working directories and plots go
  (default `./eccentricity_test_output/`).
- `--size`, `--L`, `--maxtime` — override the lattice size / box size / max
  evolution time from the input template (defaults to a small, fast 128²
  lattice).
- `--jimwlk` — keep JIMWLK evolution enabled (slower, more complete physics
  test); disabled by default for speed.
- `--max-workers N` — number of events to run in parallel (default 1; set
  this to the number of cores you have available).
- `--plot-only` — skip launching new runs and just rebuild the plot/summary
  from `seed_*/` directories already present in `--datadir` (useful after an
  interrupted run, if `--keep-logs` was used).
- `--clean` — remove stray temporary files left behind by an interrupted run
  from every `seed_*/` directory under `--datadir`, without touching the
  saved `Tmunu` snapshots, then exit.
- `--keep-logs` — keep the per-event working directories after the run
  instead of deleting them (failed events are always kept).

Output: `epsilon2_distribution.pdf` in `--datadir`, plus
`eccentricity_results.json` (or, with `--save-reference`, an updated
`eccentricity_reference.json` next to the script) recording the ε₂ values and
the exact command used to reproduce the run. The script exits non-zero if the
new distribution's mean/std or KS test drifts beyond tolerance
(`--mean-tol`, `--std-tol`) from the reference.

### Dipole correlator of the MV model

`dipole_correlator/dipole_correlator.py` (issue #123) tests the Wilson lines:
the normalization of the color charges (and that `g` drops out), the
momentum-space kernel with lattice momenta, IR regulator and UV damping,
the FFT conventions and the product over the `Ny` layers.

With a constant color-charge density (`useNucleus 0`), the layers are
Gaussian with the two-point function, summed over the layers,

G(r) = (g²μ a)²/N² Σ_k K(k)² cos(k·r),
K(k) = exp(−|k̃| u)/[(k̃ a)² + (m a)²] (K(0) = 0 for m = 0),

exactly on the periodic lattice, with the lattice momentum
k̃² a² = 4[sin²(k_x a/2) + sin²(k_y a/2)] and the UV damping u. The dipole
correlator D(r) = Re⟨tr V†(x) V(x+r)⟩/N_c therefore tends to
D_G(r) = exp{−C_F [G(0) − G(r)]} for `Ny` → ∞, with corrections of order
1/`Ny`; no continuum or infinite-volume approximation enters. The script
runs IP-Glasma (`runEvolution 0`, `writeWilsonLines 2`, input template
`dipole_correlator/input_dipole`) for each setting and seed, computes D(r)
from the two nuclei of each run (two independent configurations) along the
axes and the diagonal, averaged over all sites, and checks that

- D approaches D_G for large `Ny`,
- `Ny` (D − D_G) at the r where D_G = e^{−1/2} tends to a constant,
- the extrapolation D − D_G = d_∞ + c/`Ny` (for `Ny` ≥ 5) gives a d_∞
  consistent with 0,
- a run with `g 2` gives the same D as with `g 1` (to rounding).

It also reports Q_s/g²μ, with D(r_s) = e^{−1/2} and Q_s² = 2/r_s², to be
compared with T. Lappi, Eur. Phys. J. C 55 (2008) 285
([arXiv:0711.3039](https://arxiv.org/abs/0711.3039)).

```bash
cd dipole_correlator
# quick test: 64², g²μ = 1 GeV, m = 0.2 GeV, Ny = 1...50, 2 seeds (seconds)
python3 dipole_correlator.py --ipglasma-path ../.. --quick --max-workers 4
# the scan for the paper: 512², L = 25.6 fm, g²μ = 0.5, 1, 2 GeV,
# m = 0, 0.2, 0.4 GeV, Ny = 1...100, 4 seeds
python3 dipole_correlator.py --ipglasma-path ../.. --max-workers 8
# optionally also the UV-damped kernel, at g²μ = 1 GeV and m = 0.2 GeV
python3 dipole_correlator.py --ipglasma-path ../.. --max-workers 8 \
    --uvdamp 0.1
```

The UV damping is off by default (`UVDamp 0`, the MV model used in
practice); `--uvdamp` adds a `Ny` scan for each given value.

`--ipglasma-path` is a directory with the `ipglasma` executable and
`qs2Adj_vs_Tp_vs_Y_240.in`. The results (D(r) and D_G(r) for every setting,
the fit, the checks, Q_s/g²μ) go to `<datadir>/dipole_summary.json`
(default `./dipole_output`), with one figure per setting
(`dipole_<setting>.pdf`: D(r) against D_G(r) for several `Ny`, and
`Ny` (D − D_G) with its statistical error). Each run caches its D(r) in its
working directory, so `--analyze-only` redoes the analysis, and a rerun only
runs what is missing; `--keep-files` keeps the Wilson lines. The script exits
non-zero if a check fails. At 256² with 8 configurations per `Ny`
(`--size 256 --L 25.6 --g2mu 1 --m 0.2 --seeds 4`, 15 s with
8 workers), D − D_G at r_s falls from 0.09 (`Ny 1`) to 0.002 (`Ny 100`),
`Ny` (D − D_G) ≈ 0.1, and the extrapolation gives (1.2 ± 0.8)·10⁻³.

### Lattice symmetries

`lattice_symmetries/lattice_symmetries.py` (issue #124) tests the index
bookkeeping of the forward light cone, the evolution, T^μν, the Landau
matching, the interpolation to the output grid, the eccentricities and the
gluon spectrum. A Pb+Pb event at b = 0 (the example input without JIMWLK,
`gaussianWounding 0`, `rotateReactionPlane 0`) writes its Wilson lines and
geometry files, which are read back (`readInitialWilsonLines 2`) three
times:

- **(0)** unchanged: this must reproduce the event, to the last bit with a
  build with `-DIPGLASMA_DETERMINISTIC_FFT=ON` (the script copies the FFT
  wisdom file into every run directory);
- **(a)** transposed, x ↔ y (Wilson lines, nucleon positions, g²μ² and T_p
  maps): all outputs map onto those of (0) point by point, with
  T^xx ↔ T^yy, T^τx ↔ T^τy, T^xη ↔ T^yη, u^x ↔ u^y (and likewise π^μν),
  the same |ε_n|, Ψ_n → π/2 − Ψ_n and the same gluon spectrum;
- **(b)** projectile and target exchanged (the reflection η → −η):
  T^τη, T^xη, T^yη, u^η and π^{μη} change sign, everything else is
  unchanged.

Both hold to rounding except where the light-cone solver finds another
solution of the matching condition (issue #107): in (b) it starts from
U_B U_A instead of U_A U_B, and in (a) a link that has to restart draws
another random start, since the restart stream depends on the site and the
link direction. These links lie in the strong-field region, and their
differences spread by the final time. The script therefore reports the
fraction of T^μν points that differ beyond rounding and accepts (a) and (b)
if it is at most `--max-differing` (default 25%; a general error in the
index bookkeeping changes about half of all points, the rest being empty)
and dN/dy, dE/dy, ε₂ and ε₃ change by less than `--observable-tolerance`
(default 3%). With a restart stream that does not depend on the site, (a)
holds everywhere: T^μν to 6·10⁻¹³, the hydro fields (from the Landau
matching) to 2·10⁻⁷.

```bash
cd lattice_symmetries
# quick test: 64², L = 24 fm (about 5 s)
python3 lattice_symmetries.py --ipglasma-path ../.. --quick
# the setup of the paper: 256², L = 30 fm
python3 lattice_symmetries.py --ipglasma-path ../..
```

`--ipglasma-path` is a directory with the `ipglasma` executable,
`qs2Adj_vs_Tp_vs_Y_240.in` and `nucleusConfigurations/`; the input starts
from the example `input` (`--template`). The results go to
`<datadir>/symmetries_summary.json` (default `./symmetries_output`). At
128², L = 24 fm, `maxTime 0.3`: (0) is bit-identical; in (a) 1.4% of the
T^μν points differ and the observables by less than 0.1%; in (b) 6% and
0.3%.

### Dilute limit

`dilute_limit/dilute_limit.py` (issue #125) tests the forward light cone,
the evolution, T^μν and the gluon multiplicity in their perturbative limit.
With a constant color-charge density (`useNucleus 0`) and an IR regulator m
large compared with g²μ, both nuclei enter to second order, so the energy
density and dN/dy scale as (g²μ)⁴/g². The script runs the example input with
`useNucleus 0`, `m 1` and no JIMWLK for a scan of `g2muGeV`, all with the
same seed (the random numbers are then the same, only scaled), plus a run
with `g 2`. It checks that ε at the first time step and τε at the final time
(0.5 fm), averaged over the whole lattice, which for `useNucleus 0` is
periodic, and dN/dy scale with exponent 4 within 1% over the dilute range
(default g²μ = 0.05–0.2 GeV), and that ε g² does not depend on g. It reports
how far the larger g²μ deviate (the onset of saturation) and a fit of the
spectrum, k⁴ dN/d²k = B ln(k/k₀).

**The leading-order spectrum.** The script also compares the gluon spectrum
with its analytic leading-order value, which has no free parameter
(`leading_order.py`):
- **Formula:** to first order in each nucleus, the fields at τ = 0⁺ are the
  longitudinal fields E^η_c = (1/g) f^{abc} ∂_iΛ_A^a ∂_iΛ_B^b and
  B^η_c = (1/g) f^{abc} ε_ij ∂_iΛ_A^a ∂_jΛ_B^b, with
  (−∇² + m²)Λ = gρ. They then evolve as free fields, the E^η mode as
  E_z ∝ J₀(kτ), the B^η mode with a transverse electric field ∝ J₁(kτ).
  Averaged over the color charges, with q = k − p,

  dN/dy d²k = N_c(N_c² − 1)/g² (g²μ)⁴ /(2π)² (τ/k)
  [I_E(k) J₀(kτ)² + I_B(k) J₁(kτ)²] per area,

  I_E,B(k) = ∫ d²p/(2π)² {(p·q)², (p×q)²} / [(p² + m²)²(q² + m²)²].
- **The time of the measurement:** the Bessel functions are there because
  the code counts gluons with the electric fields, n = 2|E|²/ω, which
  oscillate with kτ. At late times J₀², J₁² → 1/(πkτ) on average, and the
  formula becomes the classical result of Kovner, McLerran and Weigert and of
  Gyulassy and McLerran, N_c(N_c² − 1)/g² (g²μ)⁴ (I_E + I_B)/(4π³k²), which
  falls like ln(k/m)/k⁴ for k ≫ m.
- **Lattice:** the formula is binned like the spectrum of the code (the
  modes of the periodic lattice in bins of the lattice momentum). It has no
  lattice artifacts, so it holds for ka ≪ 1.
- **Configurations:** the smallest g²μ is run for `--lo-seeds` (16)
  color-charge configurations, whose spectra are averaged.
- **Check:** the mean ratio of the lattice to the leading-order spectrum in
  the window 2m < k < 0.3/a must be 1 within `--lo-tolerance` (2%) plus
  two standard deviations.

```bash
cd dilute_limit
# quick test: 64², five g²μ (a few seconds)
python3 dilute_limit.py --ipglasma-path ../.. --quick --max-workers 6
# the setup of the paper: 256², L = 3.2 fm (ħc/a = 16 GeV), m = 1 GeV
python3 dilute_limit.py --ipglasma-path ../.. --max-workers 8
```

With the defaults (2 minutes with 4 workers of 2 threads), the exponents are
3.9997 ± 0.0002 (ε at τ → 0⁺), 4.014 ± 0.008 (τε) and 3.991 ± 0.002
(dN/dy). At g²μ = 0.5, 1 and 2 GeV, τε falls 2%, 6% and 15% below the
dilute scaling, and dN/dy 2.5%, 9% and 28%. In the window 2 GeV < k < 4.7 GeV
the spectrum follows B ln(k/k₀) with k₀ = 1.5 m, within 8%. Against the
leading order, with 16 configurations at g²μ = 0.05 GeV:

| k [GeV] | k a | lattice / leading order |
|---|---|---|
| 0.3–2 | < 0.13 | 1.045 ± 0.020 |
| 2–4.7 (the window) | 0.13–0.3 | 1.005 ± 0.006 |
| 4.7–9.5 | 0.3–0.6 | 1.026 ± 0.003 |
| 9.5–19 | 0.6–1.2 | 1.041 ± 0.002 |

- **The normalization agrees within 0.5% ± 0.6%** in the window, with no
  free parameter. This tests every factor from the color charges to the
  gluon number.
- **Without the Bessel functions,** i.e. with the late-time formula, the
  ratio oscillates between 0.87 and 1.32 below k ≈ 4 GeV, where kτ is only a
  few units at τ = 0.5 fm.
- **Above the window,** the lattice spectrum lies a few percent higher. These
  are lattice artifacts; e.g. the estimator divides by the lattice ω = k̃ < k,
  which alone gives about (ka)²/24.
- **dN/dy** is 1.050 ± 0.024 times the leading-order sum over the lattice
  modes; the few lowest modes, at k ≈ 2π/L, carry most of it.

Before #139 the lattice of `useNucleus 0` had clamped edges. These turned
the pure-gauge fields of the nuclei, of first order in g²μ, into physical
fields in a band near the edges. In the dilute limit that band carried 99.9%
of the energy at g²μ = 0.05 GeV, and dN/dy of the whole lattice scaled as
(g²μ)². On the periodic lattice the whole lattice follows the dilute
scaling.

The dilute range is bounded below by numerics. For g²μ a and m a small
enough (e.g. g²μ = 0.01 GeV at a = 0.0125 fm), the field produced in the
collision falls below what the light-cone solver resolves with its tolerance
(ε at τ → 0⁺ comes out 13% low), and the late energy density below the floor
of the T^μν output (10⁻¹⁶ GeV/fm³).

### Dependence on a, Ny and m at a fixed geometry

`numerical_parameters/numerical_parameters.py` (issue #127) runs the example
input without JIMWLK (b = 0, L = 22.5 fm, which suffices at b = 0, `seed 3`)
for three scans: the size (a = L/size, by default 0.234 to 0.018 fm for the
sizes 96 to 1280) at the `Ny` and `m` of the input, and `Ny` and `m` at the
size 768 (a = 0.029 fm). For another impact parameter, use L = 30 fm. With a fixed `seed` and b = 0 (or a multiple of L/128), the geometry
(nucleon positions, hot spots, T_p) does not depend on `size`, `Ny` or `m`;
only the color charges change. These are varied with the input parameter
`colorChargeSeed`, which draws them from their own stream, so that a few
seeds per setting give its statistical error. The script computes dE/dη =
τ∫T^ττ d²x at τ = 0.2 fm and at the final time (0.4 fm), dN/dy, ε₂, ε₃ and
their angles, and reports their mean and error against a, `Ny` and `m`
(`numerical_parameters_summary.json` and `numerical_parameters.pdf` in
`--datadir`, default `./numerical_parameters_output`). It checks that the
geometry is the same in all runs (N_part and N_coll lists, b) and that the
color-charge seeds give different events. The change of the physical Q_s
with `Ny` and `m` at fixed `QsMuRatio` is measured by the dipole setup.

The lattice converges once it resolves the color-charge scale g²μ =
Q_s/`QsMuRatio`, g²μ a well below 1. A run that only builds the Wilson
lines (at the size of the `Ny` and `m` scans) gives g²μ of each site of the
geometry; the script reports its median (weighted with (g²μ)², like the
multiplicity) and maximum, and the median g²μ a of each size. The hot spots
of the example input are small (`BGq 0.3` GeV⁻², a width of 0.11 fm) and
overlap in Pb to g²μ of up to 13 GeV; for `seed 3` the median is 5.2 GeV
(1/g²μ = 0.038 fm), i.e. g²μ a = 6.2, 3.1, 1.5, 0.77, 0.58 and 0.46 for the
default sizes, so that only the finest lattices resolve most of the event. As a
reference, a fourth scan runs a constant g²μ (`useNucleus 0`, on its
periodic lattice, `--uniform-g2mu`, default 5 GeV, L = 3 fm) at the `Ny`
and `m` of the input, for g²μ a = 2.4 to 0.15 (sizes 32 to 512):

| size | 32 | 64 | 128 | 256 | 512 |
|---|---|---|---|---|---|
| g²μ a | 2.38 | 1.19 | 0.59 | 0.30 | 0.15 |
| dN/dy | 3616 | 4687 | 5061 | 5260 | 5133 |
| dE/dη (0.4 fm) [GeV] | 12147 | 22022 | 28936 | 33240 | 34810 |

dN/dy is within about 4% of its limit for g²μ a ≲ 0.6, while dE/dη still
grows by 4.7% ± 1.4% from g²μ a = 0.3 to 0.15. The energy sits in modes up
to several g²μ, where the gluon spectrum falls like 1/k⁴, and converges
only like a: its deficit is about 0.3 g²μ a, so 5% needs g²μ a ≈ 0.15 (the
energy density at τ = 0 even diverges like ln² a, T. Lappi, Phys. Lett. B
643 (2006) 11). This is not an artifact of IP-Glasma: an independent
classical Yang-Mills code (a version of S. Schlichting's), run on the same
Wilson lines, gives dE/dη within 0.1–0.6% of IP-Glasma's and the same
dependence on a (see the discussion of #134). With JIMWLK, the gluon
spectrum falls less steeply (the anomalous dimension of the evolved Wilson
lines), so the energy converges even more slowly.

**Convergence check.** For each observable of the size scans, the script
compares the two finest lattices: the relative change must be below
`--convergence-tolerance` (1%) plus two standard deviations. For an
observable that fails, it prints a prominent `WARNING: NOT CONVERGED IN THE
LATTICE SPACING` with the change and an estimate of the remaining deviation
from an extrapolation linear in a, marks the panel in the figure, and lists
it under `not_converged` in the summary. With the defaults (25 minutes with
2 workers of 8 threads), the results are

| size | 96 | 192 | 384 | 768 | 1024 | 1280 |
|---|---|---|---|---|---|---|
| a [fm] | 0.234 | 0.117 | 0.059 | 0.029 | 0.022 | 0.018 |
| dN/dy | 8195 | 14375 | 18994 | 21418 | 21766 | 21931 |
| dE/dη (0.4 fm) [GeV] | 9062 | 28602 | 51825 | 68667 | 72709 | 75766 |

dN/dy changes by +0.8% ± 2.0% in the last step, ε₂ and ε₃ by less than
their errors, so they converge. dE/dη still grows by 4.2% ± 1.2%, and the
linear extrapolation puts its value at a = 0 about 16% above the finest
lattice (5% for the constant g²μ): converging the energy to the percent
level would need lattice spacings far below 0.018 fm, so its dependence on
a has to be quoted as a systematic uncertainty. The `Ny` scan (at 768²)
converges for dN/dy by `Ny` ≈ 25–50 (21418 for `Ny 50`, 21174 for
`Ny 100`); the `m` scan changes dN/dy from 28825 (`m 0.1`) to 15675
(`m 0.8`) at fixed `QsMuRatio`, which changes the physical Q_s.

```bash
cd numerical_parameters
# quick test: sizes 64 and 128, Ny and m scans at 64, the constant g²μ at
# sizes 32 and 64, 2 seeds (6 s with 8 workers)
python3 numerical_parameters.py --ipglasma-path ../.. --quick --max-workers 8
# the setup of the paper (the defaults): sizes 96-1280, Ny 1-100,
# m 0.1-0.8 GeV, the constant g²μ at sizes 32-512, 4 seeds (25 minutes)
python3 numerical_parameters.py --ipglasma-path ../.. --threads 8 \
    --max-workers 2
```

**Where the energy sits in k_T.** For each lattice of the size scans, the
script also reads the gluon spectrum (`gluonMultiplicity0.json`) and
reports, in the summary and in a table:
- **The energy split:** the energy below and above 2 and 4 times the g²μ
  scale (the median of the geometry, or the constant g²μ). The part below
  is the binned energy, and the part above is the total dE minus it, so it
  includes the modes beyond the largest bin.
- **A continuum estimate:** the binned energy below k_c plus a tail
  k⁴ dN/d²k = C beyond it, 2πC/k_c, for k_c = 0.3 and 0.4 times the largest
  k_T. It is meaningful only once k_c lies in the 1/k⁴ tail, i.e. on fine
  lattices.
- **The tail exponent:** p of k⁴ dN/d²k ∝ k^p, fitted between 6 and 12 GeV
  (`--tail-fit`); p = 0 is the 1/k⁴ tail.

The mean spectra are in the summary, and `numerical_parameters_spectra.pdf`
shows k⁴ dN/d²k of each lattice. For the constant g²μ, the energy below
2g²μ does not grow on finer lattices; all of the growth is above. The two
estimates agree on the two finest lattices, which puts the deficit at about
0.3 g²μ a (discussion of #134).

**With JIMWLK, and at another x.**
- **`--jimwlk`:** evolves the Wilson lines of all runs with the JIMWLK
  parameters of the template, i.e. the example input from x = 0.01 to 10⁻⁴.
  The g²μ scale is then that of the initial Wilson lines.
- **`--x X`:** sets `projectileX` and `targetX`. Without `--jimwlk`, this is
  the x of the Wilson lines; with it, the x the evolution ends at.
- **Separate results:** runs and summaries of these variants are tagged
  (`jimwlk_`, `x0.01_`, `numerical_parameters_jimwlk_summary.json`), so they
  never mix with the default ones.
- **Skipping scans:** `--ny` and `--m` with no values skip those scans, and
  `--uniform-sizes` with no values the constant g²μ, which does not depend
  on x.

In the discussion of #134 (τ = 0.5 fm, 2 seeds, N = 768 to 1536):
- **With JIMWLK:** dN/dy converges (+1.4% in the last step). dE/dη grows by
  about 12% per step, because the evolution hardens the tail to p ≈ 1.3.
- **At x = 0.01 without JIMWLK:** p ≈ 0.3, and dE/dη is converged within
  about 1% from N = 1024 on.

```bash
# the size scan with JIMWLK, and the Wilson lines at x = 0.01 without it
python3 numerical_parameters.py --ipglasma-path ../.. --jimwlk \
    --sizes 768 1024 1280 1536 --ny --m --uniform-sizes --color-seeds 2
python3 numerical_parameters.py --ipglasma-path ../.. --x 0.01 \
    --sizes 768 1024 1280 1536 --ny --m --uniform-sizes --color-seeds 2
```

`--ipglasma-path` is a directory with the `ipglasma` executable,
`qs2Adj_vs_Tp_vs_Y_240.in`, `nucleusConfigurations/` and `tables/`. Runs are
cached in their directories, so a rerun only runs what is missing.

### Constant-g²μ glasma

`constant_glasma/constant_glasma.py` (issue #126) tests the absolute
normalization and the time dependence of the classical Yang-Mills part
against T. Lappi, Phys. Rev. C 67 (2003) 054903 [hep-ph/0303076]: SU(3), a
constant color-charge density (`useNucleus 0`), `Ny 1` and `m 0`, on the
lattices of Lappi's Fig. 3 (N = 128, 256 and 512 at fixed g²μL, so
g²μa = g²μL/N). T^μν and the gluon multiplicity are taken from the whole
lattice, which for `useNucleus 0` is periodic, as Lappi's. The script
computes, averaged over the seeds,

- f_E = g²τε/(g²μ)³ against g²μτ; its plateau is the value at the last
  time, g²μτ = 20. After its maximum at g²μτ ≈ 1, τε decreases slowly
  (d(τε)/dτ = −P_L, and P_L is slightly positive). Lappi's Table I gives the
  asymptote of a fit A + B exp(−τ/τ₀), which his Fig. 4 reaches at
  g²μτ ≈ 20 (τ = 2 fm there);
- P_L/ε = τ²T^ηη/T^ττ and P_T/ε = (T^xx + T^yy)/(2T^ττ) against τ;
- ε(τ → 0⁺), at the first time step, against the lattice spacing;
- f_N = g²(dN/dy)/(g²μL)², the gluon number per unit rapidity and
  transverse area, at the last time. IP-Glasma computes the multiplicity in
  transverse Coulomb gauge with the free lattice dispersion, as Eqs. (27)
  and (28) of Lappi, the definition of his Table I.

It checks that P_L/ε = −1 and P_T/ε = 1 within 0.02 at the first time step
(purely longitudinal fields) and about 0 and 1/2 at the last time
(`--pressure-tolerance` 0.1), that f_E changes by less than 10% between half
the last time and the last time, that ε(0⁺) grows as a → 0 while the plateau of f_E changes by
less than 10% (beyond two standard deviations) between the two finest
lattices (T. Lappi, Phys. Lett. B 643 (2006) 11 [hep-ph/0606207]), and that
on the 512² lattices of Lappi's Table I (g²μL = 72, 120, 192) the plateau of
f_E (0.265 ± 0.005, 0.227 ± 0.003, 0.238 ± 0.005) and f_N (0.297 ± 0.006,
0.289 ± 0.003, 0.329 ± 0.006) agree with his values within three combined
standard deviations. The f_E of Krasnitz, Nara and
Venugopalan, Phys. Rev. Lett. 87 (2001) 192302 [hep-ph/0108092], as quoted
by Lappi and halved as in their erratum, Nucl. Phys. A 727 (2003) 427
[hep-ph/0305112] (0.269 at g²μL = 44.3, 0.248 at g²μL = 148.4), are printed
for information; their lattice spacing is not quoted, so they are not
checked.

```bash
cd constant_glasma
# quick test: g²μL = 36 on 128² and 256², four seeds (1 minute)
python3 constant_glasma.py --ipglasma-path ../.. --quick --max-workers 4
# the setup of the paper: g²μL = 72 and 120 on 128², 256² and 512²
python3 constant_glasma.py --ipglasma-path ../.. --max-workers 4
```

With the defaults (four seeds, 9 minutes with 2 workers of 4 threads), all
checks pass. f_E and f_N agree with Lappi's Table I on 512², and the plateau
of f_E with the points of his Fig. 3 (read off the figure) on the coarser
lattices:

| g²μL | N   | g²μa | f_E           | Lappi         | f_N           | Lappi         |
|------|-----|------|---------------|---------------|---------------|---------------|
| 72   | 128 | 0.56 | 0.246 ± 0.011 | ≈ 0.237       | 0.319 ± 0.013 |               |
| 72   | 256 | 0.28 | 0.249 ± 0.008 | ≈ 0.265       | 0.298 ± 0.011 |               |
| 72   | 512 | 0.14 | 0.264 ± 0.010 | 0.265 ± 0.005 | 0.305 ± 0.014 | 0.297 ± 0.006 |
| 120  | 128 | 0.94 | 0.199 ± 0.008 | ≈ 0.193       | 0.297 ± 0.008 |               |
| 120  | 256 | 0.47 | 0.217 ± 0.008 | ≈ 0.217       | 0.297 ± 0.009 |               |
| 120  | 512 | 0.23 | 0.231 ± 0.006 | 0.227 ± 0.003 | 0.300 ± 0.008 | 0.289 ± 0.003 |

g²ε(0⁺)/(g²μ)⁴ grows from 0.62 to 0.81 and 1.06 (g²μL = 72) and from 0.41
to 0.58 and 0.77 (g²μL = 120) as a is halved, while f_N hardly depends on
the lattice spacing. P_L/ε is between −0.995 and −0.988 at the first time
step and +0.015 to +0.025 at g²μτ = 20 (Lappi's Fig. 5 gives about +0.025),
and f_E decreases by 2.5–2.7% between g²μτ = 10 and 20. The errors are
those of the mean over four seeds; more seeds are needed to reach the
precision of Lappi's errors.

On a lattice with clamped boundaries, as `useNucleus 0` had before #139,
the pure-gauge fields of the nuclei become physical fields along the
edges; averaged over the whole lattice, f_E and f_N on 512² were then 20–43%
too large.

### Performance

`performance/performance.py` measures the time per stage of an event and
the peak memory, against the lattice size and the number of OpenMP
threads:
- **Event:** the example input (Pb+Pb with hot spots, JIMWLK from
  x = 0.01 to 10⁻⁴, the Yang–Mills evolution to τ = 0.5 fm, the gluon
  multiplicity and the eccentricities) at a fixed geometry (`seed 3`,
  b = 0), without writing Wilson lines or fields.
- **Scans:** N = 256 to 1536 at L = 30 fm with 8 threads, and 1 to 16
  threads at N = 512. With fixed L, the lattice spacing shrinks with N, so
  the Yang–Mills evolution, whose time step is a/10, needs N times more
  steps.
- **Times:** from the profile IP-Glasma writes with `IPGLASMA_PROFILE=1`
  (`ipglasma_profile_rank0.tsv`, see `OUTPUT.md`): the Wilson lines
  (`initialization.wilson_lines`), the JIMWLK evolution
  (`jimwlk.evolution`), the light-cone fields, the Yang–Mills evolution
  with the observables at the final time (`evolution.total`), and the whole
  event. The summary has every phase of the profile, with its number of
  calls.
- **Time per step:** of the JIMWLK evolution (`jimwlk.step`) and of the
  Yang–Mills evolution (`evolution.parallel_step`), the time of the phase
  divided by its calls. The number of JIMWLK steps does not depend on N
  (93 per nucleus), that of the Yang–Mills steps grows like N.
- **Memory:** the peak resident memory of each run, from the operating
  system.
- **Repeats:** the runs are sequential and repeated (`--repeats`, default
  2), and the fastest repeat is kept.

The timings depend on the machine and on what else runs on it: run the
script on an otherwise idle machine. The summary records the CPU model, the
number of CPUs and the load average, and the script warns if the load is
above 1.

```bash
cd performance
# quick test: N = 64 and 128, 1 and 2 threads, one repeat (20 s)
python3 performance.py --ipglasma-path ../.. --quick
# the setup of the paper (about 1.5 hours)
python3 performance.py --ipglasma-path ../..
```

The results go to `performance_summary.json` and `performance.pdf` in
`--datadir` (default `./performance_output`).

### JIMWLK step size

`jimwlk_step_size/jimwlk_step_size.py` (issue #128) tests the convergence of
the Langevin step of the JIMWLK evolution in `jimwlkDs`, and the error at
the step sizes recommended in the main README (0.005 with running and 0.0005
with fixed coupling). The MV model with a constant color-charge density
(`useNucleus 0`, g²μ = 2 GeV, `Ny 50`, m = 0.4 GeV, 512², L = 10 fm) is
evolved from x₀ = 0.01 with `jimwlkMass 0.4`:

- running coupling (`jimwlkAlphaS 0`, the other coupling parameters of the
  example input): `jimwlkDs` = 0.01, 0.005, 0.0025, 0.00125, to
  s = ln(x₀/x)/π² = 0.48 (y = 4.74; 48 to 384 steps);
- fixed coupling (`jimwlkAlphaS 0.2`): `jimwlkDs` = 0.002, 0.001, 0.0005,
  0.00025, to s = α_s ln(x₀/x)/π² = 0.04 (y = 1.97; 20 to 160 steps).

Q_s grows from 1.05 to 1.9 GeV (running) and 1.4 GeV (fixed), so
1/L ≪ Q_s ≪ 1/a holds (r_s ≥ 7.5 a). The evolution time and the three
snapshots in between are multiples of the largest step, so all step sizes
end and take their snapshots at exactly the same x. All step sizes use the
same seeds, so they start from the same initial Wilson lines (checked), and
both nuclei of a run are evolved, two configurations per run.

From the Wilson lines the script computes the dipole correlator D(r) (with
the functions of `dipole_correlator/`) and Q_s from D(r_s) = e^{−1/2},
Q_s² = 2/r_s². It fits the differences to the smallest step,
Q_s(ds) − Q_s(ds_min), which share the initial Wilson lines, with their
covariance over the configurations: a first-order scheme gives
c (ds − ds_min), the χ² of this fit tests it, and Q_s(0) = Q_s(ds_min) −
c ds_min is the extrapolation to zero step size. A fit c (ds^p − ds_min^p)
gives the order p. The same is done for D at the r* where D = e^{−1/2} at
the smallest step. At the final x it checks that the initial Wilson lines
agree for all step sizes, that Q_s and D converge monotonically (every
successive difference has the sign of c or is zero within 2σ), that they are
consistent with first order (p-value above 0.01) and that the deviation
c ds/Q_s(0) at the recommended step size is below 3% (`--max-deviation`).

Each step size draws its own Langevin noise, so the differences between step
sizes fluctuate by 1–2% per configuration even at 512², and many
configurations are needed: with four seeds, the target of one step size
happened to lie 2% low in three of them.

```bash
cd jimwlk_step_size
# quick test: 128² with the same lattice spacing, four seeds (1 minute)
python3 jimwlk_step_size.py --ipglasma-path ../.. --quick --max-workers 4
# the setup of the paper: 512², eight seeds (16 configurations)
python3 jimwlk_step_size.py --ipglasma-path ../.. --max-workers 5
```

With the defaults (37 minutes with 5 workers of 4 threads), all checks
pass. The step is of first order (fitted order 0.90 [0.68, 1.15] for Q_s
with running and 1.00 [0.49, 1.59] with fixed coupling, χ² p-values above
0.5), and a finite step makes the evolution too fast. At the final x, the
deviation from zero step size, c ds/Q_s(0), is

| coupling | jimwlkDs | Q_s            | D(r*)           |
|----------|----------|----------------|-----------------|
| running  | 0.01     | +3.8%          | −2.1%           |
| running  | 0.005    | +1.9% ± 0.15%  | −1.0% ± 0.08%   |
| running  | 0.0025   | +0.9%          | −0.5%           |
| fixed    | 0.002    | +1.8%          | −1.0%           |
| fixed    | 0.0005   | +0.46% ± 0.05% | −0.26% ± 0.03%  |
| fixed    | 0.00025  | +0.23%         | −0.13%          |

(the measured deviations, e.g. +2.0% ± 0.3% and +0.6% ± 0.4% at the
recommended steps, agree with the fit). The error grows about linearly with
the rapidity: at `jimwlkDs 0.005` it is +0.3%, +0.7%, +1.1% and +1.9% at
y = 1.2, 2.4, 3.6 and 4.7, about 0.4% per unit of y, or 3% of the evolution
speed d ln Q_s/dy (1.6% with fixed coupling at 0.0005). The recommended step
sizes stay below the 3% of the acceptance over the evolution of the example
input (from x = 0.01 to 10⁻⁴, y = 4.6).

The error depends strongly on the setup. Two more scripts measure it in
other setups.

**The setup of arXiv:1212.4825.** `lappi_mantysaari.py` runs the setup of
Lappi and Mäntysaari, whose curves are in `lappi_mantysaari_2012/` (see the
README there for how they were extracted):
- **Initial condition:** the MV model with g²μL = 31, `Ny 100` and no
  infrared regulator, on 1024².
- **Evolution:** the square-root running coupling, the prescription
  IP-Glasma implements (the noise prescription of the paper is #138), with
  Λ_QCD L = 6, μ₀L = 15, c = 0.2, N_f = 3 and `jimwlkMass 0`. The coupling
  of the paper has the scale 4e^{−2γ_E}/r² instead of 4/r², so
  `jimwlkLambdaQCD` and `jimwlkMu0` are e^{γ_E} times its Λ_QCD and μ₀; the
  coupling freezes at 0.76, against 0.36 in the example input.
- **Runs:** `jimwlkDs` 0.001, 0.0005 and 0.00025 from the same initial Wilson
  lines (two seeds, four configurations), to y = 8, with snapshots at every
  unit of y.
- **Analysis:** at each y, ln Q_s is fitted against the step as above and
  extrapolated to zero step; the evolution speed λ = d ln Q_s²/dy comes from
  successive snapshots.
- **BK:** `bk_sqrt.py` solves the BK equation with the same kernel, the
  mean-field limit of this JIMWLK evolution, from the mean dipole
  correlator of the initial Wilson lines (RK4 on a logarithmic grid in r,
  converged to 0.2%).
- **Checks:** the initial Wilson lines agree, the step is of first order,
  BK reproduces λ of Fig. 11 of the paper within 2%, JIMWLK at zero step is
  slower than BK, and Q_s at zero step agrees with Fig. 6 within 10%.

```bash
cd jimwlk_step_size
# quick test: 256², to y = 4, jimwlkDs 0.002, 0.001, 0.0005 (4 minutes)
python3 lappi_mantysaari.py --ipglasma-path ../.. --quick --max-workers 3
# the setup of the paper (several hours)
python3 lappi_mantysaari.py --ipglasma-path ../.. --threads 4 --max-workers 6
```

With the defaults, the step is again of first order, but its error per unit
of step and of rapidity is about 30 times larger than above:

| Q_s/Λ_QCD at y = 8 | |
|---|---|
| `jimwlkDs 0.001` | 33.2 |
| `jimwlkDs 0.0005` | 30.1 |
| `jimwlkDs 0.00025` | 28.7 |
| zero step | 27.3 |
| BK | 30.2 |
| Fig. 6 of the paper | 28.5 |

- **The paper matches a step of about 0.00025:** Q_s within 0.5% of its
  Fig. 6 at y = 4 and 8, and λ along its Fig. 8.
- **BK** reproduces λ of Fig. 11 within 0.5%.
- **At zero step,** λ of JIMWLK is 4–9% below BK for Q_s/Λ_QCD = 10–25, so
  JIMWLK is slower than BK, as in the paper (where the difference is about
  5% in Q_s at y = 7.8, here 9%).

**The example input.** `example_input_step_size.py` evolves the example
input itself (Pb+Pb with hot spots, b = 0) from x₀ = 0.0098 to 10⁻⁴, exactly
93 steps of 0.005, at `jimwlkDs` 0.005, 0.0025, 0.00125 and 0.000625, from
the same initial Wilson lines:
- **Q_s:** a nucleus is not translation invariant, so Q_s comes from the
  dipole correlator averaged over a region: the thirds of each nucleus by
  the initial g²μ (dilute, medium, dense) and all of them. The dilute
  third includes the edges and the gaps between the hot spots, where the
  coupling is strongest.
- **Observables:** dE/dη, dN/dy, ε₂ and ε₃ after the Yang–Mills evolution.
- **Analysis:** a fit linear in the step to the averages over the
  configurations. All step sizes start from the same initial Wilson lines,
  so the error comes from the spread of the same fit over the
  configurations: the fluctuations of the initial condition cancel in it,
  while those of the Langevin noise, which differs between step sizes,
  remain in the error. A quadratic fit tests the first order.
- **Checks:** the initial Wilson lines agree, the step is of first order,
  Q_s of every region and dN/dy are within 3% of zero step at 0.005, and
  dE/dη at 0.0025.

```bash
# quick test: 128², L = 30 fm, three step sizes, 2 seeds (20 s)
python3 example_input_step_size.py --ipglasma-path ../.. --quick --max-workers 3
# 768², L = 22.5 fm, 3 seeds (about an hour with 2 workers of 4 threads)
python3 example_input_step_size.py --ipglasma-path ../.. --max-workers 2
```

With the defaults, the deviation from zero step at `jimwlkDs 0.005` is:

| | 768², L = 22.5 fm | 256², L = 30 fm (the example input) |
|---|---|---|
| Q_s, all | +0.9% ± 0.6% | −0.4% ± 0.5% |
| Q_s, dilute third | +1.2% ± 0.7% | −0.4% ± 0.8% |
| Q_s, dense third | +0.5% ± 0.8% | −0.4% ± 0.4% |
| dN/dy | +2.3% ± 1.3% | +2.3% ± 3.9% |
| dE/dη | +5.4% ± 0.6% | +0.2% ± 1.6% |

- **Q_s and dN/dy:** within 1–2% of zero step at 0.005, as in the
  constant-g²μ setup above.
- **dE/dη:** its step error grows on finer lattices. It sits in the hard
  tail of the evolved spectrum, which the coarse lattice cuts off (the same
  tail makes the energy converge slowly in the lattice spacing, see
  `numerical_parameters/`). At 0.0025 it is 2.7% ± 0.3% at 768².
- **ε₂ and ε₃:** consistent with no dependence on the step, within their
  errors (about 6% and 11% at 768²).

### Vector meson production validation

```bash
cd vector_meson_production
python3 parallel_test_vector_meson_production.py \
    --maxevents 200 \
    --datadir ./jpsi_test_run \
    --ipglasma-path ../.. \
    --subnucleondiffraction-path /path/to/subnucleondiffraction \
    --max-workers 8
```

Useful flags:
- `--input-template FILE` — IP-Glasma input file to use as a template
  (default `input_vm_proton`); the Q_s table (`nucleusQsTableFileName`) is
  taken from this file rather than hardcoded.
- `--xpom X` — Bjorken x the Wilson lines are evolved to, written to the
  input as `projectileX` and `targetX` (default 0.001705).
- `--max-workers N` — number of events to run in parallel (default: the
  number of logical CPUs).
- `--plot-only` — skip the simulation and only regenerate the comparison
  plot from data files already in `--datadir`.
- `--keep-logs` — keep per-seed log files instead of deleting them.

Output: `cross_section.dat`/`.pdf` in `--datadir` (coherent/incoherent cross
section vs. reference, with a ratio panel), plus `output_1.pdf`/`output_5.pdf`
color-field (Tr V(x)) plots for a couple of representative events.

## Running on a computing cluster (SLURM)

`vector_meson_production/runvalidation_jpsi.sh` is a SLURM batch script
template for `parallel_test_vector_meson_production.py`:

```bash
#!/bin/bash
#SBATCH --job-name=validation
#SBATCH --account=<your_account>
#SBATCH --partition=small
#SBATCH --time=2:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000M

# Each event runs single-threaded; the workers run the events in parallel.
# For n threads per event, set --cpus-per-task=n (and --ntasks to the number
# of events run at a time), OMP_NUM_THREADS=n, OMP_PROC_BIND=close instead of
# spread, and OMP_WAIT_POLICY=active: with --slurm, each event then runs as
# an srun step with n CPUs (see "Threads and MPI ranks" in the README).
export OMP_NUM_THREADS=1

# Place and bind threads to single cores
# Comment the following lines if binding is not desired
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
module add python-data gsl fftw

taskset -pc $$

# Run the program
python3 parallel_test_vector_meson_production.py --slurm \
    --max-workers "${SLURM_NTASKS:-1}" \
    --datadir /path/to/datadir \
    --maxevents 1500 \
    --subnucleondiffraction-path /path/to/subnucleondiffraction \
    --ipglasma-path /path/to/ipglasma
```

Before submitting, edit the placeholders (`--account`, `--datadir`,
`--subnucleondiffraction-path`, `--ipglasma-path`, and the `module add`
line for your cluster) and adjust
`--ntasks`/`--maxevents`/`--time` to the scale of the run. Submit it from
its folder, where the job then runs:

```bash
cd vector_meson_production
sbatch runvalidation_jpsi.sh
```

Notes:
- Each `--ntasks` maps to one `--max-workers` process. With
  `--cpus-per-task=1` and `OMP_NUM_THREADS=1`, each event runs
  single-threaded and the parallelism comes from running many events at
  once; for several threads per event, see the comment in the template.
- Both scripts accept a `--slurm` flag that wraps each
  `ipglasma`/`subnucleondiffraction` invocation with
  `srun --ntasks=1 --cpus-per-task=N`, with N = `OMP_NUM_THREADS` (default
  1), which can help binding/placement on some clusters. The templates pass
  it; remove `--slurm` from the `python3` invocation if your cluster runs
  the binaries directly.
- `eccentricity_regression/runvalidation_eccentricity.sh` is the same
  template for `parallel_test_eccentricity.py` (submit it from its folder
  with `sbatch runvalidation_eccentricity.sh`); only the job name and the
  python invocation differ:
  ```bash
  python3 parallel_test_eccentricity.py --size 128 \
      --ipglasma-path /path/to/ipglasma \
      --max-workers "${SLURM_NTASKS:-1}" --maxevents 1500 \
      --datadir /path/to/datadir --input-template input_eccentricity --slurm
  ```
- Expect one `seed_<n>/` working directory per event under `--datadir` while
  the job runs; both scripts clean these up at the end unless `--keep-logs`
  is passed, so pass `--keep-logs` if you want to inspect individual events
  afterwards (or need `--plot-only`/`--clean` to resume from an interrupted
  job).

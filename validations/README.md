# Running the validation scripts

This folder contains, each in its own subfolder, two regression tests for
IP-Glasma and the small validation setups of the code paper, which test
parts of the code against exact or published results:

- **`vector_meson_production/`** (`parallel_test_vector_meson_production.py`)
  — runs many independent
  IP-Glasma + JIMWLK events, feeds the resulting Wilson lines to the external
  [`subnucleondiffraction`](https://github.com/hejajama/subnucleondiffraction)
  code, and compares the resulting coherent/incoherent J/ψ cross section to a
  stored reference (`jpsi_reference.json`) and to Fig. 1 ("CGC+shape fluct") of
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
  `g2muGeV`, `m` and `Ny`, optionally also of `UVDamp`, and the saturation
  scale against Lappi (2008) (see below).
- **`lattice_symmetries/`** — two exact symmetries of the lattice
  formulation, the transposition x ↔ y and the exchange of projectile and
  target, applied to the Wilson lines of a Pb+Pb event (see below).
- **`dilute_limit/`** — the dilute limit of the glasma: the (g²μ)⁴/g²
  scaling of the energy density and the gluon number (see below).
- **`gluon_spectrum_estimator/`** — the two definitions of the gluon number
  (`gluonSpectrumEstimator`) in the example input (see below).
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
   It reads the Wilson lines of IP-Glasma 2.0 from its commit 4e3a05c on
   (`master` since October 2026). The script computes the cross sections
   from the amplitudes itself; it does not use `python/cross_section.py`.
3. Python 3 with `numpy` and `matplotlib` (and optionally `scipy`, used for
   the two-sample KS test in the eccentricity script if available; the
   JIMWLK step-size script, the dilute limit and the `--lappi` mode of the
   dipole correlator need it).
4. Or the Docker image of `docker/Dockerfile` (see the main README), which
   has all of these: the repository with the built `ipglasma` in
   `/home/ipglasma` (the `--ipglasma-path`), subnucleondiffraction in
   `/home/subnucleondiffraction` and the Python packages. It has to be built
   from a branch that contains this folder, e.g. `devel`. A script runs in
   it with its outputs in the mounted directory (owned by the user given
   with `-u`), e.g.
   ```bash
   # from the repository root
   docker build -t ipglasma --build-arg IPGLASMA_BRANCH=devel docker/
   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" -w /work \
       ipglasma python3 \
       /home/ipglasma/validations/dipole_correlator/dipole_correlator.py \
       --ipglasma-path /home/ipglasma --quick --datadir /work/dipole
   ```

## Running interactively

From the folder of each script (paths below assume the scripts are run
from there; both scripts also accept absolute paths for `--ipglasma-path`
etc. if you invoke them from elsewhere; their input templates and
references are found next to the script):

### Eccentricity regression test (quick sanity check)

```bash
cd eccentricity_regression
# One-time: record a baseline from a known-good build (the stored
# reference has 4998 events of the seeds 0-4999)
python3 parallel_test_eccentricity.py --input-template input_eccentricity \
    --ipglasma-path ../.. --maxevents 5000 --save-reference

# Later, e.g. after changing the code: check nothing drifted (1000 events,
# about 4 minutes on 16 cores)
python3 parallel_test_eccentricity.py --input-template input_eccentricity \
    --ipglasma-path ../..
```

Useful flags:
- `--datadir DIR` — where per-event working directories and plots go
  (default `./eccentricity_test_output/`).
- `--size`, `--L`, `--maxtime` — override the lattice size / box size / max
  evolution time from the input template (defaults to a small, fast 128²
  lattice).
- `--jimwlk` — keep JIMWLK evolution enabled (slower, more complete physics
  test); disabled by default for speed.
- `--maxevents N` — number of events (default 1000).
- `--max-workers N` — number of events to run in parallel (default: the
  logical CPUs available, e.g. those of a Slurm allocation, divided by
  `OMP_NUM_THREADS`); each event runs with `OMP_NUM_THREADS` threads
  (default 1). Without `--slurm`, several events at a time run with
  `OMP_PROC_BIND=false`: binding (`OMP_PROC_BIND`, `OMP_PLACES`) would put
  the threads of all events on the same first cores.
- `--plot-only` — skip launching new runs and just rebuild the plot/summary
  from `seed_*/` directories already present in `--datadir` (useful after an
  interrupted run, which keeps them, or with `--keep-logs`).
- `--clean` — remove stray temporary files left behind by an interrupted run
  from every `seed_*/` directory under `--datadir`, without touching the
  saved `Tmunu` snapshots, then exit.
- `--keep-logs` — keep the per-event working directories after the run
  instead of deleting them (failed events are always kept); they hold the
  `Tmunu` snapshots, about 0.6 MB per event.
- `--first-seed S` — the events use the seeds S, S + 1, ...: by default 10000,
  so that a check is independent of the reference, and 0 with
  `--save-reference`.
- `--n-sigma`, `--bootstrap` — the allowed deviation in errors (default 3)
  and the bootstrap resamplings for the errors (default 1000).
- `--mean-tol`, `--std-tol` — optional fixed relative tolerances in
  addition.

Output: `epsilon2_distribution.pdf` in `--datadir`, plus
`eccentricity_results.json` (or, with `--save-reference`, an updated
`eccentricity_reference.json` next to the script) recording the ε₂ values and
the seed of each, the seeds, the version of IP-Glasma that ran the events
(from their `usedParameters*.dat`, else `ipglasma --version`), the settings
(lattice size, L, maxTime, JIMWLK) and the exact command used to reproduce
the run, and `eccentricity_check.json` with the settings of the script, the
checks and `passed`, as the summaries of the other validations.

The check: the mean and the width of the ε₂ distribution must agree with the
reference within `--n-sigma` combined errors, from a bootstrap of the run and
of the reference, and the two-sample Kolmogorov–Smirnov test must give
p > 0.01. The error of the mean is about 0.092/√N for N events, so the
default 1000 events resolve a change of the mean ε₂ of about 3% (three
combined errors); the error of the reference adds little unless the run has
a comparable number of events. The seeds of the run must not overlap those
of the reference: with
the same code, the same seeds give the same events, so the check would not be
independent. A reference saved before the seeds were recorded used the seeds
from 0, one per event. The script exits non-zero if a check fails, an event
fails or no reference is found.

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
1/`Ny`; no continuum or infinite-volume approximation enters. The
corrections are of order G(0)/`Ny`, not [G(0) − G(r)]/`Ny`: each layer also
carries the local fluctuations of the field, which cancel in D only for
`Ny` → ∞. D therefore approaches D_G as 1/`Ny` only for `Ny` ≳ G(0). With
an infrared regulator, G(0) ≈ (g²μ/m)²/(4π); without one, it grows with the
lattice as (g²μL)². The script runs IP-Glasma (`runEvolution 0`,
`writeWilsonLines 2`, input template `dipole_correlator/input_dipole`) for
each setting and seed, computes D(r) from the two nuclei of each run (two
independent configurations) along the axes and the diagonal, averaged over
all sites, and checks, at the r_s where D_G = e^{−1/2}, that

- D approaches D_G for large `Ny`: |D − D_G| at the largest `Ny` is below a
  fifth of that at the smallest,
- `Ny` (D − D_G) at the two largest `Ny` agrees within three combined
  errors (plus 5%),
- the extrapolation D − D_G = d_∞ + c/`Ny`, fitted for `Ny` ≥ 2 G(0),
  rounded to an integer (`--fit-from 2`), gives a d_∞ consistent with 0
  within three errors,
- a run with `g 2` gives the same D as with `g 1` (to rounding).

A check that cannot be made, e.g. the extrapolation with fewer than three
`Ny` in the fit range, counts as failed (`INCOMPLETE` in the printout,
`null` in the summary).

The fit starts at 2 G(0), not at G(0): there the 1/`Ny`² term still shifts
d_∞ by about 10⁻³ (from `Ny` ≥ G(0), fitted to the same runs, the default
scan gives d_∞ = 1.1, 1.6 and 1.1·10⁻³ for G(0) = 8, 65 and 260, 1.3 to 2
standard deviations), and from 2 G(0) on by a quarter of that. (`--fit-from`
also sets how far the scan is extended, so `--analyze-only --fit-from 1`
fits fewer `Ny` and gives other values.)

For G(0) ≲ 1, D − D_G changes sign near r_s and is not resolved there,
while at smaller r it is (15 to 100 standard errors for `Ny 1` and r ≤ 10a
in the two settings with G(0) = 0.5) and falls as 1/`Ny` from `Ny 1` on.
The checks are therefore made at r_s where D − D_G is resolved there at the
smallest `Ny` (more than five standard errors), and otherwise at the r ≤
r_max where it is most significant. Where fewer than three `Ny` of `--ny`
lie in the fit range (m = 0), the scan of that setting is extended by
doubling the largest `Ny`, up to `--ny-max` (3200). A setting in which D_G
does not fall to e^{−1/2} within N/4 has no r_s and is skipped. Each `Ny`
has its own seeds (1000 `Ny` + `--first-seed` + i), so the configurations
of different `Ny`, and the points of the fit, are independent. The printout
and the summary give G(0), the r of the checks, the fit range (marked `*`)
and the fit for each setting.

It also reports Q_s/g²μ, with D(r_s) = e^{−1/2} and Q_s² = 2/r_s², the
saturation scale of the fundamental representation in coordinate space.

```bash
cd dipole_correlator
# quick test: 64², g²μ = 1 GeV, m = 0.2 GeV, Ny = 1...50, 2 seeds (seconds)
python3 dipole_correlator.py --ipglasma-path ../.. --quick --max-workers 4
# the scan for the paper: 512², L = 25.6 fm, g²μ = 0.5, 1, 2 GeV,
# m = 0, 0.2, 0.4 GeV, Ny = 1...100 (for m = 0 up to 200, 800 and 3200),
# 4 seeds (about 45 minutes with 7 workers)
python3 dipole_correlator.py --ipglasma-path ../.. --max-workers 7
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
runs what is missing; `--keep-files` keeps the Wilson lines; `--threads`
sets the OpenMP threads per run (default 1). The script exits non-zero if a
check fails.

With the defaults (8 configurations per `Ny`), all checks pass; the scan
takes 43 minutes with 7 workers of one thread, most of it the four runs
with `Ny 3200` (23 minutes each; with `--ny-max 1600`, half the time, the
extrapolation for g²μ = 2 GeV, m = 0 cannot be made and the script fails):

| g²μ [GeV] | m [GeV] | G(0) | `Ny` | r of the checks | fit range | d_∞ [10⁻³] |
|---|---|---|---|---|---|---|
| 0.5 | 0 | 16.3 | 1–200 | r_s = 13a | 50–200 | 3.5 ± 3.6 |
| 0.5 | 0.2 | 0.50 | 1–100 | a | 1–100 | 0.005 ± 0.003 |
| 0.5 | 0.4 | 0.13 | | | | skipped (D_G > 0.84) |
| 1 | 0 | 65.1 | 1–800 | r_s = 6a | 200–800 | 1.7 ± 2.2 |
| 1 | 0.2 | 1.99 | 1–100 | r_s = 10a | 5–100 | −0.2 ± 0.8 |
| 1 | 0.4 | 0.50 | 1–100 | a | 1–100 | 0.021 ± 0.010 |
| 2 | 0 | 260 | 1–3200 | r_s = 3a | 800–3200 | 0.6 ± 1.3 |
| 2 | 0.2 | 7.98 | 1–100 | r_s = 4a | 25–100 | 0.5 ± 0.8 |
| 2 | 0.4 | 2.00 | 1–100 | r_s = 5a | 5–100 | −0.1 ± 0.3 |

D − D_G at r_s falls, e.g., from 0.09 (`Ny 1`) to 0.002 (`Ny 100`) for
g²μ = 1 GeV and m = 0.2 GeV, and from 0.31 to 0.006 (`Ny 3200`) for
g²μ = 2 GeV and m = 0; for G(0) = 0.5 it is at most 1.5 standard
deviations at all `Ny` at r_s, while at r = a it falls from 4·10⁻⁴ and
1.3·10⁻³ (70 and 110 standard errors) at `Ny 1` by a factor of 40 and 70.
Q_s/g²μ approaches its value for D_G accordingly, 0.529 at `Ny 100` against
0.532, and 1.029 at `Ny 3200` against 1.039. The run with `g 2` agrees with
`g 1` to 1·10⁻¹⁶. With 8 configurations, the extrapolation check is weak
where its bound, three errors of d_∞, is as large as the deviation left at
the largest `Ny` of the fit (g²μ = 0.5 GeV with m = 0, and 1 GeV with m = 0
and 0.2 GeV); more configurations tighten it.

#### Comparison with Lappi (2008)

```bash
# N = 200, g²μL = 100, g²μa = 0.5, m = 0 and 0.1 g²μ (with and without the
# zero mode), Ny = 1...150, 32 seeds (about 10 minutes with 4 workers;
# needs scipy)
python3 dipole_correlator.py --ipglasma-path ../.. --max-workers 4 --lappi
```

`--lappi` runs, instead of the scan, the setup of T. Lappi, Eur. Phys. J.
C 55 (2008) 285 ([arXiv:0711.3039](https://arxiv.org/abs/0711.3039)), that of
his Fig. 4: a 200² lattice with g²μL = 100 and g²μa = 0.5 (`g2muGeV 1`; only
these combinations enter), m = 0 and 0.1 g²μ, and `Ny` = 1, 2, 4, 10, 20,
50, 100 and 150, with 64 configurations each. It uses the definitions and
the procedure of his code and scripts (#148):
- **The adjoint correlator** is the factorized form |S(r)|² − 1, with
  S(r) = Σ_x tr V†(x) V(x+r)/N² of each configuration. His code takes
  2|S|² − 1, twice that plus a constant, which has the same maximum of
  k̃² C(k̃) but not the same C(r).
- **Q_s** is the k̃ of the maximum of k̃² C(k̃), with C averaged over the
  modes in bins of width 2√2/N in k̃ (2√2/(2N) for the fundamental
  spectrum) and each bin at the mean k̃ of its modes. The maximum of each
  configuration comes from a parabola around its largest bin, weighted
  with k̃, over a window that grows by 10% on either side until it holds at
  least 4 bins and either 10 bins or an upper limit of 1.5 times the k̃ of
  the largest bin. A fit whose maximum lies more than 50% away is dropped;
  this happens for 4–15 of the 64 fundamental spectra and none of the
  adjoint ones. Q_s and its error are the mean of these maxima and its
  standard error, as his error bars.
- **The zero mode** of the kernel is removed also for m > 0
  (`removeZeroMode 1`). The runs with m > 0 are made with and without it.
- **In coordinate space,** C(r) is the 2D momentum integral of the binned
  spectrum, Σ_bins (modes) C(k̃) J₀(k̃r)/N², with k̃ as the momentum and the
  k = 0 mode by its magnitude, as in the spectrum files of his code, and
  C(r_s) = C(0) e^{−1/2}, Q_s² = 2/r_s². The k = 0 mode of 2|S|² − 1 is
  negative (2c₀ + N² ≈ −36 000 to −39 000 for N² = 40 000), so its
  magnitude adds about 1.9 to C(r) at every r (C(0) ≈ 18.9 instead of 17).
  Fig. 7 is reproduced with this offset, a property of the procedure rather
  than of a definition of Q_s.

The script ports his scripts to Python (`lappi_procedure()`); on the same
spectra, the ports agree with his `runn.pl`, `kfit.pl` and `kfitone.pl` to all
digits they print. It checks that, for m = 0 and for m > 0 with the zero mode
removed, the factorized adjoint and the fundamental √(C_A/C_F) Q̃_s agree with
his Fig. 4 within `--lappi-n-sigma` (3) combined errors at every `Ny`, and the
coordinate-space Q_s with his Fig. 7 within `--lappi-coord-tolerance` (3%). A
check that cannot be made, e.g. at an `Ny` whose single fits are all dropped,
counts as failed. His values were read off the vector graphics of Figs. 4 and 7
(marker centres and, in Fig. 4, error bars; Fig. 7 has none) and are listed in
the script. The errors of the other columns are from the jackknife over the
configurations. The results, written to `<datadir>/dipole_lappi_summary.json`
and `<datadir>/dipole_lappi.pdf`:

| m/g²μ | `Ny` | adjoint | Lappi, Fig. 4 | fundamental √(C_A/C_F) Q̃_s | Lappi, Fig. 4 | Gaussian fit | ⟨\|tr V†V\|²⟩ − 1 | coordinate | Lappi, Fig. 7 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 0.569 ± 0.003 | 0.570 ± 0.002 | 0.571 ± 0.012 | 0.583 ± 0.005 | 0.532 | 0.689 ± 0.007 | 0.633 ± 0.003 | 0.637 |
| 0 | 2 | 0.617 ± 0.003 | 0.620 ± 0.002 | 0.602 ± 0.010 | 0.623 ± 0.006 | 0.588 | 0.652 ± 0.004 | 0.683 ± 0.003 |  |
| 0 | 4 | 0.695 ± 0.002 | 0.695 ± 0.002 | 0.687 ± 0.009 | 0.690 ± 0.010 | 0.677 | 0.735 ± 0.005 | 0.767 ± 0.002 |  |
| 0 | 10 | 0.844 ± 0.001 | 0.842 ± 0.001 | 0.853 ± 0.011 | 0.845 ± 0.006 | 0.829 | 0.899 ± 0.003 | 0.927 ± 0.001 | 0.921 |
| 0 | 20 | 0.959 ± 0.002 | 0.961 ± 0.001 | 0.950 ± 0.009 | 0.972 ± 0.007 | 0.943 | 1.017 ± 0.003 | 1.054 ± 0.001 |  |
| 0 | 50 | 1.075 ± 0.002 | 1.077 ± 0.001 | 1.074 ± 0.014 | 1.075 ± 0.008 | 1.058 | 1.148 ± 0.004 | 1.181 ± 0.001 | 1.166 |
| 0 | 100 | 1.125 ± 0.002 | 1.126 ± 0.002 | 1.102 ± 0.017 | 1.134 ± 0.013 | 1.108 | 1.198 ± 0.004 | 1.238 ± 0.001 |  |
| 0 | 150 | 1.146 ± 0.002 | 1.143 ± 0.002 | 1.147 ± 0.016 | 1.147 ± 0.014 | 1.129 | 1.219 ± 0.005 | 1.261 ± 0.002 |  |
| 0.1 | 1 | 0.521 ± 0.002 | 0.519 ± 0.002 | 0.490 ± 0.008 | 0.502 ± 0.009 | 0.476 | 0.584 ± 0.003 | 0.627 ± 0.002 | 0.623 |
| 0.1 | 2 | 0.618 ± 0.001 | 0.617 ± 0.001 | 0.590 ± 0.008 | 0.592 ± 0.004 | 0.576 | 0.661 ± 0.003 | 0.733 ± 0.001 |  |
| 0.1 | 4 | 0.705 ± 0.001 | 0.703 ± 0.002 | 0.673 ± 0.011 | 0.674 ± 0.008 | 0.675 | 0.760 ± 0.003 | 0.837 ± 0.001 |  |
| 0.1 | 10 | 0.795 ± 0.001 | 0.794 ± 0.001 | 0.776 ± 0.011 | 0.755 ± 0.012 | 0.761 | 0.860 ± 0.011 | 0.946 ± 0.001 | 0.934 |
| 0.1 | 20 | 0.832 ± 0.001 | 0.831 ± 0.002 | 0.825 ± 0.012 | 0.794 ± 0.010 | 0.790 | 0.896 ± 0.003 | 0.990 ± 0.001 |  |
| 0.1 | 50 | 0.858 ± 0.001 | 0.856 ± 0.002 | 0.837 ± 0.012 | 0.815 ± 0.010 | 0.817 | 0.919 ± 0.004 | 1.022 ± 0.001 | 1.003 |
| 0.1 | 100 | 0.864 ± 0.001 |  | 0.828 ± 0.014 |  | 0.825 | 0.926 ± 0.005 | 1.030 ± 0.001 |  |
| 0.1 | 150 | 0.864 ± 0.001 | 0.868 ± 0.002 | 0.834 ± 0.013 | 0.814 ± 0.020 | 0.823 | 0.929 ± 0.004 | 1.030 ± 0.001 |  |

(m = 0.1 g²μ with the zero mode removed.)

- **Momentum space:** the adjoint Q_s agrees with Fig. 4 within 0.5% at
  every `Ny`, within at most 1.2 (m = 0) and 1.6 (m = 0.1 g²μ) combined
  errors, and the fundamental one within at most 1.9 and 2.0 errors. For
  the adjoint spectrum, the maximum of the averaged spectrum, from the
  parabola through the largest bin and two bins on either side in bins of
  width 2π/L (`peak_momentum()`) or from his `kfit.pl`, differs from the
  mean of the single maxima by less than 1%; for the noisier fundamental
  spectrum, by up to 3.5% and 9%.
- **The zero mode:** with it kept, as IP-Glasma does by default, the
  adjoint Q_s at m = 0.1 g²μ lies 1.4–3.4% lower for `Ny` ≤ 10 (0.504, 0.597,
  0.689 at `Ny` = 1, 2, 4), 5–12 errors from Fig. 4; for large `Ny` the two
  agree.
- **Coordinate space:** Fig. 7 is reproduced within 2%. Against |S|² − 1
  on the lattice with C(r_s) = 8 e^{−1/2}, about 4.5% of Q_s comes from the
  correlator 2|S|² − 1 normalized to its own C(0) = 17, and about 8% from
  the magnitude of its k = 0 mode; the momentum integral with k̃, instead of
  the lattice Fourier transform, changes Q_s by less than 1%.
- **The definition of the paper,** ⟨|tr V†V|²⟩ − 1, gives Q_s 5–8% higher
  (21% and 12% for `Ny 1` at m = 0 and 0.1 g²μ); for `Ny` → ∞, where both
  are Gaussian, the two forms differ essentially by C_A = 3 against
  2C_F = 8/3.
- **The Gaussian fits** of his `runn.pl`, of ln C(k̃) = const − k̃²/Q_s²
  over 0.1 Q_s < k̃ < 1.3 Q_s, lie 1.5–6.5% (m = 0) and 4–9% (m = 0.1 g²μ)
  below the maxima, so Fig. 4 does not show them.

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
0.3%. In the setup of the paper (256², L = 30 fm, 15 s), T^μν and the
observables of (0) are bit-identical (the gluon spectrum to 7·10⁻¹¹); (a)
holds everywhere, T^μν to 2·10⁻¹¹; in (b) 1.2% of the T^μν points differ
and the observables by at most 0.13%.

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
  the runs count gluons with the electric fields, n = 2|E|²/ω
  (`gluonSpectrumEstimator 0`, the default of the script, unlike the code's),
  which oscillate with kτ. At late times J₀², J₁² → 1/(πkτ) on average, and the
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
- **The projection onto the free Bjorken modes:** with
  `--gluon-spectrum-estimator 1`, the runs count gluons with
  `gluonSpectrumEstimator 1`, and the spectrum is compared with the
  late-time formula, without Bessel functions. With the defaults, the ratio
  in the window is 1.008 ± 0.006 at τ = 0.5, 0.2 and 0.1 fm alike, and
  dN/dy is the same to six digits at all three times, while 2|E|²/ω gives
  15.0 ± 0.9%, 20.4 ± 0.4% and 0.6 ± 0.6% more (the mean over the 16
  configurations, `dN_dy_electric_over_projection` in the summary). The
  projection also counts the modes on the lines k_x = −π/a and k_y = −π/a,
  which 2|E|²/ω and the leading-order sum leave out. Giving each of their
  modes the mean spectrum of its k̃ bin, they carry at most 3.7 × 10⁻⁵ of
  dN/dy here (`pi_line_fraction_max` in the summary), far below the errors.

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
  ratio oscillates between 0.77 and 1.32 below k ≈ 4 GeV, where kτ is only a
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

### The two definitions of the gluon number

`gluon_spectrum_estimator/gluon_spectrum_estimator.py` (issue #140) compares
the two values of `gluonSpectrumEstimator` in the example input: Pb+Pb with
JIMWLK from x = 0.01 to 10⁻⁴, b = 0, `seed 3`, 768² at L = 22.5 fm, and two
color-charge configurations (`colorChargeSeed 1` and `2`).
- **Runs:** the gluon multiplicity is measured at the final time, so each
  time (default τ = 0.2, 0.4 and 1 fm) and configuration is a run of its
  own. The runs use `gluonSpectrumEstimator 1`, whose
  `gluonMultiplicity<id>.json` also gives dN/dy and dE/dy of `0`, from the
  same gauge-fixed fields.
- **Output:** dN/dy and ⟨k_T⟩ = (dE/dy)/(dN/dy) of both, their ratio, and
  the longitudinal part of the gauge-fixed potential from the log, a check
  of the gauge fixing.
- **Check:** 2|E|²/ω counts the gluons of a mode only on average over the
  oscillation of its fields, i.e. for kτ ≫ 1. From `--agree-from` (0.4 fm)
  on, the mean ratio of dN/dy and of ⟨k_T⟩ must therefore be 1 within
  `--tolerance` (3%).

```bash
cd gluon_spectrum_estimator
# quick test: 128², one configuration at 0.2 and 0.4 fm
python3 gluon_spectrum_estimator.py --ipglasma-path ../.. --quick
# the setup of the paper: 6 runs of about 5 minutes with 10 threads
python3 gluon_spectrum_estimator.py --ipglasma-path ../.. --threads 10
```

With the defaults, 2|E|²/ω over the projection, for the two configurations:

| τ [fm] | dN/dy | ⟨k_T⟩ |
|---|---|---|
| 0.2 | 0.860, 0.928 | 1.164, 1.078 |
| 0.4 | 1.007, 1.004 | 0.999, 1.001 |
| 1.0 | 1.013, 1.003 | 0.988, 0.997 |

- **From τ = 0.4 fm on,** the two agree within about 1%. The spectrum peaks
  near Q_s ≈ 2 GeV, where kτ ≈ 4 at 0.4 fm, and integrated over k the
  oscillation averages out. In the dilute limit, in contrast, the lowest
  modes dominate, and 2|E|²/ω is off by tens of percent at all times (see
  the dilute limit above).
- **At τ = 0.2 fm,** they differ by 7–14% in dN/dy and by 0.2–0.4 GeV in
  ⟨k_T⟩. The fields are still strongly nonlinear, and no estimator gives a
  well-defined particle number yet.
- **The longitudinal part** of the gauge-fixed potential is 2–6%.
- **The lines k = −π/a,** which only the projection counts, carry about
  10⁻⁴ of dN/dy (0.9–1.4 × 10⁻⁴, `pi_line_fraction` in the summary,
  estimated from the mean spectrum of their k̃ bins).

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
| dN/dy | 3687 | 4744 | 5154 | 5381 | 5312 |
| dE/dη (0.4 fm) [GeV] | 12147 | 22022 | 28936 | 33240 | 34810 |

dN/dy is within about 4% of its limit for g²μ a ≲ 0.6, while dE/dη still
grows by 4.7% ± 1.4% from g²μ a = 0.3 to 0.15. The energy sits in modes up
to several g²μ, where the gluon spectrum falls like 1/k⁴, and converges
only like a: its deficit is about 0.2–0.3 g²μ a, so 5% needs g²μ a ≈
0.15–0.25 (the energy density at τ = 0 even diverges like ln² a, T. Lappi,
Phys. Lett. B 643 (2006) 11). This is not an artifact of IP-Glasma: an
independent classical Yang-Mills code (a version of S. Schlichting's), run
on the same Wilson lines, gives dE/dη within 0.1–0.6% of IP-Glasma's and
the same dependence on a (see the discussion of #134). With JIMWLK, the
gluon spectrum falls less steeply (the anomalous dimension of the evolved
Wilson lines), so the energy converges even more slowly.

**Convergence check.** For each observable of the size scans, the script
compares the two finest lattices: the relative change must be below
`--convergence-tolerance` (1%) plus two standard deviations. For an
observable that fails, it prints a prominent `WARNING: NOT CONVERGED IN THE
LATTICE SPACING` with the change and an estimate of the remaining deviation
from an extrapolation linear in a, marks the panel in the figure, and lists
it under `not_converged` in the summary. With the defaults (20 minutes with
2 workers of 12 threads), the results are

| size | 96 | 192 | 384 | 768 | 1024 | 1280 |
|---|---|---|---|---|---|---|
| a [fm] | 0.234 | 0.117 | 0.059 | 0.029 | 0.022 | 0.018 |
| dN/dy | 8151 | 14170 | 18836 | 21289 | 21713 | 22046 |
| dE/dη (0.4 fm) [GeV] | 9062 | 28602 | 51825 | 68667 | 72709 | 75766 |

dN/dy changes by +1.5% ± 2.1% in the last step, ε₂ and ε₃ by less than
their errors, so they converge. dE/dη still grows by 4.2% ± 1.2%, and the
linear extrapolation puts its value at a = 0 about 16% above the finest
lattice (5% for the constant g²μ): converging the energy to the percent
level would need lattice spacings far below 0.018 fm, so its dependence on
a has to be quoted as a systematic uncertainty. The `Ny` scan (at 768²)
converges for dN/dy by `Ny` ≈ 25–50 (21289 for `Ny 50`, 20912 for
`Ny 100`); the `m` scan changes dN/dy from 28196 (`m 0.1`) to 15782
(`m 0.8`) at fixed `QsMuRatio`, which changes the physical Q_s. dN/dy is
the projection onto the free Bjorken modes (`gluonSpectrumEstimator 1`,
the default); 2|E|²/ω (`0`) gives, on the same events at 0.4 fm, between
0.5% less and 1.4% more in the size scan and 1–3.4% less for the constant
g²μ.

```bash
cd numerical_parameters
# quick test: sizes 64 and 128, Ny and m scans at 64, the constant g²μ at
# sizes 32 and 64, 2 seeds (6 s with 8 workers)
python3 numerical_parameters.py --ipglasma-path ../.. --quick --max-workers 8
# the setup of the paper (the defaults): sizes 96-1280, Ny 1-100,
# m 0.1-0.8 GeV, the constant g²μ at sizes 32-512, 4 seeds (20 minutes)
python3 numerical_parameters.py --ipglasma-path ../.. --threads 12 \
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
estimates agree within about 1% on the two finest lattices, which puts the
deficit at about 0.2 g²μ a with the projection (0.3 with 2|E|²/ω in the
discussion of #134); the extrapolation of dE/dη linear in a gives 0.3.

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

With the commands below (two seeds, N = 768 to 1536; 26 and 7 minutes with
`--threads 12 --max-workers 2`):
- **With JIMWLK:** dN/dy converges (+2.1% ± 3.1% in the last step). dE/dη
  grows by 22%, 13% and 12% in the three steps (+11.8% ± 3.1% in the
  last), because the evolution hardens the tail to p ≈ 1.4.
- **At x = 0.01 without JIMWLK:** p ≈ 0.3, and dE/dη is converged within
  about 1% from N = 1024 on (+1.0% ± 2.7% in the last step). In this dilute
  event, dN/dy differs by up to 18% between the two seeds of one lattice
  (5222 to 6512 over all eight runs), so two seeds do not resolve its
  dependence on a: the script flags its last step, +12% ± 2%, as not
  converged.
- **ε₂ and ε₃** need more than two seeds in both variants.

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
  transverse area, at the last time. The runs compute the multiplicity in
  transverse Coulomb gauge with the free lattice dispersion and the electric
  fields (`gluonSpectrumEstimator 0`), as Eqs. (27) and (28) of Lappi, the
  definition of his Table I.

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
# the setup of the paper: g²μL = 72, 120 and 192 on 128², 256² and 512²,
# 16 seeds (25 minutes with 6 workers of 4 threads)
python3 constant_glasma.py --ipglasma-path ../.. --g2muL 72 120 192 \
    --seeds 16 --threads 4 --max-workers 6
```

In the setup of the paper, all checks pass. f_E and f_N agree with Lappi's
Table I on 512² (within 1.2 combined standard deviations at g²μL = 72 and
120; at 192, f_E and f_N lie 1.9 and 1.8 standard deviations below), and
the plateau of f_E with the points of his Fig. 3 (read off the figure) on
the coarser lattices:

| g²μL | N   | g²μa | f_E           | Lappi         | f_N           | Lappi         |
|------|-----|------|---------------|---------------|---------------|---------------|
| 72   | 128 | 0.56 | 0.235 ± 0.004 | ≈ 0.237       | 0.302 ± 0.005 |               |
| 72   | 256 | 0.28 | 0.254 ± 0.004 | ≈ 0.265       | 0.303 ± 0.005 |               |
| 72   | 512 | 0.14 | 0.264 ± 0.004 | 0.265 ± 0.005 | 0.303 ± 0.005 | 0.297 ± 0.006 |
| 120  | 128 | 0.94 | 0.192 ± 0.003 | ≈ 0.193       | 0.288 ± 0.004 |               |
| 120  | 256 | 0.47 | 0.220 ± 0.003 | ≈ 0.217       | 0.298 ± 0.004 |               |
| 120  | 512 | 0.23 | 0.229 ± 0.003 | 0.227 ± 0.003 | 0.295 ± 0.004 | 0.289 ± 0.003 |
| 192  | 128 | 1.50 | 0.173 ± 0.002 |               | 0.290 ± 0.003 |               |
| 192  | 256 | 0.75 | 0.210 ± 0.003 |               | 0.316 ± 0.004 |               |
| 192  | 512 | 0.38 | 0.227 ± 0.003 | 0.238 ± 0.005 | 0.316 ± 0.004 | 0.329 ± 0.006 |

g²ε(0⁺)/(g²μ)⁴ grows from 0.59 to 0.82 and 1.07 (g²μL = 72), from 0.39
to 0.59 and 0.77 (120) and from 0.27 to 0.49 and 0.67 (192) as a is
halved, while f_N changes by less than 4% for g²μa < 1 and by less than
1% for g²μa ≤ 0.75. P_L/ε is between −0.994 and −0.985 at the first time
step and +0.016 to +0.027 at g²μτ = 20 (Lappi's Fig. 5 gives about
+0.025), and f_E decreases by 2.4–2.9% between g²μτ = 10 and 20. The
errors are those of the mean over the 16 seeds, comparable to Lappi's.
The defaults of the script (g²μL = 72 and 120, four seeds) take 9 minutes
with 2 workers of 4 threads.

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
# the setup of the paper (33 minutes)
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
# the setup of the paper (3.6 hours with 4 workers of 6 threads)
python3 lappi_mantysaari.py --ipglasma-path ../.. --threads 6 --max-workers 4
```

With the defaults, the step is again of first order, but its error per unit
of step and of rapidity is about 30 times larger than above:

| Q_s/Λ_QCD at y = 8 | |
|---|---|
| `jimwlkDs 0.001` | 33.1 |
| `jimwlkDs 0.0005` | 30.2 |
| `jimwlkDs 0.00025` | 28.7 |
| zero step | 27.3 |
| BK | 30.2 |
| Fig. 6 of the paper | 28.5 |

- **The paper matches a step of about 0.00025:** Q_s within 0.6% of its
  Fig. 6 at y = 4 and 8, and λ along its Fig. 8.
- **BK** reproduces λ of Fig. 11 within 0.6%.
- **At zero step,** λ of JIMWLK is 3–8% below BK for Q_s/Λ_QCD = 10–25, so
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
# 768², L = 22.5 fm, 3 seeds (52 minutes with 2 workers of 8 threads)
python3 example_input_step_size.py --ipglasma-path ../.. --threads 8 \
    --max-workers 2
# 256², L = 30 fm, the lattice of the example input, 16 seeds (12 minutes
# with 4 workers of 6 threads)
python3 example_input_step_size.py --ipglasma-path ../.. --size 256 --L 30 \
    --seeds 16 --threads 6 --max-workers 4
```

The deviation from zero step at `jimwlkDs 0.005` is (dN/dy with the
projection, the default):

| | 768², L = 22.5 fm, 3 seeds | 256², L = 30 fm, 16 seeds |
|---|---|---|
| Q_s, all | +0.9% ± 0.6% | −0.4% ± 0.2% |
| Q_s, dilute third | +1.2% ± 0.7% | −0.4% ± 0.4% |
| Q_s, dense third | +0.5% ± 0.8% | −0.3% ± 0.3% |
| dN/dy | +2.6% ± 1.2% | +0.9% ± 1.8% |
| dE/dη | +5.6% ± 0.6% | +0.2% ± 0.4% |

- **Q_s and dN/dy:** within about 1% (Q_s) and 3% (dN/dy) of zero step at
  0.005, as in the constant-g²μ setup above.
- **Seeds at 256²:** with four seeds, the error of dN/dy at 0.005 is 6%,
  so its 3% check can fail by chance (+4.1% ± 6.4% in one such run); 16
  seeds give 1.8%.
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
- `--max-workers N` — number of seeds to run in parallel (default: the
  logical CPUs available divided by `OMP_NUM_THREADS`); each event runs with
  `OMP_NUM_THREADS` threads (default 1), and without `--slurm` several seeds
  at a time run with `OMP_PROC_BIND=false`, since binding would put them on
  the same cores. Each seed needs about 1 GB of memory while IP-Glasma runs.
- `--plot-only` — skip the simulation and only compute the cross sections,
  the comparison plot and the check (or, with `--save-reference`, the
  reference) from the spectra already in `--datadir`; the seeds, the number
  of events, the version of IP-Glasma and the settings are taken from its
  `run_info.json` (the seeds from `--first-seed` if it is given; without
  `run_info.json`, the version from `ipglasma --version` of
  `--ipglasma-path`).
- `--keep-logs` — keep the per-seed working directories with their logs
  instead of deleting them (those of failed seeds are always kept).
- `--first-seed S` — the events use the seeds S, S + 1, ..., two events
  (projectile and target) per seed: by default 10000, so that a check is
  independent of the reference, and 0 with `--save-reference`.
- `--check-tmax`, `--n-sigma`, `--bootstrap` — the range of |t| of the
  check (default 2 GeV²), the allowed deviation in errors (default 3) and
  the bootstrap resamplings for the errors (default 1000).
- `--outlier-zmax Z` — a spectrum whose largest |A|² is more than Z robust
  standard deviations of log₁₀|A|² from the median is an outlier
  (default 7).
- `--save-reference` — save the results of the run as the new reference
  (`--reference-file`, default `jpsi_reference.json` next to the script)
  instead of comparing with it; it is not saved if a spectrum is invalid or
  an outlier, or if its seeds are not known (`--plot-only` without
  `run_info.json` or `--first-seed`). A new reference (5000 events, seeds
  0–2499):
  ```bash
  python3 parallel_test_vector_meson_production.py --maxevents 5000 \
      --save-reference --ipglasma-path ../.. \
      --subnucleondiffraction-path /path/to/subnucleondiffraction
  ```

Output: `jpsi_results.json` and `cross_section.pdf` in `--datadir`
(coherent/incoherent cross section vs. reference, with a ratio panel), plus
`output_1.pdf`/`output_5.pdf` color-field (Tr V(x)) plots for a couple of
representative events. The spectra of the events are copied to `--datadir`
as `spectra_1` to `spectra_<maxevents>`, whatever the seeds (the cross
sections and the check use only these, also if `--datadir` holds more from
an earlier run). The script computes the cross sections from the amplitudes
of the spectra as `python/cross_section.py` of subnucleondiffraction does:
the coherent one is the square of the mean of Re A, the incoherent one the
variance of Re A and Im A. `jpsi_results.json` has dσ/dt (in mb/GeV²) and its
integral over |t| ≤ `--check-tmax` (in mb), with the errors of a bootstrap
over the events, the number of events and the seeds, the versions of
IP-Glasma (`ipglasma --version`) and of subnucleondiffraction (from the
`# Git version` line of the spectra), the settings (input template and its
SHA-256, x_P, lattice, options of subnucleondiffraction) and the command;
the reference `jpsi_reference.json` has the same format.

The check (`jpsi_check.json`, with the settings of the script, the checks
and `passed`): the coherent and incoherent cross sections, integrated over
|t| ≤ `--check-tmax`, must agree with those of the reference within
`--n-sigma` errors. The error is the one expected if the two agree: the
relative bootstrap error of the reference, scaled to the events of the
reference and of the run, √(1 + N_ref/N_run). The run's own bootstrap error
(also printed) is not used: it shrinks with the cross section, so that runs
that fluctuate low would fail more often than runs that fluctuate high.
Against a reference of 2000 other events of the same code, the earlier
check (the run's own error, and the binned errors of the reference added in
quadrature) failed 1.4% (coherent) and 5% (incoherent) of the runs with 32
events beyond 3 errors, and 1.4% and 1% of those with 200 events, all on
the low side; with this error, 0.2–0.5% fail, against 0.27% expected. The
coherent cross section is the square of the mean
amplitude, which fluctuates by about 40% between events at t = 0, and the
incoherent one is its variance, so 200 events give errors of about 8% and
9%, and 1000 events about 4%. The spectra must be finite, on the same t
values, and none an outlier: its scale, the largest |A|² over t, must lie
within `--outlier-zmax` robust standard deviations (MAD/0.6745) of
log₁₀(scale) from the median, as in the outlier filter of
`utilities/compute_tdiffCrossSections.py` of IPGlasmaFramework. The check
fails instead of dropping such an event, since a blown-up event is what it
should find. The seeds of the run must not overlap those of the reference.
Settings in which the run and the reference differ are printed as warnings.
A seed for which IP-Glasma or subnucleondiffraction fails gives no spectra
and is listed in `failed_seeds.txt`; its working directory stays, with the
logs but without the Wilson-line files. The spectra `spectra_1` to
`spectra_<maxevents>` and the results (cross sections, plots, check and
`failed_seeds.txt`) of an earlier run in `--datadir` are removed when a run
starts. The script exits non-zero if a check fails or a seed fails.

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
  the job runs; both scripts remove those of the successful events at the
  end unless `--keep-logs` is passed (those of failed events stay), so pass
  `--keep-logs` if you want to inspect individual events afterwards (or
  need `--plot-only`/`--clean` to resume from an interrupted job).

# Running the validation scripts

This folder contains two parallel validation/regression scripts for IP-Glasma,
and in subfolders the small validation setups of the code paper, which test
parts of the code against exact or published results:

- **`parallel_test_vector_meson_production.py`** — runs many independent
  IP-Glasma + JIMWLK events, feeds the resulting Wilson lines to the external
  [`subnucleondiffraction`](https://github.com/hejajama/subnucleondiffraction)
  code, and compares the resulting coherent/incoherent J/ψ cross section to a
  stored reference (`validation_spectra`) and to Fig. 1 ("CGC+shape fluct") of
  [arXiv:2207.03712](https://arxiv.org/pdf/2207.03712).
- **`parallel_test_eccentricity.py`** — runs many independent IP-Glasma
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

Both scripts parallelize over random seeds using
`concurrent.futures.ProcessPoolExecutor`, so each event runs in its own
working directory (`<datadir>/seed_<n>/`) and multiple events run
simultaneously.

## Prerequisites

1. Build the `ipglasma` binary from the repository root:
   ```bash
   ./compile_IPGlasma.sh
   ```
   (dependencies: CMake, FFTW and GSL, see the main README). The scripts
   locate the binary via `--ipglasma-path` (default `.` for the eccentricity
   script and `..` for the vector-meson script), expecting it at
   `<ipglasma-path>/ipglasma`.
2. For `parallel_test_vector_meson_production.py` only: build
   [`subnucleondiffraction`](https://github.com/hejajama/subnucleondiffraction)
   separately and note its path (needed for `--subnucleondiffraction-path`).
3. Python 3 with `numpy` and `matplotlib` (and optionally `scipy`, used for
   the two-sample KS test in the eccentricity script if available).

## Running interactively

From the `validations/` directory (paths below assume the scripts are run
from there; both scripts also accept absolute paths for `--ipglasma-path`
etc. if you invoke them from elsewhere):

### Eccentricity regression test (quick sanity check)

```bash
# One-time: record a baseline from a known-good build
python3 parallel_test_eccentricity.py --input-template input_eccentricity \
    --ipglasma-path .. --maxevents 40 --max-workers 8 --save-reference

# Later, e.g. after changing the code: check nothing drifted
python3 parallel_test_eccentricity.py --input-template input_eccentricity \
    --ipglasma-path .. --maxevents 40 --max-workers 8
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

### Vector meson production validation

```bash
python3 parallel_test_vector_meson_production.py \
    --maxevents 200 \
    --datadir ./jpsi_test_run \
    --ipglasma-path .. \
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

`runvalidation_jpsi.sh` is a SLURM batch script template for
`parallel_test_vector_meson_production.py`:

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
`--ntasks`/`--maxevents`/`--time` to the scale of the run. Submit with:

```bash
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
- `runvalidation_eccentricity.sh` is the same template for
  `parallel_test_eccentricity.py` (submit it with
  `sbatch runvalidation_eccentricity.sh`); only the job name and the
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

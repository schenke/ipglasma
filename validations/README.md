# Running the validation scripts

This folder contains two parallel validation/regression scripts for IP-Glasma:

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

Both scripts parallelize over random seeds using
`concurrent.futures.ProcessPoolExecutor`, so each event runs in its own
working directory (`<datadir>/seed_<n>/`) and multiple events run
simultaneously.

## Prerequisites

1. Build the `ipglasma` binary from the repository root:
   ```bash
   ./compile_IPGlasma.sh
   ```
   (dependencies: CMake, FFTW). The scripts locate the binary via
   `--ipglasma-path` (default varies per script — see below), expecting it at
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

Output: `epsilon2_distribution.pdf`/`.png` in `--datadir`, plus
`eccentricity_results.json` (or, with `--save-reference`, an updated
`eccentricity_reference.json` next to the script) recording the ε₂ values and
the exact command used to reproduce the run. The script exits non-zero if the
new distribution's mean/std or KS test drifts beyond tolerance
(`--mean-tol`, `--std-tol`) from the reference.

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
  (default `input_vm_proton`); the lattice `NucleusQsTableFileName` is read
  from this file rather than hardcoded.
- `--plot-only` — skip the simulation and only regenerate the comparison
  plot from data files already in `--datadir`.
- `--keep-logs` — keep per-seed log files instead of deleting them.

Output: `cross_section.dat`/`.pdf` in `--datadir` (coherent/incoherent cross
section vs. reference, with a ratio panel), plus `output_1.pdf`/`output_5.pdf`
color-field (Tr V(x)) plots for a couple of representative events.

## Running on a computing cluster (SLURM)

`runvalidation.sh` is a SLURM batch script template for
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

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
module add python-data gsl fftw

python3 parallel_test_vector_meson_production.py \
    --max-workers "${SLURM_NTASKS:-1}" \
    --datadir /path/to/datadir \
    --maxevents 1500 \
    --subnucleondiffraction-path /path/to/subnucleondiffraction \
    --ipglasma-path /path/to/ipglasma_upstream
```

Before submitting, edit the placeholders (`--account`, `--datadir`,
`--subnucleondiffraction-path`, `--ipglasma-path`) and adjust
`--ntasks`/`--maxevents`/`--time` to the scale of the run. Submit with:

```bash
sbatch runvalidation.sh
```

Notes:
- Each `--ntasks` maps to one `--max-workers` process; leave
  `--cpus-per-task=1` and `OMP_NUM_THREADS=1` since each event runs
  single-threaded and parallelism comes from running many events at once.
- Both scripts accept a `--slurm` flag that wraps each
  `ipglasma`/`subnucleondiffraction` invocation with
  `srun --ntasks=1 --cpus-per-task=1`, which can help binding/placement on
  some clusters. `runvalidation.sh` above does not pass it (each worker just
  calls the binaries directly) — add `--slurm` to the `python3` invocation if
  your cluster needs `srun` per task.
- The same pattern (a batch script that calls the script with
  `--max-workers "${SLURM_NTASKS:-1}"`) works for
  `parallel_test_eccentricity.py`; just swap the python invocation, e.g.:
  ```bash
  python3 parallel_test_eccentricity.py \
      --input-template input_eccentricity \
      --ipglasma-path /path/to/ipglasma_upstream \
      --datadir /path/to/datadir \
      --maxevents 200 \
      --max-workers "${SLURM_NTASKS:-1}"
  ```
- Expect one `seed_<n>/` working directory per event under `--datadir` while
  the job runs; both scripts clean these up at the end unless `--keep-logs`
  is passed, so pass `--keep-logs` if you want to inspect individual events
  afterwards (or need `--plot-only`/`--clean` to resume from an interrupted
  job).

#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Physics regression test: run many independent IP-Glasma events in
parallel, extract the final-state (end of classical Yang-Mills evolution)
energy density from each event's Tmunu output, compute the
energy-density-weighted eccentricity epsilon_2 for every event, and plot
the resulting epsilon_2 distribution.

This does not depend on any external program (only the "ipglasma" binary
built from this repository) -- unlike parallel_test_vector_meson_production.py,
which additionally drives "subnucleondiffraction". The energy density is
read straight from the binary Tmunu-t*-*.ipgt snapshot IP-Glasma writes at
the end of the run, using the read_tmunu.py and eccentricity.py helpers
from the utilities/ directory of the IP-Glasma checkout given by
--ipglasma-path (i.e. {ipglasma-path}/utilities/*.py) -- this script
itself need not live inside that checkout.

The intended use is as a lightweight "physics output did not silently
change" check: run this once against a known-good version of the code
with --save-reference to record the epsilon_2 distribution's mean/std,
then rerun it (e.g. in CI, or after refactoring) without that flag -- the
script compares the new distribution against the saved reference and
exits non-zero if it drifts beyond its statistical errors: the mean and
the width of eps2 must agree within --n-sigma combined errors (from a
bootstrap of this run and of the reference) and the two-sample KS test
must give p > 0.01. --mean-tol and --std-tol add fixed relative tolerances.

Seeds: a run uses the seeds --first-seed, --first-seed + 1, ...; the
default, 10000, keeps a check independent of the reference, which is made
with the seeds from 0 (the default with --save-reference) and records them.
The check fails if the seeds of the run overlap those of the reference.

Usage (from validations/eccentricity_regression/):
    # one-time: record a baseline from a known-good build
    python3 parallel_test_eccentricity.py --maxevents 5000 --save-reference

    # later: check that nothing changed (1000 events, the default)
    python3 parallel_test_eccentricity.py

By default, the events run on all logical CPUs available, OMP_NUM_THREADS
(default 1) threads each, and the working directories of the successful
events are removed at the end of the run (their eps2 values are kept in
eccentricity_results.json, or in the reference with --save-reference).
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# The utilities/ helpers (read_tmunu.py, eccentricity.py) live inside the
# IP-Glasma checkout given by --ipglasma-path, which is only known once
# argparse has run -- so they can't be imported at module level. Each
# worker process (ProcessPoolExecutor re-executes this module on "spawn"
# platforms) instead imports them lazily via _load_utilities() below,
# using the ipglasma_path passed explicitly into run_seed().

def _load_utilities(ipglasma_path):
    """Import the utilities/read_tmunu.py and utilities/eccentricity.py
    helpers from inside the given IP-Glasma checkout and return
    (find_final_tmunu_file, get_energy_density, eccentricity)."""
    utilities_dir = os.path.join(ipglasma_path, "utilities")
    if utilities_dir not in sys.path:
        sys.path.insert(0, utilities_dir)
    from eccentricity import eccentricity
    from read_tmunu import find_final_tmunu_file, get_energy_density
    return find_final_tmunu_file, get_energy_density, eccentricity


# --------------------------------------------------------------------------
# input-file handling
# --------------------------------------------------------------------------

def read_input_value(template_path, key):
    """Return the (string) value of a "key value" line in an IP-Glasma
    input file, or None if the key is not present."""
    with open(template_path, "r") as f:
        for line in f:
            tokens = line.split()
            if len(tokens) >= 2 and tokens[0] == key:
                return tokens[1]
    return None


def generate_temp_input(template_path, overrides, tag, working_dir):
    """Copy an IP-Glasma input file, replacing the value of any "key
    value" line whose key is in `overrides`. Returns the path to the
    generated temp file (caller is responsible for removing it)."""
    with open(template_path, "r") as f:
        lines = f.readlines()

    remaining = dict(overrides)
    out_lines = []
    end_of_file = len(lines)
    for line in lines:
        tokens = line.split()
        if tokens[:1] == ["EndOfFile"]:
            end_of_file = len(out_lines)
        if len(tokens) >= 2 and tokens[0] in remaining:
            out_lines.append("{0} {1}\n".format(tokens[0], remaining.pop(tokens[0])))
        else:
            out_lines.append(line)
    # Override keys not already present in the template are inserted before
    # EndOfFile (IP-Glasma ignores everything after it).
    out_lines[end_of_file:end_of_file] = [
        "{0} {1}\n".format(key, value) for key, value in remaining.items()]

    fd, temp_path = tempfile.mkstemp(
        prefix="input_ecc_{0}_".format(tag), suffix=".in", dir=working_dir)
    with os.fdopen(fd, "w") as f:
        f.writelines(out_lines)
    return temp_path


# --------------------------------------------------------------------------
# running one event
# --------------------------------------------------------------------------

def threads_per_event():
    """The OpenMP threads of one event: OMP_NUM_THREADS, default 1."""
    try:
        return max(1, int(os.environ.get("OMP_NUM_THREADS", "1")))
    except ValueError:
        return 1


def available_cpus():
    """The logical CPUs this process may run on (e.g. those of a Slurm
    allocation or of taskset), else all of them."""
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:  # not available on macOS
        return cpu_count()


def run_command(cmd, cwd, log_path, slurm=False):
    if slurm:
        cmd = ["srun", "--ntasks=1",
               "--cpus-per-task={0}".format(threads_per_event())] + list(cmd)
    # one thread per event unless OMP_NUM_THREADS says otherwise: the events
    # already run in parallel
    env = dict(os.environ, OMP_NUM_THREADS=str(threads_per_event()))
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                            env=env)
    with open(log_path, "w") as f:
        f.write("$ {0}\n".format(" ".join(cmd)))
        f.write(result.stdout or "")
        f.write(result.stderr or "")
    return result


def extract_eps2(worker_dir, ipglasma_path, event_id=0):
    """Load the final-time Tmunu snapshot from a completed event's working
    directory and compute epsilon_2 from it. Returns (eps2, psi2, tau_fm,
    energy_density). Shared by run_seed() (right after an ipglasma run
    finishes) and collect_from_existing_datadir() (--plot-only, reading
    back events from a previous, possibly interrupted, run)."""
    find_final_tmunu_file, get_energy_density, eccentricity = \
        _load_utilities(ipglasma_path)
    _tau, tmunu_path = find_final_tmunu_file(worker_dir, event_id)
    energy_density, dx, dy, tau_fm = get_energy_density(tmunu_path)
    eps2, psi2 = eccentricity(energy_density, dx, dy, n=2)
    return eps2, psi2, tau_fm, energy_density


def collect_from_existing_datadir(datadir, ipglasma_path):
    """--plot-only support: rebuild the per-event result list by reading
    back already-completed seed_*/ working directories in `datadir`,
    without launching any new ipglasma runs. Useful to regenerate the
    plot (and reference comparison) after a run was interrupted or
    crashed partway through, as long as the individual working
    directories were kept (--keep-logs, or the run never reached its
    end-of-run cleanup)."""
    results = []
    for worker_dir in sorted(glob.glob(os.path.join(datadir, "seed_*"))):
        if not os.path.isdir(worker_dir):
            continue
        name = os.path.basename(worker_dir)
        try:
            seed = int(name[len("seed_"):])
        except ValueError:
            seed = name
        try:
            eps2, psi2, tau_fm, energy_density = extract_eps2(
                worker_dir, ipglasma_path)
        except Exception as exc:  # noqa: BLE001 -- report, keep going
            results.append({
                "seed": seed, "success": False, "worker_dir": worker_dir,
                "error": "failed to extract eps2: {0}".format(exc),
            })
            continue
        results.append({
            "seed": seed, "success": True, "worker_dir": worker_dir,
            "eps2": eps2, "psi2": psi2, "tau_fm": tau_fm,
            "mean_energy_density": float(np.mean(energy_density)),
        })
    return results


def run_seed(seed, ipglasma_path, template_path, qs_table_path, datadir,
             ipglasma_binary, size, L, maxtime, use_jimwlk, slurm):
    """Run a single IP-Glasma event with a distinct random seed in its own
    working directory, then extract epsilon_2 from its final energy
    density snapshot. Runs in a worker process, so it must be
    module-level / picklable and take plain (picklable) arguments."""

    print("=== Running event with seed {0} ...".format(seed), flush=True)

    worker_dir = os.path.join(datadir, "seed_{0}".format(seed))
    os.makedirs(worker_dir, exist_ok=True)

    overrides = {
        "seed": seed,
        "size": size,
        # sizeOutput is the (independent) grid the energy density / Tmunu
        # snapshot is written on; keep it equal to the evolution grid here.
        "sizeOutput": size,
        "nucleusQsTableFileName": qs_table_path,
        "nuclearConfigurationsPath": ipglasma_path + "/nucleusConfigurations/",
        "writeOutputsToHDF5": 0,
        "writeWilsonLines": 0,
        # only the final-time Tmunu snapshot: no intermediate output times
        # (the template has tau = 0.1, 0.2, 0.3, 0.4 fm/c), no hydro, Jazma
        # or initial Wilson-line snapshot
        "writeHydro": 0,
        "writeJazma": 0,
        "writeTmunu": 1,
        "writeTmunuBinary": 1,
        "outputTimes": "none",
        "writeWilsonLineSnapshot": 0,
        "useJIMWLK": 1 if use_jimwlk else 0,
        "useFluctuatingX": 0 if use_jimwlk else 1,
        "useSeedList": 0,
        "jimwlkSaveSnapshots": 0,
    }
    if L is not None:
        overrides["L"] = L
        overrides["LOutput"] = L
    if maxtime is not None:
        overrides["maxTime"] = maxtime

    temp_input = generate_temp_input(template_path, overrides, seed, worker_dir)
    log_path = os.path.join(worker_dir, "run.log")

    t0 = time.time()
    try:
        result = run_command(
            [ipglasma_binary, temp_input], cwd=worker_dir, log_path=log_path,
            slurm=slurm)
    finally:
        if os.path.exists(temp_input):
            os.remove(temp_input)
    elapsed = time.time() - t0

    #if result.returncode != 0:
    #    print("seed {0}: ipglasma failed with code {1} (see {2})".format(
    #        seed, result.returncode, log_path), flush=True)
    #    return {
    #        "seed": seed, "success": False, "worker_dir": worker_dir,
    #        "error": "ipglasma exited with code {0} (see {1})".format(
    #            result.returncode, log_path),
    #    }

    # IP-Glasma is run here as a single MPI-disabled, single-event
    # process, so it always writes rank/event index 0.
    try:
        eps2, psi2, tau_fm, energy_density = extract_eps2(
            worker_dir, ipglasma_path)
    except Exception as exc:  # noqa: BLE001 -- report, don't crash the pool
        return {
            "seed": seed, "success": False, "worker_dir": worker_dir,
            "error": "failed to extract eps2: {0} (see {1})".format(
                exc, log_path),
        }

    return {
        "seed": seed, "success": True, "worker_dir": worker_dir,
        "eps2": eps2, "psi2": psi2, "tau_fm": tau_fm,
        "mean_energy_density": float(np.mean(energy_density)),
        "elapsed_sec": elapsed,
    }


# --------------------------------------------------------------------------
# summary statistics, plotting, regression comparison
# --------------------------------------------------------------------------

def summarize(values):
    values = np.asarray(values, dtype=float)
    return {
        "n_events": int(values.size),
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "median": float(np.median(values)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
    }


def plot_distribution(eps2_values, datadir, reference=None):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    bins = np.linspace(0, max(0.6, max(eps2_values) * 1.1), 25)

    ax.hist(eps2_values, bins=bins, alpha=0.65, color="tab:blue",
            label="this run (N={0})".format(len(eps2_values)),
            density=True)
    if reference is not None and reference.get("eps2_values"):
        ax.hist(reference["eps2_values"], bins=bins, histtype="step",
                linewidth=2, color="tab:red",
                label="reference (N={0})".format(
                    len(reference["eps2_values"])),
                density=True)

    ax.set_xlabel(r"$\epsilon_2$")
    ax.set_ylabel("probability density")
    ax.set_title("Final-state $\\epsilon_2$ distribution")
    ax.legend()
    fig.tight_layout()

    out_path = os.path.join(datadir, "epsilon2_distribution.pdf")
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def load_reference(reference_file):
    if not os.path.exists(reference_file):
        return None
    with open(reference_file, "r") as f:
        return json.load(f)


def reproduce_command():
    """Reconstruct the command line this script was invoked with, so it
    can be recorded alongside a saved reference (JSON has no comment
    syntax, so this is stashed as a regular "_command" field instead)."""
    return "python3 {0}".format(" ".join([sys.argv[0]] + sys.argv[1:]))


def save_reference(reference_file, eps2_values, run_meta, seeds=None):
    # JSON doesn't support comments, so the reproduce command is stored as
    # a leading "_command" field instead -- keep it first so it's the
    # first thing a human opening the file notices.
    payload = {"_command": reproduce_command()}  # type: dict
    payload.update(summarize(eps2_values))
    payload["eps2_values"] = [float(v) for v in eps2_values]
    payload.update(run_meta)
    if seeds is not None:
        payload["first_seed"], payload["last_seed"] = seeds
    with open(reference_file, "w") as f:
        json.dump(payload, f, indent=2)


def check_reference_parameters(run_meta, reference):
    """Return a list of warning messages, one for every run parameter
    recorded in `reference` (and in this run's `run_meta`) whose value
    differs between the two -- comparing eps2 distributions obtained with
    different settings is not a meaningful regression check."""
    # Older reference files stored maxTime under the key "maxtime".
    reference = dict(reference)
    if "maxTime" not in reference and "maxtime" in reference:
        reference["maxTime"] = reference.pop("maxtime")

    warnings = []
    for key, value in run_meta.items():
        if key not in reference:
            warnings.append(
                "WARNING: parameter '{0}' is not recorded in the reference "
                "file (this run: {1})".format(key, value))
            continue
        ref_value = reference[key]
        if isinstance(value, (int, float)) and not isinstance(value, bool) \
                and isinstance(ref_value, (int, float)) \
                and not isinstance(ref_value, bool):
            same = np.isclose(value, ref_value)
        else:
            same = value == ref_value
        if not same:
            warnings.append(
                "WARNING: parameter '{0}' differs from the reference: "
                "reference {1}, this run {2}".format(key, ref_value, value))
    return warnings


def bootstrap_errors(values, n_bootstrap, rng):
    """The errors of the mean and of the standard deviation of values from
    n_bootstrap resamplings."""
    values = np.asarray(values, dtype=float)
    samples = values[rng.integers(0, len(values), (n_bootstrap, len(values)))]
    return (float(np.std(samples.mean(axis=1))),
            float(np.std(samples.std(axis=1))))


def reference_seeds(reference):
    """The seeds (first, last) of a reference; one saved before they were
    recorded used the seeds from 0, one per event."""
    first = reference.get("first_seed", 0)
    return first, reference.get("last_seed", first + reference["n_events"] - 1)


def compare_to_reference(eps2_values, reference, n_sigma, n_bootstrap,
                         mean_tol=None, std_tol=None, seeds=None):
    current = summarize(eps2_values)
    messages = []
    passed = True

    rng = np.random.default_rng(1)
    errors = bootstrap_errors(eps2_values, n_bootstrap, rng)
    if reference.get("eps2_values"):
        ref_errors = bootstrap_errors(reference["eps2_values"],
                                      n_bootstrap, rng)
    else:
        ref_errors = (0., 0.)
    for name, key, error, ref_error, tol in (
            ("mean eps2", "mean", errors[0], ref_errors[0], mean_tol),
            ("std(eps2)", "std", errors[1], ref_errors[1], std_tol)):
        diff = current[key] - reference[key]
        combined = max(float(np.hypot(error, ref_error)), 1e-12)
        ok = abs(diff) <= n_sigma * combined
        message = ("{0}: {1} {2:.5f} -> {3:.5f}, {4:+.1f} sigma (error "
                   "{5:.5f}, fail beyond {6:g} sigma)".format(
                       "{0}", name, reference[key], current[key],
                       diff / combined, combined, n_sigma))
        if tol is not None:
            relative = abs(diff) / max(abs(reference[key]), 1e-12)
            ok = ok and relative <= tol
            message += "; changed by {0:.1%} (tolerance {1:.1%})".format(
                relative, tol)
        passed = passed and ok
        messages.append(message.format("OK" if ok else "FAIL"))

    try:
        from scipy import stats as scipy_stats
        ks_stat, p_value = scipy_stats.ks_2samp(
            eps2_values, reference["eps2_values"])
        ks_ok = p_value > 0.01
        passed = passed and ks_ok
        messages.append(
            "{0}: two-sample KS test p-value = {1:.4f} (fail if <= 0.01, "
            "D={2:.4f})".format("OK" if ks_ok else "FAIL", p_value, ks_stat))
    except ImportError:
        messages.append("(scipy not available -- skipped KS distribution test)")

    if seeds is not None:
        ref_first, ref_last = reference_seeds(reference)
        independent = seeds[1] < ref_first or seeds[0] > ref_last
        passed = passed and independent
        messages.append("{0}: seeds {1}-{2}, reference {3}-{4}{5}".format(
            "OK" if independent else "FAIL", seeds[0], seeds[1], ref_first,
            ref_last, "" if independent else " (overlap: not an independent "
            "check)"))

    return passed, "\n".join(messages)


# --------------------------------------------------------------------------
# environment validation
# --------------------------------------------------------------------------

def validate_environment(ipglasma_path, ipglasma_binary, template_path,
                          qs_table_path, require_run_inputs=True):
    """Check that everything needed is in place. When require_run_inputs
    is False (--plot-only), only the utilities/ helpers are required --
    no new ipglasma runs are going to be launched, so the binary, input
    template and Qs table are not needed."""
    problems = []
    if require_run_inputs:
        if not os.path.isfile(ipglasma_binary):
            problems.append(
                "ipglasma binary not found at '{0}'. Build it first, e.g.:\n"
                "    ./compile_gnu_openmp.sh".format(ipglasma_binary))
        elif not os.access(ipglasma_binary, os.X_OK):
            problems.append(
                "'{0}' exists but is not executable.".format(ipglasma_binary))
        if not os.path.isfile(template_path):
            problems.append("input template not found at '{0}'.".format(template_path))
        if not os.path.isfile(qs_table_path):
            problems.append("Qs table file not found at '{0}'.".format(qs_table_path))
    for helper in ("read_tmunu.py", "eccentricity.py"):
        helper_path = os.path.join(ipglasma_path, "utilities", helper)
        if not os.path.isfile(helper_path):
            problems.append(
                "utilities helper not found at '{0}'. --ipglasma-path "
                "must point at an IP-Glasma checkout containing a "
                "utilities/ directory (currently '{1}').".format(
                    helper_path, ipglasma_path))
    if problems:
        raise FileNotFoundError("\n".join(problems))


# Stray output files IP-Glasma writes relative to its process cwd (see
# Init.cpp / Evolution.cpp) using "<name><eventId>.<ext>"
#  Each event's run
# is launched with cwd=worker_dir, so these normally land inside that
# event's own working directory and are removed along with it; this is a
# safety net for the current working directory the script itself was
# invoked from, in case any of them end up written there instead (e.g. an
# ipglasma run outside of this script's per-event worker dirs).
_STRAY_OUTPUT_GLOBS = [
    "NcollList*.dat", "NpartList*.dat", "NgluonEstimators*.dat",
    "usedParameters*.dat", "eccentricities*.dat", "gluonMultiplicity*.json",
    "run.log",
    # generate_temp_input()'s temp input file: normally removed by
    # run_seed()'s own finally block right after ipglasma exits, but a run
    # killed mid-event (e.g. Ctrl-C) can leave one behind -- --clean sweeps
    # those up too.
    "input_ecc_*.in",
]


def cleanup_stray_output_files(directory):
    """Remove any of _STRAY_OUTPUT_GLOBS found directly in `directory`."""
    for pattern in _STRAY_OUTPUT_GLOBS:
        for path in glob.glob(os.path.join(directory, pattern)):
            if os.path.isfile(path):
                os.remove(path)


def clean_datadir(datadir):
    """--clean support: remove stray temporary/output files from every
    seed_*/ working directory under `datadir` (and from `datadir` itself),
    keeping the Tmunu snapshots and any results/plots. Meant to tidy up
    after a run was interrupted (e.g. Ctrl-C) before it reached its own
    end-of-run cleanup. Does not touch failed-event or missing seed_*
    directories beyond removing the same stray files from them."""
    if not os.path.isdir(datadir):
        print("Nothing to clean: {0} does not exist.".format(datadir))
        return
    worker_dirs = sorted(
        d for d in glob.glob(os.path.join(datadir, "seed_*")) if os.path.isdir(d))
    for worker_dir in worker_dirs:
        cleanup_stray_output_files(worker_dir)
        print("Cleaned {0}".format(worker_dir))
    cleanup_stray_output_files(datadir)
    print("Cleaned {0} working director{1} under {2}".format(
        len(worker_dirs), "y" if len(worker_dirs) == 1 else "ies", datadir))


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxevents", type=int, default=1000,
                         help="number of independent IP-Glasma events "
                              "(random seeds) to run (default 1000)")
    parser.add_argument("--datadir", default="./eccentricity_test_output/",
                         help="directory for per-event working dirs and plots")
    parser.add_argument("--ipglasma-path", default="../..",
                         help="directory containing the ipglasma binary, "
                              "the Qs table and utilities/ (default: the "
                              "repository root, ../..)")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--input-template", default="input_eccentricity",
                         help="IP-Glasma input file to use as a template ")
    parser.add_argument("--size", type=int, default=128,
                         help="lattice size override, default 256")
    parser.add_argument("--L", type=float, default=30,
                         help="transverse box size override in fm "
                              "(default: 30 fm)")
    parser.add_argument("--maxtime", type=float, default=None,
                         help="evolution end time override in fm/c "
                              "(default: keep the template's value)")
    parser.add_argument("--jimwlk", action="store_true",
                         help="keep JIMWLK small-x evolution enabled "
                              "(slower, more complete physics test). "
                              "Default: disabled for speed.")
    parser.add_argument("--max-workers", type=int, default=None,
                         help="events run at a time (default: the logical "
                              "CPUs available divided by OMP_NUM_THREADS), "
                              "each with OMP_NUM_THREADS threads (default "
                              "1); without --slurm, several events at a "
                              "time run with OMP_PROC_BIND=false")
    parser.add_argument("--keep-logs", action="store_true",
                         help="keep the working directories of all events; "
                              "by default, those of the successful events "
                              "are removed at the end of the run (those of "
                              "failed events are always kept)")
    parser.add_argument("--slurm", action="store_true",
                         help="wrap each ipglasma invocation with "
                              "'srun --ntasks=1 --cpus-per-task=N', N = "
                              "OMP_NUM_THREADS (default 1)")
    parser.add_argument("--plot-only", action="store_true",
                         help="skip launching new ipglasma runs; instead "
                              "rebuild the eps2 distribution from the "
                              "seed_*/ working directories already present "
                              "in --datadir (e.g. left behind by an "
                              "interrupted run, or kept with --keep-logs) "
                              "and regenerate the plot / regression check")
    parser.add_argument("--clean", action="store_true",
                         help="don't run or plot anything; just remove "
                              "stray temporary/output files (run.log, "
                              "leftover temp input files, ...) from every "
                              "seed_*/ working directory under --datadir, "
                              "keeping the Tmunu snapshots, then exit -- "
                              "use this to tidy up after a run was killed "
                              "before its own end-of-run cleanup")
    parser.add_argument(
        "--reference-file",
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "eccentricity_reference.json"),
        help="JSON file with a saved baseline eps2 distribution to "
             "regression-check against, default is {0} corresponding to the case without JIMWLK".format(
                 os.path.join(
                     os.path.dirname(os.path.abspath(__file__)),
                     "eccentricity_reference.json")))
    parser.add_argument("--save-reference", action="store_true",
                         help="save this run's eps2 distribution as the "
                              "new reference instead of comparing to it")
    parser.add_argument("--first-seed", type=int, default=None,
                         help="first random seed (default 10000, beyond the "
                              "seeds of the reference; 0 with "
                              "--save-reference)")
    parser.add_argument("--n-sigma", type=float, default=3.0,
                         help="allowed change of mean(eps2) and std(eps2) "
                              "in units of the combined error")
    parser.add_argument("--bootstrap", type=int, default=1000,
                         help="bootstrap resamplings for the errors")
    parser.add_argument("--mean-tol", type=float, default=None,
                         help="optional allowed relative change in "
                              "mean(eps2), in addition")
    parser.add_argument("--std-tol", type=float, default=None,
                         help="optional allowed relative change in "
                              "std(eps2), in addition")
    args = parser.parse_args()

    first_seed = args.first_seed
    if first_seed is None:
        first_seed = 0 if args.save_reference else 10000
    seed_list = range(first_seed, first_seed + args.maxevents)

    if args.clean:
        # Just tidy up an existing --datadir and exit -- no ipglasma
        # checkout / binary / template needed for this.
        clean_datadir(os.path.abspath(args.datadir))
        sys.exit(0)

    ipglasma_path = os.path.abspath(args.ipglasma_path)
    ipglasma_binary = os.path.join(ipglasma_path, args.ipglasma_cmd)
    # the template in the working directory, else next to the script
    template_path = args.input_template
    if not os.path.isfile(template_path):
        template_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), args.input_template)

    qs_table_name = read_input_value(template_path, "nucleusQsTableFileName") \
        if os.path.isfile(template_path) else None
    qs_table_path = os.path.join(ipglasma_path, qs_table_name) \
        if qs_table_name else os.path.join(ipglasma_path, "qs2Adj_vs_Tp_vs_Y_240.in")

    validate_environment(
        ipglasma_path, ipglasma_binary, template_path, qs_table_path,
        require_run_inputs=not args.plot_only)

    # Record the values the events actually ran with
    # None means "keep the template's value", so in that case fall back to
    # reading those from the template instead of storing null.
    def effective_value(override, key):
        if override is not None:
            return override
        if not os.path.isfile(template_path):
            return None
        value = read_input_value(template_path, key)
        return float(value) if value is not None else None

    run_meta = {
        "size": args.size,
        "L": effective_value(args.L, "L"),
        "maxTime": effective_value(args.maxtime, "maxTime"),
        "jimwlk": args.jimwlk,
    }

    reference = None if args.save_reference else load_reference(args.reference_file)
    # Warn up front, before spending time on the run, if this run's
    # parameters don't match those the reference was recorded with.
    if reference is not None:
        for warning in check_reference_parameters(run_meta, reference):
            print(warning)

    datadir = os.path.abspath(args.datadir)
    os.makedirs(datadir, exist_ok=True)

    if args.plot_only:
        print("--plot-only: reading already-completed events from {0} "
              "(no new ipglasma runs) ...".format(datadir))
        results = collect_from_existing_datadir(datadir, ipglasma_path)
    else:
        max_workers = args.max_workers or max(1, min(
            available_cpus() // threads_per_event(), args.maxevents))
        if max_workers > 1 and not args.slurm:
            # binding (OMP_PROC_BIND or OMP_PLACES, as in the main README)
            # would put the threads of all events on the same first cores;
            # with --slurm, srun gives each event its own CPUs
            os.environ["OMP_PROC_BIND"] = "false"
        print("Running {0} IP-Glasma events across up to {1} parallel workers "
              "(size={2}, jimwlk={3}), lattice N={4}, L={5} fm ...".format(
                  args.maxevents, max_workers, args.size, args.jimwlk, args.size, args.L))

        results = []
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    run_seed, seed, ipglasma_path, template_path, qs_table_path,
                    datadir, ipglasma_binary, args.size, args.L, args.maxtime,
                    args.jimwlk, args.slurm): seed
                for seed in seed_list
            }
            for future in as_completed(futures):
                seed = futures[future]
                try:
                    results.append(future.result())
                except Exception as exc:  # noqa: BLE001
                    results.append({
                        "seed": seed, "success": False,
                        "error": "worker raised: {0}".format(exc),
                    })

    results.sort(key=lambda r: r["seed"])
    successes = [r for r in results if r["success"]]
    failures = [r for r in results if not r["success"]]

    if failures:
        failed_path = os.path.join(datadir, "failed_seeds.txt")
        with open(failed_path, "w") as f:
            for r in failures:
                f.write("seed {0}: {1}\n".format(r["seed"], r["error"]))
        print("{0}/{1} events FAILED (see {2}):".format(
            len(failures), len(results), failed_path))
        for r in failures:
            print("  seed {0}: {1}".format(r["seed"], r["error"]))

    if not successes:
        if args.plot_only:
            print("No successful events found under {0} -- aborting.".format(
                datadir))
        else:
            print("No successful events -- aborting.")
        sys.exit(1)

    eps2_values = [r["eps2"] for r in successes]
    seed_numbers = [r["seed"] for r in results if isinstance(r["seed"], int)]
    seeds = (min(seed_numbers), max(seed_numbers)) if seed_numbers else None
    stats = summarize(eps2_values)
    print("epsilon_2 over {0} successful events: "
          "mean={1:.5f} std={2:.5f} median={3:.5f} "
          "min={4:.5f} max={5:.5f}".format(
              stats["n_events"], stats["mean"], stats["std"],
              stats["median"], stats["min"], stats["max"]))

    plot_path = plot_distribution(eps2_values, datadir, reference=reference)
    print("Wrote epsilon_2 distribution plot to {0}".format(plot_path))

    exit_code = 1 if failures else 0

    if args.save_reference:
        save_reference(args.reference_file, eps2_values, run_meta, seeds)
        print("Saved new reference distribution to {0}".format(
            args.reference_file))
    else:
        # Not (re-)saving the reference -- still write this run's
        # distribution (same format, including the "_command" needed to
        # reproduce it) to datadir, so it isn't lost once seed_*/ working
        # directories are cleaned up.
        results_file = os.path.join(datadir, "eccentricity_results.json")
        save_reference(results_file, eps2_values, run_meta, seeds)
        print("Saved this run's eps2 distribution to {0}".format(results_file))

        if reference is not None:
            passed, message = compare_to_reference(
                eps2_values, reference, args.n_sigma, args.bootstrap,
                args.mean_tol, args.std_tol, seeds)
            print(message)
            if not passed:
                exit_code = 1
                print("REGRESSION CHECK FAILED: epsilon_2 distribution drifted "
                      "beyond its errors relative to {0}, or the seeds "
                      "overlap".format(args.reference_file))
            else:
                print("Regression check passed.")
        else:
            print("No reference file found at {0}. Rerun with --save-reference "
                  "to record a baseline once you've confirmed this run's "
                  "physics is correct.".format(args.reference_file))

    if not args.keep_logs and not args.plot_only:
        print("Removing the working directories of the {0} successful "
              "events under {1} ...".format(len(successes), datadir))
        for r in successes:
            shutil.rmtree(r["worker_dir"], ignore_errors=True)




    sys.exit(exit_code)


if __name__ == "__main__":
    main()

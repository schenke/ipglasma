#!/usr/bin/env python3
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
exits non-zero if it drifts beyond the configured tolerance.

Usage:
    # one-time: record a baseline from a known-good build
    python3 parallel_test_eccentricity.py --maxevents 40 --save-reference

    # later: check that nothing changed
    python3 parallel_test_eccentricity.py --maxevents 40
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
    for line in lines:
        tokens = line.split()
        if len(tokens) >= 2 and tokens[0] in remaining:
            out_lines.append("{0} {1}\n".format(tokens[0], remaining.pop(tokens[0])))
        else:
            out_lines.append(line)
    # Any override keys not already present in the template are appended.
    for key, value in remaining.items():
        out_lines.append("{0} {1}\n".format(key, value))

    fd, temp_path = tempfile.mkstemp(
        prefix="input_ecc_{0}_".format(tag), suffix=".in", dir=working_dir)
    with os.fdopen(fd, "w") as f:
        f.writelines(out_lines)
    return temp_path


# --------------------------------------------------------------------------
# running one event
# --------------------------------------------------------------------------

def run_command(cmd, cwd, log_path, slurm=False):
    if slurm:
        cmd = ["srun", "--ntasks=1", "--cpus-per-task=1"] + list(cmd)
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
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
        "NucleusQsTableFileName": qs_table_path,
        "writeOutputsToHDF5": 0,
        "writeWilsonLines": 0,
        "writeInitialWilsonLines": 0,
        "writeEvolution": 0,
        # writeOutputs must have its "value/4 == 1" bit set for
        # MyEigen::flowVelocity4DImpl to write any Tmunu snapshot at all
        # (see src/MyEigen.cpp); it returns immediately if writeOutputs<=0.
        # writeOutputs==5 additionally enables the extra intermediate
        # snapshots at tau=0.1,0.2,0.3,0.4 fm/c (see Evolution::run) which
        # we don't need here, so 4 is the minimal value that writes just
        # the final-time snapshot at it==itmax.
        "writeOutputs": 4,
        "writeEpsilonUHydro": 0,
        "writeTmunuBinary": 1,
        "useJIMWLK": 1 if use_jimwlk else 0,
        "useSeedList": 0,
        "useTimeForSeed": 0,
        "saveSnapshots": 0,
    }
    if L is not None:
        overrides["L"] = L
        overrides["LOutput"] = L
    if maxtime is not None:
        overrides["maxtime"] = maxtime

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
    out_path_png = os.path.join(datadir, "epsilon2_distribution.png")
    fig.savefig(out_path_png, dpi=150)
    plt.close(fig)
    return out_path_png


def load_reference(reference_file):
    if not os.path.exists(reference_file):
        return None
    with open(reference_file, "r") as f:
        return json.load(f)


def save_reference(reference_file, eps2_values, run_meta):
    payload = dict(summarize(eps2_values))
    payload["eps2_values"] = [float(v) for v in eps2_values]
    payload.update(run_meta)
    with open(reference_file, "w") as f:
        json.dump(payload, f, indent=2)


def compare_to_reference(eps2_values, reference, mean_tol, std_tol):
    current = summarize(eps2_values)
    messages = []
    passed = True

    rel_mean = abs(current["mean"] - reference["mean"]) / max(
        abs(reference["mean"]), 1e-12)
    if rel_mean > mean_tol:
        passed = False
        messages.append(
            "FAIL: mean eps2 changed by {0:.1%} (tolerance {1:.1%}): "
            "{2:.5f} -> {3:.5f}".format(
                rel_mean, mean_tol, reference["mean"], current["mean"]))
    else:
        messages.append(
            "OK: mean eps2 within tolerance ({0:.1%} <= {1:.1%}): "
            "{2:.5f} -> {3:.5f}".format(
                rel_mean, mean_tol, reference["mean"], current["mean"]))

    rel_std = abs(current["std"] - reference["std"]) / max(
        abs(reference["std"]), 1e-12)
    if rel_std > std_tol:
        passed = False
        messages.append(
            "FAIL: std(eps2) changed by {0:.1%} (tolerance {1:.1%}): "
            "{2:.5f} -> {3:.5f}".format(
                rel_std, std_tol, reference["std"], current["std"]))
    else:
        messages.append(
            "OK: std(eps2) within tolerance ({0:.1%} <= {1:.1%}): "
            "{2:.5f} -> {3:.5f}".format(
                rel_std, std_tol, reference["std"], current["std"]))

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
# Init.cpp / Evolution.cpp) using "<name><eventId>.<ext>" -- eventId is
# usually 0, but is not guaranteed to be, hence the glob. Each event's run
# is launched with cwd=worker_dir, so these normally land inside that
# event's own working directory and are removed along with it; this is a
# safety net for the current working directory the script itself was
# invoked from, in case any of them end up written there instead (e.g. an
# ipglasma run outside of this script's per-event worker dirs).
_STRAY_OUTPUT_GLOBS = [
    "NcollList*.dat", "NpartList*.dat", "NgluonEstimators*.dat",
    "usedParameters*.dat", "eccentricities*.dat", "gluonMultiplicity*.json",
    "run.log",
]


def cleanup_stray_output_files(directory):
    """Remove any of _STRAY_OUTPUT_GLOBS found directly in `directory`."""
    for pattern in _STRAY_OUTPUT_GLOBS:
        for path in glob.glob(os.path.join(directory, pattern)):
            if os.path.isfile(path):
                os.remove(path)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxevents", type=int, default=24,
                         help="number of independent IP-Glasma events "
                              "(random seeds) to run in parallel")
    parser.add_argument("--datadir", default="./eccentricity_test_output/",
                         help="directory for per-event working dirs and plots")
    parser.add_argument("--ipglasma-path", default=".",
                         help="directory containing the ipglasma binary, "
                              "input template and Qs table")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--input-template", default="input",
                         help="IP-Glasma input file to use as a template ")
    parser.add_argument("--size", type=int, default=128,
                         help="lattice size override (evolution + output), "
                              "smaller than the physics-quality default of "
                              "256/512 to keep this test fast")
    parser.add_argument("--L", type=float, default=None,
                         help="transverse box size override in fm "
                              "(default: keep the template's value)")
    parser.add_argument("--maxtime", type=float, default=None,
                         help="evolution end time override in fm/c "
                              "(default: keep the template's value)")
    parser.add_argument("--jimwlk", action="store_true",
                         help="keep JIMWLK small-x evolution enabled "
                              "(slower, more complete physics test). "
                              "Default: disabled for speed.")
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--keep-logs", action="store_true",
                         help="keep per-event working directories "
                              "(including failed ones, which are always kept)")
    parser.add_argument("--slurm", action="store_true",
                         help="wrap each ipglasma invocation with "
                              "'srun --ntasks=1 --cpus-per-task=1'")
    parser.add_argument("--plot-only", action="store_true",
                         help="skip launching new ipglasma runs; instead "
                              "rebuild the eps2 distribution from the "
                              "seed_*/ working directories already present "
                              "in --datadir (e.g. left behind by an "
                              "interrupted run with --keep-logs) and "
                              "regenerate the plot / regression check")
    parser.add_argument(
        "--reference-file",
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "eccentricity_reference.json"),
        help="JSON file with a saved baseline eps2 distribution to "
             "regression-check against")
    parser.add_argument("--save-reference", action="store_true",
                         help="save this run's eps2 distribution as the "
                              "new reference instead of comparing to it")
    parser.add_argument("--mean-tol", type=float, default=0.15,
                         help="allowed relative change in mean(eps2)")
    parser.add_argument("--std-tol", type=float, default=0.25,
                         help="allowed relative change in std(eps2)")
    args = parser.parse_args()

    ipglasma_path = os.path.abspath(args.ipglasma_path)
    ipglasma_binary = os.path.join(ipglasma_path, args.ipglasma_cmd)
    template_path = os.path.join(args.input_template)

    qs_table_name = read_input_value(template_path, "NucleusQsTableFileName") \
        if os.path.isfile(template_path) else None
    qs_table_path = os.path.join(ipglasma_path, qs_table_name) \
        if qs_table_name else os.path.join(ipglasma_path, "qs2Adj_vs_Tp_vs_Y_200.in")

    validate_environment(
        ipglasma_path, ipglasma_binary, template_path, qs_table_path,
        require_run_inputs=not args.plot_only)

    datadir = os.path.abspath(args.datadir)
    os.makedirs(datadir, exist_ok=True)

    if args.plot_only:
        print("--plot-only: reading already-completed events from {0} "
              "(no new ipglasma runs) ...".format(datadir))
        results = collect_from_existing_datadir(datadir, ipglasma_path)
    else:
        max_workers = args.max_workers or min(cpu_count(), args.maxevents)
        print("Running {0} IP-Glasma events across up to {1} parallel workers "
              "(size={2}, jimwlk={3}) ...".format(
                  args.maxevents, max_workers, args.size, args.jimwlk))

        results = []
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    run_seed, seed, ipglasma_path, template_path, qs_table_path,
                    datadir, ipglasma_binary, args.size, args.L, args.maxtime,
                    args.jimwlk, args.slurm): seed
                for seed in range(args.maxevents)
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
    stats = summarize(eps2_values)
    print("epsilon_2 over {0} successful events: "
          "mean={1:.5f} std={2:.5f} median={3:.5f} "
          "min={4:.5f} max={5:.5f}".format(
              stats["n_events"], stats["mean"], stats["std"],
              stats["median"], stats["min"], stats["max"]))

    reference = None if args.save_reference else load_reference(args.reference_file)
    plot_path = plot_distribution(eps2_values, datadir, reference=reference)
    print("Wrote epsilon_2 distribution plot to {0}".format(plot_path))

    exit_code = 1 if failures else 0

    run_meta = {
        "size": args.size, "L": args.L, "maxtime": args.maxtime,
        "jimwlk": args.jimwlk,
    }
    if args.save_reference:
        save_reference(args.reference_file, eps2_values, run_meta)
        print("Saved new reference distribution to {0}".format(
            args.reference_file))
    elif reference is not None:
        passed, message = compare_to_reference(
            eps2_values, reference, args.mean_tol, args.std_tol)
        print(message)
        if not passed:
            exit_code = 1
            print("REGRESSION CHECK FAILED: epsilon_2 distribution drifted "
                  "beyond tolerance relative to {0}".format(args.reference_file))
        else:
            print("Regression check passed.")
    else:
        print("No reference file found at {0}. Rerun with --save-reference "
              "to record a baseline once you've confirmed this run's "
              "physics is correct.".format(args.reference_file))

    if not args.keep_logs and not args.plot_only:
        for r in successes:
            shutil.rmtree(r["worker_dir"], ignore_errors=True)
        # failed events' working directories (containing run.log) are
        # always kept for debugging, regardless of --keep-logs.
        # In --plot-only mode the working directories are the input to
        # this run (not freshly produced by it), so they are left alone
        # regardless of --keep-logs -- otherwise --plot-only would delete
        # the very data it was asked to replot.
        cleanup_stray_output_files(os.getcwd())

    sys.exit(exit_code)


if __name__ == "__main__":
    main()

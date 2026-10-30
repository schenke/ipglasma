#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Parallel version of the vector meson production validation script.

This script parallelizes the loop over seeds so that multiple IP-Glasma runs execute
simultaneously. Each worker gets its own temporary input file, so the simulations are
independent and do not share configuration state.

This tool is used to test the IP-Glasma+JIMWLK+SubnucleonDiffraction framework for vector meson production.
It runs a series of simulations with different random seeds, computes the coherent and incoherent cross sections
from the amplitudes of the events, and compares them to a reference made with a large number of events
(jpsi_reference.json, made with --save-reference).

The computed cross section should also match the published results in https://arxiv.org/pdf/2207.03712 Fig. 1 ("CGC+shape fluct") curve,
except that by default the code uses a smaller lattice for performance reasons, which may lead to small deviations in the cross section.
The lattice size can be increased by modifying the input file "input_vm_proton".

The script also generates plots  color field (tr V(x)) for selected events.

Example usage (from validations/vector_meson_production/):
    python3 parallel_test_vector_meson_production.py --maxevents 100 --datadir path_for_cross_sections --subnucleondiffraction-path path_to_subnucleondiffraction_directory

Maximum number of parallel workers can be specified with --max-workers. By default, it uses the logical CPUs available, one OpenMP
thread (OMP_NUM_THREADS, default 1) per event.

The comparison plot will be saved in the specified directory as "cross_section.pdf". The color field plots will be saved as "output_<event_id>.pdf".

Seeds: the run uses the seeds --first-seed, --first-seed + 1, ... (two events,
projectile and target, per seed). The default, 10000, keeps a test run
independent of the reference, which is made with the seeds from 0 (the
default with --save-reference). The spectra are copied to the data
directory as spectra_1 ... spectra_<maxevents>, whatever the seeds.

Results: jpsi_results.json in the data directory has the cross sections
dsigma/dt and their integrals over |t| <= --check-tmax, with the errors of a
bootstrap over the events (--bootstrap resamplings), the seeds, the versions
of IP-Glasma and subnucleondiffraction and the settings; the reference has
the same format.

Check (jpsi_check.json): the integrated coherent and incoherent cross
sections must agree with those of the reference within --n-sigma errors,
where the error is the one expected if they agree: the relative error of the
reference, for the events of the reference and of the run; the spectra must
be finite and none an outlier (its largest |A|^2 more than --outlier-zmax
robust standard deviations of log10 |A|^2 from the median); and the seeds of
the run must not overlap those of the reference. Otherwise the script exits
with 1.

Use --plot-only to skip the simulation and only compute the cross sections, the plots and the check from the spectra in the specified directory.

The IP-Glasma input file used as a template can be specified with --input-template (default: input_vm_proton). The Qs table file
it points to (NucleusQsTableFileName) is read from that template rather than hardcoded, so pointing --input-template at a
differently configured input file also picks up its Qs table.

Before running this, make sure you have downloaded and built the subnucleondiffraction code from https://github.com/hejajama/subnucleondiffraction
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

import matplotlib.pyplot as plt
import numpy as np


subnucleondiffraction_path = "./subnucleondiffraction/"
ipglasma_path = "../../"
ipglasma_cmd = "ipglasma"
plot_ids = [1, 5]  # run-local event numbers, see run_seed()
default_xpom = "0.001705"
reference_file_name = "jpsi_reference.json"
default_first_seed = 10000
# dsigma/dt = |A|^2/(16 pi) in mb/GeV^2, as in python/cross_section.py of
# subnucleondiffraction
GEVSQ_TO_MB = 1.0e7 / (5.068 * 5.068) * 1e-6 / (16.0 * np.pi)
lattice_N = 512 # N*N lattice
lattice_L = 4 # fm
# the options of subnucleondiffraction, besides the wave function file and
# the Wilson lines
subnucleondiffraction_options = ["-Q2", "0", "-maxt", "4.01", "-tstep", "0.05", "-mcintpoints", "5e5"]


def _repo_root():
    return os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


def _script_dir():
    return os.path.dirname(os.path.abspath(__file__))


def read_input_value(template_path, key):
    """Return the (string) value of a "key value" line in an IP-Glasma
    input file, or None if the key is not present."""
    with open(template_path, "r") as f:
        for line in f:
            tokens = line.split()
            if len(tokens) >= 2 and tokens[0] == key:
                return tokens[1]
    return None


def generate_temp_input(source_input_file="input_vm_proton", seed=0, x_pom="0.001705", working_dir=".", qs_table_path=None):
    with open(source_input_file, "r") as f:
        lines = f.readlines()

    temp_fd, temp_path = tempfile.mkstemp(prefix="input_vm_proton_", suffix=".in", dir=working_dir)
    os.close(temp_fd)

    with open(temp_path, "w") as f:
        for line in lines:
            if line.lstrip().startswith("seed"):
                f.write(f"seed {seed}\n")
            elif line.lstrip().startswith("projectileX"):
                f.write(f"projectileX {x_pom}\n")
            elif line.lstrip().startswith("targetX"):
                f.write(f"targetX {x_pom}\n")
            elif line.lstrip().startswith("nucleusQsTableFileName") and qs_table_path:
                f.write(f"nucleusQsTableFileName {qs_table_path}\n")
            elif line.lstrip().startswith("nuclearConfigurationsPath"):
                f.write(f"nuclearConfigurationsPath {ipglasma_path}/nucleusConfigurations/\n")
            elif line.lstrip().startswith("size "):
                f.write(f"size {lattice_N}\n")
            elif line.lstrip().startswith("L "):
                f.write(f"L {lattice_L}\n")
            else:
                f.write(line)

    return temp_path


def remove_temp_input(temp_input_path):
    if temp_input_path and os.path.exists(temp_input_path):
        os.remove(temp_input_path)


def format_xvalue_for_filename(x_value):
    """Match C++ formatting: std::scientific << std::setprecision(5)."""
    return f"{float(x_value):.5e}"


def wilson_line_path(base_dir, x_value, event_id):
    return os.path.join(base_dir, f"WilsonLine_x_{format_xvalue_for_filename(x_value)}_{event_id}")


def ensure_wilson_line_file(base_dir, x_value, event_id):
    """Return canonical Wilson line file path, moving legacy Final_x output if needed."""
    target = wilson_line_path(base_dir, x_value, event_id)
    if os.path.exists(target):
        return target

    raise FileNotFoundError(
        f"Wilson line file does not exist for event {event_id}: '{target}'" 
    )


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


def run_command(cmd, cwd, stdout_path=None, stderr_path=None, slurm=False):
    print(f"Running: {' '.join(str(part) for part in cmd)}", " stdout:", stdout_path, " stderr:", stderr_path)
    srun_cmd = ["srun", "--ntasks=1", f"--cpus-per-task={threads_per_event()}"] + cmd
    if not slurm:
        srun_cmd = cmd
    # one thread per event unless OMP_NUM_THREADS says otherwise: the events
    # already run in parallel
    env = dict(os.environ, OMP_NUM_THREADS=str(threads_per_event()))
    result = subprocess.run(srun_cmd, cwd=cwd, capture_output=True, text=True, env=env)
    if stdout_path:
        with open(stdout_path, "w") as stdout_handle:
            stdout_handle.write(result.stdout)
    if stderr_path:
        with open(stderr_path, "w") as stderr_handle:
            stderr_handle.write(result.stderr)

    return result


def spectra_files(datadir):
    """The run-local numbers of the spectra_<n> files in datadir."""
    numbers = []
    for name in os.listdir(datadir):
        if name.startswith("spectra_") and name[len("spectra_"):].isdigit():
            numbers.append(int(name[len("spectra_"):]))
    return sorted(numbers)


def subnucleondiffraction_version(path):
    """The version of subnucleondiffraction that wrote a spectra file, from
    its "# Git version" line without the date of the build, or None."""
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("# Git version"):
                return line[len("# Git version"):].split(" build ")[0].strip()
    return None


def ipglasma_version(ipglasma_binary):
    """The version that "ipglasma --version" prints, or None."""
    try:
        result = subprocess.run([ipglasma_binary, "--version"], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return result.stdout.strip() or None if result.returncode == 0 else None


def read_spectra(datadir, maxconf=None):
    """The spectra_<n> files of datadir (with maxconf, only n <= maxconf):
    returns their numbers, t, the amplitudes (event, t, column; the columns
    after t: Re A and Im A of the transverse, then of the longitudinal
    polarization), the spectra that are not readable, not finite or not on
    the t values of the first (number, reason), and the versions of
    subnucleondiffraction that wrote them."""
    numbers, data, invalid, versions = [], [], [], set()
    for number in spectra_files(datadir):
        if maxconf is not None and number > maxconf:
            continue
        path = os.path.join(datadir, f"spectra_{number}")
        try:
            values = np.loadtxt(path)
        except (OSError, ValueError) as exc:  # e.g. a dangling link
            invalid.append((number, f"not readable: {exc}"))
            continue
        if values.ndim != 2 or values.shape[1] < 3:
            invalid.append((number, "no t values and amplitudes"))
        elif not np.isfinite(values).all():
            invalid.append((number, "not finite"))
        elif data and (values.shape != data[0].shape or not np.allclose(values[:, 0], data[0][:, 0], atol=1e-5)):
            invalid.append((number, "other t values than the first spectra"))
        else:
            numbers.append(number)
            data.append(values)
            version = subnucleondiffraction_version(path)
            if version:
                versions.add(version)
    if not data:
        return numbers, None, None, invalid, sorted(versions)
    data = np.array(data)
    return numbers, data[0][:, 0], data[:, :, 1:], invalid, sorted(versions)


def find_outliers(numbers, amplitudes, zmax):
    """The spectra whose scale, the largest |A|^2 (summed over the columns)
    over t, is not positive, or more than zmax robust standard deviations
    (MAD / 0.6745) of log10(scale) from its median, as the outlier filter of
    utilities/compute_tdiffCrossSections.py of IPGlasmaFramework: returns
    [(number, z)], z None for a scale that is not positive."""
    scale = np.max(np.sum(amplitudes ** 2, axis=2), axis=1)
    positive = scale > 0
    outliers = [(numbers[i], None) for i in np.where(~positive)[0]]
    x = np.log10(scale[positive])
    if len(x) >= 3:
        median = np.median(x)
        mad = np.median(np.abs(x - median))
        if mad > 0:
            z = 0.6745 * (x - median) / mad
            outliers += [(numbers[i], float(zi)) for i, zi in zip(np.where(positive)[0], z) if abs(zi) > zmax]
    return sorted(outliers)


def analyze(t, amplitudes, tmax, n_bootstrap):
    """dsigma/dt, coherent and incoherent, in mb/GeV^2, and their integrals
    over |t| <= tmax (trapezoid rule), in mb, with the errors of a bootstrap
    over the events. As python/cross_section.py of subnucleondiffraction:
    the coherent cross section is |<Re A>|^2 and the incoherent one the
    variance of Re A and Im A, of the transverse polarization."""
    real, imag = amplitudes[:, :, 0], amplitudes[:, :, 1]
    inside = t <= tmax + 1e-9
    tt = t[inside]
    weights = np.zeros(len(tt))
    weights[1:] += 0.5 * np.diff(tt)
    weights[:-1] += 0.5 * np.diff(tt)

    def cross_sections(events):
        coherent = np.mean(real[events], axis=0) ** 2
        incoherent = np.var(real[events], axis=0) + np.var(imag[events], axis=0)
        return np.array([coherent, incoherent]) * GEVSQ_TO_MB

    n_events = len(real)
    value = cross_sections(np.arange(n_events))
    rng = np.random.default_rng(1)
    samples = np.array([cross_sections(rng.integers(0, n_events, n_events)) for _ in range(n_bootstrap)])
    error = samples.std(axis=0)
    integral = value[:, inside] @ weights
    integral_error = (samples[:, :, inside] @ weights).std(axis=0)

    result = {"n_events": n_events, "bootstrap": n_bootstrap, "t": t.tolist()}
    for k, name in enumerate(("coherent", "incoherent")):
        result[f"{name}_mb_per_GeV2"] = value[k].tolist()
        result[f"{name}_error_mb_per_GeV2"] = error[k].tolist()
    result["integrated"] = {"tmax": tmax}
    for k, name in enumerate(("coherent", "incoherent")):
        result["integrated"][f"{name}_mb"] = float(integral[k])
        result["integrated"][f"{name}_error_mb"] = float(integral_error[k])
    return result


def reproduce_command():
    """The command line of this run, recorded in the results and the
    reference as their "_command"."""
    return "python3 {0}".format(" ".join([sys.argv[0]] + sys.argv[1:]))


def results_payload(result, seeds, versions, settings):
    """The results of a run as written to jpsi_results.json and, with
    --save-reference, to the reference: the command, the versions, the
    settings, the seeds and the cross sections."""
    payload = {"_command": reproduce_command(), "versions": versions, "settings": settings,
               "first_seed": seeds[0] if seeds else None, "last_seed": seeds[1] if seeds else None}
    payload.update(result)
    return payload


def load_reference(reference_path):
    if not os.path.exists(reference_path):
        return None
    with open(reference_path) as f:
        return json.load(f)


def check_against_reference(result, reference, n_sigma, seeds=None):
    """The checks of a run against the reference: its cross sections
    integrated over |t| <= tmax agree with those of the reference within
    n_sigma errors, and its seeds do not overlap those of the reference.
    The error is that of the ratio expected if the two agree: the relative
    bootstrap error of the reference, scaled to the events of the reference
    and of the run, sqrt(1 + N_ref / N_run). The run's own error is not
    used: it shrinks with its cross sections, so that runs that fluctuate
    low would fail more often (against a reference of 2000 other events of
    the same code, of the runs with 32 events 1.4% (coherent) and 5%
    (incoherent) fail beyond 3 sigma with it, all low, and 0.5% and 0.4%
    with this error). Returns (checks, summary)."""
    checks, summary = {}, {}
    run, ref = result["integrated"], reference["integrated"]
    if not np.isclose(run["tmax"], ref["tmax"]):
        print(f"  FAIL the reference has the cross sections integrated over |t| <= {ref['tmax']:g} GeV^2, "
              f"not {run['tmax']:g} (--check-tmax)", flush=True)
        checks["the reference has the integrals over the same |t|"] = False
        return checks, summary
    scale = np.sqrt(1.0 + reference["n_events"] / result["n_events"])
    print(f"Integrated over |t| <= {run['tmax']:g} GeV^2, {result['n_events']} events against the "
          f"{reference['n_events']} of the reference (bootstrap with {result['bootstrap']} resamplings):", flush=True)
    for name in ("coherent", "incoherent"):
        value, ref_value = run[f"{name}_mb"], ref[f"{name}_mb"]
        ratio = value / ref_value
        error = ref[f"{name}_error_mb"] / ref_value * scale
        own = run[f"{name}_error_mb"] / value
        ok = bool(abs(ratio - 1.0) <= n_sigma * error)
        checks[f"{name} cross section agrees with the reference within {n_sigma:g} sigma"] = ok
        summary[name] = {"sigma_mb": value, "reference_mb": ref_value, "ratio": ratio, "ratio_error": float(error),
                         "relative_error_of_the_run": own}
        print(f"  {'OK  ' if ok else 'FAIL'} {name}: ratio {ratio:.3f} +- {error:.3f} ({(ratio - 1.0) / error:+.1f} sigma; "
              f"relative error of the run alone {own:.3f})", flush=True)
    if seeds is not None:
        ref_first, ref_last = reference["first_seed"], reference["last_seed"]
        independent = bool(seeds[1] < ref_first or seeds[0] > ref_last)
        checks["seeds independent of the reference"] = independent
        print(f"  {'OK  ' if independent else 'FAIL'} seeds {seeds[0]}-{seeds[1]}, reference {ref_first}-{ref_last}", flush=True)
    return checks, summary


def settings_warnings(settings, versions, reference):
    """Warnings for the settings in which the run and the reference differ,
    and a note if their versions differ."""
    warnings = []
    for key, value in settings.items():
        ref_value = reference.get("settings", {}).get(key)
        if ref_value != value:
            warnings.append(f"WARNING: {key} differs from the reference: reference {ref_value}, this run {value}")
    for key, value in versions.items():
        ref_value = reference.get("versions", {}).get(key)
        if ref_value != value:
            warnings.append(f"note: {key} of the reference: {ref_value}; of this run: {value}")
    return warnings


def plot_comparison(result, reference, path, x_pom):
    """dsigma/dt of the run and of the reference, and their ratio."""
    t = np.array(result["t"])
    fig, (ax, ax_ratio) = plt.subplots(2, 1, sharex=True, figsize=(8, 8), gridspec_kw={"height_ratios": [3, 1]})
    for name, color, band in (("coherent", "black", "gray"), ("incoherent", "red", "pink")):
        # mb -> nb
        value = np.array(result[f"{name}_mb_per_GeV2"]) * 1e6
        error = np.array(result[f"{name}_error_mb_per_GeV2"]) * 1e6
        ax.plot(t, value, label=f"Test run {name}", color=color, linestyle="solid")
        ax.fill_between(t, value - error, value + error, color=band, alpha=0.5)
        if reference is None:
            continue
        ref_t = np.array(reference["t"])
        ref_value = np.array(reference[f"{name}_mb_per_GeV2"]) * 1e6
        ref_error = np.array(reference[f"{name}_error_mb_per_GeV2"]) * 1e6
        ax.errorbar(ref_t, ref_value, yerr=ref_error, label=f"Reference {name}", fmt="o", color="blue",
                    alpha=0.8 if name == "coherent" else 0.5)
        ref_interp = np.interp(t, ref_t, ref_value)
        ratio = value / ref_interp
        ratio_error = ratio * np.hypot(error / value, np.interp(t, ref_t, ref_error) / ref_interp)
        ax_ratio.plot(t, ratio, color=color, label=f"{name} ratio")
        ax_ratio.fill_between(t, ratio - ratio_error, ratio + ratio_error, color=band, alpha=0.3)

    ax.set_ylabel(r"$\frac{d\sigma}{dt}$ [$\mathrm{nb}/\mathrm{GeV}^2$]", fontsize=14)
    ax.tick_params(axis="both", labelsize=14)
    ax.set_ylim(bottom=0.1, top=5e2)
    ax.set_yscale("log")
    ax.set_title(r"Test run cross section, $N_\mathrm{events} = " + str(result["n_events"]) + r", x_\mathbb{P} = " + str(x_pom) + "$", fontsize=14)
    ax.legend(fontsize=14)

    ax_ratio.axhline(1.0, color="black", linestyle="--", linewidth=1)
    ax_ratio.tick_params(axis="both", labelsize=14)
    ax_ratio.set_ylabel("Ratio", fontsize=14)
    ax_ratio.set_xlabel(r"$t [\mathrm{GeV}^2]$", fontsize=14)
    ax_ratio.set_ylim(0.5, 1.5)
    ax_ratio.grid(alpha=0.3)
    if reference is not None:
        ax_ratio.legend(loc="best", fontsize=14)

    ax.set_xlim(left=0, right=t[-1])
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print("Created cross section plot:", path, "Events included:", result["n_events"], flush=True)


def evaluate(datadir, args, reference_path, seeds, maxconf, versions, settings):
    """The cross sections of the spectra_1 to spectra_<maxconf> in datadir,
    with jpsi_results.json and the plot, and the check against the
    reference, written to jpsi_check.json, or, with --save-reference, the
    new reference. Returns whether all checks passed."""
    # results of an earlier evaluation must not stand in for these
    for name in ("jpsi_results.json", "cross_section.pdf"):
        if os.path.exists(os.path.join(datadir, name)):
            os.remove(os.path.join(datadir, name))
    numbers, t, amplitudes, invalid, snd_versions = read_spectra(datadir, maxconf)
    versions = dict(versions)
    versions["subnucleondiffraction"] = snd_versions[0] if len(snd_versions) == 1 else (snd_versions or None)
    if len(snd_versions) > 1:
        print("WARNING: the spectra are from several versions of subnucleondiffraction:", snd_versions, flush=True)
    summary = {"settings": vars(args), "versions": versions, "seeds": seeds, "events": len(numbers), "checks": {}}
    checks = summary["checks"]
    for number, reason in invalid:
        print(f"  spectra_{number}: {reason}", flush=True)
    if len(numbers) < 2:
        print("Fewer than two valid spectra: no cross sections", flush=True)
        checks["at least two valid spectra"] = False
    else:
        result = analyze(t, amplitudes, args.check_tmax, args.bootstrap)
        outliers = find_outliers(numbers, amplitudes, args.outlier_zmax)
        for number, z in outliers:
            print(f"  spectra_{number}: outlier, " + ("largest |A|^2 not positive" if z is None else f"z = {z:+.1f}"), flush=True)
        checks["all spectra valid and none an outlier"] = not invalid and not outliers
        print(f"  {'OK  ' if not invalid and not outliers else 'FAIL'} {len(invalid)} invalid spectra, "
              f"{len(outliers)} outliers beyond {args.outlier_zmax:g} robust standard deviations", flush=True)
        summary["invalid_spectra"] = [{"spectra": n, "reason": r} for n, r in invalid]
        summary["outliers"] = [{"spectra": n, "z": z} for n, z in outliers]

        payload = results_payload(result, seeds, versions, settings)
        results_path = os.path.join(datadir, "jpsi_results.json")
        with open(results_path, "w") as f:
            json.dump(payload, f, indent=1)
        print("Results saved to", results_path, flush=True)

        reference = None
        if args.save_reference:
            if seeds is None:
                # a check needs the seeds of the reference, to tell whether
                # a run is independent of it
                print("No reference saved: the seeds are not known (no run_info.json in --datadir); "
                      "give them with --first-seed", flush=True)
                checks["seeds of the reference known"] = False
            elif all(checks.values()):
                with open(reference_path, "w") as f:
                    json.dump(payload, f, indent=1)
                print("Saved the new reference to", reference_path, flush=True)
            else:
                print("No reference saved: there are invalid spectra or outliers", flush=True)
        else:
            reference = load_reference(reference_path)
            if reference is None:
                print(f"No reference found at {reference_path}; make one with --save-reference", flush=True)
                checks["reference found"] = False
            else:
                for warning in settings_warnings(settings, versions, reference):
                    print(warning, flush=True)
                reference_checks, comparison = check_against_reference(result, reference, args.n_sigma, seeds)
                checks.update(reference_checks)
                summary.update(comparison)
        plot_comparison(result, reference, os.path.join(datadir, "cross_section.pdf"), settings["xpom"])
    passed = all(checks.values())
    summary["passed"] = passed
    with open(os.path.join(datadir, "jpsi_check.json"), "w") as f:
        json.dump(summary, f, indent=1)
    return passed


def run_seed(seed, args, repo_root, datadir, subnucleondiffraction_cmd, ipglasma_binary, input_file_path, reference_file_path, ipglasma_path, qs_table_path, xpom, slurm=False, first_seed=0):
    worker_dir = os.path.join(datadir, f"seed_{seed}")
    os.makedirs(worker_dir, exist_ok=True)

    print(f"===== RUNNING SEED {seed} =====", flush=True)
    # IP-Glasma numbers the Wilson-line files of a seed 2 seed + 1 and
    # 2 seed + 2; the spectra in datadir are numbered within the run, 1 to
    # maxevents, for any seeds
    id1 = 2 * seed + 1
    id2 = 2 * seed + 2
    local1 = 2 * (seed - first_seed) + 1
    local2 = 2 * (seed - first_seed) + 2

    temp_input_path = generate_temp_input(source_input_file=input_file_path, seed=seed, x_pom=xpom, working_dir=worker_dir, qs_table_path=qs_table_path)

    log_path = os.path.join(worker_dir, f"ipglasma_log_{seed}")
    ipglasma_cmd = os.path.join(ipglasma_path, ipglasma_binary)
    result = run_command([ipglasma_cmd, temp_input_path], cwd=worker_dir, stdout_path=log_path, stderr_path=log_path + ".err", slurm=slurm)
    print("Removing temp files", temp_input_path, flush=True)
    remove_temp_input(temp_input_path)
    # the Wilson-line files of the events (about 38 MB each) are removed
    # also if the seed fails: its directory is kept for the logs
    wilson_1 = wilson_line_path(worker_dir, xpom, id1)
    wilson_2 = wilson_line_path(worker_dir, xpom, id2)
    try:
        if result.returncode != 0:
            raise RuntimeError(f"ipglasma exited with code {result.returncode} (see {log_path}.err)")
        print("Done", seed, flush=True)

        # Canonicalize Wilson line file names once and reuse paths throughout.
        wilson_1 = ensure_wilson_line_file(worker_dir, xpom, id1)
        wilson_2 = ensure_wilson_line_file(worker_dir, xpom, id2)

        print(f"===== RUNNING SUBNUCLEONDIFFRACTION FOR EVENTS {id1} AND {id2} =====", flush=True)
        spectra_1_path = os.path.join(worker_dir, f"spectra_{id1}")
        spectra_2_path = os.path.join(worker_dir, f"spectra_{id2}")
        for wilson, spectra_path in ((wilson_1, spectra_1_path), (wilson_2, spectra_2_path)):
            result = run_command(subnucleondiffraction_cmd + ["-dipole", "1", "ipglasma_binary", wilson], cwd=worker_dir, stdout_path=spectra_path, stderr_path=spectra_path + ".err", slurm=slurm)
            if result.returncode != 0:
                raise RuntimeError(f"subnucleondiffraction exited with code {result.returncode} (see {spectra_path}.err)")

        for file_id, local_id in ((id1, local1), (id2, local2)):
            src = os.path.join(worker_dir, f"spectra_{file_id}")
            dst = os.path.join(datadir, f"spectra_{local_id}")
            if os.path.exists(src):
                shutil.copy2(src, dst)

        if local1 in plot_ids or local2 in plot_ids:
            plot_id = local1 if local1 in plot_ids else local2
            final_file = wilson_1 if plot_id == local1 else wilson_2
            plot_path = os.path.join(datadir, f"output_{plot_id}.pdf")
            tmp_output_path = os.path.join(worker_dir, f"tmp_output_{plot_id}")
            run_command(
                subnucleondiffraction_cmd + ["-dipole", "1", "ipglasma_binary", final_file, "-print_nucleus", "-Q2", "0"],
                cwd=worker_dir,
                stdout_path=tmp_output_path,
                stderr_path=tmp_output_path + ".err",
                slurm=slurm,
            )
            try:
                data = np.loadtxt(tmp_output_path)
                x = data[:, 0]
                y = data[:, 1]
                z = data[:, 5]

                plt.figure()
                plt.scatter(x, y, c=z, cmap="viridis", s=6, marker="s")
                plt.colorbar(label="col6")
                plt.xlabel(r"$x$")
                plt.ylabel(r"$y$")
                plt.title(fr"{final_file} $\mathrm{{Tr}} V(x)/N_c$")
                plt.tight_layout()
                plt.savefig(plot_path)
                plt.close()
                print("Created plot of V(x) for event", plot_id, "filename:", plot_path, flush=True)
            except Exception as exc:
                print(f"Plotting failed for {final_file}: {exc}", flush=True)
            if os.path.exists(tmp_output_path):
                os.remove(tmp_output_path)
            if os.path.exists(tmp_output_path + ".err"):
                os.remove(tmp_output_path + ".err")

    finally:
        print("Removing Wilson line files ", wilson_1, " and ", wilson_2, flush=True)
        for wilson in (wilson_1, wilson_2):
            if os.path.exists(wilson):
                os.remove(wilson)

    return {"seed": seed, "id1": id1, "id2": id2}


def run_settings(input_file_path, xpom):
    """The settings of a run that its results and the reference record, and
    that a check compares with those of the reference."""
    with open(input_file_path, "rb") as f:
        template_sha256 = hashlib.sha256(f.read()).hexdigest()
    return {"input_template": os.path.basename(input_file_path), "input_template_sha256": template_sha256,
            "xpom": float(xpom), "lattice_N": lattice_N, "lattice_L": lattice_L,
            "subnucleondiffraction_options": " ".join(subnucleondiffraction_options + ["-dipole", "1"])}


def main():
    global subnucleondiffraction_path
    global ipglasma_path

    parser = argparse.ArgumentParser(description="Parallel test for IP-Glasma+JIMWLK+SubnucleonDiffraction for vector meson production.")
    parser.add_argument("--plot-only", action="store_true", help="Skip the simulation: compute the cross sections, the plots and the check (or, with --save-reference, the reference) from the spectra in --datadir.")
    parser.add_argument("--maxevents", type=int, default=200, help="Maximum number of events for the runs.")
    parser.add_argument("--datadir", type=str, default="./jpsi/", help="Directory to store data files.")
    parser.add_argument("--subnucleondiffraction-path", type=str, default=subnucleondiffraction_path, help="Path to the subnucleondiffraction executable.")
    parser.add_argument("--ipglasma-path", type=str, default=ipglasma_path, help="Path to the IP-Glasma executable.")
    parser.add_argument("--input-template", type=str, default="input_vm_proton", help="IP-Glasma input file to use as a template.")
    parser.add_argument("--max-workers", type=int, default=None, help="Seeds run at a time (default: the logical CPUs available divided by OMP_NUM_THREADS), each with OMP_NUM_THREADS threads (default 1); without --slurm, several seeds at a time run with OMP_PROC_BIND=false.")
    parser.add_argument("--keep-logs", action="store_true", help="Keep the working directories of all seeds, with their logs (by default, those of the successful seeds are removed at the end of the run; those of failed seeds are always kept).")
    parser.add_argument("--slurm", action="store_true", help="Wrap each ipglasma/subnucleondiffraction invocation with 'srun --ntasks=1 --cpus-per-task=N', N = OMP_NUM_THREADS (default 1).")
    parser.add_argument("--xpom", type=float, default=default_xpom, help="Value of x_pom to use in the runs.")
    parser.add_argument("--first-seed", type=int, default=None, help=f"First random seed (default {default_first_seed}, beyond the seeds of the reference; 0 with --save-reference). With --plot-only, the seeds recorded in run_info.json.")
    parser.add_argument("--check-tmax", type=float, default=2.0, help="Largest |t| [GeV^2] of the integrated cross sections compared with the reference.")
    parser.add_argument("--n-sigma", type=float, default=3.0, help="Allowed deviation from the reference in units of the error expected if the run and the reference agree.")
    parser.add_argument("--bootstrap", type=int, default=1000, help="Bootstrap resamplings of the events for the errors.")
    parser.add_argument("--outlier-zmax", type=float, default=7.0, help="A spectrum whose largest |A|^2 is more than this many robust standard deviations of log10 |A|^2 from the median is an outlier, and the check fails.")
    parser.add_argument("--reference-file", type=str, default=os.path.join(_script_dir(), reference_file_name), help=f"The reference (default {reference_file_name} next to the script).")
    parser.add_argument("--save-reference", action="store_true", help="Save the results of this run as the new reference (--reference-file) instead of comparing with it; not if a spectrum is invalid or an outlier, or if the seeds are not known (--plot-only without run_info.json or --first-seed).")
    args = parser.parse_args()

    repo_root = _repo_root()
    script_dir = _script_dir()
    datadir = os.path.abspath(args.datadir)
    subnucleondiffraction_path = os.path.abspath(args.subnucleondiffraction_path)
    ipglasma_path = os.path.abspath(args.ipglasma_path)
    maxseed = int(args.maxevents / 2)
    max_workers = args.max_workers or max(1, min(available_cpus() // threads_per_event(), maxseed))
    first_seed = args.first_seed
    if first_seed is None:
        first_seed = 0 if args.save_reference else default_first_seed
    keep_logs = args.keep_logs
    xpom = args.xpom

    # Resolve the input template: try it as given (absolute, or relative to
    # the current working directory) first, then fall back to the script's
    # own directory so the default ("input_vm_proton") keeps working
    # regardless of where the script is invoked from.
    input_file_path = args.input_template
    if not os.path.isfile(input_file_path):
        input_file_path = os.path.join(script_dir, args.input_template)

    reference_file_path = os.path.abspath(args.reference_file)
    run_info_path = os.path.join(datadir, "run_info.json")
    ipglasma_binary = os.path.join(ipglasma_path, "ipglasma")

    if args.plot_only:
        # the seeds, the number of events, the versions and the settings of
        # the run from its run_info.json, so that, as in the run, spectra
        # of an earlier, larger run in datadir do not count
        run_info = {}
        if os.path.exists(run_info_path):
            with open(run_info_path) as f:
                run_info = json.load(f)
        seeds = tuple(run_info["seeds"]) if "seeds" in run_info else None
        maxconf = run_info.get("maxevents")
        if args.first_seed is not None:
            # the largest number, not the count: failed seeds leave gaps
            events = maxconf if maxconf is not None else max(spectra_files(datadir), default=0)
            seeds = (args.first_seed, args.first_seed + (events + 1) // 2 - 1)
        versions = run_info.get("versions")
        if versions is None:
            versions = {"ipglasma": ipglasma_version(ipglasma_binary)}
            print(f"No versions in {run_info_path}: IP-Glasma from {ipglasma_binary}: {versions['ipglasma']}", flush=True)
        settings = run_info.get("settings") or run_settings(input_file_path, xpom)
        passed = evaluate(datadir, args, reference_file_path, seeds, maxconf, versions, settings)
        if not passed:
            print("FAILED: see", os.path.join(datadir, "jpsi_check.json"), flush=True)
        raise SystemExit(0 if passed else 1)

    print("PID =", os.getpid())
    if hasattr(os, "sched_getaffinity"):
        print("Affinity =", os.sched_getaffinity(0))

    # Read the Qs table file name from the input template (NucleusQsTableFileName)
    # instead of hardcoding it, so the template controls which table is used.
    qs_table_name = read_input_value(input_file_path, "nucleusQsTableFileName") if os.path.isfile(input_file_path) else None
    if not qs_table_name:
        raise FileNotFoundError(f"nucleusQsTableFileName not found in input template '{input_file_path}'.")
    qs_table_path = os.path.join(ipglasma_path, qs_table_name)

    os.makedirs(datadir, exist_ok=True)
    shutil.copy2(input_file_path, os.path.join(datadir, "input_vm_proton"))
    seeds = (first_seed, first_seed + maxseed - 1)
    # the version of subnucleondiffraction is in the spectra
    versions = {"ipglasma": ipglasma_version(ipglasma_binary)}
    settings = run_settings(input_file_path, xpom)
    with open(run_info_path, "w") as f:
        json.dump({"seeds": seeds, "maxevents": 2 * maxseed, "versions": versions, "settings": settings}, f, indent=1)
    # the spectra and results of an earlier run in datadir must not stand in
    # for those of this one if its events fail, or do not run because it is
    # interrupted (--plot-only); spectra beyond maxevents do not count anyway
    stale_files = [f"spectra_{n}" for n in range(1, 2 * maxseed + 1)] + [f"output_{n}.pdf" for n in plot_ids] + [
        "jpsi_results.json", "cross_section.pdf", "jpsi_check.json", "failed_seeds.txt",
        # written by earlier versions of the script
        "cross_section.dat"]
    for name in stale_files:
        stale = os.path.join(datadir, name)
        if os.path.exists(stale):
            os.remove(stale)

    subnucleondiffraction_binary = os.path.join(subnucleondiffraction_path, "build", "bin", "subnucleondiffraction")
    subnucleondiffraction_wavefile = os.path.join(subnucleondiffraction_path, "gauss-boosted_mzsat.dat")

    missing = []
    non_executable = []
    for required,exec in zip([ipglasma_binary, subnucleondiffraction_binary, subnucleondiffraction_wavefile, input_file_path, qs_table_path], [True, True, False, False, False]):
        if not os.path.exists(required):
            missing.append(required)
        elif os.path.isfile(required) and not os.access(required, os.X_OK) and exec:
            non_executable.append(required)

    if missing or non_executable:
        details = []
        if missing:
            details.append("Missing required files:\n" + "\n".join(f"  - {p}" for p in missing))
        if non_executable:
            details.append("Found but not executable:\n" + "\n".join(f"  - {p}" for p in non_executable))
        raise FileNotFoundError("Required binaries/files are not available before launching workers:\n" + "\n".join(details))

    subnucleondiffraction_cmd = [subnucleondiffraction_binary, "-wavef_file", subnucleondiffraction_wavefile] + subnucleondiffraction_options

    failures = {}

    if max_workers > 1 and not args.slurm:
        # binding (OMP_PROC_BIND or OMP_PLACES, as in the main README) would
        # put the threads of all events on the same first cores; with
        # --slurm, srun gives each event its own CPUs
        os.environ["OMP_PROC_BIND"] = "false"
    print(f"Running {maxseed} seeds ({seeds[0]} to {seeds[1]}) in parallel with {max_workers} workers")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                run_seed,
                seed,
                args,
                repo_root,
                datadir,
                subnucleondiffraction_cmd,
                ipglasma_binary,
                input_file_path,
                reference_file_path,
                ipglasma_path,
                qs_table_path,
                xpom,
                args.slurm,
                first_seed,
            ): seed
            for seed in range(first_seed, first_seed + maxseed)
        }
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                seed = futures[future]
                print(f"seed {seed} failed: {exc}", flush=True)
                failures[seed] = str(exc)

    # spectra of an earlier, larger run in datadir do not count
    passed = evaluate(datadir, args, reference_file_path, seeds, 2 * maxseed, versions, settings)

    if failures:
        failure_file = os.path.join(datadir, "failed_seeds.txt")
        with open(failure_file, "w") as f:
            for seed in sorted(failures):
                f.write(f"seed {seed}: {failures[seed]}\n")
        print(f"{len(failures)} of {maxseed} seed(s) failed, without spectra; details: {failure_file}", flush=True)

    if not keep_logs:
        # the working directories of failed seeds stay, with their logs
        for seed in range(first_seed, first_seed + maxseed):
            worker_dir = os.path.join(datadir, f"seed_{seed}")
            if seed not in failures and os.path.exists(worker_dir):
                shutil.rmtree(worker_dir)
        print(f"Removed temporary directories in {datadir}", flush=True)

    if failures:
        print(f"FAILED: {len(failures)} seed(s) failed", flush=True)
    if not passed:
        print("FAILED: see", os.path.join(datadir, "jpsi_check.json"), flush=True)
    elif not failures:
        print("PASSED", flush=True)
    raise SystemExit(0 if passed and not failures else 1)

if __name__ == "__main__":
    main()

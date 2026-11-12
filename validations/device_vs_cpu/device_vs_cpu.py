#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
The device layer (GPU) of IP-Glasma against the CPU code, for a whole
event.

The script runs one Pb+Pb event of the example input (hot spots, JIMWLK
from x = 0.01 to 1.29e-4, the classical Yang-Mills evolution to
tau = 0.4 fm) at a fixed seed and b = 0

  * once on the CPU (the reference), and
  * once per requested device and precision (IPGLASMA_DEVICE,
    IPGLASMA_DEVICE_PRECISION),

and compares what the event writes:

  * the Wilson lines of both nuclei after the JIMWLK evolution (binary,
    all digits),
  * the hydrodynamic initial condition (epsilon, u^mu, pi^munu), in three
    bands of epsilon / max(epsilon),
  * dN/dy, dE/dy and <k_T> of the gluons, and the eccentricities.

The Gaussian noise of the JIMWLK evolution is drawn on the host from the
random numbers of the run, also in a device run, so both runs see the same
noise and are compared site by site, not statistically.

What to expect:

  * double precision (cuda; opencl on devices with it; host): the device
    repeats the CPU arithmetic up to rounding. The Wilson lines agree to
    about 1e-12, the integrated observables to 1e-9 or better, and the
    text output of the hydro fields in all printed digits but the last.
  * single precision (the GPUs of Apple M processors; any device with
    IPGLASMA_DEVICE_PRECISION=single): the fields carry the rounding of
    single precision, about 1e-7 per operation, which the 300 JIMWLK steps
    accumulate. The Wilson lines agree to a few 1e-6 (and are unitary to
    1e-15: the host projects what a single-precision device returns onto
    SU(3)), the integrated observables to 1e-5 - 1e-4, and epsilon to
    1e-4 or better down to 0.1% of its maximum (measured at 256^2 and
    512^2).

The tolerances below are these expectations with a margin; the script
prints PASSED or FAILED per run and writes all numbers to
<datadir>/device_vs_cpu_summary.json.

Usage (from validations/device_vs_cpu/):
    # the host backend, always built: 128^2, both precisions (1-2 minutes)
    python3 device_vs_cpu.py --ipglasma-path ../.. --quick
    # on a machine with an NVIDIA GPU (build with -DIPGLASMA_CUDA=ON)
    python3 device_vs_cpu.py --ipglasma-path ../.. --devices cuda
    # on an Apple M machine (build with -DIPGLASMA_OPENCL=ON)
    python3 device_vs_cpu.py --ipglasma-path ../.. --devices opencl
    # several runs: device[:precision], precision = single or double
    python3 device_vs_cpu.py --ipglasma-path ../.. \\
        --devices cuda:double cuda:single host:double --size 512

--ipglasma-path is a directory with the ipglasma executable,
qs2Adj_vs_Tp_vs_Y_240.in, nucleusConfigurations/ and tables/, and the
example input file `input`.
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "utilities"))
from read_wilson_lines import read_wilson_lines  # noqa: E402

# the bands of epsilon / max(epsilon) of the comparison of the hydro fields
BANDS = [("bulk", 0.1, np.inf), ("mid", 0.01, 0.1), ("edge", 0.001, 0.01)]

# tolerances: (double precision, single precision)
TOLERANCES = {
    "wilson_lines_max": (1e-10, 5e-5),
    "wilson_lines_unitarity": (1e-10, 1e-10),
    "dN": (1e-8, 5e-4),
    "dE": (1e-8, 5e-4),
    "mean_kT": (1e-8, 5e-4),
    "eccentricity_2_3": (1e-5, 1e-3),
    "epsilon_total": (1e-5, 1e-4),
    "epsilon_bulk_max_rel": (1e-5, 1e-3),
    "epsilon_mid_max_rel": (1e-4, 2e-3),
    "epsilon_edge_max_rel": (1e-3, 1e-2),
    "u_bulk_max": (2e-5, 1e-3),
}


def write_input(template, path, size, seed, steps_per_spacing):
    """The example input with the lattice size, the seed and the outputs
    of this comparison."""
    changes = {
        "size": str(size),
        "seed": str(seed),
        "bMin": "0", "bMax": "0",
        "LOutput": "30", "sizeOutput": str(size // 2),
        "timeStepsPerLatticeSpacing": str(steps_per_spacing),
        "computeEccentricities": "1",
        "computeGluonMultiplicity": "1",
        "writeWilsonLines": "2",
        "writeWilsonLineGeometry": "0",
        "writeHydro": "1",
    }
    lines = []
    seen = set()
    with open(template) as f:
        for line in f:
            key = line.split()[0] if line.split() else ""
            if key in changes:
                lines.append(f"{key} {changes[key]}\n")
                seen.add(key)
            elif key == "EndOfFile":
                break
            else:
                lines.append(line)
    for key, value in changes.items():
        if key not in seen:
            lines.append(f"{key} {value}\n")
    # a cutoff well above the rounding of the vacuum: with the default 0,
    # the columns of the eccentricity file that count the cells above the
    # cutoff depend on the sign of the rounding there
    lines.append("eccentricityCutoff 0.1\n")
    lines.append("EndOfFile\n")
    with open(path, "w") as f:
        f.writelines(lines)


def run_event(ipglasma_path, rundir, input_file, environment, threads):
    """Runs one event in rundir; returns the wall time."""
    if os.path.isdir(rundir):
        shutil.rmtree(rundir)
    os.makedirs(rundir)
    for name in ("qs2Adj_vs_Tp_vs_Y_240.in", "nucleusConfigurations",
                 "tables"):
        os.symlink(os.path.abspath(os.path.join(ipglasma_path, name)),
                   os.path.join(rundir, name))
    shutil.copy(input_file, os.path.join(rundir, "input"))
    env = dict(os.environ)
    for key in ("IPGLASMA_DEVICE", "IPGLASMA_DEVICE_PRECISION",
                "IPGLASMA_DEVICE_STAGES"):
        env.pop(key, None)
    env.update(environment)
    env["OMP_NUM_THREADS"] = str(threads)
    env["IPGLASMA_PROFILE"] = "1"
    executable = os.path.abspath(os.path.join(ipglasma_path, "ipglasma"))
    start = time.time()
    with open(os.path.join(rundir, "log.txt"), "w") as log:
        status = subprocess.call([executable, "input"], cwd=rundir, env=env,
                                 stdout=log, stderr=subprocess.STDOUT)
    elapsed = time.time() - start
    if status != 0:
        with open(os.path.join(rundir, "log.txt")) as log:
            tail = "".join(log.readlines()[-15:])
        raise RuntimeError(f"ipglasma failed in {rundir} (exit status "
                           f"{status}):\n{tail}")
    return elapsed


def read_profile(rundir):
    """{phase: seconds} of the event."""
    phases = {}
    path = os.path.join(rundir, "ipglasma_profile_rank0.tsv")
    if not os.path.exists(path):
        return phases
    with open(path) as f:
        header = [name.strip() for name in f.readline().split("\t")]
        for line in f:
            fields = line.rstrip("\n").split("\t")
            phases[fields[header.index("phase")]] = float(
                fields[header.index("seconds")])
    return phases


def read_event(rundir):
    """The outputs of the event that the comparison uses."""
    event = {}
    wilson = sorted(glob.glob(os.path.join(rundir, "WilsonLine_x_*")))
    wilson = [w for w in wilson if not w.endswith(".txt")]
    event["wilson"] = {os.path.basename(w): read_wilson_lines(w)[0]
                       for w in wilson}
    hydro = glob.glob(os.path.join(rundir, "epsilon-u-Hydro-TauHydro-*.dat"))
    event["hydro"] = np.loadtxt(hydro[0])
    with open(os.path.join(rundir, "gluonMultiplicity0.json")) as f:
        event["gluons"] = json.load(f)
    event["eccentricities"] = np.loadtxt(
        os.path.join(rundir, "eccentricities0.dat"))
    event["profile"] = read_profile(rundir)
    with open(os.path.join(rundir, "log.txt")) as f:
        event["device_line"] = next(
            (line.strip() for line in f if "[main::main]: Device " in line),
            "")
    return event


def compare(reference, event):
    """The differences of a device run from the CPU run."""
    result = {}

    # Wilson lines after JIMWLK: the files of the final x of each nucleus
    # are written last; compare every file both runs wrote
    differences = []
    unitarity = []
    for name, lines in reference["wilson"].items():
        other = event["wilson"][name]
        differences.append(float(np.max(np.abs(other - lines))))
        product = np.einsum("xyab,xycb->xyac", other, other.conj())
        unitarity.append(float(np.max(np.abs(product - np.eye(3)))))
    result["wilson_lines_max"] = max(differences)
    result["wilson_lines_unitarity"] = max(unitarity)

    # integrated observables
    for key, name in (("dN", "dN"), ("dE", "dE_GeV"),
                      ("mean_kT", "mean_kT_GeV")):
        a = reference["gluons"][name]
        b = event["gluons"][name]
        result[key] = abs(b / a - 1.)
    # columns: tau, (eps_n, Psi_n) for n = 1..6, ...
    ecc_ref = reference["eccentricities"]
    ecc = event["eccentricities"]
    result["eccentricity_2_3"] = float(max(
        abs(ecc[3] / ecc_ref[3] - 1.), abs(ecc[5] / ecc_ref[5] - 1.)))
    result["eccentricity_4_5_6"] = float(max(
        abs(ecc[i] / ecc_ref[i] - 1.) for i in (7, 9, 11)))

    # hydro fields: columns eta, x, y, epsilon, u^tau, u^x, u^y, u^eta,
    # pi^munu (10)
    hydro_ref = reference["hydro"]
    hydro = event["hydro"]
    eps_ref = hydro_ref[:, 3]
    eps = hydro[:, 3]
    result["epsilon_total"] = float(abs(eps.sum() / eps_ref.sum() - 1.))
    result["epsilon_rel_L2"] = float(
        np.linalg.norm(eps - eps_ref) / np.linalg.norm(eps_ref))
    peak = eps_ref.max()
    for name, low, high in BANDS:
        mask = (eps_ref > low * peak) & (eps_ref <= high * peak)
        if not mask.any():
            continue
        relative = np.abs(eps[mask] - eps_ref[mask]) / eps_ref[mask]
        result[f"epsilon_{name}_max_rel"] = float(relative.max())
        result[f"epsilon_{name}_rms_rel"] = float(
            np.sqrt(np.mean(relative**2)))
        result[f"u_{name}_max"] = float(
            np.abs(hydro[mask, 4:8] - hydro_ref[mask, 4:8]).max())
        result[f"cells_{name}"] = int(mask.sum())
    return result


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--ipglasma-path", required=True,
                        help="directory with the ipglasma executable, the "
                        "Qs table, nucleusConfigurations/, tables/ and "
                        "the example input")
    parser.add_argument("--devices", nargs="+", default=None,
                        metavar="DEVICE[:PRECISION]",
                        help="runs to compare with the CPU: cuda, opencl "
                        "or host, optionally with :single or :double "
                        "(default precision: that of the device; "
                        "default runs: host:double host:single)")
    parser.add_argument("--size", type=int, default=256,
                        help="lattice size (default 256; even, with prime "
                        "factors up to 32)")
    parser.add_argument("--seed", type=int, default=4242)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1,
                        help="OpenMP threads of every run")
    parser.add_argument("--steps-per-spacing", type=int, default=10,
                        help="timeStepsPerLatticeSpacing (default 10)")
    parser.add_argument("--datadir", default="device_vs_cpu_output")
    parser.add_argument("--quick", action="store_true",
                        help="a quick test: size 128")
    parser.add_argument("--no-check", action="store_true",
                        help="report the differences without applying the "
                        "tolerances")
    args = parser.parse_args()
    if args.quick:
        args.size = 128
    if args.devices is None:
        args.devices = ["host:double", "host:single"]

    datadir = os.path.abspath(args.datadir)
    os.makedirs(datadir, exist_ok=True)
    input_file = os.path.join(datadir, "input")
    write_input(os.path.join(args.ipglasma_path, "input"), input_file,
                args.size, args.seed, args.steps_per_spacing)

    print(f"Pb+Pb, b = 0, size {args.size}, seed {args.seed}, "
          f"{args.threads} threads")
    print("running on the CPU ...", flush=True)
    wall = {"cpu": run_event(args.ipglasma_path, os.path.join(datadir, "cpu"),
                             input_file, {}, args.threads)}
    reference = read_event(os.path.join(datadir, "cpu"))

    results = {}
    profiles = {"cpu": reference["profile"]}
    passed = True
    for spec in args.devices:
        device, _, precision = spec.partition(":")
        environment = {"IPGLASMA_DEVICE": device}
        if precision:
            environment["IPGLASMA_DEVICE_PRECISION"] = precision
        label = spec.replace(":", "-")
        print(f"running on {spec} ...", flush=True)
        try:
            wall[label] = run_event(
                args.ipglasma_path, os.path.join(datadir, label), input_file,
                environment, args.threads)
        except RuntimeError as error:
            print(f"  FAILED to run: {error}")
            results[label] = {"error": str(error), "passed": False}
            passed = False
            continue
        event = read_event(os.path.join(datadir, label))
        profiles[label] = event["profile"]
        result = compare(reference, event)
        single = "single precision" in event["device_line"]
        result["device"] = event["device_line"].split("[main::main]: ")[-1]
        result["single_precision"] = single
        failures = []
        for key, (tolerance_double, tolerance_single) in TOLERANCES.items():
            if key not in result:
                continue
            tolerance = tolerance_single if single else tolerance_double
            if not result[key] <= tolerance:
                failures.append(f"{key} = {result[key]:.3g} > {tolerance:g}")
        result["failures"] = failures
        result["passed"] = not failures
        results[label] = result
        if failures and not args.no_check:
            passed = False

    # report
    print()
    for label, result in results.items():
        if "error" in result:
            continue
        print(f"== {label}: {result['device']}")
        rows = [
            ("Wilson lines after JIMWLK, max |dV|", "wilson_lines_max"),
            ("  their unitarity, max |V V^+ - 1|", "wilson_lines_unitarity"),
            ("dN/dy, relative", "dN"),
            ("dE/dy, relative", "dE"),
            ("<k_T>, relative", "mean_kT"),
            ("eccentricities eps_2, eps_3, relative", "eccentricity_2_3"),
            ("eccentricities eps_4..6, relative", "eccentricity_4_5_6"),
            ("sum of epsilon, relative", "epsilon_total"),
            ("epsilon, relative L2", "epsilon_rel_L2"),
        ]
        for name, low, high in BANDS:
            band = (f"eps > {low:g} max" if not np.isfinite(high)
                    else f"{low:g} < eps/max < {high:g}")
            rows.append((f"epsilon, max rel. ({band})",
                         f"epsilon_{name}_max_rel"))
            rows.append((f"epsilon, rms rel. ({band})",
                         f"epsilon_{name}_rms_rel"))
            rows.append((f"u^mu, max |du| ({band})", f"u_{name}_max"))
        for text, key in rows:
            if key not in result:
                continue
            tolerance = ""
            if key in TOLERANCES:
                limit = TOLERANCES[key][1 if result["single_precision"] else 0]
                tolerance = f"   (tolerance {limit:g})"
            print(f"  {text:<46s} {result[key]:10.3e}{tolerance}")
        if result["failures"]:
            print("  " + ("differences above the tolerances: "
                          if args.no_check else "FAILED: ")
                  + "; ".join(result["failures"]))
        else:
            print("  PASSED")
        print()

    stages = [("jimwlk.evolution", "JIMWLK"),
              ("evolution.total", "Yang-Mills + observables"),
              ("event.total", "event")]
    print("time [s]" + "".join(f"{name:>28s}" for _, name in stages))
    for label, profile in profiles.items():
        print(f"{label:<8s}" + "".join(
            f"{profile.get(phase, float('nan')):28.2f}"
            for phase, _ in stages))

    summary = {"settings": vars(args), "results": results,
               "profiles": profiles, "wall_seconds": wall}
    with open(os.path.join(datadir, "device_vs_cpu_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{datadir}/device_vs_cpu_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

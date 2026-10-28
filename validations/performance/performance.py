#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Performance of IP-Glasma: the time per stage of an event and the peak
memory, against the lattice size and the number of OpenMP threads.

The script runs the example input (Pb+Pb with hot spots, JIMWLK from
x = 0.01 to 1e-4, the classical Yang-Mills evolution to tau = 0.5 fm, the
gluon multiplicity and the eccentricities) at a fixed geometry (--seed,
b = 0), without writing Wilson lines, JIMWLK snapshots or fields, with
IPGLASMA_PROFILE=1. Two scans:

  * the size N at --threads (by default N = 256 to 1536 at L = 30 fm with
    8 threads; the lattice spacing shrinks with N, so the Yang-Mills
    evolution, whose time step is a/10, needs N times more steps),
  * the number of OpenMP threads at --scan-size (by default 1 to 16 at
    N = 512).

The runs are sequential, and each is repeated --repeats times; the fastest
repeat is kept, since a shared machine can only make a run slower. From
the profile of the event (ipglasma_profile_rank0.tsv, see OUTPUT.md), the
script takes the time and the number of calls of every phase. The
figure shows the main stages (the Wilson lines, the JIMWLK evolution, the
light-cone fields, the Yang-Mills evolution with the observables at the
final time) and the whole event, and the time of one step of the JIMWLK
evolution (jimwlk.step) and of the Yang-Mills evolution
(evolution.parallel_step): the number of JIMWLK steps does not depend on
N, that of the Yang-Mills steps grows like N. The peak resident memory of
each run comes from the operating system (wait4).

Timings depend on the machine and on what else runs on it: run the script
on an otherwise idle machine. The summary records the CPU model, the
number of CPUs and the load average before the runs, and the script warns
if the load is above 1.

Usage (from validations/performance/):
    # quick test: N = 64 and 128, 1 and 2 threads, one repeat
    python3 performance.py --ipglasma-path ../.. --quick
    # the setup of the paper (about 1.5 hours on 16 cores)
    python3 performance.py --ipglasma-path ../..

--ipglasma-path is a directory with the ipglasma executable,
qs2Adj_vs_Tp_vs_Y_240.in, nucleusConfigurations/ and tables/. The results
go to <datadir>/performance_summary.json and performance.pdf.
"""

import argparse
import json
import os
import platform
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "numerical_parameters"))
from numerical_parameters import write_input  # noqa: E402

# the stages shown in the figure: phases of the profile that do not contain
# each other, and the whole event
STAGES = [("initialization.wilson_lines", "Wilson lines"),
          ("jimwlk.evolution", "JIMWLK"),
          ("initialization.forward_light_cone", "light-cone fields"),
          ("evolution.total", "Yang-Mills evolution and observables"),
          ("event.total", "event")]
# the phases entered once per step: the time per step is seconds/calls
STEPS = [("jimwlk.step", "JIMWLK step"),
         ("evolution.parallel_step", "Yang-Mills step")]


def read_profile(path, event=0):
    """{phase: seconds} and {phase: calls} of one event of a profile
    file."""
    phases, calls = {}, {}
    with open(path) as f:
        header = f.readline().split("\t")
        index = {name.strip(): i for i, name in enumerate(header)}
        for line in f:
            cells = line.rstrip("\n").split("\t")
            if int(cells[index["event"]]) == event:
                phase = cells[index["phase"]]
                phases[phase] = float(cells[index["seconds"]])
                calls[phase] = int(cells[index["calls"]])
    return phases, calls


def per_step(result, phase):
    """The time of one call of a phase, NaN without it."""
    calls = result.get("calls", {}).get(phase, 0)
    return result["phases"][phase] / calls if calls else float("nan")


def run_tag(cfg):
    return f"size{cfg['size']}_threads{cfg['threads']}_rep{cfg['repeat']}"


def run_one(cfg, settings):
    """Run one event; its wall time, peak memory and the times of the
    phases of its profile (cached in <workdir>/result.json)."""
    workdir = os.path.join(settings["datadir"], run_tag(cfg))
    cache = os.path.join(workdir, "result.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return json.load(f), None
    os.makedirs(workdir, exist_ok=True)
    profile = os.path.join(workdir, "ipglasma_profile_rank0.tsv")
    if os.path.exists(profile):
        os.remove(profile)
    tables = os.path.join(workdir, "tables")
    if os.path.isdir(settings["tables"]) and not os.path.exists(tables):
        os.symlink(settings["tables"], tables)
    write_input(settings["template"],
                dict(settings["overrides"], size=cfg["size"],
                     sizeOutput=cfg["size"]),
                os.path.join(workdir, "input"))
    env = dict(os.environ, OMP_NUM_THREADS=str(cfg["threads"]),
               IPGLASMA_PROFILE="1", IPGLASMA_PROFILE_DIR=workdir)
    with open(os.path.join(workdir, "run.log"), "w") as log:
        start = time.perf_counter()
        process = subprocess.Popen([settings["ipglasma"], "input"],
                                   cwd=workdir, stdout=log,
                                   stderr=subprocess.STDOUT, env=env)
        _, status, usage = os.wait4(process.pid, 0)
        wall = time.perf_counter() - start
    returncode = os.waitstatus_to_exitcode(status)
    if returncode != 0 or not os.path.exists(profile):
        return None, f"ipglasma failed (exit {returncode}), see {workdir}"
    phases, calls = read_profile(profile)
    result = dict(cfg, wall_s=wall,
                  # ru_maxrss is in kilobytes on Linux
                  max_rss_MiB=usage.ru_maxrss / 1024.,
                  phases=phases, calls=calls)
    with open(cache, "w") as f:
        json.dump(result, f, indent=1)
    for name in os.listdir(workdir):
        if name not in ("result.json", "input", "run.log", "tables",
                        "ipglasma_profile_rank0.tsv"):
            path = os.path.join(workdir, name)
            if os.path.isfile(path):
                os.remove(path)
    return result, None


def fastest(results, size, threads):
    runs = [r for r in results
            if r["size"] == size and r["threads"] == threads]
    return min(runs, key=lambda r: r["phases"].get("event.total", r["wall_s"])) \
        if runs else None


def machine():
    model = platform.processor()
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    model = line.split(":", 1)[1].strip()
                    break
    except OSError:
        pass
    return {"cpu": model, "cpus": os.cpu_count(),
            "load_average": os.getloadavg(),
            "system": platform.platform()}


def plot(summary, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#1f1f1e"]
    fig, ((left, steps), (middle, right)) = plt.subplots(2, 2,
                                                          figsize=(10, 7))
    sizes = summary["size_scan"]
    n = np.array([r["size"] for r in sizes])
    for (phase, label), color in zip(STAGES, colors):
        t = np.array([r["phases"].get(phase, np.nan) for r in sizes])
        left.plot(n, t, "o-", color=color, ms=3, label=label)
    left.set_xscale("log")
    left.set_yscale("log")
    left.set_xlabel("N")
    left.set_ylabel("time per event [s]")
    left.set_title(f"{summary['settings']['threads']} threads", fontsize=9)
    left.legend(fontsize=7)
    for (phase, label), color in zip(STEPS, colors[1::2]):
        steps.plot(n, [per_step(r, phase) for r in sizes], "o-",
                   color=color, ms=3, label=label)
    steps.set_xscale("log")
    steps.set_yscale("log")
    steps.set_xlabel("N")
    steps.set_ylabel("time per step [s]")
    steps.legend(fontsize=7)
    threads = summary["thread_scan"]
    p = np.array([r["threads"] for r in threads])
    for (phase, label), color in zip(STAGES, colors):
        t = np.array([r["phases"].get(phase, np.nan) for r in threads])
        middle.plot(p, t[0] / t * p[0], "o-", color=color, ms=3,
                    label=label)
    middle.plot(p, p / p[0], ":", color="#8a8a85", label="ideal")
    middle.set_xscale("log", base=2)
    middle.set_yscale("log", base=2)
    middle.set_xlabel("OpenMP threads")
    middle.set_ylabel(f"speed-up against {p[0]} thread(s)")
    middle.set_title(f"N = {summary['settings']['scan_size']}", fontsize=9)
    middle.legend(fontsize=7)
    right.plot(n, [r["max_rss_MiB"] / 1024. for r in sizes], "o-",
               color=colors[0], ms=3)
    right.set_xscale("log")
    right.set_yscale("log")
    right.set_xlabel("N")
    right.set_ylabel("peak memory [GiB]")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Time per stage and memory of IP-Glasma against the "
                    "lattice size and the number of threads.")
    parser.add_argument("--ipglasma-path", default="../..")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--template", default=os.path.join(HERE, "..", "..",
                                                           "input"))
    parser.add_argument("--datadir", default="./performance_output")
    parser.add_argument("--L", type=float, default=30., help="[fm]")
    parser.add_argument("--seed", type=int, default=3)
    parser.add_argument("--sizes", type=int, nargs="+",
                        default=[256, 512, 768, 1024, 1536])
    parser.add_argument("--threads", type=int, default=8,
                        help="the threads of the size scan")
    parser.add_argument("--thread-counts", type=int, nargs="+",
                        default=[1, 2, 4, 8, 16])
    parser.add_argument("--scan-size", type=int, default=512,
                        help="the size of the thread scan")
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--quick", action="store_true",
                        help="N = 64 and 128, 1 and 2 threads, one repeat")
    args = parser.parse_args()
    if args.quick:
        args.sizes, args.threads = [64, 128], 2
        args.thread_counts, args.scan_size, args.repeats = [1, 2], 128, 1

    path = os.path.abspath(args.ipglasma_path)
    settings = {
        "ipglasma": os.path.join(path, args.ipglasma_cmd),
        "template": os.path.abspath(args.template),
        "datadir": os.path.abspath(args.datadir),
        "tables": os.path.join(path, "tables"),
        "overrides": {
            "L": args.L, "LOutput": args.L, "seed": args.seed, "bMin": 0.,
            "bMax": 0., "rotateReactionPlane": 0, "writeWilsonLines": 0,
            "jimwlkSaveSnapshots": 0, "writeTmunu": 0, "writeHydro": 0,
            "computeGluonMultiplicity": 1, "computeEccentricities": 1,
            "nucleusQsTableFileName": os.path.join(
                path, "qs2Adj_vs_Tp_vs_Y_240.in"),
            "nuclearConfigurationsPath": os.path.join(
                path, "nucleusConfigurations"),
        },
    }
    os.makedirs(settings["datadir"], exist_ok=True)
    info = machine()
    print(f"{info['cpu']}, {info['cpus']} CPUs, load average "
          + ", ".join(f"{v:.2f}" for v in info["load_average"]))
    if info["load_average"][0] > 1.:
        print("WARNING: the machine is not idle; the timings will be too "
              "long.")

    configs = [{"size": n, "threads": args.threads} for n in args.sizes] \
        + [{"size": args.scan_size, "threads": t}
           for t in args.thread_counts if t != args.threads
           or args.scan_size not in args.sizes]
    runs, failures = [], []
    start = time.time()
    for repeat in range(args.repeats):
        for cfg in configs:
            result, error = run_one(dict(cfg, repeat=repeat), settings)
            if error:
                failures.append(f"{run_tag(dict(cfg, repeat=repeat))}: "
                                f"{error}")
                continue
            runs.append(result)
            print(f"  N = {cfg['size']:5d}, {cfg['threads']:2d} threads, "
                  f"repeat {repeat}: event "
                  f"{result['phases'].get('event.total', float('nan')):8.1f} s,"
                  f" peak memory {result['max_rss_MiB'] / 1024.:.2f} GiB",
                  flush=True)
    print(f"done in {time.time() - start:.0f} s")
    for failure in failures:
        print("FAILED:", failure)

    size_scan = [fastest(runs, n, args.threads) for n in args.sizes]
    thread_scan = [fastest(runs, args.scan_size, t)
                   for t in args.thread_counts]
    summary = {"settings": vars(args), "machine": info, "runs": runs,
               "size_scan": [r for r in size_scan if r],
               "thread_scan": [r for r in thread_scan if r],
               "checks": {"all runs succeeded": not failures,
                          "every event has a profile": all(
                              "event.total" in r["phases"] for r in runs)}}
    print(f"\nfastest of {args.repeats} repeat(s) [s]:")
    print(f"  {'N':>5s} {'threads':>7s}" + "".join(
        f"{label[:18]:>20s}" for _, label in STAGES + STEPS)
        + f"{'memory [GiB]':>14s}")
    for r in summary["size_scan"] + summary["thread_scan"]:
        print(f"  {r['size']:5d} {r['threads']:7d}" + "".join(
            f"{r['phases'].get(phase, float('nan')):20.2f}"
            for phase, _ in STAGES) + "".join(
            f"{per_step(r, phase):20.4f}" for phase, _ in STEPS)
            + f"{r['max_rss_MiB'] / 1024.:14.2f}")
    passed = all(summary["checks"].values())
    summary["passed"] = passed
    for name, ok in summary["checks"].items():
        print(f"{'OK  ' if ok else 'FAIL'} {name}")
    with open(os.path.join(settings["datadir"], "performance_summary.json"),
              "w") as f:
        json.dump(summary, f, indent=1)
    if summary["size_scan"] and summary["thread_scan"]:
        plot(summary, os.path.join(settings["datadir"], "performance.pdf"))
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{settings['datadir']}/performance_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

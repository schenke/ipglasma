#!/usr/bin/env python3
"""
Parallel version of the vector meson production validation script.

This script parallelizes the loop over seeds so that multiple IP-Glasma runs execute
simultaneously. Each worker gets its own temporary input file, so the simulations are
independent and do not share configuration state.

This tool is used to test the IP-Glasma+JIMWLK+SubnucleonDiffraction framework for vector meson production. 
It runs a series of simulations with different random seeds, collects the resulting cross section data, 
and compares it to reference data. The reference data has been calculated using a large number of events and is stored in the file "validation_spectra",
using a verified version of the code. 

The computed cross section should also match the published results in https://arxiv.org/pdf/2207.03712 Fig. 1 ("CGC+shape fluct") curve,
except that by default the code uses a smaller lattice for performance reasons, which may lead to small deviations in the cross section. 
The lattice size can be increased by modifying the input file "input_vm_proton".

The script also generates plots  color field (tr V(x)) for selected events.

Example usage:
    python3 test_vector_meson_production.py -maxevents 100 -datadir path_for_cross_sections -subnucleondiffraction_path path_to_subnucleondiffraction_directory

Maximum number of parallel workers can be specified with -max_workers. By default, it uses the number of logical CPUs.

The comparison plot will be saved in the specified directory as "cross_section.pdf". The color field plots will be saved as "output_<event_id>.pdf".

Use -plot_only to skip the simulation and only generate the comparison plots from existing data files in the specified directory.

Before running this, make sure you have downloaded and built the subnucleondiffraction code from https://github.com/hejajama/subnucleondiffraction
'''

"""

import argparse
import cmd
import cmd
import fnmatch
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

import matplotlib.pyplot as plt
import numpy as np


subnucleondiffraction_path = "./subnucleondiffraction/"
ipglasma_path = "../"
ipglasma_cmd = "ipglasma"
plot_ids = [1, 5]
xpom = "0.001705"
reference_cross_section_file = "validation_spectra"
lattice_N = 512 # N*N lattice
lattice_L = 4 # fm


def _repo_root():
    return os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def _script_dir():
    return os.path.dirname(os.path.abspath(__file__))


def generate_temp_input(source_input_file="input_vm_proton", seed=0, x_pom="0.001705", working_dir=".", ipglasma_path=ipglasma_path):
    with open(source_input_file, "r") as f:
        lines = f.readlines()

    temp_fd, temp_path = tempfile.mkstemp(prefix="input_vm_proton_", suffix=".in", dir=working_dir)
    os.close(temp_fd)

    with open(temp_path, "w") as f:
        for line in lines:
            if line.lstrip().startswith("seed"):
                f.write(f"seed {seed}\n")
            elif line.lstrip().startswith("x_projectile_jimwlk"):
                f.write(f"x_projectile_jimwlk {x_pom}\n")
            elif line.lstrip().startswith("x_target_jimwlk"):
                f.write(f"x_target_jimwlk {x_pom}\n")
            elif line.lstrip().startswith("NucleusQsTableFileName"):
                f.write(f"NucleusQsTableFileName {os.path.join(ipglasma_path, 'qs2Adj_vs_Tp_vs_Y_200.in')}\n")
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


def run_command(cmd, cwd, stdout_path=None, stderr_path=None, slurm=True):
    print(f"Running: {' '.join(str(part) for part in cmd)}", " stdout:", stdout_path, " stderr:", stderr_path)
    srun_cmd = ["srun", "--ntasks=1", "--cpus-per-task=1"] + cmd
    if not slurm:
        srun_cmd = cmd
    result = subprocess.run(srun_cmd, cwd=cwd, capture_output=True, text=True)
    if stdout_path:
        with open(stdout_path, "w") as stdout_handle:
            stdout_handle.write(result.stdout)
    if stderr_path:
        with open(stderr_path, "w") as stderr_handle:
            stderr_handle.write(result.stderr)

    return result


def PlotComparison(datadir="", subnucleondiffraction_path="", reference_cross_section_file="", x_pom="0.001705"):
    print("Reading data from", datadir, "and comparing to reference cross section data in", reference_cross_section_file)
    run_command(
        ["python3", os.path.join(subnucleondiffraction_path, "python", "cross_section.py"), "-maxconf", "1999", "-dir", datadir, "-ds", "0.005", "-steps", "36"],
        cwd=os.getcwd(),
        stdout_path=os.path.join(datadir, "cross_section.dat"),
        slurm=False
    )
    print("Cross section saved to", os.path.join(datadir, "cross_section.dat"), flush=True)

    events_included = None
    try:
        with open(os.path.join(datadir, "cross_section.dat"), "r") as f:
            for line in f:
                if line.startswith("#") and "events included" in line:
                    events_included = int(line.split(":")[-1].strip())
                    break
    except Exception:
        pass

    xs = np.loadtxt(os.path.join(datadir, "cross_section.dat"))
    if len(xs) == 0:
        print("No cross section data found in cross_section.dat", flush=True)
        return

    coh = xs[:, 1] * 1e6
    coherr = xs[:, 2] * 1e6
    incoh = xs[:, 3] * 1e6
    incoherr = xs[:, 4] * 1e6

    validation = np.loadtxt(reference_cross_section_file)
    ref_coh = validation[:, 1] * 1e6
    ref_coh_err = validation[:, 2] * 1e6
    ref_incoh = validation[:, 3] * 1e6
    ref_incoh_err = validation[:, 4] * 1e6

    ref_coh_interp = np.interp(xs[:, 0], validation[:, 0], ref_coh)
    ref_coh_err_interp = np.interp(xs[:, 0], validation[:, 0], ref_coh_err)
    ref_incoh_interp = np.interp(xs[:, 0], validation[:, 0], ref_incoh)
    ref_incoh_err_interp = np.interp(xs[:, 0], validation[:, 0], ref_incoh_err)

    coh_ratio = coh / ref_coh_interp
    incoh_ratio = incoh / ref_incoh_interp
    coh_ratio_err = coh_ratio * np.sqrt((coherr / coh) ** 2 + (ref_coh_err_interp / ref_coh_interp) ** 2)
    incoh_ratio_err = incoh_ratio * np.sqrt((incoherr / incoh) ** 2 + (ref_incoh_err_interp / ref_incoh_interp) ** 2)

    fig, (ax, ax_ratio) = plt.subplots(2, 1, sharex=True, figsize=(8, 8), gridspec_kw={"height_ratios": [3, 1]})

    ax.plot(xs[:, 0], coh, label="Test run coh", color="black", linestyle="solid")
    ax.fill_between(xs[:, 0], coh - coherr, coh + coherr, color="gray", alpha=0.5)
    ax.plot(xs[:, 0], incoh, label="Test run incoh", color="red", linestyle="solid")
    ax.fill_between(xs[:, 0], incoh - incoherr, incoh + incoherr, color="pink", alpha=0.5)

    ax.errorbar(validation[:, 0], ref_coh, label="Reference coh", yerr=ref_coh_err, fmt="o", color="blue", linestyle="dashed", alpha=0.8)
    ax.errorbar(validation[:, 0], ref_incoh, label="Reference incoh", yerr=ref_incoh_err, fmt="o", color="blue", linestyle="dashed", alpha=0.5)

    ax.set_ylabel(r"$\frac{d\sigma}{dt}$ [$\mu \mathrm{b}/\mathrm{GeV}^2$]", fontsize=14)
    ax.tick_params(axis="both", labelsize=14)
    ax.set_ylim(bottom=0.1, top=5e2)
    ax.set_yscale("log")
    ax.set_title(r"Test run cross section, $N_\mathrm{events} = " + str(events_included) + r", x_\mathbb{P} = " + str(x_pom) + "$", fontsize=14)
    ax.legend(fontsize=14)

    ax_ratio.axhline(1.0, color="black", linestyle="--", linewidth=1)
    ax_ratio.plot(xs[:, 0], coh_ratio, color="black", label="coh ratio")
    ax_ratio.fill_between(xs[:, 0], coh_ratio - coh_ratio_err, coh_ratio + coh_ratio_err, color="gray", alpha=0.2)
    ax_ratio.plot(xs[:, 0], incoh_ratio, color="red", label="incoh ratio", linestyle="dashed")
    ax_ratio.fill_between(xs[:, 0], incoh_ratio - incoh_ratio_err, incoh_ratio + incoh_ratio_err, color="blue", alpha=0.3)
    ax_ratio.tick_params(axis="both", labelsize=14)
    ax_ratio.set_ylabel("Ratio", fontsize=14)
    ax_ratio.set_xlabel(r"$t [\mathrm{GeV}^2]$", fontsize=14)
    ax_ratio.set_ylim(0.5, 1.5)
    ax_ratio.grid(alpha=0.3)
    ax_ratio.legend(loc="best",fontsize=14)

    ax.set_xlim(left=0, right=xs[:, 0][-1])
    fig.tight_layout()
    fig.savefig(os.path.join(datadir, "cross_section.pdf"))
    print("Created cross section plot:", os.path.join(datadir, "cross_section.pdf"), "Events included:", events_included, flush=True)
    plt.close(fig)


def run_seed(seed, args, repo_root, datadir, subnucleondiffraction_cmd, ipglasma_binary, input_file_path, reference_file_path, ipglasma_path):
    worker_dir = os.path.join(datadir, f"seed_{seed}")
    os.makedirs(worker_dir, exist_ok=True)

    print(f"===== RUNNING SEED {seed} =====", flush=True)
    temp_input_path = generate_temp_input(source_input_file=input_file_path, seed=seed, x_pom=xpom, working_dir=worker_dir, ipglasma_path=ipglasma_path)

    log_path = os.path.join(worker_dir, f"ipglasma_log_{seed}")
    ipglasma_cmd = os.path.join(ipglasma_path, ipglasma_binary)
    run_command([ipglasma_cmd, temp_input_path], cwd=worker_dir, stdout_path=log_path, stderr_path=log_path + ".err")
    print("Removing temp files", temp_input_path, flush=True)
    remove_temp_input(temp_input_path)
    print("Done", seed, flush=True)
    id1 = 2 * seed + 1
    id2 = 2 * seed + 2

    print(f"===== RUNNING SUBNUCLEONDIFFRACTION FOR EVENTS {id1} AND {id2} =====", flush=True)
    spectra_1_path = os.path.join(worker_dir, f"spectra_{id1}")
    spectra_2_path = os.path.join(worker_dir, f"spectra_{id2}")
    run_command(subnucleondiffraction_cmd + ["-dipole", "1", "ipglasma_binary", os.path.join(worker_dir, f"Final_x_{xpom}_V-{id1}")], cwd=worker_dir, stdout_path=spectra_1_path, stderr_path=spectra_1_path + ".err")
    run_command(subnucleondiffraction_cmd + ["-dipole", "1", "ipglasma_binary", os.path.join(worker_dir, f"Final_x_{xpom}_V-{id2}")], cwd=worker_dir, stdout_path=spectra_2_path, stderr_path=spectra_2_path + ".err")

    for file_name in [f"spectra_{id1}", f"spectra_{id2}"]:
        src = os.path.join(worker_dir, file_name)
        dst = os.path.join(datadir, file_name)
        if os.path.exists(src):
            shutil.copy2(src, dst)

    if id1 in plot_ids or id2 in plot_ids:
        plot_id = id1 if id1 in plot_ids else id2
        final_file = os.path.join(worker_dir, f"Final_x_{xpom}_V-{plot_id}")
        plot_path = os.path.join(datadir, f"output_{plot_id}.pdf")
        tmp_output_path = os.path.join(worker_dir, f"tmp_output_{plot_id}")
        run_command(
            subnucleondiffraction_cmd + ["-dipole", "1", "ipglasma_binary", final_file, "-print_nucleus", "-Q2", "0"],
            cwd=worker_dir,
            stdout_path=tmp_output_path,
            stderr_path=tmp_output_path + ".err",
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


    print("Removing Wilson line files ", os.path.join(worker_dir, f"Final_x_{xpom}_V-{id1}"), " and ", os.path.join(worker_dir, f"Final_x_{xpom}_V-{id2}"), flush=True)
    if os.path.exists(os.path.join(worker_dir, f"Final_x_{xpom}_V-{id1}")):
        os.remove(os.path.join(worker_dir, f"Final_x_{xpom}_V-{id1}"))
    else:
        print(f"Warning: {os.path.join(worker_dir, f'Final_x_{xpom}_V-{id1}')} does not exist.", flush=True)
    if os.path.exists(os.path.join(worker_dir, f"Final_x_{xpom}_V-{id2}")):
        os.remove(os.path.join(worker_dir, f"Final_x_{xpom}_V-{id2}"))
    else:
        print(f"Warning: {os.path.join(worker_dir, f'Final_x_{xpom}_V-{id2}')} does not exist.", flush=True)

    return {"seed": seed, "id1": id1, "id2": id2}


def main():
    global subnucleondiffraction_path
    global ipglasma_path

    parser = argparse.ArgumentParser(description="Parallel test for IP-Glasma+JIMWLK+SubnucleonDiffraction for vector meson production.")
    parser.add_argument("-plot_only", action="store_true", help="Use existing data files for plotting the comparison.")
    parser.add_argument("-maxevents", type=int, default=200, help="Maximum number of events for the runs.")
    parser.add_argument("-datadir", type=str, default="./jpsi/", help="Directory to store data files.")
    parser.add_argument("-subnucleondiffraction_path", type=str, default=subnucleondiffraction_path, help="Path to the subnucleondiffraction executable.")
    parser.add_argument("-ipglasma_path", type=str, default=ipglasma_path, help="Path to the IP-Glasma executable.")
    parser.add_argument("-max_workers", type=int, default=None, help="Maximum number of parallel workers. Defaults to the number of logical CPUs.")
    parser.add_argument("-keep_logs", type=bool, default=False, help="Keep log files for each seed. Default is False.")
    args = parser.parse_args()

    repo_root = _repo_root()
    script_dir = _script_dir()
    datadir = os.path.abspath(args.datadir)
    subnucleondiffraction_path = os.path.abspath(args.subnucleondiffraction_path)
    ipglasma_path = os.path.abspath(args.ipglasma_path)
    max_workers = args.max_workers or min(cpu_count(), max(1, int(args.maxevents / 2)))
    maxseed = int(args.maxevents / 2)
    keep_logs = args.keep_logs

    if args.plot_only:
        PlotComparison(datadir=datadir, subnucleondiffraction_path=subnucleondiffraction_path, reference_cross_section_file=os.path.join(script_dir, reference_cross_section_file), x_pom=xpom)
        raise SystemExit(0)

    print("PID =", os.getpid())
    if hasattr(os, "sched_getaffinity"):
        print("Affinity =", os.sched_getaffinity(0))

    os.makedirs(datadir, exist_ok=True)
    shutil.copy2(os.path.join(script_dir, "input_vm_proton"), os.path.join(datadir, "input_vm_proton"))
    shutil.copy2(os.path.join(ipglasma_path, "qs2Adj_vs_Tp_vs_Y_200.in"), os.path.join(script_dir, "qs2Adj_vs_Tp_vs_Y_200.in"))
    #if os.path.exists(os.path.join(repo_root, "qs2Adj_vs_Tp_vs_Y_200.in")):
    #    shutil.copy2(os.path.join(repo_root, "qs2Adj_vs_Tp_vs_Y_200.in"), os.path.join(repo_root, "qs2Adj_vs_Tp_vs_Y_200.in"))

    ipglasma_binary = os.path.join(ipglasma_path, "ipglasma")
    subnucleondiffraction_binary = os.path.join(subnucleondiffraction_path, "build", "bin", "subnucleondiffraction")
    subnucleondiffraction_wavefile = os.path.join(subnucleondiffraction_path, "gauss-boosted_mzsat.dat")

    missing = []
    non_executable = []
    for required,exec in zip([ipglasma_binary, subnucleondiffraction_binary, subnucleondiffraction_wavefile, os.path.join(script_dir, "input_vm_proton")], [True, True, False, False]):
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

    subnucleondiffraction_cmd = [
        subnucleondiffraction_binary,
        "-Q2", "0",
        "-wavef_file",
        subnucleondiffraction_wavefile,
        "-maxt", "4.01",
        "-tstep", "0.05",
        "-mcintpoints", "1e6",
    ]

    input_file_path = os.path.join(script_dir, "input_vm_proton")
    reference_file_path = os.path.join(script_dir, reference_cross_section_file)

    failures=[]

    print(f"Running {maxseed} seeds in parallel with {max_workers} workers")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_seed, seed, args, repo_root, datadir, subnucleondiffraction_cmd, ipglasma_binary, input_file_path, reference_file_path, ipglasma_path) for seed in range(maxseed)]
        for future in as_completed(futures):
            try:
                outcome = "not available"
                outcome = future.result()
            except Exception as exc:
                print("Error occurred during seed execution:", exc)
                print(outcome, flush=True)
                failures.append(outcome)
                continue

            

    PlotComparison(datadir=datadir, subnucleondiffraction_path=subnucleondiffraction_path, reference_cross_section_file=reference_file_path, x_pom=xpom)

    if failures:
        failure_file = os.path.join(datadir, "failed_seeds.txt")
        with open(failure_file, "w") as f:
            for failure in failures:
                f.write(f"{failure}\n")
        print(f"{len(failures)} seed(s) failed; details: {failure_file}", flush=True)

    if not keep_logs:
        for seed in range(maxseed):
            worker_dir = os.path.join(datadir, f"seed_{seed}")
            if os.path.exists(worker_dir):
                shutil.rmtree(worker_dir)
                print(f"Removed temporary directory for seed {seed}: {worker_dir}", flush=True)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Validation of the index bookkeeping (issue #124): two exact symmetries of
the lattice formulation.

A Pb+Pb event at b = 0 writes its Wilson lines and geometry files, which
are then read back three times (readInitialWilsonLines 2):

  (0) unchanged: reproduces the event (the reference);
  (a) transposed, x <-> y: V'(i, j) = V(j, i) for both nuclei, and the
      geometry (nucleon positions, g^2 mu^2 and T_p maps) transposed;
  (b) projectile and target exchanged.

(a) maps the lattice action onto itself, so all outputs map onto the
reference point by point (to rounding; the hydro fields and the
eccentricities, which come from the Landau matching, to about 1e-7): eps, T^tautau, T^xy, tau^2 T^etaeta, T^taueta
at (y, x) equal those at (x, y), T^xx <-> T^yy, T^taux <-> T^tauy,
T^xeta <-> T^yeta, u^x <-> u^y (and likewise pi^munu), |eps_n| is
unchanged and Psi_n -> pi/2 - Psi_n (mod 2 pi/n), and the gluon spectrum
is the same. The exception are links where the light-cone solver has to
restart: the random start of a restart depends on the site and the link
direction, which the transposition changes, so such a link can end up at
another solution of the matching condition (issue #107), in the
strong-field region. (b) is the reflection eta -> -eta: T^taueta, T^xeta,
T^yeta, u^eta and pi^{mu eta} change sign, everything else is unchanged.
Since the light-cone solver then starts from U_B U_A instead of U_A U_B,
it can find another solution in the cells where the matching condition
has several (issue #107).

So (a) and (b) agree to rounding only outside the region the differing
links reach by the final time. The script reports how many points of the
T^munu grid differ beyond rounding, and accepts (a) and (b) if that
fraction is small (a general error in the index bookkeeping would change
nearly all points) and dN/dy, dE/dy, eps_2 and eps_3 change by less than
a tolerance. (0) has to agree to rounding.

Usage (from validations/lattice_symmetries/):
    # a quick test, 64^2 (about 10 s)
    python3 lattice_symmetries.py --ipglasma-path ../.. --quick
    # the setup of the paper, 256^2
    python3 lattice_symmetries.py --ipglasma-path ../..

--ipglasma-path is a directory with the ipglasma executable,
qs2Adj_vs_Tp_vs_Y_240.in and nucleusConfigurations/. Build IP-Glasma
with -DIPGLASMA_DETERMINISTIC_FFT=ON for the reference run (0) to agree
with the original event to the last bit (the script copies the FFT
wisdom file into every run directory); otherwise it agrees to rounding.
"""

import argparse
import glob
import json
import os
import shutil
import struct
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "utilities"))
from read_tmunu import COMPONENTS, read_tmunu  # noqa: E402

WILSON_HEADER = 32  # int32 N, int32 Nc, double L, a, x
GEOMETRY_MAGIC = b"IPGGEO1\0"


# --------------------------------------------------------------------------
# transforming the files
# --------------------------------------------------------------------------

def transpose_wilson_lines(src, dst):
    """V'(i, j) = V(j, i), header unchanged."""
    with open(src, "rb") as f:
        header = f.read(WILSON_HEADER)
        data = np.frombuffer(f.read(), dtype="=f8")
    n = struct.unpack("=i", header[:4])[0]
    data = data.reshape(n, n, 18).swapaxes(0, 1)
    with open(dst, "wb") as f:
        f.write(header)
        f.write(np.ascontiguousarray(data).tobytes())


def read_geometry(path):
    with open(path, "rb") as f:
        magic = f.read(8)
        if magic != GEOMETRY_MAGIC:
            raise ValueError(f"{path}: not a geometry file")
        length = struct.unpack("<Q", f.read(8))[0]
        meta = json.loads(f.read(length).decode())
        rest = np.frombuffer(f.read(), dtype="<f8")
    nucleons, n = int(meta["nucleons"]), int(meta["N"])
    positions = rest[:4 * nucleons].reshape(nucleons, 4).copy()
    maps = rest[4 * nucleons:].reshape(-1, n, n).copy()
    return meta, positions, maps


def write_geometry(path, meta, positions, maps):
    text = json.dumps(meta, separators=(",", ":")).encode()
    with open(path, "wb") as f:
        f.write(GEOMETRY_MAGIC)
        f.write(struct.pack("<Q", len(text)))
        f.write(text)
        f.write(np.ascontiguousarray(positions, dtype="<f8").tobytes())
        f.write(np.ascontiguousarray(maps, dtype="<f8").tobytes())


def transpose_geometry(src, dst):
    meta, positions, maps = read_geometry(src)
    write_geometry(dst, meta, positions[:, [1, 0, 2, 3]],
                   maps.swapaxes(1, 2))


def exchange_geometry(src, dst):
    meta, positions, maps = read_geometry(src)
    meta["nucleus"] = {"projectile": "target",
                       "target": "projectile"}[meta["nucleus"]]
    write_geometry(dst, meta, positions, maps)


def file_pairs(directory):
    """{"wilson": (projectile, target), "geometry": (projectile, target)}:
    the projectile has the odd number n, the target n + 1."""
    out = {}
    for kind, pattern in (("wilson", "WilsonLine_x_*"),
                          ("geometry", "WilsonLineGeometry_*")):
        files = sorted(glob.glob(os.path.join(directory, pattern)),
                       key=lambda p: int(p.rsplit("_", 1)[1]))
        if len(files) != 2:
            raise RuntimeError(f"expected 2 {pattern} files in {directory}")
        out[kind] = tuple(files)
    return out


def make_files(original, target_dir, case):
    os.makedirs(target_dir, exist_ok=True)
    pairs = file_pairs(original)
    for kind in ("wilson", "geometry"):
        proj, targ = pairs[kind]
        for src, partner in ((proj, targ), (targ, proj)):
            if case == "reference":
                shutil.copy(src, os.path.join(target_dir,
                                              os.path.basename(src)))
            elif case == "transposed":
                dst = os.path.join(target_dir, os.path.basename(src))
                (transpose_wilson_lines if kind == "wilson"
                 else transpose_geometry)(src, dst)
            else:  # exchanged: the projectile's file under the target's name
                dst = os.path.join(target_dir, os.path.basename(partner))
                if kind == "wilson":
                    shutil.copy(src, dst)
                else:
                    exchange_geometry(src, dst)


# --------------------------------------------------------------------------
# running IP-Glasma
# --------------------------------------------------------------------------

def write_input(template, overrides, path):
    remaining = dict(overrides)
    lines = []
    with open(template) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) >= 2 and tokens[0] in remaining:
                lines.append(f"{tokens[0]} {remaining.pop(tokens[0])}\n")
            else:
                lines.append(line)
    if remaining:
        raise KeyError(f"keys not in the template: {sorted(remaining)}")
    with open(path, "w") as f:
        f.writelines(lines)


def run(settings, workdir, overrides, wisdom=None):
    os.makedirs(workdir, exist_ok=True)
    write_input(settings["template"], dict(settings["overrides"], **overrides),
                os.path.join(workdir, "input"))
    if wisdom and os.path.exists(wisdom):
        shutil.copy(wisdom, workdir)
    env = dict(os.environ, OMP_NUM_THREADS=str(settings["threads"]))
    with open(os.path.join(workdir, "run.log"), "w") as log:
        result = subprocess.run([settings["ipglasma"], "input"], cwd=workdir,
                                stdout=log, stderr=subprocess.STDOUT, env=env)
    if result.returncode != 0:
        sys.exit(f"ipglasma failed in {workdir} (exit {result.returncode})")


# --------------------------------------------------------------------------
# comparing the outputs
# --------------------------------------------------------------------------

def read_outputs(workdir):
    tmunu, meta = read_tmunu(glob.glob(os.path.join(workdir,
                                                    "Tmunu-t*-0.ipgt"))[0])
    hydro = np.loadtxt(os.path.join(workdir,
                                    "epsilon-u-Hydro-TauHydro-0.dat"))
    ecc = np.loadtxt(os.path.join(workdir, "eccentricities0.dat"), ndmin=2)
    with open(os.path.join(workdir, "gluonMultiplicity0.json")) as f:
        gluons = json.load(f)
    # N_part, dN/dy, T_pp, b, dE/dy (the file also holds N/A placeholders)
    with open(glob.glob(os.path.join(workdir, "NpartdNdy-t*-0.dat"))[0]) as f:
        summary = np.array([float(v) for v in f.read().split()[:5]])
    return {"tmunu": tmunu, "hydro": hydro, "ecc": ecc[-1],
            "gluons": gluons, "summary": summary}


def numbers(obj):
    """All numbers of a JSON object, in a fixed order."""
    if isinstance(obj, dict):
        return [x for k in sorted(obj) for x in numbers(obj[k])]
    if isinstance(obj, list):
        return [x for v in obj for x in numbers(v)]
    if isinstance(obj, (int, float)) and not isinstance(obj, bool):
        return [float(obj)]
    return []


def relative(a, b):
    """max |a - b| / max |a| (per component along the last axis)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    scale = np.abs(a).reshape(-1, a.shape[-1]).max(axis=0) if a.ndim > 1 \
        else np.abs(a).max()
    scale = np.where(scale > 0, scale, 1.)
    diff = np.abs(a - b).reshape(-1, a.shape[-1]).max(axis=0) if a.ndim > 1 \
        else np.abs(a - b).max()
    return float(np.max(diff / scale))


# T^munu components: index permutation and signs of the two symmetries
C = {name: i for i, name in enumerate(COMPONENTS)}
TMUNU_SWAP = list(range(len(COMPONENTS)))
for x_name, y_name in (("Txx", "Tyy"), ("neg_T0x", "neg_T0y"),
                       ("neg_tau_Txeta", "neg_tau_Tyeta")):
    TMUNU_SWAP[C[x_name]], TMUNU_SWAP[C[y_name]] = C[y_name], C[x_name]
TMUNU_ETA = np.ones(len(COMPONENTS))
for name in ("neg_tau_T0eta", "neg_tau_Txeta", "neg_tau_Tyeta"):
    TMUNU_ETA[C[name]] = -1.
# hydro columns: eta x y e utau ux uy ueta pi00 pi0x pi0y pi0eta pixx pixy
# pixeta piyy piyeta pietaeta
HYDRO_SWAP = list(range(18))
for i, j in ((1, 2), (5, 6), (9, 10), (12, 15), (14, 16)):
    HYDRO_SWAP[i], HYDRO_SWAP[j] = j, i
HYDRO_ETA = np.ones(18)
HYDRO_ETA[[7, 11, 14, 16]] = -1.


def hydro_grid(hydro):
    """The hydro rows as a (nx, ny, 18) grid, ordered by x and y."""
    xs, ys = np.unique(hydro[:, 1]), np.unique(hydro[:, 2])
    grid = np.zeros((len(xs), len(ys), hydro.shape[1]))
    grid[np.searchsorted(xs, hydro[:, 1]), np.searchsorted(ys, hydro[:, 2])] \
        = hydro
    return grid


def expected(ref, case):
    """The reference outputs transformed by the symmetry of the case."""
    out = dict(ref)
    if case == "transposed":
        # read_tmunu gives [y, x, component]
        out["tmunu"] = ref["tmunu"].swapaxes(0, 1)[..., TMUNU_SWAP]
        out["hydro_grid"] = hydro_grid(ref["hydro"]).swapaxes(0, 1)[
            ..., HYDRO_SWAP]
        ecc = ref["ecc"].copy()
        for n in range(1, 7):
            col = 2 * n  # Psi_n (column 2n + 1 of the file)
            ecc[col] = np.mod(np.pi / 2 - ref["ecc"][col], 2 * np.pi / n)
        out["ecc"] = ecc
    elif case == "exchanged":
        out["tmunu"] = ref["tmunu"] * TMUNU_ETA
        out["hydro_grid"] = hydro_grid(ref["hydro"]) * HYDRO_ETA
    else:
        out["hydro_grid"] = hydro_grid(ref["hydro"])
    return out


def angle_difference(a, b, n):
    period = 2 * np.pi / n
    d = np.mod(a - b + period / 2, period) - period / 2
    return float(abs(d))


def compare(ref, got, case):
    exp = expected(ref, case)
    result = {
        "tmunu": relative(exp["tmunu"], got["tmunu"]),
        "hydro": relative(exp["hydro_grid"], hydro_grid(got["hydro"])),
        "gluon_spectrum": relative(numbers(exp["gluons"]),
                                   numbers(got["gluons"])),
        "dN_dE": relative(exp["summary"][[1, 4]], got["summary"][[1, 4]]),
    }
    eps = [abs(exp["ecc"][2 * n - 1] - got["ecc"][2 * n - 1])
           / max(abs(exp["ecc"][2 * n - 1]), 1e-300) for n in range(1, 7)]
    psi = [angle_difference(exp["ecc"][2 * n], got["ecc"][2 * n], n)
           for n in range(1, 7)]
    result["eps_n"] = float(max(eps))
    result["Psi_n"] = float(max(psi))
    # points of the T^munu grid that differ beyond rounding, relative to the
    # largest T^tautau
    scale = np.abs(exp["tmunu"][..., C["T00"]]).max()
    differ = np.abs(exp["tmunu"] - got["tmunu"]).max(axis=-1) > 1e-8 * scale
    result["tmunu_points_differing"] = float(differ.mean())
    result["observables"] = {
        "dN/dy": float(got["summary"][1] / exp["summary"][1] - 1.),
        "dE/dy": float(got["summary"][4] / exp["summary"][4] - 1.),
        "eps_2": float(got["ecc"][3] / exp["ecc"][3] - 1.),
        "eps_3": float(got["ecc"][5] / exp["ecc"][5] - 1.),
    }
    return result


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Exact lattice symmetries: transposition x <-> y and "
                    "exchange of projectile and target (issue #124).")
    parser.add_argument("--ipglasma-path", default="../..")
    parser.add_argument("--ipglasma-cmd", default="ipglasma")
    parser.add_argument("--template", default=os.path.join(HERE, "..", "..",
                                                           "input"),
                        help="input file to start from (the example input)")
    parser.add_argument("--datadir", default="./symmetries_output")
    parser.add_argument("--size", type=int, default=256)
    parser.add_argument("--L", type=float, default=30., help="[fm]")
    parser.add_argument("--maxtime", type=float, default=0.4, help="[fm/c]")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--threads", type=int, default=4,
                        help="OpenMP threads per run")
    parser.add_argument("--tolerance", type=float, default=1e-8,
                        help="relative agreement of T^munu, the gluon "
                             "spectrum, dN/dy and dE/dy required for (0) "
                             "and (a)")
    parser.add_argument("--matching-tolerance", type=float, default=1e-5,
                        help="the same for the hydro fields and eps_n, and "
                             "for Psi_n in rad: they come from the Landau "
                             "matching, whose eigenvalue solver amplifies "
                             "rounding errors")
    parser.add_argument("--observable-tolerance", type=float, default=0.03,
                        help="largest relative change of dN/dy, dE/dy, "
                             "eps_2 and eps_3 accepted for (a) and (b)")
    parser.add_argument("--max-differing", type=float, default=0.25,
                        help="largest fraction of T^munu points that may "
                             "differ beyond rounding in (a) and (b)")
    parser.add_argument("--quick", action="store_true",
                        help="small test: size 64, L 24, maxtime 0.3")
    args = parser.parse_args()
    if args.quick:
        args.size, args.L, args.maxtime = 64, 24., 0.3

    path = os.path.abspath(args.ipglasma_path)
    datadir = os.path.abspath(args.datadir)
    settings = {
        "ipglasma": os.path.join(path, args.ipglasma_cmd),
        "template": os.path.abspath(args.template),
        "threads": args.threads,
        "overrides": {
            "size": args.size, "L": args.L, "sizeOutput": args.size,
            "LOutput": args.L, "maxTime": args.maxtime, "seed": args.seed,
            "bMin": 0, "bMax": 0, "rotateReactionPlane": 0,
            "gaussianWounding": 0, "useJIMWLK": 0, "useFluctuatingX": 0,
            "writeWilsonLines": 2, "writeWilsonLineGeometry": 1,
            "writeTmunu": 1, "writeTmunuBinary": 1, "writeHydro": 1,
            "outputTimes": "none", "computeEccentricities": 1,
            "computeGluonMultiplicity": 1, "readInitialWilsonLines": 0,
            "wilsonLinePath": "./",
            "nucleusQsTableFileName": os.path.join(
                path, "qs2Adj_vs_Tp_vs_Y_240.in"),
            "nuclearConfigurationsPath": os.path.join(
                path, "nucleusConfigurations"),
        },
    }
    if os.path.exists(datadir):
        shutil.rmtree(datadir)
    original = os.path.join(datadir, "original")
    print(f"Pb+Pb at b = 0, {args.size}^2, L = {args.L} fm: original "
          f"event ...", flush=True)
    run(settings, original, {})
    wisdom = os.path.join(original, "ipglasma_fftw_wisdom.dat")

    outputs = {}
    for case in ("reference", "transposed", "exchanged"):
        files = os.path.join(datadir, case, "wilson_lines")
        make_files(original, files, case)
        print(f"reading back the {case} Wilson lines ...", flush=True)
        run(settings, os.path.join(datadir, case),
            {"readInitialWilsonLines": 2, "wilsonLinePath": files + "/"},
            wisdom)
        outputs[case] = read_outputs(os.path.join(datadir, case))
    outputs["original"] = read_outputs(original)

    results = {
        "(0) reference reproduces the event":
            compare(outputs["original"], outputs["reference"], "reference"),
        "(a) transposition x <-> y":
            compare(outputs["reference"], outputs["transposed"],
                    "transposed"),
        "(b) exchange of projectile and target":
            compare(outputs["reference"], outputs["exchanged"],
                    "exchanged"),
    }
    passed = True
    print()
    for name, r in results.items():
        print(f"{name}:")
        print(f"  max relative difference: T^munu {r['tmunu']:.1e}, "
              f"hydro {r['hydro']:.1e}, eps_n {r['eps_n']:.1e}, "
              f"Psi_n {r['Psi_n']:.1e} rad, gluon spectrum "
              f"{r['gluon_spectrum']:.1e}, dN/dy and dE/dy {r['dN_dE']:.1e}")
        print(f"  T^munu points differing beyond 1e-8: "
              f"{100 * r['tmunu_points_differing']:.2f}%")
        if name.startswith("(0)"):
            ok = all(r[k] <= args.tolerance
                     for k in ("tmunu", "gluon_spectrum", "dN_dE")) \
                and all(r[k] <= args.matching_tolerance
                        for k in ("hydro", "eps_n", "Psi_n"))
        else:
            changes = r["observables"]
            print("  " + ", ".join(f"{k} {100 * v:+.3f}%"
                                   for k, v in changes.items()))
            ok = all(abs(v) <= args.observable_tolerance
                     for v in changes.values()) \
                and r["tmunu_points_differing"] <= args.max_differing
        r["passed"] = bool(ok)
        passed &= ok
        print(f"  {'OK' if ok else 'FAIL'}")
    summary = {"settings": vars(args), "results": results,
               "passed": bool(passed)}
    with open(os.path.join(datadir, "symmetries_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(f"\n{'PASSED' if passed else 'FAILED'}; summary in "
          f"{datadir}/symmetries_summary.json")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

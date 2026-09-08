#!/usr/bin/env python3
"""
Helpers to read the energy-momentum tensor (Tmunu) snapshots written by
IP-Glasma at the end of (or during) the classical Yang-Mills evolution.

IP-Glasma writes one file per requested proper time tau, named

    Tmunu-t{tau_fm}-{event_id}.ipgt   (binary, default when writeTmunuBinary=1)
    Tmunu-t{tau_fm}-{event_id}.dat    (legacy text format)

in the directory the run was launched from. See src/MyEigen.cpp for the
writer. This module provides pure-Python/numpy readers for both formats
plus small convenience functions to pull out the energy density (the T00
component) as a 2D grid.

Binary format ("IPGTMU01"):
    8 bytes   magic = b"IPGTMU01"
    4 bytes   little-endian uint32: length of the following JSON header
    N bytes   UTF-8 JSON metadata, e.g.
              {"format": "ipglasma-tmunu", "version": 1, "dtype": "<f4",
               "shape": [ny, nx, 10], "axis_order": ["y", "x", "component"],
               "components": ["T00", "Txx", "Tyy", "tau2_Tetaeta",
                               "neg_T0x", "neg_T0y", "neg_tau_T0eta",
                               "neg_Txy", "neg_tau_Tyeta", "neg_tau_Txeta"],
               "tau_fm": ..., "eta_points": ..., "deta": ...,
               "dx_fm": ..., "dy_fm": ..., "event_id": ...}
    rest      shape[0]*shape[1]*10 little-endian float32 values, row-major
              in [y][x][component] order.

Text format (legacy):
    header line starting with '#', containing "key= value" tokens
    (etamax=, xmax=, ymax=, deta=, dx=, dy=, ...)
    one row per grid point: "ix iy T00 Txx Tyy tau2_Tetaeta neg_T0x neg_T0y
    neg_tau_T0eta neg_Txy neg_tau_Tyeta neg_tau_Txeta" (12 whitespace
    separated columns).
"""

import glob
import json
import os
import re
import struct

import numpy as np

MAGIC = b"IPGTMU01"

# Order of the 10 independent components IP-Glasma writes out. Index 0
# (T00) is the local energy density in GeV/fm^4.
COMPONENTS = [
    "T00", "Txx", "Tyy", "tau2_Tetaeta",
    "neg_T0x", "neg_T0y", "neg_tau_T0eta",
    "neg_Txy", "neg_tau_Tyeta", "neg_tau_Txeta",
]


def find_tmunu_files(output_dir, event_id):
    """Return all Tmunu snapshot files (binary or text) for one event,
    sorted by proper time tau (ascending)."""
    pattern = re.compile(
        r"Tmunu-t([0-9.eE+-]+)-{0}\.(ipgt|dat)$".format(int(event_id)))
    hits = []
    for path in glob.glob(os.path.join(output_dir, "Tmunu-t*-{0}.*".format(
            int(event_id)))):
        m = pattern.search(os.path.basename(path))
        if m:
            hits.append((float(m.group(1)), path))
    hits.sort(key=lambda entry: entry[0])
    return hits


def find_final_tmunu_file(output_dir, event_id):
    """Return the (tau, path) of the last (largest tau) Tmunu snapshot
    written for the given event, i.e. the energy density at the end of the
    classical Yang-Mills evolution. Raises FileNotFoundError if none
    exist."""
    hits = find_tmunu_files(output_dir, event_id)
    if not hits:
        raise FileNotFoundError(
            "No Tmunu-t*-{0}.ipgt/.dat files found in {1}".format(
                event_id, output_dir))
    return hits[-1]


def read_tmunu_binary(path):
    """Read a binary .ipgt Tmunu snapshot.

    Returns (array, meta) where array has shape (ny, nx, 10) and meta is
    the JSON header dict (includes dx_fm, dy_fm, tau_fm, event_id, ...).
    """
    with open(path, "rb") as f:
        magic = f.read(len(MAGIC))
        if magic != MAGIC:
            raise ValueError(
                "{0}: bad magic {1!r}, expected {2!r}".format(
                    path, magic, MAGIC))
        (meta_len,) = struct.unpack("<I", f.read(4))
        meta = json.loads(f.read(meta_len).decode("utf-8"))
        data = np.fromfile(f, dtype="<f4")
    ny, nx, ncomp = meta["shape"]
    expected = ny * nx * ncomp
    if data.size != expected:
        raise ValueError(
            "{0}: payload has {1} floats, expected {2} from header shape "
            "{3}".format(path, data.size, expected, meta["shape"]))
    array = data.reshape(ny, nx, ncomp)
    return array, meta


_HEADER_KV_RE = re.compile(r"(\w+)\s*=\s*([-+0-9.eEnaN]+)")


def _parse_text_header(header_line):
    values = {}
    for key, val in _HEADER_KV_RE.findall(header_line):
        try:
            values[key] = float(val)
        except ValueError:
            continue
    return values


def read_tmunu_text(path):
    """Read a legacy text Tmunu-t*.dat snapshot.

    Returns (array, meta) with the same layout as read_tmunu_binary:
    array has shape (ny, nx, 10).
    """
    with open(path, "r") as f:
        header_line = f.readline()
    header = _parse_text_header(header_line)
    data = np.loadtxt(path)
    ix = data[:, 0].astype(int)
    iy = data[:, 1].astype(int)
    nx = int(ix.max()) + 1
    ny = int(iy.max()) + 1
    array = np.zeros((ny, nx, 10), dtype=np.float32)
    array[iy, ix, :] = data[:, 2:12]

    m = re.search(r"Tmunu-t([0-9.eE+-]+)-(\d+)\.dat$", os.path.basename(path))
    tau_fm = float(m.group(1)) if m else header.get("tau")
    event_id = int(m.group(2)) if m else None

    meta = {
        "format": "ipglasma-tmunu-text",
        "shape": [ny, nx, 10],
        "axis_order": ["y", "x", "component"],
        "components": list(COMPONENTS),
        "tau_fm": tau_fm,
        "dx_fm": header.get("dx"),
        "dy_fm": header.get("dy"),
        "deta": header.get("deta"),
        "event_id": event_id,
    }
    return array, meta


def read_tmunu(path):
    """Dispatch to read_tmunu_binary/read_tmunu_text based on extension."""
    if path.endswith(".ipgt"):
        return read_tmunu_binary(path)
    if path.endswith(".dat"):
        return read_tmunu_text(path)
    raise ValueError(
        "Don't know how to read Tmunu file with unknown extension: "
        "{0}".format(path))


def get_energy_density(path):
    """Load a Tmunu snapshot and return just the energy density grid.

    Returns (energy_density, dx_fm, dy_fm, tau_fm) where energy_density is
    a 2D array of shape (ny, nx) in GeV/fm^4 (the T00 component), and
    dx_fm/dy_fm are the lattice spacings in fm.
    """
    array, meta = read_tmunu(path)
    energy_density = array[..., COMPONENTS.index("T00")]
    return energy_density, meta["dx_fm"], meta["dy_fm"], meta.get("tau_fm")


def get_final_energy_density(output_dir, event_id):
    """Convenience wrapper: find and load the last-tau Tmunu snapshot for
    one event. Returns (energy_density, dx_fm, dy_fm, tau_fm)."""
    _tau, path = find_final_tmunu_file(output_dir, event_id)
    return get_energy_density(path)

#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Read the Wilson-line files IP-Glasma writes (`writeWilsonLines 1` text or
`2` binary; format in OUTPUT.md, section "Wilson lines").

    from read_wilson_lines import read_wilson_lines
    V, meta = read_wilson_lines("WilsonLine_x_1.00000e-04_1")
    # V[ix, iy] is the 3x3 complex matrix at site (ix, iy);
    # meta = {"size", "Nc", "L_fm", "a_fm", "x"} (text files: "size" only)

As a script, it prints the header and the unitarity of the matrices:

    python3 read_wilson_lines.py WilsonLine_x_1.00000e-04_1
"""

import argparse
import sys

import numpy as np

# int32 size, int32 Nc, double L, double a, double x, in the host's byte
# order (OUTPUT.md)
_HEADER = np.dtype([("size", "=i4"), ("Nc", "=i4"), ("L", "=f8"),
                    ("a", "=f8"), ("x", "=f8")])


def read_wilson_lines(path):
    """Read one Wilson-line file.

    Returns (V, meta): V is a complex array of shape (N, N, 3, 3), indexed
    [ix, iy, row, column]; meta holds the header (binary files) or the
    size (text files, which have no header)."""
    if path.endswith(".txt"):
        data = np.loadtxt(path)
        n = int(round(np.sqrt(data.shape[0])))
        if n * n != data.shape[0] or data.shape[1] != 20:
            raise ValueError(f"{path}: not a text Wilson-line file")
        values = data[:, 2::2] + 1j * data[:, 3::2]
        return values.reshape(n, n, 3, 3), {"size": n}
    with open(path, "rb") as f:
        header = np.fromfile(f, dtype=_HEADER, count=1)
        if header.size != 1:
            raise ValueError(f"{path}: file too short for the header")
        n = int(header["size"][0])
        nc = int(header["Nc"][0])
        if nc != 3 or n <= 0:
            raise ValueError(f"{path}: unexpected header size={n}, Nc={nc}")
        values = np.fromfile(f, dtype="=f8", count=n * n * 9 * 2)
    if values.size != n * n * 9 * 2:
        raise ValueError(f"{path}: expected {n * n * 9} matrix elements")
    V = (values[0::2] + 1j * values[1::2]).reshape(n, n, 3, 3)
    meta = {"size": n, "Nc": nc, "L_fm": float(header["L"][0]),
            "a_fm": float(header["a"][0]), "x": float(header["x"][0])}
    return V, meta


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()
    for path in args.files:
        V, meta = read_wilson_lines(path)
        identity = np.eye(3)
        unitarity = np.abs(
            np.einsum("xyab,xycb->xyac", V, V.conj()) - identity).max()
        det = np.abs(np.linalg.det(V) - 1.).max()
        print(f"{path}: {meta}, max |V V^+ - 1| = {unitarity:.2e}, "
              f"max |det V - 1| = {det:.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

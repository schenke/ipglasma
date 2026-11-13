#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""
Makes the posterior files of the global fits to exclusive J/psi production
(posterior_jpsi_<arXiv>_Kfactor.csv, _K1.csv and _Kfactor_K.csv in tables/)
from the MCMC chains of their Zenodo records:

  2507.14087: posterior_chain_full_fit.tar.gz of
              https://doi.org/10.5281/zenodo.15880667
  2606.20362: analysis_and_data_files.tar.gz of
              https://doi.org/10.5281/zenodo.20696856

Each chain (chain.pkl, written by pocoMC) has weighted samples of the eleven
parameters m, BG, BGq, Nq, smearingWidth, QsMuRatio, jimwlkMass,
jimwlkLambdaQCD, UVDamp, omega and K. The script draws 1000 of them with
probabilities proportional to their weights (systematic resampling: the
first uniform number of numpy's default_rng(seed) as the offset, then the
points (offset + i)/1000 on the cumulative weights), so that the sets are
equally probable (weight 1); no weight exceeds 1/1000 of their sum, so no
sample is drawn twice, and the sets keep the order of the chain.

Set 0 is the MAP parametrization, a sample of the chain, with weight 0, so
that posteriorSet -1 never draws it: for arXiv:2507.14087 the sample that
the MAP of its Table I gives rounded (the best of 25 runs of the full model
with posterior samples; checked against the rounded values), for
arXiv:2606.20362 the sample with the largest likelihood (checked against
the max_likelihood files of the record). Nq, omega and UVDamp (and K with
K = 1) were fixed in the fits, so their chain values are prior draws and
are replaced by the fixed values.

The tar file has to have the MD5 sum of the Zenodo record. The pickles are
loaded with an unpickler that only constructs numpy arrays and dtypes, and
only numeric arrays are used.

Usage (from tables/):
    python3 make_jpsi_posterior_tables.py 2507.14087 \\
        path/to/posterior_chain_full_fit.tar.gz
    python3 make_jpsi_posterior_tables.py 2606.20362 \\
        path/to/analysis_and_data_files.tar.gz
"""

import argparse
import hashlib
import io
import os
import pickle
import tarfile
import textwrap

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NSETS = 1000
# the columns of chain.pkl
CHAIN_COLUMNS = ["m", "BG", "BGq", "Nq", "smearingWidth", "QsMuRatio",
                 "jimwlkMass", "jimwlkLambdaQCD", "UVDamp", "omega", "K"]
# the fitted IP-Glasma parameters, then the parameters fixed in the fits and
# the settings of the fits that a posterior file can set
FITTED = ["m", "BG", "BGq", "smearingWidth", "QsMuRatio", "jimwlkMass",
          "jimwlkLambdaQCD"]
FIXED = {"Nq": "3", "dqMin": "0", "BGqVar": "0", "NqFluc": "0", "omega": "1",
         "UVDamp": "0", "jimwlkAlphaS": "0", "jimwlkMu0": "0.28",
         "jimwlkC": "0.2", "jimwlkInitialX": "0.01"}
UNITS = {"m": "GeV", "BG": "GeV^-2", "BGq": "GeV^-2", "dqMin": "fm",
         "jimwlkMass": "GeV", "jimwlkLambdaQCD": "GeV", "BGqVar": "GeV^-4",
         "UVDamp": "GeV^-1", "jimwlkMu0": "GeV"}

# The MAP of Table I of arXiv:2507.14087 is the run of the full model with
# the smallest chi^2/dof of 25 runs with posterior samples, given rounded.
TABLE_I = "the MAP of Table I of arXiv:2507.14087 gives it rounded (of 25 " \
          "runs of the full model with posterior samples, the one with the " \
          "smallest chi^2/dof)"
RECORDS = {
    "2507.14087": {
        "authors": "H. Mantysaari, H. Roch, F. Salazar, B. Schenke, C. Shen "
                   "and W. Zhao",
        "reference": "Phys. Rev. D 113 (2026) 014038, arXiv:2507.14087",
        "doi": "10.5281/zenodo.15880667",
        "tar": "posterior_chain_full_fit.tar.gz",
        "md5": "16c31bbeea795a0161bf2a777d56de79",
        "seed": 250714087,
        "fits": {
            "Kfactor": {
                "chain": "mcmc_full_Kfactor/chain.pkl",
                "what": "with the cross sections scaled by a free factor K "
                        "(column \"variable K\" of Table I)",
                # the paper refitted K for this run
                "map": {"sample": 39324, "K": 0.308,
                        "table": {"m": "0.51", "BG": "3.83", "BGq": "0.30",
                                  "smearingWidth": "0.88",
                                  "QsMuRatio": "0.37", "jimwlkMass": "0.16",
                                  "jimwlkLambdaQCD": "0.0348"}},
                "map_source": TABLE_I,
            },
            "K1": {
                "chain": "mcmc_full_vanilla/chain.pkl",
                "what": "with K = 1 (column \"K = 1\" of Table I)",
                "map": {"sample": 109648,
                        "table": {"m": "0.31", "BG": "2.99", "BGq": "0.07",
                                  "smearingWidth": "0.99",
                                  "QsMuRatio": "0.63", "jimwlkMass": "0.12",
                                  "jimwlkLambdaQCD": "0.0097"}},
                "map_source": TABLE_I,
            },
        },
    },
    "2606.20362": {
        "authors": "H. Mantysaari, H. Roch, B. Schenke, C. Shen and W. Zhao",
        "reference": "Phys. Rev. D 114 (2026) L031505, arXiv:2606.20362",
        "doi": "10.5281/zenodo.20696856",
        "tar": "analysis_and_data_files.tar.gz",
        "md5": "d0352c690a514f185e17652ed26b1f75",
        "seed": 260620362,
        "fits": {
            "Kfactor": {
                "chain": "MCMC_SETUP/mcmc_PCGP_full_Kfactor_tCorr/chain.pkl",
                "what": "with the gamma+Pb data corrected for electromagnetic "
                        "dissociation and the cross sections scaled by a "
                        "free factor K",
                "map": {"max_likelihood": "OSG_RUNS/max_likelihood_K_tCorr"},
                "map_source": "the sample with the largest likelihood (the "
                              "prior is flat), as in OSG_RUNS/"
                              "max_likelihood_K_tCorr of the record",
            },
            "K1": {
                "chain": "MCMC_SETUP/mcmc_PCGP_full_vanilla_tCorr/chain.pkl",
                "what": "with the gamma+Pb data corrected for electromagnetic "
                        "dissociation and K = 1",
                "map": {"max_likelihood":
                        "OSG_RUNS/max_likelihood_vanilla_tCorr"},
                "map_source": "the sample with the largest likelihood (the "
                              "prior is flat), as in OSG_RUNS/"
                              "max_likelihood_vanilla_tCorr of the record",
            },
        },
    },
}


class NumpyUnpickler(pickle.Unpickler):
    """Constructs numpy arrays, dtypes and scalars, and nothing else."""
    ALLOWED = {("numpy._core.multiarray", "_reconstruct"),
               ("numpy.core.multiarray", "_reconstruct"),
               ("numpy._core.multiarray", "scalar"),
               ("numpy.core.multiarray", "scalar"),
               ("numpy", "dtype"), ("numpy", "ndarray")}

    def find_class(self, module, name):
        if (module, name) not in self.ALLOWED:
            raise pickle.UnpicklingError(f"not allowed: {module}.{name}")
        try:
            return super().find_class(module, name)
        except ModuleNotFoundError:
            # written with numpy 2 (numpy._core), read with numpy 1
            return super().find_class(
                module.replace("numpy._core", "numpy.core"), name)


def check_md5(tar_path, expected):
    md5 = hashlib.md5()
    with open(tar_path, "rb") as tar:
        for block in iter(lambda: tar.read(1 << 20), b""):
            md5.update(block)
    if md5.hexdigest() != expected:
        raise SystemExit(f"{tar_path}: MD5 {md5.hexdigest()}, but the Zenodo "
                         f"record has {expected}")


def numeric(chain, key, member):
    """The array `key` of the chain as floats, if it is a numeric array
    (dtype kind f, i or u, so its data are numbers, not object pointers)."""
    array = chain.get(key)
    if not isinstance(array, np.ndarray) or array.dtype.kind not in "fiu":
        raise ValueError(f"{member}: {key} is no numeric numpy array")
    return array.astype(float)


def read_chain(tar_path, member):
    with tarfile.open(tar_path) as tar:
        data = tar.extractfile(member).read()
    chain = NumpyUnpickler(io.BytesIO(data)).load()
    if not isinstance(chain, dict):
        raise ValueError(f"{member}: no dictionary")
    samples, weights, logl = (numeric(chain, key, member)
                              for key in ("chain", "weights", "logl"))
    if samples.ndim != 2 or samples.shape[1] != len(CHAIN_COLUMNS):
        raise ValueError(f"{member}: expected {len(CHAIN_COLUMNS)} columns, "
                         f"got shape {samples.shape}")
    if weights.shape != samples.shape[:1] or logl.shape != samples.shape[:1]:
        raise ValueError(f"{member}: one weight and logl per sample expected")
    if (not np.all(np.isfinite(weights)) or np.any(weights < 0)
            or weights.sum() <= 0 or np.any(np.isnan(logl))):
        raise ValueError(f"{member}: weights not finite and non-negative, "
                         "or logl NaN")
    return samples, weights, logl


def read_parameter_file(tar_path, member):
    """The values of a parameter file of the record (lines "name value", in
    the order of the chain)."""
    with tarfile.open(tar_path) as tar:
        lines = tar.extractfile(member).read().decode().split("\n")
    values = [float(line.split()[1]) for line in lines if line.strip()]
    if len(values) != len(CHAIN_COLUMNS):
        raise ValueError(f"{member}: expected {len(CHAIN_COLUMNS)} values")
    return values


def find_map(tar_path, fit, spec, samples, logl):
    """The row of the chain of the MAP parametrization, checked against the
    rounded values of the paper or the max_likelihood file of the record."""
    if "sample" in spec:
        best = spec["sample"]
        for name, text in spec["table"].items():
            value = samples[best, CHAIN_COLUMNS.index(name)]
            digits = len(text.split(".")[1])
            if abs(value - float(text)) > 0.5 * 10.0**-digits + 1e-12:
                raise ValueError(f"sample {best}: {name} {value} does not "
                                 f"round to {text}")
        return best
    best = int(np.argmax(logl))
    values = read_parameter_file(tar_path, spec["max_likelihood"])
    for name in FITTED + (["K"] if fit == "Kfactor" else []):
        column = CHAIN_COLUMNS.index(name)
        if not np.isclose(samples[best, column], values[column],
                          rtol=1e-12, atol=0):
            raise ValueError(f"sample {best}: {name} differs from "
                             f"{spec['max_likelihood']}")
    return best


def systematic_resample(weights, n, rng):
    """Indices of n samples drawn with probabilities proportional to the
    weights; each sample at most once if no weight exceeds 1/n of the sum."""
    cumulative = np.cumsum(weights / weights.sum())
    cumulative[-1] = 1.
    positions = (rng.uniform() + np.arange(n)) / n
    indices = np.searchsorted(cumulative, positions)
    if len(np.unique(indices)) != n:
        raise ValueError("a weight exceeds 1/n of the sum")
    return indices


def fmt(value):
    return f"{value:.6g}"


def header(record, fit, nchain):
    """The paragraphs of the comment block that both files of a fit share."""
    r, f = RECORDS[record], RECORDS[record]["fits"][fit]
    return [
        "Posterior of the global fit to exclusive J/psi production in "
        f"gamma+p and gamma+Pb collisions of {r['authors']}, "
        f"{r['reference']}, {f['what']}.",
        f"{NSETS} samples drawn with probabilities proportional to their "
        f"weights (systematic resampling, numpy default_rng({r['seed']}); in "
        f"the order of the chain) from the {nchain} weighted samples of "
        f"{f['chain']} in {r['tar']} of https://doi.org/{r['doi']}. The sets "
        "are equally probable (weight 1), so posteriorSet -1 draws them "
        "uniformly and events of fixed sets need no weights. Made with "
        "make_jpsi_posterior_tables.py; see README.md.",
    ]


def comment(paragraphs):
    """The paragraphs as comment lines of at most 88 characters, without a
    line break before a "-1"."""
    return "".join(f"# {line.replace(chr(0xa0), ' ')}\n"
                   for paragraph in paragraphs
                   for line in textwrap.wrap(
                       paragraph.replace(" -1", "\xa0-1"), 86))


def write_tables(record, tar_path, outdir):
    r = RECORDS[record]
    check_md5(tar_path, r["md5"])
    k_column = CHAIN_COLUMNS.index("K")
    for fit in ("Kfactor", "K1"):
        f = r["fits"][fit]
        samples, weights, logl = read_chain(tar_path, f["chain"])
        indices = systematic_resample(weights, NSETS,
                                      np.random.default_rng(r["seed"]))
        best = find_map(tar_path, fit, f["map"], samples, logl)
        map_k = f["map"].get("K", samples[best, k_column])
        refitted_k = (f", with K refitted to {fmt(map_k)} "
                      f"({fmt(samples[best, k_column])} in the chain)"
                      if "K" in f["map"] else "")
        stem = f"posterior_jpsi_{record}_{fit}"
        k_file = f"posterior_jpsi_{record}_Kfactor_K.csv"
        columns = FITTED + list(FIXED) + ["weight"]
        names = [f"{c} [{UNITS[c]}]" if c in UNITS else c for c in columns]
        text = header(record, fit, len(samples)) + [
            "Values to 6 digits. Nq, omega and UVDamp were fixed in these "
            "fits (their chain values are prior draws) and are constant "
            "columns with the fixed values, as are the other settings of the "
            "fits that a posterior file can set: dqMin, BGqVar, NqFluc and "
            "jimwlkAlphaS 0 (running coupling), jimwlkMu0 0.28, jimwlkC 0.2 "
            "and jimwlkInitialX 0.01. "
            + ("The K factor, which is no IP-Glasma parameter, is in "
               f"{k_file}." if fit == "Kfactor" else "K was fixed to 1.")
            + " The fits used hot-spot nucleons (nucleonModel hotspots, "
            "smearQs 1) and JIMWLK (useJIMWLK 1) with the square-root running "
            "coupling (jimwlkRunningCoupling 1) and nFlavors 3.",
            "The first line (posteriorSet 0) is the MAP parametrization, "
            f"sample {best} of the chain: {f['map_source']}{refitted_k}. It "
            "has weight 0, so that posteriorSet -1 never draws it.",
        ]
        rows = [[fmt(samples[best, CHAIN_COLUMNS.index(c)]) for c in FITTED]
                + list(FIXED.values()) + ["0"]]
        for i in indices:
            rows.append([fmt(samples[i, CHAIN_COLUMNS.index(c)])
                         for c in FITTED] + list(FIXED.values()) + ["1"])
        with open(os.path.join(outdir, stem + ".csv"), "w") as out:
            out.write(comment(text))
            out.write(", ".join(names) + "\n")
            out.write("".join(", ".join(row) + "\n" for row in rows))
        if fit == "Kfactor":
            text = header(record, fit, len(samples)) + [
                f"The K factor of each set of {stem}.csv, in the same order: "
                "posteriorSet is the index of the set (as in "
                "usedParameters<event>.dat), sample the row of the sample in "
                "chain.pkl (counted from 0), K the factor that scales the "
                "cross sections. Set 0 is the MAP parametrization"
                + (f"; its K is the one refitted in the paper "
                   f"({fmt(samples[best, k_column])} in the chain)"
                   if "K" in f["map"] else "")
                + ". Not an IP-Glasma input file.",
            ]
            with open(os.path.join(outdir, k_file), "w") as out:
                out.write(comment(text))
                out.write("# posteriorSet, sample, K\n")
                out.write(f"0, {best}, {fmt(map_k)}\n")
                out.write("".join(
                    f"{n + 1}, {i}, {fmt(samples[i, k_column])}\n"
                    for n, i in enumerate(indices)))
        print(f"wrote {stem}.csv ({len(samples)} chain samples, MAP sample "
              f"{best})")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("record", choices=sorted(RECORDS))
    parser.add_argument("tar", help="the tar file of the Zenodo record")
    parser.add_argument("--outdir", default=HERE)
    args = parser.parse_args()
    write_tables(args.record, args.tar, args.outdir)


if __name__ == "__main__":
    main()

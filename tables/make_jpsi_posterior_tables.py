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
points (offset + i)/1000 on the cumulative weights), so that the sets are equally
probable (weight 1); no weight exceeds 1/1000 of their sum, so no sample is
drawn twice, and the sets keep the order of the chain. Set 0 is the MAP
parametrization with weight 0, so that posteriorSet -1 never draws it: for
arXiv:2507.14087 the MAP of its Table I (from full-model runs of posterior
samples), for arXiv:2606.20362 the sample of the chain with the largest
likelihood. Nq, omega and UVDamp (and K with K = 1) were fixed in the fits,
so their chain values are prior draws and are replaced by the fixed values.

The pickles are loaded with an unpickler that only allows numpy arrays, so
that loading them cannot run code.

Usage (from tables/):
    python3 make_jpsi_posterior_tables.py 2507.14087 \\
        path/to/posterior_chain_full_fit.tar.gz
    python3 make_jpsi_posterior_tables.py 2606.20362 \\
        path/to/analysis_and_data_files.tar.gz
"""

import argparse
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
FIXED = {"Nq": "3", "BGqVar": "0", "NqFluc": "0", "omega": "1", "UVDamp": "0",
         "jimwlkAlphaS": "0", "jimwlkMu0": "0.28", "jimwlkC": "0.2"}
UNITS = {"m": "GeV", "BG": "GeV^-2", "BGq": "GeV^-2", "jimwlkMass": "GeV",
         "jimwlkLambdaQCD": "GeV", "BGqVar": "GeV^-4", "UVDamp": "GeV^-1",
         "jimwlkMu0": "GeV"}

RECORDS = {
    "2507.14087": {
        "authors": "H. Mantysaari, H. Roch, F. Salazar, B. Schenke, C. Shen "
                   "and W. Zhao",
        "reference": "arXiv:2507.14087",
        "doi": "10.5281/zenodo.15880667",
        "tar": "posterior_chain_full_fit.tar.gz",
        "seed": 250714087,
        "fits": {
            "Kfactor": {
                "chain": "mcmc_full_Kfactor/chain.pkl",
                "what": "with the cross sections scaled by a free factor K "
                        "(column \"variable K\" of Table I)",
                # Table I: the best of 25 full-model runs of posterior
                # samples, with K refitted
                "map": {"m": 0.51, "BG": 3.83, "BGq": 0.30,
                        "smearingWidth": 0.88, "QsMuRatio": 0.37,
                        "jimwlkMass": 0.16, "jimwlkLambdaQCD": 0.0348,
                        "K": 0.308},
            },
            "K1": {
                "chain": "mcmc_full_vanilla/chain.pkl",
                "what": "with K = 1 (column \"K = 1\" of Table I)",
                "map": {"m": 0.31, "BG": 2.99, "BGq": 0.07,
                        "smearingWidth": 0.99, "QsMuRatio": 0.63,
                        "jimwlkMass": 0.12, "jimwlkLambdaQCD": 0.0097,
                        "K": 1.0},
            },
        },
        "map_source": "the MAP of Table I of arXiv:2507.14087 (the posterior "
                      "sample with the smallest chi^2/dof of 25 runs of the "
                      "full model, with K refitted for the free K; no sample "
                      "of the chain)",
    },
    "2606.20362": {
        "authors": "H. Mantysaari, H. Roch, B. Schenke, C. Shen and W. Zhao",
        "reference": "Phys. Rev. D 114 (2026) L031505, arXiv:2606.20362",
        "doi": "10.5281/zenodo.20696856",
        "tar": "analysis_and_data_files.tar.gz",
        "seed": 260620362,
        "fits": {
            "Kfactor": {
                "chain": "MCMC_SETUP/mcmc_PCGP_full_Kfactor_tCorr/chain.pkl",
                "what": "with the gamma+Pb data corrected for electromagnetic "
                        "dissociation and the cross sections scaled by a "
                        "free factor K",
                "map": None,
            },
            "K1": {
                "chain": "MCMC_SETUP/mcmc_PCGP_full_vanilla_tCorr/chain.pkl",
                "what": "with the gamma+Pb data corrected for electromagnetic "
                        "dissociation and K = 1",
                "map": None,
            },
        },
        "map_source": "the sample of the chain with the largest likelihood "
                      "(the prior is flat)",
    },
}


class NumpyUnpickler(pickle.Unpickler):
    """Unpickles numpy arrays and scalars, and nothing else."""
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


def read_chain(tar_path, member):
    with tarfile.open(tar_path) as tar:
        data = tar.extractfile(member).read()
    chain = NumpyUnpickler(io.BytesIO(data)).load()
    samples = np.asarray(chain["chain"], dtype=float)
    if samples.ndim != 2 or samples.shape[1] != len(CHAIN_COLUMNS):
        raise ValueError(f"{member}: expected {len(CHAIN_COLUMNS)} columns, "
                         f"got shape {samples.shape}")
    return (samples, np.asarray(chain["weights"], dtype=float),
            np.asarray(chain["logl"], dtype=float))


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
    """The paragraphs as comment lines of at most 88 characters."""
    return "".join(f"# {line}\n" for paragraph in paragraphs
                   for line in textwrap.wrap(paragraph, 86))


def write_tables(record, tar_path, outdir):
    r = RECORDS[record]
    for fit in ("Kfactor", "K1"):
        f = r["fits"][fit]
        samples, weights, logl = read_chain(tar_path, f["chain"])
        indices = systematic_resample(weights, NSETS,
                                      np.random.default_rng(r["seed"]))
        if f["map"] is None:
            best = int(np.argmax(logl))
            map_set = dict(zip(CHAIN_COLUMNS, samples[best]))
            map_sample = best
        else:
            map_set, map_sample = f["map"], -1
        if fit == "K1":
            map_set = dict(map_set, K=1.0)
        stem = f"posterior_jpsi_{record}_{fit}"
        k_file = f"posterior_jpsi_{record}_Kfactor_K.csv"
        columns = FITTED + list(FIXED) + ["weight"]
        names = [f"{c} [{UNITS[c]}]" if c in UNITS else c for c in columns]
        text = header(record, fit, len(samples)) + [
            "Values to 6 digits. Left out are the columns Nq, omega and vUV, "
            "fixed in these fits, and the K factor, which is no IP-Glasma "
            "parameter" + (f" (see {k_file})" if fit == "Kfactor" else "")
            + "; the parameters fixed in the fits are constant columns. The "
            "fits used hot-spot nucleons (nucleonModel hotspots, smearQs 1) "
            "and JIMWLK from x = 0.01 (useJIMWLK 1, jimwlkInitialX 0.01) with "
            "the square-root running coupling (jimwlkRunningCoupling 1) and "
            "nFlavors 3.",
            "The first line (posteriorSet 0) is the MAP parametrization, "
            f"{r['map_source']}, with weight 0, so that posteriorSet -1 never "
            "draws it.",
        ]
        rows = [[fmt(map_set[c]) for c in FITTED] + list(FIXED.values())
                + ["0"]]
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
                + (", which is no sample of the chain (sample -1)"
                   if map_sample < 0 else "")
                + ". Not an IP-Glasma input file.",
            ]
            k_column = CHAIN_COLUMNS.index("K")
            with open(os.path.join(outdir, k_file), "w") as out:
                out.write(comment(text))
                out.write("# posteriorSet, sample, K\n")
                out.write(f"0, {map_sample}, {fmt(map_set['K'])}\n")
                out.write("".join(
                    f"{n + 1}, {i}, {fmt(samples[i, k_column])}\n"
                    for n, i in enumerate(indices)))
        print(f"wrote {stem}.csv ({len(samples)} chain samples, MAP "
              + (f"sample {map_sample}" if map_sample >= 0 else "of the paper")
              + ")")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("record", choices=sorted(RECORDS))
    parser.add_argument("tar", help="the tar file of the Zenodo record")
    parser.add_argument("--outdir", default=HERE)
    args = parser.parse_args()
    write_tables(args.record, args.tar, args.outdir)


if __name__ == "__main__":
    main()

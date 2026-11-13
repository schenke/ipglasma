# Posterior parameter sets

The CSV files in this folder are parameter sets of Bayesian analyses for the
input parameter `posteriorFile`; the README of the repository describes the
format and how the sets are used (section "Posterior parameter sets"). In
short:

- `posteriorFile tables/<file>` takes the parameters named in the columns of
  the file from one set, instead of the values of the input file.
- `posteriorSet n` uses set `n` (counted from 0); `posteriorSet -1` draws a
  set for every event, with the probabilities of the column `weight` if the
  file has one. `usedParameters<event>.dat` records the set of each event.
- Set 0 is the maximum a posteriori (MAP) parametrization, with `weight` 0
  so that `posteriorSet -1` never draws it; the posterior samples follow with
  `weight` 1, so they are equally probable and events of fixed sets need no
  weights. `posterior5020_Nq3.csv` has no `weight` column: its set 0 is the
  parametrization of arXiv:2409.19064, and `posteriorSet -1` draws each of
  its sets with the same probability.
- The comment lines at the top of each file give its source and the
  settings of the fit.

| File | Analysis | Sets | Parameters | The input needs |
|---|---|---|---|---|
| `posterior_jpsi_2202.01998_Nq_free.csv` | exclusive J/ψ production in γ+p, number of hot spots free | MAP + 1000 | `m`, `BG`, `BGq`, `smearingWidth`, `Nq`, `QsMuRatio`, `dqMin` (`BGqVar`, `NqFluc` 0) | `nucleonModel hotspots`, `smearQs 1` |
| `posterior_jpsi_2202.01998_Nq_3.csv` | the same with `Nq` = 3 | MAP + 1000 | the same | the same |
| `posterior5020_Nq3.csv` | proton structure at x_p = 2×10⁻⁴ (LHC energies) with `Nq` = 3, from J/ψ pseudodata of JIMWLK evolution | the set of arXiv:2409.19064 + 999, no weights | the same | the same |
| `posterior_jpsi_2507.14087_Kfactor.csv` | global fit of J/ψ production in γ+p and γ+Pb, free K factor | MAP + 1000 | `m`, `BG`, `BGq`, `smearingWidth`, `QsMuRatio`, `jimwlkMass`, `jimwlkLambdaQCD` (`Nq` 3, `omega` 1, `UVDamp` 0, `jimwlkAlphaS` 0, `jimwlkMu0` 0.28, `jimwlkC` 0.2, `BGqVar`, `NqFluc` 0) | `nucleonModel hotspots`, `smearQs 1`, `useJIMWLK 1`, `jimwlkInitialX 0.01`, `jimwlkRunningCoupling 1` |
| `posterior_jpsi_2507.14087_K1.csv` | the same with K = 1 | MAP + 1000 | the same | the same |
| `posterior_jpsi_2606.20362_Kfactor.csv` | global fit of J/ψ production in γ+p and γ+Pb with the γ+Pb data corrected for electromagnetic dissociation, free K factor | MAP + 1000 | the same | the same |
| `posterior_jpsi_2606.20362_K1.csv` | the same with K = 1 | MAP + 1000 | the same | the same |

## Sources

- **arXiv:2202.01998**: H. Mäntysaari, B. Schenke, C. Shen, W. Zhao,
  Phys. Lett. B 833 (2022) 137348, the proton shape fluctuations from
  exclusive J/ψ production in γ+p.
- **posterior5020_Nq3.csv** (the sub-nucleon structure at 5 TeV): set 0 is
  the parametrization of Table 2 of arXiv:2409.19064 (H. Mäntysaari,
  B. Schenke, C. Shen, W. Zhao, Phys. Rev. C 110 (2024) 054913), where it is
  rounded to 3 digits. There, J/ψ photoproduction computed with JIMWLK
  evolution and constrained by HERA and LHC ultraperipheral data gives
  pseudodata at x_p = 2×10⁻⁴, a typical x of heavy-ion collisions at the LHC,
  and the parameters are fitted to them as in arXiv:2202.01998. Where the
  other 999 sets come from is not documented yet. The iEBE-MUSIC
  configurations of arXiv:2409.19064, arXiv:2410.22160, arXiv:2506.02930 and
  arXiv:2508.05292 (`config/arXiv_<number>/`) use set 0.
- **arXiv:2507.14087**: H. Mäntysaari, H. Roch, F. Salazar, B. Schenke,
  C. Shen, W. Zhao, the global fit to exclusive J/ψ production in γ+p and
  γ+Pb collisions; the MCMC chains are in `posterior_chain_full_fit.tar.gz`
  of [Zenodo record 15880667](https://doi.org/10.5281/zenodo.15880667)
  (`mcmc_full_Kfactor/` and `mcmc_full_vanilla/`). Set 0 is the MAP of
  Table I of the paper, the posterior sample with the smallest χ²/dof among
  25 runs of the full model (with K refitted for the free K).
- **arXiv:2606.20362**: H. Mäntysaari, H. Roch, B. Schenke, C. Shen,
  W. Zhao, Phys. Rev. D 114 (2026) L031505, the same global fit with the
  γ+Pb data corrected for electromagnetic dissociation; the chains are in
  `analysis_and_data_files.tar.gz` of
  [Zenodo record 20696856](https://doi.org/10.5281/zenodo.20696856)
  (`MCMC_SETUP/mcmc_PCGP_full_Kfactor_tCorr/` and
  `MCMC_SETUP/mcmc_PCGP_full_vanilla_tCorr/`). The paper gives no MAP
  values; set 0 is the sample of each chain with the largest likelihood
  (the prior is flat).

## The K factor

The fits with a free K factor scale the vector-meson cross sections by K,
which is not an IP-Glasma parameter. `posterior_jpsi_<arXiv>_Kfactor_K.csv`
gives it for each set of `posterior_jpsi_<arXiv>_Kfactor.csv`, in the same
order: `posteriorSet` (as in `usedParameters<event>.dat`), `sample` (the row
of the sample in `chain.pkl`, counted from 0; −1 for a MAP set that is no
sample of the chain) and `K`. It is not an input file.

## How the J/ψ files of arXiv:2507.14087 and arXiv:2606.20362 are made

`make_jpsi_posterior_tables.py` makes them from the tar files of the Zenodo
records:

```bash
cd tables
python3 make_jpsi_posterior_tables.py 2507.14087 path/to/posterior_chain_full_fit.tar.gz
python3 make_jpsi_posterior_tables.py 2606.20362 path/to/analysis_and_data_files.tar.gz
```

- **Chains:** each `chain.pkl` (written by pocoMC) has weighted samples of
  `m`, `BG`, `BGq`, `Nq`, `smearingWidth`, `QsMuRatio`, `jimwlkMass`,
  `jimwlkLambdaQCD`, `UVDamp`, `omega` and K, with their weights and
  log-likelihoods. The script reads them with an unpickler that only allows
  numpy arrays.
- **Sets:** 1000 samples drawn with probabilities proportional to their
  weights by systematic resampling (the first uniform number of numpy's
  `default_rng(seed)`, with the seeds 250714087 and 260620362, as the offset,
  then the points (offset + i)/1000 on the cumulative weights). No weight
  exceeds 1/1000 of their sum, so no sample is drawn twice, and the sets keep
  the order of the chain. The medians and the 90% intervals of the sets agree
  with the weighted ones of the chains within the sampling error of 1000
  sets.
- **Fixed parameters:** `Nq`, `omega` and `UVDamp` (and K in the fits with
  K = 1) were fixed in the fits, so their chain values are draws of the prior;
  the files have the fixed values instead (3, 1, 0). The settings of the
  fits that a posterior file can set (`jimwlkAlphaS` 0, `jimwlkMu0` 0.28 GeV,
  `jimwlkC` 0.2, `BGqVar` and `NqFluc` 0) are constant columns.
- **Values** are written to 6 digits.

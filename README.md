# nma-motivated-trial

Code and data for **"Making the most of existing evidence: valid analysis of a new trial motivated by a pairwise or network meta-analysis"** (McKeen, Wang and O'Connor; submitted to *Frontiers in Pharmacology*).

When a new two-arm trial is conducted *because* an existing pairwise or network meta-analysis (NMA) gave a "promising" result, analyzing the new trial together with that evidence inflates the type one error. This repository implements two corrections:

- **Method 1** (selection rule known): conditions on the selection event.
- **Method 2** (selection rule unknown): conditions on the observed result of the existing meta-analysis.

Both are provided as the resampling procedures described in the paper and as exact closed forms (Section 2.3.3 of the paper). A pairwise meta-analysis is the special case of a network with two treatments: the closed forms only need the existing z-statistic and its variance. The repository also includes the theoretical quantities in Section 2.3.4 (critical values of Method 1, the type one error of the unadjusted analysis, i.e. the new trial added to the meta-analysis with the selection ignored, `naive_type1()`, and conditional power and sample size for planning the new trial) and the survey of trials published in 2024 that cite meta-analyses (Sections 2.1 and 3.1).

The R code uses base R only (no packages needed); the survey scripts use Python 3. Run R scripts from the repository root.

## Contents

| File | Purpose | Manuscript |
|---|---|---|
| `R/nma_fe.R` | Fixed-effect contrast-based NMA for arm-level binary data | |
| `R/closed_form.R` | `nma_update_test()` (closed-form p-values), `crit_value()` (Method 1 critical value), `crit_value_m2()`, `naive_type1()`, `cond_power()` (conditional power, Method 2 exact), `n_for_power()` (sample size) | Sections 2.3.3–2.3.4 |
| `R/resampling.R` | Resampling versions of Methods 1 and 2 | Sections 2.3.1–2.3.2 |
| `data/BRD_Arm_Network.csv` | Bovine respiratory disease antibiotic network (O'Connor et al. 2013, *Prev Vet Med* 110:77–87; updated O'Connor et al. 2016, *Prev Vet Med* 132:130–139) | |
| `survey/` | Survey of randomized trials that cite meta-analyses: searches, samples, codebook, codes, analysis (see `survey/README.md`) | Table 1 |
| `analysis/02_simulation_closed_form.R` → `results/sim_nma_*_closed_form.csv` | Network simulation with the closed forms (about 3 min) | Tables 4–5 |
| `analysis/05_pairwise_special_case.R` → `results/sim_pairwise_*.csv` | Pairwise meta-analysis simulation, gamithromycin vs florfenicol (about 1 min) | Table 6 |
| `analysis/01_example.R` → `results/example_cefth_trim.csv` | Illustrative example (about 1.5 min) | Table 7 |
| `analysis/03_theory_tables.R` → `results/critical_values.csv`, `results/sample_size.csv` | Critical values, type one error of the unadjusted analysis, sample sizes (seconds) | Tables 8–9 |
| `analysis/04_check_resampling.R` → `results/check_resampling_vs_closed_form.csv` | Resampling vs closed form on simulated data sets (about 9 min) | Section 3.2 |

Tables 2–3 of the manuscript (resampling simulation) are from McKeen's thesis code (see below).

## Quick use

```r
source("R/nma_fe.R"); source("R/closed_form.R")
d   <- read_brd("data/BRD_Arm_Network.csv")
old <- nma_fe(d)                                # existing network: CEFTH vs TRIM
y_new <- log(.4/.6) - log(.72/.28)              # new trial: 10/25 vs 18/25 affected
v_new <- 1/(25*.4*.6) + 1/(25*.72*.28)
nma_update_test(old$z, old$se^2, y_new, v_new, criterion = "lt")   # rule: existing p < .1
# pairwise meta-analysis: pass the pooled z-statistic and its variance as the first two arguments
```

## Relation to github.com/lamckeen/Error-Control-NMA

That repository holds an earlier version of the code from McKeen's thesis work; this repository is the version that accompanies the paper. Differences:

1. **Method 1 selection rule.** As in the paper (Section 2.3.1, step 6) and in the later version of the resampling code in McKeen's notes, a resampled data set is kept when the A–B p-value of the *existing* network (the resampled network with the new trial removed) meets the selection rule. The scripts in the earlier repository applied the rule to the network including the new trial.
2. The pooled null risk for treatments A and B is computed from the data rather than from hard-coded arm totals.
3. Closed forms, theoretical quantities, planning functions, the pairwise case and the survey are added.
4. Only base R is used for the statistical code.

## License

MIT

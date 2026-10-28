# nma-trial-selection

Code for **"Methods for controlling type one error when a new study is motivated by the results of an existing network meta-analysis"** (McKeen, Wang and O'Connor).

When a new two-arm trial is conducted *because* an existing network meta-analysis (NMA) gave a "promising" result, analyzing the new trial together with that network inflates the type one error. This repository implements two corrections:

- **Method 1** (selection rule known): conditions on the selection event.
- **Method 2** (selection rule unknown): conditions on the observed result of the existing network.

Both are provided as the resampling procedures described in the paper and as exact closed forms (Section 3.4.3 of the paper). It also includes the theoretical quantities in Section 3.4.4: critical values of Method 1, the type one error of the naive analysis, and conditional power and sample size for planning the new trial.

Everything uses base R (no packages needed). Run scripts from the repository root.

## Contents

| File | Purpose |
|---|---|
| `R/nma_fe.R` | Fixed-effect contrast-based NMA for arm-level binary data |
| `R/closed_form.R` | `nma_update_test()` (closed-form p-values), `crit_value()` (Method 1 critical value), `crit_value_m2()`, `naive_type1()`, `cond_power()` (conditional power, Method 2 exact), `n_for_power()` (sample size) |
| `R/resampling.R` | Resampling versions of Methods 1 and 2 |
| `data/BRD_Arm_Network.csv` | Bovine respiratory disease antibiotic network (O'Connor et al. 2013, *Prev Vet Med* 110:77–87) |
| `analysis/01_example.R` | Illustrative example, Table 3.5 (about 1.5 min) |
| `analysis/02_simulation_closed_form.R` | Simulation with the closed forms, Tables 3.3–3.4 (about 3 min) |
| `analysis/03_theory_tables.R` | Critical values, naive type one error and sample sizes, Tables 3.6–3.7 (seconds) |
| `analysis/04_check_resampling.R` | Resampling vs closed form on simulated data sets (about 9 min) |
| `results/` | Outputs of the scripts |

## Quick use

```r
source("R/nma_fe.R"); source("R/closed_form.R")
d   <- read_brd("data/BRD_Arm_Network.csv")
old <- nma_fe(d)                                # existing network: CEFTH vs TRIM
y_new <- log(.4/.6) - log(.72/.28)              # new trial: 10/25 vs 18/25 affected
v_new <- 1/(25*.4*.6) + 1/(25*.72*.28)
nma_update_test(old$z, old$se^2, y_new, v_new, criterion = "lt")   # rule: existing p < .1
```

## Relation to github.com/lamckeen/Error-Control-NMA

This repository supersedes the original code for the paper. Differences:

1. **Method 1 selection rule.** The paper (Section 3.4.1, step 6) keeps a resampled data set when the A–B p-value of the *existing* network meets the selection rule. The original scripts applied the rule to the network *including* the new trial. Under the rule ".05 < p < .1" that version does not control the type one error. Here the rule is applied to the existing network, as in the paper.
2. The pooled null risk for treatments A and B is computed from the data rather than from hard-coded arm totals.
3. Closed forms, theoretical quantities and planning functions are added.
4. Only base R is used.

## License

MIT

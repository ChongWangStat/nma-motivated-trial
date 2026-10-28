# Codebook: how randomized trials published in 2024 use earlier meta-analyses

Each record is one trial report (PMCID) that cites at least one meta-analysis (pairwise or network). You see:
- `title`, `abstract`
- `ma_refs`: the cited references whose text mentions a meta-analysis (id -> reference text). `nma_ref_ids` lists those that are network meta-analyses (NMA).
- `contexts`: passages of the full text where a meta-analysis reference is cited (sentence before and after included), with the section heading. Passages mentioning "meta-analysis" without a citation may also appear.
- `combine_flags`: sentences anywhere in the text that might describe pooling the trial's results with earlier studies.
- `intro_last_paragraph`: the last paragraph of the Introduction (usually the aims).

Code every record with the fields below. Judge only from what is shown; do not guess. When in doubt between two levels, choose the lower one and say why in `note`.

## Fields

1. `report_type`: one of
   - `primary` — reports the main or secondary outcomes of a randomized trial (including pilot or feasibility RCTs and cluster or crossover RCTs)
   - `secondary_analysis` — a post hoc, secondary, ancillary, pooled or subgroup analysis of an earlier trial's data, or a long-term follow-up reported separately
   - `protocol` — a trial protocol or design paper without results
   - `not_rct` — not a randomized trial (for example, an observational study or a systematic review that was mis-indexed)

2. `related` (`yes`/`no`): at least one cited meta-analysis addresses the same condition or population as the trial **and** includes the trial's experimental intervention, its comparator, or the same class of intervention. A meta-analysis cited only for a method, for prevalence or for a different intervention is `no`.

3. `related_nma` (`yes`/`no`): as `related`, but the related meta-analysis is a network meta-analysis (one of `nma_ref_ids`).

4. `motivation` — the role the cited meta-analyses play in why the trial was done:
   - `2` = **motivating**. A meta-analysis result is given as a reason for conducting this trial or for choosing its comparison. Examples: a meta-analysis suggested a benefit (or a promising, inconclusive or uncertain effect) that this trial sets out to confirm; it identified a gap (for example no head-to-head trials, low-certainty evidence, too few or too small trials) that this trial fills; the aims are framed in response to its findings.
   - `1` = **background only**. A meta-analysis is cited in the Introduction or rationale for background facts (efficacy of related treatments, burden of disease, mechanisms) but is not presented as a reason for this particular trial.
   - `0` = no meta-analysis is cited in the Introduction or rationale; meta-analyses appear only in Methods or Discussion.

5. `promising_signal` (`yes`/`no`): only if `motivation` = 2 — the meta-analysis result that motivated the trial was a favourable, promising or "trend" result for the intervention (as opposed to only a gap or absence of evidence).

6. `design_use` (`yes`/`no`): a meta-analysis result was used in the design, for example the assumed effect size, control-group risk or standard deviation for the sample size, or an explicit choice of dose, comparator or outcome taken from it.

7. `contextualised` (`yes`/`no`): the Discussion compares this trial's results with the results of a cited meta-analysis (agreement or disagreement).

8. `combined` (`yes`/`no`): the report statistically combines this trial's results with the earlier meta-analysis or its trials — an updated or cumulative meta-analysis that includes this trial, or a Bayesian analysis using the meta-analysis as a prior. Merely comparing numbers is `no`.

9. `motivating_ma_type`: `nma`, `pairwise`, `both` or `none` — which kind of meta-analysis played the motivating role (use `none` unless `motivation` = 2).

10. `evidence`: a short quotation (at most 30 words) from the passages that supports the `motivation` code (and `combined` if `yes`).

11. `note`: one short sentence on any uncertainty.

## Output

Return a JSON array, one object per record, with keys `pmcid`, `report_type`, `related`, `related_nma`, `motivation`, `promising_signal`, `design_use`, `contextualised`, `combined`, `motivating_ma_type`, `evidence`, `note`. Write it to the output file named in your instructions.

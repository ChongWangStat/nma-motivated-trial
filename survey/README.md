# Survey of randomized trials that cite meta-analyses (manuscript Sections 2.1 and 3.1)

- `search_and_sample.py` — Europe PMC queries (run 28 October 2026) and the seeded random samples.
- `sample_ids.json` — sampled PMCIDs (200 citing a network meta-analysis, 100 citing only pairwise meta-analyses), the 175 + 81 with open-access full text that were coded, and the 86 coded twice.
- `extract_contexts.py` — extracts, from each JATS full text, every passage citing a meta-analysis (with section), the abstract, the last paragraph of the Introduction and sentences that may describe pooled analyses.
- `codebook.md` — the coding instructions.
- `codes/A_batch_*.json` — first coding of all 256 reports; `codes/B_rel_*.json` — independent second coding of 86 reports.
- `analyse_codes.py` — merges the codings, records the adjudication of the 10 disagreements, computes agreement (kappa) and the proportions with Wilson 95% confidence intervals; writes `survey_codes_final.csv`.

Coding was done with a large language model (Claude Opus 5.5) applying `codebook.md`; each code is accompanied by a short supporting quotation (`evidence`). Full texts are not redistributed; they can be downloaded from Europe PMC by PMCID.

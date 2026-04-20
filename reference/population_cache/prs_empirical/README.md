# Empirical PRS reference distributions

Per-PGS empirical mean + SD + per-sample scores, computed by applying the PGS to
503 EUR samples in 1000 Genomes Phase 3 (CEU + FIN + GBR + IBS + TSI).

`run_prs.py` uses these to produce calibrated z-scores, replacing the
theoretical independence-formula approximation that systematically underestimates
population SD for LDpred2 / PRS-CS / PRS-CSx style dense scores.

## File naming

Two kinds of empirical files live here. `run_prs.py` prefers subject-specific
when present; otherwise it falls back to the full-panel file.

- **`<PGS_id>.json`** — **full-panel** empirical. 1000G EUR samples scored on
  every variant in the PGS weight file. Correct when the subject has ~100%
  coverage of the PGS panel (our imputed parquets usually reach 99–100%).
- **`<PGS_id>.<subject_id>.json`** — **subject-observed** empirical. 1000G EUR
  samples scored on only the subset of weights that `<subject>` actually has
  in their parquet. Necessary when subject coverage is materially below 100%
  — otherwise the subject and 1000G references accumulate contributions over
  different variant sets, producing a systematic bias.

## When to prefer subject-observed

If `run_prs.py` reports **coverage < ~98%** and the full-panel empirical SD is
tight (e.g., `SD ~ 1` on a raw score scale), the missing-variant systematic
bias can easily push the subject z-score several SDs outside the empirical
range even for a normal individual. Symptoms:

- z-score far outside [-3, +3] despite the subject being unremarkable
- `summarize_prs.py` flagging the row as **SUSPECT** due to |z| > 5
- Subject raw score sitting below the full-panel `min` field

Remediation: run

```bash
python scripts/calibrate_prs_empirical.py <PGS_id> --subject-observed <subject_id>
```

which produces `<PGS_id>.<subject_id>.json`, and `run_prs.py` will auto-prefer
it on the next run. The general shape: for a PGS with ~10,000 variants where
the subject is missing ~8% (so the 1000G samples score over the full panel
but the subject scores over only ~92%), the systematic bias from the missing
variants can push the subject's full-panel z-score several SDs outside the
empirical distribution in either direction, purely because the missing
variants have a non-zero average weight contribution. Restricting the 1000G
calibration to the subject's observable subset eliminates that bias — the
corrected z lands well inside the plausible range.

## Files here are derived

Everything in this directory is reproducible from the PGS weight files and the
1000G Phase 3 EBI release VCFs. Safe to delete any file and re-run
`calibrate_prs_empirical.py`. Each file takes 5–15 minutes depending on panel
size and whether `--subject-observed` subsetting is applied.

## Schema

Each JSON file:

```json
{
  "pgs_id": "PGS002804",
  "reference": "1000 Genomes Phase 3 EUR (CEU + FIN + GBR + IBS + TSI, N=503)",
  "n_samples": 503,
  "mean": 177.23,
  "sd": 12.15,
  "min": 146.8,
  "max": 213.4,
  "per_sample": {"HG00096": 184.2, "HG00097": 169.7, ...}
}
```

`per_sample` is useful for future re-analysis (e.g., re-deriving distribution
quantiles without re-running the whole calibration).

# NEXT_STEPS.md — limitations, open threads, infrastructure roadmap

Living document. Update as items complete or new gaps surface. Referenced from
`CLAUDE.md` so a fresh session picks this up at the start.

**This file is subject-anonymous** (repo-level). Subject-specific open items
and completion history go in `local/NEXT_STEPS.local.md`, which is gitignored.
See `local/README.md` for the public-vs-local split.

This file is for **project state** that isn't a per-finding ledger entry and
isn't an always-on convention:
- Known limitations of the current tooling (with their impact)
- Unstarted investigation threads of general interest
- Longer-term infrastructure work whose absence blocks or degrades current analyses
- Open questions needing user input before proceeding

## Known limitations

### PRS / polygenic scores

- **Chip coverage for GWAS-scale PRS panels.** ✅ Resolved by imputation (local Beagle 5.4 + 1000G Phase 3 EUR). Imputed parquets typically reach 99–100% panel coverage on published GWAS-scale PRS. The pre-imputation ~15% coverage is no longer the bottleneck; r² ceiling of each published PRS is now the dominant source of individual-prediction uncertainty.
- **Reference distribution calibration.** ✅ Resolved by empirical calibration against 503 1000G EUR samples via `scripts/calibrate_prs_empirical.py`. Theoretical independence-based calibration is retained only as a fallback when an empirical distribution isn't cached.
- **Sub-100% coverage systematic bias.** ✅ Resolved via `calibrate_prs_empirical.py --subject-observed <id>` producing per-subject empirical distributions (`<PGS>.<subject>.json`). Required when subject coverage on a PGS is materially below 100%, otherwise the 1000G-full-panel empirical and the subject's partial-panel score aren't apples-to-apples.
- **Sex inference is from chromosome call rate**, not explicit declaration. Currently robust but fragile for edge cases (sex chromosome aneuploidies, low Y-coverage females genotyped as male). Fix: allow explicit `sex` field in `profiles/<id>.json` to override.

### Array platform (23andMe v5 specifically)

- **CYP2D6 cannot be reliably called from array data.** Copy-number variants (*2xN duplications), hybrid genes (*4.013, *68), and the nearby pseudogene CYP2D7 aren't tractable on a SNP chip. Impact: codeine metabolism, tamoxifen, some antidepressants can't be tiered. Fix: clinical CYP2D6 sequencing if ever relevant to a specific drug decision.
- **HLA-B*57:01 only via tag SNP rs2395029.** Not clinically usable for abacavir decisions. Fix: clinical HLA typing if abacavir ever on the table.
- **Prothrombin G20210A (F2 rs1799963) not on chip.** Second-most-common European thrombophilia variant; unobserved state is not reassurance. Fix: order targeted test if clinical context warrants (unexplained VTE, family history).
- **Rh factor inference is gene-presence-based**, not a clinical test. Clean RHD-region probe coverage strongly suggests at least one functional RHD copy (Rh+), but transfusion medicine requires serology. Fix: accept medical-record Rh status as authoritative if available.

### Carrier screening

- **CFTR panel on chip covers ~26 common pathogenic variants** of 100+ clinical-panel variants. Fix: clinical CFTR sequencing before reproductive decisions (standard-of-care anyway).
- **SMN1 (SMA) cannot be called from array data.** SMN1 shares >99% sequence identity with paralog SMN2; dosage requires MLPA, qPCR, or targeted long-read. `scripts/run_carrier_panel.py` writes an explicit NOT-CALLABLE finding rather than attempting a call.
- **HEXA, PAH carrier panels currently use founder-population variants** (AJ-enriched HEXA, Euro-common PAH) that are mostly absent from the 23andMe v5 chip. Fix: expand panels to imputed-parquet-coverable variants so HEXA/PAH can produce actual calls rather than NOT CALLABLE.

### Ancestry-match metadata

- **Findings with `subject_ancestry_match: unknown`.** ✅ Resolved by extending `reclassify.py` to append superseders on metadata-only changes (ancestry_match, ancestry_downgrade, tier_rule_version), not only on tier movements.

### Presentation / reporting

- **Report generator** ✅ `scripts/generate_report.py` writes `reports/YYYY-MM-DD-<subject>.md` using the active ledger and Rules 1–11.

### Ledger / supersede model

- **Supersede is single-link; tombstone accumulation.** ✅ Resolved at the schema level. Findings carry `supersede_chain_root` (earliest row in a chain, computed at write time) and `is_tombstone` (bool). `ledger_io.load_active_findings()` is the canonical reader — groups rows by chain_root, picks the latest-timestamp row per chain, skips chains whose head is a tombstone. Consumers must prefer this helper over any manual `superseded_ids = {...}` filter. One-time migration: `scripts/migrate_ledger_schema.py` backfills existing rows with a timestamped backup.

## Unstarted investigations (general-interest threads)

Ordered roughly by broad relevance. Subject-specific completed items live in
`local/NEXT_STEPS.local.md`.

1. **Personality (Big Five) PRS.** Low r² per trait; useful for curiosity, not individual prediction. Opt-in only.
2. **Income PRS** (Hill 2019). Very weak individual signal; heavy environmental confounding. Opt-in only.
3. **Additional carrier panel variants** (HEXA and PAH coverage expansion; see Known limitations above).
4. **Couple / offspring analyses.** Two-subject mechanics (kinship, combined carrier risk, offspring PRS prediction) are defined in `SKILL.md` but not yet exercised end-to-end. First real run will likely surface debugging.

## Infrastructure roadmap

High-leverage items whose absence limits current work:

- **Genotype imputation pipeline.** ✅ Landed (local Beagle 5.4 + 1000G Phase 3 EUR). See `scripts/run_imputation.py`. Remaining optional: Minimac4 as an alternative engine if Beagle's memory pressure becomes an issue on larger chromosomes (not currently a problem).
- **LD-aware PRS reference distribution.** ✅ Resolved. `scripts/calibrate_prs_empirical.py`.
- **Report generator.** ✅ Landed. `scripts/generate_report.py`.
- **`reclassify.py` v2.** ✅ Resolved. Metadata-refresh handling + tombstone skip.
- **`tier_rules.py` — PRS-specific tier handling.** Currently PRS findings use `evidence_class: well_replicated_common_variant` + `inference_confidence` — works but isn't specifically aware of PRS r², coverage fraction, or tail-regression limitations. Could promote to its own `prs` evidence class with dedicated tier logic.
- **Second-subject ingestion + couple analysis mechanics.** Schema supports two subjects; runners are all subject-parameterized. Not yet exercised on a real second subject — expect minor first-run debugging when it lands.

## Open questions for the user

- **Reproductive-planning carrier panel scope:** just CFTR expanded, or the full ACMG panel (CFTR + SMN1 + HEXA + HBB + GJB2 + others)?
- **Ethical / sensitivity filters:** any trait category to deprioritize or skip when running a PRS sweep (e.g., income, personality)?

✅ **Imputation hosting preference.** Resolved: local Beagle 5.4 + 1000G Phase 3 EUR. Data stays local.

## Space management / portability

### What's on disk

- `reference/imputation/1kg_ref_b37/` — **~16 GB**. 1000 Genomes Phase 3 EUR reference panels (bref3 + VCF + filtered VCF). Needed only during imputation; not at PRS-score time.
- `reference/imputation/1kg_ebi_release/` — **~25 GB**. EBI release VCFs for empirical PRS calibration. Needed only during calibration; not for scoring.
- `reference/imputation/jdk/` — ~200 MB portable Temurin JDK 21.
- `reference/imputation/genetic_maps/` — ~25 MB PLINK recombination maps.
- `reference/imputation/beagle/` — <500 KB Beagle + conform-gt JARs.
- `reference/haplogroups/mtdna/haplogrep3/` — ~50 MB HaploGrep3 + JAR.
- `standardized-genomes/imputed/<subject>/beagle/*.vcf.gz` — per-chromosome imputed VCFs (~few GB/subject). Source of truth before merging to parquet.
- `standardized-genomes/imputed/<subject>/vcf/*` — intermediate aligned VCFs (~300 MB). Regeneratable any time from the chip parquet.
- `standardized-genomes/imputed/<subject>.imputed.parquet` — final merged imputed data. **This is what PRS analyses read.**

### Safe-to-delete once imputation is done for all subjects you care about

- `reference/imputation/1kg_ref_b37/` — biggest win (~16 GB). Re-fetch via `scripts/imputation_download.py --chr all` and re-filter via `scripts/filter_ref_vcf.py`.
- `reference/imputation/1kg_ebi_release/` — second-biggest (~25 GB). Re-fetch via `scripts/download_1kg_canonical.py` only if re-running empirical calibrations.
- `reference/imputation/jdk/` — re-install via `scripts/install_portable_jdk.py`.
- `standardized-genomes/imputed/<subject>/vcf/` — aligned intermediates; regenerate via `scripts/parquet_to_vcf.py` + `run_imputation.py` conform-gt step.
- `standardized-genomes/imputed/<subject>/beagle/*.vcf.gz` — keep unless space-starved.

### Keep

- `standardized-genomes/imputed/<subject>.imputed.parquet` — downstream analyses depend on this.
- `reference/population_cache/1kg_eur_afs.parquet` — canonical per-locus EUR AF table (~440 MB); re-derivable via `scripts/extract_eur_afs.py` but slow to rebuild.
- `reference/population_cache/prs_empirical/*.json` — empirical PRS distributions; small, expensive to regenerate.
- `reference/prs_weights/*.gz` — PGS Catalog files; small, easy to re-fetch but convenient to cache.

### Git

`.gitignore` excludes all personal data (raw files, standardized parquets, profiles, ledger, reports, imputed data, subject-specific empirical calibrations, per-subject haplogroup outputs, and the `local/` notes directory) plus all large re-fetchable infrastructure (reference panels, EBI VCFs, JDK, genetic maps, JARs, EUR AF parquet, PGS weight files). A fresh clone gives a collaborator code + conventions + curated reference tables. See `README.md` → "What's tracked vs. gitignored" for the authoritative list.

## How to use this file

This is **living state**, not a backlog. Stale items are worse than absent ones
because they produce false confident recommendations.

- **Session start:** scan this file after reading `profiles/` and `ledger/`. Also scan `local/NEXT_STEPS.local.md` for subject-specific items. Cross-check each item against current ledger + `scripts/` + `standardized-genomes/` state before citing it.
- **On completion of a subject-specific item:** append a resolved-line to `local/NEXT_STEPS.local.md`, not this file. This file captures repo-level trajectory only.
- **On completion of a repo-level item** (infrastructure, new convention, new analysis capability that wasn't previously possible): update this file.
- **Before recommending a next step to the user:** verify it's still true — the "missing" script truly doesn't exist, the PRS hasn't already run, etc.
- **On new limitation discovery:** add under the right section with a one-line impact and a one-line fix path.
- **On new thread request:** add under "Unstarted investigations" with approximate effort.

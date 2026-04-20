---
name: genome-eval
description: Analyze personal direct-to-consumer genotype data (23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA) for one or two subjects. Build a persistent ledger of findings with full evidence metrics, propose investigation threads, and generate reports on demand. Supports pharmacogenomics, disease risk, traits, ancestry, and couple/offspring analyses. Prioritizes evidence-graded, reclassifiable findings over one-shot interpretation.
---

# genome-eval

Personal genomics workflow for the `C:\projects\genome-eval\` setup. Designed for two subjects (self + partner), multi-provider, persistent ledger, investigation-driven (not one-shot).

**Companion file:** `CLAUDE.md` at the project root carries always-on stance and structural rules (no refusals, normalize-on-ingest, append-only ledger, directory invariants). This skill assumes those rules and focuses on *how* to execute the workflow. If SKILL.md and CLAUDE.md disagree, CLAUDE.md wins — it's the project-level authority.

## Operating principles

1. **Collect continuously, report on demand.** Accumulate findings into a ledger; only produce reports when the user asks for one on a specific topic.
2. **Store metrics, not verdicts.** Every finding records N, p-value, effect size, replication count, cohort ancestry, source. Tier labels are derived from these metrics via `tier_rules.py` and can be recomputed if rules change.
3. **Propose, don't surprise.** At session start, summarize ledger state and propose 3–5 unexplored threads with effort estimates. Execute only what the user picks. Never initiate high-effort work (population API pulls, imputation, large downloads) without explicit approval.
4. **Evidence-faithful reporting.** When two sources disagree or overlap, keep them as separate findings and report both with their metrics side-by-side. Do not silently pick a winner.
5. **Be rigorous about what the array cannot tell us.** Arrays don't see indels, CNVs, repeat expansions, or rare variants absent from the chip. "Not on chip" and "genotyped as wild-type" are different claims and must be recorded differently.
6. **Normalize on ingest.** Every raw file, from any provider known or unknown, is parsed into `standardized-genomes/<id>.parquet` with the standard schema before any analysis touches it. Analysis code reads only normalized data. New provider = write one parser, not rewrite anything downstream.
7. **Enable, don't gatekeep.** If data is loaded into the system, it's legitimate — the user's act of loading is the assertion of permission. The skill analyzes; it does not audit. Caveats are epistemic ("this PRS explains 12% of variance in European samples, drops substantially outside that"), never moral. No analysis is refused on grounds of sensitivity — intelligence PRS, behavioral traits, kinship, paternity, forensic-style identity work, external-database matching — all fine when the user asks. Honest metrics go with every finding; the user decides what to do with them.

## Directory layout

```
C:\projects\genome-eval\
├── SKILL.md                        # this file
│
├── raw-source-genomes/             # provider files as received, never edited
│   ├── alice-23andme-example.txt
│   └── wife-ancestrydna.txt
│
├── standardized-genomes/           # normalized parquet (analysis reads ONLY from here)
│   ├── alice.parquet
│   └── wife.parquet
│
├── profiles/                       # per-subject metadata (flat, one file per subject)
│   ├── alice.json
│   └── wife.json
│
├── ledger/
│   ├── findings.jsonl              # one finding per line
│   ├── sources.jsonl               # one source per line
│   └── investigations.jsonl        # one investigation per line
│
├── reports/
│   └── YYYY-MM-DD-<topic>.md       # generated on request
│
├── reference/
│   ├── curated_snps.tsv            # tier 1/2 SNP table (see end of this file)
│   ├── prs_weights/                # GWAS summary statistics for polygenic scores
│   ├── ld_panels/                  # LD reference (e.g., 1000 Genomes subset) for proper PRS methods
│   └── population_cache/           # cached population allele frequencies, if fetched
│
└── scripts/
    ├── parse_23andme.py
    ├── parse_ancestry.py
    ├── parse_myheritage.py
    ├── parse_ftdna.py
    ├── normalize.py
    ├── tier_rules.py
    └── ledger_io.py
```

The subject identifier (e.g., `alice`, `bob`) is used consistently across `standardized-genomes/<id>.parquet` and `profiles/<id>.json`. Raw files keep descriptive names including provider and date.

**Hard separation between raw and standardized.** Analysis code and investigation logic read only from `standardized-genomes/` and `profiles/`. Nothing downstream of normalization should open a raw file. This means adding a new provider = write one parser that lands a parquet in `standardized-genomes/`, and everything else works unchanged.

## Initialization (first run or new subject)

Triggered when the user invokes this skill for the first time, when a new raw data file appears in `raw-source-genomes/` without a matching `standardized-genomes/` entry, or when the user explicitly says "load a new subject."

1. **Scan.** List files in `raw-source-genomes/` and `standardized-genomes/`. For standardized files, read the matching `profiles/<id>.json` if present.
2. **Ask for a display name** for any raw file not yet associated with a subject: "What ID should I use for the file `<filename>`? Short, filename-safe (e.g., `alice`, `bob`)." This ID drives all downstream filenames.
3. **Detect provider, chip version, build** from the file header (see "File formats" below).
4. **Parse and normalize** into `standardized-genomes/<id>.parquet`.
5. **Write** `profiles/<id>.json` with display name, provider, chip, build, file hash, parse stats (total SNPs, no-call rate, inferred sex from X/Y call pattern, presence of MT), and an empty `sharing_sensitivity` field (default `"unset"` — asked for only when the user requests an export, not at load time).
6. **Initialize ledger files** if missing. No findings yet.
7. **Summarize state** and propose initial investigation threads (see "Investigation workflow").

Optional subject-level metadata the user can add to `profiles/<id>.json` at any time: declared ancestry (informs PRS tier calibration), current medications (informs PGx priority ordering), known family history (informs disease-risk priority). These are helpful but never required; absent fields just mean less-informed thread proposals.

If two subjects are loaded, couple analyses become available and are proposed as threads automatically.

## File formats — provider parsing

All of the following use **GRCh37 / hg19**, 1-based coordinates, plus-strand alleles. Different column structures.

### 23andMe

TSV, `#` comments, 4 columns: `rsid, chromosome, position, genotype`. Genotype is concatenated 2-char string (`AG`, `CC`), `--` for no-call, single char for X/Y hemizygous. Chip version detectable from SNP count and presence/absence of signature rsIDs (v3 ≈ 960k, v4 ≈ 600k, v5 ≈ 640k with different content). File header includes `file_id`, `timestamp`, `signature` — capture these in metadata for provenance.

### AncestryDNA

TSV, `#` comments, 5 columns: `rsid, chromosome, position, allele1, allele2`. Alleles in separate columns. No-call is `0` (not `--`). Chromosome may use `23, 24, 25, 26` for `X, Y, XY (PAR), MT` — remap on parse. Uses GRCh37.

### MyHeritage

CSV, quoted fields, 4 columns: `"RSID","CHROMOSOME","POSITION","RESULT"`. Genotype in one 2-char field. Newer files (2022+) may use GRCh38 — **check header**; if build 38, liftover to 37 or track separately. Header contains "MyHeritage".

### FamilyTreeDNA (FamilyFinder)

CSV, similar to 23andMe but comma-separated. Columns: `RSID,CHROMOSOME,POSITION,RESULT`. Build 37. Genotype 2-char concatenated.

### Detection

Read first ~20 lines, check for provider-identifying strings in the header (`23andMe`, `AncestryDNA`, `MyHeritage`, `FamilyTreeDNA`), fall back to column structure / delimiter / header row detection. Record detected provider in metadata, along with confidence.

### Unknown provider handling

If the detector can't identify the provider (new service, old file, re-exported format), don't guess silently. Inspect the file:

1. Read the first 50 lines and the last 5 lines.
2. Identify: delimiter (tab / comma / other), column structure, whether there's a header row, how genotypes are encoded (concatenated 2-char, separate columns, quoted), how no-calls are represented (`--`, `0`, `00`, empty), and what build the file declares (search header for `GRCh37`, `GRCh38`, `hg19`, `hg38`, `build 37`, `build 38`).
3. Present findings to the user: "File appears to be [best guess]. Delimiter: X. Genotype encoding: Y. Build: Z (declared / inferred / unknown). Proceed with parser matching these assumptions?"
4. On confirmation, write a new parser following the existing `parse_*` pattern, add it to the `PARSERS` dict, and run it. Store parser source in `scripts/` with a comment noting which provider it targets.
5. If build is unclear, ask the user directly — they can usually check the provider's documentation. If the file is GRCh38, either liftover to 37 (requires chain file, medium effort) or track separately in metadata and adjust position-based lookups. Don't silently assume build 37.

### Normalization target

After parsing, write `standardized-genomes/<id>.parquet` with consistent schema regardless of source provider:

```
rsid: str
chrom: str              # '1'..'22', 'X', 'Y', 'MT'
pos: int                # GRCh37, 1-based
a1: str                 # single char, plus-strand; None if no-call
a2: str                 # single char; None for hemizygous or no-call
build: str              # 'GRCh37'
source_provider: str    # '23andMe', 'AncestryDNA', etc.
```

Alleles normalized to sorted tuple `(a1, a2)` with `a1 <= a2` for robust comparison (skip for sex chromosomes in males).

## Ledger schemas

All three ledger files are JSONL (one JSON object per line). Append-only. If a finding needs revision, append a new record with a pointer to the superseded one rather than editing in place.

### findings.jsonl

```json
{
  "finding_id": "uuid4",
  "timestamp": "2026-04-17T14:35:00Z",
  "subject_id": "alice",
  "topic": "caffeine_metabolism",
  "claim": "Fast caffeine metabolizer (CYP1A2)",
  "variants": [
    {"rsid": "rs762551", "chrom": "15", "pos": 75041917,
     "ref": "C", "alt": "A", "genotype": "AA"}
  ],
  "effect": {
    "type": "OR",
    "value": 1.8,
    "ci_low": 1.3,
    "ci_high": 2.5,
    "direction": "faster metabolism"
  },
  "study_n": 12000,
  "p_value": 3.2e-8,
  "replication_count": 3,
  "cohort_ancestry": "European",
  "subject_ancestry_match": "match",
  "source_ids": ["pmid:12345678"],
  "tier_computed": "B",
  "tier_computed_at": "2026-04-17T14:35:00Z",
  "tier_rule_version": "v1",
  "investigation_id": "uuid4",
  "supersedes": null,
  "notes": ""
}
```

Key rules:
- **One claim per record.** If two studies agree, that's two findings (both recorded, both with their metrics).
- **No-call is a valid finding** (`"claim": "not genotyped", "variants": [...with genotype: null]`) — useful to record that we looked and failed.
- **"Not on chip" is a different finding** (`"claim": "variant not present on <provider> <chip_version>"`).
- `tier_computed` is derived, not authoritative. Re-run `tier_rules.py` over all findings to update.

### sources.jsonl

```json
{
  "source_id": "pmid:12345678",
  "type": "peer-reviewed | pharmgkb | clinvar | cpic | gwas-catalog | snpedia | press | other",
  "citation": "Smith et al. 2020, Nat Genet 52:123–130",
  "url": "https://...",
  "accessed_at": "2026-04-17T14:30:00Z",
  "cohort_ancestry": "European",
  "cohort_n": 12000,
  "evidence_class": "primary GWAS | meta-analysis | review | guideline | community",
  "notes": ""
}
```

### investigations.jsonl

```json
{
  "investigation_id": "uuid4",
  "timestamp": "2026-04-17T14:20:00Z",
  "query": "Caffeine sensitivity",
  "subject_ids": ["alice"],
  "status": "proposed | in_progress | completed | abandoned",
  "effort_estimate": "low | medium | high",
  "effort_actual": "low",
  "initiated_by": "user | claude-proposal",
  "sources_consulted": ["pmid:12345678", "pharmgkb:PA166127543"],
  "findings_generated": ["finding_uuid_1", "finding_uuid_2"],
  "next_steps": [
    "Check CYP2A6 for secondary metabolism pathway",
    "Compare against wife's genotype once data available"
  ],
  "notes": ""
}
```

Effort estimate definitions:
- **low** — lookup in curated tables, single API call, <2 minutes
- **medium** — multi-source synthesis, 3–10 API calls, literature reading, 5–20 minutes
- **high** — bulk downloads, imputation, population API sweeps, >20 minutes or significant compute

## Evidence framework

### Metrics to capture (always, when available)

- `study_n` (cohort size)
- `p_value`
- `effect_size` with type (OR, β, HR, risk difference)
- `replication_count` (number of independent cohorts reproducing the finding)
- `cohort_ancestry` (e.g., "European", "East Asian", "multi-ethnic", "unknown")
- `subject_ancestry_match` (match / partial / mismatch / unknown — derived from subject's declared or inferred ancestry)

If any metric is unavailable, store `null` rather than inventing. Missing metrics drop the computed tier.

### Tier rules (`tier_rules.py`, v2)

Each finding is tiered by whichever pathway its metadata supports. The rules are checked top-down; first match wins.

| Tier | Criteria | Notes |
|---|---|---|
| **A** | CPIC Level A/B guideline | Clinical pharmacogenomics |
| **A** | ClinVar "Pathogenic/Likely Pathogenic" with ≥2-star review status | Variant-level clinical classification |
| **A** | `evidence_class ∈ {mendelian_trait, near_mendelian_trait}` with replication_count ≥ 5 | Classical single-gene traits (ABO, ABCC11 earwax, ALDH2, etc.) |
| **B** | Replicated (rep ≥ 2), combined N ≥ 10,000, p < 5e-8, ancestry-matched | GWAS-grade |
| **B** | `evidence_class = well_replicated_common_variant` with replication_count ≥ 10 | Common variants with 10+ independent replications (lactase, HERC2, TAS2R38, etc.) |
| **B** | `evidence_class ∈ {multi_snp_composite, gene_presence_inference}` with `inference_confidence = high` | High-confidence derived findings (eye color multi-SNP, Rh factor from gene-region coverage) |
| **C** | Single GWAS with N ≥ 5,000 and p < 5e-8 | |
| **C** | Multiple small studies (N 500–5000, rep ≥ 2) | |
| **C** | `evidence_class = well_replicated_common_variant` with replication_count 3–9 | |
| **C** | Inference findings with `inference_confidence = moderate` | |
| **D** | Nominal significance (p < 0.05) | Candidate-gene or heterogeneous meta-analysis |
| **D** | `evidence_class = weakly_predictive_variant` | Variants with small effect / inconsistent replication (e.g., CYP1A2 rs762551 caffeine) |
| **D** | Inference findings with `inference_confidence = low` | |
| **E** | SNPedia / press / community-class evidence | |
| **E** | Any study with N < 500 and marginal significance | |
| **unknown** | `evidence_class ∈ {not_callable_from_array, array_limitation, suspected_miscall}` | Honest "the chip can't answer this" |
| **unknown** | No classifying metadata at all | Falls through every rule |

**Ancestry mismatch downgrade:** any finding meeting higher-tier criteria but with cohort ancestry mismatched to subject drops one tier and gets a flag `ancestry_downgrade: true` in the record. Example: a European subject evaluated against an East Asian-cohort variant (rs671 ALDH2) computes A via the mendelian_trait rule but downgrades to B because the cohort ancestry doesn't match the subject.

**Evidence-class vocabulary (set by investigator at recording time):**

- `mendelian_trait` / `near_mendelian_trait` — classical single-gene traits with near-deterministic phenotype mapping
- `well_replicated_common_variant` — common variant with substantial GWAS / replication history but not Mendelian
- `weakly_predictive_variant` — real association but small effect / inconsistent replication
- `multi_snp_composite` — derived from multiple loci (set `inference_confidence`)
- `gene_presence_inference` — derived from probe-count / coverage patterns (set `inference_confidence`)
- `not_callable_from_array` / `array_limitation` — chip cannot reliably answer this question
- `suspected_miscall` — a specific variant call is believed incorrect based on cross-check evidence (LD-companion disagreement, etc.)

Set `inference_confidence` to `high` / `moderate` / `low` for multi-SNP and gene-presence findings. Leave unset for other evidence classes.

### Reclassification

`tier_rules.py` includes a version field. If rules change, bump the version and re-run over all findings. Re-derivation is cheap; no re-research needed. Record which rule version produced each tier.

### Tier discipline: evidence only, not clinical actionability

The tier label records **evidence quality for the claim the finding makes.** It does not record clinical actionability, guideline recommendations, regulatory status, social sensitivity, or anything else. If the evidence for the claim is real, compute the tier from the metrics — regardless of whether a clinical body recommends acting on it.

**Don't:**
- Leave a well-studied variant at `unknown` tier because ACMG (or any guideline body) advises against clinical action (e.g., MTHFR C677T). That's gatekeeping through the tier system.
- Omit `study_n`, `p_value`, or `replication_count` so the tier rule falls through to `unknown` — that buries a finding via the schema.
- Use `unknown` when what's actually meant is "clinically not actionable," "contested," or "controversial."
- Attach the ClinVar significance annotation only when the subject carries the variant. ClinVar evidence is a property of the locus, not of the subject's genotype — a negative result at a Pathogenic locus has the same evidence weight as a positive one, just a different claim.

**Do:**
- Record the real metrics for the specific claim the finding makes (biochemical effect, population frequency, penetrance, whatever the locus is about). Let the tier rule compute what it computes.
- Record the clinical-actionability position as a prominent note, clearly labeled (e.g., `ACTIONABILITY CAVEAT: ...`), carrying the guideline source. Evidence tier and clinical actionability are two separate axes — both belong in the record.
- If the evidence is genuinely thin, a low tier or `unknown` is fine — but the reason must be evidence, not policy.

**Why:** The project stance explicitly rejects gatekeeping. "Caveats are epistemic, never moral." A reader is capable of integrating "the enzyme runs at 65% of normal *and* ACMG says don't act on it in healthy people." Collapsing one of those into the other reduces the information the reader receives.

### ClinVar annotation rule

ClinVar significance and review stars are properties of the **locus and its evidence**, not of the subject's genotype. Attach them to every finding about a ClinVar-annotated locus, regardless of whether the subject carries the variant:

- Subject is a carrier → claim says "heterozygous / homozygous"; evidence annotation unchanged.
- Subject is wild-type at a Pathogenic locus → claim says "reference / non-carrier"; evidence annotation unchanged (still Pathogenic, same star count).
- Subject's locus was not genotyped → claim says "not genotyped"; evidence annotation still applies to the locus and is recorded; the note makes clear that subject status is *unknown*, not *absent*.

This keeps the evidence axis independent of the genotype-call axis, which is the whole point of storing metrics rather than verdicts.

## Investigation workflow

### Session start

1. Read `profiles/*.json` — which subjects are loaded?
2. Read `ledger/*.jsonl` — summarize: N findings, latest investigation, topics covered.
3. Identify **unexplored threads**: topics in the curated SNP table with no corresponding finding yet; recent PharmGKB guideline updates; couple analyses (if second subject exists) that haven't been run.
4. Present to user:
   - 2-sentence summary of ledger state
   - 3–5 proposed investigation threads, each with: topic, effort estimate, what might be found, why it's worth doing
   - Ask which to pursue, or wait for a specific request

Example opener:

> Ledger has 18 findings for `alice` across pharmacogenomics and traits. Nothing yet on cardiovascular risk variants or ancestry markers. No second subject loaded.
>
> Unexplored threads:
> - **Cardiovascular risk panel** (low effort) — Factor V Leiden, prothrombin, APOE, LPA. All on-chip.
> - **HFE / iron overload** (low effort) — two SNPs, quick.
> - **Fluoropyrimidine (DPYD) PGx** (low effort) — high-value if you ever need 5-FU or capecitabine.
> - **Caffeine & alcohol metabolism panel** (low effort) — combines CYP1A2, ADH1B, ALDH2.
> - **Polygenic height prediction** (medium effort, ancestry-caveats apply) — requires external PRS weights.
>
> Which of these, or something else?

### During investigation

- Use curated tables first (`reference/curated_snps.tsv`). If the topic is covered by the curated set, one file lookup answers most of the question.
- For topics beyond the curated set: web search PharmGKB, ClinVar, GWAS Catalog, OMIM, specific papers. Prefer primary sources over SNPedia.
- Record every source consulted in `sources.jsonl` before referencing it in a finding.
- When two sources disagree (different effect sizes, different direction), create **two separate findings**, flag as `conflict_with: <other_finding_id>`.

### End of investigation

- Write investigation summary to `investigations.jsonl`.
- Ask the user: "Want a report now, or keep accumulating?"
- If report requested, write to `reports/YYYY-MM-DD-<topic>.md` (see "Reports").

## Couple analyses (both subjects loaded)

Enabled only when `profiles/` has two subjects with parsed data. Skill asks for consent from both where possible ("are both subjects aware their data is being jointly analyzed?") and records in metadata.

### Capabilities, by reliability

**High reliability** (single-locus Mendelian):
- ABO blood type — from rs8176719 (and optionally rs8176746, rs8176747). Clean inheritance.
- Earwax / axillary odor — rs17822931 (ABCC11). Fully Mendelian.
- Lactase persistence — rs4988235 in Europeans; different variants in other populations.
- RhD — not well-covered by arrays (RHD gene deletion is the common negative allele; arrays can't see deletions reliably).

**Medium reliability** (oligogenic, dominant locus + modifiers):
- Eye color — HERC2 rs12913832 is the dominant locus; ~75% of blue/brown variance. Secondary loci (OCA2, TYR, SLC24A4, SLC45A2) add resolution. Published models (IrisPlex, HIrisPlex-S) hit ~85% accuracy for blue vs brown in Europeans, degrade for green/hazel and in admixed populations.
- Hair color — HIrisPlex-S model, 11 SNPs. Accuracy: ~80% for red, ~70% for black, worse for brown/blonde distinction.
- Skin pigmentation — multi-locus; models exist but performance is ancestry-dependent.

**Low reliability** (highly polygenic):
- Height — ~10,000 contributing loci, array captures some. Can predict parent-midparent correlation but specific offspring prediction has wide CI.
- Educational attainment, personality, etc. — do not attempt from this data. PRS for these are weak and ethically fraught.

**Carrier-carrier screening** (high clinical value):
- For each recessive condition with pathogenic variants on both subjects' chips, compute whether compound het / homozygous affected is possible in offspring.
- Main conditions worth screening: CFTR (cystic fibrosis), HBB (sickle cell, beta-thalassemia), HEXA (Tay-Sachs, limited array coverage — caveat), SMN1 (SMA, not well-detected by arrays), GJB2 (connexin-26 deafness), plus founder-population-specific variants if relevant ancestry.
- ACMG recommends expanded carrier screening for ~100+ conditions in couples planning pregnancy; array data covers a fraction of these. State clearly what was and wasn't screened.

**Ancestry overlap & IBD** (if both subjects consent):
- Approximate IBD (identity by descent) — can estimate relatedness from shared haplotype segments. Requires both files in same build, phasing. Moderate effort.
- Shared ancestry inference — PCA against 1000 Genomes reference. Medium effort.
- Actual relatedness calculation uses KING or similar tools; mention as an option, don't auto-run.

### Offspring prediction output format

For each trait, produce a probability distribution, not a point estimate. Example:

```
Eye color prediction (based on HIrisPlex-S, 11 SNPs):
  Brown: 62%
  Blue:  24%
  Intermediate/green: 14%
  Model accuracy for European-ancestry parents: ~85% (Walsh 2017)
  Caveats: model trained on European cohort; performance degrades
  outside that ancestry. Non-HERC2 loci contribute substantial uncertainty.
```

Never give a single-label prediction without distribution + confidence.

## Population reference data (explicit-consent gated)

**This is opt-in.** Do not fetch population data without explicit user approval.

### When to propose it

- User asks "is this common?" or "what's typical for [ancestry]?"
- Second subject data not available; user wants couple analysis with a placeholder
- Ancestry inference needed for tier calibration

Prompt format: "This needs population allele frequencies for [ancestry X]. I can fetch from [source]. Effort: [low/medium/high]. Proceed?"

### Sources (in order of preference)

1. **Ensembl REST API** — per-SNP, free, fast. `GET https://rest.ensembl.org/variation/human/<rsid>?pops=1`. Returns allele freqs across 1000 Genomes populations. Use for one-off queries and targeted sets (<500 SNPs).
2. **gnomAD** — larger cohort, more granular populations. GraphQL API at `https://gnomad.broadinstitute.org/api`. Use for modern high-quality population freqs. Good for individual loci and bulk if paginated.
3. **ALFA (NCBI)** — Allele Frequency Aggregator. Use if Ensembl/gnomAD missing a variant.
4. **1000 Genomes VCF direct** — only for bulk work (whole-genome population profiles). High effort; download is substantial.

### "Placeholder partner" workflow (requires explicit request)

If user wants to run couple analysis before partner's data is available:
1. Ask partner's declared ancestry (population codes: EUR, AFR, EAS, SAS, AMR for 1000 Genomes).
2. Fetch allele frequencies for all curated SNPs from Ensembl (medium effort — ~60 API calls with rate limiting).
3. Cache to `reference/population_cache/<population>.tsv`.
4. For offspring prediction, simulate partner genotype as a random draw from the population; report distributions across Monte Carlo iterations (10,000+).
5. Clearly label all output as "simulated partner, ancestry-based" and contrast with what real data would give.

This is interesting-but-rough. Real partner data is dramatically more informative.

## Tier 1 — Pharmacogenomics (start here, curated)

All alleles given as **forward-strand / plus-strand**, matching what appears in raw data files directly.

### CYP2C19 — clopidogrel, SSRIs, PPIs, voriconazole

| Star allele | rsID | Forward ref>alt | Effect of alt |
|---|---|---|---|
| *2 | rs4244285 | G>A | Loss of function (splicing) |
| *3 | rs4986893 | G>A | Loss of function; mostly East Asian |
| *17 | rs12248560 | C>T | Increased expression (ultra-rapid) |

**Diplotype → phenotype:**
- *1/*1 → Normal Metabolizer
- *1/*17, *17/*17 → Rapid / Ultra-rapid
- *1/*2, *1/*3, *2/*17 → Intermediate
- *2/*2, *2/*3, *3/*3 → Poor

Clopidogrel: PM/IM → use alternative (prasugrel/ticagrelor). CPIC Level A.

**Phasing note:** array data is unphased. *2 + *17 combinations require care. Flag as "unphased; most probable diplotype."

### CYP2C9 — warfarin, phenytoin, NSAIDs

| Star | rsID | Ref>alt | Effect |
|---|---|---|---|
| *2 | rs1799853 | C>T | Decreased function |
| *3 | rs1057910 | A>C | Strongly decreased |

Combine with VKORC1 for warfarin (IWPC / Gage algorithm).

### VKORC1 — warfarin sensitivity

| rsID | Ref>alt | Effect |
|---|---|---|
| rs9923231 | C>T | T allele → increased sensitivity, lower dose |

**Strand note:** Literature says "−1639G>A" (gene coding strand). Forward strand is C>T. 23andMe shows C or T, not G or A.

### DPYD — 5-FU / capecitabine (toxicity severe)

| Star | rsID | Ref>alt | Effect |
|---|---|---|---|
| *2A | rs3918290 | C>T | No function (splice) |
| *13 | rs55886062 | A>C | No function |
| c.2846A>T | rs67376798 | T>A | Decreased function (forward-strand; c.2846A>T on gene-strand — DPYD is on chr1 minus strand) |
| HapB3 tag | rs75017182 | G>C | Decreased (forward-strand; gene-strand c.1129-5923C>G — flipped because DPYD is minus-strand). May be absent on v4 chip. |

Any decreased/no-function → dose reduction or alternative. CPIC Level A. Missed DPYD variants cause fatal toxicity; treat positives as high-priority.

### TPMT — thiopurines (azathioprine, 6-MP, thioguanine)

| Star | rsID | Ref>alt |
|---|---|---|
| *2 | rs1800462 | C>G |
| *3B | rs1800460 | C>T |
| *3C | rs1142345 | T>C |

*3A = *3B + *3C in cis. Het → dose reduction; hom → severe myelosuppression risk.

### NUDT15 — thiopurines (complements TPMT, higher impact in East/South Asian)

| rsID | Ref>alt | Effect |
|---|---|---|
| rs116855232 | C>T | T = decreased function |

### SLCO1B1 — statin myopathy (especially simvastatin)

| rsID | Ref>alt | Effect |
|---|---|---|
| rs4149056 | T>C | C = *5 allele; CC → meaningful myopathy risk |

### UGT1A1 — irinotecan, atazanavir

*28 is a TA repeat (arrays can't genotype repeats). Tag SNP:

| rsID | Ref>alt | Effect |
|---|---|---|
| rs887829 | C>T | T tags *28 (reduced glucuronidation) |

### HLA-B*57:01 — abacavir hypersensitivity

Arrays can't type HLA directly. rs2395029 is a reasonable tag in Europeans but unreliable across ancestries. **Never issue an abacavir-safe verdict from array data.** Recommend formal HLA typing.

### CYP2D6 — DO NOT CALL PHENOTYPE FROM ARRAY DATA

Copy number variation and CYP2D7 hybrids drive most CYP2D6 variability, and arrays can't detect them. Key SNPs (rs3892097 *4, rs5030655 *6, rs1065852 *10, rs35742686 *3) can be reported as raw genotypes but do NOT assign star alleles or phenotypes. Clinical PGx testing required.

## Tier 2 — Disease risk / carrier variants

| Gene / condition | rsID(s) | Ref>alt | Interpretation |
|---|---|---|---|
| **APOE** (Alzheimer's, CV) | rs429358 + rs7412 | See below | ε2/ε3/ε4 haplotypes |
| **HFE** C282Y (hemochromatosis) | rs1800562 | G>A | A/A → high iron overload risk |
| **HFE** H63D | rs1799945 | C>G | Lower penetrance; compound het with C282Y matters |
| **Factor V Leiden** | rs6025 | C>T | T → VTE risk (~4–8× het, higher hom) |
| **Prothrombin G20210A** | rs1799963 | G>A | A → VTE risk |
| **SERPINA1** PI*S | rs17580 | T>A | A = S deficiency |
| **SERPINA1** PI*Z | rs28929474 | C>T | T = Z deficiency (more severe) |
| **LPA** (lipoprotein(a), CV) | rs10455872, rs3798220 | varies | Elevated Lp(a) → CV risk |
| **MTHFR** C677T | rs1801133 | G>A | A = reduced enzyme; **ACMG: routine testing not recommended** |
| **MTHFR** A1298C | rs1801131 | T>G | Same caveat |

### APOE (special handling)

- **rs429358:** T = Cys112, C = Arg112
- **rs7412:** C = Arg158, T = Cys158

Haplotypes (forward strand):
- ε2: T at rs429358, T at rs7412
- ε3: T at rs429358, C at rs7412
- ε4: C at rs429358, C at rs7412

Array data is unphased. Most diplotypes are unambiguous from the two genotype pairs; one case (het at both) is formally ambiguous between ε1/ε3 and ε2/ε4. ε1 is vanishingly rare — call it ε2/ε4 and state the assumption.

**Treat APOE as sensitive.** ε4/ε4 is associated with substantially elevated Alzheimer's risk. Present findings neutrally, recommend genetic counseling, do not editorialize.

## Tier 3 — Traits (do only if asked or for couple analysis)

Each row lists the plain-English phrasing to use when surfacing results to the reader, alongside the gene/rsID mapping. **Never report a trait result using only the jargon label** (e.g., "lactase non-persistent", "PAV/AVI", "*1F/*1A"); always translate. See "Presentation conventions → Rule 7" for the full do/don't table.

Every row lists a **confidence level** for the genotype → phenotype mapping. Confidence is the inherent biological predictivity of the locus (how deterministic the relationship is in the general population); it's distinct from the per-subject genotype-call confidence (array error rate, probe quality) and must be surfaced alongside the result. See Rule 10 in "Presentation conventions."

| Trait (plain language) | Jargon form (do NOT use alone) | rsID(s) | Coding | Confidence (genotype → phenotype) |
|---|---|---|---|---|
| Eye color (dominant locus) | — | rs12913832 (HERC2) | G/G → brown; A/A → blue; A/G → mixed. European-calibrated. | **Low–moderate** alone: rs12913832 explains ~74% of blue/brown variance; single-SNP calls are also subject to array miscall (~0.1–0.5%). Always cross-check with linkage-companion SNPs rs1129038 and rs1667394 (tight LD in Europeans). |
| Eye color (multi-SNP composite) | HIrisPlex-S | 11 SNPs incl. HERC2, OCA2, SLC24A4, SLC45A2, TYR, IRF4, TYRP1 | Weighted composite across loci; direction-of-effect per allele published. | **Moderate–high**. Near-HIrisPlex-S accuracy (~85% for blue/brown in Europeans) when the full model is used. Lets miscalls at any single probe be detected via internal inconsistency. |
| Hair color | HIrisPlex-S | same 11 SNPs | Model output. | Moderate in Europeans; lower outside that ancestry. |
| Adult lactose tolerance | "lactase persistence / non-persistence" | rs4988235 (MCM6) | T allele → lactose-tolerant adult (European-derived persistence). C/C (ancestral) → reduced adult lactase activity. Different causal variants in African / Arabian populations. | **Moderate**. Genotype-call is high-confidence; penetrance is probabilistic — about 70% of G/G Europeans are symptomatic with typical dairy doses; ~30% are asymptomatic or mildly symptomatic due to microbiome adaptation, dose, product type, regular consumption. Never report C/C as "lactose intolerant" without the probabilistic framing. |
| Alcohol flush reaction (East Asian variant) | "ALDH2 deficiency" | rs671 (ALDH2) | A allele → strong flushing + nausea after alcohol; rare outside East Asian ancestry. | **Very high** (near-Mendelian). A-allele carriers reliably flush; G/G reliably do not. |
| Alcohol metabolism speed | "ADH1B Arg47His" | rs1229984 (ADH1B) | A allele → faster ethanol→acetaldehyde (buildup may worsen perceived hangover; weakly protective against alcoholism). | **Moderate**. Genetic effect is real and replicated, but behavioral variance (drinking habits, weight, timing with food) dominates perceived phenotype. |
| Earwax type + underarm body odor | "ABCC11 dry-earwax variant" | rs17822931 (ABCC11) | T/T → dry flaky earwax + reduced axillary odor (typical of East Asian ancestry). C-carrier → wet earwax + typical body odor. | **Very high** (Mendelian). Dominant C allele produces wet phenotype; T/T homozygotes produce dry phenotype. |
| Caffeine sensitivity | "CYP1A2 *1F/*1A" | rs762551 (CYP1A2) | A/A fast metabolizer; C/C slow; A/C intermediate. | **Low**. Single-variant explains a small fraction of caffeine-response variance; published effect sizes inconsistent. Treat as a weak signal; sleep/caffeine behavior dominates. |
| Supertaster for bitter compounds (PTC) | "TAS2R38 PAV/AVI haplotype" | rs713598 + rs1726866 + rs10246939 | Three-SNP haplotype. More "taster" alleles → stronger perception of bitterness (brussels sprouts, broccoli, grapefruit, tonic water). | **Moderate–high** at extremes (PAV/PAV ≈ strong taster; AVI/AVI ≈ non-taster). **Moderate** for intermediate haplotype combinations (3–5 taster alleles) where individual variance in self-reported sensitivity is large. |
| Blood type — **ABO component only** | — | rs8176719 (O-deletion: D=O allele, I=A-or-B allele), rs8176746, rs8176747 | Combine for A / B / AB / O typing. Full blood type requires Rh (next row). | **Very high** (Mendelian). ABO diplotype maps directly to phenotype. |
| Blood type — **Rh factor component** | "RhD antigen" | Direct RHD typing from tag SNPs is unreliable on most arrays (see "Rh factor gotcha" below). | Infer Rh+ vs Rh- from RHD gene-region coverage: if ~15+ RHD-region probes return clean calls → at least one functional RHD copy → Rh positive. High RHD-region no-call rate → likely homozygous deletion → Rh negative. | **High for Rh+ call** (gene-presence inference is robust). **Moderate for Rh- call** (requires interpretation of no-call / cross-hybridization patterns; serological confirmation recommended). |

### Blood-type reporting — always report both components

A user asking about blood type typically wants the full answer (e.g., A+, O-). **ABO alone is not sufficient.** When reporting blood type:

1. Give ABO from the three ABO SNPs above.
2. Give Rh (+/-) from RHD gene-region coverage (not from single tag SNPs — see below).
3. Combine into full blood type ("A+", "O-", etc.).
4. If Rh cannot be confidently inferred, say so and ask whether the subject knows their Rh status from medical records (blood donation cards, pregnancy records, surgery paperwork).

### Rh factor gotcha — array caveat

The Rh-negative phenotype in European populations is caused almost entirely by **complete deletion of the RHD gene** (not by SNP-level variation within an intact gene). Arrays cannot directly genotype a gene-scale deletion. **Single-SNP "Rh tag SNPs" in the literature (rs590787, rs676785, rs17418085, rs586178, rs675072) are poorly covered on v5 and the allele-direction mappings vary by strand and source — do not rely on a single tag SNP.**

Recommended approach: count probes in the RHD region (GRCh37 chr1:25,598,884–25,656,935) that returned clean calls vs. no-calls. A fully genotyped RHD region strongly indicates at least one functional RHD allele (Rh+). A high no-call rate across the RHD region is consistent with homozygous deletion (Rh-). This is not a clinical test — transfusion medicine requires serology.

## CFTR — variant encoding gotcha (23andMe i-probe IDs)

**23andMe uses internal Illumina probe IDs (prefixed `i`, e.g., `i3000001`) for variants that don't have standard rsIDs on the chip — primarily indels and custom-designed carrier-screen assays.** A lookup that only queries `rs`-prefixed IDs will miss these entirely.

The F508del deletion — the single most common cystic fibrosis variant — is the canonical example: dbSNP assigns it `rs113993960`, but that rsID is typically **not on the 23andMe v5 chip** as a standard SNP probe. F508del is assayed via the custom probe **`i3000001`** at chr7:117,199,646 (GRCh37). Allele encoding is **D = deletion present (pathogenic)** and **I = insertion / reference**. A D/I genotype is heterozygous F508del — i.e., a CF carrier.

### Rule: always query both rsID and i-probe ID for any CFTR carrier screen

Use this i-probe mapping for 23andMe chips (positions GRCh37, subject to confirmation against `reference/curated_snps.tsv` when available):

| Probe | chr7 pos | Variant | Protein |
|---|---|---|---|
| `i3000001` | 117,199,646 | F508del (c.1521_1523delCTT) | p.Phe508del |
| `i4000292` | 117,199,644 | F508del flanking | — |
| `i4000295` | 117,171,029 | R117H region (paired with rs78655421) | p.Arg117His |
| `i4000296` | 117,180,284 | 621+1G>T flanking | splice |
| `i4000300–i4000302` | 117,227,832–117,227,855 | G542X region | p.Gly542* |
| `i4000305–i4000306` | 117,227,860–117,227,865 | G551D region | p.Gly551Asp |
| `i4000307` | 117,227,887 | R553X region | p.Arg553* |
| `i4000317` | 117,227,792 | 1717-1G>A flanking | splice |
| `i4000318–i4000319` | 117,230,494 / 117,232,273 | 2184delA region | frameshift |
| `i4000320` | 117,242,922 | 3120+1G>A flanking | splice |
| `i4000321` | 117,246,808 | 3659delC flanking | frameshift |
| `i4000311` | 117,292,931 | W1282X region | p.Trp1282* |
| `i5012121` | 117,304,824 | N1303K region | p.Asn1303Lys |

The same principle extends to other carrier-screen panels encoded as i-probes (HEXA, HBB, SMN1 and others may use custom probes when the variant is an indel or when 23andMe added it specifically for carrier reports). Whenever a user reports a carrier-status result from the commercial report that the `rs`-based query fails to reproduce, **search the region by genomic coordinates and inspect i-probe calls** before concluding the variant isn't on the chip.

## Tier 4 — Polygenic scores (complex traits)

Polygenic scores (PRS) aggregate effects across hundreds to millions of SNPs. Different class of analysis from single-variant lookup: requires external weights files (usually from published GWAS summary statistics), chip-to-weights SNP overlap check, and careful ancestry calibration.

**Store these under the same ledger schema**, but with `variants: [...]` populated as an aggregate (number of weight SNPs, number overlapping the chip, missingness rate) rather than a single variant row. Effect size is the PRS z-score or percentile relative to a reference cohort.

### Reasonably well-supported (Tier B–C under normal rules)

| Trait | Source GWAS | N | Variance explained | Portability |
|---|---|---|---|---|
| Height | Yengo et al. 2022 (GIANT) | ~5.4M | ~40% in Europeans | Degrades significantly in other ancestries |
| BMI | Yengo et al. 2018 / Loh et al. | ~700k+ | ~6–9% | European-calibrated |
| Educational attainment (EA4) | Okbay et al. 2022 | ~3M | ~12–16% in Europeans | Poor outside European ancestry; within-family effects ~50% of between-family |
| Coronary artery disease | Inouye et al. / Khera et al. | ~500k+ | Top 5% PRS ≈ 3× CAD risk | European-calibrated; clinically validated for CAD risk stratification in some settings |
| Type 2 diabetes | Mahajan et al. 2022 | ~1.4M | Modest; useful combined with clinical risk factors | Multi-ancestry GWAS available |
| LDL cholesterol | GLGC | ~1.6M | ~10–15% | Good |
| Schizophrenia | PGC3 | ~150k | ~2.6% on liability scale | European |
| Major depressive disorder | PGC / Howard et al. | ~500k+ | ~2% | European |

### Weaker evidence (Tier D typically)

- Big Five personality dimensions (N ~100k per trait; variance explained single digits)
- Subjective well-being
- Risk tolerance / risk-taking
- Political orientation (weak effects)
- Specific cognitive subdomains beyond general factor

### How to run a PRS (implemented pipeline)

The project has a working end-to-end PRS pipeline. Files:

- `scripts/prs_download.py` — fetch PGS Catalog scoring files (harmonized to GRCh37/38). Usage: `python scripts/prs_download.py PGS000297 --build GRCh37`.
- `scripts/prs_pipeline.py` — reusable core. `load_weights()` parses the PGS Catalog format; `compute_score()` performs strand-aware allele matching and computes the weighted sum; `zscore_and_percentile()` converts raw score to percentile given a reference distribution.
- `scripts/prs_fetch_afs.py` — fetches 1000 Genomes Phase 3 EUR allele frequencies for the SNPs overlapping a subject's chip (uses Ensembl REST single-variant endpoint; batch endpoint ignores `pops` parameter). Caches to `reference/population_cache/1kg_eur_afs.json`. Usage: `python scripts/prs_fetch_afs.py PGS000297 <subject_id>`.
- `scripts/run_prs.py` — generic runner. Dispatches on `scripts/prs_traits.py`, which is the registry of per-trait calibration metadata (r², anthropometric anchors, claim/notes templates, self-report key for Rule-10/11 cross-checks). Add a new trait by adding an entry there; no new per-trait runner needed.

**Pipeline stages:**

1. **Download weights** — from PGS Catalog. Harmonized files include chr, pos, effect_allele, other_allele, effect_weight at minimum.
2. **Match by (chr, pos)** — build-aware. The current pipeline is GRCh37. Watch out for chr-column formatting (`"1"` vs `"chr1"`).
3. **Strand-aware allele matching:**
   - If subject alleles ⊆ {effect, other}: direct match. Dosage = count of effect in subject.
   - Else if subject alleles ⊆ {complement(effect), complement(other)}: strand flip. Dosage uses the complement.
   - For palindromic SNPs (A/T or C/G): only the direct-strand match is safe; if direct matches, accept (both PGS Catalog and 23andMe default to forward strand). If only the flip matches, skip — unresolvable from genotype alone.
   - If neither strand fits: genuine allele mismatch. Skip and count.
   - **Critical:** a subject homozygous for the OTHER allele is NOT a mismatch — it's dosage 0. Early pipeline versions got this wrong by requiring the effect allele to appear in the subject's genotype; fixed by switching to subset-checking against {effect, other}.
4. **Fetch per-SNP EUR allele frequencies** — Ensembl REST `/variation/human/{rsid}?pops=1`, filter for `1000GENOMES:phase_3:EUR`. Cache per-rsID. ~0.1 s/call × variants on overlap. Skip variants not in overlap.
5. **Compute reference distribution.** Two methods, in preference order:

   **5a. Empirical calibration from 1000G EUR (preferred).** Apply the PRS to the 503 EUR samples in the 1000 Genomes Phase 3 reference panel; take the empirical mean and SD across those scores. This is produced once per PGS via `scripts/calibrate_prs_empirical.py`, which writes `reference/population_cache/prs_empirical/<PGS_id>.json` with the mean, SD, and per-sample distribution. `run_prs.py` automatically uses this when present. Correctly captures LD-induced variance (which the theoretical method ignores). This is the published-standard approach for PRS z-score calibration.

   Caveat: 1000G EUR is a specific sample of 503 individuals across five sub-populations (CEU/FIN/GBR/IBS/TSI), weighted toward Northern/Western European. It's the best publicly-available proxy for "the European population distribution" but isn't identical to any specific training cohort (most PRS trained on UK Biobank don't have public individual-level calibration data). Sampling error on the empirical SD is ~3% at N=503. Record the reference as `1000G EUR N=503` in the finding's `calibration_method` field.

   **5b. Theoretical calibration from AFs (fallback, only when empirical not available).**
   - Mean = Σ 2 × p_i × β_i (over contributing SNPs with AF available)
   - Variance = Σ 2 × p_i × (1 − p_i) × β_i² (**assumes independence — ignores LD**; severely underestimates SD for LDpred2-style scores where correlated SNPs produce inflated true population variance)
   - Restrict the subject's score to the same subset of SNPs so scales match.
   - Fallback-to-fallback: if fewer than 50% of contributing SNPs have AF data, use uniform p=0.35.
   - Known failure mode: for dense LDpred2 scores the theoretical SD can be 3–4× too small, producing implausible z-scores. Always prefer empirical calibration when feasible.
6. **z-score, percentile, predicted trait value** — standardize subject score; convert to phenotype using published r² and population anchors (mean + SD for the target demographic, e.g., "European adult male height: 177 cm ± 7 cm"). Residual SD = trait_SD × √(1 − r²) sets the prediction uncertainty.
7. **Log to ledger** — aggregate `variants` row with panel size, chip overlap, contributing count, strand flips, palindromics skipped; `effect` carries raw_score, z_score, percentile, predicted value, CI, r², calibration method; `evidence_class: well_replicated_common_variant`; `inference_confidence` tier keyed off r² and ancestry match.

**Critical calibration pitfall:** when comparing the subject's partial-panel score against the theoretical reference distribution, compute both over the **same set of contributing SNPs**. Computing subject score over 333 SNPs and reference over all 3290 gives an apples-to-oranges comparison that yields absurd z-scores (the first version of the pipeline predicted 68 cm height; the fix was restricting reference to contributing SNPs and using real per-SNP AFs).

**Chip-coverage limits.** 23andMe v5 typically overlaps only ~15% of GWAS-scale PRS SNPs directly. Imputation (see "Genotype imputation pipeline" below) bumps coverage to 80–95%.

### Genotype imputation pipeline (local Beagle 5.4 + 1000G Phase 3 EUR)

Implemented end-to-end in `scripts/`. Converts a subject's chip parquet → imputed parquet with tens of millions of additional genotype calls inferred from linkage-disequilibrium patterns in the reference panel.

### Fresh-install / bootstrap instructions

`.gitignore` keeps large / re-fetchable / personal data out of the repo, so cloning gives you code + docs + curated tables but no runtime state. To rebuild:

1. Install Python deps: `pip install pandas pyarrow requests` (or your environment of choice).
2. `python scripts/install_portable_jdk.py` — ~200 MB, no admin rights.
3. `python scripts/imputation_download.py --beagle --maps` — Beagle + conform-gt JARs + genetic maps (~25 MB total).
4. `python scripts/imputation_download.py --chr all` — 1000G Phase 3 EUR reference panels, ~16 GB across 22 autosomes. Supports retry + resume on connection drops. Budget 30–60 min depending on your bandwidth.
5. `python scripts/filter_ref_vcf.py` — parallel-filter reference VCFs to biallelic SNPs with unique IDs (~5-10 min with default 6 workers). Produces `chr*.filtered.vcf.gz` alongside the originals.
6. Place your chip parquet at `standardized-genomes/<subject>.parquet` and a matching profile at `profiles/<subject>.json` (see existing format; fields include `subject_id`, `provider`, `raw_file`, `parse_stats`, `declared_ancestry`, `self_reported_phenotypes`).
7. `python scripts/parquet_to_vcf.py <subject> --chr all` — produce input VCFs per chromosome.
8. `python scripts/run_imputation.py <subject> --all --mem 8g --threads 4` — conform-gt + Beagle per chromosome; ~5–30 min depending on chip density.
9. `python scripts/vcf_to_parquet.py <subject> --chr all` — merge per-chromosome output into `standardized-genomes/imputed/<subject>.imputed.parquet`.

After bootstrap, the 16 GB `reference/imputation/1kg_ref_b37/` directory is safe to delete unless you plan to impute more subjects later. All PRS / downstream analyses read from `alice.imputed.parquet`, not from the reference panel.

**One-time setup:**

1. `python scripts/install_portable_jdk.py` — downloads Temurin OpenJDK 21 ZIP (~200 MB) and extracts to `reference/imputation/jdk/`. No admin rights needed.
2. `python scripts/imputation_download.py --beagle --maps` — Beagle JAR, conform-gt JAR, PLINK genetic maps.
3. `python scripts/imputation_download.py --chr all` — 1000G Phase 3 EUR reference panels (both bref3 format for Beagle and VCF format for conform-gt). Total ~10–12 GB across 22 autosomes. Downloader has retry + resume logic because the host occasionally drops connections mid-stream.

**Per-subject run:**

1. `python scripts/parquet_to_vcf.py <subject> --all` — reads `standardized-genomes/<subject>.parquet`, looks up canonical REF/ALT from the reference VCF per position, writes properly biallelic VCFs per chromosome. Handles strand flips (complement the subject's alleles if they don't match direct), dedupes multi-probe positions (23andMe sometimes has both `rs*` and `i*` probes at one coordinate), skips indels / no-calls / non-canonical bases.
2. `python scripts/run_imputation.py <subject> --all --mem 8g --threads 4` — per-chromosome: runs conform-gt to align alleles → runs Beagle with reference panel (bref3) + genetic map. Output: one imputed VCF per chromosome under `standardized-genomes/imputed/<subject>/beagle/chrN.vcf.gz`.
3. `python scripts/vcf_to_parquet.py <subject> --all` — merges per-chromosome VCFs into `standardized-genomes/imputed/<subject>.imputed.parquet`. Schema: `rsid, chrom, pos, ref, alt, a1, a2, dosage, gp_hom_ref, gp_het, gp_hom_alt, imputed` (bool; True for imputed variants, False for directly measured).

**What PRS runners should do.** Prefer the imputed parquet when present, fall back to the chip parquet. See `run_prs.py` (which does this) for the pattern. Imputed rows carry a `dosage` field in `[0, 2]`; PRS scoring should use the dosage (probability-weighted count of ALT) rather than the discrete 0/1/2 genotype for imputed SNPs, because dosages capture imputation uncertainty.

**Quality expectations.**
- Beagle study-marker acceptance: 98%+ (conform-gt handles allele alignment; chip positions not in the reference panel ~1-2% are legitimately skipped).
- Per-chromosome Beagle runtime: seconds (chr22) to a minute (chr1) on 4 threads.
- Imputation accuracy: >98% for MAF>5% in European ancestry; degrades for MAF<1% and for non-European subjects.
- Typical coverage expansion: ~50× per chromosome. A v5 chip (~640k SNPs) → imputed parquet (~30–40M SNPs).

**Known limitations.**
- Autosomes only. X-chromosome imputation uses a different reference panel format and script handling; not implemented.
- Palindromic SNPs (A/T or C/G) rely on the assumption that both input and reference use forward strand — correct for PGS Catalog harmonized files + 23andMe, but worth flagging if a different chip vendor is added.
- No error-correction for mis-genotyped sites on the chip. If the chip has a systematic miscall, it feeds through to imputation and can mislead inference at nearby sites. (Example pattern observed during validation: an rs12913832 HERC2 miscall was detected via LD-companion disagreement with rs1129038 and rs1667394 — multi-SNP haplotype consistency is a useful sanity check for single-probe miscalls at LD-tight loci.)

### Why naive PRS underperforms — the LD problem

SNPs near each other on a chromosome are inherited together. This is **linkage disequilibrium (LD)**: neighboring SNPs are correlated, so their genotypes aren't independent. When a GWAS reports 50 significant SNPs in a single LD block, those 50 SNPs are largely redundant — one or two are probably causal and the rest are tag SNPs riding along.

A naive weighted sum double-counts that signal. Every tag SNP contributes its effect size to the score, even though they're all measuring the same underlying genetic signal. In dense regions the score gets inflated; in sparse regions it's underweighted. Effect-size estimates from GWAS are also noisy (winner's curse: significant hits are systematically over-estimated).

Proper PRS methods address this:

- **Clumping + thresholding** (`PLINK --clump` then `PLINK --score`): within each LD block, keep the SNP with the lowest p-value and discard the rest. Then sum effects of the survivors. Simple, fast, works reasonably well for a baseline. Needs an LD reference panel to identify the blocks.
- **LDpred2** (R package `bigsnpr`): Bayesian. Treats effect sizes as drawn from a prior distribution and shrinks them based on the observed LD structure. Handles correlated SNPs by modeling the joint distribution rather than summing them independently.
- **PRS-CS** (Python/C++, Ge et al. 2019): similar Bayesian approach with a continuous shrinkage prior. Widely considered state-of-the-art for polygenic traits. Needs: GWAS summary statistics, LD reference panel (1000 Genomes subset by ancestry, ~several GB), Python environment with specific deps.

**LD reference panel** = a population dataset used to estimate which SNPs are correlated with which. 1000 Genomes Phase 3 (1KG) is standard; use the population subset matching your subject's ancestry (EUR, EAS, AFR, SAS, AMR). Pre-computed LD matrices for common panels are downloadable from the PRS-CS and LDpred2 authors, which saves the user from recomputing from raw VCFs.

**Practical recommendation:** for a first pass on a single trait, naive weighted sum or PLINK clump+score is fine — gets you 70–80% of the way with an hour of setup. For anything the user takes seriously, or when comparing across multiple traits, move to PRS-CS or LDpred2. Budget 2–4 hours for first-time setup (downloading LD panels, installing R or Python deps, getting weights in the right format), then subsequent runs are fast.

### Offspring PRS prediction

Given both parents' PRS, offspring PRS distribution is approximately Normal with:
- Mean = (PRS_parent1 + PRS_parent2) / 2 (midparent value)
- Variance ≈ (1/2) × σ²_PRS_population (due to Mendelian segregation)

Report as distribution with CI, not point estimate.

## Tier 5 — Kinship and relatedness (two subjects)

Requires both subjects loaded.

- **KING kinship coefficient** — robust to population structure. Values: 0.5 = self/MZ twin, 0.25 = first-degree (parent-child, full siblings), 0.125 = second-degree, 0.0625 = third-degree, etc. Reliable to ~5th degree.
- **PLINK IBD** — identity-by-descent estimation. Fast, works well for close relatives.
- **Parent-child check** — at every autosomal SNP where both parent and child are genotyped, child must carry at least one allele from each parent. Concordance should be essentially 100% (allowing for genotyping error ~0.1%). A handful of mismatches is normal; dozens indicate non-parentage or sample mix-up.
- **Shared ancestry inference** — PCA against 1000 Genomes reference panel. Places both subjects on a global ancestry map; shared components indicate common ancestry.

Tools: KING, PLINK, ADMIXTURE. Not trivial to set up on Windows; budget medium-to-high effort for first run.

**External database matching.** GEDmatch-style uploads, One-to-One / One-to-Many kinship searches, genealogy databases — all doable. Document which database, what was uploaded, and what came back in `investigations.jsonl`. User manages the external accounts; skill can prep files in the right format (typically 23andMe-style TSV from the parquet) and interpret returned matches.

## Critical caveats (every report must include)

1. **Array error rate.** Per-genotype error ~0.1% overall, higher for rare variants. Anything actionable → confirm via CLIA-certified clinical test.
2. **CYP2D6 unreliable.** Copy number variants and gene conversions invisible to arrays.
3. **BRCA non-comprehensive.** 23andMe's FDA-cleared report covers 3 Ashkenazi founder mutations out of thousands of pathogenic variants. Negative = nearly meaningless.
4. **No indel / CNV detection.** Arrays genotype SNPs. Deletions, insertions, repeat expansions, CNVs not represented.
5. **MTHFR clinically not actionable.** ACMG explicit recommendation against routine testing.
6. **Ancestry bias.** Most PRS and GWAS were trained on European cohorts. Non-European subjects: predictions degrade, sometimes severely.
7. **Phasing.** Array data is unphased. Some diplotype assignments involve assumption; flag explicitly.
8. **"Not on chip" ≠ "wild-type."** A variant absent from the chip is unobserved, not reassuring.

## Presentation conventions

These apply to both live conversation responses (when surfacing findings to the user) **and** to the `reports/` files. Ledger rows remain metrics-first; presentation is where the ledger data gets translated into something a curious non-specialist can read.

**Assume the reader doesn't know gene symbols, star-allele notation, drug compound names, or biochemical jargon.** Translate. Gene names, rsIDs, star alleles, and technical designations still appear — but next to plain-language context, not as the primary label.

### Rule 1 — Headers state the category AND the result, not just the gene

Header format:

> **[Plain-language category / dimension]: [plain-language result for this subject] [(inline metric if probabilistic)] *([technical identifiers: gene, star allele, rsID])*

The category names the dimension being measured (e.g., "Adult lactose tolerance," "Cystic fibrosis carrier status," "Blood type"). The result states the subject's specific finding on that dimension in plain English. The inline metric appears when phenotype is probabilistic — put the actual number (or a calibrated adjective like "mild," "strong") in the header so the reader sees the confidence without having to parse the body. Technical identifiers in parens at the end for search/traceability.

| Don't | Do |
|---|---|
| `APOE` | `Alzheimer's / cardiovascular risk: neutral baseline *(APOE ε3/ε3)*` |
| `UGT1A1 *1/*28` | `Gilbert's syndrome: carrier (mild benign bilirubin elevation expected) *(UGT1A1 *1/*28, rs887829 C/T)*` |
| `SLCO1B1 rs4149056` | `Statin muscle-pain risk: normal *(SLCO1B1 rs4149056 T/T — no *5 allele)*` |
| `CYP2C19` | `Blood-thinner / antidepressant / heartburn drug metabolism: normal metabolizer *(CYP2C19 *1/*1)*` |
| `Factor V Leiden rs6025` | `Blood-clot risk (Factor V Leiden): not a carrier *(rs6025 C/C)*` |
| `HFE C282Y` | `Iron overload / hemochromatosis risk: not a carrier *(HFE C282Y rs1800562 G/G)*` |
| `Adult lactose tolerance` *(category alone; reader can't tell what was found)* | `Adult lactose tolerance: intolerant genotype (~70% symptomatic expression in Europeans) *(rs4988235 G/G, MCM6)*` |
| `Eye color (dominant locus)` *(category alone)* | `Eye color: blue *(multi-SNP HERC2 inference; rs12913832 GG is a suspected miscall disagreeing with rs1129038 TT + rs1667394 TT)*` |
| `Cystic fibrosis carrier screen` *(category alone)* | `Cystic fibrosis carrier: F508del heterozygous carrier *(CFTR, 23andMe probe i3000001 D/I)*` |
| `Blood type` *(category alone)* | `Blood type: A+ *(ABO diplotype AO; Rh factor inferred positive from RHD gene-region coverage)*` |
| `Caffeine metabolism` *(category alone)* | `Caffeine metabolism: intermediate — weak signal *(CYP1A2 *1F/*1A, rs762551 A/C; variant explains small fraction of variance)*` |

**When the phenotype is probabilistic, put the number in the header.** Examples:

- "Adult lactose tolerance: **intolerant genotype (~70% symptomatic expression)**" — the 70% goes with the result so the reader sees the uncertainty up front, not buried in notes.
- "CYP1A2 caffeine metabolism: **intermediate — weak signal**" — where a number isn't published, a calibrated adjective ("weak signal," "moderate penetrance," "strong deterministic effect") serves the same role.
- "Rh factor: **very likely positive**" — explicit hedging in the header when the call is an inference, not a direct read.

The plain-language half tells the reader *what was found and how certain it is*. The parenthesized half gives the searchable identifiers so a curious reader can look deeper (dbSNP, CPIC, ClinVar, etc.).

Clarifications, caveats, self-report cross-checks, mechanism, and extended context go in the body below the header — not in the header itself.

### Rule 2 — Lead with one plain-English sentence

State in one sentence what the finding means in practical terms (what it affects, whether it's common, whether it requires action). Then the technical detail: genotype, diplotype, tier, effect size, mechanism, caveats.

**Template:**

> #### [Plain-language header] *([gene / star allele / rsID])*
>
> **[Plain-English one-liner — what it means for them in practice.]** Follow-up sentences with frequency context ("about 1 in 3 Europeans carry this"), what it does/doesn't imply, what to do or watch for.
>
> *Technical:* `rs... XY`, diplotype notation, tier, CI, source (CPIC/ClinVar/PubMed ID).

### Rule 3 — Proportion response to signal

- **Deviations from reference, Tier A/B findings, actionable items:** full plain-language treatment (header + one-liner + context + technical block).
- **Reference / unremarkable / negative results:** one-liner and move on. Example: `Blood-clot risk variant (Factor V Leiden) — not a carrier.` Don't force the reader through paragraphs to discover nothing notable.
- **"Not genotyped" at a locus the reader would reasonably expect to have been tested:** flag it explicitly. Absence of data ≠ absence of risk.

### Rule 4 — Technical rigor stays; it just moves

The ledger continues to carry study N, p-value, effect size, replication count, cohort ancestry, tier, and source IDs per the evidence framework. Those numbers appear in the `*Technical:*` block or in the methodology section of a report, not in the plain-language lede. Nothing is omitted — it's repositioned.

### Rule 5 — Never use an unexplained acronym or compound name

First use of any acronym, gene symbol, drug name, or biochemical term gets a gloss:

- "irinotecan (a colorectal-cancer chemo drug)"
- "atazanavir (an HIV medication)"
- "alpha-1 antitrypsin deficiency (a genetic cause of early lung and liver disease)"
- "Lp(a) (a cholesterol-like particle linked to heart disease)"

### Rule 6 — Reference-table notes must be specialized to the subject

When inline-quoting a reference-table note (e.g., from `reference/curated_snps.tsv` or this file's SNP tables), **annotate whether the caveat applies to the subject's actual chip, or strip it.** Reference-table notes describe generic coverage across chip versions; they're not claims about a loaded subject.

- Don't echo "may be absent on v4 chip" unchanged in a finding about a v5 subject.
- If the caveat is about a different chip version than the subject's, either drop it or prefix explicitly: "(re v4 chips only; this subject is on v5, which genotypes this SNP — see variant row)."
- In ledger `notes` fields, prefer concrete per-subject statements ("genotyped as GG, forward-strand reference") over generic reference-table boilerplate.

### Rule 7 — Translate trait-result labels, not just gene names

Rule 1 covers gene-name headers. Rule 7 covers **result labels**: the short phrases describing *which state* the subject is in for a trait. Many such labels are jargon lifted from the literature — intelligible only if you already know the biology. Always pair the technical label with the plain-English consequence.

| Don't (unexplained technical label) | Do (plain-English + technical) |
|---|---|
| "Lactase non-persistent" | "**Adult lactose intolerance** (technical: lactase non-persistent, G/G at rs4988235)" |
| "Lactase persistent" | "**Lactose-tolerant as an adult** (technical: lactase persistence allele present at rs4988235)" |
| "ALDH2 flush variant present" | "**Strong alcohol flush reaction** — nausea, redness, racing heart after even small amounts of alcohol (technical: rs671 A allele present)" |
| "Fast caffeine metabolizer (*1F/*1F)" | "**Caffeine clears quickly** — less per-dose effect, shorter jitters (technical: CYP1A2 *1F/*1F at rs762551)" |
| "Slow caffeine metabolizer (*1A/*1A)" | "**Caffeine lingers** — more jittery per dose, longer half-life (technical: CYP1A2 *1A/*1A at rs762551)" |
| "3/6 taster alleles" | "**Moderate sensitivity to bitter compounds** like brussels sprouts, grapefruit, tonic water (technical: TAS2R38 heterozygous PAV/AVI, 3 of 6 taster alleles)" |
| "Wet earwax, typical odor" | "**Wet, sticky earwax and typical body odor** (technical: C-carrier at rs17822931 ABCC11)" |
| "CYP2C19 *1/*1 Normal Metabolizer" | "**Typical metabolism** for clopidogrel, some antidepressants, and some proton-pump inhibitors (technical: CYP2C19 *1/*1)" |
| "Blood type A (AO)" | "**Blood type A** (one A-allele, one O-allele; ABO genotype AO). Rh factor reported separately below." |

**Don't report a trait result using only the technical label.** A reader who already understood the label wouldn't need the analysis; the whole point of the project is translation.

### Rule 8 — Compound traits: don't partially-answer

Some everyday concepts require multiple independent loci. A partial answer is worse than acknowledging the gap, because the reader may assume the partial answer is complete.

- **Blood type** = ABO type × Rh factor. Reporting "Blood type A" without Rh is incomplete; the reader expects "A+" or "A-". Always address both, or explicitly flag the missing half. See Tier 3 table and the Rh factor gotcha.
- **Eye color** at the full-resolution level needs HIrisPlex-S (11 SNPs), not just rs12913832. Single-SNP reporting is fine for the dominant blue/brown axis but should flag "green / hazel / intermediate shades require the multi-SNP model."
- **PGx drug-response profiles** usually involve multiple variants forming a diplotype; single-SNP reports can miss the call entirely (e.g., a CYP2C19 *17 allele changes the result dramatically). Always compute the diplotype across all relevant SNPs the chip covers.
- **Carrier screens** for recessive conditions are multi-variant by nature; reporting "no pathogenic variant found" is only meaningful if accompanied by which variants were actually screened — reference `reference/curated_snps.tsv` / the relevant i-probe tables explicitly.

### Rule 9 — Before sending, run the presentation checklist

For every finding surfaced to the user (in conversation or in a report), verify:

- [ ] Header is a plain-English phenotype, not a bare gene symbol (Rule 1).
- [ ] One plain-English sentence leads, before technical detail (Rule 2).
- [ ] Signal is proportioned: reference/unremarkable findings are one-liners; deviations get full treatment (Rule 3).
- [ ] First use of every acronym, gene, drug, or biochemical term has a gloss (Rule 5).
- [ ] Result labels are translated, not just gene names (Rule 7).
- [ ] If the concept is compound (blood type, eye color at resolution, PGx diplotypes), all components are addressed or missing-components are explicitly flagged (Rule 8).
- [ ] Confidence is attached and visible (Rule 10).
- [ ] Genomic evidence and self-reported phenotype are kept epistemically separate (Rule 11).

If a finding fails any of these, revise before sending. This is the cheapest quality gate in the project — a 30-second review that prevents the reader from having to decode terminology.

### Rule 10 — Every finding surfaces a confidence phrase, not just a tier

Tier (A/B/C/D/E) captures evidence quality for the **general claim** (how replicated is the association, how large the cohort, etc.). It does **not** capture how reliably the genotype predicts phenotype in a specific subject. Those are different axes, and both matter to the reader.

**Every trait or disease finding must display, in the normally visible body of the output (not just in ledger `effect` metadata):**

1. **Genotype-call confidence** — how reliable the subject's genotype call is. Usually "high" for array-genotyped common SNPs with clean calls; "moderate" for heterozygous calls on noisy probes; "low" or "suspected miscall" when cross-check evidence (LD companions, internal consistency) disagrees.
2. **Genotype → phenotype confidence** — how deterministic the mapping is in the general population. See the Tier 3 table's confidence column for per-trait values; for PGx / disease variants, use published penetrance where available.
3. **A one-line reasoning note** — why the confidence is what it is. Examples: "Mendelian, essentially deterministic" / "explains ~74% of variance, needs multi-SNP cross-check" / "penetrance ~70% for symptoms; dose and microbiome modulate."

**Format example:**

> #### Adult lactose tolerance *(rs4988235, MCM6)*
>
> **Genotype predicts reduced adult lactase activity, but phenotype is probabilistic.** ~70% of G/G Europeans have clinically noticeable symptoms with typical dairy doses; ~30% are asymptomatic.
>
> **Confidence:** genotype call — high (homozygous C/C, clean probe). Genotype → phenotype — moderate (penetrance ~70% for symptoms; dose, frequency, microbiome, product type all modulate).
>
> *Technical:* rs4988235 C/C; Tier A evidence for the persistence variant; cohort European (Enattah 2002, Itan 2010); N > 100,000 across replication cohorts.

If confidence is low or an array miscall is suspected, **lead with that** — don't bury it in a technical block. The reader's interpretation changes significantly when "likely miscall" enters the picture.

#### Rule 10.1 — PRS results across multiple traits use a cross-trait table

PRS findings are numerically dense (PRS z-score, percentile, phenotype anchor, predicted value + CI, residual SD, r²) in a way that single-SNP findings aren't. When multiple PRS are surfaced together — session reviews, `reports/`, cross-trait comparisons — use a **single table with rows = traits** so the reader can scan across traits without re-reading the same metric names.

This convention applies specifically to PRS findings. Single-SNP, PGx, carrier, and other finding types keep the prose-based Rule-1 / Rule-10 format; they don't need a table.

**Required columns, in this order:**

| Column | What it holds |
|---|---|
| **Trait (PGS id)** | plain-language trait name + the PGS Catalog ID |
| **PRS z-score** | signed, e.g., `+0.58 SD` |
| **Percentile** | on the PRS's own distribution |
| **Phenotype anchor** | population mean ± SD in the trait's units (e.g., `177 cm ± 7 cm`) |
| **Predicted phenotype** | the conditional-expectation value in trait units |
| **95% CI** | on the predicted phenotype |
| **r²** | variance captured, as a fraction (`0.40`) |
| **Naive (if r²=1)** | the value that PRS z-score would imply if PRS explained all variance. Shows the regression-to-mean shrinkage explicitly. |

**Example (synthetic three-PRS table for illustration — values are made up, not anyone's actual results):**

| Trait (PGS id) | PRS z | %ile | Phenotype anchor | Predicted | 95% CI | r² | Naive (if r²=1) |
|---|---|---|---|---|---|---|---|
| Height — Yengo 2022 (PGS002804) | +0.62 | 73 | 177 cm ± 7 cm (♂ EUR) | 180.7 cm | 170.1 – 191.4 cm | 0.40 | 181.3 cm |
| Educational attainment (PGS002231) | −0.25 | 40 | 14 yrs ± 3 yrs | 13.5 yrs | 8.0 – 19.0 yrs | 0.12 | 13.3 yrs |
| Cognitive ability, g-like (PGS002135) | +0.40 | 66 | 100 IQ ± 15 IQ | 101.6 IQ | 73.2 – 130.0 IQ | 0.07 | 106.0 IQ |

**Why the "Naive (if r²=1)" column matters.** Readers who haven't thought about this before tend to expect "+0.40 PRS SD = 106.0 IQ." The cross-trait table lets them see the un-shrunk naive value alongside the regression-to-mean-adjusted predicted value, making the effect of r² concrete. A low-r² trait shows naive ≫ predicted (in absolute deviation from the mean); a high-r² trait shows the two converging.

**When a single PRS is surfaced alone** (e.g., `run_prs.py` console output, or a report discussing only one trait), the cross-trait table collapses to a trivial one-row table. Prefer a **key-value list** in that case — it reads better for a single trait:

> Height — Yengo 2022 (PGS002804)
> - PRS z-score: +0.62 SD (73rd percentile)
> - Phenotype anchor: 177 cm ± 7 cm (European adult male)
> - Predicted: 180.7 cm (95% CI 170.1 – 191.4 cm)
> - r²: 0.40 (PRS explains ~40% of variance)
> - Naive (if r²=1): 181.3 cm — regression-to-mean attenuates to 180.7 given r²=0.40

**When r² is low (say < 0.15), lead with the PRS z-score / percentile in any prose framing.** The predicted-value row still appears in the table, but the sentence introducing the finding should frame it as PRS position ("The subject is at +0.40 SD on the cognitive-ability polygenic score"), not as a phenotype point estimate ("predicted 101.6 IQ"). Reserve leading with the phenotype number for high-r² PRS (height Yengo 2022 at r²=0.40).

**US-customary display by default.** Height in feet-and-inches; weight in pounds; temperature in °F. The ledger itself stores values in the canonical unit (typically SI: cm, kg, °C) so the data is interoperable, but the display layer converts to US-customary for any reader-facing output (conversation, `run_prs.py` console, `summarize_prs.py` tables, reports). This is configurable per-subject: set `display_preferences.units` to `"si"` in `profiles/<subject>.json` to override, or set the environment variable `GENOME_EVAL_UNITS=si` for a one-off. Helpers live in `scripts/run_prs.py` as `fmt_value(...)` and `fmt_anchor(...)`. Default applies to: cm → feet/inches, kg → lb, °C → °F. Traits already in US-natural units (years of schooling, IQ points) are unchanged.

**Only one active PRS per trait.** When a methodologically superior PRS becomes available for a trait already in the ledger (higher r², larger training cohort, better method — e.g., Yengo 2022 vs. an older 2020 score), supersede the older finding with a tombstone pointing to the new canonical. The older run stays in ledger history for provenance but is not an active finding. `scripts/summarize_prs.py` only reads active findings, so this automatically keeps cross-trait summaries showing the current best. Users reviewing the project should never see an inferior height PRS listed alongside the state-of-the-art one — the comparison is confusing and the older number is a worse answer for the same question.

### Rule 11 — Separate genomic evidence from self-reported phenotype; never conflate

This project is a genomic analysis tool. The primary claims in findings are **genomic inferences derived from the data**. Self-reported phenotype (eye color, lactose tolerance, drug response history, medical diagnoses) is valuable context but is a **different category of evidence** and must be tracked separately.

**Why separation matters.** If a subject self-reports "blue eyes" and a single SNP says "brown," it's easy to let the phenotype report silently override the genomic call ("must be a miscall"). Sometimes that's right, but the reasoning should rest on independent genomic evidence (linkage-companion SNPs disagreeing with the outlier call), not on the phenotype report. Otherwise the tool degrades into "user tells us what they have, we retrofit the genome to match" — which has zero analytical value. Genomic predictions must stand on genomic evidence; self-reports are a separate observation that can corroborate or contradict.

**Where self-reported phenotype lives.**
- `profiles/<subject>.json` carries a `self_reported_phenotypes` field. Each entry is `{value, reported_at, source, confidence}` where `source` is one of `subject_self_report`, `medical_record`, `clinical_test`, `family_observation`. Medical records and clinical tests get high confidence; self-reports get whatever the subject's certainty is.
- **Findings in `findings.jsonl` record genomic conclusions only.** The `claim`, `variants`, `effect`, and `tier_computed` fields reflect what the genotype data says, reasoned from the genotype data.
- When a finding's genomic conclusion and the subject's self-report either corroborate or diverge, record the cross-check as a **separate, clearly labeled line** in `notes` — e.g., `"Self-report cross-check: subject reports blue eyes (2026-04-17); matches genomic inference."` — never fold the self-report into the genomic reasoning.

**Workflow.**
1. Start from genotype data. Compute the best genomic inference using only genomic evidence (SNP calls, linkage, published weights, chip coverage).
2. Record the finding based on that inference, with its own confidence (Rule 10).
3. If the subject has reported phenotype for this trait, cross-check in a separate `notes` line. Flag **matches** and **mismatches** explicitly; a mismatch is valuable information (may indicate miscall, incomplete penetrance, misclassification, or genuine genotype-phenotype divergence — all worth investigating further).
4. When reporting to the user, surface both (genomic inference + self-report cross-check) as **separate statements** — not as one sentence that blends them.

**Do:**
> Genomic inference: blue eyes, high confidence. Two HERC2 linkage-companion SNPs (rs1129038 TT, rs1667394 TT) are homozygous for the blue-eye haplotype; the disagreeing rs12913832 GG call is the outlier and most likely a single-probe miscall.
>
> Self-report cross-check: subject reports blue eyes. Consistent with the genomic inference.

**Don't:**
> Because the subject has blue eyes, the rs12913832 GG call is a miscall.  ← phenotype driving the genomic conclusion

**Don't:**
> Eye color: blue (subject reports blue eyes; genomic data GG at HERC2).  ← blends the two, reader can't tell what the data alone says

## Reports

Generated on request to `reports/YYYY-MM-DD-<topic>.md`. Reports follow the presentation conventions above — plain-language headers and ledes, technical detail in sub-blocks. Structure:

```
# <Topic> — <Subject(s)>
Generated: <timestamp>
Ledger state: N findings at time of report

## Summary
<1 paragraph. State scope and main takeaways.>

## Findings
<For each finding in scope, sorted by tier then by clinical actionability:
  - Claim (plain-English; observed effect in practical terms)
  - Variants + genotypes
  - Effect on this subject (what it means for them)
  - Confidence: genotype call + genotype→phenotype + one-line reasoning (Rule 10)
  - Metrics: N, p, replication, ancestry
  - Tier
  - Sources
  - Notes / caveats
  - Self-report cross-check (if subject has reported phenotype for this trait):
    - State the subject's reported phenotype separately from the genomic inference (Rule 11)
    - Flag match / partial match / mismatch explicitly
    - A mismatch is a finding in itself — do not silently retrofit>

## Conflicts or gaps
<Findings that disagree; topics searched but not found; variants not on chip;
 genomic-inference vs. self-report mismatches with their candidate explanations
 (miscall, incomplete penetrance, misclassified self-report, genuine divergence).>

## Methodology
<Parse stats, tier rule version, which sources consulted.>

## Recommendations
<Where relevant: "confirm clinically," "discuss with physician/genetic counselor." No medical advice.>
```

Reports reference ledger finding IDs so claims are traceable.

### Shareable exports

If subject metadata flags `sharing_sensitivity: shareable-anonymized`, on request generate an anonymized export: replace subject display names with `subject_1` / `subject_2`, strip raw data file paths, keep findings and sources. Raw genotype data is **never** included in shared exports unless user overrides.

## Implementation starter — Python

Minimum deps: `pandas`, `pyarrow`. Optional: `requests` (Ensembl/gnomAD), `click` (CLI).

### Multi-provider parser dispatch

```python
from pathlib import Path
import gzip
import pandas as pd

def detect_provider(path: Path) -> str:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        head = "\n".join(fh.readline() for _ in range(30))
    if "23andMe" in head or "23andme" in head.lower():
        return "23andMe"
    if "AncestryDNA" in head or "ancestry.com" in head.lower():
        return "AncestryDNA"
    if "MyHeritage" in head:
        return "MyHeritage"
    if "FamilyTreeDNA" in head or "FTDNA" in head:
        return "FamilyTreeDNA"
    # Fallback: structural detection
    raise ValueError(f"Could not detect provider from header: {head[:500]}")

def parse_23andme(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                continue
            rsid, chrom, pos, geno = parts
            a1 = geno[0] if len(geno) >= 1 and geno[0] != "-" else None
            a2 = geno[1] if len(geno) >= 2 and geno[1] != "-" else a1 if chrom in ("Y", "MT") else None
            rows.append((rsid, chrom, int(pos), a1, a2))
    return pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "a1", "a2"])

def parse_ancestry(path: Path) -> pd.DataFrame:
    chrom_map = {"23": "X", "24": "Y", "25": "X", "26": "MT"}  # 25 = PAR
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("rsid"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 5:
                continue
            rsid, chrom, pos, a1, a2 = parts
            chrom = chrom_map.get(chrom, chrom)
            a1 = None if a1 == "0" else a1
            a2 = None if a2 == "0" else a2
            rows.append((rsid, chrom, int(pos), a1, a2))
    return pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "a1", "a2"])

def parse_myheritage(path: Path) -> pd.DataFrame:
    # CSV, quoted. Check header for build.
    df = pd.read_csv(path, comment="#")
    df.columns = [c.lower() for c in df.columns]
    # Expect rsid, chromosome, position, result
    out = pd.DataFrame({
        "rsid": df["rsid"],
        "chrom": df["chromosome"].astype(str),
        "pos": df["position"].astype(int),
    })
    out["a1"] = df["result"].str[0].where(df["result"].str[0] != "-", None)
    out["a2"] = df["result"].str[1].where(df["result"].str[1] != "-", None)
    return out

def parse_ftdna(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#")
    df.columns = [c.lower() for c in df.columns]
    out = pd.DataFrame({
        "rsid": df["rsid"],
        "chrom": df["chromosome"].astype(str),
        "pos": df["position"].astype(int),
    })
    out["a1"] = df["result"].str[0].where(df["result"].str[0] != "-", None)
    out["a2"] = df["result"].str[1].where(df["result"].str[1] != "-", None)
    return out

PARSERS = {
    "23andMe": parse_23andme,
    "AncestryDNA": parse_ancestry,
    "MyHeritage": parse_myheritage,
    "FamilyTreeDNA": parse_ftdna,
}
```

### Lookup helper

```python
def lookup(df: pd.DataFrame, rsid: str) -> tuple[str | None, str | None]:
    hit = df.loc[df["rsid"] == rsid]
    if hit.empty:
        return (None, None)  # not on chip
    row = hit.iloc[0]
    return (row["a1"], row["a2"])

def genotype_sorted(a1, a2) -> tuple[str, str] | None:
    if a1 is None and a2 is None:
        return None  # no-call
    if a2 is None:  # hemizygous
        return (a1, a1)
    return tuple(sorted([a1, a2]))
```

### Ledger writers

```python
import json, uuid
from datetime import datetime

def append_finding(ledger_path: Path, **kwargs):
    rec = {
        "finding_id": str(uuid.uuid4()),
        "timestamp": datetime.utcnow().isoformat() + "Z",
        **kwargs,
    }
    with open(ledger_path, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    return rec["finding_id"]
```

### VCF export

For downstream tooling, emit minimal VCFv4.2 from a standardized parquet. Use plus-strand alleles, GRCh37 reference. Needs `ref` column (not in raw data) — look up from curated table, from dbSNP, or skip variants where ref is unknown. Mark `##reference=GRCh37` explicitly. VCF is the lingua franca for external tools (bcftools, ANNOVAR, snpEff, KING, PLINK for some operations); it's worth producing once per subject and caching.

### Parquet schema (reference)

All standardized files have the same columns:

```
rsid: str
chrom: str              # '1'..'22', 'X', 'Y', 'MT'
pos: int                # GRCh37, 1-based
a1: str | None          # single char, plus-strand; None if no-call
a2: str | None          # single char; None for hemizygous (male X/Y, MT) or no-call
```

Provider and build are metadata properties, not per-row columns — stored in `profiles/<id>.json`.

## Framing sensitive or contested topics

Caveats to include in findings and reports, not reasons to skip analyses. Every analysis runs when the user asks; honesty comes through evidence metrics attached to the output.

- **Cognitive / intelligence PRS.** Largest current sources: Okbay et al. 2022 (EA4, N ≈ 3M, explains ~12–16% of variance in educational attainment in Europeans) and predecessor EA3 (Lee et al. 2018). IQ-proper GWAS are smaller and weaker. Under the standard tier rules these land in Tier B for European subjects, Tier C with ancestry downgrade. Record with every finding: (a) within-family effects are roughly half of population effects (direct vs indirect/assortative genetic components), (b) educational attainment is a proxy for cognitive ability with substantial environmental contribution, (c) portability outside the training ancestry is poor.

- **Personality / behavioral PRS.** Big Five GWAS up to N ~100k per dimension; variance explained typically single-digit percent. Political orientation, risk tolerance, subjective well-being have varying evidence. Most land in Tier D under standard rules. Run when asked, report with the actual metrics.

- **Kinship, relatedness, paternity.** Runs between any loaded subjects. KING and PLINK IBD give reliable coefficients to ~5th degree; parent-child verification is technically robust with a few hundred concordant SNPs.

- **Forensic-style identity analyses.** CODIS-style short tandem repeat loci aren't on SNP arrays, but SNP-based identification panels (Pakstis et al., Kidd et al.) can be computed. Ancestry-informative marker panels likewise. Run when asked.

- **External database matching.** Out-of-system queries (GEDmatch, genealogy platforms) are the user's to initiate; skill prepares export files and interprets returned results.

- **Medical findings.** Report with evidence grades; cite guidelines (CPIC, ACMG); note when clinical confirmation makes sense before acting. This is a factual point about array error rates and the gap between research and clinical testing — it travels with the finding, not as a prelude to it.

### Default answer

Run the analysis. Attach the metrics. Let the user judge.

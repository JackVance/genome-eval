# genome-eval

**Explore your own DNA data, honestly and in depth, on your own machine.**

If you've done a consumer DNA test (23andMe, AncestryDNA, MyHeritage, or
FamilyTreeDNA) and downloaded your raw data, this project reads that file and
produces a growing, evidence-graded notebook of findings about your genome:
drug-response predictions (pharmacogenomics), carrier status for recessive
diseases, ancestral origins, trait probabilities, and polygenic risk scores
for common conditions like heart disease and type 2 diabetes.

Everything runs locally. Your raw DNA data never leaves your machine. The
project is designed to **show its work** — every finding records the evidence
(study size, effect size, ancestry cohort, confidence), so you can see *why*
it's saying what it's saying, and re-evaluate later when the science
improves.

## What this does, step by step

1. **Reads your raw data file** (from 23andMe etc.) and converts it to a
   standard format the rest of the code can work with.
2. **Imputes** — fills in the genetic variants the chip didn't directly
   measure, using a public reference dataset of ~2,500 well-characterized
   genomes (see "What imputation means" below).
3. **Runs analyses**:
   - **Pharmacogenomics (PGx)**: which drugs you metabolize normally, fast,
     or slow, based on well-established genetic variants with FDA drug-label
     backing.
   - **Carrier screening**: do you carry a hidden copy of any serious
     recessive disease mutations (relevant if you plan to have children).
   - **Traits**: things like lactose tolerance, eye color, bitter taste
     sensitivity, athletic fast-twitch vs. endurance lean, caffeine metabolism.
   - **Polygenic scores (PRS)**: for complex conditions like coronary
     artery disease, type 2 diabetes, cholesterol, blood pressure, height —
     based on the combined effect of thousands to millions of small-effect
     variants.
   - **Haplogroups**: deep ancestry — your maternal line (mtDNA) and paternal
     line (Y-chromosome, males only) traced back tens of thousands of years.
4. **Logs findings to a ledger** — a chronological, append-only record with
   every piece of evidence attached. The ledger is designed so findings can
   be revised or superseded as the science evolves, without losing the history.
5. **Generates a report** on request — a human-readable markdown document
   summarizing the active findings in plain language.

## Key concepts, explained

**Genotype chip:** A consumer DNA test (like 23andMe v5) directly measures
roughly **640,000 specific positions** in your genome — a tiny sample of the
~3 billion total bases. Think of it as sampling the book of your DNA at
650,000 marked spots, not reading the whole thing.

**Imputation:** Because the chip only measures a sample, statistical methods
can "fill in" the unmeasured positions by comparing your pattern against a
reference database of **fully-sequenced genomes** (the 1000 Genomes Project,
2,500+ individuals). If your chip positions match a specific pattern seen in
the reference database, the unmeasured positions very likely match the
reference's nearby positions too. Imputation expands ~640K measured variants
to tens of millions of "known" variants with high confidence. **This is what
makes polygenic scores possible** — those scores need millions of variants,
far more than any chip measures directly.

**Polygenic score (PRS):** A single number that combines the effects of many
genetic variants (thousands to millions) to estimate your genetic
predisposition to a trait or condition. It's a *relative* number — where you
rank in the population — not a diagnosis. The predictive power varies a lot
by trait: height PRS explains ~40% of variance, coronary artery disease ~10%,
IQ/cognitive ability ~7%. Low variance explained = wide uncertainty around
any individual prediction.

**Carrier screening:** For recessive diseases, you can be a "carrier" —
you have one copy of a mutation but aren't sick yourself. If both you and
a partner are carriers for the same disease, your children have a 1-in-4
chance of being affected. The major ACMG-recommended genes to screen are
CFTR (cystic fibrosis), SMN1 (spinal muscular atrophy), HEXA (Tay-Sachs),
HBB (sickle cell / β-thalassemia), and GJB2 (hereditary hearing loss).

**Pharmacogenomics (PGx):** How your genetic variants affect your response
to specific drugs. For example, a variant in UGT1A1 predicts how you handle
certain chemo drugs; VKORC1 predicts warfarin dosing; CYP2C19 predicts
clopidogrel response. These are the genetic findings with the strongest
clinical backing — many have FDA-recognized drug labels.

**Haplogroup:** Your deep-ancestry maternal (mtDNA) or paternal (Y) lineage
— a label that traces back through specific branches of a human phylogenetic
tree covering 50,000+ years of migration. E.g., "R1b-M269" is the most common
Western European paternal line; "U5a" is one of the oldest surviving
European maternal lines, associated with pre-Neolithic hunter-gatherers.

**Evidence tier:** Every finding in the ledger has a letter grade (A–E) for
evidence quality: A for ClinVar/ACMG-level pathogenic variants, B for
well-replicated common variants, C for moderately-supported associations,
D for weaker signals, E for single-study claims. This lets you distinguish
"this is genuine medical-grade information" from "this is a curious
statistical association I wouldn't act on."

**Confidence:** Separate from evidence tier, every trait/disease finding
also reports two axes of confidence: (1) how reliable the genetic call is
on this chip for this variant, and (2) how reliably the genotype predicts
the phenotype. A variant can be A-tier evidence but still have
moderate-to-low confidence for a specific person (e.g., the genotype call
is crisp, but the phenotype is only 70% predicted by the genotype).

## What you get — example output structure

The final output is a markdown report grouped by category, for example:

```
### Cystic fibrosis carrier: F508del heterozygous carrier (CFTR, probe i3000001 D/I)

You carry one copy of F508del — the single most common CF-causing mutation.
You don't have CF (that requires two damaged copies), but any child with
a partner who is also a CF carrier has a 1-in-4 chance of being affected.

- Evidence tier: A (pathogenic, extensively replicated).
- Genotype-call confidence: high.
- Genotype→phenotype confidence: high — F508del causes CF with
  near-complete penetrance when in trans with another pathogenic CFTR variant.
```

...plus pharmacogenomics tables, carrier status across multiple genes,
polygenic scores with 95% confidence intervals, haplogroup background, and
"not callable from this chip" sections for honest limitations (some things,
like CYP2D6 copy number, cannot be determined from chip data alone).

## Privacy and safety

**Your data never leaves your machine.** All analyses run locally. Nothing
is uploaded to any service. The project's `.gitignore` is designed for
public-repo safety — if you push the code to GitHub, everything
subject-identifying stays behind on your machine (your raw file, the
standardized parquet, your profile, the ledger, reports, anything under
`local/`). Only the generic code and curated reference tables are tracked.

The reference data this project uses (1000 Genomes Project, PGS Catalog,
ClinVar, PhyloTree, ISOGG) is **all public** research data — none of it is
personal to you or anyone else you'd know.

## Prerequisites

- **Python 3.13** (developed and tested on 3.13.12; 3.11+ probably works).
- **Git** (one dependency installs via `pip install git+...`).
- **~40 GB free disk** for the full reference-panel + 1000G VCF set.
  Reclaimable after your imputation runs (see
  `NEXT_STEPS.md` → "Space management").
- **Windows, macOS, or Linux.** Tested on Windows 11.
  - On Windows, use Git Bash or WSL.
- **No system Java required** — the project bootstraps its own portable
  Java runtime.
- **Command-line comfort.** You don't need to write code, but you should
  be comfortable running `pip install`, copying files, and running
  `python scripts/foo.py` from a terminal.

## One-time setup

```bash
# 1. Clone the repo
git clone <repo-url> genome-eval && cd genome-eval

# 2. Create a Python virtual environment and install dependencies.
#    A venv keeps this project's Python packages isolated from the rest
#    of your system.
python -m venv .venv
source .venv/Scripts/activate           # Windows Git Bash
# or: source .venv/bin/activate         # macOS / Linux
pip install -r requirements.txt

# 3. Download the external tools and reference data.
#    Each script is idempotent — re-running is safe, nothing is duplicated.
#    Expect ~45 minutes total wall time here, mostly downloading.

python scripts/install_portable_jdk.py       # Java runtime (~200 MB)
python scripts/imputation_download.py        # Beagle imputer + 1000G reference panels (~16 GB)
python scripts/download_1kg_canonical.py     # 1000G full VCFs, for PRS calibration (~25 GB)
python scripts/haplogrep_download.py         # HaploGrep3 (mtDNA haplogroup caller)
python scripts/extract_eur_afs.py            # Build the allele-frequency table (~20 min)

# 4. Drop your raw DNA file into raw-source-genomes/.
#    This is the file you downloaded from 23andMe / AncestryDNA / MyHeritage /
#    FamilyTreeDNA. Leave it untouched — the project treats this directory as
#    a read-only boundary.

cp /path/to/genome_YourName_v5_Full_*.txt raw-source-genomes/

# 5. Ingest: detect provider, parse, and produce the standardized internal
#    format. You pick a short ID for yourself (any lowercase word — "alice",
#    "me", a nickname). This ID is used in filenames downstream.

python scripts/normalize.py <subject_id>     # short filename-safe id, e.g. "alice"
```

## Running analyses on your data

After setup, every analysis is a single command. The long one is imputation
— expect 30–90 minutes of mostly-background compute for that step. After
that, everything else takes seconds to a few minutes.

```bash
# 6. Imputation (the long one — 30–90 min per subject; runs in background-like fashion).
#    This is the step that fills in the ~97% of variants the chip didn't
#    directly measure, by cross-referencing against the 1000 Genomes panel.
python scripts/parquet_to_vcf.py <subject_id> --all
python scripts/run_imputation.py <subject_id> --chr all
python scripts/vcf_to_parquet.py <subject_id>

# 7. Run the analyses — each is independent and can be run in any order.
python scripts/run_mtdna_haplogroup.py <subject_id>         # maternal-line deep ancestry
python scripts/run_y_haplogroup.py <subject_id>             # paternal-line deep ancestry (males only)
python scripts/run_carrier_panel.py <subject_id>            # recessive-disease carrier status
python scripts/run_prs.py <subject_id> <trait>              # polygenic score for one trait
python scripts/investigate_pgx_cardio.py <subject_id>       # cardiovascular drug-metabolism panel
# ... there are more investigate_* scripts for specific trait families;
# see the full catalogue in scripts/ or the "Common operations" table below.

# 8. Generate a human-readable report of everything in your ledger.
python scripts/generate_report.py <subject_id>
# Writes reports/YYYY-MM-DD-<subject>.md — open it in any markdown viewer.
```

## What each command does (quick reference)

| Command | What it does |
|---|---|
| `normalize.py <id>` | Detect your DNA provider, parse the raw file, produce the standard internal parquet |
| `parquet_to_vcf.py <id> --all` | Convert chip data to VCF format for the imputer |
| `run_imputation.py <id> --chr all` | Fill in the unmeasured variants using Beagle + 1000 Genomes |
| `vcf_to_parquet.py <id>` | Merge imputed per-chromosome files into one parquet |
| `run_mtdna_haplogroup.py <id>` | Determine maternal deep-ancestry haplogroup |
| `run_y_haplogroup.py <id>` | Determine paternal deep-ancestry haplogroup (males) |
| `run_carrier_panel.py <id>` | Screen for recessive-disease carrier status across 5 major genes |
| `run_prs.py <id> <trait>` | Compute a polygenic score for one trait (see `scripts/prs_traits.py` for the list) |
| `investigate_pgx_cardio.py <id>` | Cardiovascular pharmacogenomics (warfarin, clopidogrel, statins, etc.) |
| `summarize_findings.py <id>` | Print a concise list of all active findings in the ledger |
| `summarize_prs.py <id>` | Print a cross-trait table of all polygenic scores |
| `generate_report.py <id>` | Write a full markdown report to `reports/` |
| `reclassify.py <id>` | Re-derive evidence tiers after a rules update |
| `calibrate_prs_empirical.py <PGS_id> --subject-observed <id>` | Recalibrate a PRS if it gives an implausible result — see `reference/population_cache/prs_empirical/README.md` |

## Dependency summary

### Python packages (installed by `pip install -r requirements.txt`)

| Package | Role |
|---|---|
| `pandas` | Data wrangling; the standardized genotype format is a pandas DataFrame backed by parquet. |
| `pyarrow` | Parquet I/O engine for pandas. |
| `requests` | HTTP fetch in a handful of bootstrap scripts (`prs_download.py`, others). |
| `yhaplo` | Y-chromosome haplogroup caller (ISOGG tree). Installed from GitHub. Non-commercial-use license. |

### External tools (installed by the bootstrap scripts — no manual install)

| Tool | Purpose | Bootstrap script |
|---|---|---|
| Temurin JDK 21 | Runs HaploGrep3 and Beagle | `scripts/install_portable_jdk.py` |
| Beagle 5.4 + conform-gt | Genotype imputation | `scripts/imputation_download.py` |
| HaploGrep3 3.2.2 | mtDNA haplogroup classifier (PhyloTree 17.2) | `scripts/haplogrep_download.py` |
| 1000G Phase 3 EUR reference panels | Imputation target (bref3 + VCF) | `scripts/imputation_download.py` |
| 1000G Phase 3 EBI release VCFs | Empirical PRS calibration (503 EUR samples) | `scripts/download_1kg_canonical.py` |
| rCRS (NC_012920.1) mtDNA reference | mtDNA reference sequence | Downloaded on first mtDNA run |
| PhyloTree 17.2 | mtDNA haplogroup tree | Bundled with HaploGrep3 |
| ISOGG Y-tree 2016.01.04 | Y-haplogroup tree | Bundled with the yhaplo package |

## Directory layout

```
genome-eval/
├── README.md                       # This file
├── CLAUDE.md                       # Project stance, invariants, working notes (for Claude Code sessions)
├── MAINTENANCE.md                  # Governance for convention evolution
├── NEXT_STEPS.md                   # Living roadmap and open items
├── requirements.txt                # Python dependencies
├── .gitignore                      # Excludes personal data + re-fetchable infrastructure
│
├── raw-source-genomes/             # Your provider DNA file — never edited, gitignored
├── standardized-genomes/           # Normalized parquet (analyses read from here, gitignored)
│   ├── <id>.parquet                # Chip-only data
│   ├── imputed/<id>.imputed.parquet  # After imputation — this is what PRS uses
│   └── haplogroups/<id>.*          # mtDNA / Y outputs
├── profiles/                       # Per-subject metadata (name, sex, self-reports) — gitignored
│
├── ledger/                         # Append-only findings log — gitignored
│   ├── findings.jsonl
│   ├── sources.jsonl
│   └── investigations.jsonl
│
├── reports/                        # Generated markdown reports — gitignored
├── local/                          # Per-subject working notes — gitignored
│
├── reference/
│   ├── curated_snps.tsv            # Curated variant table (tracked)
│   ├── carrier_panels/*.tsv        # ACMG per-gene variant panels (tracked)
│   ├── imputation/                 # Beagle + 1000G panels (gitignored, re-fetchable)
│   ├── haplogroups/
│   │   ├── mtdna/rCRS.fasta        # Reference sequence (tracked)
│   │   └── mtdna/haplogrep3/       # HaploGrep3 (gitignored, re-fetchable)
│   ├── population_cache/           # 1000G EUR frequencies + PRS calibrations (mostly gitignored)
│   └── prs_weights/                # PGS Catalog files (gitignored, re-fetchable)
│
└── scripts/                        # All runners and library code (tracked)
```

## What's tracked vs. gitignored

**Tracked** (safe to push to public remotes):

- Code, `README.md`, `CLAUDE.md`, `MAINTENANCE.md`, `NEXT_STEPS.md`, `requirements.txt`
- `reference/curated_snps.tsv` and `reference/carrier_panels/*.tsv` — curated variant tables
- `reference/haplogroups/mtdna/rCRS.fasta` — public reference sequence
- `reference/population_cache/prs_empirical/<PGS>.json` — full-panel empirical distributions (content is 1000G data, not yours)

**Gitignored** (never leaves your machine):

- All personal data: `raw-source-genomes/`, standardized parquets, profiles, ledger, reports, imputed data, haplogroup outputs, `local/` notes
- Subject-specific PRS calibrations (`<PGS>.<subject>.json`)
- All large re-fetchable infrastructure: 1000G reference panels, EBI VCFs, JDK, genetic maps, JARs, PGS weight files

A fresh clone gives a collaborator code + conventions + curated reference tables. Each user bootstraps their own local state via the one-time setup.

## Troubleshooting

- **`java: command not found` / HaploGrep3 complains about JRE**: the bundled JDK
  works but Windows `.exe` wrappers may not find it. The scripts invoke the
  JAR directly via the project's own JDK, so this is usually only a problem if
  you're invoking HaploGrep3 by hand.
- **`pip install yhaplo` fails with "no matching distribution"**: a name-squatter
  occupies `yhaplo` on PyPI. The real package is the GitHub repo — `requirements.txt`
  already points at it correctly.
- **yhaplo reports "root haplogroup assigned"**: the `.genos.txt` input format is
  tab-delimited starting with `ID`, despite what `--help` implies. `run_y_haplogroup.py`
  already handles this.
- **PRS produces an extreme z-score (|z| > 5)**: usually either a coordinate-column
  mismatch (already fixed in `prs_pipeline.py` for harmonized PGS Catalog files) or
  a missing-variant bias from sub-100% coverage (fix: use
  `calibrate_prs_empirical.py --subject-observed <id>`). See
  `reference/population_cache/prs_empirical/README.md` for details.

## Further reading

- **`CLAUDE.md`** — project stance (no refusals, metrics-first, epistemic caveats),
  hard invariants, presentation conventions, and a dense log of working-notes
  about specific bugs and gotchas. If you're using this with Claude Code, this
  is your always-on context.
- **`NEXT_STEPS.md`** — living roadmap. What's open, what's resolved, what's
  intentionally out of scope.
- **`MAINTENANCE.md`** — how conventions evolve. Read if you want to contribute
  or fork and extend.
- **`.claude/skills/genome-eval/SKILL.md`** — detailed workflow, SNP tables by
  tier, PRS method docs, presentation rules with examples. This is the
  deep-dive reference for anyone actively running the analyses.
- **`local/README.md`** — describes the gitignored per-subject notes convention.

## License

Project code is MIT-licensed — see `LICENSE` for the full text.

**Important:** the MIT license covers *this project's* code and curated data
tables. Third-party runtime tools (yhaplo, Beagle, HaploGrep3, Temurin JDK)
are downloaded at runtime and carry their own licenses:

- **yhaplo** — non-commercial academic research use only.
- **Beagle 5.4** — GPL-3.
- **HaploGrep3** — MIT.
- **Temurin JDK 21** — GPLv2 with Classpath Exception.

If your use case is commercial or you're redistributing derivative tools,
check the upstream licenses directly.

## Contributing / extending

See `MAINTENANCE.md`. Short version: small additions go in `scripts/`, new
conventions go in `CLAUDE.md` working notes, detailed workflow guidance goes
in `.claude/skills/genome-eval/SKILL.md`. Before pushing to a public remote,
run the pre-push audit checklist in `MAINTENANCE.md` → "Public-repo safety"
to verify no subject-identifying content has drifted into tracked files.

# Carrier-screening panels

Per-gene TSVs of pathogenic / likely-pathogenic variants for ACMG tier-1
carrier-screening genes. Queried by `scripts/run_carrier_panel.py`.

## Schema

Every panel TSV uses:

| Column | Description |
|---|---|
| `gene` | HGNC gene symbol |
| `rsid` | dbSNP identifier (or blank if indel-only / probe-only) |
| `hgvs_c` | HGVS cDNA nomenclature |
| `hgvs_p` | HGVS protein nomenclature (or blank for splice/UTR) |
| `chrom` | GRCh37 chromosome (1–22, X, Y, MT) |
| `pos_b37` | GRCh37 position (integer) |
| `ref` | REF allele on the forward strand |
| `alt` | ALT allele on the forward strand |
| `probe_type` | `rsid` (SNP-tractable on array), `iprobe` (needs i-probe), `cnv` (not callable from array) |
| `condition` | Disease / phenotype name |
| `notes` | One-line annotation — founder population, clinical flag, verify note |

## What's here

- `cftr.tsv` — cystic fibrosis, ACMG-23 + CFTR2 common variants
- `hexa.tsv` — Tay-Sachs (HEXA gene)
- `hbb.tsv` — sickle cell + β-thalassemia (HBB gene)
- `gjb2.tsv` — nonsyndromic hearing loss (connexin 26 / GJB2)
- `pah.tsv` — phenylketonuria (PAH gene)

## Not included (not reliably callable from SNP array / imputation)

- **SMN1** (SMA) — requires MLPA or qPCR dosage; paralog SMN2 confounds.
- **FMR1** (Fragile X) — CGG triplet expansion, not SNP-based.
- **CYP21A2** (CAH) — pseudogene CYP21A1P has >98% homology.
- **PKD1** — 6 pseudogenes on chr16; private variants; no common panel.

## Sources

- CFTR2.org quarterly releases (CFTR panel)
- ClinGen Variant Curation Expert Panels (CFTR, Hearing Loss, Hemoglobinopathy)
- ACMG Gregg 2021 *Genet Med* (tier-3 pan-ethnic panel)
- ClinVar (`clinvar_GRCh37.vcf.gz`) for variant-level authoritative coordinates
- gnomAD v2.1.1 for population AF sanity-check

## Adding a new panel

1. Create `<gene>.tsv` with the schema above.
2. Add the gene to the default panel list in `run_carrier_panel.py`.
3. Positions must be GRCh37 — any hg38-only variant needs lift-over first.

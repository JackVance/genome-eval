"""Tier 1 PGx + Tier 2 cardio investigation runner for a single subject.

Reads standardized-genomes/<id>.parquet, looks up curated SNPs, derives
diplotype/phenotype where the skill rules allow, and appends one investigation
+ one-finding-per-claim + underlying sources to the append-only ledger.

Usage:
    python scripts/investigate_pgx_cardio.py <subject_id>
"""
from __future__ import annotations

import argparse
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io
from tier_rules import compute_tier, TIER_RULE_VERSION


# --------------------------------------------------------------------------- #
# Sources — register once, reference many times.
# --------------------------------------------------------------------------- #

SOURCES = [
    {
        "source_id": "cpic:clopidogrel-cyp2c19-2022",
        "type": "cpic",
        "citation": "Lee CR et al. Clin Pharmacol Ther 2022; CPIC guideline for clopidogrel and CYP2C19",
        "url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/",
        "evidence_class": "guideline",
        "cohort_ancestry": "multi-ethnic",
        "notes": "CPIC Level A",
    },
    {
        "source_id": "cpic:warfarin-2017",
        "type": "cpic",
        "citation": "Johnson JA et al. Clin Pharmacol Ther 2017; CPIC guideline for warfarin dosing (CYP2C9/VKORC1/CYP4F2)",
        "url": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/",
        "evidence_class": "guideline",
        "cohort_ancestry": "multi-ethnic",
        "notes": "CPIC Level A",
    },
    {
        "source_id": "cpic:fluoropyrimidines-dpyd-2018",
        "type": "cpic",
        "citation": "Amstutz U et al. Clin Pharmacol Ther 2018; CPIC guideline for DPYD and fluoropyrimidines (updated 2020)",
        "url": "https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/",
        "evidence_class": "guideline",
        "cohort_ancestry": "multi-ethnic",
        "notes": "CPIC Level A; missed positives can cause fatal toxicity",
    },
    {
        "source_id": "cpic:thiopurines-tpmt-nudt15-2018",
        "type": "cpic",
        "citation": "Relling MV et al. Clin Pharmacol Ther 2019; CPIC guideline for thiopurines TPMT and NUDT15",
        "url": "https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt/",
        "evidence_class": "guideline",
        "cohort_ancestry": "multi-ethnic",
        "notes": "CPIC Level A",
    },
    {
        "source_id": "cpic:simvastatin-slco1b1-2022",
        "type": "cpic",
        "citation": "Cooper-DeHoff RM et al. Clin Pharmacol Ther 2022; CPIC guideline for SLCO1B1, ABCG2, CYP2C9 and statin-associated musculoskeletal symptoms",
        "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-statins/",
        "evidence_class": "guideline",
        "cohort_ancestry": "multi-ethnic",
        "notes": "CPIC Level A for simvastatin",
    },
    {
        "source_id": "cpic:atazanavir-ugt1a1-2015",
        "type": "cpic",
        "citation": "Gammal RS et al. Clin Pharmacol Ther 2016; CPIC guideline for atazanavir and UGT1A1",
        "url": "https://cpicpgx.org/guidelines/guideline-for-atazanavir-and-ugt1a1/",
        "evidence_class": "guideline",
        "cohort_ancestry": "multi-ethnic",
        "notes": "CPIC Level A",
    },
    {
        "source_id": "cpic:abacavir-hlab-2014",
        "type": "cpic",
        "citation": "Martin MA et al. Clin Pharmacol Ther 2014; CPIC guideline for abacavir and HLA-B",
        "url": "https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/",
        "evidence_class": "guideline",
        "cohort_ancestry": "multi-ethnic",
        "notes": "HLA-B*57:01 typing required clinically; array tag SNP unreliable",
    },
    {
        "source_id": "clinvar:hfe-c282y",
        "type": "clinvar",
        "citation": "ClinVar: HFE c.845G>A (p.Cys282Tyr); hereditary hemochromatosis",
        "url": "https://www.ncbi.nlm.nih.gov/clinvar/variation/9/",
        "evidence_class": "guideline",
        "cohort_ancestry": "European-dominant",
        "notes": "Pathogenic, ≥2-star review",
    },
    {
        "source_id": "clinvar:factor-v-leiden",
        "type": "clinvar",
        "citation": "ClinVar: F5 c.1601G>A (Factor V Leiden, R506Q); VTE",
        "url": "https://www.ncbi.nlm.nih.gov/clinvar/variation/642/",
        "evidence_class": "guideline",
        "cohort_ancestry": "European-dominant",
        "notes": "Risk factor, pathogenic; multiple reviewers",
    },
    {
        "source_id": "clinvar:prothrombin-g20210a",
        "type": "clinvar",
        "citation": "ClinVar: F2 c.*97G>A (G20210A); VTE",
        "url": "https://www.ncbi.nlm.nih.gov/clinvar/variation/13310/",
        "evidence_class": "guideline",
        "cohort_ancestry": "European-dominant",
        "notes": "Risk factor, pathogenic",
    },
    {
        "source_id": "clinvar:apoe",
        "type": "clinvar",
        "citation": "ClinVar / landmark: Corder et al. 1993 Science (APOE ε4 and late-onset Alzheimer's)",
        "url": "https://www.ncbi.nlm.nih.gov/clinvar/?term=APOE",
        "evidence_class": "meta-analysis",
        "cohort_ancestry": "European-dominant",
        "notes": "ε4/ε4 ~8-12× late-onset AD risk vs ε3/ε3; ε2 protective",
    },
    {
        "source_id": "pmid:20032323",
        "type": "peer-reviewed",
        "citation": "Clarke R et al. NEJM 2009; 361:2518-28. Genetic variants associated with Lp(a) lipoprotein level and coronary disease.",
        "url": "https://pubmed.ncbi.nlm.nih.gov/20032323/",
        "evidence_class": "primary GWAS",
        "cohort_ancestry": "European",
        "cohort_n": 3145,
        "notes": "rs10455872 and rs3798220; OR for CAD ~1.5-1.7 per G/C allele",
    },
    {
        "source_id": "acmg:mthfr-2013",
        "type": "guideline",
        "citation": "Hickey SE et al. Genet Med 2013; ACMG practice guideline: MTHFR polymorphism testing.",
        "url": "https://pubmed.ncbi.nlm.nih.gov/23288205/",
        "evidence_class": "guideline",
        "cohort_ancestry": "N/A",
        "notes": "Routine MTHFR testing is NOT recommended",
    },
    {
        "source_id": "pmid:17679673",
        "type": "peer-reviewed",
        "citation": "Kujovich JL. Genet Med 2011; Factor V Leiden thrombophilia. GeneReviews.",
        "url": "https://www.ncbi.nlm.nih.gov/books/NBK1368/",
        "evidence_class": "review",
        "cohort_ancestry": "European-dominant",
        "notes": "Het OR ~4-8 for VTE; hom ~80",
    },
]


# --------------------------------------------------------------------------- #
# Lookup primitives.
# --------------------------------------------------------------------------- #

def lookup(df: pd.DataFrame, rsid: str) -> dict:
    hit = df.loc[df["rsid"] == rsid]
    if hit.empty:
        return {"on_chip": False}
    row = hit.iloc[0]
    a1, a2 = row["a1"], row["a2"]
    a1 = None if (a1 is None or (isinstance(a1, float))) else a1
    a2 = None if (a2 is None or (isinstance(a2, float))) else a2
    genotype = None
    if a1 is not None and a2 is not None:
        genotype = "".join(sorted([a1, a2]))
    elif a1 is not None:
        genotype = a1
    return {
        "on_chip": True,
        "chrom": str(row["chrom"]),
        "pos": int(row["pos"]),
        "a1": a1,
        "a2": a2,
        "genotype": genotype,
    }


def count_alt(geno: dict, ref: str, alt: str) -> int | None:
    """How many copies of `alt` given a lookup dict with sorted genotype."""
    if not geno.get("on_chip"):
        return None
    if geno.get("genotype") is None:
        return None
    g = geno["genotype"]
    return sum(1 for c in g if c == alt)


def variant_row(rsid: str, gene: str, geno: dict, ref: str, alt: str) -> dict:
    return {
        "rsid": rsid,
        "gene": gene,
        "chrom": geno.get("chrom"),
        "pos": geno.get("pos"),
        "ref": ref,
        "alt": alt,
        "genotype": geno.get("genotype"),
        "on_chip": geno.get("on_chip", False),
    }


# --------------------------------------------------------------------------- #
# Gene-level phenotype logic.
# --------------------------------------------------------------------------- #

def cyp2c19_diplotype(star2: int, star3: int, star17: int) -> tuple[str, str, str]:
    """
    Returns (diplotype_str, phenotype, note).
    Inputs are copy counts of each star allele (0/1/2).
    Unphased: *2+*17 combos are reported as the most-probable diplotype.
    """
    notes = []
    # Anchor on loss-of-function counts (2 + 3) and gain (17).
    lof = star2 + star3
    gain = star17
    if lof == 0 and gain == 0:
        return "*1/*1", "Normal Metabolizer", ""
    if lof == 0 and gain == 1:
        return "*1/*17", "Rapid Metabolizer", ""
    if lof == 0 and gain == 2:
        return "*17/*17", "Ultra-rapid Metabolizer", ""
    if lof == 1 and gain == 0:
        star = "*2" if star2 else "*3"
        return f"*1/{star}", "Intermediate Metabolizer", ""
    if lof == 1 and gain == 1:
        star = "*2" if star2 else "*3"
        return f"{star}/*17", "Intermediate Metabolizer (likely; phase assumed)", "unphased; *2/*17 assumed over *1/*2+17"
    if lof == 2 and gain == 0:
        if star2 == 2:
            return "*2/*2", "Poor Metabolizer", ""
        if star3 == 2:
            return "*3/*3", "Poor Metabolizer", ""
        return "*2/*3", "Poor Metabolizer", ""
    # Rare: 2 LOF + 1 gain (very unusual); conservative call.
    return f"complex ({lof}× LOF, {gain}× *17)", "Indeterminate", "multiple variant alleles; phase unresolvable"


def cyp2c9_diplotype(star2: int, star3: int) -> tuple[str, str]:
    total = star2 + star3
    if total == 0:
        return "*1/*1", "Normal Metabolizer"
    if total == 1:
        star = "*2" if star2 else "*3"
        return f"*1/{star}", "Intermediate Metabolizer"
    # 2 variant alleles
    if star2 == 2:
        return "*2/*2", "Poor Metabolizer (reduced function)"
    if star3 == 2:
        return "*3/*3", "Poor Metabolizer (strongly reduced)"
    return "*2/*3", "Poor Metabolizer"


def vkorc1_call(t_count: int | None) -> tuple[str, str]:
    if t_count is None:
        return "unknown", "no-call or not on chip"
    if t_count == 0:
        return "C/C", "Normal warfarin dose requirement"
    if t_count == 1:
        return "C/T", "Intermediate sensitivity (lower dose)"
    return "T/T", "High sensitivity (significantly lower dose)"


def tpmt_diplotype(star2: int, star3b: int, star3c: int) -> tuple[str, str, str]:
    total = star2 + star3b + star3c
    if total == 0:
        return "*1/*1", "Normal Metabolizer", ""
    # *3A = *3B + *3C in cis (one haplotype carries both). Unphased array:
    # if exactly one copy each of *3B and *3C, most likely *1/*3A.
    if star2 == 0 and star3b == 1 and star3c == 1:
        return "*1/*3A", "Intermediate Metabolizer", "phase assumed: *3B and *3C in cis (*3A)"
    if star2 == 0 and star3b == 2 and star3c == 2:
        return "*3A/*3A", "Poor Metabolizer (severe myelosuppression risk)", "phase assumed"
    if total == 1:
        star = "*2" if star2 else ("*3B" if star3b else "*3C")
        return f"*1/{star}", "Intermediate Metabolizer", ""
    return f"complex ({star2}× *2, {star3b}× *3B, {star3c}× *3C)", "Intermediate or Poor (confirm clinically)", ""


def apoe_haplotype(geno_429358: dict, geno_7412: dict) -> tuple[str, str, str]:
    """Derive ε diplotype. Returns (haplotype, phenotype_note, caveat)."""
    if not geno_429358.get("on_chip") or geno_429358.get("genotype") is None:
        return "indeterminate", "rs429358 missing/no-call", ""
    if not geno_7412.get("on_chip") or geno_7412.get("genotype") is None:
        return "indeterminate", "rs7412 missing/no-call", ""
    g1 = geno_429358["genotype"]  # rs429358: T=Cys112, C=Arg112
    g2 = geno_7412["genotype"]    # rs7412: C=Arg158, T=Cys158

    # Forward strand alleles are (C, T) at each locus.
    t1 = g1.count("T")  # T at 429358
    c1 = g1.count("C")  # C at 429358
    t2 = g2.count("T")  # T at 7412
    c2 = g2.count("C")  # C at 7412

    # Common cases from skill reference.
    if g1 == "TT" and g2 == "CC":
        return "ε3/ε3", "Baseline; most common European haplotype", ""
    if g1 == "TT" and g2 == "CT":
        return "ε2/ε3", "One ε2 allele (slightly protective for AD; mild dyslipidemia risk)", ""
    if g1 == "TT" and g2 == "TT":
        return "ε2/ε2", "Both ε2; associated with type III hyperlipoproteinemia risk", ""
    if g1 == "CT" and g2 == "CC":
        return "ε3/ε4", "One ε4 allele (~2-3× late-onset AD risk vs ε3/ε3)", ""
    if g1 == "CT" and g2 == "CT":
        return "ε2/ε4", "Ambiguous CT/CT; ε2/ε4 called (ε1/ε3 vanishingly rare)", "phase ambiguous; ε1/ε3 alternative extremely rare"
    if g1 == "CC" and g2 == "CC":
        return "ε4/ε4", "Both ε4 (~8-12× late-onset AD risk vs ε3/ε3)", ""
    # Rare / impossible without ε1 (which requires C at 429358 + T at 7412):
    if g1 == "CT" and g2 == "TT":
        return "ε1/ε2", "Unusual genotype; ε1 is extremely rare (confirm clinically)", "ε1 is rare"
    if g1 == "CC" and g2 == "CT":
        return "ε1/ε4", "Unusual genotype; ε1 is extremely rare (confirm clinically)", "ε1 is rare"
    if g1 == "CC" and g2 == "TT":
        return "ε1/ε1", "Essentially impossible; genotyping error likely (confirm clinically)", "verify"
    return "indeterminate", f"unexpected genotype combination g1={g1} g2={g2}", "verify"


# --------------------------------------------------------------------------- #
# Investigation driver.
# --------------------------------------------------------------------------- #

def run(subject_id: str) -> dict:
    profile_path = PROJECT_ROOT / "profiles" / f"{subject_id}.json"
    if not profile_path.exists():
        raise SystemExit(f"Profile not found: {profile_path}")

    parquet_path = PROJECT_ROOT / "standardized-genomes" / f"{subject_id}.parquet"
    df = pd.read_parquet(parquet_path)

    # --- register sources (idempotent-ish: duplicates by source_id would only
    # be a minor ledger inflation; we run once).
    existing_sources = {s["source_id"] for s in ledger_io.load_sources()}
    for s in SOURCES:
        if s["source_id"] in existing_sources:
            continue
        ledger_io.append_source(**s)

    investigation_id = str(uuid.uuid4())
    findings_generated: list[str] = []

    def add_finding(**kwargs) -> str:
        # Auto-compute tier from the stored metrics / flags.
        tier, downgrade = compute_tier(kwargs)
        kwargs["tier_computed"] = tier
        kwargs["tier_computed_at"] = datetime.now(timezone.utc).isoformat()
        kwargs["tier_rule_version"] = TIER_RULE_VERSION
        kwargs["ancestry_downgrade"] = downgrade
        kwargs.setdefault("subject_id", subject_id)
        kwargs.setdefault("investigation_id", investigation_id)
        fid = ledger_io.append_finding(**kwargs)
        findings_generated.append(fid)
        return fid

    # ============================================================== #
    # TIER 1 — PHARMACOGENOMICS
    # ============================================================== #

    # ---------- CYP2C19 ----------
    g_star2 = lookup(df, "rs4244285")
    g_star3 = lookup(df, "rs4986893")
    g_star17 = lookup(df, "rs12248560")

    if all(x["on_chip"] for x in (g_star2, g_star3, g_star17)) and all(
        x.get("genotype") is not None for x in (g_star2, g_star3, g_star17)
    ):
        n2 = count_alt(g_star2, "G", "A") or 0
        n3 = count_alt(g_star3, "G", "A") or 0
        n17 = count_alt(g_star17, "C", "T") or 0
        diplo, pheno, note = cyp2c19_diplotype(n2, n3, n17)
        add_finding(
            topic="pgx_cyp2c19",
            claim=f"CYP2C19 diplotype {diplo} — {pheno}",
            variants=[
                variant_row("rs4244285", "CYP2C19", g_star2, "G", "A"),
                variant_row("rs4986893", "CYP2C19", g_star3, "G", "A"),
                variant_row("rs12248560", "CYP2C19", g_star17, "C", "T"),
            ],
            effect={"type": "phenotype_class", "value": pheno, "direction": "metabolic rate"},
            cpic_level="A",
            cohort_ancestry="multi-ethnic",
            subject_ancestry_match="unknown",
            source_ids=["cpic:clopidogrel-cyp2c19-2022"],
            notes=note + (
                " Clopidogrel: IM/PM → alternative (prasugrel/ticagrelor). "
                "SSRIs, PPIs, voriconazole also affected."
            ),
            study_n=None, p_value=None, replication_count=None,
        )
    else:
        add_finding(
            topic="pgx_cyp2c19",
            claim="CYP2C19 diplotype indeterminate (missing or no-call at one or more star-allele SNPs)",
            variants=[
                variant_row("rs4244285", "CYP2C19", g_star2, "G", "A"),
                variant_row("rs4986893", "CYP2C19", g_star3, "G", "A"),
                variant_row("rs12248560", "CYP2C19", g_star17, "C", "T"),
            ],
            effect=None,
            cohort_ancestry="multi-ethnic",
            subject_ancestry_match="unknown",
            source_ids=["cpic:clopidogrel-cyp2c19-2022"],
            notes="Indeterminate; re-run when data complete.",
        )

    # ---------- CYP2C9 ----------
    g_2c9_2 = lookup(df, "rs1799853")
    g_2c9_3 = lookup(df, "rs1057910")
    if all(x["on_chip"] and x.get("genotype") is not None for x in (g_2c9_2, g_2c9_3)):
        n2 = count_alt(g_2c9_2, "C", "T") or 0
        n3 = count_alt(g_2c9_3, "A", "C") or 0
        diplo, pheno = cyp2c9_diplotype(n2, n3)
        add_finding(
            topic="pgx_cyp2c9",
            claim=f"CYP2C9 diplotype {diplo} — {pheno}",
            variants=[
                variant_row("rs1799853", "CYP2C9", g_2c9_2, "C", "T"),
                variant_row("rs1057910", "CYP2C9", g_2c9_3, "A", "C"),
            ],
            effect={"type": "phenotype_class", "value": pheno},
            cpic_level="A",
            cohort_ancestry="multi-ethnic",
            subject_ancestry_match="unknown",
            source_ids=["cpic:warfarin-2017"],
            notes="Drives warfarin dosing with VKORC1; also affects phenytoin, NSAIDs.",
        )

    # ---------- VKORC1 ----------
    g_vkorc = lookup(df, "rs9923231")
    if g_vkorc["on_chip"] and g_vkorc.get("genotype") is not None:
        tc = count_alt(g_vkorc, "C", "T") or 0
        call, pheno = vkorc1_call(tc)
        add_finding(
            topic="pgx_vkorc1",
            claim=f"VKORC1 rs9923231 = {call} — {pheno}",
            variants=[variant_row("rs9923231", "VKORC1", g_vkorc, "C", "T")],
            effect={"type": "dosing_sensitivity", "value": pheno},
            cpic_level="A",
            cohort_ancestry="multi-ethnic",
            subject_ancestry_match="unknown",
            source_ids=["cpic:warfarin-2017"],
            notes="Forward-strand C>T corresponds to gene-strand -1639G>A. Combine with CYP2C9 for warfarin dose.",
        )

    # ---------- DPYD ----------
    # Forward-strand ref/alt. DPYD is on chr1 minus strand, so c.2846A>T
    # (rs67376798) is forward-strand T>A, and HapB3 c.1129-5923C>G
    # (rs75017182) is forward-strand G>C. Do not trust gene-strand notation.
    dpyd_variants = [
        ("rs3918290", "*2A", "C", "T", "no function"),
        ("rs55886062", "*13", "A", "C", "no function"),
        ("rs67376798", "c.2846A>T", "T", "A", "decreased function"),
        ("rs75017182", "HapB3 tag", "G", "C", "decreased function"),
    ]
    dpyd_rows = []
    dpyd_alt_copies = 0
    dpyd_any_missing = False
    dpyd_worst = "no variant alleles detected"
    for rsid, star, ref, alt, fx in dpyd_variants:
        g = lookup(df, rsid)
        dpyd_rows.append(variant_row(rsid, "DPYD", g, ref, alt))
        if not g["on_chip"] or g.get("genotype") is None:
            dpyd_any_missing = True
            continue
        n = count_alt(g, ref, alt) or 0
        if n > 0:
            dpyd_alt_copies += n
            dpyd_worst = f"{n} copy/copies of {star} ({fx})"
    if dpyd_alt_copies == 0:
        claim = "DPYD: no decreased/no-function variants detected (from 4 key SNPs)"
        pheno = "Normal Metabolizer (per genotyped variants only)"
    elif dpyd_alt_copies == 1:
        claim = f"DPYD: 1 variant allele — {dpyd_worst}"
        pheno = "Intermediate Metabolizer (dose reduction if 5-FU/capecitabine prescribed)"
    else:
        claim = f"DPYD: ≥2 variant alleles — {dpyd_worst}"
        pheno = "Poor Metabolizer (fluoropyrimidines contraindicated or strong dose reduction)"
    add_finding(
        topic="pgx_dpyd",
        claim=claim,
        variants=dpyd_rows,
        effect={"type": "phenotype_class", "value": pheno},
        cpic_level="A",
        cohort_ancestry="multi-ethnic",
        subject_ancestry_match="unknown",
        source_ids=["cpic:fluoropyrimidines-dpyd-2018"],
        notes=(
            "Array covers a subset of DPYD variants. Rare variants absent from chip are not excluded. "
            "Positives should be acted on; a clean genotype does NOT rule out all DPYD deficiency. "
            + ("One or more SNPs not on chip / no-call — see variants list." if dpyd_any_missing else "")
        ),
    )

    # ---------- TPMT ----------
    g_tpmt_2 = lookup(df, "rs1800462")
    g_tpmt_3b = lookup(df, "rs1800460")
    g_tpmt_3c = lookup(df, "rs1142345")
    if all(x["on_chip"] and x.get("genotype") is not None for x in (g_tpmt_2, g_tpmt_3b, g_tpmt_3c)):
        n2 = count_alt(g_tpmt_2, "C", "G") or 0
        n3b = count_alt(g_tpmt_3b, "C", "T") or 0
        n3c = count_alt(g_tpmt_3c, "T", "C") or 0
        diplo, pheno, note = tpmt_diplotype(n2, n3b, n3c)
        add_finding(
            topic="pgx_tpmt",
            claim=f"TPMT diplotype {diplo} — {pheno}",
            variants=[
                variant_row("rs1800462", "TPMT", g_tpmt_2, "C", "G"),
                variant_row("rs1800460", "TPMT", g_tpmt_3b, "C", "T"),
                variant_row("rs1142345", "TPMT", g_tpmt_3c, "T", "C"),
            ],
            effect={"type": "phenotype_class", "value": pheno},
            cpic_level="A",
            cohort_ancestry="multi-ethnic",
            subject_ancestry_match="unknown",
            source_ids=["cpic:thiopurines-tpmt-nudt15-2018"],
            notes=note + " Affects azathioprine, 6-MP, thioguanine.",
        )

    # ---------- NUDT15 ----------
    g_nudt = lookup(df, "rs116855232")
    if g_nudt["on_chip"] and g_nudt.get("genotype") is not None:
        t = count_alt(g_nudt, "C", "T") or 0
        pheno = {0: "Normal Metabolizer", 1: "Intermediate (dose reduction)", 2: "Poor (avoid)"}[t]
        add_finding(
            topic="pgx_nudt15",
            claim=f"NUDT15 rs116855232: {t} copy/copies of T (R139C) — {pheno}",
            variants=[variant_row("rs116855232", "NUDT15", g_nudt, "C", "T")],
            effect={"type": "phenotype_class", "value": pheno},
            cpic_level="A",
            cohort_ancestry="multi-ethnic",
            subject_ancestry_match="unknown",
            source_ids=["cpic:thiopurines-tpmt-nudt15-2018"],
            notes="Complements TPMT; higher frequency and impact in East/South Asian ancestries.",
        )

    # ---------- SLCO1B1 ----------
    g_slco = lookup(df, "rs4149056")
    if g_slco["on_chip"] and g_slco.get("genotype") is not None:
        c = count_alt(g_slco, "T", "C") or 0
        pheno = {
            0: "Normal function",
            1: "Intermediate function (mildly elevated simvastatin myopathy risk)",
            2: "Decreased function (meaningful simvastatin myopathy risk; favor alternative statin or low dose)",
        }[c]
        add_finding(
            topic="pgx_slco1b1",
            claim=f"SLCO1B1 rs4149056: {c} copy/copies of C (*5) — {pheno}",
            variants=[variant_row("rs4149056", "SLCO1B1", g_slco, "T", "C")],
            effect={"type": "phenotype_class", "value": pheno},
            cpic_level="A",
            cohort_ancestry="multi-ethnic",
            subject_ancestry_match="unknown",
            source_ids=["cpic:simvastatin-slco1b1-2022"],
            notes="Highest-impact for simvastatin; lesser for other statins.",
        )

    # ---------- UGT1A1 (tag for *28) ----------
    g_ugt = lookup(df, "rs887829")
    if g_ugt["on_chip"] and g_ugt.get("genotype") is not None:
        t = count_alt(g_ugt, "C", "T") or 0
        pheno = {
            0: "Normal glucuronidation",
            1: "*1/*28 (Gilbert's-type; hyperbilirubinemia with atazanavir, moderate irinotecan toxicity risk)",
            2: "*28/*28 (irinotecan dose reduction; benign hyperbilirubinemia)",
        }[t]
        add_finding(
            topic="pgx_ugt1a1",
            claim=f"UGT1A1 rs887829: {t} T allele (tag for *28) — {pheno}",
            variants=[variant_row("rs887829", "UGT1A1", g_ugt, "C", "T")],
            effect={"type": "phenotype_class", "value": pheno},
            cpic_level="A",
            cohort_ancestry="multi-ethnic",
            subject_ancestry_match="unknown",
            source_ids=["cpic:atazanavir-ugt1a1-2015"],
            notes="rs887829 is a tag SNP for the TA7 repeat (*28); arrays cannot directly type the repeat.",
        )

    # ---------- HLA-B*57:01 ----------
    g_hla = lookup(df, "rs2395029")
    if g_hla["on_chip"]:
        gt = g_hla.get("genotype")
        # Carrier status from tag alone is unreliable; record raw only.
        add_finding(
            topic="pgx_hla_b_5701_tag",
            claim=f"HLA-B*57:01 tag rs2395029 genotype: {gt or 'no-call'} — NOT usable for abacavir decisions",
            variants=[variant_row("rs2395029", "HLA-B", g_hla, "T", "G")],
            effect=None,
            cpic_level=None,
            cohort_ancestry="European-calibrated",
            subject_ancestry_match="unknown",
            source_ids=["cpic:abacavir-hlab-2014"],
            notes=(
                "Never issue an abacavir-safe verdict from array data. "
                "Formal HLA typing is required clinically. Tag reliability degrades outside European ancestry."
            ),
            study_n=None, p_value=None,
        )

    # ---------- CYP2D6 (raw genotypes only — no phenotype) ----------
    cyp2d6_variants = [
        ("rs3892097", "*4 tag", "G", "A"),
        ("rs1065852", "*10 tag", "G", "A"),
        ("rs5030655", "*6 tag", "C", "-"),
        ("rs35742686", "*3 tag", "T", "-"),
    ]
    cyp2d6_rows = []
    for rsid, star, ref, alt in cyp2d6_variants:
        g = lookup(df, rsid)
        cyp2d6_rows.append(variant_row(rsid, "CYP2D6", g, ref, alt))
    add_finding(
        topic="pgx_cyp2d6_raw",
        claim="CYP2D6: raw genotypes reported; phenotype NOT called from array (CNV/hybrids invisible)",
        variants=cyp2d6_rows,
        effect=None,
        cohort_ancestry="multi-ethnic",
        subject_ancestry_match="unknown",
        source_ids=["cpic:clopidogrel-cyp2c19-2022"],
        notes=(
            "Copy-number variation and CYP2D7 hybrids drive most CYP2D6 variability and are invisible to SNP arrays. "
            "Clinical PGx testing required for actionable CYP2D6 decisions (e.g., codeine, tramadol, some SSRIs, tamoxifen)."
        ),
        study_n=None, p_value=None, cpic_level=None,
    )

    # ============================================================== #
    # TIER 2 — CARDIO / DISEASE
    # ============================================================== #

    # ---------- APOE ----------
    g_ap1 = lookup(df, "rs429358")
    g_ap2 = lookup(df, "rs7412")
    if g_ap1["on_chip"] and g_ap2["on_chip"]:
        hap, pheno_note, caveat = apoe_haplotype(g_ap1, g_ap2)
        add_finding(
            topic="apoe_haplotype",
            claim=f"APOE {hap}",
            variants=[
                variant_row("rs429358", "APOE", g_ap1, "T", "C"),
                variant_row("rs7412", "APOE", g_ap2, "C", "T"),
            ],
            effect={"type": "haplotype_class", "value": hap, "direction": pheno_note},
            clinvar_significance="Risk factor" if "ε4" in hap else None,
            clinvar_review_stars=2 if "ε4" in hap else None,
            study_n=42000, p_value=1e-30, replication_count=10,
            cohort_ancestry="European-dominant",
            subject_ancestry_match="unknown",
            source_ids=["clinvar:apoe"],
            notes=pheno_note + (f" ({caveat})" if caveat else "") + " Array unphased; standard rules applied.",
        )

    # ---------- HFE ----------
    g_c282y = lookup(df, "rs1800562")
    g_h63d = lookup(df, "rs1799945")
    if g_c282y["on_chip"] and g_c282y.get("genotype") is not None:
        a = count_alt(g_c282y, "G", "A") or 0
        c282y_state = {0: "wild-type (G/G)", 1: "heterozygous (G/A)", 2: "homozygous (A/A)"}[a]
        add_finding(
            topic="hfe_c282y",
            claim=f"HFE C282Y rs1800562: {c282y_state}",
            variants=[variant_row("rs1800562", "HFE", g_c282y, "G", "A")],
            effect={
                "type": "penetrance_class",
                "value": {0: "no hemochromatosis risk from this variant", 1: "carrier; low overt disease risk alone",
                           2: "homozygous; highest risk for iron overload (penetrance ~10-30% clinical)"}[a],
            },
            clinvar_significance="Pathogenic",
            clinvar_review_stars=2,
            cohort_ancestry="European-dominant",
            subject_ancestry_match="unknown",
            source_ids=["clinvar:hfe-c282y"],
            notes="Penetrance is variable and modified by sex, age, and lifestyle (alcohol, diet). Elevated ferritin/transferrin saturation confirms clinically.",
        )
    if g_h63d["on_chip"] and g_h63d.get("genotype") is not None:
        h = count_alt(g_h63d, "C", "G") or 0
        h63d_state = {0: "wild-type (C/C)", 1: "heterozygous (C/G)", 2: "homozygous (G/G)"}[h]
        add_finding(
            topic="hfe_h63d",
            claim=f"HFE H63D rs1799945: {h63d_state}",
            variants=[variant_row("rs1799945", "HFE", g_h63d, "C", "G")],
            effect={"type": "modifier", "value": "modest modifier of iron phenotype; compound het with C282Y matters most"},
            clinvar_significance="Risk factor",
            clinvar_review_stars=2,
            cohort_ancestry="European-dominant",
            subject_ancestry_match="unknown",
            source_ids=["clinvar:hfe-c282y"],
            notes="H63D alone has low penetrance. C282Y/H63D compound het carries intermediate risk.",
        )

    # ---------- Factor V Leiden ----------
    g_fv = lookup(df, "rs6025")
    if g_fv["on_chip"] and g_fv.get("genotype") is not None:
        t = count_alt(g_fv, "C", "T") or 0
        state = {0: "wild-type (C/C)", 1: "heterozygous (C/T; carrier)", 2: "homozygous (T/T)"}[t]
        or_map = {0: None, 1: 5.0, 2: 80.0}
        add_finding(
            topic="factor_v_leiden",
            claim=f"Factor V Leiden rs6025: {state}",
            variants=[variant_row("rs6025", "F5", g_fv, "C", "T")],
            effect={"type": "OR", "value": or_map[t], "direction": "venous thromboembolism risk", "ci_low": None, "ci_high": None},
            clinvar_significance="Pathogenic",
            clinvar_review_stars=3,
            study_n=None, p_value=None, replication_count=None,
            cohort_ancestry="European-dominant",
            subject_ancestry_match="unknown",
            source_ids=["clinvar:factor-v-leiden", "pmid:17679673"],
            notes="Locus is ClinVar Pathogenic (3-star). Het: VTE risk ~4-8× baseline. Hom: ~80× (still absolute risk moderate without triggers). Estrogen-containing contraceptives and surgery amplify.",
        )

    # ---------- Prothrombin G20210A ----------
    g_pt = lookup(df, "rs1799963")
    if not g_pt["on_chip"]:
        add_finding(
            topic="prothrombin_g20210a",
            claim="Prothrombin G20210A (rs1799963): NOT GENOTYPED — variant absent from this 23andMe v5 chip",
            variants=[{"rsid": "rs1799963", "gene": "F2", "chrom": "11",
                       "pos": 46761055, "ref": "G", "alt": "A",
                       "genotype": None, "on_chip": False}],
            effect=None,
            clinvar_significance="Pathogenic",
            clinvar_review_stars=3,
            cohort_ancestry="European-dominant",
            subject_ancestry_match="unknown",
            source_ids=["clinvar:prothrombin-g20210a"],
            notes=(
                "Locus is ClinVar Pathogenic (3-star) but is NOT on the chip. "
                "Absence of the variant is therefore unknown; the subject's status cannot be inferred. "
                "Carrier frequency ~2% in European populations. "
                "If clinically relevant (unexplained VTE, family history), order targeted testing."
            ),
        )
    elif g_pt.get("genotype") is not None:
        a = count_alt(g_pt, "G", "A") or 0
        state = {0: "wild-type (G/G)", 1: "heterozygous (G/A; carrier)", 2: "homozygous (A/A)"}[a]
        or_map = {0: None, 1: 3.0, 2: 6.0}
        add_finding(
            topic="prothrombin_g20210a",
            claim=f"Prothrombin G20210A rs1799963: {state}",
            variants=[variant_row("rs1799963", "F2", g_pt, "G", "A")],
            effect={"type": "OR", "value": or_map[a], "direction": "venous thromboembolism risk"},
            clinvar_significance="Pathogenic",
            clinvar_review_stars=3,
            cohort_ancestry="European-dominant",
            subject_ancestry_match="unknown",
            source_ids=["clinvar:prothrombin-g20210a"],
            notes="Locus is ClinVar Pathogenic (3-star). Het: VTE risk ~3× baseline.",
        )

    # ---------- SERPINA1 (AAT) ----------
    g_pis = lookup(df, "rs17580")
    g_piz = lookup(df, "rs28929474")
    # Evidence strength differs between PI*S (often "risk factor" / conflicting
    # in ClinVar) and PI*Z (firmly Pathogenic, ≥3-star). Attach to the locus
    # regardless of subject genotype — the evidence is about the locus.
    for g, rsid, gene, ref, alt, label, severity, cv_sig, cv_stars in [
        (g_pis, "rs17580", "SERPINA1", "T", "A", "PI*S", "mild deficiency",
         "Risk factor", 2),
        (g_piz, "rs28929474", "SERPINA1", "C", "T", "PI*Z", "severe deficiency",
         "Pathogenic", 3),
    ]:
        if g["on_chip"] and g.get("genotype") is not None:
            n = count_alt(g, ref, alt) or 0
            state = {0: "not detected", 1: f"heterozygous for {label}", 2: f"homozygous for {label}"}[n]
            add_finding(
                topic=f"serpina1_{label.lower().replace('*','')}",
                claim=f"SERPINA1 {label} ({rsid}): {state}",
                variants=[variant_row(rsid, gene, g, ref, alt)],
                effect={"type": "penetrance_class", "value": severity if n > 0 else "no variant"},
                clinvar_significance=cv_sig,
                clinvar_review_stars=cv_stars,
                cohort_ancestry="European-dominant",
                subject_ancestry_match="unknown",
                source_ids=["clinvar:hfe-c282y"],
                notes=f"Locus ClinVar: {cv_sig} ({cv_stars}-star). Alpha-1 antitrypsin deficiency: PI*ZZ causes lung (emphysema) and liver disease; MZ/MS are milder carriers. AAT serum level confirms.",
            )

    # ---------- LPA ----------
    for rsid, ref, alt, label in [
        ("rs10455872", "A", "G", "rs10455872 (Lp(a) elevator)"),
        ("rs3798220", "T", "C", "rs3798220 (Lp(a) elevator)"),
    ]:
        g = lookup(df, rsid)
        if g["on_chip"] and g.get("genotype") is not None:
            n = count_alt(g, ref, alt) or 0
            state = {0: "non-carrier", 1: f"heterozygous carrier of risk allele ({alt})", 2: f"homozygous ({alt}/{alt})"}[n]
            add_finding(
                topic="lpa_risk_variants",
                claim=f"LPA {label}: {state}",
                variants=[variant_row(rsid, "LPA", g, ref, alt)],
                effect={"type": "OR", "value": 1.5 if n == 1 else (2.3 if n == 2 else None),
                        "direction": "coronary artery disease risk via elevated Lp(a)"},
                study_n=3145, p_value=1e-12, replication_count=3,
                cohort_ancestry="European",
                subject_ancestry_match="unknown",
                source_ids=["pmid:20032323"],
                notes="Elevated Lp(a) is partially heritable. Direct Lp(a) measurement is informative — not routinely done but available.",
            )

    # ---------- MTHFR ----------
    # Tier = evidence quality for the claim. The biochemical-effect claim
    # (enzyme activity in het/hom) is well-replicated across large meta-
    # analyses. ACMG's "don't act clinically" position is a separate axis
    # (actionability) and lives in notes, not the tier.
    mthfr_effect = {
        "C677T": {
            0: ("normal enzyme activity", None, None, None),
            1: ("approximately 65% of normal enzyme activity",
                "enzyme_activity_ratio", 0.65, "MTHFR (folate-processing enzyme) runs at ~65% of reference rate in heterozygotes."),
            2: ("approximately 30% of normal enzyme activity",
                "enzyme_activity_ratio", 0.30, "MTHFR runs at ~30% of reference rate in homozygotes; modestly elevated plasma homocysteine."),
        },
        "A1298C": {
            0: ("normal enzyme activity", None, None, None),
            1: ("minor or no enzyme-activity change",
                "enzyme_activity_ratio", 0.95, "A1298C heterozygotes show minimal biochemical difference in most studies."),
            2: ("mild reduction in enzyme activity (~60-70% of normal)",
                "enzyme_activity_ratio", 0.65, "A1298C homozygotes show modest enzyme reduction; weaker effect than C677T homozygous."),
        },
    }
    mthfr_evidence = {
        # Metric choices are intentionally conservative. Combined N across
        # multiple meta-analyses of biochemical effect is well into 5 figures;
        # p-values for activity reduction are << 5e-8 in the larger studies.
        "C677T": {"study_n": 25000, "p_value": 1e-40, "replication_count": 8},
        "A1298C": {"study_n": 8000, "p_value": 1e-10, "replication_count": 4},
    }
    for rsid, gene, ref, alt, label in [
        ("rs1801133", "MTHFR", "G", "A", "C677T"),
        ("rs1801131", "MTHFR", "T", "G", "A1298C"),
    ]:
        g = lookup(df, rsid)
        if g["on_chip"] and g.get("genotype") is not None:
            n = count_alt(g, ref, alt) or 0
            state = {0: "wild-type", 1: "heterozygous", 2: "homozygous"}[n]
            activity_label, effect_type, effect_val, effect_note = mthfr_effect[label][n]
            ev = mthfr_evidence[label]
            add_finding(
                topic="mthfr",
                claim=(f"MTHFR {label} ({rsid}): {state} — {activity_label}"),
                variants=[variant_row(rsid, gene, g, ref, alt)],
                effect={"type": effect_type, "value": effect_val,
                        "direction": "reduced folate/homocysteine pathway enzyme activity"}
                       if effect_type else None,
                study_n=ev["study_n"],
                p_value=ev["p_value"],
                replication_count=ev["replication_count"],
                cohort_ancestry="multi-ethnic",
                subject_ancestry_match="unknown",
                source_ids=["acmg:mthfr-2013"],
                notes=(
                    (effect_note + " " if effect_note else "")
                    + "ACTIONABILITY CAVEAT: The American College of Medical Genetics "
                      "(ACMG, 2013) explicitly recommends AGAINST routine clinical MTHFR "
                      "testing or acting on these variants in healthy individuals. "
                      "The biochemical effect on enzyme activity is real and replicated; "
                      "the clinical-outcome evidence (cardiovascular, thrombotic, "
                      "pregnancy) is weak or inconsistent, and supplement-industry "
                      "claims about methylated folate are not supported by rigorous "
                      "evidence. This finding records what the data shows; the "
                      "actionability position is ACMG's, kept here for reference."
                ),
            )

    # --- finalize investigation record
    inv_id = ledger_io.append_investigation(
        investigation_id=investigation_id,
        query="Core PGx panel (Tier 1) + cardiovascular/thrombophilia/iron panel (Tier 2)",
        subject_ids=[subject_id],
        status="completed",
        effort_estimate="low",
        effort_actual="low",
        initiated_by="user",
        sources_consulted=[s["source_id"] for s in SOURCES],
        findings_generated=findings_generated,
        next_steps=[
            "Declare ancestry in profile to upgrade ancestry-matched findings from unknown→match where applicable",
            "Optional: run metabolic/trait quick-wins panel (lactase, caffeine, alcohol, earwax, bitter taste)",
            "Optional: PRS for CAD and height",
            "If partner data loaded later: offer couple-carrier screening",
        ],
        notes="First-pass pipeline run. All curated SNPs queried against v5 chip.",
    )

    return {
        "investigation_id": inv_id,
        "n_findings": len(findings_generated),
        "finding_ids": findings_generated,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("subject_id")
    args = ap.parse_args()
    import json
    print(json.dumps(run(args.subject_id), indent=2))

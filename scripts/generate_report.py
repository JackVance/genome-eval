"""Generate a markdown report from the active ledger findings.

Applies Rules 1-11 from CLAUDE.md / SKILL.md: headers state category AND
result, plain-English leads, confidence visible, deviations get full
treatment, reference findings get one-liners, compound traits in full,
genomic inference kept separate from self-report, US-customary units.

Usage:
    python scripts/generate_report.py alice
    python scripts/generate_report.py alice --out reports/2026-04-19.md
    python scripts/generate_report.py alice --topic prs    # filter by topic prefix
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io


def _cat_for(topic: str, effect: dict, evidence_class: str | None) -> str:
    """Assign each finding to a presentation category so related items group.

    Order (by category key) determines section order in the report.
    """
    t = topic or ""
    e_type = (effect or {}).get("type", "")
    if t.startswith("carrier_"):
        return "2_carrier"
    if t.startswith("pgx_"):
        return "1_pgx"
    if t in ("cftr_carrier_i3000001", "cftr_carrier_panel_summary"):
        return "2_carrier"
    if "haplogroup" in t:
        return "6_haplogroup"
    if e_type == "prs" or t.startswith("prs_"):
        return "7_prs"
    if t.startswith("apoe") or t.startswith("factor_v") or t.startswith("hfe") \
            or t.startswith("prothrombin") or t.startswith("mthfr") \
            or t.startswith("serpina") or t.startswith("lpa"):
        return "3_disease_risk"
    if t.startswith("trait_"):
        return "4_trait"
    if evidence_class == "not_callable_from_array":
        return "9_not_callable"
    return "8_other"


SECTION_TITLES = {
    "1_pgx": "Pharmacogenomics",
    "2_carrier": "Carrier screening",
    "3_disease_risk": "Disease risk (single/few variant)",
    "4_trait": "Traits",
    "6_haplogroup": "Deep ancestry (haplogroups)",
    "7_prs": "Polygenic scores",
    "8_other": "Other",
    "9_not_callable": "Not callable from this chip",
}

CATEGORY_ORDER = [
    "2_carrier", "1_pgx", "3_disease_risk", "4_trait",
    "6_haplogroup", "7_prs", "9_not_callable", "8_other",
]


def _format_header(row: dict) -> str:
    """Produce a Rule-1 style header: 'Category: result *(technical ids)*'.

    The ledger's `claim` already contains an informative one-line summary;
    use that as the header body when the topic doesn't map to a better
    human-readable category label.
    """
    topic = row.get("topic", "")
    claim = (row.get("claim") or "").strip()
    return f"### {claim}"


def _confidence_lines(row: dict) -> list[str]:
    """Emit evidence tier + genotype-call / phenotype-confidence lines."""
    tier = row.get("tier_computed", "?")
    inf = row.get("inference_confidence", "?")
    evidence_class = row.get("evidence_class", "?")
    out = [
        f"- **Evidence tier:** {tier} *(evidence class: {evidence_class})*",
        f"- **Inference confidence:** {inf}",
    ]
    # For PRS, surface the self-report cross-check as its own line.
    ccheck = (row.get("effect") or {}).get("self_report_crosscheck") or {}
    if ccheck:
        status = ccheck.get("status", "?").upper()
        out.append(
            f"- **Self-report cross-check:** {status}"
            + (f" — gap {ccheck['gap']:+.2f} ({ccheck.get('gap_residual_sd', 0):+.2f} residual SDs)"
               if ccheck.get("gap") is not None else "")
        )
    return out


def _variant_summary(row: dict) -> str:
    """Compact rsid/genotype list for the footer of each finding block."""
    variants = row.get("variants") or []
    if not variants:
        return ""
    parts = []
    for v in variants[:10]:
        rsid = v.get("rsid") or v.get("pgs_id") or v.get("gene", "")
        call = (v.get("genotype") or v.get("call") or "")
        if v.get("dosage") is not None:
            call = f"{call} dosage={v['dosage']}"
        if rsid and call:
            parts.append(f"`{rsid}`={call}")
        elif rsid:
            parts.append(f"`{rsid}`")
    suffix = " …" if len(variants) > 10 else ""
    return ", ".join(parts) + suffix


def _render_finding(row: dict, is_reference: bool) -> str:
    header = _format_header(row)
    lines: list[str] = [header, ""]

    # Reference / unremarkable findings get a one-liner per Rule 6.
    if is_reference:
        tier = row.get("tier_computed", "?")
        vs = _variant_summary(row)
        return f"- **[{tier}]** {(row.get('claim') or '').strip()} " + (f"*({vs})*" if vs else "")

    # Full treatment for deviations/positive findings.
    notes = (row.get("notes") or "").strip()
    if notes:
        lines.append(notes)
        lines.append("")

    lines.extend(_confidence_lines(row))
    vs = _variant_summary(row)
    if vs:
        lines.append(f"- **Variants:** {vs}")
    lines.append("")
    return "\n".join(lines)


def _is_reference_finding(row: dict) -> bool:
    """Heuristic: a finding is 'reference/unremarkable' if its claim says so.

    Wild-type, non-carrier, normal-metabolizer etc. get one-liner treatment.
    Carrier/affected/heterozygous/mismatch findings get full treatment.
    """
    claim = (row.get("claim") or "").lower()
    if any(tok in claim for tok in [
        "wild-type", "wild type", "non-carrier", "non carrier",
        "not detected", "normal metabolizer", "no pathogenic variants",
        "not genotyped",  # variant absent from chip
    ]):
        return True
    if any(tok in claim for tok in [
        "carrier", "heterozygous", "homozygous", "intermediate",
        "elevated", "reduced", "mismatch",
    ]):
        return False
    return False  # default to full treatment


def generate(subject_id: str, out_path: Path, topic_filter: str | None = None) -> None:
    active = ledger_io.load_active_findings(subject_id=subject_id)
    if topic_filter:
        active = [r for r in active if topic_filter in r.get("topic", "")]

    by_cat: dict[str, list[dict]] = {}
    for r in active:
        cat = _cat_for(r.get("topic", ""), r.get("effect") or {}, r.get("evidence_class"))
        by_cat.setdefault(cat, []).append(r)

    # Subject profile context
    profile_path = PROJECT_ROOT / "profiles" / f"{subject_id}.json"
    profile = json.loads(profile_path.read_text()) if profile_path.exists() else {}
    ancestry = profile.get("declared_ancestry") or {}
    sr = profile.get("self_reported_phenotypes") or {}

    today = date.today().isoformat()
    md: list[str] = [
        f"# Genome report — {profile.get('display_name', subject_id)}",
        "",
        f"_Report generated {today}. Subject ID: `{subject_id}`._",
        "",
        "## Subject profile",
        "",
        f"- **Declared ancestry:** {ancestry.get('top_level', 'unknown')}",
    ]
    if ancestry.get("breakdown"):
        for k, v in ancestry["breakdown"].items():
            pct = v.get("percent_of_total") if isinstance(v, dict) else None
            if pct is not None:
                md.append(f"  - {k.replace('_', ' ')}: {pct}%")
    md.append(f"- **Provider/chip:** {profile.get('provider', '?')} {profile.get('chip_version', '')}".rstrip())
    md.append(f"- **Active findings:** {len(active)}")
    if topic_filter:
        md.append(f"- **Topic filter:** `{topic_filter}`")
    md.append("")

    if sr:
        md.extend([
            "## Self-reported phenotypes on file",
            "",
            "These are subject self-reports, stored separately from genomic inference per Rule 11.",
            "",
        ])
        for key, val in sr.items():
            disp = val.get("value_imperial") or val.get("value") or val.get("value_approx") \
                   or val.get("value_years") or val.get("value_cm")
            md.append(f"- **{key.replace('_', ' ').title()}:** {disp} "
                      f"(source: {val.get('source', '?')}, confidence: {val.get('confidence', '?')})")
        md.append("")

    # Emit sections in canonical order
    for cat in CATEGORY_ORDER:
        rows = by_cat.get(cat, [])
        if not rows:
            continue
        md.append(f"## {SECTION_TITLES[cat]}")
        md.append("")

        # Split reference vs deviation within each category
        references: list[dict] = []
        deviations: list[dict] = []
        for r in rows:
            (references if _is_reference_finding(r) else deviations).append(r)

        for r in deviations:
            md.append(_render_finding(r, is_reference=False))
        if references:
            md.append("")
            md.append("### Reference / unremarkable (one-liners)")
            md.append("")
            for r in references:
                md.append(_render_finding(r, is_reference=True))
            md.append("")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(md), encoding="utf-8")
    print(f"Wrote {out_path} ({out_path.stat().st_size/1024:.1f} KB)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("subject_id")
    parser.add_argument("--topic", help="Filter by topic substring (e.g., 'prs' or 'carrier')")
    parser.add_argument("--out", help="Output path (default: reports/YYYY-MM-DD-<subject>.md)")
    args = parser.parse_args()

    out = Path(args.out) if args.out else (
        PROJECT_ROOT / "reports" / f"{date.today().isoformat()}-{args.subject_id}.md"
    )
    generate(args.subject_id, out, topic_filter=args.topic)


if __name__ == "__main__":
    main()

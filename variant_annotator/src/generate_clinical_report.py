"""
generate_report.py

Section 4 of the VCF Annotator project: clinical filtering + automated
HTML summary report.

Pipeline:
    1. Parse VCF -> Variant objects              (custom_vcf_parser)
    2. Split multiallelic sites into single alleles
    3. Batch-query Ensembl VEP for consequence + gene/transcript
    4. Query MyVariant.info for ClinVar / dbSNP / gnomAD data
                                                   (pathogenicity_scoring)
    5. Filter: gnomAD population AF < threshold AND predicted
       High/Moderate impact consequence
    6. Render a self-contained HTML report (Jinja2) with outbound
       links to PubMed and ClinVar for each candidate variant

Usage:
    python generate_report.py path/to/variants.vcf
    python generate_report.py path/to/variants.vcf --out report.html --af-threshold 0.01

Requires: requests, jinja2  (pip install requests jinja2 --break-system-packages)
"""

import argparse
import sys
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import requests
from jinja2 import Template

from custom_vcf_parser import Variant, parse_vcf, split_multiallelic
from variant_annotation import HEADERS, SERVER, build_vep_region_string, parse_vep_consequence
from pathogenicity_scoring import annotate_variant_pathogenicity

# ---------------------------------------------------------------------------
# Filtering criteria
# ---------------------------------------------------------------------------

HIGH_MODERATE_IMPACTS = {"Nonsense", "Frameshift/Indel", "Missense"}
AF_THRESHOLD = 0.01


# ---------------------------------------------------------------------------
# Data model for one report row
# ---------------------------------------------------------------------------

@dataclass
class ReportRow:
    chrom: str
    pos: int
    ref: str
    alt: str
    rsid: str = "N/A"
    gene: str = "N/A"
    transcript: str = "N/A"
    consequence: str = "Unknown"
    gnomad_af: float = 0.0
    qual: Optional[float] = None
    coverage: Optional[str] = None
    clinvar_sig: str = "None"
    conditions: List[str] = field(default_factory=list)
    pathogenicity: str = "Uncertain / Unannotated"
    clinvar_id: Optional[str] = None

    @property
    def variant_label(self) -> str:
        return f"{self.chrom}:{self.pos} {self.ref}>{self.alt}"

    @property
    def pathogenicity_class(self) -> str:
        """CSS class hook for badge coloring."""
        key = self.pathogenicity.lower()
        if key.startswith("pathogenic"):
            return "badge-path"
        if key.startswith("likely pathogenic"):
            return "badge-likely-path"
        if "uncertain" in key or "vus" in key:
            return "badge-vus"
        if key.startswith("likely benign"):
            return "badge-likely-benign"
        if key.startswith("benign"):
            return "badge-benign"
        return "badge-unknown"

    @property
    def pubmed_url(self) -> str:
        term = self.gene if self.gene not in (None, "N/A") else self.variant_label
        return f"https://pubmed.ncbi.nlm.nih.gov/?term={term}+pathogenic+variant"

    @property
    def clinvar_url(self) -> str:
        if self.clinvar_id:
            return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{self.clinvar_id}/"
        if self.rsid and self.rsid != "N/A":
            return f"https://www.ncbi.nlm.nih.gov/clinvar/?term={self.rsid}"
        return f"https://www.ncbi.nlm.nih.gov/clinvar/?term={self.chrom}%3A{self.pos}"


# ---------------------------------------------------------------------------
# VEP lookup (keeps gene/transcript detail, not just the collapsed impact
# string that variant_annotation.determine_batch_variant_impacts returns)
# ---------------------------------------------------------------------------

def _extract_vep_detail(vep_entry: Dict[str, Any]) -> Dict[str, Any]:
    consequence = parse_vep_consequence(vep_entry)
    transcripts = vep_entry.get("transcript_consequences", [])

    gene, transcript = "N/A", "N/A"
    if transcripts:
        canonical = next(
            (t for t in transcripts if str(t.get("canonical")) == "1"), transcripts[0]
        )
        gene = canonical.get("gene_symbol", "N/A") or "N/A"
        transcript = canonical.get("transcript_id", "N/A") or "N/A"

    return {"consequence": consequence, "gene": gene, "transcript": transcript}


def batch_vep_lookup(variants: List[Variant]) -> Dict[str, Dict[str, Any]]:
    """Batch POST to Ensembl VEP. Returns lookup_key -> vep detail dict."""
    if not variants:
        return {}

    post_hgvs_list, key_mapping = [], {}
    for v in variants:
        region_str, lookup_key = build_vep_region_string(v.chrom, v.pos, v.ref, v.alt)
        post_hgvs_list.append(region_str)
        key_mapping[region_str] = lookup_key

    results: Dict[str, Dict[str, Any]] = {}
    try:
        response = requests.post(
            SERVER + "/vep/homo_sapiens/region",
            headers=HEADERS,
            json={"variants": post_hgvs_list},
            timeout=15,
        )
        if response.ok:
            for item in response.json():
                lookup_key = key_mapping.get(item.get("input"))
                if lookup_key:
                    results[lookup_key] = _extract_vep_detail(item)
        else:
            print(f"[warning] Batch VEP call failed: HTTP {response.status_code}", file=sys.stderr)
    except requests.exceptions.RequestException as e:
        print(f"[warning] Batch VEP connection error: {e}", file=sys.stderr)

    return results


def single_vep_lookup(chrom: str, pos: int, ref: str, alt: str, retries: int = 2) -> Dict[str, Any]:
    """Fallback single-variant GET, used when a variant is missing from the batch response."""
    region_str, _ = build_vep_region_string(chrom, pos, ref, alt)
    ext = f"/vep/homo_sapiens/region/{region_str}?"

    for attempt in range(retries + 1):
        try:
            response = requests.get(SERVER + ext, headers=HEADERS, timeout=10)
            if response.ok:
                data = response.json()
                return _extract_vep_detail(data[0]) if data else {
                    "consequence": "No_Transcript_Overlap", "gene": "N/A", "transcript": "N/A"
                }
            if response.status_code >= 500 and attempt < retries:
                time.sleep(1 + attempt)
                continue
            return {"consequence": f"API_Error_{response.status_code}", "gene": "N/A", "transcript": "N/A"}
        except requests.exceptions.RequestException:
            if attempt == retries:
                return {"consequence": "Connection_Error", "gene": "N/A", "transcript": "N/A"}
            time.sleep(1)

    return {"consequence": "Unknown_Failure", "gene": "N/A", "transcript": "N/A"}


# ---------------------------------------------------------------------------
# Row building
# ---------------------------------------------------------------------------

def build_report_rows(variants: List[Variant]) -> List[ReportRow]:
    """Runs the full VEP + MyVariant.info annotation pipeline and returns
    one ReportRow per variant."""
    vep_results = batch_vep_lookup(variants)
    rows: List[ReportRow] = []

    for v in variants:
        _, lookup_key = build_vep_region_string(v.chrom, v.pos, v.ref, v.alt)
        vep_detail = vep_results.get(lookup_key) or single_vep_lookup(v.chrom, v.pos, v.ref, v.alt)
        print(f"[debug] {v.chrom}:{v.pos} -> {vep_detail}", file=sys.stderr)

        pathogenicity_record = annotate_variant_pathogenicity(
            chrom=v.chrom,
            pos=v.pos,
            ref=v.ref,
            alt=v.alt,
            functional_impact=vep_detail["consequence"],
        )

        coverage = v.info.get("DP")
        if coverage is None and v.sample_data:
            coverage = next(iter(v.sample_data.values())).get("DP")

        gene = vep_detail["gene"]
        if gene == "N/A":
            gene = v.info.get("GENE", "N/A")

        rows.append(
            ReportRow(
                chrom=v.chrom,
                pos=v.pos,
                ref=v.ref,
                alt=v.alt,
                rsid=pathogenicity_record["rsid"],
                gene=gene,
                transcript=vep_detail["transcript"],
                consequence=vep_detail["consequence"],
                gnomad_af=pathogenicity_record["gnomad_af"],
                qual=v.qual,
                coverage=coverage,
                clinvar_sig=pathogenicity_record["clinvar_sig"],
                conditions=pathogenicity_record["conditions"],
                pathogenicity=pathogenicity_record["pathogenicity_summary"],
                clinvar_id=pathogenicity_record["clinvar_id"],
            )
        )

    return rows


def filter_candidate_variants(rows: List[ReportRow], af_threshold: float = AF_THRESHOLD) -> List[ReportRow]:
    """Clinical filtering strategy: population AF < threshold AND
    predicted High/Moderate functional impact."""
    return [r for r in rows if r.gnomad_af < af_threshold and r.consequence in HIGH_MODERATE_IMPACTS]


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

REPORT_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Variant Pathogenicity Report</title>
<style>
  :root {
    --bg: #0f172a; --panel: #1e293b; --text: #e2e8f0; --muted: #94a3b8;
    --border: #334155; --accent: #38bdf8;
  }
  * { box-sizing: border-box; }
  body {
    font-family: "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    background: var(--bg); color: var(--text); margin: 0; padding: 32px;
  }
  h1 { font-size: 22px; margin-bottom: 4px; }
  .subtitle { color: var(--muted); font-size: 13px; margin-bottom: 24px; }
  .stats { display: flex; gap: 16px; margin-bottom: 24px; flex-wrap: wrap; }
  .stat-card {
    background: var(--panel); border: 1px solid var(--border); border-radius: 10px;
    padding: 14px 18px; min-width: 160px;
  }
  .stat-card .num { font-size: 24px; font-weight: 700; color: var(--accent); }
  .stat-card .label { font-size: 12px; color: var(--muted); text-transform: uppercase; letter-spacing: .04em; }
  table { width: 100%; border-collapse: collapse; background: var(--panel); border-radius: 10px; overflow: hidden; }
  th, td { padding: 10px 12px; text-align: left; font-size: 13px; border-bottom: 1px solid var(--border); }
  th { background: #16233a; color: var(--muted); text-transform: uppercase; font-size: 11px; letter-spacing: .04em; }
  tr:hover { background: #24344f; }
  .mono { font-family: "SF Mono", Consolas, monospace; }
  .badge {
    display: inline-block; padding: 3px 9px; border-radius: 999px; font-size: 11px; font-weight: 600;
  }
  .badge-path { background: #7f1d1d; color: #fecaca; }
  .badge-likely-path { background: #7c2d12; color: #fed7aa; }
  .badge-vus { background: #713f12; color: #fde68a; }
  .badge-likely-benign { background: #14532d; color: #bbf7d0; }
  .badge-benign { background: #14532d; color: #bbf7d0; }
  .badge-unknown { background: #334155; color: #cbd5e1; }
  .links a { color: var(--accent); text-decoration: none; margin-right: 10px; font-size: 12px; }
  .links a:hover { text-decoration: underline; }
  .conditions { color: var(--muted); font-size: 12px; }
  .empty { padding: 40px; text-align: center; color: var(--muted); }
  footer { margin-top: 24px; color: var(--muted); font-size: 11px; }
</style>
</head>
<body>
  <h1>Variant Pathogenicity Report</h1>
  <div class="subtitle">
    Source: {{ vcf_path }} &nbsp;|&nbsp; Generated: {{ generated_at }} &nbsp;|&nbsp;
    Filter: gnomAD AF &lt; {{ af_threshold }} AND predicted High/Moderate impact
  </div>

  <div class="stats">
    <div class="stat-card"><div class="num">{{ total_variants }}</div><div class="label">Total variants parsed</div></div>
    <div class="stat-card"><div class="num">{{ candidate_count }}</div><div class="label">Candidate variants</div></div>
    <div class="stat-card"><div class="num">{{ af_threshold }}</div><div class="label">AF threshold</div></div>
  </div>

  {% if rows %}
  <table>
    <thead>
      <tr>
        <th>Variant</th>
        <th>Gene</th>
        <th>Transcript</th>
        <th>Consequence</th>
        <th>gnomAD AF</th>
        <th>Qual</th>
        <th>Coverage</th>
        <th>Clinical significance</th>
        <th>Pathogenicity</th>
        <th>Conditions</th>
        <th>Links</th>
      </tr>
    </thead>
    <tbody>
      {% for r in rows %}
      <tr>
        <td class="mono">{{ r.variant_label }}<br><span class="conditions">{{ r.rsid }}</span></td>
        <td>{{ r.gene }}</td>
        <td class="mono">{{ r.transcript }}</td>
        <td>{{ r.consequence }}</td>
        <td>{{ "%.5f"|format(r.gnomad_af) }}</td>
        <td>{{ r.qual if r.qual is not none else "N/A" }}</td>
        <td>{{ r.coverage if r.coverage else "N/A" }}</td>
        <td>{{ r.clinvar_sig }}</td>
        <td><span class="badge {{ r.pathogenicity_class }}">{{ r.pathogenicity }}</span></td>
        <td class="conditions">{{ r.conditions|join(", ") if r.conditions else "—" }}</td>
        <td class="links">
          <a href="{{ r.pubmed_url }}" target="_blank" rel="noopener">PubMed</a>
          <a href="{{ r.clinvar_url }}" target="_blank" rel="noopener">ClinVar</a>
        </td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
  {% else %}
  <div class="empty">No variants passed the filtering criteria.</div>
  {% endif %}

  <footer>Generated by generate_report.py &middot; Data: Ensembl VEP, MyVariant.info (dbSNP, gnomAD, ClinVar)</footer>
</body>
</html>
"""


def render_html_report(candidate_rows: List[ReportRow], total_variants: int, vcf_path: str, out_path: str, af_threshold: float) -> None:
    template = Template(REPORT_TEMPLATE)
    html = template.render(
        vcf_path=vcf_path,
        generated_at=time.strftime("%Y-%m-%d %H:%M"),
        total_variants=total_variants,
        candidate_count=len(candidate_rows),
        af_threshold=af_threshold,
        rows=candidate_rows,
    )
    with open(out_path, "w") as f:
        f.write(html)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter VCF variants by pathogenicity and generate an HTML summary report."
    )
    parser.add_argument(
        "vcf_path",
        nargs="?",
        default="../data/real_grch38_test.vcf",
        help="Path to input .vcf file (default: ../data/real_grch38_test.vcf)",
    )
    parser.add_argument("--out", default="../output/report.html", help="Output HTML file path")
    parser.add_argument(
        "--af-threshold", type=float, default=AF_THRESHOLD,
        help="Max gnomAD population allele frequency to keep (default: 0.01)",
    )
    args = parser.parse_args()

    raw_variants = parse_vcf(args.vcf_path)
    variants = [sv for v in raw_variants for sv in split_multiallelic(v)]
    print(f"Parsed {len(variants)} single-allele variants from {args.vcf_path}")

    rows = build_report_rows(variants)
    candidates = filter_candidate_variants(rows, args.af_threshold)
    print(
        f"{len(candidates)} candidate variants passed filtering "
        f"(gnomAD AF < {args.af_threshold}, High/Moderate impact)"
    )

    render_html_report(candidates, len(variants), args.vcf_path, args.out, args.af_threshold)
    print(f"Report written to {args.out}")


if __name__ == "__main__":
    main()
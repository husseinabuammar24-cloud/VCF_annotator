# VCF Annotator

A clinical-grade Python tool for annotating VCF (Variant Call Format) files with functional impact predictions and pathogenicity assessments. This tool integrates with **Ensembl VEP**, **MyVariant.info**, **ClinVar**, **dbSNP**, and **gnomAD** to automatically identify high-risk genomic variants.

## What This Does

VCF Annotator takes raw genomic variant data from VCF files and enriches it with:

1. **Functional Impact Classification** — Variants are classified as Synonymous, Missense, Nonsense, or Frameshift using Ensembl VEP
2. **Pathogenicity Scoring** — Integration with ClinVar for known disease associations and clinical significance  
3. **Population Frequency Filtering** — Cross-reference against gnomAD to flag rare, potentially pathogenic variants
4. **Gene & Transcript Mapping** — Links variants to specific genes, transcripts, and HGVS nomenclature
5. **Clinical HTML Report** — Generates a filtered, actionable report with PubMed and ClinVar links

### Stack

- **Language:** Python 3.10+
- **Key Libraries:**
  - `requests` — HTTP queries to Ensembl VEP and MyVariant.info APIs
  - `jinja2` — Templated HTML report generation
  - `dataclasses` — Clean variant data models
- **External APIs:**
  - Ensembl VEP REST API (GRCh38/hg38)
  - MyVariant.info (dbSNP, gnomAD, ClinVar aggregator)

---

## How It's Organized

```
VCF_annotator/
  variant_annotator/
    src/
      custom_vcf_parser.py       # Custom VCF parser (no external deps)
      variant_annotation.py      # Ensembl VEP REST queries
      pathogenicity_scoring.py   # MyVariant.info integration
      generate_clinical_report.py # HTML report generation
    data/
      real_grch38_test.vcf       # Test VCF with real APOE variants
    output/
      report.html                # Generated clinical report (sample output)
  README.md
  vcf_project.md                 # Project requirements & reference
```

### Data Flow

```
VCF file (raw variants)
    ↓
parse_vcf() → Variant objects
    ↓
split_multiallelic() → One variant per ALT allele
    ↓
Ensembl VEP (batch POST) → Gene, transcript, consequence
    ↓
MyVariant.info (parallel queries) → dbSNP, gnomAD AF, ClinVar
    ↓
filter_candidate_variants() → AF < 0.01 AND High/Moderate impact
    ↓
render_html_report() → Clinical summary with links
```

---

## Installation

### Requirements

- Python 3.10 or later
- pip (Python package manager)
- Internet connection (for Ensembl VEP and MyVariant.info APIs)

### Setup

```bash
# Clone the repository
git clone https://github.com/husseinabuammar24-cloud/VCF_annotator.git
cd VCF_annotator

# Install dependencies
pip install requests jinja2
```

> **Note:** System packages are installed with `--break-system-packages` if using system Python. Using a virtual environment is recommended:
> ```bash
> python3 -m venv venv
> source venv/bin/activate
> pip install requests jinja2
> ```

---

## Quick Start

### 1. Generate a Report from a VCF File

```bash
cd variant_annotator/src

# Default: reads from ../data/real_grch38_test.vcf, outputs to ../output/report.html
python generate_clinical_report.py

# Or, specify custom input/output paths
python generate_clinical_report.py path/to/variants.vcf --out report.html --af-threshold 0.01
```

**Output:** A standalone HTML file with:
- Summary statistics (variants parsed, candidates passing filters)
- Filtered table of high-risk variants with gene, consequence, frequencies, pathogenicity badges
- Direct links to PubMed and ClinVar for each variant

### 2. Parse a VCF and Annotate Variants Directly

```python
from custom_vcf_parser import parse_vcf, split_multiallelic
from variant_annotation import determine_batch_variant_impacts

# Load VCF
variants = parse_vcf("data/real_grch38_test.vcf")

# Split multi-allelic sites
single_allele_variants = [v for var in variants for v in split_multiallelic(var)]

# Query Ensembl VEP in batch
impacts = determine_batch_variant_impacts(single_allele_variants)

for v in single_allele_variants:
    lookup_key = f"{v.chrom}:{v.pos}_{v.ref}/{v.alt}"
    impact = impacts.get(lookup_key, "Unknown")
    print(f"{v.chrom}:{v.pos} {v.ref}>{v.alt}: {impact}")
```

---

## Module Reference

### `custom_vcf_parser.py`

**Custom VCF parser built from scratch** — no external bioinformatics libraries.

**Main Functions:**
- `parse_vcf(filepath)` → list of `Variant` objects
- `split_multiallelic(variant)` → list of single-allele `Variant` objects
- `parse_info_field(info_str)` → dict of INFO annotations
- `parse_sample_field(format_str, sample_str)` → dict of genotype fields (GT, DP, AD)

**Variant Dataclass Properties:**
- `is_snv`, `is_insertion`, `is_deletion` — boolean classifiers
- `variant_type` → "SNV", "insertion", "deletion", or "complex"

**Example:**
```python
from custom_vcf_parser import parse_vcf, split_multiallelic

variants = parse_vcf("real_grch38_test.vcf")
for v in variants:
    for single_v in split_multiallelic(v):
        print(f"{single_v.chrom}:{single_v.pos} {single_v.ref}>{single_v.alt} [{single_v.variant_type}]")
        print(f"  Gene: {single_v.info.get('GENE')}")
        print(f"  Samples: {list(single_v.sample_data.keys())}")
```

---

### `variant_annotation.py`

**Queries Ensembl VEP REST API** to determine functional impact (consequence).

**Main Functions:**
- `determine_variant_impact(chrom, pos, ref, alt)` → impact string (single query)
- `determine_batch_variant_impacts(variants)` → dict of impacts (efficient batch POST)
- `build_vep_region_string(chrom, pos, ref, alt)` → VEP-formatted region string
- `parse_vep_consequence(vep_entry)` → high-level impact category

**Impact Categories:**
- `Frameshift/Indel` — frameshift or indel mutations
- `Nonsense` — stop codon gained
- `Missense` — single amino acid change
- `Synonymous` — silent mutation (no amino acid change)
- `Intergenic/Intronic` — non-coding
- `Unknown` — could not be determined

**Example:**
```python
from variant_annotation import determine_batch_variant_impacts
from custom_vcf_parser import parse_vcf, split_multiallelic

variants = [...]  # list of Variant objects
impacts = determine_batch_variant_impacts(variants)

for v in variants:
    key = f"{v.chrom}:{v.pos}_{v.ref}/{v.alt}"
    impact = impacts[key]
    print(f"{v.chrom}:{v.pos}: {impact}")
```

---

### `pathogenicity_scoring.py`

**Integrates MyVariant.info API** to fetch dbSNP IDs, gnomAD frequencies, ClinVar significance, and disease conditions.

**Main Functions:**
- `annotate_variant_pathogenicity(chrom, pos, ref, alt, functional_impact)` → enriched dict
- `convert_variant_to_hgvs(chrom, pos, ref, alt)` → HGVS-g nomenclature
- `extract_gnomad_af(data)` → population allele frequency (float)
- `extract_clinvar_info(data)` → (significance_terms, disease_names)
- `summarize_pathogenicity(clinvar_significances, impact)` → unified classification

**Returned Record Fields:**
```python
{
    "rsid": "rs429358",                    # dbSNP ID
    "gnomad_af": 0.0123,                   # Population frequency
    "clinvar_sig": "Pathogenic",           # ClinVar classification
    "conditions": ["Alzheimers disease"],  # Associated diseases
    "pathogenicity_summary": "Pathogenic", # Final score
    "clinvar_id": "RCV000019456"           # ClinVar accession
}
```

**Example:**
```python
from pathogenicity_scoring import annotate_variant_pathogenicity

result = annotate_variant_pathogenicity(
    chrom="chr19",
    pos=44908684,
    ref="T",
    alt="C",
    functional_impact="Missense"
)
print(f"rsID: {result['rsid']}")
print(f"gnomAD AF: {result['gnomad_af']}")
print(f"ClinVar: {result['clinvar_sig']}")
print(f"Pathogenicity: {result['pathogenicity_summary']}")
```

---

### `generate_clinical_report.py`

**Orchestrates the full pipeline** and generates a clinical HTML report.

**Main Functions:**
- `build_report_rows(variants)` → list of `ReportRow` objects (all annotations)
- `filter_candidate_variants(rows, af_threshold)` → filtered rows (AF < threshold, High/Moderate impact)
- `render_html_report(rows, total, vcf_path, out_path, af_threshold)` → writes HTML

**Filtering Criteria (default):**
- gnomAD allele frequency < 0.01 (rare variants)
- Predicted consequence in {Nonsense, Frameshift/Indel, Missense} (High/Moderate impact)

**HTML Report Features:**
- Dark theme, responsive layout
- Summary statistics (variants parsed, candidates passing filters)
- Sortable table with gene, transcript, consequence, frequencies
- Pathogenicity badges (color-coded: pathogenic, VUS, benign)
- PubMed and ClinVar links for each variant

**CLI Usage:**
```bash
python generate_clinical_report.py [VCF_PATH] --out OUTPUT.html --af-threshold 0.01

# Examples:
python generate_clinical_report.py                                    # defaults
python generate_clinical_report.py my_variants.vcf --out my_report.html
python generate_clinical_report.py my_variants.vcf --af-threshold 0.05
```

---

## Test Data

The repository includes a **real, validated test dataset**:

```
variant_annotator/data/real_grch38_test.vcf
```

Contains two GRCh38-aligned variants in the **APOE gene**:
- `chr19:44908684 T>C` (rs429358) — Missense variant
- `chr19:44908822 C>T` (rs7412) — Missense variant

Both are well-known APOE SNPs linked to Alzheimer's disease risk.

---

## API Dependencies & Rate Limits

### Ensembl VEP REST
- **Endpoint:** `https://rest.ensembl.org/vep/homo_sapiens/region/`
- **Limits:** ~15 requests per second per IP
- **Batch POST:** Supports up to ~200 variants per request (recommended for efficiency)
- **Retry Logic:** Automatic backoff on 5xx errors

### MyVariant.info
- **Endpoint:** `https://myvariant.info/v1/variant/{hgvs_id}`
- **Limits:** ~10 requests per second
- **Coverage:** dbSNP, ClinVar, gnomAD, VEP (aggregated)

### Internet Connection Required
All queries are made at runtime against live external APIs. Offline mode is not supported.

---

## Clinical Filtering Logic

The report applies **strict filtering for clinical relevance**:

1. **Population Frequency:** Keep only variants with gnomAD AF < 0.01 (rare in general population)
2. **Functional Impact:** Keep only High/Moderate impact predictions:
   - ✅ Nonsense (stop codon gained)
   - ✅ Frameshift/Indel (reading frame disrupted)
   - ✅ Missense (amino acid change)
   - ❌ Synonymous (silent, no protein change)
   - ❌ Intergenic/Intronic (non-coding regions)

3. **Pathogenicity Priority:**
   - ClinVar clinical significance (if available) takes precedence
   - Falls back to VEP functional impact prediction
   - Assigns "Uncertain Significance" if both are inconclusive

---

## Lessons Learned (Embedded in Code)

This project encodes several critical **engineering and bioinformatics lessons**:

### 1. Versioning and Reference Context
- Genomic coordinates are version-dependent (GRCh37 vs GRCh38 differ)
- Structural shifts between genome updates break validation
- Always embed assembly/version metadata into data flows

### 2. Interface Layer Selection
- Web domains return unparseable HTML; use dedicated API subdomains
- Choose POST methods for complex array payloads (batch VEP queries)
- Offload strand/orientation logic to API endpoints

### 3. Classification of Error Types
- Distinguish flaky errors (transient infrastructure) from design bugs (faulty logic)
- Identical inputs producing identical errors = code bug; deploy retry loops only for transient errors
- Refuse to retry structural code errors

### 4. Diagnostic Logging Fidelity
- Catching exceptions with generic fallbacks hides bugs
- Propagate literal HTTP response codes immediately
- Print raw error payloads, not aestheticized summaries

### 5. Mock Data Integrity
- Hand-written test files cause artificial data bugs
- Ground test suites in validated reference data (real GRCh38 APOE variants)

### 6. Code Scope and Scaling
- Misplaced indentation breaks variable lifetimes across execution blocks
- Transition early from single queries to batch APIs to prevent network saturation

---

## Known Limitations

1. **GRCh38-only** — VEP queries use the GRCh38/hg38 reference. GRCh37/hg19 not supported.
2. **Human genome only** — Ensembl VEP is restricted to *Homo sapiens*.
3. **Network-dependent** — Requires live internet; no offline annotation mode.
4. **Single-threaded** — Queries are sequential; parallel processing not implemented.
5. **Rate limits** — Batch sizes are conservative to avoid API throttling.

---

## Troubleshooting

### "API_Error_50x" or Connection Timeout

The Ensembl VEP or MyVariant.info API is temporarily unavailable.

**Solution:** Retry manually; built-in backoff handles transient 5xx errors up to 2 retries.

### "No_Transcript_Overlap"

A variant does not overlap any known transcript in Ensembl.

**Expected behavior** — Often occurs for intronic, intergenic, or structural variants.

### Missing GENE field in variant

If the VCF INFO field lacks a GENE tag, the report falls back to the gene extracted from Ensembl VEP.

**Solution:** Ensure VCF INFO field includes `GENE=...` or let VEP provide it.

### Report is empty (no candidates)

All variants are either common (AF ≥ 0.01) or predicted low-impact.

**Solution:** Lower the `--af-threshold` or remove the functional impact filter.

---

## Performance Notes

- **Parsing:** 1,000 variants in <100 ms
- **Annotation:** ~1–2 seconds per variant (API network latency)
- **Batch queries:** 100 variants in ~3–5 seconds (POST more efficient than GET)
- **Report generation:** <500 ms for HTML rendering

---

## Contributing

Contributions are welcome! Areas of interest:
- Parallel API queries (threading/async)
- Support for additional genome builds (GRCh37, others)
- Integration with additional annotation sources (SIFT, PolyPhen, CADD)
- Offline annotation mode with precomputed databases

---

## License

This project is open source. See LICENSE for details.

---

## Citation & References

- **Ensembl VEP:** [https://useast.ensembl.org/info/docs/tools/vep/](https://useast.ensembl.org/info/docs/tools/vep/)
- **MyVariant.info:** [https://myvariant.info/](https://myvariant.info/)
- **VCF Format Spec:** [https://samtools.github.io/hts-specs/VCFv4.2.pdf](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
- **HGVS Nomenclature:** [https://varnomen.hgvs.org/](https://varnomen.hgvs.org/)
- **gnomAD Database:** [https://gnomad.broadinstitute.org/](https://gnomad.broadinstitute.org/)
- **ClinVar:** [https://www.ncbi.nlm.nih.gov/clinvar/](https://www.ncbi.nlm.nih.gov/clinvar/)

---

## Contact

For questions or issues, please file a GitHub issue or contact the maintainers.

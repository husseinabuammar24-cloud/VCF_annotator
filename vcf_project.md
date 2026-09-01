# Variant Effect Predictor & Variant Call Format (VCF) Annotator

This is a core project used in **clinical genetics and cancer genomics**. You will build a Python tool that reads raw genomic variant files, determines how mutations affect proteins, and pulls live biological data to rank variants by pathogenicity.

## Key Requirements & Modules

### 1. VCF Parser (Custom & Library-Based)

* Read `.vcf` files containing genomic coordinates:

  * Chromosome
  * Position
  * Reference allele
  * Alternate allele
* Parse the `INFO` and `FORMAT` fields to extract:

  * Read coverage
  * Allele frequency
  * Quality scores

### 2. Gene Coordinates & Variant Annotation Engine

* Fetch human genome annotations (e.g., **GRCh38/hg38**) using:

  * `Biopython`
  * **Ensembl REST API**
* Determine the impact of a variant:

  * **Synonymous:** Coding mutation with no amino acid change.
  * **Missense:** Single amino acid substitution.
  * **Nonsense:** Premature stop codon created.
  * **Frameshift/Indel:** Insertion or deletion that alters the codon reading frame.

### 3. API Integration for Pathogenicity Scoring

Query external APIs such as:

* **ClinVar** — known disease associations and clinical significance.
* **dbSNP** — known variant identifiers and population information.
* **MyVariant.info** — aggregated variant annotations.
* **gnomAD** — population allele frequencies.

Use the retrieved information to rank variants based on their potential pathogenicity.

### 4. Filtering & HTML Summary Report

Implement a clinical filtering strategy, for example:

> Keep variants with **population allele frequency < 0.01** AND **predicted High/Moderate impact**.

Generate an automated **HTML report** summarizing candidate disease-causing variants.

The report should include:

* Variant coordinates
* Reference and alternate alleles
* Gene and transcript information
* Predicted consequence
* Allele frequency
* Quality and coverage metrics
* Clinical significance
* Pathogenicity ranking
* Interactive links to:

  * PubMed
  * ClinVar

## Suggested Python Stack

| Tool / Library             | Purpose                                                  |
| -------------------------- | -------------------------------------------------------- |
| **`pysam`** / **`cyvcf2`** | High-performance VCF parsing                             |
| **`biopython`**            | Sequence handling, codon translation, and PubMed queries |
| **`pandas`**               | Data manipulation and tabular reporting                  |
| **`requests`**             | Querying Ensembl and MyVariant REST APIs                 |
| **`jinja2`**               | Generating the final HTML report                         |

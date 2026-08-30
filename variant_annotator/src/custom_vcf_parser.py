"""
custom_vcf_parser.py

A from-scratch VCF (Variant Call Format) parser, built without any
bioinformatics libraries -- just string splitting -- so the internal
structure of a VCF file is fully visible rather than hidden inside a
library call.

VCF line anatomy (tab-separated):
    CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE1 [SAMPLE2 ...]

    CHROM   Chromosome (e.g. "chr17")
    POS     1-based genomic position of the variant
    ID      Known variant ID if any (e.g. an rsID like "rs80357382"), else "."
    REF     Reference allele -- what the genome normally has here
    ALT     Alternate allele -- what this sample has instead
    QUAL    Phred-scaled confidence that the variant call is real
    FILTER  "PASS" if the variant passed quality filters, else a fail reason
    INFO    Semicolon-separated site-level annotations (see below)
    FORMAT  Colon-separated list of per-sample field NAMES (see below)
    SAMPLE* One column per sample, holding per-sample field VALUES

INFO field:    semicolon-separated key=value pairs,
               e.g. "DP=45;AF=0.002;GENE=BRCA1"
FORMAT field:  colon-separated field NAMES, e.g. "GT:DP:AD"
SAMPLE field:  colon-separated VALUES matching those names,
               e.g. "0/1:45:32,13" -> GT=0/1, DP=45, AD=32,13

--------------------------------------------------------------------------
HOW THE PIECES FIT TOGETHER (read this before diving into each function)
--------------------------------------------------------------------------
1. parse_vcf()             opens the file and turns every data line into
                            a Variant object. It is the main entry point --
                            in normal use, this is the only function you
                            call directly.
2. parse_info_field()      a helper used INSIDE parse_vcf() to turn one
                            variant's INFO column into a dict.
3. parse_sample_field()    a helper used INSIDE parse_vcf() to turn one
                            sample's genotype column into a dict, using
                            FORMAT as the list of field names.
4. Variant                 the data object parse_vcf() hands back to you,
                            one per variant call, with convenience
                            properties (is_snv, variant_type, etc.) for
                            classifying it.
5. split_multiallelic()    an optional POST-processing step you run on
                            each Variant after parsing, if a site has more
                            than one ALT allele and you want one Variant
                            per allele instead of one Variant per line.

Typical usage:
    variants = parse_vcf("real_grch38_test.vcf")
    for v in variants:
        for single_allele_variant in split_multiallelic(v):
            print(single_allele_variant.variant_type)
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Variant:
    """
    Represents a single variant call: one row of a VCF file, for one
    genomic position, plus the per-sample genotype data attached to it.

    This is a plain data container (a @dataclass, so __init__ is
    generated automatically from the fields below) -- you don't create
    Variant objects yourself in normal use. parse_vcf() builds them for
    you; you just read their fields and properties.

    Fields (set once, when the Variant is created):
        chrom       Chromosome name, e.g. "chr17"
        pos         1-based position on that chromosome
        id          Known variant ID (rsID), or "." if unknown
        ref         Reference allele, e.g. "G"
        alt         Alternate allele, e.g. "A" (may be comma-separated
                    if multi-allelic, until split_multiallelic() is run)
        qual        Confidence score as a float, or None if the file
                    had "." (missing) in the QUAL column
        filter      "PASS", or a reason the variant failed filtering
        info        dict of INFO annotations, e.g. {"DP": "45", "GENE": "BRCA1"}
        sample_data dict mapping sample name -> that sample's fields, e.g.
                    {"SAMPLE1": {"GT": "0/1", "DP": "45", "AD": "32,13"}}

    Properties (computed on the fly from ref/alt -- see below):
        is_snv, is_insertion, is_deletion, variant_type
    """
    
    chrom: str          # e.g. "chr17"
    pos: int             # 1-based position on that chromosome
    id: str              # known variant ID (rsID) or "."
    ref: str              # reference allele, e.g. "G"
    alt: str              # alternate (mutant) allele, e.g. "A"
    qual: Optional[float]  # confidence score, or None if "." in the file
    filter: str            # "PASS" or a fail reason like "LowQual"
    info: dict = field(default_factory=dict)
    # sample_data maps sample name -> that sample's GT/DP/AD fields, e.g.:
    #   {"SAMPLE1": {"GT": "0/1", "DP": "45", "AD": "32,13"}, "SAMPLE2": {...}}
    sample_data: dict = field(default_factory=dict)




    @property
    def is_snv(self) -> bool:
        """
        SNV = Single Nucleotide Variant: exactly one base swapped for
        another, e.g. G -> A.

        Use: call this (or check variant_type) when you want to filter
        a variant list down to simple single-base substitutions and
        skip indels.

        Returns True only when both REF and ALT are a single character.
        """
        return len(self.ref) == 1 and len(self.alt) == 1

    @property
    def is_insertion(self) -> bool:
        """
        Use: check this when you specifically want to find variants
        where extra DNA was added at this position.

        Returns True when ALT is longer than REF, e.g. REF="A", ALT="ATG".
        """
        return len(self.alt) > len(self.ref)

    @property
    def is_deletion(self) -> bool:
        """
        Use: check this when you specifically want to find variants
        where DNA was removed at this position.

        Returns True when REF is longer than ALT, e.g. REF="ATG", ALT="A".
        """
        return len(self.ref) > len(self.alt)

    @property
    def variant_type(self) -> str:
        """
        Human-readable classification of this variant: one of
        "SNV", "insertion", "deletion", or "complex".

        Use: this is usually the property you actually want in reports
        or print statements, rather than checking is_snv / is_insertion
        / is_deletion individually -- it's the single-call summary of
        the other three properties.

        It's computed from is_snv / is_insertion / is_deletion each
        time it's accessed (rather than stored as a fixed field) so it
        can never disagree with the current ref/alt values.

        "complex" covers anything that isn't a clean single-base swap,
        pure insertion, or pure deletion -- e.g. a multi-base
        substitution like REF="AG", ALT="CT" (same length, >1 base).
        """
        if self.is_snv:
            return "SNV"
        elif self.is_insertion:
            return "insertion"
        elif self.is_deletion:
            return "deletion"
        return "complex"  # e.g. multi-base substitution like "AG" -> "CT"


def parse_info_field(info_str: str) -> dict:
    """
    Parse a VCF INFO string into a dictionary.

    Use: called internally by parse_vcf() once per data line, to fill
    in Variant.info. You generally won't call this yourself unless
    you're parsing an INFO string you got from somewhere other than
    parse_vcf() (e.g. one you built or edited by hand).

    INFO entries are separated by ';'. Most entries are key=value pairs
    (e.g. "DP=45"), but VCF also allows bare flags with no value at all
    (e.g. "SOMATIC" on its own) -- those are stored as True, meaning
    "this flag is present."

    Args:
        info_str: the raw INFO column text from one VCF line, e.g.
                  "DP=45;AF=0.002;AC=1;GENE=BRCA1"

    Returns:
        A dict of {field_name: value}. Values are strings for key=value
        pairs, or the literal True for bare flags.

    Example:
        parse_info_field("DP=45;AF=0.002;AC=1;GENE=BRCA1")
        -> {"DP": "45", "AF": "0.002", "AC": "1", "GENE": "BRCA1"}

        parse_info_field("DP=45;SOMATIC")
        -> {"DP": "45", "SOMATIC": True}
    """
    info = {}
    for entry in info_str.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)  # split only on the FIRST '='
            info[key] = value
        else:
            info[entry] = True  # bare flag, e.g. "SOMATIC"
    return info


def parse_sample_field(format_str: str, sample_str: str) -> dict:
    """
    Parse one sample's genotype data using the FORMAT column as a key.

    Use: called internally by parse_vcf(), once per sample per data
    line, to fill in Variant.sample_data[sample_name]. You generally
    won't call this yourself unless you're parsing a FORMAT/sample pair
    you obtained separately from parse_vcf().

    FORMAT lists the *names* of the per-sample fields, in a fixed order,
    e.g. "GT:DP:AD". Each SAMPLE column then holds the *values* for that
    sample, in that same order, e.g. "0/1:45:32,13". zip() pairs them up
    position by position: 1st name with 1st value, 2nd with 2nd, etc.

    Args:
        format_str: the FORMAT column text, e.g. "GT:DP:AD" -- the
                    field names, in order.
        sample_str: one sample's column text, e.g. "0/1:45:32,13" --
                    the values, in the same order as format_str.

    Returns:
        A dict of {field_name: value}, e.g.
        {"GT": "0/1", "DP": "45", "AD": "32,13"}
        Note: if format_str and sample_str have different numbers of
        entries, zip() silently truncates to the shorter one.

    Example:
        parse_sample_field("GT:DP:AD", "0/1:45:32,13")
        -> {"GT": "0/1", "DP": "45", "AD": "32,13"}
    """
    keys = format_str.split(":")
    values = sample_str.split(":")
    return dict(zip(keys, values))


def parse_vcf(filepath: str) -> list[Variant]:
    """
    Read a VCF file from disk and return a list of Variant objects,
    one per data row, each carrying its own per-sample genotype data.

    Use: this is the main entry point of the module -- the function you
    call directly. It reads the whole file, skips the metadata lines,
    reads the sample names from the header, and builds one Variant per
    variant-call line, using parse_info_field() and parse_sample_field()
    as helpers along the way.

    VCF files have three kinds of lines:
      1. "##..."    metadata (file format version, field definitions) -- skipped
      2. "#CHROM..." the single column header line -- used to learn sample names
      3. everything else -- one actual variant call per line

    Args:
        filepath: path to a VCF file on disk.

    Returns:
        A list of Variant objects, in the same order as the file's
        data lines. If a variant is multi-allelic (ALT has commas,
        e.g. "G,T"), it is still returned as ONE Variant at this stage
        -- call split_multiallelic() afterwards if you want it broken
        into one Variant per allele.

    Example:
        variants = parse_vcf("real_grch38_test.vcf")
        variants[0].chrom      -> "chr17"
        variants[0].info       -> {"DP": "45", "GENE": "BRCA1"}
        variants[0].sample_data["SAMPLE1"]["GT"]  -> "0/1"
    """
    variants = []
    # Filled in once we hit the "#CHROM" header line. Holds sample column
    # names in file order, e.g. ["SAMPLE1", "SAMPLE2"], so we know which
    # name to attach to each value column later on.
    sample_names: list[str] = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()  # remove the trailing newline

            # The #CHROM line is the column header. Columns 0-8 are
            # always the same 9 fixed fields (CHROM..FORMAT); everything
            # from column 9 onward is a sample name, and there can be
            # any number of them.
            if line.startswith("#CHROM"):
                header_fields = line.split("\t")
                sample_names = header_fields[9:]
                continue

            # "##" lines are metadata/header declarations -- not needed
            # for basic parsing, so skip them.
            if line.startswith("#"):
                continue
            if not line:
                continue  # skip any blank lines

            # Every data row is tab-separated. The first 8 fields are
            # always in this fixed order, per the VCF spec.
            fields = line.split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info_str = fields[:8]

            variant = Variant(
                chrom=chrom,
                pos=int(pos),
                id=vid,
                ref=ref,
                alt=alt,
                # QUAL is sometimes "." (missing) rather than a number
                qual=float(qual) if qual != "." else None,
                filter=filt,
                info=parse_info_field(info_str),
            )

            # FORMAT + sample columns (fields 9 onward) are optional --
            # some VCFs are site-only, with no per-sample genotype data.
            if len(fields) > 9:
                format_str = fields[8]        # field NAMES, e.g. "GT:DP:AD"
                sample_columns = fields[9:]   # one column of VALUES per sample
                # Pair each sample name with its matching value column,
                # in order, and parse each one individually.
                for name, sample_str in zip(sample_names, sample_columns):
                    variant.sample_data[name] = parse_sample_field(format_str, sample_str)

            variants.append(variant)
    return variants


def split_multiallelic(variant: Variant) -> list[Variant]:
    """
    Split one multi-allelic Variant (alt like "G,T") into multiple
    single-allele Variant objects (one for A>G, one for A>T).

    Use: call this on each Variant AFTER parse_vcf(), typically in a
    loop, when you want every Variant you work with to represent
    exactly one allele. This matters because a single VCF line can
    report more than one alternate allele at the same position (e.g. a
    site that is G in the reference but appears as either A or T across
    different samples) -- most downstream analysis expects one
    ref->alt pair per Variant, not a comma-joined list.

    It is safe to call on EVERY variant, not just the multi-allelic
    ones: if ALT has no comma, the variant is returned unchanged,
    wrapped in a list of length 1.

    Args:
        variant: a Variant, as produced by parse_vcf(). Its own
                 sample_data and info are read but not modified.

    Returns:
        A list of Variant objects, one per allele in the original ALT:
          - length 1, unchanged, if the variant was already single-allelic
          - length N, one new Variant per allele, if ALT had N comma-
            separated alleles

        Each new Variant copies chrom/pos/id/ref/qual/filter from the
        original, gets its own single ALT allele, and gets its own
        info dict (a shallow copy, so editing one variant's info won't
        affect another's).

        Per-sample AD (allele depth) is split too: AD normally holds
        one depth for REF plus one depth per ALT allele, in order
        (e.g. "32,10,5" for REF plus two ALTs). For ALT index i, the
        new Variant keeps only the REF depth and that allele's own
        depth, e.g. "32,10" for the first ALT, "32,5" for the second.

    Example:
        # variant.alt == "G,T", variant.sample_data["S1"]["AD"] == "32,10,5"
        split_multiallelic(variant)
        # -> [Variant(alt="G", sample_data["S1"]["AD"]="32,10"),
        #     Variant(alt="T", sample_data["S1"]["AD"]="32,5")]

        # variant.alt == "G" (already single-allelic)
        split_multiallelic(variant)
        # -> [variant]   (the same object, unchanged)
    """
    alt_alleles = variant.alt.split(",")

    # Already single-allelic -- nothing to split.
    if len(alt_alleles) == 1:
        return [variant]

    split_variants = []
    for i, single_alt in enumerate(alt_alleles):
        # Build a new Variant, identical to the original except for
        # a single ALT allele instead of the comma-joined list.
        new_variant = Variant(
            chrom=variant.chrom,
            pos=variant.pos,
            id=variant.id,
            ref=variant.ref,
            alt=single_alt,
            qual=variant.qual,
            filter=variant.filter,
            info=dict(variant.info),  # shallow copy so edits don't leak back
        )

        # Per-sample fields need slicing too: AD has one value for REF
        # plus one value per ALT allele, in order. For ALT index i, we
        # want AD[0] (ref depth) and AD[i+1] (this allele's depth).
        for sample_name, fields in variant.sample_data.items():
            new_fields = dict(fields)  # copy so we don't mutate the original
            if "AD" in fields:
                ad_values = fields["AD"].split(",")
                ref_depth = ad_values[0]
                this_allele_depth = ad_values[i + 1]
                new_fields["AD"] = f"{ref_depth},{this_allele_depth}"
            new_variant.sample_data[sample_name] = new_fields

        split_variants.append(new_variant)

    return split_variants

if __name__ == "__main__":
    import os
    script_dir = os.path.dirname(os.path.abspath(__file__))
    vcf_path = os.path.join(script_dir, "..", "data", "real_grch38_test.vcf")
    raw_variants = parse_vcf(vcf_path)

    # Expand every multi-allelic row into separate single-allele variants
    variants = []
    for v in raw_variants:
        variants.extend(split_multiallelic(v))

    for v in variants:
        print(f"{v.chrom}:{v.pos} {v.ref}>{v.alt} [{v.variant_type}] gene={v.info.get('GENE')}")
        for sample_name, fields in v.sample_data.items():
            print(f"    {sample_name}: GT={fields.get('GT')} DP={fields.get('DP')} AD={fields.get('AD')}")
import requests
import time
from custom_vcf_parser import parse_vcf, split_multiallelic


################################
# Step 1: The Automated Blueprint (Using Ensembl VEP)
################################
SERVER = "https://rest.ensembl.org"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

def determine_variant_impact(chrom: str, pos: int, ref: str, alt: str, retries: int = 2) -> str:
    """
    Queries the Ensembl VEP REST API to retrieve the exact biological impact
    of a genomic variant on human transcripts (GRCh38).
    """
    # Clean chromosome names for Ensembl formatting rules
    chrom_clean = chrom.replace("chr", "")

    # VCF encodes indels with a shared "anchor" base at the start of REF/ALT
    # (e.g. REF=ATCT, ALT=A means "delete TCT, keep the leading A").
    # VEP's region-based GET endpoint wants just the deleted/inserted bases
    # and a position that excludes that anchor -- passing the raw VCF-style
    # REF/ALT straight through gives coordinates that don't reliably line
    # up with the actual variant, which is why deletions were flip-flopping
    # between categories across identical reruns.
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt[0]:
        # Deletion: drop the anchor base from REF, query starts one base later.
        deleted = ref[1:]
        start = pos + 1
        end = pos + len(deleted)
        vep_allele = "-"
    elif len(alt) > 1 and len(ref) == 1 and alt[0] == ref[0]:
        # Insertion: VEP wants start > end for pure insertions, allele = inserted bases.
        inserted = alt[1:]
        start = pos + 1
        end = pos
        vep_allele = inserted
    else:
        # SNV or complex substitution -- no anchor base to strip.
        start = pos
        end = pos + len(ref) - 1
        vep_allele = alt

    # Build VEP REST API notation path string
    ext = f"/vep/homo_sapiens/region/{chrom_clean}:{start}-{end}/{vep_allele}?"

    try:
        response = None
        for attempt in range(retries + 1):  # try up to retries+1 times total
            response = requests.get(SERVER + ext, headers=HEADERS)
            if response.ok:
                break
            if response.status_code >= 500 and attempt < retries:
                time.sleep(1 + attempt)   # wait 1s, then 2s, before retrying
                continue
            break   # not a retryable 5xx, or out of retries -- stop trying

        if not response.ok:
            # TEMPORARY DEBUG: print the response body so we can see exactly
            # why VEP is rejecting the request on a 400 (client-side error).
            # Remove this print once the 400s are understood/fixed.
            if response.status_code == 400:
                print(f"[debug] 400 body for {chrom}:{pos} {ref}>{alt} "
                      f"(queried as {chrom_clean}:{start}-{end}/{vep_allele}): "
                      f"{response.text[:300]}")
            return f"API_Error_{response.status_code}"

        vep_data = response.json()
        if not vep_data:
            return "No_Transcript_Overlap"

        # Inspect consequences predicted for the variant entry
        consequences = vep_data[0].get("transcript_consequences", [])
        if not consequences:
            return "Intergenic/Intronic"

        # Pick the highest-impact consequence terms across transcripts
        # Ensembl Sequence Ontology (SO) terms mapping to your target categories:
        impact_terms = set()
        for csq in consequences:
            for term in csq.get("consequence_terms", []):
                impact_terms.add(term)

        # Classify the primary functional impact profile
        if "frameshift_variant" in impact_terms:
            return "Frameshift/Indel"
        elif "stop_gained" in impact_terms:
            return "Nonsense"
        elif "missense_variant" in impact_terms:
            return "Missense"
        elif "synonymous_variant" in impact_terms:
            return "Synonymous"

        # Return fallback if it's a different coding/splicing consequence type
        return list(impact_terms)[0] if impact_terms else "Unknown"

    except Exception as e:
        print(f"[warning] {chrom}:{pos} failed: {e}")
        return "Connection_Error"


###########################
# Step 2: Integration Loop
###########################

# Assuming you already loaded your file using: raw_variants = parse_vcf("input.vcf")
raw_variants = parse_vcf("variant_annotator/data/real_grch38_test.vcf")

print(f"{'VARIANT':<20} | {'TYPE':<12} | {'FUNCTIONAL IMPACT'}")
print("-" * 55)

for raw_v in raw_variants:
    for v in split_multiallelic(raw_v):

        # Call the automated evaluation function
        functional_impact = determine_variant_impact(v.chrom, v.pos, v.ref, v.alt)

        variant_str = f"{v.chrom}:{v.pos} {v.ref}>{v.alt}"
        print(f"{variant_str:<20} | {v.variant_type:<12} | {functional_impact}")

        # Protect against Ensembl endpoint rate throttling (Limit < 15 req/sec)
        time.sleep(0.1)
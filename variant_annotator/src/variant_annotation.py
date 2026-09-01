import time
from typing import Any, Dict, List
import requests
from custom_vcf_parser import parse_vcf, split_multiallelic

SERVER = "https://rest.ensembl.org"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}


def build_vep_region_string(
    chrom: str, pos: int, ref: str, alt: str
) -> tuple[str, str]:
    """Translates VCF variants into Ensembl VEP REST region format.

    Handles anchor base trimming for insertions and deletions.
    """
    chrom_clean = chrom.replace("chr", "")

    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt[0]:
        # Deletion
        deleted = ref[1:]
        start = pos + 1
        end = pos + len(deleted)
        vep_allele = "-"
    elif len(alt) > 1 and len(ref) == 1 and alt[0] == ref[0]:
        # Insertion
        inserted = alt[1:]
        start = pos + 1
        end = pos
        vep_allele = inserted
    else:
        # SNV / Substitution
        start = pos
        end = pos + len(ref) - 1
        vep_allele = alt

    region_str = f"{chrom_clean}:{start}-{end}/{vep_allele}"
    lookup_key = f"{chrom}:{pos}_{ref}/{alt}"

    return region_str, lookup_key


def parse_vep_consequence(vep_entry: Dict[str, Any]) -> str:
    """Parses raw VEP JSON output into a single high-level impact category."""
    consequences = vep_entry.get("transcript_consequences", [])
    if not consequences:
        return "Intergenic/Intronic"

    impact_terms = set()
    for csq in consequences:
        for term in csq.get("consequence_terms", []):
            impact_terms.add(term)

    # Priority ranking of consequence terms
    if "frameshift_variant" in impact_terms:
        return "Frameshift/Indel"
    elif "stop_gained" in impact_terms:
        return "Nonsense"
    elif "missense_variant" in impact_terms:
        return "Missense"
    elif "synonymous_variant" in impact_terms:
        return "Synonymous"

    return list(impact_terms)[0] if impact_terms else "Unknown"


def determine_variant_impact(
    chrom: str, pos: int, ref: str, alt: str, retries: int = 2
) -> str:
    """Queries single variant against Ensembl VEP REST GET endpoint."""
    region_str, _ = build_vep_region_string(chrom, pos, ref, alt)
    ext = f"/vep/homo_sapiens/region/{region_str}?"

    for attempt in range(retries + 1):
        try:
            response = requests.get(
                SERVER + ext, headers=HEADERS, timeout=10
            )
            if response.ok:
                vep_data = response.json()
                return (
                    parse_vep_consequence(vep_data[0])
                    if vep_data
                    else "No_Transcript_Overlap"
                )
            if response.status_code >= 500 and attempt < retries:
                time.sleep(1 + attempt)
                continue
            return f"API_Error_{response.status_code}"
        except requests.exceptions.RequestException as e:
            if attempt == retries:
                return f"Connection_Error: {e}"
            time.sleep(1)

    return "Unknown_Failure"


def determine_batch_variant_impacts(variants: List[Any]) -> Dict[str, str]:
    """Queries a batch of variants in a SINGLE POST request to Ensembl VEP.

    Significantly reduces network overhead for larger VCFs.
    """
    post_hgvs_list = []
    key_mapping = {}

    for v in variants:
        region_str, lookup_key = build_vep_region_string(
            v.chrom, v.pos, v.ref, v.alt
        )
        post_hgvs_list.append(region_str)
        key_mapping[region_str] = lookup_key

    ext = "/vep/homo_sapiens/region"
    payload = {"variants": post_hgvs_list}

    results = {}
    try:
        response = requests.post(
            SERVER + ext, headers=HEADERS, json=payload, timeout=15
        )
        if response.ok:
            for item in response.json():
                input_str = item.get("input")
                lookup_key = key_mapping.get(input_str)
                if lookup_key:
                    results[lookup_key] = parse_vep_consequence(item)
        else:
            print(f"[warning] Batch VEP call failed: HTTP {response.status_code}")
    except requests.exceptions.RequestException as e:
        print(f"[warning] Batch VEP connection error: {e}")

    return results


def main():
    vcf_path = "variant_annotator/data/real_grch38_test.vcf"
    raw_variants = parse_vcf(vcf_path)

    # Flatten multiallelic sites into distinct variant instances
    split_variants = [
        v for raw_v in raw_variants for v in split_multiallelic(raw_v)
    ]

    print(f"{'VARIANT':<20} | {'TYPE':<12} | {'FUNCTIONAL IMPACT'}")
    print("-" * 55)

    # Option A: Fast Batch POST API query
    batch_impacts = determine_batch_variant_impacts(split_variants)

    for v in split_variants:
        lookup_key = f"{v.chrom}:{v.pos}_{v.ref}/{v.alt}"
        impact = batch_impacts.get(lookup_key)

        # Fallback to single GET if batch returned no data for this variant
        if not impact:
            impact = determine_variant_impact(v.chrom, v.pos, v.ref, v.alt)

        variant_str = f"{v.chrom}:{v.pos} {v.ref}>{v.alt}"
        print(f"{variant_str:<20} | {v.variant_type:<12} | {impact}")


if __name__ == "__main__":
    main()
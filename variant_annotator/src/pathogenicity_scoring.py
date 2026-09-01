import json
from typing import Any, Dict, List, Optional
import requests


def convert_variant_to_hgvs(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Converts VCF coordinates into an HGVS genomic variant identifier

    compatible with the MyVariant.info REST API.
    """
    clean_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return f"{clean_chrom}:g.{pos}{ref}>{alt}"


def extract_dbsnp_rsid(data: Dict[str, Any]) -> str:
    """Extracts dbSNP rsID identifier if present in MyVariant payload."""
    dbsnp_data = data.get("dbsnp", {})
    if isinstance(dbsnp_data, dict):
        rsid = dbsnp_data.get("rsid")
        if rsid:
            return rsid if str(rsid).startswith("rs") else f"rs{rsid}"
    elif isinstance(dbsnp_data, list) and dbsnp_data:
        rsid = dbsnp_data[0].get("rsid")
        if rsid:
            return rsid if str(rsid).startswith("rs") else f"rs{rsid}"
    return "N/A"


def extract_gnomad_af(data: Dict[str, Any]) -> float:
    """Safely extracts maximum population allele frequency across gnomAD exomes/genomes.

    Returns 0.0 if variant is absent from gnomAD (novel/rare).
    """
    gnomad_exome = data.get("gnomad_exome", {})
    gnomad_genome = data.get("gnomad_genome", {})

    af_exome = gnomad_exome.get("af", {}).get("af") if isinstance(gnomad_exome, dict) else None
    af_genome = gnomad_genome.get("af", {}).get("af") if isinstance(gnomad_genome, dict) else None

    frequencies = []
    for val in (af_exome, af_genome):
        if val is not None:
            try:
                frequencies.append(float(val))
            except (ValueError, TypeError):
                continue

    return max(frequencies) if frequencies else 0.0


def extract_clinvar_info(data: Dict[str, Any]) -> tuple[set[str], list[str]]:
    """Extracts unique clinical significance terms and associated condition names

    from ClinVar records inside MyVariant payload.
    """
    clinvar_data = data.get("clinvar")
    significances = set()
    conditions = set()

    if not clinvar_data or not isinstance(clinvar_data, dict):
        return significances, []

    rcv_records = clinvar_data.get("rcv", [])
    if isinstance(rcv_records, dict):
        rcv_records = [rcv_records]

    for record in rcv_records:
        sig = record.get("clinical_significance")
        if sig:
            significances.add(sig)
        
        # Extract disease conditions
        conditions_data = record.get("conditions", {})
        if isinstance(conditions_data, dict):
            name = conditions_data.get("name")
            if name and name.lower() != "not specified":
                conditions.add(name)
        elif isinstance(conditions_data, list):
            for c in conditions_data:
                name = c.get("name")
                if name and name.lower() != "not specified":
                    conditions.add(name)

    return significances, list(conditions)


def summarize_pathogenicity(clinvar_significances: set[str], impact: str) -> str:
    """Assigns a unified clinical classification by prioritizing ClinVar entries

    and falling back to VEP functional impact predictions.
    """
    terms = {s.lower() for s in clinvar_significances}

    # Priority 1: Direct ClinVar Evidence
    if any("pathogenic" in t and "likely" not in t for t in terms):
        return "Pathogenic"
    elif any("likely pathogenic" in t for t in terms):
        return "Likely Pathogenic"
    elif any("uncertain significance" in t for t in terms):
        return "Uncertain Significance (VUS)"
    elif any("likely benign" in t for t in terms):
        return "Likely Benign"
    elif any("benign" in t for t in terms):
        return "Benign"
    elif any("risk factor" in t for t in terms):
        return "Risk Factor"

    # Priority 2: VEP Consequence Fallback when ClinVar is unannotated
    if impact in ["Nonsense", "Frameshift/Indel"]:
        return "Predicted High Impact (VUS)"
    elif impact == "Missense":
        return "Predicted Moderate Impact (VUS)"
    
    return "Uncertain / Unannotated"


def annotate_variant_pathogenicity(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    functional_impact: str = "Unknown",
    timeout: int = 10
) -> Dict[str, Any]:
    """Master API coordinator: queries MyVariant.info and consolidates dbSNP,

    gnomAD, ClinVar, and Pathogenicity scoring into a unified dict.
    """
    hgvs_id = convert_variant_to_hgvs(chrom, pos, ref, alt)
    url = f"https://myvariant.info/v1/variant/{hgvs_id}"
    params = {"assembly": "hg38"}

    record = {
        "hgvs": hgvs_id,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "rsid": "N/A",
        "functional_impact": functional_impact,
        "clinvar_sig": "None",
        "pathogenicity_summary": "Uncertain / Unannotated",
        "gnomad_af": 0.0,
        "conditions": [],
        "clinvar_id": None
    }

    try:
        response = requests.get(url, params=params, timeout=timeout)
        if response.status_code == 200:
            data = response.json()

            # Extract component annotations
            record["rsid"] = extract_dbsnp_rsid(data)
            record["gnomad_af"] = extract_gnomad_af(data)
            
            significances, conditions = extract_clinvar_info(data)
            if significances:
                record["clinvar_sig"] = ", ".join(significances)
            record["conditions"] = conditions
            
            # Extract ClinVar Variation ID for linking
            clinvar_data = data.get("clinvar", {})
            if isinstance(clinvar_data, dict):
                record["clinvar_id"] = clinvar_data.get("variant_id") or clinvar_data.get("rcv", [{}])[0].get("accession") if isinstance(clinvar_data.get("rcv"), list) else None

            # Calculate score
            record["pathogenicity_summary"] = summarize_pathogenicity(significances, functional_impact)

    except requests.exceptions.RequestException as e:
        print(f"[warning] API request failed for {hgvs_id}: {e}")

    return record


def extract_all_sources(data: dict) -> dict:
    """Extracts ClinVar, dbSNP, MyVariant, and gnomAD parameters

    from a MyVariant.info JSON response payload.
    """
    # 1. MyVariant.info Top-Level HGVS ID
    myvariant_id = data.get("_id", "N/A")

    # 2. dbSNP rsID
    dbsnp_id = data.get("dbsnp", {}).get("rsid") or data.get("rsid", "N/A")

    # 3. ClinVar Annotations & Disease Conditions
    clinvar_info = data.get("clinvar", {})
    clinvar_id = clinvar_info.get("variant_id", "N/A")

    rcv = clinvar_info.get("rcv", [])
    if isinstance(rcv, dict):
        rcv = [rcv]

    conditions = list(
        {item.get("conditions", {}).get("name") for item in rcv if item.get("conditions")}
    )

    # 4. gnomAD Population Allele Frequency
    gnomad_af = (
        data.get("gnomad_genome", {}).get("af", {}).get("af")
        or data.get("gnomad_exome", {}).get("af", {}).get("af")
        or 0.0
    )

    return {
        "myvariant_id": myvariant_id,
        "dbsnp_id": dbsnp_id,
        "clinvar_id": clinvar_id,
        "conditions": conditions,
        "gnomad_af": float(gnomad_af),
    }
    
if __name__ == "__main__":
    # Test dataset containing both variants
    test_variants = [
        {
            "chrom": "chr19",
            "pos": 44908684,
            "ref": "T",
            "alt": "C",
            "functional_impact": "Missense",
        },
        {
            "chrom": "chr19",
            "pos": 44908822,
            "ref": "C",
            "alt": "T",
            "functional_impact": "Missense",
        },
    ]

    # Process and print results for all variants
    annotated_results = [
        annotate_variant_pathogenicity(**v) for v in test_variants
    ]
    print(json.dumps(annotated_results, indent=2))
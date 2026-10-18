#!/usr/bin/env python3
"""Collect and normalize the first public data batch for the ADC ABC-substrate project."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import openpyxl


ROOT = Path(__file__).resolve().parents[1]
RAW_ROTH = ROOT / "raw" / "roth_2026"
RAW_DAOOD = ROOT / "raw" / "daood_2025"
PROCESSED = ROOT / "processed"
PROCESSED.mkdir(parents=True, exist_ok=True)

ROTH_XML = RAW_ROTH / "PMC12883344_fulltext.xml"
DAOOD_XLSX = RAW_DAOOD / "mp5c01065_si_002.xlsx"

ROTH_DOI = "10.20517/cdr.2025.151"
ROTH_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC12883344/"
DAOOD_DOI = "10.1021/acs.molpharmaceut.5c01065"
DAOOD_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC12587445/"


PUBCHEM_QUERY_ALIASES = {
    "α-Amanitin": ["alpha-Amanitin"],
    "Ansamitocin P3": ["Ansamitocin P3"],
    "Auristatin E": ["Auristatin E"],
    "Auristatin F": ["Auristatin F"],
    "Calicheamicin γ1": ["Calicheamicin gamma1", "Calicheamicin gamma 1"],
    "Camptothecin": ["Camptothecin"],
    "Dasatinib": ["Dasatinib"],
    "DM1": ["Mertansine", "DM1"],
    "DM4": ["Ravtansine", "DM4"],
    "Dolastatin 10": ["Dolastatin 10"],
    "Doxorubicin": ["Doxorubicin"],
    "Duocarmycin DM": ["Duocarmycin DM"],
    "Duocarmycin MA": ["Duocarmycin MA"],
    "Duocarmycin SA": ["Duocarmycin SA"],
    "Duocarmycin TA": ["Duocarmycin TA"],
    "Maytansinol": ["Maytansinol"],
    "MMAE": ["Monomethyl auristatin E", "MMAE"],
    "MMAF": ["Monomethyl auristatin F", "MMAF"],
    "Monomethyl dolastatin 10": ["Monomethyl dolastatin 10", "Demethyldolastatin 10"],
    "N-acetyl-calicheamicin γ1": ["N-acetyl-calicheamicin gamma1"],
    "Nemorubicin": ["Nemorubicin"],
    "PBD Dimer (SGD-1882)": ["SGD-1882"],
    "PNU-159682": ["PNU-159682"],
    "SN-38": ["SN-38"],
    "Tubulysin IM-2": ["Tubulysin IM-2"],
    "Tubulysin M": ["Tubulysin M"],
    "Vinblastine sulfate": ["Vinblastine sulfate", "Vinblastine"],
    "Exatecan": ["Exatecan"],
    "Deruxtecan": ["Deruxtecan", "DXd"],
    "SJG136": ["SJG136", "SJG-136"],
    "SGD-1882": ["SGD-1882"],
    "SG2057": ["SG2057"],
    "SG3199": ["SG3199"],
}


STRUCTURE_REVIEW_NOTES = {
    "Deruxtecan": (
        "manual_review",
        "PubChem name resolution may return a linker-payload/maleimide record rather than the released DXd species. "
        "Do not use this structure for modeling until the exact material tested in Roth 2026 is verified.",
    ),
    "Vinblastine sulfate": (
        "salt_form_review",
        "The paper specifies vinblastine sulfate. Preserve the experimental salt identity and standardize a neutral parent only in a derived table.",
    ),
}

PAPER_NAME_CANONICAL = {
    "PBD Dimer (SGD-1882)": "SGD-1882",
}


def text_content(node) -> str:
    return " ".join("".join(node.itertext()).split())


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_mean_sd(value: str) -> tuple[float | None, float | None]:
    match = re.match(r"^\s*([0-9.]+)\s*±\s*([0-9.]+)\s*$", value)
    if match:
        return float(match.group(1)), float(match.group(2))
    try:
        return float(value), None
    except ValueError:
        return None, None


def parse_roth_tables() -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    root = ET.parse(ROTH_XML).getroot()
    tables = {}
    for table_wrap in root.findall(".//table-wrap"):
        label = table_wrap.findtext("label") or ""
        rows = []
        for tr in table_wrap.findall(".//tr"):
            rows.append([text_content(cell) for cell in list(tr)])
        tables[label] = rows

    table1 = [
        {
            "adc_name": row[0],
            "brand_name": row[1],
            "payload_reported": row[2],
            "linker_class": row[3],
            "fda_status_as_reported": row[4],
            "source_doi": ROTH_DOI,
            "source_table": "Table 1",
            "source_url": ROTH_URL,
        }
        for row in tables["Table 1"][1:]
    ]

    table2 = [
        {
            "panel_id": f"ROTH_T2_{index:02d}",
            "paper_compound_name": row[0],
            "vendor": row[1],
            "mechanism_of_action": row[2],
            "screen_panel_n": 27,
            "screen_result_available_in_article_table": "no",
            "source_doi": ROTH_DOI,
            "source_table": "Table 2",
            "source_url": ROTH_URL,
        }
        for index, row in enumerate(tables["Table 2"][1:], start=1)
    ]

    wide = []
    long = []
    transporter_defs = [
        ("ABCB1", "P-gp", "MDR19", 2, 3),
        ("ABCG2", "BCRP", "R5", 4, 5),
        ("ABCC1", "MRP1", "MRP1", 6, 7),
    ]
    for index, row in enumerate(tables["Table 3"][2:], start=1):
        match = re.match(r"^(.*?)\s*\((nM|μM|pM)\)$", row[0])
        compound_name = match.group(1) if match else row[0]
        unit = match.group(2) if match else ""
        control_mean, control_sd = parse_mean_sd(row[1])
        record = {
            "confirmatory_id": f"ROTH_T3_{index:02d}",
            "paper_compound_name": compound_name,
            "concentration_unit": unit,
            "pcDNA_ic50_mean": control_mean,
            "pcDNA_ic50_sd": control_sd,
            "ABCB1_MDR19_ic50_raw": row[2],
            "ABCB1_RR": float(row[3]),
            "ABCG2_R5_ic50_raw": row[4],
            "ABCG2_RR": float(row[5]),
            "ABCC1_MRP1_ic50_raw": row[6],
            "ABCC1_RR": float(row[7]),
            "source_doi": ROTH_DOI,
            "source_table": "Table 3",
            "source_url": ROTH_URL,
        }
        wide.append(record)

        for gene, protein, cell_line, value_idx, rr_idx in transporter_defs:
            test_mean, test_sd = parse_mean_sd(row[value_idx])
            long.append(
                {
                    "assay_record_id": f"ROTH_T3_{index:02d}_{gene}",
                    "paper_compound_name": compound_name,
                    "transporter_gene": gene,
                    "transporter_protein": protein,
                    "species": "Homo sapiens",
                    "assay_system": "HEK293 stable transfectant cytotoxicity",
                    "control_cell_line": "HEK293-pcDNA",
                    "transporter_cell_line": f"HEK293-{cell_line}",
                    "readout": "IC50",
                    "concentration_unit": unit,
                    "control_ic50_mean": control_mean,
                    "control_ic50_sd": control_sd,
                    "transporter_ic50_mean": test_mean,
                    "transporter_ic50_sd": test_sd,
                    "resistance_ratio": float(row[rr_idx]),
                    "binary_substrate_label": "",
                    "binary_label_status": "not_assigned_continuous_rr_preserved",
                    "specific_inhibitor_rescue_in_same_assay": "no",
                    "evidence_tier": "B",
                    "source_doi": ROTH_DOI,
                    "source_table": "Table 3",
                    "source_url": ROTH_URL,
                    "curation_note": "Confirmatory transfected-cell cytotoxicity result; functional RR is not an intrinsic transport rate.",
                }
            )
    return table1, table2, wide, long


def fetch_pubchem_property(query_name: str) -> tuple[dict | None, str, str]:
    properties = (
        "Title,CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,MolecularWeight,"
        "XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,Charge"
    )
    encoded = urllib.parse.quote(query_name, safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/{properties}/JSON"
    request = urllib.request.Request(url, headers={"User-Agent": "ADC-substrate-data-curation/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            data = json.load(response)
        items = data.get("PropertyTable", {}).get("Properties", [])
        return (items[0] if items else None), url, "200"
    except urllib.error.HTTPError as exc:
        return None, url, str(exc.code)
    except Exception as exc:  # noqa: BLE001 - error text belongs in the query log
        return None, url, f"error:{type(exc).__name__}:{exc}"


def build_compound_registry(table2: list[dict], table3_wide: list[dict]) -> tuple[list[dict], list[dict]]:
    panel_lookup = {row["paper_compound_name"]: row for row in table2}
    reported_names = list(panel_lookup)
    for row in table3_wide:
        if row["paper_compound_name"] not in reported_names:
            reported_names.append(row["paper_compound_name"])

    canonical_order = []
    canonical_to_reported: dict[str, list[str]] = {}
    for reported_name in reported_names:
        canonical_name = PAPER_NAME_CANONICAL.get(reported_name, reported_name)
        if canonical_name not in canonical_order:
            canonical_order.append(canonical_name)
        canonical_to_reported.setdefault(canonical_name, []).append(reported_name)

    registry = []
    query_log = []
    for index, paper_name in enumerate(canonical_order, start=1):
        names_as_reported = canonical_to_reported[paper_name]
        aliases = PUBCHEM_QUERY_ALIASES.get(paper_name, [paper_name])
        result = None
        used_alias = ""
        final_url = ""
        final_status = ""
        attempts = []
        for alias in aliases:
            result, final_url, final_status = fetch_pubchem_property(alias)
            attempts.append(f"{alias}:{final_status}")
            time.sleep(0.27)
            if result:
                used_alias = alias
                break

        status = "resolved" if result else "unresolved"
        note = ""
        if paper_name in STRUCTURE_REVIEW_NOTES:
            status, note = STRUCTURE_REVIEW_NOTES[paper_name]
        elif not result:
            note = "No PubChem CID was resolved with the recorded aliases; manual structure curation is required."

        panel_name = next((name for name in names_as_reported if name in panel_lookup), "")
        panel_row = panel_lookup.get(panel_name, {})
        registry.append(
            {
                "compound_id": f"ADC_PAYLOAD_{index:03d}",
                "paper_compound_name": paper_name,
                "paper_names_as_reported": " | ".join(names_as_reported),
                "in_roth_27_panel": "yes" if panel_name else "no",
                "in_roth_table3_confirmatory": "yes" if any(r["paper_compound_name"] in names_as_reported for r in table3_wide) else "no",
                "vendor_as_reported": panel_row.get("vendor", ""),
                "mechanism_as_reported": panel_row.get("mechanism_of_action", ""),
                "pubchem_query_name": used_alias or aliases[-1],
                "pubchem_cid": result.get("CID", "") if result else "",
                "pubchem_title": result.get("Title", "") if result else "",
                "isomeric_smiles": result.get("SMILES", result.get("IsomericSMILES", "")) if result else "",
                "connectivity_smiles": result.get("ConnectivitySMILES", result.get("CanonicalSMILES", "")) if result else "",
                "inchikey": result.get("InChIKey", "") if result else "",
                "molecular_formula": result.get("MolecularFormula", "") if result else "",
                "molecular_weight": result.get("MolecularWeight", "") if result else "",
                "xlogp": result.get("XLogP", "") if result else "",
                "tpsa": result.get("TPSA", "") if result else "",
                "hbond_donor_count": result.get("HBondDonorCount", "") if result else "",
                "hbond_acceptor_count": result.get("HBondAcceptorCount", "") if result else "",
                "rotatable_bond_count": result.get("RotatableBondCount", "") if result else "",
                "formal_charge": result.get("Charge", "") if result else "",
                "structure_status": status,
                "structure_note": note,
                "structure_source_url": final_url,
                "source_doi": ROTH_DOI,
            }
        )
        query_log.append(
            {
                "paper_compound_name": paper_name,
                "aliases_attempted": " | ".join(attempts),
                "resolved_query_name": used_alias,
                "http_status_final": final_status,
                "pubchem_cid": result.get("CID", "") if result else "",
                "query_url_final": final_url,
                "queried_at": "2026-07-26",
            }
        )
    return registry, query_log


def extract_daood() -> tuple[list[dict], list[dict], list[dict]]:
    workbook = openpyxl.load_workbook(DAOOD_XLSX, read_only=True, data_only=True)
    sheet = workbook["Table S2"]
    headers = [cell.value for cell in sheet[3]]
    raw_rows = []
    long_rows = []
    endpoint_counts: dict[str, Counter] = {name: Counter() for name in headers[4:]}
    source_counts = Counter()

    endpoint_map = {"Pgp_substrate": ("ABCB1", "P-gp"), "BCRP_substrate": ("ABCG2", "BCRP")}
    record_index = 0
    for values in sheet.iter_rows(min_row=4, values_only=True):
        row = dict(zip(headers, values))
        raw_rows.append(row)
        source_counts[row["Source"]] += 1
        for endpoint in headers[4:]:
            if row[endpoint] is not None:
                endpoint_counts[endpoint][int(row[endpoint])] += 1
        for endpoint, (gene, protein) in endpoint_map.items():
            if row[endpoint] is None:
                continue
            record_index += 1
            long_rows.append(
                {
                    "training_record_id": f"DAOOD2025_SUB_{record_index:05d}",
                    "pid": row["PID"],
                    "source_database": row["Source"],
                    "preferred_name": row["Preferred_name"],
                    "smiles_as_published": row["SMILES"],
                    "transporter_gene": gene,
                    "transporter_protein": protein,
                    "endpoint": "substrate",
                    "published_binary_label": int(row[endpoint]),
                    "label_source": "Daood 2025 Table S2 curated training set",
                    "processing_status": "raw_as_published_not_yet_deduplicated_against_adc_test_set",
                    "source_doi": DAOOD_DOI,
                    "source_url": DAOOD_URL,
                }
            )

    count_rows = []
    for endpoint, counts in endpoint_counts.items():
        count_rows.append(
            {
                "dataset": "Daood 2025 Table S2",
                "endpoint": endpoint,
                "negative_n": counts[0],
                "positive_n": counts[1],
                "labeled_total_n": counts[0] + counts[1],
                "unique_table_rows_total": len(raw_rows),
                "source_doi": DAOOD_DOI,
            }
        )
    for source, count in sorted(source_counts.items()):
        count_rows.append(
            {
                "dataset": "Daood 2025 Table S2 source composition",
                "endpoint": source,
                "negative_n": "",
                "positive_n": "",
                "labeled_total_n": count,
                "unique_table_rows_total": len(raw_rows),
                "source_doi": DAOOD_DOI,
            }
        )
    return raw_rows, long_rows, count_rows


def normalize_literal_name(value: str | None) -> str:
    """Normalize names for a conservative literal overlap screen, not identity resolution."""
    ascii_text = unicodedata.normalize("NFKD", value or "").encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", "", ascii_text.lower())


def build_preliminary_overlap_audit(registry: list[dict], daood_rows: list[dict]) -> list[dict]:
    """Flag exact published-SMILES or normalized-name overlap with the ADC registry.

    This deliberately does not canonicalize salts, tautomers, stereochemistry, or SMILES.
    A negative result here therefore cannot be interpreted as proof that a compound is absent.
    """
    audit_rows = []
    overlap_index = 0
    for compound in registry:
        reported_names = compound["paper_names_as_reported"].split(" | ")
        candidate_names = {
            normalize_literal_name(name)
            for name in [
                compound["paper_compound_name"],
                *reported_names,
                compound["pubchem_query_name"],
                compound["pubchem_title"],
            ]
            if name
        }
        candidate_smiles = {
            value
            for value in [compound["isomeric_smiles"], compound["connectivity_smiles"]]
            if value
        }
        compound_hits = []
        for row in daood_rows:
            match_types = []
            if normalize_literal_name(row["Preferred_name"]) in candidate_names:
                match_types.append("normalized_name")
            if row["SMILES"] and row["SMILES"] in candidate_smiles:
                match_types.append("exact_smiles")
            if match_types:
                compound_hits.append((row, "+".join(match_types)))

        for hit_number, (row, match_kind) in enumerate(compound_hits, start=1):
            overlap_index += 1
            substrate_labels = [row["Pgp_substrate"], row["BCRP_substrate"]]
            audit_rows.append(
                {
                    "overlap_record_id": f"OVERLAP_{overlap_index:03d}",
                    "compound_id": compound["compound_id"],
                    "adc_compound_name": compound["paper_compound_name"],
                    "unique_adc_entity_marker": "yes" if hit_number == 1 else "",
                    "adc_structure_status": compound["structure_status"],
                    "daood_pid": row["PID"],
                    "daood_source": row["Source"],
                    "daood_preferred_name": row["Preferred_name"],
                    "match_kind": match_kind,
                    "pgp_substrate_label": row["Pgp_substrate"],
                    "bcrp_substrate_label": row["BCRP_substrate"],
                    "has_pgp_or_bcrp_substrate_label": (
                        "yes" if any(value not in ("", None) for value in substrate_labels) else "no"
                    ),
                    "audit_scope": "normalized literal name and exact published SMILES only",
                    "leakage_action": "remove_or_flag_before_external_evaluation",
                    "audit_note": "Preliminary screen; canonical InChIKey, salt, tautomer, stereochemistry and similarity audit remains required.",
                }
            )
    return audit_rows


def check_pubchem_assays() -> list[dict]:
    aids = [2061167, 2061304, 2061305, 2061306, 2061317, 2061318, 2061319, 2061320, 2061321, 2061322, 2061323, 2061324, 2202239, 2202241]
    rows = []
    for aid in aids:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/CSV"
        request = urllib.request.Request(url, headers={"User-Agent": "ADC-substrate-data-curation/1.0"})
        status = ""
        detail = ""
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                payload = response.read()
                status = str(response.status)
                detail = f"downloadable_csv_bytes={len(payload)}"
        except urllib.error.HTTPError as exc:
            status = str(exc.code)
            detail = exc.read().decode("utf-8", errors="replace")[:240]
        except Exception as exc:  # noqa: BLE001
            status = f"error:{type(exc).__name__}"
            detail = str(exc)[:240]
        rows.append(
            {
                "aid": aid,
                "role": "project_summary" if aid == 2061167 else ("paper_declared_second_aid" if aid == 2202241 else "related_assay"),
                "pug_rest_csv_http_status": status,
                "detail": detail,
                "checked_at": "2026-07-26",
                "url": url,
            }
        )
        time.sleep(0.27)
    return rows


def source_inventory() -> list[dict]:
    rows = []
    for path in sorted((RAW_ROTH).glob("*")) + sorted((RAW_DAOOD).glob("*")):
        if not path.is_file():
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(
            {
                "relative_path": str(path.relative_to(ROOT)),
                "size_bytes": path.stat().st_size,
                "sha256": digest,
                "source_group": "Roth 2026" if "roth_2026" in str(path) else "Daood 2025",
                "acquired_at": "2026-07-26",
            }
        )
    return rows


def main() -> None:
    table1, table2, table3_wide, table3_long = parse_roth_tables()
    registry, pubchem_log = build_compound_registry(table2, table3_wide)
    name_to_id = {}
    for compound in registry:
        for reported_name in compound["paper_names_as_reported"].split(" | "):
            name_to_id[reported_name] = compound["compound_id"]
    for row in table2:
        row["compound_id"] = name_to_id[row["paper_compound_name"]]
    for row in table3_wide:
        row["compound_id"] = name_to_id[row["paper_compound_name"]]
    for row in table3_long:
        row["compound_id"] = name_to_id[row["paper_compound_name"]]
    daood_raw, daood_long, endpoint_counts = extract_daood()
    overlap_audit = build_preliminary_overlap_audit(registry, daood_raw)
    overlap_pids = {row["daood_pid"] for row in overlap_audit}
    for row in daood_long:
        if row["pid"] in overlap_pids:
            row["processing_status"] = "preliminary_exact_overlap_flag_remove_or_holdout"
    assay_availability = check_pubchem_assays()

    write_csv(PROCESSED / "roth_2026_table1_adc_context.csv", table1)
    write_csv(PROCESSED / "roth_2026_table2_payload_panel_27.csv", table2)
    write_csv(PROCESSED / "roth_2026_table3_confirmatory_wide_13.csv", table3_wide)
    write_csv(PROCESSED / "roth_2026_table3_transporter_records_long_39.csv", table3_long)
    write_csv(
        PROCESSED / "external_truth_set_v1_pgp_bcrp_continuous_26.csv",
        [row for row in table3_long if row["transporter_gene"] in {"ABCB1", "ABCG2"}],
    )
    write_csv(PROCESSED / "compound_registry_v1.csv", registry)
    write_csv(PROCESSED / "pubchem_structure_query_log.csv", pubchem_log)
    write_csv(PROCESSED / "pubchem_assay_availability_log.csv", assay_availability)
    write_csv(PROCESSED / "daood_2025_tableS2_all_endpoints_8794.csv", daood_raw)
    write_csv(PROCESSED / "daood_2025_pgp_bcrp_substrate_training_long.csv", daood_long)
    write_csv(PROCESSED / "daood_2025_endpoint_counts.csv", endpoint_counts)
    write_csv(PROCESSED / "preliminary_exact_overlap_audit.csv", overlap_audit)
    write_csv(PROCESSED / "source_inventory_sha256.csv", source_inventory())

    summary = {
        "generated_at": "2026-07-26",
        "roth_table1_adc_rows": len(table1),
        "roth_table2_panel_rows": len(table2),
        "roth_table3_compounds": len(table3_wide),
        "roth_table3_all_transporter_records": len(table3_long),
        "roth_external_pgp_bcrp_records": sum(r["transporter_gene"] in {"ABCB1", "ABCG2"} for r in table3_long),
        "compound_registry_rows": len(registry),
        "compound_structures_resolved": sum(r["pubchem_cid"] != "" for r in registry),
        "compound_structures_unresolved": sum(r["pubchem_cid"] == "" for r in registry),
        "daood_unique_compound_rows": len(daood_raw),
        "daood_pgp_bcrp_substrate_labeled_records": len(daood_long),
        "preliminary_exact_overlap_rows": len(overlap_audit),
        "preliminary_exact_overlap_adc_entities": sum(r["unique_adc_entity_marker"] == "yes" for r in overlap_audit),
        "preliminary_exact_overlap_labeled_training_records": sum(
            r["processing_status"].startswith("preliminary_exact_overlap") for r in daood_long
        ),
        "warning": "No binary substrate labels were assigned to Roth Table 3. Continuous RR is preserved pending a preregistered rule and author-interpretation audit.",
    }
    (PROCESSED / "first_batch_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build an ABCB1 regulatory panel by reusing locally archived PMDA reviews.

The source documents were originally archived for an ABCG2 confirmation task,
but many contain direct human P-gp substrate experiments. This script records
the ABCB1 rows separately, resolves structures, and performs the same prior-
identity audit used by the project. No model scoring occurs here.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from abcb1_screener_core_v0_1 import standardize_smiles  # noqa: E402


PROCESSED = PROJECT / "数据/processed"
RAW_DIR = PROJECT / "数据/raw/human_abcg2_prospective_confirmation"
PUBCHEM_DIR = PROJECT / "数据/raw/abcb1_regulatory_reuse/pubchem"

PRIOR_DATASETS = {
    "frozen_443_training": (PROCESSED / "abcb1_screener_v0_1_training_audit.csv", ["standardized_inchikey"]),
    "soskuti2024_development": (PROCESSED / "soskuti_2024_mdr1_species_panel_v1.csv", ["standardized_inchikey"]),
    "q1_external_development": (
        PROCESSED / "abcb1_q1_direct_external_entity_predictions_v1.csv",
        ["model_standardized_inchikey", "standardized_inchikey"],
    ),
    "kelly2022_postfreeze_stress": (
        PROCESSED / "abcb1_kelly2022_table7_external_stress_v1.csv",
        ["standardized_inchikey"],
    ),
    "broccatelli2012_new_candidate_pool": (
        PROCESSED / "broccatelli_2012_abcb1_entities_v1.csv",
        ["standardized_inchikey"],
    ),
}


SOURCES = {
    "PMDA_APREMILAST": {
        "url": "https://www.pmda.go.jp/files/000230195.pdf",
        "local_pdf": "apremilast_pmda_primary.pdf",
        "local_text": "apremilast_pmda_primary.txt",
        "location": "Section 4.5.2 Transporters; PDF review page 19; extracted text lines 1291-1303",
    },
    "PMDA_BICTEGRAVIR_FDC": {
        "url": "https://www.pmda.go.jp/drugs/2019/P20190423002/230867000_23100AMX00302_I100_1.pdf",
        "local_pdf": "bictegravir_pmda_primary.pdf",
        "local_text": "bictegravir_pmda_primary.txt",
        "location": "Section 2.6.5 Tables AD-141-2278, AD-120-2018, and AD-236-2004; extracted text lines 11175-11191, 11473-11492, 12272-12296",
    },
    "PMDA_CRIZOTINIB": {
        "url": "https://www.pmda.go.jp/files/000153949.pdf",
        "local_pdf": "crizotinib_pmda_primary.pdf",
        "local_text": "crizotinib_pmda_primary.txt",
        "location": "Section 3(ii).A.(5).3)a; extracted text lines 1568-1573",
    },
    "PMDA_DORAVIRINE": {
        "url": "https://www.pmda.go.jp/drugs/2019/P20191223002/170050000_30200AMX00002_I100_1.pdf",
        "local_pdf": "doravirine_pmda_primary.pdf",
        "local_text": "doravirine_pmda_primary.txt",
        "location": "Section Efflux of Doravirine by BCRP and P-gp; extracted text lines 1932-1936",
    },
    "PMDA_OLODATEROL": {
        "url": "https://www.pmda.go.jp/files/000153757.pdf",
        "local_pdf": "olodaterol_pmda_primary.pdf",
        "local_text": "olodaterol_pmda_primary.txt",
        "location": "Section 3(ii).A.(5).2 Transporters; extracted text lines 1040-1047",
    },
    "PMDA_QUIZARTINIB": {
        "url": "https://www.pmda.go.jp/files/000240288.pdf",
        "local_pdf": "quizartinib_pmda_primary.pdf",
        "local_text": "quizartinib_pmda_primary.txt",
        "location": "Section 4.5.3 Transporters; extracted text lines 917-923",
    },
    "PMDA_ROXADUSTAT": {
        "url": "https://www.pmda.go.jp/drugs/2019/P20191007001/800126000_30100AMX00239_K100_1.pdf",
        "local_pdf": "roxadustat_pmda_primary.pdf",
        "local_text": "roxadustat_pmda_primary.txt",
        "location": "Section 2.7.2.2.1.4.2; extracted text lines 3892-3908",
    },
    "PMDA_SEPETAPROST": {
        "url": "https://www.pmda.go.jp/files/000280580.pdf",
        "local_pdf": "sepetaprost_pmda_primary.pdf",
        "local_text": "sepetaprost_pmda_primary.txt",
        "location": "Section (5) Transport by drug transporters; extracted text lines 1336-1365",
    },
    "PMDA_DEUCRAVACITINIB": {
        "url": "https://www.pmda.go.jp/drugs/2022/P20221006002/670605000_30400AMX00412_I100_1.pdf",
        "local_pdf": "deucravacitinib_pmda_primary.pdf",
        "local_text": "deucravacitinib_pmda_primary.txt",
        "location": "Table 2.6.5.16F BMS-986165 MDR1 transporter kinetics; extracted text lines 5644-5669",
    },
    "PMDA_OSIMERTINIB": {
        "url": "https://www.pmda.go.jp/drugs/2016/P20160411001/670227000_22800AMX00385_I100_1.pdf",
        "local_pdf": "osimertinib_pmda_primary.pdf",
        "local_text": "osimertinib_pmda_primary.txt",
        "location": "Table 35 and detailed MDCK substrate tables; extracted text lines 2741-2773 and 5915-5957",
    },
    "PMDA_ASUNAPREVIR": {
        "url": "https://www.pmda.go.jp/files/000209021.pdf",
        "local_pdf": "asunaprevir_pmda_primary.pdf",
        "local_text": "asunaprevir_pmda_primary.txt",
        "location": "Transporter section; extracted text lines 2018-2024",
    },
}


ENTITIES = [
    {
        "entity_name": "Apremilast", "source_key": "PMDA_APREMILAST", "label": 1,
        "assay_system": "human P-gp-transfected LLC-PK1 versus control LLC-PK1", "concentration_uM": 10,
        "er_transporter": 19, "er_parental": 1.1, "net_er": 19 / 1.1, "inhibitor": "ketoconazole 30 uM",
        "inhibitor_er": "", "inhibitor_relative_reduction": 0.92,
        "strict_matched_parent_candidate": "yes", "qc_note": "direct parent comparison and 92% inhibition",
    },
    {
        "entity_name": "Bictegravir", "source_key": "PMDA_BICTEGRAVIR_FDC", "label": 1,
        "assay_system": "human P-gp-overexpressing MDCKII versus wild-type MDCKII", "concentration_uM": 10,
        "er_transporter": 7.5, "er_parental": 1.3, "net_er": 7.5 / 1.3, "inhibitor": "cyclosporin A 10 uM",
        "inhibitor_er": 2.4, "inhibitor_relative_reduction": (7.5 - 2.4) / 7.5,
        "strict_matched_parent_candidate": "yes", "qc_note": "full Papp table and inhibitor condition",
    },
    {
        "entity_name": "Tenofovir alafenamide", "source_key": "PMDA_BICTEGRAVIR_FDC", "label": 1,
        "assay_system": "human MDR1-transfected MDCKII versus wild-type MDCKII", "concentration_uM": 8,
        "er_transporter": 66.2, "er_parental": 4.8, "net_er": 66.2 / 4.8, "inhibitor": "cyclosporin A 10 uM",
        "inhibitor_er": 5.6, "inhibitor_relative_reduction": (66.2 - 5.6) / 66.2,
        "strict_matched_parent_candidate": "yes", "qc_note": "full bidirectional Papp, recovery, parent and inhibitor rows",
    },
    {
        "entity_name": "Tenofovir", "source_key": "PMDA_BICTEGRAVIR_FDC", "label": 0,
        "assay_system": "human MDR1-transfected MDCKII versus wild-type MDCKII", "concentration_uM": 13,
        "er_transporter": 1.0, "er_parental": 1.5, "net_er": 1.0 / 1.5, "inhibitor": "cyclosporin A 10 uM",
        "inhibitor_er": 1.5, "inhibitor_relative_reduction": 0,
        "strict_matched_parent_candidate": "yes", "qc_note": "full bidirectional Papp, recovery, parent and inhibitor rows",
    },
    {
        "entity_name": "Crizotinib", "source_key": "PMDA_CRIZOTINIB", "label": 1,
        "assay_system": "human P-gp-overexpressing MDCK versus non-P-gp MDCK", "concentration_uM": "0.1-20",
        "er_transporter": "3.9-33.6", "er_parental": "reported lower without exact range", "net_er": "",
        "inhibitor": "none reported in extracted substrate result", "inhibitor_er": "", "inhibitor_relative_reduction": "",
        "strict_matched_parent_candidate": "yes", "qc_note": "multi-concentration direct transfectant range; source explicit substrate call",
    },
    {
        "entity_name": "Doravirine", "source_key": "PMDA_DORAVIRINE", "label": 1,
        "assay_system": "human MDR1-transfected LLC-PK1 versus non-transfected LLC-PK1", "concentration_uM": "0.1,0.5,1",
        "er_transporter": "3.5-4.7", "er_parental": "0.9-1.0", "net_er": "3.5-5.22", "inhibitor": "none reported",
        "inhibitor_er": "", "inhibitor_relative_reduction": "", "strict_matched_parent_candidate": "yes",
        "qc_note": "multi-concentration direct matched-parent range",
    },
    {
        "entity_name": "Olodaterol", "source_key": "PMDA_OLODATEROL", "label": 1,
        "assay_system": "Caco-2 endogenous efflux", "concentration_uM": 10, "er_transporter": 43.5,
        "er_parental": "", "net_er": "", "inhibitor": "cyclosporin or verapamil", "inhibitor_er": "approximately 1",
        "inhibitor_relative_reduction": "approximately 0.98", "strict_matched_parent_candidate": "no",
        "qc_note": "strong inhibitor rescue but endogenous multi-transporter epithelium",
    },
    {
        "entity_name": "Quizartinib", "source_key": "PMDA_QUIZARTINIB", "label": 1,
        "assay_system": "human P-gp-expressing MDCKII versus non-expressing MDCKII", "concentration_uM": 3,
        "er_transporter": 5.8, "er_parental": "implied by net ratio", "net_er": 3.6,
        "inhibitor": "verapamil 200 uM or ketoconazole 50 uM", "inhibitor_er": "1.5 or 2.5",
        "inhibitor_relative_reduction": (5.8 - 1.5) / 5.8, "strict_matched_parent_candidate": "yes",
        "qc_note": "matched non-expressing comparison and two inhibitor conditions",
    },
    {
        "entity_name": "Roxadustat", "source_key": "PMDA_ROXADUSTAT", "label": 0,
        "assay_system": "human MDR1-expressing LLC-PK1 versus control LLC-PK1", "concentration_uM": "1,10,100",
        "er_transporter": "", "er_parental": "", "net_er": "<2 at almost all concentrations", "inhibitor": "verapamil control",
        "inhibitor_er": "", "inhibitor_relative_reduction": "", "strict_matched_parent_candidate": "yes",
        "qc_note": "source explicit non-substrate call; positive control behaved as expected",
    },
    {
        "entity_name": "Sepetaprost", "source_key": "PMDA_SEPETAPROST", "label": 0,
        "assay_system": "human P-gp-expressing MDCKII versus non-expressing MDCKII", "concentration_uM": "0.2-20",
        "er_transporter": "", "er_parental": "", "net_er": "0.34-1.36", "inhibitor": "none needed for negative range",
        "inhibitor_er": "", "inhibitor_relative_reduction": "", "strict_matched_parent_candidate": "yes",
        "qc_note": "multi-concentration net ER range entirely below 2",
    },
    {
        "entity_name": "ONO-AG-367", "source_key": "PMDA_SEPETAPROST", "label": 0,
        "assay_system": "human P-gp-expressing MDCKII versus non-expressing MDCKII", "concentration_uM": "0.2-20",
        "er_transporter": "", "er_parental": "", "net_er": "0.56-0.68", "inhibitor": "none needed for negative range",
        "inhibitor_er": "", "inhibitor_relative_reduction": "", "strict_matched_parent_candidate": "yes",
        "qc_note": "multi-concentration net ER range entirely below 2",
    },
    {
        "entity_name": "Deucravacitinib", "source_key": "PMDA_DEUCRAVACITINIB", "label": 1,
        "assay_system": "human MDR1-expressing MDCKII versus wild-type MDCKII", "concentration_uM": "0.07-150",
        "er_transporter": "9.02-96.2", "er_parental": "", "net_er": "2.59-42.2", "inhibitor": "none in kinetic table",
        "inhibitor_er": "", "inhibitor_relative_reduction": "", "strict_matched_parent_candidate": "yes",
        "qc_note": "eight concentrations, triplicates, direct net ER and positive control",
    },
    {
        "entity_name": "Osimertinib", "source_key": "PMDA_OSIMERTINIB", "label": 1,
        "assay_system": "human MDR1-MDCK versus non-expressing MDCK", "concentration_uM": 1,
        "er_transporter": 13.4, "er_parental": 1.55, "net_er": 8.69, "inhibitor": "valspodar 1 uM",
        "inhibitor_er": 1.01, "inhibitor_relative_reduction": (13.4 - 1.01) / 13.4,
        "strict_matched_parent_candidate": "yes", "qc_note": "full Papp, recovery, parent and inhibitor table",
    },
    {
        "entity_name": "AZ5104", "source_key": "PMDA_OSIMERTINIB", "label": 1,
        "assay_system": "human MDR1-MDCK versus non-expressing MDCK", "concentration_uM": 1,
        "er_transporter": 186, "er_parental": "", "net_er": 57.4, "inhibitor": "valspodar 1 uM",
        "inhibitor_er": 0.875, "inhibitor_relative_reduction": (186 - 0.875) / 186,
        "strict_matched_parent_candidate": "yes", "qc_note": "parent-corrected net ER and inhibitor rescue",
    },
    {
        "entity_name": "AZ7550", "source_key": "PMDA_OSIMERTINIB", "label": 1,
        "assay_system": "human MDR1-MDCK versus non-expressing MDCK", "concentration_uM": "10,30",
        "er_transporter": "11.9,2.38", "er_parental": "4.49,1.15", "net_er": "2.65,2.06",
        "inhibitor": "valspodar condition not quantifiable at 1 uM", "inhibitor_er": "ND", "inhibitor_relative_reduction": "",
        "strict_matched_parent_candidate": "yes", "qc_note": "positive at 10 and 30 uM; 1 uM A-to-B below LLOQ",
    },
    {
        "entity_name": "Asunaprevir", "source_key": "PMDA_ASUNAPREVIR", "label": 1,
        "assay_system": "wild-type MDCK plus mouse P-gp knockout context", "concentration_uM": 1,
        "er_transporter": 12.1, "er_parental": "", "net_er": "", "inhibitor": "GF120918/MK-571 context",
        "inhibitor_er": "", "inhibitor_relative_reduction": "", "strict_matched_parent_candidate": "no",
        "qc_note": "P-gp/MRP2 attribution but not a human P-gp matched-transfectant substrate assay",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def identity_index(path: Path, columns: list[str]) -> tuple[set[str], set[str]]:
    exact: set[str] = set()
    connectivity: set[str] = set()
    for row in read_csv(path):
        key = next((row.get(column, "") for column in columns if row.get(column, "")), "")
        if key:
            exact.add(key)
            connectivity.add(key.split("-")[0])
    return exact, connectivity


def safe_cache_name(name: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "compound"
    return f"{stem}_{hashlib.sha1(name.encode('utf-8')).hexdigest()[:10]}.json"


def fetch_pubchem(name: str) -> tuple[dict, str, str]:
    PUBCHEM_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = PUBCHEM_DIR / safe_cache_name(name)
    if cache_path.exists():
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
        status = "cached_pubchem_exact_name"
    else:
        encoded = urllib.parse.quote(name, safe="")
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{encoded}/property/CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,MolecularWeight,XLogP,TPSA/JSON"
        )
        request = urllib.request.Request(url, headers={"User-Agent": "ABCB1-research-audit/1.0"})
        payload = {}
        last_error = ""
        for attempt in range(4):
            try:
                with urllib.request.urlopen(request, timeout=40) as response:
                    payload = json.loads(response.read().decode("utf-8"))
                cache_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
                status = "live_pubchem_exact_name"
                time.sleep(0.12)
                break
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
                last_error = f"{type(exc).__name__}: {exc}"
                if attempt == 3:
                    return {}, "pubchem_query_failed", last_error
                time.sleep(1.5 * (attempt + 1))
    properties = payload.get("PropertyTable", {}).get("Properties", [])
    if not properties:
        return {}, "pubchem_no_property_record", "no PropertyTable.Properties result"
    prop = properties[0]
    smiles = prop.get("SMILES") or prop.get("IsomericSMILES") or prop.get("CanonicalSMILES") or ""
    return {
        "pubchem_cid": prop.get("CID", ""),
        "pubchem_inchikey": prop.get("InChIKey", ""),
        "pubchem_smiles": smiles,
        "molecular_formula": prop.get("MolecularFormula", ""),
        "molecular_weight": prop.get("MolecularWeight", ""),
        "xlogp": prop.get("XLogP", ""),
        "tpsa": prop.get("TPSA", ""),
        "pubchem_cache_file": str(cache_path.relative_to(PROJECT)),
    }, status, ""


def main() -> None:
    prior_indexes = {
        name: identity_index(path, columns)
        for name, (path, columns) in PRIOR_DATASETS.items()
    }
    entity_rows: list[dict] = []
    for index, manual in enumerate(ENTITIES, start=1):
        source = SOURCES[manual["source_key"]]
        pubchem, pubchem_status, pubchem_error = fetch_pubchem(manual["entity_name"])
        standardized = standardize_smiles(pubchem.get("pubchem_smiles", ""))
        key = standardized.get("standardized_inchikey", "")
        block = key.split("-")[0] if key else ""
        exact_hits = [name for name, (exact, _) in prior_indexes.items() if key and key in exact]
        connectivity_hits = [name for name, (_, conn) in prior_indexes.items() if block and block in conn]
        if standardized.get("standardization_status") != "ok":
            identity_status = "exclude_unresolved_structure"
        elif exact_hits:
            identity_status = "prior_exact_identity_overlap"
        elif connectivity_hits:
            identity_status = "prior_connectivity_identity_overlap"
        else:
            identity_status = "no_exact_or_connectivity_prior_overlap"

        strict = manual["strict_matched_parent_candidate"] == "yes"
        candidate_eligible = strict and identity_status == "no_exact_or_connectivity_prior_overlap"
        local_pdf = RAW_DIR / source["local_pdf"]
        entity_rows.append(
            {
                "regulatory_entity_id": f"ABCB1_REG_REUSE_{index:03d}",
                **manual,
                **pubchem,
                "pubchem_resolution_status": pubchem_status,
                "pubchem_resolution_error": pubchem_error,
                **standardized,
                "inchikey_connectivity_block": block,
                "overlap_exact_prior_sets": " | ".join(exact_hits),
                "overlap_connectivity_prior_sets": " | ".join(connectivity_hits),
                "identity_status": identity_status,
                "strict_external_candidate_eligible": "yes" if candidate_eligible else "no",
                "source_id": manual["source_key"],
                "source_url": source["url"],
                "source_location": source["location"],
                "local_pdf": str(local_pdf.relative_to(PROJECT)),
                "local_pdf_sha256": sha256(local_pdf) if local_pdf.exists() else "missing",
                "local_text": str((RAW_DIR / source["local_text"]).relative_to(PROJECT)),
                "evidence_tier": (
                    "A1_regulatory_direct_human_transfectant_matched_control"
                    if strict
                    else "B1_regulatory_supporting_nonmatched_context"
                ),
                "claim_boundary": "regulatory source-call candidate; not model validation until prelocked scoring and source-independence checks",
            }
        )

    write_csv(PROCESSED / "abcb1_regulatory_reuse_entities_v1.csv", entity_rows)
    write_csv(PROCESSED / "abcb1_regulatory_reuse_observations_v1.csv", entity_rows)

    source_rows = []
    for source_id, source in SOURCES.items():
        local_pdf = RAW_DIR / source["local_pdf"]
        rows = [row for row in entity_rows if row["source_id"] == source_id]
        source_rows.append(
            {
                "source_id": source_id,
                "source_url": source["url"],
                "source_location": source["location"],
                "local_pdf": str(local_pdf.relative_to(PROJECT)),
                "local_pdf_sha256": sha256(local_pdf) if local_pdf.exists() else "missing",
                "entity_count": len(rows),
                "strict_candidate_count": sum(row["strict_external_candidate_eligible"] == "yes" for row in rows),
                "allowed_use": "ABCB1 regulatory direct-evidence candidate collection and future frozen scoring",
            }
        )
    write_csv(PROCESSED / "abcb1_regulatory_reuse_source_registry_v1.csv", source_rows)

    strict_rows = [row for row in entity_rows if row["strict_external_candidate_eligible"] == "yes"]
    summary = {
        "source_documents": len(source_rows),
        "entity_rows": len(entity_rows),
        "standardization_ok": sum(row["standardization_status"] == "ok" for row in entity_rows),
        "strict_nonoverlap_candidates": len(strict_rows),
        "strict_candidate_class_counts": dict(Counter(str(row["label"]) for row in strict_rows)),
        "secondary_context_rows": sum(row["strict_matched_parent_candidate"] == "no" for row in entity_rows),
        "formal_status": "candidate_collection_complete_no_model_scoring_in_this_script",
    }
    (PROCESSED / "abcb1_regulatory_reuse_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

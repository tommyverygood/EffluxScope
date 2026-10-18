#!/usr/bin/env python3
"""Collect and tier the open ASAP MDCK-MDR1 ChEMBL datasets.

This collector preserves activity-level source values and reconstructs one
wide observation per assay/molecule/source-record/replicate.  It deliberately
keeps three evidence roles separate:

* bidirectional MDCK-MDR1 with inhibitor: weak/continuous ABCB1 evidence;
* bidirectional MDCK-MDR1 without inhibitor: high/low-efflux continuous data;
* unidirectional Papp only: permeability/QC data, not an ABCB1 substrate label.

No predictor is loaded or scored.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from abcb1_screener_core_v0_1 import FP_GENERATOR, standardize_smiles  # noqa: E402
from collect_abcb1_regulatory_reuse_panel_v1 import identity_index  # noqa: E402


PROCESSED = PROJECT / "数据/processed"
RAW_DIR = PROJECT / "数据/raw/abcb1_asap_mdck_mdr1_v1"

DOCUMENTS = {
    "CHEMBL5442303": {
        "release": "CHEMBL_35",
        "year": 2024,
        "doi": "10.6019/CHEMBL5442303",
        "predeclared_role": "auxiliary_development_or_source_heldout_sensitivity_only",
    },
    "CHEMBL5674852": {
        "release": "CHEMBL_36",
        "year": 2025,
        "doi": "10.6019/CHEMBL5674852",
        "predeclared_role": "reserve_for_source_heldout_continuous_evaluation_before_any_scoring",
    },
}

PRIOR_DATASETS = {
    "frozen_443_training": (
        PROCESSED / "abcb1_screener_v0_1_training_audit.csv",
        ["standardized_inchikey", "model_standardized_inchikey"],
    ),
    "soskuti2024_direct_development": (
        PROCESSED / "soskuti_2024_mdr1_species_panel_v1.csv",
        ["standardized_inchikey"],
    ),
    "q1_historical_external": (
        PROCESSED / "abcb1_q1_direct_external_entity_predictions_v1.csv",
        ["model_standardized_inchikey", "standardized_inchikey"],
    ),
    "regulatory_strict_v4_scored": (
        PROCESSED / "abcb1_regulatory_strict_external_pool_v4.csv",
        ["standardized_inchikey"],
    ),
    "regulatory_batch5_reserved": (
        PROCESSED / "abcb1_regulatory_strict_increment_batch5_reserved_unscored_v1.csv",
        ["standardized_inchikey"],
    ),
    "npmmp_weak_pairs": (
        PROCESSED / "abcb1_npmmp_inhibitor_pair_reconstruction_v1.csv",
        ["standardized_inchikey"],
    ),
}

ACTIVITY_FIELDS = [
    "document_chembl_id",
    "assay_chembl_id",
    "record_id",
    "activity_id",
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
    "molecule_pref_name",
    "canonical_smiles",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "data_validity_comment",
    "potential_duplicate",
]

ENDPOINT_MAP = {
    "Papp A to B (mean)": "papp_ab",
    "Papp B to A (mean)": "papp_ba",
    "Recovery A to B (mean)": "recovery_ab_pct",
    "Recovery B to A (mean)": "recovery_ba_pct",
    "Efflux ratio": "efflux_ratio_source",
    "Papp A to B plus inhibitor (mean)": "papp_ab_inhibitor",
    "Papp B to A plus inhibitor (mean)": "papp_ba_inhibitor",
    "Recovery A to B plus inhibitor (mean)": "recovery_ab_inhibitor_pct",
    "Recovery B to A plus inhibitor (mean)": "recovery_ba_inhibitor_pct",
    "Efflux ratio plus inhibitor": "efflux_ratio_inhibitor_source",
}


def fetch_json(url: str, params: dict | None = None) -> dict:
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(url, headers={"User-Agent": "ABCB1-data-audit/1.0"})
    with urllib.request.urlopen(request, timeout=90) as response:
        return json.load(response)


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def number(value: object) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def relation_value(group: pd.DataFrame, endpoint: str) -> tuple[str, float | None, str, str]:
    row = group[group["standard_type"].eq(endpoint)]
    if row.empty:
        return "", None, "", ""
    source = row.iloc[0]
    return (
        str(source.get("standard_relation") or ""),
        number(source.get("standard_value")),
        str(source.get("standard_units") or ""),
        str(source.get("activity_id") or ""),
    )


def load_prior_fingerprints(path: Path, smiles_columns: list[str]) -> list:
    if not path.exists():
        return []
    frame = pd.read_csv(path, low_memory=False)
    fps = []
    for col in smiles_columns:
        if col not in frame.columns:
            continue
        for smiles in frame[col].dropna().astype(str):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                fps.append(FP_GENERATOR.GetFingerprint(mol))
        if fps:
            break
    return fps


def max_similarity(smiles: str, fingerprints: list) -> float | str:
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None or not fingerprints:
        return ""
    fp = FP_GENERATOR.GetFingerprint(mol)
    return round(max(DataStructs.BulkTanimotoSimilarity(fp, fingerprints)), 8)


def assay_inhibitor(assay: dict) -> str:
    for parameter in assay.get("assay_parameters") or []:
        if str(parameter.get("standard_type", "")).lower() == "inhibitor name":
            return str(parameter.get("standard_text_value") or parameter.get("text_value") or "")
    return ""


def classify_observation(row: dict) -> tuple[str, str, str]:
    er = number(row.get("efflux_ratio"))
    er_inhibitor = number(row.get("efflux_ratio_with_inhibitor"))
    rescue = number(row.get("inhibitor_rescue_fraction"))
    has_bidir = number(row.get("papp_ab")) is not None and number(row.get("papp_ba")) is not None
    if er is not None and er_inhibitor is not None:
        if er >= 2 and rescue is not None and rescue >= 0.5:
            return (
                "positive_like_ER_ge_2_and_rescue_ge_50pct",
                "1",
                "B1_direct_MDR1_MDCK_bidirectional_paired_inhibitor_no_parent",
            )
        if er < 2:
            return (
                "negative_like_ER_lt_2_in_MDR1_MDCK",
                "0",
                "B1_direct_MDR1_MDCK_bidirectional_paired_inhibitor_no_parent",
            )
        return (
            "high_ER_not_rescued_ge_50pct_or_ambiguous",
            "",
            "B1_direct_MDR1_MDCK_bidirectional_paired_inhibitor_no_parent",
        )
    if has_bidir and er is not None:
        return (
            "high_efflux_ER_ge_2_no_inhibitor_confirmation" if er >= 2 else "low_efflux_ER_lt_2_no_inhibitor_confirmation",
            "1" if er >= 2 else "0",
            "B2_direct_MDR1_MDCK_bidirectional_no_parent_no_inhibitor",
        )
    if has_bidir:
        return (
            "bidirectional_ER_nonestimable_zero_or_censored_papp_ab",
            "",
            "B2_direct_MDR1_MDCK_bidirectional_no_parent_no_inhibitor",
        )
    return (
        "unidirectional_permeability_only_no_substrate_call",
        "",
        "C1_MDR1_MDCK_unidirectional_Papp_and_recovery_only",
    )


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    all_activities: list[dict] = []
    assay_lookup: dict[str, dict] = {}
    source_rows: list[dict] = []

    for document_id, metadata in DOCUMENTS.items():
        document = fetch_json(f"https://www.ebi.ac.uk/chembl/api/data/document/{document_id}.json")
        assays = fetch_json(
            "https://www.ebi.ac.uk/chembl/api/data/assay.json",
            {"document_chembl_id": document_id, "limit": 1000},
        )
        activities = fetch_json(
            "https://www.ebi.ac.uk/chembl/api/data/activity.json",
            {"document_chembl_id": document_id, "limit": 1000},
        )
        document_path = RAW_DIR / f"{document_id}_document.json"
        assay_path = RAW_DIR / f"{document_id}_assays.json"
        activity_path = RAW_DIR / f"{document_id}_activities.json"
        write_json(document_path, document)
        write_json(assay_path, assays)
        write_json(activity_path, activities)
        for assay in assays.get("assays", []):
            assay_lookup[assay["assay_chembl_id"]] = assay
        all_activities.extend(activities.get("activities", []))
        source_rows.append(
            {
                "document_chembl_id": document_id,
                "title": document.get("title", ""),
                "year": document.get("year", metadata["year"]),
                "doi": metadata["doi"],
                "document_url": f"https://www.ebi.ac.uk/chembl/explore/document/{document_id}",
                "document_api_url": f"https://www.ebi.ac.uk/chembl/api/data/document/{document_id}.json",
                "assay_api_url": f"https://www.ebi.ac.uk/chembl/api/data/assay.json?document_chembl_id={document_id}&limit=1000",
                "activity_api_url": f"https://www.ebi.ac.uk/chembl/api/data/activity.json?document_chembl_id={document_id}&limit=1000",
                "chembl_release": metadata["release"],
                "activity_rows": len(activities.get("activities", [])),
                "unique_molecules": len(
                    {row.get("molecule_chembl_id") for row in activities.get("activities", [])}
                ),
                "assay_count": len(assays.get("assays", [])),
                "predeclared_role": metadata["predeclared_role"],
                "document_raw_sha256": sha256(document_path),
                "assay_raw_sha256": sha256(assay_path),
                "activity_raw_sha256": sha256(activity_path),
                "allowed_use": "continuous endpoint modeling, weak supervision, source-heldout sensitivity analysis, and assay-QC research",
                "prohibited_use": "not strict matched-parent ABCB1 truth; not intact-ADC or release-species validation",
            }
        )

    activity_rows = [{field: row.get(field, "") for field in ACTIVITY_FIELDS} for row in all_activities]
    activity_output = PROCESSED / "abcb1_asap_mdck_mdr1_activity_rows_v1.csv"
    write_csv(activity_output, activity_rows)

    frame = pd.DataFrame(activity_rows)
    frame["activity_id_numeric"] = pd.to_numeric(frame["activity_id"], errors="coerce")
    frame = frame.sort_values(
        ["document_chembl_id", "assay_chembl_id", "molecule_chembl_id", "record_id", "activity_id_numeric"]
    )
    frame["endpoint_occurrence"] = frame.groupby(
        ["document_chembl_id", "assay_chembl_id", "molecule_chembl_id", "record_id", "standard_type"],
        dropna=False,
    ).cumcount()

    prior_indexes = {
        name: identity_index(path, columns)
        for name, (path, columns) in PRIOR_DATASETS.items()
        if path.exists()
    }
    training_fps = load_prior_fingerprints(
        PROCESSED / "abcb1_screener_v0_1_training_audit.csv",
        ["standardized_isomeric_smiles", "derived_isomeric_smiles"],
    )
    direct_fps = load_prior_fingerprints(
        PROCESSED / "soskuti_2024_mdr1_species_panel_v1.csv",
        ["standardized_isomeric_smiles"],
    )

    observations: list[dict] = []
    group_keys = [
        "document_chembl_id",
        "assay_chembl_id",
        "molecule_chembl_id",
        "record_id",
        "endpoint_occurrence",
    ]
    for key, group in frame.groupby(group_keys, dropna=False):
        document_id, assay_id, molecule_id, record_id, occurrence = key
        source = group.iloc[0]
        structure = standardize_smiles(str(source.get("canonical_smiles") or ""))
        inchikey = structure.get("standardized_inchikey", "")
        connectivity = inchikey.split("-")[0] if inchikey else ""
        exact_hits = [name for name, (exact, _) in prior_indexes.items() if inchikey and inchikey in exact]
        connectivity_hits = [
            name for name, (_, blocks) in prior_indexes.items() if connectivity and connectivity in blocks
        ]
        identity_status = (
            "exclude_unresolved_structure"
            if structure.get("standardization_status") != "ok"
            else "prior_exact_identity_overlap"
            if exact_hits
            else "prior_connectivity_identity_overlap"
            if connectivity_hits
            else "no_exact_or_connectivity_prior_overlap"
        )
        endpoint_values: dict[str, object] = {}
        endpoint_activity_ids: list[str] = []
        endpoint_relations: list[str] = []
        source_units: list[str] = []
        for source_endpoint, output_endpoint in ENDPOINT_MAP.items():
            relation, value, units, activity_id = relation_value(group, source_endpoint)
            endpoint_values[output_endpoint] = value if value is not None else ""
            endpoint_values[f"{output_endpoint}_relation"] = relation
            endpoint_values[f"{output_endpoint}_source_units"] = units
            if activity_id:
                endpoint_activity_ids.append(activity_id)
            if relation:
                endpoint_relations.append(f"{output_endpoint}:{relation}")
            if units:
                source_units.append(f"{output_endpoint}:{units}")

        papp_ab = number(endpoint_values.get("papp_ab"))
        papp_ba = number(endpoint_values.get("papp_ba"))
        er_source = number(endpoint_values.get("efflux_ratio_source"))
        er = er_source if er_source is not None else (papp_ba / papp_ab if papp_ab and papp_ba is not None else None)
        papp_ab_inh = number(endpoint_values.get("papp_ab_inhibitor"))
        papp_ba_inh = number(endpoint_values.get("papp_ba_inhibitor"))
        er_inh_source = number(endpoint_values.get("efflux_ratio_inhibitor_source"))
        er_inhibitor = (
            er_inh_source
            if er_inh_source is not None
            else (papp_ba_inh / papp_ab_inh if papp_ab_inh and papp_ba_inh is not None else None)
        )
        rescue = 1 - er_inhibitor / er if er and er_inhibitor is not None else None
        assay = assay_lookup[str(assay_id)]
        row = {
            "observation_id": f"ASAP_MDR1::{document_id}::{assay_id}::{molecule_id}::{record_id}::R{int(occurrence)+1}",
            "document_chembl_id": document_id,
            "document_year": DOCUMENTS[str(document_id)]["year"],
            "document_doi": DOCUMENTS[str(document_id)]["doi"],
            "assay_chembl_id": assay_id,
            "assay_group": assay.get("assay_group", ""),
            "assay_cell_type": assay.get("assay_cell_type", ""),
            "assay_description": assay.get("description", ""),
            "inhibitor": assay_inhibitor(assay),
            "molecule_chembl_id": molecule_id,
            "parent_molecule_chembl_id": source.get("parent_molecule_chembl_id", ""),
            "molecule_pref_name": source.get("molecule_pref_name", ""),
            "record_id": record_id,
            "replicate_index_within_record": int(occurrence) + 1,
            "source_canonical_smiles": source.get("canonical_smiles", ""),
            **structure,
            "inchikey_connectivity_block": connectivity,
            **endpoint_values,
            "efflux_ratio": round(er, 8) if er is not None else "",
            "efflux_ratio_with_inhibitor": round(er_inhibitor, 8) if er_inhibitor is not None else "",
            "inhibitor_rescue_fraction": round(rescue, 8) if rescue is not None else "",
            "endpoint_activity_ids": " | ".join(endpoint_activity_ids),
            "endpoint_relations": " | ".join(endpoint_relations),
            "source_units": " | ".join(source_units),
            "identity_status": identity_status,
            "overlap_exact_prior_sets": " | ".join(exact_hits),
            "overlap_connectivity_prior_sets": " | ".join(connectivity_hits),
            "max_morgan_similarity_to_frozen_443": max_similarity(
                structure.get("standardized_isomeric_smiles", ""), training_fps
            ),
            "max_morgan_similarity_to_soskuti_direct_panel": max_similarity(
                structure.get("standardized_isomeric_smiles", ""), direct_fps
            ),
            "upstream_experiment_cluster": f"ASAP::{document_id}::{assay_id}",
            "source_url": f"https://www.ebi.ac.uk/chembl/explore/assay/{assay_id}",
            "predeclared_source_role": DOCUMENTS[str(document_id)]["predeclared_role"],
            "strict_binary_truth_eligible": "no",
            "model_scoring_performed": "no",
        }
        reconstructed_class, weak_label, tier = classify_observation(row)
        row.update(
            {
                "reconstructed_class": reconstructed_class,
                "weak_binary_label": weak_label,
                "evidence_tier": tier,
                "development_candidate_eligible": (
                    "yes_weak_or_continuous_only"
                    if identity_status == "no_exact_or_connectivity_prior_overlap"
                    and str(document_id) == "CHEMBL5442303"
                    and reconstructed_class != "unidirectional_permeability_only_no_substrate_call"
                    else "reserved_source_heldout_continuous_only"
                    if identity_status == "no_exact_or_connectivity_prior_overlap"
                    and str(document_id) == "CHEMBL5674852"
                    and reconstructed_class != "unidirectional_permeability_only_no_substrate_call"
                    else "permeability_only_or_prior_overlap"
                ),
                "allowed_use": "assay-aware continuous modeling, weak supervision, source-grouped sensitivity, and historical source-heldout evaluation",
                "prohibited_use": "not strict matched-parent ABCB1 truth; not universal binary substrate label; not ADC validation",
            }
        )
        observations.append(row)

    collision_map: dict[tuple[str, str], set[str]] = {}
    for row in observations:
        collision_map.setdefault(
            (row["document_chembl_id"], row["standardized_inchikey"]), set()
        ).add(row["molecule_chembl_id"])
    for row in observations:
        source_ids = collision_map[(row["document_chembl_id"], row["standardized_inchikey"])]
        row["standardization_collision_source_molecule_count"] = len(source_ids)
        row["standardization_collision_source_molecule_ids"] = " | ".join(sorted(source_ids))
        row["standardization_collision_flag"] = "yes" if len(source_ids) > 1 else "no"

    observations.sort(
        key=lambda row: (
            row["document_chembl_id"],
            row["assay_chembl_id"],
            row["molecule_chembl_id"],
            int(row["record_id"]),
            row["replicate_index_within_record"],
        )
    )
    observation_output = PROCESSED / "abcb1_asap_mdck_mdr1_observations_v1.csv"
    source_output = PROCESSED / "abcb1_asap_mdck_mdr1_source_registry_v1.csv"
    write_csv(observation_output, observations)
    write_csv(source_output, source_rows)

    identity_clean = [row for row in observations if row["identity_status"] == "no_exact_or_connectivity_prior_overlap"]
    bidirectional_measured = [
        row
        for row in observations
        if number(row.get("papp_ab")) is not None and number(row.get("papp_ba")) is not None
    ]
    estimable_er = [row for row in observations if number(row.get("efflux_ratio")) is not None]
    paired_inhibitor = [
        row for row in observations if number(row.get("efflux_ratio_with_inhibitor")) is not None
    ]
    summary = {
        "created_by": Path(__file__).name,
        "source_documents": len(DOCUMENTS),
        "activity_rows": len(activity_rows),
        "observation_rows": len(observations),
        "unique_molecules": len({row["standardized_inchikey"] for row in observations if row["standardized_inchikey"]}),
        "source_molecule_chembl_ids": len({row["molecule_chembl_id"] for row in observations}),
        "document_source_molecule_chembl_ids": {
            document_id: len(
                {
                    row["molecule_chembl_id"]
                    for row in observations
                    if row["document_chembl_id"] == document_id
                }
            )
            for document_id in DOCUMENTS
        },
        "document_unique_molecules": {
            document_id: len(
                {
                    row["standardized_inchikey"]
                    for row in observations
                    if row["document_chembl_id"] == document_id and row["standardized_inchikey"]
                }
            )
            for document_id in DOCUMENTS
        },
        "bidirectional_measured_observation_rows": len(bidirectional_measured),
        "estimable_efflux_ratio_observation_rows": len(estimable_er),
        "bidirectional_nonestimable_ER_rows": len(bidirectional_measured) - len(estimable_er),
        "paired_inhibitor_observation_rows": len(paired_inhibitor),
        "paired_inhibitor_unique_molecules": len(
            {row["standardized_inchikey"] for row in paired_inhibitor if row["standardized_inchikey"]}
        ),
        "identity_status_counts": dict(Counter(row["identity_status"] for row in observations)),
        "identity_clean_unique_molecules": len(
            {row["standardized_inchikey"] for row in identity_clean if row["standardized_inchikey"]}
        ),
        "standardization_collision_groups": sum(len(source_ids) > 1 for source_ids in collision_map.values()),
        "standardization_collision_source_molecule_ids": sum(
            len(source_ids) for source_ids in collision_map.values() if len(source_ids) > 1
        ),
        "reconstructed_class_counts": dict(Counter(row["reconstructed_class"] for row in observations)),
        "evidence_tier_counts": dict(Counter(row["evidence_tier"] for row in observations)),
        "development_role_counts": dict(
            Counter(row["development_candidate_eligible"] for row in observations)
        ),
        "strict_binary_truth_rows_created": 0,
        "model_scoring_performed": False,
        "formal_status": "open ASAP MDCK-MDR1 endpoints collected and evidence-tiered; 2025 source remains unscored",
        "output_hashes": {
            str(activity_output.relative_to(PROJECT)): sha256(activity_output),
            str(observation_output.relative_to(PROJECT)): sha256(observation_output),
            str(source_output.relative_to(PROJECT)): sha256(source_output),
        },
    }
    summary_output = PROCESSED / "abcb1_asap_mdck_mdr1_summary_v1.json"
    write_json(summary_output, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

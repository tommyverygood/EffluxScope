#!/usr/bin/env python3
"""Build an observation-level MDR1 panel from Soskuti et al. 2024 Table S1.

The source reports efflux-ratio mean +/- SD values for 83 analysis compounds
plus a separately marked famotidine permeability control in a
matched Mock Abcb1KO-MDCKII line, four MDR1 ortholog lines, and hBCRP/MDCKII.
This script deliberately creates *screening-evidence* labels only.  The paper
does not report compound-level inhibitor rescue in Table S1, so a net efflux
ratio >= 2 is not promoted to a universally confirmed substrate label.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from abcb1_screener_core_v0_1 import standardize_smiles  # noqa: E402


SOURCE_ID = "SRC_SOSKUTI_2024_MDR1_SPECIES_TABLE_S1"
SOURCE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11207571/"
DOI = "10.3390/pharmaceutics16060736"
PMCID = "PMC11207571"
REGULATORY_THRESHOLD = 2.0
PAIR_RE = re.compile(r"(?P<mean>\d+(?:\.\d+)?)\s*±\s*(?P<sd>\d+(?:\.\d+)?)")

SYSTEMS = ["Mock", "hMDR1", "rMDR1", "mMDR1", "cyMDR1", "BCRP"]
SYSTEM_META = {
    "Mock": ("parental_control", "", "canine", "Abcb1KO-MDCKII"),
    "hMDR1": ("transporter_expressing", "ABCB1", "human", "Abcb1KO-MDCKII"),
    "rMDR1": ("transporter_expressing", "Mdr1 ortholog", "rat", "Abcb1KO-MDCKII"),
    "mMDR1": ("transporter_expressing", "Mdr1 ortholog", "mouse", "Abcb1KO-MDCKII"),
    "cyMDR1": ("transporter_expressing", "MDR1 ortholog", "cynomolgus monkey", "Abcb1KO-MDCKII"),
    "BCRP": ("transporter_expressing", "ABCG2", "human", "MDCKII"),
}

NAME_ALIASES = {
    "Amitriptylline": "Amitriptyline",
}


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
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def extract_pdf_text(pdf_path: Path) -> str:
    completed = subprocess.run(
        ["pdftotext", "-layout", str(pdf_path), "-"],
        check=True,
        text=True,
        capture_output=True,
    )
    return completed.stdout


def parse_table_s1(text: str) -> list[dict]:
    start_marker = "Table S1. Efflux Ratios across the six cell lines for 83 compounds"
    end_marker = "*Excluded from data analysis due to low recovery"
    start = text.find(start_marker)
    end = text.find(end_marker, start)
    if start < 0 or end < 0:
        raise ValueError("Table S1 boundaries were not found in extracted PDF text")
    block = text[start:end]
    rows: list[dict] = []
    for line in block.splitlines():
        pairs = list(PAIR_RE.finditer(line))
        if len(pairs) != 6:
            continue
        raw_name = line[: pairs[0].start()].replace("\f", " ").strip()
        if not raw_name or raw_name.lower().startswith("compounds"):
            continue
        low_recovery = "*" in raw_name
        permeability_control = "#" in raw_name
        name_source = re.sub(r"\s*[*#]\s*", "", raw_name).strip()
        name_normalized = NAME_ALIASES.get(name_source, name_source)
        values = {}
        for system, match in zip(SYSTEMS, pairs):
            values[f"{system}_efflux_ratio_mean"] = float(match.group("mean"))
            values[f"{system}_efflux_ratio_sd"] = float(match.group("sd"))
        rows.append(
            {
                "entity_name_source": name_source,
                "entity_name_normalized": name_normalized,
                "low_recovery_flag": "yes" if low_recovery else "no",
                "permeability_control_flag": "yes" if permeability_control else "no",
                **values,
            }
        )
    # The title says 83 compounds, while the body contains those 83 analysis
    # compounds plus famotidine marked with # as a permeability control.
    if len(rows) != 84:
        raise ValueError(f"expected 84 Table S1 body rows (83 compounds + 1 control), parsed {len(rows)}")
    if len({row['entity_name_source'] for row in rows}) != 84:
        raise ValueError("Table S1 entity names are not unique")
    return rows


def safe_cache_name(name: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "compound"
    suffix = hashlib.sha1(name.encode("utf-8")).hexdigest()[:10]
    return f"{stem}_{suffix}.json"


def fetch_pubchem(name: str, cache_dir: Path, delay: float) -> tuple[dict, str, str]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / safe_cache_name(name)
    if cache_path.exists():
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
        status = "cached_pubchem_name_query"
    else:
        encoded = urllib.parse.quote(name, safe="")
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{encoded}/property/CanonicalSMILES,IsomericSMILES,InChIKey,"
            "MolecularFormula,MolecularWeight,XLogP,TPSA/JSON"
        )
        request = urllib.request.Request(url, headers={"User-Agent": "ABCB1-research-audit/1.0"})
        last_error = ""
        payload = {}
        for attempt in range(4):
            try:
                with urllib.request.urlopen(request, timeout=40) as response:
                    payload = json.loads(response.read().decode("utf-8"))
                cache_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
                status = "live_pubchem_name_query"
                time.sleep(delay)
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
    smiles = prop.get("SMILES") or prop.get("IsomericSMILES") or prop.get("ConnectivitySMILES") or prop.get("CanonicalSMILES") or ""
    record = {
        "pubchem_cid": prop.get("CID", ""),
        "pubchem_inchikey": prop.get("InChIKey", ""),
        "pubchem_isomeric_smiles": smiles,
        "pubchem_connectivity_smiles": prop.get("ConnectivitySMILES") or prop.get("CanonicalSMILES") or "",
        "pubchem_molecular_formula": prop.get("MolecularFormula", ""),
        "pubchem_molecular_weight": prop.get("MolecularWeight", ""),
        "pubchem_xlogp": prop.get("XLogP", ""),
        "pubchem_tpsa": prop.get("TPSA", ""),
        "pubchem_cache_file": str(cache_path.relative_to(PROJECT)),
    }
    return record, status, ""


def approximate_ratio_interval(numerator: float, numerator_sd: float, denominator: float, denominator_sd: float) -> tuple[float, float, float, float]:
    ratio = numerator / denominator
    # First-order independent-error propagation. This is a sensitivity flag,
    # not a replacement for replicate-level statistical inference.
    relative_variance = (numerator_sd / numerator) ** 2 + (denominator_sd / denominator) ** 2
    ratio_sd = ratio * math.sqrt(relative_variance)
    lower = max(0.0, ratio - 1.96 * ratio_sd)
    upper = ratio + 1.96 * ratio_sd
    return ratio, ratio_sd, lower, upper


def build_identity_indexes(training_rows: list[dict], external_rows: list[dict]) -> tuple[dict, dict, set[str], set[str]]:
    exact_labels: dict[str, set[str]] = defaultdict(set)
    connectivity_labels: dict[str, set[str]] = defaultdict(set)
    for row in training_rows:
        # Only the 443 rows actually used by the frozen supervised model are
        # admissible here.  Blank PU labels are not a third class and must
        # never be interpreted as discordance.
        if row.get("target_label", "") not in {"0", "1"}:
            continue
        key = row.get("standardized_inchikey", "")
        block = key.split("-")[0] if key else ""
        label = row.get("target_label", "")
        if key:
            exact_labels[key].add(label)
        if block:
            connectivity_labels[block].add(label)
    external_exact: set[str] = set()
    external_connectivity: set[str] = set()
    for row in external_rows:
        key = row.get("model_standardized_inchikey") or row.get("standardized_inchikey", "")
        block = row.get("model_standardized_connectivity_block") or (key.split("-")[0] if key else "")
        if key:
            external_exact.add(key)
        if block:
            external_connectivity.add(block)
    return exact_labels, connectivity_labels, external_exact, external_connectivity


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdf",
        type=Path,
        default=PROJECT / "数据/raw/pharmaceutics_2024_mdr1_species/s001/pharmaceutics-3015208-supplementary.pdf",
    )
    parser.add_argument("--pubchem-delay", type=float, default=0.12)
    args = parser.parse_args()

    processed = PROJECT / "数据/processed"
    raw_dir = PROJECT / "数据/raw/pharmaceutics_2024_mdr1_species"
    cache_dir = raw_dir / "pubchem"
    training_path = processed / "abcb1_screener_v0_1_training_audit.csv"
    external_path = PROJECT / "outputs/20260803_abcb1_q1_external_benchmark_v1/ABCB1_v0.2_triage_external_overlay_diagnostic.csv"
    training_rows = read_csv(training_path)
    external_rows = read_csv(external_path)
    exact_labels, connectivity_labels, external_exact, external_connectivity = build_identity_indexes(training_rows, external_rows)

    pdf_text = extract_pdf_text(args.pdf)
    parsed = parse_table_s1(pdf_text)
    panel_rows: list[dict] = []
    long_rows: list[dict] = []
    structure_failures: list[dict] = []

    for index, row in enumerate(parsed, start=1):
        pubchem, query_status, query_error = fetch_pubchem(row["entity_name_normalized"], cache_dir, args.pubchem_delay)
        standardized = standardize_smiles(pubchem.get("pubchem_isomeric_smiles", ""))
        key = standardized.get("standardized_inchikey", "")
        block = key.split("-")[0] if key else ""
        exact_overlap = key in exact_labels if key else False
        connectivity_overlap = block in connectivity_labels if block else False
        external_exact_overlap = key in external_exact if key else False
        external_connectivity_overlap = block in external_connectivity if block else False
        weak_labels = sorted(exact_labels.get(key, set()) or connectivity_labels.get(block, set()))

        net, net_sd, net_low, net_high = approximate_ratio_interval(
            row["hMDR1_efflux_ratio_mean"],
            row["hMDR1_efflux_ratio_sd"],
            row["Mock_efflux_ratio_mean"],
            row["Mock_efflux_ratio_sd"],
        )
        low_recovery = row["low_recovery_flag"] == "yes"
        control = row["permeability_control_flag"] == "yes"
        if low_recovery:
            derived_class = "excluded_low_recovery_le_60pct"
            derived_binary: int | str = ""
            admission = "exclude_source_qc"
        elif control:
            derived_class = "excluded_permeability_control"
            derived_binary = ""
            admission = "exclude_source_control"
        elif net >= REGULATORY_THRESHOLD:
            derived_class = "positive_like_net_ER_ge_2_no_inhibitor_confirmation"
            derived_binary = 1
            admission = "development_screen_label"
        else:
            derived_class = "negative_like_net_ER_lt_2_no_inhibitor_confirmation"
            derived_binary = 0
            admission = "development_screen_label"

        if low_recovery or control:
            certainty = "not_applicable_excluded"
            conservative_binary: int | str = ""
        elif net_low >= REGULATORY_THRESHOLD:
            certainty = "interval_entirely_positive_like"
            conservative_binary = 1
        elif net_high < REGULATORY_THRESHOLD:
            certainty = "interval_entirely_negative_like"
            conservative_binary = 0
        else:
            certainty = "borderline_interval_crosses_threshold"
            conservative_binary = ""

        weak_conflict = "not_comparable_no_weak_identity_overlap"
        if weak_labels and derived_binary != "":
            weak_conflict = "concordant" if str(derived_binary) in weak_labels else "discordant"

        identity_overlap = (
            "exact_model_standardized_identity_overlap"
            if exact_overlap
            else "connectivity_only_model_standardized_overlap"
            if connectivity_overlap
            else "no_exact_or_connectivity_overlap"
        )
        development_external_overlap = (
            "exact_model_standardized_identity_overlap"
            if external_exact_overlap
            else "connectivity_only_model_standardized_overlap"
            if external_connectivity_overlap
            else "no_exact_or_connectivity_overlap"
        )
        base = {
            "panel_entity_id": f"SOSKUTI2024_{index:03d}",
            **row,
            **pubchem,
            "structure_resolution_status": query_status,
            "structure_resolution_error": query_error,
            **standardized,
            "model_standardized_connectivity_block": block,
            "overlap_with_frozen_training": identity_overlap,
            "overlap_with_current_development_external": development_external_overlap,
            "weak_training_labels_at_identity": "|".join(weak_labels),
            "weak_vs_direct_screen_label_status": weak_conflict,
            "human_mdr1_net_efflux_ratio_vs_mock": round(net, 6),
            "human_mdr1_net_efflux_ratio_approx_sd": round(net_sd, 6),
            "human_mdr1_net_efflux_ratio_approx_95ci_low": round(net_low, 6),
            "human_mdr1_net_efflux_ratio_approx_95ci_high": round(net_high, 6),
            "regulatory_screen_threshold": REGULATORY_THRESHOLD,
            "derived_human_mdr1_screen_class": derived_class,
            "derived_human_mdr1_screen_binary": derived_binary,
            "threshold_certainty_status": certainty,
            "conservative_human_mdr1_screen_binary": conservative_binary,
            "training_admission_role": admission,
            "inhibitor_confirmation_available": "no_table_S1_compound_level_rescue",
            "evidence_tier": "A_direct_matched_parental_numeric_no_inhibitor",
            "source_id": SOURCE_ID,
            "source_doi": DOI,
            "source_pmcid": PMCID,
            "source_url": SOURCE_URL,
            "source_table": "Supplementary Table S1",
            "source_pdf_pages": "3-4",
            "concentration_value": 1,
            "concentration_unit": "uM",
            "duration_value": 120,
            "duration_unit": "min",
            "claim_boundary": "matched-parental regulatory-style screening evidence; no compound-level inhibitor confirmation; not universal substrate truth",
        }
        panel_rows.append(base)
        if standardized.get("standardization_status") != "ok":
            structure_failures.append(
                {
                    "panel_entity_id": base["panel_entity_id"],
                    "entity_name_source": row["entity_name_source"],
                    "entity_name_normalized": row["entity_name_normalized"],
                    "structure_resolution_status": query_status,
                    "structure_resolution_error": query_error,
                    "standardization_error": standardized.get("standardization_error", ""),
                }
            )

        for system in SYSTEMS:
            role, gene, species, background = SYSTEM_META[system]
            er_mean = row[f"{system}_efflux_ratio_mean"]
            er_sd = row[f"{system}_efflux_ratio_sd"]
            if system == "Mock":
                ratio, ratio_sd, ratio_low, ratio_high = 1.0, "", "", ""
            else:
                ratio, ratio_sd, ratio_low, ratio_high = approximate_ratio_interval(
                    er_mean, er_sd, row["Mock_efflux_ratio_mean"], row["Mock_efflux_ratio_sd"]
                )
            long_rows.append(
                {
                    "observation_id": f"SOSKUTI2024_{index:03d}_{system}",
                    "panel_entity_id": base["panel_entity_id"],
                    "entity_name_source": row["entity_name_source"],
                    "entity_name_normalized": row["entity_name_normalized"],
                    "cell_line_code": system,
                    "cell_line_role": role,
                    "transporter_gene": gene,
                    "transporter_species": species,
                    "cell_background": background,
                    "efflux_ratio_mean": er_mean,
                    "efflux_ratio_sd": er_sd,
                    "parent_mock_efflux_ratio_mean": row["Mock_efflux_ratio_mean"],
                    "net_efflux_ratio_vs_mock": round(float(ratio), 6) if ratio != "" else "",
                    "net_efflux_ratio_approx_sd": round(float(ratio_sd), 6) if ratio_sd != "" else "",
                    "net_efflux_ratio_approx_95ci_low": round(float(ratio_low), 6) if ratio_low != "" else "",
                    "net_efflux_ratio_approx_95ci_high": round(float(ratio_high), 6) if ratio_high != "" else "",
                    "derived_human_mdr1_screen_class": derived_class if system == "hMDR1" else "not_target_human_mdr1_row",
                    "derived_human_mdr1_screen_binary": derived_binary if system == "hMDR1" else "",
                    "threshold_certainty_status": certainty if system == "hMDR1" else "not_target_human_mdr1_row",
                    "low_recovery_flag": row["low_recovery_flag"],
                    "permeability_control_flag": row["permeability_control_flag"],
                    "concentration_value": 1,
                    "concentration_unit": "uM",
                    "duration_value": 120,
                    "duration_unit": "min",
                    "standardized_inchikey": key,
                    "standardized_isomeric_smiles": standardized.get("standardized_isomeric_smiles", ""),
                    "source_id": SOURCE_ID,
                    "source_doi": DOI,
                    "source_pmcid": PMCID,
                    "source_url": SOURCE_URL,
                    "source_table": "Supplementary Table S1",
                    "source_pdf_pages": "3-4",
                    "evidence_tier": "A_direct_matched_parental_numeric_no_inhibitor",
                    "claim_boundary": "numeric observation; derived net ER uses matched Mock mean; approximate interval assumes independent mean errors",
                }
            )

    panel_path = processed / "soskuti_2024_mdr1_species_panel_v1.csv"
    observations_path = processed / "soskuti_2024_mdr1_species_observations_v1.csv"
    failures_path = processed / "soskuti_2024_mdr1_species_structure_resolution_issues_v1.csv"
    validation_path = processed / "soskuti_2024_mdr1_species_validation_v1.json"
    source_registry_path = processed / "soskuti_2024_mdr1_species_source_registry_v1.csv"
    write_csv(panel_path, panel_rows)
    write_csv(observations_path, long_rows)
    write_csv(failures_path, structure_failures)
    write_csv(
        source_registry_path,
        [
            {
                "source_id": SOURCE_ID,
                "citation": "Soskuti et al. Pharmaceutics 2024;16:736",
                "doi": DOI,
                "pmcid": PMCID,
                "source_url": SOURCE_URL,
                "local_source_pdf": str(args.pdf.relative_to(PROJECT)),
                "local_source_pdf_sha256": sha256(args.pdf),
                "source_location": "Supplementary Table S1, PDF pages 3-4",
                "extraction_method": "pdftotext -layout plus deterministic six-pair row parser",
                "data_scope": "83 analysis compounds plus 1 permeability control x 6 cell lines; ER mean and SD; 1 uM; 120 min",
                "label_derivation": "human hMDR1 ER / matched Mock ER; point threshold 2; approximate propagated interval is a sensitivity flag",
                "claim_boundary": "development screening evidence; low-recovery and permeability-control rows excluded; inhibitor rescue absent",
            }
        ],
    )

    eligible = [row for row in panel_rows if row["derived_human_mdr1_screen_binary"] != ""]
    conservative = [row for row in panel_rows if row["conservative_human_mdr1_screen_binary"] != ""]
    validation = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_pdf_sha256": sha256(args.pdf),
        "frozen_443_training_audit_sha256": sha256(training_path),
        "development_external_overlay_sha256": sha256(external_path),
        "parsed_entity_count": len(panel_rows),
        "long_observation_count": len(long_rows),
        "expected_long_observation_count": 84 * 6,
        "low_recovery_entity_count": sum(row["low_recovery_flag"] == "yes" for row in panel_rows),
        "permeability_control_entity_count": sum(row["permeability_control_flag"] == "yes" for row in panel_rows),
        "point_label_eligible_entity_count": len(eligible),
        "point_label_counts": Counter(str(row["derived_human_mdr1_screen_binary"]) for row in eligible),
        "conservative_interval_label_entity_count": len(conservative),
        "conservative_label_counts": Counter(str(row["conservative_human_mdr1_screen_binary"]) for row in conservative),
        "borderline_interval_entity_count": sum(row["threshold_certainty_status"] == "borderline_interval_crosses_threshold" for row in panel_rows),
        "structure_standardization_ok_count": sum(row["standardization_status"] == "ok" for row in panel_rows),
        "structure_resolution_issue_count": len(structure_failures),
        "frozen_training_exact_overlap_count": sum(row["overlap_with_frozen_training"] == "exact_model_standardized_identity_overlap" for row in panel_rows),
        "frozen_training_connectivity_only_overlap_count": sum(row["overlap_with_frozen_training"] == "connectivity_only_model_standardized_overlap" for row in panel_rows),
        "current_development_external_exact_overlap_count": sum(row["overlap_with_current_development_external"] == "exact_model_standardized_identity_overlap" for row in panel_rows),
        "current_development_external_connectivity_only_overlap_count": sum(row["overlap_with_current_development_external"] == "connectivity_only_model_standardized_overlap" for row in panel_rows),
        "weak_vs_direct_screen_status_counts": Counter(row["weak_vs_direct_screen_label_status"] for row in panel_rows),
        "checks": {
            "parsed_83_compounds_plus_1_control": len(panel_rows) == 84,
            "generated_504_long_observations": len(long_rows) == 504,
            "unique_panel_entity_ids": len({row["panel_entity_id"] for row in panel_rows}) == 84,
            "all_six_systems_per_entity": all(count == 6 for count in Counter(row["panel_entity_id"] for row in long_rows).values()),
            "all_numeric_er_positive": all(float(row["efflux_ratio_mean"]) > 0 for row in long_rows),
            "no_low_recovery_or_control_admitted": all(row["low_recovery_flag"] == "no" and row["permeability_control_flag"] == "no" for row in eligible),
        },
        "limitations": [
            "PubChem name-to-structure matches are automated and provisional until blinded human identity review.",
            "The net-ER approximate interval assumes independent errors and is only a sensitivity screen.",
            "Table S1 has no compound-level selective-inhibitor rescue; derived calls are not confirmed universal substrate labels.",
            "This 2024 panel is assigned to development and must never be used as the future untouched temporal/prospective validation set.",
        ],
        "outputs": {
            "panel": str(panel_path.relative_to(PROJECT)),
            "observations": str(observations_path.relative_to(PROJECT)),
            "structure_issues": str(failures_path.relative_to(PROJECT)),
            "source_registry": str(source_registry_path.relative_to(PROJECT)),
        },
    }
    validation_path.write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

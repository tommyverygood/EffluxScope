#!/usr/bin/env python3
"""Validate the ABCB1 external evidence-channel ablation package."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"
PROTOCOL = PROCESSED / "abcb1_external_evidence_channel_ablation_protocol_v1.json"
SUMMARY = PROCESSED / "abcb1_external_evidence_channel_ablation_summary_v1.json"
VALIDATION = PROCESSED / "abcb1_external_evidence_channel_ablation_validation_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text())
    summary = json.loads(SUMMARY.read_text())
    metrics = pd.read_csv(PROCESSED / "abcb1_external_evidence_channel_ablation_metrics_v1.csv")
    contrasts = pd.read_csv(PROCESSED / "abcb1_external_evidence_channel_set_contrasts_v1.csv")
    pairwise = pd.read_csv(PROCESSED / "abcb1_external_evidence_channel_pairwise_bootstrap_v1.csv")
    figure_manifest = pd.read_csv(PROCESSED / "abcb1_external_evidence_channel_ablation_figure_manifest_v1.csv")

    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, evidence: object) -> None:
        checks.append({"check": name, "pass": bool(condition), "evidence": evidence})

    check("protocol_status_truthful_posthoc", protocol["status"].startswith("locked_posthoc_component_ablation"), protocol["status"])
    check("summary_completed", summary["status"] == "completed_frozen_model_posthoc_evidence_channel_ablation", summary["status"])
    check("summary_protocol_hash", summary["protocol_sha256"] == sha256(PROTOCOL), summary["protocol_sha256"])
    for label, entry in protocol["frozen_inputs"].items():
        path = PROJECT / entry["path"]
        check(f"frozen_input_exists_{label}", path.exists(), entry["path"])
        check(f"frozen_input_hash_{label}", path.exists() and sha256(path) == entry["sha256"], entry["sha256"])
    for relative, expected_hash in summary["output_hashes"].items():
        path = PROJECT / relative
        check(f"output_exists_{path.name}", path.exists(), relative)
        check(f"output_hash_{path.name}", path.exists() and sha256(path) == expected_hash, expected_hash)

    variant_ids = [variant["variant_id"] for variant in protocol["locked_variants"]]
    check("locked_variant_count", len(variant_ids) == 12 and len(set(variant_ids)) == 12, variant_ids)
    check("metrics_row_count", len(metrics) == 24, len(metrics))
    for population in ["identity_clean", "scaffold_clean"]:
        frame = metrics[metrics["population"].eq(population)]
        check(f"all_variants_present_{population}", set(frame["variant_id"]) == set(variant_ids), sorted(frame["variant_id"].tolist()))
        check(f"metric_bounds_{population}", bool(frame["coverage"].between(0, 1).all() and frame["balanced_accuracy"].between(0, 1).all() and frame["mcc"].between(-1, 1).all()), "coverage/BAcc/MCC bounded")
        check(f"bootstrap_complete_{population}", bool((frame["bootstrap_valid_replicates"] == 2000).all()), sorted(frame["bootstrap_valid_replicates"].unique().tolist()))

    reproduction = summary["full_v0_4_reconstruction"]
    check("full_reconstruction_exact", reproduction["exactly_reproduced"] and reproduction["symmetric_difference_n"] == 0, reproduction)
    check("full_reconstruction_anchor", reproduction["frozen_accepted_n"] == 453 and np.isclose(reproduction["balanced_accuracy"], 0.73510345) and np.isclose(reproduction["mcc"], 0.49857011), reproduction)

    identity = metrics[metrics["population"].eq("identity_clean")].set_index("variant_id")
    scaffold = metrics[metrics["population"].eq("scaffold_clean")].set_index("variant_id")
    check("direct_no_local_exceeds_primary_no_local", identity.loc["V05_direct_support_no_local", "balanced_accuracy"] > identity.loc["V03_primary_support_no_local", "balanced_accuracy"] and scaffold.loc["V05_direct_support_no_local", "balanced_accuracy"] > scaffold.loc["V03_primary_support_no_local", "balanced_accuracy"], identity.loc[["V03_primary_support_no_local", "V05_direct_support_no_local"], ["coverage", "balanced_accuracy", "mcc"]].to_dict(orient="index"))
    check("direct_complete_exceeds_primary_complete", identity.loc["V06_direct_channel_complete", "balanced_accuracy"] > identity.loc["V04_primary_channel_complete", "balanced_accuracy"] and scaffold.loc["V06_direct_channel_complete", "balanced_accuracy"] > scaffold.loc["V04_primary_channel_complete", "balanced_accuracy"], identity.loc[["V04_primary_channel_complete", "V06_direct_channel_complete"], ["coverage", "balanced_accuracy", "mcc"]].to_dict(orient="index"))
    check("cross_model_guard_adds_value", identity.loc["V08_dual_cross_model_no_local", "balanced_accuracy"] > identity.loc["V07_dual_support_only", "balanced_accuracy"] and scaffold.loc["V08_dual_cross_model_no_local", "balanced_accuracy"] > scaffold.loc["V07_dual_support_only", "balanced_accuracy"], identity.loc[["V07_dual_support_only", "V08_dual_cross_model_no_local"], ["coverage", "balanced_accuracy", "mcc"]].to_dict(orient="index"))
    check("direct_local_guard_harms_external_bacc", identity.loc["V10_dual_cross_direct_local", "balanced_accuracy"] < identity.loc["V08_dual_cross_model_no_local", "balanced_accuracy"] and scaffold.loc["V10_dual_cross_direct_local", "balanced_accuracy"] < scaffold.loc["V08_dual_cross_model_no_local", "balanced_accuracy"], identity.loc[["V08_dual_cross_model_no_local", "V10_dual_cross_direct_local"], ["coverage", "balanced_accuracy", "mcc"]].to_dict(orient="index"))
    primary_local_delta = identity.loc["V09_dual_cross_primary_local", "balanced_accuracy"] - identity.loc["V08_dual_cross_model_no_local", "balanced_accuracy"]
    direct_local_delta = identity.loc["V10_dual_cross_direct_local", "balanced_accuracy"] - identity.loc["V08_dual_cross_model_no_local", "balanced_accuracy"]
    check("direct_local_loss_larger_than_primary_local_loss", direct_local_delta < primary_local_delta, {"primary_local_delta": primary_local_delta, "direct_local_delta": direct_local_delta})
    check("descriptive_maximum_flagged_not_selection", summary["descriptive_maximum_not_for_model_selection"]["variant_id"] == "V05_direct_support_no_local" and "cannot be selected" in summary["descriptive_maximum_not_for_model_selection"]["warning"], summary["descriptive_maximum_not_for_model_selection"])

    check("set_contrast_row_count", len(contrasts) == 48, len(contrasts))
    check("set_contrast_regions_complete", set(contrasts["region"]) == {"intersection", "full_only", "comparator_only", "neither"}, sorted(contrasts["region"].unique().tolist()))
    check("set_jaccard_bounds", bool(contrasts["accepted_set_jaccard"].between(0, 1).all()), contrasts["accepted_set_jaccard"].describe().to_dict())
    check("pairwise_row_count", len(pairwise) == 36, len(pairwise))
    check("pairwise_bootstrap_complete", bool((pairwise["valid_bootstrap_replicates"] == 2000).all()), sorted(pairwise["valid_bootstrap_replicates"].unique().tolist()))

    def pair(comparator: str, metric: str) -> pd.Series:
        frame = pairwise[(pairwise["population"] == "identity_clean") & (pairwise["comparator_variant"] == comparator) & (pairwise["metric"] == metric)]
        assert len(frame) == 1
        return frame.iloc[0]

    full_vs_primary = pair("V04_primary_channel_complete", "balanced_accuracy")
    full_vs_direct = pair("V06_direct_channel_complete", "balanced_accuracy")
    full_vs_no_local = pair("V08_dual_cross_model_no_local", "balanced_accuracy")
    check("full_exceeds_primary_channel_ci", full_vs_primary["ci95_low"] > 0, full_vs_primary.to_dict())
    check("direct_complete_exceeds_full_ci", full_vs_direct["ci95_high"] < 0, full_vs_direct.to_dict())
    check("no_local_exceeds_full_ci", full_vs_no_local["ci95_high"] < 0, full_vs_no_local.to_dict())

    check("figure_manifest_complete", len(figure_manifest) == 3 and set(figure_manifest["format"]) == {"png", "svg", "pdf"}, figure_manifest[["format", "path"]].to_dict(orient="records"))
    for _, row in figure_manifest.iterrows():
        path = PROJECT / row["path"]
        check(f"figure_exists_{path.name}", path.exists(), row["path"])
        check(f"figure_hash_{path.name}", path.exists() and sha256(path) == row["sha256"], row["sha256"])
        check(f"figure_nontrivial_{path.name}", path.exists() and path.stat().st_size >= 20_000, int(path.stat().st_size) if path.exists() else 0)

    result = {
        "validation_version": "ABCB1_EXTERNAL_EVIDENCE_CHANNEL_ABLATION_VALIDATION_V1",
        "validation_date": "2026-08-05",
        "overall_status": "pass" if all(item["pass"] for item in checks) else "fail",
        "checks_passed": int(sum(item["pass"] for item in checks)),
        "checks_total": int(len(checks)),
        "checks": checks,
    }
    VALIDATION.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ["overall_status", "checks_passed", "checks_total"]}, ensure_ascii=False, indent=2))
    if result["overall_status"] != "pass":
        raise SystemExit(json.dumps([item for item in checks if not item["pass"]], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

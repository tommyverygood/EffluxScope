#!/usr/bin/env python3
"""Backfill row-level ChEMBL assay and primary-document provenance.

The input is the 512-row Tier A/B candidate set created in batch 10.  This
script downloads all current human ABCB1/ABCG2 target activities, intersects
them with the candidate ChEMBL molecule identifiers, resolves assay and
document records, and creates a prioritized manual-review queue.

No published binary label is changed and no API hit is automatically promoted
to Tier A.  Direct-assay candidates still require primary-full-text review of
the tested entity, controls, assay validity, and result interpretation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


DATA = Path(__file__).resolve().parents[1]
PROCESSED = DATA / "processed"
RAW = DATA / "raw" / "eleventh_batch_chembl_assay_backfill"
INPUT = PROCESSED / "transporter_truthset_AB_candidate_v1.csv"
RUN_DATE = "2026-07-28"
BASE = "https://www.ebi.ac.uk/chembl/api/data"
PAGE_LIMIT = 1000
USER_AGENT = "ADC-transporter-evidence-audit/1.0 (research data provenance)"

MANUAL_PRIMARY_SOURCES = [
    {
        "source_id": "T11_SRC_MANUAL_001",
        "url": "https://pure.uva.nl/ws/files/1325385/100585_09.pdf",
        "filename": "cefadroxil_ABCG2_fulltext_10.1124_dmd.111.041731.pdf",
        "source_type": "author_institutional_fulltext",
        "evidence_use": "cefadroxil human ABCG2 direct vesicle negative with mock and ATP controls",
    },
    {
        "source_id": "T11_SRC_MANUAL_002",
        "url": "https://ira.lib.polyu.edu.hk/bitstream/10397/95467/1/Wong_Flavonoid_Monomers_As.pdf",
        "filename": "Ac32Az19_ABCG2_accepted_manuscript_10.1021_acs.jmedchem.1c00779.pdf",
        "source_type": "author_institutional_accepted_manuscript",
        "evidence_use": "Ac32Az19 human ABCG2 direct parent-versus-overexpressing cell negative",
    },
    {
        "source_id": "T11_SRC_MANUAL_003",
        "url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=15290871&retmode=xml",
        "filename": "pubmed_15290871_sex_steroid_ABCB1.xml",
        "source_type": "primary_article_abstract_xml",
        "evidence_use": "ethynyl estradiol, estrone, estriol and norethindrone ABCB1 direct bidirectional assay results",
    },
    {
        "source_id": "T11_SRC_MANUAL_004",
        "url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=22415933&retmode=xml",
        "filename": "pubmed_22415933_7HC_glucuronide_transport.xml",
        "source_type": "primary_article_abstract_xml",
        "evidence_use": "7-hydroxycoumarin glucuronide direct vesicle transporter selectivity result",
    },
]

TARGETS = {
    "ABCB1": {
        "target_chembl_id": "CHEMBL4302",
        "pref_name": "ATP-dependent translocase ABCB1",
        "organism": "Homo sapiens",
    },
    "ABCG2": {
        "target_chembl_id": "CHEMBL5393",
        "pref_name": "Broad substrate specificity ATP-binding cassette transporter ABCG2",
        "organism": "Homo sapiens",
    },
}

DIRECT_PATTERNS = [
    r"\bbidirectional\b",
    r"\bapical\s*(?:to|[-→])\s*basolateral\b",
    r"\bbasolateral\s*(?:to|[-→])\s*apical\b",
    r"\bpapp\b",
    r"\befflux ratio\b",
    r"\bnet efflux\b",
    r"\btransport(?:ed|ation)?\b",
    r"\bpermeab(?:ility|le)\b",
    r"\bsubstrate\b",
]
INHIBITION_PATTERNS = [r"\binhibit(?:ion|or|ory|ed|s)?\b", r"\bic50\b", r"\bki\b"]
FUNCTIONAL_PATTERNS = [
    r"\baccumulation\b",
    r"\bcytotoxic(?:ity)?\b",
    r"\bresistan(?:ce|t)\b",
    r"\batpase\b",
    r"\bcalcein\b",
    r"\brhodamine\b",
    r"\bmitoxantrone\b",
]
DIRECT_STANDARD_TYPES = {
    "papp",
    "permeability",
    "efflux ratio",
    "transport",
    "uptake",
    "substrate",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    if fields is None:
        if not rows:
            raise ValueError(f"No rows and no explicit fields for {path}")
        fields = list(rows[0])
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch_json(url: str, cache_path: Path, retries: int = 5) -> dict[str, Any]:
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read()
            parsed = json.loads(payload)
            temporary = cache_path.with_suffix(cache_path.suffix + ".tmp")
            temporary.write_bytes(payload)
            temporary.replace(cache_path)
            return parsed
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError):
            if attempt + 1 == retries:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def fetch_binary(url: str, cache_path: Path, retries: int = 5) -> None:
    if cache_path.exists():
        return
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read()
            temporary = cache_path.with_suffix(cache_path.suffix + ".tmp")
            temporary.write_bytes(payload)
            temporary.replace(cache_path)
            return
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError):
            if attempt + 1 == retries:
                raise
            time.sleep(2 ** attempt)


def manual_verified_binary_rows() -> list[dict[str, Any]]:
    common = {
        "manual_review_date": RUN_DATE,
        "published_label_retained": "yes",
        "claim_boundary": "binary call is transporter- and assay-context-specific; it is not a complete-ADC label",
    }
    rows = [
        {
            "manual_label_id": "T11_LABEL_001",
            "audit_record_id": "T10_AUD_00699",
            "transporter_gene": "ABCG2",
            "preferred_name": "CEFADROXIL",
            "source_chembl_id": "CHEMBL1644",
            "binary_label": 0,
            "evidence_grade": "A1_fulltext_direct_assay",
            "experimental_system": "human ABCG2 recombinant baculovirus-infected Sf21 plasma membrane vesicles",
            "transporter_comparison_control": "mock-transduced vesicles; with versus without ATP; ABCC1/2/3/4 positive transport comparisons",
            "test_concentration": "10 uM [3H]-cefadroxil",
            "duration": "6 min",
            "measurement": "radiolabeled vesicular uptake; average plus SD",
            "direct_result": "primary full text states ABCG2 did not transport cefadroxil",
            "replicate_qc": "Figure 2 reports one of two independent experiments and mock/ATP controls",
            "chembl_activity_id": "15450646",
            "assay_chembl_id": "CHEMBL3537693",
            "primary_doi": "10.1124/dmd.111.041731",
            "primary_pmid": "22166395",
            "primary_source_url": MANUAL_PRIMARY_SOURCES[0]["url"],
            "local_primary_file": f"raw/eleventh_batch_chembl_assay_backfill/manual_primary_sources/{MANUAL_PRIMARY_SOURCES[0]['filename']}",
            "primary_source_location": "full text Figure 2 and Results: ABCG2 did not transport cefadroxil",
            "chembl_primary_discrepancy": "ChEMBL description says 6 hrs; primary Figure 2 says 6 min; primary value retained",
            "training_use": "yes_strict_binary_core",
            "weight_recommendation": 1.0,
            **common,
        },
        {
            "manual_label_id": "T11_LABEL_002",
            "audit_record_id": "T10_AUD_00865",
            "transporter_gene": "ABCG2",
            "preferred_name": "Ac32Az19",
            "source_chembl_id": "CHEMBL4877206",
            "binary_label": 0,
            "evidence_grade": "A1_fulltext_direct_assay",
            "experimental_system": "human BCRP-overexpressing HEK293/R2 versus wild-type HEK293/pcDNA3.1 cells",
            "transporter_comparison_control": "matched parental/wild-type cells; known BCRP substrates and Ko143/GF120918 controls in the study",
            "test_concentration": "2 and 10 uM",
            "duration": "2 h",
            "measurement": "intracellular Ac32Az19 by UPLC-MS/MS",
            "direct_result": "same intracellular level in HEK293/R2 and parental cells; authors conclude not a BCRP transport substrate",
            "replicate_qc": "Figure 4D, n=3 independent experiments; mechanistic controls reported",
            "chembl_activity_id": "23238003",
            "assay_chembl_id": "CHEMBL4828495",
            "primary_doi": "10.1021/acs.jmedchem.1c00779",
            "primary_pmid": "34606270",
            "primary_source_url": MANUAL_PRIMARY_SOURCES[1]["url"],
            "local_primary_file": f"raw/eleventh_batch_chembl_assay_backfill/manual_primary_sources/{MANUAL_PRIMARY_SOURCES[1]['filename']}",
            "primary_source_location": "accepted manuscript Figure 4D and Results/Discussion",
            "chembl_primary_discrepancy": "none identified",
            "training_use": "yes_strict_binary_core",
            "weight_recommendation": 1.0,
            **common,
        },
        {
            "manual_label_id": "T11_LABEL_003",
            "audit_record_id": "T10_AUD_01051",
            "transporter_gene": "ABCG2",
            "preferred_name": "7-hydroxycoumarin glucuronide",
            "source_chembl_id": "CHEMBL1528022",
            "binary_label": 0,
            "evidence_grade": "A2_primary_abstract_plus_traceable_assay",
            "experimental_system": "human ABCG2 recombinant membrane vesicle drug-transport assay",
            "transporter_comparison_control": "ATP versus AMP; multiple efflux transporter panel",
            "test_concentration": "400 uM",
            "duration": "10 min",
            "measurement": "reversed-phase HPLC; AMP result below detection",
            "direct_result": "primary abstract reports 7-HC-G transported by MRP3/4, not BCRP; ChEMBL direct ABCG2 call is Not active",
            "replicate_qc": "full replicate and vesicle QC not exposed in the abstract record",
            "chembl_activity_id": "15462557 | 15462565",
            "assay_chembl_id": "CHEMBL3530760 | CHEMBL3533970",
            "primary_doi": "10.1124/dmd.111.044438",
            "primary_pmid": "22415933",
            "primary_source_url": MANUAL_PRIMARY_SOURCES[3]["url"],
            "local_primary_file": f"raw/eleventh_batch_chembl_assay_backfill/manual_primary_sources/{MANUAL_PRIMARY_SOURCES[3]['filename']}",
            "primary_source_location": "primary abstract plus ChEMBL assay descriptions",
            "chembl_primary_discrepancy": "none identified; full text still required for A1",
            "training_use": "yes_sensitivity_lower_weight",
            "weight_recommendation": 0.75,
            **common,
        },
    ]
    steroid_rows = [
        ("T11_LABEL_004", "T10_AUD_00002", "ETHINYL ESTRADIOL", "CHEMBL691", 1, "10.3", "10997684"),
        ("T11_LABEL_005", "T10_AUD_00125", "ESTRONE", "CHEMBL1405", 1, "6.9", "10997683"),
        ("T11_LABEL_006", "T10_AUD_00124", "ESTRIOL", "CHEMBL193482", 1, "2.8", "10997682"),
        ("T11_LABEL_007", "T10_AUD_00001", "NORETHINDRONE", "CHEMBL1162", 0, "not reported; explicit non-substrate", "11001531"),
    ]
    for label_id, audit_id, name, chembl_id, label, result, activity_id in steroid_rows:
        rows.append(
            {
                "manual_label_id": label_id,
                "audit_record_id": audit_id,
                "transporter_gene": "ABCB1",
                "preferred_name": name,
                "source_chembl_id": chembl_id,
                "binary_label": label,
                "evidence_grade": "A2_primary_abstract_plus_traceable_assay",
                "experimental_system": "MDR1-transfected MDCK versus parental MDCK bidirectional transport",
                "transporter_comparison_control": "with and without P-gp inhibitor GG918; human MDR1 ATPase and control membranes also studied",
                "test_concentration": "not exposed in abstract",
                "duration": "not exposed in abstract",
                "measurement": "LC/MS bidirectional transport",
                "direct_result": f"primary abstract efflux ratio/result: {result}",
                "replicate_qc": "direct method and quantitative result/author call reported; full row-level methods still desirable",
                "chembl_activity_id": activity_id,
                "assay_chembl_id": "CHEMBL2078577",
                "primary_doi": "10.1023/b:pham.0000033017.52484.81",
                "primary_pmid": "15290871",
                "primary_source_url": MANUAL_PRIMARY_SOURCES[2]["url"],
                "local_primary_file": f"raw/eleventh_batch_chembl_assay_backfill/manual_primary_sources/{MANUAL_PRIMARY_SOURCES[2]['filename']}",
                "primary_source_location": "primary abstract Methods and Results plus ChEMBL TP_TRANSPORTER row",
                "chembl_primary_discrepancy": "none identified; test concentration unavailable in abstract",
                "training_use": "yes_sensitivity_lower_weight",
                "weight_recommendation": 0.75,
                **common,
            }
        )
    return rows


def chunks(values: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def download_target_activities(gene: str, target_id: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    manifest: list[dict[str, Any]] = []
    offset = 0
    while True:
        query = urllib.parse.urlencode(
            {"target_chembl_id": target_id, "limit": PAGE_LIMIT, "offset": offset}
        )
        url = f"{BASE}/activity.json?{query}"
        cache = RAW / "activity_pages" / f"{gene}_{target_id}_offset_{offset:06d}.json"
        page = fetch_json(url, cache)
        page_rows = page.get("activities", [])
        rows.extend(page_rows)
        manifest.append(
            {
                "raw_file": str(cache.relative_to(DATA)),
                "request_url": url,
                "retrieved_date": RUN_DATE,
                "response_rows": len(page_rows),
                "file_size_bytes": cache.stat().st_size,
                "sha256": sha256(cache),
            }
        )
        total = int(page.get("page_meta", {}).get("total_count", len(rows)))
        offset += len(page_rows)
        if not page_rows or offset >= total:
            break
    return rows, manifest


def batch_download(
    resource: str,
    id_field: str,
    ids: list[str],
    filter_field: str | None = None,
    batch_size: int = 40,
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    output: dict[str, dict[str, Any]] = {}
    manifest: list[dict[str, Any]] = []
    plural = f"{resource}s"
    for batch_index, batch in enumerate(chunks(sorted(set(ids)), batch_size), start=1):
        offset = 0
        while True:
            query = urllib.parse.urlencode(
                {f"{filter_field or id_field}__in": ",".join(batch), "limit": 1000, "offset": offset}
            )
            url = f"{BASE}/{resource}.json?{query}"
            cache = RAW / f"{resource}_batches" / f"{resource}_{batch_index:04d}_offset_{offset:06d}.json"
            page = fetch_json(url, cache)
            rows = page.get(plural, [])
            for row in rows:
                output[str(row[id_field])] = row
            manifest.append(
                {
                    "raw_file": str(cache.relative_to(DATA)),
                    "request_url": url,
                    "retrieved_date": RUN_DATE,
                    "response_rows": len(rows),
                    "file_size_bytes": cache.stat().st_size,
                    "sha256": sha256(cache),
                }
            )
            total = int(page.get("page_meta", {}).get("total_count", len(rows)))
            offset += len(rows)
            if not rows or offset >= total:
                break
    return output, manifest


def candidate_chembl_ids(row: dict[str, str]) -> list[str]:
    return sorted(set(re.findall(r"CHEMBL\d+", row.get("source_pids", ""), flags=re.IGNORECASE)))


def combined_text(activity: dict[str, Any]) -> str:
    values = [
        activity.get("assay_description"),
        activity.get("activity_comment"),
        activity.get("standard_type"),
        activity.get("type"),
        activity.get("standard_text_value"),
        activity.get("text_value"),
    ]
    return " ".join(str(value) for value in values if value).lower()


def activity_evidence_stratum(activity: dict[str, Any]) -> tuple[str, str]:
    text = combined_text(activity)
    description = str(activity.get("assay_description") or "").lower()
    activity_comment = str(activity.get("activity_comment") or "").lower()
    standard_type = str(activity.get("standard_type") or "").strip().lower()
    direct_hits = [pattern for pattern in DIRECT_PATTERNS if re.search(pattern, text)]
    inhibitor_hits = [pattern for pattern in INHIBITION_PATTERNS if re.search(pattern, text)]
    functional_hits = [pattern for pattern in FUNCTIONAL_PATTERNS if re.search(pattern, text)]
    explicit_substrate_assay = "substrate activity" in description or "substrate [" in activity_comment or "not a substrate" in activity_comment
    if inhibitor_hits and not explicit_substrate_assay:
        return "inhibitor_or_probe_activity", "inhibition/probe endpoint cannot establish transport of the test compound"
    if standard_type in DIRECT_STANDARD_TYPES or direct_hits or explicit_substrate_assay:
        return "direct_transport_candidate", "direct transport/permeability/substrate terminology; primary full text still required"
    if inhibitor_hits:
        return "inhibitor_or_probe_activity", "inhibition endpoint cannot establish test compound as transporter substrate"
    if functional_hits:
        return "functional_indirect_activity", "accumulation/resistance/ATPase or probe readout is indirect for substrate labeling"
    return "other_target_activity", "target-mapped activity lacks direct substrate terminology"


def explicit_direct_substrate_call(activity: dict[str, Any], stratum: str) -> tuple[str, str]:
    """Return 0/1 only for an explicit call in a direct-substrate assay context."""
    if stratum != "direct_transport_candidate":
        return "", "not_a_direct_transport_candidate"
    description = str(activity.get("assay_description") or "").strip().lower()
    call_text = " ".join(
        str(activity.get(field) or "").strip().lower()
        for field in ("activity_comment", "standard_text_value", "text_value")
    )
    if "not a substrate" in call_text or "substrate [-]" in call_text:
        return "0", "explicit_not_substrate_call"
    if "substrate [+]" in call_text:
        return "1", "explicit_substrate_positive_call"
    direct_call_assay = any(
        term in description
        for term in (
            "substrate activity",
            "transepithelial transport",
            "transporter-mediated drug uptake",
            "assessed as drug transport",
        )
    )
    if direct_call_assay and "not active" in call_text:
        return "0", "not_active_in_explicit_direct_substrate_assay"
    if direct_call_assay and re.search(r"(?:^|\s)active(?:\s|$)", call_text) and "not active" not in call_text:
        return "1", "active_in_explicit_direct_substrate_assay"
    return "", "numeric_or_ambiguous_direct_result_requires_fulltext_rule"


def value_or_blank(value: Any) -> Any:
    return "" if value is None else value


def main() -> None:
    candidates = read_csv(INPUT)
    if len(candidates) != 512:
        raise ValueError(f"Expected 512 Tier A/B candidates, found {len(candidates)}")

    candidate_by_gene_and_id: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in candidates:
        for chembl_id in candidate_chembl_ids(row):
            candidate_by_gene_and_id[(row["transporter_gene"], chembl_id.upper())].append(row)

    raw_manifest: list[dict[str, Any]] = []
    target_activity_rows: dict[str, list[dict[str, Any]]] = {}
    for gene, metadata in TARGETS.items():
        rows, manifest = download_target_activities(gene, metadata["target_chembl_id"])
        target_activity_rows[gene] = rows
        raw_manifest.extend(manifest)

    matched: list[tuple[dict[str, str], dict[str, Any], str]] = []
    for gene, activities in target_activity_rows.items():
        for activity in activities:
            exact_id = str(activity.get("molecule_chembl_id") or "").upper()
            parent_id = str(activity.get("parent_molecule_chembl_id") or "").upper()
            matches = candidate_by_gene_and_id.get((gene, exact_id), [])
            match_type = "exact_molecule_chembl_id"
            if not matches and parent_id:
                matches = candidate_by_gene_and_id.get((gene, parent_id), [])
                match_type = "parent_molecule_chembl_id"
            for candidate in matches:
                matched.append((candidate, activity, match_type))

    assay_ids = sorted({str(activity.get("assay_chembl_id")) for _, activity, _ in matched if activity.get("assay_chembl_id")})
    document_ids = sorted({str(activity.get("document_chembl_id")) for _, activity, _ in matched if activity.get("document_chembl_id")})
    assays, assay_manifest = batch_download("assay", "assay_chembl_id", assay_ids)
    documents, document_manifest = batch_download("document", "document_chembl_id", document_ids)
    compound_records, compound_record_manifest = batch_download(
        "compound_record",
        "record_id",
        sorted({chembl_id for _, chembl_id in candidate_by_gene_and_id}),
        filter_field="molecule_chembl_id",
        batch_size=10,
    )
    raw_manifest.extend(assay_manifest)
    raw_manifest.extend(document_manifest)
    raw_manifest.extend(compound_record_manifest)

    manual_source_registry: list[dict[str, Any]] = []
    manual_source_dir = RAW / "manual_primary_sources"
    for source in MANUAL_PRIMARY_SOURCES:
        local_path = manual_source_dir / source["filename"]
        fetch_binary(source["url"], local_path)
        registry_row = {
            "source_id": source["source_id"],
            "source_type": source["source_type"],
            "evidence_use": source["evidence_use"],
            "source_url": source["url"],
            "local_file": str(local_path.relative_to(DATA)),
            "retrieved_date": RUN_DATE,
            "file_size_bytes": local_path.stat().st_size,
            "sha256": sha256(local_path),
        }
        manual_source_registry.append(registry_row)
        raw_manifest.append(
            {
                "raw_file": registry_row["local_file"],
                "request_url": registry_row["source_url"],
                "retrieved_date": registry_row["retrieved_date"],
                "response_rows": 1,
                "file_size_bytes": registry_row["file_size_bytes"],
                "sha256": registry_row["sha256"],
            }
        )

    activity_fields = [
        "backfill_record_id", "audit_record_id", "transporter_gene", "published_binary_label",
        "preferred_names", "candidate_source_pids", "candidate_standardized_inchikey",
        "candidate_match_type", "activity_id", "activity_record_id", "activity_molecule_chembl_id", "parent_molecule_chembl_id",
        "source_compound_key", "source_compound_name",
        "assay_chembl_id", "document_chembl_id", "evidence_stratum", "evidence_reason",
        "explicit_direct_substrate_call", "explicit_call_basis",
        "assay_description", "activity_comment", "assay_type", "bao_label", "standard_type",
        "standard_relation", "standard_value", "standard_units", "standard_text_value", "pchembl_value",
        "data_validity_comment", "potential_duplicate", "target_chembl_id", "target_organism",
        "assay_organism", "assay_confidence_score", "assay_relationship_type", "assay_bao_format",
        "document_doi", "document_pubmed_id", "document_title", "document_year", "source_activity_url",
        "tier_A_auto_upgrade", "manual_review_requirement",
    ]
    activity_output: list[dict[str, Any]] = []
    for candidate, activity, match_type in sorted(
        matched,
        key=lambda item: (
            item[0]["transporter_gene"], item[0]["audit_record_id"], int(item[1].get("activity_id") or 0)
        ),
    ):
        stratum, reason = activity_evidence_stratum(activity)
        explicit_call, explicit_call_basis = explicit_direct_substrate_call(activity, stratum)
        assay = assays.get(str(activity.get("assay_chembl_id") or ""), {})
        document = documents.get(str(activity.get("document_chembl_id") or ""), {})
        compound_record = compound_records.get(str(activity.get("record_id") or ""), {})
        activity_output.append(
            {
                "backfill_record_id": f"T11_ACT_{len(activity_output) + 1:06d}",
                "audit_record_id": candidate["audit_record_id"],
                "transporter_gene": candidate["transporter_gene"],
                "published_binary_label": candidate["published_binary_label"],
                "preferred_names": candidate["preferred_names"],
                "candidate_source_pids": candidate["source_pids"],
                "candidate_standardized_inchikey": candidate["standardized_inchikey"],
                "candidate_match_type": match_type,
                "activity_id": value_or_blank(activity.get("activity_id")),
                "activity_record_id": value_or_blank(activity.get("record_id")),
                "activity_molecule_chembl_id": value_or_blank(activity.get("molecule_chembl_id")),
                "parent_molecule_chembl_id": value_or_blank(activity.get("parent_molecule_chembl_id")),
                "source_compound_key": value_or_blank(compound_record.get("compound_key")),
                "source_compound_name": value_or_blank(compound_record.get("compound_name")),
                "assay_chembl_id": value_or_blank(activity.get("assay_chembl_id")),
                "document_chembl_id": value_or_blank(activity.get("document_chembl_id")),
                "evidence_stratum": stratum,
                "evidence_reason": reason,
                "explicit_direct_substrate_call": explicit_call,
                "explicit_call_basis": explicit_call_basis,
                "assay_description": value_or_blank(activity.get("assay_description")),
                "activity_comment": value_or_blank(activity.get("activity_comment")),
                "assay_type": value_or_blank(activity.get("assay_type")),
                "bao_label": value_or_blank(activity.get("bao_label")),
                "standard_type": value_or_blank(activity.get("standard_type")),
                "standard_relation": value_or_blank(activity.get("standard_relation")),
                "standard_value": value_or_blank(activity.get("standard_value")),
                "standard_units": value_or_blank(activity.get("standard_units")),
                "standard_text_value": value_or_blank(activity.get("standard_text_value")),
                "pchembl_value": value_or_blank(activity.get("pchembl_value")),
                "data_validity_comment": value_or_blank(activity.get("data_validity_comment")),
                "potential_duplicate": value_or_blank(activity.get("potential_duplicate")),
                "target_chembl_id": value_or_blank(activity.get("target_chembl_id")),
                "target_organism": value_or_blank(activity.get("target_organism")),
                "assay_organism": value_or_blank(assay.get("assay_organism")),
                "assay_confidence_score": value_or_blank(assay.get("confidence_score")),
                "assay_relationship_type": value_or_blank(assay.get("relationship_type")),
                "assay_bao_format": value_or_blank(assay.get("bao_format")),
                "document_doi": value_or_blank(document.get("doi")),
                "document_pubmed_id": value_or_blank(document.get("pubmed_id")),
                "document_title": value_or_blank(document.get("title")),
                "document_year": value_or_blank(document.get("year")),
                "source_activity_url": f"https://www.ebi.ac.uk/chembl/explore/activity/{activity.get('activity_id')}",
                "tier_A_auto_upgrade": "no",
                "manual_review_requirement": "verify exact tested entity, direct endpoint, controls, assay validity, numeric result and author interpretation in primary full text",
            }
        )

    assay_output: list[dict[str, Any]] = []
    for assay_id in sorted(assays):
        assay = assays[assay_id]
        assay_output.append(
            {
                "assay_chembl_id": assay_id,
                "description": value_or_blank(assay.get("description")),
                "assay_type": value_or_blank(assay.get("assay_type")),
                "assay_organism": value_or_blank(assay.get("assay_organism")),
                "target_chembl_id": value_or_blank(assay.get("target_chembl_id")),
                "confidence_score": value_or_blank(assay.get("confidence_score")),
                "relationship_type": value_or_blank(assay.get("relationship_type")),
                "bao_format": value_or_blank(assay.get("bao_format")),
                "cell_chembl_id": value_or_blank(assay.get("cell_chembl_id")),
                "tissue_chembl_id": value_or_blank(assay.get("tissue_chembl_id")),
                "document_chembl_id": value_or_blank(assay.get("document_chembl_id")),
                "source_url": f"https://www.ebi.ac.uk/chembl/api/data/assay/{assay_id}.json",
            }
        )

    document_output: list[dict[str, Any]] = []
    for document_id in sorted(documents):
        document = documents[document_id]
        doi = value_or_blank(document.get("doi"))
        document_output.append(
            {
                "document_chembl_id": document_id,
                "doi": doi,
                "pubmed_id": value_or_blank(document.get("pubmed_id")),
                "title": value_or_blank(document.get("title")),
                "year": value_or_blank(document.get("year")),
                "journal": value_or_blank(document.get("journal")),
                "volume": value_or_blank(document.get("volume")),
                "first_page": value_or_blank(document.get("first_page")),
                "last_page": value_or_blank(document.get("last_page")),
                "primary_source_url": f"https://doi.org/{doi}" if doi else f"https://www.ebi.ac.uk/chembl/api/data/document/{document_id}.json",
            }
        )

    by_candidate: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in activity_output:
        by_candidate[str(row["audit_record_id"])].append(row)

    triage_output: list[dict[str, Any]] = []
    for candidate in candidates:
        rows = by_candidate.get(candidate["audit_record_id"], [])
        strata = Counter(str(row["evidence_stratum"]) for row in rows)
        dois = sorted({str(row["document_doi"]) for row in rows if row["document_doi"]})
        direct_rows = [row for row in rows if row["evidence_stratum"] == "direct_transport_candidate"]
        direct_positive = [row for row in direct_rows if row["explicit_direct_substrate_call"] == "1"]
        direct_negative = [row for row in direct_rows if row["explicit_direct_substrate_call"] == "0"]
        direct_ambiguous = [row for row in direct_rows if row["explicit_direct_substrate_call"] == ""]
        resolved_names = sorted({str(row["source_compound_name"]) for row in direct_rows if row["source_compound_name"]})
        explicit_calls = {row["explicit_direct_substrate_call"] for row in direct_rows if row["explicit_direct_substrate_call"]}
        reconstructed_label = next(iter(explicit_calls)) if len(explicit_calls) == 1 else ""
        reconstructed_concordant = reconstructed_label != "" and reconstructed_label == candidate["published_binary_label"]
        high_confidence_traceable = any(
            row["explicit_direct_substrate_call"] == reconstructed_label
            and row["candidate_match_type"] == "exact_molecule_chembl_id"
            and row["target_organism"] == "Homo sapiens"
            and row["document_doi"]
            and str(row["assay_confidence_score"]) in {"8", "9"}
            and not row["data_validity_comment"]
            for row in direct_rows
        ) if reconstructed_label else False
        gap = "P0_ABCG2_confirmed_negative_or_weak_substrate" if (
            candidate["transporter_gene"] == "ABCG2" and candidate["published_binary_label"] == "0"
        ) else "P0_ABCB1_independent_positive" if (
            candidate["transporter_gene"] == "ABCB1" and candidate["published_binary_label"] == "1"
        ) else "P1_other_class_provenance_backfill"
        priority_score = 0
        if gap.startswith("P0"):
            priority_score += 40
        priority_score += min(25, 10 * len(direct_rows))
        if dois:
            priority_score += 10
        if any(row["standard_value"] != "" for row in direct_rows):
            priority_score += 10
        if not rows:
            priority_score -= 5
        triage_output.append(
            {
                "audit_record_id": candidate["audit_record_id"],
                "transporter_gene": candidate["transporter_gene"],
                "published_binary_label": candidate["published_binary_label"],
                "preferred_names": candidate["preferred_names"],
                "resolved_source_compound_names": " | ".join(resolved_names),
                "source_pids": candidate["source_pids"],
                "standardized_inchikey": candidate["standardized_inchikey"],
                "label_gap_priority": gap,
                "n_target_activity_hits": len(rows),
                "n_unique_assays": len({row["assay_chembl_id"] for row in rows if row["assay_chembl_id"]}),
                "n_unique_documents": len({row["document_chembl_id"] for row in rows if row["document_chembl_id"]}),
                "n_direct_transport_candidates": strata["direct_transport_candidate"],
                "n_explicit_direct_positive_calls": len(direct_positive),
                "n_explicit_direct_negative_calls": len(direct_negative),
                "n_ambiguous_direct_results": len(direct_ambiguous),
                "n_inhibitor_or_probe_activities": strata["inhibitor_or_probe_activity"],
                "n_functional_indirect_activities": strata["functional_indirect_activity"],
                "n_other_target_activities": strata["other_target_activity"],
                "document_dois": " | ".join(dois),
                "api_backfill_status": "direct_candidate_found" if direct_rows else "target_activity_found_no_direct_candidate" if rows else "no_human_target_activity_hit",
                "reconstructed_binary_label": reconstructed_label,
                "reconstructed_label_matches_published": "yes" if reconstructed_concordant else "no" if reconstructed_label else "not_determined",
                "explicit_direct_call_conflict": "yes" if len(explicit_calls) > 1 else "no",
                "tier_A_upgrade_status": (
                    "traceable_direct_call_concordant_pending_fulltext_QC"
                    if reconstructed_concordant and high_confidence_traceable
                    else "direct_call_conflict_requires_resolution"
                    if len(explicit_calls) > 1
                    else "direct_call_disagrees_with_published_requires_resolution"
                    if reconstructed_label and not reconstructed_concordant
                    else "direct_numeric_or_ambiguous_result_requires_fulltext_rule"
                    if direct_rows
                    else "not_upgradable_from_current_target_activity_hit"
                ),
                "priority_score": priority_score,
                "published_label_retained": "yes",
                "model_training_change": "none_until_manual_evidence_review",
            }
        )

    triage_output.sort(key=lambda row: (-int(row["priority_score"]), row["transporter_gene"], row["audit_record_id"]))
    for rank, row in enumerate(triage_output, start=1):
        row["manual_review_rank"] = rank

    reconstructed_output = [
        {
            "reconstruction_id": f"T11_RECON_{index:04d}",
            "audit_record_id": row["audit_record_id"],
            "transporter_gene": row["transporter_gene"],
            "preferred_names": row["preferred_names"],
            "resolved_source_compound_names": row["resolved_source_compound_names"],
            "source_pids": row["source_pids"],
            "published_binary_label": row["published_binary_label"],
            "reconstructed_binary_label": row["reconstructed_binary_label"],
            "reconstructed_label_matches_published": row["reconstructed_label_matches_published"],
            "n_explicit_direct_positive_calls": row["n_explicit_direct_positive_calls"],
            "n_explicit_direct_negative_calls": row["n_explicit_direct_negative_calls"],
            "n_ambiguous_direct_results": row["n_ambiguous_direct_results"],
            "document_dois": row["document_dois"],
            "tier_A_upgrade_status": row["tier_A_upgrade_status"],
            "strict_tier_A_training_use": "no_pending_primary_fulltext_QC",
        }
        for index, row in enumerate(
            [row for row in triage_output if row["tier_A_upgrade_status"] == "traceable_direct_call_concordant_pending_fulltext_QC"],
            start=1,
        )
    ]
    conflict_output = [
        {
            "conflict_id": f"T11_CONFLICT_{index:04d}",
            "audit_record_id": row["audit_record_id"],
            "transporter_gene": row["transporter_gene"],
            "preferred_names": row["preferred_names"],
            "resolved_source_compound_names": row["resolved_source_compound_names"],
            "source_pids": row["source_pids"],
            "published_binary_label": row["published_binary_label"],
            "reconstructed_binary_label": row["reconstructed_binary_label"],
            "n_explicit_direct_positive_calls": row["n_explicit_direct_positive_calls"],
            "n_explicit_direct_negative_calls": row["n_explicit_direct_negative_calls"],
            "document_dois": row["document_dois"],
            "conflict_type": "opposing_explicit_direct_calls" if row["explicit_direct_call_conflict"] == "yes" else "reconstructed_call_disagrees_with_published",
            "required_action": "review primary methods, exact tested form, concentration, assay controls and independent replication; retain all contexts",
            "model_use": "exclude_from_strict_binary_training_until_resolved",
        }
        for index, row in enumerate(
            [
                row for row in triage_output
                if row["explicit_direct_call_conflict"] == "yes"
                or row["tier_A_upgrade_status"] == "direct_call_disagrees_with_published_requires_resolution"
            ],
            start=1,
        )
    ]

    queue_output: list[dict[str, Any]] = []
    for row in triage_output:
        if row["n_direct_transport_candidates"] or (
            str(row["label_gap_priority"]).startswith("P0") and len(queue_output) < 150
        ):
            queue_output.append(
                {
                    "queue_id": f"T11_QUEUE_{len(queue_output) + 1:04d}",
                    "manual_review_rank": row["manual_review_rank"],
                    "audit_record_id": row["audit_record_id"],
                    "transporter_gene": row["transporter_gene"],
                    "published_binary_label": row["published_binary_label"],
                    "preferred_names": row["preferred_names"],
                    "source_pids": row["source_pids"],
                    "label_gap_priority": row["label_gap_priority"],
                    "n_direct_transport_candidates": row["n_direct_transport_candidates"],
                    "n_target_activity_hits": row["n_target_activity_hits"],
                    "document_dois": row["document_dois"],
                    "tier_A_upgrade_status": row["tier_A_upgrade_status"],
                    "required_primary_review": "exact identity; human transporter; direct substrate endpoint; parental/empty-vector or inhibitor control; monolayer/assay QC; numeric result; author interpretation",
                    "reject_as_binary_truth_if": "inhibitor-only; probe accumulation only; cytotoxicity-only; ATPase-only; species mismatch; no tested negative; unresolved assay conflict",
                }
            )

    write_csv(PROCESSED / "chembl_candidate_activity_hits_v1.csv", activity_output, activity_fields)
    write_csv(PROCESSED / "chembl_candidate_assay_registry_v1.csv", assay_output)
    write_csv(PROCESSED / "chembl_candidate_document_registry_v1.csv", document_output)
    write_csv(PROCESSED / "chembl_candidate_provenance_triage_v1.csv", triage_output)
    write_csv(
        PROCESSED / "chembl_reconstructed_binary_label_candidates_v1.csv",
        reconstructed_output,
        fields=[
            "reconstruction_id", "audit_record_id", "transporter_gene", "preferred_names", "resolved_source_compound_names", "source_pids",
            "published_binary_label", "reconstructed_binary_label", "reconstructed_label_matches_published",
            "n_explicit_direct_positive_calls", "n_explicit_direct_negative_calls", "n_ambiguous_direct_results",
            "document_dois", "tier_A_upgrade_status", "strict_tier_A_training_use",
        ],
    )
    write_csv(
        PROCESSED / "chembl_reconstructed_label_conflicts_v1.csv",
        conflict_output,
        fields=[
            "conflict_id", "audit_record_id", "transporter_gene", "preferred_names", "resolved_source_compound_names", "source_pids",
            "published_binary_label", "reconstructed_binary_label", "n_explicit_direct_positive_calls",
            "n_explicit_direct_negative_calls", "document_dois", "conflict_type", "required_action", "model_use",
        ],
    )
    write_csv(PROCESSED / "eleventh_batch_primary_source_queue_v1.csv", queue_output)
    manual_labels = manual_verified_binary_rows()
    write_csv(PROCESSED / "eleventh_batch_manual_verified_binary_labels_v1.csv", manual_labels)
    write_csv(PROCESSED / "eleventh_batch_binary_training_supplement_v1.csv", manual_labels)
    write_csv(PROCESSED / "eleventh_batch_manual_primary_source_registry.csv", manual_source_registry)
    write_csv(PROCESSED / "eleventh_batch_raw_download_log.csv", raw_manifest)

    priority_counts = Counter(row["label_gap_priority"] for row in triage_output)
    direct_by_gene_label: dict[str, dict[str, int]] = defaultdict(dict)
    for gene in ("ABCB1", "ABCG2"):
        for label in ("0", "1"):
            direct_by_gene_label[gene][label] = sum(
                row["transporter_gene"] == gene
                and row["published_binary_label"] == label
                and int(row["n_direct_transport_candidates"]) > 0
                for row in triage_output
            )
    summary = {
        "batch": 11,
        "date": RUN_DATE,
        "purpose": "chembl_row_level_activity_assay_document_backfill",
        "input_tier_AB_candidates": len(candidates),
        "candidate_chembl_identifiers": len(candidate_by_gene_and_id),
        "downloaded_target_activities": {gene: len(rows) for gene, rows in target_activity_rows.items()},
        "matched_activity_rows": len(activity_output),
        "matched_candidates": sum(int(row["n_target_activity_hits"]) > 0 for row in triage_output),
        "candidates_with_direct_transport_activity_candidates": sum(int(row["n_direct_transport_candidates"]) > 0 for row in triage_output),
        "direct_candidate_counts_by_gene_and_published_label": direct_by_gene_label,
        "unique_assays": len(assay_output),
        "unique_documents": len(document_output),
        "compound_record_rows": len(compound_records),
        "documents_with_doi": sum(bool(row["doi"]) for row in document_output),
        "traceable_concordant_direct_call_candidates": len(reconstructed_output),
        "explicit_direct_call_conflicts_or_disagreements": len(conflict_output),
        "priority_queue_rows": len(queue_output),
        "priority_class_counts": dict(sorted(priority_counts.items())),
        "manual_verified_binary_labels": len(manual_labels),
        "manual_verified_by_gene_and_label": {
            gene: {
                label: sum(
                    row["transporter_gene"] == gene and str(row["binary_label"]) == label
                    for row in manual_labels
                )
                for label in ("0", "1")
            }
            for gene in ("ABCB1", "ABCG2")
        },
        "strict_A1_binary_labels": sum(
            str(row["evidence_grade"]).startswith("A1_") for row in manual_labels
        ),
        "A2_binary_labels": sum(
            str(row["evidence_grade"]).startswith("A2_") for row in manual_labels
        ),
        "manual_primary_source_files": len(manual_source_registry),
        "primary_source_discrepancies": sum(
            str(row["chembl_primary_discrepancy"]).lower().startswith("chembl description")
            for row in manual_labels
        ),
        "automatic_tier_A_upgrades": 0,
        "claim_boundary": "ChEMBL target activity hits are retrieval leads, not automatic substrate truth. Tier A requires primary-full-text confirmation of exact entity, direct endpoint, controls, assay validity, and result.",
    }
    (PROCESSED / "eleventh_batch_chembl_backfill_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

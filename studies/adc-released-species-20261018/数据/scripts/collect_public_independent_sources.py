#!/usr/bin/env python3
"""Collect two independent, openly downloadable ABC transporter label sources.

The outputs intentionally keep the original multiclass call and provenance.  The
binary labels are a *derived* analysis view:

* Hazai et al. (2013, PMC3641962): the table section is the curated ground-truth
  class (human wild-type ABCG2); the PDF's ``Predicted property`` column is kept
  separately and is not used as the label.
* Lavecchia et al. (2019, PMC6571636): ``substrate`` is positive; the source's
  ``non_substrate`` and ``Pgp_nonInhibitor`` classes are merged into a negative
  ``non-active`` view, exactly as described by the paper.  ``Pgp_Inhibitor`` is
  excluded from the binary view instead of being silently treated as negative.

No record is silently upgraded to a Tier-A direct assay truth.  These datasets
are independent-source candidates / sensitivity or external-holdout data until
row-level primary-assay review and leakage checks are complete.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from datetime import date
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw" / "public_independent"
PROCESSED = ROOT / "processed"
ACQUIRED_AT = date.today().isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download_if_missing(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return
    tmp = path.with_suffix(path.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(path)


def ensure_sources(download: bool = False) -> None:
    """Download only when requested; parsing works offline from tracked raw files."""
    RAW.mkdir(parents=True, exist_ok=True)
    specs = [
        (
            "PMC3641962",
            "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3641962/fullTextXML",
            RAW / "PMC3641962_fullText.xml",
        ),
        (
            "PMC3641962_supplementaryFiles",
            "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3641962/supplementaryFiles",
            RAW / "PMC3641962_supplementaryFiles.zip",
        ),
        (
            "PMC6571636",
            "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6571636/fullTextXML",
            RAW / "PMC6571636_fullText.xml",
        ),
        (
            "PMC6571636_supplementaryFiles",
            "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6571636/supplementaryFiles",
            RAW / "PMC6571636_supplementaryFiles.zip",
        ),
    ]
    if download:
        for _, url, path in specs:
            print(f"download {url} -> {path}", file=sys.stderr)
            download_if_missing(url, path)

    # Extracted files are kept alongside the archive; extraction is idempotent.
    for pmc in ("PMC3641962", "PMC6571636"):
        archive = RAW / f"{pmc}_supplementaryFiles.zip"
        out = RAW / f"{pmc}_supplementary"
        if not archive.exists():
            continue
        out.mkdir(exist_ok=True)
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(out)


def _clean_text(value: str) -> str:
    value = value.replace("\f", " ")
    value = re.sub(r"\s+", " ", value)
    return value.strip()


def parse_hazai_table(pdf_path: Path) -> List[Dict[str, object]]:
    """Parse Additional file 1's table while retaining source predicted calls.

    The PDF has a few names wrapped across lines (e.g. Trp-P-1).  We join lines
    up to the terminal ``S/NS + tr/te/e`` token and remove only the final numeric
    reference bracket, leaving chemical-name brackets intact.
    """
    txt_path = pdf_path.with_suffix(".txt")
    if not txt_path.exists() or txt_path.stat().st_mtime < pdf_path.stat().st_mtime:
        subprocess.run(
            ["pdftotext", "-layout", str(pdf_path), str(txt_path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    lines = txt_path.read_text(errors="replace").splitlines()
    records: List[Dict[str, object]] = []
    section: Optional[str] = None
    acc: List[str] = []
    skip_headers = {
        "Compound",
        "Reference Predicted",
        "property set",
        "Predicted property",
        "Data",
    }
    for raw in lines:
        s = _clean_text(raw)
        if not s or re.fullmatch(r"\d+", s):
            continue
        if s.startswith("Substrates of wild-type BCRP"):
            section = "substrate"
            acc = []
            continue
        if s.startswith("Non-substrates of wild-type BCRP"):
            section = "non_substrate"
            acc = []
            continue
        if section is None or any(s.startswith(h) for h in skip_headers):
            continue
        # The table ends before the numbered reference list.  No reference-list
        # line has the terminal class tokens, so simply stop after all 263 rows.
        acc.append(s)
        m_label = re.search(r"\b(NS|S)\s+(tr|te|e)\s*$", s)
        if not m_label:
            continue
        rec = _clean_text(" ".join(acc))
        acc = []
        label_start = rec.rfind(m_label.group(0))
        prefix = rec[:label_start].strip()
        # Pick the last numeric reference bracket; brackets in names such as
        # [4,5-f] are not numeric-only and are therefore preserved.
        refs = list(re.finditer(r"\[\s*\d+(?:\s*,\s*\d+)*\s*\]", prefix))
        if not refs:
            raise ValueError(f"Cannot parse Hazai references from row: {rec}")
        ref_match = refs[-1]
        name = _clean_text(prefix[: ref_match.start()] + " " + prefix[ref_match.end() :])
        records.append(
            {
                "source_section_class": 1 if section == "substrate" else 0,
                "compound_name": name,
                "source_predicted_property": m_label.group(1),
                "source_data_set": m_label.group(2),
                "source_reference_ids": ref_match.group(0).strip("[] "),
            }
        )
        if section == "non_substrate" and len(records) >= 263:
            break
    if len(records) != 263:
        raise ValueError(f"Expected 263 Hazai rows, parsed {len(records)}")
    if sum(int(r["source_section_class"]) for r in records) != 164:
        raise ValueError("Hazai substrate section count is not 164")
    return records


def read_sdf_records(path: Path) -> List[Dict[str, object]]:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    out: List[Dict[str, object]] = []
    for idx, mol in enumerate(supplier):
        if mol is None:
            out.append({"sdf_record_index": idx, "sdf_parse_status": "failed"})
            continue
        # First line is generally a PubChem CID but a few records are named
        # OpenBabel; preserve both the raw name and the parsed structure.
        raw_name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        cid = int(raw_name) if raw_name.strip().isdigit() else None
        try:
            # OpenBabel SDFs contain explicit hydrogens.  Remove them before
            # emitting the structure so downstream standardization sees the
            # same implicit-H representation as PubChem/Daood.
            mol_props = Chem.RemoveHs(mol)
            smiles = Chem.MolToSmiles(mol_props, canonical=True, isomericSmiles=True)
        except Exception:
            smiles = ""
        try:
            inchikey = inchi.MolToInchiKey(mol_props) if smiles else ""
        except Exception:
            inchikey = ""
        out.append(
            {
                "sdf_record_index": idx,
                "sdf_mol_name": raw_name,
                "source_pubchem_cid": cid,
                "source_smiles_from_sdf": smiles,
                "source_inchikey": inchikey,
                "sdf_parse_status": "ok" if smiles else "failed",
            }
        )
    return out


def make_hazai_records() -> pd.DataFrame:
    pdf = RAW / "PMC3641962_supplementary" / "1471-2105-14-130-S1.pdf"
    rows = parse_hazai_table(pdf)
    s2 = read_sdf_records(RAW / "PMC3641962_supplementary" / "1471-2105-14-130-S2.sdf")
    s3 = read_sdf_records(RAW / "PMC3641962_supplementary" / "1471-2105-14-130-S3.sdf")
    if len(s2) != 164 or len(s3) != 99:
        raise ValueError(f"Hazai SDF counts unexpected: S2={len(s2)}, S3={len(s3)}")
    out: List[Dict[str, object]] = []
    for i, row in enumerate(rows):
        mol = (s2 if i < 164 else s3)[i if i < 164 else i - 164]
        out.append(
            {
                "independent_record_id": f"HAZAI2013_ABCG2_{i + 1:04d}",
                "transporter_gene": "ABCG2",
                "transporter_protein": "human wild-type BCRP",
                "binary_label": int(row["source_section_class"]),
                "strict_binary_label": int(row["source_section_class"]),
                "strict_binary_basis": "explicit_Hazai_table_section_substrate_or_non_substrate_class",
                "source_class": "substrate" if row["source_section_class"] else "non_substrate",
                "compound_name": row["compound_name"],
                "source_dataset_split": row["source_data_set"],
                "source_predicted_property": row["source_predicted_property"],
                "source_reference_ids": row["source_reference_ids"],
                "source_pubchem_cid": mol.get("source_pubchem_cid"),
                "source_sdf_mol_name": mol.get("sdf_mol_name", ""),
                "source_smiles": mol.get("source_smiles_from_sdf", ""),
                "standardized_smiles": mol.get("source_smiles_from_sdf", ""),
                "standardized_inchikey": mol.get("source_inchikey", ""),
                "inchikey_connectivity_block": str(mol.get("source_inchikey", "")).split("-")[0] if mol.get("source_inchikey") else "",
                "preferred_name": row["compound_name"],
                "source_structure_status": mol.get("sdf_parse_status", "failed"),
                "source_pmc": "PMC3641962",
                "source_doi": "10.1186/1471-2105-14-130",
                "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3641962/",
                "evidence_grade": "B2_curated_literature_consensus",
                "label_basis": "table section class; source used literature transport/resistance assays and accepted conflicting calls only when confirmed by >=2 independent studies",
                "negative_evidence_basis": "source_curated_non_substrate_section" if not row["source_section_class"] else "not_applicable",
                "recommended_use": "independent_external_or_sensitivity; primary-assay review before strict training",
                "strict_training_default": "no",
                "uncertainty_note": "curated literature class; row-level direct assay controls and exact tested entity are not exposed in this table",
            }
        )
    return pd.DataFrame(out)


def canonicalize_smiles(smiles: object) -> Tuple[str, str]:
    if not isinstance(smiles, str) or not smiles.strip():
        return "", ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "", ""
        can = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        return can, inchi.MolToInchiKey(mol)
    except Exception:
        return "", ""


def make_lavecchia_records() -> pd.DataFrame:
    xlsx = RAW / "PMC6571636_supplementary" / "molecules-503177-2nd-revised-supple" / "S1.xlsx"
    df = pd.read_excel(xlsx)
    df = df.drop(columns=[c for c in df.columns if str(c).startswith("Unnamed")], errors="ignore")
    required = {"ID", "Name", "Smiles", "Set", "Class", "Prediction target 1", "Prediction target 2", "Prediction target 3"}
    if not required.issubset(df.columns):
        raise ValueError(f"Lavecchia S1 missing columns: {required - set(df.columns)}")
    expected = {"Pgp_Inhibitor": 1178, "substrate": 477, "Pgp_nonInhibitor": 787, "non_substrate": 70}
    if df["Class"].value_counts().to_dict() != expected:
        raise ValueError(f"Unexpected Lavecchia class counts: {df['Class'].value_counts().to_dict()}")
    rows: List[Dict[str, object]] = []
    for _, r in df.iterrows():
        source_class = str(r["Class"])
        if source_class == "substrate":
            binary_label: Optional[int] = 1
            binary_view = "substrate"
            negative_basis = "not_applicable"
        elif source_class in {"non_substrate", "Pgp_nonInhibitor"}:
            binary_label = 0
            binary_view = "non_active"
            negative_basis = "explicit_non_substrate_class" if source_class == "non_substrate" else "source_non_inhibitor_class_merged_into_non_active"
        else:
            binary_label = None
            binary_view = "excluded_inhibitor"
            negative_basis = "not_a_non_substrate; inhibitor class excluded"
        can, ik = canonicalize_smiles(r["Smiles"])
        rows.append(
            {
                "independent_record_id": f"LAVECCHIA2019_ABCB1_{int(r['ID']):04d}",
                "transporter_gene": "ABCB1",
                "transporter_protein": "human P-gp",
                "binary_label": binary_label,
                "strict_binary_label": binary_label if source_class in {"substrate", "non_substrate"} else None,
                "strict_binary_basis": (
                    "explicit_source_substrate_or_non_substrate_class"
                    if source_class in {"substrate", "non_substrate"}
                    else "not_assigned_for_inhibitor_or_non_inhibitor_class"
                ),
                "binary_view": binary_view,
                "source_class": source_class,
                "compound_name": r["Name"],
                "source_iupac_name": r["IUPAC Name"],
                "source_dataset_split": r["Set"],
                "source_prediction_target_1_inhibitor": r["Prediction target 1"],
                "source_prediction_target_2_substrate": r["Prediction target 2"],
                "source_prediction_target_3_nonactive": r["Prediction target 3"],
                "source_record_id": int(r["ID"]),
                "source_smiles": r["Smiles"],
                "standardized_smiles": can,
                "standardized_isomeric_smiles": can,
                "standardized_inchikey": ik,
                "inchikey_connectivity_block": ik.split("-")[0] if ik else "",
                "preferred_name": r["Name"],
                "source_structure_status": "ok" if can and ik else "failed",
                "source_pmc": "PMC6571636",
                "source_doi": "10.3390/molecules24102006",
                "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6571636/",
                "evidence_grade": "B2_curated_multiclass_literature_dataset",
                "label_basis": "source Class; source merged Pgp_nonInhibitor and non_substrate into non-active for its classifier",
                "negative_evidence_basis": negative_basis,
                "recommended_use": "binary_view sensitivity/external holdout; exclude inhibitor class; primary-assay review before strict training",
                "strict_training_default": "no",
                "uncertainty_note": "class is curated from heterogeneous P-gp assays; inhibitors are not negatives and are retained only as excluded records",
            }
        )
    return pd.DataFrame(rows)


def build_overlap_audit(independent: pd.DataFrame) -> pd.DataFrame:
    daood_path = PROCESSED / "daood_training_records_standardized_v1.csv"
    daood = pd.read_csv(daood_path, low_memory=False)
    daood["connectivity_block"] = daood["standardized_inchikey"].fillna("").astype(str).str.split("-").str[0]
    # Exact stereochemical InChIKey is preferred; connectivity block is reported
    # separately because salt/tautomer/stereo choices can differ across sources.
    rows = []
    for _, r in independent.iterrows():
        gene = r["transporter_gene"]
        key = str(r.get("standardized_inchikey") or "")
        conn = key.split("-")[0] if key else ""
        q = daood[daood["transporter_gene"].astype(str) == str(gene)]
        exact = q[q["standardized_inchikey"].astype(str) == key] if key else q.iloc[0:0]
        conn_hits = q[q["connectivity_block"].astype(str) == conn] if conn else q.iloc[0:0]
        hit = exact if len(exact) else conn_hits
        if len(exact): overlap_type = "exact_inchikey"
        elif len(conn_hits): overlap_type = "connectivity_block_only"
        else: overlap_type = "no_daood_overlap"
        labels = sorted({int(x) for x in hit["published_binary_label"].dropna().tolist()})
        src_label = r.get("binary_label")
        if pd.isna(src_label) or src_label is None or not labels:
            agreement = "not_comparable_or_no_hit"
        elif len(labels) > 1:
            agreement = "daood_internal_label_conflict"
        elif int(src_label) == labels[0]:
            agreement = "agree"
        else:
            agreement = "disagree"
        rows.append(
            {
                "independent_record_id": r["independent_record_id"],
                "transporter_gene": gene,
                "independent_binary_label": src_label,
                "independent_source_class": r.get("source_class", ""),
                "independent_compound_name": r.get("compound_name", ""),
                "standardized_inchikey": key,
                "overlap_type": overlap_type,
                "daood_match_count": len(hit),
                "daood_record_ids": " | ".join(hit["training_record_id"].astype(str).tolist()),
                "daood_preferred_names": " | ".join(hit["preferred_name"].astype(str).tolist()),
                "daood_labels": " | ".join(map(str, labels)),
                "cross_source_label_agreement": agreement,
                "action": (
                    "hold_for_conflict_review"
                    if agreement in {"disagree", "daood_internal_label_conflict"}
                    else "holdout_or_sensitivity_only_until_primary_assay_review"
                ),
            }
        )
    return pd.DataFrame(rows)


def build_binary_candidate_view(independent: pd.DataFrame, audit: pd.DataFrame) -> pd.DataFrame:
    """Conservative sensitivity/external view; never silently calls it truth."""
    cols = ["independent_record_id", "overlap_type", "cross_source_label_agreement", "action"]
    out = independent.merge(audit[cols], on="independent_record_id", how="left")
    keep = (
        out["binary_label"].notna()
        & (out["source_structure_status"].astype(str) == "ok")
        & (out["overlap_type"] == "no_daood_overlap")
    )
    out = out.loc[keep].copy()
    out["candidate_use"] = "independent_external_or_sensitivity_only"
    out["strict_training_default"] = "no"
    out["candidate_reason"] = "parseable structure + derived binary label + no Daood exact/connectivity overlap; primary assay review still required"
    return out


def build_strict_binary_candidate_view(independent: pd.DataFrame, audit: pd.DataFrame) -> pd.DataFrame:
    """Use only explicit substrate/non-substrate source classes."""
    cols = ["independent_record_id", "overlap_type", "cross_source_label_agreement", "action"]
    out = independent.merge(audit[cols], on="independent_record_id", how="left")
    keep = (
        out["strict_binary_label"].notna()
        & (out["source_structure_status"].astype(str) == "ok")
        & (out["overlap_type"] == "no_daood_overlap")
    )
    out = out.loc[keep].copy()
    out["candidate_use"] = "strict_external_or_sensitivity_only"
    out["strict_training_default"] = "no"
    out["candidate_reason"] = "explicit source substrate/non-substrate class + parseable structure + no Daood exact/connectivity overlap; primary assay review still required"
    return out


def build_eleventh_batch_queue_intersection(independent: pd.DataFrame) -> pd.DataFrame:
    """Map new-source structures onto the 512-row ChEMBL review queue."""
    triage_path = PROCESSED / "chembl_candidate_provenance_triage_v1.csv"
    triage = pd.read_csv(triage_path, low_memory=False)
    triage["connectivity_block"] = triage["standardized_inchikey"].fillna("").astype(str).str.split("-").str[0]
    triage = triage[triage["connectivity_block"].astype(bool)].copy()
    triage = triage.drop_duplicates(["transporter_gene", "connectivity_block", "audit_record_id"])
    ind = independent.copy()
    ind["connectivity_block"] = ind["standardized_inchikey"].fillna("").astype(str).str.split("-").str[0]
    ind = ind[ind["connectivity_block"].astype(bool)].copy()
    keep_cols = [
        "audit_record_id", "transporter_gene", "preferred_names", "published_binary_label",
        "standardized_inchikey", "connectivity_block", "label_gap_priority",
    ]
    q = ind.merge(triage[keep_cols], on=["transporter_gene", "connectivity_block"], suffixes=("_independent", "_queue"))
    if q.empty:
        return pd.DataFrame(columns=["audit_record_id"])
    q["match_level"] = q.apply(
        lambda r: "exact_inchikey" if str(r["standardized_inchikey_independent"]) == str(r["standardized_inchikey_queue"]) else "connectivity_block_only",
        axis=1,
    )
    def agreement(r: pd.Series) -> str:
        if pd.isna(r.get("binary_label")):
            return "independent_label_not_binary"
        if pd.isna(r.get("published_binary_label")):
            return "queue_label_missing"
        return "agree" if int(r["binary_label"]) == int(r["published_binary_label"]) else "disagree"
    q["queue_label_agreement"] = q.apply(agreement, axis=1)
    q["action"] = q["queue_label_agreement"].map(
        {"disagree": "hold_for_primary_assay_conflict_review", "agree": "independent_support_but_not_Tier_A", "independent_label_not_binary": "no_binary_support", "queue_label_missing": "manual_queue_label_review"}
    ).fillna("manual_review")
    cols = [
        "audit_record_id", "transporter_gene", "preferred_names", "published_binary_label",
        "label_gap_priority", "independent_record_id", "compound_name", "source_class",
        "binary_label", "strict_binary_label", "evidence_grade", "source_pmc", "source_doi",
        "match_level", "queue_label_agreement", "action",
    ]
    return q[cols].sort_values(["transporter_gene", "audit_record_id", "independent_record_id"]).reset_index(drop=True)


def source_registry() -> pd.DataFrame:
    files = [
        ("SRC_B12_HAZAI_XML", "PMC3641962", "PMC3641962_fullText.xml", "full_text_xml", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3641962/fullTextXML"),
        ("SRC_B12_HAZAI_ARCHIVE", "PMC3641962", "PMC3641962_supplementaryFiles.zip", "supplementary_archive", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3641962/supplementaryFiles"),
        ("SRC_B12_HAZAI_TABLE", "PMC3641962", "PMC3641962_supplementary/1471-2105-14-130-S1.pdf", "supplementary_table_pdf", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3641962/bin/1471-2105-14-130-S1.pdf"),
        ("SRC_B12_HAZAI_S2", "PMC3641962", "PMC3641962_supplementary/1471-2105-14-130-S2.sdf", "substrate_structure_sdf", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3641962/bin/1471-2105-14-130-S2.sdf"),
        ("SRC_B12_HAZAI_S3", "PMC3641962", "PMC3641962_supplementary/1471-2105-14-130-S3.sdf", "nonsubstrate_structure_sdf", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3641962/bin/1471-2105-14-130-S3.sdf"),
        ("SRC_B12_LAVECCHIA_XML", "PMC6571636", "PMC6571636_fullText.xml", "full_text_xml", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6571636/fullTextXML"),
        ("SRC_B12_LAVECCHIA_ARCHIVE", "PMC6571636", "PMC6571636_supplementaryFiles.zip", "supplementary_archive", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6571636/supplementaryFiles"),
        ("SRC_B12_LAVECCHIA_S1", "PMC6571636", "PMC6571636_supplementary/molecules-503177-2nd-revised-supple/S1.xlsx", "multiclass_label_table", "https://pmc.ncbi.nlm.nih.gov/articles/PMC6571636/bin/molecules-24-02006-s001.zip"),
    ]
    rows = []
    for source_id, paper_pmc, rel, typ, url in files:
        p = RAW / rel
        rows.append(
            {
                "source_id": source_id,
                "paper_pmc": paper_pmc,
                "relative_path": f"raw/public_independent/{rel}",
                "source_url": url,
                "source_type": typ,
                "acquired_at": ACQUIRED_AT,
                "exists": p.exists(),
                "size_bytes": p.stat().st_size if p.exists() else None,
                "sha256": sha256(p) if p.exists() else "",
                "allowed_use": "independent_source_candidate; binary sensitivity/external holdout; not automatic Tier-A truth",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true", help="download Europe PMC archives/XML if absent")
    args = parser.parse_args()
    ensure_sources(download=args.download)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    hazai = make_hazai_records()
    lavecchia = make_lavecchia_records()
    hazai.to_csv(PROCESSED / "public_independent_hazai_2013_abcg2_labels_v1.csv", index=False, encoding="utf-8-sig")
    lavecchia.to_csv(PROCESSED / "public_independent_lavecchia_2019_abcb1_multiclass_v1.csv", index=False, encoding="utf-8-sig")
    independent = pd.concat([hazai, lavecchia], ignore_index=True, sort=False)
    independent.to_csv(PROCESSED / "public_independent_transporter_labels_long_v1.csv", index=False, encoding="utf-8-sig")
    audit = build_overlap_audit(independent)
    audit.to_csv(PROCESSED / "public_independent_cross_source_overlap_audit_v1.csv", index=False, encoding="utf-8-sig")
    candidates = build_binary_candidate_view(independent, audit)
    candidates.to_csv(PROCESSED / "public_independent_binary_candidate_view_v1.csv", index=False, encoding="utf-8-sig")
    strict_candidates = build_strict_binary_candidate_view(independent, audit)
    strict_candidates.to_csv(PROCESSED / "public_independent_strict_binary_candidate_view_v1.csv", index=False, encoding="utf-8-sig")
    queue_intersection = build_eleventh_batch_queue_intersection(independent)
    queue_intersection.to_csv(PROCESSED / "public_independent_eleventh_queue_intersection_v1.csv", index=False, encoding="utf-8-sig")
    registry = source_registry()
    registry.to_csv(PROCESSED / "public_independent_source_registry_v1.csv", index=False, encoding="utf-8-sig")
    summary = {
        "batch": "public_independent_sources_v1",
        "acquired_at": ACQUIRED_AT,
        "sources": {
            "hazai_2013_abcg2": {
                "pmc": "PMC3641962",
                "doi": "10.1186/1471-2105-14-130",
                "rows": int(len(hazai)),
                "substrates": int((hazai.binary_label == 1).sum()),
                "non_substrates": int((hazai.binary_label == 0).sum()),
                "structure_parse_ok": int((hazai.source_structure_status == "ok").sum()),
            },
            "lavecchia_2019_abcb1": {
                "pmc": "PMC6571636",
                "doi": "10.3390/molecules24102006",
                "rows_all_classes": int(len(lavecchia)),
                "substrates": int((lavecchia.source_class == "substrate").sum()),
                "non_active_binary_negatives": int(lavecchia.binary_label.eq(0).sum()),
                "inhibitors_excluded": int((lavecchia.source_class == "Pgp_Inhibitor").sum()),
                "binary_view_rows": int(lavecchia.binary_label.notna().sum()),
                "structure_parse_ok_all_classes": int((lavecchia.source_structure_status == "ok").sum()),
                "structure_parse_ok_binary_view": int(((lavecchia.binary_label.notna()) & (lavecchia.source_structure_status == "ok")).sum()),
            },
        },
        "cross_source_overlap": {
            "rows_audited": int(len(audit)),
            "exact_inchikey": int((audit.overlap_type == "exact_inchikey").sum()),
            "connectivity_block_only": int((audit.overlap_type == "connectivity_block_only").sum()),
            "no_daood_overlap": int((audit.overlap_type == "no_daood_overlap").sum()),
            "agreements": int((audit.cross_source_label_agreement == "agree").sum()),
            "disagreements": int((audit.cross_source_label_agreement == "disagree").sum()),
            "conflicts": int((audit.cross_source_label_agreement == "daood_internal_label_conflict").sum()),
        },
        "conservative_binary_candidate_view": {
            "rows": int(len(candidates)),
            "by_gene_and_label": {
                f"{gene}_{int(label)}": int(len(g))
                for (gene, label), g in candidates.groupby(["transporter_gene", "binary_label"])
            },
        },
        "strict_explicit_class_candidate_view": {
            "rows": int(len(strict_candidates)),
            "by_gene_and_label": {
                f"{gene}_{int(label)}": int(len(g))
                for (gene, label), g in strict_candidates.groupby(["transporter_gene", "strict_binary_label"])
            },
        },
        "eleventh_batch_queue_intersection": {
            "rows": int(len(queue_intersection)),
            "unique_queue_records": int(queue_intersection["audit_record_id"].nunique()) if len(queue_intersection) else 0,
            "by_gene_unique_queue_records": {
                gene: int(g["audit_record_id"].nunique()) for gene, g in queue_intersection.groupby("transporter_gene")
            } if len(queue_intersection) else {},
            "label_agreement_rows": int((queue_intersection["queue_label_agreement"] == "agree").sum()) if len(queue_intersection) else 0,
            "label_disagreement_rows": int((queue_intersection["queue_label_agreement"] == "disagree").sum()) if len(queue_intersection) else 0,
        },
        "claim_boundary": "independent curated literature datasets; retain source classes, exclude inhibitor-as-negative, and require primary-assay review before strict binary training",
    }
    (PROCESSED / "public_independent_sources_summary_v1.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

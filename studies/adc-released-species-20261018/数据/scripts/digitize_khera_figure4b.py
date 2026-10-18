#!/usr/bin/env python3
"""Reproducibly digitize the mean PAMPA bars in Khera et al. Figure 4B.

The source is page 6 of the downloaded article PDF rendered at 220 dpi to
``tmp/pdfs/eighth_batch/khera-06.png``.  Only bar means are extracted.  Error
bars are deliberately not converted to numerical uncertainty.  The two axis
calibrations below are algorithmic checks on the same image, not independent
curation, so every output row remains exploratory and truth-set ineligible.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
IMAGE = ROOT.parent / "tmp" / "pdfs" / "eighth_batch" / "khera-06.png"
OUTPUT = ROOT / "processed" / "khera_figure4b_pampa_digitization_v1.csv"


@dataclass(frozen=True)
class Bar:
    analyte_as_figure: str
    normalized_entity_name: str
    release_species_id: str
    payload_class: str
    x1: int
    x2: int


BARS = [
    Bar("MMAE", "MMAE", "ADC_REL_007", "microtubule_inhibitor", 881, 900),
    Bar("Ahx-May", "Ahx-maytansine", "", "microtubule_inhibitor", 919, 937),
    Bar("Lys-MCC-DM1*", "Lys-MCC-DM1", "ADC_REL_001", "microtubule_inhibitor", 955, 974),
    Bar("CaliD", "calicheamicin D", "", "DNA_interacting_agent", 992, 1011),
    Bar("PBD", "SG3199", "ADC_REL_013", "DNA_interacting_agent", 1030, 1048),
    Bar("CaliG", "calicheamicin G", "", "DNA_interacting_agent", 1067, 1085),
    Bar("Exatecan", "Exatecan", "ADC_REL_011", "topoisomerase_I_inhibitor", 1104, 1123),
    Bar("DXd", "DXd", "ADC_REL_010", "topoisomerase_I_inhibitor", 1142, 1160),
    Bar("SN38", "SN-38", "ADC_REL_009", "topoisomerase_I_inhibitor", 1179, 1197),
]

# Tick centers detected from the Figure 4B y-axis in the full-page 220 dpi PNG.
# They span 3.0 to 0.0 in steps of 0.5, even though the 3.0 label is near the
# upper chart boundary.
TICK_Y = np.array([1700.0, 1740.5, 1781.5, 1821.5, 1862.5, 1903.5, 1944.5])
TICK_VALUE = np.array([3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0])


def first_full_width_dark_row(gray: np.ndarray, x1: int, x2: int) -> int:
    """Detect the top horizontal bar border while rejecting error-bar caps."""
    interior = gray[1700:1945, x1 + 1 : x2]
    dark_fraction = (interior < 120).mean(axis=1)
    candidates = np.flatnonzero(dark_fraction >= 0.90)
    if candidates.size == 0:
        raise ValueError(f"No bar top found for x={x1}:{x2}")
    return int(candidates[0] + 1700)


def digitize() -> list[dict[str, object]]:
    gray = np.asarray(Image.open(IMAGE).convert("L"))
    if gray.shape != (2566, 1961):
        raise ValueError(f"Unexpected rendered page dimensions: {gray.shape}")

    # Pass 1: least-squares calibration across all seven detected tick centers.
    slope_1, intercept_1 = np.polyfit(TICK_Y, TICK_VALUE, 1)
    # Pass 2: endpoint-only calibration using 0.0 and 3.0 ticks.
    slope_2 = (3.0 - 0.0) / (1700.0 - 1944.5)
    intercept_2 = -slope_2 * 1944.5

    rows: list[dict[str, object]] = []
    for index, bar in enumerate(BARS, start=1):
        top_first = first_full_width_dark_row(gray, bar.x1, bar.x2)
        top_center = top_first + 0.5
        value_1 = slope_1 * top_center + intercept_1
        value_2 = slope_2 * top_center + intercept_2
        consensus = (value_1 + value_2) / 2.0
        rows.append(
            {
                "digitization_record_id": f"KHERA_F4B_{index:02d}",
                "analyte_as_figure": bar.analyte_as_figure,
                "normalized_entity_name": bar.normalized_entity_name,
                "release_species_id": bar.release_species_id,
                "payload_class": bar.payload_class,
                "assay": "Gentest pre-coated PAMPA",
                "working_concentration_value": 200,
                "working_concentration_unit": "uM",
                "donor_solution": "20% methanol in PBS",
                "donor_volume_uL": 300,
                "acceptor_solution": "20% methanol in PBS",
                "acceptor_volume_uL": 200,
                "incubation_time": "2.5 h or 5 h; analyte-specific time not reported",
                "analysis_method": "C18 HPLC against standard curve",
                "x_interval_pixels": f"{bar.x1}:{bar.x2}",
                "bar_top_first_row_pixel": top_first,
                "bar_top_center_y_pixel": round(top_center, 3),
                "axis_pass1_slope": round(float(slope_1), 10),
                "axis_pass1_intercept": round(float(intercept_1), 10),
                "pampa_pass1_1e6_cm_s": round(float(value_1), 4),
                "axis_pass2_slope": round(float(slope_2), 10),
                "axis_pass2_intercept": round(float(intercept_2), 10),
                "pampa_pass2_1e6_cm_s": round(float(value_2), 4),
                "pampa_consensus_1e6_cm_s": round(float(consensus), 4),
                "one_pixel_value_uncertainty_1e6_cm_s": round(abs(float(slope_1)), 4),
                "pass_difference_1e6_cm_s": round(abs(float(value_1 - value_2)), 6),
                "error_bar_numeric_value": "",
                "numeric_origin": "algorithmic_digitization_from_220dpi_pdf_render_bar_mean",
                "digitization_status": "exploratory_single_figure_algorithmic_no_independent_curator",
                "truth_set_eligible": "no",
                "model_binary_label": "",
                "source_id": "SRC_B8_001",
                "source_location": "Figure 4B; article page 315; PDF physical page 6",
                "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8832507/",
                "note": "Mean bar only; error bars not digitized. Figure star on Lys-MCC-DM1 points to Ogitani et al. 2016. Figure label PBD is mapped to SG3199-PBD from the Methods.",
            }
        )
    return rows


def main() -> None:
    rows = digitize()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows -> {OUTPUT}")
    for row in rows:
        print(row["normalized_entity_name"], row["pampa_consensus_1e6_cm_s"])


if __name__ == "__main__":
    main()

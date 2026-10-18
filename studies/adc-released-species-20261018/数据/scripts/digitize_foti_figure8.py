#!/usr/bin/env python3
"""Reproducibly digitize mean efflux-ratio bars from Foti et al. Figure 8.

The source figure is a JPEG with a linear y axis.  This script deliberately
extracts bar means only; error bars are not converted into numerical SD values.
Two algorithmic calibration passes are reported, but they are not a substitute
for a second independent curator.  Consequently, the output remains
exploratory and is not eligible for the model truth set.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
IMAGE = ROOT / "raw" / "second_batch_release_species" / "Foti_2015_Figure8_LysMCCDM1_Pgp.jpg"
OUTPUT = ROOT / "processed" / "foti_figure8_digitization_v1.csv"


@dataclass(frozen=True)
class Bar:
    group: str
    analyte: str
    concentration_value: str
    concentration_unit: str
    inhibitor: str
    cell_system: str
    x1: int
    x2: int


GROUPS = [
    ("14C-Mann", "14C-mannitol", "", "", ""),
    ("3H-CsA", "3H-cyclosporin A", "", "", ""),
    ("3H-CsA + GG", "3H-cyclosporin A", "", "", "elacridar; concentration not reported in caption"),
    ("LMD - 0.1 uM", "Lys-MCC-DM1", "0.1", "uM", ""),
    ("LMD - 1 uM", "Lys-MCC-DM1", "1", "uM", ""),
    ("LMD - 5 uM", "Lys-MCC-DM1", "5", "uM", ""),
    ("LMD - 10 uM", "Lys-MCC-DM1", "10", "uM", ""),
    ("LMD - 5 uM + GG", "Lys-MCC-DM1", "5", "uM", "elacridar; concentration not reported in caption"),
]

X_INTERVALS = [
    ((158, 180), (183, 205)),
    ((227, 249), (252, 274)),
    ((297, 319), (322, 344)),
    ((366, 388), (392, 414)),
    ((436, 458), (461, 483)),
    ((506, 528), (531, 553)),
    ((575, 597), (600, 622)),
    ((645, 667), (670, 692)),
]

# Tick centers were detected from the 13-pixel horizontal tick marks at x=130:149.
TICK_Y = np.array([68.5, 130.5, 192.5, 254.5, 316.5, 378.5, 439.5], dtype=float)
TICK_VALUE = np.array([12, 10, 8, 6, 4, 2, 0], dtype=float)


def build_bars() -> list[Bar]:
    bars: list[Bar] = []
    for (group, analyte, concentration_value, concentration_unit, inhibitor), pair in zip(GROUPS, X_INTERVALS):
        for cell_system, (x1, x2) in zip(("LLC-PK1 parental", "LLC-PK1 human P-gp"), pair):
            bars.append(
                Bar(
                    group=group,
                    analyte=analyte,
                    concentration_value=concentration_value,
                    concentration_unit=concentration_unit,
                    inhibitor=inhibitor,
                    cell_system=cell_system,
                    x1=x1,
                    x2=x2,
                )
            )
    return bars


def first_full_width_dark_row(gray: np.ndarray, x1: int, x2: int) -> int:
    """Find the top horizontal bar edge while rejecting narrower error-bar caps."""
    interior = gray[50:442, x1 + 2 : x2 - 1]
    fraction_dark = (interior < 80).mean(axis=1)
    candidates = np.flatnonzero(fraction_dark >= 0.90)
    if candidates.size == 0:
        raise ValueError(f"No full-width bar edge detected for x={x1}:{x2}")
    return int(candidates[0] + 50)


def digitize() -> list[dict[str, object]]:
    gray = np.asarray(Image.open(IMAGE).convert("L"))
    if gray.shape != (719, 708):
        raise ValueError(f"Unexpected source image dimensions: {gray.shape}")

    # Pass 1: least-squares calibration across all seven y-axis ticks.
    slope_1, intercept_1 = np.polyfit(TICK_Y, TICK_VALUE, 1)
    # Pass 2: independent endpoint calibration using only 0 and 12 ticks.
    slope_2 = (12.0 - 0.0) / (68.5 - 439.5)
    intercept_2 = -slope_2 * 439.5

    rows: list[dict[str, object]] = []
    for index, bar in enumerate(build_bars(), start=1):
        top_first_row = first_full_width_dark_row(gray, bar.x1, bar.x2)
        top_center_y = top_first_row + 0.5
        pass_1 = slope_1 * top_center_y + intercept_1
        pass_2 = slope_2 * top_center_y + intercept_2
        consensus = (pass_1 + pass_2) / 2.0
        pixel_uncertainty = abs(slope_1)  # +/- one source pixel
        rows.append(
            {
                "digitization_record_id": f"FOTI_F8_{index:02d}",
                "figure_group": bar.group,
                "analyte": bar.analyte,
                "concentration_value": bar.concentration_value,
                "concentration_unit": bar.concentration_unit,
                "inhibitor": bar.inhibitor,
                "cell_system": bar.cell_system,
                "transporter_gene": "ABCB1" if "human P-gp" in bar.cell_system else "",
                "x_interval_pixels": f"{bar.x1}:{bar.x2}",
                "bar_top_first_row_pixel": top_first_row,
                "bar_top_center_y_pixel": round(top_center_y, 3),
                "axis_pass1_slope": round(float(slope_1), 9),
                "axis_pass1_intercept": round(float(intercept_1), 9),
                "efflux_ratio_pass1": round(float(pass_1), 4),
                "axis_pass2_slope": round(float(slope_2), 9),
                "axis_pass2_intercept": round(float(intercept_2), 9),
                "efflux_ratio_pass2": round(float(pass_2), 4),
                "efflux_ratio_consensus": round(float(consensus), 4),
                "one_pixel_value_uncertainty": round(float(pixel_uncertainty), 4),
                "pass_difference": round(float(abs(pass_1 - pass_2)), 6),
                "numeric_origin": "algorithmic_digitization_from_jpeg_bar_mean",
                "digitization_status": "exploratory_single_figure_algorithmic_no_independent_curator",
                "truth_set_eligible": "no",
                "source_id": "SRC_B5_002",
                "source_location": "Figure 8",
                "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11024933/",
                "note": "Mean bar only. Error bars were not digitized. Two calibration passes are algorithmic, not two independent curators.",
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


if __name__ == "__main__":
    main()

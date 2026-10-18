"""Reproduce current ADC Figure 5 and current ADC table from the frozen 14-species CSV.

Run from any directory with an environment containing matplotlib, numpy, Pillow,
and RDKit: python make_figure5.py. Source structures, scores and decision reasons
are read together; no historical example, inferred transport label, or hand-edited
score is used. The figures retain the original Songti SC font and project palette.
"""
from pathlib import Path
import csv
import hashlib
import io
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import palette

DATA = ROOT / "outputs/20261012_adc_completion/scoring/species_predictions.csv"
OUT = ROOT / "outputs/20261018_manuscript_cleanup/figure5"
FIGURES = ROOT / "论文主稿/latex_jcim/figures"
TABLE = ROOT / "论文主稿/latex_jcim/tables/table5_adc.tex"
C = palette.C
S = (9.5, 8.5, 7.5)
BLUE, GREY = C["blue"], C["grey"]
DMH = "N-acetyl-gamma-calicheamicin dimethylhydrazide"
CARD_NAMES = [("MMAE", "cys-VC-MMAE", "DXd"),
              ("DM4", "S-methyl-DM4", "Lys-MCC-DM1")]
ORDER = ["MMAE", "cys-VC-MMAE", "MMAF", "cys-mcMMAF", "DM1", "Lys-MCC-DM1",
         "DM4", "S-methyl-DM4", "Lys-SPDB-DM4", "SN-38", "DXd", "Exatecan", "SG3199", DMH]
REASONS = {
    frozenset(): "None",
    frozenset(["probability_gray_zone"]): "Score gray zone",
    frozenset(["outside_both_reference_panels"]): "Reference support",
    frozenset(["outside_physicochemical_envelope"]): "Physicochemical support",
    frozenset(["outside_both_reference_panels", "outside_physicochemical_envelope"]): "Both support criteria",
}


def display_name(name):
    if name == "Lys-MCC-DM1":
        return "Lys-MCC-DM1*"
    return "Calicheamicin DMH" if name == DMH else name


def family(name):
    if name in ORDER[:4]:
        return "Auristatin"
    if name in ORDER[4:9]:
        return "Maytansinoid"
    if name in ORDER[9:12]:
        return "TOP1 inhibitor"
    return "PBD" if name == "SG3199" else "Calicheamicin"


def reason_color(row):
    rs = row["reasons"]
    if not rs:
        return BLUE
    if rs == {"probability_gray_zone"}:
        return GREY
    if len(rs) == 2:
        return C["purple"]
    if "outside_physicochemical_envelope" in rs:
        return C["red"]
    return C["green"]


def load_and_verify():
    rows = list(csv.DictReader(DATA.open(encoding="utf-8-sig")))
    assert len(rows) == 14 and {r["entity_name"] for r in rows} == set(ORDER)
    rows = {r["entity_name"]: r for r in rows}
    checks = []
    for name, row in rows.items():
        mol = Chem.MolFromSmiles(row["source_smiles"])
        assert mol is not None, name
        key = Chem.MolToInchiKey(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        assert key == row["source_inchikey"], (name, key, row["source_inchikey"])
        assert formula == row["source_formula"], (name, formula, row["source_formula"])
        row["mol"] = mol
        row["score"] = float(row["primary_calibrated_probability"])
        row["low"] = float(row["bootstrap_probability_p10"])
        row["high"] = float(row["bootstrap_probability_p90"])
        row["z"] = float(row["physchem_max_abs_robust_z"])
        row["reasons"] = frozenset(json.loads(row["abstention_reasons_json"]))
        row["accepted"] = row["accepted_prediction"] == "yes"
        assert row["reasons"] in REASONS
        assert 0 <= row["low"] <= row["high"] <= 1
        assert float(row["physchem_threshold_q99"]) == 7
        assert row["accepted"] == (row["risk_call"] == "high_ABCB1_efflux_risk")
        assert (not row["reasons"]) == row["accepted"]
        checks.append({"name": name, "source_inchikey": key, "source_formula": formula,
                       "sulfur_atoms": sum(a.GetAtomicNum() == 16 for a in mol.GetAtoms()),
                       "primary_score": row["score"], "direct_score": float(row["direct_panel_probability"]),
                       "p10": row["low"], "p90": row["high"], "robust_z_max": row["z"],
                       "decision": row["risk_call"], "reasons": sorted(row["reasons"])})
    assert sum(r["accepted"] for r in rows.values()) == 5
    assert rows["MMAF"]["reasons"] == {"probability_gray_zone"}
    assert sum(a.GetAtomicNum() == 16 for a in rows["S-methyl-DM4"]["mol"].GetAtoms()) == 1
    return rows, checks


def molecule_image(row):
    drawer = rdMolDraw2D.MolDraw2DCairo(1200, 430)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.8
    opts.padding = 0.06
    opts.fixedFontSize = 23
    opts.addStereoAnnotation = False
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, row["mol"])
    drawer.FinishDrawing()
    image_bytes = drawer.GetDrawingText()
    (OUT / (row["entity_name"] + "_source_structure.png")).write_bytes(image_bytes)
    (HERE / "fig5_structures" / (row["entity_name"] + ".png")).write_bytes(image_bytes)
    return mpimg.imread(io.BytesIO(image_bytes), format="png")


def card(ax, row):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.add_patch(Rectangle((0.006, 0.01), 0.988, 0.98, transform=ax.transAxes,
                           facecolor="none", edgecolor="#777777", lw=0.6))
    img = molecule_image(row)
    ax.figure.canvas.draw()
    bbox = ax.get_window_extent(ax.figure.canvas.get_renderer())
    box_ratio = bbox.width / bbox.height
    image_ratio = img.shape[1] / img.shape[0]
    width, height = 0.95, 0.64
    if image_ratio > width / height * box_ratio:
        height = width * box_ratio / image_ratio
    else:
        width = height * image_ratio / box_ratio
    ax.imshow(img, extent=(0.5 - width / 2, 0.5 + width / 2,
                           0.49 - height / 2, 0.49 + height / 2), aspect="auto")
    ax.text(0.05, 0.955, display_name(row["entity_name"]), fontsize=S[1], va="top", ha="left",
            color="#1a1a1a", fontweight=700)
    ax.text(0.5, 0.045, "High efflux risk" if row["accepted"] else "Classification withheld",
            fontsize=S[2], va="bottom", ha="center", color="#444444")


def strip_panel(ax, rows):
    for y, row in enumerate(rows):
        col = BLUE if row["accepted"] else GREY
        tint = palette.tint("blue" if row["accepted"] else "grey", 0.55)
        ax.plot([row["low"], row["high"]], [y, y], color=tint, lw=1.5, solid_capstyle="butt")
        ax.scatter(row["score"], y, s=25, facecolor=col if row["accepted"] else "white",
                   edgecolor=col, linewidths=0.9, zorder=3)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([display_name(r["entity_name"]) for r in rows], fontsize=S[2])
    for tick, row in zip(ax.get_yticklabels(), rows):
        tick.set_color(BLUE if row["accepted"] else GREY)
    ax.set_xlim(0.54, 1.035); ax.set_xticks([0.6, 0.8, 1.0])
    ax.set_xlabel("ABCB1 model score", fontsize=S[1], labelpad=3)
    ax.set_title("Scores and P10–P90 ranges", fontsize=S[0], loc="left", pad=8)
    ax.tick_params(labelsize=S[2]); ax.margins(y=0.035)


def support_panel(ax, rows):
    ax.axvspan(0, 7, color=palette.tint("grey", 0.13), zorder=0)
    ax.axvline(7, color=GREY, lw=0.9, ls=(0, (3, 2)))
    for y, row in enumerate(rows):
        col = reason_color(row)
        ax.plot([0, row["z"]], [y, y], color="#dedede", lw=0.6, zorder=1)
        ax.scatter(row["z"], y, s=25, facecolor=col if row["accepted"] else "white",
                   edgecolor=col, linewidths=0.9, zorder=3)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([display_name(r["entity_name"]) for r in rows], fontsize=S[2])
    for tick, row in zip(ax.get_yticklabels(), rows):
        tick.set_color(reason_color(row))
    ax.set_xlim(0, 14); ax.set_xticks([0, 4, 8, 12])
    ax.set_xlabel("Maximum absolute robust z", fontsize=S[1], labelpad=3)
    ax.set_title("Physicochemical support", fontsize=S[0], loc="left", pad=8)
    ax.tick_params(labelsize=S[2]); ax.margins(y=0.035)
    ax.annotate("Cutoff: 7", xy=(7, 1.7), xytext=(5, 0), textcoords="offset points",
                fontsize=S[2], color="#555555", ha="left", va="center")


def key(fig, x, y, color, label, filled=False):
    fig.patches.append(Rectangle((x, y - 0.0045), 0.012, 0.009,
                                facecolor=color if filled else "none", edgecolor=color,
                                lw=0.9, transform=fig.transFigure))
    fig.text(x + 0.018, y, label, fontsize=S[2], color="#444444", va="center")


def draw(rows):
    fig = plt.figure(figsize=(6.5, 6.3))
    grid = fig.add_gridspec(2, 3, left=0.055, right=0.975, top=0.950,
                            bottom=0.615, hspace=0.12, wspace=0.16)
    for ri, names in enumerate(CARD_NAMES):
        for ci, name in enumerate(names):
            card(fig.add_subplot(grid[ri, ci]), rows[name])
    ordered = sorted(rows.values(), key=lambda r: r["score"])
    grid2 = fig.add_gridspec(1, 2, left=0.190, right=0.975, top=0.455,
                             bottom=0.105, wspace=0.80)
    ax_b, ax_c = fig.add_subplot(grid2[0]), fig.add_subplot(grid2[1])
    strip_panel(ax_b, ordered); support_panel(ax_c, ordered)
    key(fig, 0.190, 0.579, BLUE, "Accepted high", True)
    key(fig, 0.420, 0.579, GREY, "Abstained")
    fig.text(0.650, 0.579, "14 structures: 5 high, 9 abstain", fontsize=S[2], color="#555555", va="center")
    fig.text(0.055, 0.545, "Abstention criteria (c):", fontsize=S[2], color="#444444", va="center")
    key(fig, 0.300, 0.545, C["green"], "Reference support")
    key(fig, 0.550, 0.545, C["red"], "Physicochemical support")
    key(fig, 0.300, 0.514, C["purple"], "Both support criteria")
    key(fig, 0.550, 0.514, GREY, "Score gray zone")
    fig.text(0.055, 0.963, "a", fontsize=12, fontweight=700, va="bottom", color="#1a1a1a")
    for ax, label in ((ax_b, "b"), (ax_c, "c")):
        ax.text(-0.48, 1.04, label, transform=ax.transAxes, fontsize=12,
                fontweight=700, va="bottom", color="#1a1a1a")
    return fig


def write_table(rows):
    lines = [r"% Generated by figures_src/make_figure5.py from the frozen 14-species CSV.",
             r"\begin{table}[htbp]", r"\centering", r"\fontsize{9}{11}\selectfont",
             r"\begin{threeparttable}",
             r"\caption{Current scores and decisions for 14 ADC-related small-molecule structures.}",
             r"\label{tab:adc}",
             r"\setlength{\tabcolsep}{4pt}",
             r"\begin{tabular}{@{}p{0.22\textwidth}p{0.14\textwidth}rrp{0.11\textwidth}p{0.25\textwidth}@{}}",
             r"\toprule", r"Species & Family & $p_{0.1}$ & $p_{0.3}$ & Decision & Abstention reason \\", r"\midrule"]
    for name in ORDER:
        row = rows[name]
        label = display_name(name)
        reason = REASONS[row["reasons"]] if row["reasons"] else "--"
        lines.append(f"{label} & {family(name)} & {row['score']:.3f} & {float(row['direct_panel_probability']):.3f} & "
                     + ("High risk" if row["accepted"] else "Abstain") + f" & {reason} " + r"\\")
    lines.extend([r"\bottomrule", r"\end{tabular}",
                  r"\begin{tablenotes}[flushleft]\fontsize{9}{11}\selectfont",
                  r"\item $p_{0.1}$ and $p_{0.3}$ are the frozen weak-label and direct-panel model scores, respectively. The v0.4 rule retains 5/14 high-risk calls (35.7\%) and withholds 9/14 classifications; there are no accepted low-risk calls. These scores are not efflux ratios or validated ADC-specific probabilities.",
                  r"\item Reference support: similarity falls below the frozen boundary in both reference panels. Physicochemical support: maximum absolute robust $z$ exceeds the frozen cutoff of 7. Both support criteria: both conditions hold. Score gray zone: the primary score lies between the low- and high-risk thresholds; this does not imply failed chemical support.",
                  r"\item cys-VC-MMAE is shorthand for Cys-MC-VC-PABC-MMAE, a processing intermediate. Calicheamicin DMH abbreviates $N$-acetyl-$\gamma$-calicheamicin dimethylhydrazide. Source chemical identities and model-standardized representations remain distinct; scores refer to the standardized inputs, while the structures in Figure~\ref{fig:adc} depict source identities.",
                  r"\item Lys-MCC-DM1* depicts the admitted representative stereoisomer; reported source material is an epimeric mixture, not a stereochemically resolved match.",
                  r"\item No ADC-specific classification accuracy is estimated. The scope is free or post-release small molecules, not intact-ADC transport, intracellular exposure or efficacy.",
                  r"\end{tablenotes}", r"\end{threeparttable}", r"\end{table}"])
    TABLE.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    palette.apply_base_style()
    palette.assert_palette([BLUE, GREY, C["green"], C["red"], C["purple"]])
    rows, checks = load_and_verify()
    fig = draw(rows)
    destinations = [FIGURES / "Figure5_adc_triage.pdf", FIGURES / "Figure5_adc_triage.eps",
                    OUT / "Figure5_adc_triage.png"]
    for target in destinations:
        fig.savefig(target, dpi=600 if target.suffix == ".png" else 300)
    plt.close(fig)
    write_table(rows)
    outputs = destinations + [TABLE]
    report = {"source_csv": str(DATA.relative_to(ROOT)), "source_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
              "rdkit_version": rdBase.rdkitVersion, "matplotlib_version": matplotlib.__version__,
              "n_scored": 14, "n_accepted_high": 5, "n_abstained": 9, "n_accepted_low": 0,
              "source_identity_and_formula_checks_passed": len(checks),
              "depictions_from": "source_smiles", "chemical_arrows": 0,
              "support_axis": "physchem_max_abs_robust_z", "support_cutoff": 7,
              "uncertainty_range": "bootstrap_probability_p10 to bootstrap_probability_p90; not confidence intervals",
              "entities": checks,
              "outputs": [{"path": str(p.relative_to(ROOT)), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in outputs]}
    (OUT / "verification.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k not in ["entities", "outputs"]}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

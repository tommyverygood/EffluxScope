"""Figure 4: errors and abstention on model-score and measured-ER scales.

Blue series and pale amber regions use the author-approved reduced palette.
Panel b uses paired bars; centered bold panel headings and frozen geometry remain.
All plotted values come from the rebuilt handoff/fig4_data.json, with explicit
frozen-input checks and a reproducible scaffold-cluster bootstrap provenance.
The fixed 0.49 outcome-aware band is descriptive. Mean error exceeding threshold
distance does not imply that every individual classification is incorrect.
"""

from pathlib import Path
import json
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = HERE.parent / "latex_jcim/tools"
OUTPUT = HERE.parent / "latex_jcim/figures"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(TOOLS))
import palette
FS = palette.FIGURE_STYLE
from jcim_fig34_style import apply_jcim_readability, export_eps_from_pdf

HANDOFF = HERE / "../../handoff"
FROZEN = HERE / "../../数据/processed/abcb1_public_predictor_metrics_v1.csv"

# 端点对齐模型（输出即预测 log10 ER，可进 MAE 轴）
ALIGNED = ["astrazeneca_gnn_mtl_mdck_er_2025", "astrazeneca_gnn_mtl_nih_mdck_er_2025"]
# b 面板的模型顺序：两个对齐 → 本项目 → 两个不同端点
MODEL_ORDER = ALIGNED + ["internal_v0_4_primary_score",
                         "deeppk_pgp_substrate_2026_08_05",
                         "admet_ai_v2_pgp_broccatelli_inhibition"]
MLAB = {
    "astrazeneca_gnn_mtl_mdck_er_2025": "MDCK-ER (aligned)",
    "astrazeneca_gnn_mtl_nih_mdck_er_2025": "NIH-MDCK-ER (aligned)",
    "internal_v0_4_primary_score": "This work (weak label)",
    "deeppk_pgp_substrate_2026_08_05": "Public substrate score",
    "admet_ai_v2_pgp_broccatelli_inhibition": "P-gp inhibition (mism.)",
}
# Reduced palette: model rows are named directly; light/dark blue denotes inside/outside.
MCOL = {cid: palette.C["blue"] for cid in MODEL_ORDER}

MARK = {
    "astrazeneca_gnn_mtl_mdck_er_2025": "s",
    "astrazeneca_gnn_mtl_nih_mdck_er_2025": "D",
    "internal_v0_4_primary_score": "o",
    "deeppk_pgp_substrate_2026_08_05": "^",
    "admet_ai_v2_pgp_broccatelli_inhibition": "v",
}
FOCUS = {"astrazeneca_gnn_mtl_mdck_er_2025", "internal_v0_4_primary_score"}


def load():
    return json.load(open(HANDOFF / "fig4_data.json", encoding="utf-8"))


def verify_palette():
    """出图前守卫：所有语义用色必须来自用户规定的色板或其 tint 派生。

    This figure uses blue for data series and pale amber for diagnostic regions.
    Model identity is given by direct labels and solid/dashed lines in panel a.
    """
    palette.assert_palette(list(MCOL.values()) +
                           [palette.C["amber"], palette.C["blue"]])
    return True


def verify_against_frozen(D):
    """Verify all BAcc points from confusion counts, not an invalid average bound.

    A pooled BAcc need not lie between subgroup BAcc values when class mixes
    differ. Reconstruct TP/TN and denominators directly from frozen rows instead.
    """
    import hashlib
    import pandas as pd
    source = ROOT / D["provenance"]["source"]
    if not source.is_file():
        return verify_aggregate_counts(D)
    assert hashlib.sha256(source.read_bytes()).hexdigest() == D["provenance"]["source_sha256"]
    long = pd.read_csv(source)
    full = pd.read_csv(FROZEN)
    for cid in MODEL_ORDER:
        rows = long[long.comparator_id == cid]
        assert len(rows) == 2624 and rows.external_row_id.is_unique
        y = rows.observed_er2_label.to_numpy(int)
        pred = rows.predicted_label.to_numpy(int)
        inside = abs(rows.observed_log10_mdr1_mdck_er.to_numpy() - D["thr"]) <= D["band_x"]
        for side, mask in [("in", inside), ("out", ~inside)]:
            expected = D["b"][f"{cid}|{side}"]
            tp = int(((y == 1) & (pred == 1) & mask).sum())
            tn = int(((y == 0) & (pred == 0) & mask).sum())
            n_pos, n_neg = int(((y == 1) & mask).sum()), int(((y == 0) & mask).sum())
            assert (tp, tn, n_pos, n_neg) == (expected["tp"], expected["tn"], expected["n_pos"], expected["n_neg"])
            assert abs(.5 * (tp / n_pos + tn / n_neg) - expected["bacc"]) < 1e-14
        actual_full = .5 * ((pred[y == 1] == 1).mean() + (pred[y == 0] == 0).mean())
        expected_full = full[(full.comparator_id == cid) & (full.population == "all_identity_clean")].balanced_accuracy.iloc[0]
        assert abs(actual_full - expected_full) < 1e-14
    print("Verified: 2624 frozen entities; all 10 subgroup and 5 full-set BAcc points match.")
    return 2624


def verify_aggregate_counts(D):
    """Public-archive redraw from audited aggregates without restricted rows.

    This validates aggregate arithmetic, not the original per-molecule inputs
    or the bootstrap calculation. The latter require the CSVs identified by
    hashes in handoff provenance and are not included in the public overlay.
    """
    import csv
    with open(FROZEN, encoding="utf-8") as handle:
        full = {r["comparator_id"]: float(r["balanced_accuracy"])
                for r in csv.DictReader(handle) if r["population"] == "all_identity_clean"}
    totals = []
    for cid in MODEL_ORDER:
        records = [D["b"][f"{cid}|{side}"] for side in ["in", "out"]]
        for record in records:
            assert record["n"] == record["n_pos"] + record["n_neg"]
            assert 0 <= record["tp"] <= record["n_pos"]
            assert 0 <= record["tn"] <= record["n_neg"]
            actual = .5 * (record["tp"] / record["n_pos"] + record["tn"] / record["n_neg"])
            assert abs(actual - record["bacc"]) < 1e-14
        total = {key: sum(r[key] for r in records) for key in ["n", "n_pos", "n_neg", "tp", "tn"]}
        actual_full = .5 * (total["tp"] / total["n_pos"] + total["tn"] / total["n_neg"])
        assert abs(actual_full - full[cid]) < 1e-14
        assert total["n"] == 2624
        totals.append((total["n"], total["n_pos"], total["n_neg"]))
    assert len(set(totals)) == 1
    assert D["band"]["n_in"] + D["band"]["n_out"] == 2624
    print("Public redraw: aggregate confusion-count arithmetic verified; "
          "per-molecule inputs and bootstrap were NOT rerun.")
    return 2624


def density_strip(ax, vals, thr_x, S):
    """顶部密度带：外部集到判定线的带符号距离分布。"""
    ax.hist(vals, bins=70, range=(-1.25, 1.65), color=FS["blue_fill"],
            edgecolor="none", zorder=2)
    ax.axvline(0.0, color=palette.C["grey"], lw=0.8, ls=(0, (3, 2)), zorder=3)
    ax.set_xlim(-1.25, 1.65)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)
    return ax


def triangle_panel(ax, D, S):
    """Rolling MAE relative to threshold distance; shading is a diagnostic aid."""
    xs = np.linspace(-1.25, 1.65, 400)
    # 三角：y = |x|。y > |x| 的区域内，平均误差本身就大于到判定线的距离
    ax.fill_between(xs, np.abs(xs), 1.05, color=palette.tint("amber", 0.10),
                    lw=0, zorder=1)
    ax.plot(xs, np.abs(xs), color=palette.tint("amber", .70), lw=0.9, alpha=0.9, zorder=2)
    ax.axvline(0.0, color=palette.C["grey"], lw=0.8, ls=(0, (3, 2)), zorder=2)

    for cid in ALIGNED:
        a = D["a"][cid]
        ax.plot(a["x"], a["mae"], color=(MCOL[cid] if cid in FOCUS else palette.tint("blue", .55)),
                lw=1.9 if cid in FOCUS else 1.3,
                ls="-" if cid in FOCUS else (0, (4, 2)), zorder=5,
                solid_capstyle="round")
        # 直接标注在线的右端（优于图例框）
        # 标在曲线右段的内侧上方（标在末端外侧会伸进 b 面板）
        # 标在曲线上方/下方留足间距，且避开三角斜边
        # 焦点线标在右端上方；对照线标在其自身最低点下方（该处两线间距最大，
        # 且下方有 0.33 log 单位空白，不会被任一曲线穿过）
        if cid in FOCUS:
            k = int(0.90 * len(a["x"])); dy, va, ha = 16, "bottom", "right"
            xoff = -6
        else:
            k = min(range(len(a["mae"])), key=lambda i: a["mae"][i])
            dy, va, ha = -14, "top", "center"
            xoff = 0
        ax.annotate(MLAB[cid], xy=(a["x"][k], a["mae"][k]),
                    xytext=(xoff, dy), textcoords="offset points",
                    fontsize=S[2] - 0.5, color=MCOL[cid], va=va, ha=ha)

    # Historical fixed descriptive band; not asserted to equal both curves' crossings.
    xb = D["band_x"]
    ax.axvspan(-xb, xb, color=FS["panel"], lw=0, zorder=0)
    # 三条注记纵向分层，避免互撞：带宽在轴外顶部，阴影说明在左上，窗宽在右下
    ax.text(0.985, 0.135,
            f"band |dist| $\\leq$ {xb:.2f} (ER {D['band']['er_lo']:.2f}"
            f"$-${D['band']['er_hi']:.1f})\n"
            f"{D['band']['n_in']:,} of "
            f"{D['band']['n_in'] + D['band']['n_out']:,} compounds",
            transform=ax.transAxes, fontsize=S[2] - 0.5, color="#444444",
            va="bottom", ha="right", linespacing=1.35)
    ax.text(0.02, 0.96, "shaded: mean error exceeds\nthe distance to the line",
            transform=ax.transAxes, fontsize=S[2] - 0.5, color=palette.shade("amber", 1.45),
            va="top", ha="left", linespacing=1.3)


    ax.set_xlim(-1.25, 1.65)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel("Measured $\\log_{10}$ ER $-$ decision line "
                  "($\\log_{10} 2$)", fontsize=S[1])
    ax.set_ylabel("Rolling mean absolute error\n($\\log_{10}$ ER units)",
                  fontsize=S[1])
    ax.tick_params(labelsize=S[2])
    return ax


def coverage_panel(axes, D, S):
    """c 面板：拒判在【自身分数】标尺上是硬边界，映射到实测标尺后几乎消失。

    两条曲线的横轴不是同一个量，双轴叠在一格里读不出各自范围，
    因此改为上下两格：上格自身分数标尺（硬边界所在），下格实测标尺（与 a 同标尺）。
    """
    ax_own, ax_obs = axes
    c2, cov, bd = D["c2"], D["cov"], D["band"]
    assert c2["coordinate"] == "signed_score_minus_0.5"

    # 上格：接受率对分数到自身判定线的距离
    ax_own.plot(c2["x"], [100 * v for v in c2["rate"]], color=palette.C["blue"],
                lw=1.7, zorder=4, solid_capstyle="round")
    ax_own.axvspan(c2["lo_t"] - .5, c2["hi_t"] - .5, color=palette.tint("amber", 0.10), lw=0, zorder=0)
    ax_own.set_xlim(-0.50, 0.50)
    ax_own.set_ylim(0, 82)
    ax_own.set_yticks([0, 40, 80])
    ax_own.set_xticks([-0.4, -0.2, 0.0, 0.2, 0.4])
    ax_own.set_xlabel("Score $-$ own decision line ($p - 0.5$)", fontsize=S[1],
                      labelpad=1.5)
    ax_own.tick_params(labelsize=S[2])
    ax_own.text(0.50, 0.97,
                f"abstain: $p \\in$ ({c2['lo_t']:.2f}, {c2['hi_t']:.2f})\n"
                f"{c2['n_dead']:,}; none accepted",
                transform=ax_own.transAxes, fontsize=S[2] - 1,
                color=palette.shade("amber", 1.45), va="top", ha="center", linespacing=1.3)

    # 下格：接受率对实测距离（与 a 面板同标尺）
    ax_obs.plot(cov["x"], [100 * v for v in cov["rate"]], color=palette.C["blue"],
                lw=1.7, zorder=4, solid_capstyle="round")
    ax_obs.axvspan(-D["band_x"], D["band_x"], color=FS["panel"], lw=0, zorder=0)
    ax_obs.axvline(0.0, color=palette.C["grey"], lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax_obs.set_xlim(-1.25, 1.65)
    ax_obs.set_ylim(0, 58)
    ax_obs.set_yticks([0, 20, 40])
    ax_obs.set_xlabel("Measured $\\log_{10}$ ER $-$ decision line "
                      "($\\log_{10} 2$)", fontsize=S[1], labelpad=1.5)
    ax_obs.tick_params(labelsize=S[2])
    ax_obs.text(0.985, 1.08,
                f"inside $-$ outside retention: {100 * bd['acc_diff']:+.1f} pp "
                f"[{100 * bd['acc_lo']:+.1f}, {100 * bd['acc_hi']:+.1f}]\n"
                f"score-to-endpoint $\\rho$ = {c2['rho_own']:.2f} "
                f"(aligned model {c2['rho_az']:.2f})",
                transform=ax_obs.transAxes, fontsize=S[2] - 0.5, color="#444444",
                va="bottom", ha="right", linespacing=1.35, clip_on=False)

    ax_own.set_ylabel("Accepted (%)", fontsize=S[1], labelpad=2)
    ax_obs.set_ylabel("Accepted (%)", fontsize=S[1], labelpad=2)
    return ax_own


def paired_panel(ax, D, S):
    """Paired horizontal bars with the same frozen BAcc values and intervals."""
    from matplotlib.colors import to_rgb
    from matplotlib.patches import Patch
    ys = np.arange(len(MODEL_ORDER))[::-1]
    for y, cid in zip(ys, MODEL_ORDER):
        inside, outside = D["b"][f"{cid}|in"], D["b"][f"{cid}|out"]
        col = MCOL[cid]
        light = tuple(1 - .30 * (1 - v) for v in to_rgb(col))
        for record, yy, fill in [(inside, y + .14, light), (outside, y - .14, col)]:
            ax.barh(yy, record["bacc"], height=.20, left=0, color=fill,
                    edgecolor="none", linewidth=0, zorder=3)
            ax.errorbar(record["bacc"], yy,
                        xerr=[[record["bacc"] - record["lo"]], [record["hi"] - record["bacc"]]],
                        fmt="none", ecolor=FS["ink"], elinewidth=.85,
                        capsize=1.4, capthick=.85, zorder=5)
        ax.text(0, y + .42, MLAB[cid], fontsize=7.2, color=FS["ink"],
                va="center", ha="left")
        ax.text(1, y + .42, f"{outside['bacc'] - inside['bacc']:+.3f}",
                fontsize=7.2, color=FS["ink"], va="center", ha="right")

    ax.axvline(.5, ymax=.885, color=palette.C["grey"], lw=.85, ls=":", zorder=4)
    ax.set_yticks([])
    ax.set_ylim(-.46, 4.95)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, .25, .5, .75, 1])
    ax.set_xticklabels(["0", "0.25", "0.50", "0.75", "1.00"])
    ax.set_xlabel("Balanced accuracy (forced binary)", fontsize=S[1], labelpad=2)
    ax.tick_params(labelsize=S[2])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.legend(handles=[Patch(facecolor=FS["blue_fill"], edgecolor="none", label="Inside band"),
                       Patch(facecolor=palette.C["blue"], edgecolor="none", label="Outside band")],
              ncol=2, fontsize=7.2, frameon=False, loc="upper center",
              bbox_to_anchor=(.5, 1.015), handlelength=1.4, handleheight=.7,
              handletextpad=.4, columnspacing=1.4, borderpad=0, borderaxespad=0)
    return ax


def draw():
    D = load()
    verify_palette()
    n_ext = verify_against_frozen(D)
    S = (8.5, 7.5, 6.5)

    fig = plt.figure(figsize=(8.5, 5.6))
    fig.suptitle("Prediction error and abstention near the decision threshold",
                 x=0.5, y=0.985, ha="center", va="top", fontsize=14.17,
                 fontweight="bold", color="#000000")
    gs = fig.add_gridspec(2, 2, left=0.083, right=0.940, top=0.835, bottom=0.120,
                          width_ratios=[1.0, 0.86], hspace=0.60, wspace=0.46)
    # a 跨左列两行（三角图需要高度才看得出交点）；b、c 各占右列一行
    gsa = gs[:, 0].subgridspec(2, 1, height_ratios=[0.14, 1.0], hspace=0.04)
    ax_top = fig.add_subplot(gsa[0])
    ax_a = fig.add_subplot(gsa[1])
    ax_b = fig.add_subplot(gs[0, 1])
    bb = ax_b.get_position()
    ax_b.set_position([bb.x0, bb.y0 - .045, bb.width, bb.height + .045])
    gsc = gs[1, 1].subgridspec(2, 1, hspace=1.55)
    ax_c_own = fig.add_subplot(gsc[0])
    ax_c_obs = fig.add_subplot(gsc[1])

    density_strip(ax_top, D["hist"]["vals"], D["thr"], S)
    triangle_panel(ax_a, D, S)

    paired_panel(ax_b, D, S)

    coverage_panel((ax_c_own, ax_c_obs), D, S)
    # Match Figure 5: a bold panel letter and bold title share one baseline.
    renderer = fig.canvas.get_renderer()
    for ax, letter, title in [
        (ax_top, "a", "Error structure relative to the decision line"),
        (ax_b, "b", "Lower accuracy near the decision threshold"),
        (ax_c_own, "c", "Abstention on model and measured scales"),
    ]:
        bb = ax.get_position(); cx = (bb.x0 + bb.x1) / 2
        y = bb.y1 + .035
        lf = FontProperties(family="Arial", size=10, weight="bold")
        tf = FontProperties(family="Arial", size=8.4, weight="bold")
        lw = renderer.get_text_width_height_descent(letter, lf, False)[0]
        tw = renderer.get_text_width_height_descent(title, tf, False)[0]
        gap = 5 * fig.dpi / 72
        left = cx - (lw + gap + tw) / (2 * fig.bbox.width)
        fig.text(left, y, letter, fontproperties=lf, ha="left", va="baseline", color=FS["ink"])
        fig.text(left + (lw + gap) / fig.bbox.width, y, title,
                 fontproperties=tf, ha="left", va="baseline", color=FS["ink"])
    for ax in fig.axes:
        ax.set_facecolor("white")
        for side in ("left", "bottom"):
            ax.spines[side].set_color(FS["border"])

    return fig


def main():
    palette.apply_base_style({"savefig.bbox": None})
    fig = draw()
    stem = "Figure4_error_distance"
    apply_jcim_readability(fig)
    from matplotlib.text import Text
    for text in fig.findobj(match=Text):
        text.set_fontfamily("Arial")
        text.set_color(FS["ink"])
        text.set_alpha(1.0)
    for ext in ("pdf", "png"):
        fig.savefig(OUTPUT / f"{stem}.{ext}", dpi=600 if ext == "png" else None, bbox_inches=None)
    export_eps_from_pdf(OUTPUT / f"{stem}.pdf")
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()

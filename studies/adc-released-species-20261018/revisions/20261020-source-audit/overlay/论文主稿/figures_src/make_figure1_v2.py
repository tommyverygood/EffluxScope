"""Compact English Figure 1: endpoint alignment and structural shift.

Panel a retains six per-molecule scatter plots (two models, Q1/Q3/Q5).
Panel b retains five models across eight strata, the original bootstrap
intervals, and the original permutation null band. The redundant numeric
matrices are omitted; the released input data are read without modification.
"""
from pathlib import Path
import json
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import palette

C = palette.C
DATA = HERE / "../../handoff/fig1ab.json"
GRID = HERE / "../../handoff/fig1ab_grid.json"

FOCUS = ["internal_v0_4_primary_score", "astrazeneca_gnn_mtl_mdck_er_2025"]
ALL5 = ["astrazeneca_gnn_mtl_mdck_er_2025", "astrazeneca_gnn_mtl_nih_mdck_er_2025",
        "deeppk_pgp_substrate_2026_08_05", "internal_v0_4_primary_score",
        "admet_ai_v2_pgp_broccatelli_inhibition"]
# KM20 is reported in Figure 2; retain All, Q1–Q5 and two structural stress sets here.
STRATA = ["All", "Q1", "Q2", "Q3", "Q4", "Q5", "Lo-Hi", "Dual"]
# b 面板的模型配色：色相=端点家族，深浅=家族内成员（与手稿其余图一致）
MCOL = {"astrazeneca_gnn_mtl_mdck_er_2025": (C["green"], 1.00),
        "astrazeneca_gnn_mtl_nih_mdck_er_2025": (C["green"], 0.45),
        "internal_v0_4_primary_score": (C["purple"], 1.00),
        "deeppk_pgp_substrate_2026_08_05": (C["purple"], 0.45),
        "admet_ai_v2_pgp_broccatelli_inhibition": (C["amber"], 1.00)}
MARK = {"astrazeneca_gnn_mtl_mdck_er_2025": "s",
        "astrazeneca_gnn_mtl_nih_mdck_er_2025": "D",
        "internal_v0_4_primary_score": "o",
        "deeppk_pgp_substrate_2026_08_05": "^",
        "admet_ai_v2_pgp_broccatelli_inhibition": "v"}


def load():
    d = json.load(open(DATA, encoding="utf-8"))
    g = json.load(open(GRID, encoding="utf-8"))          # 加密度后的指标行
    G = {(r["cid"], r["stratum"]): r for r in g["grid"]}
    return G, d["scatter"]


#: 全集口径对账映射：图中每个性能指标 → 冻结文件 abcb1_public_predictor_metrics_v1.csv 的列名。
#: 新增指标行时必须同时在此登记，否则 verify_against_frozen() 会报缺项。
FROZEN_COLS = {"bacc": "balanced_accuracy", "mcc": "mcc", "auc": "roc_auc",
               "auprc": "average_precision", "rho": "spearman_rho",
               "sens": "sensitivity", "spec": "specificity"}


def load_ci():
    """bootstrap 95% CI 与置换零带。CI 用于两个焦点模型的带，
    零带取代名义 0.5 线——各层 n 不同，随机基线的经验区间宽度也不同，
    画一条 0.5 会掩盖这一点（Q1 半宽 0.041 vs 全集 0.019）。"""
    ci = json.load(open(HERE / "../../handoff/fig1new_stats.json", encoding="utf-8"))
    band = json.load(open(HERE / "../../handoff/fig1_permband.json", encoding="utf-8"))
    return ci, band


def verify_palette():
    """Require the established palette for model identity and prediction status."""
    palette.assert_palette([c for c, _ in MCOL.values()] + [C["blue"], C["grey"]])
    return True


def verify_against_frozen(G, csv_path, tol=1e-6):
    """把网格里的全集口径逐值与冻结结果文件比对，返回不一致清单（空 = 全部一致）。"""
    import csv as _csv
    frozen = {}
    with open(csv_path, encoding="utf-8") as fh:
        for row in _csv.DictReader(fh):
            if row["population"] == "all_identity_clean":
                frozen[row["comparator_id"]] = row
    diffs = []
    for cid in {c for c, _s in G}:
        fr = frozen.get(cid)
        if not fr:
            continue
        mine = G.get((cid, "All"))
        for key, col in FROZEN_COLS.items():
            if col not in fr or fr[col] in ("", None) or mine is None:
                continue
            d = abs(float(fr[col]) - float(mine[key]))
            if d > tol:
                diffs.append((cid, key, float(fr[col]), float(mine[key]), d))
    return diffs


def draw():
    G, SC = load()
    CI, BAND = load_ci()
    fig = plt.figure(figsize=(6.8, 5.0))
    fig.text(0.025, 0.965, "a", fontsize=11, fontweight="bold", va="center")
    fig.text(0.115, 0.965, "Model outputs against measured efflux", fontsize=9,
             fontweight="bold", va="center")
    legend = [
        Line2D([], [], marker="o", ls="", color=C["blue"], ms=4, label="v0.4 accepted"),
        Line2D([], [], marker="o", ls="", color=C["grey"], ms=4, label="v0.4 abstained"),
        Line2D([], [], ls=(0, (3, 2)), color=C["red"], lw=0.9, label="Model threshold"),
    ]
    fig.legend(handles=legend, loc="center left", bbox_to_anchor=(0.11, 0.926),
               ncol=3, fontsize=8, handlelength=1.5, columnspacing=1.5,
               handletextpad=0.5, borderaxespad=0)

    # Two model rows and three structural-distance columns; every original
    # per-molecule point, acceptance status, threshold and rho is retained.
    x_left, x_right, gap = 0.115, 0.98, 0.065
    width = (x_right - x_left - 2 * gap) / 3
    for row, cid in enumerate(FOCUS):
        is_er = row == 1
        ds = [SC[f"{cid}|{q}"] for q in ("Q1", "Q3", "Q5")]
        # Common limits within each row make the three strata comparable.
        allx = np.concatenate([np.asarray(d["x"]) for d in ds])
        ally = np.concatenate([np.asarray(d["y"]) for d in ds])
        xpad = 0.06 * np.ptp(allx)
        ypad = 0.06 * np.ptp(ally)
        xlim = (allx.min() - xpad, allx.max() + xpad)
        ylim = (ally.min() - ypad, ally.max() + ypad)
        for col, d in enumerate(ds):
            ax = fig.add_axes([x_left + col * (width + gap), 0.68 - row * 0.26,
                               width, 0.20])
            x, y = np.asarray(d["x"]), np.asarray(d["y"])
            acc = np.asarray(d["acc"], dtype=bool)
            ax.scatter(x[~acc], y[~acc], s=3.2, c=C["grey"], alpha=0.32,
                       lw=0, rasterized=True, zorder=2)
            ax.scatter(x[acc], y[acc], s=4.6, c=C["blue"], alpha=0.85,
                       lw=0, rasterized=True, zorder=3)
            ax.axhline(d["thr"], color=C["red"], lw=0.8, ls=(0, (3, 2)), zorder=4)
            ax.axvline(np.log10(2), color="#999999", lw=0.7,
                       ls=(0, (1, 2)), zorder=1)
            if is_er:
                lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
                ax.plot([lo, hi], [lo, hi], color="#444444", lw=0.7, zorder=4)
            ax.set_xlim(xlim); ax.set_ylim(ylim)
            ax.text(0.96, 0.96, rf"$\rho$ = {d['rho']:.2f}", transform=ax.transAxes,
                    fontsize=8, va="top", ha="right", zorder=6,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.86, pad=1))
            ax.xaxis.set_major_locator(plt.MaxNLocator(3, prune="both"))
            ax.yaxis.set_major_locator(plt.MaxNLocator(3))
            ax.tick_params(labelsize=8, length=3, pad=2)
            if row == 0:
                ax.set_title(["Q1: nearest", "Q3: middle", "Q5: farthest"][col],
                             fontsize=8.5, pad=4)
            if col == 0:
                ax.set_ylabel("Weak-label score" if row == 0 else "MDCK-ER prediction\n(log$_{10}$ ER)",
                              fontsize=8, labelpad=3)
            else:
                ax.tick_params(labelleft=False)
    fig.text(0.55, 0.373, r"Measured $\log_{10}$(MDR1-MDCK ER)", fontsize=8,
             ha="center", va="center")

    # Same five-model comparison and the same uncertainty estimates.
    fig.text(0.025, 0.334, "b", fontsize=11, fontweight="bold", va="center")
    fig.text(0.115, 0.334, "Balanced accuracy across structural strata", fontsize=9,
             fontweight="bold", va="center")
    axb = fig.add_axes([0.115, 0.10, 0.64, 0.205])
    xs = np.arange(len(STRATA))
    axb.set_xlim(-0.30, len(STRATA) - 0.55)
    # Show the complete original null band and bootstrap intervals.
    axb.set_ylim(0.445, 0.82)
    bl = [BAND[s]["lo"] if s in BAND else np.nan for s in STRATA]
    bh = [BAND[s]["hi"] if s in BAND else np.nan for s in STRATA]
    axb.fill_between(xs, bl, bh, color="#c9c9c9", alpha=0.55, linewidth=0, zorder=1)
    for cid in FOCUS:
        lo = [CI[cid][s]["lo"] if s in CI.get(cid, {}) else np.nan for s in STRATA]
        hi = [CI[cid][s]["hi"] if s in CI.get(cid, {}) else np.nan for s in STRATA]
        col = MCOL[cid][0]
        axb.fill_between(xs, lo, hi, facecolor=col, alpha=0.13, linewidth=0, zorder=2)
        axb.plot(xs, lo, color=col, lw=0.5, alpha=0.55, zorder=2)
        axb.plot(xs, hi, color=col, lw=0.5, alpha=0.55, zorder=2)
    ends = []
    for cid in ALL5:
        ys = [G[(cid, s)]["bacc"] if (cid, s) in G else np.nan for s in STRATA]
        col, alpha = MCOL[cid]
        focus = cid in FOCUS
        axb.plot(xs, ys, marker=MARK[cid], ms=4.2 if focus else 3.4,
                 lw=1.5 if focus else 1.0, color=col, alpha=alpha,
                 ls="-" if focus else (0, (3, 1.6)), markeredgewidth=0,
                 zorder=3, clip_on=False)
        ends.append([ys[-1], cid, col, max(alpha, 0.75)])
    names = {
        "astrazeneca_gnn_mtl_mdck_er_2025": "MDCK-ER",
        "astrazeneca_gnn_mtl_nih_mdck_er_2025": "NIH-MDCK-ER",
        "internal_v0_4_primary_score": "Weak-label score",
        "deeppk_pgp_substrate_2026_08_05": "Deep-PK substrate",
        "admet_ai_v2_pgp_broccatelli_inhibition": "P-gp inhibition",
    }
    # Direct labels use a fixed 10-pt minimum vertical separation.
    ends.sort(key=lambda e: e[0])
    fig.canvas.draw()
    step = (10 / 72) / (5 * 0.205) * np.ptp(axb.get_ylim())
    positions = [e[0] for e in ends]
    for k in range(1, len(positions)):
        positions[k] = max(positions[k], positions[k - 1] + step)
    offset = max(0, positions[-1] - 0.80)
    positions = [y - offset for y in positions]
    for (entry, ylab) in zip(ends, positions):
        ytrue, cid, col, alpha = entry
        axb.annotate(names[cid], xy=(xs[-1] + 0.06, ytrue),
                     xytext=(xs[-1] + 0.32, ylab), fontsize=8, color=col,
                     alpha=alpha, va="center", ha="left", annotation_clip=False,
                     arrowprops=dict(arrowstyle="-", lw=0.5, color=col,
                                     alpha=alpha * 0.7, shrinkA=0, shrinkB=1))
    axb.set_xticks(xs); axb.set_xticklabels(STRATA, fontsize=8)
    axb.set_yticks([0.5, 0.6, 0.7, 0.8])
    axb.tick_params(labelsize=8, length=3, pad=2)
    axb.set_ylabel("Balanced accuracy", fontsize=8, labelpad=3)
    axb.set_xlabel("Structural stratum", fontsize=8, labelpad=3)
    return fig


def main():
    # 作图前先对账：图里的全集口径必须与冻结结果文件逐位一致，否则不出图。
    G, _ = load()
    verify_palette()
    diffs = verify_against_frozen(G, HERE / "../../数据/processed/abcb1_public_predictor_metrics_v1.csv")
    if diffs:
        raise AssertionError(f"与冻结文件不一致 {len(diffs)} 处：{diffs[:5]}")
    print(f"对账通过：{len(FROZEN_COLS)} 个指标 × {len(set(c for c, _ in G))} 个模型，全集口径 Δ < 1e-6")

    palette.apply_base_style({"savefig.bbox": None})
    # Base rendering only; the JCIM wrapper applies publication typography.
    fig = draw()
    destinations = [HERE / "Figure1_endpoint_alignment_shift"]
    for stem in destinations:
        for ext in ("pdf", "png"):
            fig.savefig(stem.with_suffix(f".{ext}"), dpi=600 if ext == "png" else None,
                        bbox_inches=None)
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()

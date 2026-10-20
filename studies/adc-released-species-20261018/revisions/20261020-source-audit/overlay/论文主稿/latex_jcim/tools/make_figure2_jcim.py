"""JCIM-specific Figure 2 typography and stroke revision.

Numeric inputs and plotting semantics are unchanged from figures_src/make_figure2.py.
Outputs are isolated to latex_jcim/figures.
"""
"""Figure 2：两种划分策略下的化学空间分离度与模型表现。

版式参照 DataSAIL（Joeres et al., Nat Commun 2025, 10.1038/s41467-025-58606-8）Fig 3：
每行 = 两幅 t-SNE 投影（划分前 / 划分后）+ 一幅双纵轴性能对比。

本文的两行不是两个数据集，而是**两种把外部集判为域外的方式**：
  行 1（a–c）Murcko 骨架清洁 —— 常被当作 OOD 划分，但只清掉 62/2,624；
  行 2（d–f）KMeans-OOD20  —— 真正把化学空间分开的划分。
对照两行即本文论点：骨架排除不构成结构域外，性能几乎不动；
只有真正分开化学空间时，异质弱标签分数才掉到与随机不可区分。

所有数值现场重算，全集口径与冻结的
数据/processed/abcb1_public_predictor_metrics_v1.csv 逐位一致（出图前强制校验）。
嵌入用 Morgan(r=2,2048) 指纹的 Jaccard 距离 + t-SNE(perplexity=30, seed=20260817)。
"""
from pathlib import Path
import json
import sys
import shutil
import subprocess

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parents[2] / "figures_src"
OUTPUT = Path(__file__).resolve().parents[1] / "figures"
sys.path.insert(0, str(HERE))
import palette
from palette import C

HANDOFF = HERE / "../../handoff"
FROZEN = HERE / "../../数据/processed/abcb1_public_predictor_metrics_v1.csv"

#: 行定义：(行键, 面板字母, 划分名, 三个切片的显示名)
ROWS = [
    ("row1", ("_", "c", "d"), "Murcko scaffold cleaning",
     ["All (100%)", "Scaffold-clean", "Scaffold-removed"]),
    ("row2", ("_", "e", "f"), "KMeans structural split",
     ["All (100%)", "KMeans in-domain", "KMeans-OOD20"]),
]

#: 模型显示顺序（与 Figure 1 的色彩语义一致：端点对齐 = 绿，弱标签 = 紫，错配 = 琥珀）
MODEL_ORDER = ["astrazeneca_gnn_mtl_mdck_er_2025", "astrazeneca_gnn_mtl_nih_mdck_er_2025",
               "internal_v0_4_primary_score", "deeppk_pgp_substrate_2026_08_05",
               "admet_ai_v2_pgp_broccatelli_inhibition"]
MODEL_LABEL = {"astrazeneca_gnn_mtl_mdck_er_2025": "MDCK-ER (aligned)",
               "astrazeneca_gnn_mtl_nih_mdck_er_2025": "NIH-MDCK-ER (aligned)",
               "internal_v0_4_primary_score": "This work (weak label)",
               "deeppk_pgp_substrate_2026_08_05": "Public substrate prob.",
               "admet_ai_v2_pgp_broccatelli_inhibition": "P-gp inhibition (mism.)"}


def model_colors():
    """端点家族定色（外层），家族内成员用深浅（内层）—— 与 Figure 1 同一方案。"""
    def light(hexv, f=0.45):
        r, g, b = (int(hexv[i:i + 2], 16) for i in (1, 3, 5))
        return "#%02X%02X%02X" % tuple(int(v + (255 - v) * f) for v in (r, g, b))
    return {"astrazeneca_gnn_mtl_mdck_er_2025": C["green"],
            "astrazeneca_gnn_mtl_nih_mdck_er_2025": palette.tint("green", 0.55),
            "internal_v0_4_primary_score": C["purple"],
            "deeppk_pgp_substrate_2026_08_05": palette.tint("purple", 0.55),
            "admet_ai_v2_pgp_broccatelli_inhibition": C["amber"]}


#: 全集口径对账映射（新增指标须在此登记）
FROZEN_COLS = {"bacc": "balanced_accuracy", "auc": "roc_auc"}


def verify_palette(colmap):
    """出图前守卫：语义用色必须来自色板或其 tint 派生。"""
    palette.assert_palette(list(colmap.values()) + [C["blue"], C["grey"]])
    return True


def verify_against_frozen(perf, tol=1e-6):
    """图中的 All (100%) 列必须与冻结结果文件逐位一致，否则不出图。"""
    fr = pd.read_csv(FROZEN)
    fr = fr[fr.population == "all_identity_clean"].set_index("comparator_id")
    diffs = []
    for cid in MODEL_ORDER:
        mine = perf[(perf.cid == cid) & (perf.split == "All (100%)")]
        if not len(mine) or cid not in fr.index:
            continue
        for key, col in FROZEN_COLS.items():
            a, b = float(fr.loc[cid, col]), float(mine[key].iloc[0])
            if abs(a - b) > tol:
                diffs.append((cid, key, a, b, abs(a - b)))
    return diffs


def load():
    emb = pd.read_csv(HANDOFF / "fig2_emb.csv")
    ref = np.load(HANDOFF / "fig2_ref_xy.npy")
    perf = pd.read_csv(HANDOFF / "fig2_perf.csv")
    leak = json.load(open(HANDOFF / "fig2_leak.json", encoding="utf-8"))
    hist = json.load(open(HANDOFF / "fig2_hist.json", encoding="utf-8"))
    ci = json.load(open(HANDOFF / "fig2_ci.json", encoding="utf-8"))
    return emb, ref, perf, leak, hist, ci


def scatter_panel(ax, ref_xy, ext, *, split_col, split_val, title, S, show_axis_cue):
    """一幅 t-SNE 投影。§6.6：去刻度，用角标注明轴。

    参考面板始终画在底层（灰），外部集按该划分着色：
    保留侧 = 蓝，被划出侧 = 红（本文的告警色，专用于"被排除/拒判"）。
    """
    ax.scatter(ref_xy[:, 0], ref_xy[:, 1], s=5, c=C["grey"], alpha=0.35,
               linewidths=0, zorder=1, rasterized=True)
    if split_val is None:                     # 划分前：外部集不分色
        ax.scatter(ext.tx, ext.ty, s=4, c=C["blue"], alpha=0.45,
                   linewidths=0, zorder=2, rasterized=True)
    else:
        keep = ext[~ext[split_col]]
        out = ext[ext[split_col]]
        ax.scatter(keep.tx, keep.ty, s=4, c=C["blue"], alpha=0.40,
                   linewidths=0, zorder=2, rasterized=True)
        ax.scatter(out.tx, out.ty, s=9, facecolors="none", edgecolors=C["red"],
                   linewidths=0.85, zorder=3, rasterized=True)
    ax.set_title(title, fontsize=9.0 if split_val is not None else S[0], loc="left", pad=3)
    ax.set_aspect("equal", adjustable="box")   # 等距嵌入：纵横比锁 1，否则形状失真
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.85); sp.set_edgecolor("#999999")
    if show_axis_cue:                          # §6.6 角标代替轴标签
        ax.annotate("", xy=(0.14, 0.035), xytext=(0.035, 0.035),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", lw=0.85, color="#666666"))
        ax.annotate("", xy=(0.035, 0.145), xytext=(0.035, 0.035),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", lw=0.85, color="#666666"))
        ax.text(0.155, 0.030, "t-SNE 1", fontsize=S[2] - 1, color="#666666",
                transform=ax.transAxes, va="center", ha="left")
        ax.text(0.030, 0.160, "t-SNE 2", fontsize=S[2] - 1, color="#666666",
                transform=ax.transAxes, va="bottom", ha="left", rotation=90)


def hist_panel(ax, emb, hist, S):
    """左列下格：D5 分布 + 参考 LOO 分层线 + 两种划分各自挑出的子集分布。

    要点：骨架清洁挑出的 62 个分子 D5 中位与保留侧几乎相同（0.755 vs 0.773），
    而 KMeans-OOD20 明显偏远（0.800）—— 即骨架排除没有在结构距离上分开任何东西。
    """
    d5 = emb.joint_model_d5.values
    bins = np.linspace(0.25, 0.90, 46)
    ax.hist(d5, bins=bins, color=C["grey"], alpha=0.30, linewidth=0,
            label="All (2,624)", zorder=2)
    ax.hist(emb[emb.removed_scaffold].joint_model_d5, bins=bins, histtype="step",
            color=C["red"], linewidth=0.9, label="Scaffold-rm. (62)", zorder=4)
    ax.hist(emb[emb.is_km_ood].joint_model_d5, bins=bins, histtype="step",
            color=C["blue"], linewidth=0.9, label="KMeans-OOD (544)", zorder=3)
    for name, xv in hist["loo_q"].items():
        ax.axvline(xv, color="#777777", lw=0.85, ls=(0, (2, 1.8)), zorder=5)
        ax.annotate(name, xy=(xv, 0), xytext=({"Q50": -2, "Q80": -5, "Q95": 6}.get(name, 0), -20 if name != "Q95" else -31), textcoords="offset points",
                    fontsize=S[2] - 1.5, color="#666666", ha="center", va="top",
                    annotation_clip=False)
    ax.set_xlabel(r"Structural distance $D_5$", fontsize=S[1], labelpad=24)
    ax.set_ylabel("External compounds", fontsize=S[1])
    ax.tick_params(labelsize=S[2])
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=1,
              fontsize=S[2] - 1, loc="upper right",
              handlelength=1.1, handletextpad=0.5, labelspacing=0.55, borderpad=0.15,
              bbox_to_anchor=(1.02, 1.02))
    ax.margins(x=0.02)


def perf_panel(ax, perf, leak, row, splits, S, ci, *, show_ylab, show_leg):
    """双纵轴性能对比（照 DataSAIL Fig 3 的 c/f 面板）：左轴分组柱 = BAcc，右轴点 = 到参考面板的接近度。"""
    cols = model_colors()
    verify_palette(cols)
    nm, ns = len(MODEL_ORDER), len(splits)
    w = 0.8 / nm
    for k, cid in enumerate(MODEL_ORDER):
        vals, xs = [], []
        for j, sp in enumerate(splits):
            r = perf[(perf.row == row) & (perf.split == sp) & (perf.cid == cid)]
            vals.append(float(r.bacc.iloc[0]) if len(r) else np.nan)
            xs.append(j - 0.4 + w * (k + 0.5))
        ax.bar(xs, vals, width=w * 0.92, color=cols[cid], linewidth=0,
               zorder=3, label=MODEL_LABEL[cid] if show_leg else None)
        # bootstrap 95% CI：n=62 的骨架层半宽约 0.13，没有误差条时那一组
        # 看起来像"对齐模型下降、弱标签上升"，实为纯噪声。
        er = np.array([[vals[j] - ci[f"{cid}|{row}|{sp}"]["lo"],
                        ci[f"{cid}|{row}|{sp}"]["hi"] - vals[j]]
                       for j, sp in enumerate(splits)]).T
        ax.errorbar(xs, vals, yerr=er, fmt="none", ecolor="#3a3a3a",
                    elinewidth=0.85, capsize=1.5, capthick=0.85, zorder=5)
    ax.axhline(0.5, color="#888888", lw=0.85, ls=(0, (3, 2)), zorder=2)
    # Chance reference remains explicit as the dashed line and 0.50 tick.
    ax.set_ylim(0.48, 0.80)
    ax.set_yticks(np.arange(0.50, 0.801, 0.05))
    ax.set_xticks(range(ns))
    NL = chr(10)   # f-string 表达式内不能含反斜杠（Py<3.12）
    labs = []
    for s in splits:
        head = s.replace(" (100%)", "").replace("KMeans ", "KMeans" + NL).replace("Scaffold-", "Scaffold" + NL).replace("KMeans-OOD20", "KMeans" + NL + "OOD20")
        labs.append(head + NL + "n=" + format(leak[row + "|" + s]["n"], ","))
    ax.set_xticklabels(labs, fontsize=S[2] - 0.5)
    if show_ylab:
        ax.set_ylabel("Balanced accuracy", fontsize=S[1])
    ax.tick_params(labelsize=S[2])

    # 右轴：该切片与参考面板的平均 top-5 相似度（越高 = 越接近参考 = 越像"内插"）
    ax2 = ax.twinx()
    sims = [leak[f"{row}|{sp}"]["mean_sim"] for sp in splits]
    ax2.plot(range(ns), sims, marker="D", ms=4.2, lw=1.0, color="#3a3a3a",
             mfc="white", mew=0.9, zorder=6, clip_on=False)
    ax2.set_ylim(0.185, 0.265)
    ax2.set_ylabel("Mean top-5 similarity\nto reference panel", fontsize=S[1] - 0.5,
                   color="#3a3a3a")
    ax2.tick_params(labelsize=S[2], colors="#3a3a3a")
    for sp in ("top",):
        ax2.spines[sp].set_visible(False)
    return ax2


def draw():
    emb, ref_xy, perf, leak, hist, ci = load()
    S = (10.0, 9.5, 8.0)

    fig = plt.figure(figsize=(7.4, 5.65))
    gs = fig.add_gridspec(2, 3, left=0.080, right=0.790, top=0.900, bottom=0.205,
                          width_ratios=[1.0, 1.0, 1.30], hspace=0.72, wspace=0.48)

    # (a) 划分前只画一次：两行用的是同一批外部分子与同一嵌入，
    #     重复两幅会让读者误以为是两个数据集。它跨两行居中。
    axA = fig.add_subplot(gs[0, 0])
    scatter_panel(axA, ref_xy, emb, split_col="removed_scaffold", split_val=None,
                  title="Before any split", S=S, show_axis_cue=True)

    leg_handles = None
    for i, (row, letters, split_name, splits) in enumerate(ROWS):
        col = "removed_scaffold" if row == "row1" else "is_km_ood"
        n_out = int(emb[col].sum())
        ax1 = fig.add_subplot(gs[i, 1])
        scatter_panel(ax1, ref_xy, emb, split_col=col, split_val="yes",
                      title=f"{split_name}\n{n_out:,} set aside "
                            f"({100 * n_out / len(emb):.1f}%)",
                      S=S, show_axis_cue=False)
        ax2 = fig.add_subplot(gs[i, 2])
        perf_panel(ax2, perf, leak, row, splits, S, ci, show_ylab=True, show_leg=(i == 0))
        if i == 0:
            leg_handles = ax2.get_legend_handles_labels()
        for letter, ax in zip(letters[1:], (ax1, ax2)):
            bb = ax.get_position()
            fig.text(bb.x0 - 0.040, bb.y1 + 0.016, letter, fontsize=S[0] + 1.5,
                     fontweight=700, va="bottom", ha="left")
    axH = fig.add_subplot(gs[1, 0])
    hist_panel(axH, emb, hist, S)
    bbH = axH.get_position()
    fig.text(bbH.x0 - 0.046, bbH.y1 + 0.016, "b", fontsize=S[0] + 1.5,
             fontweight=700, va="bottom", ha="left")

    bbA = axA.get_position()
    fig.text(bbA.x0 - 0.046, bbA.y1 + 0.016, "a", fontsize=S[0] + 1.5,
             fontweight=700, va="bottom", ha="left")

    # 图例（§7.3 置于白处）
    h, l = leg_handles
    scat = [Line2D([], [], marker="o", ls="", ms=3.6, mfc=C["grey"], mec="none",
                   label="Reference panel (443)"),
            Line2D([], [], marker="o", ls="", ms=3.6, mfc=C["blue"], mec="none",
                   label="External, retained"),
            Line2D([], [], marker="o", ls="", ms=4.2, mfc="none", mec=C["red"],
                   mew=0.8, label="Set aside by split"),
            Line2D([], [], marker="D", ls="-", lw=1.0, ms=4.0, color="#3a3a3a",
                   mfc="white", mew=0.9, label="Similarity to reference (right axis)")]
    fig.legend(handles=scat, loc="lower left", bbox_to_anchor=(0.080, 0.001),
               ncol=2, frameon=False, fontsize=S[2] - 0.5,
               handletextpad=0.4, columnspacing=1.3)
    fig.legend(handles=h, labels=l, loc="center left", bbox_to_anchor=(0.836, 0.55),
               frameon=False, fontsize=S[2] - 0.5, handletextpad=0.5,
               labelspacing=0.42)
    fig.text(0.787, 0.952, "higher = better \u2191", fontsize=S[2] - 1,
             color="#666666", va="bottom", ha="right")
    return fig


def main():
    emb, ref_xy, perf, leak, hist, ci = load()
    diffs = verify_against_frozen(perf)
    if diffs:
        raise AssertionError(f"与冻结文件不一致 {len(diffs)} 处：{diffs[:4]}")
    print(f"对账通过：{len(FROZEN_COLS)} 个指标 × {len(MODEL_ORDER)} 个模型，"
          f"All (100%) 口径 Δ < 1e-6")
    palette.apply_base_style({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "axes.linewidth": 0.90, "xtick.major.width": 0.85, "ytick.major.width": 0.85}, cjk=False)
    fig = draw()
    stem = "Figure2_split_separation"
    for ext in ("pdf", "png"):
        fig.savefig(OUTPUT / f"{stem}.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    gs = shutil.which("gs")
    if gs is None:
        raise RuntimeError("Ghostscript is required to regenerate the EPS from the PDF")
    subprocess.run([gs, "-dNOPAUSE", "-dBATCH", "-sDEVICE=eps2write",
                    "-dEPSCrop", "-dLanguageLevel=3",
                    "-sOutputFile=" + str(OUTPUT / f"{stem}.eps"),
                    str(OUTPUT / f"{stem}.pdf")], check=True)
    print("done")


if __name__ == "__main__":
    main()

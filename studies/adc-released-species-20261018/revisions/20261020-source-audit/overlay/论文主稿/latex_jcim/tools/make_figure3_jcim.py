"""Figure 3：预先声明门槛下的独立验证。

版式参照三篇近两年论文：
  · Matbench Discovery (Nat Mach Intell 2025, 10.1038/s42256-025-01055-1)
      —— 面板上方叠加样本密度；把误差与判定阈值画进同一坐标；
         指标画成对"已验证条数"的累积曲线而非单点。
  · UQ under data shifts (JCIM 2026, 10.1021/acs.jcim.5c02381)
      —— 可靠性差异落成一个可报的数，bootstrap 显著性画进主图。
  · TDC leaderboard audit (JCIM 2026, 10.1021/acs.jcim.6c00819)
      —— "阶段 × 对象"的通过/不通过矩阵；阴性结果照实呈现。

四个面板：
  a  ASAP 2024 来源留出：冻结模型概率 vs 实测 log2 外排比，事前门槛
     （骨架 bootstrap 95% 下界 > 0）标出，上方叠加样本密度。
  b  三个来源的端点诊断，区分 GI50 耐药倍数与直接运输外排比。
  c  外部集（n = 2,624）覆盖率—可靠性曲线：拒判的价值依赖基础模型。
  d  法规严格池盲测的 8 条门槛判定 + 混淆结构，fail 如实呈现。

ASAP 来源留出的运输关联满足事前门槛；Roth 为回顾性身份核验后的耐药
端点诊断。不同模型及端点不构成同一个运输结论的三次重复验证。

所有数值由 数据/processed 下的冻结结果文件现场读取；出图前对账，
不一致则拒绝出图。
"""
from pathlib import Path
import json
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

TOOLS = Path(__file__).resolve().parent
HERE = TOOLS.parent.parent / "figures_src"
OUTPUT = TOOLS.parent / "figures"
sys.path.insert(0, str(HERE))
import palette
from jcim_fig34_style import apply_jcim_readability, export_eps_from_pdf

DATA = HERE / "../../handoff/fig3_data.json"
PROC = HERE / "../../数据/processed"

C = palette.C
S = (8.5, 7.0, 6.0)
SFOOT = 5.5          # §5.2 三级字号阶梯

# 与 Fig 1/2 一致的模型配色（§4.1 颜色即交叉引用）
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
MLAB = {"astrazeneca_gnn_mtl_mdck_er_2025": "MDCK-ER (aligned)",
        "astrazeneca_gnn_mtl_nih_mdck_er_2025": "NIH-MDCK-ER (aligned)",
        "internal_v0_4_primary_score": "This work (weak label)",
        "deeppk_pgp_substrate_2026_08_05": "Public substrate prob.",
        "admet_ai_v2_pgp_broccatelli_inhibition": "P-gp inhibition (mism.)"}
FOCUS = ["astrazeneca_gnn_mtl_mdck_er_2025", "internal_v0_4_primary_score"]

# a 面板两个冻结模型：同一族（本项目自有冻结模型），用深浅区分
# a 面板画的是本项目的冻结筛查器 → 用 purple（色板语义 weak_label），与 c 面板
# 的「本项目弱标签」同色；blue 在本图专供门槛判定状态，避免同色两义。
ACOL = {"frozen_ABCB1_v0_1": (C["purple"], 1.00),
        "frozen_ABCB1_assay_aware_v0_3_dev": (C["purple"], 0.45)}
ALAB = {"frozen_ABCB1_v0_1": "fixed v0.1",
        "frozen_ABCB1_assay_aware_v0_3_dev": "fixed v0.3-dev (assay-aware)"}
STATUS_LABEL = {"pass": "MET", "fail": "NOT MET"}


def load():
    return json.load(open(DATA, encoding="utf-8"))


def verify_palette():
    """出图前守卫：所有语义用色必须来自用户规定的色板。

    色相语义（palette.SEMANTIC）：green=端点对齐、purple=异质弱标签、
    amber=端点错配对照、blue=接受/通过、red=拒判/不通过。
    族内成员由元组第二项的 alpha 表示深浅，色相本身不变。
    """
    palette.assert_palette([c for c, _ in MCOL.values()] +
                           [c for c, _ in ACOL.values()] +
                           [gate_color(s) for s in ("pass", "fail", "not supported")])
    return True


def verify_against_frozen(D):
    """Check all panels, including the identity-audited Roth endpoint."""
    import csv
    bad = []
    add = json.load(open(PROC / "abcb1_asap2025_sourceheldout_postscoring_addendum_v1.json",
                         encoding="utf-8"))
    for k, g in D["a"]["gate"].items():
        fz = add["gate_outcome"][k]
        for f in ("spearman_rho", "ci95_lower", "ci95_upper", "n", "status"):
            if isinstance(fz[f], str):
                if str(g[f]) != fz[f]:
                    bad.append(f"a.{k}.{f}: {g[f]} vs {fz[f]}")
            elif abs(float(g[f]) - float(fz[f])) > 1e-8:
                bad.append(f"a.{k}.{f}: {g[f]} vs {fz[f]}")
    with open(PROC / "baseline_roth_external_metrics_v2_identity_safe.csv",
              encoding="utf-8-sig") as fh:
        roth = next(r for r in csv.DictReader(fh)
                    if r["transporter_gene"] == "ABCB1" and r["model_name"] == "morgan_logistic")
    rb = next(b for b in D["b"] if b["label"].startswith("Roth"))
    for field, source in [("n", "n_compounds"),
                          ("rho", "spearman_probability_vs_continuous_rr"),
                          ("p", "spearman_p_value")]:
        if abs(float(rb[field]) - float(roth[source])) > 1e-8:
            bad.append(f"b.Roth.{field}: {rb[field]} vs {roth[source]}")
    if rb.get("endpoint") != "GI50 resistance ratio":
        bad.append("b.Roth.endpoint: must identify GI50 resistance ratio, not direct transport ER")
    if rb.get("excluded_identity_ambiguous_compounds") != roth["excluded_identity_ambiguous_compounds"]:
        bad.append("b.Roth.identity: exclusion must match the identity-audited source")
    asap_b = next(b for b in D["b"] if b["label"].startswith("ASAP"))
    asap_a = D["a"]["gate"]["frozen_ABCB1_v0_1"]
    for field, source in [("n", "n"), ("rho", "spearman_rho"),
                          ("lo", "ci95_lower"), ("hi", "ci95_upper")]:
        if abs(float(asap_b[field]) - float(asap_a[source])) > 1e-8:
            bad.append(f"b.ASAP.{field}: does not match panel a")
    bs = json.load(open(PROC / "abcb1_regulatory_v4_blind_summary_v1.json", encoding="utf-8"))
    for f, key in [("bacc", "balanced_accuracy"), ("mcc", "mcc"), ("auc", "roc_auc"),
                   ("sens", "sensitivity"), ("spec", "specificity")]:
        if abs(float(D["d"][f]) - float(bs["frozen_primary_metrics"][key])) > 1e-8:
            bad.append(f"d.{f}: {D['d'][f]} vs {bs['frozen_primary_metrics'][key]}")
    if [g["status"] for g in D["d"]["gates"]] != [g["status"] for g in bs["claim_gates"]]:
        bad.append("d.gates: 门槛判定序列与冻结文件不一致")
    # c 面板全覆盖点必须等于冻结指标文件的全集值
    with open(PROC / "abcb1_public_predictor_metrics_v1.csv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for cid, s in D["c"].items():
        m = [r for r in rows if r["comparator_id"] == cid
             and r["population"] == "all_identity_clean"]
        if m and abs(float(s["bacc"][-1]) - float(m[0]["balanced_accuracy"])) > 1e-6:
            bad.append(f"c.{cid}: 全覆盖 {s['bacc'][-1]} vs 冻结 {m[0]['balanced_accuracy']}")
    return bad


def gate_color(status):
    """门槛判定的统一取色（b、d 面板共用，登记在 palette.SEMANTIC）。

    pass          有区间且满足门槛           → blue
    fail          有区间且明确跌破门槛        → red
    not supported 仅点估计、无区间，无法评价  → grey（不等于失败）
    """
    return {"pass": C["blue"], "fail": C["red"]}.get(status, C["grey"])


def panel_a(ax_main, ax_top, ax_stats, ax_legend, D):
    """Keep every annotation outside the measured data area."""
    from matplotlib.lines import Line2D
    a = D["a"]
    x = np.asarray(a["x"], float)
    open_mask = np.asarray(a["band"]) == "outside_empirical_support"
    for key, col in [("frozen_ABCB1_v0_1", "p01"),
                     ("frozen_ABCB1_assay_aware_v0_3_dev", "p03")]:
        y = np.asarray(a[col], float)
        c, al = ACOL[key]
        ax_main.scatter(x[~open_mask], y[~open_mask], s=18, color=c, alpha=al,
                        linewidths=0, zorder=3, label=ALAB[key])
        ax_main.scatter(x[open_mask], y[open_mask], s=18, facecolors="none",
                        edgecolors=c, alpha=al, linewidths=0.85, zorder=3)
    ax_main.axvline(1.0, color="#999999", lw=0.85, ls=(0, (2, 2)), zorder=1)
    ax_main.set_xlabel("Measured log2 (MDR1-MDCK efflux ratio)\nDashed line: ER = 2", fontsize=9.0)
    ax_main.set_ylabel("Fixed model score", fontsize=9.0)
    ax_main.tick_params(labelsize=9.0)
    ax_main.set_xticks([0, 2, 4, 6])
    ax_main.set_ylim(0.26, 0.97)
    ax_main.margins(x=0.05)

    ax_top.hist(x, bins=18, color="#c9c9c9", linewidth=0)
    ax_top.set_xlim(ax_main.get_xlim())
    ax_top.set_yticks([]); ax_top.set_xticks([])
    for sp in ("top", "right", "left"):
        ax_top.spines[sp].set_visible(False)

    ax_stats.axis("off")
    ax_stats.text(0, 1.0, f"n = {a['n']} rows, {a['n_scaf']} scaffolds",
                  fontsize=9.0, color="#555555", va="top")
    ax_stats.text(0, 0.72, "Criterion: scaffold-bootstrap 95% lower bound > 0",
                  fontsize=9.0, color="#333333", va="top")
    for yy, key, label in [(0.44, "frozen_ABCB1_v0_1", "v0.1"),
                            (0.16, "frozen_ABCB1_assay_aware_v0_3_dev", "v0.3")]:
        g = a["gate"][key]
        ax_stats.text(0, yy,
            f"{label}  $\\rho$ = {g['spearman_rho']:.3f} "
            f"[{g['ci95_lower']:.3f}, {g['ci95_upper']:.3f}]  "
            f"{STATUS_LABEL.get(g['status'], 'NOT ASSESSED')}",
            fontsize=9.0, color="#333333", va="top")
    ax_legend.axis("off")
    h, labels = ax_main.get_legend_handles_labels()
    h.append(Line2D([], [], marker="o", ls="none", ms=4.2,
                    markerfacecolor="none", markeredgecolor="#555555", markeredgewidth=0.85))
    labels.append("open = outside empirical support")
    ax_legend.legend(h, labels, fontsize=9.0, frameon=False,
                     loc="upper left", bbox_to_anchor=(0, 1),
                     borderaxespad=0, handletextpad=0.4, labelspacing=0.35)


def panel_b(ax, ax_labels, D):
    """Separate a fixed-width study column from the numerical forest plot."""
    B = D["b"]
    ys = np.arange(len(B))[::-1] * 1.5
    for yy, b in zip(ys, B):
        col = gate_color(b["status"])
        if b["lo"] is not None:
            ax.plot([b["lo"], b["hi"]], [yy, yy], color=col, lw=1.5,
                    solid_capstyle="butt", zorder=3)
            for cap in (b["lo"], b["hi"]):
                ax.plot([cap, cap], [yy - 0.10, yy + 0.10], color=col, lw=1, zorder=3)
        ax.scatter([b["rho"]], [yy], s=34,
                   marker="o" if b["status"] == "pass" else "s",
                   color=col, zorder=4, linewidths=0)
    ax.axvline(0.0, color="#777777", lw=0.85, zorder=1)
    ax.set_xlim(-0.90, 0.98)
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.set_ylim(-0.63, ys[0] + 0.63)
    ax.set_yticks([])
    ax.tick_params(labelsize=9.0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Spearman $\\rho$", fontsize=9.0)
    ax_labels.set_ylim(ax.get_ylim()); ax_labels.set_xlim(0, 1); ax_labels.axis("off")
    for yy, b in zip(ys, B):
        label = b["label"].replace(" panel", "")
        if "ASAP" in label:
            title = f"ASAP 2024 · n = {b['n']}"
            obj = "fixed v0.1; MDR1-MDCK ER"
        elif "Roth" in label:
            title = f"Roth 2026 · n = {b['n']}*"
            obj = "Morgan baseline; GI50 RR"
        else:
            title = f"{label} · n = {b['n']}"
            obj = "fixed v0.1; MDCK ER"
        ax_labels.text(0, yy + 0.47, title, fontsize=9.0, color="#1a1a1a",
                       ha="left", va="top", linespacing=1.2)
        obj_y = yy + 0.20
        ax_labels.text(0, obj_y, obj, fontsize=8.6, color="#777777",
                       style="italic", ha="left", va="top")
        val = (f"$\\rho$ = {b['rho']:+.3f}; 95% CI" if b["lo"] is not None
               else f"$\\rho$ = {b['rho']:+.3f}; p = {b['p']:.3f}")
        ax_labels.text(0, yy - 0.20, val, fontsize=8.6,
                       color=gate_color(b["status"]), ha="left", va="top")


def panel_c(ax, D):
    """覆盖率—可靠性：拒判换来的可靠性依赖基础模型（照 Matbench Discovery 的累积曲线）。"""
    for cid, s in D["c"].items():
        col, al = MCOL[cid]
        is_focus = cid in FOCUS
        ax.plot(s["cov"], s["bacc"], marker=MARK[cid], ms=4.0 if is_focus else 3.2,
                lw=1.5 if is_focus else 1.0, color=col, alpha=al,
                ls="-" if is_focus else (0, (3, 1.6)), markeredgewidth=0, zorder=3)
    ax.axhline(0.5, color="#999999", lw=0.7, ls=(0, (1, 2)), zorder=1)
    ax.text(2, 0.507, "0.5", fontsize=S[2] - 1, color="#777777", va="bottom", ha="left")
    ax.set_xlabel("Coverage retained after abstention (% of external set)\nhigher = better $\\uparrow$", fontsize=S[2])
    ax.set_ylabel("Balanced accuracy", fontsize=S[2])
    ax.set_xlim(-3, 118)
    ax.set_ylim(0.485, 1.12)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.tick_params(labelsize=S[2])
    # 头条数字（§2.5）：对齐模型峰值→全覆盖的落差，与弱标签的落差对比
    az = D["c"]["astrazeneca_gnn_mtl_mdck_er_2025"]
    wk = D["c"]["internal_v0_4_primary_score"]
    ia = max(range(len(az["bacc"])), key=lambda i: az["bacc"][i])
    iw = max(range(len(wk["bacc"])), key=lambda i: wk["bacc"][i])
    ax.text(0.30, 0.995,
            f"peak $\\rightarrow$ full coverage:  aligned "
            f"{az['bacc'][ia]:.3f} $\\rightarrow$ {az['bacc'][-1]:.3f} "
            f"({az['bacc'][-1] - az['bacc'][ia]:+.3f})\n"
            f"weak label {wk['bacc'][iw]:.3f} $\\rightarrow$ {wk['bacc'][-1]:.3f} "
            f"({wk['bacc'][-1] - wk['bacc'][iw]:+.3f})",
            transform=ax.transAxes, fontsize=S[2] - 1.5, color="#444444",
            va="top", ha="left", linespacing=1.5)
    # 直接标注（§6.3）优于图例框
    ends = sorted([(s["bacc"][-1], cid) for cid, s in D["c"].items()])
    # 间距按实际渲染行高换算成数据单位，再对相撞簇整体居中展开
    # （单向抬升会把整簇推离真值，且末端会溢出轴顶）
    fig = ax.figure
    fig.canvas.draw()
    _r = fig.canvas.get_renderer()
    _pr = ax.text(0, 0, "Ag", fontsize=7.2)
    _h = _pr.get_window_extent(_r).height
    _pr.remove()
    y0, y1 = ax.get_ylim()
    step = 1.45 * _h * (y1 - y0) / ax.get_window_extent(_r).height
    ys = [e[0] for e in ends]
    i = 0
    while i < len(ys):
        j = i
        while j + 1 < len(ys) and ys[j + 1] - ys[j] < step:
            j += 1
        if j > i:
            mid = sum(ys[i:j + 1]) / (j - i + 1)
            k = j - i + 1
            for m in range(k):
                ys[i + m] = mid + (m - (k - 1) / 2) * step
        i = j + 1
    for (yv, (y_true, cid)) in zip(ys, ends):
        col, al = MCOL[cid]
        ax.annotate(MLAB[cid], xy=(101, y_true), xytext=(107, yv),
                    fontsize=S[2] - 1, color=col, alpha=max(al, 0.75),
                    va="center", ha="left", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=col,
                                    alpha=max(al, 0.75) * 0.7, shrinkA=0, shrinkB=1))


def panel_d(axes, D):
    """d 面板：8 条事前门槛的判定，符号网格 + 分组汇总条。

    版式参照 robvis（McGuinness & Higgins, Res Synth Methods 2020）的
    traffic-light plot 与 summary plot 配对做法，并按 PROBAST+AI 的三层聚合
    思路把门槛分组：判定不用文字而用统一符号，同名门槛按重抽样口径并列，
    使「按分子计通过、按实验簇计不通过」的翻转一眼可见。

    左格为判定网格（行=组，列=组内门槛），右格为每组通过比例的堆叠条。
    所有判定与数值直接取自冻结的评分摘要，不重算。
    """
    ax_g, ax_s = axes
    d2 = D["d2"]
    ROWH = 1.55          # 每组占的行高（数据单位）
    rows = d2["rows"]
    nrow = len(rows)
    ncol = max(r["n"] for r in rows)

    # ---- 左格：符号网格 ----
    ax_g.set_xlim(-1.55, ncol - 0.28)
    ax_g.set_ylim(-0.55, (nrow - 1) * 1.55 + 2.18)
    ax_g.invert_yaxis()
    for sp in ax_g.spines.values():
        sp.set_visible(False)
    ax_g.set_xticks([]); ax_g.set_yticks([])

    for ri, r in enumerate(rows):
        y = ri * 1.55 + 0.55
        # 组名（行首，左对齐于轴外）
        ax_g.text(-0.10, y - 0.13, r["group"], fontsize=S[1] - 0.5, va="center",
                  ha="right", color="#1a1a1a")
        # 组内列按该组门槛数均分行宽：第三组只有 2 条，不该挤在 3 列的格位上
        step = ncol / r["n"]
        for ci, it in enumerate(r["items"]):
            ok = it["status"] == "pass"
            col = gate_color(it["status"])
            cx = ci * step + 0.28
            # 判定符号：实心圆 + 白色勾/叉（线段绘制，衬线体无这两个字形）
            # 用 scatter 而非 Circle：Circle 半径按数据单位，非等比坐标下会变成椭圆
            ax_g.scatter([cx], [y - 0.13], s=118, c=col, marker="o",
                         linewidths=0, zorder=3)
            # 勾/叉用 marker 叠加（衬线体无这两个字形，且 marker 不随坐标比例变形）
            # marker "P" 是加号不是勾：通过则用折线画勾，不通过用 X marker
            if ok:
                ax_g.plot([cx - 0.115, cx - 0.035, cx + 0.125],
                          [y - 0.135, y - 0.075, y - 0.215],
                          color="white", lw=1.15, solid_capstyle="round",
                          transform=ax_g.transData, zorder=4)
            else:
                ax_g.scatter([cx], [y - 0.13], s=26, c="white", marker="X",
                             linewidths=0, zorder=4)
            # 门槛名、规则、实际值
            # 名称与门槛规则合成一段，避免两段文字随行名换行而贴住
            # 偏移按行高的比例给（数据单位），画布尺寸变化时相对间距不变
            ax_g.text(cx, y - 0.13 - 0.16 * ROWH, f"{it['short']}\n({it['rule']})",
                      fontsize=S[2] - 0.5, va="bottom", ha="center",
                      color="#1a1a1a", linespacing=1.2)
            ax_g.text(cx, y - 0.13 + 0.13 * ROWH, it["val"], fontsize=S[2] - 1,
                      va="top", ha="center", color=col)

        # 组内空列留白（不画符号，也不画占位）
        if ri < nrow - 1:
            ax_g.plot([-1.5, ncol - 0.35], [ri * 1.55 + 1.30] * 2,
                      color="#ececec", lw=0.6, zorder=0,
                      solid_capstyle="butt")

    # ---- 右格：分组通过比例堆叠条 ----
    ax_s.set_xlim(0, 1.0)
    ax_s.set_ylim(-0.55, (nrow - 1) * 1.55 + 2.18)
    ax_s.invert_yaxis()
    for sp in ("top", "right", "left"):
        ax_s.spines[sp].set_visible(False)
    ax_s.set_yticks([])
    ax_s.set_xticks([0, 0.5, 1.0])
    ax_s.set_xticklabels(["0", "50", "100"], fontsize=S[2])
    ax_s.set_xlabel("Criteria met (%)", fontsize=S[1], labelpad=1.5)
    ax_s.tick_params(labelsize=S[2])

    for ri, r in enumerate(rows):
        y = ri * 1.55 + 0.55
        f = r["npass"] / r["n"]
        ax_s.barh(y - 0.13, f, height=0.30, color=C["blue"], zorder=3)
        ax_s.barh(y - 0.13, 1 - f, left=f, height=0.30, color=C["red"], zorder=3)
        ax_s.text(0.5, y + 0.24, f"{r['npass']} / {r['n']}", fontsize=S[2] - 0.5,
                  va="top", ha="center", color="#444444")

    # 总体判定：单独一条，加粗区分（照 robvis summary plot 的 Overall 行）
    yo = (nrow - 1) * 1.55 + 1.82
    ax_g.plot([-1.5, ncol - 0.35], [yo - 0.42] * 2, color="#c8c8c8", lw=0.8,
              zorder=0, solid_capstyle="butt")
    fo = d2["total_pass"] / d2["total"]
    ax_s.barh(yo, fo, height=0.34, color=C["blue"], zorder=3)
    ax_s.barh(yo, 1 - fo, left=fo, height=0.34, color=C["red"], zorder=3)
    ax_g.text(-0.10, yo, f"Overall assessment:\n{STATUS_LABEL.get(d2['overall'], 'NOT ASSESSED')}",
              fontsize=S[1] - 0.5, va="center", ha="right", color=C["red"],
              fontweight=700)
    ax_s.text(0.5, yo - 0.26, f"{d2['total_pass']} / {d2['total']}",
              fontsize=S[2] - 0.5, va="bottom", ha="center", color="#444444",
              fontweight=700)

    return ax_g


def draw():
    D = load()
    fig = plt.figure(figsize=(8.5, 7.8))
    # The upper row has dedicated annotation areas: no text overlays the data.
    ax_stats = fig.add_axes([0.078, 0.805, 0.405, 0.110])
    ax_top = fig.add_axes([0.078, 0.745, 0.405, 0.043])
    ax_a = fig.add_axes([0.078, 0.510, 0.405, 0.230])
    ax_legend = fig.add_axes([0.078, 0.377, 0.405, 0.070])
    ax_labels = fig.add_axes([0.555, 0.510, 0.205, 0.310])
    ax_b = fig.add_axes([0.775, 0.510, 0.200, 0.310])
    ax_c = fig.add_axes([0.078, 0.120, 0.290, 0.200])
    ax_d = fig.add_axes([0.620, 0.100, 0.235, 0.220])
    ax_d2 = fig.add_axes([0.872, 0.100, 0.103, 0.220])
    panel_a(ax_a, ax_top, ax_stats, ax_legend, D)
    panel_b(ax_b, ax_labels, D)
    panel_c(ax_c, D)
    panel_d((ax_d, ax_d2), D)
    titles = [(0.078, 0.940, "a", "Source-held-out endpoint association\n(ASAP 2024)"),
              (0.555, 0.940, "b", "Endpoint diagnostics across sources\nand models"),
              (0.078, 0.340, "c", "Abstention performance depends\non the model"),
              (0.620, 0.340, "d", "Predefined criteria in the scoring-locked\nregulatory set")]
    for xx, yy, letter, title in titles:
        fig.text(xx - 0.052, yy, letter, fontsize=11, fontweight="bold", va="bottom")
        fig.text(xx, yy, title, fontsize=8.5, va="bottom", color="#1a1a1a")
    fig.text(0.555, 0.910, "ASAP criterion: point estimate > 0 and\nbootstrap 95% lower bound > 0",
             fontsize=8.3, color="#555555", va="top", linespacing=1.3)
    fig.text(0.555, 0.442, "* Roth: retrospective identity audit; Deruxtecan excluded.\n"
             "GI50 resistance ratio is not direct transport ER.\n"
             "Grey: point estimates only; blue: ASAP bootstrap 95% CI.",
             fontsize=8.0, color="#777777", va="top", linespacing=1.3)
    fig.text(0.078, 0.006,
             "Archived analyses; models were fixed before evaluation and were not refitted.\n"
             f"Panel (d): {D['d2']['n']} compounds, {D['d2']['clusters']} source experimental clusters; "
             "symbols follow the risk-of-bias traffic-light convention.",
             fontsize=7.2, color="#555555", va="bottom", ha="left", linespacing=1.5)
    return fig


def main():
    palette.apply_base_style()
    D = load()
    verify_palette()
    bad = verify_against_frozen(D)
    if bad:
        raise AssertionError("与冻结文件不一致：\n  " + "\n  ".join(bad))
    print(f"对账通过：a 面板 2 个门槛、b 面板 Roth 身份核验与 ASAP 区间、d 面板 8 条判定与 5 个指标、"
          f"c 面板 {len(D['c'])} 个模型全覆盖口径，均与冻结文件一致")
    fig = draw()
    stem = "Figure3_validation_gates"
    apply_jcim_readability(fig)
    for ext in ("pdf", "png"):
        fig.savefig(OUTPUT / f"{stem}.{ext}", dpi=600 if ext == "png" else None)
    export_eps_from_pdf(OUTPUT / f"{stem}.pdf")
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()

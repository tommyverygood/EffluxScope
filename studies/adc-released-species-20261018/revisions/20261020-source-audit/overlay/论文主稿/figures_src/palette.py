"""论文一图表统一配色。

来源：用户 2026-08-17 提供的配色卡，按卡面标注的 RGB/HEX 原样录入并交叉校验。
每张图都从这里取色，不要在作图脚本里手写十六进制值。

用法：
    from palette import C, cycle, apply_base_style
    fig, ax = plt.subplots()
    apply_base_style()
    ax.plot(x, y, color=C["blue"])
"""
from pathlib import Path
import json
import re

_HERE = Path(__file__).resolve().parent
_CFG = json.load(open(_HERE / "palette.json", encoding="utf-8"))

#: 名称 → HEX
C = dict(_CFG["hex_by_name"])
#: 卡面从左到右的顺序
ORDER = list(_CFG["order"])
#: 按卡面顺序的 HEX 列表，可直接喂给 color cycle
CYCLE = [C[n] for n in ORDER]

#: 名称 → (r, g, b)，0-255
RGB = {p["name"]: tuple(p["rgb"]) for p in _CFG["palette"]}


def rgb01(name):
    """返回 0-1 归一化的 RGB 三元组（matplotlib 需要时用）。"""
    r, g, b = RGB[name]
    return (r / 255, g / 255, b / 255)


def cycle(*names):
    """按给定名称取一组颜色；不给名称时返回卡面完整顺序。

    cycle("blue", "grey")  → ['#2B6CE2', '#878787']
    """
    return [C[n] for n in names] if names else list(CYCLE)


#: 中文字体候选（本机 MacTeX/系统已装，按优先级）
CJK_FONTS = ["Songti SC", "PingFang SC", "Heiti SC", "STHeiti", "Arial Unicode MS"]


def pick_cjk_font():
    """返回本机可用的第一个中文字体名；都不可用时返回 None。"""
    import matplotlib.font_manager as fm
    have = {f.name for f in fm.fontManager.ttflist}
    for name in CJK_FONTS:
        if name in have:
            return name
    return None


def apply_base_style(rc=None, cjk=True):
    """设定本项目图表的基础 rcParams（字号、线宽、去顶右框、色环、中文字体）。

    只做与配色和可读性相关的最小设定，不锁定图幅——图幅由每张图自己决定。
    cjk=True 时把可用中文字体插到 sans-serif 首位，使中文标签不再显示为方框。
    """
    import matplotlib as mpl
    from cycler import cycler

    families = list(mpl.rcParams["font.sans-serif"])
    if cjk:
        f = pick_cjk_font()
        if f and f not in families:
            families.insert(0, f)

    base = {
        "font.sans-serif": families,
        "axes.unicode_minus": False,   # 中文字体缺 U+2212，负号用 ASCII hyphen
        "axes.prop_cycle": cycler(color=CYCLE),
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "axes.labelsize": 9,
        "axes.titlesize": 9.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.fontsize": 8,
        "legend.frameon": False,
        "lines.linewidth": 1.4,
        "lines.markersize": 5,
        "grid.color": C["grey"],
        "grid.alpha": 0.25,
        "grid.linewidth": 0.6,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,   # 保证 PDF 里字体可编辑、可被期刊排版系统识别
        "ps.fonttype": 42,
    }
    if rc:
        base.update(rc)
    mpl.rcParams.update(base)
    return base


#: 本项目的语义映射建议（每张图可覆盖，但同一含义在全篇应保持同色）
SEMANTIC = {
    "accepted": "blue",        # v0.4 接受层
    "forced": "grey",          # 全覆盖强制二分类
    # 注意与下面的 screener_abstain 区分，二者不是同义：
    #   forced_on_rejected = 对被拒判分子仍强制输出二分类（「不听拒判会怎样」的对照层），
    #                        它衡量的是错误代价，用 red 合理
    #   screener_abstain   = 拒判这个动作本身（模型不给判定），用 grey
    "forced_on_rejected": "red",   # 拒判后强制分类（对照层，非拒判本身）
    "endpoint_aligned": "green",   # 端点对齐的连续 ER 模型
    "weak_label": "purple",    # 异质弱标签分数
    "mismatched": "amber",     # 端点错配负对照
    # teal 目前全图册未使用：Figure 5 的 13 个实体全是 ADC 释放实体、无非 ADC 对照，
    # 色相不区分任何东西；且该图的结构保留 CPK 化学惯例，RDKit 的氟色 #33CCCC
    # 与 teal #11B8C5 距离仅 40.1，同图并置会混。若将来出现 ADC 与非 ADC 并置的
    # 图，可启用此条——但那张图不应同时画分子结构。
    "adc": "teal",             # ADC 释放物（保留待用）
    # 三态筛查器的输出状态（Figure 1a 与 Figure 5 共用一套编码，因为是同一冻结规则）：
    #   scored   = 模型给出风险判定
    #   abstain  = 模型自认不可靠、不给判定（这不是失败，故不用 red——
    #              red 在 Figure 3 已表示「明确跌破预先声明的门槛」）
    "screener_scored": "blue",
    "screener_abstain": "grey",
    # Figure 5 把「拒判」进一步按判据分色，因为两个判据是独立的，
    # 合成单一灰色会把一条实质结论藏起来（较大的分子可能只触发参考距离
    # 而不触发理化包络）。这三色在 Figure 5 不与模型身份冲突——该图没有
    # 任何模型对照，只有一个冻结筛查器作用于实体；而 red 不用于此，
    # 它在 Figure 3 表示「跌破预先声明的门槛」。
    # 色距实测：teal 与 blue 只差 85.4，二者在 Figure 5c 并列时偏近；
    # green 与 blue 差 164.0，故用 green 承担「仅参考面板之外」。
    # green 在 Figure 1–4 是端点对齐模型的身份色，但 Figure 5 没有任何模型
    # 对照（只有一个冻结筛查器作用于实体），故色相不会被误读——这与 purple
    # 在本图承担「两判据兼有」是同一条理由。
    "abstain_reference": "green",    # 仅「位于两个参考面板之外」
    # red 是本色板中语义最重的一色，全图册已承担五种用途，各图注分别声明：
    #   Fig 1  判定阈线
    #   Fig 2  骨架剔除的 62 个分子
    #   Fig 3  gate_fail（跌破预先声明的门槛）
    #   Fig 4  误差 > 判定距离的失效区、与规则自身的拒判窗
    #   Fig 5c 「超出理化包络」这一判据
    # 其中 Fig 3 与 Fig 5c 的冲突最直接（同为「判定结果」语义），因此 Fig 5 的
    # 图注与正文都显式声明此处的红不是门槛失败，否则读者会把拒判读成失败。
    # 新图若再用 red，先检查是否会与同版面的 Fig 3 / Fig 5 并读。
    "abstain_physchem": "red",      # 仅「超出理化包络」
    "abstain_both": "purple",       # 两个判据同时触发
    # 门槛判定三态（Figure 3 的 b、d 面板共用一套编码）：
    #   pass      = 有区间且满足门槛
    #   fail      = 有区间且明确跌破门槛
    #   inconcl.  = 只有点估计、无区间，门槛无法评价（不等于失败）
    "gate_pass": "blue",
    "gate_fail": "red",
    "gate_inconclusive": "grey",
}

def tint(name, alpha):
    """色板色与白底按 alpha 混合，返回等效不透明 hex。

    用于「族内成员」的浅色变体：色相由 SEMANTIC 决定，深浅表示族内位置。
    必须由此派生，不得手写近似 hex —— 否则改色板时浅色不会跟着变。
    """
    h = C[name].lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    f = lambda v: round(v * alpha + 255 * (1 - alpha))
    return "#%02X%02X%02X" % (f(r), f(g), f(b))


def assert_palette(colors, allow_neutral=True):
    """出图前守卫：所有承载语义的颜色必须来自色板或其 tint 派生。

    中性灰（R=G=B，含黑白）视为版式用色，不受色板管辖。
    """
    ok = set(C.values())
    for a in (0.30, 0.40, 0.42, 0.45, 0.50, 0.55, 0.60, 0.70):
        ok |= {tint(n, a) for n in C}
    bad = []
    for c in colors:
        u = str(c).upper()
        if u in {x.upper() for x in ok}:
            continue
        if allow_neutral and re.fullmatch(r"#([0-9A-F]{2})\1\1", u):
            continue
        bad.append(c)
    if bad:
        raise AssertionError(f"非色板色: {sorted(set(bad))}")
    return True


def shade(name, k):
    """色板色的明度变体：k<1 与白混合（更浅），k>1 与黑混合（更深）。

    用于连续色阶：色相锁定为色板色，只改明度。k=1 返回色板色本身。
    """
    h = C[name].lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    if k <= 1.0:
        f = lambda v: round(v * k + 255 * (1 - k))
    else:
        f = lambda v: round(v * (2.0 - k))
    return "#%02X%02X%02X" % (f(r), f(g), f(b))

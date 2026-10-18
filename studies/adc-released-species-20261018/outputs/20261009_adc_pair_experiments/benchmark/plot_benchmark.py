#!/usr/bin/env python3
"""Render already completed benchmark metrics; no recalculation or tuning."""
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
with (HERE / "metrics.csv").open(encoding="utf-8-sig", newline="") as handle:
    rows = list(csv.DictReader(handle))
main = [r for r in rows if r["version"] == "identity_audited_v2"]
POLICIES = ["v01_margin", "v03_margin", "full_v04_pool_primary_margin", "no_local_pool_primary_margin"]
LABELS = ["v0.1 own margin", "v0.3 own margin", "Full v0.4", "No local + margin"]
COLORS = ["#7F8C8D", "#2171B5", "#B35D2B", "#765A9D"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 200, "pdf.fonttype": 42})


def error(r, metric):
    p = float(r[metric])
    return [[max(0, p - float(r[metric + "_ci95_low"]))], [max(0, float(r[metric + "_ci95_high"]) - p)]]


fixed = [next(r for r in main if r["target"] == "fixed_full_coverage" and r["policy"] == p) for p in POLICIES]
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
x = np.arange(4)
for i, r in enumerate(fixed):
    axes[0].bar(i, float(r["balanced_accuracy"]), width=.64, color=COLORS[i], alpha=.90)
    axes[0].errorbar(i, float(r["balanced_accuracy"]), yerr=error(r, "balanced_accuracy"), color="#222222", capsize=4, fmt="none")
    axes[0].text(i, float(r["balanced_accuracy_ci95_high"]) + .018, f"{float(r['balanced_accuracy']):.3f}", ha="center", fontsize=10)
    axes[2].bar(i, float(r["npv"]), width=.64, color=COLORS[i], alpha=.90)
    axes[2].errorbar(i, float(r["npv"]), yerr=error(r, "npv"), color="#222222", capsize=4, fmt="none")
    axes[2].text(i, float(r["npv_ci95_high"]) + .018, f"{float(r['npv']):.3f}", ha="center", fontsize=10)
for j, (metric, color, label) in enumerate([
    ("positive_class_coverage", "#CB6A4C", "ER >= 2"),
    ("negative_class_coverage", "#4C85AA", "ER < 2"),
]):
    offset = -.18 if j == 0 else .18
    for i, r in enumerate(fixed):
        axes[1].bar(i + offset, float(r[metric]) * 100, width=.32, color=color,
                    label=label if i == 0 else None)
        axes[1].errorbar(i + offset, float(r[metric]) * 100, yerr=np.array(error(r, metric))*100,
                        color="#222222", capsize=3, fmt="none")
axes[0].axhline(.5, color="#999999", ls=":", lw=1)
axes[0].set(ylim=(0, 1), ylabel="Balanced accuracy", title="a  Accuracy at equal retained count")
axes[1].set(ylim=(0, 30), ylabel="True-class coverage (%)", title="b  Retention within each true class")
axes[1].legend(frameon=False, loc="upper right")
axes[2].set(ylim=(0, 1.04), ylabel="Negative predictive value", title="c  Reliability of low-risk calls")
for ax in axes:
    ax.set_xticks(x, LABELS, rotation=25, ha="right")
    ax.grid(axis="y", alpha=.15)
    ax.set_axisbelow(True)
fig.suptitle("Frozen scores on Biogen/Fang: all four policies retain 453 / 2624 (17.2637%)\n95% whole-cluster bootstrap intervals; retrospective comparison", fontsize=12)
fig.savefig(HERE / "fixed_coverage_comparison.png")
fig.savefig(HERE / "fixed_coverage_comparison.pdf")
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), constrained_layout=True)
for policy, label, color in zip(POLICIES, LABELS, COLORS):
    rr = sorted([r for r in main if r["policy"] == policy and r["target"] != "native_working_point"], key=lambda r: int(r["selected_n"]))
    # Retain exact fixed-full count rather than duplicate the adjacent nominal17.3% point in plot.
    rr = [r for r in rr if r["target"] != "coverage_17_3pct"]
    coverage = np.array([float(r["coverage"]) * 100 for r in rr])
    for ax, metric in zip(axes, ["balanced_accuracy", "npv"]):
        values = np.array([float(r[metric]) for r in rr])
        low = np.array([float(r[metric + "_ci95_low"]) for r in rr])
        high = np.array([float(r[metric + "_ci95_high"]) for r in rr])
        ax.plot(coverage, values, marker="o", color=color, label=label, lw=1.6)
        ax.fill_between(coverage, low, high, color=color, alpha=.10)
for ax in axes:
    ax.axvline(453 / 2624 * 100, ls=":", color="#555555", lw=1)
    ax.set(xlabel="Overall coverage (%)", xlim=(3, 27), ylim=(.45, 1))
    ax.grid(alpha=.15)
axes[0].set(ylabel="Balanced accuracy", title="a  Accuracy versus coverage")
axes[1].set(ylabel="Negative predictive value", title="b  Low-risk NPV versus coverage")
axes[0].legend(frameon=False, fontsize=9, loc="lower right")
fig.suptitle("Full / no-local policies rank only inside their accepted pool\nTargets above the pool capacity are unavailable; abstentions are not force filled", fontsize=11)
fig.savefig(HERE / "coverage_curves.png")
fig.savefig(HERE / "coverage_curves.pdf")
plt.close(fig)

manifest = {
    "status": "rendered_from_completed_metrics",
    "matplotlib_version": matplotlib.__version__,
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "input_sha256": hashlib.sha256((HERE / "metrics.csv").read_bytes()).hexdigest(),
    "fixed_figure": "all four policies retain exactly453; full uses original mask; no-local is reranked within551 accepted pool",
    "curve_figure": "same targets except the nominal17.3% adjacent point omitted visually; exact table retains every target",
    "outputs": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob("*.png"))} | {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob("*.pdf"))},
}
(HERE / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
completion_path = HERE / "completion.json"
if completion_path.exists():
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    completion["figure_rendered_utc"] = datetime.now(timezone.utc).isoformat()
    completion["output_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(HERE.iterdir()) if p.is_file() and p.name != "completion.json"
    }
    completion_path.write_text(json.dumps(completion, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"status": "rendered", "figures": list(manifest["outputs"])}, ensure_ascii=False))

"""Rebuild Figure 4 from frozen per-molecule predictions; no model fitting.

Run from any directory with numpy and pandas. The confidence intervals are a
new, explicitly documented reconstruction, not a claim about the unknown old
interval generator. The outcome-aware band is held fixed at its historical
0.49 log10-unit half-width during resampling.
"""
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs/20261018_manuscript_cleanup/figure4"
SOURCE = ROOT / "数据/processed/abcb1_public_predictor_predictions_long_v1.csv"
INPUT = ROOT / "数据/processed/abcb1_public_predictor_benchmark_input_v1.csv"
METRICS = ROOT / "数据/processed/abcb1_public_predictor_metrics_v1.csv"
HANDOFF = ROOT / "handoff/fig4_data.json"
ALIGNED = ["astrazeneca_gnn_mtl_mdck_er_2025", "astrazeneca_gnn_mtl_nih_mdck_er_2025"]
MODELS = ALIGNED + ["internal_v0_4_primary_score", "deeppk_pgp_substrate_2026_08_05",
                   "admet_ai_v2_pgp_broccatelli_inhibition"]
SEED = 20261018
N_BOOT = 2000


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bacc(y, predicted):
    return float((np.mean(predicted[y == 1] == 1) + np.mean(predicted[y == 0] == 0)) / 2)


def rolling(x, values, grid, width, minimum=20):
    result = {"x": [], "value": [], "n": []}
    for z in grid:
        keep = np.abs(x - z) <= width / 2
        if keep.sum() >= minimum:
            result["x"].append(float(z))
            result["value"].append(float(np.mean(values[keep])))
            result["n"].append(int(keep.sum()))
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    backup = OUT / "fig4_data.before.json"
    old = json.loads(backup.read_text()) if backup.exists() else None
    long = pd.read_csv(SOURCE)
    base = pd.read_csv(INPUT)
    assert len(base) == 2624 and base.external_row_id.is_unique
    ids = base.external_row_id.to_list()
    models = {m: long[long.comparator_id == m].set_index("external_row_id").loc[ids] for m in MODELS}
    for data in models.values():
        assert len(data) == len(base) and data.index.is_unique
        for key in ["observed_log10_mdr1_mdck_er", "observed_er2_label", "v0_4_accepted", "bootstrap_cluster"]:
            assert np.array_equal(data[key].to_numpy(), base[key].to_numpy()), key
    scaffold = base.murcko_scaffold_smiles.fillna("")
    nonempty = scaffold != ""
    assert (base.loc[nonempty, "bootstrap_cluster"] == scaffold[nonempty]).all()
    acyclic = base.loc[~nonempty]
    assert (acyclic.bootstrap_cluster == "ACYCLIC::" + acyclic.standardized_inchikey.str.split("-").str[0]).all()
    y = base.observed_er2_label.to_numpy(int)
    measured = base.observed_log10_mdr1_mdck_er.to_numpy(float)
    accepted = base.v0_4_accepted.to_numpy(bool)
    score = models[MODELS[2]].score.to_numpy(float)
    x = measured - np.log10(2)
    inside = np.abs(x) <= .49
    assert np.array_equal(y, (measured >= np.log10(2)).astype(int))
    if old is not None:
        assert np.allclose(x, old["hist"]["vals"], rtol=0, atol=0)
    assert inside.sum() == 1349
    assert accepted.sum() == 453
    points = []
    D = {"thr": float(np.log10(2)), "window": .15, "band_x": .49,
         "a": {}, "hist": {"vals": x.tolist()}, "b": {}}
    for m in ALIGNED:
        err = np.abs(models[m].score.to_numpy(float) - measured)
        r = rolling(x, err, np.arange(-1.2, 1.620001, .02), .15, minimum=25)
        # Restrict the historical plotting grid to <= 1.60 (floating tolerance).
        r = {k: [v for v, xx in zip(vals, r["x"]) if xx <= 1.600001] for k, vals in r.items()}
        if old is not None:
            np.testing.assert_allclose(r["x"], old["a"][m]["x"], rtol=0, atol=1e-14)
            np.testing.assert_allclose(r["value"], old["a"][m]["mae"], rtol=0, atol=1e-14)
            assert r["n"] == old["a"][m]["n"]
        D["a"][m] = {"x": r["x"], "mae": r["value"], "n": r["n"]}
        points += [{"panel": "a", "model": m, "x": xx, "value": v, "n": n}
                   for xx, v, n in zip(r["x"], r["value"], r["n"])]
    r = rolling(x, accepted, D["a"][ALIGNED[0]]["x"], .15, minimum=25)
    if old is not None:
        np.testing.assert_allclose(r["value"], old["cov"]["rate"], rtol=0, atol=1e-14)
        assert r["n"] == old["cov"]["n"]
    D["cov"] = {"x": r["x"], "rate": r["value"], "n": r["n"]}
    points += [{"panel": "c_measured", "model": MODELS[2], "x": xx, "value": v, "n": n}
               for xx, v, n in zip(r["x"], r["value"], r["n"])]
    # Validate the old absolute-margin curve before replacing its misleading
    # symmetric-window display with the signed model-score scale.
    if old is not None:
        old_margin = rolling(abs(score - .5), accepted, old["c2"]["x"], .05)
        np.testing.assert_allclose(old_margin["value"], old["c2"]["rate"], rtol=0, atol=1e-14)
        assert old_margin["n"] == old["c2"]["n"]
    r = rolling(score - .5, accepted, np.arange(-.5, .50001, .008), .05)
    dead = (score > .35) & (score < .70)
    assert int(dead.sum()) == 1689 and not accepted[dead].any()
    D["c2"] = {"x": r["x"], "rate": r["value"], "n": r["n"],
               "coordinate": "signed_score_minus_0.5", "window": .05,
               "lo_t": .35, "hi_t": .70, "n_dead": int(dead.sum()),
               "acc_dead": float(accepted[dead].mean()), "interval": "open",
               "rho_own": float(pd.Series(score).rank().corr(pd.Series(measured).rank())),
               "rho_az": float(pd.Series(models[MODELS[0]].score.to_numpy()).rank().corr(pd.Series(measured).rank()))}
    points += [{"panel": "c_score", "model": MODELS[2], "x": xx, "value": v, "n": n}
               for xx, v, n in zip(r["x"], r["value"], r["n"])]
    frozen = pd.read_csv(METRICS)
    for m in MODELS:
        pred = models[m].predicted_label.to_numpy(int)
        assert np.array_equal(pred, (models[m].score.to_numpy() >= models[m].binary_threshold.to_numpy()).astype(int))
        expected = frozen[(frozen.comparator_id == m) & (frozen.population == "all_identity_clean")].balanced_accuracy.iloc[0]
        assert abs(bacc(y, pred) - expected) < 1e-14
        for name, mask in [("in", inside), ("out", ~inside)]:
            val = bacc(y[mask], pred[mask])
            if old is not None:
                assert abs(val - old["b"][f"{m}|{name}"]["bacc"]) < 1e-14
            D["b"][f"{m}|{name}"] = {"n": int(mask.sum()), "n_pos": int(y[mask].sum()),
                "n_neg": int((y[mask] == 0).sum()), "tp": int(((pred == 1) & (y == 1) & mask).sum()),
                "tn": int(((pred == 0) & (y == 0) & mask).sum()), "bacc": val}
            points.append({"panel": "b", "model": m, "x": name, "value": val, "n": int(mask.sum())})
    # Sample exact Murcko scaffold clusters from the full external set, with
    # replacement, preserving all members; use each draw across both bands and
    # every comparator. No molecule-level or experimental-source resampling.
    groups, codes = np.unique(base.bootstrap_cluster.to_numpy(str), return_inverse=True)
    ng = len(groups)
    cols, colnames = [], []
    for side, mask in [("in", inside), ("out", ~inside)]:
        for nm, values in [("n", mask), ("pos", mask & (y == 1)), ("neg", mask & (y == 0)), ("accepted", mask & accepted)]:
            cols.append(np.bincount(codes, weights=values.astype(int), minlength=ng))
            colnames.append(f"{side}_{nm}")
        for m in MODELS:
            pred = models[m].predicted_label.to_numpy(int)
            for nm, values in [("tp", mask & (y == 1) & (pred == 1)), ("tn", mask & (y == 0) & (pred == 0))]:
                cols.append(np.bincount(codes, weights=values.astype(int), minlength=ng))
                colnames.append(f"{side}_{m}_{nm}")
    counts = np.column_stack(cols)
    rng = np.random.default_rng(SEED)
    boot = np.empty((N_BOOT, len(cols)))
    for draw in range(N_BOOT):
        boot[draw] = counts[rng.integers(0, ng, size=ng)].sum(axis=0)
    B = {name: boot[:, i] for i, name in enumerate(colnames)}
    for m in MODELS:
        for side in ["in", "out"]:
            valid = (B[f"{side}_pos"] > 0) & (B[f"{side}_neg"] > 0)
            values = .5 * (B[f"{side}_{m}_tp"][valid] / B[f"{side}_pos"][valid]
                          + B[f"{side}_{m}_tn"][valid] / B[f"{side}_neg"][valid])
            lo, hi = np.quantile(values, [.025, .975])
            D["b"][f"{m}|{side}"].update(lo=float(lo), hi=float(hi), valid_draws=int(valid.sum()))
    delta = B["in_accepted"] / B["in_n"] - B["out_accepted"] / B["out_n"]
    lo, hi = np.quantile(delta, [.025, .975])
    D["band"] = {"n_in": int(inside.sum()), "n_out": int((~inside).sum()),
        "pos_in": float(y[inside].mean()), "pos_out": float(y[~inside].mean()),
        "acc_in": float(accepted[inside].mean()), "acc_out": float(accepted[~inside].mean()),
        "acc_diff": float(accepted[inside].mean() - accepted[~inside].mean()),
        "acc_lo": float(lo), "acc_hi": float(hi), "er_lo": float(2 * 10**(-.49)),
        "er_hi": float(2 * 10**.49), "n_clusters": ng}
    if old is not None:
        assert abs(D["band"]["acc_diff"] - old["band"]["acc_diff"]) < 1e-14
    D["provenance"] = {"analysis_date": "2026-10-18", "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": digest(SOURCE), "benchmark_input_sha256": digest(INPUT),
        "generator": str(Path(__file__).relative_to(ROOT)), "generator_sha256": digest(Path(__file__)),
        "bootstrap": {"unit": "exact Murcko scaffold group; acyclic fallback to connectivity block", "column": "bootstrap_cluster", "clusters": ng,
            "acyclic_fallback_entities": int((~nonempty).sum()),
            "replicates": N_BOOT, "seed": SEED, "ci": "2.5 and 97.5 percentile",
            "sampling": "full-set clusters with replacement; same draw for all models and both fixed bands",
            "single_class_handling": "exclude that BAcc replicate", "conditional_on_fixed_posthoc_band": True,
            "historical_interval_provenance": "unavailable; old intervals replaced, not reproduced"},
        "fixed_band": "absolute measured log10(ER/2) <= 0.49; descriptive, outcome-aware, not a deployment rule",
        "smoothing": "full-width 0.15 measured units or 0.05 score units; closed centered windows; measured-scale grid -1.20:0.02:1.60 with n >= 25; score grid -0.5:0.008:0.5 with n >= 20",
        "old_point_estimates": "all historical MAE, measured coverage, absolute-score coverage and BAcc points verified at atol 1e-14 in the 2026-10-18 cleanup audit",
        "changed_panel": "score axis changed from absolute to signed margin so asymmetric 0.35/0.70 window is shown correctly"}
    HANDOFF.write_text(json.dumps(D, ensure_ascii=False, indent=2) + "\n")
    pd.DataFrame(points).to_csv(OUT / "verified_plot_points.csv", index=False)
    report = {"status": "PASS", "old_handoff_sha256": digest(backup) if backup.exists() else None, "new_handoff_sha256": digest(HANDOFF),
        "point_rows": len(points), "n_entities": len(base), "n_clusters": ng,
        "new_coverage_delta": D["band"], "bootstrap": D["provenance"]["bootstrap"],
        "bacc": D["b"], "signed_score_points": len(D["c2"]["x"])}
    (OUT / "data_verification.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "n": len(base), "clusters": ng, "delta": D["band"]["acc_diff"],
                      "ci": [float(lo), float(hi)], "points": len(points)}, indent=2))


if __name__ == "__main__":
    main()

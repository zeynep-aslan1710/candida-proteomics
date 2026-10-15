"""Compare the two stresses with each other and with the published results.

1. Strength of the moderated t-test prior (d0) for each comparison.
2. Gre2 and Prn1: our fold changes next to the values in the paper.
3. A robust "core" response per stress: proteins significant at both doses
   in the same direction.
4. Shared vs stress-specific proteins, a fold-change scatter plot
   (H2O2 5 mM vs acetic acid 40 mM) and tables of the top hits.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import special

from differential import COMPARISONS, MARKERS, is_marker, trigamma_inverse

LOG2 = Path("data/processed/protein_log2.tsv")
DE_DIR = Path("data/processed")
FIG_DIR = Path("figures")

# Published log2 ratios (paper, Table 1) for the positive controls
PAPER = {("GRE2", "5mM"): 4.97, ("GRE2", "10mM"): 4.60,
         ("PRN1", "5mM"): 4.90, ("PRN1", "10mM"): 5.15}

COLOURS = {"both stresses": "#eb6834", "H2O2 only": "#2a78d6",
           "acetic acid only": "#1baf7a", "neither": "#c3c2b7"}


def prior_df(diffs: np.ndarray) -> float:
    """d0 of the moderated t-test: how many 'extra' degrees of freedom the
    shared variance adds to each protein (same estimate as differential.py)."""
    d = diffs.shape[1] - 1
    s2 = diffs.var(axis=1, ddof=1)
    s2 = np.maximum(s2, 1e-5 * np.median(s2))
    e = np.log(s2) - special.digamma(d / 2) + np.log(d / 2)
    extra_var = e.var(ddof=1) - special.polygamma(1, d / 2)
    return 2 * trigamma_inverse(extra_var) if extra_var > 0 else np.inf


def core(de: dict, a: str, b: str) -> pd.Series:
    """Direction ('up'/'down') of proteins significant in both a and b, same way."""
    x = de[a].set_index("Protein.Group")["direction"]
    y = de[b].set_index("Protein.Group")["direction"]
    both = x.index.intersection(y.index)
    x, y = x[both], y[both]
    return x[(x != "ns") & (x == y)]


def main() -> None:
    log2 = pd.read_csv(LOG2, sep="\t")
    names = [c for c in ("Protein.Names", "First.Protein.Description") if c in log2.columns]
    info = log2.set_index("Protein.Group")[["Genes"] + names]
    de = {name: pd.read_csv(DE_DIR / f"de_{name}.tsv", sep="\t")
          for name in list(COMPARISONS)[:4]}

    # 1) Prior strength
    print("Moderated t-test prior (d0 = extra degrees of freedom from all proteins)")
    for name in de:
        treated, control, reps = COMPARISONS[name]
        t = log2[[f"{treated}_R{r}" for r in reps]].to_numpy()
        c = log2[[f"{control}_R{r}" for r in reps]].to_numpy()
        diffs = (t - c)[~np.isnan(t - c).any(axis=1)]
        print(f"  {name:20s} own df = {len(reps) - 1}, d0 = {prior_df(diffs):5.1f}")

    # 2) Positive controls vs the paper
    print("\nPositive controls: this reanalysis vs the paper (log2 fold change)")
    for marker in MARKERS:
        for dose in ("5mM", "10mM"):
            res = de[f"H2O2_{dose}_vs_ctrl"]
            hit = res[is_marker(res["Genes"], marker)].iloc[0]
            print(f"  {marker} {dose:5s} ours {hit['log2FC']:5.2f} (q {hit['q']:.1e})"
                  f"   paper {PAPER[(marker, dose)]:5.2f}")

    # 3) Core responses and their overlap
    h2o2 = core(de, "H2O2_5mM_vs_ctrl", "H2O2_10mM_vs_ctrl")
    aa = core(de, "AA_40mM_vs_ctrl", "AA_60mM_vs_ctrl")
    shared = h2o2.index.intersection(aa.index)
    same_way = (h2o2[shared] == aa[shared]).sum()
    print(f"\nCore H2O2 response: {len(h2o2)} proteins "
          f"({(h2o2 == 'up').sum()} up, {(h2o2 == 'down').sum()} down)")
    print(f"Core acetic acid response: {len(aa)} proteins "
          f"({(aa == 'up').sum()} up, {(aa == 'down').sum()} down)")
    print(f"Shared: {len(shared)} proteins, {same_way} change in the same direction")

    # 4) Top hits per stress
    rows = []
    for stress, members, a, b in [("H2O2", h2o2, "H2O2_5mM_vs_ctrl", "H2O2_10mM_vs_ctrl"),
                                  ("acetic acid", aa, "AA_40mM_vs_ctrl", "AA_60mM_vs_ctrl")]:
        fc = (de[a].set_index("Protein.Group")["log2FC"]
              + de[b].set_index("Protein.Group")["log2FC"]) / 2
        table = info.loc[members.index].assign(stress=stress, direction=members,
                                               mean_log2FC=fc[members.index].round(2),
                                               in_other_stress=members.index.isin(shared))
        rows.append(table.sort_values("mean_log2FC", key=abs, ascending=False))
    top = pd.concat(rows).reset_index()
    top.to_csv(DE_DIR / "core_responses.tsv", sep="\t", index=False)
    show = ["stress", "Genes"] + names[:1] + ["mean_log2FC", "in_other_stress"]
    for stress in ("H2O2", "acetic acid"):
        for direction in ("up", "down"):
            sub = top[(top["stress"] == stress) & (top["direction"] == direction)]
            print(f"\nTop 10 {direction} with {stress}")
            print(sub[show].head(10).to_string(index=False))

    # 5) Fold-change scatter: H2O2 5 mM vs acetic acid 40 mM
    x = de["H2O2_5mM_vs_ctrl"].set_index("Protein.Group")
    y = de["AA_40mM_vs_ctrl"].set_index("Protein.Group")
    both = x.index.intersection(y.index)
    sx, sy = x.loc[both, "direction"] != "ns", y.loc[both, "direction"] != "ns"
    group = np.select([sx & sy, sx, sy], list(COLOURS)[:3], "neither")
    r = np.corrcoef(x.loc[both, "log2FC"], y.loc[both, "log2FC"])[0, 1]

    sns.set_theme(style="ticks", font_scale=0.8)
    fig, ax = plt.subplots(figsize=(6.5, 6))
    for label in ["neither", "acetic acid only", "H2O2 only", "both stresses"]:
        pts = group == label
        ax.scatter(x.loc[both, "log2FC"][pts], y.loc[both, "log2FC"][pts], s=8,
                   color=COLOURS[label], linewidth=0, alpha=0.8,
                   label=f"{label} ({pts.sum()})")
    for marker in MARKERS:
        hit = x.loc[both][is_marker(x.loc[both, "Genes"], marker)]
        if len(hit):
            px, py = hit["log2FC"].iloc[0], y.loc[hit.index[0], "log2FC"]
            ax.scatter(px, py, s=40, facecolor="none", edgecolor="black", linewidth=1)
            ax.annotate(marker.capitalize(), (px, py), xytext=(6, 4),
                        textcoords="offset points", fontsize=8, weight="bold")
    ax.axhline(0, color="0.6", linewidth=0.8)
    ax.axvline(0, color="0.6", linewidth=0.8)
    ax.set(xlabel="log2 fold change, H2O2 5 mM", ylabel="log2 fold change, acetic acid 40 mM",
           title=f"Two stresses, one proteome ({len(both)} proteins, r = {r:.2f})")
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    sns.despine(fig=fig)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stress_comparison.png", dpi=200)
    print(f"\nFold-change correlation H2O2 vs acetic acid: r = {r:.2f}")
    print("saved figures/stress_comparison.png and data/processed/core_responses.tsv")


if __name__ == "__main__":
    main()
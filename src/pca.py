"""PCA of the protein table: do the samples group by treatment?

Makes three PCA plots (all samples, H2O2 experiment, acetic acid experiment)
using only proteins quantified in every sample of the plot.
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from qc import COLOURS, condition, group_name

LOG2 = Path("data/processed/protein_log2.tsv")
SUMMARY = Path("data/processed/sample_summary.tsv")
FIG_DIR = Path("figures")
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # same colour-blind safe set


def run_pca(log2: pd.DataFrame, labels: list[str]):
    """PCA with samples as points and proteins as measurements (numpy SVD)."""
    data = log2[labels].dropna()                 # proteins found in every sample
    x = data.T.to_numpy()                        # one row per sample
    x = x - x.mean(axis=0)                       # centre: subtract each protein's mean
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    scores = pd.DataFrame(u[:, :2] * s[:2], index=labels, columns=["PC1", "PC2"])
    variance = s**2 / (s**2).sum() * 100         # % of variation per component
    return scores, variance[:2], len(data)


def main() -> None:
    log2 = pd.read_csv(LOG2, sep="\t")
    summary = pd.read_csv(SUMMARY, sep="\t")
    labels = [l for l in summary["label"] if "reinjection" not in l]

    cond = pd.Series({l: condition(l) for l in labels})                    # e.g. "AA treated"
    group = pd.Series({l: group_name(l).split("_", 1)[1] for l in labels})  # e.g. "40mM"
    panels = [
        ("All samples", labels, cond, list(COLOURS)),
        ("H2O2 experiment", [l for l in labels if l.startswith("H2O2")], group,
         ["ctrl", "5mM", "10mM"]),
        ("Acetic acid experiment", [l for l in labels if l.startswith("AA")], group,
         ["ctrl_40", "40mM", "ctrl_60", "60mM"]),
    ]

    sns.set_theme(style="ticks", font_scale=0.8)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax, (title, labs, groups, order) in zip(axes, panels):
        scores, var, n = run_pca(log2, labs)
        for colour, name in zip(PALETTE, order):
            pts = scores[groups[labs] == name]
            ax.scatter(pts["PC1"], pts["PC2"], s=70, color=colour,
                       edgecolor="white", linewidth=1.5, label=name)
            for label, row in pts.iterrows():   # write R1..R4 next to each point
                rep = re.search(r"_(R\d+)", label).group(1)
                ax.annotate(rep, (row["PC1"], row["PC2"]), xytext=(5, 4),
                            textcoords="offset points", fontsize=7, color="0.35")
        ax.set(title=f"{title} ({n} proteins)",
               xlabel=f"PC1 ({var[0]:.0f}% of variation)",
               ylabel=f"PC2 ({var[1]:.0f}% of variation)")
        ax.legend(title="", fontsize=7)
        print(f"{title}: {n} proteins, PC1 {var[0]:.1f}%, PC2 {var[1]:.1f}%")

    sns.despine(fig=fig)
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    fig.savefig(FIG_DIR / "pca.png", dpi=200)
    print("saved figures/pca.png")


if __name__ == "__main__":
    main()
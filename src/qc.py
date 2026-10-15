"""Quality control of the DIA-NN protein table.

Plots proteins per sample and intensity distributions, and prints the
coefficient of variation (CV) per group, data completeness and the
agreement between the AA40R3 measurement and its re-injection.
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

LOG2 = Path("data/processed/protein_log2.tsv")
SUMMARY = Path("data/processed/sample_summary.tsv")
FIG_DIR = Path("figures")

# Colour = experiment x control/treated (palette checked for colour blindness)
COLOURS = {
    "H2O2 control": "#2a78d6",
    "H2O2 treated": "#eb6834",
    "AA control": "#1baf7a",
    "AA treated": "#eda100",
}


def condition(label: str) -> str:
    """H2O2_ctrl_R1 -> 'H2O2 control', AA_40mM_R2 -> 'AA treated'."""
    experiment = label.split("_")[0]
    kind = "control" if "_ctrl" in label else "treated"
    return f"{experiment} {kind}"


def group_name(label: str) -> str:
    """H2O2_5mM_R3 -> H2O2_5mM (drop the replicate part)."""
    return re.sub(r"_R\d+.*$", "", label)


def save(fig: plt.Figure, name: str) -> None:
    sns.despine(fig=fig)
    fig.tight_layout()
    fig.savefig(FIG_DIR / name, dpi=200)
    print(f"saved figures/{name}")


def main() -> None:
    data = pd.read_csv(LOG2, sep="\t")
    summary = pd.read_csv(SUMMARY, sep="\t")
    labels = summary["label"].tolist()
    log2 = data[labels]
    summary["condition"] = summary["label"].map(condition)
    FIG_DIR.mkdir(exist_ok=True)
    sns.set_theme(style="ticks", font_scale=0.8)

    # 1) Proteins quantified per sample
    fig, ax = plt.subplots(figsize=(7, 8))
    sns.barplot(data=summary, y="label", x="proteins_quantified", hue="condition",
                palette=COLOURS, dodge=False, ax=ax)
    ax.set(xlabel="Protein groups quantified", ylabel="",
           title="Proteins quantified per sample")
    ax.legend(title="", loc="lower right")
    save(fig, "qc_proteins_per_sample.png")

    # 2) Intensity distribution of every sample
    long = log2.melt(var_name="label", value_name="log2").dropna()
    long["condition"] = long["label"].map(condition)
    fig, ax = plt.subplots(figsize=(7, 8))
    sns.boxplot(data=long, y="label", x="log2", hue="condition", palette=COLOURS,
                dodge=False, fliersize=1, linewidth=0.8, ax=ax)
    ax.set(xlabel="log2 protein quantity", ylabel="",
           title="Intensity distributions")
    ax.legend(title="", loc="lower right")
    save(fig, "qc_intensity_distributions.png")

    # 3) CV within each group (replicates should give similar values)
    main_labels = summary.loc[~summary["label"].str.contains("reinjection"), "label"]
    rows = []
    for group, cols in main_labels.groupby(main_labels.map(group_name)):
        linear = 2 ** log2[cols.tolist()]          # CV is computed on raw scale
        linear = linear.dropna(thresh=3)           # proteins seen in >= 3 replicates
        cv = linear.std(axis=1) / linear.mean(axis=1) * 100
        rows.append({"group": group, "replicates": len(cols),
                     "proteins": len(linear), "median_CV_%": round(cv.median(), 1)})
    print("\nCoefficient of variation per group")
    print(pd.DataFrame(rows).to_string(index=False))

    # 4) Completeness and technical reproducibility
    main = log2[main_labels]
    n = main.shape[1]
    found_in = main.notna().sum(axis=1)
    print(f"\nMissing values overall: {main.isna().mean().mean():.1%}")
    print(f"Proteins quantified in all {n} samples: {(found_in == n).sum()}")
    print(f"Proteins quantified in at least half: {(found_in >= n / 2).sum()}")

    pair = log2[["AA_40mM_R3", "AA_40mM_R3_reinjection"]].dropna()
    r = pair.corr().iloc[0, 1]
    print(f"\nAA40R3 vs re-injection: {len(pair)} shared proteins, Pearson r = {r:.3f}")


if __name__ == "__main__":
    main()
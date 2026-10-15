"""Differential abundance: which proteins change with H2O2 or acetic acid?

Each treated culture is compared with the control from the same replicate.
The paired differences are tested with a moderated t-test (the empirical
Bayes approach of limma): each protein's variance is shrunk towards the
typical variance of all proteins, which gives more power with few replicates.
P-values are corrected for testing thousands of proteins with the
Benjamini-Hochberg method. Results are saved per comparison and shown as
volcano plots.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import special, stats

LOG2 = Path("data/processed/protein_log2.tsv")
OUT_DIR = Path("data/processed")
FIG_DIR = Path("figures")

MIN_LOG2FC = 0.5   # at least ~1.4-fold change (same cut-off as the paper)
MAX_Q = 0.05       # at most 5 % false discoveries
UP, DOWN, NS = "#eb6834", "#2a78d6", "#c3c2b7"
MARKERS = ["GRE2", "PRN1"]  # positive controls from the paper

# name: (treated group, control group, replicates to pair)
COMPARISONS = {
    "H2O2_5mM_vs_ctrl": ("H2O2_5mM", "H2O2_ctrl", [1, 2, 3, 4]),
    "H2O2_10mM_vs_ctrl": ("H2O2_10mM", "H2O2_ctrl", [1, 2, 3, 4]),
    "AA_40mM_vs_ctrl": ("AA_40mM", "AA_ctrl_40", [1, 2, 3, 4]),
    "AA_60mM_vs_ctrl": ("AA_60mM", "AA_ctrl_60", [2, 3, 4]),  # CC60R1 is missing
    # sensitivity check: H2O2 without the outlying replicate 4 (see PCA)
    "H2O2_5mM_vs_ctrl_noR4": ("H2O2_5mM", "H2O2_ctrl", [1, 2, 3]),
    "H2O2_10mM_vs_ctrl_noR4": ("H2O2_10mM", "H2O2_ctrl", [1, 2, 3]),
}


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg: turn p-values into q-values (false discovery rate)."""
    n = len(p)
    order = np.argsort(p)                              # smallest p first
    q = p[order] * n / np.arange(1, n + 1)             # p x n / rank
    q = np.minimum.accumulate(q[::-1])[::-1]           # keep q non-decreasing
    out = np.empty(n)
    out[order] = np.minimum(q, 1)
    return out


def trigamma_inverse(x: float) -> float:
    """Solve trigamma(y) = x for y with Newton's method (as in limma)."""
    if x > 1e7:
        return 1 / np.sqrt(x)
    if x < 1e-6:
        return 1 / x
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        step = tri * (1 - tri / x) / special.polygamma(2, y)
        y += step
        if -step / y < 1e-8:
            break
    return y


def moderated_paired_t(diffs: np.ndarray) -> np.ndarray:
    """Moderated t-test on paired differences (one row per protein).

    1. Every protein gets its own variance s2 from its few differences.
    2. A prior variance s0 and its weight d0 are estimated from the
       distribution of all s2 values (limma's fitFDist).
    3. Each s2 is pulled towards s0: s2_post = (d0*s0 + d*s2) / (d0 + d).
    Returns two-sided p-values.
    """
    n = diffs.shape[1]
    d = n - 1                                          # degrees of freedom per protein
    mean = diffs.mean(axis=1)
    s2 = diffs.var(axis=1, ddof=1)
    s2 = np.maximum(s2, 1e-5 * np.median(s2))          # avoid log(0)

    e = np.log(s2) - special.digamma(d / 2) + np.log(d / 2)
    extra_var = e.var(ddof=1) - special.polygamma(1, d / 2)
    if extra_var > 0:
        d0 = 2 * trigamma_inverse(extra_var)
        s0 = np.exp(e.mean() + special.digamma(d0 / 2) - np.log(d0 / 2))
        s2_post = (d0 * s0 + d * s2) / (d0 + d)
    else:                                              # all proteins equally variable
        d0, s2_post = np.inf, np.full_like(s2, np.exp(e.mean()))

    t = mean / np.sqrt(s2_post / n)
    return 2 * stats.t.sf(np.abs(t), d + d0)


def compare(data: pd.DataFrame, treated: str, control: str, reps: list[int]) -> pd.DataFrame:
    """Moderated paired t-test for every protein quantified in all pairs."""
    t_cols = [f"{treated}_R{r}" for r in reps]
    c_cols = [f"{control}_R{r}" for r in reps]
    complete = data[t_cols + c_cols].notna().all(axis=1)
    sub = data[complete]
    t_vals, c_vals = sub[t_cols].to_numpy(), sub[c_cols].to_numpy()

    res = sub[["Protein.Group", "Genes"]].copy()
    res["log2FC"] = (t_vals - c_vals).mean(axis=1)     # mean of paired differences
    res["p"] = moderated_paired_t(t_vals - c_vals)
    res["q"] = bh_adjust(res["p"].to_numpy())
    res["n_pairs"] = len(reps)
    big = res["log2FC"].abs() > MIN_LOG2FC
    res["direction"] = np.where(big & (res["q"] < MAX_Q),
                                np.where(res["log2FC"] > 0, "up", "down"), "ns")
    return res.sort_values("p")


def is_marker(genes: pd.Series, name: str) -> pd.Series:
    """Exact gene-name match (GRE2 must not match e.g. GRE22)."""
    return genes.fillna("").str.upper().str.split(";").apply(lambda g: name in g)


def volcano(ax: plt.Axes, res: pd.DataFrame, title: str) -> None:
    y = -np.log10(res["p"])
    for direction, colour in [("ns", NS), ("down", DOWN), ("up", UP)]:
        pts = res["direction"] == direction
        n = pts.sum()
        label = {"ns": f"not significant ({n})", "up": f"up ({n})", "down": f"down ({n})"}
        ax.scatter(res.loc[pts, "log2FC"], y[pts], s=8, color=colour,
                   linewidth=0, alpha=0.8, label=label[direction])
    for name in MARKERS:                                # label Gre2 and Prn1
        hit = res[is_marker(res["Genes"], name)]
        if len(hit):
            x0, y0 = hit["log2FC"].iloc[0], -np.log10(hit["p"].iloc[0])
            ax.scatter(x0, y0, s=40, facecolor="none", edgecolor="black", linewidth=1)
            ax.annotate(name.capitalize(), (x0, y0), xytext=(6, 4),
                        textcoords="offset points", fontsize=8, weight="bold")
    for x in (-MIN_LOG2FC, MIN_LOG2FC):
        ax.axvline(x, color="0.6", linewidth=0.8, linestyle="--")
    ax.set(title=title, xlabel="log2 fold change (treated / control)",
           ylabel="-log10 p-value")
    ax.legend(fontsize=7, loc="upper left", frameon=False)


def main() -> None:
    data = pd.read_csv(LOG2, sep="\t")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(exist_ok=True)

    results, rows = {}, []
    for name, (treated, control, reps) in COMPARISONS.items():
        res = compare(data, treated, control, reps)
        res.to_csv(OUT_DIR / f"de_{name}.tsv", sep="\t", index=False)
        results[name] = res
        rows.append({"comparison": name, "pairs": len(reps), "tested": len(res),
                     "p<0.05": (res["p"] < 0.05).sum(),
                     "up": (res["direction"] == "up").sum(),
                     "down": (res["direction"] == "down").sum()})
    print("Significant proteins (q < 0.05 and |log2FC| > 0.5)")
    print(pd.DataFrame(rows).to_string(index=False))

    print("\nPositive controls")
    for name, res in results.items():
        for marker in MARKERS:
            hit = res[is_marker(res["Genes"], marker)]
            if len(hit):
                h = hit.iloc[0]
                print(f"{name:24s} {marker}: log2FC {h['log2FC']:5.2f}  "
                      f"p {h['p']:.1e}  q {h['q']:.1e}  {h['direction']}")
            else:
                print(f"{name:24s} {marker}: not quantified in all pairs")

    # Overlap between the full and the no-R4 analysis
    for dose in ("5mM", "10mM"):
        full = results[f"H2O2_{dose}_vs_ctrl"]
        no_r4 = results[f"H2O2_{dose}_vs_ctrl_noR4"]
        a = set(full.loc[full["direction"] != "ns", "Protein.Group"])
        b = set(no_r4.loc[no_r4["direction"] != "ns", "Protein.Group"])
        print(f"\nH2O2 {dose}: {len(a)} significant with R4, {len(b)} without R4, "
              f"{len(a & b)} in both")

    sns.set_theme(style="ticks", font_scale=0.8)
    main_four = list(COMPARISONS)[:4]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for ax, name in zip(axes.flat, main_four):
        volcano(ax, results[name], name.replace("_", " "))
    sns.despine(fig=fig)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "volcano.png", dpi=200)
    print("\nsaved figures/volcano.png")


if __name__ == "__main__":
    main()
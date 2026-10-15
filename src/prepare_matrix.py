"""Combine the DIA-NN protein table with the sample annotation.

Reads the DIA-NN protein group matrix, renames the run columns to readable
sample labels (e.g. H2O2_5mM_R1), log2-transforms the quantities and saves
the result together with a per-sample summary.
"""
from pathlib import Path, PureWindowsPath

import numpy as np
import pandas as pd

PG_MATRIX = Path("C:/proteomics_data/diann_output/report.pg_matrix.tsv")
ANNOTATION = Path("data/sample_annotation.tsv")
OUT_DIR = Path("data/processed")


def make_label(row: pd.Series) -> str:
    """Build a readable label, e.g. H2O2_5mM_R1 or AA_ctrl_40_R2."""
    label = f"{row['experiment']}_{row['group']}_R{row['replicate']}"
    if isinstance(row["note"], str) and row["note"]:
        label += "_reinjection"
    return label


def main() -> None:
    pg = pd.read_csv(PG_MATRIX, sep="\t")
    ann = pd.read_csv(ANNOTATION, sep="\t", dtype={"replicate": str})
    ann["label"] = ann.apply(make_label, axis=1)
    file_to_label = dict(zip(ann["file"], ann["label"]))

    # Run columns are the full paths of the .raw files; the rest is protein info
    run_cols = [c for c in pg.columns if c.lower().endswith(".raw")]
    info = pg[[c for c in pg.columns if c not in run_cols]]
    quant = pg[run_cols].rename(columns=lambda c: file_to_label[PureWindowsPath(c).name])

    # Missing values can appear as empty cells or zeros; log2 needs NaN for both
    log2 = np.log2(quant.replace(0, np.nan))
    log2 = log2[ann["label"]]  # same order as the annotation table

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.concat([info, log2], axis=1).to_csv(OUT_DIR / "protein_log2.tsv", sep="\t", index=False)

    summary = ann[["label", "file", "sample_name"]].copy()
    summary["proteins_quantified"] = log2.notna().sum().values
    summary["median_log2"] = log2.median().round(2).values
    summary.to_csv(OUT_DIR / "sample_summary.tsv", sep="\t", index=False)

    print(f"{len(log2)} protein groups x {log2.shape[1]} samples\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
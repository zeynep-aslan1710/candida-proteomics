"""First look at the DIA-NN protein table: how many proteins, and Gre2/Prn1."""
import pandas as pd

PG_MATRIX = "C:/proteomics_data/diann_output/report.pg_matrix.tsv"

pg = pd.read_csv(PG_MATRIX, sep="\t")
print(f"{len(pg)} protein groups quantified\n")

# Find Gre2 and Prn1 by gene or protein name
names = pg["Genes"].fillna("") + " " + pg["Protein.Names"].fillna("")
hits = pg[names.str.contains("GRE2|PRN1", case=False)]

# Show replicate 1: control (S1704_143) vs 5 mM H2O2 (S1704_144)
run_cols = [c for c in pg.columns if "S1704_143" in c or "S1704_144" in c]
table = hits[["Genes"] + run_cols].copy()
table.columns = ["Gene", "ctrl R1 (143)", "5mM H2O2 R1 (144)"]
print(table.to_string(index=False))
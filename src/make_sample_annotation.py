"""Build a sample annotation table from the .raw file headers.

Reads data/raw_headers.txt (created by read_raw_headers.py), takes the
sample name stored in each file header and decodes it into experiment,
group and replicate. Writes the result to data/sample_annotation.tsv.

Sample name patterns used by the lab:
    CH2O2R1, 5mMH2O2R1, 10mMH2O2R1   -> hydrogen peroxide experiment
    CC40R1, AA40R1, CC60R2, AA60R2   -> acetic acid experiment
"""
import csv
import re
from pathlib import Path

HEADERS_FILE = Path("data/raw_headers.txt")
OUT_FILE = Path("data/sample_annotation.tsv")

H2O2_PATTERN = re.compile(r"^(C|5mM|10mM)H2O2R(\d+)$")
AA_PATTERN = re.compile(r"^(CC|AA)(40|60)R(\d+)$")


def read_sample_names(path: Path) -> dict[str, str]:
    """Return {raw file name: sample name}. The sample name is the line after 'r_System'."""
    lines = path.read_text(encoding="utf-8").splitlines()
    sample_names = {}
    current_file = None
    for i, line in enumerate(lines):
        header = re.match(r"^=== (.+\.raw) ===$", line)
        if header:
            current_file = header.group(1)
        elif line == "r_System" and current_file and i + 1 < len(lines):
            sample_names[current_file] = lines[i + 1].strip()
    return sample_names


def decode_sample_name(name: str) -> dict[str, str]:
    """Turn a lab sample name into experiment, group and replicate."""
    match = H2O2_PATTERN.match(name)
    if match:
        dose, replicate = match.groups()
        group = "ctrl" if dose == "C" else dose
        return {"experiment": "H2O2", "group": group, "replicate": replicate}

    match = AA_PATTERN.match(name)
    if match:
        kind, dose, replicate = match.groups()
        group = f"ctrl_{dose}" if kind == "CC" else f"{dose}mM"
        return {"experiment": "AA", "group": group, "replicate": replicate}

    return {"experiment": "unknown", "group": "unknown", "replicate": ""}


def main() -> None:
    sample_names = read_sample_names(HEADERS_FILE)

    rows = []
    seen = set()
    for file_name in sorted(sample_names):
        name = sample_names[file_name]
        info = decode_sample_name(name)
        note = "repeated injection of an earlier sample" if name in seen else ""
        seen.add(name)
        rows.append({"file": file_name, "sample_name": name, **info, "note": note})

    with OUT_FILE.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["file", "sample_name", "experiment", "group", "replicate", "note"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(f"{row['file']:<26} {row['sample_name']:<12} {row['experiment']:<6} "
              f"{row['group']:<8} R{row['replicate']}  {row['note']}")
    print(f"\n{len(rows)} samples written to {OUT_FILE}")


if __name__ == "__main__":
    main()
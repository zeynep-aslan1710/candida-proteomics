"""Read sample metadata from the headers of Thermo .raw files.

Only the first part of each file is downloaded, so the full
4 GB files are not needed to find out which sample each run is.
"""
import re
import urllib.request
from pathlib import Path

URL_FILE = Path("raw_urls.txt")
OUT_FILE = Path("data/raw_headers.txt")
N_BYTES = 200_000          # the header is at the start of the file
MAX_STRINGS = 40           # readable text pieces to keep per file

# Thermo stores header text as UTF-16 (each letter followed by a zero byte)
UTF16_TEXT = re.compile(rb"(?:[\x20-\x7e]\x00){3,}")


def read_header_strings(url: str) -> list[str]:
    """Download the first bytes of a file and return the readable text in it."""
    request = urllib.request.Request(url, headers={"Range": f"bytes=0-{N_BYTES - 1}"})
    with urllib.request.urlopen(request, timeout=120) as response:
        chunk = response.read(N_BYTES)
    strings = []
    for match in UTF16_TEXT.finditer(chunk):
        text = match.group().decode("utf-16-le").strip()
        if text and text not in strings:
            strings.append(text)
    return strings[:MAX_STRINGS]


def main() -> None:
    urls = [line.strip() for line in URL_FILE.read_text(encoding="utf-8").splitlines()
            if line.strip().endswith(".raw")]
    print(f"{len(urls)} raw files found")

    lines = []
    for url in urls:
        name = url.rsplit("/", 1)[-1]
        strings = read_header_strings(url)
        print(f"done: {name}")
        lines.append(f"=== {name} ===")
        lines.extend(strings)
        lines.append("")

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FILE.write_text("\n".join(lines), encoding="utf-8")
    print(f"Saved to {OUT_FILE}")


if __name__ == "__main__":
    main()
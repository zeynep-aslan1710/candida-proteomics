"""Download the raw files of PRIDE dataset PXD020283.

The files (~85 GB) are stored outside the project folder, because they are
too large for GitHub and OneDrive. Partly downloaded files are resumed and
complete files are skipped, so the script can simply be run again.
"""
import time
import urllib.request
from pathlib import Path

URL_FILE = Path("raw_urls.txt")
TARGET_DIR = Path("C:/proteomics_data/raw")
CHUNK = 8 * 1024 * 1024      # read 8 MB at a time
MAX_TRIES = 5


def remote_size(url: str) -> int:
    """Ask the server how large the file is, without downloading it."""
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=120) as response:
        return int(response.headers["Content-Length"])


def download(url: str, target: Path) -> None:
    """Download one file, resuming a partial download if there is one."""
    total = remote_size(url)
    done = target.stat().st_size if target.exists() else 0
    if done == total:
        print(f"{target.name}: already complete, skipped")
        return
    if done > total:
        done = 0

    headers = {"Range": f"bytes={done}-"} if done else {}
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=300) as response:
        resumed = done > 0 and response.status == 206
        if not resumed:
            done = 0
        with target.open("ab" if resumed else "wb") as f:
            while chunk := response.read(CHUNK):
                f.write(chunk)
                done += len(chunk)
                print(f"\r{target.name}: {done / 1e9:.2f} / {total / 1e9:.2f} GB",
                      end="", flush=True)
    print()
    if done != total:
        raise IOError(f"incomplete download ({done} of {total} bytes)")


def main() -> None:
    TARGET_DIR.mkdir(parents=True, exist_ok=True)
    urls = [line.strip() for line in URL_FILE.read_text(encoding="utf-8").splitlines()
            if line.strip().endswith(".raw")]
    print(f"{len(urls)} raw files to download into {TARGET_DIR}\n")

    failed = []
    for i, url in enumerate(urls, start=1):
        target = TARGET_DIR / url.rsplit("/", 1)[-1]
        print(f"[{i}/{len(urls)}]")
        for attempt in range(1, MAX_TRIES + 1):
            try:
                download(url, target)
                break
            except Exception as error:
                print(f"\n  attempt {attempt} failed: {error}")
                time.sleep(30)
        else:
            failed.append(target.name)

    print("\nFinished.")
    if failed:
        print("These files could not be downloaded, run the script again:")
        for name in failed:
            print(f"  {name}")


if __name__ == "__main__":
    main()
"""
Download aligned TCGA data via gdc-client.

Reads data types from config.toml and downloads each dtype's aligned manifest.

Usage:
    uv run python scripts/02_download.py
    uv run python scripts/02_download.py --types expression mutation
    uv run python scripts/02_download.py --gdc-client /path/to/gdc-client
    uv run python scripts/02_download.py --token ~/tcga-token.txt
"""

import argparse
import shutil
import subprocess
import tomllib
from pathlib import Path

CONFIG_PATH = Path("config.toml")
MANIFEST_DIR = Path("data/manifests/aligned")
RAW_DIR = Path("data/raw")
DEFAULT_TOKEN = Path.home() / "tcga-token.txt"


def load_config():
    with open(CONFIG_PATH, "rb") as f:
        return tomllib.load(f)


def find_gdc_client(override: str = None) -> str:
    if override:
        return override
    found = shutil.which("gdc-client")
    if found:
        return found
    # common install locations
    for candidate in [
        Path.home() / ".local/bin/gdc-client",
        Path("/usr/local/bin/gdc-client"),
    ]:
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError(
        "gdc-client not found. Install it or pass --gdc-client /path/to/gdc-client"
    )


def main():
    parser = argparse.ArgumentParser(description="Download aligned TCGA data via gdc-client")
    parser.add_argument("--types", nargs="+", help="Data types to download (default: all from config)")
    parser.add_argument("--gdc-client", dest="gdc_client", help="Path to gdc-client binary")
    parser.add_argument("--token", default=str(DEFAULT_TOKEN), help="GDC token file path")
    parser.add_argument("--threads", type=int, default=8, help="Download threads (default: 8)")
    parser.add_argument("--retries", type=int, default=5, help="Retry attempts (default: 5)")
    args = parser.parse_args()

    cfg = load_config()
    dtypes = args.types or list(cfg["data_types"].keys())
    gdc_client = find_gdc_client(args.gdc_client)
    token_path = Path(args.token)

    print(f"gdc-client : {gdc_client}")
    print(f"token      : {token_path} ({'found' if token_path.exists() else 'NOT FOUND'})")
    print(f"data types : {dtypes}\n")

    for dtype in dtypes:
        manifest = MANIFEST_DIR / f"manifest_{dtype}.tsv"
        if not manifest.exists():
            print(f"[{dtype}] manifest not found: {manifest} — skipping")
            continue

        out_dir = RAW_DIR / dtype
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"{'='*50}")
        print(f"  Downloading: {dtype}")
        print(f"{'='*50}")

        cmd = [
            gdc_client, "download",
            "-m", str(manifest),
            "-d", str(out_dir),
            "-n", str(args.threads),
            "--retry-amount", str(args.retries),
        ]
        if token_path.exists():
            cmd += ["-t", str(token_path)]

        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"  WARNING: gdc-client exited with code {result.returncode} for {dtype}")
        else:
            print(f"  Done: {dtype}\n")

    print("All downloads complete!")
    for d in sorted(RAW_DIR.iterdir()):
        if d.is_dir():
            size = sum(f.stat().st_size for f in d.rglob("*") if f.is_file())
            print(f"  {d.name}: {size / 1e9:.2f} GB")


if __name__ == "__main__":
    main()

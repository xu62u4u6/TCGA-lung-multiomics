"""
Download TCGA data with resume support and progress tracking.

Features:
- TSV-based download log for resume (skips already-downloaded files)
- Auto-detects files from previous gdc-client runs
- MD5 verification
- Progress bar (tqdm)

Usage:
    uv run python scripts/download_with_resume.py --type cnv
    uv run python scripts/download_with_resume.py --type all
    uv run python scripts/download_with_resume.py --status
"""

import argparse
import hashlib
import io
import tarfile
import time
from pathlib import Path

import requests
from tqdm import tqdm

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
MANIFEST_DIR = DATA_DIR / "manifests"
RAW_DIR = DATA_DIR / "raw"
LOG_PATH = DATA_DIR / "download_log.tsv"

GDC_DATA = "https://api.gdc.cancer.gov/data"
DATA_TYPES = ["expression", "mutation", "cnv", "cnv_gene", "methylation", "clinical"]
BATCH_SIZE = 50
MAX_RETRIES = 3


def load_manifest(data_type: str) -> list[dict]:
    path = MANIFEST_DIR / f"manifest_{data_type}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")
    rows = []
    with open(path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                rows.append({
                    "file_id": parts[0],
                    "filename": parts[1],
                    "md5": parts[2],
                    "size": int(parts[3]),
                })
    return rows


def load_done_ids() -> set[str]:
    """Load completed file_ids from download log."""
    done = set()
    if LOG_PATH.exists():
        with open(LOG_PATH) as f:
            f.readline()  # skip header
            for line in f:
                parts = line.strip().split("\t")
                if parts:
                    done.add(parts[0])
    return done


def append_log(file_id: str, data_type: str, filename: str, status: str):
    """Append one entry to the download log."""
    if not LOG_PATH.exists():
        LOG_PATH.write_text("file_id\tdata_type\tfilename\tstatus\ttime\n")
    with open(LOG_PATH, "a") as f:
        f.write(f"{file_id}\t{data_type}\t{filename}\t{status}\t{time.strftime('%Y-%m-%d %H:%M:%S')}\n")


def scan_existing(data_type: str, manifest: list[dict]) -> set[str]:
    """Detect files already on disk from previous gdc-client downloads."""
    out_dir = RAW_DIR / data_type
    if not out_dir.exists():
        return set()

    on_disk = set()
    existing_dirs = {d.name for d in out_dir.iterdir() if d.is_dir()}
    existing_files = {f.name for f in out_dir.iterdir() if f.is_file()}

    for r in manifest:
        if r["file_id"] in existing_dirs or r["filename"] in existing_files:
            on_disk.add(r["file_id"])
    return on_disk


def md5_check(filepath: Path, expected: str) -> bool:
    if not expected:
        return True
    h = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest() == expected


def download_batch(file_ids: list[str], output_dir: Path, token: str = None) -> dict[str, Path]:
    headers = {"Content-Type": "application/json"}
    if token:
        headers["X-Auth-Token"] = token

    r = requests.post(GDC_DATA, json={"ids": file_ids}, headers=headers, stream=True)
    r.raise_for_status()

    results = {}
    content_type = r.headers.get("Content-Type", "")

    if len(file_ids) == 1 and "application/x-tar" not in content_type:
        cd = r.headers.get("Content-Disposition", "")
        fname = cd.split("filename=")[-1].strip('"') if "filename=" in cd else f"{file_ids[0]}.dat"
        out = output_dir / fname
        out.write_bytes(r.content)
        results[file_ids[0]] = out
    else:
        buf = io.BytesIO(r.content)
        with tarfile.open(fileobj=buf) as tar:
            for member in tar.getmembers():
                if not member.isfile() or member.name.endswith("MANIFEST.txt"):
                    continue
                fname = Path(member.name).name
                parent = Path(member.name).parent.name
                f = tar.extractfile(member)
                if f:
                    out = output_dir / fname
                    out.write_bytes(f.read())
                    results[parent] = out

    return results


def download_data_type(data_type: str, done_ids: set[str], token: str = None):
    manifest = load_manifest(data_type)

    # Scan disk for files from previous gdc-client runs, log them
    on_disk = scan_existing(data_type, manifest) - done_ids
    if on_disk:
        id_to_row = {r["file_id"]: r for r in manifest}
        for fid in on_disk:
            append_log(fid, data_type, id_to_row[fid]["filename"], "done")
            done_ids.add(fid)
        tqdm.write(f"  found {len(on_disk)} existing files on disk")

    pending = [r for r in manifest if r["file_id"] not in done_ids]
    if not pending:
        print(f"  {data_type}: all done")
        return

    total_size = sum(r["size"] for r in pending)
    print(f"  {len(pending)} files remaining ({total_size / 1e9:.2f} GB)")

    out_dir = RAW_DIR / data_type
    out_dir.mkdir(parents=True, exist_ok=True)

    id_to_row = {r["file_id"]: r for r in pending}
    all_ids = [r["file_id"] for r in pending]

    pbar = tqdm(total=len(all_ids), desc=data_type, unit="file")

    for i in range(0, len(all_ids), BATCH_SIZE):
        batch = all_ids[i: i + BATCH_SIZE]
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                results = download_batch(batch, out_dir, token=token)
                for fid in batch:
                    row = id_to_row[fid]
                    filepath = results.get(fid)
                    if filepath and filepath.exists():
                        status = "done" if md5_check(filepath, row["md5"]) else "md5_fail"
                    elif (out_dir / row["filename"]).exists():
                        status = "done"
                    else:
                        status = "error"
                    append_log(fid, data_type, row["filename"], status)
                    done_ids.add(fid)
                pbar.update(len(batch))
                break
            except Exception as e:
                if attempt == MAX_RETRIES:
                    for fid in batch:
                        append_log(fid, data_type, id_to_row[fid]["filename"], "error")
                        done_ids.add(fid)
                    pbar.update(len(batch))
                    tqdm.write(f"  batch failed after {MAX_RETRIES} retries: {e}")
                else:
                    tqdm.write(f"  retry {attempt}/{MAX_RETRIES}: {e}")
                    time.sleep(2 ** attempt)

    pbar.close()


def show_status():
    done_ids = load_done_ids()
    print("\n Download Progress")
    print("=" * 65)
    for dt in DATA_TYPES:
        try:
            manifest = load_manifest(dt)
        except FileNotFoundError:
            print(f"  {dt:15s}  no manifest")
            continue

        total = len(manifest)
        on_disk = scan_existing(dt, manifest)
        done_in_log = {r["file_id"] for r in manifest if r["file_id"] in done_ids}
        done = len(on_disk | done_in_log)
        remaining = total - done

        pct = done / total * 100 if total else 0
        bar = "█" * int(pct // 5) + "░" * (20 - int(pct // 5))
        print(f"  {dt:15s} {bar} {pct:5.1f}%  {done}/{total}  ({remaining} remaining)")
    print()


def main():
    parser = argparse.ArgumentParser(description="Download TCGA data with resume support")
    parser.add_argument("--type", type=str, default="all",
                        help="expression, mutation, cnv, methylation, clinical, or all")
    parser.add_argument("--token", type=str, default=str(Path.home() / "tcga-token.txt"))
    parser.add_argument("--status", action="store_true", help="Show download progress only")
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    args = parser.parse_args()

    if args.status:
        show_status()
        return

    types = DATA_TYPES if args.type == "all" else [args.type]

    token = None
    token_path = Path(args.token)
    if token_path.exists():
        token = token_path.read_text().strip()
        print(f"Using token: {token_path}")
    else:
        print("No token file found, downloading without auth")

    done_ids = load_done_ids()

    for dt in types:
        print(f"\n--- {dt} ---")
        download_data_type(dt, done_ids, token=token)

    show_status()


if __name__ == "__main__":
    main()

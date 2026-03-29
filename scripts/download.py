"""
Download aligned TCGA data via gdc-client with resume and MD5 verification.

Python manages resume logic (skip completed files), calls gdc-client for actual
downloading, then verifies MD5 and records results to data/download_log.tsv.

Usage:
    uv run python scripts/02_download.py                      # download all dtypes
    uv run python scripts/02_download.py --types expression   # specific dtype
    uv run python scripts/02_download.py --status             # show progress
    uv run python scripts/02_download.py --verify             # re-verify MD5 only
"""

import argparse
import hashlib
import shutil
import subprocess
import tempfile
import time
import tomllib
from pathlib import Path

import pandas as pd

CONFIG_PATH    = Path("config.toml")
MANIFEST_DIR   = Path("data/manifests/aligned")
RAW_DIR        = Path("data/raw")
LOG_PATH       = Path("data/download_log.tsv")


# ── Config ───────────────────────────────────────────────────────────────────

def load_config() -> dict:
    with open(CONFIG_PATH, "rb") as f:
        return tomllib.load(f)


def find_gdc_client(cfg: dict) -> str:
    # 1. config [download].gdc_client
    configured = cfg.get("download", {}).get("gdc_client")
    if configured:
        p = Path(configured)
        if p.exists():
            return str(p)
    # 2. PATH
    found = shutil.which("gdc-client")
    if found:
        return found
    # 3. common locations
    for candidate in [Path("tools/gdc-client"), Path.home() / ".local/bin/gdc-client"]:
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError(
        "gdc-client not found. Download it to tools/gdc-client or set [download].gdc_client in config.toml"
    )


# ── Manifest ─────────────────────────────────────────────────────────────────

def load_manifest(dtype: str) -> pd.DataFrame:
    path = MANIFEST_DIR / f"manifest_{dtype}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Aligned manifest not found: {path}\nRun 01_align_manifests.py first.")
    return pd.read_csv(path, sep="\t").rename(columns={"id": "file_id"})


def write_temp_manifest(df: pd.DataFrame, path: Path):
    with open(path, "w") as f:
        f.write("id\tfilename\tmd5\tsize\tstate\n")
        for _, row in df.iterrows():
            f.write(f"{row['file_id']}\t{row['filename']}\t{row.get('md5','')}\t{row.get('size',0)}\t{row.get('state','')}\n")


# ── Log ──────────────────────────────────────────────────────────────────────

def load_log() -> pd.DataFrame:
    if LOG_PATH.exists():
        return pd.read_csv(LOG_PATH, sep="\t")
    return pd.DataFrame(columns=["file_id", "dtype", "filename", "md5_ok", "timestamp"])


def save_log(log_df: pd.DataFrame):
    log_df.to_csv(LOG_PATH, sep="\t", index=False)


def append_log(log_df: pd.DataFrame, records: list[dict]) -> pd.DataFrame:
    new_rows = pd.DataFrame(records)
    return pd.concat([log_df, new_rows], ignore_index=True)


# ── Disk scan ────────────────────────────────────────────────────────────────

def scan_existing(dtype: str, manifest: pd.DataFrame) -> set[str]:
    """Find file_ids already on disk from previous downloads."""
    out_dir = RAW_DIR / dtype
    if not out_dir.exists():
        return set()

    existing_dirs  = {d.name for d in out_dir.iterdir() if d.is_dir()}
    existing_files = {f.name for f in out_dir.iterdir() if f.is_file()}

    found = set()
    for _, row in manifest.iterrows():
        if row["file_id"] in existing_dirs or row["filename"] in existing_files:
            found.add(row["file_id"])
    return found


def find_downloaded_file(file_id: str, filename: str, dtype: str) -> Path | None:
    out_dir = RAW_DIR / dtype
    # gdc-client stores as {out_dir}/{file_id}/{filename}
    by_uuid = out_dir / file_id / filename
    if by_uuid.exists():
        return by_uuid
    # flat layout fallback
    flat = out_dir / filename
    if flat.exists():
        return flat
    return None


# ── MD5 ──────────────────────────────────────────────────────────────────────

def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_files(dtype: str, manifest: pd.DataFrame, file_ids: set[str]) -> list[dict]:
    records = []
    for _, row in manifest[manifest["file_id"].isin(file_ids)].iterrows():
        fpath = find_downloaded_file(row["file_id"], row["filename"], dtype)
        if fpath is None:
            ok = False
        elif row.get("md5"):
            ok = md5sum(fpath) == row["md5"]
        else:
            ok = True  # no md5 in manifest, assume ok
        records.append({
            "file_id":   row["file_id"],
            "dtype":     dtype,
            "filename":  row["filename"],
            "md5_ok":    ok,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        })
    return records


# ── gdc-client ───────────────────────────────────────────────────────────────

def run_gdc_client(manifest_path: Path, out_dir: Path, token: str | None,
                   threads: int, retries: int, cfg: dict) -> int:
    cmd = [
        find_gdc_client(cfg),
        "download",
        "-m", str(manifest_path),
        "-d", str(out_dir),
        "-n", str(threads),
        "--retry-amount", str(retries),
    ]
    if token and Path(token).exists():
        cmd += ["-t", token]
    return subprocess.run(cmd).returncode





# ── Download ─────────────────────────────────────────────────────────────────

def download_dtype(dtype: str, cfg: dict, log_df: pd.DataFrame) -> pd.DataFrame:
    dl = cfg.get("download", {})
    token   = dl.get("token")
    threads = dl.get("threads", 8)
    retries = dl.get("retries", 5)

    manifest = load_manifest(dtype)
    done_ids = set(log_df[log_df["dtype"] == dtype]["file_id"])

    # also count files on disk not yet in log
    on_disk = scan_existing(dtype, manifest) - done_ids
    if on_disk:
        print(f"  Found {len(on_disk)} existing files not in log — verifying...")
        new_records = verify_files(dtype, manifest, on_disk)
        log_df = append_log(log_df, new_records)
        save_log(log_df)
        done_ids |= on_disk

    remaining = manifest[~manifest["file_id"].isin(done_ids)]
    total = len(manifest)
    print(f"  {dtype}: {total - len(remaining)}/{total} done, {len(remaining)} remaining")

    if remaining.empty:
        return log_df

    out_dir = RAW_DIR / dtype
    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False, mode="w") as tmp:
        tmp_path = Path(tmp.name)
    write_temp_manifest(remaining, tmp_path)

    print(f"  Calling gdc-client ({len(remaining)} files, {threads} threads)...")
    rc = run_gdc_client(tmp_path, out_dir, token, threads, retries, cfg)
    tmp_path.unlink(missing_ok=True)

    if rc != 0:
        print(f"  WARNING: gdc-client exited with code {rc}")

    # verify newly downloaded files
    newly_done = scan_existing(dtype, manifest) - done_ids
    if newly_done:
        print(f"  Verifying {len(newly_done)} downloaded files...")
        new_records = verify_files(dtype, manifest, newly_done)
        failed = [r for r in new_records if not r["md5_ok"]]
        if failed:
            print(f"  WARNING: {len(failed)} MD5 failures: {[r['filename'] for r in failed]}")
        log_df = append_log(log_df, new_records)
        save_log(log_df)

    return log_df


# ── Status ───────────────────────────────────────────────────────────────────

def show_status(dtypes: list[str]):
    log_df = load_log()
    print(f"\n{'dtype':<15} {'progress':>22}  {'md5_ok':>8}  {'md5_fail':>8}")
    print("-" * 60)
    for dtype in dtypes:
        try:
            manifest = load_manifest(dtype)
        except FileNotFoundError:
            print(f"  {dtype:<13}  no manifest")
            continue
        total    = len(manifest)
        done_ids = set(log_df[log_df["dtype"] == dtype]["file_id"])
        on_disk  = scan_existing(dtype, manifest)
        done     = len(done_ids | on_disk)
        ok       = len(log_df[(log_df["dtype"] == dtype) & (log_df["md5_ok"] == True)])
        fail     = len(log_df[(log_df["dtype"] == dtype) & (log_df["md5_ok"] == False)])
        pct      = done / total * 100 if total else 0
        bar      = "█" * int(pct // 5) + "░" * (20 - int(pct // 5))
        print(f"  {dtype:<13} {bar} {pct:5.1f}%  {ok:>8}  {fail:>8}")
    print()


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Download aligned TCGA data via gdc-client")
    parser.add_argument("--types",   nargs="+", help="Data types (default: all from config)")
    parser.add_argument("--status",  action="store_true", help="Show download progress and exit")
    parser.add_argument("--verify",  action="store_true", help="Re-verify MD5 of existing files only")
    parser.add_argument("--token",   help="Override token path from config")
    parser.add_argument("--threads", type=int, help="Override thread count from config")
    args = parser.parse_args()

    cfg    = load_config()
    dtypes = args.types or list(cfg["data_types"].keys())

    if args.token:
        cfg.setdefault("download", {})["token"] = args.token
    if args.threads:
        cfg.setdefault("download", {})["threads"] = args.threads

    if args.status:
        show_status(dtypes)
        return

    log_df = load_log()

    if args.verify:
        for dtype in dtypes:
            manifest = load_manifest(dtype)
            on_disk  = scan_existing(dtype, manifest)
            print(f"{dtype}: verifying {len(on_disk)} files...")
            records  = verify_files(dtype, manifest, on_disk)
            log_df   = append_log(log_df, records)
        save_log(log_df)
        show_status(dtypes)
        return

    for dtype in dtypes:
        print(f"\n{'='*50}\n  {dtype}\n{'='*50}")
        log_df = download_dtype(dtype, cfg, log_df)

    print("\nAll done!")
    show_status(dtypes)


if __name__ == "__main__":
    main()

"""
Align manifests across data types — keep only cases with all omics.

Supports two input modes:
  pipeline  (default)  Read per-dtype manifests from data/manifests/full/
                        Uses data/file_to_case.tsv built by 00_generate_manifests.py
  portal               Read a GDC portal manifest (--input FILE)
                        Queries GDC API to classify files and map to cases

Output:
    data/manifests/aligned/manifest_{dtype}.tsv  (gdc-client ready)
    data/aligned_cases.tsv

Usage:
    uv run python scripts/01_align_manifests.py
    uv run python scripts/01_align_manifests.py --mode portal --input data/gdc_manifest.2026-03-29.txt
"""

import argparse
import json
import time
import tomllib
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

GDC_FILES = "https://api.gdc.cancer.gov/files"
CONFIG_PATH = Path("config.toml")
FULL_DIR = Path("data/manifests/full")
ALIGNED_DIR = Path("data/manifests/aligned")
MAPPING_PATH = Path("data/file_to_case.tsv")
ALIGNED_CASES_PATH = Path("data/aligned_cases.tsv")
BATCH_SIZE = 200
MAX_RETRIES = 3


def load_config():
    with open(CONFIG_PATH, "rb") as f:
        return tomllib.load(f)


def query_file_metadata(file_ids: list[str]) -> list[dict]:
    filters = {"op": "in", "content": {"field": "file_id", "value": file_ids}}
    payload = {
        "filters": json.dumps(filters),
        "fields": (
            "file_id,file_name,file_size,md5sum,state,"
            "data_type,analysis.workflow_type,"
            "cases.submitter_id,cases.samples.sample_type,cases.project.project_id"
        ),
        "size": len(file_ids),
        "format": "json",
    }
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.post(GDC_FILES, json=payload, timeout=60)
            r.raise_for_status()
            return r.json()["data"]["hits"]
        except Exception as e:
            if attempt == MAX_RETRIES:
                raise
            tqdm.write(f"  retry {attempt}: {e}")
            time.sleep(2 ** attempt)


# ── Pipeline mode ────────────────────────────────────────────────────────────

def load_pipeline(data_types: dict) -> tuple[dict[str, pd.DataFrame], dict]:
    """Read full manifests + file_to_case.tsv."""
    mapping_df = pd.read_csv(MAPPING_PATH, sep="\t")
    mapping = {
        row["file_id"]: {
            "dtype":       row["dtype"],
            "case_id":     row["case_id"],
            "sample_type": row["sample_type"],
            "project":     row.get("project"),
        }
        for _, row in mapping_df.iterrows()
    }

    manifests = {}
    for label in data_types:
        path = FULL_DIR / f"manifest_{label}.tsv"
        if not path.exists():
            raise FileNotFoundError(f"Full manifest not found: {path}\nRun 00_generate_manifests.py first.")
        manifests[label] = pd.read_csv(path, sep="\t").rename(columns={"id": "file_id"})

    return manifests, mapping


# ── Portal mode ──────────────────────────────────────────────────────────────

def load_portal(portal_path: Path, data_types: dict) -> tuple[dict[str, pd.DataFrame], dict]:
    """Read portal manifest, query GDC for metadata, classify by dtype."""
    dtype_to_label = {dt["data_type"]: label for label, dt in data_types.items()}

    portal_df = pd.read_csv(portal_path, sep="\t")
    file_ids = portal_df["id"].tolist()
    print(f"Portal manifest: {len(file_ids)} files\n")

    # Check which dtypes are present by filename pattern (fast preview)
    missing = _preview_dtypes(portal_df, dtype_to_label)
    if missing:
        raise SystemExit(
            f"Aborted: add {missing} to your portal cart and re-download the manifest."
        )

    print(f"\nQuerying GDC metadata ({len(file_ids)} files)...")
    meta = {}
    for i in tqdm(range(0, len(file_ids), BATCH_SIZE), desc="Fetching", unit="batch"):
        batch = file_ids[i: i + BATCH_SIZE]
        for hit in query_file_metadata(batch):
            fid = hit["file_id"]
            case_id = sample_type = project = None
            if hit.get("cases"):
                c = hit["cases"][0]
                case_id = c.get("submitter_id")
                project = c.get("project", {}).get("project_id")
                samples = c.get("samples", [])
                if samples:
                    sample_type = samples[0].get("sample_type")
            meta[fid] = {
                "filename":    hit.get("file_name", ""),
                "md5":         hit.get("md5sum", ""),
                "size":        hit.get("file_size", 0),
                "state":       hit.get("state", ""),
                "data_type":   hit.get("data_type", ""),
                "case_id":     case_id,
                "sample_type": sample_type,
                "project":     project,
            }

    manifests_rows: dict[str, list] = {label: [] for label in data_types}
    mapping: dict[str, dict] = {}

    for fid in file_ids:
        m = meta.get(fid)
        if not m:
            continue
        label = dtype_to_label.get(m["data_type"])
        if not label:
            continue
        manifests_rows[label].append({
            "file_id":  fid,
            "filename": m["filename"],
            "md5":      m["md5"],
            "size":     m["size"],
            "state":    m["state"],
        })
        mapping[fid] = {
            "dtype":       label,
            "case_id":     m["case_id"],
            "sample_type": m["sample_type"],
            "project":     m["project"],
        }

    manifests = {
        label: pd.DataFrame(rows)
        for label, rows in manifests_rows.items()
    }
    return manifests, mapping


def _preview_dtypes(portal_df: pd.DataFrame, dtype_to_label: dict) -> list[str]:
    """Quick filename-based preview of what's in the portal manifest.
    Returns list of missing dtype labels."""
    patterns = {
        "aliquot_ensemble_masked": "mutation",
        "ascat3.gene_level":       "cnv_gene",
        "ascat3.allelic_specific": "cnv_seg",
        "sesame.level3betas":      "methylation",
        "star_gene_counts":        "expression",
    }
    counts = {}
    for fname in portal_df["filename"]:
        for pat, label in patterns.items():
            if pat in fname:
                counts[label] = counts.get(label, 0) + 1
                break
        else:
            counts.setdefault("other", 0)
            counts["other"] += 1

    print("Portal manifest contents (by filename pattern):")
    for label, n in sorted(counts.items()):
        mark = "" if label in dtype_to_label.values() else " ← not in config"
        print(f"  {label}: {n}{mark}")

    missing = [l for l in dtype_to_label.values() if l not in counts]
    if missing:
        print(f"\n  WARNING: missing data types: {missing}")
        print("  Add them to your portal cart before continuing.")
    return missing


# ── Alignment ────────────────────────────────────────────────────────────────

def align_and_write(manifests: dict, mapping: dict, data_types: dict):
    ALIGNED_DIR.mkdir(parents=True, exist_ok=True)

    # Cases per dtype
    dtype_cases: dict[str, set] = {}
    for label in data_types:
        df = manifests.get(label)
        if df is None or df.empty:
            raise ValueError(f"No files for dtype '{label}'. Check your manifest / config.")
        cases = {
            mapping[fid]["case_id"]
            for fid in df["file_id"]
            if fid in mapping and mapping[fid].get("case_id")
        }
        dtype_cases[label] = cases
        print(f"  {label}: {len(cases)} cases")

    aligned = set.intersection(*dtype_cases.values())
    n = len(data_types)
    print(f"\n  → aligned ({n}-way): {len(aligned)} cases\n")

    # Write per-dtype aligned manifests
    for label in data_types:
        df = manifests[label].copy()
        keep = {
            fid for fid, info in mapping.items()
            if info.get("dtype") == label and info.get("case_id") in aligned
        }
        filtered = df[df["file_id"].isin(keep)]
        out = ALIGNED_DIR / f"manifest_{label}.tsv"
        with open(out, "w") as f:
            f.write("id\tfilename\tmd5\tsize\tstate\n")
            for _, row in filtered.iterrows():
                f.write(f"{row['file_id']}\t{row['filename']}\t{row.get('md5','')}\t{row.get('size',0)}\t{row.get('state','')}\n")
        print(f"  {label}: {len(filtered)} files → {out}")

    # aligned_cases.tsv
    project_map = {
        info["case_id"]: info["project"]
        for info in mapping.values()
        if info.get("case_id") in aligned and info.get("project")
    }
    rows = [{"submitter_id": c, "project": project_map.get(c, "unknown")} for c in sorted(aligned)]
    aligned_df = pd.DataFrame(rows)
    aligned_df.to_csv(ALIGNED_CASES_PATH, sep="\t", index=False)
    print(f"\nSaved: {ALIGNED_CASES_PATH}")
    print(aligned_df["project"].value_counts().to_string())


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Align manifests across data types")
    parser.add_argument("--mode", choices=["pipeline", "portal"], default="pipeline",
                        help="Input source (default: pipeline)")
    parser.add_argument("--input", type=str,
                        help="Path to portal manifest file (required for --mode portal)")
    args = parser.parse_args()

    cfg = load_config()
    data_types = cfg["data_types"]

    print(f"=== Aligning manifests (mode: {args.mode}) ===\n")

    if args.mode == "pipeline":
        manifests, mapping = load_pipeline(data_types)
    else:
        if not args.input:
            parser.error("--input is required for portal mode")
        manifests, mapping = load_portal(Path(args.input), data_types)

    print("Cases per data type:")
    align_and_write(manifests, mapping, data_types)


if __name__ == "__main__":
    main()

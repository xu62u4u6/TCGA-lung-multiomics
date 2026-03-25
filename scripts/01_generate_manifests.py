"""
Generate GDC manifests and file-to-case mapping for aligned multi-omics cases.

1. Query GDC API for all files matching aligned cases
2. Write per-datatype manifests (for gdc-client download)
3. Build file_id → case_id + sample_type mapping table

Usage:
    uv run python scripts/01_generate_manifests.py
"""

import json
import time
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

GDC_FILES = "https://api.gdc.cancer.gov/files"

DATA_TYPES = {
    "expression": {
        "data_type": "Gene Expression Quantification",
        "workflow_type": "STAR - Counts",
    },
    "mutation": {
        "data_type": "Masked Somatic Mutation",
    },
    "cnv": {
        "data_type": "Masked Copy Number Segment",
    },
    "cnv_gene": {
        "data_type": "Gene Level Copy Number",
    },
    "methylation": {
        "data_type": "Methylation Beta Value",
    },
    "clinical": {
        "data_type": "Clinical Supplement",
    },
}

PROJECTS = ["TCGA-LUAD", "TCGA-LUSC"]
ALIGNED_CASES = Path("data/aligned_cases.tsv")
MANIFEST_DIR = Path("data/manifests")
MAPPING_PATH = Path("data/file_to_case.tsv")

BATCH_SIZE = 200
MAX_RETRIES = 3


def load_aligned_cases() -> dict[str, list[str]]:
    df = pd.read_csv(ALIGNED_CASES, sep="\t")
    return {proj: group["submitter_id"].tolist() for proj, group in df.groupby("project")}


def get_files(project_id: str, data_type: str, cases: list[str],
              workflow_type: str = None) -> list[dict]:
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": project_id}},
            {"op": "=", "content": {"field": "data_type", "value": data_type}},
            {"op": "in", "content": {"field": "cases.submitter_id", "value": cases}},
        ],
    }
    if workflow_type:
        filters["content"].append(
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": workflow_type}}
        )

    payload = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,file_size,md5sum,state",
        "size": 5000,
        "format": "json",
    }
    r = requests.post(GDC_FILES, json=payload)
    r.raise_for_status()
    return r.json()["data"]["hits"]


def write_manifest(hits: list[dict], output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("id\tfilename\tmd5\tsize\tstate\n")
        for h in hits:
            f.write(f"{h['file_id']}\t{h['file_name']}\t{h.get('md5sum','')}\t{h['file_size']}\t{h.get('state','')}\n")


def query_case_mapping(file_ids: list[str]) -> list[dict]:
    """Query GDC API for file_id → case_id + sample_type."""
    filters = {
        "op": "in",
        "content": {"field": "file_id", "value": file_ids},
    }
    payload = {
        "filters": json.dumps(filters),
        "fields": "file_id,cases.submitter_id,cases.samples.sample_type",
        "size": len(file_ids),
        "format": "json",
    }
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.post(GDC_FILES, json=payload)
            r.raise_for_status()
            return r.json()["data"]["hits"]
        except Exception as e:
            if attempt == MAX_RETRIES:
                raise
            tqdm.write(f"  retry {attempt}: {e}")
            time.sleep(2 ** attempt)


def build_file_mapping(all_files: pd.DataFrame) -> pd.DataFrame:
    """Build file_id → case_id + sample_type mapping via GDC API."""
    unique_ids = all_files["file_id"].unique().tolist()
    print(f"\nBuilding file mapping ({len(unique_ids)} unique files)...")

    mapping = {}
    pbar = tqdm(total=len(unique_ids), desc="Querying GDC", unit="file")

    for i in range(0, len(unique_ids), BATCH_SIZE):
        batch = unique_ids[i: i + BATCH_SIZE]
        hits = query_case_mapping(batch)
        for hit in hits:
            fid = hit["file_id"]
            case_id = None
            sample_type = None
            if hit.get("cases"):
                case_id = hit["cases"][0].get("submitter_id")
                samples = hit["cases"][0].get("samples", [])
                if samples:
                    sample_type = samples[0].get("sample_type")
            mapping[fid] = (case_id, sample_type)
        pbar.update(len(batch))

    pbar.close()

    all_files["case_id"] = all_files["file_id"].map(lambda x: mapping.get(x, (None, None))[0])
    all_files["sample_type"] = all_files["file_id"].map(lambda x: mapping.get(x, (None, None))[1])
    return all_files


def main():
    aligned = load_aligned_cases()

    # Stage 1: Generate manifests
    print("=== Generating manifests ===\n")
    all_file_records = []
    all_hits = []

    for label, dt in DATA_TYPES.items():
        label_hits = []
        for proj in PROJECTS:
            hits = get_files(proj, dt["data_type"], aligned[proj], dt.get("workflow_type"))
            label_hits.extend(hits)
            print(f"  {proj} / {label}: {len(hits)} files")

        write_manifest(label_hits, MANIFEST_DIR / f"manifest_{label}.tsv")
        all_hits.extend(label_hits)

        # Collect for mapping
        for h in label_hits:
            all_file_records.append({
                "file_id": h["file_id"],
                "filename": h["file_name"],
                "data_type": label,
            })

        total_size = sum(h["file_size"] for h in label_hits)
        print(f"  → manifest_{label}.tsv ({len(label_hits)} files, {total_size / 1e9:.1f} GB)\n")

    write_manifest(all_hits, MANIFEST_DIR / "manifest_all.tsv")
    total = sum(h["file_size"] for h in all_hits)
    print(f"Total: {len(all_hits)} files, {total / 1e9:.1f} GB")

    # Stage 2: Build file-to-case mapping
    print("\n=== Building file-to-case mapping ===")
    all_files_df = pd.DataFrame(all_file_records)
    mapping_df = build_file_mapping(all_files_df)
    mapping_df.to_csv(MAPPING_PATH, sep="\t", index=False)

    print(f"\nSaved: {MAPPING_PATH}")
    print(f"  {len(mapping_df)} rows, {mapping_df['case_id'].notna().sum()} with case_id")
    print(f"  sample_types: {mapping_df['sample_type'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()

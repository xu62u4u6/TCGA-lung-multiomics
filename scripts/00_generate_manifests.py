"""
Generate full GDC manifests for all TCGA-LUAD/LUSC files matching config data types.

No case-level filtering — produces complete per-dtype manifests and a
file_id → case_id / sample_type / project mapping table.

Output:
    data/manifests/full/manifest_{dtype}.tsv
    data/file_to_case.tsv

Usage:
    uv run python scripts/00_generate_manifests.py
"""

import json
import time
import tomllib
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

GDC_FILES = "https://api.gdc.cancer.gov/files"
CONFIG_PATH = Path("config.toml")
MANIFEST_DIR = Path("data/manifests/full")
MAPPING_PATH = Path("data/file_to_case.tsv")
BATCH_SIZE = 200
MAX_RETRIES = 3


def load_config():
    with open(CONFIG_PATH, "rb") as f:
        return tomllib.load(f)


def get_files(project_id: str, data_type: str, workflow_type: str = None) -> list[dict]:
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": project_id}},
            {"op": "=", "content": {"field": "data_type", "value": data_type}},
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
    r = requests.post(GDC_FILES, json=payload, timeout=60)
    r.raise_for_status()
    return r.json()["data"]["hits"]


def write_manifest(hits: list[dict], path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("id\tfilename\tmd5\tsize\tstate\n")
        for h in hits:
            f.write(f"{h['file_id']}\t{h['file_name']}\t{h.get('md5sum','')}\t{h['file_size']}\t{h.get('state','')}\n")


def query_case_mapping(file_ids: list[str]) -> list[dict]:
    filters = {"op": "in", "content": {"field": "file_id", "value": file_ids}}
    payload = {
        "filters": json.dumps(filters),
        "fields": "file_id,cases.submitter_id,cases.samples.sample_type,cases.project.project_id",
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


def build_file_mapping(records: list[dict]) -> pd.DataFrame:
    unique_ids = list({r["file_id"] for r in records})
    print(f"\nBuilding file mapping ({len(unique_ids)} files)...")

    mapping = {}
    for i in tqdm(range(0, len(unique_ids), BATCH_SIZE), desc="Querying GDC", unit="batch"):
        batch = unique_ids[i: i + BATCH_SIZE]
        for hit in query_case_mapping(batch):
            fid = hit["file_id"]
            case_id = sample_type = project = None
            if hit.get("cases"):
                c = hit["cases"][0]
                case_id = c.get("submitter_id")
                project = c.get("project", {}).get("project_id")
                samples = c.get("samples", [])
                if samples:
                    sample_type = samples[0].get("sample_type")
            mapping[fid] = (case_id, sample_type, project)

    df = pd.DataFrame(records)
    df["case_id"]     = df["file_id"].map(lambda x: mapping.get(x, (None, None, None))[0])
    df["sample_type"] = df["file_id"].map(lambda x: mapping.get(x, (None, None, None))[1])
    df["project"]     = df["file_id"].map(lambda x: mapping.get(x, (None, None, None))[2])
    return df


def main():
    cfg = load_config()
    projects = cfg["projects"]
    data_types = cfg["data_types"]

    print("=== Generating full manifests ===\n")
    records = []

    for label, dt in data_types.items():
        label_hits = []
        for proj in projects:
            hits = get_files(proj, dt["data_type"], dt.get("workflow_type"))
            label_hits.extend(hits)
            print(f"  {proj} / {label}: {len(hits)} files")

        write_manifest(label_hits, MANIFEST_DIR / f"manifest_{label}.tsv")
        for h in label_hits:
            records.append({"file_id": h["file_id"], "filename": h["file_name"], "dtype": label})

        total_gb = sum(h["file_size"] for h in label_hits) / 1e9
        print(f"  → manifest_{label}.tsv ({len(label_hits)} files, {total_gb:.1f} GB)\n")

    mapping_df = build_file_mapping(records)
    mapping_df.to_csv(MAPPING_PATH, sep="\t", index=False)
    print(f"Saved: {MAPPING_PATH} ({len(mapping_df)} rows)")


if __name__ == "__main__":
    main()

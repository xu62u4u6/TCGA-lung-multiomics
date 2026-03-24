"""Generate GDC manifests for aligned multi-omics cases."""

import json
import os
from pathlib import Path

import requests

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
    "methylation": {
        "data_type": "Methylation Beta Value",
    },
    "clinical": {
        "data_type": "Clinical Supplement",
    },
}

PROJECTS = ["TCGA-LUAD", "TCGA-LUSC"]


def load_aligned_cases(path: str) -> dict[str, set]:
    result = {}
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            proj, sid = line.strip().split("\t")
            result.setdefault(proj, set()).add(sid)
    return result


def get_files(project_id: str, data_type: str, aligned_cases: set,
              workflow_type: str = None) -> list[dict]:
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": project_id}},
            {"op": "=", "content": {"field": "data_type", "value": data_type}},
            {"op": "in", "content": {"field": "cases.submitter_id", "value": list(aligned_cases)}},
        ],
    }
    if workflow_type:
        filters["content"].append(
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": workflow_type}}
        )

    payload = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,file_size,md5sum,state,cases.submitter_id",
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


def main():
    aligned = load_aligned_cases("data/aligned_cases.tsv")
    manifest_dir = Path("data/manifests")

    all_hits = []
    for label, dt in DATA_TYPES.items():
        label_hits = []
        for proj in PROJECTS:
            hits = get_files(proj, dt["data_type"], aligned[proj], dt.get("workflow_type"))
            label_hits.extend(hits)
            print(f"  {proj} / {label}: {len(hits)} files")

        # Per data-type manifest
        write_manifest(label_hits, manifest_dir / f"manifest_{label}.tsv")
        all_hits.extend(label_hits)

        total_size = sum(h["file_size"] for h in label_hits)
        print(f"  → manifest_{label}.tsv ({len(label_hits)} files, {total_size / 1e9:.1f} GB)\n")

    # Combined manifest
    write_manifest(all_hits, manifest_dir / "manifest_all.tsv")
    total = sum(h["file_size"] for h in all_hits)
    print(f"Total: {len(all_hits)} files, {total / 1e9:.1f} GB")
    print(f"Manifests saved to {manifest_dir}/")


if __name__ == "__main__":
    main()

"""
Download aligned multi-omics data for TCGA-LUAD and TCGA-LUSC.

5 data types: Expression, Mutation (MAF), CNV, Methylation, Clinical
Aligns samples at case level (submitter_id) then downloads only aligned cases.

Usage:
    uv run python scripts/download_aligned.py --token /path/to/gdc-token.txt
    uv run python scripts/download_aligned.py --token /path/to/gdc-token.txt --align-only  # just show counts
"""

import argparse
import gzip
import io
import json
import os
import tarfile
from pathlib import Path

import requests

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data"

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


def get_cases_by_datatype(project_id: str, data_type: str, workflow_type: str = None) -> set:
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

    params = {
        "filters": json.dumps(filters),
        "fields": "cases.submitter_id",
        "size": 5000,
        "format": "json",
    }
    r = requests.get(GDC_FILES, params=params)
    r.raise_for_status()
    cases = set()
    for hit in r.json()["data"]["hits"]:
        for c in hit.get("cases", []):
            cases.add(c["submitter_id"])
    return cases


def align_cases(projects: list[str], data_types: dict) -> dict[str, set]:
    """Return {project: aligned_case_ids}"""
    result = {}
    for proj in projects:
        print(f"\n{'='*50}")
        print(f"  {proj}")
        print(f"{'='*50}")
        case_sets = {}
        for label, dt in data_types.items():
            cases = get_cases_by_datatype(proj, dt["data_type"], dt.get("workflow_type"))
            case_sets[label] = cases
            print(f"  {label}: {len(cases)} cases")

        aligned = set.intersection(*case_sets.values())
        print(f"  >>> Aligned ({len(case_sets)}-way): {len(aligned)} cases")
        result[proj] = aligned
    return result


def get_file_ids(project_id: str, data_type: str, aligned_cases: set,
                 workflow_type: str = None) -> list[dict]:
    """Get file_id + metadata for aligned cases."""
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
        "fields": "file_id,file_name,file_size,cases.submitter_id,cases.samples.sample_type",
        "size": 5000,
        "format": "json",
    }
    r = requests.post(GDC_FILES, json=payload)
    r.raise_for_status()
    return r.json()["data"]["hits"]


def download_files(file_ids: list[str], output_dir: Path, token: str = None,
                   batch_size: int = 50):
    """Download files from GDC in batches."""
    output_dir.mkdir(parents=True, exist_ok=True)
    headers = {"Content-Type": "application/json"}
    if token:
        headers["X-Auth-Token"] = token

    for i in range(0, len(file_ids), batch_size):
        batch = file_ids[i : i + batch_size]
        print(f"    Downloading batch {i // batch_size + 1} ({len(batch)} files)...")

        payload = {"ids": batch}
        r = requests.post(GDC_DATA, json=payload, headers=headers, stream=True)
        r.raise_for_status()

        content_type = r.headers.get("Content-Type", "")

        if "application/x-tar" in content_type or len(batch) > 1:
            # Multi-file: tar.gz response
            buf = io.BytesIO(r.content)
            try:
                with tarfile.open(fileobj=buf) as tar:
                    for member in tar.getmembers():
                        if member.isfile():
                            fname = os.path.basename(member.name)
                            if fname == "MANIFEST.txt":
                                continue
                            f = tar.extractfile(member)
                            if f:
                                (output_dir / fname).write_bytes(f.read())
            except tarfile.ReadError:
                # Single file returned despite batch
                (output_dir / "batch_{i}.dat").write_bytes(buf.getvalue())
        else:
            # Single file
            cd = r.headers.get("Content-Disposition", "")
            fname = cd.split("filename=")[-1].strip('"') if "filename=" in cd else f"{batch[0]}.dat"
            (output_dir / fname).write_bytes(r.content)

    print(f"    Done: {len(file_ids)} files → {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="Download aligned TCGA lung multi-omics data")
    parser.add_argument("--token", type=str, help="Path to GDC token file")
    parser.add_argument("--outdir", type=str, default="data/raw", help="Output base directory")
    parser.add_argument("--align-only", action="store_true", help="Only show alignment counts, no download")
    args = parser.parse_args()

    # Align
    aligned = align_cases(PROJECTS, DATA_TYPES)
    total = sum(len(v) for v in aligned.values())
    print(f"\nTotal aligned cases: {total}")

    # Save aligned case IDs
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    aligned_file = outdir.parent / "aligned_cases.tsv"
    with open(aligned_file, "w") as f:
        f.write("project\tsubmitter_id\n")
        for proj, cases in aligned.items():
            for c in sorted(cases):
                f.write(f"{proj}\t{c}\n")
    print(f"Saved aligned case IDs → {aligned_file}")

    if args.align_only:
        return

    # Read token
    token = None
    if args.token:
        token = Path(args.token).read_text().strip()

    # Download each data type
    for label, dt in DATA_TYPES.items():
        print(f"\n--- {label} ---")
        for proj in PROJECTS:
            proj_short = proj.split("-")[1].lower()
            out = outdir / label / proj_short
            hits = get_file_ids(proj, dt["data_type"], aligned[proj], dt.get("workflow_type"))
            fids = [h["file_id"] for h in hits]
            print(f"  {proj}: {len(fids)} files")
            download_files(fids, out, token=token)


if __name__ == "__main__":
    main()

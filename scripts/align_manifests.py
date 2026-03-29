"""
Align manifests across data types — keep only cases with all omics.

Usage:
    uv run python scripts/align_manifests.py
    uv run python scripts/align_manifests.py --mode portal --input data/gdc_manifest.txt
"""

import argparse
from pathlib import Path
from pymaftools.io.tcga import GDCClient

parser = argparse.ArgumentParser(description="Align manifests across data types")
parser.add_argument("--mode",  choices=["pipeline", "portal"], default="pipeline")
parser.add_argument("--input", help="Portal manifest file (--mode portal)")
args = parser.parse_args()

client = GDCClient.from_config("config.toml")
client.align_manifests(
    mode               = args.mode,
    full_manifest_dir  = Path("data/manifests/full"),
    portal_path        = Path(args.input) if args.input else None,
    mapping_path       = Path("data/file_to_case.tsv"),
    outdir             = Path("data/manifests/aligned"),
    aligned_cases_path = Path("data/aligned_cases.tsv"),
)

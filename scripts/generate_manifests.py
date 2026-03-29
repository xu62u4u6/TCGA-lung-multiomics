"""
Generate full GDC manifests for all projects and data types in config.toml.

Usage:
    uv run python scripts/generate_manifests.py
"""

from pathlib import Path
from pymaftools.io.tcga import GDCClient

client = GDCClient.from_config("config.toml")
client.generate_full_manifests(
    outdir       = Path("data/manifests/full"),
    mapping_path = Path("data/file_to_case.tsv"),
)

"""
Build raw tables from downloaded TCGA data.

Reads downloaded files via pymaftools builders and saves each table as HDF5.
No filtering — tables contain all sample types.

Output:
    data/processed/tables/{dtype}.h5
    data/processed/tables/cnv_seg.parquet   (segment-level, from cnv_seg)
    data/processed/tables/cnv_cytoband.h5   (arm-level, derived from cnv_seg)
    data/processed/tables/mutation.maf      (raw MAF)

Usage:
    uv run python scripts/build_tables.py
    uv run python scripts/build_tables.py --types expression mutation
"""

from __future__ import annotations

import argparse
import tomllib
from pathlib import Path

import pandas as pd

from pymaftools.io.tcga import (
    TCGACNVGeneBuilder,
    TCGACNVSegmentBuilder,
    TCGAExpressionBuilder,
    TCGAMethylationBuilder,
    TCGAMutationBuilder,
    load_file_mapping,
)

CONFIG_PATH = Path("config.toml")
RAW_DIR     = Path("data/raw")
MAPPING     = Path("data/file_to_case.tsv")
TABLE_DIR   = Path("data/processed/tables")


def load_config():
    with open(CONFIG_PATH, "rb") as f:
        return tomllib.load(f)


def save_pivot(table, name: str):
    out = TABLE_DIR / f"{name}.h5"
    with pd.HDFStore(str(out), mode="w") as store:
        store.put("data", pd.DataFrame(table))
        store.put("feature_metadata", table.feature_metadata)
        store.put("sample_metadata", table.sample_metadata)
    return out


def main():
    cfg    = load_config()
    dtypes = list(cfg["data_types"].keys())

    parser = argparse.ArgumentParser(description="Build raw tables from TCGA downloads")
    parser.add_argument("--types", nargs="+", default=dtypes)
    args = parser.parse_args()

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    mapping = load_file_mapping(MAPPING)
    print(f"Mapping: {len(mapping)} files\n")

    builders = {
        "expression": lambda: TCGAExpressionBuilder(RAW_DIR / "expression", mapping).build(),
        "mutation":   lambda: TCGAMutationBuilder(RAW_DIR / "mutation", mapping).build(),
        "cnv_seg":    lambda: TCGACNVSegmentBuilder(RAW_DIR / "cnv_seg", mapping),
        "cnv_gene":   lambda: TCGACNVGeneBuilder(RAW_DIR / "cnv_gene", mapping).build(),
        "methylation":lambda: TCGAMethylationBuilder(RAW_DIR / "methylation", mapping).build(),
    }

    for dtype in args.types:
        if dtype not in builders:
            print(f"No builder for '{dtype}', skipping\n")
            continue

        print(f"=== {dtype} ===")
        try:
            result = builders[dtype]()

            if dtype == "cnv_seg":
                seg_df = result.build()
                seg_df.to_parquet(TABLE_DIR / "cnv_seg.parquet")
                print(f"  saved segments: {TABLE_DIR / 'cnv_seg.parquet'}")

                cyto = result.build_cytoband_table(seg_df)
                save_pivot(cyto, "cnv_cytoband")
                print(f"  saved cytoband: {TABLE_DIR / 'cnv_cytoband.h5'}")

            elif dtype == "mutation":
                result.to_MAF(TABLE_DIR / "mutation.maf")
                print(f"  saved MAF: {TABLE_DIR / 'mutation.maf'}")
                pivot = result.to_pivot_table()
                out = save_pivot(pivot, "mutation")
                print(f"  saved pivot: {out} ({pivot.shape})")

            else:
                out = save_pivot(result, dtype)
                print(f"  saved: {out} ({result.shape})")

        except Exception as e:
            print(f"  ERROR: {e}")
        print()

    print("Tables:")
    for f in sorted(TABLE_DIR.iterdir()):
        print(f"  {f.name} ({f.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()

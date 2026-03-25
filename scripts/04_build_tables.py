"""
Build raw tables from downloaded TCGA data.

Reads all downloaded files via TCGA builders and saves each table
independently as HDF5. No filtering — tables contain all samples
(all sample types).

Output: data/processed/tables/{type}.h5

Usage:
    uv run python scripts/04_build_tables.py
    uv run python scripts/04_build_tables.py --types expression mutation
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pymaftools.io.tcga import (
    TCGAClinicalBuilder,
    TCGACNVGeneBuilder,
    TCGACNVSegmentBuilder,
    TCGAExpressionBuilder,
    TCGAMethylationBuilder,
    TCGAMutationBuilder,
    load_file_mapping,
)

RAW_DIR = Path("data/raw")
MAPPING = Path("data/file_to_case.tsv")
TABLE_DIR = Path("data/processed/tables")

ALL_TYPES = ["expression", "mutation", "cnv_seg", "cnv_gene", "methylation", "clinical"]


def main():
    parser = argparse.ArgumentParser(description="Build raw tables from TCGA downloads")
    parser.add_argument("--types", nargs="+", default=ALL_TYPES)
    args = parser.parse_args()

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    mapping = load_file_mapping(MAPPING)
    print(f"Mapping: {len(mapping)} files\n")

    builders = {
        "expression": lambda: TCGAExpressionBuilder(RAW_DIR / "expression", mapping).build(),
        "mutation": lambda: TCGAMutationBuilder(RAW_DIR / "mutation", mapping).build(),
        "cnv_seg": lambda: TCGACNVSegmentBuilder(RAW_DIR / "cnv", mapping),
        "cnv_gene": lambda: TCGACNVGeneBuilder(RAW_DIR / "cnv_gene", mapping).build(),
        "methylation": lambda: TCGAMethylationBuilder(RAW_DIR / "methylation", mapping).build(),
        "clinical": lambda: TCGAClinicalBuilder(RAW_DIR / "clinical", mapping).build(),
    }

    for dtype in args.types:
        if dtype not in builders:
            print(f"Unknown type: {dtype}, skipping\n")
            continue

        print(f"=== {dtype} ===")
        try:
            result = builders[dtype]()

            if dtype == "cnv_seg":
                # CNV segment builder: save both raw segments and cytoband table
                seg_builder = result
                seg_df = seg_builder.build()
                seg_df.to_parquet(TABLE_DIR / "cnv_seg.parquet")
                print(f"  saved segments: {TABLE_DIR / 'cnv_seg.parquet'}")

                cyto = seg_builder.build_cytoband_table(seg_df)
                store = pd.HDFStore(str(TABLE_DIR / "cnv_cytoband.h5"), mode="w")
                store.put("data", pd.DataFrame(cyto))
                store.put("feature_metadata", cyto.feature_metadata)
                store.put("sample_metadata", cyto.sample_metadata)
                store.close()
                print(f"  saved cytoband: {TABLE_DIR / 'cnv_cytoband.h5'}")

            elif dtype == "clinical":
                # Clinical is a plain DataFrame
                result.to_parquet(TABLE_DIR / "clinical.parquet")
                print(f"  saved: {TABLE_DIR / 'clinical.parquet'} ({result.shape})")

            elif dtype == "mutation":
                # MAF: save as .maf + pivot table as h5
                result.to_MAF(TABLE_DIR / "mutation.maf")
                print(f"  saved MAF: {TABLE_DIR / 'mutation.maf'}")

                pivot = result.to_pivot_table()
                store = pd.HDFStore(str(TABLE_DIR / "mutation.h5"), mode="w")
                store.put("data", pd.DataFrame(pivot))
                store.put("feature_metadata", pivot.feature_metadata)
                store.put("sample_metadata", pivot.sample_metadata)
                store.close()
                print(f"  saved pivot: {TABLE_DIR / 'mutation.h5'} ({pivot.shape})")

            else:
                # PivotTable types: save as HDF5 with metadata
                out = TABLE_DIR / f"{dtype}.h5"
                store = pd.HDFStore(str(out), mode="w")
                store.put("data", pd.DataFrame(result))
                store.put("feature_metadata", result.feature_metadata)
                store.put("sample_metadata", result.sample_metadata)
                store.close()
                print(f"  saved: {out} ({result.shape})")

        except Exception as e:
            print(f"  ERROR: {e}")

        print()

    print("Done! Tables:")
    for f in sorted(TABLE_DIR.iterdir()):
        size_mb = f.stat().st_size / 1e6
        print(f"  {f.name} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()

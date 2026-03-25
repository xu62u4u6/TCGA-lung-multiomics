"""
Assemble filtered tables into a Cohort and save as HDF5.

Reads filtered tables from data/processed/filtered/,
builds a Cohort object, and saves to data/processed/tcga_lung.h5.

Usage:
    uv run python scripts/06_build_cohort.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pymaftools.core.Cohort import Cohort
from pymaftools.core.CopyNumberVariationTable import CopyNumberVariationTable
from pymaftools.core.ExpressionTable import ExpressionTable
from pymaftools.core.PivotTable import PivotTable

FILTERED_DIR = Path("data/processed/filtered")
COHORT_PATH = Path("data/processed/tcga_lung.h5")

TABLE_TYPES = ["expression", "mutation", "cnv_cytoband", "cnv_gene", "methylation"]


def load_pivot_table(name: str, cls=PivotTable):
    path = FILTERED_DIR / f"{name}.h5"
    with pd.HDFStore(str(path), mode="r") as store:
        data = store.get("data")
        feature_meta = store.get("feature_metadata")
        sample_meta = store.get("sample_metadata")

    table = cls(data)
    table.feature_metadata = feature_meta
    table.sample_metadata = sample_meta
    return table


def main():
    parser = argparse.ArgumentParser(description="Assemble Cohort from filtered tables")
    parser.add_argument("--types", nargs="+", default=TABLE_TYPES)
    parser.add_argument("--output", type=str, default=str(COHORT_PATH))
    args = parser.parse_args()

    cohort = Cohort(name="TCGA-Lung", description="TCGA LUAD+LUSC multi-omics cohort")

    # Clinical sets sample_IDs
    print("=== Clinical ===")
    clin = pd.read_parquet(FILTERED_DIR / "clinical.parquet")
    cohort.add_sample_metadata(clin)
    print(f"  sample_IDs: {len(cohort.sample_IDs)}")

    # Add tables
    type_to_cls = {
        "expression": ExpressionTable,
        "mutation": PivotTable,
        "cnv_cytoband": CopyNumberVariationTable,
        "cnv_gene": CopyNumberVariationTable,
        "methylation": PivotTable,
    }

    for dtype in args.types:
        print(f"\n=== {dtype} ===")
        try:
            cls = type_to_cls.get(dtype, PivotTable)
            table = load_pivot_table(dtype, cls)
            cohort.add_table(table, dtype)
            print(f"  {cohort.tables[dtype].shape}")
        except Exception as e:
            print(f"  ERROR: {e}")

    # Save
    cohort.to_hdf5(args.output)
    print(f"\nCohort saved: {args.output}")
    for name, t in cohort.tables.items():
        fm = list(t.feature_metadata.columns) if hasattr(t, "feature_metadata") else []
        print(f"  {name}: {t.shape}, feature_metadata: {fm}")


if __name__ == "__main__":
    main()

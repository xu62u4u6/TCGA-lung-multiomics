"""
Filter raw tables by sample_type and aligned cases.

Reads tables from data/processed/tables/, applies filters,
and saves filtered tables to data/processed/filtered/.

Usage:
    uv run python scripts/05_filter_tables.py
    uv run python scripts/05_filter_tables.py --sample-type "Primary Tumor"
    uv run python scripts/05_filter_tables.py --types expression mutation
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pymaftools.core.CopyNumberVariationTable import CopyNumberVariationTable
from pymaftools.core.ExpressionTable import ExpressionTable
from pymaftools.core.MAF import MAF
from pymaftools.core.PivotTable import PivotTable

TABLE_DIR = Path("data/processed/tables")
FILTERED_DIR = Path("data/processed/filtered")
ALIGNED_CASES = Path("data/aligned_cases.tsv")

ALL_TYPES = ["expression", "mutation", "cnv_cytoband", "cnv_gene", "methylation", "clinical"]


def load_aligned_ids() -> set[str]:
    df = pd.read_csv(ALIGNED_CASES, sep="\t")
    return set(df["submitter_id"])


def load_pivot_table(name: str, cls=PivotTable):
    """Load a PivotTable from HDF5."""
    path = TABLE_DIR / f"{name}.h5"
    with pd.HDFStore(str(path), mode="r") as store:
        data = store.get("data")
        feature_meta = store.get("feature_metadata")
        sample_meta = store.get("sample_metadata")

    table = cls(data)
    table.feature_metadata = feature_meta
    table.sample_metadata = sample_meta
    return table


def save_pivot_table(table, name: str, out_dir: Path):
    """Save a PivotTable to HDF5."""
    out = out_dir / f"{name}.h5"
    with pd.HDFStore(str(out), mode="w") as store:
        store.put("data", pd.DataFrame(table))
        store.put("feature_metadata", table.feature_metadata)
        store.put("sample_metadata", table.sample_metadata)
    return out


def filter_by_sample_type_and_aligned(table, aligned: set[str], sample_type: str):
    """Filter PivotTable columns: sample_type match AND in aligned set."""
    mask = table.sample_metadata["sample_type"] == sample_type
    keep = [c for c in mask[mask].index if c in aligned]
    return table.subset(samples=keep), len(table.columns), len(keep)


def main():
    parser = argparse.ArgumentParser(description="Filter raw tables")
    parser.add_argument("--types", nargs="+", default=ALL_TYPES)
    parser.add_argument("--sample-type", type=str, default="Primary Tumor")
    args = parser.parse_args()

    FILTERED_DIR.mkdir(parents=True, exist_ok=True)
    aligned = load_aligned_ids()
    print(f"Aligned cases: {len(aligned)}")
    print(f"Sample type filter: {args.sample_type!r}\n")

    for dtype in args.types:
        print(f"=== {dtype} ===")
        try:
            if dtype == "expression":
                table = load_pivot_table("expression", ExpressionTable)
                filtered, before, after = filter_by_sample_type_and_aligned(
                    table, aligned, args.sample_type
                )
                out = save_pivot_table(filtered, "expression", FILTERED_DIR)
                print(f"  {before} → {after} samples, saved: {out}")

            elif dtype == "mutation":
                # Filter MAF rows, then save both MAF and pivot
                maf_df = pd.read_csv(TABLE_DIR / "mutation.maf", sep="\t", comment="#", low_memory=False)
                maf_df = maf_df[maf_df["sample_type"] == args.sample_type]
                maf_df = maf_df[maf_df["sample_ID"].isin(aligned)]
                before_cases = pd.read_csv(TABLE_DIR / "mutation.maf", sep="\t", comment="#", usecols=["sample_ID"])["sample_ID"].nunique()
                after_cases = maf_df["sample_ID"].nunique()

                maf = MAF(maf_df)
                maf.to_MAF(FILTERED_DIR / "mutation.maf")

                pivot = maf.to_pivot_table()
                save_pivot_table(pivot, "mutation", FILTERED_DIR)
                print(f"  {before_cases} → {after_cases} cases, saved MAF + pivot")

            elif dtype == "cnv_cytoband":
                table = load_pivot_table("cnv_cytoband", CopyNumberVariationTable)
                filtered, before, after = filter_by_sample_type_and_aligned(
                    table, aligned, args.sample_type
                )
                out = save_pivot_table(filtered, "cnv_cytoband", FILTERED_DIR)
                print(f"  {before} → {after} samples, saved: {out}")

            elif dtype == "cnv_gene":
                table = load_pivot_table("cnv_gene", CopyNumberVariationTable)
                filtered, before, after = filter_by_sample_type_and_aligned(
                    table, aligned, args.sample_type
                )
                out = save_pivot_table(filtered, "cnv_gene", FILTERED_DIR)
                print(f"  {before} → {after} samples, saved: {out}")

            elif dtype == "methylation":
                table = load_pivot_table("methylation")
                filtered, before, after = filter_by_sample_type_and_aligned(
                    table, aligned, args.sample_type
                )
                out = save_pivot_table(filtered, "methylation", FILTERED_DIR)
                print(f"  {before} → {after} samples, saved: {out}")

            elif dtype == "clinical":
                clin = pd.read_parquet(TABLE_DIR / "clinical.parquet")
                keep = [c for c in clin.index if c in aligned]
                filtered = clin.loc[keep]
                filtered.to_parquet(FILTERED_DIR / "clinical.parquet")
                print(f"  {len(clin)} → {len(filtered)} patients")

        except Exception as e:
            print(f"  ERROR: {e}")

        print()

    print("Done! Filtered tables:")
    for f in sorted(FILTERED_DIR.iterdir()):
        size_mb = f.stat().st_size / 1e6
        print(f"  {f.name} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()

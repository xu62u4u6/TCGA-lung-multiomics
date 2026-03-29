"""
Filter raw tables to aligned Primary Tumor samples.

Reads tables from data/processed/tables/, keeps only samples that are
in aligned_cases.tsv and match the target sample type.

Output:
    data/processed/filtered/{dtype}.h5
    data/processed/filtered/mutation.maf

Usage:
    uv run python scripts/filter_tables.py
    uv run python scripts/filter_tables.py --sample-type "Primary Tumor"
    uv run python scripts/filter_tables.py --types expression mutation
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pymaftools.core.CopyNumberVariationTable import CopyNumberVariationTable
from pymaftools.core.ExpressionTable import ExpressionTable
from pymaftools.core.MAF import MAF
from pymaftools.core.PivotTable import PivotTable

TABLE_DIR    = Path("data/processed/tables")
FILTERED_DIR = Path("data/processed/filtered")
ALIGNED_CASES = Path("data/aligned_cases.tsv")

# cnv_cytoband is derived from cnv_seg during build_tables — filter it too
ALL_TYPES = ["expression", "mutation", "cnv_cytoband", "cnv_gene", "methylation"]


def load_aligned_ids() -> set[str]:
    return set(pd.read_csv(ALIGNED_CASES, sep="\t")["submitter_id"])


def load_pivot(name: str, cls=PivotTable):
    path = TABLE_DIR / f"{name}.h5"
    with pd.HDFStore(str(path), mode="r") as store:
        data        = store.get("data")
        feature_meta = store.get("feature_metadata")
        sample_meta  = store.get("sample_metadata")
    table = cls(data)
    table.feature_metadata = feature_meta
    table.sample_metadata  = sample_meta
    return table


def save_pivot(table, name: str):
    out = FILTERED_DIR / f"{name}.h5"
    with pd.HDFStore(str(out), mode="w") as store:
        store.put("data", pd.DataFrame(table))
        store.put("feature_metadata", table.feature_metadata)
        store.put("sample_metadata", table.sample_metadata)
    return out


def filter_pivot(table, aligned: set[str], sample_type: str):
    mask = table.sample_metadata["sample_type"] == sample_type
    keep = [c for c in mask[mask].index if c in aligned]
    return table.subset(samples=keep), len(table.columns), len(keep)


def main():
    parser = argparse.ArgumentParser(description="Filter raw tables to aligned Primary Tumor samples")
    parser.add_argument("--types", nargs="+", default=ALL_TYPES)
    parser.add_argument("--sample-type", default="Primary Tumor")
    args = parser.parse_args()

    FILTERED_DIR.mkdir(parents=True, exist_ok=True)
    aligned = load_aligned_ids()
    print(f"Aligned cases : {len(aligned)}")
    print(f"Sample type   : {args.sample_type!r}\n")

    cls_map = {
        "expression":   ExpressionTable,
        "cnv_cytoband": CopyNumberVariationTable,
        "cnv_gene":     CopyNumberVariationTable,
        "methylation":  PivotTable,
    }

    for dtype in args.types:
        print(f"=== {dtype} ===")
        try:
            if dtype == "mutation":
                maf_df = pd.read_csv(TABLE_DIR / "mutation.maf", sep="\t",
                                     comment="#", low_memory=False)
                before = maf_df["sample_ID"].nunique()
                maf_df = maf_df[
                    (maf_df["sample_type"] == args.sample_type) &
                    (maf_df["sample_ID"].isin(aligned))
                ]
                after = maf_df["sample_ID"].nunique()
                maf = MAF(maf_df)
                maf.to_MAF(FILTERED_DIR / "mutation.maf")
                save_pivot(maf.to_pivot_table(), "mutation")
                print(f"  {before} → {after} cases, saved MAF + pivot")

            else:
                cls = cls_map.get(dtype, PivotTable)
                table = load_pivot(dtype, cls)
                filtered, before, after = filter_pivot(table, aligned, args.sample_type)
                save_pivot(filtered, dtype)
                print(f"  {before} → {after} samples")

        except Exception as e:
            print(f"  ERROR: {e}")
        print()

    print("Filtered tables:")
    for f in sorted(FILTERED_DIR.iterdir()):
        print(f"  {f.name} ({f.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()

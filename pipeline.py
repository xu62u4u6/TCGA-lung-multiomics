"""
TCGA lung multi-omics pipeline.

Steps (in order):
    generate_manifests   Query GDC API → full per-dtype manifests
    align_manifests      Find case intersection → aligned manifests
    download             Download via gdc-client (resume-safe)
    build_tables         Raw files → HDF5 tables
    filter_tables        Filter to aligned Primary Tumor samples
    build_cohort         Assemble Cohort → HDF5

Usage:
    uv run python pipeline.py run
    uv run python pipeline.py run --from download
    uv run python pipeline.py run --steps build_tables filter_tables build_cohort
    uv run python pipeline.py status
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path

STEPS = [
    "generate_manifests",
    "align_manifests",
    "download",
    "build_tables",
    "filter_tables",
    "build_cohort",
]

SCRIPTS = {step: Path(f"scripts/{step}.py") for step in STEPS}

# Output paths used to determine if a step is complete
STEP_OUTPUTS = {
    "generate_manifests": Path("data/manifests/full"),
    "align_manifests":    Path("data/manifests/aligned"),
    "download":           Path("data/download_log.tsv"),
    "build_tables":       Path("data/processed/tables"),
    "filter_tables":      Path("data/processed/filtered"),
    "build_cohort":       Path("data/processed/tcga_lung.h5"),
}


def step_done(step: str) -> bool:
    out = STEP_OUTPUTS.get(step)
    if out is None:
        return False
    if out.is_dir():
        return any(out.iterdir())
    return out.exists()


def run_step(step: str, extra_args: list[str] = None) -> bool:
    script = SCRIPTS[step]
    cmd = [sys.executable, str(script)] + (extra_args or [])
    print(f"\n{'='*60}")
    print(f"  {step}")
    print(f"{'='*60}")
    t0 = time.time()
    result = subprocess.run(cmd)
    elapsed = time.time() - t0
    if result.returncode != 0:
        print(f"\n  FAILED (exit {result.returncode}) after {elapsed:.0f}s")
        return False
    print(f"\n  done in {elapsed:.0f}s")
    return True


def show_status():
    print(f"\n{'step':<22} {'status'}")
    print("-" * 40)
    for step in STEPS:
        mark = "✓" if step_done(step) else "·"
        print(f"  {mark}  {step}")
    print()


def main():
    parser = argparse.ArgumentParser(description="TCGA lung multi-omics pipeline")
    sub = parser.add_subparsers(dest="command")

    run_p = sub.add_parser("run", help="Run pipeline steps")
    group = run_p.add_mutually_exclusive_group()
    group.add_argument("--from", dest="from_step", metavar="STEP",
                       choices=STEPS, help="Start from this step")
    group.add_argument("--steps", nargs="+", metavar="STEP",
                       choices=STEPS, help="Run only these steps")

    sub.add_parser("status", help="Show pipeline progress")

    args = parser.parse_args()

    if args.command == "status" or args.command is None:
        show_status()
        if args.command is None:
            parser.print_help()
        return

    if args.command == "run":
        if args.steps:
            to_run = args.steps
        elif args.from_step:
            idx = STEPS.index(args.from_step)
            to_run = STEPS[idx:]
        else:
            to_run = STEPS

        print(f"Steps to run: {', '.join(to_run)}")

        for step in to_run:
            ok = run_step(step)
            if not ok:
                print(f"\nPipeline stopped at '{step}'.")
                sys.exit(1)

        print("\nPipeline complete.")
        show_status()


if __name__ == "__main__":
    main()

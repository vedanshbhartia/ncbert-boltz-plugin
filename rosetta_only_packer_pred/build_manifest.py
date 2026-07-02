#!/usr/bin/env python3
"""Build a jobs.tsv for the Rosetta-only packer de novo design protocol.

Unlike the threading workflow (build_1a85_target_manifest.py /
build_explicit_manifest.py), this protocol does NOT thread a fixed sequence —
rosetta_only_packer_pred/design_script.xml designs all 18 chain-A positions
from scratch with the composition + net-charge constraints baked in. So the
manifest only needs to map each input backbone to a job; the ``pepseq`` column
is a harmless poly-glycine placeholder (the XML never references ``%%pepseq%%``,
but run_rosetta_batch_threadseq.sh always passes a pepseq script var).

One row per backbone. Run several designs per backbone with the runner's
``--nstruct`` flag (each nstruct is an independent stochastic design
trajectory).

Inputs:
  --backbones-dir   Directory of input backbone PDBs (chain A = peptide,
                    chain B = MMP8 target).
  --glob            Filename glob within that dir (default: *.pdb).
  --out             Output jobs.tsv path.

Output columns (consumed by scripts/run_rosetta_batch_threadseq.sh):
  job_id, source_file, line_no, pepseq, input_pdb
"""

import argparse
import csv
import re
from pathlib import Path

PLACEHOLDER_PEPSEQ = "G" * 18  # not used by the XML; keeps the runner happy.


def natural_key(path: Path):
    """Sort PerturbN backbones numerically (Perturb2 before Perturb10)."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", path.name)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--backbones-dir",
        type=Path,
        default=Path("perturbed_backbones/perturbed_backbones"),
        help="Directory of input backbone PDBs (default: %(default)s).",
    )
    parser.add_argument(
        "--glob", default="*.pdb", help="Filename glob (default: %(default)s)."
    )
    parser.add_argument("--out", type=Path, required=True, help="Output jobs.tsv path.")
    parser.add_argument(
        "--source-file",
        default="rosetta_only_packer_design",
        help="Value for the source_file column (default: %(default)s).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.backbones_dir.is_dir():
        raise SystemExit(f"Backbones dir not found: {args.backbones_dir}")

    backbones = sorted(args.backbones_dir.glob(args.glob), key=natural_key)
    if not backbones:
        raise SystemExit(
            f"No backbones matched {args.glob!r} in {args.backbones_dir}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["job_id", "source_file", "line_no", "pepseq", "input_pdb"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for i, pdb in enumerate(backbones, start=1):
            writer.writerow(
                {
                    "job_id": pdb.stem,
                    "source_file": args.source_file,
                    "line_no": str(i),
                    "pepseq": PLACEHOLDER_PEPSEQ,
                    "input_pdb": str(pdb),
                }
            )

    print(f"Wrote {len(backbones)} job(s) to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

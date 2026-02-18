#!/usr/bin/env python3
"""Build a Rosetta threading job manifest from peptide sequence text files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Input sequence text files",
    )
    parser.add_argument("--out", required=True, help="Output jobs TSV path")
    parser.add_argument(
        "--expected-length",
        type=int,
        default=18,
        help="Expected peptide length counting bracket blocks as one residue",
    )
    return parser.parse_args()


def count_residues(peptide: str) -> int:
    count = 0
    i = 0
    while i < len(peptide):
        char = peptide[i]
        if char == "[":
            close = peptide.find("]", i + 1)
            if close == -1:
                raise ValueError(f"Unclosed bracket in sequence: {peptide}")
            token = peptide[i + 1 : close]
            if not token:
                raise ValueError(f"Empty bracket token in sequence: {peptide}")
            count += 1
            i = close + 1
            continue
        if not char.isspace():
            count += 1
        i += 1
    return count


def manifest_prefix(path: Path) -> str:
    stem = path.stem
    suffix = "_backbone_design"
    if stem.endswith(suffix):
        stem = stem[: -len(suffix)]
    return stem


def main() -> int:
    args = parse_args()
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    for input_name in args.inputs:
        input_path = Path(input_name)
        prefix = manifest_prefix(input_path)
        with input_path.open("r", encoding="utf-8") as handle:
            for line_no, raw in enumerate(handle, start=1):
                pepseq = raw.strip()
                if not pepseq:
                    continue
                residue_count = count_residues(pepseq)
                if residue_count != args.expected_length:
                    raise ValueError(
                        f"{input_path}:{line_no} has length {residue_count}, "
                        f"expected {args.expected_length}: {pepseq}"
                    )
                job_id = f"{prefix}_L{line_no:02d}"
                rows.append(
                    {
                        "job_id": job_id,
                        "source_file": input_path.name,
                        "line_no": str(line_no),
                        "pepseq": pepseq,
                    }
                )

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["job_id", "source_file", "line_no", "pepseq"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} jobs to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

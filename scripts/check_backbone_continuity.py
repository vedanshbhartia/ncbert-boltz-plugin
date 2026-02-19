#!/usr/bin/env python3
"""Check peptide backbone continuity from PDB atom coordinates."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ResidueRecord:
    resname: str
    resseq: int
    icode: str
    n_xyz: tuple[float, float, float] | None = None
    c_xyz: tuple[float, float, float] | None = None

    @property
    def label(self) -> str:
        if self.icode.strip():
            return f"{self.resname} {self.resseq}{self.icode.strip()}"
        return f"{self.resname} {self.resseq}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, help="Input PDB path")
    parser.add_argument("--chain", default="A", help="Chain ID to check (default: A)")
    parser.add_argument(
        "--max-cn-distance",
        type=float,
        default=1.8,
        help="Maximum allowed C(i)-N(i+1) distance in Angstrom (default: 1.8)",
    )
    parser.add_argument(
        "--cyclic",
        action="store_true",
        help="Also check closure C(last)-N(first)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Print only failing pairs and summary",
    )
    parser.add_argument(
        "--include-hetatm",
        action="store_true",
        help="Include HETATM records in addition to ATOM when building the chain",
    )
    return parser.parse_args()


def parse_chain_residues(
    pdb_path: Path, chain_id: str, include_hetatm: bool
) -> list[ResidueRecord]:
    residues: list[ResidueRecord] = []
    idx_by_key: dict[tuple[int, str], int] = {}

    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM") and not (
            include_hetatm and line.startswith("HETATM")
        ):
            continue
        if len(line) < 54:
            continue
        if line[21] != chain_id:
            continue

        atom_name = line[12:16].strip()
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue

        resname = line[17:20].strip()
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        icode = line[26]
        key = (resseq, icode)

        if key not in idx_by_key:
            idx_by_key[key] = len(residues)
            residues.append(ResidueRecord(resname=resname, resseq=resseq, icode=icode))
        rec = residues[idx_by_key[key]]

        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue

        if atom_name == "N" and rec.n_xyz is None:
            rec.n_xyz = xyz
        elif atom_name == "C" and rec.c_xyz is None:
            rec.c_xyz = xyz

    return residues


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.dist(a, b)


def main() -> int:
    args = parse_args()
    pdb_path = Path(args.pdb)
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    residues = parse_chain_residues(pdb_path, args.chain, args.include_hetatm)
    if not residues:
        record_types = "ATOM/HETATM" if args.include_hetatm else "ATOM"
        raise ValueError(
            f"No {record_types} records found for chain {args.chain} in {pdb_path}"
        )
    if len(residues) < 2:
        raise ValueError(f"Need at least 2 residues in chain {args.chain} to check continuity")

    failures: list[str] = []

    def check_pair(left: ResidueRecord, right: ResidueRecord, tag: str) -> None:
        if left.c_xyz is None:
            failures.append(f"{tag}: missing C atom on {left.label}")
            return
        if right.n_xyz is None:
            failures.append(f"{tag}: missing N atom on {right.label}")
            return
        dist = distance(left.c_xyz, right.n_xyz)
        line = (
            f"{tag}: {left.label} C -> {right.label} N = {dist:.3f} A "
            f"(limit {args.max_cn_distance:.3f})"
        )
        if not args.quiet:
            print(line)
        if dist > args.max_cn_distance:
            failures.append(line)

    for i in range(len(residues) - 1):
        check_pair(residues[i], residues[i + 1], "linear")

    if args.cyclic:
        check_pair(residues[-1], residues[0], "closure")

    if failures:
        if args.quiet:
            for line in failures:
                print(line)
        print(
            f"Backbone continuity check FAILED for chain {args.chain}: "
            f"{len(failures)} issue(s)."
        )
        return 1

    print(
        f"Backbone continuity check PASSED for chain {args.chain}: "
        f"{len(residues)} residues."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

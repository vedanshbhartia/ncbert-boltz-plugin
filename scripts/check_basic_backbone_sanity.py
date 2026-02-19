#!/usr/bin/env python3
"""Basic backbone sanity checks excluding omega and closure-specific checks."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, field
from pathlib import Path


IDEAL_N_CA = 1.458
IDEAL_CA_C = 1.525
IDEAL_C_O = 1.231
IDEAL_C_N_CA = 121.7
IDEAL_N_CA_C = 111.2
IDEAL_CA_C_N = 116.2
IDEAL_CA_C_O = 120.8


@dataclass
class ResidueRecord:
    resname: str
    resseq: int
    icode: str
    atoms: dict[str, tuple[float, float, float]] = field(default_factory=dict)

    @property
    def resid(self) -> str:
        ins = self.icode.strip()
        return f"{self.resseq}{ins}" if ins else str(self.resseq)

    @property
    def label(self) -> str:
        return f"{self.resname} {self.resid}"


@dataclass
class AtomRecord:
    residue_index: int
    residue_label: str
    atom_name: str
    xyz: tuple[float, float, float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, help="Input PDB path")
    parser.add_argument("--chain", default="A", help="Chain ID to check (default: A)")
    parser.add_argument(
        "--include-hetatm",
        action="store_true",
        help="Include HETATM records in addition to ATOM",
    )
    parser.add_argument(
        "--cyclic",
        action="store_true",
        help="Treat first and last residues as covalently connected for clash exclusions",
    )
    parser.add_argument(
        "--max-nca-dev",
        type=float,
        default=0.15,
        help="Max allowed |N-CA - 1.458| in Angstrom (default: 0.15)",
    )
    parser.add_argument(
        "--max-cac-dev",
        type=float,
        default=0.15,
        help="Max allowed |CA-C - 1.525| in Angstrom (default: 0.15)",
    )
    parser.add_argument(
        "--max-co-dev",
        type=float,
        default=0.12,
        help="Max allowed |C-O - 1.231| in Angstrom (default: 0.12)",
    )
    parser.add_argument(
        "--max-cnca-dev",
        type=float,
        default=10.0,
        help="Max allowed |C-N-CA - 121.7| in degrees (default: 10.0)",
    )
    parser.add_argument(
        "--max-ncac-dev",
        type=float,
        default=12.0,
        help="Max allowed |N-CA-C - 111.2| in degrees (default: 12.0)",
    )
    parser.add_argument(
        "--max-cacn-dev",
        type=float,
        default=10.0,
        help="Max allowed |CA-C-N - 116.2| in degrees (default: 10.0)",
    )
    parser.add_argument(
        "--max-caco-dev",
        type=float,
        default=15.0,
        help="Max allowed |CA-C-O - 120.8| in degrees (default: 15.0)",
    )
    parser.add_argument(
        "--min-caca",
        type=float,
        default=2.9,
        help="Minimum allowed CA(i)-CA(i+1) distance (default: 2.9)",
    )
    parser.add_argument(
        "--max-caca",
        type=float,
        default=4.3,
        help="Maximum allowed CA(i)-CA(i+1) distance (default: 4.3)",
    )
    parser.add_argument(
        "--min-heavy-distance",
        type=float,
        default=2.0,
        help=(
            "Flag non-neighbor heavy-atom pairs below this distance in Angstrom "
            "(default: 2.0)"
        ),
    )
    parser.add_argument(
        "--skip-clash-check",
        action="store_true",
        help="Skip heavy-atom close-contact scan",
    )
    return parser.parse_args()


def infer_element(line: str, atom_name: str) -> str:
    if len(line) >= 78:
        element = line[76:78].strip()
        if element:
            return element.upper()
    for ch in atom_name:
        if ch.isalpha():
            return ch.upper()
    return ""


def parse_chain(
    pdb_path: Path, chain_id: str, include_hetatm: bool
) -> tuple[list[ResidueRecord], list[AtomRecord]]:
    residues: list[ResidueRecord] = []
    heavy_atoms: list[AtomRecord] = []
    idx_by_key: dict[tuple[int, str], int] = {}

    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        rec = line[:6].strip()
        if rec == "ATOM":
            pass
        elif rec == "HETATM" and include_hetatm:
            pass
        else:
            continue
        if len(line) < 54:
            continue
        if line[21] != chain_id:
            continue

        atom_name = line[12:16].strip()
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue

        try:
            resseq = int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        icode = line[26]
        resname = line[17:20].strip()
        key = (resseq, icode)

        if key not in idx_by_key:
            idx_by_key[key] = len(residues)
            residues.append(ResidueRecord(resname=resname, resseq=resseq, icode=icode))
        res_idx = idx_by_key[key]
        residue = residues[res_idx]

        if atom_name in residue.atoms:
            continue
        residue.atoms[atom_name] = xyz

        element = infer_element(line, atom_name)
        if element != "H":
            heavy_atoms.append(
                AtomRecord(
                    residue_index=res_idx,
                    residue_label=residue.label,
                    atom_name=atom_name,
                    xyz=xyz,
                )
            )

    return residues, heavy_atoms


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.dist(a, b)


def angle_deg(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> float:
    ba = (a[0] - b[0], a[1] - b[1], a[2] - b[2])
    bc = (c[0] - b[0], c[1] - b[1], c[2] - b[2])
    nba = math.sqrt(ba[0] ** 2 + ba[1] ** 2 + ba[2] ** 2)
    nbc = math.sqrt(bc[0] ** 2 + bc[1] ** 2 + bc[2] ** 2)
    if nba == 0 or nbc == 0:
        return float("nan")
    cos_theta = (ba[0] * bc[0] + ba[1] * bc[1] + ba[2] * bc[2]) / (nba * nbc)
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return math.degrees(math.acos(cos_theta))


def main() -> int:
    args = parse_args()
    pdb_path = Path(args.pdb)
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    residues, heavy_atoms = parse_chain(pdb_path, args.chain, args.include_hetatm)
    if not residues:
        record_types = "ATOM/HETATM" if args.include_hetatm else "ATOM"
        raise ValueError(
            f"No {record_types} records found for chain {args.chain} in {pdb_path}"
        )
    if len(residues) < 2:
        raise ValueError(f"Need at least 2 residues in chain {args.chain} to run checks")

    issues: list[str] = []

    for residue in residues:
        missing = [atom for atom in ("N", "CA", "C", "O") if atom not in residue.atoms]
        if missing:
            issues.append(f"Missing backbone atom(s) on {residue.label}: {','.join(missing)}")
            continue

        n = residue.atoms["N"]
        ca = residue.atoms["CA"]
        c = residue.atoms["C"]
        o = residue.atoms["O"]

        nca = distance(n, ca)
        nca_dev = abs(nca - IDEAL_N_CA)
        if nca_dev > args.max_nca_dev:
            issues.append(
                f"N-CA deviation on {residue.label}: {nca_dev:.3f} A (distance {nca:.3f})"
            )

        cac = distance(ca, c)
        cac_dev = abs(cac - IDEAL_CA_C)
        if cac_dev > args.max_cac_dev:
            issues.append(
                f"CA-C deviation on {residue.label}: {cac_dev:.3f} A (distance {cac:.3f})"
            )

        co = distance(c, o)
        co_dev = abs(co - IDEAL_C_O)
        if co_dev > args.max_co_dev:
            issues.append(
                f"C-O deviation on {residue.label}: {co_dev:.3f} A (distance {co:.3f})"
            )

        ncac = angle_deg(n, ca, c)
        ncac_dev = abs(ncac - IDEAL_N_CA_C)
        if ncac_dev > args.max_ncac_dev:
            issues.append(
                f"N-CA-C angle deviation on {residue.label}: {ncac_dev:.3f} deg "
                f"(angle {ncac:.3f})"
            )

        caco = angle_deg(ca, c, o)
        caco_dev = abs(caco - IDEAL_CA_C_O)
        if caco_dev > args.max_caco_dev:
            issues.append(
                f"CA-C-O angle deviation on {residue.label}: {caco_dev:.3f} deg "
                f"(angle {caco:.3f})"
            )

    for i, residue in enumerate(residues):
        prev_residue: ResidueRecord | None = None
        next_residue: ResidueRecord | None = None
        if i > 0:
            prev_residue = residues[i - 1]
        elif args.cyclic:
            prev_residue = residues[-1]

        if i < len(residues) - 1:
            next_residue = residues[i + 1]
        elif args.cyclic:
            next_residue = residues[0]

        if (
            prev_residue is not None
            and "C" in prev_residue.atoms
            and "N" in residue.atoms
            and "CA" in residue.atoms
        ):
            cnca = angle_deg(
                prev_residue.atoms["C"], residue.atoms["N"], residue.atoms["CA"]
            )
            cnca_dev = abs(cnca - IDEAL_C_N_CA)
            if cnca_dev > args.max_cnca_dev:
                issues.append(
                    f"C-N-CA angle deviation on {residue.label}: {cnca_dev:.3f} deg "
                    f"(angle {cnca:.3f})"
                )

        if (
            next_residue is not None
            and "CA" in residue.atoms
            and "C" in residue.atoms
            and "N" in next_residue.atoms
        ):
            cacn = angle_deg(
                residue.atoms["CA"], residue.atoms["C"], next_residue.atoms["N"]
            )
            cacn_dev = abs(cacn - IDEAL_CA_C_N)
            if cacn_dev > args.max_cacn_dev:
                issues.append(
                    f"CA-C-N angle deviation on {residue.label}: {cacn_dev:.3f} deg "
                    f"(angle {cacn:.3f})"
                )

    for left, right in zip(residues[:-1], residues[1:]):
        if "CA" not in left.atoms or "CA" not in right.atoms:
            continue
        d = distance(left.atoms["CA"], right.atoms["CA"])
        if d < args.min_caca or d > args.max_caca:
            issues.append(
                f"CA-CA distance out of range: {left.label} -> {right.label} = {d:.3f} A "
                f"(allowed {args.min_caca:.3f}..{args.max_caca:.3f})"
            )

    if not args.skip_clash_check:
        last_residue_index = len(residues) - 1
        for i in range(len(heavy_atoms) - 1):
            left = heavy_atoms[i]
            for j in range(i + 1, len(heavy_atoms)):
                right = heavy_atoms[j]
                if abs(left.residue_index - right.residue_index) <= 1:
                    continue
                if args.cyclic and (
                    (left.residue_index == 0 and right.residue_index == last_residue_index)
                    or (
                        right.residue_index == 0
                        and left.residue_index == last_residue_index
                    )
                ):
                    continue
                d = distance(left.xyz, right.xyz)
                if d < args.min_heavy_distance:
                    issues.append(
                        f"Possible heavy-atom clash: {left.residue_label}:{left.atom_name} "
                        f"vs {right.residue_label}:{right.atom_name} = {d:.3f} A"
                    )

    if issues:
        for issue in issues:
            print(issue)
        print(
            f"Basic backbone sanity FAILED for chain {args.chain}: "
            f"{len(issues)} issue(s)."
        )
        return 1

    print(
        f"Basic backbone sanity PASSED for chain {args.chain}: "
        f"{len(residues)} residues, {len(heavy_atoms)} heavy atoms checked."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Replicate the bond-geometry checks performed by geometry_check.m."""

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


LENGTH_THRESH = 0.2
ANGLE_THRESH = 10.0
OMEGA_THRESH = 20.0
IDEAL_LENGTHS = (1.458, 1.524, 1.329)
IDEAL_ANGLES = (121.7, 111.2, 116.2)
BOND_TYPES = ("N-CA", "CA-C", "C-N")
ANGLE_TYPES = ("C-N-CA", "N-CA-C", "CA-C-N")


@dataclass
class Residue:
    resname: str
    resseq: int
    icode: str
    n_xyz: tuple[float, float, float] | None = None
    ca_xyz: tuple[float, float, float] | None = None
    c_xyz: tuple[float, float, float] | None = None
    o_xyz: tuple[float, float, float] | None = None


@dataclass
class Issue:
    issue_type: str
    residue_index: int
    detail_type: str
    deviation: float
    message: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        required=True,
        help="PDB file or directory to scan recursively for PDB files",
    )
    parser.add_argument(
        "--chain",
        default="B",
        help="Peptide chain to evaluate (default: B, matching geometry_check.m)",
    )
    parser.add_argument(
        "--summary-out",
        default=None,
        help="Optional TSV summary path",
    )
    parser.add_argument(
        "--issues-out",
        default=None,
        help="Optional TSV listing individual issues",
    )
    parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        help=(
            "Path component to exclude from the scan (matches any directory or "
            "file name along the path). May be passed multiple times."
        ),
    )
    return parser.parse_args()


def find_pdbs(root: Path, exclude_names: set[str]) -> list[Path]:
    if root.is_file():
        return [] if exclude_names & set(root.parts) else [root]
    return sorted(
        path
        for path in root.rglob("*.pdb")
        if path.is_file() and not (exclude_names & set(path.parts))
    )


def structure_name(pdb_path: Path) -> str:
    return str(Path(pdb_path.parent.name) / pdb_path.stem)


def parse_chain_residues(pdb_path: Path, chain_id: str) -> list[Residue]:
    residues: list[Residue] = []
    idx_by_key: dict[tuple[int, str], int] = {}

    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 54 or line[21] != chain_id:
                continue

            atom_name = line[12:16].strip()
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue

            try:
                resseq = int(line[22:26])
                xyz = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            except ValueError:
                continue

            key = (resseq, line[26])
            if key not in idx_by_key:
                idx_by_key[key] = len(residues)
                residues.append(
                    Residue(
                        resname=line[17:20].strip(),
                        resseq=resseq,
                        icode=line[26],
                    )
                )
            residue = residues[idx_by_key[key]]
            if atom_name == "N" and residue.n_xyz is None:
                residue.n_xyz = xyz
            elif atom_name == "CA" and residue.ca_xyz is None:
                residue.ca_xyz = xyz
            elif atom_name == "C" and residue.c_xyz is None:
                residue.c_xyz = xyz
            elif atom_name == "O" and residue.o_xyz is None:
                residue.o_xyz = xyz

    return residues


def vector(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm(a: tuple[float, float, float]) -> float:
    return math.sqrt(dot(a, a))


def cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def angle3(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> float:
    ba = vector(a, b)
    bc = vector(c, b)
    denom = norm(ba) * norm(bc)
    if denom == 0.0:
        return float("nan")
    cos_theta = max(-1.0, min(1.0, dot(ba, bc) / denom))
    return math.degrees(math.acos(cos_theta))


def dihedral(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
    d: tuple[float, float, float],
) -> float:
    b0 = (a[0] - b[0], a[1] - b[1], a[2] - b[2])
    b1 = (c[0] - b[0], c[1] - b[1], c[2] - b[2])
    b2 = (d[0] - c[0], d[1] - c[1], d[2] - c[2])

    b1_norm = norm(b1)
    if b1_norm == 0.0:
        return float("nan")
    b1_unit = (b1[0] / b1_norm, b1[1] / b1_norm, b1[2] / b1_norm)

    v = (
        b0[0] - dot(b0, b1_unit) * b1_unit[0],
        b0[1] - dot(b0, b1_unit) * b1_unit[1],
        b0[2] - dot(b0, b1_unit) * b1_unit[2],
    )
    w = (
        b2[0] - dot(b2, b1_unit) * b1_unit[0],
        b2[1] - dot(b2, b1_unit) * b1_unit[1],
        b2[2] - dot(b2, b1_unit) * b1_unit[2],
    )
    x = dot(v, w)
    y = dot(cross(b1_unit, v), w)
    return math.degrees(math.atan2(y, x))


def require_backbone(residue: Residue, chain_id: str) -> None:
    missing = []
    if residue.n_xyz is None:
        missing.append("N")
    if residue.ca_xyz is None:
        missing.append("CA")
    if residue.c_xyz is None:
        missing.append("C")
    if residue.o_xyz is None:
        missing.append("O")
    if missing:
        joined = ", ".join(missing)
        raise ValueError(
            f"Missing backbone atom(s) at residue {residue.resseq}{residue.icode.strip()} "
            f"(chain {chain_id}): {joined}"
        )


def check_structure(pdb_path: Path, chain_id: str) -> list[Issue]:
    residues = parse_chain_residues(pdb_path, chain_id)
    if not residues:
        raise ValueError(f"No residues found for chain {chain_id} in {pdb_path}")
    if len(residues) < 2:
        raise ValueError(f"Need at least 2 residues in chain {chain_id}: {pdb_path}")

    for residue in residues:
        require_backbone(residue, chain_id)

    issues: list[Issue] = []
    complex_name = structure_name(pdb_path)
    n_res = len(residues)

    for i, residue in enumerate(residues):
        prev_residue = residues[i - 1]
        next_residue = residues[(i + 1) % n_res]

        assert residue.n_xyz and residue.ca_xyz and residue.c_xyz
        assert prev_residue.c_xyz and next_residue.n_xyz and next_residue.ca_xyz

        omega = dihedral(residue.ca_xyz, residue.c_xyz, next_residue.n_xyz, next_residue.ca_xyz)
        omega_diff = abs(omega - 180.0)
        omega_diff = min(omega_diff, 360.0 - omega_diff)
        residue_index = i + 1
        if omega_diff > OMEGA_THRESH:
            issues.append(
                Issue(
                    issue_type="omega",
                    residue_index=residue_index,
                    detail_type="omega",
                    deviation=omega_diff,
                    message=(
                        f"Omega angle deviation: structure {complex_name}, residue {residue_index}, "
                        f"deviation {omega_diff:.3f} degrees"
                    ),
                )
            )

        bond_lengths = (
            math.dist(residue.n_xyz, residue.ca_xyz),
            math.dist(residue.ca_xyz, residue.c_xyz),
            math.dist(residue.c_xyz, next_residue.n_xyz),
        )
        bond_angles = (
            angle3(prev_residue.c_xyz, residue.n_xyz, residue.ca_xyz),
            angle3(residue.n_xyz, residue.ca_xyz, residue.c_xyz),
            angle3(residue.ca_xyz, residue.c_xyz, next_residue.n_xyz),
        )

        for bond_type, actual, ideal in zip(BOND_TYPES, bond_lengths, IDEAL_LENGTHS):
            deviation = actual - ideal
            if deviation > LENGTH_THRESH:
                issues.append(
                    Issue(
                        issue_type="bond_length",
                        residue_index=residue_index,
                        detail_type=bond_type,
                        deviation=deviation,
                        message=(
                            f"Bond length deviation: structure {complex_name}, residue {residue_index}, "
                            f"bond {bond_type}, deviation {deviation:.3f} A"
                        ),
                    )
                )

        for angle_type, actual, ideal in zip(ANGLE_TYPES, bond_angles, IDEAL_ANGLES):
            deviation = actual - ideal
            if deviation > ANGLE_THRESH:
                issues.append(
                    Issue(
                        issue_type="bond_angle",
                        residue_index=residue_index,
                        detail_type=angle_type,
                        deviation=deviation,
                        message=(
                            f"Bond angle deviation: structure {complex_name}, residue {residue_index} "
                            f"angle {angle_type} deviation {deviation:.3f} degree"
                        ),
                    )
                )

    return issues


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pdb_path",
                "structure",
                "chain",
                "status",
                "issue_count",
                "omega_issues",
                "bond_length_issues",
                "bond_angle_issues",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_issues(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pdb_path",
                "structure",
                "chain",
                "issue_type",
                "residue_index",
                "detail_type",
                "deviation",
                "message",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    root = Path(args.root)
    pdb_paths = find_pdbs(root, set(args.exclude))
    if not pdb_paths:
        raise FileNotFoundError(f"No PDB files found under {root}")

    print("\n\nStart bond geometry check ...")

    summary_rows: list[dict[str, str]] = []
    issue_rows: list[dict[str, str]] = []

    for pdb_path in pdb_paths:
        issues = check_structure(pdb_path, args.chain)
        for issue in issues:
            print(issue.message)
            issue_rows.append(
                {
                    "pdb_path": str(pdb_path),
                    "structure": structure_name(pdb_path),
                    "chain": args.chain,
                    "issue_type": issue.issue_type,
                    "residue_index": str(issue.residue_index),
                    "detail_type": issue.detail_type,
                    "deviation": f"{issue.deviation:.6f}",
                    "message": issue.message,
                }
            )

        omega_issues = sum(issue.issue_type == "omega" for issue in issues)
        length_issues = sum(issue.issue_type == "bond_length" for issue in issues)
        angle_issues = sum(issue.issue_type == "bond_angle" for issue in issues)
        summary_rows.append(
            {
                "pdb_path": str(pdb_path),
                "structure": structure_name(pdb_path),
                "chain": args.chain,
                "status": "pass" if not issues else "fail",
                "issue_count": str(len(issues)),
                "omega_issues": str(omega_issues),
                "bond_length_issues": str(length_issues),
                "bond_angle_issues": str(angle_issues),
            }
        )

    if args.summary_out:
        write_summary(Path(args.summary_out), summary_rows)
    if args.issues_out:
        write_issues(Path(args.issues_out), issue_rows)

    passed = sum(row["status"] == "pass" for row in summary_rows)
    failed = len(summary_rows) - passed
    print(f"Checked {len(summary_rows)} structure(s): pass={passed}, fail={failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

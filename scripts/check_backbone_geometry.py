#!/usr/bin/env python3
"""Chain-aware backbone geometry checks for PDB files."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, field
from pathlib import Path


IDEAL_CN_LENGTH = 1.329
IDEAL_N_CA_C_ANGLE = 111.2


@dataclass
class ResidueRecord:
    chain: str
    resname: str
    resseq: int
    icode: str
    n_xyz: tuple[float, float, float] | None = None
    ca_xyz: tuple[float, float, float] | None = None
    c_xyz: tuple[float, float, float] | None = None

    @property
    def resid(self) -> str:
        ins = self.icode.strip()
        return f"{self.resseq}{ins}" if ins else str(self.resseq)

    @property
    def label(self) -> str:
        return f"{self.chain}{self.resid}"


@dataclass
class ChainSegment:
    chain: str
    residues: list[ResidueRecord] = field(default_factory=list)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, help="Input PDB path")
    parser.add_argument(
        "--chains",
        default="",
        help="Comma-separated chains to check (default: all chains found)",
    )
    parser.add_argument(
        "--cyclic-chains",
        default="",
        help="Comma-separated chains for closure checks (default: none)",
    )
    parser.add_argument(
        "--max-cn-dev",
        type=float,
        default=0.20,
        help="Max allowed |C(i)-N(i+1) - 1.329| in Angstrom (default: 0.20)",
    )
    parser.add_argument(
        "--max-ncac-dev",
        type=float,
        default=12.0,
        help="Max allowed |N-CA-C - 111.2| in degrees (default: 12.0)",
    )
    parser.add_argument(
        "--max-omega-dev",
        type=float,
        default=30.0,
        help="Max allowed omega deviation from trans/cis nearest state (default: 30.0)",
    )
    parser.add_argument(
        "--structure-id",
        default="",
        help="Structure ID used in output lines (default: PDB basename)",
    )
    return parser.parse_args()


def parse_chain_set(raw: str) -> set[str]:
    if not raw.strip():
        return set()
    return {item.strip() for item in raw.split(",") if item.strip()}


def parse_segments(pdb_path: Path) -> list[ChainSegment]:
    segments: list[ChainSegment] = []
    current: ChainSegment | None = None
    current_key: tuple[str, int, str] | None = None

    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        rec = line[:6].strip()
        if rec == "TER":
            current = None
            current_key = None
            continue
        if rec != "ATOM":
            continue
        if len(line) < 54:
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue

        chain = line[21].strip() or "_"
        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        try:
            resseq = int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        icode = line[26]
        key = (chain, resseq, icode)

        if current is None or current.chain != chain:
            current = ChainSegment(chain=chain)
            segments.append(current)
            current_key = None

        if key != current_key:
            current.residues.append(
                ResidueRecord(chain=chain, resname=resname, resseq=resseq, icode=icode)
            )
            current_key = key

        residue = current.residues[-1]
        if atom_name == "N" and residue.n_xyz is None:
            residue.n_xyz = xyz
        elif atom_name == "CA" and residue.ca_xyz is None:
            residue.ca_xyz = xyz
        elif atom_name == "C" and residue.c_xyz is None:
            residue.c_xyz = xyz

    return segments


def vector(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(
    a: tuple[float, float, float], b: tuple[float, float, float]
) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(v: tuple[float, float, float]) -> float:
    return math.sqrt(dot(v, v))


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.dist(a, b)


def angle_deg(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> float:
    ba = vector(a, b)
    bc = vector(c, b)
    cos_theta = dot(ba, bc) / (norm(ba) * norm(bc))
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return math.degrees(math.acos(cos_theta))


def dihedral_deg(
    p1: tuple[float, float, float],
    p2: tuple[float, float, float],
    p3: tuple[float, float, float],
    p4: tuple[float, float, float],
) -> float:
    b0 = vector(p2, p1)
    b1 = vector(p3, p2)
    b2 = vector(p4, p3)
    b1n = tuple(x / norm(b1) for x in b1)
    v = tuple(b0[i] - dot(b0, b1n) * b1n[i] for i in range(3))
    w = tuple(b2[i] - dot(b2, b1n) * b1n[i] for i in range(3))
    x = dot(v, w)
    y = dot(cross(b1n, v), w)
    return math.degrees(math.atan2(y, x))


def normalize_180(angle: float) -> float:
    out = ((angle + 180.0) % 360.0) - 180.0
    if out == -180.0:
        return 180.0
    return out


def omega_deviation_nearest(omega: float) -> float:
    omega_abs = abs(normalize_180(omega))
    return min(omega_abs, abs(180.0 - omega_abs))


def main() -> int:
    args = parse_args()
    pdb_path = Path(args.pdb)
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    include_chains = parse_chain_set(args.chains)
    cyclic_chains = parse_chain_set(args.cyclic_chains)
    structure_id = args.structure_id or pdb_path.stem

    segments = parse_segments(pdb_path)
    if not segments:
        raise ValueError(f"No ATOM records found in {pdb_path}")

    issues: list[str] = []

    def check_residue_angle(residue: ResidueRecord) -> None:
        if residue.n_xyz is None or residue.ca_xyz is None or residue.c_xyz is None:
            return
        angle = angle_deg(residue.n_xyz, residue.ca_xyz, residue.c_xyz)
        dev = abs(angle - IDEAL_N_CA_C_ANGLE)
        if dev > args.max_ncac_dev:
            issues.append(
                "Bond angle deviation: "
                f"structure {structure_id}, residue {residue.chain}{residue.resid}, "
                f"angle N-CA-C deviation {dev:.3f} degree "
                f"(value {angle:.3f})"
            )

    def check_bond_and_omega(left: ResidueRecord, right: ResidueRecord) -> None:
        if left.c_xyz is not None and right.n_xyz is not None:
            cn_dist = distance(left.c_xyz, right.n_xyz)
            cn_dev = abs(cn_dist - IDEAL_CN_LENGTH)
            if cn_dev > args.max_cn_dev:
                issues.append(
                    "Bond length deviation: "
                    f"structure {structure_id}, residue {left.chain}{left.resid}, "
                    f"bond C-N, deviation {cn_dev:.3f} A (distance {cn_dist:.3f})"
                )

        if (
            left.ca_xyz is not None
            and left.c_xyz is not None
            and right.n_xyz is not None
            and right.ca_xyz is not None
        ):
            omega = dihedral_deg(left.ca_xyz, left.c_xyz, right.n_xyz, right.ca_xyz)
            omega_dev = omega_deviation_nearest(omega)
            if omega_dev > args.max_omega_dev:
                issues.append(
                    "Omega angle deviation: "
                    f"structure {structure_id}, residue {left.chain}{left.resid}, "
                    f"deviation {omega_dev:.3f} degrees (omega {omega:.3f})"
                )

    for seg in segments:
        if include_chains and seg.chain not in include_chains:
            continue
        if len(seg.residues) < 1:
            continue

        for residue in seg.residues:
            check_residue_angle(residue)

        for i in range(len(seg.residues) - 1):
            check_bond_and_omega(seg.residues[i], seg.residues[i + 1])

        if seg.chain in cyclic_chains and len(seg.residues) > 1:
            check_bond_and_omega(seg.residues[-1], seg.residues[0])

    if issues:
        for line in issues:
            print(line)
        print(f"Geometry check FAILED: {len(issues)} issue(s).")
        return 1

    print("Geometry check PASSED: no issues above thresholds.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

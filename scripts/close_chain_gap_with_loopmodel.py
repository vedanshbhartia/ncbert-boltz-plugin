#!/usr/bin/env python3
"""Close a specific backbone gap in a PDB chain using Rosetta loopmodel (KIC refine)."""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


DEFAULT_ROOT = (
    Path.home()
    / "Downloads"
    / "rosetta_binary_ubuntu_3.15_bundle"
    / "rosetta.binary.ubuntu.release-408"
    / "main"
)
DEFAULT_ROSETTA_BIN = DEFAULT_ROOT / "source" / "bin" / "loopmodel.static.linuxgccrelease"
DEFAULT_ROSETTA_DB = DEFAULT_ROOT / "database"


@dataclass
class ResidueRecord:
    resname: str
    resseq: int
    icode: str
    n_xyz: tuple[float, float, float] | None = None
    c_xyz: tuple[float, float, float] | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-pdb", required=True, help="Input PDB path")
    parser.add_argument("--out-pdb", required=True, help="Output PDB path")
    parser.add_argument("--chain", default="A", help="Chain ID to repair (default: A)")
    parser.add_argument(
        "--gap-left",
        type=int,
        default=3,
        help="Residue index on C-terminal side of the gap (default: 3)",
    )
    parser.add_argument(
        "--gap-right",
        type=int,
        default=4,
        help="Residue index on N-terminal side of the gap (default: 4)",
    )
    parser.add_argument(
        "--loop-start",
        type=int,
        default=None,
        help="Loop start residue for loopmodel (default: gap-left-1, bounded to >=1)",
    )
    parser.add_argument(
        "--loop-end",
        type=int,
        default=None,
        help="Loop end residue for loopmodel (default: gap-right+1, bounded to <=N)",
    )
    parser.add_argument(
        "--max-cn-distance",
        type=float,
        default=1.8,
        help="Maximum accepted C-N distance after closure (default: 1.8)",
    )
    parser.add_argument(
        "--work-dir",
        default="runs/template_gap_closure",
        help="Working directory for loopmodel outputs",
    )
    parser.add_argument(
        "--nstruct",
        type=int,
        default=3,
        help="Number of loopmodel decoys to generate (default: 3)",
    )
    parser.add_argument(
        "--test-cycles",
        action="store_true",
        help="Run Rosetta with -run:test_cycles (faster, lower sampling)",
    )
    parser.add_argument(
        "--rosetta-bin",
        default=str(DEFAULT_ROSETTA_BIN),
        help=f"loopmodel executable (default: {DEFAULT_ROSETTA_BIN})",
    )
    parser.add_argument(
        "--rosetta-db",
        default=str(DEFAULT_ROSETTA_DB),
        help=f"Rosetta database path (default: {DEFAULT_ROSETTA_DB})",
    )
    return parser.parse_args()


def parse_chain_residues(pdb_path: Path, chain_id: str) -> list[ResidueRecord]:
    residues: list[ResidueRecord] = []
    idx_by_key: dict[tuple[int, str], int] = {}

    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM"):
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
        except ValueError:
            continue
        icode = line[26]
        key = (resseq, icode)

        if key not in idx_by_key:
            idx_by_key[key] = len(residues)
            residues.append(
                ResidueRecord(
                    resname=line[17:20].strip(),
                    resseq=resseq,
                    icode=icode,
                )
            )

        rec = residues[idx_by_key[key]]
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if atom_name == "N" and rec.n_xyz is None:
            rec.n_xyz = xyz
        elif atom_name == "C" and rec.c_xyz is None:
            rec.c_xyz = xyz

    return residues


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.dist(a, b)


def get_gap_distance(pdb_path: Path, chain: str, left: int, right: int) -> float:
    residues = parse_chain_residues(pdb_path, chain)
    by_resseq = {r.resseq: r for r in residues}
    if left not in by_resseq:
        raise ValueError(f"Residue {left} not found in chain {chain} ({pdb_path})")
    if right not in by_resseq:
        raise ValueError(f"Residue {right} not found in chain {chain} ({pdb_path})")
    left_rec = by_resseq[left]
    right_rec = by_resseq[right]
    if left_rec.c_xyz is None:
        raise ValueError(f"Missing C atom on residue {left} ({pdb_path})")
    if right_rec.n_xyz is None:
        raise ValueError(f"Missing N atom on residue {right} ({pdb_path})")
    return distance(left_rec.c_xyz, right_rec.n_xyz)


def collect_linear_distances(pdb_path: Path, chain: str) -> list[tuple[int, int, float]]:
    residues = parse_chain_residues(pdb_path, chain)
    out: list[tuple[int, int, float]] = []
    for left, right in zip(residues[:-1], residues[1:]):
        if left.c_xyz is None or right.n_xyz is None:
            continue
        out.append((left.resseq, right.resseq, distance(left.c_xyz, right.n_xyz)))
    return out


def select_best_model(
    pdbs: list[Path], chain: str, gap_left: int, gap_right: int
) -> tuple[Path, float]:
    scored: list[tuple[float, Path]] = []
    for pdb_path in pdbs:
        scored.append((get_gap_distance(pdb_path, chain, gap_left, gap_right), pdb_path))
    scored.sort(key=lambda x: x[0])
    return scored[0][1], scored[0][0]


def main() -> int:
    args = parse_args()

    in_pdb = Path(args.in_pdb)
    out_pdb = Path(args.out_pdb)
    rosetta_bin = Path(args.rosetta_bin)
    rosetta_db = Path(args.rosetta_db)
    work_dir = Path(args.work_dir)
    score_path = work_dir / "score.sc"
    log_path = work_dir / "loopmodel.log"
    loop_path = work_dir / "gap.loop"

    if not in_pdb.is_file():
        raise FileNotFoundError(f"Input PDB not found: {in_pdb}")
    if not rosetta_bin.is_file():
        raise FileNotFoundError(f"Rosetta loopmodel binary not found: {rosetta_bin}")
    if not rosetta_db.is_dir():
        raise FileNotFoundError(f"Rosetta database directory not found: {rosetta_db}")
    if args.nstruct < 1:
        raise ValueError("--nstruct must be >= 1")

    residues = parse_chain_residues(in_pdb, args.chain)
    if not residues:
        raise ValueError(f"No residues found for chain {args.chain} in {in_pdb}")
    n_res = len(residues)

    if args.loop_start is None:
        loop_start = max(1, args.gap_left - 1)
    else:
        loop_start = args.loop_start
    if args.loop_end is None:
        loop_end = min(n_res, args.gap_right + 1)
    else:
        loop_end = args.loop_end

    if not (1 <= loop_start <= args.gap_left <= loop_end):
        raise ValueError("Require loop_start <= gap_left <= loop_end")
    if not (1 <= loop_start <= args.gap_right <= loop_end):
        raise ValueError("Require loop_start <= gap_right <= loop_end")

    pre_gap = get_gap_distance(in_pdb, args.chain, args.gap_left, args.gap_right)
    print(
        f"Input gap distance ({args.gap_left}->{args.gap_right}) on chain {args.chain}: "
        f"{pre_gap:.3f} A"
    )

    work_dir.mkdir(parents=True, exist_ok=True)
    loop_path.write_text(f"LOOP {loop_start} {loop_end} {args.gap_left} 0 0\n", encoding="utf-8")

    cmd = [
        str(rosetta_bin),
        "-database",
        str(rosetta_db),
        "-in:file:s",
        str(in_pdb),
        "-in:file:fullatom",
        "-loops:loop_file",
        str(loop_path),
        "-loops:remodel",
        "no",
        "-loops:refine",
        "refine_kic",
        "-nstruct",
        str(args.nstruct),
        "-out:path:all",
        str(work_dir),
        "-out:file:scorefile",
        str(score_path),
        "-overwrite",
    ]
    if args.test_cycles:
        cmd.append("-run:test_cycles")

    with log_path.open("w", encoding="utf-8") as handle:
        result = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"loopmodel failed with exit code {result.returncode}. See {log_path}")

    models = sorted(p for p in work_dir.glob("*.pdb") if p.is_file())
    if not models:
        raise RuntimeError(f"No PDB outputs from loopmodel in {work_dir}")

    best_pdb, best_gap = select_best_model(models, args.chain, args.gap_left, args.gap_right)
    print(f"Selected model: {best_pdb.name} (gap {best_gap:.3f} A)")

    if best_gap > args.max_cn_distance:
        raise RuntimeError(
            f"Gap remains too large after loopmodel: {best_gap:.3f} A "
            f"(limit {args.max_cn_distance:.3f} A)."
        )

    linear = collect_linear_distances(best_pdb, args.chain)
    linear_bad = [x for x in linear if x[2] > args.max_cn_distance]
    if linear_bad:
        joined = ", ".join(f"{a}->{b}:{d:.3f}" for a, b, d in linear_bad)
        raise RuntimeError(
            f"Backbone continuity still broken after loopmodel (>{args.max_cn_distance:.3f} A): {joined}"
        )

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(best_pdb, out_pdb)
    print(f"Wrote repaired scaffold: {out_pdb}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

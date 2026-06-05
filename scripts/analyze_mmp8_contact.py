#!/usr/bin/env python3
"""Per-peptide-residue min distance to MMP8 for Rosetta-relaxed outputs.

For each Rosetta job directory under --run-dir, locate the relaxed model PDB
(the non-input_scaffold .pdb), parse chain A (peptide) and chain B (MMP8), and
emit a TSV with two distance metrics per peptide residue:

  - min_dist_any : min distance from any heavy atom in the residue to any
                   heavy atom in MMP8. Captures "is this residue in contact
                   with the target at all?"
  - min_dist_sc  : same but restricted to side-chain heavy atoms (everything
                   past CA except backbone C/N/O; CB included). Captures
                   "is the side chain doing binding work?" Useful for
                   deciding whether a swap to PRO/HPR/AIB sacrifices contacts.
                   GLY has no side-chain heavy atoms; min_dist_sc is NaN.

A residue is "non-contacting" by the user's experiment-1 heuristic when
min_dist_any exceeds --contact-threshold (default 5.0 A).
"""

import argparse
import csv
import math
from pathlib import Path

import numpy as np


BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT", "H", "HA", "1H", "2H", "3H"}


def parse_pdb_chain(pdb_path: Path, chain_id: str):
    """Return list of (resseq, resname, [(atom_name, xyz), ...]) for one chain.

    Heavy atoms only. Reads both ATOM and HETATM (Rosetta writes NCAAs as
    HETATM). Skips altloc != ' '/'A'.
    """
    residues: dict[tuple[int, str], dict] = {}
    order: list[tuple[int, str]] = []

    with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 54 or line[21] != chain_id:
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip() if len(line) >= 78 else atom_name[0]
            if element == "H" or atom_name.startswith("H"):
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
            icode = line[26]
            key = (resseq, icode)
            if key not in residues:
                residues[key] = {
                    "resname": line[17:20].strip(),
                    "atoms": [],
                }
                order.append(key)
            residues[key]["atoms"].append((atom_name, xyz))

    return [
        (key[0], residues[key]["resname"], residues[key]["atoms"])
        for key in order
    ]


def is_sidechain_atom(atom_name: str, resname: str) -> bool:
    if atom_name in BACKBONE_ATOMS:
        return False
    return True


def find_design_pdb(job_dir: Path) -> Path | None:
    candidates = sorted(
        p for p in job_dir.glob("*.pdb") if p.name != "input_scaffold.pdb"
    )
    return candidates[0] if candidates else None


def per_residue_distances(pdb_path: Path, peptide_chain: str, target_chain: str):
    peptide = parse_pdb_chain(pdb_path, peptide_chain)
    target = parse_pdb_chain(pdb_path, target_chain)

    if not peptide:
        raise ValueError(f"No peptide residues (chain {peptide_chain}) in {pdb_path}")
    if not target:
        raise ValueError(f"No target residues (chain {target_chain}) in {pdb_path}")

    target_coords = np.array(
        [xyz for _, _, atoms in target for _, xyz in atoms],
        dtype=float,
    )

    rows = []
    for resseq, resname, atoms in peptide:
        all_xyz = np.array([xyz for _, xyz in atoms], dtype=float)
        sc_xyz = np.array(
            [xyz for name, xyz in atoms if is_sidechain_atom(name, resname)],
            dtype=float,
        )

        all_d = np.linalg.norm(all_xyz[:, None, :] - target_coords[None, :, :], axis=2)
        min_any = float(all_d.min())

        if sc_xyz.size:
            sc_d = np.linalg.norm(sc_xyz[:, None, :] - target_coords[None, :, :], axis=2)
            min_sc = float(sc_d.min())
        else:
            min_sc = float("nan")

        rows.append(
            {
                "resseq": resseq,
                "resname": resname,
                "min_dist_any": min_any,
                "min_dist_sc": min_sc,
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Rosetta run directory containing per-design subdirectories.",
    )
    parser.add_argument(
        "--peptide-chain",
        default="A",
        help="Chain id for the peptide (default: A, matching explicit-batch outputs).",
    )
    parser.add_argument(
        "--target-chain",
        default="B",
        help="Chain id for the target protein, MMP8 (default: B).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output TSV with per-residue distances.",
    )
    parser.add_argument(
        "--contact-threshold",
        type=float,
        default=5.0,
        help="Residues with min_dist_any > this are flagged non_contact=1 in the TSV (default: 5.0 A).",
    )
    parser.add_argument(
        "--designs",
        nargs="*",
        default=None,
        help="Subset of design ids to process (default: every subdirectory under --run-dir).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.designs:
        design_dirs = [args.run_dir / d for d in args.designs]
    else:
        design_dirs = sorted(
            p for p in args.run_dir.iterdir()
            if p.is_dir() and p.name not in {"reference"} and not p.name.startswith(".")
        )

    rows_out = []
    skipped = []

    for job_dir in design_dirs:
        pdb_path = find_design_pdb(job_dir)
        if pdb_path is None:
            skipped.append((job_dir.name, "no PDB"))
            continue
        try:
            residues = per_residue_distances(pdb_path, args.peptide_chain, args.target_chain)
        except ValueError as e:
            skipped.append((job_dir.name, str(e)))
            continue
        for r in residues:
            rows_out.append(
                {
                    "design": job_dir.name,
                    "pdb_path": str(pdb_path),
                    "resseq": r["resseq"],
                    "resname": r["resname"],
                    "min_dist_any": f"{r['min_dist_any']:.3f}",
                    "min_dist_sc": "nan" if math.isnan(r["min_dist_sc"]) else f"{r['min_dist_sc']:.3f}",
                    "non_contact": "1" if r["min_dist_any"] > args.contact_threshold else "0",
                }
            )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["design", "pdb_path", "resseq", "resname", "min_dist_any", "min_dist_sc", "non_contact"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows_out)

    n_designs_done = len({r["design"] for r in rows_out})
    n_flagged = sum(1 for r in rows_out if r["non_contact"] == "1")
    print(f"Wrote {len(rows_out)} rows across {n_designs_done} designs to {args.out}")
    print(f"Flagged {n_flagged} residue(s) with min_dist_any > {args.contact_threshold} A")
    if skipped:
        print(f"Skipped {len(skipped)} design(s):")
        for name, reason in skipped:
            print(f"  {name}: {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

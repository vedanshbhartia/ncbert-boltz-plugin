#!/usr/bin/env python3
"""Per-residue Cα deviation between Boltz peptide-only prediction and the
input perturbed backbone scaffold.

For each design under --boltz-root, locate the best peptide-only Boltz model
(input_model_0.pdb by default — Boltz orders samples by confidence), find the
matching perturbed backbone PDB by extracting the Perturb<N> token from the
design id, Kabsch-align Cαs (peptide chain only), and emit per-residue
deviations.

Both structures are read from the peptide chain only. The Boltz peptide-only
output uses chain B (matching the YAML id), while the perturbed backbone has
the peptide on chain A. Residue ordering is by resseq, which matches between
the two for these 18-residue cyclic peptides.
"""

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np


PERTURB_RE = re.compile(r"(Perturb\d+)")


def parse_ca_coords(pdb_path: Path, chain_id: str) -> tuple[list[int], list[str], np.ndarray]:
    """Return parallel lists (resseq, resname) and an (N, 3) array of Cα xyz.

    Residues are returned in ascending resseq order.
    """
    by_res: dict[int, dict] = {}
    with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 54 or line[21] != chain_id:
                continue
            if line[12:16].strip() != "CA":
                continue
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
            if resseq in by_res:
                continue
            by_res[resseq] = {"resname": line[17:20].strip(), "xyz": xyz}

    ordered = sorted(by_res.items())
    resseqs = [k for k, _ in ordered]
    resnames = [v["resname"] for _, v in ordered]
    coords = np.array([v["xyz"] for _, v in ordered], dtype=float)
    return resseqs, resnames, coords


def kabsch_align(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (aligned_mobile, R) where aligned_mobile = (mobile - mobile_c) @ R + target_c.

    Standard Kabsch via SVD with reflection guard.
    """
    if mobile.shape != target.shape:
        raise ValueError(f"shape mismatch: mobile {mobile.shape} vs target {target.shape}")

    mobile_c = mobile.mean(axis=0)
    target_c = target.mean(axis=0)
    p = mobile - mobile_c
    q = target - target_c

    h = p.T @ q
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    s = np.diag([1.0, 1.0, d])
    r = u @ s @ vt
    aligned = p @ r + target_c
    return aligned, r


def find_best_boltz_model(design_dir: Path) -> tuple[Path, dict] | None:
    """Return (pdb_path, confidence_dict) for the best peptide-only Boltz model.

    Boltz writes input_model_<k>.pdb ordered by confidence (model_0 is best),
    but we re-rank by confidence_score from the JSONs to be explicit.
    """
    pred_dir = design_dir / "boltz_out" / "boltz_results_input" / "predictions" / "input"
    if not pred_dir.is_dir():
        return None
    best = None
    best_score = -float("inf")
    for conf_path in sorted(pred_dir.glob("confidence_input_model_*.json")):
        with conf_path.open() as fh:
            conf = json.load(fh)
        model_stem = conf_path.stem.replace("confidence_", "")  # e.g. input_model_0
        pdb_path = pred_dir / f"{model_stem}.pdb"
        if not pdb_path.is_file():
            continue
        score = conf.get("confidence_score", -float("inf"))
        if score > best_score:
            best_score = score
            best = (pdb_path, conf)
    return best


def find_backbone_for_design(design_id: str, backbones_dir: Path) -> Path | None:
    m = PERTURB_RE.search(design_id)
    if not m:
        return None
    candidate = backbones_dir / f"alpha_01_1a85_cyclic_001_{m.group(1)}.pdb"
    return candidate if candidate.is_file() else None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--boltz-root",
        type=Path,
        required=True,
        help="Root directory of per-design Boltz peptide-only runs (default layout: boltz/runs_peptide_only/).",
    )
    parser.add_argument(
        "--backbones-dir",
        type=Path,
        required=True,
        help="Directory containing perturbed backbone PDBs (alpha_01_1a85_cyclic_001_Perturb<N>.pdb).",
    )
    parser.add_argument(
        "--boltz-chain",
        default="B",
        help="Chain id of the peptide in Boltz output (default: B).",
    )
    parser.add_argument(
        "--backbone-chain",
        default="A",
        help="Chain id of the peptide in the perturbed backbone PDB (default: A).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output TSV with per-residue Cα deviations.",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=None,
        help="Optional second TSV with one row per design (overall RMSD, max deviation, confidence).",
    )
    parser.add_argument(
        "--deviation-threshold",
        type=float,
        default=2.0,
        help="Residues with Cα deviation > this are flagged high_dev=1 (default: 2.0 A).",
    )
    parser.add_argument(
        "--designs",
        nargs="*",
        default=None,
        help="Subset of design ids to process (default: every subdirectory under --boltz-root).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.designs:
        design_dirs = [args.boltz_root / d for d in args.designs]
    else:
        design_dirs = sorted(
            p for p in args.boltz_root.iterdir()
            if p.is_dir() and not p.name.startswith(".")
        )

    rows: list[dict] = []
    summary_rows: list[dict] = []
    skipped: list[tuple[str, str]] = []

    for design_dir in design_dirs:
        design = design_dir.name

        best = find_best_boltz_model(design_dir)
        if best is None:
            skipped.append((design, "no Boltz predictions"))
            continue
        boltz_pdb, conf = best

        backbone_pdb = find_backbone_for_design(design, args.backbones_dir)
        if backbone_pdb is None:
            skipped.append((design, "no matching perturbed backbone"))
            continue

        b_resseqs, b_resnames, b_ca = parse_ca_coords(boltz_pdb, args.boltz_chain)
        s_resseqs, s_resnames, s_ca = parse_ca_coords(backbone_pdb, args.backbone_chain)

        if b_ca.shape != s_ca.shape:
            skipped.append(
                (design, f"Cα count mismatch boltz={b_ca.shape[0]} backbone={s_ca.shape[0]}")
            )
            continue

        aligned, _ = kabsch_align(b_ca, s_ca)
        per_res_dev = np.linalg.norm(aligned - s_ca, axis=1)
        overall_rmsd = float(np.sqrt((per_res_dev ** 2).mean()))

        for i, dev in enumerate(per_res_dev):
            rows.append(
                {
                    "design": design,
                    "boltz_pdb": str(boltz_pdb),
                    "backbone_pdb": str(backbone_pdb),
                    "resseq": s_resseqs[i],
                    "backbone_resname": s_resnames[i],
                    "boltz_resname": b_resnames[i],
                    "ca_deviation": f"{float(dev):.3f}",
                    "high_dev": "1" if dev > args.deviation_threshold else "0",
                }
            )

        summary_rows.append(
            {
                "design": design,
                "overall_rmsd": f"{overall_rmsd:.3f}",
                "max_deviation": f"{float(per_res_dev.max()):.3f}",
                "n_high_dev": str(int((per_res_dev > args.deviation_threshold).sum())),
                "confidence_score": f"{conf.get('confidence_score', float('nan')):.4f}",
                "complex_plddt": f"{conf.get('complex_plddt', float('nan')):.4f}",
                "ptm": f"{conf.get('ptm', float('nan')):.4f}",
                "boltz_pdb": str(boltz_pdb),
                "backbone_pdb": str(backbone_pdb),
            }
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "design", "boltz_pdb", "backbone_pdb", "resseq",
                "backbone_resname", "boltz_resname", "ca_deviation", "high_dev",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    if args.summary_out is not None:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        with args.summary_out.open("w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "design", "overall_rmsd", "max_deviation", "n_high_dev",
                    "confidence_score", "complex_plddt", "ptm",
                    "boltz_pdb", "backbone_pdb",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(summary_rows)

    n_designs_done = len(summary_rows)
    print(f"Wrote {len(rows)} residue rows across {n_designs_done} designs to {args.out}")
    if args.summary_out is not None:
        print(f"Wrote design summary to {args.summary_out}")
    if skipped:
        print(f"Skipped {len(skipped)} design(s):")
        for name, reason in skipped:
            print(f"  {name}: {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

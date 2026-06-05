#!/usr/bin/env python3
"""Build a v4 manifest from the v3 mutation plan.

v4 uses the same candidate-position set as v3 (non-contact residues with
phi in the PRO/DPRO range, after the adjacency filter), but instead of
pre-picking the residue type the planner writes the *positions* into the
jobs.tsv's optional extra_script_vars column. The XML's design step then
lets Rosetta choose among {PRO, HYP, DPRO, DHYP, AIB} at each position.

Inputs:
  --plan-tsv        runs/proline_mutations_v3/mutation_plan.tsv
  --baseline-sel    runs/explicit_20260506_014900/selected.tsv  (provides
                    baseline pepseq + per-job input_pdb)
  --out-dir         New v4 run directory.

Output:
  <out-dir>/jobs.tsv         job_id, source_file, line_no, pepseq, input_pdb,
                             extra_script_vars  (extra column carries
                             "design_positions=7A,12A")
  <out-dir>/selected.tsv     for provenance; mirrors baseline selected.tsv
                             with a new design_positions column.

Designs with zero accepted positions in the plan are skipped — running them
would be identical to baseline.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def load_plan(path: Path) -> dict[str, list[int]]:
    """Return {design_id: [accepted positions]}."""
    by_design: dict[str, list[int]] = defaultdict(list)
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("applied") != "1":
                continue
            by_design[r["design_id"]].append(int(r["position"]))
    return {d: sorted(set(ps)) for d, ps in by_design.items()}


def load_selected(path: Path) -> dict[str, dict]:
    with path.open() as fh:
        return {r["job_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plan-tsv", type=Path, required=True)
    parser.add_argument("--baseline-sel", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--peptide-chain", default="A",
                        help="Chain id to suffix to each position (default: A).")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    positions_by_design = load_plan(args.plan_tsv)
    baseline = load_selected(args.baseline_sel)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    jobs_path = args.out_dir / "jobs.tsv"
    sel_path = args.out_dir / "selected.tsv"

    n_jobs = 0
    n_skipped = 0
    with jobs_path.open("w", newline="") as fh, sel_path.open("w", newline="") as sf:
        jw = csv.DictWriter(
            fh,
            fieldnames=["job_id", "source_file", "line_no", "pepseq", "input_pdb", "extra_script_vars"],
            delimiter="\t",
            lineterminator="\n",
        )
        jw.writeheader()
        sw = csv.DictWriter(
            sf,
            fieldnames=[
                "job_id", "source_file", "line_no", "raw_pepseq", "pepseq",
                "score", "input_pdb", "design_positions",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        sw.writeheader()

        for i, design_id in enumerate(sorted(baseline), start=1):
            positions = positions_by_design.get(design_id, [])
            if not positions:
                n_skipped += 1
                continue
            pos_str = ",".join(f"{p}{args.peptide_chain}" for p in positions)
            extra_sv = f"design_positions={pos_str}"
            sel = baseline[design_id]
            jw.writerow(
                {
                    "job_id": design_id,
                    "source_file": "v4_prolike_design.tsv",
                    "line_no": str(i),
                    "pepseq": sel["pepseq"],
                    "input_pdb": sel["input_pdb"],
                    "extra_script_vars": extra_sv,
                }
            )
            sw.writerow(
                {
                    "job_id": design_id,
                    "source_file": "v4_prolike_design.tsv",
                    "line_no": str(i),
                    "raw_pepseq": sel.get("raw_pepseq", ""),
                    "pepseq": sel["pepseq"],
                    "score": sel.get("score", ""),
                    "input_pdb": sel["input_pdb"],
                    "design_positions": pos_str,
                }
            )
            n_jobs += 1

    print(f"Wrote {n_jobs} v4 jobs to {jobs_path} (skipped {n_skipped} design(s) with no accepted positions).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

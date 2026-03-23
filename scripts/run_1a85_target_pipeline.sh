#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

LIMIT="${LIMIT:-50}"
PARALLEL="${PARALLEL:-1}"
NSTRUCT="${NSTRUCT:-1}"
RUN_TAG="${RUN_TAG:-1a85_targeted_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${RUN_DIR:-runs/$RUN_TAG}"
JOBS_TSV="$RUN_DIR/jobs.tsv"
SELECTED_TSV="$RUN_DIR/selected.tsv"
SUMMARY_TSV="$RUN_DIR/summary.tsv"
GEOM_SUMMARY_TSV="$RUN_DIR/geometry_summary.tsv"
GEOM_ISSUES_TSV="$RUN_DIR/geometry_issues.tsv"
PASS_TSV="$RUN_DIR/geometry_pass.tsv"
FAIL_TSV="$RUN_DIR/geometry_fail.tsv"

SCAFFOLD_PDB="${SCAFFOLD_PDB:-/home/vedansh/Downloads/shasha/02_12_rosetta_refine_only/SA474935_cand1_L09/alpha_001_chainA_del4_5_18mer_0001.pdb}"
ROSETTA_CONFIG="${ROSETTA_CONFIG:-$REPO_ROOT/config/rosetta.env}"

mkdir -p "$RUN_DIR"

python3 scripts/build_1a85_target_manifest.py \
  --design-dir design \
  --out "$JOBS_TSV" \
  --selected-out "$SELECTED_TSV" \
  --limit "$LIMIT"

scripts/run_rosetta_batch_threadseq.sh \
  --config "$ROSETTA_CONFIG" \
  --jobs "$JOBS_TSV" \
  --input-pdb "$SCAFFOLD_PDB" \
  --xml relax_script_thread_1a85_chainA.xml \
  --run-dir "$RUN_DIR" \
  --peptide-chain A \
  --include-hetatm-backbone-check \
  --parallel "$PARALLEL" \
  --nstruct "$NSTRUCT"

python3 scripts/collect_scores.py \
  --jobs "$JOBS_TSV" \
  --run-dir "$RUN_DIR" \
  --out "$SUMMARY_TSV"

python3 scripts/geometry_check.py \
  --root "$RUN_DIR" \
  --chain A \
  --exclude input_scaffold.pdb \
  --summary-out "$GEOM_SUMMARY_TSV" \
  --issues-out "$GEOM_ISSUES_TSV"

python3 - <<'PY' "$GEOM_SUMMARY_TSV" "$PASS_TSV" "$FAIL_TSV"
import csv
import sys
from pathlib import Path

summary_path = Path(sys.argv[1])
pass_path = Path(sys.argv[2])
fail_path = Path(sys.argv[3])

with summary_path.open("r", encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

pass_rows = [row for row in rows if row["status"] == "pass"]
fail_rows = [row for row in rows if row["status"] != "pass"]

for path, subset in ((pass_path, pass_rows), (fail_path, fail_rows)):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(subset)

print(f"Geometry pass: {len(pass_rows)}")
print(f"Geometry fail: {len(fail_rows)}")
PY

echo "Run directory: $RUN_DIR"
echo "Selected sequences: $SELECTED_TSV"
echo "Score summary: $SUMMARY_TSV"
echo "Geometry summary: $GEOM_SUMMARY_TSV"

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

LIMIT="${LIMIT:-6}"
PARALLEL="${PARALLEL:-3}"
NSTRUCT="${NSTRUCT:-1}"
RUN_TAG="${RUN_TAG:-design_relax_variant_sa474935_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${RUN_DIR:-runs/$RUN_TAG}"
ROSETTA_CONFIG="${ROSETTA_CONFIG:-$REPO_ROOT/config/rosetta.env}"
SCAFFOLD_PDB="${SCAFFOLD_PDB:-$REPO_ROOT/SA474935_cand1_backbone.pdb}"
JOBS_TSV="$RUN_DIR/jobs.tsv"
SELECTED_TSV="$RUN_DIR/selected.tsv"
SUMMARY_TSV="$RUN_DIR/variant_summary.tsv"

VARIANTS=(
  "baseline:relax_script_thread_pep.xml"
  "rigid_jump:relax_script_thread_pep_rigid_jump.xml"
  "staged_backbone:relax_script_thread_pep_staged.xml"
  "staged_backbone_rigid_jump:relax_script_thread_pep_staged_rigid_jump.xml"
)

mkdir -p "$RUN_DIR"

python3 scripts/build_1a85_target_manifest.py \
  --design-dir design \
  --out "$JOBS_TSV" \
  --selected-out "$SELECTED_TSV" \
  --limit "$LIMIT"

write_pass_fail_tables() {
  local geom_summary_tsv="$1"
  local pass_tsv="$2"
  local fail_tsv="$3"

  python3 - <<'PY' "$geom_summary_tsv" "$pass_tsv" "$fail_tsv"
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
PY
}

summarize_variant() {
  local variant_name="$1"
  local variant_dir="$2"
  local out_tsv="$3"

  python3 - <<'PY' "$variant_name" "$variant_dir" "$out_tsv"
import csv
import sys
from pathlib import Path

variant_name = sys.argv[1]
variant_dir = Path(sys.argv[2])
out_tsv = Path(sys.argv[3])
summary_path = variant_dir / "summary.tsv"
geom_path = variant_dir / "geometry_summary.tsv"

with summary_path.open("r", encoding="utf-8", newline="") as handle:
    summary_rows = list(csv.DictReader(handle, delimiter="\t"))
with geom_path.open("r", encoding="utf-8", newline="") as handle:
    geom_rows = list(csv.DictReader(handle, delimiter="\t"))

row = {
    "variant": variant_name,
    "total_jobs": str(len(summary_rows)),
    "rosetta_success": str(sum(r["status"] == "success" for r in summary_rows)),
    "rosetta_failed": str(sum(r["status"] != "success" for r in summary_rows)),
    "geometry_checked": str(len(geom_rows)),
    "geometry_pass": str(sum(r["status"] == "pass" for r in geom_rows)),
    "geometry_fail": str(sum(r["status"] != "pass" for r in geom_rows)),
}

header = [
    "variant",
    "total_jobs",
    "rosetta_success",
    "rosetta_failed",
    "geometry_checked",
    "geometry_pass",
    "geometry_fail",
]

write_header = not out_tsv.exists()
out_tsv.parent.mkdir(parents=True, exist_ok=True)
with out_tsv.open("a", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
    if write_header:
        writer.writeheader()
    writer.writerow(row)
PY
}

for variant_spec in "${VARIANTS[@]}"; do
  variant_name="${variant_spec%%:*}"
  variant_xml="${variant_spec#*:}"
  variant_dir="$RUN_DIR/$variant_name"
  geom_summary_tsv="$variant_dir/geometry_summary.tsv"
  geom_issues_tsv="$variant_dir/geometry_issues.tsv"
  pass_tsv="$variant_dir/geometry_pass.tsv"
  fail_tsv="$variant_dir/geometry_fail.tsv"

  mkdir -p "$variant_dir"
  echo "=== Running variant: $variant_name ($variant_xml) ==="

  scripts/run_rosetta_batch_threadseq.sh \
    --config "$ROSETTA_CONFIG" \
    --jobs "$JOBS_TSV" \
    --input-pdb "$SCAFFOLD_PDB" \
    --xml "$variant_xml" \
    --run-dir "$variant_dir" \
    --peptide-chain B \
    --include-hetatm-backbone-check \
    --parallel "$PARALLEL" \
    --nstruct "$NSTRUCT"

  python3 scripts/collect_scores.py \
    --jobs "$JOBS_TSV" \
    --run-dir "$variant_dir" \
    --out "$variant_dir/summary.tsv"

  python3 scripts/geometry_check.py \
    --root "$variant_dir" \
    --chain B \
    --exclude input_scaffold.pdb \
    --summary-out "$geom_summary_tsv" \
    --issues-out "$geom_issues_tsv"

  write_pass_fail_tables "$geom_summary_tsv" "$pass_tsv" "$fail_tsv"
  summarize_variant "$variant_name" "$variant_dir" "$SUMMARY_TSV"
done

echo
echo "Selected sequences: $SELECTED_TSV"
echo "Pilot summary: $SUMMARY_TSV"
cat "$SUMMARY_TSV"

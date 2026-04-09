#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

ROSETTA_CONFIG="${ROSETTA_CONFIG:-$REPO_ROOT/config/rosetta.env}"
SEQS_PER_BACKBONE="${SEQS_PER_BACKBONE:-3}"
PARALLEL="${PARALLEL:-3}"
NSTRUCT="${NSTRUCT:-1}"
RUN_TAG="${RUN_TAG:-relax_variant_pilot_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${RUN_DIR:-runs/$RUN_TAG}"
SUMMARY_TSV="$RUN_DIR/variant_summary.tsv"

PAIR_NAMES=(
  "SA474935_cand1"
  "SA2989311_cand2"
)
PAIR_DESIGNS=(
  "SA474935_cand1_backbone_design.txt"
  "SA2989311_cand2_backbone_design.txt"
)
PAIR_PDBS=(
  "SA474935_cand1_backbone.pdb"
  "SA2989311_cand2_backbone.pdb"
)

VARIANTS=(
  "baseline:relax_script_thread_pep.xml"
  "rigid_jump:relax_script_thread_pep_rigid_jump.xml"
  "staged_backbone:relax_script_thread_pep_staged.xml"
  "staged_backbone_rigid_jump:relax_script_thread_pep_staged_rigid_jump.xml"
)

mkdir -p "$RUN_DIR"

limit_manifest() {
  local input_tsv="$1"
  local output_tsv="$2"
  local limit="$3"

  python3 - <<'PY' "$input_tsv" "$output_tsv" "$limit"
import csv
import sys
from pathlib import Path

input_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
limit = int(sys.argv[3])

with input_path.open("r", encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

output_path.parent.mkdir(parents=True, exist_ok=True)
with output_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["job_id", "source_file", "line_no", "pepseq"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows[:limit])
PY
}

summarize_variant() {
  local variant_dir="$1"
  local variant_name="$2"
  local out_tsv="$3"

  python3 - <<'PY' "$variant_dir" "$variant_name" "$out_tsv"
import csv
import sys
from pathlib import Path

variant_dir = Path(sys.argv[1])
variant_name = sys.argv[2]
out_tsv = Path(sys.argv[3])

row = {
    "variant": variant_name,
    "pair_count": "0",
    "total_jobs": "0",
    "rosetta_success": "0",
    "rosetta_failed": "0",
    "geometry_checked": "0",
    "geometry_pass": "0",
    "geometry_fail": "0",
}

pair_count = 0
for pair_dir in sorted(path for path in variant_dir.iterdir() if path.is_dir()):
    summary_path = pair_dir / "summary.tsv"
    geom_path = pair_dir / "geometry_summary.tsv"
    if not summary_path.exists():
      continue
    pair_count += 1
    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    row["total_jobs"] = str(int(row["total_jobs"]) + len(summary_rows))
    row["rosetta_success"] = str(
        int(row["rosetta_success"]) + sum(r["status"] == "success" for r in summary_rows)
    )
    row["rosetta_failed"] = str(
        int(row["rosetta_failed"]) + sum(r["status"] != "success" for r in summary_rows)
    )

    if geom_path.exists():
        with geom_path.open("r", encoding="utf-8", newline="") as handle:
            geom_rows = list(csv.DictReader(handle, delimiter="\t"))
        row["geometry_checked"] = str(int(row["geometry_checked"]) + len(geom_rows))
        row["geometry_pass"] = str(
            int(row["geometry_pass"]) + sum(r["status"] == "pass" for r in geom_rows)
        )
        row["geometry_fail"] = str(
            int(row["geometry_fail"]) + sum(r["status"] != "pass" for r in geom_rows)
        )

row["pair_count"] = str(pair_count)

header = [
    "variant",
    "pair_count",
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
  mkdir -p "$variant_dir"

  echo "=== Running variant: $variant_name ($variant_xml) ==="

  for idx in "${!PAIR_NAMES[@]}"; do
    pair_name="${PAIR_NAMES[$idx]}"
    design_file="${PAIR_DESIGNS[$idx]}"
    backbone_pdb="${PAIR_PDBS[$idx]}"
    pair_dir="$variant_dir/$pair_name"
    manifest_dir="$variant_dir/manifests"
    jobs_all_tsv="$manifest_dir/${pair_name}_jobs_all.tsv"
    jobs_tsv="$manifest_dir/${pair_name}_jobs.tsv"

    mkdir -p "$pair_dir"
    mkdir -p "$manifest_dir"

    python3 scripts/build_job_manifest.py \
      --inputs "$design_file" \
      --out "$jobs_all_tsv" \
      --expected-length 18

    limit_manifest "$jobs_all_tsv" "$jobs_tsv" "$SEQS_PER_BACKBONE"

    scripts/run_rosetta_batch.sh \
      --config "$ROSETTA_CONFIG" \
      --jobs "$jobs_tsv" \
      --input-pdb "$backbone_pdb" \
      --xml "$variant_xml" \
      --run-dir "$pair_dir" \
      --peptide-chain B \
      --include-hetatm-backbone-check \
      --parallel "$PARALLEL" \
      --nstruct "$NSTRUCT"

    python3 scripts/collect_scores.py \
      --jobs "$jobs_tsv" \
      --run-dir "$pair_dir" \
      --out "$pair_dir/summary.tsv"

    python3 scripts/geometry_check.py \
      --root "$pair_dir" \
      --chain B \
      --exclude input_scaffold.pdb \
      --summary-out "$pair_dir/geometry_summary.tsv" \
      --issues-out "$pair_dir/geometry_issues.tsv"
  done

  summarize_variant "$variant_dir" "$variant_name" "$SUMMARY_TSV"
done

echo
echo "Pilot summary: $SUMMARY_TSV"
cat "$SUMMARY_TSV"

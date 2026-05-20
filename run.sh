#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

LIMIT="${LIMIT:-6}"
PARALLEL="${PARALLEL:-3}"
NSTRUCT="${NSTRUCT:-1}"
RUN_TAG="${RUN_TAG:-sa474935_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${RUN_DIR:-runs/$RUN_TAG}"
ROSETTA_CONFIG="${ROSETTA_CONFIG:-$SCRIPT_DIR/config/rosetta.env}"
SCAFFOLD_PDB="${SCAFFOLD_PDB:-$SCRIPT_DIR/SA474935_cand1_backbone.pdb}"
JOBS_TSV="$RUN_DIR/jobs.tsv"
SELECTED_TSV="$RUN_DIR/selected.tsv"
REF_DIR="$RUN_DIR/reference"
REF_JOBS_TSV="$REF_DIR/jobs.tsv"

mkdir -p "$RUN_DIR" "$REF_DIR"

# --- Reference backbone run (no threading) ---
printf 'job_id\tsource_file\tline_no\tpepseq\nreference\tbackbone\t0\t\n' > "$REF_JOBS_TSV"

scripts/run_rosetta_batch.sh \
  --config "$ROSETTA_CONFIG" \
  --jobs "$REF_JOBS_TSV" \
  --input-pdb "$SCAFFOLD_PDB" \
  --xml relax_script_reference_pep.xml \
  --run-dir "$REF_DIR" \
  --peptide-chain B \
  --include-hetatm-backbone-check \
  --parallel 1 \
  --nstruct 1

REF_PDB="$(ls "$REF_DIR/reference/"*.pdb 2>/dev/null | head -1)"
REFERENCE_DG="$(python3 scripts/extract_dg_separated.py "${REF_PDB:-}")"

# --- Design jobs ---
python3 scripts/build_1a85_target_manifest.py \
  --design-dir design \
  --out "$JOBS_TSV" \
  --selected-out "$SELECTED_TSV" \
  --limit "$LIMIT"

scripts/run_rosetta_batch_threadseq.sh \
  --config "$ROSETTA_CONFIG" \
  --jobs "$JOBS_TSV" \
  --xml relax_script_thread_1a85_chainA_staged_rigid_jump.xml \
  --run-dir "$RUN_DIR" \
  --peptide-chain A \
  --include-hetatm-backbone-check \
  --parallel "$PARALLEL" \
  --nstruct "$NSTRUCT"

COLLECT_ARGS=(
  --jobs "$JOBS_TSV"
  --run-dir "$RUN_DIR"
  --out "$RUN_DIR/summary.tsv"
)
if [[ -n "${REFERENCE_DG:-}" ]]; then
  COLLECT_ARGS+=(--reference-dg "$REFERENCE_DG")
fi
python3 scripts/collect_scores.py "${COLLECT_ARGS[@]}"

python3 scripts/geometry_check.py \
  --root "$RUN_DIR" \
  --chain A \
  --exclude input_scaffold.pdb \
  --exclude reference \
  --summary-out "$RUN_DIR/geometry_summary.tsv" \
  --issues-out "$RUN_DIR/geometry_issues.tsv"

echo
echo "Selected sequences: $SELECTED_TSV"
echo "Summary:            $RUN_DIR/summary.tsv"
echo "Geometry:           $RUN_DIR/geometry_summary.tsv"

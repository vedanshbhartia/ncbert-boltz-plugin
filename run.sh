#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROSETTA_CONFIG="${ROSETTA_CONFIG:-$SCRIPT_DIR/config/rosetta.env}"

run_pair() {
  local design_txt="$1"
  local backbone_pdb="$2"
  local run_tag="$3"
  local jobs_tsv="runs/jobs_${run_tag}.tsv"
  local run_dir="runs/full_run_${run_tag}"

  python3 scripts/build_job_manifest.py \
    --inputs "$design_txt" \
    --out "$jobs_tsv" \
    --expected-length 18

  python3 scripts/check_backbone_continuity.py \
    --pdb "$backbone_pdb" \
    --chain B \
    --cyclic \
    --include-hetatm \
    --max-cn-distance 1.8

  local -a batch_cmd=(
    scripts/run_rosetta_batch.sh
    --config "$ROSETTA_CONFIG"
    --jobs "$jobs_tsv"
    --input-pdb "$backbone_pdb"
    --xml relax_script_thread_pep.xml
    --run-dir "$run_dir"
    --peptide-chain B
    --include-hetatm-backbone-check
    --parallel 8
    --nstruct 1
  )
  if [[ -n "${EXTRA_RES_FA_DIR:-}" ]]; then
    batch_cmd+=(--extra-res-fa-dir "$EXTRA_RES_FA_DIR")
  fi
  "${batch_cmd[@]}"

  python3 scripts/collect_scores.py \
    --jobs "$jobs_tsv" \
    --run-dir "$run_dir" \
    --out "$run_dir/summary.tsv"
}

run_pair \
  SA474935_cand1_backbone_design.txt \
  SA474935_cand1_backbone.pdb \
  SA474935_cand1

run_pair \
  SA2989311_cand2_backbone_design.txt \
  SA2989311_cand2_backbone.pdb \
  SA2989311_cand2

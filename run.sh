set -euo pipefail

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

  scripts/run_rosetta_batch.sh \
    --jobs "$jobs_tsv" \
    --input-pdb "$backbone_pdb" \
    --xml relax_script_thread_pep.xml \
    --run-dir "$run_dir" \
    --peptide-chain B \
    --include-hetatm-backbone-check \
    --parallel 8 \
    --nstruct 1

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

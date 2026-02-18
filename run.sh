python3 scripts/build_job_manifest.py \
  --inputs SA474935_cand1_backbone_design.txt SA2989311_cand2_backbone_design.txt \
  --out runs/jobs.tsv \
  --expected-length 18

python3 scripts/prepare_trimmed_input.py \
  --in-cif alpha_001_1a85_cyclic_199.cif \
  --out-pdb inputs/alpha_001_chainA_del4_5_18mer_raw.pdb \
  --chain A \
  --delete 4,5 \
  --expect-old-seq DHSSSGGGTFTDSSTGNIAY \
  --expect-new-seq DHSGGGTFTDSSTGNIAY

python3 scripts/close_chain_gap_with_loopmodel.py \
  --in-pdb inputs/alpha_001_chainA_del4_5_18mer_raw.pdb \
  --out-pdb inputs/alpha_001_chainA_del4_5_18mer.pdb \
  --chain A \
  --gap-left 3 \
  --gap-right 4 \
  --loop-start 2 \
  --loop-end 5 \
  --nstruct 3 \
  --test-cycles \
  --work-dir runs/template_gap_closure

python3 scripts/check_backbone_continuity.py \
  --pdb inputs/alpha_001_chainA_del4_5_18mer.pdb \
  --chain A \
  --cyclic \
  --max-cn-distance 1.8

scripts/run_rosetta_batch.sh \
  --jobs runs/jobs.tsv \
  --input-pdb inputs/alpha_001_chainA_del4_5_18mer.pdb \
  --xml relax_script_thread_pep.xml \
  --run-dir runs/full_run \
  --parallel 8 \
  --nstruct 1

python3 scripts/collect_scores.py \
  --jobs runs/jobs.tsv \
  --run-dir runs/full_run \
  --out runs/full_run/summary.tsv

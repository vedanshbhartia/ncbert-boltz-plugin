#!/bin/bash
# SLURM array submit for the rosetta_only_packer_pred de novo design protocol.
# One array task per backbone; each task runs --nstruct 5 (5 designs).
# Array size 13 x nstruct 5 = 65 designs (~"about 64").
#
# Submit with the wrapper at the bottom of this file's header comment, which
# pre-builds $RUN_DIR/jobs64.tsv and exports RUN_DIR into the job.
#
#   cd <repo>/lib
#   mkdir -p slurm_logs
#   RUN_DIR=runs/packer_design_torch_$(date +%Y%m%d_%H%M%S); mkdir -p "$RUN_DIR"
#   python3 rosetta_only_packer_pred/build_manifest.py \
#       --backbones-dir perturbed_backbones/perturbed_backbones --out "$RUN_DIR/jobs.tsv"
#   head -n 14 "$RUN_DIR/jobs.tsv" > "$RUN_DIR/jobs64.tsv"   # header + 13 backbones
#   N=$(($(wc -l < "$RUN_DIR/jobs64.tsv") - 1))
#   sbatch --array=1-$N --export=ALL,RUN_DIR="$RUN_DIR" rosetta_only_packer_pred/run_torch_array.sbatch

#SBATCH --job-name=pkr_design
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1          # Rosetta packer/relax is single-threaded per structure
#SBATCH --mem=8G
#SBATCH --time=12:00:00            # 5 nstruct, several min..~1h each -> headroom
#SBATCH --output=slurm_logs/pkr_%A_%a.out
#SBATCH --error=slurm_logs/pkr_%A_%a.err
### ---- EDIT for your torch allocation ----
# --account is REQUIRED on torch: the default `users` account is rejected. Get a
# project account (torch_pr_xxx_yyy) from your PI via projects.hpc.nyu.edu, then
# pass it here or as `sbatch --account=...`. Check yours with:
#   sacctmgr -n -P show assoc user=$USER format=account
###SBATCH --account=torch_pr_xxx_yyy
#SBATCH --partition=cpu_short

set -euo pipefail

# Run from the repo root (lib/) so the XML's repo-root-relative .comp/.charge
# paths resolve. SLURM_SUBMIT_DIR is wherever you ran `sbatch`.
cd "${SLURM_SUBMIT_DIR:?submit from the repo lib/ directory}"

# ---- Rosetta on the cluster (EDIT) -------------------------------------------
# Option A: copy the static binary + database to the cluster and point here:
export ROSETTA_BIN="${ROSETTA_BIN:-$HOME/rosetta/rosetta_scripts.static.linuxgccrelease}"
export ROSETTA_DB="${ROSETTA_DB:-$HOME/rosetta/database}"
# Option B: if torch provides a Rosetta module, instead do e.g.:
#   module purge && module load rosetta/3.14
#   export ROSETTA_BIN=$(command -v rosetta_scripts.static.linuxgccrelease)
#   export ROSETTA_DB=$ROSETTA_ROOT/main/database
# ------------------------------------------------------------------------------

: "${RUN_DIR:?export RUN_DIR=... when submitting (see header)}"

# This task's job row = data line SLURM_ARRAY_TASK_ID of the pre-built manifest.
TASK_JOBS="$RUN_DIR/jobs_task${SLURM_ARRAY_TASK_ID}.tsv"
{ head -n 1 "$RUN_DIR/jobs64.tsv"
  sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$RUN_DIR/jobs64.tsv"; } > "$TASK_JOBS"

scripts/run_rosetta_batch_threadseq.sh \
  --jobs "$TASK_JOBS" \
  --xml rosetta_only_packer_pred/design_script.xml \
  --run-dir "$RUN_DIR" \
  --peptide-chain A \
  --rosetta-bin "$ROSETTA_BIN" \
  --rosetta-db "$ROSETTA_DB" \
  --extra-res-fa-dir ncaa_params \
  --parallel 1 --nstruct 5

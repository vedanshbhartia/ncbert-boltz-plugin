# Running the rosetta-only packer design

Two ways to run [design_script.xml](design_script.xml): locally with the runner's
own process parallelism, or on the NYU **torch** cluster as a SLURM job array.
Both designs all 18 chain-A positions de novo and emit only designs that pass the
hard composition gates (≥1 proline-like, ≥1 acidic, ≤2 Arg, net charge ≤ −1).

> Run everything from the repo root (`lib/`). The XML resolves its `.comp` /
> `.charge` paths against the launch cwd.

---

## Local run (single machine, N jobs in parallel)

Parallelism here is the runner's `--parallel` (that many backbone jobs run as
concurrent processes on this node); each job runs `--nstruct` designs
sequentially. Example: 10 backbones × `--nstruct 5` = 50 designs.

```bash
cd <repo>/lib

RUN_DIR=runs/packer_design_local_$(date +%Y%m%d_%H%M%S)
mkdir -p "$RUN_DIR"

# Manifest of the first 10 backbones (header + 10 rows).
python3 rosetta_only_packer_pred/build_manifest.py \
  --backbones-dir perturbed_backbones/perturbed_backbones --out "$RUN_DIR/jobs.tsv"
head -n 11 "$RUN_DIR/jobs.tsv" > "$RUN_DIR/jobs10.tsv"

# Run them, 10 at a time, 5 designs each.
scripts/run_rosetta_batch_threadseq.sh \
  --config config/rosetta.env \
  --jobs "$RUN_DIR/jobs10.tsv" \
  --xml rosetta_only_packer_pred/design_script.xml \
  --run-dir "$RUN_DIR" \
  --peptide-chain A \
  --parallel 10 --nstruct 5
```

Collect afterward:

```bash
python3 scripts/collect_scores.py --jobs "$RUN_DIR/jobs10.tsv" --run-dir "$RUN_DIR" --out "$RUN_DIR/summary.tsv"
python3 scripts/geometry_check.py --root "$RUN_DIR" --chain A --exclude input_scaffold.pdb \
  --summary-out "$RUN_DIR/geometry_summary.tsv" --issues-out "$RUN_DIR/geometry_issues.tsv"
```

`--parallel` should be ≲ the number of free cores; each Rosetta process is
single-threaded but memory-hungry (a few GB), so watch RAM when fanning out wide.

---

## NYU torch cluster (SLURM job array)

Cluster parallelism comes from the **scheduler**, not `--parallel`: the array
runs one task per backbone, each doing `--nstruct 5`. 13 tasks × 5 = 65 designs.
Submit script: [run_torch_array.sbatch](run_torch_array.sbatch).

### One-time setup on torch

Everything below is staged under `/scratch/vb1467` (5 TB quota; `$HOME` is only
50 GB). torch has **no Rosetta module**, so the static binary is copied in.

1. Get the repo there (clone or `rsync` `lib/`, including `scripts/`,
   `ncaa_params/`, `perturbed_backbones/`, `config/`). Landing spot:
   `/scratch/vb1467/git-install/lib`.
2. Get Rosetta there. The local binary is **statically linked** (confirmed to run
   on a torch compute-node's glibc), so just copy it (deref the symlink with
   `-L`) and the database — no module needed:
   ```bash
   rsync -aL <local>/.../source/bin/rosetta_scripts.static.linuxgccrelease  torch:/scratch/vb1467/rosetta/
   rsync -a  <local>/.../database                                           torch:/scratch/vb1467/rosetta/
   ```
3. The scratch `run_torch_array.sbatch` already has `--partition=cpu_short` and
   the scratch `ROSETTA_BIN` / `ROSETTA_DB` filled in. torch's CPU partitions are
   `cpu_short` and `cpu_prem`.

> **Slurm account required.** torch now rejects the default `users` account:
> submissions must pass `--account=torch_pr_xxx_yyy`. Your PI registers the
> project at <https://projects.hpc.nyu.edu>. Check what you have with
> `sacctmgr -n -P show assoc user=$USER format=account`, then add
> `--account=<acct>` to every `sbatch` below.

### Submit (from `/scratch/vb1467/git-install/lib` on torch)

```bash
mkdir -p slurm_logs
RUN_DIR=runs/packer_design_torch_$(date +%Y%m%d_%H%M%S); mkdir -p "$RUN_DIR"

python3 rosetta_only_packer_pred/build_manifest.py \
  --backbones-dir perturbed_backbones/perturbed_backbones --out "$RUN_DIR/jobs.tsv"
head -n 14 "$RUN_DIR/jobs.tsv" > "$RUN_DIR/jobs64.tsv"        # header + 13 backbones
N=$(($(wc -l < "$RUN_DIR/jobs64.tsv") - 1))

sbatch --account=torch_pr_xxx_yyy \
  --array=1-$N --export=ALL,RUN_DIR="$RUN_DIR" rosetta_only_packer_pred/run_torch_array.sbatch
```

Smoke-test first with one backbone / one design before the full array:

```bash
head -n 2 "$RUN_DIR/jobs.tsv" > "$RUN_DIR/jobs64.tsv"         # header + 1 backbone
sbatch --account=torch_pr_xxx_yyy --array=1-1 --time=01:00:00 \
  --export=ALL,RUN_DIR="$RUN_DIR" rosetta_only_packer_pred/run_torch_array.sbatch
```

Collect after the array finishes (or as a `--dependency=afterok:<jobid>` job):

```bash
python3 scripts/collect_scores.py --jobs "$RUN_DIR/jobs64.tsv" --run-dir "$RUN_DIR" --out "$RUN_DIR/summary.tsv"
python3 scripts/geometry_check.py --root "$RUN_DIR" --chain A --exclude input_scaffold.pdb \
  --summary-out "$RUN_DIR/geometry_summary.tsv" --issues-out "$RUN_DIR/geometry_issues.tsv"
```

### Notes

- Each task runs its 5 `nstruct` sequentially (~minutes–~1h each); all tasks run
  concurrently, so wall time ≈ one task. For lower wall time, use more tasks with
  smaller `--nstruct`.
- The hard gates discard non-conforming designs, so *attempts* ≥ *emitted PDBs*.
  For margin, use all 17 backbones (`head -n 18`, `--array=1-17` → 85 attempts).
- Outputs land in `$RUN_DIR/<backbone_stem>/`; SLURM logs in `slurm_logs/`.
- SLURM on torch requires a valid project `--account` (see setup note above); the
  default `users` account is rejected. Single-threaded CPU work → `cpu_short`.

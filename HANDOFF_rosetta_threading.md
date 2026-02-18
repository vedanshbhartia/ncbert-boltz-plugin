# Handoff: 18-mer Chain-A Threading + Rosetta Batch Refinement

## Date / Environment
- Date of this run: **2026-02-12**
- Repo: `/home/vedansh/projects/random/git-install/lib`
- Rosetta binary used:
  - `~/Downloads/rosetta_binary_ubuntu_3.15_bundle/rosetta.binary.ubuntu.release-408/main/source/bin/rosetta_scripts.static.linuxgccrelease`
- Rosetta database used:
  - `~/Downloads/rosetta_binary_ubuntu_3.15_bundle/rosetta.binary.ubuntu.release-408/main/database`

## Goal Completed in This Chat
- Trim chain A from 20 to 18 residues by deleting residues 4 and 5.
- Thread sequences from:
  - `SA474935_cand1_backbone_design.txt`
  - `SA2989311_cand2_backbone_design.txt`
- Run Rosetta relax workflow for each sequence (`nstruct=1`, parallel=8).
- Continue on per-sequence failure and collect scores/status in one table.

## Important Files Added / Modified

### Added scripts
1. `scripts/prepare_trimmed_input.py`
   - Parses mmCIF atom_site loop.
   - Deletes specified chain residues (here `A:4,5`).
   - Renumbers trimmed chain to 1..N.
   - Writes PDB for Rosetta input.
   - Has optional sequence assertions (`--expect-old-seq`, `--expect-new-seq`).

2. `scripts/build_job_manifest.py`
   - Reads sequence text files.
   - Validates effective peptide length (treats bracketed tokens as single residues).
   - Writes `jobs.tsv` columns:
     - `job_id`, `source_file`, `line_no`, `pepseq`

3. `scripts/run_rosetta_batch.sh`
   - Runs one Rosetta job per manifest row.
   - Uses:
     - `-parser:protocol relax_script_thread_pep.xml`
     - `-parser:script_vars pepseq=...`
     - `-load_PDB_components`
     - `-nstruct` (default 1)
   - Parallel execution with background jobs and `wait -n`.
   - Writes per-job:
     - `rosetta.log`
     - `status.env` (exit code, timestamps, error line)
     - `score.sc` (native Rosetta scorefile if present)
   - Fallback behavior:
     - If Rosetta exits 0 but no `score.sc`, script extracts final pose total from output PDB `pose` line and writes synthetic `score.sc`.

4. `scripts/collect_scores.py`
   - Aggregates run outputs into summary TSV:
     - `job_id, source_file, line_no, pepseq, status, total_score, model_path, scorefile_path, log_path, error`
   - Reads score from `score.sc` if present.
   - Fallback: parse `pose` total from PDB energies table.
   - Marks success when `exit_code==0` and model PDB exists.

### Modified protocol
5. `relax_script_thread_pep.xml`
   - Changed:
     - `pep_last` selector from negative index slice to explicit 18:
       - from `from="-1" to="-1"`
       - to `from="18" to="18"`
   - Why: in this Rosetta build, `Slice ... from="-1"` returned zero residues, causing:
     - `DeclareBond expected a residue selector to select exactly one residues. Instead, it selected 0`

### Modified ignore rules
6. `.gitignore`
   - Added:
     - `runs/`
     - `inputs/*.pdb`
     - `ROSETTA_CRASH.log`

## Generated Outputs
- Trimmed input PDB:
  - `inputs/alpha_001_chainA_del4_5_18mer.pdb`
- Full manifest:
  - `runs/jobs.tsv` (20 jobs)
- Full run directory:
  - `runs/full_run/`
- Consolidated summary:
  - `runs/full_run/summary.tsv`

## Validation Performed
1. Trim validation
   - Old chain A sequence: `DHSSSGGGTFTDSSTGNIAY` (20)
   - New chain A sequence: `DHSGGGTFTDSSTGNIAY` (18)
   - Chain B residue count preserved (158).

2. Pilot tests
   - Canonical sequence succeeded.
   - Bracket-containing sequence failed cleanly, logged parse error.

3. Full batch integrity
   - 20 rows in summary.
   - Every row has status.

## Full Batch Result (Current State)
- Total: 20
- Success: 4
- Failed: 16

Successful jobs (sorted by total_score ascending):
1. `SA2989311_cand2_L10` -> `-580.892000`
2. `SA2989311_cand2_L08` -> `-574.923000`
3. `SA474935_cand1_L06` -> `-559.398000`
4. `SA2989311_cand2_L07` -> `-449.332000`

All 16 failures share the same root cause:
- `SimpleThreadingMover::determine_mutations_oneletter(): Could not parse character "[" in sequence ...`

## Root Cause of Noncanonical Failure
- Current XML uses:
  - `SimpleThreadingMover ... sequence_mode="oneletter" thread_sequence="%%pepseq%%"`
- Input sequence files encode noncanonicals with bracket tokens (e.g., `[MLY]`, `[DTH]`, `[CME]`).
- `SimpleThreadingMover` in current mode/build does **not** accept `[` bracket token syntax in this field.

## Repro Commands

### 1) Build manifest
```bash
python3 scripts/build_job_manifest.py \
  --inputs SA474935_cand1_backbone_design.txt SA2989311_cand2_backbone_design.txt \
  --out runs/jobs.tsv \
  --expected-length 18
```

### 2) Build trimmed input
```bash
python3 scripts/prepare_trimmed_input.py \
  --in-cif alpha_001_1a85_cyclic_199.cif \
  --out-pdb inputs/alpha_001_chainA_del4_5_18mer.pdb \
  --chain A \
  --delete 4,5 \
  --expect-old-seq DHSSSGGGTFTDSSTGNIAY \
  --expect-new-seq DHSGGGTFTDSSTGNIAY
```

### 3) Run batch
```bash
scripts/run_rosetta_batch.sh \
  --jobs runs/jobs.tsv \
  --input-pdb inputs/alpha_001_chainA_del4_5_18mer.pdb \
  --xml relax_script_thread_pep.xml \
  --run-dir runs/full_run \
  --parallel 8 \
  --nstruct 1
```

### 4) Aggregate summary
```bash
python3 scripts/collect_scores.py \
  --jobs runs/jobs.tsv \
  --run-dir runs/full_run \
  --out runs/full_run/summary.tsv
```

## Recommended Next Steps
1. Add a **noncanonical threading strategy** that is compatible with Rosetta (current blocker).
2. Keep canonical sequence path unchanged (already works).
3. For noncanonical rows, implement a preprocessing + mutation path, e.g.:
   - Parse bracket tokens into a tokenized representation.
   - Thread a compatible base sequence first.
   - Apply post-thread residue mutations using Rosetta-supported residue names/types for modified residues.
4. Re-run only the 16 failed jobs first (incremental retry mode), then regenerate `summary.tsv`.

## Known Caveats
- `relax_script_thread_pep.xml` currently assumes trimmed 18-mer chain A (explicit selector `from=18 to=18`).
- `runs/` and trimmed `inputs/*.pdb` are ignored by git (by design).
- `score.sc` may be synthetic for successful jobs when Rosetta doesn't emit one; in that case score is taken from the output PDB `pose` energy table.

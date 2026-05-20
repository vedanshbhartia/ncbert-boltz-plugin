# Rosetta-based design evaluation workflow

This documents the pipeline driven by [run.sh](run.sh), which threads MMP8-targeted
peptide designs onto cyclic backbones, relaxes them with RosettaScripts, and
collects scores plus geometry diagnostics. It also captures the recent change
that made the threading step pick the correct perturbed backbone per design.

## Entry point: `run.sh`

The script orchestrates four stages:

1. **Reference backbone run** — runs the relax/scoring protocol on a fixed
   reference scaffold without threading, to obtain a baseline `dG_separated`.
2. **Manifest construction** — selects design candidates from
   [design/](design/) and emits a Rosetta thread-sequence manifest.
3. **Design batch** — runs the threading + relax protocol for each selected
   design.
4. **Score collection + geometry check** — aggregates scorefiles into a
   summary TSV and runs bond-geometry diagnostics on output PDBs.

Tunable env vars (with defaults):

| Variable        | Default                                            | Notes                                              |
|-----------------|----------------------------------------------------|----------------------------------------------------|
| `LIMIT`         | `6`                                                | Number of design candidates to select.             |
| `PARALLEL`      | `3`                                                | Concurrent Rosetta jobs.                           |
| `NSTRUCT`       | `1`                                                | `nstruct` per Rosetta invocation.                  |
| `RUN_TAG`       | `sa474935_<timestamp>`                             | Run tag used in `RUN_DIR`.                         |
| `RUN_DIR`       | `runs/$RUN_TAG`                                    | Output directory.                                  |
| `ROSETTA_CONFIG`| `config/rosetta.env`                               | Path to Rosetta defaults file.                     |
| `SCAFFOLD_PDB`  | `SA474935_cand1_backbone.pdb`                      | Used by the **reference** run only (see caveats). |

## Inputs

### Design candidates: `design/*.jsonl`

Files are named like `alpha_01_1a85_cyclic_001_PerturbN_design.jsonl` and
contain a `samples_b` array. Each sample provides:

- `aa3_tokens` — three-letter codes including NCAA tokens such as `HYP`, `AIB`,
  `DPRO`, `DHYP`.
- `aa1_string` — best-effort one-letter rendering.
- `score` — model score used for ranking inside a category.
- `ncaa_list` — NCAA residues present in the sample.

Selection logic in [scripts/build_1a85_target_manifest.py](scripts/build_1a85_target_manifest.py):

- A sample is eligible only if all three-letter tokens are in
  `{20 canonical AAs, HYP, AIB, DPRO/DPR, DHYP}`.
- It is then categorized as one of `HYP`, `AIB`, `DPRO_LIKE`, `PRO` (first match
  wins). Samples without any of those go into none of the buckets and are
  skipped.
- Selection is round-robin across `(HYP, AIB, DPRO_LIKE, PRO)` in order, until
  `--limit` candidates are picked.

### Perturbed backbones: `perturbed_backbones/perturbed_backbones/`

PDB files named `alpha_01_1a85_cyclic_001_PerturbN.pdb` for `N = 1..17`.
Layout:

- **Chain A**: 18-residue cyclic peptide (residues 1..18).
- **Chain B**: MMP8 target (residues 1..158).

This is the **opposite** chain assignment from the legacy `SA474935_cand1_backbone.pdb`,
where the peptide is on chain B and MMP8 on chain A.

### Threaded thread sequence format

`build_1a85_target_manifest.py` converts `aa3_tokens` to Rosetta one-letter
threading syntax:

| Token        | Threaded form |
|--------------|---------------|
| Canonical AA | one-letter    |
| `HYP`        | `X[HPR]`      |
| `AIB`        | `X[AIB]`      |
| `DPRO`/`DPR` | `X[DPRO]`     |
| `DHYP`       | `X[DHYP]`     |

[scripts/normalize_thread_sequence.py](scripts/normalize_thread_sequence.py) is
the more general normalizer used by the non-threadseq batch script — it maps
extra aliases (e.g. `NLE`→`NLU`, `CRO`/`GYS`/`CR2`/`PIA`→`CRO`, etc.) and
falls back to CCD one-letter codes via the Rosetta DB. The threadseq path
skips it because the manifest is already pre-normalized.

## Pipeline detail

### 1. Reference run

```
scripts/run_rosetta_batch.sh \
  --jobs reference/jobs.tsv \
  --input-pdb $SCAFFOLD_PDB \
  --xml relax_script_reference_pep.xml \
  --peptide-chain B \
  --include-hetatm-backbone-check
```

`relax_script_reference_pep.xml` ([relax_script_reference_pep.xml](relax_script_reference_pep.xml))
expects peptide=B, MMP8=A — it does **not** thread; it just relaxes the
pre-existing peptide. The output PDB is parsed by
[scripts/extract_dg_separated.py](scripts/extract_dg_separated.py) to extract
`dG_separated`, which is fed to `collect_scores.py` as `--reference-dg` so the
summary gets a per-design `ddg`.

### 2. Manifest construction

```
python3 scripts/build_1a85_target_manifest.py \
  --design-dir design \
  --backbones-dir perturbed_backbones/perturbed_backbones \
  --out $RUN_DIR/jobs.tsv \
  --selected-out $RUN_DIR/selected.tsv \
  --limit $LIMIT
```

Outputs:

- `jobs.tsv` (5 columns): `job_id`, `source_file`, `line_no`, `pepseq`, `input_pdb`.
- `selected.tsv` (richer): adds `category`, `score`, `aa1_string`,
  `aa3_sequence`, `ncaa_list`, `input_pdb`.

The new `input_pdb` column is derived from `source_file`:
strip the trailing `_design` from the stem, append `.pdb`, look up in
`--backbones-dir`. If the file is missing, the script aborts.

#### Explicit list builder: `scripts/build_explicit_manifest.py`

For ad-hoc runs of a hand-curated set of designs that may not exist in the
`design/*.jsonl` files (or use different sample indices), use
[scripts/build_explicit_manifest.py](scripts/build_explicit_manifest.py).

Input format — whitespace-separated, one row per design:

```
<job_id>[_<nstruct_idx>]   <score>   <one_letter_pepseq>
```

Example:

```
alpha_01_1a85_cyclic_001_Perturb5_design_S0167_0001 -24.526 ALNIAHEDVVpDTPAYAL
alpha_01_1a85_cyclic_001_Perturb12_design_S0007_0001 -15.606 lIEHSEDIITDSNPESLT
```

Behavior:

- Strips any trailing `_NNNN` nstruct suffix from the job id (so `_0001` is
  removed). The canonical id matches the design jsonl stem convention used by
  `build_1a85_target_manifest.py`.
- Resolves the matching `_PerturbN.pdb` from `--backbones-dir` by extracting
  the `Perturb<digits>` token from the job id.
- Converts **lowercase letters to D-amino acids** in Rosetta `X[DXXX]` bracket
  form. Mapping (matches the residue types in the design palette):

  | char | type   | char | type   | char | type   | char | type   |
  |------|--------|------|--------|------|--------|------|--------|
  | a    | DALA   | f    | DPHE   | k    | DLYS   | r    | DARG   |
  | d    | DASP   | h    | DHIS   | l    | DLEU   | s    | DSER   |
  | e    | DGLU   | i    | DILE   | n    | DASN   | t    | DTHR   |
  |      |        |      |        | p    | DPRO   | v    | DVAL   |
  |      |        |      |        | q    | DGLN   | w    | DTRP   |
  |      |        |      |        |      |        | y    | DTYR   |

  Uppercase letters pass through unchanged. Any character outside this set
  raises an error.
- Validates the threaded length against `--expected-length` (default 18,
  counting each bracket token as one residue).
- Emits the same 5-column `jobs.tsv` the threadseq batch script consumes, and
  optionally a richer `selected.tsv` that also keeps the original raw
  sequence and the user-supplied score for traceability.

Typical usage (skipping the round-robin selector and the reference run):

```
RUN_DIR=runs/explicit_$(date +%Y%m%d_%H%M%S)
mkdir -p "$RUN_DIR"
# put your <job_id> <score> <pepseq> rows in $RUN_DIR/sequences.txt

python3 scripts/build_explicit_manifest.py \
  --input "$RUN_DIR/sequences.txt" \
  --backbones-dir perturbed_backbones/perturbed_backbones \
  --out "$RUN_DIR/jobs.tsv" \
  --selected-out "$RUN_DIR/selected.tsv"

scripts/run_rosetta_batch_threadseq.sh \
  --config config/rosetta.env \
  --jobs "$RUN_DIR/jobs.tsv" \
  --xml relax_script_thread_1a85_chainA_staged_rigid_jump.xml \
  --run-dir "$RUN_DIR" \
  --peptide-chain A \
  --include-hetatm-backbone-check \
  --parallel 4 --nstruct 1

python3 scripts/collect_scores.py \
  --jobs "$RUN_DIR/jobs.tsv" \
  --run-dir "$RUN_DIR" \
  --out "$RUN_DIR/summary.tsv"

python3 scripts/geometry_check.py \
  --root "$RUN_DIR" \
  --chain A \
  --exclude input_scaffold.pdb --exclude reference \
  --summary-out "$RUN_DIR/geometry_summary.tsv" \
  --issues-out "$RUN_DIR/geometry_issues.tsv"
```

Without a reference run, `summary.tsv` will not have a `ddg` column.

### 3. Design batch (threaded sequences)

```
scripts/run_rosetta_batch_threadseq.sh \
  --jobs $RUN_DIR/jobs.tsv \
  --xml relax_script_thread_1a85_chainA_staged_rigid_jump.xml \
  --peptide-chain A \
  --include-hetatm-backbone-check \
  --parallel $PARALLEL --nstruct $NSTRUCT
```

Per-job behavior in [scripts/run_rosetta_batch_threadseq.sh](scripts/run_rosetta_batch_threadseq.sh):

- Reads the optional 5th column `input_pdb`. If present, that scaffold is used
  for that job; otherwise the global `--input-pdb` (now optional) is used as
  fallback. If neither resolves, the script aborts before running anything.
- Backbone-continuity preflight runs once per **unique** scaffold collected
  from the manifest, instead of a single check on a global scaffold. The
  preflight calls [scripts/check_backbone_continuity.py](scripts/check_backbone_continuity.py)
  with `--cyclic` and the configured peptide chain; failures abort the batch.
- Each job snapshots its scaffold to `<job_dir>/input_scaffold.pdb` for
  provenance, alongside `rosetta.log`, `score.sc`, `status.env`, and any model
  PDBs.
- `status.env` records `input_pdb` so downstream consumers can recover which
  scaffold was used.
- `run_config.txt` lists every unique `job_input_pdb=` and the count.

The other batch script [scripts/run_rosetta_batch.sh](scripts/run_rosetta_batch.sh)
is the non-threadseq variant used by the reference run. It still expects a
single `--input-pdb` and runs the bracket-token normalizer on each `pepseq`.

#### Rosetta invocation specifics

Both batch scripts call `rosetta_scripts.static.linuxgccrelease` with:

- `-database $ROSETTA_DB`
- `-s <input scaffold>`
- `-parser:protocol <xml>`
- `-parser:script_vars pepseq=<thread sequence>`
- `-load_PDB_components`
- `-in:file:extra_res_fa <ncaa params...>` — params files in
  [ncaa_params/](ncaa_params/) are auto-discovered, but ones whose `NAME`
  collides with a Rosetta `fa_standard` residue type are skipped (those names
  are recorded as `extra_res_fa_skipped=` in `run_config.txt`).
- `-nstruct $NSTRUCT -overwrite -out:path:all <job dir> -out:file:scorefile <job_dir>/score.sc`

If Rosetta returns success but writes no scorefile (which can happen on some
builds), the script reconstructs a minimal scorefile (just `total_score
description`) by parsing the `pose` total score line out of the output PDB.
This fallback is in fact what triggers in our environment — the per-job
`score.sc` typically has only those two columns. The richer per-term metrics
(`dG_separated`, `shape`, `dSASA_int`, etc.) are still recovered, because
`collect_scores.py` also parses the `BEGIN_POSE_ENERGIES_TABLE` and the
trailing key/value lines that `InterfaceAnalyzerMover` writes inside the
output PDB itself.

If `rosetta.log` contains "no jobs were attempted", the exit code is bumped to
`2` so the failure surfaces in `status.env`.

### 4. Score collection + geometry

```
python3 scripts/collect_scores.py \
  --jobs $RUN_DIR/jobs.tsv \
  --run-dir $RUN_DIR \
  --out $RUN_DIR/summary.tsv \
  [--reference-dg $REFERENCE_DG]
```

[scripts/collect_scores.py](scripts/collect_scores.py) walks each `<job>/score.sc`,
picks the best (lowest `total_score`) row, and joins it with the job manifest
and `status.env`. With `--reference-dg`, it computes `ddg = dG_separated - reference_dg`.

```
python3 scripts/geometry_check.py \
  --root $RUN_DIR \
  --chain A \
  --exclude input_scaffold.pdb \
  --summary-out $RUN_DIR/geometry_summary.tsv \
  --issues-out $RUN_DIR/geometry_issues.tsv
```

[scripts/geometry_check.py](scripts/geometry_check.py) replicates the bond-geometry
audits from the matlab reference: bond lengths vs ideal `(1.458, 1.524, 1.329)`,
bond angles vs ideal `(121.7°, 111.2°, 116.2°)`, and omega torsion deviations,
with thresholds `0.2 Å / 10° / 20°`. It runs over every PDB it finds under
`--root`, excluding `input_scaffold.pdb` snapshots so only generated models are
audited. `--chain A` matches the new perturbed-backbone layout.

`--exclude` matches **any path component**, not just basenames. That lets you
skip whole subtrees by name. In `run.sh` we pass both
`--exclude input_scaffold.pdb` (the per-job scaffold snapshots) and
`--exclude reference` (the entire reference subtree, whose output PDB has the
peptide on chain B, so scanning chain A there would treat MMP8 as the
peptide and trip on residue 158 missing backbone atoms).

## RosettaScripts protocols

Three XML files, all referencing `ref2015` with re-weighted backbone hbonds,
chainbreak, omega, and a small `voids_penalty`:

- [relax_script_reference_pep.xml](relax_script_reference_pep.xml) — peptide=B,
  MMP8=A. No threading. Used by the reference run.
- [relax_script_thread_pep_staged_rigid_jump.xml](relax_script_thread_pep_staged_rigid_jump.xml) —
  legacy: peptide=B, MMP8=A, `start_position="162B"`, `interface="A_B"`.
  Designed for `SA474935_cand1_backbone.pdb`. **No longer used by `run.sh`**
  after the per-perturbation switch, but kept for reference / older scaffolds.
- [relax_script_thread_1a85_chainA_staged_rigid_jump.xml](relax_script_thread_1a85_chainA_staged_rigid_jump.xml) —
  current: peptide=A, MMP8=B (named `select_target` in this file),
  `start_position="1A"`, `interface="B_A"`. This is the file `run.sh` now uses
  for the design step.

Common protocol structure (threading variants):

1. `SimpleThreadingMover` writes `pepseq` onto the peptide chain starting at
   the configured residue, with `pack_neighbors=true`, neighbor distance 8 Å.
2. Termini variants are stripped from the first/last peptide residues, a
   `DeclareBond` closes head-to-tail, and `CUTPOINT_UPPER`/`CUTPOINT_LOWER`
   variants are added so cyclization is treated as a real bond.
3. Near-peptide MMP8 residues are stored as a static `subset_near_peptide` so
   the relax does not see the selection drift after threading.
4. Stage 1 `FastRelax` (`frlx_fixed`) keeps backbone fixed, allows jumps and
   sidechain repacking on the peptide and the stored MMP8 neighborhood.
5. `OversaturatedHbondAcceptorFilter` rejects designs with >0 backbone hbond
   acceptors with more than 2 donors.
6. Stage 2 `FastRelax` (`frlx_move`) reintroduces peptide backbone movement.
7. `ShapeComplementarity` and `InterfaceAnalyzerMover` produce the final
   metrics. `pack_separated=true` is what gives `dG_separated`.

The custom `CustomBaseTypePackerPalette` lists the D-amino-acid and NCAA
types the design palette is allowed to use (e.g. `HYP`, `AIB`, `ORN`, `NLU`,
`HPR`, etc.) — relevant if any future protocol does design rather than just
threading + relax.

## Per-job backbone selection (recent change)

The original pipeline used a single global scaffold for every design. With
many perturbed cyclic backbones available, each design must be relaxed on the
backbone it was generated against, otherwise the threaded peptide is being
forced onto unrelated geometry.

Mapping:

```
design/alpha_01_1a85_cyclic_001_Perturb12_design.jsonl
  -> perturbed_backbones/perturbed_backbones/alpha_01_1a85_cyclic_001_Perturb12.pdb
```

Implementation summary:

- [scripts/build_1a85_target_manifest.py](scripts/build_1a85_target_manifest.py)
  resolves the matching backbone for every selected sample, emits an
  `input_pdb` column, and aborts if any expected backbone PDB is missing.
- [scripts/run_rosetta_batch_threadseq.sh](scripts/run_rosetta_batch_threadseq.sh)
  reads the new column, makes `--input-pdb` optional (used only as fallback),
  validates each per-job scaffold up front, runs the backbone-continuity
  preflight per **unique** scaffold, snapshots per-job inputs, and records
  `input_pdb` in `status.env` and `run_config.txt`.
- [run.sh](run.sh) switched the design step to the chain-A XML, set
  `--peptide-chain A`, dropped its global `--input-pdb` for that step, and
  changed the geometry check to `--chain A`.

## Config

[config/rosetta.env](config/rosetta.env) sets the defaults shared by both
batch scripts:

```
ROSETTA_ROOT=$HOME/Downloads/rosetta_binary_ubuntu_3.15_bundle/rosetta.binary.ubuntu.release-408/main
ROSETTA_BIN=$ROSETTA_ROOT/source/bin/rosetta_scripts.static.linuxgccrelease
LOOPMODEL_BIN=$ROSETTA_ROOT/source/bin/loopmodel.static.linuxgccrelease
ROSETTA_DB=$ROSETTA_ROOT/database
EXTRA_RES_FA_DIR=ncaa_params
```

CLI flags (`--rosetta-bin`, `--rosetta-db`, `--extra-res-fa-dir`) override the
config values when passed.

## Reading `summary.tsv`

Useful columns:

| column                        | meaning                                                                                  |
|-------------------------------|------------------------------------------------------------------------------------------|
| `total_score`                 | Rosetta total score (lower is better).                                                   |
| `dg_separated` (`score_dG_separated`) | Interface binding score from `InterfaceAnalyzerMover` with `pack_separated=true`. |
| `score_shape`                 | Shape complementarity from the `<ShapeComplementarity>` filter (0..1, higher is better). |
| `score_shape_int_area`        | Interface area (Å²) reported by the same filter when `write_int_area=1`.                 |
| `score_dSASA_int`             | Total interface ΔSASA from `InterfaceAnalyzerMover`.                                     |
| `score_hbonds_int`            | Interface hbond count.                                                                   |
| `score_oversat`               | `OversaturatedHbondAcceptorFilter` score; 0 means the filter passed.                     |

**Watch out** — there are two columns with similar names because Rosetta
writes the same metric from two different machines:

- `score_shape` is the `<ShapeComplementarity>` filter's output. **Use this.**
- `score_sc_value` comes from `InterfaceAnalyzerMover` and is `0` for every
  job in our outputs (the IAM doesn't compute SC under our setup). It is
  *not* an SC failure; just ignore it.

## Viewing output PDBs in PyMOL

The output models look "missing sidechains" in default PyMOL views — they are
not. Reasons:

- **Cartoon hides sidechains.** Default `cartoon` only draws the backbone
  ribbon. Show sidechains with `show sticks, chain A` (or restrict to heavy
  atoms: `show sticks, chain A and not name N+CA+C+O+H*`).
- **Non-canonicals are HETATM records.** `HPR`, `DPRO`, `DTYR`, `DGLN`,
  `DALA`, `DLEU`, `DASP`, `AIB`, etc. are written as `HETATM`. PyMOL treats
  HETATM as ligand by default, so the cartoon ribbon "skips" those positions
  and only the NCAA's sidechain pops out as sticks. After
  `show sticks, chain A` it looks consistent.
- **The cyclic peptide bond is not drawn automatically.** The head-to-tail
  amide between residue 1's `N` and residue 18's `C` is real (declared via
  `DeclareBond` and reflected in the energy terms) but PyMOL won't infer the
  bond across an ATOM/HETATM gap. Force it with:
  ```
  bond chain A and resi 1 and name N, chain A and resi 18 and name C
  ```
- **Some MMP8 residues genuinely have no sidechain past CB** because they
  are GLY (no CB) or ALA (only CB). About 13/158 in these scaffolds. That is
  chemistry, not data loss.
- **The peptide ring really is closed.** The output PDB writes an explicit
  `LINK N <res1> ... C <res18>` record at distance ~1.35 Å, and measuring
  N(res1)→C(res18) across every model in the explicit batch gives 1.32–1.36 Å
  (canonical peptide bond). If PyMOL shows a "gap" it is purely a rendering
  default, not a broken structure.
- **Extended K/R sidechains can look like an out-of-place "ladder."** The
  movemap factories only allow chi minimization on the peptide and on MMP8
  residues within 8 Å of the peptide, so surface-exposed Lys/Arg residues on
  the peptide that don't pack into a binding-site pocket end up in extended
  conformations sticking into solvent. With explicit hydrogens shown as
  sticks, that looks like a long ladder with an NH₃⁺/guanidinium cap. That
  is expected for the stage-2 relax with `disable_design="true"` — the
  protocol doesn't search alternate rotamers for these solvent-exposed
  residues. To declutter:
  ```
  hide sticks, chain A and resn LYS+ARG and not (name N+CA+C+O+CB)
  ```

A reasonable PyMOL visualization recipe for one of these models:

```
load alpha_..._Perturb14_0001.pdb, design
hide everything
show cartoon, design
color grey80, design and chain B
show sticks, design and chain A
color cyan, design and chain A and elem C
hide everything, hydro
bond design and chain A and resi 1 and name N, design and chain A and resi 18 and name C
zoom design and chain A, 6
```

## Boltz-2 evaluation of Rosetta outputs

[boltz/run_boltz.py](boltz/run_boltz.py) feeds the Rosetta-relaxed PDBs into
Boltz-2 to get an independent structural/binding-confidence read. It is CPU
only by default (the existing setup runs `--accelerator cpu`).

### Source layout assumption

Rosetta outputs in `runs/explicit_<tag>/<design>/<perturb>_0001.pdb` use
**chain A = peptide, chain B = MMP8**. The Boltz YAML keeps the project
convention of **chain A = MMP8, chain B = peptide**, so the script remaps
template chain ids: YAML chain A ← source chain B, YAML chain B ← source
chain A.

### D-amino acids and NCAAs

Rosetta writes D-residues and NCAAs as `HETATM` records with 3-letter codes
like `DAL, DAS, DGN, DLE, DPR, DSN, DTH, DTR, DTY, DVA, HPR, AIB`. The
script detects modifications **by residue name**, not by phi-angle geometry
— the phi heuristic gave many false positives on relaxed cyclic peptides
because lots of L residues legitimately have phi > 30° in cyclic
conformations. The mapping is in
`ROSETTA_NCAA_INFO` (Rosetta 3-letter → canonical L 1-letter +
CCD code). For each modified residue, the YAML carries the canonical L
1-letter in the `sequence` and a `modifications: [{position, ccd}]` entry
that tells Boltz the actual chemistry.

The template CIF is built by `make_clean_template_cif`:

- Keeps every polymer residue (canonical L + Rosetta NCAAs) so the residue
  count matches the YAML polymer.
- Renames NCAA residues to their canonical L 3-letter (`DLE` → `LEU`, `HPR`
  → `PRO`, etc.) so the template residue sequence aligns with the YAML
  sequence boltz constructs from the modifications list.
- Drops single-residue HETATM ions/cofactors (metals, buffer molecules) that
  cannot be polymer template residues.
- Drops `LINK` / connection records (the cyclic bond is conveyed via
  `cyclic: true` on the peptide entity in YAML; LINK records in templates
  trip Boltz's parser).

Two `read_pdb_atoms` gotchas worth recording:

1. The reader **must** parse both `ATOM` and `HETATM` records. The original
   `ATOM`-only version silently dropped every D-residue, producing a
   17-residue peptide and zero modifications.
2. Position indices in modifications are **1-indexed within the chain**, not
   PDB residue numbers. The script enumerates `sorted(atoms.keys())` so it
   does the right thing as long as the chain is contiguous.

### CLI

```
python3 boltz/run_boltz.py \
  [--rosetta-run-dir runs/explicit_<tag>] \
  [--out-root boltz/runs] \
  [--designs DESIGN_ID [DESIGN_ID ...]] \
  [--n-samples 3] \
  [--skip-existing]
```

Defaults run on every design subdirectory under
`runs/explicit_20260506_014900` with 3 diffusion samples per design.
`--skip-existing` lets you resume after an interrupted batch — designs
whose `boltz_out/.../predictions/input/` directory already exists are
skipped.

### What Boltz reports

Per-design output goes to `boltz/runs/<design>/boltz_out/boltz_results_input/predictions/input/`:

- `input_model_<k>.pdb` — k-th sampled structure, ordered by confidence.
- `confidence_input_model_<k>.json` — per-sample confidence scores.

`run_boltz.py` prints a summary per design with `confidence_score`,
`complex_plddt`, `ptm`, `iptm`, `protein_iptm`, `complex_iplddt` and
interprets them at coarse thresholds. Boltz **affinity is unavailable for
these designs** because the affinity head requires a small-molecule ligand
(<128 heavy atoms) and an 18-residue cyclic peptide exceeds that. Use
`protein_iptm` and `confidence_score` as binding-quality proxies.

### Working around colabfold MSA-server outages

`--use_msa_server` queries the colabfold mmseqs2 endpoint at
`https://api.colabfold.com`. That server intermittently returns errors like
`Server didn't reply with json: mkdir ...: no space left on device`. When
that happens, Boltz **prints "Failed to process … Skipping" and still exits
0** — without producing a `predictions/input/` directory. The hardened
`run_boltz` checks for that directory after the subprocess returns and
raises `RuntimeError` (printed as `[fail] <design>: ...`) so the failure is
not silent.

Workaround when the server is misbehaving: reuse the MSA from a previously
successful design. The MMP8 sequence is identical across all 18 designs, so
its `input_0.csv` is interchangeable. The peptide MSA is essentially
"self-only" anyway (cyclic NCAA peptides have no useful homologs), so a
2-row CSV with the design's own sequence is fine.

```bash
# pick any successful design's MMP8 MSA (~15k sequences for MMP8)
SRC=boltz/runs/<good_design>/boltz_out/boltz_results_input/msa/input_0.csv

# for each failed design, stage cached MSAs and patch the YAML to point at
# them, then run boltz WITHOUT --use_msa_server:
for design in <failed_designs>; do
  DST=boltz/runs/$design
  rm -rf "$DST/boltz_out"
  cp "$SRC" "$DST/cached_mmp8.csv"
  pep=$(python3 -c "import yaml; print(yaml.safe_load(open('$DST/input.yaml'))['sequences'][1]['protein']['sequence'])")
  printf 'key,sequence\n0,%s\n-1,%s\n' "$pep" "$pep" > "$DST/cached_pep.csv"
  python3 -c "
import yaml; p='$DST/input.yaml'
d=yaml.safe_load(open(p))
d['sequences'][0]['protein']['msa']='$DST/cached_mmp8.csv'
d['sequences'][1]['protein']['msa']='$DST/cached_pep.csv'
yaml.dump(d, open(p,'w'), default_flow_style=False, sort_keys=False)
"
  venv/bin/boltz predict "$DST/input.yaml" --out_dir "$DST/boltz_out" \
    --use_potentials --diffusion_samples 3 --recycling_steps 3 \
    --output_format pdb --accelerator cpu --override
done
```

The `msa:` field in the YAML overrides the server lookup. Multi-chain inputs
require CSV format (rows with the same `key` are mutually paired across
chains; `key=-1, -2, ...` rows are unpaired entries).

## Verification runs to date

### Smoke run on a single design

`runs/verify_one_220249/` — `LIMIT=1 PARALLEL=1 NSTRUCT=1`. One HYP-bucket
design (`Perturb11_S0056`, sequence `AYEVSX[HPR]EVVVDDTPAKYI`). Verified that
the per-job `input_scaffold.pdb` md5 equals the source `Perturb11.pdb`,
chain A residue 6 is the threaded `HPR` HETATM, and the geometry check passes
(after the `--exclude reference` change described above).

### Explicit 18-design batch

`runs/explicit_20260506_014900/` — built from
[scripts/build_explicit_manifest.py](scripts/build_explicit_manifest.py) on a
hand-curated list with both canonical and D-amino-acid sequences. All 18 jobs
succeeded; geometry 17/18 pass. The single failure (`Perturb17_S0543`,
sequence `pIQEAAYLFYNNREVYFI` with D-Pro at position 1) had a 23.3° omega
deviation at residue 18 — i.e. at the head-to-tail cyclization seam right
next to the D-Pro head, which is the most strained position. Top hits by
`total_score`:

| rank | job                | total   | dG_sep  | shape | int_area |
|------|--------------------|---------|---------|-------|----------|
| 1    | Perturb14_S0662    | −664.07 | −76.24  | 0.68  | 1027     |
| 2    | Perturb5_S0167     | −650.47 | −79.99  | 0.67  | 1088     |
| 3    | Perturb11_S0231    | −644.34 | −77.68  | 0.59  | 1136     |
| 4    | Perturb15_S0368    | −641.49 | −68.09  | 0.60  | 1113     |
| 5    | Perturb17_S0304    | −634.27 | −70.78  | 0.69  | 1059     |

(Full table in `summary.tsv`. No reference run, so no `ddg` column.)

## Caveats and rough edges

- **Reference vs design scaffolds are heterogeneous.** The reference run uses
  `SA474935_cand1_backbone.pdb` (peptide=B, MMP8=A), while design jobs now use
  per-perturbation `1a85`-derived backbones (peptide=A, MMP8=B). The single
  `dG_separated` produced by the reference is fed into `collect_scores.py` as
  the baseline for `ddg`, but the comparison is no longer apples-to-apples.
  A more principled fix is a per-perturbation reference (relax each
  perturbed backbone with its native sequence and use that backbone's
  `dG_separated` as the baseline for designs derived from the same
  perturbation), but that is not implemented.
- **`run_rosetta_batch.sh` was not updated** — it still expects a single
  `--input-pdb`. Only the threadseq variant supports per-job scaffolds. If the
  reference workflow ever needs per-job scaffolds, the same change must be
  ported.
- **Backbone preflight chain.** The preflight uses `--peptide-chain` for the
  cyclic-closure check. If a perturbed backbone ever swaps chains, this would
  fail; the manifest assumes chain A is the peptide for every file in
  `perturbed_backbones/perturbed_backbones/`.
- **NCAA params name collisions.** Custom params in [ncaa_params/](ncaa_params/)
  whose `NAME` matches a Rosetta `fa_standard` residue are silently skipped to
  avoid double-registration; check `extra_res_fa_skipped=` in
  `run_config.txt` if a residue type is unexpectedly absent.
- **Jobs TSV header.** `run.sh` writes a header `job_id\tsource_file\t...`;
  both batch scripts skip rows where the first field equals `job_id` (or is
  empty), so manual edits should preserve a header line.

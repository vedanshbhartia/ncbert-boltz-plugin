# Rosetta-packer de novo design (constraint-driven)

This directory holds a **design-from-scratch** alternative to the threading
pipeline documented in [../rosetta_workflow.md](../rosetta_workflow.md).

## Idea

The threading workflow takes a sequence proposed elsewhere, threads it onto a
backbone, repacks, relaxes, and then **filters** designs after the fact for the
properties we want (composition, net charge, etc.). Because sidechains are
repacked by the Rosetta packer afterwards anyway, we can let the **packer do the
sequence design itself**, and instead of selecting good designs out of a pool,
we push the optimizer toward the features we want by encoding them as
**constraints** that are active during design:

- **Amino-acid composition** — at least 1 proline-like residue (PRO/DPRO/HYP/
  DHYP/AIB) and at most 2 Arg, via the `aa_composition` score term and
  [desired_makeup_MMP8.comp](desired_makeup_MMP8.comp).
- **Net charge** — net charge of −1 or lower, via the `netcharge` score term and
  [desired_peptide.charge](desired_peptide.charge).

The `aa_composition` / `netcharge` terms are **soft** — they bias the packer
during design but cannot guarantee the result. The non-negotiable requirements
(proline-like count, an acidic cyclization handle) are therefore *also* enforced
as **hard gates** (see *Filters*), so a design that slips past the soft guidance
is discarded rather than emitted.

All 18 chain-A positions are designed de novo from the input backbone; no
sequence is threaded.

## System

These scaffolds are the `1a85`/MMP8 cyclic-peptide complexes:

- **Chain A** — 18-residue cyclic peptide (the design target), polyGly backbone.
- **Chain B** — MMP8 target (158 residues). No metal.

[design_script.xml](design_script.xml) is adapted from the Mulligan trimmed
cyclic-peptide design demo. The demo used absolute residue numbers
(1–161 = MMP8, 162–179 = peptide) and a catalytic Zn; this version instead:

- expresses the peptide / target split and the head-to-tail cyclization with the
  same `Chain` + `Slice` selectors as `../relax_script_thread_1a85_chainA_*.xml`
  (peptide = chain A, residues 1 and 18 are the cyclization seam);
- drops the `SetupMetalsMover`, the `metalbinding_constraint` reweights, and the
  Zn fold-tree file (`foldtree_MMP8.txt` was deleted — it encoded MMP8 residue
  numbers and a Zn jump that do not exist here; the script uses the default fold
  tree + `DeclareBond`, like the rest of the repo);
- `RamaMutationSelector` is intersected with the peptide chain because this
  Rosetta build's selector has no sub-selector attribute;
- adds an `InterfaceAnalyzerMover` so `dG_separated`/shape are produced for
  [../scripts/collect_scores.py](../scripts/collect_scores.py), comparable with
  the relax scripts.

The protocol runs three rounds of `FastRelax` (full palette repack + minimize)
alternating with `FastDesign` (layer-based sequence design), ramping the
backbone-hbond emphasis down, with CA–CA pocket constraints holding the peptide
in place until the final round.

### Filters

- `oversat` (oversaturated hbond acceptors) — **hard gate**, blocks a
  scorefunction artifact.
- `min_proline_like` (count of PRO/DPRO/AIB/HYP/DHYP in the peptide) — **hard
  gate**, requires at least 1. A `ResidueName` selector
  (`residue_name3="PRO,DPR,AIB,HYP,DHY"`) picks the set; a `ResidueCount` filter
  with `confidence=1` rejects designs below the floor.
- `min_acidic` (count of ASP/GLU, L or D, in the peptide) — **hard gate**,
  requires at least 1 (the cyclization handle). Counts via the `NEGATIVE_CHARGE`
  residue property (which the D-forms inherit), so it needs no D-form name3 codes.
- `internal_hbonds` (peptide internal hbonds) — **report-only** (`confidence=0`).
  A target-bound peptide trades internal hbonds for interface hbonds, so this is
  recorded as a metric rather than used to discard designs. Raise its confidence
  back to `1` (and tune `hbond_cutoff`) if you want it as a gate.
- `shape` (shape complementarity / interface area) — reported.

## Running

> **Run from the repo root** (`lib/`). The constraint files are referenced from
> `design_script.xml` by repo-root-relative path
> (`rosetta_only_packer_pred/...`), and Rosetta resolves them against the
> working directory it is launched from.

```bash
cd /home/vedansh/projects/random/git-install/lib

RUN_DIR=runs/packer_design_$(date +%Y%m%d_%H%M%S)

# 1. Build a one-row-per-backbone manifest. pepseq is a poly-G placeholder
#    (the XML designs from scratch and never reads it).
python3 rosetta_only_packer_pred/build_manifest.py \
  --backbones-dir perturbed_backbones/perturbed_backbones \
  --out "$RUN_DIR/jobs.tsv"

# 2. Run the packer design. Use --nstruct > 1 to draw several independent
#    designs per backbone (each is a separate stochastic trajectory).
scripts/run_rosetta_batch_threadseq.sh \
  --config config/rosetta.env \
  --jobs "$RUN_DIR/jobs.tsv" \
  --xml rosetta_only_packer_pred/design_script.xml \
  --run-dir "$RUN_DIR" \
  --peptide-chain A \
  --parallel 4 --nstruct 4

# 3. Collect scores (dG_separated, shape, etc.).
python3 scripts/collect_scores.py \
  --jobs "$RUN_DIR/jobs.tsv" \
  --run-dir "$RUN_DIR" \
  --out "$RUN_DIR/summary.tsv"

# 4. Bond-geometry audit on the designed peptides.
python3 scripts/geometry_check.py \
  --root "$RUN_DIR" \
  --chain A \
  --exclude input_scaffold.pdb \
  --summary-out "$RUN_DIR/geometry_summary.tsv" \
  --issues-out "$RUN_DIR/geometry_issues.tsv"
```

The batch runner reuses the existing
[../scripts/run_rosetta_batch_threadseq.sh](../scripts/run_rosetta_batch_threadseq.sh):
it auto-discovers the NCAA params in [../ncaa_params/](../ncaa_params/), runs the
cyclic backbone-continuity preflight on chain A, snapshots each scaffold, and
writes `score.sc` / `status.env` per job. The design step is compute-heavy
(three relax + three design rounds with `ex1/ex2` rotamers, `voids_penalty`, and
`buried_unsatisfied_penalty`), so expect several minutes per `nstruct`.

## Tuning the design requirements

- Edit [desired_makeup_MMP8.comp](desired_makeup_MMP8.comp) to change the
  composition targets (the `ABSOLUTE` count and `PENALTIES` curve per residue
  group). `TYPE` codes are residue `name3` values (`DPR` = D-Pro, `DHY` =
  D-Hyp).
- Edit [desired_peptide.charge](desired_peptide.charge) to change the target net
  charge.
- The **hard gates** (`min_proline_like`, `min_acidic`) hold the non-negotiable
  counts. If you retarget the proline-like requirement, update *both* the `.comp`
  `ABSOLUTE` (soft guidance) and the gate's `min_residue_count` (hard floor) so
  they agree — otherwise the packer is steered to one count but judged at another.
- The per-layer allowed residue sets live in the `*_restrictions`
  `RestrictToSpecifiedBaseResidueTypes` task operations in the XML; keep the
  proline-like types in them so the composition constraint stays satisfiable.

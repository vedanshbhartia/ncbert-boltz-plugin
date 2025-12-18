# Chirality Correction Module

Detects and corrects chirality violations in protein structure predictions by inverting residue stereochemistry and re-running predictions with corrected structures as templates.

## Quick Start

```python
from chirality_correction import ChiralityCorrector

corrector = ChiralityCorrector(cache_dir="/cache")

results = corrector.check_and_correct(
    prediction_file=Path("prediction.cif"),
    yaml_input=Path("input.yaml"),
    expected_sequence="DWWPLAFEALLR",
    output_dir=Path("corrected"),
    chain_id='B',
    expected_chirality='D'
)
```

## How It Works

1. Analyzes chirality using both CCD codes and N-C-CA-CB dihedral angles
2. Corrects violations by reflecting CB atom across CA: `CB' = CA - (CB - CA)`
3. Re-runs Boltz with corrected structure as template
4. Iterates until violations < 10% or max iterations (default 3)

## Key Functions

**`analyze_chirality_violations(structure_file, expected_sequence, chain_id, expected_chirality)`**

Returns `(ccd_violation_rate, geom_violation_rate, detailed_analysis)`. Use for analysis without correction.

**`correct_structure_chirality(input_file, output_file, expected_sequence, chain_id, expected_chirality)`**

Returns `(num_corrections, output_path)`. Corrects violations in-place.

**`ChiralityCorrector.check_and_correct(...)`**

Full iterative pipeline. Returns dict with per-iteration results and convergence status.

## Integration Example

```python
from main import ModifiedSequenceBoltzPipeline
from chirality_correction import ChiralityCorrector

# Run prediction
pipeline = ModifiedSequenceBoltzPipeline(work_dir=Path("results"))
output_dir = pipeline.process_single_sequence("DWWP[DLE]AFEALLR", "peptide")

# Correct chirality
corrector = ChiralityCorrector()
best_model = sorted(output_dir.rglob("*_model_*.cif"))[0]

results = corrector.check_and_correct(
    prediction_file=best_model,
    yaml_input=pipeline.yaml_dir / "peptide.yaml",
    expected_sequence="DWWPLAFEALLR",
    output_dir=Path("corrected"),
    expected_chirality='D'
)
```

## Notes

- Both CCD and geometric methods must agree for correction to trigger
- Only corrects residues with CA and CB atoms (excludes glycine)
- Template-based re-prediction doesn't always improve results
- Requires BioPython, NumPy, PyYAML, and Boltz

## Configuration

```python
ChiralityCorrector(
    cache_dir=None,              # Boltz cache directory
    use_msa_server=True,         # Use MSA server for re-runs
    diffusion_samples=5,         # Samples per re-run
    recycling_steps=10           # Recycling steps per re-run
)

corrector.check_and_correct(
    ...,
    max_iterations=3,            # Max correction iterations
    violation_threshold=0.1      # Acceptable violation rate (10%)
)
```

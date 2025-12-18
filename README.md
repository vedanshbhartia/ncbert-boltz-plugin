# Boltz Modified Sequence Predictor

Pipeline for running Boltz structure predictions on peptide sequences with non-canonical amino acid modifications.

## Configuration

Edit `config.yaml` to set paths and parameters:
- `cache_dir`: Cache directory for Boltz (null for default)
- `target_pdb_path`: Path to target PDB file for docking (null to skip)
- `use_msa_server`: Use MSA server for alignment generation
- `diffusion_samples`: Number of diffusion samples
- `recycling_steps`: Number of recycling steps
- `work_dir`: Output directory

## Installation

```bash
pip install -r requirements.txt
```

For installing boltz, check [https://github.com/jwohlwend/boltz](https://github.com/jwohlwend/boltz)
Installing the latest version:
```
git clone https://github.com/jwohlwend/boltz.git
cd boltz; pip install -e .[cuda]
```

## Usage

### As a script

```bash
python main.py
```

### As a module

```python
from main import ModifiedSequenceBoltzPipeline, load_config

config = load_config()
pipeline = ModifiedSequenceBoltzPipeline(
    work_dir=config['work_dir'],
    boltz_cache=config['cache_dir'],
    use_msa_server=config['use_msa_server']
)

# Process a single sequence
result = pipeline.process_single_sequence(
    seq_string="RGDG[MSE]GCGV",
    name="peptide_1",
    run_boltz=True,
    target_pdb_path=Path("target.pdb")
)
```

## Sequence Format

Sequences use bracket notation for modified residues:
```
RGDG[MSE]GCGV  # Selenomethionine at position 5
```

Supported modifications are defined in the `CCD_TO_BASE_AA` dictionary.

## Functions

### `load_config(config_path=None)`
Loads configuration from YAML file. Returns dict with config values.

### `ModifiedSequenceParser`
Parses sequences with CCD modifications in bracket notation.
- `parse_sequence(seq_string)`: Returns (base_sequence, modifications)
- `parse_file(file_path)`: Parses multiple sequences from file
- `parse_string(multi_seq_string)`: Parses multiple sequences from string

### `PDBParser`
Handles PDB file operations using Gemmi.
- `extract_sequence(pdb_path, chain_id=None)`: Extracts amino acid sequence from PDB
- `ensure_pdb_has_seqres(pdb_path, convert_to_standard=True)`: Ensures PDB has SEQRES records and optionally converts non-canonical residues

### `BoltzInputGenerator`
Generates Boltz-compatible YAML input files.
- `create_single_protein_yaml(sequence, modifications, ...)`: Creates YAML for single peptide
- `create_multi_protein_yaml(sequences_and_mods, ...)`: Creates YAML for multiple peptides
- `save_yaml(yaml_dict, output_path)`: Saves YAML to file

### `BoltzRunner`
Runs Boltz predictions.
- `run_prediction(input_yaml, output_dir, ...)`: Runs Boltz on single input
- `batch_predict(input_yamls, output_base_dir, ...)`: Runs Boltz on multiple inputs

### `ModifiedSequenceBoltzPipeline`
Complete pipeline from sequence parsing to prediction.
- `process_single_sequence(seq_string, name, ...)`: Processes one sequence through full pipeline
- `process_multiple_sequences(seq_strings, names, ...)`: Processes multiple sequences
- `process_file(file_path, ...)`: Processes sequences from file

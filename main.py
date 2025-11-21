#!/usr/bin/env python3
"""
Modified Sequence to Boltz Pipeline
Parses sequences with CCD modifications in bracket notation and runs Boltz predictions

Format: SEQUENCE[CCD]SEQUENCE[CCD]
Example: FGPLPNGEILDTYGHDT[CME]

Author: Modular design for reusability
"""

import re
import yaml
import subprocess
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# ============================================================================
# MODULE 1: CCD Code Mappings
# ============================================================================

# Mapping of CCD codes to their base amino acids
CCD_TO_BASE_AA = {
    # Modified cysteines
    'CME': 'C',  # S-methylcysteine
    'CSD': 'C',  # Cysteine sulfinic acid
    'CSO': 'C',  # S-hydroxycysteine
    'CSP': 'C',  # Cysteine persulfide
    'CSS': 'C',  # S-sulfocysteine
    'CSW': 'C',  # Cysteine-S-dioxide
    'CSX': 'C',  # S-oxy cysteine
    'OCS': 'C',  # Cysteinesulfonic acid
    
    # Modified lysines
    'MLY': 'K',  # Dimethyllysine
    'M3L': 'K',  # Trimethyllysine
    'ALY': 'K',  # Acetyllysine
    'LYZ': 'K',  # 5-hydroxylysine
    
    # Modified methionines
    'FME': 'M',  # N-formylmethionine
    'MSE': 'M',  # Selenomethionine
    'DME': 'M',  # N,N-dimethylmethionine
    
    # Modified threonines
    'TPO': 'T',  # Phosphothreonine
    'DTH': 'T',  # D-Threonine (different chirality)
    
    # Modified serines
    'SEP': 'S',  # Phosphoserine
    'DSE': 'S',  # D-Serine
    
    # Modified tyrosines
    'PTR': 'Y',  # Phosphotyrosine
    'TYS': 'Y',  # O-sulfo-L-tyrosine
    'DTY': 'Y',  # D-Tyrosine
    
    # Modified prolines
    'HYP': 'P',  # Hydroxyproline
    'DPR': 'P',  # D-Proline
    
    # Modified histidines
    'HIE': 'H',  # ND1-protonated histidine
    'HID': 'H',  # NE2-protonated histidine
    'HIP': 'H',  # Doubly protonated histidine
    'DHI': 'H',  # D-Histidine
    
    # Other modified amino acids
    'ABA': 'A',  # Alpha-aminobutyric acid (Alanine-like)
    'NLE': 'L',  # Norleucine
    'NVA': 'V',  # Norvaline
    'TRO': 'W',  # 2-carboxy-tryptophan
    'DTR': 'W',  # D-Tryptophan
    
    # D-amino acids
    'DAL': 'A', 'DAR': 'R', 'DAS': 'N', 'DSP': 'D', 'DCY': 'C',
    'DGL': 'E', 'DGN': 'Q', 'DIL': 'I', 'DLE': 'L', 'DLY': 'K',
    'DPN': 'F', 'DVA': 'V',
    
    # Ambiguous
    'ASX': 'B',  # Asparagine or Aspartic acid
    'GLX': 'Z',  # Glutamine or Glutamic acid
    'UNK': 'X',  # Unknown
    'XAA': 'X',  # Unknown amino acid
}

# ============================================================================
# MODULE 2: Sequence Parsing
# ============================================================================

class ModifiedSequenceParser:
    """Parse sequences with CCD modifications in bracket notation"""
    
    def __init__(self, ccd_mapping: Dict[str, str] = None):
        self.ccd_mapping = ccd_mapping or CCD_TO_BASE_AA
        self.pattern = re.compile(r'\[([A-Z0-9]+)\]')
    
    def parse_sequence(self, seq_string: str) -> Tuple[str, List[Dict]]:
        """
        Parse a sequence with CCD modifications
        
        Args:
            seq_string: Sequence like "SEQUENCE[CCD]SEQUENCE[CCD]"
        
        Returns:
            Tuple of (base_sequence, modifications_list)
            - base_sequence: Standard amino acid sequence
            - modifications_list: List of {position: int, ccd: str}
        """
        base_sequence = ""
        modifications = []
        position = 1
        
        # Split by CCD codes but keep them
        parts = self.pattern.split(seq_string)
        
        for i, part in enumerate(parts):
            if i % 2 == 0:
                # Regular sequence
                base_sequence += part
                position += len(part)
            else:
                # CCD code in brackets
                ccd_code = part
                base_aa = self.ccd_mapping.get(ccd_code, 'X')
                
                # Add base amino acid to sequence
                base_sequence += base_aa
                
                # Record modification at current position
                modifications.append({
                    'position': position,
                    'ccd': ccd_code
                })
                position += 1
        
        return base_sequence, modifications
    
    def parse_file(self, file_path: Path) -> List[Tuple[str, List[Dict]]]:
        """
        Parse a file with multiple sequences (one per line)
        
        Returns:
            List of (base_sequence, modifications) tuples
        """
        sequences = []
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    seq, mods = self.parse_sequence(line)
                    sequences.append((seq, mods))
        
        return sequences
    
    def parse_string(self, multi_seq_string: str) -> List[Tuple[str, List[Dict]]]:
        """Parse multiple sequences from a string (one per line)"""
        sequences = []
        
        for line in multi_seq_string.strip().split('\n'):
            line = line.strip()
            if line and not line.startswith('#'):
                seq, mods = self.parse_sequence(line)
                sequences.append((seq, mods))
        
        return sequences

# ============================================================================
# MODULE 3: Boltz YAML Generation
# ============================================================================

class BoltzInputGenerator:
    """Generate Boltz-compatible YAML input files"""
    
    @staticmethod
    def create_single_protein_yaml(
        sequence: str,
        modifications: List[Dict],
        chain_id: str = 'A',
        msa_path: Optional[str] = None
    ) -> Dict:
        """
        Create YAML structure for a single protein
        
        Args:
            sequence: Base amino acid sequence
            modifications: List of {position: int, ccd: str}
            chain_id: Chain identifier
            msa_path: Optional path to pre-computed MSA
        
        Returns:
            Dictionary representing YAML structure
        """
        protein_entry = {
            'id': chain_id,
            'sequence': sequence
        }
        
        # Add modifications if present
        if modifications:
            protein_entry['modifications'] = modifications
        
        # Add MSA path if provided
        if msa_path:
            protein_entry['msa'] = msa_path
        
        return {
            'sequences': [
                {'protein': protein_entry}
            ]
        }
    
    @staticmethod
    def create_multi_protein_yaml(
        sequences_and_mods: List[Tuple[str, List[Dict]]],
        chain_ids: Optional[List[str]] = None,
        msa_paths: Optional[List[str]] = None
    ) -> Dict:
        """
        Create YAML structure for multiple proteins
        
        Args:
            sequences_and_mods: List of (sequence, modifications) tuples
            chain_ids: Optional list of chain IDs (default: A, B, C, ...)
            msa_paths: Optional list of MSA paths
        
        Returns:
            Dictionary representing YAML structure
        """
        if chain_ids is None:
            chain_ids = [chr(65 + i) for i in range(len(sequences_and_mods))]
        
        if msa_paths is None:
            msa_paths = [None] * len(sequences_and_mods)
        
        sequences = []
        
        for i, (seq, mods) in enumerate(sequences_and_mods):
            protein_entry = {
                'id': chain_ids[i],
                'sequence': seq
            }
            
            if mods:
                protein_entry['modifications'] = mods
            
            if msa_paths[i]:
                protein_entry['msa'] = msa_paths[i]
            
            sequences.append({'protein': protein_entry})
        
        return {'sequences': sequences}
    
    @staticmethod
    def save_yaml(yaml_dict: Dict, output_path: Path):
        """Save YAML dictionary to file"""
        with open(output_path, 'w') as f:
            yaml.dump(yaml_dict, f, default_flow_style=False, sort_keys=False)
        
        logging.info(f"Created YAML input: {output_path}")

# ============================================================================
# MODULE 4: Boltz Execution
# ============================================================================

class BoltzRunner:
    """Run Boltz predictions"""
    
    def __init__(
        self,
        cache_dir: str = "/scratch/vb1467/.cache",
        use_msa_server: bool = True,
        diffusion_samples: int = 5,
        recycling_steps: int = 10
    ):
        self.cache_dir = cache_dir
        self.use_msa_server = use_msa_server
        self.diffusion_samples = diffusion_samples
        self.recycling_steps = recycling_steps
    
    def run_prediction(
        self,
        input_yaml: Path,
        output_dir: Path,
        additional_args: Optional[List[str]] = None
    ) -> bool:
        """
        Run Boltz prediction
        
        Args:
            input_yaml: Path to input YAML file
            output_dir: Output directory for results
            additional_args: Optional additional command-line arguments
        
        Returns:
            True if successful, False otherwise
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            'boltz', 'predict', str(input_yaml),
            # '--cache', self.cache_dir,
            '--accelerator', 'cpu',
            '--use_potentials',
            '--out_dir', str(output_dir),
            '--diffusion_samples', str(self.diffusion_samples),
            '--recycling_steps', str(self.recycling_steps)
        ]
        
        if self.use_msa_server:
            cmd.append('--use_msa_server')
        
        if additional_args:
            cmd.extend(additional_args)
        
        logging.info(f"Running Boltz: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False  # Don't raise on non-zero exit
            )
            
            # Check if prediction files were created
            prediction_files = list(output_dir.rglob("*_model_*.cif"))
            
            if prediction_files:
                logging.info(f"Boltz prediction completed: {len(prediction_files)} models generated")
                return True
            else:
                logging.error(f"No prediction files generated")
                if result.stderr:
                    logging.error(f"Error output: {result.stderr[:500]}")
                return False
                
        except Exception as e:
            logging.error(f"Error running Boltz: {e}")
            return False
    
    def batch_predict(
        self,
        input_yamls: List[Path],
        output_base_dir: Path,
        name_prefix: str = "prediction"
    ) -> Dict[str, bool]:
        """
        Run Boltz predictions on multiple inputs
        
        Returns:
            Dictionary mapping input file names to success status
        """
        results = {}
        
        for i, yaml_file in enumerate(input_yamls):
            name = f"{name_prefix}_{i}"
            output_dir = output_base_dir / name
            
            logging.info(f"[{i+1}/{len(input_yamls)}] Processing {yaml_file.name}")
            success = self.run_prediction(yaml_file, output_dir)
            results[yaml_file.name] = success
        
        return results

# ============================================================================
# MODULE 5: Complete Pipeline
# ============================================================================

class ModifiedSequenceBoltzPipeline:
    """Complete pipeline: Parse → Generate YAML → Run Boltz"""
    
    def __init__(
        self,
        work_dir: Path = Path("boltz_modified_sequences"),
        ccd_mapping: Optional[Dict[str, str]] = None,
        boltz_cache: str = "/scratch/vb1467/.cache",
        use_msa_server: bool = True
    ):
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(exist_ok=True)
        
        # Create subdirectories
        self.yaml_dir = self.work_dir / "yaml_inputs"
        self.predictions_dir = self.work_dir / "predictions"
        self.yaml_dir.mkdir(exist_ok=True)
        self.predictions_dir.mkdir(exist_ok=True)
        
        # Initialize modules
        self.parser = ModifiedSequenceParser(ccd_mapping)
        self.yaml_gen = BoltzInputGenerator()
        self.runner = BoltzRunner(
            cache_dir=boltz_cache,
            use_msa_server=use_msa_server
        )
    
    def process_single_sequence(
        self,
        seq_string: str,
        name: str,
        run_boltz: bool = True
    ) -> Optional[Path]:
        """
        Process a single sequence through the complete pipeline
        
        Args:
            seq_string: Sequence with CCD modifications
            name: Name for output files
            run_boltz: Whether to run Boltz prediction
        
        Returns:
            Path to output directory if successful
        """
        logging.info(f"Processing sequence: {name}")
        
        # Parse sequence
        sequence, modifications = self.parser.parse_sequence(seq_string)
        
        logging.info(f"  Base sequence: {sequence} ({len(sequence)} residues)")
        logging.info(f"  Modifications: {len(modifications)}")
        for mod in modifications:
            logging.info(f"    Position {mod['position']}: {mod['ccd']}")
        
        # Generate YAML
        yaml_dict = self.yaml_gen.create_single_protein_yaml(
            sequence, modifications
        )
        
        yaml_file = self.yaml_dir / f"{name}.yaml"
        self.yaml_gen.save_yaml(yaml_dict, yaml_file)
        
        # Run Boltz if requested
        if run_boltz:
            output_dir = self.predictions_dir / name
            success = self.runner.run_prediction(yaml_file, output_dir)
            
            if success:
                return output_dir
            else:
                return None
        
        return yaml_file
    
    def process_multiple_sequences(
        self,
        seq_strings: List[str],
        names: Optional[List[str]] = None,
        run_boltz: bool = True
    ) -> Dict[str, Optional[Path]]:
        """
        Process multiple sequences
        
        Args:
            seq_strings: List of sequences with CCD modifications
            names: Optional list of names (default: seq_0, seq_1, ...)
            run_boltz: Whether to run Boltz predictions
        
        Returns:
            Dictionary mapping names to output directories
        """
        if names is None:
            names = [f"seq_{i}" for i in range(len(seq_strings))]
        
        results = {}
        
        for seq_string, name in zip(seq_strings, names):
            result = self.process_single_sequence(seq_string, name, run_boltz)
            results[name] = result
        
        return results
    
    def process_file(
        self,
        file_path: Path,
        run_boltz: bool = True
    ) -> Dict[str, Optional[Path]]:
        """
        Process sequences from a file
        
        Args:
            file_path: Path to file with sequences (one per line)
            run_boltz: Whether to run Boltz predictions
        
        Returns:
            Dictionary mapping sequence indices to output directories
        """
        logging.info(f"Processing file: {file_path}")
        
        sequences_and_mods = self.parser.parse_file(file_path)
        
        logging.info(f"Found {len(sequences_and_mods)} sequences")
        
        results = {}
        
        for i, (sequence, modifications) in enumerate(sequences_and_mods):
            name = f"{file_path.stem}_seq_{i}"
            
            yaml_dict = self.yaml_gen.create_single_protein_yaml(
                sequence, modifications
            )
            
            yaml_file = self.yaml_dir / f"{name}.yaml"
            self.yaml_gen.save_yaml(yaml_dict, yaml_file)
            
            if run_boltz:
                output_dir = self.predictions_dir / name
                success = self.runner.run_prediction(yaml_file, output_dir)
                results[name] = output_dir if success else None
            else:
                results[name] = yaml_file
        
        return results

# ============================================================================
# MAIN / EXAMPLES
# ============================================================================

def main():
    """Main function with example usage"""
    
    # Example sequence from the user
    example_sequence = """
FGPLPNGEILDTYGHDT[CME]
CG[MLY]LPGGWNSITFGVGNR
GGPMADGEL[FME]DPRG[ABA]GVW
HGPVAEGKIDF[TPO]AGHGGD
CGDME[MLY]GEVA[OCS]SLGWGFY
[CSD]GAMPGGMDPDPSGYG[DTH]H
NGYPPGGMWQFPMGYGAH
NGPFKSSEFVLTGGAGGF
WGPLDKGVSCDPNG[MLY]GCQ
FGYLNNYHWFNGAGRGGY
""".strip()
    
    # Initialize pipeline
    pipeline = ModifiedSequenceBoltzPipeline(
        work_dir=Path("boltz_modified_results"),
        use_msa_server=True
    )
    
    # Parse sequences
    sequences = pipeline.parser.parse_string(example_sequence)
    
    print("="*80)
    print("PARSED SEQUENCES WITH MODIFICATIONS")
    print("="*80)
    
    for i, (seq, mods) in enumerate(sequences):
        print(f"\nSequence {i+1}:")
        print(f"  Length: {len(seq)}")
        print(f"  Sequence: {seq}")
        print(f"  Modifications ({len(mods)}):")
        for mod in mods:
            base_aa = seq[mod['position']-1]
            print(f"    Position {mod['position']:3d}: {base_aa} → [{mod['ccd']}]")
    
    print("\n" + "="*80)
    print("RUNNING BOLTZ PREDICTIONS")
    print("="*80)
    
    # Process all sequences
    seq_strings = example_sequence.strip().split('\n')
    results = pipeline.process_multiple_sequences(
        seq_strings,
        names=[f"SA2989311_cand2_pred_{i+1}" for i in range(len(seq_strings))],
        run_boltz=True
    )
    
    # Print summary
    print("\n" + "="*80)
    print("RESULTS SUMMARY")
    print("="*80)
    
    success_count = sum(1 for v in results.values() if v is not None)
    print(f"Successful predictions: {success_count}/{len(results)}")
    print(f"Results directory: {pipeline.work_dir}")
    print(f"YAML inputs: {pipeline.yaml_dir}")
    print(f"Predictions: {pipeline.predictions_dir}")
    
    for name, result in results.items():
        if result:
            print(f"  ✓ {name}: {result}")
        else:
            print(f"  ✗ {name}: Failed")

if __name__ == "__main__":
    main()
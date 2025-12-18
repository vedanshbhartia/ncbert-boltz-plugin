#!/usr/bin/env python3

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
# Configuration Loading
# ============================================================================

def load_config(config_path: Optional[Path] = None) -> Dict:
    """
    Load configuration from YAML file

    Args:
        config_path: Path to config file (default: config.yaml in same directory as script)

    Returns:
        Dictionary with configuration values
    """
    if config_path is None:
        config_path = Path(__file__).parent / "config.yaml"

    config = {
        'cache_dir': None,
        'target_pdb_path': None,
        'use_msa_server': True,
        'diffusion_samples': 5,
        'recycling_steps': 10,
        'work_dir': 'boltz_modified_results'
    }

    if config_path.exists():
        with open(config_path, 'r') as f:
            loaded_config = yaml.safe_load(f)
            if loaded_config:
                config.update(loaded_config)

    return config

# ============================================================================
# Configuration Constants
# ============================================================================

# Default values (can be overridden by config file)
DEFAULT_CACHE_DIR = None
DEFAULT_TARGET_PDB = None
DEFAULT_USE_MSA_SERVER = True
DEFAULT_DIFFUSION_SAMPLES = 5
DEFAULT_RECYCLING_STEPS = 10
DEFAULT_WORK_DIR = "boltz_modified_results"

# Mapping from 1-letter to 3-letter standard amino acid codes
ONE_TO_THREE_AA = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
    'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
    'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
    'X': 'UNK', 'B': 'ASX', 'Z': 'GLX'
}

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
    'DPN': 'F', 'DVA': 'V', 'DGU': 'E',  # D-Glutamic acid

    # Additional non-standard residues commonly found in PDBs
    'NLU': 'L',  # N-methyl-leucine or similar leucine variant
    'DCS': 'C',  # D-Cysteine
    'B30': 'X',  # Unknown/modified residue - placeholder
    'HYP': 'P',  # Already listed above but confirming

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
# MODULE 2.5: PDB Parsing
# ============================================================================

class PDBParser:
    """PDB parser to extract sequences using Gemmi"""

    @staticmethod
    def extract_sequence(pdb_path: Path, chain_id: Optional[str] = None) -> str:
        """
        Extract amino acid sequence from PDB file using Gemmi
        
        Args:
            pdb_path: Path to PDB file
            chain_id: Optional chain ID to extract. If None, extracts first chain found.
            
        Returns:
            Sequence string
        """
        try:
            import gemmi
        except ImportError:
            logging.error("Gemmi not found. Please install gemmi or run in an environment with Boltz installed.")
            raise

        try:
            st = gemmi.read_structure(str(pdb_path))
            st.setup_entities()
            
            target_chain_id = chain_id
            
            # Find the entity corresponding to the chain
            # gemmi entities are not directly mapped by chain ID in a simple dict
            # We iterate entities and check their subchains
            
            target_entity = None
            
            if target_chain_id:
                for ent in st.entities:
                    if ent.entity_type == gemmi.EntityType.Polymer:
                        # Check if this entity includes our chain
                        # subchains are named like "A" or "A1" etc.
                        for sub in ent.subchains:
                            # subchain name usually starts with chain ID
                            if sub.startswith(target_chain_id):
                                target_entity = ent
                                break
                    if target_entity:
                        break
                
                if not target_entity:
                    # Fallback: try to find any polymer if specific chain not found/matched
                    # Or maybe the chain ID is different in subchains
                    logging.warning(f"Could not find polymer entity for chain {target_chain_id} using setup_entities. Falling back to manual extraction.")
                    # Fallback logic below
            else:
                # Use first polymer entity
                for ent in st.entities:
                    if ent.entity_type == gemmi.EntityType.Polymer:
                        target_entity = ent
                        break
            
            if target_entity:
                # gemmi returns list of strings (residue names)
                # These can be 1-letter or 3-letter codes depending on source
                seq_list = target_entity.full_sequence
                
                # If full_sequence is empty (no SEQRES), extract from residues in subchains
                if not seq_list:
                    logging.info("Entity full_sequence is empty (no SEQRES). Extracting from residues.")
                    seq_list = []
                    for sub in target_entity.subchains:
                        # Simpler: iterate subchains of entity and find them in structure
                        for model in st:
                            for chain in model:
                                for res in chain:
                                    if res.subchain == sub:
                                        if res.is_water(): continue
                                        code = gemmi.find_tabulated_residue(res.name).one_letter_code
                                        if code.strip():
                                            seq_list.append(code.upper())
                
                # Convert to 1-letter codes if necessary
                final_seq = []
                for item in seq_list:
                    if len(item) > 1:
                        code = gemmi.find_tabulated_residue(item).one_letter_code
                        if code.strip():
                            final_seq.append(code.upper())
                        else:
                            # Fallback for unknown? Use X?
                            final_seq.append('X')
                    else:
                        final_seq.append(item.upper())
                        
                sequence = "".join(final_seq)
                logging.info(f"Extracted sequence for chain {target_chain_id or 'first'}: {sequence}")
                return sequence
            
            # FALLBACK: Manual extraction if setup_entities fails to give us what we want
            # (e.g. if connectivity is bad and gemmi splits it or fails to make polymer)
            logging.info("Falling back to manual extraction (ignoring connectivity).")
            model = st[0]
            if target_chain_id and target_chain_id in model:
                chain = model[target_chain_id]
            else:
                chain = model[0]
                
            seq = []
            for res in chain:
                if res.is_water(): continue
                code = gemmi.find_tabulated_residue(res.name).one_letter_code
                if code.strip():
                    seq.append(code.upper())
            return "".join(seq)

        except Exception as e:
            logging.error(f"Error parsing PDB {pdb_path}: {e}")
            raise

    @staticmethod
    def ensure_pdb_has_seqres(pdb_path: Path, convert_to_standard: bool = True) -> Path:
        """
        Ensure PDB file has SEQRES records (populated entities).
        If not, create a fixed version with SEQRES derived from residues.
        Optionally converts non-canonical amino acids to standard ones.

        Args:
            pdb_path: Path to original PDB
            convert_to_standard: If True, convert non-canonical residues to standard amino acids

        Returns:
            Path to PDB file (original or fixed)
        """
        try:
            import gemmi

            st = gemmi.read_structure(str(pdb_path))
            st.setup_entities()

            # Check if any polymer entity has empty sequence or has non-standard residues
            needs_fix = False

            for ent in st.entities:
                if ent.entity_type == gemmi.EntityType.Polymer:
                    if not ent.full_sequence:
                        needs_fix = True
                        break
                    # Check for non-standard residues in existing SEQRES
                    if convert_to_standard:
                        for res_name in ent.full_sequence:
                            if res_name in CCD_TO_BASE_AA:
                                needs_fix = True
                                break

            if not needs_fix:
                return pdb_path

            logging.info(f"PDB {pdb_path} needs fixing (missing SEQRES or non-standard residues). Creating fixed version...")

            # Fix entities
            for ent in st.entities:
                if ent.entity_type == gemmi.EntityType.Polymer:
                    seq_list = []

                    # If entity has existing sequence, process it
                    if ent.full_sequence and convert_to_standard:
                        for res_name in ent.full_sequence:
                            if res_name in CCD_TO_BASE_AA:
                                # Convert to standard 1-letter code
                                standard_aa = CCD_TO_BASE_AA[res_name]
                                # Convert to 3-letter code for SEQRES
                                standard_3letter = ONE_TO_THREE_AA.get(standard_aa, 'UNK')
                                seq_list.append(standard_3letter)
                                logging.info(f"  Converted {res_name} -> {standard_3letter} ({standard_aa})")
                            else:
                                seq_list.append(res_name)

                        if seq_list:
                            ent.full_sequence = seq_list
                            logging.info(f"Converted entity {ent.name} with {len(seq_list)} residues")
                            continue

                    # Otherwise, extract from residues (original behavior)
                    if not ent.full_sequence:
                        for model in st:
                            for chain in model:
                                for res in chain:
                                    if res.subchain in ent.subchains:
                                        if res.is_water(): continue

                                        res_name = res.name

                                        # Convert non-canonical to standard if requested
                                        if convert_to_standard and res_name in CCD_TO_BASE_AA:
                                            standard_aa = CCD_TO_BASE_AA[res_name]
                                            standard_3letter = ONE_TO_THREE_AA.get(standard_aa, 'UNK')
                                            seq_list.append(standard_3letter)
                                            logging.info(f"  Converted {res_name} -> {standard_3letter} ({standard_aa})")
                                        else:
                                            seq_list.append(res_name)

                        ent.full_sequence = seq_list
                        logging.info(f"Populated entity {ent.name} with {len(seq_list)} residues")

            # Also convert ATOM records if converting to standard
            if convert_to_standard:
                # Convert ATOM records
                for model in st:
                    for chain in model:
                        for res in chain:
                            if res.is_water(): continue
                            if res.name in CCD_TO_BASE_AA:
                                standard_aa = CCD_TO_BASE_AA[res.name]
                                standard_3letter = ONE_TO_THREE_AA.get(standard_aa, 'UNK')
                                old_name = res.name
                                res.name = standard_3letter
                                logging.info(f"  Converted ATOM record: {old_name} -> {standard_3letter}")

            # Save fixed PDB
            fixed_path = pdb_path.parent / f"{pdb_path.stem}_fixed{pdb_path.suffix}"
            st.write_pdb(str(fixed_path))
            logging.info(f"Saved fixed PDB to {fixed_path}")

            return fixed_path

        except Exception as e:
            logging.error(f"Error ensuring SEQRES for {pdb_path}: {e}")
            # Return original if fix fails, let downstream handle it (or crash)
            return pdb_path

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
        msa_path: Optional[str] = None,
        target_pdb_path: Optional[Path] = None,
        target_chain_id: str = 'A',
        cyclic: bool = True
    ) -> Dict:
        """
        Create YAML structure for a single protein (peptide) potentially docked to a target
        
        Args:
            sequence: Base amino acid sequence of the peptide
            modifications: List of {position: int, ccd: str} for the peptide
            chain_id: Chain identifier for the peptide (default 'A' if no target, else 'B')
            msa_path: Optional path to pre-computed MSA for peptide
            target_pdb_path: Optional path to target PDB file
            target_chain_id: Chain identifier for the target (default 'A')
        
        Returns:
            Dictionary representing YAML structure
        """
        sequences_list = []
        
        # 1. Handle Target Protein (if provided)
        if target_pdb_path:
            # Extract sequence from PDB
            target_seq = PDBParser.extract_sequence(target_pdb_path)
            
            target_entry = {
                'protein': {
                    'id': target_chain_id,
                    'sequence': target_seq,
                    'msa': 'empty' # Usually empty for target if we use template, or could be auto
                }
            }
            sequences_list.append(target_entry)
            
            # If target is A, peptide should be B
            if chain_id == target_chain_id:
                chain_id = chr(ord(target_chain_id) + 1)
        
        # 2. Handle Peptide
        peptide_entry = {
            'id': chain_id,
            'sequence': sequence
        }
        
        # Add modifications if present
        if modifications:
            peptide_entry['modifications'] = modifications
        
        # Add MSA path if provided
        if msa_path:
            peptide_entry['msa'] = msa_path
        else:
             # If no MSA provided and we are using server, it's fine. 
             # If we want to force empty for peptide (often good for short peptides), we can do that.
             # For now leaving it to default (which might require MSA or use server)
             pass
        
        if cyclic:
            peptide_entry['cyclic'] = True

        sequences_list.append({'protein': peptide_entry})
        
        yaml_output = {
            'sequences': sequences_list
        }
        
        # 3. Add Templates if target provided
        if target_pdb_path:
            yaml_output['templates'] = [
                {
                    'pdb': str(target_pdb_path),
                    'chain_id': [target_chain_id],
                    # We assume the PDB has the chain we want. 
                    # If we extracted from chain A, we map to chain A.
                }
            ]
            
        return yaml_output
    
    @staticmethod
    def create_multi_protein_yaml(
        sequences_and_mods: List[Tuple[str, List[Dict]]],
        chain_ids: Optional[List[str]] = None,
        msa_paths: Optional[List[str]] = None,
        target_pdb_path: Optional[Path] = None
    ) -> Dict:
        """
        Create YAML structure for multiple proteins (e.g. multiple peptides or complex)
        Note: This method might need adjustment if we want to dock multiple peptides to ONE target.
        For now, assuming this is just for batch processing where each is separate, 
        BUT the current implementation puts them all in one YAML (complex).
        
        If target_pdb_path is provided, it's added as the first chain.
        """
        sequences_list = []
        current_chain_idx = 0
        
        # 1. Handle Target
        if target_pdb_path:
            target_seq = PDBParser.extract_sequence(target_pdb_path)
            target_chain = 'A'
            sequences_list.append({
                'protein': {
                    'id': target_chain,
                    'sequence': target_seq,
                    'msa': 'empty'
                }
            })
            current_chain_idx = 1 # Next chain starts at B
            
        if chain_ids is None:
            chain_ids = [chr(65 + current_chain_idx + i) for i in range(len(sequences_and_mods))]
        
        if msa_paths is None:
            msa_paths = [None] * len(sequences_and_mods)
        
        for i, (seq, mods) in enumerate(sequences_and_mods):
            protein_entry = {
                'id': chain_ids[i],
                'sequence': seq
            }
            
            if mods:
                protein_entry['modifications'] = mods
            
            if msa_paths[i]:
                protein_entry['msa'] = msa_paths[i]
            
            sequences_list.append({'protein': protein_entry})
        
        yaml_output = {'sequences': sequences_list}
        
        if target_pdb_path:
            yaml_output['templates'] = [{
                'pdb': str(target_pdb_path),
                'chain_id': ['A']
            }]
            
        return yaml_output
    
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
        cache_dir: Optional[str] = None,
        use_msa_server: bool = True,
        diffusion_samples: int = 5,
        recycling_steps: int = 10
    ):
        self.cache_dir = cache_dir if cache_dir is not None else DEFAULT_CACHE_DIR
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
        boltz_cache: Optional[str] = None,
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
            cache_dir=boltz_cache if boltz_cache is not None else DEFAULT_CACHE_DIR,
            use_msa_server=use_msa_server
        )
    
    def process_single_sequence(
        self,
        seq_string: str,
        name: str,
        run_boltz: bool = True,
        target_pdb_path: Optional[Path] = None
    ) -> Optional[Path]:
        """
        Process a single sequence through the complete pipeline
        
        Args:
            seq_string: Sequence with CCD modifications
            name: Name for output files
            run_boltz: Whether to run Boltz prediction
            target_pdb_path: Optional path to target PDB file
        
        Returns:
            Path to output directory if successful
        """
        logging.info(f"Processing sequence: {name}")
        
        # Ensure target PDB has SEQRES if provided
        fixed_target_pdb = None
        if target_pdb_path:
            fixed_target_pdb = PDBParser.ensure_pdb_has_seqres(target_pdb_path)
            if fixed_target_pdb != target_pdb_path:
                logging.info(f"Using fixed PDB: {fixed_target_pdb}")
        
        # Parse sequence
        sequence, modifications = self.parser.parse_sequence(seq_string)
        
        logging.info(f"  Base sequence: {sequence} ({len(sequence)} residues)")
        logging.info(f"  Modifications: {len(modifications)}")
        for mod in modifications:
            logging.info(f"    Position {mod['position']}: {mod['ccd']}")
        
        # Generate YAML
        yaml_dict = self.yaml_gen.create_single_protein_yaml(
            sequence, modifications, target_pdb_path=fixed_target_pdb or target_pdb_path
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
        run_boltz: bool = True,
        target_pdb_path: Optional[Path] = None
    ) -> Dict[str, Optional[Path]]:
        """
        Process multiple sequences
        
        Args:
            seq_strings: List of sequences with CCD modifications
            names: Optional list of names (default: seq_0, seq_1, ...)
            run_boltz: Whether to run Boltz predictions
            target_pdb_path: Optional path to target PDB file (applied to all)
        
        Returns:
            Dictionary mapping names to output directories
        """
        if names is None:
            names = [f"seq_{i}" for i in range(len(seq_strings))]
        
        # Ensure target PDB has SEQRES if provided (once for all)
        fixed_target_pdb = None
        if target_pdb_path:
            fixed_target_pdb = PDBParser.ensure_pdb_has_seqres(target_pdb_path)
            if fixed_target_pdb != target_pdb_path:
                logging.info(f"Using fixed PDB for batch: {fixed_target_pdb}")
        
        results = {}
        
        for seq_string, name in zip(seq_strings, names):
            # Pass fixed path directly to avoid re-checking
            result = self.process_single_sequence(
                seq_string, name, run_boltz, target_pdb_path=fixed_target_pdb or target_pdb_path
            )
            results[name] = result
        
        return results
    
    def process_file(
        self,
        file_path: Path,
        run_boltz: bool = True,
        target_pdb_path: Optional[Path] = None
    ) -> Dict[str, Optional[Path]]:
        """
        Process sequences from a file
        
        Args:
            file_path: Path to file with sequences (one per line)
            run_boltz: Whether to run Boltz predictions
            target_pdb_path: Optional path to target PDB file
        
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
                sequence, modifications, target_pdb_path=target_pdb_path
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

    # Load configuration
    config = load_config()

    # Example sequence from the user
    example_sequence = """
RGDG[DIL]GCGVSFKKYHGWA
CGSKFLGGHAHYTGKN[CSS]A
AGDNQGHGSMDLTN[CGU]CCV
CGDTKGMGK[DLY]GQSVDCCT
HGNGKCQGAA[DTY]GGTVNGW
GGPNYVMGAGTPHAVWNF
AYDDKGHGC[PCA]GKRDWHHC
AYSGQGTGRSG[ARG]DVVLHD
FGDRRGYGIGYDQN[YCM]NEF
GGPFQGGGR[DTH]HQYYVA[CSO]T
""".strip()

    # Initialize pipeline
    pipeline = ModifiedSequenceBoltzPipeline(
        work_dir=Path(config['work_dir']),
        boltz_cache=config['cache_dir'],
        use_msa_server=config['use_msa_server']
    )

    # Get target PDB from config
    target_pdb = None
    if config['target_pdb_path']:
        target_pdb = Path(config['target_pdb_path'])
        if not target_pdb.exists():
            logging.warning(f"Target PDB not found: {target_pdb}. Running without target.")
            target_pdb = None
    
    print("="*80)
    print("PARSED SEQUENCES WITH MODIFICATIONS")
    print("="*80)
    
    # Parse sequences just for display
    sequences = pipeline.parser.parse_string(example_sequence)
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
        names=[f"SA474935_cand1_pred_{i+1}" for i in range(len(seq_strings))],
        run_boltz=True,
        target_pdb_path=target_pdb
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

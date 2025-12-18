#!/usr/bin/env python3

import numpy as np
import subprocess
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# ============================================================================
# Chirality Detection and Correction Module
# ============================================================================

D_AMINO_ACID_CCD = {
    'A': 'DAL', 'R': 'DAR', 'N': 'DAS', 'D': 'DSP', 'C': 'DCY',
    'E': 'DGL', 'Q': 'DGN', 'G': 'GLY', 'H': 'DHI', 'I': 'DIL',
    'L': 'DLE', 'K': 'DLY', 'M': 'DME', 'F': 'DPN', 'P': 'DPR',
    'S': 'DSE', 'T': 'DTH', 'W': 'DTR', 'Y': 'DTY', 'V': 'DVA',
}

L_AMINO_ACID_CCD = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
    'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

# ============================================================================
# Chirality Analysis Functions
# ============================================================================

def determine_chirality_from_residue_name(residue_name: str) -> str:
    """Determine chirality from residue CCD code.

    Returns 'L', 'D', 'Achiral', or 'Unknown'.
    """
    d_amino_ccds = set(D_AMINO_ACID_CCD.values()) - {'GLY'}

    if residue_name in d_amino_ccds:
        return 'D'
    elif residue_name in L_AMINO_ACID_CCD:
        return 'L'
    elif residue_name == 'GLY':
        return 'Achiral'
    else:
        return 'Unknown'


def calculate_chirality_from_coords(ca_atom, n_atom, c_atom, cb_atom) -> Tuple[str, float]:
    """Calculate chirality from atomic coordinates using dihedral angle.

    Returns (chirality, dihedral_angle) where chirality is 'L', 'D', or 'Planar'
    and dihedral_angle is the N-CA-CB-C dihedral in degrees.

    Positive dihedral = L-amino acid
    Negative dihedral = D-amino acid
    """
    try:
        ca_coord = np.array(ca_atom.get_coord())
        n_coord = np.array(n_atom.get_coord())
        c_coord = np.array(c_atom.get_coord())
        cb_coord = np.array(cb_atom.get_coord()) if cb_atom else None

        if cb_coord is None:
            return 'Planar', 0.0

        def torsion_angle(p1, p2, p3, p4):
            """Calculate torsion angle between four points (N-C-CA-CB)."""
            b1 = p2 - p1
            b2 = p3 - p2
            b3 = p4 - p3

            n1 = np.cross(b1, b2)
            n1_norm = np.linalg.norm(n1)
            if n1_norm == 0:
                return 0.0
            n1 = n1 / n1_norm

            n2 = np.cross(b2, b3)
            n2_norm = np.linalg.norm(n2)
            if n2_norm == 0:
                return 0.0
            n2 = n2 / n2_norm

            x = np.dot(n1, n2)
            b2_normalized = b2 / np.linalg.norm(b2)
            y = np.dot(np.cross(n1, n2), b2_normalized)

            chi = np.arctan2(y, x)
            return np.degrees(chi)

        dihedral = torsion_angle(n_coord, c_coord, ca_coord, cb_coord)

        if abs(dihedral) < 10:
            chirality = 'Planar'
        elif dihedral > 0:
            chirality = 'L'
        else:
            chirality = 'D'

        return chirality, dihedral

    except Exception as e:
        logging.warning(f"Chirality calculation failed: {e}")
        return 'Planar', 0.0


def analyze_chirality_violations(
    structure_file: Path,
    expected_sequence: str,
    chain_id: str = 'B',
    expected_chirality: str = 'D'
) -> Tuple[float, float, List[Dict]]:
    """Analyze chirality violations in a structure file.

    Returns (ccd_violation_rate, geom_violation_rate, detailed_analysis).
    """
    try:
        from Bio.PDB import MMCIFParser
    except ImportError:
        logging.error("BioPython not found. Please install biopython.")
        raise

    parser = MMCIFParser(QUIET=True)

    try:
        structure = parser.get_structure('pred', structure_file)
        model = structure[0]

        if chain_id not in model:
            logging.error(f"Chain {chain_id} not found in structure")
            return None, None, []

        chain = model[chain_id]

        residues = []
        for residue in chain.get_residues():
            if residue.get_id()[0] == ' ' or residue.get_id()[0].startswith('H_'):
                residues.append(residue)

        ccd_violations = []
        geom_violations = []
        detailed_analysis = []

        for i, residue in enumerate(residues):
            if i >= len(expected_sequence):
                break

            expected_aa = expected_sequence[i]
            residue_name = residue.get_resname()

            ccd_chirality = determine_chirality_from_residue_name(residue_name)

            try:
                ca_atom = residue['CA'] if residue.has_id('CA') else None
                n_atom = residue['N'] if residue.has_id('N') else None
                c_atom = residue['C'] if residue.has_id('C') else None
                cb_atom = residue['CB'] if residue.has_id('CB') else None

                if all([ca_atom, n_atom, c_atom]):
                    geom_chirality, dihedral = calculate_chirality_from_coords(
                        ca_atom, n_atom, c_atom, cb_atom
                    )
                else:
                    geom_chirality = 'Missing_atoms'
                    dihedral = 0.0
            except Exception as e:
                geom_chirality = 'Error'
                dihedral = 0.0

            ccd_is_violation = (
                ccd_chirality != expected_chirality and
                ccd_chirality != 'Achiral'
            )
            geom_is_violation = (
                geom_chirality != expected_chirality and
                geom_chirality not in ['Achiral', 'Planar', 'Missing_atoms', 'Error']
            )

            if ccd_chirality != 'Achiral':
                ccd_violations.append(ccd_is_violation)
            if geom_chirality not in ['Achiral', 'Planar', 'Missing_atoms', 'Error']:
                geom_violations.append(geom_is_violation)

            detailed_analysis.append({
                'position': i+1,
                'expected_aa': expected_aa,
                'residue_name': residue_name,
                'ccd_chirality': ccd_chirality,
                'geom_chirality': geom_chirality,
                'dihedral_angle': dihedral,
                'expected_chirality': expected_chirality,
                'ccd_violation': ccd_is_violation,
                'geom_violation': geom_is_violation
            })

        ccd_violation_rate = sum(ccd_violations) / len(ccd_violations) if ccd_violations else 0
        geom_violation_rate = sum(geom_violations) / len(geom_violations) if geom_violations else 0

        return ccd_violation_rate, geom_violation_rate, detailed_analysis

    except Exception as e:
        logging.error(f"Error analyzing structure: {e}")
        return None, None, []


# ============================================================================
# Chirality Correction Functions
# ============================================================================

def invert_residue_chirality(residue) -> bool:
    """Invert the chirality of a residue by reflecting CB atom.

    Returns True if successful, False otherwise.
    """
    try:
        if not residue.has_id('CA') or not residue.has_id('CB'):
            return False

        ca_coord = np.array(residue['CA'].get_coord())
        cb_coord = np.array(residue['CB'].get_coord())

        ca_to_cb = cb_coord - ca_coord
        reflected_cb = ca_coord - ca_to_cb

        residue['CB'].set_coord(reflected_cb)

        return True

    except Exception as e:
        logging.warning(f"Failed to invert residue chirality: {e}")
        return False


def correct_structure_chirality(
    input_structure_file: Path,
    output_structure_file: Path,
    expected_sequence: str,
    chain_id: str = 'B',
    expected_chirality: str = 'D'
) -> Tuple[int, Path]:
    """Correct chirality violations in a structure.

    Returns (num_corrections, output_file_path).
    """
    try:
        from Bio.PDB import MMCIFParser, MMCIFIO
    except ImportError:
        logging.error("BioPython not found. Please install biopython.")
        raise

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('pred', input_structure_file)

    model = structure[0]
    if chain_id not in model:
        logging.error(f"Chain {chain_id} not found")
        return 0, None

    chain = model[chain_id]

    residues = []
    for residue in chain.get_residues():
        if residue.get_id()[0] == ' ' or residue.get_id()[0].startswith('H_'):
            residues.append(residue)

    num_corrections = 0

    for i, residue in enumerate(residues):
        if i >= len(expected_sequence):
            break

        expected_aa = expected_sequence[i]
        if expected_aa == 'G':
            continue

        residue_name = residue.get_resname()
        ccd_chirality = determine_chirality_from_residue_name(residue_name)

        if ccd_chirality != expected_chirality and ccd_chirality != 'Achiral':
            try:
                ca_atom = residue['CA'] if residue.has_id('CA') else None
                n_atom = residue['N'] if residue.has_id('N') else None
                c_atom = residue['C'] if residue.has_id('C') else None
                cb_atom = residue['CB'] if residue.has_id('CB') else None

                if all([ca_atom, n_atom, c_atom, cb_atom]):
                    geom_chirality, dihedral = calculate_chirality_from_coords(
                        ca_atom, n_atom, c_atom, cb_atom
                    )

                    if geom_chirality != expected_chirality and geom_chirality != 'Planar':
                        if invert_residue_chirality(residue):
                            num_corrections += 1
                            logging.info(
                                f"Corrected chirality at position {i+1} "
                                f"({expected_aa}): {geom_chirality} → {expected_chirality}"
                            )
            except Exception as e:
                logging.warning(f"Failed to correct position {i+1}: {e}")

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(output_structure_file))

    logging.info(f"Saved corrected structure to {output_structure_file}")
    logging.info(f"Total corrections made: {num_corrections}")

    return num_corrections, output_structure_file


# ============================================================================
# Iterative Correction Pipeline
# ============================================================================

class ChiralityCorrector:
    """Pipeline for checking and correcting chirality in predictions."""

    def __init__(
        self,
        cache_dir: Optional[str] = None,
        use_msa_server: bool = True,
        diffusion_samples: int = 5,
        recycling_steps: int = 10
    ):
        self.cache_dir = cache_dir
        self.use_msa_server = use_msa_server
        self.diffusion_samples = diffusion_samples
        self.recycling_steps = recycling_steps

    def check_and_correct(
        self,
        prediction_file: Path,
        yaml_input: Path,
        expected_sequence: str,
        output_dir: Path,
        chain_id: str = 'B',
        expected_chirality: str = 'D',
        max_iterations: int = 3,
        violation_threshold: float = 0.1
    ) -> Dict:
        """Check chirality and re-run prediction with corrections if needed.

        Returns dictionary with results from each iteration.
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        results = {
            'iterations': [],
            'final_violation_rate': None,
            'converged': False
        }

        current_prediction = prediction_file

        for iteration in range(max_iterations):
            logging.info(f"\n=== Iteration {iteration + 1}/{max_iterations} ===")

            ccd_rate, geom_rate, analysis = analyze_chirality_violations(
                current_prediction,
                expected_sequence,
                chain_id,
                expected_chirality
            )

            if ccd_rate is None or geom_rate is None:
                logging.error("Chirality analysis failed")
                break

            logging.info(f"CCD violation rate: {ccd_rate:.1%}")
            logging.info(f"Geometric violation rate: {geom_rate:.1%}")

            iteration_result = {
                'iteration': iteration + 1,
                'ccd_violation_rate': ccd_rate,
                'geom_violation_rate': geom_rate,
                'detailed_analysis': analysis
            }

            results['iterations'].append(iteration_result)

            if geom_rate <= violation_threshold:
                logging.info(
                    f"Chirality acceptable (≤{violation_threshold:.1%}). "
                    f"No correction needed."
                )
                results['final_violation_rate'] = geom_rate
                results['converged'] = True
                break

            logging.info(
                f"Chirality violations detected ({geom_rate:.1%}). "
                f"Applying corrections..."
            )

            corrected_file = output_dir / f"corrected_iter_{iteration + 1}.cif"
            num_corrections, corrected_path = correct_structure_chirality(
                current_prediction,
                corrected_file,
                expected_sequence,
                chain_id,
                expected_chirality
            )

            iteration_result['num_corrections'] = num_corrections
            iteration_result['corrected_file'] = str(corrected_path)

            if num_corrections == 0:
                logging.warning("No corrections made. Stopping.")
                results['final_violation_rate'] = geom_rate
                break

            if iteration < max_iterations - 1:
                logging.info("Re-running Boltz with corrected structure as template...")

                new_prediction = self._run_boltz_with_template(
                    yaml_input,
                    corrected_path,
                    output_dir / f"rerun_iter_{iteration + 2}",
                    chain_id
                )

                if new_prediction:
                    current_prediction = new_prediction
                    iteration_result['new_prediction'] = str(new_prediction)
                else:
                    logging.error("Re-run failed. Stopping.")
                    results['final_violation_rate'] = geom_rate
                    break
            else:
                results['final_violation_rate'] = geom_rate

        return results

    def _run_boltz_with_template(
        self,
        yaml_input: Path,
        template_structure: Path,
        output_dir: Path,
        template_chain_id: str
    ) -> Optional[Path]:
        """Run Boltz prediction using a corrected structure as template."""
        import yaml

        with open(yaml_input, 'r') as f:
            yaml_dict = yaml.safe_load(f)

        if 'templates' not in yaml_dict:
            yaml_dict['templates'] = []

        yaml_dict['templates'].append({
            'pdb': str(template_structure),
            'chain_id': [template_chain_id]
        })

        templated_yaml = output_dir / 'input_with_template.yaml'
        with open(templated_yaml, 'w') as f:
            yaml.dump(yaml_dict, f, default_flow_style=False, sort_keys=False)

        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            'boltz', 'predict', str(templated_yaml),
            '--accelerator', 'cpu',
            '--use_potentials',
            '--out_dir', str(output_dir),
            '--diffusion_samples', str(self.diffusion_samples),
            '--recycling_steps', str(self.recycling_steps)
        ]

        if self.cache_dir:
            cmd.extend(['--cache', self.cache_dir])

        if self.use_msa_server:
            cmd.append('--use_msa_server')

        logging.info(f"Running Boltz: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode != 0:
                logging.error(f"Boltz failed with return code {result.returncode}")
                if result.stderr:
                    logging.error(f"Error: {result.stderr[:500]}")
                return None

            prediction_files = list(output_dir.rglob("*_model_*.cif"))

            if prediction_files:
                best_model = sorted(prediction_files)[0]
                logging.info(f"Re-run completed: {best_model}")
                return best_model
            else:
                logging.error("No prediction files generated")
                return None

        except Exception as e:
            logging.error(f"Error running Boltz: {e}")
            return None


# ============================================================================
# Convenience Functions
# ============================================================================

def check_and_correct_prediction(
    prediction_file: Path,
    yaml_input: Path,
    expected_sequence: str,
    output_dir: Path,
    chain_id: str = 'B',
    expected_chirality: str = 'D',
    cache_dir: Optional[str] = None,
    max_iterations: int = 3
) -> Dict:
    """Convenience function to check and correct a single prediction."""
    corrector = ChiralityCorrector(cache_dir=cache_dir)

    return corrector.check_and_correct(
        prediction_file,
        yaml_input,
        expected_sequence,
        output_dir,
        chain_id,
        expected_chirality,
        max_iterations
    )


# ============================================================================
# Main / Example Usage
# ============================================================================

def main():
    """Example usage of chirality correction pipeline."""

    prediction_file = Path("path/to/prediction.cif")
    yaml_input = Path("path/to/input.yaml")
    expected_sequence = "DWWPLAFEALLR"
    output_dir = Path("chirality_correction_output")

    corrector = ChiralityCorrector(
        cache_dir=None,
        use_msa_server=True,
        diffusion_samples=5,
        recycling_steps=10
    )

    results = corrector.check_and_correct(
        prediction_file=prediction_file,
        yaml_input=yaml_input,
        expected_sequence=expected_sequence,
        output_dir=output_dir,
        chain_id='B',
        expected_chirality='D',
        max_iterations=3,
        violation_threshold=0.1
    )

    print("\n=== Chirality Correction Results ===")
    print(f"Number of iterations: {len(results['iterations'])}")
    print(f"Converged: {results['converged']}")
    print(f"Final violation rate: {results['final_violation_rate']:.1%}")

    for iter_result in results['iterations']:
        print(f"\nIteration {iter_result['iteration']}:")
        print(f"  CCD violation rate: {iter_result['ccd_violation_rate']:.1%}")
        print(f"  Geometric violation rate: {iter_result['geom_violation_rate']:.1%}")
        if 'num_corrections' in iter_result:
            print(f"  Corrections made: {iter_result['num_corrections']}")


if __name__ == "__main__":
    main()

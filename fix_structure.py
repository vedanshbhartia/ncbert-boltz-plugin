
import argparse
import sys
from pathlib import Path
import logging

# Ensure we can import from local lib
sys.path.append(str(Path(__file__).parent.parent))

from lib.structure_check import StructureValidator
from lib.optimize_structure import run_optimization, init_pyrosetta, check_pyrosetta

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)

def fix_structure(input_path: Path, output_path: Path, max_retries: int = 3) -> bool:
    """
    Attempt to fix a structure by running checking and optimization in a loop.
    """
    validator = StructureValidator(clash_threshold=0.6)
    
    current_input = input_path
    current_output = output_path
    
    # Check if initial structure is already good
    print(f"Checking {current_input}...")
    if validator.validate_file(current_input):
        print("Structure is already valid! No fixing needed.")
        # If output path is different, copy it over? Or just leave it.
        # Let's save a copy to output for consistency if requested
        if current_input != current_output:
            import shutil
            shutil.copy(current_input, current_output)
            print(f"Copied to {current_output}")
        return True

    # Initialize PyRosetta once
    if not check_pyrosetta():
        logging.error("PyRosetta is required but not found/installable.")
        return False
        
    init_pyrosetta()

    for i in range(1, max_retries + 1):
        print(f"\n--- Attempt {i}/{max_retries} ---")
        
        # Run Optimization
        # Use 'relax' mode as it is best for fixing clashes
        print("Running FastRelax...")
        run_optimization(current_input, current_output, mode='relax', constrain=True)
        
        # Check result
        print(f"Validating result {current_output}...")
        if validator.validate_file(current_output):
            print(f"\nSUCCESS: Structure fixed after {i} iteration(s)!")
            return True
            
        # Prepare for next iteration
        # Use the "fixed" (but still failing) structure as input for next round?
        # Sometimes iterative relax helps.
        current_input = current_output
        
    logging.error(f"Failed to fix structure after {max_retries} attempts.")
    return False

def main():
    parser = argparse.ArgumentParser(description="Automated Structure Fixer using PyRosetta")
    parser.add_argument("input_file", help="Input PDB/CIF file")
    parser.add_argument("--output", "-o", help="Output PDB file path")
    parser.add_argument("--retries", type=int, default=1, help="Max optimization retries (default 1)")
    
    args = parser.parse_args()
    
    input_path = Path(args.input_file)
    if not input_path.exists():
        logging.error(f"Input file not found: {input_path}")
        sys.exit(1)
        
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.with_name(f"{input_path.stem}_fixed{input_path.suffix}")
        
    success = fix_structure(input_path, output_path, max_retries=args.retries)
    
    if success:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()

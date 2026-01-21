import os
import sys
import subprocess
import argparse
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('batch_optimization.log')
    ]
)

def main():
    parser = argparse.ArgumentParser(description="Batch optimize protein structures.")
    parser.add_argument("target_dir", help="Directory containing structure files")
    parser.add_argument("--script_path", default="./optimize_structure.py", help="Path to optimize_structure.py")
    
    args = parser.parse_args()
    
    target_dir = Path(args.target_dir).resolve()
    script_path = Path(args.script_path).resolve()
    
    if not target_dir.exists():
        logging.error(f"Target directory does not exist: {target_dir}")
        sys.exit(1)
        
    if not script_path.exists():
        logging.error(f"Optimization script not found: {script_path}")
        sys.exit(1)

    cif_files = []
    for root, _, files in os.walk(target_dir):
        for file in files:
            if file.endswith(".cif") and "_optimized" not in file:
                cif_files.append(Path(root) / file)
    
    total_files = len(cif_files)
    logging.info(f"Found {total_files} CIF files to process.")
    
    for i, cif_file in enumerate(cif_files, 1):
        logging.info(f"Processing ({i}/{total_files}): {cif_file.name}")
        
        # Define output filename
        # Based on user request: "name the new structures appropriately"
        # optimize_structure.py defaults to {stem}_optimized{suffix} if output is not specified
        # We will let the script handle naming or pass it explicitly if needed.
        # Let's trust the default behavior of optimize_structure.py to keep it simple,
        # but check if we want to place it in the same dir (default).
        
        cmd = [
            sys.executable,
            str(script_path),
            str(cif_file),
            "--mode", "cartesian", # Use cartesian to fix geometry issues
            "--cyclic"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logging.info(f"Successfully processed {cif_file.name}")
            logging.info(result.stdout)
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to process {cif_file.name}")
            logging.error(f"Error output: {e.stderr}")
            # Continue to next file instead of crashing
            continue

if __name__ == "__main__":
    main()

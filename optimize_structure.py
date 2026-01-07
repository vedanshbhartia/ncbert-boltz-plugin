
import argparse
import sys
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)

def check_pyrosetta():
    """Check if pyrosetta is importable."""
    try:
        import pyrosetta
        return True
    except ImportError:
        return False

def init_pyrosetta():
    """Initialize PyRosetta with standard options."""
    import pyrosetta
    # mute core to reduce verbosity, enable standard options
    cmd_line_options = "-mute core.chemical -mute core.scoring -corrections::beta_nov16"
    pyrosetta.init(cmd_line_options)

def run_optimization(
    input_path: Path, 
    output_path: Path, 
    mode: str = 'relax',
    constrain: bool = True
):
    import pyrosetta
    from pyrosetta import rosetta

    print(f"Loading {input_path}...")
    
    # Load Pose
    try:
        pose = pyrosetta.pose_from_file(str(input_path))
    except Exception as e:
        logging.error(f"Failed to load structure: {e}")
        sys.exit(1)

    # Setup ScoreFunction
    # beta_nov16 is widely used for protein structure prediction/refinement
    # Use cartesian version if we are doing cartesian relax
    score_name = "beta_nov16_cart" if mode == 'cartesian' else "beta_nov16"
    sf = pyrosetta.create_score_function(score_name)

    # Setup Constraints (highly recommended to keep structure close to input)
    if constrain:
        print("Applying coordinate constraints...")
        # Add constraints to the pose
        # Constraint all heavy atoms to their starting coordinates with a harmonic potential
        
        # Easier way in PyRosetta: use CoordinateConstraintGenerator
        cg = rosetta.protocols.constraint_generator.CoordinateConstraintGenerator()
        cg.set_sd(0.5) # Standard deviation of the harmonic potential (lower = tighter)
        cg.set_bounded(False)
        
        # Apply constraints
        ac = rosetta.protocols.constraint_generator.AddConstraints()
        ac.add_generator(cg)
        ac.apply(pose)
        
        # Add constraint weight to score function
        sf.set_weight(rosetta.core.scoring.coordinate_constraint, 1.0)

    # IDEALIZATION STEP
    # If geometry is bad (bond lengths/angles), we MUST idealize first.
    # We do this before minimization/relax.
    print("Running IdealizeMover to fix bond lengths and angles...")
    idealizer = rosetta.protocols.idealize.IdealizeMover()
    idealizer.fast(True) # Fast mode
    idealizer.apply(pose)

    initial_energy = sf(pose)
    print(f"Initial Energy: {initial_energy:.2f}")

    if mode == 'minimize':
        print("Running MinMover (Gradient Descent)...")
        # MinMover: strictly minimizes energy in the local basin
        min_mover = rosetta.protocols.minimization_packing.MinMover()
        min_mover.movemap(rosetta.core.kinematics.MoveMap()) # Update all dofs
        # Customize movemap
        mm = rosetta.core.kinematics.MoveMap()
        mm.set_bb(True)
        mm.set_chi(True)
        mm.set_jump(True)
        min_mover.movemap(mm)
        min_mover.score_function(sf)
        min_mover.min_type('lbfgs_armijo_nonmonotone')
        min_mover.tolerance(0.001)
        
        min_mover.apply(pose)

    elif mode == 'relax':
        print("Running FastRelax (Iterative Packing/Minimization)...")
        # FastRelax: Iterative cycles of packing and minimization
        relax = pyrosetta.rosetta.protocols.relax.FastRelax()
        relax.set_scorefxn(sf)
        
        # If we added constraints, we typically want the relax to respect them
        if constrain:
            relax.constrain_relax_to_start_coords(True)
            relax.coord_constrain_sidechains(True)
        
        relax.apply(pose)

    elif mode == 'cartesian':
        print("Running Cartesian FastRelax...")
        relax = pyrosetta.rosetta.protocols.relax.FastRelax()
        relax.set_scorefxn(sf)
        relax.cartesian(True) # ENABLE CARTESIAN MINIMIZATION
        
        if constrain:
            relax.constrain_relax_to_start_coords(True)
            relax.coord_constrain_sidechains(True)
            
            # For Cartesian, we might need to be careful with prolines or constraints
            # but usually FastRelax handles it well.
        
        relax.apply(pose)

    final_energy = sf(pose)
    print(f"Final Energy:   {final_energy:.2f}")
    print(f"Delta Energy:   {final_energy - initial_energy:.2f}")

    print(f"Saving to {output_path}...")
    if output_path.suffix.lower() == '.cif':
        # PyRosetta mostly dumps PDB. Use Gemmi to convert if CIF is requested.
        temp_pdb = output_path.with_suffix('.pdb.tmp')
        pose.dump_pdb(str(temp_pdb))
        try:
            import gemmi
            st = gemmi.read_pdb(str(temp_pdb))
            st.setup_entities() # Helps with CIF structure
            st.make_mmcif_document().write_file(str(output_path))
            temp_pdb.unlink()
        except ImportError:
            logging.error("Gemmi not found. Cannot convert to CIF. Saving as PDB instead.")
            pose.dump_pdb(str(output_path.with_suffix('.pdb')))
    else:
        pose.dump_pdb(str(output_path))
    # Note: PyRosetta dump_pdb writes PDB format. If CIF is strictly required, 
    # we might need extra conversion, but PDB is usually fine for these sizes.

def main():
    parser = argparse.ArgumentParser(description="Optimize protein structure using PyRosetta")
    parser.add_argument("input_file", help="Input PDB/CIF file")
    parser.add_argument("--output", "-o", help="Output PDB file path")
    parser.add_argument("--mode", choices=['minimize', 'relax', 'cartesian'], default='cartesian',
                      help="Optimization mode: 'minimize', 'relax' (torsion), or 'cartesian' (fixes bond lengths)")
    parser.add_argument("--no-constraints", action="store_true", 
                      help="Disable coordinate constraints (warning: structure may drift largely)")
    
    args = parser.parse_args()
    
    input_path = Path(args.input_file)
    if not input_path.exists():
        logging.error(f"Input file not found: {input_path}")
        sys.exit(1)
        
    if args.output:
        output_path = Path(args.output)
    else:
        # Auto-name: input_minimized.pdb
        # Auto-name: input_optimized.ext
        output_path = input_path.with_name(f"{input_path.stem}_optimized{input_path.suffix}")

    if not check_pyrosetta():
        logging.error("PyRosetta not found!")
        logging.error("Please install it (e.g., via Conda with your license key).")
        logging.error("Command: conda install -c https://conda.rosetta.commons.org pyrosetta")
        sys.exit(1)

    init_pyrosetta()
    
    run_optimization(
        input_path, 
        output_path, 
        mode=args.mode, 
        constrain=not args.no_constraints
    )

    print("Done.")

if __name__ == "__main__":
    main()

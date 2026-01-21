
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

def run_cyclic_optimization(
    input_path: Path, 
    output_path: Path, 
    constrain: bool = True
):
    import pyrosetta
    from pyrosetta import rosetta
    from pyrosetta.rosetta.core.select.residue_selector import ChainSelector, NotResidueSelector, NeighborhoodResidueSelector, AndResidueSelector
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    
    print(f"Loading {input_path} for CYCLIC optimization...")
    try:
        pose = pyrosetta.pose_from_file(str(input_path))
    except Exception as e:
        logging.error(f"Failed to load structure: {e}")
        sys.exit(1)

    # 1. IDENTIFY CHAINS
    # Assume last chain is peptide
    n_chains = pose.conformation().num_chains()
    peptide_chain_id = n_chains
    print(f"Assuming Chain {peptide_chain_id} is the peptide.")
    
    # 2. SELECTORS
    pep_sel = ChainSelector(peptide_chain_id)
    prot_sel = NotResidueSelector(pep_sel)
    # Interface: Protein residues within 8A of peptide
    interface_sel = AndResidueSelector(prot_sel, NeighborhoodResidueSelector(pep_sel, 8.0, False))
    
    # Get residue indices (1-based)
    # Note: apply returns a vector1_bool, which is 1-indexed in internal access but mapped to list here
    # We iterate 1..total_residue
    total_res = pose.total_residue()
    pep_mask = pep_sel.apply(pose)
    interface_mask = interface_sel.apply(pose)
    
    pep_resi = [i for i in range(1, total_res+1) if pep_mask[i]]
    interface_resi = [i for i in range(1, total_res+1) if interface_mask[i]]
    
    if not pep_resi:
        logging.error("No residues found in the assumed peptide chain!")
        sys.exit(1)
        
    pep_start = pep_resi[0]
    pep_end = pep_resi[-1]
    print(f"Peptide residues: {pep_start}-{pep_end}")
    print(f"Interface residues: {len(interface_resi)} residues")

    # 3. SETUP CYCLIZATION
    print("Applying cyclic constraints (DeclareBond)...")
    # Remove termini variants
    try:
        rosetta.core.pose.remove_variant_type_from_pose_residue(pose, rosetta.core.chemical.LOWER_TERMINUS_VARIANT, pep_start)
        rosetta.core.pose.remove_variant_type_from_pose_residue(pose, rosetta.core.chemical.UPPER_TERMINUS_VARIANT, pep_end)
    except:
        logging.warning("Could not remove termini variants (might not exist).")

    # Declare Bond (Connect C of end to N of start)
    # Expert XML: res1=start atom1=N res2=end atom2=C
    print(f"Bonding {pep_start}:N to {pep_end}:C")
    try:
        pose.conformation().declare_chemical_bond(pep_start, "N", pep_end, "C")
    except Exception as e:
        logging.error(f"Manual bond declaration failed: {e}")
        # Proceeding might be fatal, but let's see.
        
    # Add Cutpoints for Relax (Constraint scoring)
    # Expert XML: Start -> CUTPOINT_UPPER, End -> CUTPOINT_LOWER
    # This matches the logic that the "gap" is between End(Lower) and Start(Upper) in the cyclic sense effectively?
    # Actually, if we connect End->Start. End is residues i. Start is i+1.
    # So End is Lower, Start is Upper.
    try:
        rosetta.core.pose.add_variant_type_to_pose_residue(pose, rosetta.core.chemical.CUTPOINT_UPPER, pep_start)
        rosetta.core.pose.add_variant_type_to_pose_residue(pose, rosetta.core.chemical.CUTPOINT_LOWER, pep_end)
    except Exception as e:
        logging.warning(f"Failed to add cutpoint variants: {e}")

    # 4. MOVEMAP
    mm = MoveMap()
    mm.set_bb(False); mm.set_chi(False); mm.set_jump(False)
    
    # Allow Peptide BB/Chi
    for r in pep_resi:
        mm.set_bb(r, True)
        mm.set_chi(r, True)
        
    # Allow Interface Chi
    for r in interface_resi:
        mm.set_chi(r, True)
        
    # 5. SCORE FUNCTION
    sf = pyrosetta.create_score_function("beta_nov16")
    sf.set_weight(rosetta.core.scoring.chainbreak, 100.0)
    sf.set_weight(rosetta.core.scoring.hbond_lr_bb, 6.0)
    sf.set_weight(rosetta.core.scoring.hbond_sr_bb, 6.0)
    
    # 6. RUN FASTRELAX
    print("Running FastRelax with Cyclic Constraints...")
    fr = FastRelax()
    fr.set_scorefxn(sf)
    fr.set_movemap(mm)
    if constrain:
        fr.constrain_relax_to_start_coords(True)
        fr.coord_constrain_sidechains(True)
    
    fr.apply(pose)
    
    # 7. CLEANUP & SAVE
    # Remove cutpoints
    try:
        rosetta.core.pose.remove_variant_type_from_pose_residue(pose, rosetta.core.chemical.CUTPOINT_UPPER, pep_start)
        rosetta.core.pose.remove_variant_type_from_pose_residue(pose, rosetta.core.chemical.CUTPOINT_LOWER, pep_end)
    except:
        pass
        
    final_energy = sf(pose)
    print(f"Final Score: {final_energy:.2f}")
    
    print(f"Saving to {output_path}...")
    if output_path.suffix.lower() == '.cif':
        pose.dump_pdb(str(output_path.with_suffix('.pdb')))
        # Conversion logic same as above if needed, but keeping it simple for now
    else:
        pose.dump_pdb(str(output_path))

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
    parser.add_argument("--cyclic", action="store_true", help="Enable cyclic peptide optimization (expert mode)")
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
    
    if args.cyclic:
        run_cyclic_optimization(
            input_path,
            output_path,
            constrain=not args.no_constraints
        )
    else:
        run_optimization(
            input_path, 
            output_path, 
            mode=args.mode, 
            constrain=not args.no_constraints
        )

    print("Done.")

if __name__ == "__main__":
    main()

"""PyRosetta-based cyclic peptide optimization.

Optimizes cyclic peptide-protein complexes by declaring N->C bonds and
applying appropriate constraints during FastRelax.
"""

import argparse
import logging
import sys
from pathlib import Path

PYROSETTA_INIT_FLAGS = "-mute core.chemical"
CONSTRAINT_SD = 0.5
INTERFACE_CUTOFF = 8.0  # angstroms

CHAINBREAK_WEIGHT = 20.0
HBOND_LR_WEIGHT = 6.0
HBOND_SR_WEIGHT = 6.0
CART_BONDED_WEIGHT = 1

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def init_pyrosetta(flags: str = PYROSETTA_INIT_FLAGS):
    import pyrosetta
    pyrosetta.init(flags)


def _manage_cyclic_variants(pose, pep_start: int, pep_end: int, add_cutpoints: bool = False):
    from pyrosetta import rosetta
    
    if not add_cutpoints:
        # Remove terminus variants
        try:
            rosetta.core.pose.remove_variant_type_from_pose_residue(
                pose, rosetta.core.chemical.LOWER_TERMINUS_VARIANT, pep_start
            )
            rosetta.core.pose.remove_variant_type_from_pose_residue(
                pose, rosetta.core.chemical.UPPER_TERMINUS_VARIANT, pep_end
            )
        except Exception:
            logger.warning("Couldn't remove terminus variants")
    else:
        # Add cutpoint variants for relax
        try:
            rosetta.core.pose.add_variant_type_to_pose_residue(
                pose, rosetta.core.chemical.CUTPOINT_UPPER, pep_start
            )
            rosetta.core.pose.add_variant_type_to_pose_residue(
                pose, rosetta.core.chemical.CUTPOINT_LOWER, pep_end
            )
        except Exception as e:
            logger.warning(f"Cutpoint variant setup failed: {e}")


def _cleanup_cutpoints(pose, pep_start: int, pep_end: int):
    """Remove cutpoint variants after optimization."""
    from pyrosetta import rosetta
    try:
        rosetta.core.pose.remove_variant_type_from_pose_residue(
            pose, rosetta.core.chemical.CUTPOINT_UPPER, pep_start
        )
        rosetta.core.pose.remove_variant_type_from_pose_residue(
            pose, rosetta.core.chemical.CUTPOINT_LOWER, pep_end
        )
    except Exception:
        pass


def _save_pose(pose, output_path: Path):
    """Save pose as both PDB and CIF formats."""
    # Always establish paths for both formats
    pdb_path = output_path.with_suffix('.pdb')
    cif_path = output_path.with_suffix('.cif')
    
    # 1. Always save PDB (PyRosetta native)
    logger.info(f"Saving PDB to {pdb_path}...")
    pose.dump_pdb(str(pdb_path))
    
    # 2. Try to save CIF using Gemmi (needs the PDB we just saved)
    logger.info(f"Saving CIF to {cif_path}...")
    try:
        import gemmi
        st = gemmi.read_pdb(str(pdb_path))
        st.setup_entities()
        st.make_mmcif_document().write_file(str(cif_path))
    except ImportError:
        logger.warning(f"Gemmi not installed. Could not generate {cif_path}")
    except Exception as e:
        logger.error(f"Failed to convert PDB to CIF: {e}")


def optimize_cyclic_peptide(input_path: Path, output_path: Path, constrain: bool = True):
    """    
    Assumes last chain is the peptide. Declares N->C bond, sets up cutpoints,
    and optimizes with appropriate constraints.
    """
    import pyrosetta
    from pyrosetta import rosetta
    from pyrosetta.rosetta.core.select.residue_selector import (
        ChainSelector, NotResidueSelector, NeighborhoodResidueSelector, AndResidueSelector
    )
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    
    logger.info(f"Loading {input_path} for cyclic optimization...")
    
    try:
        pose = pyrosetta.pose_from_file(str(input_path))
    except Exception as e:
        logger.error(f"Failed to load: {e}")
        sys.exit(1)

    # Assume last chain = peptide
    n_chains = pose.conformation().num_chains()
    peptide_chain = n_chains
    logger.info(f"{n_chains} chains detected, using chain {peptide_chain} as peptide")
    
    pep_sel = ChainSelector(peptide_chain)
    prot_sel = NotResidueSelector(pep_sel)
    interface_sel = AndResidueSelector(
        prot_sel, 
        NeighborhoodResidueSelector(pep_sel, INTERFACE_CUTOFF, False)
    )
    
    total_res = pose.total_residue()
    pep_mask = pep_sel.apply(pose)
    interface_mask = interface_sel.apply(pose)
    
    pep_resi = [i for i in range(1, total_res + 1) if pep_mask[i]]
    interface_resi = [i for i in range(1, total_res + 1) if interface_mask[i]]
    
    if not pep_resi:
        logger.error("No peptide residues found!")
        sys.exit(1)
        
    pep_start, pep_end = pep_resi[0], pep_resi[-1]
    logger.info(f"Peptide: {pep_start}-{pep_end}, Interface: {len(interface_resi)} residues")

    logger.info("Setting up cyclic bond...")
    # _manage_cyclic_variants(pose, pep_start, pep_end, add_cutpoints=False)
    
    # logger.info(f"Declaring bond: {pep_start}:N <-> {pep_end}:C")
    # try:
    #     pose.conformation().declare_chemical_bond(pep_start, "N", pep_end, "C")
    # except Exception as e:
    #     logger.error(f"Bond declaration failed: {e}")
        
    # _manage_cyclic_variants(pose, pep_start, pep_end, add_cutpoints=True)

    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    mm.set_jump(False)
    
    for r in pep_resi:
        mm.set_bb(r, True)
        mm.set_chi(r, True)
    for r in interface_resi:
        mm.set_chi(r, True)

    sf = pyrosetta.create_score_function("ref2015_cart")
    # sf.set_weight(rosetta.core.scoring.chainbreak, CHAINBREAK_WEIGHT)
    sf.set_weight(rosetta.core.scoring.hbond_lr_bb, HBOND_LR_WEIGHT)
    sf.set_weight(rosetta.core.scoring.hbond_sr_bb, HBOND_SR_WEIGHT)
    sf.set_weight(rosetta.core.scoring.cart_bonded, CART_BONDED_WEIGHT)

    sf.show(pose)
    
    logger.info("Running Cartesian FastRelax...")
    relax = FastRelax()
    relax.set_scorefxn(sf)
    relax.set_movemap(mm)
    relax.cartesian(True)
    relax.minimize_bond_angles(True)
    relax.minimize_bond_lengths(True)
    if constrain:
        relax.constrain_relax_to_start_coords(True)
        relax.coord_constrain_sidechains(True)
    
    relax.apply(pose)
    
    # _cleanup_cutpoints(pose, pep_start, pep_end)
    
    final_score = sf(pose)
    logger.info(f"Final score: {final_score:.2f}")
    sf.show(pose)
    
    _save_pose(pose, output_path)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_file", help="Input PDB/CIF file")
    parser.add_argument("--output", "-o", help="Output path (default: <input>_optimized.<ext>)")
    parser.add_argument("--no-constraints", action="store_true",help="Disable coordinate constraints")
    args = parser.parse_args()
    
    input_path = Path(args.input_file)
    if not input_path.exists():
        logger.error(f"File not found: {input_path}")
        sys.exit(1)
    
    output_path = Path(args.output) if args.output else \
                  input_path.with_name(f"{input_path.stem}_optimized{input_path.suffix}")
    
    init_pyrosetta()
    
    constrain = not args.no_constraints
    optimize_cyclic_peptide(input_path, output_path, constrain)

    logger.info("Done")


if __name__ == "__main__":
    main()


import argparse
import sys
import math
from pathlib import Path
from typing import List, Dict, Optional, Set
import gemmi
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)

class StructureValidator:
    def __init__(self, clash_threshold: float = 0.6):
        """
        Initialize validator.
        
        Args:
            clash_threshold: Tolerance factor for VDW radii sum. 
                             If distance < (r1 + r2) * recall_threshold, it's a clash.
                             Default 0.6 is a common strictness for severe clashes.
        """
        self.clash_threshold = clash_threshold
        # Approximate VDW radii (in Angstroms) for common elements
        self.vdw_radii = {
            'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
            'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75,
            'Br': 1.85, 'I': 1.98
        }
        self.default_radius = 1.70

    def get_radius(self, element: str) -> float:
        return self.vdw_radii.get(element, self.default_radius)

    def check_clashes(self, structure: gemmi.Structure) -> List[Dict]:
        """Check for steric clashes between non-bonded atoms using NeighborSearch."""
        clashes = []
        
        search_radius = 5.0
        
        # Ensure we have a valid cell for neighbor search
        # If no cell info, create a large dummy box to allow NeighborSearch to work
        if not structure.cell.is_crystal() or structure.cell.a < 10.0:
            min_x = min_y = min_z = float('inf')
            max_x = max_y = max_z = float('-inf')
            
            for chain in structure[0]:
                 for res in chain:
                     for atom in res:
                         pos = atom.pos
                         min_x = min(min_x, pos.x)
                         max_x = max(max_x, pos.x)
                         min_y = min(min_y, pos.y)
                         max_y = max(max_y, pos.y)
                         min_z = min(min_z, pos.z)
                         max_z = max(max_z, pos.z)

            extent = max(max_x - min_x, max_y - min_y, max_z - min_z) + 50.0
            box_size = max(1000.0, extent * 2)
            structure.cell = gemmi.UnitCell(box_size, box_size, box_size, 90, 90, 90)

        ns = gemmi.NeighborSearch(structure[0], structure.cell, search_radius).populate()
        
        # We need to iterate atoms to query neighbors
        seen_pairs = set()

        for chain in structure[0]:
            for res in chain:
                for atom in res:
                    if atom.element.name == "H": continue
                    
                    r1 = self.get_radius(atom.element.name)
                    
                    # Find neighbors
                    marks = ns.find_atoms(atom.pos, '\0', radius=search_radius)
                    
                    for mark in marks:
                        # mark is a NeighborSearch.Mark, has .atom, .pos, etc
                        other_atom = mark.to_cra(structure[0]).atom
                        
                        if other_atom.element.name == "H": continue
                        
                        # Avoid self-clash
                        if atom is other_atom:  # Identity check
                            continue
                            
                        # Avoid duplicates (A-B vs B-A)
                        # Use memory address for stable comparison
                        id1, id2 = id(atom), id(other_atom)
                        if id1 >= id2: 
                            continue
                            
                        # Use image distance (handling periodicity if applicable)
                        dist = atom.pos.dist(mark.pos)
                        
                        # 1. Check exclusions
                        res1 = res
                        res2 = mark.to_cra(structure[0]).residue
                        chain1 = chain
                        chain2 = mark.to_cra(structure[0]).chain

                        # Same residue
                        if res1.seqid == res2.seqid and chain1.name == chain2.name:
                            continue

                        # Peptide bond (i to i+1)
                        if chain1.name == chain2.name:
                            diff = res1.seqid.num - res2.seqid.num
                            if abs(diff) == 1:
                                names = {atom.name, other_atom.name}
                                if 'C' in names and 'N' in names:
                                    continue
                        
                        # Disulfide bonds
                        if atom.name == 'SG' and other_atom.name == 'SG':
                             if res1.name == 'CYS' and res2.name == 'CYS':
                                 if dist < 2.3:
                                     continue

                        # 2. Check overlap
                        r2 = self.get_radius(other_atom.element.name)
                        limit = (r1 + r2) * self.clash_threshold
                        
                        if dist < limit:
                             clashes.append({
                                'atom1': f"{chain1.name}/{res1.name}{res1.seqid}/{atom.name}",
                                'atom2': f"{chain2.name}/{res2.name}{res2.seqid}/{other_atom.name}",
                                'distance': dist,
                                'limit': limit,
                                'overlap': limit - dist
                            })

        return clashes

    def check_backbone_geometry(self, structure: gemmi.Structure) -> List[Dict]:
        """
        Check standard protein bond lengths and angles.
        """
        issues = []
        
        for model in structure:
            for chain in model:
                residues = list(chain)
                for i in range(len(residues) - 1):
                    res1 = residues[i]
                    res2 = residues[i+1]
                    
                    # 1. Peptide Bond Length (C-N)
                    c_atom = res1.find_atom('C', '*')
                    n_atom = res2.find_atom('N', '*')
                    
                    if c_atom and n_atom:
                        pos_c = c_atom.pos
                        pos_n = n_atom.pos
                        dist = pos_c.dist(pos_n)
                        
                        if dist > 1.5 or dist < 1.1:
                            issues.append({
                                'type': 'Bond Length',
                                'label': 'Peptide C-N',
                                'res1': f"{res1.name}{res1.seqid}",
                                'res2': f"{res2.name}{res2.seqid}",
                                'value': dist,
                                'expected': '1.33 ± 0.2'
                            })

                    # 2. Backbone Angle (N-CA-C)
                    # We need N, CA, C within the SAME residue (res1)
                    n_local = res1.find_atom('N', '*')
                    ca_local = res1.find_atom('CA', '*')
                    c_local = res1.find_atom('C', '*')
                    
                    if n_local and ca_local and c_local:
                        # Vectors: CA->N and CA->C
                        v1 = gemmi.Position(n_local.pos.x - ca_local.pos.x, 
                                          n_local.pos.y - ca_local.pos.y, 
                                          n_local.pos.z - ca_local.pos.z)
                        v2 = gemmi.Position(c_local.pos.x - ca_local.pos.x, 
                                          c_local.pos.y - ca_local.pos.y, 
                                          c_local.pos.z - ca_local.pos.z)
                                          
                        # Calculate angle
                        # gemmi.calculate_angle uses 3 positions: A-B-C returns angle at B
                        angle_rad = gemmi.calculate_angle(n_local.pos, ca_local.pos, c_local.pos)
                        angle_deg = math.degrees(angle_rad)
                        
                        # Expected ~111 degrees. Allow range 95-125 safely.
                        if angle_deg < 90.0 or angle_deg > 130.0:
                             issues.append({
                                'type': 'Bond Angle',
                                'label': 'N-CA-C',
                                'res1': f"{res1.name}{res1.seqid}",
                                'res2': "-",
                                'value': angle_deg,
                                'expected': '111 ± 20'
                            })

        return issues

    def validate_file(self, file_path: Path) -> bool:
        """Run validation on a file. Returns True if passed (no severe issues)."""
        print(f"Validating {file_path}...")
        try:
            st = gemmi.read_structure(str(file_path))
        except Exception as e:
            logging.error(f"Failed to read structure: {e}")
            return False
            
        clashes = self.check_clashes(st)
        geometry_issues = self.check_backbone_geometry(st)
        
        passed = True
        
        if clashes:
            print(f"  FOUND {len(clashes)} CLASHES:")
            for c in clashes[:10]: # Show top 10
                 print(f"    {c['atom1']} - {c['atom2']}: {c['distance']:.3f} A (Limit {c['limit']:.3f})")
            if len(clashes) > 10:
                print(f"    ... and {len(clashes) - 10} more.")
            passed = False
        else:
            print("  No severe clashes found.")
            
        if geometry_issues:
            print(f"  FOUND {len(geometry_issues)} GEOMETRY ISSUES:")
            for b in geometry_issues[:10]:
                val_str = f"{b['value']:.1f} deg" if 'Angle' in b['type'] else f"{b['value']:.3f} A"
                print(f"    {b['res1']} {b['label']} ({b['type']}): {val_str} (Expected {b['expected']})")
            if len(geometry_issues) > 10:
                print(f"    ... and {len(geometry_issues) - 10} more.")
            passed = False
        else:
            print("  Backbone geometry seems reasonable.")
            
        return passed

def main():
    parser = argparse.ArgumentParser(description="Check structure for clashes and geometry issues.")
    parser.add_argument("files", nargs='+', help="CIF/PDB files to check")
    parser.add_argument("--threshold", type=float, default=0.6, help="Clash threshold factor (default 0.6)")
    
    args = parser.parse_args()
    
    validator = StructureValidator(clash_threshold=args.threshold)
    
    issues_found = False
    for f in args.files:
        path = Path(f)
        if not path.exists():
            logging.error(f"File not found: {path}")
            continue
            
        if not validator.validate_file(path):
            issues_found = True
            
    if issues_found:
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Run Boltz-2 on Rosetta MMP8-cyclic peptide designs.
Handles D-amino acids via phi-angle detection and CCD modifications.
NOTE: Boltz affinity is only for small molecule ligands (<128 heavy atoms);
      cyclic peptides are too large, so affinity metrics are skipped.
      We report confidence_score, pLDDT, ptm, iptm, protein_iptm instead.

Source layout: the new Rosetta outputs in runs/explicit_*/<design>/ have
chain A = 18-residue cyclic peptide and chain B = MMP8 target. The Boltz YAML
keeps the project convention of chain A = MMP8 and chain B = peptide, so the
template chain_ids in the YAML are remapped: YAML A <- source B, YAML B <- source A.
"""
import argparse
import math
import json
import re
import subprocess
import sys
import yaml
from pathlib import Path

BASE = Path("/home/vedansh/projects/random/git-install")
DEFAULT_RUN_DIR = BASE / "lib" / "runs" / "explicit_20260506_014900"
DEFAULT_BOLTZ_OUT = BASE / "lib" / "boltz" / "runs"
DEFAULT_PEPTIDE_OUT = BASE / "lib" / "boltz" / "runs_peptide_only"
VENV_BOLTZ = BASE / "venv" / "bin" / "boltz"
PERTURB_RE = re.compile(r"(alpha_\d+_1a85_cyclic_\d+_Perturb\d+)")

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C",
}

# Rosetta 3-letter residue code -> (canonical L one-letter for the sequence
# string, CCD code to register as a Boltz modification). Boltz takes a
# canonical L sequence plus per-position modifications.
ROSETTA_NCAA_INFO = {
    # D-amino acids: Rosetta's 3-letter codes happen to be valid CCD codes too.
    "DAL": ("A", "DAL"),
    "DAR": ("R", "DAR"),
    "DSG": ("N", "DSG"),  # D-Asn (not seen in current outputs but kept for completeness)
    "DAS": ("D", "DAS"),
    "DCY": ("C", "DCY"),
    "DGN": ("Q", "DGN"),
    "DGL": ("E", "DGL"),
    "DHI": ("H", "DHI"),
    "DIL": ("I", "DIL"),
    "DLE": ("L", "DLE"),
    "DLY": ("K", "DLY"),
    "MED": ("M", "MED"),
    "DPN": ("F", "DPN"),
    "DPR": ("P", "DPR"),
    "DSN": ("S", "DSN"),  # D-Ser
    "DTH": ("T", "DTH"),
    "DTR": ("W", "DTR"),
    "DTY": ("Y", "DTY"),
    "DVA": ("V", "DVA"),
    # NCAAs that may appear in 1a85 designs:
    "HPR": ("P", "HYP"),  # trans-4-hydroxy-L-proline
    "AIB": ("A", "AIB"),  # alpha-aminoisobutyric acid (achiral; sequence position uses A)
}


def read_pdb_atoms(pdb_path, chain):
    """Return {resi: {resn, atomname: [x,y,z]}} for ATOM/HETATM records in chain.

    Rosetta writes D-amino acids and other NCAA residues as HETATM records, so
    both record types must be parsed to reconstruct the full peptide chain.
    """
    atoms = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[21] == chain:
                resi = int(line[22:26])
                name = line[12:16].strip()
                resn = line[17:20].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                if resi not in atoms:
                    atoms[resi] = {"resn": resn}
                atoms[resi][name] = [x, y, z]
    return atoms


def dihedral(p1, p2, p3, p4):
    def cross(a, b):
        return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    def dot(a, b):
        return sum(a[i]*b[i] for i in range(3))
    def norm(v):
        n = math.sqrt(dot(v, v))
        return [x/n for x in v]
    b1 = [p2[i]-p1[i] for i in range(3)]
    b2 = [p3[i]-p2[i] for i in range(3)]
    b3 = [p4[i]-p3[i] for i in range(3)]
    n1, n2 = cross(b1, b2), cross(b2, b3)
    m1 = cross(n1, norm(b2))
    return math.degrees(math.atan2(dot(m1, n2), dot(n1, n2)))


def extract_sequence_and_modifications(atoms):
    """Walk chain residues in order; return (one_letter_seq, modifications).

    Modifications come from Rosetta NCAA residue names (DAL, DLE, DPR, HPR,
    AIB, ...) registered in ROSETTA_NCAA_INFO. The sequence string uses the
    canonical L one-letter code for each modified position so Boltz sees a
    standard polymer with per-position modifications applied.
    """
    seq_chars: list[str] = []
    mods: list[dict] = []
    for position, resi in enumerate(sorted(atoms.keys()), start=1):
        resn = atoms[resi]["resn"]
        if resn in THREE_TO_ONE:
            seq_chars.append(THREE_TO_ONE[resn])
        elif resn in ROSETTA_NCAA_INFO:
            one_letter, ccd = ROSETTA_NCAA_INFO[resn]
            seq_chars.append(one_letter)
            mods.append({"position": position, "ccd": ccd, "from_resn": resn})
        else:
            seq_chars.append("X")
            print(f"  [warn] unknown residue {resn} at position {position}; sequence char = X, no modification emitted")
    return "".join(seq_chars), mods


def extract_sequence(atoms):
    """Backwards-compat shim: just the one-letter sequence."""
    return extract_sequence_and_modifications(atoms)[0]


STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

# Rosetta NCAA 3-letter -> canonical L 3-letter equivalent for template
# residue renaming. Boltz applies the actual D/NCAA chemistry via the
# modifications list in YAML, so the template just needs L names at the
# correct positions.
NCAA_TO_CANONICAL_L_3 = {
    "DAL": "ALA", "DAR": "ARG", "DSG": "ASN", "DAS": "ASP", "DCY": "CYS",
    "DGN": "GLN", "DGL": "GLU", "DHI": "HIS", "DIL": "ILE", "DLE": "LEU",
    "DLY": "LYS", "MED": "MET", "DPN": "PHE", "DPR": "PRO", "DSN": "SER",
    "DTH": "THR", "DTR": "TRP", "DTY": "TYR", "DVA": "VAL",
    "HPR": "PRO", "AIB": "ALA",
}


def make_clean_template_cif(pdb_path: Path, out_path: Path) -> None:
    """
    Convert Rosetta PDB to a clean mmCIF template for boltz.
    - Keeps every polymer residue (canonical + D-amino acids + NCAAs); these
      are renamed to their canonical L 3-letter so the template residue
      sequence matches the YAML polymer sequence boltz constructs from the
      `modifications` list.
    - Drops single-residue HETATM ions/cofactors (metals etc.) which boltz
      cannot consume as part of a polymer template.
    - Sets entity sequences explicitly so gemmi/boltz reads them correctly.
    - Writes mmCIF format (avoids SEQRES parsing issues with gemmi+PDB).
    """
    import gemmi

    st = gemmi.read_structure(str(pdb_path))
    model = st[0]

    st.entities.clear()

    for entity_id, ch_name in enumerate(["A", "B"], start=1):
        ch = model[ch_name]

        # Drop residues that are neither canonical L AAs nor recognised NCAAs.
        # These are typically metal ions or buffer molecules that came in as
        # single-residue HETATMs and have no polymer position.
        keep_set = STD_AA | set(NCAA_TO_CANONICAL_L_3)
        non_polymer_indices = [i for i, res in enumerate(ch) if res.name not in keep_set]
        for i in reversed(non_polymer_indices):
            del ch[i]

        # Rename NCAA residues to their canonical L name so the template
        # residue sequence matches what boltz expects from the YAML.
        for res in ch:
            if res.name in NCAA_TO_CANONICAL_L_3:
                res.name = NCAA_TO_CANONICAL_L_3[res.name]

        seq = [res.name for res in ch]

        ent = gemmi.Entity(str(entity_id))
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = gemmi.PolymerType.PeptideL
        ent.full_sequence = seq
        ent.subchains = [ch_name]  # maps label_asym_id -> entity in _struct_asym
        st.entities.append(ent)

        for res in ch:
            res.entity_id = str(entity_id)
            res.entity_type = gemmi.EntityType.Polymer
            res.subchain = ch_name  # sets label_asym_id in CIF

    st.connections.clear()  # cyclic bond conveyed via YAML cyclic:true; LINK records cause KeyError in boltz template parser
    doc = st.make_mmcif_document()
    doc.write_file(str(out_path))


def create_yaml(design_name, pdb_path, run_dir, peptide_only=False):
    """Build Boltz-2 YAML input and return path.

    Source PDB convention: chain A = peptide (18 residues), chain B = MMP8.
    Bound mode YAML convention (kept for parity with prior runs): id A = MMP8,
    id B = peptide. In peptide_only mode the MMP8 chain and templates are
    dropped and the peptide is folded de novo with msa: empty.
    """
    atoms_pep = read_pdb_atoms(pdb_path, "A")    # source chain A = peptide
    pep_seq, pep_mods = extract_sequence_and_modifications(atoms_pep)
    modifications = [{"position": m["position"], "ccd": m["ccd"]} for m in pep_mods]

    print(f"\n{'='*60}")
    print(f"Design: {design_name}  (mode: {'peptide-only' if peptide_only else 'bound'})")
    print(f"  Peptide (source chain A): {len(pep_seq)} residues  [{pep_seq}]")
    print(f"  Modifications detected ({len(pep_mods)}):")
    for m in pep_mods:
        print(f"    pos {m['position']:2d} {m['from_resn']:3s} -> ccd {m['ccd']}")

    run_dir.mkdir(parents=True, exist_ok=True)

    if peptide_only:
        boltz_input = {
            "version": 1,
            "sequences": [
                {
                    "protein": {
                        "id": "B",
                        "sequence": pep_seq,
                        "msa": "empty",
                        "cyclic": True,
                        "modifications": modifications,
                    }
                },
            ],
        }
    else:
        atoms_mmp8 = read_pdb_atoms(pdb_path, "B")   # source chain B = MMP8
        mmp8_seq, _ = extract_sequence_and_modifications(atoms_mmp8)
        print(f"  MMP8 (source chain B): {len(mmp8_seq)} residues")

        template_cif = run_dir / "template.cif"
        make_clean_template_cif(pdb_path, template_cif)

        boltz_input = {
            "version": 1,
            "sequences": [
                {"protein": {"id": "A", "sequence": mmp8_seq}},
                {
                    "protein": {
                        "id": "B",
                        "sequence": pep_seq,
                        "cyclic": True,
                        "modifications": modifications,
                    }
                },
            ],
            "templates": [
                # Source CIF chain B = MMP8, chain A = peptide. YAML chain A = MMP8,
                # chain B = peptide. So template_id is the swapped source chain.
                {
                    "cif": str(template_cif),
                    "chain_id": "A",
                    "template_id": "B",
                    "force": True,
                    "threshold": 2.0,
                },
                {
                    "cif": str(template_cif),
                    "chain_id": "B",
                    "template_id": "A",
                },
            ],
        }
    yaml_path = run_dir / "input.yaml"
    with open(yaml_path, "w") as f:
        yaml.dump(boltz_input, f, default_flow_style=False, sort_keys=False)
    print(f"  YAML written: {yaml_path}")
    return yaml_path


def run_boltz(yaml_path, out_dir, n_samples=3):
    """Run boltz predict; raise if predictions are missing.

    Boltz exits 0 even when an input is skipped (e.g. when the colabfold MSA
    server returns "no space left on device"). We post-check that the
    predictions directory exists; if not, surface the tail of stderr/stdout
    as a RuntimeError so the caller can mark the design as failed instead of
    silently moving on.
    """
    cmd = [
        str(VENV_BOLTZ), "predict", str(yaml_path),
        "--out_dir", str(out_dir),
        "--use_msa_server",
        "--use_potentials",
        "--diffusion_samples", str(n_samples),
        "--recycling_steps", "3",
        "--output_format", "pdb",
        "--accelerator", "cpu",
        "--override",
    ]
    print(f"\nRunning: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("STDERR:", result.stderr[-4000:])
        raise RuntimeError(f"boltz failed (exit {result.returncode})")

    pred_dir = Path(out_dir) / "boltz_results_input" / "predictions" / "input"
    if not pred_dir.is_dir():
        # Boltz exited 0 but produced nothing. Most often this is the colabfold
        # MSA server returning a transient error.
        diag = []
        for stream_name, content in (("STDOUT", result.stdout), ("STDERR", result.stderr)):
            tail = content.splitlines()[-30:]
            diag.append(f"--- {stream_name} (last 30 lines) ---\n" + "\n".join(tail))
        print("\n".join(diag))
        raise RuntimeError(
            f"boltz exited 0 but produced no predictions at {pred_dir}; "
            "check the dump above (commonly an MSA-server outage)."
        )
    return result


def parse_results(out_dir, design_name):
    """Parse confidence JSON files and report metrics."""
    pred_dir = out_dir / "boltz_results_input" / "predictions" / "input"
    if not pred_dir.exists():
        print(f"  No predictions directory found at {pred_dir}")
        return

    conf_files = sorted(pred_dir.glob("confidence_input_model_*.json"))
    if not conf_files:
        print("  No confidence JSON files found.")
        return

    print(f"\n{'='*60}")
    print(f"Results for {design_name}:")
    print(f"{'='*60}")

    best_score = -1
    best_model = None

    def fmt(val):
        return f"{val:.4f}" if isinstance(val, (int, float)) else "N/A"

    for cf in conf_files:
        with open(cf) as f:
            conf = json.load(f)
        model_name = cf.stem.replace("confidence_", "")
        score = conf.get("confidence_score", 0)
        if score > best_score:
            best_score = score
            best_model = model_name
        print(f"\n  Model: {model_name}")
        print(f"    confidence_score : {fmt(conf.get('confidence_score'))}")
        print(f"    complex_plddt    : {fmt(conf.get('complex_plddt'))}")
        print(f"    ptm              : {fmt(conf.get('ptm'))}")
        if "iptm" in conf and conf["iptm"]:
            print(f"    iptm             : {fmt(conf.get('iptm'))}")
            print(f"    protein_iptm     : {fmt(conf.get('protein_iptm'))}  <- peptide interface")
        print(f"    complex_iplddt   : {fmt(conf.get('complex_iplddt'))}")

    pdb_files = sorted(pred_dir.glob("input_model_*.pdb"))
    print(f"\n  Best model: {best_model}  (confidence_score={best_score:.4f})")
    print(f"  Structure files: {[p.name for p in pdb_files]}")

    interpret_results(best_score,
                      conf.get("complex_plddt", 0),
                      conf.get("protein_iptm"))


def interpret_results(confidence, plddt, protein_iptm):
    print("\n  --- INTERPRETATION ---")
    # Confidence score
    if confidence >= 0.8:
        print(f"  Overall quality : EXCELLENT ({confidence:.2f}) - very well folded complex")
    elif confidence >= 0.6:
        print(f"  Overall quality : GOOD ({confidence:.2f}) - likely reasonable complex")
    elif confidence >= 0.4:
        print(f"  Overall quality : MODERATE ({confidence:.2f}) - uncertain structure")
    else:
        print(f"  Overall quality : LOW ({confidence:.2f}) - poor prediction, likely mis-folded")

    # pLDDT
    if plddt >= 0.8:
        print(f"  pLDDT           : HIGH ({plddt:.2f}) - confident residue-level prediction")
    elif plddt >= 0.6:
        print(f"  pLDDT           : MEDIUM ({plddt:.2f}) - moderate per-residue confidence")
    else:
        print(f"  pLDDT           : LOW ({plddt:.2f}) - uncertain residue positions")

    if protein_iptm:
        # protein_iptm - binding interface (bound mode only)
        if protein_iptm >= 0.6:
            print(f"  Binding pose    : CONFIDENT ({protein_iptm:.2f}) - interface well-predicted")
        elif protein_iptm >= 0.4:
            print(f"  Binding pose    : MODERATE ({protein_iptm:.2f}) - interface uncertain")
        else:
            print(f"  Binding pose    : LOW ({protein_iptm:.2f}) - poor interface prediction")
        print("\n  NOTE: Boltz-2 affinity (IC50) requires a small-molecule ligand chain (<128 heavy")
        print("  atoms). An 18-residue cyclic peptide exceeds this limit, so affinity values are")
        print("  not available. Use protein_iptm and confidence_score as binding quality proxies.")
    else:
        print("  (peptide-only run: no interface metrics; ptm/plddt reflect intrinsic foldability)")


def find_design_pdb(rosetta_run_dir: Path, design: str) -> Path | None:
    """Locate the Rosetta output PDB for a given design id.

    The Rosetta job dir is <run>/<design>/, and the output PDB is named after
    the per-job perturbed scaffold (e.g. alpha_..._Perturb14_0001.pdb), not
    the design id itself. We pick the first non-input_scaffold PDB.
    """
    job_dir = rosetta_run_dir / design
    if not job_dir.is_dir():
        return None
    candidates = sorted(p for p in job_dir.glob("*.pdb") if p.name != "input_scaffold.pdb")
    return candidates[0] if candidates else None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rosetta-run-dir",
        type=Path,
        default=DEFAULT_RUN_DIR,
        help=f"Rosetta batch run directory containing per-design output PDBs "
             f"(default: {DEFAULT_RUN_DIR})",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=DEFAULT_BOLTZ_OUT,
        help=f"Root directory under which per-design boltz workspaces are created "
             f"(default: {DEFAULT_BOLTZ_OUT})",
    )
    parser.add_argument(
        "--designs",
        nargs="*",
        default=None,
        help="Subset of design ids to run. Defaults to every design subdirectory "
             "found under --rosetta-run-dir.",
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=3,
        help="Number of diffusion samples per design (default: 3).",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip designs whose boltz_out/.../predictions/input dir already exists.",
    )
    parser.add_argument(
        "--peptide-only",
        action="store_true",
        help="Drop the MMP8 chain and templates; fold the cyclic peptide alone "
             "(msa: empty). When set without --out-root, output goes to "
             f"{DEFAULT_PEPTIDE_OUT} instead of {DEFAULT_BOLTZ_OUT}.",
    )
    args = parser.parse_args()
    if args.peptide_only and args.out_root == DEFAULT_BOLTZ_OUT:
        args.out_root = DEFAULT_PEPTIDE_OUT
    return args


def main():
    args = parse_args()

    if args.designs:
        designs = args.designs
    else:
        designs = sorted(
            p.name for p in args.rosetta_run_dir.iterdir()
            if p.is_dir() and p.name not in {"reference"} and not p.name.startswith(".")
        )

    print(f"Rosetta run dir : {args.rosetta_run_dir}")
    print(f"Boltz output    : {args.out_root}")
    print(f"Designs to run  : {len(designs)}")

    completed = 0
    for design in designs:
        pdb_path = find_design_pdb(args.rosetta_run_dir, design)
        if pdb_path is None:
            print(f"[skip] {design}: no PDB found under {args.rosetta_run_dir / design}")
            continue
        run_dir = args.out_root / design
        out_dir = run_dir / "boltz_out"
        if args.skip_existing and (out_dir / "boltz_results_input" / "predictions" / "input").is_dir():
            print(f"[skip] {design}: boltz output already exists at {out_dir}")
            continue

        yaml_path = create_yaml(design, pdb_path, run_dir, peptide_only=args.peptide_only)
        try:
            run_boltz(yaml_path, out_dir, n_samples=args.n_samples)
        except RuntimeError as e:
            print(f"[fail] {design}: {e}")
            continue
        parse_results(out_dir, design)
        completed += 1

    print(f"\nDone. Completed {completed}/{len(designs)} designs.")


if __name__ == "__main__":
    main()

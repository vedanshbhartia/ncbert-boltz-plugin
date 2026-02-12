#!/usr/bin/env python3
import argparse
import json
import math
import random
import subprocess
import sys
from typing import Optional
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

try:
    import gemmi
except ImportError as exc:
    raise SystemExit(
        "Gemmi is required for this script. Install it with: pip install gemmi"
    ) from exc

try:
    from main import ModifiedSequenceParser, CCD_TO_BASE_AA, ONE_TO_THREE_AA
except Exception:
    import importlib.util

    main_path = REPO_ROOT / "main.py"
    if not main_path.exists():
        raise SystemExit("main.py not found in repo root; cannot load sequence parser")
    spec = importlib.util.spec_from_file_location("main", main_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    ModifiedSequenceParser = module.ModifiedSequenceParser
    CCD_TO_BASE_AA = module.CCD_TO_BASE_AA
    ONE_TO_THREE_AA = module.ONE_TO_THREE_AA


def die(msg: str) -> None:
    raise SystemExit(msg)


def resolve_rosetta_main(path_str: str) -> Path:
    root = Path(path_str).expanduser().resolve()
    if (root / "source").exists() and (root / "database").exists():
        return root
    if (root / "main").exists():
        main = root / "main"
        if (main / "source").exists() and (main / "database").exists():
            return main
    die(f"Rosetta main directory not found under: {root}")


def find_binary(bin_dir: Path, base: str) -> Path:
    candidates = [
        f"{base}.static.linuxgccrelease",
        f"{base}.linuxgccrelease",
        f"{base}.static.macosclangrelease",
        f"{base}.macosclangrelease",
    ]
    for name in candidates:
        path = bin_dir / name
        if path.exists():
            return path
    die(f"Could not find {base} binary in {bin_dir}")


def parse_sequence(seq: str):
    parser = ModifiedSequenceParser()
    base_seq, mods = parser.parse_sequence(seq)
    return base_seq, mods


def d_residue_name(one_letter: str) -> Optional[str]:
    if one_letter == "G":
        return None
    base_three = ONE_TO_THREE_AA.get(one_letter)
    if not base_three or base_three in {"UNK", "ASX", "GLX"}:
        return None
    return f"D{base_three}"


def d_mutations(mods):
    muts = []
    for mod in mods:
        ccd = mod.get("ccd", "")
        if not ccd.startswith("D"):
            continue
        base = CCD_TO_BASE_AA.get(ccd)
        if not base:
            continue
        dname = d_residue_name(base)
        if not dname:
            continue
        muts.append((mod["position"], dname, ccd))
    return muts


def read_sequence_from_file(path: Path) -> str:
    seq = ""
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq = line
            break
    if not seq:
        die(f"No sequence found in {path}")
    return seq


def build_peptide(build_bin: Path, db_dir: Path, fasta_path: Path, out_pdb: Path):
    cmd = [
        str(build_bin),
        "-database",
        str(db_dir),
        "-in:file:fasta",
        str(fasta_path),
        "-out:file:o",
        str(out_pdb),
    ]
    subprocess.run(cmd, check=True)


def clone_chain(chain: gemmi.Chain, new_id: str) -> gemmi.Chain:
    new_chain = chain.clone()
    new_chain.name = new_id
    return new_chain


def extract_chain(structure: gemmi.Structure, chain_id: str) -> gemmi.Chain:
    model = structure[0]
    for chain in model:
        if chain.name == chain_id:
            return chain
    die(f"Chain {chain_id} not found in structure")


def get_residue_by_num(chain: gemmi.Chain, resnum: int) -> Optional[gemmi.Residue]:
    for res in chain:
        if res.seqid.num == resnum:
            return res
    return None




def find_atom(res: gemmi.Residue, name: str) -> Optional[gemmi.Atom]:
    try:
        atom = res.find_atom(name, '\0')
    except TypeError:
        atom = None
    if atom:
        return atom
    for a in res:
        if a.name == name:
            return a
    return None


def residue_atoms(res: gemmi.Residue):
    return [atom for atom in res if atom.element.name != "H"]


def chain_centroid(chain: gemmi.Chain) -> gemmi.Position:
    atoms = []
    for res in chain:
        if res.is_water():
            continue
        atoms.extend(residue_atoms(res))
    if not atoms:
        die("No atoms found for peptide chain")
    x = sum(a.pos.x for a in atoms) / len(atoms)
    y = sum(a.pos.y for a in atoms) / len(atoms)
    z = sum(a.pos.z for a in atoms) / len(atoms)
    return gemmi.Position(x, y, z)


def random_rotation(rng: random.Random):
    u1 = rng.random()
    u2 = rng.random()
    u3 = rng.random()
    q1 = math.sqrt(1 - u1) * math.sin(2 * math.pi * u2)
    q2 = math.sqrt(1 - u1) * math.cos(2 * math.pi * u2)
    q3 = math.sqrt(u1) * math.sin(2 * math.pi * u3)
    q4 = math.sqrt(u1) * math.cos(2 * math.pi * u3)
    # Rotation matrix
    return [
        [1 - 2 * (q3 * q3 + q4 * q4), 2 * (q2 * q3 - q1 * q4), 2 * (q2 * q4 + q1 * q3)],
        [2 * (q2 * q3 + q1 * q4), 1 - 2 * (q2 * q2 + q4 * q4), 2 * (q3 * q4 - q1 * q2)],
        [2 * (q2 * q4 - q1 * q3), 2 * (q3 * q4 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3)],
    ]


def apply_rotation(pos: gemmi.Position, rot):
    x = pos.x * rot[0][0] + pos.y * rot[0][1] + pos.z * rot[0][2]
    y = pos.x * rot[1][0] + pos.y * rot[1][1] + pos.z * rot[1][2]
    z = pos.x * rot[2][0] + pos.y * rot[2][1] + pos.z * rot[2][2]
    return gemmi.Position(x, y, z)


def random_offset(rng: random.Random, radius: float) -> gemmi.Position:
    # Uniform within a sphere
    u = rng.random()
    v = rng.random()
    theta = 2 * math.pi * u
    phi = math.acos(2 * v - 1)
    r = radius * (rng.random() ** (1.0 / 3.0))
    x = r * math.sin(phi) * math.cos(theta)
    y = r * math.sin(phi) * math.sin(theta)
    z = r * math.cos(phi)
    return gemmi.Position(x, y, z)


def place_peptide(chain: gemmi.Chain, center: gemmi.Position, rng: random.Random, radius: float):
    centroid = chain_centroid(chain)
    rot = random_rotation(rng)
    offset = random_offset(rng, radius)
    for res in chain:
        if res.is_water():
            continue
        for atom in res:
            rel = gemmi.Position(atom.pos.x - centroid.x, atom.pos.y - centroid.y, atom.pos.z - centroid.z)
            rotated = apply_rotation(rel, rot)
            atom.pos = gemmi.Position(
                rotated.x + center.x + offset.x,
                rotated.y + center.y + offset.y,
                rotated.z + center.z + offset.z,
            )


def compute_pocket_from_complex(complex_path: Path, receptor_chain: str, peptide_chain: str, cutoff: float):
    st = gemmi.read_structure(str(complex_path))
    model = st[0]
    rec = extract_chain(st, receptor_chain)
    pep = extract_chain(st, peptide_chain)

    pep_atoms = []
    for res in pep:
        if res.is_water():
            continue
        pep_atoms.extend(residue_atoms(res))
    if not pep_atoms:
        die("No peptide atoms found in pocket template")

    pocket = []
    for res in rec:
        if res.is_water():
            continue
        min_dist = None
        for atom in residue_atoms(res):
            for patom in pep_atoms:
                dist = atom.pos.dist(patom.pos)
                if min_dist is None or dist < min_dist:
                    min_dist = dist
        if min_dist is not None and min_dist <= cutoff:
            pocket.append({
                "chain": receptor_chain,
                "resnum": res.seqid.num,
                "resname": res.name,
                "min_dist": float(min_dist),
            })

    if not pocket:
        die("No pocket residues found. Try increasing cutoff.")

    # Pocket center from CA atoms
    coords = []
    for entry in pocket:
        res = get_residue_by_num(rec, entry["resnum"])
        if res is None:
            continue
        ca = find_atom(res, "CA")
        if ca:
            coords.append(ca.pos)
        else:
            for atom in residue_atoms(res):
                coords.append(atom.pos)
    cx = sum(p.x for p in coords) / len(coords)
    cy = sum(p.y for p in coords) / len(coords)
    cz = sum(p.z for p in coords) / len(coords)

    return {
        "receptor_chain": receptor_chain,
        "peptide_chain": peptide_chain,
        "cutoff": cutoff,
        "center": [cx, cy, cz],
        "residues": pocket,
    }


def write_pocket_json(info, path: Path):
    with path.open("w") as f:
        json.dump(info, f, indent=2)


def load_pocket_json(path: Path):
    with path.open() as f:
        return json.load(f)


def format_resnum_chain(resnum: int, chain: str) -> str:
    return f"{resnum}{chain}"


def constraint_lines_from_pose(pose_path: Path, pocket_info: dict, pep_chain: str, anchor_index: int, top_n: int, sd: float):
    st = gemmi.read_structure(str(pose_path))
    model = st[0]
    pep = extract_chain(st, pep_chain)

    pep_residues = [res for res in pep if not res.is_water()]
    if not pep_residues:
        die("No peptide residues found in pose")
    if anchor_index < 1 or anchor_index > len(pep_residues):
        die(f"Anchor residue {anchor_index} is out of range for peptide length {len(pep_residues)}")

    anchor_res = pep_residues[anchor_index - 1]
    anchor_atom = find_atom(anchor_res, "CA") or next(iter(anchor_res))

    rec_chain = pocket_info["receptor_chain"]
    rec = extract_chain(st, rec_chain)

    distances = []
    for entry in pocket_info["residues"]:
        resnum = entry["resnum"]
        res = get_residue_by_num(rec, resnum)
        if res is None:
            continue
        ca = find_atom(res, "CA") or next(iter(res))
        dist = anchor_atom.pos.dist(ca.pos)
        distances.append((dist, resnum))

    distances.sort(key=lambda x: x[0])
    distances = distances[:top_n]

    lines = []
    for dist, resnum in distances:
        lines.append(
            "AtomPair CA {pep} CA {rec} HARMONIC {d:.3f} {sd:.3f}".format(
                pep=format_resnum_chain(anchor_res.seqid.num, pep_chain),
                rec=format_resnum_chain(resnum, rec_chain),
                d=dist,
                sd=sd,
            )
        )
    return lines


def write_constraints(lines, path: Path):
    with path.open("w") as f:
        for line in lines:
            f.write(line + "\n")


def write_mutate_xml(mutations, out_path: Path):
    parts = [
        "<ROSETTASCRIPTS>",
        "  <SCOREFXNS>",
        "    <ScoreFunction name=\"ref2015\" weights=\"ref2015\" />",
        "  </SCOREFXNS>",
        "  <MOVERS>",
    ]
    for i, (pos, resname, chain) in enumerate(mutations, start=1):
        parts.append(
            f"    <MutateResidue name=\"mut_{i}\" target=\"{pos}{chain}\" new_res=\"{resname}\" preserve_atom_coords=\"true\" />"
        )
    parts.append("  </MOVERS>")
    parts.append("  <PROTOCOLS>")
    for i in range(1, len(mutations) + 1):
        parts.append(f"    <Add mover=\"mut_{i}\" />")
    parts.append("  </PROTOCOLS>")
    parts.append("  <OUTPUT scorefxn=\"ref2015\" />")
    parts.append("</ROSETTASCRIPTS>")
    out_path.write_text("\n".join(parts) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Rosetta-only de novo peptide docking prep.")
    parser.add_argument("--rosetta-dir", required=True, help="Path to Rosetta main directory")
    parser.add_argument("--receptor", required=True, help="Receptor PDB path (MMP8 backbone)")
    parser.add_argument("--receptor-chain", default="A", help="Receptor chain ID")
    parser.add_argument("--peptide-chain", default="B", help="Peptide chain ID")
    parser.add_argument("--sequence", help="Peptide sequence (bracket mods allowed)")
    parser.add_argument("--sequence-file", help="FASTA or plain sequence file")
    parser.add_argument("--pocket-from", help="Complex (PDB/CIF) to derive pocket residues")
    parser.add_argument("--pocket-json", help="Precomputed pocket JSON")
    parser.add_argument("--pocket-cutoff", type=float, default=6.0, help="Distance cutoff for pocket residues")
    parser.add_argument("--start-radius", type=float, default=5.0, help="Random placement radius (Angstrom)")
    parser.add_argument("--nstart", type=int, default=10, help="Number of randomized starting poses")
    parser.add_argument("--anchor", type=int, help="Anchor peptide residue index (1-based). Defaults to middle residue.")
    parser.add_argument("--constraints-top-n", type=int, default=10, help="Number of pocket residues to constrain")
    parser.add_argument("--constraints-sd", type=float, default=2.0, help="Std dev for harmonic constraints")
    parser.add_argument("--outdir", default="rosetta_de_novo", help="Output directory")
    parser.add_argument("--seed", type=int, default=1234, help="Random seed")
    parser.add_argument("--emit-run-script", action="store_true", help="Write run_rosetta.sh")

    args = parser.parse_args()

    if not args.sequence and not args.sequence_file:
        die("Provide --sequence or --sequence-file")

    rosetta_main = resolve_rosetta_main(args.rosetta_dir)
    bin_dir = rosetta_main / "source" / "bin"
    db_dir = rosetta_main / "database"

    build_bin = find_binary(bin_dir, "BuildPeptide")
    flexpep_bin = find_binary(bin_dir, "FlexPepDocking")
    rs_bin = find_binary(bin_dir, "rosetta_scripts")

    receptor_path = Path(args.receptor).expanduser().resolve()
    if not receptor_path.exists():
        die(f"Receptor PDB not found: {receptor_path}")

    if args.sequence:
        seq_input = args.sequence.strip()
    else:
        seq_input = read_sequence_from_file(Path(args.sequence_file))

    base_seq, mods = parse_sequence(seq_input)
    pep_len = len(base_seq)

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    starts_dir = outdir / "starts"
    starts_dir.mkdir(exist_ok=True)
    constraints_dir = outdir / "constraints"
    constraints_dir.mkdir(exist_ok=True)

    fasta_path = outdir / "peptide.fasta"
    fasta_path.write_text(f">peptide\n{base_seq}\n")

    peptide_raw = outdir / "peptide_raw.pdb"
    build_peptide(build_bin, db_dir, fasta_path, peptide_raw)

    # Rename peptide chain
    pep_st = gemmi.read_structure(str(peptide_raw))
    pep_chain = pep_st[0][0]
    pep_chain.name = args.peptide_chain
    # Renumber residues starting at 1
    for i, res in enumerate(pep_chain, start=1):
        res.seqid.num = i
    peptide_pdb = outdir / "peptide.pdb"
    pep_st.write_pdb(str(peptide_pdb))

    pocket_info = None
    if args.pocket_json:
        pocket_info = load_pocket_json(Path(args.pocket_json))
    elif args.pocket_from:
        pocket_info = compute_pocket_from_complex(
            Path(args.pocket_from), args.receptor_chain, args.peptide_chain, args.pocket_cutoff
        )
        write_pocket_json(pocket_info, outdir / "pocket.json")

    if not pocket_info:
        die("Provide --pocket-from or --pocket-json to define the pocket center.")

    center = gemmi.Position(*pocket_info["center"])
    rng = random.Random(args.seed)

    # Load receptor chain
    rec_st = gemmi.read_structure(str(receptor_path))
    rec_chain = extract_chain(rec_st, args.receptor_chain)

    for i in range(1, args.nstart + 1):
        # Clone chains
        pep_clone = clone_chain(pep_chain, args.peptide_chain)
        place_peptide(pep_clone, center, rng, args.start_radius)

        # Load a fresh receptor structure and keep only the requested chain
        st = gemmi.read_structure(str(receptor_path))
        model = st[0]
        for ch in list(model):
            if ch.name != args.receptor_chain:
                model.remove_chain(ch.name)
        model.add_chain(pep_clone)

        start_path = starts_dir / f"start_{i:04d}.pdb"
        st.write_pdb(str(start_path))

        anchor = args.anchor or ((pep_len + 1) // 2)
        lines = constraint_lines_from_pose(
            start_path,
            pocket_info,
            args.peptide_chain,
            anchor,
            args.constraints_top_n,
            args.constraints_sd,
        )
        cst_path = constraints_dir / f"constraints_{i:04d}.cst"
        write_constraints(lines, cst_path)

    # D-residue mutations
    d_muts = d_mutations(mods)
    if d_muts:
        # Prepare mutation script
        mutate_xml = outdir / "mutate_peptide.xml"
        # Include chain for each mutation (peptide chain)
        muts = [(pos, resname, args.peptide_chain) for pos, resname, _ in d_muts]
        write_mutate_xml(muts, mutate_xml)

    # Emit run script
    if args.emit_run_script:
        run_path = outdir / "run_rosetta.sh"
        mmp8_len = len([res for res in rec_chain if not res.is_water()])
        pep_start = mmp8_len + 1
        pep_end = mmp8_len + pep_len
        relax_script = REPO_ROOT / "relax_script.xml"

        lines = [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "",
            f"ROSETTA_MAIN=\"{rosetta_main}\"",
            f"ROSETTA_DB=\"{db_dir}\"",
            f"FLEXPEP=\"{flexpep_bin}\"",
            f"ROSETTA_SCRIPTS=\"{rs_bin}\"",
            "",
            f"RELAX_XML=\"{relax_script}\"",
            f"MMP8_END={mmp8_len}",
            f"PEP_START={pep_start}",
            f"PEP_END={pep_end}",
            "",
            "mkdir -p mutated prepack refine relax",
            "",
        ]

        if d_muts:
            lines.append("# Apply D-residue mutations to each start pose")
            lines.append("for f in starts/start_*.pdb; do")
            lines.append("  name=$(basename \"$f\" .pdb)")
            lines.append("  $ROSETTA_SCRIPTS -database \"$ROSETTA_DB\" -parser:protocol mutate_peptide.xml -s \"$f\" -out:path:all mutated -out:suffix _mut")
            lines.append("done")
            lines.append("")

        lines.append("# Prepack and refine (adjust -nstruct as needed)")
        lines.append("for f in starts/start_*.pdb; do")
        lines.append("  name=$(basename \"$f\" .pdb)")
        lines.append("  cst=constraints/constraints_${name#start_}.cst")
        lines.append("  input=$f")
        if d_muts:
            lines.append("  if [ -f mutated/${name}_mut_0001.pdb ]; then input=mutated/${name}_mut_0001.pdb; fi")
        lines.append("  input_base=$(basename \"$input\" .pdb)")
        lines.append("  $FLEXPEP -database \"$ROSETTA_DB\" -s \"$input\" -flexPepDocking:flexpep_prepack -out:path:all prepack -out:suffix _prepack")
        lines.append("  prepack=prepack/${input_base}_prepack_0001.pdb")
        lines.append("  $FLEXPEP -database \"$ROSETTA_DB\" -s \"$prepack\" -flexPepDocking:pep_refine -nstruct 200 -constraints:cst_file \"$cst\" -constraints:cst_fa_weight 1.0 -out:path:all refine -out:suffix _refine")
        lines.append("done")
        lines.append("")
        lines.append("# Relax best models (edit input glob as needed)")
        lines.append("for f in refine/*_refine_*.pdb; do")
        lines.append("  $ROSETTA_SCRIPTS -database \"$ROSETTA_DB\" -s \"$f\" -parser:protocol \"$RELAX_XML\" -parser:script_vars MMP8_END=$MMP8_END PEP_START=$PEP_START PEP_END=$PEP_END -out:path:all relax -out:suffix _relax")
        lines.append("done")

        run_path.write_text("\n".join(lines) + "\n")
        run_path.chmod(0o755)

    print("Done. Outputs:")
    print(f"  {outdir}")


if __name__ == "__main__":
    main()

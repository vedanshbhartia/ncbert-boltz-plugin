#!/usr/bin/env python3
"""Validate NCAA side-chain geometry in PDB outputs against PDB CCD references.

For each selected residue in one or more PDB files, this script checks:
1) heavy side-chain atom presence versus CCD atom IDs
2) side-chain-involving bond lengths versus CCD ideal coordinates
"""

from __future__ import annotations

import argparse
import math
import shlex
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


CANONICAL_AA = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

ROSETTA_TO_CCD = {
    "A06": "HIC",
    "A82": "MHS",
    "B35": "FTR",
    "DILE": "DIL",
    "DLYS": "DLY",
    "DTHR": "DTH",
    "DTYR": "DTY",
    "HP2": "NEP",
    "HPR": "HYP",
    "NLU": "NLE",
    "V03": "3FG",
}

TWO_LETTER_ELEMENTS = {
    "BR",
    "CA",
    "CL",
    "CO",
    "CU",
    "FE",
    "HG",
    "MG",
    "MN",
    "NA",
    "NI",
    "SE",
    "SI",
    "ZN",
}


@dataclass(frozen=True)
class ObservedAtom:
    element: str
    xyz: tuple[float, float, float]


@dataclass(frozen=True)
class ObservedResidue:
    pdb_path: Path
    chain: str
    resseq: int
    icode: str
    resname: str
    atoms: dict[str, ObservedAtom]

    @property
    def label(self) -> str:
        chain = self.chain or "_"
        icode = self.icode or ""
        return f"{self.pdb_path}:{chain}:{self.resseq}{icode}:{self.resname}"


@dataclass(frozen=True)
class CCDAtom:
    element: str
    ideal_xyz: tuple[float, float, float] | None


@dataclass(frozen=True)
class CCDComponent:
    code: str
    atoms: dict[str, CCDAtom]
    bonds: frozenset[tuple[str, str]]


@dataclass(frozen=True)
class ValidationResult:
    residue: ObservedResidue
    ccd_code: str
    expected_sidechain_heavy_atoms: int
    observed_sidechain_heavy_atoms: int
    missing_sidechain_atoms: tuple[str, ...]
    extra_sidechain_atoms: tuple[str, ...]
    bond_checks: int
    max_bond_dev: float
    worst_bond: tuple[str, str, float, float, float] | None

    @property
    def passed(self) -> bool:
        if self.observed_sidechain_heavy_atoms == 0:
            return False
        if self.missing_sidechain_atoms:
            return False
        if self.extra_sidechain_atoms:
            return False
        if self.bond_checks == 0:
            return False
        return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdb",
        action="append",
        default=[],
        help="PDB file to validate (repeat for multiple files)",
    )
    parser.add_argument(
        "--run-dir",
        action="append",
        default=[],
        help="Directory to recursively scan for *.pdb outputs (repeatable)",
    )
    parser.add_argument(
        "--rosetta-db",
        required=True,
        help="Rosetta database root (expects chemical/pdb_components)",
    )
    parser.add_argument(
        "--chain",
        default="",
        help="Optional chain ID filter (single-character chain in PDB column 22)",
    )
    parser.add_argument(
        "--resname",
        action="append",
        default=[],
        help="Residue name to validate (repeatable). Default: noncanonical polymer residues",
    )
    parser.add_argument(
        "--max-bond-dev",
        type=float,
        default=0.12,
        help="Maximum allowed absolute bond-length deviation in Angstroms (default: 0.12)",
    )
    parser.add_argument(
        "--min-checked",
        type=int,
        default=1,
        help="Require at least this many residues to be validated (default: 1)",
    )
    parser.add_argument(
        "--fail-on-missing-ccd",
        action="store_true",
        help="Treat missing CCD references as failures instead of skips",
    )
    parser.add_argument(
        "--fail-on-empty-sidechain",
        action="store_true",
        help="Fail when a selected residue has no heavy side-chain atoms in the model",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Only print summary and failures",
    )
    return parser.parse_args()


def infer_element(atom_name: str) -> str:
    letters = "".join(ch for ch in atom_name if ch.isalpha()).upper()
    if not letters:
        return ""
    if len(letters) >= 2 and letters[:2] in TWO_LETTER_ELEMENTS:
        return letters[:2]
    return letters[:1]


def is_hydrogen(element: str) -> bool:
    return element.upper() in {"H", "D", "T"}


def parse_float(value: str) -> float | None:
    if value in {"?", "."}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(
        (a[0] - b[0]) * (a[0] - b[0])
        + (a[1] - b[1]) * (a[1] - b[1])
        + (a[2] - b[2]) * (a[2] - b[2])
    )


def collect_pdb_paths(pdb_args: list[str], run_dirs: list[str]) -> list[Path]:
    out: set[Path] = set()
    for value in pdb_args:
        path = Path(value)
        if not path.is_file():
            raise FileNotFoundError(f"PDB file not found: {path}")
        out.add(path)
    for run_dir in run_dirs:
        root = Path(run_dir)
        if not root.is_dir():
            raise FileNotFoundError(f"Run directory not found: {root}")
        for path in root.rglob("*.pdb"):
            if path.is_file():
                out.add(path)
    return sorted(out)


def parse_pdb_residues(pdb_path: Path, chain_filter: str) -> list[ObservedResidue]:
    residues: dict[tuple[str, int, str, str], dict[str, ObservedAtom]] = {}
    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.startswith(("ATOM  ", "HETATM")):
                continue
            if len(raw) < 54:
                continue
            altloc = raw[16]
            if altloc not in {" ", "A"}:
                continue
            chain = raw[21].strip()
            if chain_filter and chain != chain_filter:
                continue
            resname = raw[17:20].strip().upper()
            atom_name = raw[12:16].strip()
            if not atom_name:
                continue
            try:
                resseq = int(raw[22:26])
                x = float(raw[30:38])
                y = float(raw[38:46])
                z = float(raw[46:54])
            except ValueError:
                continue
            icode = raw[26].strip()
            element = raw[76:78].strip().upper() if len(raw) >= 78 else ""
            if not element:
                element = infer_element(atom_name)

            key = (chain, resseq, icode, resname)
            atom_map = residues.setdefault(key, {})
            atom_map.setdefault(atom_name, ObservedAtom(element=element, xyz=(x, y, z)))

    out: list[ObservedResidue] = []
    for (chain, resseq, icode, resname), atoms in sorted(
        residues.items(), key=lambda item: (item[0][0], item[0][1], item[0][2], item[0][3])
    ):
        out.append(
            ObservedResidue(
                pdb_path=pdb_path,
                chain=chain,
                resseq=resseq,
                icode=icode,
                resname=resname,
                atoms=atoms,
            )
        )
    return out


def tokenize_cif_line(line: str) -> list[str]:
    lexer = shlex.shlex(line, posix=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    return list(lexer)


def parse_cif_loops(block_lines: list[str]) -> list[tuple[list[str], list[list[str]]]]:
    loops: list[tuple[list[str], list[list[str]]]] = []
    i = 0
    n = len(block_lines)
    while i < n:
        if block_lines[i].strip() != "loop_":
            i += 1
            continue
        i += 1
        headers: list[str] = []
        while i < n and block_lines[i].lstrip().startswith("_"):
            headers.append(block_lines[i].strip())
            i += 1
        if not headers:
            continue

        rows: list[list[str]] = []
        buf: list[str] = []
        while i < n:
            stripped = block_lines[i].strip()
            if not stripped:
                i += 1
                continue
            if stripped == "#":
                i += 1
                break
            if stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
                break
            buf.extend(tokenize_cif_line(block_lines[i]))
            while len(buf) >= len(headers):
                rows.append(buf[: len(headers)])
                buf = buf[len(headers) :]
            i += 1
        loops.append((headers, rows))
    return loops


def read_component_block_lines(cif_path: Path, code: str) -> list[str]:
    block_tag = f"data_{code.upper()}"
    in_block = False
    block_lines: list[str] = []
    with cif_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if stripped.startswith("data_"):
                if in_block:
                    break
                in_block = stripped == block_tag
                continue
            if in_block:
                block_lines.append(line)
    if not in_block:
        raise ValueError(f"CCD block {block_tag} not found in {cif_path}")
    return block_lines


def components_file(rosetta_db: Path, ccd_code: str) -> Path:
    letter = ccd_code[0].upper()
    return rosetta_db / "chemical" / "pdb_components" / f"components.{letter}.cif"


@lru_cache(maxsize=None)
def load_ccd_component(rosetta_db: str, ccd_code: str) -> CCDComponent:
    code = ccd_code.upper()
    rosetta_db_path = Path(rosetta_db)
    cif_path = components_file(rosetta_db_path, code)
    if not cif_path.is_file():
        raise FileNotFoundError(f"CCD file not found: {cif_path}")

    block_lines = read_component_block_lines(cif_path, code)
    loops = parse_cif_loops(block_lines)

    atom_rows: list[dict[str, str]] | None = None
    bond_rows: list[dict[str, str]] | None = None
    for headers, rows in loops:
        header_set = set(headers)
        if (
            "_chem_comp_atom.atom_id" in header_set
            and "_chem_comp_atom.type_symbol" in header_set
        ):
            atom_rows = [dict(zip(headers, row)) for row in rows]
        elif (
            "_chem_comp_bond.atom_id_1" in header_set
            and "_chem_comp_bond.atom_id_2" in header_set
        ):
            bond_rows = [dict(zip(headers, row)) for row in rows]

    if not atom_rows:
        raise ValueError(f"No _chem_comp_atom loop found for {code}")
    if not bond_rows:
        raise ValueError(f"No _chem_comp_bond loop found for {code}")

    atoms: dict[str, CCDAtom] = {}
    for row in atom_rows:
        atom_id = row.get("_chem_comp_atom.atom_id", "").strip()
        if not atom_id:
            continue
        element = row.get("_chem_comp_atom.type_symbol", "").strip().upper()
        if not element:
            element = infer_element(atom_id)
        x = parse_float(row.get("_chem_comp_atom.pdbx_model_Cartn_x_ideal", "?"))
        y = parse_float(row.get("_chem_comp_atom.pdbx_model_Cartn_y_ideal", "?"))
        z = parse_float(row.get("_chem_comp_atom.pdbx_model_Cartn_z_ideal", "?"))
        if x is None or y is None or z is None:
            x = parse_float(row.get("_chem_comp_atom.model_Cartn_x", "?"))
            y = parse_float(row.get("_chem_comp_atom.model_Cartn_y", "?"))
            z = parse_float(row.get("_chem_comp_atom.model_Cartn_z", "?"))
        xyz = None if x is None or y is None or z is None else (x, y, z)
        atoms[atom_id] = CCDAtom(element=element, ideal_xyz=xyz)

    bonds: set[tuple[str, str]] = set()
    for row in bond_rows:
        a1 = row.get("_chem_comp_bond.atom_id_1", "").strip()
        a2 = row.get("_chem_comp_bond.atom_id_2", "").strip()
        if not a1 or not a2 or a1 == a2:
            continue
        bonds.add((a1, a2) if a1 < a2 else (a2, a1))

    if not atoms:
        raise ValueError(f"No atoms parsed for CCD {code}")
    if not bonds:
        raise ValueError(f"No bonds parsed for CCD {code}")
    return CCDComponent(code=code, atoms=atoms, bonds=frozenset(bonds))


def should_check_residue(residue: ObservedResidue, selected_resnames: set[str]) -> bool:
    if selected_resnames:
        return residue.resname in selected_resnames
    if residue.resname in CANONICAL_AA:
        return False
    return {"N", "CA", "C"}.issubset(residue.atoms)


def sidechain_heavy_atom_names_from_ccd(component: CCDComponent) -> set[str]:
    out = set()
    for atom_name, atom in component.atoms.items():
        if atom_name in BACKBONE_ATOMS:
            continue
        if is_hydrogen(atom.element):
            continue
        out.add(atom_name)
    return out


def sidechain_heavy_atom_names_from_observed(residue: ObservedResidue) -> set[str]:
    out = set()
    for atom_name, atom in residue.atoms.items():
        if atom_name in BACKBONE_ATOMS:
            continue
        element = atom.element or infer_element(atom_name)
        if is_hydrogen(element):
            continue
        out.add(atom_name)
    return out


def validate_residue(
    residue: ObservedResidue,
    component: CCDComponent,
) -> ValidationResult:
    ref_sc_atoms = sidechain_heavy_atom_names_from_ccd(component)
    obs_sc_atoms = sidechain_heavy_atom_names_from_observed(residue)

    missing = tuple(sorted(ref_sc_atoms - obs_sc_atoms))
    extra = tuple(sorted(obs_sc_atoms - ref_sc_atoms))

    bond_checks = 0
    max_dev = 0.0
    worst_bond: tuple[str, str, float, float, float] | None = None

    common_atoms = set(residue.atoms) & set(component.atoms)
    for a1, a2 in component.bonds:
        if a1 not in common_atoms or a2 not in common_atoms:
            continue
        if not (a1 in ref_sc_atoms or a2 in ref_sc_atoms):
            continue
        ref1 = component.atoms[a1].ideal_xyz
        ref2 = component.atoms[a2].ideal_xyz
        if ref1 is None or ref2 is None:
            continue

        expected = dist(ref1, ref2)
        observed = dist(residue.atoms[a1].xyz, residue.atoms[a2].xyz)
        dev = abs(observed - expected)
        bond_checks += 1
        if dev > max_dev:
            max_dev = dev
            worst_bond = (a1, a2, observed, expected, dev)

    return ValidationResult(
        residue=residue,
        ccd_code=component.code,
        expected_sidechain_heavy_atoms=len(ref_sc_atoms),
        observed_sidechain_heavy_atoms=len(obs_sc_atoms),
        missing_sidechain_atoms=missing,
        extra_sidechain_atoms=extra,
        bond_checks=bond_checks,
        max_bond_dev=max_dev,
        worst_bond=worst_bond,
    )


def describe_result(
    result: ValidationResult,
    max_bond_dev: float,
    fail_on_empty_sidechain: bool,
) -> tuple[str, str]:
    if (
        result.expected_sidechain_heavy_atoms > 0
        and result.observed_sidechain_heavy_atoms == 0
    ):
        status = "FAIL" if fail_on_empty_sidechain else "SKIP"
        return (status, "no_heavy_sidechain_atoms_present_in_model")
    if result.missing_sidechain_atoms:
        return (
            "FAIL",
            f"missing_sidechain_atoms={','.join(result.missing_sidechain_atoms)}",
        )
    if result.extra_sidechain_atoms:
        return (
            "FAIL",
            f"extra_sidechain_atoms={','.join(result.extra_sidechain_atoms)}",
        )
    if result.bond_checks == 0:
        return ("FAIL", "no_sidechain_bond_checks")
    if result.max_bond_dev > max_bond_dev:
        if result.worst_bond is None:
            return ("FAIL", f"max_bond_dev={result.max_bond_dev:.3f}A")
        a1, a2, observed, expected, dev = result.worst_bond
        return (
            "FAIL",
            (
                f"max_bond_dev={dev:.3f}A bond={a1}-{a2} "
                f"observed={observed:.3f}A expected={expected:.3f}A"
            ),
        )
    return (
        "PASS",
        f"sidechain_atoms_ok bond_checks={result.bond_checks} max_bond_dev={result.max_bond_dev:.3f}A",
    )


def main() -> int:
    args = parse_args()
    pdb_paths = collect_pdb_paths(args.pdb, args.run_dir)
    if not pdb_paths:
        raise ValueError("No input PDB files found. Use --pdb and/or --run-dir.")

    rosetta_db = Path(args.rosetta_db)
    if not rosetta_db.is_dir():
        raise FileNotFoundError(f"Rosetta DB not found: {rosetta_db}")

    selected_resnames = {name.strip().upper() for name in args.resname if name.strip()}
    selected_count = 0
    checked_count = 0
    pass_count = 0
    fail_count = 0
    skipped_ccd = 0
    skipped_empty_sidechain = 0
    missing_ccd_failures = 0

    for pdb_path in pdb_paths:
        residues = parse_pdb_residues(pdb_path, args.chain.strip())
        for residue in residues:
            if not should_check_residue(residue, selected_resnames):
                continue
            selected_count += 1
            ccd_code = ROSETTA_TO_CCD.get(residue.resname, residue.resname)
            try:
                component = load_ccd_component(str(rosetta_db), ccd_code)
            except (FileNotFoundError, ValueError) as exc:
                skipped_ccd += 1
                status = "FAIL" if args.fail_on_missing_ccd else "SKIP"
                if args.fail_on_missing_ccd:
                    missing_ccd_failures += 1
                if not args.quiet or status == "FAIL":
                    print(f"{status}\t{residue.label}\tresname={residue.resname}\tccd={ccd_code}\t{exc}")
                continue

            result = validate_residue(residue, component)
            status, detail = describe_result(
                result,
                args.max_bond_dev,
                args.fail_on_empty_sidechain,
            )
            if status == "PASS":
                checked_count += 1
                pass_count += 1
                if not args.quiet:
                    print(
                        f"PASS\t{residue.label}\tresname={residue.resname}\tccd={result.ccd_code}\t{detail}"
                    )
            elif status == "FAIL":
                checked_count += 1
                fail_count += 1
                print(
                    f"FAIL\t{residue.label}\tresname={residue.resname}\tccd={result.ccd_code}\t{detail}"
                )
            else:
                skipped_empty_sidechain += 1
                if not args.quiet:
                    print(
                        f"SKIP\t{residue.label}\tresname={residue.resname}\tccd={result.ccd_code}\t{detail}"
                    )

    print(
        "SUMMARY\t"
        f"files={len(pdb_paths)} "
        f"selected={selected_count} "
        f"checked={checked_count} "
        f"pass={pass_count} "
        f"fail={fail_count} "
        f"skipped_missing_ccd={skipped_ccd} "
        f"skipped_empty_sidechain={skipped_empty_sidechain}"
    )

    exit_code = 0
    if checked_count < args.min_checked:
        print(
            f"ERROR\tvalidated {checked_count} residues, below --min-checked={args.min_checked}"
        )
        exit_code = 1
    if fail_count > 0:
        exit_code = 1
    if missing_ccd_failures > 0:
        exit_code = 1
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())

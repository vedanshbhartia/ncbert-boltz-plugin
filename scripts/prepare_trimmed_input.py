#!/usr/bin/env python3
"""Trim residues from one chain in an mmCIF and write a PDB for Rosetta input."""

from __future__ import annotations

import argparse
import shlex
from pathlib import Path


THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-cif", required=True, help="Input mmCIF file")
    parser.add_argument("--out-pdb", required=True, help="Output PDB file")
    parser.add_argument("--chain", default="A", help="Chain to trim (default: A)")
    parser.add_argument(
        "--delete",
        required=True,
        help="Comma-separated residue numbers (label_seq_id in mmCIF), e.g. 4,5",
    )
    parser.add_argument(
        "--expect-old-seq",
        default=None,
        help="Optional expected sequence for the chain before deletion",
    )
    parser.add_argument(
        "--expect-new-seq",
        default=None,
        help="Optional expected sequence for the chain after deletion",
    )
    return parser.parse_args()


def parse_delete_set(raw: str) -> set[int]:
    values: set[int] = set()
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        values.add(int(token))
    if not values:
        raise ValueError("No residue numbers provided in --delete")
    return values


def parse_mmcif_atom_site(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    for i, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        j = i + 1
        cols: list[str] = []
        while j < len(lines) and lines[j].startswith("_"):
            cols.append(lines[j].strip())
            j += 1
        if not cols or not any(col.startswith("_atom_site.") for col in cols):
            continue
        colnames = [c.split(".", 1)[1] for c in cols]
        rows: list[dict[str, str]] = []
        k = j
        while k < len(lines):
            raw = lines[k].strip()
            if not raw or raw == "#" or raw == "loop_" or raw.startswith("_"):
                break
            fields = shlex.split(raw)
            if len(fields) != len(colnames):
                raise ValueError(
                    f"Unexpected atom_site field count on line {k + 1}: "
                    f"expected {len(colnames)}, got {len(fields)}"
                )
            rows.append(dict(zip(colnames, fields)))
            k += 1
        if rows:
            return rows
    raise ValueError("Could not find _atom_site loop in input mmCIF")


def safe_int(value: str) -> int | None:
    if value in (".", "?"):
        return None
    try:
        return int(value)
    except ValueError:
        return int(float(value))


def safe_float(value: str, default: float) -> float:
    if value in (".", "?"):
        return default
    return float(value)


def format_atom_name(atom_name: str, element: str) -> str:
    if len(atom_name) >= 4:
        return atom_name[:4]
    if atom_name and atom_name[0].isdigit():
        return f"{atom_name:<4}"
    if len(element.strip()) == 1:
        return f" {atom_name:<3}"
    return f"{atom_name:<4}"


def collect_chain_sequence(rows: list[dict[str, str]], chain: str) -> tuple[str, list[int]]:
    residues: dict[int, str] = {}
    for row in rows:
        if row["label_asym_id"] != chain:
            continue
        resid = safe_int(row["label_seq_id"])
        if resid is None:
            continue
        residues.setdefault(resid, row["label_comp_id"])
    ordered = sorted(residues)
    seq = "".join(THREE_TO_ONE.get(residues[r], "X") for r in ordered)
    return seq, ordered


def write_trimmed_pdb(
    rows: list[dict[str, str]],
    out_pdb: Path,
    trim_chain: str,
    delete_set: set[int],
) -> tuple[str, str, int, int]:
    original_seq, original_resids = collect_chain_sequence(rows, trim_chain)
    missing = sorted(delete_set.difference(original_resids))
    if missing:
        raise ValueError(f"Residues not found in chain {trim_chain}: {missing}")

    kept_resids = [r for r in original_resids if r not in delete_set]
    renumber = {old: new for new, old in enumerate(kept_resids, start=1)}

    kept_seq = "".join(
        THREE_TO_ONE.get(
            next(
                row["label_comp_id"]
                for row in rows
                if row["label_asym_id"] == trim_chain and safe_int(row["label_seq_id"]) == resid
            ),
            "X",
        )
        for resid in kept_resids
    )

    out_pdb.parent.mkdir(parents=True, exist_ok=True)

    sorted_rows = sorted(rows, key=lambda r: safe_int(r["id"]) or 0)
    serial = 1
    prev_chain = None
    prev_resname = None
    prev_resseq = None

    with out_pdb.open("w", encoding="utf-8") as handle:
        for row in sorted_rows:
            record = row.get("group_PDB", "ATOM")
            if record not in ("ATOM", "HETATM"):
                continue

            label_chain = row["label_asym_id"]
            label_resid = safe_int(row["label_seq_id"])
            if label_chain == trim_chain and label_resid in delete_set:
                continue

            chain_id = row.get("auth_asym_id", ".")
            if chain_id in (".", "?"):
                chain_id = label_chain

            if label_chain == trim_chain:
                if label_resid is None:
                    raise ValueError("Encountered chain residue with missing label_seq_id")
                resseq = renumber[label_resid]
            else:
                auth_seq = safe_int(row.get("auth_seq_id", "."))
                resseq = auth_seq if auth_seq is not None else (label_resid or 0)

            resname = row["label_comp_id"]
            atom_name = format_atom_name(row["label_atom_id"], row["type_symbol"])
            altloc = row["label_alt_id"]
            altloc = " " if altloc in (".", "?") else altloc
            icode = row.get("pdbx_PDB_ins_code", ".")
            icode = " " if icode in (".", "?") else icode

            x = safe_float(row["Cartn_x"], 0.0)
            y = safe_float(row["Cartn_y"], 0.0)
            z = safe_float(row["Cartn_z"], 0.0)
            occ = safe_float(row.get("occupancy", "1.0"), 1.0)
            bfac = safe_float(row.get("B_iso_or_equiv", "0.0"), 0.0)
            element = row["type_symbol"].strip()[:2]

            if prev_chain is not None and chain_id != prev_chain:
                handle.write(
                    f"TER   {serial:>5}      {prev_resname:>3} {prev_chain:1}{prev_resseq:>4}\n"
                )
                serial += 1

            handle.write(
                f"{record:<6}{serial:>5} {atom_name:4s}{altloc:1s}{resname:>3s} {chain_id:1s}"
                f"{resseq:>4d}{icode:1s}   "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{bfac:>6.2f}          {element:>2s}\n"
            )
            serial += 1
            prev_chain = chain_id
            prev_resname = resname
            prev_resseq = resseq

        if prev_chain is not None:
            handle.write(
                f"TER   {serial:>5}      {prev_resname:>3} {prev_chain:1}{prev_resseq:>4}\n"
            )
            handle.write("END\n")

    return original_seq, kept_seq, len(original_resids), len(kept_resids)


def main() -> int:
    args = parse_args()
    delete_set = parse_delete_set(args.delete)
    rows = parse_mmcif_atom_site(Path(args.in_cif))
    old_seq, new_seq, old_count, new_count = write_trimmed_pdb(
        rows=rows,
        out_pdb=Path(args.out_pdb),
        trim_chain=args.chain,
        delete_set=delete_set,
    )

    print(f"Input chain {args.chain} residues: {old_count}")
    print(f"Output chain {args.chain} residues: {new_count}")
    print(f"Old sequence: {old_seq}")
    print(f"New sequence: {new_seq}")

    if args.expect_old_seq and old_seq != args.expect_old_seq:
        raise ValueError(
            f"Old sequence mismatch: expected {args.expect_old_seq}, got {old_seq}"
        )
    if args.expect_new_seq and new_seq != args.expect_new_seq:
        raise ValueError(
            f"New sequence mismatch: expected {args.expect_new_seq}, got {new_seq}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

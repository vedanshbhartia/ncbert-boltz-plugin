#!/usr/bin/env python3
"""Normalize bracketed peptide sequences for Rosetta one-letter threading.

This utility preserves canonical one-letter inputs and converts bracketed
three-letter residue tokens into Rosetta-safe one-letter thread strings.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from pathlib import Path


# Preserve D-amino-acid stereochemistry for residues supported by this build.
SUPPORTED_D_TOKEN_TO_ROSETTA = {
    "DIL": "DILE",
    "DLY": "DLYS",
    "DTY": "DTYR",
    "DTH": "DTHR",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pepseq", required=True, help="Input peptide sequence")
    parser.add_argument("--rosetta-db", required=True, help="Rosetta database path")
    return parser.parse_args()


def tokenize_peptide(pepseq: str) -> list[str]:
    tokens: list[str] = []
    i = 0
    while i < len(pepseq):
        char = pepseq[i]
        if char.isspace():
            i += 1
            continue
        if char == "[":
            close = pepseq.find("]", i + 1)
            if close == -1:
                raise ValueError(f"Unclosed bracket token in peptide: {pepseq}")
            token = pepseq[i + 1 : close].strip()
            if not token:
                raise ValueError(f"Empty bracket token in peptide: {pepseq}")
            tokens.append(f"[{token}]")
            i = close + 1
            continue
        tokens.append(char)
        i += 1
    return tokens


@lru_cache(maxsize=None)
def lookup_one_letter_from_ccd(rosetta_db: str, token: str) -> str | None:
    upper = token.upper()
    pdb_components_dir = Path(rosetta_db) / "chemical" / "pdb_components"
    cif_path = pdb_components_dir / f"components.{upper[0]}.cif"
    if not cif_path.exists():
        return None

    block_tag = f"data_{upper}"
    in_block = False
    with cif_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("data_"):
                if in_block:
                    break
                in_block = line == block_tag
                continue
            if not in_block:
                continue
            if not line.startswith("_chem_comp.one_letter_code"):
                continue
            parts = line.split()
            if len(parts) < 2:
                return None
            value = parts[-1].strip().strip("'\"")
            if value in {"?", "."}:
                return None
            if len(value) != 1:
                return None
            return value.upper()
    return None


def normalize_peptide(pepseq: str, rosetta_db: str) -> str:
    out: list[str] = []
    for token in tokenize_peptide(pepseq):
        if token.startswith("[") and token.endswith("]"):
            residue = token[1:-1].strip().upper()
            if residue in SUPPORTED_D_TOKEN_TO_ROSETTA:
                out.append(f"X[{SUPPORTED_D_TOKEN_TO_ROSETTA[residue]}]")
                continue

            one_letter = lookup_one_letter_from_ccd(rosetta_db, residue)
            if one_letter is None:
                raise ValueError(
                    f"No one-letter mapping found for bracket token [{residue}]"
                )
            out.append(one_letter)
            continue

        out.append(token.upper())
    return "".join(out)


def main() -> int:
    args = parse_args()
    normalized = normalize_peptide(args.pepseq, args.rosetta_db)
    print(normalized)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

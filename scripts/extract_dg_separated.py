#!/usr/bin/env python3
"""Extract dG_separated from a Rosetta output PDB and print it.

Reads the bare key-value lines written by InterfaceAnalyzerMover.
Exits silently (empty output) if the file does not exist or the term
is absent — the caller treats an empty result as "no reference available".
"""

import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) < 2:
        return 0
    pdb_path = Path(sys.argv[1])
    if not pdb_path.exists():
        return 0

    with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) == 2 and parts[0] == "dG_separated":
                try:
                    print(f"{float(parts[1]):.6f}")
                    return 0
                except ValueError:
                    pass
    return 0


if __name__ == "__main__":
    main()

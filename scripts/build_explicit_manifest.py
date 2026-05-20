#!/usr/bin/env python3
"""Build a Rosetta thread-sequence manifest from an explicit sequence list.

Each input line is whitespace-separated: ``<job_id> <score> <pepseq>``.
``job_id`` may include a trailing ``_<nstruct_index>`` suffix (e.g. ``_0001``);
it is stripped so the canonical id matches a ``*_design.jsonl`` stem.
``pepseq`` is one-letter; lowercase characters are interpreted as D-amino
acids and rewritten as ``X[DXXX]`` bracket tokens that Rosetta's
SimpleThreadingMover accepts.

The matching perturbed backbone PDB is resolved from the ``PerturbN`` token
in the job id (same convention as build_1a85_target_manifest.py).
"""

import argparse
import csv
import re
from pathlib import Path


LOWERCASE_TO_D = {
    "a": "DALA",
    "d": "DASP",
    "e": "DGLU",
    "f": "DPHE",
    "h": "DHIS",
    "i": "DILE",
    "k": "DLYS",
    "l": "DLEU",
    "n": "DASN",
    "p": "DPRO",
    "q": "DGLN",
    "r": "DARG",
    "s": "DSER",
    "t": "DTHR",
    "v": "DVAL",
    "w": "DTRP",
    "y": "DTYR",
}

NSTRUCT_SUFFIX_RE = re.compile(r"_(\d{4,})$")
PERTURB_RE = re.compile(r"(alpha_\d+_1a85_cyclic_\d+_Perturb\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        required=True,
        help="Text file with whitespace-separated <job_id> <score> <pepseq> rows.",
    )
    parser.add_argument(
        "--backbones-dir",
        default="perturbed_backbones/perturbed_backbones",
        help="Directory holding alpha_..._PerturbN.pdb scaffolds.",
    )
    parser.add_argument("--out", required=True, help="Output jobs TSV path")
    parser.add_argument(
        "--selected-out",
        default=None,
        help="Optional richer TSV including the original score column.",
    )
    parser.add_argument(
        "--expected-length",
        type=int,
        default=18,
        help="Expected peptide length, counting bracket blocks as one residue.",
    )
    return parser.parse_args()


def strip_nstruct(job_id: str) -> str:
    return NSTRUCT_SUFFIX_RE.sub("", job_id)


def resolve_backbone(job_id: str, backbones_dir: Path) -> Path:
    match = PERTURB_RE.search(job_id)
    if not match:
        raise ValueError(f"Cannot extract Perturb token from job_id: {job_id}")
    backbone = backbones_dir / f"{match.group(1)}.pdb"
    if not backbone.is_file():
        raise FileNotFoundError(f"No backbone for {job_id}: expected {backbone}")
    return backbone


def to_thread_sequence(pepseq: str) -> str:
    out: list[str] = []
    for ch in pepseq:
        if ch.isupper():
            out.append(ch)
        elif ch in LOWERCASE_TO_D:
            out.append(f"X[{LOWERCASE_TO_D[ch]}]")
        elif ch.isspace():
            continue
        else:
            raise ValueError(f"Unsupported character in sequence: {ch!r} ({pepseq})")
    return "".join(out)


def count_residues(thread_seq: str) -> int:
    count = 0
    i = 0
    while i < len(thread_seq):
        ch = thread_seq[i]
        if ch == "[":
            close = thread_seq.find("]", i + 1)
            if close == -1:
                raise ValueError(f"Unclosed bracket: {thread_seq}")
            count += 1
            i = close + 1
        elif ch == "X":
            i += 1
        else:
            count += 1
            i += 1
    return count


def main() -> int:
    args = parse_args()
    backbones_dir = Path(args.backbones_dir)
    if not backbones_dir.is_dir():
        raise FileNotFoundError(f"Backbones dir not found: {backbones_dir}")

    rows: list[dict[str, str]] = []
    input_path = Path(args.input)
    with input_path.open("r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"{input_path}:{line_no} expected 3 fields, got: {line}")
            raw_job_id, score, pepseq_raw = parts[0], parts[1], parts[2]
            job_id = strip_nstruct(raw_job_id)
            backbone = resolve_backbone(job_id, backbones_dir)
            thread_seq = to_thread_sequence(pepseq_raw)
            residues = count_residues(thread_seq)
            if residues != args.expected_length:
                raise ValueError(
                    f"{input_path}:{line_no} length {residues} != "
                    f"{args.expected_length}: {pepseq_raw} -> {thread_seq}"
                )
            rows.append(
                {
                    "job_id": job_id,
                    "source_file": input_path.name,
                    "line_no": str(line_no),
                    "pepseq": thread_seq,
                    "input_pdb": str(backbone),
                    "raw_pepseq": pepseq_raw,
                    "score": score,
                }
            )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["job_id", "source_file", "line_no", "pepseq", "input_pdb"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in writer.fieldnames})

    if args.selected_out:
        selected_path = Path(args.selected_out)
        selected_path.parent.mkdir(parents=True, exist_ok=True)
        with selected_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "job_id",
                    "source_file",
                    "line_no",
                    "raw_pepseq",
                    "pepseq",
                    "score",
                    "input_pdb",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

    print(f"Wrote {len(rows)} Rosetta jobs to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

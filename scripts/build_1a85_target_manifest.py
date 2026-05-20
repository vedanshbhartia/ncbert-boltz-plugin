#!/usr/bin/env python3
"""Select 1a85 design candidates and emit Rosetta thread-sequence jobs."""

import argparse
import csv
import json
from collections import deque
from pathlib import Path


CANONICAL_AA3_TO_AA1 = {
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

TARGET_CATEGORY_ORDER = ("HYP", "AIB", "DPRO_LIKE", "PRO")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--design-dir",
        default="design",
        help="Directory containing *_design.jsonl files (default: design)",
    )
    parser.add_argument("--out", required=True, help="Output jobs TSV path")
    parser.add_argument(
        "--selected-out",
        required=True,
        help="Output TSV describing the selected designs",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=50,
        help="Number of sequences to select (default: 50)",
    )
    parser.add_argument(
        "--backbones-dir",
        default="perturbed_backbones/perturbed_backbones",
        help=(
            "Directory containing per-perturbation scaffold PDBs. "
            "For a design source file alpha_..._PerturbN_design.jsonl, "
            "the matching backbone is alpha_..._PerturbN.pdb in this directory."
        ),
    )
    return parser.parse_args()


def resolve_backbone(source_file: str, backbones_dir: Path) -> Path:
    stem = Path(source_file).stem
    suffix = "_design"
    if not stem.endswith(suffix):
        raise ValueError(
            f"Design source file does not end in '_design.jsonl': {source_file}"
        )
    backbone_stem = stem[: -len(suffix)]
    backbone_path = backbones_dir / f"{backbone_stem}.pdb"
    if not backbone_path.is_file():
        raise FileNotFoundError(
            f"No matching backbone for {source_file}: expected {backbone_path}"
        )
    return backbone_path


def classify_tokens(tokens: list[str]) -> str | None:
    if "HYP" in tokens or "DHYP" in tokens:
        return "HYP"
    if "AIB" in tokens:
        return "AIB"
    if "DPRO" in tokens or "DPR" in tokens:
        return "DPRO_LIKE"
    if "PRO" in tokens:
        return "PRO"
    return None


def has_only_supported_tokens(tokens: list[str]) -> bool:
    supported = set(CANONICAL_AA3_TO_AA1) | {"HYP", "AIB", "DPRO", "DHYP", "DPR"}
    return set(tokens).issubset(supported)


def to_thread_sequence(tokens: list[str]) -> str:
    out: list[str] = []
    for token in tokens:
        if token in CANONICAL_AA3_TO_AA1:
            out.append(CANONICAL_AA3_TO_AA1[token])
        elif token == "HYP":
            out.append("X[HPR]")
        elif token == "AIB":
            out.append("X[AIB]")
        elif token in {"DPRO", "DPR"}:
            out.append("X[DPRO]")
        elif token == "DHYP":
            out.append("X[DHYP]")
        else:
            raise ValueError(f"Unsupported token in selected sequence: {token}")
    return "".join(out)


def load_candidates(design_dir: Path, backbones_dir: Path) -> list[dict[str, str]]:
    candidates: list[dict[str, str]] = []
    for json_path in sorted(design_dir.glob("*_design.jsonl")):
        backbone_path = resolve_backbone(json_path.name, backbones_dir)
        data = json.loads(json_path.read_text(encoding="utf-8"))
        for sample_index, sample in enumerate(data["samples_b"], start=1):
            tokens = sample["aa3_tokens"]
            if not has_only_supported_tokens(tokens):
                continue
            category = classify_tokens(tokens)
            if category is None:
                continue
            candidates.append(
                {
                    "source_file": json_path.name,
                    "line_no": str(sample_index),
                    "job_id": f"{json_path.stem}_S{sample_index:04d}",
                    "category": category,
                    "pepseq": to_thread_sequence(tokens),
                    "aa3_sequence": " ".join(tokens),
                    "aa1_string": sample["aa1_string"],
                    "score": f"{sample['score']:.6f}",
                    "ncaa_list": ",".join(sample["ncaa_list"]),
                    "input_pdb": str(backbone_path),
                }
            )
    return candidates


def round_robin_select(candidates: list[dict[str, str]], limit: int) -> list[dict[str, str]]:
    buckets = {
        category: deque(candidate for candidate in candidates if candidate["category"] == category)
        for category in TARGET_CATEGORY_ORDER
    }

    selected: list[dict[str, str]] = []
    while len(selected) < limit:
        progress = False
        for category in TARGET_CATEGORY_ORDER:
            bucket = buckets[category]
            if not bucket:
                continue
            selected.append(bucket.popleft())
            progress = True
            if len(selected) >= limit:
                break
        if not progress:
            break
    return selected


def write_jobs(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["job_id", "source_file", "line_no", "pepseq", "input_pdb"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "job_id": row["job_id"],
                    "source_file": row["source_file"],
                    "line_no": row["line_no"],
                    "pepseq": row["pepseq"],
                    "input_pdb": row["input_pdb"],
                }
            )


def write_selected(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "job_id",
                "source_file",
                "line_no",
                "category",
                "score",
                "aa1_string",
                "aa3_sequence",
                "pepseq",
                "ncaa_list",
                "input_pdb",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    design_dir = Path(args.design_dir)
    backbones_dir = Path(args.backbones_dir)
    if not backbones_dir.is_dir():
        raise FileNotFoundError(f"Backbones directory not found: {backbones_dir}")
    candidates = load_candidates(design_dir, backbones_dir)
    selected = round_robin_select(candidates, args.limit)
    if len(selected) != args.limit:
        raise ValueError(f"Only found {len(selected)} eligible sequences; expected {args.limit}")

    out_path = Path(args.out)
    selected_path = Path(args.selected_out)
    write_jobs(out_path, selected)
    write_selected(selected_path, selected)

    counts: dict[str, int] = {}
    for row in selected:
        counts[row["category"]] = counts.get(row["category"], 0) + 1

    print(f"Wrote {len(selected)} Rosetta jobs to {out_path}")
    print(f"Wrote selection detail to {selected_path}")
    for category in TARGET_CATEGORY_ORDER:
        print(f"{category}: {counts.get(category, 0)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

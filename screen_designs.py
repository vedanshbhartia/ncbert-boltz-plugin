#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path


ALIASES = {
    "DMET": "DME",
    "DPRO": "DPR",
    "DCYS": "DCY",
    "DHYP": "DHYP",
}


def normalize_residue(residue):
    token = str(residue).strip().upper()
    return ALIASES.get(token, token)


def load_design_file(path):
    text = Path(path).read_text(encoding="utf-8")
    text = text.strip()
    if not text:
        return []

    # Preferred format in this repo: one pretty JSON object with "samples_b".
    try:
        obj = json.loads(text)
        if isinstance(obj, dict):
            return obj.get("samples_b", [])
    except json.JSONDecodeError:
        pass

    # Fallback: true JSONL (one JSON object per line).
    samples = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        row = json.loads(line)
        if isinstance(row, dict) and "aa3_tokens" in row:
            samples.append(row)
        elif isinstance(row, dict) and "samples_b" in row:
            samples.extend(row.get("samples_b", []))
    return samples


def estimate_net_charge(aa3_tokens):
    # Simple side-chain estimate (rough pH ~7):
    # +1: ARG/LYS (and D-forms), -1: ASP/GLU (and D-forms), HIS ignored (0).
    positive = {"ARG", "LYS", "DAR", "DLY"}
    negative = {"ASP", "GLU", "DAN", "DGL"}
    charge = 0
    for residue in aa3_tokens:
        r = normalize_residue(residue)
        if r in positive:
            charge += 1
        elif r in negative:
            charge -= 1
    return charge


def evaluate_sample(aa3_tokens):
    residues = [normalize_residue(token) for token in aa3_tokens]
    n = len(residues)

    banned = {"MET", "DME", "GLY", "CYS", "DCY"}
    proline_like = {"PRO", "DPR", "AIB", "HYP", "DHYP"}

    net_charge = estimate_net_charge(residues)
    banned_present = sorted({r for r in residues if r in banned})
    proline_like_count = sum(1 for r in residues if r in proline_like)
    min_proline_like = 1
    has_asp_or_glu = any(r in {"ASP", "GLU", "DAN", "DGL"} for r in residues)
    arg_count = sum(1 for r in residues if r in {"ARG", "DAR"})

    failures = []
    if net_charge > 0:
        failures.append("net_charge_gt_0")
    if banned_present:
        failures.append("contains_banned_residues")
    if proline_like_count < min_proline_like:
        failures.append("insufficient_proline_like")
    if not has_asp_or_glu:
        failures.append("missing_ASP_or_GLU")
    if arg_count > 2:
        failures.append("too_many_ARG")

    return {
        "length": n,
        "net_charge": net_charge,
        "banned_present": banned_present,
        "proline_like_count": proline_like_count,
        "min_proline_like": min_proline_like,
        "has_asp_or_glu": has_asp_or_glu,
        "arg_count": arg_count,
        "pass": len(failures) == 0,
        "failures": failures,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Screen ncAA ProteinMPNN design outputs against peptide rules."
    )
    parser.add_argument(
        "--input_dir",
        type=Path,
        default=Path("design"),
        help="Directory containing *_design.jsonl files.",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="*_design.jsonl",
        help="Glob pattern for design files inside --input_dir.",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=Path("screened"),
        help="Directory to write screening reports.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob(args.pattern))
    if not files:
        raise FileNotFoundError(f"No files matched {args.pattern} in {input_dir}")

    all_rows = []
    passed_rows = []
    per_file_summary = []

    for file_path in files:
        samples = load_design_file(file_path)
        total = len(samples)
        passed = 0

        for idx, sample in enumerate(samples, start=1):
            aa3_tokens = sample.get("aa3_tokens", [])
            stats = evaluate_sample(aa3_tokens)
            row = {
                "source_file": file_path.name,
                "sample_index": idx,
                "aa3_sequence": " ".join(aa3_tokens),
                "aa1_string": sample.get("aa1_string", ""),
                "length": stats["length"],
                "net_charge": stats["net_charge"],
                "banned_present": ";".join(stats["banned_present"]),
                "proline_like_count": stats["proline_like_count"],
                "min_proline_like": stats["min_proline_like"],
                "has_asp_or_glu": stats["has_asp_or_glu"],
                "arg_count": stats["arg_count"],
                "pass": stats["pass"],
                "failures": ";".join(stats["failures"]),
            }
            all_rows.append(row)
            if stats["pass"]:
                passed_rows.append(row)
                passed += 1

        per_file_summary.append(
            {
                "source_file": file_path.name,
                "total_samples": total,
                "passed_samples": passed,
                "pass_rate": (passed / total) if total > 0 else 0.0,
            }
        )

    # Write detailed per-sample CSV.
    details_csv = output_dir / "screening_details.csv"
    if all_rows:
        fieldnames = list(all_rows[0].keys())
        with details_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows)

    # Write passed-only CSV.
    passed_csv = output_dir / "screening_passed.csv"
    if passed_rows:
        fieldnames = list(passed_rows[0].keys())
        with passed_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(passed_rows)
    else:
        passed_csv.write_text("", encoding="utf-8")

    # Write per-file summary JSON.
    summary_json = output_dir / "screening_summary.json"
    summary_json.write_text(json.dumps(per_file_summary, indent=2), encoding="utf-8")

    total_all = len(all_rows)
    total_pass = len(passed_rows)
    print(f"Processed {len(files)} files, {total_all} samples total.")
    print(f"Passed: {total_pass} ({(100.0 * total_pass / total_all) if total_all else 0.0:.2f}%)")
    print(f"Details: {details_csv}")
    print(f"Passed : {passed_csv}")
    print(f"Summary: {summary_json}")


if __name__ == "__main__":
    main()

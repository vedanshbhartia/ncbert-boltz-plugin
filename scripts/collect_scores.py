#!/usr/bin/env python3
"""Collect Rosetta job status and scores into a single TSV summary."""

import argparse
import csv
from collections import deque
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", required=True, help="Input jobs TSV")
    parser.add_argument("--run-dir", required=True, help="Rosetta batch run directory")
    parser.add_argument(
        "--out",
        default=None,
        help="Output summary TSV (default: <run-dir>/summary.tsv)",
    )
    parser.add_argument(
        "--reference-dg",
        type=float,
        default=None,
        help="Reference dG_separated (REU) from backbone run; when provided adds a ddg column",
    )
    return parser.parse_args()


def read_status_env(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    if not path.exists():
        return data
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[key] = value
    return data


def parse_best_score_terms(path: Path) -> tuple[dict[str, float] | None, str | None]:
    if not path.exists():
        return None, None
    header: list[str] | None = None
    best_terms: dict[str, float] | None = None
    best_desc: str | None = None

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("SCORE:"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            if parts[1] == "total_score":
                header = parts[1:]
                continue
            if header is None:
                continue
            values = parts[1:]
            if len(values) != len(header):
                continue
            row = dict(zip(header, values))
            raw_total = row.get("total_score")
            if raw_total is None:
                continue
            try:
                total = float(raw_total)
            except ValueError:
                continue
            description = row.get("description")
            if best_terms is None or total < best_terms.get("total_score", float("inf")):
                numeric_terms: dict[str, float] = {}
                for key, raw in row.items():
                    if key == "description":
                        continue
                    try:
                        numeric_terms[key] = float(raw)
                    except ValueError:
                        continue
                best_terms = numeric_terms
                best_desc = description

    return best_terms, best_desc


def tail_error_line(log_path: Path) -> str:
    if not log_path.exists():
        return ""
    keep = deque(maxlen=40)
    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            keep.append(line.rstrip("\n"))
    for line in reversed(keep):
        stripped = line.strip()
        if not stripped:
            continue
        if any(token in stripped for token in ("ERROR", "Exception", "EXCN", "FATAL", "failed")):
            return stripped
    for line in reversed(keep):
        stripped = line.strip()
        if stripped:
            return stripped
    return ""


def resolve_model_path(job_dir: Path, description: str | None) -> str:
    if description:
        candidate = job_dir / f"{description}.pdb"
        if candidate.exists():
            return str(candidate)
    pdbs = sorted(job_dir.glob("*.pdb"))
    if pdbs:
        return str(pdbs[0])
    return ""


def parse_pose_terms_from_pdb(pdb_path: Path) -> dict[str, float] | None:
    if not pdb_path.exists():
        return None
    header_terms: list[str] | None = None
    out: dict[str, float] = {}
    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("label "):
                parts = stripped.split()
                if len(parts) > 1:
                    header_terms = parts[1:]
                continue
            if stripped.startswith("pose "):
                if not header_terms:
                    continue
                values = stripped.split()[1:]
                if len(values) != len(header_terms):
                    continue
                for term, raw in zip(header_terms, values):
                    key = "total_score" if term == "total" else term
                    try:
                        out[key] = float(raw)
                    except ValueError:
                        continue
                continue
            # Bare "key value" lines written by InterfaceAnalyzerMover and similar
            parts = stripped.split()
            if len(parts) == 2:
                key, raw = parts
                if key[0].isalpha() or key[0] == "_":
                    try:
                        out[key] = float(raw)
                    except ValueError:
                        pass
    return out if out else None


def main() -> int:
    args = parse_args()
    jobs_path = Path(args.jobs)
    run_dir = Path(args.run_dir)
    out_path = Path(args.out) if args.out else run_dir / "summary.tsv"
    reference_dg: float | None = args.reference_dg
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows_out: list[dict[str, str]] = []
    score_term_columns: set[str] = set()
    success = 0
    failed = 0

    with jobs_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            job_id = row["job_id"]
            job_dir = run_dir / job_id
            status_file = job_dir / "status.env"
            score_file = job_dir / "score.sc"
            log_file = job_dir / "rosetta.log"

            status_data = read_status_env(status_file)
            raw_exit = status_data.get("exit_code", "")
            exit_code = int(raw_exit) if raw_exit.isdigit() else None
            score_terms, description = parse_best_score_terms(score_file)
            model_path = resolve_model_path(job_dir, description)
            pose_terms = parse_pose_terms_from_pdb(Path(model_path)) if model_path else None

            combined_terms: dict[str, float] = {}
            if pose_terms:
                combined_terms.update(pose_terms)
            if score_terms:
                combined_terms.update(score_terms)
            total_score = combined_terms.get("total_score")

            status = "success" if (exit_code == 0 and model_path) else "failed"
            if status == "success":
                success += 1
            else:
                failed += 1

            error_msg = status_data.get("error", "").strip()
            if status == "failed" and not error_msg:
                error_msg = tail_error_line(log_file)

            dg_sep = combined_terms.get("dG_separated")
            if dg_sep is None:
                dg_sep = combined_terms.get("dg_separated")
            ddg: float | None = None
            if reference_dg is not None and dg_sep is not None:
                ddg = dg_sep - reference_dg

            out_row = {
                "job_id": job_id,
                "source_file": row["source_file"],
                "line_no": row["line_no"],
                "pepseq": row["pepseq"],
                "status": status,
                "total_score": "" if total_score is None else f"{total_score:.6f}",
                "dg_separated": "" if dg_sep is None else f"{dg_sep:.6f}",
                "ddg": "" if ddg is None else f"{ddg:.6f}",
                "model_path": model_path,
                "scorefile_path": str(score_file),
                "log_path": str(log_file),
                "error": error_msg,
            }
            for term, value in combined_terms.items():
                if term == "total_score":
                    continue
                col = f"score_{term}"
                out_row[col] = f"{value:.6f}"
                score_term_columns.add(col)
            rows_out.append(out_row)

    base_fields = [
        "job_id",
        "source_file",
        "line_no",
        "pepseq",
        "status",
        "total_score",
        "dg_separated",
        "ddg",
        "model_path",
        "scorefile_path",
        "log_path",
        "error",
    ]
    score_fields = sorted(score_term_columns)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=base_fields + score_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote summary: {out_path}")
    print(f"Success: {success}")
    print(f"Failed: {failed}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

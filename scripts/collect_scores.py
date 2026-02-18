#!/usr/bin/env python3
"""Collect Rosetta job status and scores into a single TSV summary."""

from __future__ import annotations

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
    return parser.parse_args()


def read_status_env(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    if not path.exists():
        return data
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[key] = value
    return data


def parse_scorefile(path: Path) -> tuple[float | None, str | None]:
    if not path.exists():
        return None, None
    header: list[str] | None = None
    best_total: float | None = None
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
            if best_total is None or total < best_total:
                best_total = total
                best_desc = description

    return best_total, best_desc


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


def parse_pose_total_from_pdb(pdb_path: Path) -> float | None:
    if not pdb_path.exists():
        return None
    pose_total: float | None = None
    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("pose "):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        pose_total = float(parts[-1])
                    except ValueError:
                        pass
                break
    return pose_total


def main() -> int:
    args = parse_args()
    jobs_path = Path(args.jobs)
    run_dir = Path(args.run_dir)
    out_path = Path(args.out) if args.out else run_dir / "summary.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows_out: list[dict[str, str]] = []
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
            total_score, description = parse_scorefile(score_file)
            model_path = resolve_model_path(job_dir, description)
            if total_score is None and model_path:
                total_score = parse_pose_total_from_pdb(Path(model_path))

            status = "success" if (exit_code == 0 and model_path) else "failed"
            if status == "success":
                success += 1
            else:
                failed += 1

            error_msg = status_data.get("error", "").strip()
            if status == "failed" and not error_msg:
                error_msg = tail_error_line(log_file)

            rows_out.append(
                {
                    "job_id": job_id,
                    "source_file": row["source_file"],
                    "line_no": row["line_no"],
                    "pepseq": row["pepseq"],
                    "status": status,
                    "total_score": "" if total_score is None else f"{total_score:.6f}",
                    "model_path": model_path,
                    "scorefile_path": str(score_file),
                    "log_path": str(log_file),
                    "error": error_msg,
                }
            )

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "job_id",
                "source_file",
                "line_no",
                "pepseq",
                "status",
                "total_score",
                "model_path",
                "scorefile_path",
                "log_path",
                "error",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote summary: {out_path}")
    print(f"Success: {success}")
    print(f"Failed: {failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

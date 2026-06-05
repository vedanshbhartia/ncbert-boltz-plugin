#!/usr/bin/env python3
"""Plan proline-like substitutions at non-MMP8-contacting peptide residues, and
write a new Rosetta threading manifest that the explicit threadseq batch
script can consume.

Inputs:
  --contact-tsv     Output of analyze_mmp8_contact.py
  --selected-tsv    Baseline run's selected.tsv (carries pepseq, input_pdb)
  --run-dir         Baseline run dir; per-design relaxed PDB is read from here
                    to measure phi at each candidate position.
  --out-dir         New run dir to write into (mutation_plan.tsv, jobs.tsv,
                    selected.tsv). Created if missing.

Outputs:
  <out-dir>/mutation_plan.tsv  Auditable per-substitution plan.
  <out-dir>/jobs.tsv           Manifest for run_rosetta_batch_threadseq.sh.
  <out-dir>/selected.tsv       Selected-format file with raw_pepseq, score,
                               input_pdb columns preserved + mutated pepseq.

Substitution rule (phi taken from the baseline relaxed PDB chain A):
  -100 <= phi <= -30  -> L-PRO  (token: 'P')   or HPR if CB near polar MMP8
   +30 <= phi <= +100 -> D-PRO  (token: 'X[DPRO]') or DHYP if CB near polar MMP8
  otherwise           -> no substitution (applied=0)

Use --allow-aib-fallback to substitute X[AIB] when phi is outside PRO/DPRO
ranges. AIB's preferred phi is roughly +/-55 deg; at extended phi
(~+/-120-180 deg) it pays a Rosetta intra-chain penalty during relax even
though the substitution looks innocuous on paper.

Use --no-hydroxyproline to disable the PRO->HPR / DPRO->DHYP swap and always
use plain proline.

Adjacency filter (--skip-adjacent-proline-like, default on): a candidate
position is skipped if either neighbor (i-1 or i+1, cyclic) is already
proline-like in the baseline OR is being substituted with proline-like in
this plan. Adjacent prolines/PRO-DPRO pairs invite cis-omega and ring-ring
clashes that show up as +100..+800 fa_rep penalties post-relax. When two
candidates conflict, we keep the one with the larger min_dist_any (i.e. the
more clearly non-contacting residue).

Residues whose current token is already in {P, X[HPR], X[AIB], X[DPRO],
X[DHYP]} are left untouched.

If --from-plan is provided, the plan is loaded instead of being recomputed,
which lets you hand-edit substitutions before generating the manifest.
"""

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path


PROLINE_LIKE_TOKENS = {"P", "X[HPR]", "X[AIB]", "X[DPRO]", "X[DHYP]"}

PHI_PRO_RANGE = (-100.0, -30.0)
PHI_DPRO_RANGE = (30.0, 100.0)

# Heavy-atom MMP8 "polar" atom names — backbone N, O, OXT and sidechain
# heteroatoms common to ASN/ASP/GLN/GLU/SER/THR/TYR/HIS/ARG/LYS/TRP/CYS.
POLAR_ATOMS = {
    "N", "O", "OXT",
    "OD1", "OD2", "ND2",  # Asp/Asn
    "OE1", "OE2", "NE2",  # Glu/Gln
    "OG", "OG1", "OH",    # Ser/Thr/Tyr
    "ND1", "NE2",         # His
    "NE", "NH1", "NH2",   # Arg
    "NZ",                 # Lys
    "NE1",                # Trp
    "SG",                 # Cys (counted; donor capability is weak but treat as polar)
}

HPR_REACH_THRESHOLD = 4.5  # CB-to-nearest-polar-MMP8-atom A for PRO->HPR swap

TOKEN_RE = re.compile(r"X\[[A-Za-z0-9]+\]|[A-Z]")


def tokenize_pepseq(pepseq: str) -> list[str]:
    """Split a Rosetta thread sequence into per-residue tokens."""
    tokens = TOKEN_RE.findall(pepseq)
    rebuilt = "".join(tokens)
    if rebuilt != pepseq.strip():
        raise ValueError(f"Unparseable pepseq: {pepseq!r} (rebuilt {rebuilt!r})")
    return tokens


def pick_substitution(
    phi: float,
    cb_to_polar: float,
    allow_aib: bool,
    allow_hydroxyproline: bool,
) -> tuple[str | None, str]:
    """Return (token_or_None, rationale) given backbone phi and CB-to-polar
    MMP8-atom distance at the position.

    Token=None means "leave residue unchanged" (applied=0 in the plan).
    When CB is near (<HPR_REACH_THRESHOLD A) a polar MMP8 atom and HPR is
    permitted, the L-PRO / D-PRO substitution upgrades to HPR / DHYP so the
    hydroxyl can H-bond.
    """
    if PHI_PRO_RANGE[0] <= phi <= PHI_PRO_RANGE[1]:
        if allow_hydroxyproline and cb_to_polar < HPR_REACH_THRESHOLD:
            return (
                "X[HPR]",
                f"phi={phi:.1f} in PRO range; CB->polar={cb_to_polar:.2f}A (<{HPR_REACH_THRESHOLD}A) -> HPR",
            )
        return "P", f"phi={phi:.1f} in PRO range; CB->polar={cb_to_polar:.2f}A"
    if PHI_DPRO_RANGE[0] <= phi <= PHI_DPRO_RANGE[1]:
        if allow_hydroxyproline and cb_to_polar < HPR_REACH_THRESHOLD:
            return (
                "X[DHYP]",
                f"phi={phi:.1f} in DPRO range; CB->polar={cb_to_polar:.2f}A (<{HPR_REACH_THRESHOLD}A) -> DHYP",
            )
        return "X[DPRO]", f"phi={phi:.1f} in DPRO range; CB->polar={cb_to_polar:.2f}A"
    if allow_aib:
        return "X[AIB]", f"phi={phi:.1f} outside PRO/DPRO ranges; AIB fallback"
    return None, f"phi={phi:.1f} outside PRO/DPRO ranges; no substitution"


@dataclass
class CaResidue:
    resseq: int
    resname: str
    n: tuple[float, float, float] | None = None
    ca: tuple[float, float, float] | None = None
    c: tuple[float, float, float] | None = None
    cb: tuple[float, float, float] | None = None


def read_backbone(pdb_path: Path, chain: str) -> list[CaResidue]:
    """Read N, CA, C, CB of every residue in chain (in ascending resseq order).

    CB is optional (GLY has none) and is only used for HPR-reach heuristics.
    """
    by_res: dict[int, CaResidue] = {}
    order: list[int] = []
    with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 54 or line[21] != chain:
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            atom = line[12:16].strip()
            if atom not in {"N", "CA", "C", "CB"}:
                continue
            try:
                resseq = int(line[22:26])
                xyz = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            except ValueError:
                continue
            res = by_res.get(resseq)
            if res is None:
                res = CaResidue(resseq=resseq, resname=line[17:20].strip())
                by_res[resseq] = res
                order.append(resseq)
            if atom == "N" and res.n is None:
                res.n = xyz
            elif atom == "CA" and res.ca is None:
                res.ca = xyz
            elif atom == "C" and res.c is None:
                res.c = xyz
            elif atom == "CB" and res.cb is None:
                res.cb = xyz
    return [by_res[k] for k in sorted(order)]


def read_chain_polar_atoms(pdb_path: Path, chain: str) -> list[tuple[float, float, float]]:
    """Return list of polar heavy-atom xyz from chain (uses POLAR_ATOMS set)."""
    coords: list[tuple[float, float, float]] = []
    with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 54 or line[21] != chain:
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            atom = line[12:16].strip()
            if atom not in POLAR_ATOMS:
                continue
            try:
                coords.append(
                    (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                )
            except ValueError:
                continue
    return coords


def min_distance(point: tuple[float, float, float], cloud: list[tuple[float, float, float]]) -> float:
    if not cloud:
        return float("inf")
    return min(
        math.sqrt(
            (point[0] - c[0]) ** 2 + (point[1] - c[1]) ** 2 + (point[2] - c[2]) ** 2
        )
        for c in cloud
    )


def dihedral(p1, p2, p3, p4) -> float:
    def sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    def norm(a):
        n = math.sqrt(dot(a, a))
        return (a[0] / n, a[1] / n, a[2] / n)

    b1 = sub(p2, p1)
    b2 = sub(p3, p2)
    b3 = sub(p4, p3)
    n1 = cross(b1, b2)
    n2 = cross(b2, b3)
    m1 = cross(n1, norm(b2))
    return math.degrees(math.atan2(dot(m1, n2), dot(n1, n2)))


def compute_phi(residues: list[CaResidue], idx_one_based: int) -> float:
    """Phi at position `idx_one_based` (1..N). Cyclic: residue 1 uses residue N's C."""
    i = idx_one_based - 1
    n = len(residues)
    prev_res = residues[(i - 1) % n]
    curr = residues[i]
    for atoms_owner, label in [(prev_res, "C"), (curr, "N"), (curr, "CA"), (curr, "C")]:
        atom = getattr(atoms_owner, label.lower())
        if atom is None:
            raise ValueError(f"Missing {label} for residue {atoms_owner.resseq}")
    assert prev_res.c and curr.n and curr.ca and curr.c
    return dihedral(prev_res.c, curr.n, curr.ca, curr.c)


def find_relaxed_pdb(run_dir: Path, design_id: str) -> Path:
    """Same logic as analyze_mmp8_contact.find_design_pdb."""
    job_dir = run_dir / design_id
    candidates = sorted(p for p in job_dir.glob("*.pdb") if p.name != "input_scaffold.pdb")
    if not candidates:
        raise FileNotFoundError(f"No relaxed PDB in {job_dir}")
    return candidates[0]


def load_contact_tsv(path: Path, contact_threshold: float) -> dict[str, list[dict]]:
    """Return {design_id: [{position, resname, min_dist_any, min_dist_sc}, ...]}.

    Only rows with min_dist_any > contact_threshold are kept.
    """
    out: dict[str, list[dict]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            min_any = float(row["min_dist_any"])
            if min_any <= contact_threshold:
                continue
            out.setdefault(row["design"], []).append(
                {
                    "position": int(row["resseq"]),
                    "resname": row["resname"],
                    "min_dist_any": min_any,
                    "min_dist_sc": (
                        float("nan") if row["min_dist_sc"] == "nan"
                        else float(row["min_dist_sc"])
                    ),
                }
            )
    return out


def load_selected_tsv(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out[row["job_id"]] = row
    return out


def build_plan(
    contact_by_design: dict[str, list[dict]],
    selected_by_design: dict[str, dict],
    run_dir: Path,
    positions_filter: set[int] | None,
    peptide_chain: str,
    target_chain: str,
    allow_aib: bool,
    allow_hydroxyproline: bool,
    skip_adjacent_proline_like: bool,
) -> list[dict]:
    plan_rows: list[dict] = []
    for design_id, contacts in contact_by_design.items():
        sel = selected_by_design.get(design_id)
        if sel is None:
            print(f"  [skip] {design_id}: not in selected.tsv")
            continue
        tokens = tokenize_pepseq(sel["pepseq"])
        relaxed = find_relaxed_pdb(run_dir, design_id)
        residues = read_backbone(relaxed, peptide_chain)
        if len(residues) != len(tokens):
            print(
                f"  [skip] {design_id}: token count {len(tokens)} != residue count {len(residues)} in {relaxed}"
            )
            continue
        polar_cloud = read_chain_polar_atoms(relaxed, target_chain)
        n_res = len(residues)

        candidates: list[dict] = []
        for hit in sorted(contacts, key=lambda h: h["position"]):
            pos = hit["position"]
            if positions_filter is not None and pos not in positions_filter:
                continue
            current_token = tokens[pos - 1]
            if current_token in PROLINE_LIKE_TOKENS:
                plan_rows.append(
                    {
                        "design_id": design_id,
                        "position": pos,
                        "current_token": current_token,
                        "current_resname": hit["resname"],
                        "phi": float("nan"),
                        "new_token": current_token,
                        "min_dist_any": hit["min_dist_any"],
                        "min_dist_sc": hit["min_dist_sc"],
                        "cb_polar": float("nan"),
                        "rationale": "already proline-like; no change",
                        "applied": "0",
                    }
                )
                continue
            phi = compute_phi(residues, pos)
            cb = residues[pos - 1].cb
            cb_polar = min_distance(cb, polar_cloud) if cb is not None else float("inf")
            new_token, rationale = pick_substitution(
                phi, cb_polar, allow_aib, allow_hydroxyproline
            )
            candidates.append(
                {
                    "design_id": design_id,
                    "position": pos,
                    "current_token": current_token,
                    "current_resname": hit["resname"],
                    "phi": phi,
                    "new_token": new_token,
                    "min_dist_any": hit["min_dist_any"],
                    "min_dist_sc": hit["min_dist_sc"],
                    "cb_polar": cb_polar,
                    "rationale": rationale,
                    "applied": "1" if new_token is not None else "0",
                }
            )

        # Adjacency filter: greedy-keep candidates sorted by min_dist_any
        # descending. We only care about candidates that would actually
        # substitute (applied==1); rejected ones (applied==0 from
        # pick_substitution) stay as audit rows.
        applies = [c for c in candidates if c["applied"] == "1"]
        non_applies = [c for c in candidates if c["applied"] != "1"]
        applies.sort(key=lambda c: c["min_dist_any"], reverse=True)

        accepted_positions: set[int] = set()
        # Positions already proline-like in the baseline pepseq count as
        # blockers from the start.
        baseline_blockers = {
            i + 1 for i, t in enumerate(tokens) if t in PROLINE_LIKE_TOKENS
        }
        for cand in applies:
            pos = cand["position"]
            if skip_adjacent_proline_like:
                # cyclic neighbors
                left = ((pos - 2) % n_res) + 1
                right = (pos % n_res) + 1
                neighbor_blocked = (
                    left in baseline_blockers
                    or right in baseline_blockers
                    or left in accepted_positions
                    or right in accepted_positions
                )
                if neighbor_blocked:
                    blockers = []
                    if left in baseline_blockers:
                        blockers.append(f"baseline pos {left}={tokens[left - 1]}")
                    if right in baseline_blockers:
                        blockers.append(f"baseline pos {right}={tokens[right - 1]}")
                    if left in accepted_positions:
                        blockers.append(f"planned pos {left}")
                    if right in accepted_positions:
                        blockers.append(f"planned pos {right}")
                    cand["new_token"] = cand["current_token"]
                    cand["applied"] = "0"
                    cand["rationale"] = (
                        cand["rationale"] + "; skipped: adjacent proline-like ("
                        + ", ".join(blockers) + ")"
                    )
                    plan_rows.append(cand)
                    continue
            accepted_positions.add(pos)
            plan_rows.append(cand)

        plan_rows.extend(non_applies)

    return plan_rows


def write_plan(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "design_id", "position", "current_token", "current_resname",
                "phi", "new_token", "min_dist_any", "min_dist_sc", "cb_polar",
                "rationale", "applied",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for r in rows:
            out = dict(r)
            out["phi"] = "nan" if isinstance(r["phi"], float) and math.isnan(r["phi"]) else f"{r['phi']:.2f}"
            out["min_dist_any"] = f"{r['min_dist_any']:.3f}"
            out["min_dist_sc"] = "nan" if math.isnan(r["min_dist_sc"]) else f"{r['min_dist_sc']:.3f}"
            cb = r.get("cb_polar", float("nan"))
            if isinstance(cb, float) and (math.isnan(cb) or math.isinf(cb)):
                out["cb_polar"] = "inf" if math.isinf(cb) else "nan"
            else:
                out["cb_polar"] = f"{cb:.3f}"
            writer.writerow(out)


def load_plan(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            row["position"] = int(row["position"])
            rows.append(row)
    return rows


def apply_plan_to_selected(
    selected_by_design: dict[str, dict],
    plan_rows: list[dict],
) -> list[dict]:
    """Return mutated selected.tsv rows. Designs without applied mutations are
    still emitted unchanged so the re-run captures all baseline designs."""
    by_design: dict[str, list[dict]] = {}
    for r in plan_rows:
        if r.get("applied", "1") != "1":
            continue
        if r["new_token"] == r["current_token"]:
            continue
        by_design.setdefault(r["design_id"], []).append(r)

    out_rows: list[dict] = []
    for design_id, sel in selected_by_design.items():
        tokens = tokenize_pepseq(sel["pepseq"])
        applied = by_design.get(design_id, [])
        for r in applied:
            pos = int(r["position"])
            tokens[pos - 1] = r["new_token"]
        new_pepseq = "".join(tokens)
        new_row = dict(sel)
        new_row["pepseq"] = new_pepseq
        new_row["n_mutations"] = str(len(applied))
        out_rows.append(new_row)
    return out_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--contact-tsv", type=Path, required=True)
    parser.add_argument("--selected-tsv", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True,
                        help="Baseline run directory containing per-design relaxed PDBs.")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--contact-threshold", type=float, default=5.0)
    parser.add_argument("--peptide-chain", default="A")
    parser.add_argument("--positions", nargs="*", type=int, default=None,
                        help="Restrict substitutions to these residue positions (default: any non-contact).")
    parser.add_argument("--designs", nargs="*", default=None,
                        help="Restrict to these design ids (default: all in contact-tsv).")
    parser.add_argument("--from-plan", type=Path, default=None,
                        help="Skip planning; load this plan TSV (e.g. after hand-editing).")
    parser.add_argument("--allow-aib-fallback", action="store_true",
                        help="Substitute X[AIB] when phi is outside PRO/DPRO ranges. "
                             "Default leaves such positions unchanged because AIB at "
                             "extended phi typically pays a Rosetta intra-chain penalty.")
    parser.add_argument("--no-hydroxyproline", action="store_true",
                        help="Disable the PRO->HPR / DPRO->DHYP upgrade when CB sits near "
                             "a polar MMP8 atom; always use plain proline.")
    parser.add_argument("--no-adjacency-filter", action="store_true",
                        help="Disable the rule that skips a candidate if any neighbor (i+/-1) "
                             "is already proline-like in the baseline or planned.")
    parser.add_argument("--target-chain", default="B",
                        help="Chain id of MMP8 in the relaxed PDB (default: B). Used to find "
                             "polar atoms for the HPR-reach heuristic.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    selected_by_design = load_selected_tsv(args.selected_tsv)
    if args.designs:
        selected_by_design = {k: v for k, v in selected_by_design.items() if k in set(args.designs)}

    if args.from_plan is not None:
        plan_rows = load_plan(args.from_plan)
    else:
        contact_by_design = load_contact_tsv(args.contact_tsv, args.contact_threshold)
        if args.designs:
            contact_by_design = {k: v for k, v in contact_by_design.items() if k in set(args.designs)}
        positions_filter = set(args.positions) if args.positions else None
        plan_rows = build_plan(
            contact_by_design,
            selected_by_design,
            args.run_dir,
            positions_filter,
            args.peptide_chain,
            args.target_chain,
            args.allow_aib_fallback,
            not args.no_hydroxyproline,
            not args.no_adjacency_filter,
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    plan_path = args.out_dir / "mutation_plan.tsv"
    write_plan(plan_path, plan_rows)

    mutated_rows = apply_plan_to_selected(selected_by_design, plan_rows)

    jobs_path = args.out_dir / "jobs.tsv"
    with jobs_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["job_id", "source_file", "line_no", "pepseq", "input_pdb"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for i, row in enumerate(mutated_rows, start=1):
            writer.writerow(
                {
                    "job_id": row["job_id"],
                    "source_file": "proline_mutations.tsv",
                    "line_no": str(i),
                    "pepseq": row["pepseq"],
                    "input_pdb": row["input_pdb"],
                }
            )

    sel_out_path = args.out_dir / "selected.tsv"
    sel_fields = [
        "job_id", "source_file", "line_no", "raw_pepseq", "pepseq",
        "score", "input_pdb", "n_mutations",
    ]
    with sel_out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=sel_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in mutated_rows:
            writer.writerow({k: row.get(k, "") for k in sel_fields})

    n_applied = sum(1 for r in plan_rows if r.get("applied", "1") == "1" and r["new_token"] != r["current_token"])
    n_designs_touched = len({r["design_id"] for r in plan_rows if r.get("applied", "1") == "1" and r["new_token"] != r["current_token"]})
    print(f"Plan: {len(plan_rows)} candidate position(s), {n_applied} applied across {n_designs_touched} design(s).")
    print(f"  plan : {plan_path}")
    print(f"  jobs : {jobs_path}")
    print(f"  sel  : {sel_out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

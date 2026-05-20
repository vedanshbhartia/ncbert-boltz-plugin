#!/usr/bin/env bash
set -u -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEFAULT_CONFIG_FILE="$REPO_ROOT/config/rosetta.env"

usage() {
  cat <<'USAGE'
Usage:
  scripts/run_rosetta_batch_threadseq.sh \
    [--config config/rosetta.env] \
    --jobs jobs.tsv \
    [--input-pdb input.pdb] \
    --xml relax_script.xml \
    [--run-dir runs/rosetta_batch_YYYYmmdd_HHMMSS] \
    [--rosetta-bin /path/to/rosetta_scripts.static.linuxgccrelease] \
    [--rosetta-db /path/to/database] \
    [--extra-res-fa-dir ncaa_params] \
    [--peptide-chain A] \
    [--include-hetatm-backbone-check] \
    [--parallel 8] \
    [--nstruct 1] \
    [--skip-backbone-check]

Required:
  --jobs         Job manifest TSV (job_id, source_file, line_no, pepseq[, input_pdb])
  --xml          RosettaScripts XML protocol

Optional:
  --input-pdb    Fallback scaffold PDB used when a job row has no input_pdb column.

If the jobs TSV contains a 5th input_pdb column, each job uses its per-row
scaffold; --input-pdb is only used when that column is empty. At least one of
the two must resolve to a real PDB for every job.

This variant expects pepseq to already contain the final Rosetta
thread_sequence string (for example AX[AIB]P...).
USAGE
}

DEFAULT_ROSETTA_ROOT="$HOME/Downloads/rosetta_binary_ubuntu_3.15_bundle/rosetta.binary.ubuntu.release-408/main"
CONFIG_FILE="${ROSETTA_CONFIG:-$DEFAULT_CONFIG_FILE}"
ROSETTA_ROOT="${ROSETTA_ROOT:-}"
ROSETTA_BIN="${ROSETTA_BIN:-}"
ROSETTA_DB="${ROSETTA_DB:-}"
BACKBONE_CHECK_SCRIPT="$SCRIPT_DIR/check_backbone_continuity.py"
JOBS_TSV=""
INPUT_PDB=""
XML_FILE=""
RUN_DIR=""
PARALLEL=1
NSTRUCT=1
SKIP_BACKBONE_CHECK=0
PEPTIDE_CHAIN="A"
INCLUDE_HETATM_BACKBONE_CHECK=0
EXTRA_RES_FA_DIR="${EXTRA_RES_FA_DIR:-}"
ROSETTA_BIN_SET=0
ROSETTA_DB_SET=0
EXTRA_RES_FA_DIR_SET=0
declare -a EXTRA_RES_FA_FILES=()
declare -a EXTRA_RES_FA_SKIPPED=()
declare -A DEFAULT_FA_STANDARD_NAMES=()

trim_whitespace() {
  local value="$1"
  value="${value#"${value%%[![:space:]]*}"}"
  value="${value%"${value##*[![:space:]]}"}"
  printf '%s' "$value"
}

expand_config_value() {
  local value="$1"
  value="${value//\$\{HOME\}/$HOME}"
  value="${value//\$HOME/$HOME}"
  value="${value//\$\{ROSETTA_ROOT\}/${ROSETTA_ROOT:-}}"
  value="${value//\$ROSETTA_ROOT/${ROSETTA_ROOT:-}}"
  printf '%s' "$value"
}

load_config_file() {
  local path="$1"
  local line key value

  if [[ -z "$path" ]]; then
    return 0
  fi
  if [[ ! -f "$path" ]]; then
    echo "Config file not found: $path" >&2
    exit 1
  fi

  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%#*}"
    line="$(trim_whitespace "$line")"
    if [[ -z "$line" ]]; then
      continue
    fi
    if [[ "$line" != *=* ]]; then
      echo "Ignoring malformed config line: $line" >&2
      continue
    fi
    key="$(trim_whitespace "${line%%=*}")"
    value="$(trim_whitespace "${line#*=}")"

    if [[ ${#value} -ge 2 && "${value:0:1}" == '"' && "${value: -1}" == '"' ]]; then
      value="${value:1:${#value}-2}"
    elif [[ ${#value} -ge 2 && "${value:0:1}" == "'" && "${value: -1}" == "'" ]]; then
      value="${value:1:${#value}-2}"
    fi
    value="${value//\$\{HOME\}/$HOME}"
    value="${value//\$HOME/$HOME}"

    case "$key" in
      ROSETTA_ROOT)
        if [[ -z "$ROSETTA_ROOT" ]]; then
          ROSETTA_ROOT="$value"
        fi
        ;;
      ROSETTA_BIN)
        if [[ $ROSETTA_BIN_SET -eq 0 && -z "$ROSETTA_BIN" ]]; then
          ROSETTA_BIN="$value"
        fi
        ;;
      ROSETTA_DB)
        if [[ $ROSETTA_DB_SET -eq 0 && -z "$ROSETTA_DB" ]]; then
          ROSETTA_DB="$value"
        fi
        ;;
      EXTRA_RES_FA_DIR)
        if [[ $EXTRA_RES_FA_DIR_SET -eq 0 && -z "$EXTRA_RES_FA_DIR" ]]; then
          EXTRA_RES_FA_DIR="$value"
        fi
        ;;
    esac
  done < "$path"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_FILE="$2"
      shift 2
      ;;
    --jobs)
      JOBS_TSV="$2"
      shift 2
      ;;
    --input-pdb)
      INPUT_PDB="$2"
      shift 2
      ;;
    --xml)
      XML_FILE="$2"
      shift 2
      ;;
    --run-dir)
      RUN_DIR="$2"
      shift 2
      ;;
    --rosetta-bin)
      ROSETTA_BIN="$2"
      ROSETTA_BIN_SET=1
      shift 2
      ;;
    --rosetta-db)
      ROSETTA_DB="$2"
      ROSETTA_DB_SET=1
      shift 2
      ;;
    --extra-res-fa-dir)
      EXTRA_RES_FA_DIR="$2"
      EXTRA_RES_FA_DIR_SET=1
      shift 2
      ;;
    --parallel)
      PARALLEL="$2"
      shift 2
      ;;
    --nstruct)
      NSTRUCT="$2"
      shift 2
      ;;
    --peptide-chain)
      PEPTIDE_CHAIN="$2"
      shift 2
      ;;
    --include-hetatm-backbone-check)
      INCLUDE_HETATM_BACKBONE_CHECK=1
      shift
      ;;
    --skip-backbone-check)
      SKIP_BACKBONE_CHECK=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ "$CONFIG_FILE" != /* ]]; then
  CONFIG_FILE="$REPO_ROOT/$CONFIG_FILE"
fi

load_config_file "$CONFIG_FILE"

if [[ -z "$ROSETTA_ROOT" ]]; then
  ROSETTA_ROOT="$DEFAULT_ROSETTA_ROOT"
fi
if [[ $ROSETTA_BIN_SET -eq 0 && -z "$ROSETTA_BIN" ]]; then
  ROSETTA_BIN="$ROSETTA_ROOT/source/bin/rosetta_scripts.static.linuxgccrelease"
fi
if [[ $ROSETTA_DB_SET -eq 0 && -z "$ROSETTA_DB" ]]; then
  ROSETTA_DB="$ROSETTA_ROOT/database"
fi
if [[ $EXTRA_RES_FA_DIR_SET -eq 0 && -z "$EXTRA_RES_FA_DIR" ]]; then
  EXTRA_RES_FA_DIR="ncaa_params"
fi

ROSETTA_BIN="$(expand_config_value "$ROSETTA_BIN")"
ROSETTA_DB="$(expand_config_value "$ROSETTA_DB")"
EXTRA_RES_FA_DIR="$(expand_config_value "$EXTRA_RES_FA_DIR")"

if [[ -n "$EXTRA_RES_FA_DIR" && "$EXTRA_RES_FA_DIR" != /* ]]; then
  EXTRA_RES_FA_DIR="$REPO_ROOT/$EXTRA_RES_FA_DIR"
fi

if [[ -z "$JOBS_TSV" || -z "$XML_FILE" ]]; then
  echo "Missing required arguments." >&2
  usage
  exit 1
fi

if [[ -z "$RUN_DIR" ]]; then
  RUN_DIR="runs/rosetta_batch_$(date +%Y%m%d_%H%M%S)"
fi

if [[ ! -f "$JOBS_TSV" ]]; then
  echo "Jobs TSV not found: $JOBS_TSV" >&2
  exit 1
fi
if [[ -n "$INPUT_PDB" && ! -f "$INPUT_PDB" ]]; then
  echo "Input PDB not found: $INPUT_PDB" >&2
  exit 1
fi
if [[ ! -f "$XML_FILE" ]]; then
  echo "XML file not found: $XML_FILE" >&2
  exit 1
fi
if [[ ! -x "$ROSETTA_BIN" ]]; then
  echo "Rosetta executable not found or not executable: $ROSETTA_BIN" >&2
  exit 1
fi
if [[ ! -d "$ROSETTA_DB" ]]; then
  echo "Rosetta database directory not found: $ROSETTA_DB" >&2
  exit 1
fi
if [[ "$SKIP_BACKBONE_CHECK" -eq 0 && ! -f "$BACKBONE_CHECK_SCRIPT" ]]; then
  echo "Backbone check script not found: $BACKBONE_CHECK_SCRIPT" >&2
  exit 1
fi
if ! [[ "$PARALLEL" =~ ^[0-9]+$ ]] || [[ "$PARALLEL" -lt 1 ]]; then
  echo "Invalid --parallel value: $PARALLEL" >&2
  exit 1
fi
if ! [[ "$NSTRUCT" =~ ^[0-9]+$ ]] || [[ "$NSTRUCT" -lt 1 ]]; then
  echo "Invalid --nstruct value: $NSTRUCT" >&2
  exit 1
fi
if [[ -z "$PEPTIDE_CHAIN" ]]; then
  echo "Invalid --peptide-chain value: must be non-empty." >&2
  exit 1
fi

residue_types_index="$ROSETTA_DB/chemical/residue_type_sets/fa_standard/residue_types.txt"
residue_types_root="$ROSETTA_DB/chemical/residue_type_sets/fa_standard"
if [[ -f "$residue_types_index" ]]; then
  while IFS= read -r rel_path; do
    rel_path="${rel_path%%#*}"
    rel_path="${rel_path#"${rel_path%%[![:space:]]*}"}"
    rel_path="${rel_path%"${rel_path##*[![:space:]]}"}"
    if [[ -z "$rel_path" ]]; then
      continue
    fi
    if [[ "$rel_path" != residue_types/*.params ]]; then
      continue
    fi
    full_path="$residue_types_root/$rel_path"
    if [[ ! -f "$full_path" ]]; then
      continue
    fi
    res_name="$(awk '$1=="NAME"{print $2; exit}' "$full_path")"
    if [[ -n "$res_name" ]]; then
      DEFAULT_FA_STANDARD_NAMES["$res_name"]=1
    fi
  done < "$residue_types_index"
fi

if [[ -n "$EXTRA_RES_FA_DIR" && -d "$EXTRA_RES_FA_DIR" ]]; then
  while IFS= read -r -d '' param_file; do
    extra_name="$(awk '$1=="NAME"{print $2; exit}' "$param_file")"
    if [[ -n "$extra_name" && -n "${DEFAULT_FA_STANDARD_NAMES[$extra_name]:-}" ]]; then
      EXTRA_RES_FA_SKIPPED+=("$param_file:$extra_name")
      continue
    fi
    EXTRA_RES_FA_FILES+=("$param_file")
  done < <(find "$EXTRA_RES_FA_DIR" -maxdepth 1 -type f -name '*.params' -print0 | sort -z)
fi

declare -A UNIQUE_INPUT_PDBS=()
declare -a JOB_INPUT_PDB_LIST=()
HEADER_LINE=""
while IFS= read -r raw_line || [[ -n "$raw_line" ]]; do
  if [[ -z "$raw_line" ]]; then
    continue
  fi
  IFS=$'\t' read -r f1 f2 f3 f4 f5 <<<"$raw_line"
  f1="${f1%$'\r'}"
  f5="${f5%$'\r'}"
  if [[ -z "$HEADER_LINE" ]]; then
    HEADER_LINE="$raw_line"
    if [[ "$f1" == "job_id" ]]; then
      continue
    fi
  fi
  if [[ -z "$f1" ]]; then
    continue
  fi
  job_pdb="$f5"
  if [[ -z "$job_pdb" ]]; then
    job_pdb="$INPUT_PDB"
  fi
  if [[ -z "$job_pdb" ]]; then
    echo "Job $f1 has no input_pdb column and no --input-pdb fallback was provided." >&2
    exit 1
  fi
  if [[ ! -f "$job_pdb" ]]; then
    echo "Input PDB for job $f1 not found: $job_pdb" >&2
    exit 1
  fi
  JOB_INPUT_PDB_LIST+=("$f1=$job_pdb")
  UNIQUE_INPUT_PDBS["$job_pdb"]=1
done < "$JOBS_TSV"

if [[ "$SKIP_BACKBONE_CHECK" -eq 0 ]]; then
  for pdb_path in "${!UNIQUE_INPUT_PDBS[@]}"; do
    check_cmd=(
      python3 "$BACKBONE_CHECK_SCRIPT"
      --pdb "$pdb_path"
      --chain "$PEPTIDE_CHAIN"
      --cyclic
      --max-cn-distance 1.8
      --quiet
    )
    if [[ "$INCLUDE_HETATM_BACKBONE_CHECK" -eq 1 ]]; then
      check_cmd+=(--include-hetatm)
    fi
    if ! "${check_cmd[@]}"; then
      echo "Backbone continuity preflight failed for $pdb_path; refusing to run Rosetta batch." >&2
      echo "Use --skip-backbone-check only if you intentionally want to run on a broken scaffold." >&2
      exit 1
    fi
  done
fi

mkdir -p "$RUN_DIR"
if [[ "$(realpath "$JOBS_TSV")" != "$(realpath -m "$RUN_DIR/jobs.tsv")" ]]; then
  cp "$JOBS_TSV" "$RUN_DIR/jobs.tsv"
fi
if [[ -n "$INPUT_PDB" ]] && \
   [[ "$(realpath "$INPUT_PDB")" != "$(realpath -m "$RUN_DIR/input_scaffold.pdb")" ]]; then
  cp "$INPUT_PDB" "$RUN_DIR/input_scaffold.pdb"
fi
if [[ "$(realpath "$XML_FILE")" != "$(realpath -m "$RUN_DIR/protocol.xml")" ]]; then
  cp "$XML_FILE" "$RUN_DIR/protocol.xml"
fi

cat > "$RUN_DIR/run_config.txt" <<EOF
config_file=$CONFIG_FILE
rosetta_root=$ROSETTA_ROOT
jobs_tsv=$JOBS_TSV
input_pdb=$INPUT_PDB
xml_file=$XML_FILE
rosetta_bin=$ROSETTA_BIN
rosetta_db=$ROSETTA_DB
peptide_chain=$PEPTIDE_CHAIN
include_hetatm_backbone_check=$INCLUDE_HETATM_BACKBONE_CHECK
parallel=$PARALLEL
nstruct=$NSTRUCT
extra_res_fa_dir=$EXTRA_RES_FA_DIR
extra_res_fa_count=${#EXTRA_RES_FA_FILES[@]}
extra_res_fa_skipped_count=${#EXTRA_RES_FA_SKIPPED[@]}
unique_job_input_pdbs=${#UNIQUE_INPUT_PDBS[@]}
start_time=$(date -Iseconds)
EOF
for pdb_path in "${!UNIQUE_INPUT_PDBS[@]}"; do
  echo "job_input_pdb=$pdb_path" >> "$RUN_DIR/run_config.txt"
done
for skipped in "${EXTRA_RES_FA_SKIPPED[@]}"; do
  echo "extra_res_fa_skipped=$skipped" >> "$RUN_DIR/run_config.txt"
done

run_one() {
  local job_id="$1"
  local source_file="$2"
  local line_no="$3"
  local pepseq="$4"
  local job_input_pdb="$5"
  local job_dir="$RUN_DIR/$job_id"
  local log_file="$job_dir/rosetta.log"
  local score_file="$job_dir/score.sc"
  local status_file="$job_dir/status.env"
  local start_time end_time exit_code error_line model_file model_desc pose_total

  if [[ -z "$job_input_pdb" ]]; then
    job_input_pdb="$INPUT_PDB"
  fi

  mkdir -p "$job_dir"
  find "$job_dir" -maxdepth 1 -type f \
    \( -name '*.pdb' -o -name '*.cif' -o -name '*.sc' -o -name 'status.env' -o -name 'rosetta.log' \) \
    -delete
  if [[ "$(realpath "$job_input_pdb")" != "$(realpath -m "$job_dir/input_scaffold.pdb")" ]]; then
    cp "$job_input_pdb" "$job_dir/input_scaffold.pdb"
  fi
  start_time="$(date -Iseconds)"
  : > "$log_file"

  "$ROSETTA_BIN" \
    -database "$ROSETTA_DB" \
    -s "$job_input_pdb" \
    -parser:protocol "$XML_FILE" \
    -parser:script_vars "pepseq=$pepseq" \
    -load_PDB_components \
    "${EXTRA_RES_FA_FILES[@]:+-in:file:extra_res_fa}" "${EXTRA_RES_FA_FILES[@]}" \
    -nstruct "$NSTRUCT" \
    -overwrite \
    -out:path:all "$job_dir" \
    -out:file:scorefile "$score_file" \
    >>"$log_file" 2>&1
  exit_code=$?
  if grep -q "no jobs were attempted" "$log_file"; then
    exit_code=2
  fi
  end_time="$(date -Iseconds)"

  if [[ "$exit_code" -eq 0 && ! -s "$score_file" ]]; then
    model_file="$(find "$job_dir" -maxdepth 1 -type f -name '*.pdb' | head -n 1)"
    if [[ -n "$model_file" ]]; then
      pose_total="$(awk '$1=="pose"{val=$NF} END{print val}' "$model_file")"
      if [[ -n "$pose_total" ]]; then
        model_desc="$(basename "$model_file" .pdb)"
        {
          echo "SCORE: total_score description"
          echo "SCORE: $pose_total $model_desc"
        } > "$score_file"
      fi
    fi
  fi

  error_line=""
  if [[ "$exit_code" -ne 0 ]]; then
    error_line="$(grep -m1 -E 'ERROR|Exception|EXCN|FATAL|unrecognized|failed' "$log_file" || true)"
    if [[ -z "$error_line" ]]; then
      error_line="$(tail -n 1 "$log_file" 2>/dev/null || true)"
    fi
  fi
  error_line="${error_line//$'\t'/ }"
  error_line="${error_line//$'\r'/ }"
  error_line="${error_line//$'\n'/ }"

  cat > "$status_file" <<EOF
job_id=$job_id
source_file=$source_file
line_no=$line_no
pepseq=$pepseq
input_pdb=$job_input_pdb
start_time=$start_time
end_time=$end_time
exit_code=$exit_code
scorefile=$score_file
log=$log_file
error=$error_line
EOF
}

running=0
total_jobs=0

while IFS=$'\t' read -r job_id source_file line_no pepseq input_pdb_col; do
  if [[ "$job_id" == "job_id" || -z "$job_id" ]]; then
    continue
  fi
  job_id="${job_id%$'\r'}"
  source_file="${source_file%$'\r'}"
  line_no="${line_no%$'\r'}"
  pepseq="${pepseq%$'\r'}"
  input_pdb_col="${input_pdb_col%$'\r'}"

  run_one "$job_id" "$source_file" "$line_no" "$pepseq" "$input_pdb_col" &
  running=$((running + 1))
  total_jobs=$((total_jobs + 1))

  if (( running >= PARALLEL )); then
    wait -n || true
    running=$((running - 1))
  fi
done < "$JOBS_TSV"

while (( running > 0 )); do
  wait -n || true
  running=$((running - 1))
done

success=0
failed=0

while IFS= read -r -d '' status_file; do
  exit_code="$(awk -F= '$1=="exit_code"{print $2}' "$status_file")"
  if [[ "$exit_code" == "0" ]]; then
    success=$((success + 1))
  else
    failed=$((failed + 1))
  fi
done < <(find "$RUN_DIR" -mindepth 2 -maxdepth 2 -name status.env -print0)

cat >> "$RUN_DIR/run_config.txt" <<EOF
end_time=$(date -Iseconds)
total_jobs=$total_jobs
success_jobs=$success
failed_jobs=$failed
EOF

echo "Run directory: $RUN_DIR"
echo "Total jobs: $total_jobs"
echo "Succeeded: $success"
echo "Failed: $failed"

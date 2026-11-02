#!/usr/bin/env bash
# Run onsight on a list of signal files, optionally with profiling.
#
# Usage: scripts/hw_prof/run_profile.sh [--prof] [--perf] [--keys] [--attr] [--pore r9|r10] [--neighbor FILE]
#                                       OUTDIR INDEX FILE_LIST [THREADS] [BATCH]
#   --prof      use bin/onsight-prof (make prof); writes OUTDIR/prof.*
#   --perf      wrap the run in `perf stat` (OUTDIR/perf_stat.txt) and `perf record`
#               at 499 Hz with call graphs (OUTDIR/perf.data)
#   --keys      count lookups per index key (ONS_PROF_KEYS=1)
#   --attr      use bin/onsight-attr (make attr): hot helpers out of line, for perf attribution
#   --pore      r9 (default) or r10
#   --neighbor  regions file for the near-neighbor check (a panel's regions.lcr); off without it
#   --panel     panel.tsv (target names); without it each index taxon is a target
#   FILE_LIST   text file, one signal file per line, in sequencing order
#   THREADS / BATCH  -t and -K (default 4 and 4M)
#
# Always writes OUTDIR/{time.txt,stderr.log,run_info.txt,read_records.jsonl,read_calls.tsv,target_calls.tsv,
# sig_summary.tsv}. Run from anywhere, with the conda environment active.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
BIN=$ROOT/bin/onsight; PERF=0; PORE=r9; NEIGHBOR=(); PANEL=()
while [[ $# -gt 0 && $1 == --* ]]; do
  case $1 in
    --prof) BIN=$ROOT/bin/onsight-prof ;;
    --perf) PERF=1 ;;
    --keys) export ONS_PROF_KEYS=1 ;;
    --attr) BIN=$ROOT/bin/onsight-attr ;;
    --pore) PORE=$2; shift ;;
    --neighbor) NEIGHBOR=(--neighbor-check "$(realpath "$2")"); shift ;;
    --panel) PANEL=(--panel "$(realpath "$2")"); shift ;;
    *) echo "unknown option $1" >&2; exit 1 ;;
  esac
  shift
done
OUT=$1; IDX=$(realpath "$2"); LIST=$3; T=${4:-4}; K=${5:-4M}
[[ -x /usr/bin/time ]] || { echo "needs GNU time: sudo apt install time" >&2; exit 1; }
[[ -x $BIN ]] || { echo "$BIN not found: run make, make prof or make attr" >&2; exit 1; }
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
export ONS_PROF_OUT="$OUT/prof"

WRAP=()
if [[ $PERF == 1 ]]; then
  WRAP=(perf stat -o "$OUT/perf_stat.txt"
        -e cycles:u,instructions:u,cache-references:u,cache-misses:u,branch-misses:u,L1-dcache-loads:u,L1-dcache-load-misses:u,dTLB-load-misses:u
        -- perf record -q -F 499 -g -o "$OUT/perf.data" --)
fi

{ echo "binary: $BIN"; echo "git: $(git -C "$ROOT" rev-parse --short HEAD) $(git -C "$ROOT" status --porcelain | wc -l) uncommitted"
  echo "start: $(date -Is)"; echo "pore: $PORE threads: $T batch: $K keys: ${ONS_PROF_KEYS:-0} perf: $PERF"
  echo "neighbor check: ${NEIGHBOR[1]:-off} panel: ${PANEL[1]:-index taxa}"
  echo "files: $(wc -l < "$LIST") from $LIST"
  echo "host: $(uname -m) $(nproc) cpus, $(lscpu 2>/dev/null | sed -n 's/^Model name: *//p' | head -1)"
  # Jetson: the power mode and clocks set the timing (sudo nvpmodel -m 0; sudo jetson_clocks)
  if command -v nvpmodel >/dev/null; then echo "jetson power mode: $(nvpmodel -q 2>/dev/null | head -1)"; fi
  if [[ -r /etc/nv_tegra_release ]]; then echo "jetson L4T: $(head -1 /etc/nv_tegra_release)"; fi
} > "$OUT/run_info.txt"

"${WRAP[@]}" /usr/bin/time -v -o "$OUT/time.txt" \
"$BIN" --pore "$PORE" -t "$T" -K "$K" "${PANEL[@]}" "${NEIGHBOR[@]}" \
  --read-calls "$OUT/read_calls.tsv" --target-calls "$OUT/target_calls.tsv" \
  --sig-summary "$OUT/sig_summary.tsv" --stream-emit-fd 3 \
  "$IDX" $(cat "$LIST") \
  1>/dev/null 2>"$OUT/stderr.log" 3>"$OUT/read_records.jsonl"

echo "end: $(date -Is)" >> "$OUT/run_info.txt"

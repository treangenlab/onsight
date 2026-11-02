#!/usr/bin/env bash
# Build one index of a panel from its signatures (signatures.fasta, signatures.csv) with strobe-event seeds (two
# strobes of 6 quantized events, the second 8-20 events after the first), as every current index. Writes <index> and
# <index>.sigmeta.tsv (the engine); the build log goes to stdout. R9.4.1: 6-mer pore model; R10.4.1: 9-mer model and
# the R10 event segmentation (--pore).
#   conda activate s2s
#   build_index.sh <signatures folder> <r9|r10> <index file> [onsight binary] [threads]
set -euo pipefail
[ $# -ge 3 ] || { sed -n 2,7p "$0"; exit 1; }
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MODELS=$HERE/../../extern/local_kmer_models
sigs=$1 chem=$2 idx=$3 bin=${4:-$(command -v onsight)} threads=${5:-10}
case $chem in
  r9)  model=$MODELS/r94_means_only.txt ;;
  r10) model=$MODELS/uncalled_r1041_model_only_means.txt ;;
  *)   echo "chemistry must be r9 or r10" >&2; exit 2 ;;
esac
mkdir -p "$(dirname "$idx")"
cmd=("$bin" --pore "$chem" -p "$model" -t "$threads" --metadata "$sigs/signatures.csv"
     --seed-mode strobe --strobe-n 2 --strobe-ks 6 --strobe-wmin 8 --strobe-wmax 20
     -d "$idx" "$sigs/signatures.fasta")
echo "${cmd[*]}"
"${cmd[@]}"

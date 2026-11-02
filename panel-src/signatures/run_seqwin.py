#!/usr/bin/env python3
"""Seqwin on one set's genome lists (targets.txt, negatives.txt of select_genomes.py), with the paper's settings.

Settings (Seqwin 0.3.1): k-mer length 21, minimizer window 200, stringency 3, minimum signature length 200 bp, Mash
threshold estimation and BLAST check on, penalty threshold cap 0.2, seed 42.

Writes Seqwin's folder --out/<set>/ (<set> = the name of the --lists folder): signatures.fasta, signatures.csv,
assemblies.csv, config.json, seqwin.log and Seqwin's work files. A set that already has signatures.fasta is skipped.

    conda activate seqwin
    python3 run_seqwin.py --lists <dir>/<set> --out <dir> [--threads 12]
"""
import argparse
from pathlib import Path

from seqwin import Config, run

SETTINGS = dict(kmerlen=21, windowsize=200, stringency=3, min_len=200, run_mash=True, run_blast=True,
                penalty_th_cap=0.2, seed=42)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lists", type=Path, required=True, help="folder with targets.txt and negatives.txt")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=12)
    a = ap.parse_args()
    name = a.lists.name
    if (a.out / name / "signatures.fasta").exists():
        print(f"{name}: built before, skipped")
        return 0
    a.out.mkdir(parents=True, exist_ok=True)
    run(Config(tar_paths=a.lists / "targets.txt", neg_paths=a.lists / "negatives.txt", prefix=a.out, title=name,
               overwrite=False, n_cpu=a.threads, **SETTINGS))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

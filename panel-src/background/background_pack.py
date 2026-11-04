#!/usr/bin/env python3
"""The background genome pack: one genome per species from a local RefSeq library, the genomes every build screens
its signatures against and draws near neighbors from (the paper's pack has 27,592 genomes).

One genome per species (NCBI species taxid) in each library folder (bacteria, archaea, fungi, viral): NCBI's reference
genome if there is one, else the most complete assembly (complete, chromosome, scaffold, contig), ties broken by
accession. Genomes are shuffled (--seed) and dealt into batches of about --batch-gb (the batch column; the panel build
does not use it).

--library: a folder per kingdom holding one FASTA per assembly, named <accession>_<assembly name>_genomic.fa (as the
RefSeq FTP site and kraken2-build --download-library name them). --summary: RefSeq's assembly_summary_refseq.txt.

Writes <out>.tsv: genome, kingdom, file, species_taxid, organism, genus, refseq_category, assembly_level, genome_size,
batch (genus: the organism name's first word; empty for viruses).

    conda activate onsight-panel
    python3 background_pack.py --library <dir> --summary <assembly_summary_refseq.txt> --out <prefix>
"""
import argparse
import csv
import os
from pathlib import Path

import numpy as np

KINGDOMS = ("bacteria", "archaea", "fungi", "viral")
FASTA = ("_genomic.fa", "_genomic.fna", "_genomic.fa.gz", "_genomic.fna.gz")   # other files (e.g. .fai) are skipped
LEVEL = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}


def accession(file_name: str) -> str:
    return "_".join(file_name.split("_")[:2])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--batch-gb", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    library = {}
    for kingdom in (k for k in KINGDOMS if (a.library / k).is_dir()):     # a kingdom folder may be missing
        for f in sorted(os.listdir(a.library / kingdom)):
            if f.endswith(FASTA):
                library[accession(f)] = (kingdom, a.library / kingdom / f)
    best = {}
    for line in open(a.summary):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] not in library:
            continue
        kingdom, path = library[f[0]]
        rank = (0 if f[4] == "reference genome" else 1, LEVEL.get(f[11], 4), f[0])
        key = (kingdom, f[6])
        if key not in best or rank < best[key][0]:
            name = f[7].replace("[", "").replace("]", "")
            best[key] = (rank, dict(genome=f[0], kingdom=kingdom, file=str(path), species_taxid=f[6], organism=f[7],
                                    genus="" if kingdom == "viral" else name.split()[0], refseq_category=f[4],
                                    assembly_level=f[11], genome_size=int(f[25])))
    rows = [r for _, r in sorted(best.values(), key=lambda v: v[0][2])]
    order = np.random.default_rng(a.seed).permutation(len(rows))
    batch, filled = 0, 0
    for i in order:
        if filled >= a.batch_gb * 1e9:
            batch, filled = batch + 1, 0
        rows[i]["batch"] = batch
        filled += rows[i]["genome_size"]
    rows = [rows[i] for i in order]
    with open(f"{a.out}.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    size = sum(r["genome_size"] for r in rows)
    print(f"{len(rows):,} genomes, {size / 1e9:.2f} Gb -> {a.out}.tsv")
    for kingdom in KINGDOMS:
        sub = [r for r in rows if r["kingdom"] == kingdom]
        print(f"  {kingdom}: {len(sub):,} genomes, {sum(r['genome_size'] for r in sub) / 1e9:.2f} Gb")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

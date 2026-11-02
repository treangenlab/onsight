#!/usr/bin/env python3
"""More genomes per near-neighbor species of one panel target, by the own-strain rule (own_strains.py): for each
near-neighbor species of the target (near_neighbors.tsv, one representative per group at >= 95% ANI), up to 10 RefSeq
complete genomes (latest) of the representative's species taxid from the local RefSeq libraries, in a seeded random
order (seed 42), skipping any at >= 99.9% ANI (skani) to a genome the test reads came from and any genome already in
the near-neighbor folder. One genome per species can lack a homolog that other strains of the species carry; these
genomes add them to the near-neighbor side.

--read-sources: a text file of genome files the test reads came from. Only those of the target's near-neighbor species
can be >= 99.9% ANI to a candidate, so the list may hold just those.

Writes <out>/genomes.tsv (genome, file), <out>/taxids.tsv (genome, species_taxid, near_neighbor_of = the target's
signature sets, organism, representative, ani_to_read_sources) and <out>/considered.tsv, for blast_homologs.py.

    conda activate squiggle
    python3 near_neighbor_strains.py --near <near-neighbor folder> --target Enterobacter --read-sources <list> \
        --summary <assembly_summary_refseq.txt> --libraries <dir> [<dir> ...] --out <dir>
"""
import argparse
import csv
import random
from collections import defaultdict
from pathlib import Path

import sys
sys.path[:0] = [str(d) for d in (Path(__file__).resolve().parents[1]).iterdir() if d.is_dir()]   # panel-src folders
from near_neighbors import read_tsv, write_tsv
from own_strains import BATCH, N_KEEP, SKIP_ANI, max_ani, on_disk as library_files


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--near", type=Path, required=True, help="near-neighbor folder (targets.tsv, near_neighbors.tsv)")
    ap.add_argument("--target", required=True)
    ap.add_argument("--read-sources", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True, help="RefSeq assembly summary")
    ap.add_argument("--libraries", type=Path, nargs="+", required=True, help="folders of RefSeq genomes (*_genomic.fa)")
    ap.add_argument("--skani", default="skani")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    sets = next(t for t in read_tsv(a.near / "targets.tsv") if t["target"] == a.target)["sets"]
    reps = [r for r in read_tsv(a.near / "near_neighbors.tsv") if r["target"] == a.target]
    listed = {l.split("\t")[0] for l in open(a.near / "genomes.tsv") if l.strip()}
    refs = [l.strip() for l in open(a.read_sources) if l.strip()]
    complete = defaultdict(dict)
    for l in open(a.summary):
        if not l.startswith("#"):
            r = l.rstrip("\n").split("\t")
            if r[11] == "Complete Genome" and r[10] == "latest":
                complete[r[6]][r[0]] = r
    on_disk = library_files(a.libraries)
    kept, considered = [], []
    for rep in reps:
        info = complete.get(rep["species_taxid"], {})
        accs = sorted(info)
        random.Random(a.seed).shuffle(accs)
        n, i = 0, 0
        while n < N_KEEP and i < len(accs):
            batch = accs[i:i + BATCH]
            i += BATCH
            files = {g: on_disk[g] for g in batch if g in on_disk and g not in listed}
            ani = max_ani(list(files.values()), refs, a.threads, a.skani) if files else {}
            for g in batch:
                v = ani.get(files.get(g), "")
                if g in listed:
                    status = "already a near neighbor"
                elif g not in on_disk:
                    status = "not on disk"
                elif v >= SKIP_ANI:
                    status = "skipped: >= 99.9% ANI to a read source"
                else:
                    status = "kept" if n < N_KEEP else "not needed"
                considered.append(dict(representative=rep["genome"], species_taxid=rep["species_taxid"], genome=g,
                                       ani_to_read_sources=v, status=status))
                if status == "kept":
                    kept.append(dict(genome=g, file=files[g], species_taxid=rep["species_taxid"],
                                     near_neighbor_of=sets, organism=f"{info[g][7]} {info[g][8]}".strip(),
                                     representative=rep["genome"], ani_to_read_sources=v))
                    n += 1
        print(rep["organism"], rep["species_taxid"], "complete genomes", len(accs), "kept", n, flush=True)
    with open(a.out / "genomes.tsv", "w") as f:
        f.writelines(f"{k['genome']}\t{k['file']}\n" for k in kept)
    write_tsv(a.out / "taxids.tsv", [{k: v for k, v in r.items() if k != "file"} for r in kept])
    write_tsv(a.out / "considered.tsv", considered)
    print(a.target, "near-neighbor species", len(reps), "with a genome kept",
          len({k["representative"] for k in kept}), "genomes kept", len(kept), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

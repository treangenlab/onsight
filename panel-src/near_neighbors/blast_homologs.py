#!/usr/bin/env python3
"""Each signature's homologous region in each genome of a folder (near neighbors or own strains), for the near-neighbor check.

blastn -task blastn, e-value 1e-5, one database per genome; every hit is written as PAF, the format
neighbor_regions.py reads (query = signature, 0-based half-open coordinates, strand, target contig, column 10 =
identical bases), so a signature that aligns to a genome in several pieces keeps each piece. Per signature, the hit
with the most identical bases (then the higher bit score) comes first; homolog_counts.tsv counts that one.

The folder holds genomes.tsv (genome, FASTA file) and taxids.tsv (genome, species_taxid, and the panel species whose
signatures go against it: column near_neighbor_of or species). Writes <folder>/paf/<genome>.paf,
<folder>/blast/<genome>.tsv.gz (all hits) and <folder>/homolog_counts.tsv; <folder>/tmp/ (query files, BLAST
databases) is removed at the end.
--reuse <folder>: a genome already searched there with the same panel species keeps that folder's blastn hits
(copied) instead of a new search; homolog_counts.tsv then also says which genomes were reused.

    conda activate onsight-panel
    python3 blast_homologs.py --folder <dir> --signatures <signatures.fasta> --metadata <signatures.csv> \
      --blastn <blastn> --makeblastdb <makeblastdb> --processes 20
"""
import argparse
import csv
import gzip
import multiprocessing as mp
import shutil
import subprocess
import tempfile
from pathlib import Path

import sys
sys.path[:0] = [str(d) for d in (Path(__file__).resolve().parents[1]).iterdir() if d.is_dir()]   # panel-src folders

FIELDS = "qseqid qlen qstart qend sseqid slen sstart send nident length pident evalue bitscore"
A = {}



def read_fasta(path: str):
    name, chunks = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)

def one(job):
    genome, fasta, species = job
    folder, out = A["folder"], {}
    old = A.get("reuse")
    reused = old is not None and A["old_species"].get(genome) == species
    if reused:
        shutil.copy(old / "blast" / f"{genome}.tsv.gz", folder / "blast" / f"{genome}.tsv.gz")
        with gzip.open(folder / "blast" / f"{genome}.tsv.gz", "rt") as f:
            hits = [l.rstrip("\n").split("\t") for l in f][1:]
    else:
        with tempfile.TemporaryDirectory(dir=folder / "tmp") as tmp:
            db = Path(tmp) / "db"
            src = fasta
            if fasta.endswith(".gz"):                   # makeblastdb reads plain FASTA only
                src = str(Path(tmp) / "genome.fa")
                with gzip.open(fasta, "rb") as fi, open(src, "wb") as fo:
                    shutil.copyfileobj(fi, fo)
            subprocess.run([A["makeblastdb"], "-in", src, "-dbtype", "nucl", "-out", str(db)], check=True,
                           capture_output=True)
            hits = []
            for sp in species:
                query = folder / "tmp" / "queries" / f"{sp}.fa"
                res = subprocess.run([A["blastn"], "-task", "blastn", "-query", str(query),
                                      "-db", str(db), "-evalue", "1e-5", "-max_target_seqs", "500", "-max_hsps", "3",
                                      "-num_threads", "1", "-outfmt", "6 " + FIELDS], check=True,
                                     capture_output=True, text=True)
                hits += [l.split("\t") for l in res.stdout.splitlines() if l]
        with gzip.open(folder / "blast" / f"{genome}.tsv.gz", "wt") as f:
            f.write(FIELDS.replace(" ", "\t") + "\n")
            f.writelines("\t".join(h) + "\n" for h in hits)
    for h in hits:
        key = (int(h[8]), float(h[12]))
        if h[0] not in out or key > out[h[0]][0]:
            out[h[0]] = (key, h)
    with open(folder / "paf" / f"{genome}.paf", "w") as f:
        for h in sorted(hits, key=lambda h: (h[0], -int(h[8]), -float(h[12]))):
            q = h[0]
            qlen, qs, qe, slen, ss, se = (int(v) for v in (h[1], h[2], h[3], h[5], h[6], h[7]))
            strand, ts, te = ("+", ss - 1, se) if ss <= se else ("-", se - 1, ss)
            f.write(f"{q}\t{qlen}\t{qs - 1}\t{qe}\t{strand}\t{h[4]}\t{slen}\t{ts}\t{te}\t{h[8]}\t{h[9]}\t255\t"
                    f"pi:f:{h[10]}\tev:f:{h[11]}\n")
    rows = []
    for sp in species:
        names = A["names"][sp]
        good = lambda h: float(h[10]) >= 75 and int(h[9]) >= 0.5 * int(h[1])
        rows.append(dict(genome=genome, species=sp, signatures=len(names),
                         blastn=sum(1 for q in names if q in out),
                         blastn_75pct_half=sum(1 for q in names if q in out and good(out[q][1])),
                         **({"reused": int(reused)} if old else {})))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", type=Path, required=True)
    ap.add_argument("--signatures", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--blastn", required=True)
    ap.add_argument("--makeblastdb", required=True)
    ap.add_argument("--processes", type=int, default=20)
    ap.add_argument("--reuse", type=Path, help="a folder searched before: its genomes keep their hits")
    a = ap.parse_args()
    for sub in ("paf", "blast", "tmp/queries"):
        (a.folder / sub).mkdir(parents=True, exist_ok=True)
    taxon = {r["fasta_header"]: r["source_taxon"].replace("_vs_genus", "") for r in csv.DictReader(open(a.metadata))}
    names = {}
    seqs = list(read_fasta(str(a.signatures)))
    for sp in sorted(set(taxon.values())):
        names[sp] = [n for n, _ in seqs if taxon[n] == sp]
        with open(a.folder / "tmp" / "queries" / f"{sp}.fa", "w") as f:
            f.writelines(f">{n}\n{s}\n" for n, s in seqs if taxon[n] == sp)
    files = dict(l.rstrip("\n").split("\t") for l in open(a.folder / "genomes.tsv") if l.strip())
    tax = list(csv.DictReader(open(a.folder / "taxids.tsv"), delimiter="\t"))
    col = "near_neighbor_of" if "near_neighbor_of" in tax[0] else "species"
    jobs = sorted(((r["genome"], files[r["genome"]], r[col].split(",")) for r in tax),
                  key=lambda j: -Path(j[1]).stat().st_size)
    A.update(folder=a.folder, blastn=a.blastn, makeblastdb=a.makeblastdb, names=names)
    if a.reuse:
        old = list(csv.DictReader(open(a.reuse / "taxids.tsv"), delimiter="\t"))
        c = "near_neighbor_of" if "near_neighbor_of" in old[0] else "species"
        A.update(reuse=a.reuse, old_species={r["genome"]: r[c].split(",") for r in old
                                              if (a.reuse / "blast" / f"{r['genome']}.tsv.gz").exists()})
        print(f"{sum(1 for g, _, s in jobs if A['old_species'].get(g) == s)} of {len(jobs)} genomes reused from "
              f"{a.reuse}", flush=True)
    rows = []
    with mp.Pool(a.processes) as pool:
        for k, part in enumerate(pool.imap_unordered(one, jobs), 1):
            rows += part
            if k % 25 == 0 or k == len(jobs):
                print(f"{k} / {len(jobs)} genomes", flush=True)
    rows.sort(key=lambda r: (r["species"], r["genome"]))
    with open(a.folder / "homolog_counts.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for sp in sorted({r["species"] for r in rows}):
        rs = [r for r in rows if r["species"] == sp]
        print(sp, "genomes", len(rs), "signature-genome pairs with a hit", sum(r["blastn"] for r in rs), flush=True)
    shutil.rmtree(a.folder / "tmp")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Own strains of the other panel targets on a target's near-neighbor side, for the near-neighbor check.

near_neighbors.py candidates drops a candidate genome that is another panel target ("other panel species <target>"),
so on its own the check never compares a target's read calls with another target's genomes, and a read of one target
that hits a signature of another is kept (24-target panel, 10-28: E. coli reads kept as Enterobacter and
S. enterica). Here, when target s dropped a candidate as target t, t's own strains (own_strains.py) go on s's
near-neighbor side: each is searched against s's signature sets. Their species_taxid stays t's taxid, so
neighbor_regions.py puts their regions on the near side of s's signatures. A target of the keep-off panel that is not
a target here has no own strains; it is reported and left out.

Writes <out>/genomes.tsv (genome, file; no header, as the other homolog folders) and <out>/taxids.tsv (genome,
species_taxid, near_neighbor_of = the sets of every such s, ncbi_species_taxid, organism) for blast_homologs.py; one
line per target pair goes to the log.

    conda activate squiggle
    python3 other_targets.py --near <near folder> --own <own folder> --out <dir>
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

DROPPED = "other panel species "


def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--near", type=Path, required=True, help="near_neighbors.py folder: targets.tsv, candidates.tsv")
    ap.add_argument("--own", type=Path, required=True, help="own_strains.py folder: genomes.tsv, taxids.tsv")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    sets = {t["target"]: t["sets"].split(",") for t in read_tsv(a.near / "targets.tsv")}
    target_of = {s: t for t, ss in sets.items() for s in ss}
    dropped = defaultdict(int)                       # (s, t): candidates of target s dropped as target t
    for r in read_tsv(a.near / "candidates.tsv"):
        if r["dropped"].startswith(DROPPED):
            for t in r["dropped"][len(DROPPED):].split("+"):
                dropped[r["target"], t] += 1
    against = defaultdict(set)                       # t -> the targets whose near side gets t's own strains
    for s, t in sorted(dropped):
        if t in sets:
            against[t].add(s)
        else:
            print(f"{s}: {dropped[s, t]} candidates kept off as {t}, not a target of this panel (no own strains)")

    files = dict(l.rstrip("\n").split("\t") for l in open(a.own / "genomes.tsv") if l.strip())
    rows, added = [], defaultdict(int)
    for r in read_tsv(a.own / "taxids.tsv"):
        t = target_of[r["species"].split(",")[0]]
        if t not in against:
            continue
        for s in against[t]:
            added[s, t] += 1
        rows.append(dict(genome=r["genome"], species_taxid=r["species_taxid"],
                         near_neighbor_of=",".join(x for s in sorted(against[t]) for x in sets[s]),
                         ncbi_species_taxid=r["ncbi_species_taxid"], organism=r["organism"]))

    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "genomes.tsv").write_text("".join(f"{r['genome']}\t{files[r['genome']]}\n" for r in rows))
    with open(a.out / "taxids.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["genome", "species_taxid", "near_neighbor_of", "ncbi_species_taxid",
                                          "organism"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for s, t in sorted(dropped):
        if t in sets:
            print(f"{s}: {added[s, t]} own strains of {t} ({dropped[s, t]} candidates dropped as it)", flush=True)
    print(f"{len(rows)} own strains of {len(against)} targets, on the near side of "
          f"{len({s for s, t in dropped if t in sets})} targets", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

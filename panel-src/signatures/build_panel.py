#!/usr/bin/env python3
"""Combine Seqwin outputs (one folder per species) into one signature panel for onsight.

Writes in --out:
  signatures.fasta  every signature, unfiltered, as Seqwin wrote it (eskap5 was built the same way)
  signatures.csv    Seqwin's per-signature table with source_taxon in front; the index reads
                    fasta_header and source_taxon from it (--metadata)
  sources.tsv       one row per species: Seqwin folder, version, run date, settings, signatures, bp

Usage:
    python3 build_panel.py --out <dir> <seqwin-out/Genus_species_vs_genus> [...]
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    recs, name, seq = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line.strip())
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("seqwin_dirs", type=Path, nargs="+")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    seen, fasta, rows, sources, header = {}, [], [], [], None
    for d in a.seqwin_dirs:
        taxon = d.name
        recs = read_fasta(d / "signatures.fasta")
        with open(d / "signatures.csv") as fh:
            table = list(csv.DictReader(fh))
        if header is None:
            header = ["source_taxon"] + list(table[0])
        if {r["fasta_header"] for r in table} != {n for n, _ in recs}:
            sys.exit(f"{d}: signatures.csv and signatures.fasta list different signatures")
        for n, _ in recs:
            if n in seen:
                sys.exit(f"duplicate signature name {n} in {taxon} and {seen[n]}")
            seen[n] = taxon
        fasta += recs
        rows += [{"source_taxon": taxon, **r} for r in table]
        log = (d / "seqwin.log").read_text()
        cfg = json.loads((d / "config.json").read_text())
        version = re.search(r"Running Seqwin (v[\d.]+)", log)
        sources.append(dict(source_taxon=taxon, seqwin_dir=str(d.resolve()),
                            seqwin_version=version.group(1) if version else "", run_date=log[:10],
                            kmerlen=cfg["kmerlen"], windowsize=cfg["windowsize"], stringency=cfg["stringency"],
                            min_len=cfg["min_len"], signatures=len(recs), bp=sum(len(s) for _, s in recs)))

    with open(a.out / "signatures.fasta", "w") as fh:
        for n, s in fasta:
            fh.write(f">{n}\n{s}\n")
    with open(a.out / "signatures.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        w.writerows(rows)
    with open(a.out / "sources.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sources[0]), delimiter="\t")
        w.writeheader()
        w.writerows(sources)
    print(f"{len(fasta):,} signatures, {sum(len(s) for _, s in fasta):,} bp, {len(sources)} species -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

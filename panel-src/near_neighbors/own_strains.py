#!/usr/bin/env python3
"""Own strains of each panel target, for the near-neighbor check's own side.

Per NCBI species of the target's signature sets (targets.tsv of near_neighbors.py sources; a lineage set shares its
species' strains): up to 10 RefSeq complete genomes (latest, from --summary) of the species taxid that are in the
local RefSeq libraries, in a seeded random order (seed 42), in batches of 12. A genome at >= 99.9% ANI (skani) to one
of the species' test strains (--test-strains: the genomes the test reads came from) is skipped. A genus target
(Enterobacter) gets the rule per species that has a set. Its strains carry the target's taxid in taxids.tsv
(species_taxid, the taxid the near-neighbor check uses for the target), their own NCBI species in ncbi_species_taxid,
and go against the sets of their own species only.

Writes <out>/genomes.tsv (genome, file), <out>/taxids.tsv (genome, species_taxid, species = the sets whose signatures
go against it, ncbi_species_taxid, organism, ani_to_read_sources) and <out>/considered.tsv.

    conda activate onsight-panel
    python3 own_strains.py --targets <near folder>/targets.tsv --test-strains <tsv: species, file> \
        --summary <assembly_summary_refseq.txt> --libraries <dir> [<dir> ...] --taxdump <dir> --out <dir>
"""
import argparse
import csv
import random
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import sys
sys.path[:0] = [str(d) for d in (Path(__file__).resolve().parents[1]).iterdir() if d.is_dir()]   # panel-src folders
from near_neighbors import Taxonomy, acc

N_KEEP, SKIP_ANI, BATCH = 10, 99.9, 12


def max_ani(files: list, refs: list, threads: int, skani: str = "skani") -> dict:
    """file -> largest skani ANI to any of refs (0 when skani reports nothing)."""
    with tempfile.TemporaryDirectory() as tmp:
        q, r = Path(tmp) / "q.txt", Path(tmp) / "r.txt"
        q.write_text("".join(f + "\n" for f in files))
        r.write_text("".join(f + "\n" for f in refs))
        out = subprocess.run([skani, "dist", "--ql", str(q), "--rl", str(r), "-t", str(threads)],
                             capture_output=True, text=True, check=True).stdout
    best = {f: 0.0 for f in files}
    for row in csv.DictReader(out.splitlines(), delimiter="\t"):
        best[row["Query_file"]] = max(best[row["Query_file"]], float(row["ANI"]))
    return best


def complete_genomes(summary: Path) -> list:
    return [l.rstrip("\n").split("\t") for l in open(summary) if not l.startswith("#")]


FASTA = ("_genomic.fa", "_genomic.fna", "_genomic.fa.gz", "_genomic.fna.gz")


def on_disk(libraries: list) -> dict:
    """accession -> genome file in the RefSeq library folders (<accession>_<name>_genomic.fa, .fna, gzipped or not)."""
    return {acc(p.name): str(p) for d in libraries for p in sorted(Path(d).glob("GC*_genomic.*"))
            if p.name.endswith(FASTA)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", type=Path, required=True)
    ap.add_argument("--test-strains", type=Path, help="species (scientific name) and file of each test strain genome")
    ap.add_argument("--summary", type=Path, required=True, help="RefSeq assembly summary")
    ap.add_argument("--libraries", type=Path, nargs="+", required=True, help="folders of RefSeq genomes (*_genomic.fa)")
    ap.add_argument("--taxdump", type=Path, required=True)
    ap.add_argument("--skani", default="skani")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    targets = list(csv.DictReader(open(a.targets), delimiter="\t"))
    test = defaultdict(list)
    for r in (csv.DictReader(open(a.test_strains), delimiter="\t") if a.test_strains else []):
        test[r["species"]].append(r["file"])
    # NCBI species -> (panel taxid, the sets whose signatures go against its strains)
    species = defaultdict(lambda: [0, []])
    for e in targets:
        for s in e["sets"].split(","):
            sci = " ".join(s.split("_")[:2])
            species[sci][0] = int(e["taxid"])
            species[sci][1].append(s)
    ids = Taxonomy(a.taxdump).ids(set(species))
    rows = complete_genomes(a.summary)
    files_of = on_disk(a.libraries)
    kept, considered = [], []
    for sci, (panel_taxid, sets) in sorted(species.items()):
        taxid = ids[sci]
        refs = test.get(sci, [])
        info = {r[0]: r for r in rows if r[6] == str(taxid) and r[11] == "Complete Genome" and r[10] == "latest"}
        accs = sorted(info)
        random.Random(a.seed).shuffle(accs)
        n, i = 0, 0
        while n < N_KEEP and i < len(accs):
            batch = accs[i:i + BATCH]
            i += BATCH
            files = {g: files_of[g] for g in batch if g in files_of}
            ani = (max_ani(list(files.values()), refs, a.threads, a.skani) if files and refs
                   else {f: 0.0 for f in files.values()})
            for g in batch:
                if g not in files:
                    considered.append(dict(species=sci, genome=g, ani_to_read_sources="", status="not on disk"))
                    continue
                v = ani[files[g]]
                status = "kept" if v < SKIP_ANI and n < N_KEEP else ("skipped: >= 99.9% ANI to a read source"
                                                                     if v >= SKIP_ANI else "not needed")
                considered.append(dict(species=sci, genome=g, ani_to_read_sources=v if refs else "", status=status))
                if status == "kept":
                    kept.append(dict(genome=g, file=files[g], species_taxid=panel_taxid, species=",".join(sets),
                                     ncbi_species_taxid=taxid, organism=f"{info[g][7]} {info[g][8]}".strip(),
                                     ani_to_read_sources=v if refs else ""))
                    n += 1
        print(sci, taxid, "complete genomes", len(accs), "kept", n, "sets", ",".join(sets), flush=True)
    with open(a.out / "genomes.tsv", "w") as f:
        f.writelines(f"{k['genome']}\t{k['file']}\n" for k in kept)
    with open(a.out / "taxids.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[k for k in kept[0] if k != "file"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(kept)
    with open(a.out / "considered.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(considered[0]), delimiter="\t")
        w.writeheader()
        w.writerows(considered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

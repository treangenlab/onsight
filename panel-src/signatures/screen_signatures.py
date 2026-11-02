#!/usr/bin/env python3
"""Which genomes carry a Seqwin build's signatures: background genomes (signatures to drop), held-out genomes (how
well the signatures cover strains the build did not see), and the mock communities' genomes (check only).

For each build (signatures in --root/seqwin/<build>/signatures.fasta, lists in --root/lists/<build>/):
- background genomes: one RefSeq genome per species of the same family in other genera (the all-genome background
  pack, --background). A background genome that carries the species' genome by the species rule (>= 95% ANI over
  >= 50% of the NCBI reference genome the ANI rule used; e.g. Shigella for E. coli) is the species itself and left
  out. A signature
  carried by any other background genome would call another species; final_sets.py drops it.
- held-out genomes (held_out/): the fraction of signatures each carries, median and 10th percentile: same_lineage,
  other_unnamed (other lineages under the name), other_named_species (assemblies under the name that the species rule
  gives to another named species; carrying means a false match).
- check only: the other species' genomes of the mock communities (--mock-genomes; the paper's: Zymo D6331 and D6322),
  and the fraction carried by the species' own mock strain genomes (--mock-strains; the lowest if several).
A genome carries a signature when minimap2 -x asm10 aligns it over >= 95% of its length in matching bases and >= 90%
of its length in the alignment.

Writes in --out: <build>.tsv (one row: signatures, bp, held-out carriage, background and mock signatures and genomes)
and <build>.carried.tsv (each background or mock genome that carries a signature, and which).

    conda activate squiggle     # minimap2, skani
    python3 screen_signatures.py --root <dir> --builds <build,...> --background <genomes.tsv> --taxdump <dir> \
        [--mock-genomes <dir> ...] [--mock-strains panels/mock_strains.tsv] --out <dir>/screen
"""
import argparse
import csv
import multiprocessing as mp
import subprocess
from pathlib import Path

import numpy as np

HELD_OUT = ("same_lineage", "other_unnamed", "other_named_species")


def taxonomy(taxdump: Path):
    """up(taxid, rank) -> the taxid's ancestor of that rank; ids: scientific name -> taxid."""
    parent, rank, ids = {}, {}, {}
    for line in open(taxdump / "nodes.dmp"):
        f = line.split("\t|\t", 3)
        parent[int(f[0])], rank[int(f[0])] = int(f[1]), f[2]
    for line in open(taxdump / "names.dmp"):
        f = line.split("\t|\t")
        if f[3].startswith("scientific name"):
            ids.setdefault(f[1], int(f[0]))
    merged = {int(f[0]): int(f[1].split("\t")[0]) for f in (l.split("\t|\t") for l in open(taxdump / "merged.dmp"))}

    def up(t, want):
        t = merged.get(t, t)
        while t in parent and t != 1:
            if rank[t] == want:
                return t
            t = parent[t]
        return None
    return up, ids


def carry_species_genome(reference: str, genomes: list, out: Path) -> set:
    """The genomes that carry the species' genome by the species rule (>= 95% ANI over >= 50% of its reference
    genome, skani)."""
    if not reference or not genomes:
        return set()
    out.with_suffix(".list").write_text("\n".join(genomes) + "\n")
    subprocess.run(["skani", "dist", "-t", "4", "-q", reference, "--rl", str(out.with_suffix(".list")), "-o",
                    str(out)], check=True, stderr=subprocess.DEVNULL)
    return {r["Ref_file"] for r in csv.DictReader(open(out), delimiter="\t")
            if float(r["ANI"]) >= 95 and float(r["Align_fraction_query"]) >= 50}


def carried(job):
    """The signatures (query names) that one genome carries."""
    genome, query = job
    out = subprocess.run(["minimap2", "-c", "-x", "asm10", "--secondary=no", "-t", "1", genome, query],
                         capture_output=True, text=True).stdout
    got = set()
    for line in out.splitlines():
        f = line.split("\t")
        q = int(f[1])
        if int(f[9]) >= 0.95 * q and int(f[10]) >= 0.9 * q:
            got.add(f[0])
    return genome, got


def write_tsv(path: Path, rows: list, fields: list) -> None:
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True, help="folder with lists/<build>/ and seqwin/<build>/")
    ap.add_argument("--builds", required=True, help="comma separated")
    ap.add_argument("--background", type=Path, required=True,
                    help="genomes.tsv of the background pack (kingdom, species_taxid, file per genome)")
    ap.add_argument("--taxdump", type=Path, required=True, help="NCBI taxonomy dump folder (nodes, names, merged)")
    ap.add_argument("--mock-genomes", type=Path, nargs="*", default=[], help="folders of mock genome FASTA (check only)")
    ap.add_argument("--mock-strains", type=Path, help="species and file of each species' own mock strain genomes")
    ap.add_argument("--processes", type=int, default=8)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    up, ids = taxonomy(a.taxdump)
    background = list(csv.DictReader(open(a.background), delimiter="\t"))
    mocks = [p for d in a.mock_genomes for p in sorted(d.glob("*.fasta"))]
    strain_files = {}
    for r in (csv.DictReader(open(a.mock_strains), delimiter="\t") if a.mock_strains else []):
        strain_files.setdefault(r["species"], []).append(r["file"])
    for build in a.builds.split(","):
        lists, signatures = a.root / "lists" / build, a.root / "seqwin" / build / "signatures.fasta"
        species = next(csv.DictReader(open(lists / "targets.tsv"), delimiter="\t"))["species"]
        reference = next((r["path"] for r in csv.DictReader(open(lists / "reference_genome.tsv"), delimiter="\t")
                          if r["species"] == species and r["used_for_ani"] == "yes"), "")
        tax = ids[species]
        family, genus = up(tax, "family"), up(tax, "genus")
        others = [r["file"] for r in background if r["kingdom"] in ("bacteria", "archaea", "fungi")
                  and up(int(r["species_taxid"]), "family") == family and up(int(r["species_taxid"]), "genus") != genus]
        species_itself = carry_species_genome(reference, others, a.out / f"{build}.background_species_itself.tsv")
        others = [g for g in others if g not in species_itself]
        own = set(strain_files.get(species, [])) | {species.replace(" ", "_")}
        mock_others = [str(p) for p in mocks if p.name not in own and not p.name.startswith(species.replace(" ", "_"))]
        mock_strains = [str(p) for p in mocks if p.name in own]
        held = {k: [l.strip() for l in open(lists / "held_out" / f"{k}.txt") if l.strip()]
                if (lists / "held_out" / f"{k}.txt").exists() else [] for k in HELD_OUT}
        n_signatures = bp = 0
        for line in open(signatures):
            if line.startswith(">"):
                n_signatures += 1
            else:
                bp += len(line.strip())
        genomes = {g: "held_out" for k in HELD_OUT for g in held[k]} | {g: "background" for g in others} | \
                  {g: "mock" for g in mock_others} | {g: "mock_strain" for g in mock_strains}
        print(f"{build}: {n_signatures} signatures, {len(held['same_lineage'])} + {len(held['other_unnamed'])} + "
              f"{len(held['other_named_species'])} held-out genomes, {len(others)} background genomes of the family "
              f"({len(species_itself)} more are the species itself), {len(mock_others)} mock genomes", flush=True)
        with mp.Pool(a.processes) as pool:
            got = dict(pool.imap_unordered(carried, [(g, str(signatures)) for g in genomes]))
        fraction = lambda g: len(got[g]) / max(n_signatures, 1)
        row = dict(build=build, species=species, signatures=n_signatures, bp=bp)
        for k in HELD_OUT:
            v = [fraction(g) for g in held[k]]
            row[f"{k}_median"] = round(float(np.median(v)), 3) if v else ""
            row[f"{k}_p10"] = round(float(np.percentile(v, 10)), 3) if v else ""
        for kind in ("background", "mock"):
            row[f"{kind}_signatures"] = len({s for g, kd in genomes.items() if kd == kind for s in got[g]})
            row[f"{kind}_genomes"] = sum(1 for g, kd in genomes.items() if kd == kind and got[g])
        row["mock_strain_carried"] = round(min(fraction(g) for g in mock_strains), 3) if mock_strains else ""
        write_tsv(a.out / f"{build}.tsv", [row], list(row))
        hits = [dict(build=build, kind=kind, genome=Path(g).name, carried=len(got[g]),
                     signatures=";".join(sorted(got[g])))
                for g, kind in genomes.items() if kind in ("background", "mock") and got[g]]
        write_tsv(a.out / f"{build}.carried.tsv", hits, ["build", "kind", "genome", "carried", "signatures"])
        print(row, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

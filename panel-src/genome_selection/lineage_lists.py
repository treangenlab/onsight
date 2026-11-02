#!/usr/bin/env python3
"""One Seqwin set per lineage when a species name covers several lineages.

Input: a set folder of select_genomes.py. Its high-quality target-labelled assemblies that the species rule leaves
under the name (targets.tsv: high_quality, no closer named species) are clustered at 95% ANI (skani triangle),
greedily from the best assembly down (CheckM completeness - 5 x contamination; assemblies with no CheckM value after
those with one), as GTDB clusters species around a representative genome. With --other-lineages-only, only the
assemblies outside the reference genome's lineage (< 95% ANI to it) are clustered: the lineages next to the set itself,
which its targets leave out (e.g. L. monocytogenes lineage I next to the reference genome's lineage II).

Each lineage with at least --min-size members is one set: its members as targets, a quarter of them (at least one,
random, --seed) held out; the negatives of the input set. The name is detected when any of its sets is.

Writes, for lineage i (in the order of the clustering), --lists/<set>_lineage<i>/ (with --other-lineages-only:
<set>_other_lineage<i>/): targets.txt, negatives.txt, held_out/same_lineage.txt (other held-out lists empty),
targets.tsv (members: accession, species, held_out, target, path), reference_genome.tsv (the input set's). Work files
in --lists/<set>/lineages/: skani_triangle.tsv, lineages.tsv (representative and members per lineage), members.tsv.

    conda activate seqwin       # skani 0.4.0 by --skani
    python3 lineage_lists.py --lists <dir> --set <name> [--other-lineages-only]
"""
import argparse
import csv
import random
import shutil
import subprocess
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lists", type=Path, required=True, help="folder of select_genomes.py's set folders")
    ap.add_argument("--set", required=True)
    ap.add_argument("--other-lineages-only", action="store_true",
                    help="cluster only the assemblies outside the reference genome's lineage (< 95%% ANI to it)")
    ap.add_argument("--min-size", type=int, default=5)
    ap.add_argument("--skani", default="skani")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    src = a.lists / a.set
    table = [r for r in csv.DictReader(open(src / "targets.tsv"), delimiter="\t")
             if r["high_quality"] == "True" and not r["closer_named_species"]
             and not (a.other_lineages_only and r["reference_lineage"] == "True")]
    quality = lambda r: (r["completeness"] != "",
                         float(r["completeness"] or 0) - 5 * float(r["contamination"] or 0))
    table.sort(key=quality, reverse=True)
    paths = [r["path"] for r in table]
    by_path = {r["path"]: r for r in table}
    work = src / "lineages"
    work.mkdir(exist_ok=True)
    listed = "\n".join(paths) + "\n"
    triangle = work / "skani_triangle.tsv"
    if not (triangle.exists() and (work / "clustered.list").exists() and (work / "clustered.list").read_text() == listed):
        (work / "clustered.list").write_text(listed)      # a triangle is never reused for another list of genomes
        subprocess.run([a.skani, "triangle", "-t", str(a.threads), "-l", str(work / "clustered.list"), "-E",
                        "-o", str(triangle)], check=True, stderr=subprocess.DEVNULL)
    ani = {}
    for r in csv.DictReader(open(triangle), delimiter="\t"):
        ani[(r["Ref_file"], r["Query_file"])] = ani[(r["Query_file"], r["Ref_file"])] = float(r["ANI"])
    representatives, member_of = [], {}
    for p in paths:            # best first; join the first representative at >= 95%, else start a new lineage
        rep = next((q for q in representatives if ani.get((p, q), 0) >= 95), None)
        if rep is None:
            representatives.append(p)
            rep = p
        member_of[p] = rep

    kind = "other_lineage" if a.other_lineages_only else "lineage"
    rng = random.Random(a.seed)
    negatives = (src / "negatives.txt").read_text()
    rows = []
    for i, rep in enumerate(representatives):
        members = [p for p in paths if member_of[p] == rep]
        name = f"{a.set}_{kind}{i}"
        rows.append(dict(lineage=name, representative=Path(rep).name, members=len(members)))
        if len(members) < a.min_size:
            continue
        held = rng.sample(members, max(1, len(members) // 4))
        out = a.lists / name
        (out / "held_out").mkdir(parents=True, exist_ok=True)
        (out / "targets.txt").write_text("".join(p + "\n" for p in members if p not in held))
        (out / "negatives.txt").write_text(negatives)
        (out / "held_out" / "same_lineage.txt").write_text("".join(p + "\n" for p in held))
        for empty in ("other_lineages.txt", "other_named_species.txt", "other_unnamed.txt"):
            (out / "held_out" / empty).write_text("")
        with open(out / "targets.tsv", "w") as f:
            f.write("accession\tspecies\theld_out\ttarget\tpath\n")
            f.writelines(f"{by_path[p]['accession']}\t{by_path[p]['species']}\t{'same_lineage' if p in held else ''}\t"
                         f"{int(p not in held)}\t{p}\n" for p in members)
        shutil.copy(src / "reference_genome.tsv", out / "reference_genome.tsv")
    with open(work / "lineages.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(work / "members.tsv", "w") as f:
        f.write("genome\tlineage\n")
        f.writelines(f"{p}\t{a.set}_{kind}{representatives.index(member_of[p])}\n" for p in paths)
    for r in rows:
        print(r, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

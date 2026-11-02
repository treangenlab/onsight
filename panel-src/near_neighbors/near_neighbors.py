#!/usr/bin/env python3
"""Near-neighbor genomes of each panel target, for the near-neighbor check (taxonomy and ANI only, nothing fitted to
reads).

Targets: config rows that share a target name are one target. A name that is the genus of its sets (Enterobacter) is
a genus target, otherwise a species target. A target's own genomes are its sets' NCBI reference genomes (the config's
first; final_sets.py manifests), plus the representative of each lineage set.

Rule, per target:
  candidates  background genomes (non-viral; one RefSeq genome per NCBI species taxid, the all-genome background pack)
              of the target's NCBI family, plus any background genome at >= 80% ANI (skani) to one of its own genomes;
  dropped     the target itself: the same species taxid (a genus target: any species of the genus) or >= 95% ANI to
              one of its own genomes; and every other target of the keep-off panel, judged the same way (the config column
              keep_off names that panel; empty: this panel);
  grouped     candidates at >= 95% ANI with >= 50% aligned (skani triangle, joined transitively) are one species. Each
              group has one representative: a named species before an "sp.", then RefSeq "reference genome", then
              Complete Genome, then accession order.
skani reports nothing below about 80% ANI; such genomes get "< 80".

Subcommands, in order:
  sources     targets.tsv (target, name, taxid, rank, sets), keep_off.tsv, sources.tsv (each target's own genomes),
              sources.txt
  ani         skani_sources.tsv: skani dist of every background genome against the own genomes (-s 80)
  candidates  candidates.tsv and candidates_<target>.txt
  triangles   triangle_<target>.tsv: skani triangle of each target's candidates
  groups      near_neighbors.tsv (one row per representative), near_neighbor_groups.tsv (every member), genomes.tsv
              and taxids.tsv (the homolog folder lists for blast_homologs.py; near_neighbor_of = the target's sets)
A skani table is reused only when it was made from the same genome lists (saved next to it).

    conda activate squiggle
    python3 near_neighbors.py sources --config panels/<panel>.tsv --sets <root>/sets --taxdump <dir> --out <dir>
    python3 near_neighbors.py ani --background <genomes.tsv> --out <dir> [--skani skani --threads 4]
    python3 near_neighbors.py candidates --background <genomes.tsv> --taxdump <dir> --out <dir>
    python3 near_neighbors.py triangles --out <dir> [--skani skani --threads 4]
    python3 near_neighbors.py groups --background <genomes.tsv> --out <dir>
"""
import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

NEAR, SAME = 80.0, 95.0
ACC = re.compile(r"(GC[AF]_\d+\.\d+)")


def acc(path: str) -> str:
    return ACC.search(Path(path).name).group(1)


def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def write_tsv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def rank(g: dict) -> tuple:
    return (" sp." in g["organism"], g["refseq_category"] != "reference genome",
            g["assembly_level"] != "Complete Genome", g["genome"])


def made_from(out: Path, inputs: dict) -> bool:
    """True when the skani table out exists and was made from these inputs (suffix -> text, saved next to it).
    Otherwise the inputs are saved and False is returned: a table is never reused for another list of genomes."""
    if out.exists() and all(out.with_suffix(k).exists() and out.with_suffix(k).read_text() == t for k, t in inputs.items()):
        return True
    for k, t in inputs.items():
        out.with_suffix(k).write_text(t)
    return False


class Taxonomy:
    """NCBI taxonomy dump: the ancestor of a given rank, and taxids of scientific names."""

    def __init__(self, taxdump: Path):
        self.parent, self.rank = {}, {}
        for line in open(taxdump / "nodes.dmp"):
            f = line.split("\t|\t", 3)
            self.parent[int(f[0])], self.rank[int(f[0])] = int(f[1]), f[2]
        self.merged = {int(f[0]): int(f[1].split("\t")[0])
                       for f in (l.split("\t|\t") for l in open(taxdump / "merged.dmp"))}
        self.taxdump = taxdump

    def up(self, t: int, want: str) -> int:
        t = self.merged.get(t, t)
        while t in self.parent and t != 1:
            if self.rank[t] == want:
                return t
            t = self.parent[t]
        return 0

    def ids(self, names: set) -> dict:
        out = {}
        for line in open(self.taxdump / "names.dmp"):
            f = line.split("\t|\t")
            if f[1] in names and f[3].startswith("scientific name"):
                out.setdefault(f[1], int(f[0]))
        return out


def targets_of(rows: list) -> dict:
    """target -> (name, scientific name, rank, sets) from config rows: rows that share a name are one target."""
    out = {}
    for r in rows:
        genus, sp = r["set"].split("_")[:2]
        if r["name"] == genus:                       # a genus target (Enterobacter)
            key, sci, rk = genus, genus, "genus"
        else:
            key, sci, rk = f"{genus}_{sp}", f"{genus} {sp}", "species"
        out.setdefault(key, [r["name"], sci, rk, []])[3].append(r["set"])
    return out


def sources(a) -> None:
    rows_cfg = read_tsv(a.config)
    ent = targets_of(rows_cfg)
    tax = Taxonomy(a.taxdump)
    ids = tax.ids({v[1] for v in ent.values()})
    rows, erows = [], []
    for key, (name, sci, rk, sets) in ent.items():
        erows.append(dict(target=key, name=name, scientific_name=sci, taxid=ids[sci], rank=rk, sets=",".join(sets)))
        for s in sets:
            m = next(csv.DictReader(open(a.sets / s / "manifest.tsv"), delimiter="\t"))
            lists = Path(m["lists"])            # the set's genome lists (select_genomes.py, lineage_lists.py)
            first = m["reference_genome"].split(";")[0]      # the config's first reference genome of the set
            refs = list(csv.DictReader(open(lists / "reference_genome.tsv"), delimiter="\t"))
            # empty in the config: the reference genome the genome selection used for its ANI rule
            r = next((r for r in refs if r["accession"] == first) if first else
                     (r for r in refs if r.get("used_for_ani") == "yes"), None) or (None if first else refs[0])
            if r is None:
                raise SystemExit(f"{s}: reference genome {first} is not in {lists / 'reference_genome.tsv'}")
            rows.append(dict(target=key, set=s, assembly=r["accession"], file=r["path"], why="reference genome"))
            if m["lineage_representative"]:          # lineage set: its representative too
                g = m["lineage_representative"].split()[0]
                hit = [t["path"] for t in csv.DictReader(open(lists / "targets.tsv"), delimiter="\t")
                       if t["accession"] == g]
                rows.append(dict(target=key, set=s, assembly=g, file=hit[0], why="lineage representative"))
    # the other targets kept off each target's near-neighbor side: those of the panel its config row names
    keep = []
    for key, (name, sci, rk, sets) in ent.items():
        panel = next(r["keep_off"] for r in rows_cfg if r["set"] == sets[0])
        others = targets_of(read_tsv(a.config.parent / f"{panel}.tsv")) if panel else ent
        missing = set(others) - set(ent)
        if missing:
            raise SystemExit(f"{key}: keep-off targets {sorted(missing)} are not targets of this panel")
        keep.append(dict(target=key, keep_off=",".join(k for k in others if k != key)))
    write_tsv(a.out / "targets.tsv", erows)
    write_tsv(a.out / "keep_off.tsv", keep)
    write_tsv(a.out / "sources.tsv", rows)
    (a.out / "sources.txt").write_text("".join(sorted({r["file"] + "\n" for r in rows})))
    for e in erows:
        print(e["target"], e["taxid"], e["rank"], len(e["sets"].split(",")), "sets",
              sum(1 for r in rows if r["target"] == e["target"]), "own genomes", flush=True)


def ani(a) -> None:
    queries = "".join(g["file"] + "\n" for g in read_tsv(a.background) if g["kingdom"] != "viral")
    out = a.out / "skani_sources.tsv"
    if not made_from(out, {".queries": queries, ".sources": (a.out / "sources.txt").read_text()}):
        subprocess.run([a.skani, "dist", "--ql", str(out.with_suffix(".queries")), "--rl",
                        str(out.with_suffix(".sources")), "-t", str(a.threads), "-s", "80", "-o", str(out)], check=True)
    print("skani table", out, sum(1 for _ in open(out)) - 1, "rows", flush=True)


def candidates(a) -> None:
    ent = read_tsv(a.out / "targets.tsv")
    keep_off = {r["target"]: set(filter(None, r["keep_off"].split(","))) for r in read_tsv(a.out / "keep_off.tsv")}
    src = read_tsv(a.out / "sources.tsv")
    of = defaultdict(set)
    for r in src:
        of[r["file"]].add(r["target"])
    ani_of = {}                              # (genome, target) -> largest ANI to the target's own genomes
    for r in read_tsv(a.out / "skani_sources.tsv"):
        for e in of[r["Ref_file"]]:
            k = (acc(r["Query_file"]), e)
            ani_of[k] = max(ani_of.get(k, 0.0), float(r["ANI"]))
    tax = Taxonomy(a.taxdump)
    genomes = [g for g in read_tsv(a.background) if g["kingdom"] != "viral"]
    fam = {t: tax.up(t, "family") for t in {int(g["species_taxid"]) for g in genomes}}
    gen = {t: tax.up(t, "genus") for t in fam}

    def own(e, t, g):
        v = ani_of.get((g, e["target"]))
        same = gen[t] == int(e["taxid"]) if e["rank"] == "genus" else t == int(e["taxid"])
        return same or (v is not None and v >= SAME)

    rows = []
    for e in ent:
        efam = tax.up(int(e["taxid"]), "family")
        files = []
        for g in genomes:
            t = int(g["species_taxid"])
            v = ani_of.get((g["genome"], e["target"]))
            by_family, by_ani = fam[t] == efam != 0, v is not None and v >= NEAR
            if not (by_family or by_ani):
                continue
            other = [o["target"] for o in ent if o["target"] in keep_off[e["target"]] and own(o, t, g["genome"])]
            dropped = "same species" if own(e, t, g["genome"]) else ("other panel species " + "+".join(other)
                                                                     if other else "")
            rows.append(dict(target=e["target"], genome=g["genome"], organism=g["organism"], genus_taxid=gen[t],
                             family_taxid=fam[t], species_taxid=t, ani=v if v is not None else "< 80",
                             why="+".join(w for w, ok in (("family", by_family), ("ani", by_ani)) if ok),
                             dropped=dropped))
            if not dropped:
                files.append(g["file"])
        (a.out / f"candidates_{e['target']}.txt").write_text("".join(f + "\n" for f in files))
        mine = [r for r in rows if r["target"] == e["target"]]
        print(e["target"], "family", efam, "candidates", len(files), "dropped as the target itself",
              sum(1 for r in mine if r["dropped"] == "same species"), "as another panel target",
              sum(1 for r in mine if r["dropped"].startswith("other")), flush=True)
    write_tsv(a.out / "candidates.tsv", rows)


def triangles(a) -> None:
    for e in read_tsv(a.out / "targets.tsv"):
        listed = (a.out / f"candidates_{e['target']}.txt").read_text()
        out = a.out / f"triangle_{e['target']}.tsv"
        if listed.count("\n") < 2:                   # nothing to pair
            out.write_text("")
            continue
        if not made_from(out, {".list": listed}):
            subprocess.run([a.skani, "triangle", "-l", str(out.with_suffix(".list")), "-E", "-t", str(a.threads),
                            "-o", str(out)], check=True)
        print(e["target"], "candidates", listed.count("\n"), flush=True)


def groups(a) -> None:
    ent = {e["target"]: e for e in read_tsv(a.out / "targets.tsv")}
    meta = {g["genome"]: g for g in read_tsv(a.background)}
    cand = [r for r in read_tsv(a.out / "candidates.tsv") if not r["dropped"]]
    reps, members, listed = [], [], {}
    for sp in ent:
        mine = [r for r in cand if r["target"] == sp]
        link = defaultdict(set)
        tri = a.out / f"triangle_{sp}.tsv"
        for r in (read_tsv(tri) if tri.exists() and tri.stat().st_size else []):
            if float(r["ANI"]) >= SAME and max(float(r["Align_fraction_ref"]), float(r["Align_fraction_query"])) >= 50:
                x, y = acc(r["Ref_file"]), acc(r["Query_file"])
                link[x].add(y)
                link[y].add(x)
        seen = set()
        row = {x["genome"]: x for x in mine}
        for r in sorted(mine, key=lambda r: rank(meta[r["genome"]])):
            if r["genome"] in seen:
                continue
            group, stack = set(), [r["genome"]]
            while stack:
                x = stack.pop()
                if x not in group:
                    group.add(x)
                    stack.extend((link[x] & set(row)) - group)
            seen |= group
            rep = min(group, key=lambda g: rank(meta[g]))
            reps.append(dict(target=sp, genome=rep, organism=meta[rep]["organism"],
                             species_taxid=meta[rep]["species_taxid"], ani=row[rep]["ani"], why=row[rep]["why"],
                             group_size=len(group),
                             group_members=";".join(f"{g} {meta[g]['organism']}" for g in sorted(group) if g != rep)))
            members += [dict(target=sp, representative=rep, genome=g, organism=meta[g]["organism"],
                             species_taxid=meta[g]["species_taxid"], ani=row[g]["ani"]) for g in sorted(group)]
            listed.setdefault(rep, (meta[rep]["file"], meta[rep]["species_taxid"], []))[2].extend(
                ent[sp]["sets"].split(","))
        print(sp, "candidates", len(mine), "species groups", sum(1 for x in reps if x["target"] == sp), flush=True)
    write_tsv(a.out / "near_neighbors.tsv", reps)
    write_tsv(a.out / "near_neighbor_groups.tsv", members)
    with open(a.out / "genomes.tsv", "w") as f:
        f.writelines(f"{g}\t{v[0]}\n" for g, v in sorted(listed.items()))
    with open(a.out / "taxids.tsv", "w") as f:
        f.write("genome\tspecies_taxid\tnear_neighbor_of\n")
        f.writelines(f"{g}\t{v[1]}\t{','.join(v[2])}\n" for g, v in sorted(listed.items()))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=("sources", "ani", "candidates", "triangles", "groups"))
    ap.add_argument("--config", type=Path, help="panel config (panels/<panel>.tsv)")
    ap.add_argument("--sets", type=Path, help="final sets folder (final/<set>/manifest.tsv, final_sets.py)")
    ap.add_argument("--background", type=Path, help="genomes.tsv of the background pack")
    ap.add_argument("--taxdump", type=Path)
    ap.add_argument("--skani", default="skani")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    {"sources": sources, "ani": ani, "candidates": candidates, "triangles": triangles, "groups": groups}[a.what](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

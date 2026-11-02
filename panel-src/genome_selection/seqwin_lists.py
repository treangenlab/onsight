#!/usr/bin/env python3
"""Seqwin target and negative lists by NCBI label from one NCBI query, and the genome download.

Assemblies: every GenBank assembly under --scope, with Seqwin's own download filters (atypical genomes and large
multi-isolate projects excluded, MAGs kept, all assembly levels). Each assembly's species and genus come from NCBI
Taxonomy (the local dump TAXDUMP; NCBI's API only for tax ids it does not know).
- Targets: the assemblies under --target (a species or a genus).
- Negatives: the assemblies of every other named species under --scope ("Genus sp." and genus-only labels are left
  out, as in the March and May negative lists), at most --neg-cap per species (random, --seed) when given.

Writes in --out:
  assemblies.tsv  one row per assembly: labels, NCBI reference genome or not, role, genome path
  species.tsv     assemblies per species and role
  targets.txt, negatives.txt   one genome path per line (Seqwin's tar_paths and neg_paths)
Genomes go to --genomes as one NCBI datasets package (gzip FASTA).

    conda activate seqwin
    python3 seqwin_lists.py --target "Akkermansia muciniphila" --scope Akkermansia --neg-cap 50 --taxdump <dir> \
        --out <dir> --genomes <dir>
"""
import argparse
import csv
import json
import random
import subprocess
import tempfile
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

FILTERS = ["--assembly-source", "GenBank", "--exclude-atypical", "--exclude-multi-isolate", "--mag", "all"]
TAXDUMP = Path("taxdump")       # NCBI taxonomy dump folder (--taxdump); without one, NCBI's API is asked for every tax id


def reports(*args):
    p = subprocess.Popen(["datasets", *args, "--as-json-lines"], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        if line.strip():
            yield json.loads(line)
    if p.wait():
        raise RuntimeError(f"datasets {' '.join(args)} exited {p.returncode}")


def classify_local(taxids) -> dict:
    """classify() from the local NCBI Taxonomy dump (TAXDUMP); tax ids it does not know are left out."""
    parent, rank, name = {}, {}, {}
    for line in open(TAXDUMP / "nodes.dmp"):
        f = line.split("\t|\t", 3)
        parent[int(f[0])], rank[int(f[0])] = int(f[1]), f[2]
    for line in open(TAXDUMP / "names.dmp"):
        f = line.split("\t|\t")
        if f[3].startswith("scientific name"):
            name[int(f[0])] = f[1]
    merged = {int(f[0]): int(f[1].split("\t")[0]) for f in (line.split("\t|\t") for line in open(TAXDUMP / "merged.dmp"))}
    out = {}
    for tid in set(taxids):
        t, lab = merged.get(tid, tid), {}
        while t in rank and "genus" not in lab:
            if rank[t] in ("species", "genus"):
                lab.setdefault(rank[t], (t, name.get(t, "")))
            if parent[t] == t:
                break
            t = parent[t]
        if t in rank:
            out[tid] = lab
    return out


def classify(taxids) -> dict:
    """Tax id -> {"species": (id, name), "genus": (id, name)}, as far as NCBI Taxonomy has them. Uses the local dump
    when present and asks NCBI only for tax ids the dump does not know."""
    taxids = set(taxids)
    out = classify_local(taxids) if (TAXDUMP / "nodes.dmp").exists() else {}
    ids = sorted(taxids - set(out))
    for i in range(0, len(ids), 1000):
        with tempfile.NamedTemporaryFile("w", suffix=".txt") as f:
            f.write("\n".join(map(str, ids[i:i + 1000])) + "\n")
            f.flush()
            try:
                for r in reports("summary", "taxonomy", "taxon", "--inputfile", f.name):
                    t = r.get("taxonomy") or {}
                    c = t.get("classification") or {}
                    if "tax_id" in t:
                        out[t["tax_id"]] = {k: (c[k]["id"], c[k]["name"]) for k in ("species", "genus") if k in c}
            except RuntimeError as e:
                print(f"skipped tax ids NCBI does not know: {e}", flush=True)
    unknown = sorted(taxids - set(out))
    if unknown:
        print(f"{len(unknown)} tax ids without a label, their assemblies are left out: {unknown[:10]}", flush=True)
    return out


def named(species: str) -> bool:
    return " sp." not in species and not species.startswith(("uncultured", "unclassified"))


def download(accessions: list, folder: Path, workers: int) -> dict:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "accessions.txt").write_text("\n".join(accessions) + "\n")
    package = folder / "package.zip"
    package.unlink(missing_ok=True)
    subprocess.run(["datasets", "download", "genome", "accession", "--inputfile", str(folder / "accessions.txt"),
                    "--include", "genome", "--dehydrated", "--no-progressbar", "--filename", str(package)], check=True)
    with zipfile.ZipFile(package) as z:
        z.extractall(folder)
    subprocess.run(["datasets", "rehydrate", "--directory", str(folder), "--gzip", "--max-workers", str(workers),
                    "--no-progressbar"], check=True)
    paths = {}
    for acc in accessions:
        files = sorted((folder / "ncbi_dataset" / "data" / acc).glob("*.fna*"))
        if files:
            paths[acc] = files[0]
    return paths


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--scope", required=True)
    ap.add_argument("--neg-cap", type=int)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--genomes", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--taxdump", type=Path, help="NCBI taxonomy dump folder (nodes.dmp, names.dmp, merged.dmp)")
    a = ap.parse_args()
    global TAXDUMP
    TAXDUMP = a.taxdump or TAXDUMP
    a.out.mkdir(parents=True, exist_ok=True)

    target = next(reports("summary", "taxonomy", "taxon", a.target))["taxonomy"]
    rank, target_id = target["rank"].lower(), target["tax_id"]
    if rank not in ("species", "genus"):
        raise SystemExit(f"{a.target}: rank {rank}, expected species or genus")
    rows = [dict(accession=r["accession"], tax_id=r["organism"]["tax_id"], organism=r["organism"]["organism_name"],
                 level=r.get("assembly_info", {}).get("assembly_level", ""))
            for r in reports("summary", "genome", "taxon", a.scope, *FILTERS)]
    refs = set()
    for r in reports("summary", "genome", "taxon", a.scope, "--reference"):
        refs.update(x for x in (r.get("accession"), r.get("paired_accession")) if x)
    labels = classify(r["tax_id"] for r in rows)
    print(f"{a.target} ({rank} {target_id}) in {a.scope}: {len(rows)} assemblies, {len(labels)} tax ids classified", flush=True)

    by_species = defaultdict(list)
    for r in rows:
        lab = labels.get(r["tax_id"], {})
        r["species_id"], r["species"] = lab.get("species", ("", ""))
        r["genus_id"], r["genus"] = lab.get("genus", ("", ""))
        r["reference"] = "yes" if r["accession"] in refs else "no"
        if (r["species_id"] if rank == "species" else r["genus_id"]) == target_id:
            r["role"] = "target"
        elif not r["species"] or not named(r["species"]):
            r["role"] = "dropped_unnamed"
        else:
            r["role"] = "negative"
            by_species[r["species_id"]].append(r)
    rng = random.Random(a.seed)
    for sid in sorted(by_species):
        group = sorted(by_species[sid], key=lambda r: r["accession"])
        if a.neg_cap and len(group) > a.neg_cap:
            keep = {r["accession"] for r in rng.sample(group, a.neg_cap)}
            for r in group:
                if r["accession"] not in keep:
                    r["role"] = "dropped_cap"

    chosen = [r["accession"] for r in rows if r["role"] in ("target", "negative")]
    print("to download:", dict(Counter(r["role"] for r in rows)), flush=True)
    paths = download(chosen, a.genomes, a.workers)
    for r in rows:
        r["path"] = str(paths.get(r["accession"], ""))
    missing = [acc for acc in chosen if acc not in paths]
    if missing:
        print(f"{len(missing)} chosen assemblies have no genome file, e.g. {missing[:5]}", flush=True)

    cols = ["accession", "tax_id", "organism", "species_id", "species", "genus_id", "genus", "level", "reference",
            "role", "path"]
    with open(a.out / "assemblies.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    counts = Counter((r["species_id"], r["species"], r["genus"], r["role"]) for r in rows)
    with open(a.out / "species.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["species_id", "species", "genus", "role", "assemblies"])
        w.writerows([*k, n] for k, n in sorted(counts.items(), key=lambda kv: (kv[0][3], -kv[1])))
    for role, name in (("target", "targets.txt"), ("negative", "negatives.txt")):
        (a.out / name).write_text("".join(f"{r['path']}\n" for r in rows if r["role"] == role and r["path"]))
    roles = Counter(r["role"] for r in rows)
    print("roles:", dict(roles), "| species in negatives:", len({r["species_id"] for r in rows if r["role"] == "negative"}))
    if not roles["target"] or not roles["negative"]:
        raise SystemExit("no targets or no negatives")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

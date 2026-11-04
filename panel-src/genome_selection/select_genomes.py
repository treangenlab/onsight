#!/usr/bin/env python3
"""Seqwin's target and negative genomes for one signature set: genome selection with the 95% ANI species rule.

Input: assemblies.tsv of seqwin_lists.py (one row per assembly: NCBI species label, assembly level, NCBI reference
genome or not, role, genome path).
- Targets: the assemblies labelled with the target species that are high quality and at >= 95% ANI (skani) to the
  species' NCBI reference genome. High quality is MIMAG high quality by NCBI's CheckM (completeness > 90%,
  contamination < 5%); an assembly with no CheckM value at NCBI counts when it is at scaffold level or better.
- Species rule: >= 95% ANI with >= 50% of both genomes aligned to a named species' NCBI reference genome.
  - A target-labelled assembly closer to another named species' reference genome than to its own is that species:
    it becomes a negative (at most --neg-cap per species, random, seed + 1).
  - A negative that carries the target's genome (>= 95% ANI over >= 50% of the target's reference genome, also a
    hybrid) is left out.
- Negatives: every other named species ("Genus sp." and genus-only labels left out), at most --neg-cap per species
  (random, --seed). With --target, only assemblies with a genome path (downloaded).
- Held-out genomes, put aside at random (--seed) before the targets are chosen, never targets:
  - reference lineage: 25% of the high-quality assemblies at >= 95% ANI to the reference genome, at least 3 and at
    most 20 (none when fewer than 6);
  - other lineages: 25% of the high-quality assemblies under the label at < 95% ANI, at least 1 and at most 20.
- At most --target-cap targets per species (random, sorted by accession).

Random draws, in this order from one generator: the negative cap, then the held-out genomes, then the target cap.

Writes in --out/<set>/:
  targets.txt, negatives.txt     one genome path per line (Seqwin's tar_paths and neg_paths)
  held_out/same_lineage.txt      the held-out genomes of the reference lineage; held_out/other_lineages.txt those of
                                 other lineages (other_named_species.txt: those the species rule gives to another
                                 named species; other_unnamed.txt: the rest)
  targets.tsv                    every target-labelled assembly: CheckM, ANI to the reference genome, quality,
                                 closer named species, held out, target, genome path
  reference_genome.tsv           every NCBI reference genome of each target species (species, accession, path);
                                 used_for_ani marks the one the ANI rule uses (the first in assemblies.tsv)
  checkm.tsv                     the CheckM values used; skani/ the skani tables (kept, so a rerun reuses them)

    conda activate seqwin       # NCBI datasets; skani 0.4.0 by --skani
    python3 select_genomes.py --assemblies <assemblies.tsv> --set <name> [--target "<species>"] \
        [--neg-cap 50] [--target-cap 1000] [--checkm <saved checkm.tsv>] --out <dir>
"""
import argparse
import csv
import json
import random
import subprocess
from pathlib import Path

SCAFFOLD_OR_BETTER = {"Complete Genome", "Chromosome", "Scaffold"}


def named(species: str) -> bool:
    return bool(species) and " sp." not in species and not species.startswith("uncultured") and len(species.split()) >= 2


def checkm_values(accessions: list, saved: Path, out: Path) -> dict:
    """accession -> (completeness, contamination), None when NCBI has no CheckM value. Values in --checkm are used as
    they are; the rest come from NCBI datasets. The values used are written to out."""
    have = {r["accession"]: r for r in csv.DictReader(open(saved), delimiter="\t")} if saved else {}
    todo = [a for a in accessions if a not in have]
    for i in range(0, len(todo), 1000):
        query = out.parent / "checkm_query.txt"
        query.write_text("\n".join(todo[i:i + 1000]) + "\n")
        report = subprocess.run(["datasets", "summary", "genome", "accession", "--inputfile", str(query),
                                 "--as-json-lines"], capture_output=True, text=True, check=True).stdout
        for line in report.splitlines():
            if line.startswith("{"):
                d = json.loads(line)
                c = d.get("checkm_info")          # NCBI leaves out fields that are zero
                have[d["accession"]] = dict(accession=d["accession"],
                                            completeness=c.get("completeness", 0.0) if c else "",
                                            contamination=c.get("contamination", 0.0) if c else "")
        for a in todo[i:i + 1000]:
            have.setdefault(a, dict(accession=a, completeness="", contamination=""))
        query.unlink()
    with open(out, "w") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "completeness", "contamination"], delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(have[a] for a in accessions)
    return {a: (float(have[a]["completeness"]) if have[a]["completeness"] != "" else None,
                float(have[a]["contamination"]) if have[a]["contamination"] != "" else None) for a in accessions}


def made_from(out: Path, inputs: dict) -> bool:
    """True when the skani table out exists and was made from these inputs (suffix -> text, saved next to it).
    Otherwise the inputs are saved and False is returned: a table is never reused for another list of genomes."""
    if out.exists() and all(out.with_suffix(k).exists() and out.with_suffix(k).read_text() == t for k, t in inputs.items()):
        return True
    for k, t in inputs.items():
        out.with_suffix(k).write_text(t)
    return False


def ani_to_reference(skani: str, reference: str, paths: list, out: Path, threads: int) -> dict:
    """genome path -> ANI to the reference genome (skani reports nothing below about 80%, taken as 0)."""
    if not made_from(out, {".list": "\n".join(paths) + "\n", ".reference": reference + "\n"}):
        subprocess.run([skani, "dist", "-t", str(threads), "-q", reference, "--rl", str(out.with_suffix(".list")),
                        "-o", str(out)], check=True, stderr=subprocess.DEVNULL)
    return {r["Ref_file"]: float(r["ANI"]) for r in csv.DictReader(open(out), delimiter="\t")}


def closest_named_species(skani: str, genomes: list, references: dict, out: Path, threads: int,
                          reference_side_only: bool = False) -> dict:
    """genome path -> (species, ANI) of the named species' reference genome it is closest to, by the species rule
    (>= 95% ANI, >= 50% of both genomes aligned). reference_side_only: >= 50% of the reference genome aligned,
    whatever the share of the genome itself (a hybrid or contaminated assembly that carries the species' genome)."""
    if not genomes or not references:
        return {}
    if not made_from(out, {".genomes": "\n".join(genomes) + "\n", ".references": "\n".join(references.values()) + "\n"}):
        subprocess.run([skani, "dist", "-t", str(threads), "--ql", str(out.with_suffix(".references")), "--rl",
                        str(out.with_suffix(".genomes")), "-o", str(out)], check=True, stderr=subprocess.DEVNULL)
    species_of = {p: sp for sp, p in references.items()}
    best = {}
    for r in csv.DictReader(open(out), delimiter="\t"):        # query: a reference genome; ref: a genome
        ani = float(r["ANI"])
        aligned = float(r["Align_fraction_query"]) if reference_side_only else \
            min(float(r["Align_fraction_ref"]), float(r["Align_fraction_query"]))
        if ani >= 95 and aligned >= 50 and ani > best.get(r["Ref_file"], ("", 0))[1]:
            best[r["Ref_file"]] = (species_of[r["Query_file"]], ani)
    return best


def cap_per_species(rows: list, cap: int, rng: random.Random, key: str = "species") -> list:
    """At most cap rows per species (random), species in sorted order."""
    by = {}
    for r in rows:
        by.setdefault(r[key], []).append(r)
    return [r for sp in sorted(by) for r in (rng.sample(by[sp], cap) if len(by[sp]) > cap else by[sp])]


def write_paths(path: Path, rows: list) -> None:
    path.write_text("".join(r["path"] + "\n" for r in rows))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assemblies", type=Path, required=True, help="assemblies.tsv of seqwin_lists.py")
    ap.add_argument("--set", required=True, help="set name, the output folder")
    ap.add_argument("--target", help="target species, comma separated (default: the assemblies with role target); "
                                     "negatives are then every other named species")
    ap.add_argument("--neg-cap", type=int, help="at most this many negatives per species (random)")
    ap.add_argument("--target-cap", type=int, help="at most this many targets per species (random)")
    ap.add_argument("--checkm", type=Path, help="saved CheckM table (accession, completeness, contamination); "
                                                "accessions not in it are looked up at NCBI")
    ap.add_argument("--skani", default="skani")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    out = a.out / a.set
    (out / "skani").mkdir(parents=True, exist_ok=True)
    (out / "held_out").mkdir(exist_ok=True)

    rows = list(csv.DictReader(open(a.assemblies), delimiter="\t"))
    if a.target:
        want = set(a.target.split(","))
        targets = [r for r in rows if r["species"] in want]
        # downloaded genomes only: seqwin_lists.py capped each species and left the rest without a path
        negatives = [r for r in rows if r["species"] not in want and named(r["species"]) and r["path"]]
    else:
        targets = [r for r in rows if r["role"] == "target"]
        negatives = [r for r in rows if r["role"] == "negative"]
    rng = random.Random(a.seed)
    if a.neg_cap:
        negatives = cap_per_species(negatives, a.neg_cap, rng)

    checkm = checkm_values([r["accession"] for r in targets], a.checkm, out / "checkm.tsv")
    ani, reference_rows = {}, []
    for sp in sorted({r["species"] for r in targets}):
        mine = [r for r in targets if r["species"] == sp]
        reference = [r for r in rows if r["species"] == sp and r["reference"] == "yes"]
        if not reference:
            print(f"{sp}: no NCBI reference genome; none of its assemblies can be a target")
            continue
        reference_rows += [(sp, r["accession"], r["path"], "yes" if i == 0 else "no") for i, r in enumerate(reference)]
        got = ani_to_reference(a.skani, reference[0]["path"], [r["path"] for r in mine],
                               out / "skani" / f"ani_to_reference_{sp.replace(' ', '_')}.tsv", a.threads)
        ani.update({r["accession"]: got.get(r["path"], 0.0) for r in mine})
    with open(out / "reference_genome.tsv", "w") as f:
        f.write("species\taccession\tpath\tused_for_ani\n" + "".join("\t".join(r) + "\n" for r in reference_rows))
    for r in targets:
        r["completeness"], r["contamination"] = checkm[r["accession"]]
        r["ani_to_reference"] = ani.get(r["accession"], 0.0)
        r["high_quality"] = (r["completeness"] > 90 and r["contamination"] < 5 if r["completeness"] is not None
                             else r["level"] in SCAFFOLD_OR_BETTER)
        r["reference_lineage"] = r["ani_to_reference"] >= 95

    # the species rule, against every named species' NCBI reference genome
    references = {}
    for r in rows:
        if r["reference"] == "yes" and named(r["species"]):
            references.setdefault(r["species"], r["path"])
    target_species = {r["species"] for r in targets}
    closer = closest_named_species(a.skani, [r["path"] for r in targets],
                                   {sp: p for sp, p in references.items() if sp not in target_species},
                                   out / "skani" / "targets_to_other_references.tsv", a.threads)
    for r in targets:
        sp, ani_other = closer.get(r["path"], ("", 0))
        r["closer_named_species"] = sp if sp and ani_other > r["ani_to_reference"] else ""
    carry_target = closest_named_species(a.skani, [r["path"] for r in negatives],
                                         {sp: p for sp, p in references.items() if sp in target_species},
                                         out / "skani" / "negatives_to_target_references.tsv", a.threads,
                                         reference_side_only=True)

    # held-out genomes, put aside before the targets are chosen
    lineage_pool = [r for r in targets if r["high_quality"] and r["reference_lineage"] and not r["closer_named_species"]]
    other_pool = [r for r in targets if r["high_quality"] and not r["reference_lineage"]]
    held_lineage = (rng.sample(lineage_pool, min(len(lineage_pool), max(3, min(20, len(lineage_pool) // 4))))
                    if len(lineage_pool) >= 6 else [])
    held_other = rng.sample(other_pool, min(len(other_pool), 20, max(1, len(other_pool) // 4))) if other_pool else []
    held = {r["accession"] for r in held_lineage + held_other}
    pool = [r for r in targets if r["accession"] not in held]
    chosen = [r for r in pool if r["high_quality"] and r["reference_lineage"] and not r["closer_named_species"]]
    if a.target_cap:
        by = {}
        for r in chosen:
            by.setdefault(r["species"], []).append(r)
        chosen = [r for sp in sorted(by) for r in (sorted(rng.sample(by[sp], a.target_cap), key=lambda x: x["accession"])
                                                   if len(by[sp]) > a.target_cap else by[sp])]
    moved = [dict(r, species_moved_to=r["closer_named_species"]) for r in pool if r["closer_named_species"]]
    if a.neg_cap:          # the moved assemblies are negatives of the species they match: same cap
        moved = cap_per_species(moved, a.neg_cap, random.Random(a.seed + 1), key="species_moved_to")
    negatives = [r for r in negatives if r["path"] not in carry_target] + moved

    write_paths(out / "targets.txt", chosen)
    write_paths(out / "negatives.txt", negatives)
    write_paths(out / "held_out" / "same_lineage.txt", held_lineage)
    write_paths(out / "held_out" / "other_lineages.txt", held_other)
    write_paths(out / "held_out" / "other_named_species.txt", [r for r in held_other if r["closer_named_species"]])
    write_paths(out / "held_out" / "other_unnamed.txt", [r for r in held_other if not r["closer_named_species"]])
    chosen_ids = {r["accession"] for r in chosen}
    cols = ["accession", "species", "level", "completeness", "contamination", "ani_to_reference", "high_quality",
            "reference_lineage", "closer_named_species"]
    with open(out / "targets.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(cols + ["held_out", "target", "path"])
        for r in targets:
            kind = "same_lineage" if r in held_lineage else "other_lineage" if r in held_other else ""
            w.writerow([r[c] for c in cols] + [kind, int(r["accession"] in chosen_ids), r["path"]])
    print(f"{a.set}: {len(chosen)} targets, {len(negatives)} negatives; {len(targets)} labelled, high quality "
          f"{sum(r['high_quality'] for r in targets)}, reference lineage {sum(r['reference_lineage'] for r in targets)}; "
          f"held out {len(held_lineage)} of the reference lineage, {len(held_other)} of other lineages "
          f"({sum(bool(r['closer_named_species']) for r in held_other)} another named species); "
          f"closer to another named species {sum(bool(r['closer_named_species']) for r in targets)}; "
          f"negatives carrying the target's genome {len(carry_target)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build a panel from its config, step by step, ending in the panel folder onsight reads.

--config (panels/<panel>.tsv), one row per final set:
  set, build        the set's name in the panel; the build folder under sets/lists/ and sets/seqwin/
  short, name       the panel table's species and the target name a call is reported under
  input, target, neg_cap, target_cap
                    select_genomes.py: <inputs>/<input>/assemblies.tsv, the target species when the input holds
                    several, and the caps (empty: none)
  lineages_of, lineages
                    for a lineage set: the build it is split from, and "other" (lineages outside the reference
                    genome's lineage) or "all" (lineage_lists.py)
  report_as, panel, mock_label, reference_genome, lineage_representative
                    carried into the final sets' manifests
  keep_off          the panel config whose other targets are kept off this target's near-neighbor side (empty: this
                    panel)
  test_strains      genomes the test reads of this species came from, as <paths key>/<file>, comma separated; an own
                    strain at >= 99.9% ANI to one is skipped
  extra_near_strains
                    a list of read-source genomes (file in inputs): this target gets more strains per near-neighbor
                    species (near_neighbor_strains.py)
--paths (panels/paths.tsv): inputs (folders with assemblies.tsv, checkm.tsv and read-source lists), background,
taxdump, mock_genomes (comma separated), mock_strains, skani, seqwin_python (Seqwin's conda env), refseq_summary,
refseq_libraries (comma separated), the folders that test_strains name, blastn, makeblastdb, onsight (the binary
that builds the indexes); relative paths are read from the paths file's folder, and a bare name (skani, blastn) is a
tool on PATH. mock_genomes and mock_strains are checks only and may be left empty (panels/example_paths.tsv).

Steps (--steps, in this order): lists (select_genomes.py, lineage_lists.py), seqwin (run_seqwin.py), screen
(screen_signatures.py), sets (final_sets.py), panel (build_panel.py, and the panel table panel.tsv), near
(near_neighbors.py), own (own_strains.py), extra (near_neighbor_strains.py), others (other_targets.py: the other
targets' own strains on each target's near-neighbor side), blast (blast_homologs.py on each homolog folder), index
(index/build_index.sh, R9.4.1 and R10.4.1), regions (neighbor_regions.py). Work already there is skipped. Leaving out
the others step (with no work/other_targets/) builds the region file the way the paper's panels were built.

Layout of --root:
  sets/                          per signature set, shared by every panel built in this root
    lists/<build>/ seqwin/<build>/ screen/<build>.* final/<set>/
  <panel>/onsight/               the panel folder onsight reads (--panel-dir): panel.tsv, regions.lcr,
                                 r9/index.idx and r10/index.idx (each with its .sigmeta.tsv)
  <panel>/work/                  signatures/ (the panel's signatures), near_neighbors/, own_strains/,
                                 near_strains/<target>/, other_targets/ (homolog folders: genome lists, blast hits,
                                 paf/)
  <panel>/logs/<step>.log        each step's output

    conda activate squiggle     # minimap2 and skani on PATH for the screen
    python3 run_panel.py --config panels/eskape.tsv --paths panels/paths.tsv --root <dir> [--steps lists,seqwin]
"""
import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = ("lists", "seqwin", "screen", "sets", "panel", "near", "own", "extra", "others", "blast", "index",
         "regions")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--paths", type=Path, required=True)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--steps", default=",".join(STEPS))
    ap.add_argument("--threads", type=int, default=8, help="skani, Seqwin and index threads")
    ap.add_argument("--processes", type=int, default=8, help="minimap2 processes of the screen, blastn processes")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.config), delimiter="\t"))
    paths = {r["key"]: r["value"] for r in csv.DictReader(open(a.paths), delimiter="\t")}
    def at(v: str) -> str:
        """A paths.tsv value: absolute as is; relative to the paths file's folder; a bare tool name from PATH."""
        if v.startswith("/"):
            return v
        p = a.paths.parent / v
        return str(p) if "/" in v or p.exists() else (shutil.which(v) or v)
    inputs = Path(at(paths["inputs"]))
    steps = a.steps.split(",")
    builds = list(dict.fromkeys(r["build"] for r in rows))
    py = sys.executable

    shared = a.root / "sets"
    lists, seqwin, screen, final = (shared / d for d in ("lists", "seqwin", "screen", "final"))
    base = a.root / a.config.stem
    out, work, logs = base / "onsight", base / "work", base / "logs"
    for d in (shared, out, work, logs):
        d.mkdir(parents=True, exist_ok=True)

    def run(step: str, cmd: list) -> None:
        """Runs one command; its output goes to logs/<step>.log."""
        print(f"{step}: {' '.join(str(c) for c in cmd)}", flush=True)
        with open(logs / f"{step}.log", "a") as log:
            log.write("+ " + " ".join(str(c) for c in cmd) + "\n")
            log.flush()
            subprocess.run([str(c) for c in cmd], check=True, stdout=log, stderr=subprocess.STDOUT)

    if "lists" in steps:
        for r in rows:
            if r["lineages_of"] or (lists / r["build"] / "targets.txt").exists():
                continue
            cmd = [py, HERE / "genome_selection/select_genomes.py", "--assemblies",
                   inputs / r["input"] / "assemblies.tsv", "--set", r["build"],
                   "--skani", at(paths["skani"]), "--threads", a.threads, "--out", lists]
            if (inputs / "checkm.tsv").exists():          # saved CheckM values; otherwise NCBI's are fetched
                cmd += ["--checkm", inputs / "checkm.tsv"]
            for opt, col in (("--target", "target"), ("--neg-cap", "neg_cap"), ("--target-cap", "target_cap")):
                if r[col]:
                    cmd += [opt, r[col]]
            run("lists", cmd)
        for parent in dict.fromkeys(r["lineages_of"] for r in rows if r["lineages_of"]):
            mine = [r for r in rows if r["lineages_of"] == parent]
            if all((lists / r["build"] / "targets.txt").exists() for r in mine):
                continue
            cmd = [py, HERE / "genome_selection/lineage_lists.py", "--lists", lists, "--set", parent,
                   "--skani", at(paths["skani"]), "--threads", a.threads]
            if mine[0]["lineages"] == "other":
                cmd.append("--other-lineages-only")
            run("lists", cmd)
    if "seqwin" in steps:
        for b in builds:
            if not (seqwin / b / "signatures.fasta").exists():
                run("seqwin", [at(paths["seqwin_python"]), HERE / "signatures/run_seqwin.py", "--lists", lists / b,
                               "--out", seqwin, "--threads", a.threads])
    if "screen" in steps:
        todo = [b for b in builds if not (screen / f"{b}.tsv").exists()]
        if todo:
            mocks = [at(p) for p in paths.get("mock_genomes", "").split(",") if p]     # checks only; optional
            run("screen", [py, HERE / "signatures/screen_signatures.py", "--root", shared, "--builds", ",".join(todo),
                           "--background", at(paths["background"]), "--taxdump", at(paths["taxdump"]),
                           *(["--mock-genomes", *mocks] if mocks else []),
                           *(["--mock-strains", at(paths["mock_strains"])] if paths.get("mock_strains") else []),
                           "--processes", a.processes, "--out", screen])
    if "sets" in steps:
        run("sets", [py, HERE / "signatures/final_sets.py", "--root", shared, "--config", a.config])
    signatures = work / "signatures"
    if "panel" in steps:
        run("panel", [py, HERE / "signatures/build_panel.py", "--out", signatures, *[final / r["set"] for r in rows]])
        with open(out / "panel.tsv", "w") as f:
            f.write("short\ttaxon\tname\tgenome\n")
            f.writelines(f"{r['short']}\t{r['set']}\t{r['name']}\t-\n" for r in rows)

    near, own, extra, others = (work / d for d in ("near_neighbors", "own_strains", "near_strains", "other_targets"))
    nn = HERE / "near_neighbors"
    skani, libraries = at(paths["skani"]), [at(p) for p in paths.get("refseq_libraries", "").split(",") if p]
    if "near" in steps and not (near / "genomes.tsv").exists():
        common = ["--out", near]
        run("near", [py, nn / "near_neighbors.py", "sources", "--config", a.config, "--sets", final,
                     "--taxdump", at(paths["taxdump"]), *common])
        run("near", [py, nn / "near_neighbors.py", "ani", "--background", at(paths["background"]), "--skani", skani,
                     "--threads", a.threads, *common])
        run("near", [py, nn / "near_neighbors.py", "candidates", "--background", at(paths["background"]),
                     "--taxdump", at(paths["taxdump"]), *common])
        run("near", [py, nn / "near_neighbors.py", "triangles", "--skani", skani, "--threads", a.threads, *common])
        run("near", [py, nn / "near_neighbors.py", "groups", "--background", at(paths["background"]), *common])
    target_of = {}                                  # set -> panel target (targets.tsv of the near step)
    if (near / "targets.tsv").exists():
        for t in csv.DictReader(open(near / "targets.tsv"), delimiter="\t"):
            target_of.update({s: t["target"] for s in t["sets"].split(",")})
    if "own" in steps and not (own / "genomes.tsv").exists():
        own.mkdir(parents=True, exist_ok=True)
        with open(own / "test_strains.tsv", "w") as f:       # species and file of each test strain genome
            f.write("species\tfile\n")
            done = set()
            for r in rows:
                sp = " ".join(r["set"].split("_")[:2])
                for t in filter(None, r["test_strains"].split(",")):
                    key, name = t.split("/", 1)
                    if (sp, t) not in done:
                        f.write(f"{sp}\t{at(paths[key])}/{name}\n")
                        done.add((sp, t))
        run("own", [py, nn / "own_strains.py", "--targets", near / "targets.tsv", "--test-strains",
                    own / "test_strains.tsv", "--summary", at(paths["refseq_summary"]), "--libraries", *libraries,
                    "--taxdump", at(paths["taxdump"]), "--skani", skani, "--threads", a.threads, "--out", own])
    extras = list(dict.fromkeys((target_of[r["set"]], r["extra_near_strains"]) for r in rows if r["extra_near_strains"]))
    if "extra" in steps:
        for target, sources in extras:
            if not (extra / target / "genomes.tsv").exists():
                run("extra", [py, nn / "near_neighbor_strains.py", "--near", near, "--target", target,
                              "--read-sources", inputs / sources, "--summary", at(paths["refseq_summary"]),
                              "--libraries", *libraries, "--skani", skani, "--threads", a.threads,
                              "--out", extra / target])
    if "others" in steps and not (others / "taxids.tsv").exists():
        run("others", [py, nn / "other_targets.py", "--near", near, "--own", own, "--out", others])
    has = lambda d: (d / "genomes.tsv").exists() and (d / "genomes.tsv").stat().st_size > 0
    folders = [near, *[extra / t for t, _ in extras], *([others] if has(others) else []), own]
    if "blast" in steps:
        for d in folders:
            genomes = [l.split("\t")[0] for l in open(d / "genomes.tsv") if l.strip()]
            if all((d / "paf" / f"{g}.paf").exists() for g in genomes):
                continue
            run("blast", [py, nn / "blast_homologs.py", "--folder", d, "--signatures", signatures / "signatures.fasta",
                          "--metadata", signatures / "signatures.csv", "--blastn", at(paths["blastn"]),
                          "--makeblastdb", at(paths["makeblastdb"]), "--processes", a.processes])
    if "index" in steps:
        for chem in ("r9", "r10"):
            if not (out / chem / "index.idx").exists():
                run("index", ["bash", HERE / "index/build_index.sh", signatures, chem, out / chem / "index.idx",
                              at(paths["onsight"]), a.threads])
    if "regions" in steps and not (out / "regions.lcr").exists():
        # one region file serves both chemistries: the two indexes list the same signatures in the same order
        sig = {c: [l.split("\t")[:2] for l in open(out / c / "index.idx.sigmeta.tsv")] for c in ("r9", "r10")}
        if sig["r9"] != sig["r10"]:
            raise SystemExit("the R9 and R10 indexes list different signatures")
        run("regions", [py, nn / "neighbor_regions.py", "build", "--panel", out / "panel.tsv",
                        "--targets", near / "targets.tsv", "--sigmeta", out / "r9" / "index.idx.sigmeta.tsv",
                        "--signatures", signatures / "signatures.fasta", "--near", *folders[:-1], "--own", own,
                        "--max-pad", 0, "--processes", a.processes, "--out", out / "regions.lcr"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

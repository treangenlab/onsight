#!/usr/bin/env python3
"""Is the sample's strain the species its Seqwin set targets? skani ANI of the strain genome against every assembly
in the set's lists (seqwin_lists.py), per NCBI species.

Species rule: ANI >= 95% with >= 50% of the strain genome aligned (skani 0.4.0), the rule of the all-genome background
and the B. spizizenii check.
- PASS: the target species' NCBI reference genome passes the rule and no other species' reference genome passes it
  with a higher ANI. When no reference genome of the target is in the lists: the closest assembly passing the rule
  carries the target label.
- CHECK otherwise; the set gets no Seqwin run until decided.

Writes in --out: skani.tsv (raw), strain_check.tsv (per species: assemblies, passing the rule, best ANI and aligned
share, the reference genome's ANI and aligned share), verdict.txt.

    conda activate onsight-panel
    python3 seqwin_strain_check.py --genome <fasta> --lists <dir> --target "<species>" --out <dir> --threads 16
"""
import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path

MIN_ANI, MIN_AF = 95.0, 50.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", type=Path, required=True)
    ap.add_argument("--lists", type=Path, required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    rows = [r for r in csv.DictReader(open(a.lists / "assemblies.tsv"), delimiter="\t") if r["path"]]
    by_path = {r["path"]: r for r in rows}
    (a.out / "references.txt").write_text("".join(f"{p}\n" for p in by_path))
    subprocess.run(["skani", "dist", "-q", str(a.genome), "--rl", str(a.out / "references.txt"), "-t", str(a.threads),
                    "-o", str(a.out / "skani.tsv")], check=True)
    hits = {}
    for h in csv.DictReader(open(a.out / "skani.tsv"), delimiter="\t"):
        hits[h["Ref_file"]] = (float(h["ANI"]), float(h["Align_fraction_query"]))
    passes = lambda x: x is not None and x[0] >= MIN_ANI and x[1] >= MIN_AF

    per = defaultdict(lambda: dict(assemblies=0, passing=0, best=None, reference=None))
    for p, r in by_path.items():
        s = per[r["species"] or "(no species)"]
        s["assemblies"] += 1
        h = hits.get(p)
        s["passing"] += passes(h)
        if h and (s["best"] is None or h > s["best"]):
            s["best"] = h
        if r["reference"] == "yes" and h and (s["reference"] is None or h > s["reference"]):
            s["reference"] = h
    fmt = lambda x: ("", "") if x is None else (f"{x[0]:.2f}", f"{x[1]:.1f}")
    table = sorted(per.items(), key=lambda kv: kv[1]["best"] or (0, 0), reverse=True)
    with open(a.out / "strain_check.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "assemblies", "passing_rule", "best_ani", "best_af_strain", "reference_ani",
                    "reference_af_strain"])
        for sp, s in table:
            w.writerow([sp, s["assemblies"], s["passing"], *fmt(s["best"]), *fmt(s["reference"])])

    own = per.get(a.target, {}).get("reference")
    others = [(s["reference"], sp) for sp, s in per.items() if sp != a.target and passes(s["reference"])]
    if any(r["reference"] == "yes" and r["species"] == a.target for r in rows):
        ok = passes(own) and not any(ref > own for ref, _ in others)
        why = f"reference genome of {a.target}: " + \
              ("no ANI reported (below skani's range)" if own is None else f"ANI {own[0]:.2f}, aligned {own[1]:.1f}%")
    else:
        best = max(((h, by_path[p]["species"]) for p, h in hits.items() if passes(h)), default=None)
        ok = best is not None and best[1] == a.target
        why = f"no reference genome of {a.target} in the lists; closest passing assembly: " + \
              ("none" if best is None else f"{best[1]} (ANI {best[0][0]:.2f})")
    if others:
        why += "; other species' reference genomes passing the rule: " + \
               ", ".join(f"{sp} ({ref[0]:.2f})" for ref, sp in sorted(others, reverse=True))
    verdict = "PASS" if ok else "CHECK"
    (a.out / "verdict.txt").write_text(f"{verdict}\t{why}\n")
    print(verdict, why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

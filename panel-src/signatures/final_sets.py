#!/usr/bin/env python3
"""The final signature sets of a panel: each Seqwin build without the signatures that background genomes carry.

--config (one row per final set): set (the name in the panel), build (the folder under lists/ and seqwin/),
report_as (the name a call is reported under), panel, mock_label, reference_genome (accession),
lineage_representative (lineage sets only).

For each row, --root/seqwin/<build>/signatures.fasta and signatures.csv are copied to --root/final/<set>/ without the
signatures that screen_signatures.py found in a background genome (--root/screen/<build>.carried.tsv): those would
call another species. Seqwin's seqwin.log and config.json are copied along (build_panel.py reads both). manifest.tsv
holds the config row, the build's target and negative counts, signatures built, removed and kept, held-out carriage,
and as checks only the mock genomes. --root/final/manifest.tsv indexes every set written so far. build_panel.py then
makes the panel from the final/<set>/ folders.

    conda activate onsight-panel
    python3 final_sets.py --root <dir> --config <sets.tsv> [--sets <set,...>]
"""
import argparse
import csv
import shutil
from pathlib import Path


def count_lines(path: Path) -> int:
    return sum(1 for line in open(path) if line.strip())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True, help="folder with lists/, seqwin/ and screen/")
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--sets", help="only these sets (comma separated)")
    a = ap.parse_args()
    only = set(a.sets.split(",")) if a.sets else None
    for s in csv.DictReader(open(a.config), delimiter="\t"):
        lists, build = a.root / "lists" / s["build"], a.root / "seqwin" / s["build"]
        screen = a.root / "screen" / f"{s['build']}.tsv"
        if (only and s["set"] not in only) or not screen.exists():
            continue                     # not asked for, or not screened yet
        e = next(csv.DictReader(open(screen), delimiter="\t"))
        drop = {name for r in csv.DictReader(open(a.root / "screen" / f"{s['build']}.carried.tsv"), delimiter="\t")
                if r["kind"] == "background" for name in r["signatures"].split(";")}
        out = a.root / "final" / s["set"]
        out.mkdir(parents=True, exist_ok=True)
        n_built = n_kept = bp = 0
        keep = False
        with open(out / "signatures.fasta", "w") as f:
            for line in open(build / "signatures.fasta"):
                if line.startswith(">"):
                    n_built += 1
                    keep = line[1:].split()[0] not in drop
                    n_kept += keep
                elif keep:
                    bp += len(line.strip())
                if keep:
                    f.write(line)
        with open(build / "signatures.csv") as src, open(out / "signatures.csv", "w") as dst:
            dst.write(next(src))
            dst.writelines(line for line in src if line.split(",", 1)[0] not in drop)
        for f in ("seqwin.log", "config.json"):
            shutil.copyfile(build / f, out / f)
        row = dict(s, targets=count_lines(lists / "targets.txt"), negatives=count_lines(lists / "negatives.txt"),
                   signatures_built=n_built, removed_by_background=n_built - n_kept, signatures=n_kept,
                   kb=round(bp / 1000), held_out_same_lineage_median=e["same_lineage_median"],
                   held_out_same_lineage_p10=e["same_lineage_p10"],
                   held_out_other_named_species_median=e["other_named_species_median"],
                   mock_signatures_check=e["mock_signatures"], mock_genomes_check=e["mock_genomes"],
                   mock_strain_check=e["mock_strain_carried"], lists=str(lists.resolve()), build=str(build.resolve()))
        with open(out / "manifest.tsv", "w") as f:
            w = csv.DictWriter(f, fieldnames=list(row), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerow(row)
        print(row["set"], row["signatures"], f"(removed {row['removed_by_background']})", flush=True)
    rows = [r for m in sorted(a.root.glob("final/*/manifest.tsv")) for r in csv.DictReader(open(m), delimiter="\t")]
    if rows:
        with open(a.root / "final" / "manifest.tsv", "w") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

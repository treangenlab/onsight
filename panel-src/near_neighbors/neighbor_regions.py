#!/usr/bin/env python3
"""Region file for the engine's near-neighbor check (src/neighbor_check.c, --neighbor-check).

Per signature of the index (sigmeta rid): its sequence and the regions the check compares it with, as 2-bit DNA.
Regions: every blastn hit of the signature in each own-strain and near-neighbor genome (each folder's paf/), so a
signature that aligns to a genome in several pieces has a region for each. A region is on the own side when its
genome's species taxid is the taxid of the signature's panel target. For each hit the file keeps every genome base a band on the signature can reach with a margin up to --max-pad:
the whole signature carried through the hit, plus the margin. Windows that overlap on one contig are stored once.

    conda activate squiggle
    python3 neighbor_regions.py build --panel <panel.tsv> --targets <near folder>/targets.tsv --sigmeta <sigmeta.tsv> \
      --signatures <signatures.fasta> --near <near folder> [<near folder> ...] --own <own folder> --max-pad 0 \
      --out <regions.lcr>
(--panel defaults to $S2S_PANEL.) The reads and cuts subcommands, test inputs for the C check against the Python
check, are in s2s/archive/python-scripts/near_neighbors/neighbor_regions_with_tests.py.

File layout (little-endian; onsight's src/neighbor_check.c reads it):
  header   magic "LCR1", u32 max_pad, u32 signatures, u32 genomes, u64 regions, u64 bases, u64 N runs, u64 name bytes
  sigs     per rid: u32 base, u32 length (0: no sequence), u32 first region, u32 regions
  regions  u32 base, i32 off, u16 length, u16 qs, u16 qe, u16 flags (bits 0-13 genome, 14 hit on '-', 15 own side)
  N runs   u32 start, u32 length: bases that are not A, C, G or T (stored as A in the 2-bit array)
LCR2 (magic "LCR2"), for files of 2^32 bases or more: each sig and region record ends with one more u32, the high 32
bits of base, and each N run is u32 start, u32 length, u32 high 32 bits of start. --format auto (default) writes LCR2
only when LCR1 cannot hold the bases, so smaller files stay byte-identical to earlier builds. onsight reads both.
  names    genome accessions, each ending in a zero byte
  bases    2 bits per base (A 0, C 1, G 2, T 3), base i in byte i / 4 at bit 2 * (i % 4)
A region's window covers genome bases [w0, w0 + length) of the hit's contig; off is the genome coordinate, minus w0,
of signature base 0 for a '+' hit (ts - qs) or of the base after the last for a '-' hit (te + qs).
"""
import argparse
import csv
import gzip
import multiprocessing as mp
import os
import struct
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

MAGIC = b"LCR1"
HEADER = struct.Struct("<4sIIIQQQQ")
SIG = np.dtype([("base", "<u4"), ("len", "<u4"), ("reg0", "<u4"), ("n_reg", "<u4")])
REG = np.dtype([("base", "<u4"), ("off", "<i4"), ("len", "<u2"), ("qs", "<u2"), ("qe", "<u2"), ("flags", "<u2")])
RUN = np.dtype([("start", "<u4"), ("len", "<u4")])
SIG2 = np.dtype(SIG.descr + [("base_hi", "<u4")])           # LCR2: the high 32 bits of base after the LCR1 fields
REG2 = np.dtype(REG.descr + [("base_hi", "<u4")])
RUN2 = np.dtype(RUN.descr + [("start_hi", "<u4")])
BASE4 = np.full(256, 4, np.uint8)
for _i, _b in enumerate("ACGT"):
    BASE4[ord(_b)] = BASE4[ord(_b.lower())] = _i


def read_fasta(path: str):
    name, chunks = None, []
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def homologs(folders: list) -> tuple:
    """Every hit per (signature, genome) from each folder's paf/, the genome files and taxids. A genome listed in two
    folders keeps the hits of the later one for a signature both searched."""
    hits, files, taxid = defaultdict(dict), {}, {}
    for folder in folders:
        files.update(l.rstrip("\n").split("\t") for l in open(folder / "genomes.tsv") if l.strip())
        taxid.update((r["genome"], int(r["species_taxid"])) for r in
                     csv.DictReader(open(folder / "taxids.tsv"), delimiter="\t"))
        for genome in (l.split("\t")[0] for l in open(folder / "genomes.tsv") if l.strip()):
            per = defaultdict(list)
            for line in open(folder / "paf" / f"{genome}.paf"):
                f = line.split("\t")
                per[f[0]].append((int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8])))
            for sig, hs in per.items():
                hits[sig][genome] = hs
    return dict(hits), files, taxid


def pack(seq: str) -> tuple:
    """2-bit bytes of seq (length padded to a multiple of 4) and its N runs [(start, length)]."""
    c = BASE4[np.frombuffer(seq.encode(), np.uint8)]
    runs = []
    n = np.flatnonzero(c == 4)
    if len(n):
        cut = np.flatnonzero(np.diff(n) > 1)
        starts, ends = np.r_[n[0], n[cut + 1]], np.r_[n[cut], n[-1]] + 1
        runs = list(zip(starts.tolist(), (ends - starts).tolist()))
        c[n] = 0
    c = np.concatenate([c, np.zeros(-len(c) % 4, np.uint8)])
    return (c[0::4] | (c[1::4] << 2) | (c[2::4] << 4) | (c[3::4] << 6)).astype(np.uint8).tobytes(), runs


def windows(path: str, wins: list) -> tuple:
    """One genome: each window (contig, a, b) clipped to its contig; windows that overlap on a contig are packed once.
    Returns (2-bit bytes, N runs, per window (w0, w1, base)), bases counted from the start of this genome's bytes;
    base -1 for an empty window."""
    seqs = {n: q.upper() for n, q in read_fasta(path)}
    by = defaultdict(list)
    for i, (contig, a, b) in enumerate(wins):
        by[contig].append(i)
    out, parts, runs, pos = [None] * len(wins), [], [], 0
    for contig, idx in by.items():
        s = seqs[contig]
        blocks = []
        for w0, w1, i in sorted((max(wins[i][1], 0), min(wins[i][2], len(s)), i) for i in idx):
            if w1 <= w0:
                out[i] = (0, 0, -1)
                continue
            if blocks and w0 <= blocks[-1][1]:
                blocks[-1][1] = max(blocks[-1][1], w1)
                blocks[-1][2].append((i, w0, w1))
            else:
                blocks.append([w0, w1, [(i, w0, w1)]])
        for b0, b1, members in blocks:
            p, r = pack(s[b0:b1])
            parts.append(p)
            runs += [(pos + x, n) for x, n in r]
            for i, w0, w1 in members:
                out[i] = (w0, w1, pos + w0 - b0)
            pos += 4 * len(p)
    return b"".join(parts), runs, out


def target_taxids(sigmeta: Path, panel: Path, targets: Path) -> dict:
    """rid -> taxid of its signature's panel target (panel.tsv: taxon -> target name; targets.tsv: name -> taxid, the
    check's same-species rule)."""
    name = {r["taxon"]: r["name"] for r in csv.DictReader(open(panel), delimiter="\t")}
    taxid = {r["name"]: int(r["taxid"]) for r in csv.DictReader(open(targets), delimiter="\t")}
    return {int(r["rid"]): taxid.get(name.get(r["source_taxon"])) for r in csv.DictReader(open(sigmeta), delimiter="\t")}


def build(a) -> int:
    t = time.time()
    meta = list(csv.DictReader(open(a.sigmeta), delimiter="\t"))
    rid_name = {int(r["rid"]): r["signature"] for r in meta}
    seqs = {n: s.upper() for n, s in read_fasta(str(a.signatures))}
    hits, files, taxid = homologs([*a.near, a.own])
    target = target_taxids(a.sigmeta, a.panel, a.targets)
    genomes = list(files)
    gi = {g: i for i, g in enumerate(genomes)}
    assert len(genomes) < 1 << 14, "genome index needs more than 14 bits"
    no_target = sorted({r["source_taxon"] for r in meta if target[int(r["rid"])] is None})
    print(f"signatures {len(rid_name)} (with a sequence {sum(n in seqs for n in rid_name.values())}); genomes "
          f"{len(genomes)}; signatures with a hit {len(hits)}; no panel taxid: {no_target or 'none'}", flush=True)

    # windows per genome: (signature, genome, hit) -> the genome bases any band can reach (margin up to max_pad)
    P = a.max_pad
    per_genome, where, risky = defaultdict(list), {}, 0
    for name, per in hits.items():
        if name not in seqs:
            continue
        L = len(seqs[name])
        assert L + 2 * P < 1 << 16, (name, L)
        for g, hs in per.items():
            for k, (qs, qe, ori, contig, ts, te) in enumerate(hs):
                anchor = ts - qs if ori == "+" else te + qs
                a0, b0 = (anchor - P, anchor + L + P) if ori == "+" else (anchor - L - P, anchor + P)
                if ori == "-" and te - (qe - qs) + 1 <= 0:     # near_neighbor_region's end could reach 0 or below
                    risky += 1
                where[(name, g, k)] = (len(per_genome[g]), anchor)
                per_genome[g].append((contig, a0, b0))
    print(f"hits {len(where)}; '-' hits whose cut could end at or before the contig start: {risky}", flush=True)

    # signatures first, then each genome's windows
    blob, runs, pos = [], [], 0
    sig_base = {}
    for name in sorted({n for n in rid_name.values() if n in seqs}):
        p, r = pack(seqs[name])
        sig_base[name] = pos
        blob.append(p)
        runs += [(pos + x, n) for x, n in r]
        pos += 4 * len(p)
    win = {}
    jobs = [(files[g], per_genome[g]) for g in genomes if per_genome[g]]
    names = [g for g in genomes if per_genome[g]]
    with mp.Pool(a.processes) as pool:
        for g, (p, r, out) in zip(names, pool.imap(job_windows, jobs)):
            blob.append(p)
            runs += [(pos + x, n) for x, n in r]
            win[g] = [(w0, w1, base + pos if base >= 0 else -1) for w0, w1, base in out]
            pos += 4 * len(p)
    wide = a.format == "2" or (a.format == "auto" and pos >= 1 << 32)
    assert wide or pos < 1 << 32, "bases need more than 32 bits: use --format 2 or auto"
    print(f"bases {pos:,}; N runs {len(runs)}; format {'LCR2' if wide else 'LCR1'}; {time.time() - t:.0f} s", flush=True)
    lo = lambda v: v & 0xFFFFFFFF

    # regions per rid, in the check's order (hits[name], near folder then own folder)
    n_sigs = max(rid_name) + 1
    sigs = np.zeros(n_sigs, SIG2 if wide else SIG)
    regs, n_own, n_empty = [], 0, 0
    for rid in range(n_sigs):
        name = rid_name.get(rid)
        if name not in seqs:
            continue
        b = sig_base[name]
        sigs[rid] = (lo(b), len(seqs[name]), len(regs), 0, b >> 32) if wide else (b, len(seqs[name]), len(regs), 0)
        for g, hs in hits.get(name, {}).items():
            for k, (qs, qe, ori, contig, ts, te) in enumerate(hs):
                j, anchor = where[(name, g, k)]
                w0, w1, base = win[g][j]
                if base < 0:
                    n_empty += 1
                    continue
                own = taxid.get(g) == target[rid] and target[rid] is not None
                n_own += own
                flags = gi[g] | (ori == "-") << 14 | own << 15
                regs.append((lo(base), anchor - w0, w1 - w0, qs, qe, flags, base >> 32) if wide
                            else (base, anchor - w0, w1 - w0, qs, qe, flags))
        sigs["n_reg"][rid] = len(regs) - sigs["reg0"][rid]
    regs = np.array(regs, REG2 if wide else REG)
    runs = np.array([(lo(x), n, x >> 32) for x, n in runs] if wide else runs, RUN2 if wide else RUN)
    names_b = b"".join(g.encode() + b"\0" for g in genomes)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "wb") as f:
        f.write(HEADER.pack(b"LCR2" if wide else MAGIC, P, n_sigs, len(genomes), len(regs), pos, len(runs), len(names_b)))
        f.write(sigs.tobytes())
        f.write(regs.tobytes())
        f.write(runs.tobytes())
        f.write(names_b)
        for p in blob:
            f.write(p)
    size = a.out.stat().st_size
    print(f"wrote {a.out}: {size / 1e6:.1f} MB (bases {pos / 4e6:.1f} MB, regions {len(regs):,} = "
          f"{regs.nbytes / 1e6:.1f} MB, own side {n_own:,}; empty windows dropped {n_empty}); "
          f"{time.time() - t:.0f} s", flush=True)
    return 0


def job_windows(job):
    return windows(*job)


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="what", required=True)
    b = sub.add_parser("build")
    b.add_argument("--panel", type=Path, default=os.environ.get("S2S_PANEL"), help="panel.tsv (default: $S2S_PANEL)")
    b.add_argument("--targets", type=Path, required=True, help="targets.tsv of near_neighbors.py sources")
    b.add_argument("--sigmeta", type=Path, required=True)
    b.add_argument("--signatures", type=Path, required=True)
    b.add_argument("--near", type=Path, nargs="+", required=True)
    b.add_argument("--own", type=Path, required=True)
    b.add_argument("--max-pad", type=int, default=0, help="largest margin the engine may use (default 0, adopted)")
    b.add_argument("--processes", type=int, default=6)
    b.add_argument("--format", choices=("auto", "1", "2"), default="auto",
                   help="region file format: LCR1 (32-bit base positions), LCR2 (64-bit), auto (LCR2 only when needed)")
    b.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.panel is None:
        ap.error("give --panel or set S2S_PANEL")
    return build(a)


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Summarize an ONS_PROF run (bin/onsight-prof) for hardware design.

Usage: prof_report.py PREFIX [--sample-rate HZ] [--out REPORT.md]

PREFIX is the ONS_PROF_OUT value of the run (e.g. runs/x/prof); the script reads
PREFIX.{reads,batches,index,summary,hist,keys}.tsv and, if present, time.txt and perf_stat.txt
next to them. Standard library only.
"""
import argparse
import csv
import math
import os
import sys

# Per-read stage timers (ns), by pipeline step. Step 2 runs per read on all threads (prepare), then in read order.
STEP0 = [("t_read", "read signal file (decompress, pA, 30-200 pA filter)")]
STEP1 = [("t_norm", "normalize (z-score, |z| < 3)"), ("t_prefix", "prefix sums"),
         ("t_tstat", "t-statistic (2 windows)"), ("t_peaks", "peak detection"),
         ("t_events", "segment -> event (sort, IQR mean)"), ("t_seeds", "4-bit quantization + strobe-event hash"),
         ("t_lookup", "index lookup"), ("t_accum", "seed hits per signature (grouped by signature, hit cap)")]
STEP2 = [("t_summary", "signature scores (--sig-summary)"), ("t_bands", "diagonal bands over slopes, band p-values"),
         ("t_pvalue", "read p-value, read call"), ("t_slopes", "slope range: chain of a confident call"),
         ("t_neighbor", "near-neighbor check (subsequence DTW)")]

PLATFORMS = [("MinION (512 ch)", 512), ("PromethION flow cell (2,675 ch)", 2675)]


def read_kv(path):
    d = {}
    if not os.path.exists(path):
        return d
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("key\t"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) == 2:
                d[p[0]] = p[1]
    return d


def pct(sorted_vals, q):
    if not sorted_vals:
        return 0
    i = min(len(sorted_vals) - 1, max(0, int(math.ceil(q * len(sorted_vals))) - 1))
    return sorted_vals[i]


def hist_quantile(h, q):
    tot = sum(h.values())
    acc = 0
    for b in sorted(h):
        acc += h[b]
        if acc >= q * tot:
            return b
    return max(h) if h else 0


def fmt(x, nd=1):
    if isinstance(x, int) or (isinstance(x, float) and x.is_integer() and abs(x) >= 1000):
        return f"{int(x):,}"
    if abs(x) >= 100:
        return f"{x:,.0f}"
    return f"{x:,.{nd}f}"


def fbytes(b):
    for u in ["B", "KB", "MB", "GB"]:
        if abs(b) < 1024 or u == "GB":
            return f"{b:,.1f} {u}" if u != "B" else f"{int(b)} B"
        b /= 1024.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--sample-rate", type=float, default=4000.0,
                    help="signal sample rate in Hz (R9.4.1: 4000, R10.4.1: 5000)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    P = a.prefix
    out = open(a.out, "w") if a.out else sys.stdout
    w = lambda s="": print(s, file=out)

    summ = read_kv(P + ".summary.tsv")
    idx = read_kv(P + ".index.tsv")

    # ---- per-read table: totals and distributions ----
    tot = {}
    dist_cols = ["l_sig", "n_events", "n_seeds", "n_found", "post_found", "post_visited", "n_sigs",
                 "n_hits_kept", "km_cap", "hits_bytes", "t_step1", "t_prepare", "t_commit"]
    dists = {c: [] for c in dist_cols}
    n = n_call = n_checked = n_kept = 0
    with open(P + ".reads.tsv") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        ci = {h: i for i, h in enumerate(hdr)}
        num_cols = [h for h in hdr if h not in ("name", "call", "checked")]
        for h in num_cols:
            tot[h] = 0
        for row in rd:
            n += 1
            for h in num_cols:
                tot[h] += int(row[ci[h]])
            for c in dist_cols:
                dists[c].append(int(row[ci[c]]))
            n_call += int(row[ci["call"]])
            ck = int(row[ci["checked"]])
            n_checked += ck >= 0
            n_kept += ck == 1
    for c in dist_cols:
        dists[c].sort()

    batches = []
    with open(P + ".batches.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            batches.append(r)
    bsum = lambda k: sum(int(b[k]) for b in batches)
    b_step0, b_step1, b_busy, b_busy_max = bsum("step0_ns"), bsum("step1_ns"), bsum("busy_sum_ns"), bsum("busy_max_ns")
    b_par, b_ser, b_end, b_step2 = bsum("step2_par_ns"), bsum("step2_ser_ns"), bsum("end_batch_ns"), bsum("step2_ns")
    max_rss = max(int(b["rss_bytes"]) for b in batches) if batches else 0

    hist = {}
    if os.path.exists(P + ".hist.tsv"):
        with open(P + ".hist.tsv") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                hist.setdefault(r["hist"], {})[int(r["bin"])] = int(r["count"])

    S = tot["l_sig"]; E = tot["n_events"]; Q = tot["n_seeds"]
    wall_mono = float(summ.get("wall_s", 0))
    # Clock correction. Under WSL2 the kernel clock (CLOCK_MONOTONIC) can run fast against real time;
    # /usr/bin/time's elapsed time is real time. Scale every measured duration by real / monotonic.
    cf, real = 1.0, None
    tp = os.path.join(os.path.dirname(P), "time.txt")
    if os.path.exists(tp):
        for line in open(tp):
            if "Elapsed (wall clock)" in line:
                v = [float(x) for x in line.rsplit(" ", 1)[1].strip().split(":")]
                real = v[-1] + (v[-2] * 60 if len(v) > 1 else 0) + (v[-3] * 3600 if len(v) > 2 else 0)
    if real and wall_mono > 0:
        cf = real / wall_mono
        if abs(cf - 1.0) < 0.02:
            cf = 1.0  # within measurement noise: leave uncorrected
    for k in list(tot):
        if k.startswith("t_") and k != "t_step1_cpu":  # CPU-time accounting runs at real speed
            tot[k] = int(tot[k] * cf)
    for c in ("t_step1", "t_prepare", "t_commit"):
        dists[c] = [int(x * cf) for x in dists[c]]
    b_step0 *= cf; b_step1 *= cf; b_busy *= cf; b_busy_max *= cf
    b_par *= cf; b_ser *= cf; b_end *= cf; b_step2 *= cf
    wall = wall_mono * cf
    workers = int(summ.get("workers", 1))

    w(f"# ONS_PROF report: `{os.path.basename(P)}`")
    w()
    w("## Run")
    w()
    w("| | |\n|---|---|")
    w(f"| Reads | {fmt(n)} |")
    w(f"| Signal samples (after pA filter) | {fmt(S)} |")
    w(f"| Events | {fmt(E)} ({S / max(1, E):.2f} samples/event) |")
    w(f"| Seeds = strobe-events (= index lookups) | {fmt(Q)} ({Q / max(1, E):.2f} seeds/event) |")
    w(f"| Read calls (p ≤ top cutoff) | {fmt(n_call)}; near-neighbor checked {fmt(n_checked)}, kept {fmt(n_kept)} |")
    w(f"| Wall time | {wall:,.1f} s |")
    w(f"| CPU time (user + sys) | {float(summ.get('user_cpu_s', 0)) + float(summ.get('sys_cpu_s', 0)):,.1f} s |")
    if cf != 1.0:
        w(f"| Clock correction applied | ×{cf:.4f} to stage timers (CLOCK_MONOTONIC {wall_mono:,.1f} s vs real elapsed {real:,.1f} s); CPU times are not affected |")
    w(f"| Peak RSS | {fbytes(int(summ.get('peak_rss_bytes', 0)))} (profiler's own: {fbytes(int(summ.get('profiler_self_bytes', 0)))}) |")
    w(f"| Worker threads / I/O threads / pipeline threads | {workers} / {summ.get('io_threads')} / {summ.get('pipeline_threads')} |")
    w(f"| Minibatch `-K` (samples) | {fmt(int(summ.get('opt_batch_samples', 0)))} ({fmt(len(batches))} minibatches, {n / max(1, len(batches)):.1f} reads each) |")
    w()

    # ---- pipeline ----
    w("## Pipeline: where the wall time goes")
    w()
    w("Each minibatch goes through three steps: read signal files → events, seeds and seed hits (all workers) →")
    w("output (per read on all workers, then in read order). One minibatch at a time is in each step; the steps of")
    w("different minibatches overlap, so the rows add up to more than 100%")
    w(f"(at most {int(summ.get('pipeline_threads', 1)) * 100}% with {summ.get('pipeline_threads')} pipeline threads). A step near 100% is the bottleneck.")
    w()
    w("| Step | Busy (s) | Busy, % of wall time | Notes |\n|---|---|---|---|")
    w(f"| 0 read signal files | {b_step0 / 1e9:,.1f} | {b_step0 / 1e9 / max(1e-9, wall) * 100:.0f}% | serial |")
    w(f"| 1 events, seeds, seed hits ({workers} workers) | {b_step1 / 1e9:,.1f} | {b_step1 / 1e9 / max(1e-9, wall) * 100:.0f}% | worker busy {b_busy / 1e9:,.1f} s; parallel efficiency {b_busy / max(1, b_step1 * workers) * 100:.0f}% (slowest worker {b_busy_max / max(1, b_step1) * 100:.0f}% of step) |")
    w(f"| 2 output and decisions | {b_step2 / 1e9:,.1f} | {b_step2 / 1e9 / max(1e-9, wall) * 100:.0f}% | per read on {workers} workers {b_par / 1e9:,.1f} s; in read order {b_ser / 1e9:,.1f} s; slope grid update {b_end / 1e9:,.2f} s |")
    w()
    w("Wait and work per pipeline thread (s):")
    w()
    w("| Thread | Step | Work | Wait |\n|---|---|---|---|")
    for k in sorted(summ):
        if k.startswith("pipeline_thread") and k.endswith("_work_s"):
            base = k[:-len("_work_s")]
            th, st = base.split("_")[1], base.split("_")[2]
            w(f"| {th} | {st} | {float(summ[k]) * cf:,.1f} | {float(summ.get(base + '_wait_s', 0)) * cf:,.1f} |")
    w()

    # ---- stage breakdown ----
    w("## Per-stage time")
    w()
    w("Wall time measured inside each stage on the thread doing it. Stages that run on several workers add up")
    w("across them; per-unit costs are what a single core spends.")
    w()
    w("| Step | Stage | Total (s) | Share of step | ns / read | ns / sample | ns / event | ns / seed |")
    w("|---|---|---|---|---|---|---|---|")
    row = lambda step, name, t, denom: w(
        f"| {step} | {name} | {t / 1e9:,.2f} | {t / max(1, denom) * 100:.1f}% | {t / max(1, n):,.0f} | "
        f"{t / max(1, S):.2f} | {t / max(1, E):.1f} | {t / max(1, Q):.1f} |")
    row("0", STEP0[0][1], tot["t_read"], tot["t_read"])
    acc = 0
    for col, name in STEP1:
        acc += tot[col]
        row("1", name, tot[col], tot["t_step1"])
    row("1", "other (read setup)", tot["t_step1"] - acc, tot["t_step1"])
    row("1", "**total step 1**", tot["t_step1"], tot["t_step1"])
    d2 = tot["t_summary"] + tot["t_prepare"] + tot["t_commit"]
    acc = 0
    for col, name in STEP2:
        acc += tot[col]
        row("2", name, tot[col], d2)
    row("2", "other per-read scoring", tot["t_summary"] + tot["t_prepare"] - acc, d2)
    row("2", "in read order: rows, slope ring, read record, read/target calls", tot["t_commit"], d2)
    row("2", "**total step 2 (per read)**", d2, d2)
    w()
    w(f"Step 1 CPU time measured per thread: {tot['t_step1_cpu'] / 1e9:,.1f} s vs {tot['t_step1'] / 1e9:,.1f} s wall in-thread "
      f"({tot['t_step1_cpu'] / max(1, tot['t_step1']) * 100:.0f}%; the gap is time the thread was not running).")
    w()

    # ---- workload ----
    w("## Workload per read")
    w()
    w("| Quantity | Mean | p50 | p99 | Max |\n|---|---|---|---|---|")
    labels = {"l_sig": "signal samples", "n_events": "events", "n_seeds": "seeds / lookups",
              "n_found": "lookups that hit", "post_found": "index positions returned",
              "post_visited": "positions added as seed hits (after repeat filters)",
              "n_sigs": "signatures with hits", "n_hits_kept": "seed hits kept (all signatures)",
              "km_cap": "memory pool bytes (signal processing)", "hits_bytes": "seed-hit array bytes",
              "t_step1": "step 1 time (ns)", "t_prepare": "stream scoring time (ns)",
              "t_commit": "in-order output time (ns)"}
    for c in dist_cols:
        v = dists[c]
        w(f"| {labels[c]} | {fmt(sum(v) / max(1, n))} | {fmt(pct(v, .5))} | {fmt(pct(v, .99))} | {fmt(v[-1] if v else 0)} |")
    w()
    w("| Ratio | Value |\n|---|---|")
    w(f"| Lookups that hit | {tot['n_found'] / max(1, Q) * 100:.1f}% |")
    w(f"| Hits with a single position | {tot['n_single'] / max(1, tot['n_found']) * 100:.1f}% |")
    w(f"| Positions per hit | {tot['post_found'] / max(1, tot['n_found']):.1f} |")
    w(f"| Extra hash probes per lookup | {tot['n_probe'] / max(1, Q):.2f} |")
    w(f"| Positions dropped by the repeat cutoff `mid_occ` ({summ.get('opt_mid_occ')}) | {tot['post_flt_occ'] / max(1, tot['post_found']) * 100:.1f}% |")
    w(f"| Positions dropped by `--sig-max-hash-occ` ({summ.get('opt_sig_max_hash_occ')}) | {tot['post_skip_adm'] / max(1, tot['post_found']) * 100:.1f}% |")
    w(f"| New (signature, position) per seed hit | {tot['n_pos_new'] / max(1, tot['post_visited']) * 100:.1f}% |")
    w(f"| Hit cap full: lightest-hit scans / replacements | {fmt(tot['n_cap_scan'])} / {fmt(tot['n_cap_evict'])} |")
    w(f"| Targets per read (stream scoring) | {tot['n_taxa'] / max(1, n):.1f} |")
    w(f"| Reads below `--min-events` | {tot['below_min_events'] / max(1, n) * 100:.1f}% |")
    w()

    # ---- real-time requirement ----
    sps = S / max(1, n)
    w(f"## Real-time requirement (sample rate {a.sample_rate:,.0f} Hz)")
    w()
    w("Rates the chip must sustain when every channel streams, and the CPU cores this code would need for it.")
    w()
    w("| Quantity | per sample | " + " | ".join(p for p, _ in PLATFORMS) + " |")
    w("|---|---|" + "---|" * len(PLATFORMS))
    rows = [("samples / s", 1.0), ("events / s", E / max(1, S)), ("index lookups / s", Q / max(1, S)),
            ("index positions read / s", tot["post_found"] / max(1, S)),
            ("seed-hit table updates / s", tot["post_visited"] / max(1, S)), ("reads finished / s", n / max(1, S))]
    for name, per in rows:
        w(f"| {name} | {per:.3g} | " + " | ".join(fmt(per * ch * a.sample_rate) for _, ch in PLATFORMS) + " |")
    core_rows = [("step 1 (events, seeds, seed hits)", tot["t_step1"] / max(1, S)),
                 ("step 2 (scoring and decisions)", d2 / max(1, S)),
                 ("step 0 (file reading, not on chip)", tot["t_read"] / max(1, S))]
    for name, nsps in core_rows:
        w(f"| CPU cores for {name} | {nsps:.2f} ns | " + " | ".join(f"{nsps * ch * a.sample_rate / 1e9:,.1f}" for _, ch in PLATFORMS) + " |")
    w()
    w(f"Mean read is {sps / a.sample_rate:,.1f} s of signal ({fmt(sps)} samples).")
    w()

    # ---- memory ----
    w("## Memory")
    w()
    ib = sum(int(idx.get(k, 0)) for k in ["bytes_bucket_array", "bytes_hash_keys", "bytes_hash_vals",
                                         "bytes_hash_flags", "bytes_postings_array"])
    w("| Structure | Size |\n|---|---|")
    w(f"| Index, as laid out now (buckets + hash keys/values/flags + positions) | {fbytes(ib)} |")
    w(f"| … keys / positions | {fmt(int(idx.get('keys', 0)))} keys ({idx.get('bits_key_used')} key bits after the bucket), {fmt(int(idx.get('total_postings', 0)))} positions |")
    w(f"| … minimum payload (positions × 64 bit) | {fbytes(int(idx.get('total_postings', 0)) * 8)} |")
    w(f"| … hash slots / load factor | {fmt(int(idx.get('hash_slots', 0)))} / {float(idx.get('load_factor', 0)):.2f} |")
    if dists["km_cap"]:
        w(f"| Memory pool of a thread (kept between reads), p99 / max | {fbytes(pct(dists['km_cap'], .99))} / {fbytes(dists['km_cap'][-1])} |")
        w(f"| Seed-hit arrays per read, p99 / max | {fbytes(pct(dists['hits_bytes'], .99))} / {fbytes(dists['hits_bytes'][-1])} |")
    w(f"| Peak process RSS (minibatch sampling) | {fbytes(max_rss)} |")
    w()

    # ---- histograms ----
    w("## Distributions")
    w()
    if "segment_len" in hist:
        h = hist["segment_len"]
        w(f"- Segment length (samples per event): p50 {hist_quantile(h, .5)}, p99 {hist_quantile(h, .99)}, "
          f"p99.9 {hist_quantile(h, .999)}, longest {max(h)}{'+' if max(h) >= 1023 else ''}. Segments outside the accepted length range are dropped.")
    if "positions_per_found_seed" in hist:
        h = hist["positions_per_found_seed"]
        w(f"- Positions per hit: p50 {hist_quantile(h, .5)}, p90 {hist_quantile(h, .9)}, p99 {hist_quantile(h, .99)}, max {max(h)}.")
    if "probe_steps_per_lookup" in hist:
        h = hist["probe_steps_per_lookup"]
        w(f"- Extra probe steps per lookup: p50 {hist_quantile(h, .5)}, p99 {hist_quantile(h, .99)}, max {max(h)}.")
    if "hits_per_read_sig" in hist:
        h = hist["hits_per_read_sig"]
        w(f"- Seed hits kept per (read, signature): p50 {hist_quantile(h, .5)}, p99 {hist_quantile(h, .99)}, max {max(h)}.")
    w()
    kp = P + ".keys.tsv"
    if os.path.exists(kp):
        meta, rows = {}, []
        with open(kp) as f:
            for line in f:
                if line.startswith("#"):
                    k, v = line[2:].strip().split("\t"); meta[k] = int(v)
                elif not line.startswith("top_n"):
                    rows.append(line.split())
        w("## Index key locality")
        w()
        w(f"{fmt(meta.get('distinct_keys_touched', 0))} distinct keys were looked up and hit, out of {fmt(int(idx.get('keys', 0)))}. "
          "Share of hits absorbed by the hottest N keys:")
        w()
        w("| Hottest keys | Share of hits |\n|---|---|")
        for r in rows:
            if int(r[0]) >= 1024:
                w(f"| {fmt(int(r[0]))} | {float(r[2]) * 100:.1f}% |")
        w()
    ps = os.path.join(os.path.dirname(P), "perf_stat.txt")
    if os.path.exists(ps):
        w("## Hardware counters (perf stat, whole process)")
        w()
        w("```")
        with open(ps) as f:
            for line in f:
                if any(k in line for k in ["cycles", "instructions", "cache", "branch", "TLB", "elapsed"]):
                    w(line.rstrip())
        w("```")


if __name__ == "__main__":
    main()

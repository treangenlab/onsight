# Profiling onsight for hardware co-design

This is the same profiling we ran on Squig2Sig (`hw-prof` branch, `scripts/hw_prof/`), ported to onsight. It
measures per-stage time and per-read work counts, the busy time of each pipeline step, the index layout and its
lookup locality, and optionally hardware counters with `perf`.

Profiling never changes the results. All profiling code is behind `#ifdef ONS_PROF`, the normal `bin/onsight`
does not contain it, and the profiling binary writes the same output files as the normal one (checked in step 3).

## What was added

| File | Change |
|---|---|
| `src/onsprof.h`, `src/onsprof.c` | new: the per-read record, histograms and output files |
| `src/pipeline.c` | timers for the 3 pipeline steps, per minibatch and per read; one profiling record per read, written in read order |
| `src/event_detection.c` | timers for normalize / prefix sums / t-statistic / peaks / events; segment counts and lengths |
| `src/seed_lookup.c`, `src/index.c` | lookup counts, positions per found seed, hash probe steps, lookups per bucket and per key; index layout dump |
| `src/seed_hits.c` | lookup and seed-hit timers; repeat-filter, `--sig-max-hash-occ` and hit-cap counts |
| `src/streaming.c` | timers for the diagonal bands, read p-value, slope range and near-neighbor check |
| `src/third_party/kthread.c` | wait and work time of each pipeline thread per step |
| `src/signal_input.h`, `src/seed_hits.h` | one profiling field each (time to read the signal; the read's record) |
| `Makefile`, `src/Makefile` | `make prof` → `bin/onsight-prof`, `make attr` → `bin/onsight-attr` (objects in `src/obj-prof/`, `src/obj-attr/`) |
| `scripts/hw_prof/run_profile.sh` | runs onsight on a file list, with or without profiling and `perf` |
| `scripts/hw_prof/prof_report.py` | turns the profiling files into a Markdown report |
| `scripts/hw_prof/perf_split.py` | splits a multi-event `perf report` into one table per event (unchanged from Squig2Sig) |

## 1. Build

```bash
cd onsight
git checkout hw-prof
git submodule update --init extern/slow5lib
conda activate onsight      # or the s2s environment from the Squig2Sig handoff; it has the same libraries
make -j 8                   # bin/onsight        (normal, unchanged)
make prof                   # bin/onsight-prof   (profiling)
make attr                   # bin/onsight-attr   (only for per-function perf attribution)
```

`make prof` and `make attr` build the normal binary first (for the external libraries), then the profiling one.

You also need `python3` for the report (standard library only) and optionally Linux `perf` for hardware counters
(`kernel.perf_event_paranoid` must be ≤ 2).

## 2. Index and file list

The Zymo run from the Squig2Sig handoff kit (`2026-10-18_hardware_zymo-ff1`: same FAST5 files, same signatures).
onsight has its own index format, so build the index with onsight (2 s):

```bash
KIT=<path to 2026-10-18_hardware_zymo-ff1>
FAST5=<folder with the 201 .fast5 files>
WORK=<an output folder>; mkdir -p $WORK

bin/onsight --pore r9 -t 4 -d $WORK/zymo.idx \
  --metadata $KIT/reference/signatures.csv $KIT/reference/signatures.fasta

sed "s#^#$FAST5/#" $KIT/fast5_in_run_order.txt > $WORK/fast5_in_run_order.txt
head -2 $WORK/fast5_in_run_order.txt > $WORK/first2.txt
```

`--pore r9` sets everything else; the defaults are the settings of the paper (strobe-events 2 × 6 events, gap 8–20).

Keep the FAST5 files on a local disk. Reading them from a network or Windows-mounted drive (`/mnt/d` on WSL2) slows
step 0 about 3× and changes the timing picture.

## 3. Check: 2 files with profiling (about 1 min)

```bash
scripts/hw_prof/run_profile.sh --prof $WORK/check_prof $WORK/zymo.idx $WORK/first2.txt
md5sum $WORK/check_prof/{read_records.jsonl,read_calls.tsv,sig_summary.tsv}
```

Expected (the same with `bin/onsight` and `bin/onsight-prof`):

| File | md5 |
|---|---|
| `read_records.jsonl` | `93822f9506fb0a9c8645dc2f8fc4a7fd` |
| `read_calls.tsv` | `6736fceffca1f95296a063fdb24fc862` |
| `sig_summary.tsv` | `1922fc8654852fe1662de9abda3149a6` (also the same as Squig2Sig's on these 2 files) |

`target_calls.tsv` has a wall-clock column (`first_call_time_s`), so compare it without that column:
`cut -f1-7,9- target_calls.tsv`. All 8 species are called within the first 2 files.

### `run_profile.sh`

```
scripts/hw_prof/run_profile.sh [--prof] [--perf] [--keys] [--attr] [--pore r9|r10] [--neighbor FILE] [--panel FILE]
                               OUTDIR INDEX FILE_LIST [THREADS] [BATCH]
```

| Option | What it does |
|---|---|
| `--prof` | use `bin/onsight-prof`; writes `OUTDIR/prof.*.tsv` |
| `--perf` | `perf stat` (`perf_stat.txt`) and `perf record` with call graphs at 499 Hz (`perf.data`) |
| `--keys` | count lookups per index key (key locality; +35 MB memory) |
| `--attr` | use `bin/onsight-attr` (hot helpers not inlined, for per-function attribution) |
| `--pore` | `r9` (default) or `r10` |
| `--neighbor FILE` | turn on the near-neighbor check with a panel's `regions.lcr` (off without it) |
| `--panel FILE` | a panel's `panel.tsv` (target names); without it each taxon of the index is a target |
| `THREADS`, `BATCH` | `-t` and `-K`; default `4` and `4M` |

Every run writes `time.txt`, `stderr.log`, `run_info.txt` (binary, git commit, settings) and onsight's outputs
(`read_calls.tsv`, `target_calls.tsv`, `sig_summary.tsv`, and the read records from fd 3 in `read_records.jsonl`).

**Threads.** With `-t 4` onsight uses 1 I/O thread and 3 workers (`-t` minus `--io-thread`), and 2 pipeline threads.
Steps 1 and 2 each run on all 3 workers. **`-K` changes the results** (the slope range is refitted between
minibatches), so keep `4M` except for an explicit sweep.

## 4. Whole run (802,603 reads)

```bash
scripts/hw_prof/run_profile.sh --prof --perf --keys $WORK/whole_run $WORK/zymo.idx $WORK/fast5_in_run_order.txt
python3 scripts/hw_prof/prof_report.py $WORK/whole_run/prof --out $WORK/whole_run/report.md
```

The 2-file check took 65 s for 8,000 reads, so the whole run should take about 1.5–2 h on 4 threads. Run it under
`nohup` or `tmux`. Leave out `--perf` if `perf` is not available. Use `--sample-rate 5000` in the report for R10.4.1.

### Report sections

| Section | What to look at |
|---|---|
| Run | reads, samples, events, seeds, read calls, wall/CPU time, peak memory |
| Pipeline: where the wall time goes | busy % of each step (read files / events, seeds, seed hits / output and decisions); a step near 100% is the bottleneck |
| Per-stage time | ns per signal sample for each stage of step 1 and step 2 |
| Workload per read | events, lookups, positions, seed hits, signatures touched per read (mean, p50, p99) |
| Real-time requirement | the same rates scaled to a MinION (512 channels × 4 kHz) and a PromethION flow cell |
| Memory | index layout and size, per-read memory pool and seed-hit arrays, peak RSS |
| Distributions, Index key locality | segment lengths, positions per hit, probe steps, share of lookups absorbed by the hottest keys |
| Hardware counters | IPC, cache and TLB misses (with `--perf`) |

With `--perf`: `perf report -i $WORK/whole_run/perf.data --no-children --sort symbol` shows time per function. For
per-stage cache misses, run with `--attr --perf` and split the events with
`perf report --stdio ... > full.txt; python3 scripts/hw_prof/perf_split.py full.txt OUTDIR`.

## 5. The stages, and how they map to Squig2Sig's

onsight's pipeline differs from Squig2Sig's, so the stage names differ. Steps 0 and 1 are the same work. Step 2 was
the serial vote/checkpoint gate in Squig2Sig. In onsight it is the decision path: diagonal bands over the slope grid,
read p-value, near-neighbor check and target call, and it runs per read on all workers.

| Step | Timer (`prof.reads.tsv`) | Stage | Squig2Sig timer |
|---|---|---|---|
| 0 | `t_read` | read the signal (FAST5 decompress, pA, 30–200 pA filter) | `t_rd_attr + t_rd_data + t_rd_conv` |
| 1 | `t_norm` | normalize (z-score, keep \|z\| < 3) | `t_norm` |
| 1 | `t_prefix` | prefix sums | `t_prefix` |
| 1 | `t_tstat` | t-statistic, 2 windows | `t_tstat` |
| 1 | `t_peaks` | peak detection | `t_peaks` |
| 1 | `t_events` | segment → event (sort, IQR mean) | `t_events` |
| 1 | `t_seeds` | 4-bit quantization + strobe-event hashing | `t_sketch` |
| 1 | `t_lookup` | index lookup | `t_lookup` |
| 1 | `t_accum` | seed hits per signature (grouped by signature, hit cap) | `t_accum` |
| 1 | `t_step1` | whole step 1 for the read (contains the above), `t_step1_cpu` its CPU time | `t_map`, `t_map_cpu` |
| 2 | `t_summary` | signature ranking, only with `--sig-summary` or `--sig-profile` | (in `t_rank`) |
| 2 | `t_bands` | diagonal bands over the slope grid, non-overlapping hit counts, band p-values | `t_diag` |
| 2 | `t_pvalue` | read p-value and read call | — |
| 2 | `t_slopes` | sliding-window slope range (chain of a confident call) | — |
| 2 | `t_neighbor` | near-neighbor check (subsequence DTW), only with `--neighbor` | — |
| 2 | `t_prepare` | whole stream scoring (contains `t_bands` … `t_neighbor`) | `t_stream` |
| 2 | `t_commit` | in read order: rows, slope ring, read record, read and target calls | `t_ckpt` + output |

Other columns: `l_sig`, `n_events`, `n_seeds` (= lookups), `n_found`, `n_single`, `n_probe`, `post_found` (index
positions returned), `n_flt_occ`/`post_flt_occ` (repeat cutoff `mid_occ`), `n_skip_adm`/`post_skip_adm`
(`--sig-max-hash-occ`), `post_visited` (positions added as seed hits), `n_pos_new`, `n_cap_scan`/`n_cap_evict` (the
1,024-hit cap per signature), `n_sigs` (signatures with hits), `n_hits_kept`, `hits_bytes`, `km_cap` (memory pool),
`n_taxa`, `call` (1: read call), `checked` (near-neighbor check: −2 not checked, −1 no band, 0 dropped, 1 kept).

Other files:

| File | Contents |
|---|---|
| `prof.batches.tsv` | per minibatch: reads, samples, step 0 / 1 time, worker busy sum and max, step 2 split into per-read (`step2_par_ns`), in read order (`step2_ser_ns`) and the slope-grid update (`end_batch_ns`), RSS |
| `prof.summary.tsv` | totals, CPU time, peak memory, work/wait per pipeline thread and step |
| `prof.index.tsv` | index layout: buckets, keys, hash slots, positions, bytes |
| `prof.hist.tsv` | segment length, positions per found seed, probe steps per lookup, seed hits per (read, signature) |
| `prof.buckets.tsv`, `prof.keys.tsv` | lookups per bucket; lookups absorbed by the hottest N keys (`--keys`) |

```python
import pandas as pd
r = pd.read_csv("whole_run/prof.reads.tsv", sep="\t")
r.filter(regex="^t_").sum() / r.l_sig.sum()      # ns per signal sample, by timer (some timers contain others, see above)
```

## 6. First numbers (2-file check, 8,000 reads, `-t 4 -K 4M`)

Ryzen 7 9800X3D under WSL2, data on ext4. These are from the check run only; the whole run will be more reliable.

| | |
|---|---|
| Species called | 8 / 8 (all within the first 2 files) |
| Wall time, CPU, peak memory | 65 s, 4.4 cores, 0.54 GB |
| Pipeline step busy | read 19%, step 1 61%, **step 2 95% (the bottleneck)** |
| Step 1 (ns per sample) | 103: event detection 45, quantization + hashing 5, lookup 11, seed hits 38 |
| Step 2 (ns per sample) | 158: **diagonal bands 147**, signature scores 9, the rest < 2 |
| Per read | 9,022 signatures with hits (of 14,578); 32,000 seed hits kept |
| Lookups that hit / positions per hit | 41% / 11.3 |
| `perf` (user cycles) | `ri_band_nonoverlap` 57%, `ri_collect_seed_hits` 9%, `detect_events` 7%, sorts 8%, `ri_idx_get` 4%; IPC 1.11 |
| Index key locality | 4.1 M of 5.6 M keys hit; the hottest 64 K keys take 21% of lookups, 1 M keys 74% |

The main difference from Squig2Sig is that the decision step now dominates. Most of it is the diagonal-band count
over 151 slopes for every signature a read touches, not the read p-value or the target call. The near-neighbor
check is not in these numbers: the Zymo kit has no `regions.lcr`. To profile it, run with a panel that has one
(`--neighbor panel/regions.lcr`).

## Running on an NVIDIA Jetson Nano (ARM64)

The profiling code is plain C/C++ (`clock_gettime`, `/proc/self/statm`, GCC atomics) with nothing x86-specific, and
the build already handles aarch64 (`-march=native`, `-fsigned-char`, POD5's linux-arm64 libraries). `make`, `make
prof` and `make attr` all build for aarch64, and the scripts run on JetPack 4's Python 3.6.

**Setup (once).**

```bash
# Miniforge for aarch64 (JetPack's Ubuntu 18.04 is too old for the build's GCC 11 requirement)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh
bash Miniforge3-Linux-aarch64.sh
sudo apt install time                    # /usr/bin/time, used by run_profile.sh

cd onsight && git checkout hw-prof
git submodule update --init extern/slow5lib
sed 's/sysroot_linux-64=/sysroot_linux-aarch64=/' environment.yml > environment-aarch64.yml
conda env create -f environment-aarch64.yml && conda activate onsight
make -j 2 && make prof && make attr      # -j 2: the Nano's 4 GB is shared with the GPU
```

**Before each timed run** set the maximum power mode and fixed clocks. Otherwise DVFS and thermal throttling change
the stage times from run to run:

```bash
sudo nvpmodel -m 0        # MAXN, 4 cores at full speed
sudo jetson_clocks        # fix the clocks
```

`run_profile.sh` writes the architecture, CPU, power mode and L4T release to `run_info.txt`. Watch the temperature on
long runs (`tegrastats`).

**Data.** Put the FAST5 (or BLOW5) files on a USB 3 SSD or NVMe drive, not the SD card; an SD card makes step 0 the
bottleneck. The whole Zymo run (157 GB) does not fit on most SD cards anyway. Start with the 2-file check (1.4 GB) to
calibrate the run time. FAST5 needs the VBZ plugin, which the conda environment provides for aarch64 as well.

**Threads.** The Nano has 4 cores, so the default `-t 4` (1 I/O thread + 3 workers) matches the x86 runs. Peak memory
in the check was 0.54 GB, which fits.

**perf.** JetPack does not ship `perf`. Build it from the L4T kernel sources (`tools/perf`), or try the Ubuntu
`linux-tools` package that matches the kernel version. Then check `perf stat -e cycles:u,instructions:u true` and
`cat /proc/sys/kernel/perf_event_paranoid` (must be ≤ 2). The Cortex-A57 PMU has the generic events
`run_profile.sh --perf` asks for; any it lacks show as `<not supported>`, and the run is not affected. Without `perf`,
leave out `--perf`. The `ONS_PROF` timers and counts do not need it.

**Comparing with x86.** On the Jetson, check that `bin/onsight` and `bin/onsight-prof` give the same md5s, as in
step 3. Do not expect the md5s from the x86 run. The aarch64 build lets GCC fuse multiply-adds
(`-ffp-contract=off` is set only for x86 in `src/Makefile`), and the two glibc maths libraries can differ in the last
bit, so read p-values near a cutoff can come out differently. Compare the target calls and per-read workload counts
with the x86 run instead. Index files built on the two machines differ only in 14 header bytes that hold memory
addresses.

## Note on timing under WSL2 or VMs

On WSL2 the kernel clock can run fast against real time. `prof_report.py` compares the profiler's wall time with the
real elapsed time in `time.txt` and scales the stage timers when they differ by more than 2%. The factor is printed in
the Run table. On a normal Linux machine nothing is changed.

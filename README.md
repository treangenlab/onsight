<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="doc/logo_dark.png">
    <img src="doc/logo.png" alt="onsight" width="480">
  </picture>
</p>

# onsight

**Pathogen detection on nanopore raw signal, during the run, on CPU only. No basecalling, no GPU.**

onsight reads the raw signal of nanopore reads and matches it against a panel: short DNA pieces (signatures) that
are specific to each target species. It reports which targets are in the sample and after how many reads each one
was called. It is built to run next to the sequencer on a small device.

```
raw signal → events → 4-bit quantization → strobe-events → seed hits in the index
           → diagonal bands → non-overlapping hit count → read p-value, read call
           → near-neighbor check → target call
```

Reads **FAST5, POD5, SLOW5 and BLOW5**, from R9.4.1 or R10.4.1 flow cells.

---

## Install

```bash
git clone https://github.com/treangenlab/onsight.git
cd onsight
git submodule update --init extern/slow5lib
conda env create -f environment.yml
conda activate onsight
make -j 8
make install                  # puts onsight on PATH, and the pore models next to it
onsight --version
```

The build takes under a minute. `make` downloads the POD5 library from GitHub, so it needs internet access once.

`make install` installs into the active conda environment (`$CONDA_PREFIX`), otherwise into `/usr/local`; change it
with `make install PREFIX=/your/path` and undo it with `make uninstall`. Without `make install`, run `./bin/onsight`.

`environment.yml` supplies everything the build needs:

* **Compiler.** POD5's libraries need GCC 11 or newer; many servers ship GCC 8, which fails at link time.
* **glibc.** `sysroot_linux-64=2.17` sets the oldest glibc the binary needs, so one build runs on RHEL 7 through
  current Ubuntu. This line is for x86-64; on ARM64 (for example NVIDIA Jetson) use `sysroot_linux-aarch64=2.17`.
* **HDF5 1.14.3 and ONT's VBZ plugin.** FAST5 files written by MinKNOW since about 2020 store the signal
  VBZ-compressed, and HDF5 reads them only with the plugin. The environment sets `HDF5_PLUGIN_PATH` when activated.
  POD5 and SLOW5/BLOW5 files do not need it.

### Without conda

```bash
git clone --recursive https://github.com/treangenlab/onsight.git
cd onsight
make -j 8
```

Builds HDF5 and zstd from `extern/` (several minutes more). Needs GCC 11+ and GNU Make. To read VBZ-compressed FAST5,
install ONT's [vbz_compression](https://github.com/nanoporetech/vbz_compression) plugin and set `HDF5_PLUGIN_PATH` to
its folder.

---

## Try the example

`example/` holds a small panel and 160 simulated reads, to check that onsight runs:

```bash
conda activate onsight
onsight --pore r10 --panel-dir example/panel -o example_out/ -t 4 example/reads.pod5
```

It takes a few seconds and about 0.3 GB of memory. In `example_out/target_calls.tsv`, the rows with rule
`read_pvalue+neighbor_check` should show *S. aureus* called at read 4 and *E. faecalis* at read 9
(`first_call_reads`); `first_call_time_s` depends on the machine.

| File | What it is |
|---|---|
| `example/panel/` | a two-target panel: the signatures and near-neighbor regions of *S. aureus* and *E. faecalis* in the 24-target panel, with an R10.4.1 index only |
| `example/reads.pod5` | 160 R10.4.1 reads simulated with squigulator, in run order: 20 each of *S. aureus*, *E. faecalis* and their near neighbors *S. epidermidis* and *E. faecium*, and 80 of two gut species (*Bacteroides difficilis*, *Faecalibacterium prausnitzii*) |
| `example/truth.tsv` | the species of each read, to compare with `example_out/read_calls.tsv` |

---

## Run

### 1. Get a panel

A panel is one folder. Download the 24-target panel from [Zenodo](https://zenodo.org/records/23967694), or build your
own (below), and keep it anywhere; you only pass its path. `example/panel` is a small one to try.

```bash
wget https://zenodo.org/records/23967694/files/panel24.tar.gz      # 2.7 GB
tar -xzf panel24.tar.gz                                             # makes panel24/
```

```
panel24/
├── panel.tsv                  the targets: which signature set belongs to which target name
├── regions.lcr                near-neighbor regions (the check that keeps close relatives from being called)
├── r9/index.idx               index for R9.4.1 reads, with r9/index.idx.sigmeta.tsv
└── r10/index.idx              index for R10.4.1 reads, with r10/index.idx.sigmeta.tsv
```

The download also has `signatures/` (the signatures as FASTA and table, which onsight does not read) and `md5sums.txt`
(`md5sum -c md5sums.txt` in the folder checks the files). The same record has the simulated reads of the paper,
`sim_r9.tar.gz` (27 GB) and `sim_r10.tar.gz` (38 GB): three sets per pore, each with its POD5, squigulator's
alignments and the species of every read.

### 2. Detect

```bash
conda activate onsight
onsight --pore r10 --panel-dir panel24/ -o results/ -t 4 reads/
```

Activate the environment first: it tells HDF5 where the VBZ plugin is, and most FAST5 files need it.

| Option | Meaning |
|---|---|
| `--pore r9` or `--pore r10` | the flow cell the reads come from (R9.4.1 or R10.4.1); **required**. It sets the pore model, the event detection and which index of the panel is used |
| `--panel-dir DIR` | the panel folder |
| `-o DIR` | where the results go (made if missing) |
| `-t INT` | threads (default 3) |
| reads | one or more FAST5, POD5, SLOW5 or BLOW5 files, or folders of them |

Reads are processed in the order given; a folder's files are read in the order the file system lists them. To get
first-call times that follow the run, give the files in sequencing order, for example as a sorted list:
`onsight ... $(cat files_in_order.txt)`. Very large single files (hundreds of thousands of reads) slow the run down;
split them, for example into 4,000 reads per file.

`onsight -h` lists every option. The defaults are the settings used in the paper; leave them unless you know why.

### 3. Read the results

**`results/target_calls.tsv`** is the answer: one row per target and rule. Use the rows with rule
`read_pvalue+neighbor_check`.

| Column | Meaning |
|---|---|
| `entry` | the target |
| `rule` | `read_pvalue`: read calls as they are; `read_pvalue+neighbor_check`: each read call counts with the larger of its read p-value and its near-neighbor p-value (the result to use) |
| `called` | 1: the target is in the sample; 0: not called |
| `first_call_reads` | reads processed when the target was called |
| `first_call_time_s` | seconds of processing when the target was called |
| `calls` | reads called for this target (read p-value ≤ 0.001); with `read_pvalue+neighbor_check`, those the check kept, so it stops growing soon after the target call, at the end of that batch of reads (see `skipped` below) |
| `reads` | reads processed in the whole run |
| `log_lr`, `threshold` | the test statistic and the value it must reach to call the target; set so that the chance of any false target call in a run is at most 1% |
| `test`, `alpha_read`, `t_min` | the test's settings (averaged over read p-value cutoffs from `t_min` to 0.001) |

**`results/read_calls.tsv`**: one row per read that matched a target.

| Column | Meaning |
|---|---|
| `read` | read name |
| `position` | the read's number in processing order, from 0 |
| `entry` | the target the read matched best |
| `p` | read p-value: the chance of a match this good by chance; a read call is p ≤ 0.001 |
| `nonoverlap_hits`, `hits` | seed hits on the best signature, without and with overlapping ones |
| `read_events`, `sig_len` | events in the read; length of the signature in bases |
| `n_sigs` | signatures of this target in the panel |
| `slopes` | slopes searched (read speed against the signature) |
| `rid`, `signature` | the signature matched |
| `neighbor_check` | `kept` when the larger of `p` and `neighbor_p` is at most 0.001, else `dropped`; `no band` when the match gave no region to check; `skipped` when the target was already called, since the check can no longer change a target call (`--neighbor-check-all` checks these too, more slowly); `-` when not checked (p > 0.001, or the panel has no `regions.lcr`) |
| `neighbor_p` | near-neighbor p-value: the chance that a read of a near neighbor fits the target this much better than every near neighbor; 0 when no near neighbor has the region. Worked out only where it can change the call, from the larger of `p` and 1e-8 up to 0.001: at the lower end it means "at most", above 0.001 "at least"; `-` when not checked or skipped |

**`results/sig_summary.tsv`** (only with `--debug`): one row per read: its best signature, its top five signatures and
their scores. Use it to look into what a read matched. `--debug` sorts every read's signatures, so the run is slower;
`read_calls.tsv` and `target_calls.tsv` do not change.

### Memory and speed

With the 24-target panel and 4 threads, peak memory was 2.5–3.1 GB on an Intel Xeon Gold 6230 server and 2.3–3.3 GB on
an NVIDIA Jetson Orin Nano; most of it is the index and `regions.lcr`, loaded whole. On the server, keeping pace with
the sequencing of three Zymo mock runs (R9.4.1) took 0.31–0.72 CPU cores. `-K` (minibatch size, default 4M signal
samples) also sets memory, but it changes the results, because the slope grid is updated between minibatches; keep
the default.

### Building a panel

A panel is built offline from genomes, once, with the Python pipeline in [`panel-src/`](panel-src/README.md). The
signatures come from [Seqwin](https://github.com/treangenlab/Seqwin): for each target it finds sequences that are
frequent in the target's genomes and rare in its relatives' genomes, from a pan-genome minimizer graph. `panel-src`
picks the genomes, runs Seqwin, removes signatures that other species carry, finds each target's near neighbors and
their homologous regions, and writes the panel folder above:

```bash
conda env create -f panel-src/environment.yml && conda env create -f panel-src/seqwin.yml
conda activate onsight-panel
python3 panel-src/run_panel.py --config mypanel.tsv --paths mypaths.tsv --root build/
onsight --pore r10 --panel-dir build/mypanel/onsight -o results/ reads/
```

[`panel-src/README.md`](panel-src/README.md) walks through a two-target example step by step (inputs, genome
downloads, the example config `panel-src/panels/example.tsv` and paths file `example_paths.tsv`), then lists every
config column and step.

---

## Build options

All signal formats are enabled by default.

| Flag | Default | Disables |
|---|---|---|
| `ENABLE_POD5` | `ON` | `.pod5` |
| `ENABLE_HDF5` | `ON` | `.fast5` |
| `ENABLE_SLOW5` | `ON` | `.slow5`, `.blow5` |

`FP_CONTRACT=off` on ARM (for example NVIDIA Jetson) turns off fused multiply-add, as x86 builds always do. The
default keeps the compiler's setting, as in the paper's Jetson runs. Fused multiply-add changes the last digits of
p-values, and through the slope grid some read calls: an x86 build with it (`FP_CONTRACT=fast`) gives the Jetson's
read calls. With `off`, ARM should give x86's results; this has not been run on ARM yet.

```bash
make ENABLE_HDF5=OFF                          # drop .fast5
make ENABLE_POD5=OFF ENABLE_HDF5=OFF          # pure C99, SLOW5/BLOW5 only
```

Disabling all three is a build error. When `CONDA_PREFIX` is set, the build uses conda's HDF5
and zstd rather than building `extern/hdf5` and `extern/zstd`.

Debug and sanitizer builds: `make DEBUG=1`, `make asan=1`, `make tsan=1`. Per-stage timers and work counts: `make prof`
(`bin/onsight-prof`, same results; see `src/onsprof.h` and `scripts/hw_prof/`).

The version is the number in `VERSION` (`onsight --version`).

> **Platform note:** POD5 ships prebuilt libraries only for linux-x64, linux-arm64, and
> osx-arm64. On Windows, Intel macOS, or musl, build with `ENABLE_POD5=OFF`.

---

## License

GPL-3.0, see [LICENSE](LICENSE). onsight is a modified version of [RawHash2](https://github.com/CMU-SAFARI/RawHash)
(GPL-3.0), changed in 2025-2026. It includes klib, minimap2's FASTA reader, thread and memory-pool code, and
`hdf5_tools.hpp`, under the MIT license ([src/third_party/LICENSE](src/third_party/LICENSE)). If you use onsight,
please also cite RawHash and RawHash2 ([doc/CITATION.md](doc/CITATION.md)).

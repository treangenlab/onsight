# panel-src: build a panel for onsight

`run_panel.py` turns a list of targets into the panel folder onsight reads (`onsight --panel-dir`):

```
<root>/<panel>/onsight/
├── panel.tsv           which signature set belongs to which target
├── regions.lcr         near-neighbor regions for the near-neighbor check
├── r9/index.idx        R9.4.1 index (and r9/index.idx.sigmeta.tsv)
└── r10/index.idx       R10.4.1 index (and r10/index.idx.sigmeta.tsv)
```

**Signatures.** A target's signatures are short sequences that are frequent in its genomes and rare in the genomes of
its relatives. They come from [Seqwin](https://github.com/treangenlab/Seqwin), which builds a weighted pan-genome
minimizer graph of the target and non-target genomes and walks it to find such sequences ("Seqwin: Ultrafast
identification of signature sequences in microbial genomes", bioRxiv 2025, doi:10.1101/2025.11.07.687294). We run
Seqwin 0.3.1 with stringency 3 (`signatures/run_seqwin.py`).

**Near neighbors.** Each target's near neighbors are genomes of its family, or at 80% ANI or more to its own
genomes, that are not the target. The other targets of the panel compete too, through their own strains. For every
signature, blastn finds the matching regions in these genomes and in the target's own strains. onsight keeps a read
call only when the read's event levels fit the target better than every near neighbor.

## Tutorial: a two-target panel

This builds a panel of *Klebsiella pneumoniae* and *Staphylococcus aureus* from the example files in `panels/`
(`example.tsv`, `example_paths.tsv`). Everything goes in one working folder; `<onsight>` is your clone of this repo.

**1. Environments.** onsight itself is built as in the main README (`make && make install`); the index step calls it.

```bash
conda env create -f <onsight>/panel-src/environment.yml   # onsight-panel: Python, numpy, skani, minimap2, BLAST, NCBI datasets
conda env create -f <onsight>/panel-src/seqwin.yml        # seqwin 0.3.1, used by the seqwin step only
conda activate onsight-panel
mkdir work && cd work
```

**2. NCBI taxonomy** (species names and ranks):

```bash
mkdir taxdump && wget -qO- https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | tar -xz -C taxdump
```

**3. RefSeq genomes and the background pack.** The build screens signatures against one genome per species of
RefSeq (the background pack) and takes each target's own strains from RefSeq complete genomes. It needs
`refseq/assembly_summary_refseq.txt` (from `https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt`)
and the genomes as one FASTA per assembly, named as on the RefSeq FTP site (`<accession>_<name>_genomic.fna.gz`,
under the `ftp_path` of each assembly in that table), in `refseq/library/bacteria/` (and `archaea/`, `fungi/`, `viral/`
for a pack across all of them, as ours: 27,592 genomes, 66 GB). Then:

```bash
python3 <onsight>/panel-src/background/background_pack.py --library refseq/library \
    --summary refseq/assembly_summary_refseq.txt --out background            # writes background.tsv
```

**4. Genomes of each target and its relatives.** One command per signature set makes both lists: it asks NCBI for
every GenBank assembly under `--scope`; the assemblies of the `--target` species are the targets, and those of every
other named species under `--scope` are the non-targets (at most `--neg-cap` per species, picked at random). It
downloads them and writes `inputs/<build>/assemblies.tsv`:

```bash
for t in "Klebsiella pneumoniae:Klebsiella" "Staphylococcus aureus:Staphylococcus"; do
  species=${t%%:*}; scope=${t##*:}; build=${species// /_}_vs_genus
  python3 <onsight>/panel-src/genome_selection/seqwin_lists.py --target "$species" --scope $scope --neg-cap 50 \
      --taxdump taxdump --out inputs/$build --genomes genomes/$build
done
```

Use the parent taxon as `--scope`: the genus for a species target (here *K. variicola*, *K. quasipneumoniae* and
the other *Klebsiella* species are the non-targets of *K. pneumoniae*), the family for a genus target. A genus
target such as *Enterobacter* is built one species at a time, each against the family, so the other species of the
genus are non-targets too (rows that share a `name` in the config are one target; see `panels/eskape.tsv`).

**5. Config and paths.** Copy the examples; the config needs no change for this panel, and the paths file only needs
the seqwin environment's Python (`conda run -n seqwin which python3`):

```bash
cp <onsight>/panel-src/panels/example.tsv mypanel.tsv
cp <onsight>/panel-src/panels/example_paths.tsv mypaths.tsv     # edit seqwin_python
```

**6. Build and use.** The panel is named after the config file (`mypanel`):

```bash
python3 <onsight>/panel-src/run_panel.py --config mypanel.tsv --paths mypaths.tsv --root build/ --threads 8 --processes 8
onsight --pore r10 --panel-dir build/mypanel/onsight -o results/ reads/
```

A stopped build continues where it was: work already there is skipped. Each step's output is in
`build/mypanel/logs/<step>.log`.

## Config

One row per signature set (tab-separated). `panels/example.tsv` is the smallest; `panels/eskape.tsv` and
`panels/zymo_log.tsv` are the paper's panels, with genus targets, lineage sets and the near-neighbor options.

| Column | Meaning | Example |
|---|---|---|
| `set` | the set's name in the panel | `Klebsiella_pneumoniae` |
| `build` | its build folder; one Seqwin run | `Klebsiella_pneumoniae_vs_genus` |
| `short` | the species, `Genus_species` | `Klebsiella_pneumoniae` |
| `name` | the target name calls are reported under; rows that share it are one target | `K. pneumoniae` |
| `input` | the folder under `inputs` with `assemblies.tsv` (`seqwin_lists.py --out`) | `Klebsiella_pneumoniae_vs_genus` |
| `target` | the target species, only when the input holds several (a genus input) | empty |
| `neg_cap`, `target_cap` | at most this many non-target genomes per species, and target genomes; empty: no cap | `50`, `1000` |
| `lineages_of`, `lineages` | for an extra set of a species with several lineages (95% ANI clusters): the build it is split from, and `other` or `all`; empty otherwise | empty |
| `report_as`, `panel`, `mock_label` | labels kept in the set's manifest | `Klebsiella pneumoniae`, `example`, empty |
| `reference_genome` | the species' NCBI reference genome accession; empty: the one the genome selection used | empty |
| `lineage_representative` | a lineage set's representative genome; empty otherwise | empty |
| `keep_off` | empty (the other targets of this panel are not near neighbors; they compete through their own strains) | empty |
| `test_strains` | genomes your test reads came from (`<paths key>/<file>`, comma separated): own strains at 99.9% ANI or more to one are skipped, so tests stay honest | empty |
| `extra_near_strains` | a list (file in `inputs`) of genomes test reads came from: this target gets up to 10 genomes per near-neighbor species, none at 99.9% ANI to one of them | empty |

## Paths

Key and value per line (tab-separated); `panels/example_paths.tsv` is a filled-in example. Relative paths are read
from the paths file's folder, and a bare name (`skani`, `onsight`) is a program on PATH.

| Key | Value |
|---|---|
| `inputs` | the folder of `seqwin_lists.py` outputs (step 4) |
| `background` | `background.tsv` (step 3) |
| `taxdump` | the NCBI taxonomy folder (step 2) |
| `refseq_summary`, `refseq_libraries` | `assembly_summary_refseq.txt`, and the genome folders own strains come from (comma separated) |
| `seqwin_python` | the seqwin environment's `python3` |
| `onsight` | the onsight binary |
| `skani`, `blastn`, `makeblastdb` | the tools |
| `mock_genomes`, `mock_strains` | mock community genomes, as a check only; may be empty |

## Steps

| Step | Script | Writes |
|---|---|---|
| lists | `genome_selection/select_genomes.py`, `lineage_lists.py` | `sets/lists/<build>/`: targets (high quality, 95% ANI to the reference), non-targets, held-out genomes |
| seqwin | `signatures/run_seqwin.py` | `sets/seqwin/<build>/`: candidate signatures |
| screen | `signatures/screen_signatures.py` | `sets/screen/<build>.*`: background genomes that carry each signature (minimap2) |
| sets | `signatures/final_sets.py` | `sets/final/<set>/`: signatures no background genome carries |
| panel | `signatures/build_panel.py` | `<panel>/work/signatures/`, `<panel>/onsight/panel.tsv` |
| near | `near_neighbors/near_neighbors.py` | `<panel>/work/near_neighbors/`: `targets.tsv`, candidates, one genome per species group |
| own | `near_neighbors/own_strains.py` | `<panel>/work/own_strains/`: up to 10 RefSeq complete genomes per target species |
| extra | `near_neighbors/near_neighbor_strains.py` | `<panel>/work/near_strains/<target>/` (targets with `extra_near_strains`) |
| others | `near_neighbors/other_targets.py` | `<panel>/work/other_targets/`: other targets' own strains on each target's near-neighbor side |
| blast | `near_neighbors/blast_homologs.py` | `paf/` in each folder above: every blastn hit of every signature |
| index | `index/build_index.sh` | `<panel>/onsight/r9/`, `<panel>/onsight/r10/` |
| regions | `near_neighbors/neighbor_regions.py` | `<panel>/onsight/regions.lcr` |

`sets/` is shared by every panel built in the same root, so panels that share targets reuse their Seqwin runs.
`--steps` runs only some steps (comma separated).

**Time and space.** Seqwin and the blast step take most of the time: from minutes for a small genus to several hours
for a species with thousands of genomes (our 24-target panel took about a day on 48 cores). The 24-target panel
folder is 3.5 GB; the build root with genomes is tens of GB.

**Region file format.** `regions.lcr` is written as LCR1 (32-bit base positions) when it fits, and as LCR2 (64-bit)
past 4.29 billion stored bases, as the 24-target panel's is. onsight 0.0.2 reads both; 0.0.1 reads LCR1 only.

## Layout of `panel-src/`

| Folder | Holds |
|---|---|
| `background/` | `background_pack.py`: the background genome pack |
| `genome_selection/` | genome lists and downloads (`seqwin_lists.py`), genome selection with the 95% ANI species rule (`select_genomes.py`), lineage sets (`lineage_lists.py`), and a check that a sample's strain is the species its set targets (`seqwin_strain_check.py`) |
| `signatures/` | Seqwin (`run_seqwin.py`), the screen, the final sets, the panel's signature files |
| `near_neighbors/` | near neighbors, own strains, extra strains, other targets, blastn, the region file |
| `index/` | `build_index.sh`: one index (calls onsight) |
| `panels/` | panel configs (`example.tsv`, `eskape.tsv`, `zymo_log.tsv`), `example_paths.tsv`, `mock_strains.tsv` |

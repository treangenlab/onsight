# Citing

onsight builds on RawHash and RawHash2 (<https://github.com/CMU-SAFARI/RawHash>). If you use onsight, please also
cite them:

```bibtex
@article{firtina_rawhash_2023,
	title = {{RawHash}: enabling fast and accurate real-time analysis of raw nanopore signals for large genomes},
	author = {Firtina, Can and Mansouri Ghiasi, Nika and Lindegger, Joel and Singh, Gagandeep and Cavlak, Meryem Banu and Mao, Haiyu and Mutlu, Onur},
	journal = {Bioinformatics},
	volume = {39},
	number = {Supplement_1},
	pages = {i297-i307},
	month = jun,
	year = {2023},
	doi = {10.1093/bioinformatics/btad272},
	issn = {1367-4811},
	url = {https://doi.org/10.1093/bioinformatics/btad272},
}

@article{firtina_rawhash2_2024,
	title = {{RawHash2}: mapping raw nanopore signals using hash-based seeding and adaptive quantization},
	volume = {40},
	issn = {1367-4811},
	url = {https://doi.org/10.1093/bioinformatics/btae478},
	doi = {10.1093/bioinformatics/btae478},
	number = {8},
	journal = {Bioinformatics},
	author = {Firtina, Can and Soysal, Melina and Lindegger, Joël and Mutlu, Onur},
	month = aug,
	year = {2024},
	pages = {btae478},
}
```

## Acknowledgement

onsight keeps RawHash2's signal input, event detection (with the R9.4.1 segmentation parameters of
[Sigmap](https://github.com/haowenz/sigmap)), quantization, index and seed lookup, and uses
[klib](https://github.com/attractivechaos/klib) and code from [minimap2](https://github.com/lh3/minimap2) (the FASTA
reader, the thread pipeline, the memory pool); their MIT notices are in `src/third_party/LICENSE`.

The pore models in `extern/local_kmer_models/` come with RawHash2: `r94_means_only.txt` holds the level means of ONT's
R9.4 6-mer model (`legacy_r9.4_180mv_450bps_6mer`, [kmer_models](https://github.com/nanoporetech/kmer_models)), and
`uncalled_r1041_model_only_means.txt` the R10.4.1 9-mer model of [Uncalled4](https://github.com/skovaka/uncalled4).

The panel build (`panel-src/`) finds signatures with [Seqwin](https://github.com/treangenlab/Seqwin).

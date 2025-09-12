# AdapterRemoval v3 - Dang fast(Q) processing

[![build](https://github.com/MikkelSchubert/adapterremoval/actions/workflows/build-and-test.yaml/badge.svg)](https://github.com/MikkelSchubert/adapterremoval/actions/workflows/build-and-test.yml) [![docs](https://readthedocs.org/projects/adapterremoval/badge/?version=latest)](https://adapterremoval.readthedocs.io/)

AdapterRemoval trims adapter sequences and low quality bases from High-Throughput Sequencing (HTS) data in FASTQ format. For paired-end data, AdapterRemoval can merge overlapping paired-ended reads into (longer) consensus sequences. Additionally, AdapterRemoval can demultiplex FASTQ reads, and construct a consensus adapter sequence for paired-ended reads, if this information is not available.

AdapterRemoval v3 is a major revision of AdapterRemoval v2 that aims at simplifying usage via sensible default settings. In addition, v3 adds support for detailed QC reports, trimming of poly-X tails and improved trimming of low quality bases, and greatly increased throughput. See below for more information.

AdapterRemoval v3 is still a work in progress, but alpha release 3 is [available for download](https://github.com/MikkelSchubert/adapterremoval/releases/tag/v3.0.0-alpha3/). Documentation is available at [Read the Docs](https://adapterremoval.readthedocs.io/en/v3.0.0-alpha3/), including a guide on how to migrate from v2. For questions, bug reports, and/or suggestions, please use the [GitHub tracker](https://github.com/MikkelSchubert/adapterremoval/issues/).

## Major features

- Trimming of adapters sequences from single-end and paired-end FASTQ reads
  - Trimming of multiple, different adapters or adapter pairs
- Detailed [human-readable](https://mikkelschubert.github.io/adapterremoval/examples/3.0.0-alpha3.html) and [machine-readable](https://mikkelschubert.github.io/adapterremoval/examples/3.0.0-alpha3.json) QC reports **(v3)**
  - The ability to perform QC-only runs with or without read processing **(v3)**
- Barcode based demultiplexing with or without trimming of adapter sequences
  - Support for samples identified by multiple barcode pairs **(v3)**
  - Support for mixed orientation barcodes **(v3)**
- Support for multiple methods for trimming low quality bases/reads
  - Quality trimming using windows or constants thresholds
  - Quality trimming using the modified Mott algorithm **(v3)**
  - Poly-X tail trimming, supporting any combination of trailing bases **(v3)**
- Filtering of reads based on complexity **(v3)**
- Reconstruction of adapter sequences by pair-wise alignment of paired-end reads
- Merging of overlapping read-pairs into higher-quality consensus sequences
- Support for reading interleaved FASTQ files
- Support for arbitrary splitting/interleaving of output files **(v3)**
- Support for writing BGZF compressed SAM and BAM files **(v3)**
- Support for SSE2, AVX2, AVX512, and NEON accelerated alignments **(v3)**

## Performance

AdapterRemoval v3 features greatly increased throughput compared to AdapterRemoval v2. This is accomplished through support for additional SIMD instruction sets, improved parallelization of I/O and computationally expensive tasks, including block based compression of output files, as well as defaulting to a lower compression level.

Please note that these results are preliminary:

![Throughput for ARv2, ARv3, and fastp](https://raw.githubusercontent.com/MikkelSchubert/adapterremoval/main/docs/images/throughput.svg)

Point labels indicate the number of worker threads configured for each program, while the X-axis indicates observed CPU-usage for a given number of worker threads. The Y-axis indicates millions of 150bp paired-end reads processed per second, for gzipped input and output, with merging enabled and duplication estimation disabled.

Benchmarking was performed on an Intel i9-11900K with 8 physical cores, and plotting is therefore limited to ~8 CPUs. The CPU usage of fastp being higher than the number of worker threads, is due to additional threads being used (compressed) I/O. ARv2 does not scale past 4 threads.

## Documentation

For a detailed description of program installation and usage, please refer to the [online documentation](https://adapterremoval.readthedocs.io/en/latest). A summary of command-line options may also be found in the [manual page](https://adapterremoval.readthedocs.io/en/latest/manpage.html), accessible via the command `man adapterremoval3` once AdapterRemoval has been installed.

## Citation

If you use AdapterRemoval v3, then please cite the paper::

    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 12;9(1):88 <http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2>

AdapterRemoval was originally published in Lindgreen 2012:

    Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation Sequencing Reads, BMC Research Notes, 5:337 <http://www.biomedcentral.com/1756-0500/5/337/>

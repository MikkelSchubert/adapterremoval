# AdapterRemoval

[![build](https://github.com/MikkelSchubert/adapterremoval/actions/workflows/build-and-test.yaml/badge.svg)](https://github.com/MikkelSchubert/adapterremoval/actions/workflows/build-and-test.yml) [![coverage](https://coveralls.io/repos/github/MikkelSchubert/adapterremoval/badge.svg?branch=master)](https://coveralls.io/github/MikkelSchubert/adapterremoval) [![docs](https://readthedocs.org/projects/adapterremoval/badge/?version=latest)](https://adapterremoval.readthedocs.io/)

AdapterRemoval trims adapter sequences and low quality bases from High-Throughput Sequencing (HTS) data in FASTQ format. For paired-end data, AdapterRemoval can merge overlapping paired-ended reads into (longer) consensus sequences. Additionally, AdapterRemoval can demultiplex FASTQ reads, and construct a consensus adapter sequence for paired-ended reads, if this information is not available.

For questions, bug reports, and/or suggestions, please use the [GitHub tracker](https://github.com/MikkelSchubert/adapterremoval/issues/).

## AdapterRemoval v3 - Dang fast(Q) processing

AdapterRemoval v3 is a major revision of AdapterRemoval, that aims to simplify usage by picking a sensible set of default settings, providing [human-readable](https://mikkelschubert.github.io/adapterremoval/examples/3.0.0-alpha3.html) and [machine-readable](https://mikkelschubert.github.io/adapterremoval/examples/3.0.0-alpha3.json) QC reports, as well as greatly improving overall throughput, and more.

AdapterRemoval v3 is still a work in progress, but alpha release 3 is [available for download](https://github.com/MikkelSchubert/adapterremoval/releases/tag/v3.0.0-alpha3/). Documentation is available at [Read the Docs](https://adapterremoval.readthedocs.io/en/v3.0.0-alpha3/), including a guide on how to migrate from v2. Bug reports, feature requests, and other feedback is greatly appreciated.

For more information, please see [README-v3.md](README-v3.md).

## AdapterRemoval v2

If you use AdapterRemoval v2, then please cite the paper:

    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 12;9(1):88 <http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2>

AdapterRemoval was originally published in Lindgreen 2012:

    Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation Sequencing Reads, BMC Research Notes, 5:337 <http://www.biomedcentral.com/1756-0500/5/337/>

### Overview of major features

- Trimming of adapters sequences from single-end and paired-end FASTQ reads.
- Trimming of multiple, different adapters or adapter pairs.
- Demultiplexing of single or double indexed reads, with or without trimming
      of adapter sequences.
- Reconstruction of adapter sequences from paired-end reads, by the pairwise
      alignment of reads in the absence of a known adapter sequence.
- Merging of overlapping read-pairs into higher-quality consensus sequences.
- Multi-threading of all operations for increased throughput.
- Reading and writing of gzip and bzip2 compressed files.
- Reading and writing of interleaved FASTQ files.

### Documentation

For a detailed description of program installation and usage, please refer to the [online documentation](https://adapterremoval.readthedocs.io/). A summary of command-line options may also be found in the manual page, accessible via the command "man AdapterRemoval" once AdapterRemoval has been installed.

### Installation

#### Installation with Conda

If you have Conda installed on your system:

   conda install -c bioconda adapterremoval

#### Installing from sources

Installing AdapterRemoval from sources requires libz and libbz2.

To compile AdapterRemoval, download the latest release, unpack the archive and then simply run "make" in the resulting folder:

    wget -O adapterremoval-2.3.4.tar.gz <https://github.com/MikkelSchubert/adapterremoval/archive/v2.3.4.tar.gz> tar xvzf adapterremoval-2.3.4.tar.gz cd adapterremoval-2.3.4 make

The resulting 'AdapterRemoval' executable is located in the 'build' subdirectory and may be installed by running "make install":

    sudo make install

### Getting started

To run AdapterRemoval, specify the location of pair 1 and (optionally) pair 2 FASTQ using the --file1 and --file2 command-line options:

    AdapterRemoval --file1 reads_1.fastq.gz --file2 reads_2.fastq.gz

By default, AdapterRemoval will save the trimmed reads in the current working directly, using filenames starting with 'your_output'.

More examples of common usage may be found in the [Examples](https://adapterremoval.readthedocs.io/en/latest/examples.html) section of the online documentation:

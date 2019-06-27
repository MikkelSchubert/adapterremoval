# AdapterRemoval [![Travis](https://img.shields.io/travis/MikkelSchubert/adapterremoval/master.svg)](https://travis-ci.org/MikkelSchubert/adapterremoval) [![Coveralls](https://img.shields.io/coveralls/MikkelSchubert/adapterremoval.svg)](https://coveralls.io/github/MikkelSchubert/adapterremoval)

AdapterRemoval searches for and removes adapter sequences from High-Throughput
Sequencing (HTS) data and (optionally) trims low quality bases from the 3' end
of reads following adapter removal. AdapterRemoval can analyze both single end
and paired end data, and can be used to merge overlapping paired-ended reads
into (longer) consensus sequences. Additionally, AdapterRemoval can construct a
consensus adapter sequence for paired-ended reads, if which this information is
not available.

For comments, suggestions  and feedback please contact Mikkel Schubert
(MikkelSch@gmail.com) and Stinus Lindgreen (stinus@binf.ku.dk).

If you use AdapterRemoval v2, then please cite the paper:

    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter
    trimming, identification, and read merging. BMC Research Notes, 12;9(1):88
    http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2

AdapterRemoval was originally published in Lindgreen 2012:

    Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation
    Sequencing Reads, BMC Research Notes, 5:337
    http://www.biomedcentral.com/1756-0500/5/337/


## Overview of major features

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


## Documentation

For a detailed description of program installation and usage, please refer to
the [online documentation](https://adapterremoval.readthedocs.io/). A summary
of command-line options may also be found in the manual page, accessible via
the command "man AdapterRemoval" once AdapterRemoval has been installed.


## Installation

### Installation with Conda

If you have `Conda`_ installed on your system:

    conda install -c bioconda adapterremoval


### Installing on OSX

MacOSX users may install AdapterRemoval using [Homebrew](https://brew.sh/):

    brew install homebrew/science/adapterremoval


### Installing from sources

Installing AdapterRemoval from sources requires libz and libbz2.

To compile AdapterRemoval, download the latest release, unpack the archive and
then simply run "make" in the resulting folder:

    wget -O adapterremoval-2.3.1.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v2.3.1.tar.gz
    tar xvzf adapterremoval-2.3.1.tar.gz
    cd adapterremoval-2.3.1
    make

The resulting 'AdapterRemoval' executable is located in the 'build'
subdirectory and may be installed by running "make install":

    sudo make install


## Getting started

To run AdapterRemoval, specify the location of pair 1 and (optionally) pair 2
FASTQ using the --file1 and --file2 command-line options:

    AdapterRemoval --file1 myreads_1.fastq.gz --file2 myreads_2.fastq.gz

By default, AdapterRemoval will save the trimmed reads in the current working
directly, using filenames starting with 'your_output'. 


More examples of common usage may be found in the
[Examples](https://adapterremoval.readthedocs.io/en/latest/examples.html)
section of the online documentation:


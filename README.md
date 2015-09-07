==============
AdapterRemoval
==============

This program searches for and removes remnant adapter sequences from your read data and (optionally) trims low quality bases from the 3' end of reads following adapter removal.  The program can analyze both single end and paired end data, and can be used to merge overlapping paired-ended reads into (longer) consensus sequences. Additionally, the the program may to recover a consensus adapter sequence for paired-ended data.

For detailed explanation of the parameters, please refer to the man page.  For comments, suggestions  and feedback please contact Mikkel Schubert (MikkelSch@gmail.com) and Stinus Lindgreen (stinus@binf.ku.dk).

If you use the program, please cite the paper:
    S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation Sequencing Reads, BMC Research Notes, 5:337
    http://www.biomedcentral.com/1756-0500/5/337/


Installation
============

Note that AdapterRemoval requires that the zlib library and headers (www.zlib.net) are installed, and that the pthread library and headers are installed. Please refer to your operating system documentation for installation instructions.


Download and unpack the newest release from GitHub:

    $ wget -O adapterremoval-2.0.0.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v2.0.0.tar.gz
    $ tar xvzf adapterremoval-2.0.0.tar.gz
    $ cd adapterremoval-2.0.0

or

    $ git clone https://github.com/MikkelSchubert/adapterremoval.git
    $ cd adapterremoval

To compile, run

    $ make

The resulting binary and man page is located in the "build" folder.

To install, run

    $ sudo make install


Usage
=====

For program usage, please refer to the manual page.

If AdapterRemoval has been installed, this may be accessed using the command "man AdapterRemoval". If AdapterRemoval has not been installed, the manual page may be read using the command "man build/AdapterRemoval.1" in the source folder once "make" has been run. Alternatively, the manual may be read on GitHub:

https://github.com/MikkelSchubert/adapterremoval/blob/master/AdapterRemoval.pod


Demultiplexing
==============

As of version 2.1, AdapterRemoval supports simultanious demultiplexing and adapter trimming; demultiplexing is carried out using a simple comparison between the specified barcode sequences and the first N bases of the reads, corresponding to the length of the barcodes. Reads identified as containing a specific barcode or pair of barcodes are then trimmed using adapter sequences including these barcodes. Demultiplexing is accomplished by creating a table of barcodes, the first column of which species the sample name (using characters [a-zA-Z0-9_]) and the second and (optional) third columns specifies the mate 1 and mate 2 barcode sequences.

For example, a table of barcodes from a double-indexed run might be

    sample_1 ATGCGGA TGAATCT
    sample_2 ATGGATT ATAGTGA
    sample_7 CAAAACT TCGCTGC

AdapterRemoval is invoked with the --barcode-list option, specifying the path to this table, and generates a set of output files for each sample thus specified, using the basename (--basename) as the prefix, followed by a dot and the sample name, followed by a dot and the default name for a given file type. For example, the output files for sample_2 could be

    your_output.sample_2.discarded
    your_output.sample_2.pair1.truncated
    your_output.sample_2.pair2.truncated
    your_output.sample_2.settings
    your_output.sample_2.singleton.truncated

The settings files generated for each sample summarizes the reads for that sample only; in addition, a basename.settings file is generated which summarizes the number and proportion of reads identified as belonging to each sample.

The maximum number of mismatches allowed when comparing barocdes is controlled using the options --barcode-mm, --barcode-mm-r1, and --barcode-mm-r2, which specify the maximum number of mismatches total, and the maximum number of mismatches for the mate 1 and mate 2 barcodes respectively. Thus, if mm_1(i) and mm_2(i) represents the number of mismatches observed for barcode-pair i for a given pair of reads, these options require that

   1. mm_1(i) <= --barcode-mm-r1
   2. mm_2(i) <= --barcode-mm-r2
   3. mm_1(i) + mm_2(i) <= --barcode-mm


A note on specifying adapter sequences
======================================

Please note that the --pcr1 and --pcr2 options used with AdapterRemoval v1.x have been deprecated in favor of options --adapter1 and --adapter2. For both --adapter1 and --adapter2 the adapter sequence are expected to be observed in the raw mate 1 and mate 2 reads respectively, exactly as specified on the command-line, which corresponds to the behavior of most adapter trimming programs.

Default adapter #1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG

Default adapter #2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

Assuming these were the adapters used to generate our data, we should therefore see these in the FASTQ files (typically followed by a low-quality A-tail), when ignoring any difference in case and treating Ns as wildcards:

    $ grep -i "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC......ATCTCGTATGCCGTCTTCTGCTTG" file1.fq
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAACAAGAAT
    CTGGAGTTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGAATCTCGTATGCCGTCTTCTGCTTGAAAAAAA
    GGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGAATCTCGTATGCCGTCTTCTGCTTGCAAATTGAAAACAC
    ...

    $ grep -i "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" file2.fq
    CAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTCAAAAAAAGAAAAACATCTTG
    GAACTCCAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTCAAAAAAAATAGA
    GAACTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTCAAAAACATAAGACCTA
    ...

The options --pcr1 and --adapter1 are functionally equivalent, while the option --pcr2 expects the reverse complement of the --adapter2 sequence. Thus, the default for --pcr2 is AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, the reverse complement of the default for --adapter2.

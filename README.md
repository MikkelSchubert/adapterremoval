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

# AdapterRemoval [![Travis](https://img.shields.io/travis/MikkelSchubert/adapterremoval/master.svg)](https://travis-ci.org/MikkelSchubert/adapterremoval) [![Coveralls](https://img.shields.io/coveralls/MikkelSchubert/adapterremoval.svg)](https://coveralls.io/github/MikkelSchubert/adapterremoval)

This program searches for and removes remnant adapter sequences from High-Throughput Sequencing (HTS) data and (optionally) trims low quality bases from the 3' end of reads following adapter removal. AdapterRemoval can analyze both single end and paired end data, and can be used to merge overlapping paired-ended reads into (longer) consensus sequences. Additionally, the AdapterRemoval may be used to recover a consensus adapter sequence for paired-ended data, for which this information is not available.

For comments, suggestions  and feedback please contact Mikkel Schubert (MikkelSch@gmail.com) and Stinus Lindgreen (stinus@binf.ku.dk). If you use AdapterRemoval v2, then please cite the paper:

    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 12;9(1):88
    http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2

AdapterRemoval was originally published in Lindgreen 2012:

    Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation Sequencing Reads, BMC Research Notes, 5:337
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


## Installation

To install, first download and unpack the newest release from GitHub:

    $ wget -O adapterremoval-2.1.7.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v2.1.7.tar.gz
    $ tar xvzf adapterremoval-2.1.7.tar.gz
    $ cd adapterremoval-2.1.7

or

    $ git clone https://github.com/MikkelSchubert/adapterremoval.git
    $ cd adapterremoval

To compile, run

    $ make

The resulting binary and man page is located in the "build" folder.

To install, run

    $ sudo make install


It is also possible to compile AdapterRemoval as a static library:

    $ sudo make static


Note that AdapterRemoval requires that the zlib library and headers (www.zlib.net) are installed, that the bzlib2 library and headers are installed, and that the compiler used supports c++11. Please refer to your operating system documentation for installation instructions.


## Documentation

For detailed program usage, please refer to the manual page.

If AdapterRemoval has been installed, this may be accessed using the command "man AdapterRemoval". If AdapterRemoval has not been installed, the manual page may be read using the command "man build/AdapterRemoval.1" in the source folder once "make" has been run. Alternatively, the manual may be read online:

https://github.com/MikkelSchubert/adapterremoval/blob/master/AdapterRemoval.pod



## Examples

The following examples make use of the data included in the 'examples' folder:



### Trimming single-end reads

The following command removes adapters from the file 'reads\_1.fq' trims both Ns and low quality bases from the reads, and gzip compresses the resulting files. The --basename option is used to specify the prefix for output files.

    $ AdapterRemoval --file1 reads_1.fq --basename output_single --trimns --trimqualities --gzip

Since --gzip and --basename is specified, the trimmed FASTQ reads are written to 'output_single.truncated.gz', the discarded FASTQ reads are written to 'output_single.discarded.gz', and settings and summary statistics are written to 'output_single.settings'.

Note that by default, AdapterRemoval does not require a minimum number of bases overlapping with the adapter sequence, before reads are trimmed. This may result in an excess of very short (1 - 3 bp) 3' fragments being falsely identified as adapter sequences, and trimmed. This behavior may be changed using the --minadapteroverlap option, which allows the specification of a minimum number of bases (excluding Ns) that must be aligned to carry trimming. For example, use --minadapteroverlap 3 to require an overlap of at least 3 bp.


### Trimming paired-end reads

The following command removes adapters from a paired-end reads, where the mate 1 and mate 2 reads are kept in files 'reads\_1.fq' and 'reads\_2.fq', respectively. The reads are trimmed for both Ns and low quality bases, and overlapping reads (at least 11 nucleotides, per default) are merged (collapsed):

    $ AdapterRemoval --file1 reads_1.fq --file2 reads_2.fq --basename output_paired --trimns --trimqualities --collapse

This command generates the files 'output_paired.pair1.truncated' and 'output_paired.pair2.truncated', which contain trimmed pairs of reads which were not collapsed, 'output_paired.singleton.truncated' containing reads where one mate was discarded, 'output_paired.collapsed' containing merged reads, and 'output_paired.collapsed.truncated' containing merged reads that have been trimmed due to the --trimns or --trimqualities options. Finally, the 'output_paired.discarded' and 'output_paired.settings' files correspond to those of the single-end run.


### Multiple input FASTQ files

More than one input file may be specified for mate 1 and mate 2 reads. This is accomplished simply by listing more than one file after the --file1 and the --file2 options.

For single-end reads:

    $ AdapterRemoval --file1 reads_1a.fq reads_1b.fq reads_1c.fq

And for paired-end reads:

    $ AdapterRemoval --file1 reads_1a.fq reads_1b.fq reads_1c.fq --file2 reads_2a.fq reads_2b.fq reads_2c.fq

AdapterRemoval will process these files as if they had been concatenated into a single file or pair of files prior to invoking AdapterRemoval. For paired reads, the files must be specified in the same order for --file1 and --file2.


### Interleaved FASTQ reads

AdapterRemoval is able to read and write paired-end reads stored in a single, so-called interleaved FASTQ file (one pair at a time, first mate 1, then mate 2). This is accomplished by specifying the location of the file using --file1 and *also* setting the --interleaved command-line option:

    $ AdapterRemoval --interleaved --file1 interleaved.fq --basename output_interleaved

Other than taking just a single input file, this mode operates almost exactly like paired end trimming (as described above); the mode differs only in that paired reads are not written to a 'pair1' and a 'pair2' file, but instead these are instead written to a single, interleaved file, named 'paired'. The location of this file is controlled using the --output1 option. Enabling either reading or writing of interleaved FASTQ files, both not both, can be accomplished by specifying the either of the --interleaved-input and --interleaved-output options, both of which are enabled by the --interleaved option.


### Combining FASTQ output

By default, AdapterRemoval will create one output file for each mate, one file for discarded reads, and (in PE mode) one file paired reads where one mate has been discarded, and (optionally) two files for collapsed reads. Alternatively, these files may be combined using the --combined-output, in which case all output is directed to the mate 1 and (in PE mode) to the mate 2 file. In cases where reads are discarded due to trimming to due to being collapsed into a single sequence, the sequence and quality scores of the discarded read is replaced with a single 'N' with base-quality 0. This option may be combined with --interleaved / --interleaved-output, to write a single, interleaved file in paired-end mode.


### Different quality score encodings

By default, AdapterRemoval expects the quality scores in FASTQ reads to be Phred+33 encoded, meaning that the error probabilities are encoded as (char)('!' - 10 * log10(p)). Most data will be encoded using Phred+33, but Phred+64 and 'Solexa' encoded quality scores are also supported. These are selected by specifying the --qualitybase command-line option (specifying either '33', '64', or 'solexa'):

    $ AdapterRemoval --qualitybase 64 --file1 reads_q64.fq --basename output_phred_64

By default, reads are written using the *same* encoding as the input. If a different encoding is desired, this may be accomplished using the --qualitybase-output option::

    $ AdapterRemoval --qualitybase 64 --qualitybase-output 33 --file1 reads_q64.fq --basename output_phred_33

Note furthermore that AdapterRemoval by default only expects quality scores in the range 0 - 41 (or -5 to 41 in the case of Solexa encoded scores). If input data using a different maximum quality score is to be processed, or if the desired maximum quality score of collapsed reads is greater than 41, then this limit may be increased using the --qualitymax option::

    $ AdapterRemoval --qualitymax 50 --file1 reads_1.fq --file2 reads_2.fq --collapse --basename output_collapsed_q50

For a detailed overview of Phred encoding schemes currently and previously in use, see e.g. the Wikipedia article on the subject:
https://en.wikipedia.org/wiki/FASTQ_format#Encoding


### Trimming paired-end reads with multiple adapter pairs

It is possible to trim data that contains multiple adapter pairs, by providing a one or two-column table containing possible adapter combinations (for single-end and paired-end trimming, respectively; see e.g. examples/adapters.txt):

    $ cat adapters.txt
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACCTAATCTCGTATGCCGTCTTCTGCTTG    AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
    AAACTTGCTCTGTGCCCGCTCCGTATGTCACAACAGTGCGTGTATCACCTCAATGCAGGACTCA    GATCGGGAGTAATTTGGAGGCAGTAGTTCGTCGAAACTCGGAGCGTCTTTAGCAGGAG
    CTAATTTGCCGTAGCGACGTACTTCAGCCTCCAGGAATTGGACCCTTACGCACACGCATTCATG    TACCGTGAAAGGTGCGCTTAGTGGCATATGCGTTAAGAGCTAGGTAACGGTCTGGAGG
    GTTCATACGACGACGACCAATGGCACACTTATCCGGTACTTGCGTTTCAATGCGCATGCCCCAT    TAAGAAACTCGGAGTTTGGCCTGCGAGGTAGCTTGGGTGTTATGAAGAACGGCATGCG
    CCATGCCCCGAAGATTCCTATACCCTTAAGGTCGCAATTGTTCGAGTAAGCTGTACGCGCCCAT    GTTGCATTGACCCGAAGGGCTCGATGTTTAGGGAGGTCAGAAGTTGAGCGGGTTCAAA

This table is then specified using the --adapter-list option:

    $ AdapterRemoval --file1 reads_1.fq --file2 reads_2.fq --basename output_multi --trimns --trimqualities --collapse --adapter-list adapters.txt

The resulting .summary file contains an overview of how frequently each adapter (pair) was used.

Note that in the case of paired-end adapters, AdapterRemoval considers only the combinations of adapters specified in the table, one combination per row. For single-end trimming, only the first column of the table file is required, and the list may therefore take the form of a file containing one sequence per line.


### Identifying adapter sequences from paired-ended reads

If we did not know the adapter sequences for the 'reads\_*.fq' files, AdapterRemoval may be used to generate a consensus adapter sequence based on fragments identified as belonging to the adapters through pairwise alignments of the reads, provided that the data set contains only a single adapter sequence (not counting differences in index sequences).

In the following example, the identified adapters corresponds to the default adapter sequences with a poly-A tail resulting from sequencing past the end of the insert + templates. It is not necessary to specify this tail when using the --adapter1 or --adapter2 command-line options. The characters shown under each of the consensus sequences represented the phred-encoded fraction of bases identical to the consensus base, with adapter 1 containing the index CACCTA:

    $ AdapterRemoval --identify-adapters --file1 reads_1.fq --file2 reads_2.fq

    Attemping to identify adapter sequences ...
    Processed a total of 1,000 reads in 0.0s; 129,000 reads per second on average ...
       Found 394 overlapping pairs ...
       Of which 119 contained adapter sequence(s) ...

    Printing adapter sequences, including poly-A tails:
      --adapter1:  AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
                   ||||||||||||||||||||||||||||||||||******||||||||||||||||||||||||
       Consensus:  AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACCTAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAAAA
         Quality:  55200522544444/4411330333330222222/1.1.1.1111100-00000///..+....--*-)),,+++++++**(('%%%$

        Top 5 most common 9-bp 5'-kmers:
                1: AGATCGGAA = 96.00% (96)
                2: AGATGGGAA =  1.00% (1)
                3: AGCTCGGAA =  1.00% (1)
                4: AGAGCGAAA =  1.00% (1)
                5: AGATCGGGA =  1.00% (1)


      --adapter2:  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
                   ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
       Consensus:  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
         Quality:  525555555144141441430333303.2/22-2/-1..11111110--00000///..+....--*-),,,+++++++**(%'%%%$

        Top 5 most common 9-bp 5'-kmers:
                1: AGATCGGAA = 100.00% (100)

No files are generated from running the adapter identification step.

The consensus sequences inferred are compared to those specified using the --adapter1 and --adapter2 command-line options, or with the default values for these if no values have been given (as in this case). Pipes (|) indicate matches between the provided sequences and the consensus sequence, and "*" indicate the presence of unspecified bases (Ns).


### Demultiplexing and adapter-trimming

As of version 2.1, AdapterRemoval supports simultaneous demultiplexing and adapter trimming; demultiplexing is carried out using a simple comparison between the specified barcode (a sequence of A, C, G, and T) and the first N bases of the mate 1 read, where N is the length of the barcode. Demultiplexing of double-indexed sequences is also supported, in which case two barcodes must be specified for each sample. The first barcode is then compared to first N_1 bases of the mate 1 read, and the second barcode is compared to the first N_2 bases of the mate 2 read. By default, this comparison requires a perfect match. Reads identified as containing a specific barcode(s) are then trimmed using adapter sequences including the barcode(s) as necessary. Reads for which no (pair of) barcodes matched are written to a separate file or pair of files (for paired end reads).

Demultiplexing is enabled by creating a table of barcodes, the first column of which species the sample name (using characters [a-zA-Z0-9_]) and the second and (optional) third columns specifies the barcode sequences expected at the 5' termini of mate 1 and mate 2 reads, respectively.

For example, a table of barcodes from a double-indexed run might be as follows (see examples/barcodes.txt):

    $ cat barcodes.txt
    sample_1 ATGCGGA TGAATCT
    sample_2 ATGGATT ATAGTGA
    sample_7 CAAAACT TCGCTGC

In the case of single-read reads, only the first two columns are required. AdapterRemoval is invoked with the --barcode-list option, specifying the path to this table:

    $ AdapterRemoval --file1 demux_1.fq --file2 demux_2.fq --basename output_demux --barcode-list barcodes.txt

This generates a set of output files for each sample specified in the barcode table, using the basename (--basename) as the prefix, followed by a dot and the sample name, followed by a dot and the default name for a given file type. For example, the output files for sample_2 would be

    output_demux.sample_2.discarded
    output_demux.sample_2.pair1.truncated
    output_demux.sample_2.pair2.truncated
    output_demux.sample_2.settings
    output_demux.sample_2.singleton.truncated

The settings files generated for each sample summarizes the reads for that sample only; in addition, a basename.settings file is generated which summarizes the number and proportion of reads identified as belonging to each sample.

The maximum number of mismatches allowed when comparing barocdes is controlled using the options --barcode-mm, --barcode-mm-r1, and --barcode-mm-r2, which specify the maximum number of mismatches total, and the maximum number of mismatches for the mate 1 and mate 2 barcodes respectively. Thus, if mm_1(i) and mm_2(i) represents the number of mismatches observed for barcode-pair i for a given pair of reads, these options require that

   1. mm_1(i) <= --barcode-mm-r1
   2. mm_2(i) <= --barcode-mm-r2
   3. mm_1(i) + mm_2(i) <= --barcode-mm


### Demultiplexing mode

As of version 2.2, AdapterRemoval can furthermore be used to demultiplex reads, without carrying out other forms of adapter trimming. This is accomplished by specifying the --demultiplex-only option:

    $ AdapterRemoval --file1 demux_1.fq --file2 demux_2.fq --basename output_only_demux --barcode-list barcodes.txt --demultiplex-only

Options listed under "TRIMMING SETTINGS" (see 'AdapterRemoval --help') do not apply to this mode, but compression (--gzip, --bzip2), multi-threading (--threads), interleaving (--interleaved, etc.) and other such options may be used in conjunction with --demultiplex-only.

AdapterRemoval will generate a '.settings' file for each sample listed in the --barcode-list file, along with the adapter-sequences that should be used when trimming reads for a given sample. These adapters correspond to the adapters that were specified when running AdapterRemoval in demultiplexing mode, with the barcode prefixed as appropriate. An underscore is used to demarcate the location at which the barcode ends and the adapter beings.

It is important to use these, updated, adapter sequences when trimming the demultiplexed reads, to avoid the inclusion of barcode sequences in reads extending past the 3' termini of the DNA template sequence.


## A note on specifying adapter sequences

Please note that the --pcr1 and --pcr2 options used with AdapterRemoval v1.x have been deprecated in favor of options --adapter1 and --adapter2. For both --adapter1 and --adapter2 the adapter sequence are expected to be observed in the raw mate 1 and mate 2 reads respectively, exactly as specified on the command-line, which corresponds to the behavior of most adapter trimming programs.

Default adapter #1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG

Default adapter #2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

Assuming these were the adapters used to generate our data, we should therefore see these in the FASTQ files (assuming that the read lengths are sufficiently long and that insert sizes are sufficently short), typically followed by a low-quality A-tail, when ignoring any difference in case and treating Ns as wildcards:

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

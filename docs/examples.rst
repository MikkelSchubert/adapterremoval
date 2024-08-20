.. highlight:: Bash


Example usage
=============

The following examples all make use of the data included in the 'examples' folder.


Basic usage
-----------

To trim single-end data, simplify specify one or more files using the --in-file1 option::

    adapterremoval3 --in-file1 reads_1.fastq --out-prefix output_single

To trim paired-end FASTQ reads, instead specify two input files using the ``--in-file1`` and ``--in-file2`` options::

    adapterremoval3 --in-file1 reads_1.fastq --in-file2 reads_2.fastq --out-prefix output_paired --merge

The ``--out-prefix`` option ensure that all output files start with the specified prefix. For the paired reads, the ``--merge`` option enables merging of reads overlapping at least 11bp (by default).

Reading from STDIN and Writing to STDOUT
----------------------------------------

Reading from STDIN and writing to STDOUT can be accomplished *either* by using the special ``/dev/stdin`` and ``/dev/stdout`` files, *or* by using the filename ``-`` for either ``--in-file*`` or ``--out-*`` options::

    some-command | adapterremoval3 --in-file1 - --out-file1 - | some-other-command

Input from STDIN and output to STDOUT can freely be interleaved as described in the *Interleaved input and output* below.

Multiple input FASTQ files
--------------------------

More than one input file may be specified listed after the ``--in-file1`` and ``--in-file2`` options. Files are processed in the specified order, as if they had been concatenated using ``cat`` or ``zcat``::

    adapterremoval3 --in-file1 reads_1a.fastq reads_1b.fastq reads_1c.fastq
    adapterremoval3 --in-file1 reads_1a.fastq reads_1b.fastq reads_1c.fastq --in-file2 reads_2a.fastq reads_2b.fastq reads_2c.fastq


Interleaved input and output
----------------------------

AdapterRemoval is able to read and write paired-end reads stored in a single, so-called interleaved FASTQ file (one pair at a time, first mate 1, then mate 2). This is accomplished by specifying the location of the file using ``--in-file1`` and *also* setting the ``--interleaved`` command-line option::

    adapterremoval3 --interleaved --in-file1 interleaved.fastq --out-prefix output_interleaved

Other than taking just a single input file, this mode operates almost exactly like paired end trimming (as described above); the mode differs only in that paired reads are not written to a 'r1' and a 'r2' file, but instead these are instead written to a single file. The location of this file is controlled using the ``--out-file1`` option.

Enabling either reading or writing of interleaved FASTQ files, both not both, can be accomplished by specifying the either of the ``--interleaved-input`` and ``--interleaved-output`` options, both of which are enabled by the ``--interleaved`` option.

Alternatively, you can simply specify the same file for both ``--in-file1`` and ``--in-file2``::

    adapterremoval3 --in-file1 interleaved.fastq --in-file2 interleaved.fastq --out-file1 output_interleaved.fastq.gz --out-file2 output_interleaved.fastq.gz

The ability to interleave input extends to all read types, and one could for example write discarded and singleton reads to the same file using the following command::

   adapterremoval3 --in-file1 interleaved.fastq --in-file2 interleaved.fastq --out-prefix output_interleaved --out-discarded output_interleaved.discarded.fastq.gz --out-singleton output_interleaved.discarded.fastq.gz

Different quality score encodings
---------------------------------

By default, AdapterRemoval expects the quality scores in FASTQ reads to be Phred+33 encoded, meaning that the error probabilities are encoded as ``(char)('!' - 10 * log10(p))``. Most data will be encoded using Phred+33, but Phred+64 and 'Solexa' encoded quality scores are also supported. These are selected by specifying the ``--quality-format`` command-line option (specifying either '33', '64', or 'solexa')::

    adapterremoval3 --quality-format 64 --in-file1 reads_q64.fastq --out-prefix output_phred_64

Output is always saved as Phred+33. See `this Wikipedia article`_ for a detailed overview of Phred encoding schemes currently and previously in use.


Trimming paired-end reads with multiple adapter pairs
-----------------------------------------------------

It is possible to trim data that contains multiple adapter pairs, by providing a one (SE) or two-column (PE) table containing possible adapter combinations::

    cat adapters.txt
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA    AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    AAACTTGCTCTGTGCCCGCTCCGTATGTCACAA    GATCGGGAGTAATTTGGAGGCAGTAGTTCGTCG
    CTAATTTGCCGTAGCGACGTACTTCAGCCTCCA    TACCGTGAAAGGTGCGCTTAGTGGCATATGCGT
    GTTCATACGACGACGACCAATGGCACACTTATC    TAAGAAACTCGGAGTTTGGCCTGCGAGGTAGCT
    CCATGCCCCGAAGATTCCTATACCCTTAAGGTC    GTTGCATTGACCCGAAGGGCTCGATGTTTAGGG

This table is then specified using the ``--adapter-list`` option::

    adapterremoval3 --in-file1 reads_1.fastq --in-file2 reads_2.fastq --out-prefix output_multi --merge --adapter-list adapters.txt

AdapterRemoval uses the exact pairs of adapters listed in the table. The resulting reports contains overviews of how frequently each adapter (pair) was used.


Identifying adapter sequences from paired-ended reads
-----------------------------------------------------

If we did not know the adapter sequences for the ``reads_*.fastq`` files, AdapterRemoval may be used to generate a consensus adapter sequence based on fragments identified as belonging to the adapters through pairwise alignments of the reads, provided that the data set contains only a single adapter sequence (not counting differences in index sequences).

In the following example, the identified adapters corresponds to the default adapter sequences with a poly-A tail resulting from sequencing past the end of the insert + templates. It is not necessary to specify this tail when using the ``--adapter1`` or ``--adapter2`` command-line options. The characters shown under each of the consensus sequences represented the Phred-encoded fraction of bases identical to the consensus base::

    adapterremoval3 --identify-adapters --in-file1 reads_1.fastq --in-file2 reads_2.fastq

    Attempting to identify adapter sequences ...
    Processed a total of 1,000 reads in 0.0s; 129,000 reads per second on average ...
       Found 394 overlapping pairs ...
       Of which 119 contained adapter sequence(s) ...

    Printing adapter sequences, including poly-A tails:
      --adapter1:  AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
                   |||||||||||||||||||||||||||||||||
       Consensus:  AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACCTAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAAAA
         Quality:  55200522544444/4411330333330222222/1.1.1.1111100-00000///..+....--*-)),,+++++++**(('%%%$

        Top 5 most common 9-bp 5'-kmers:
                1: AGATCGGAA = 96.00% (96)
                2: AGATGGGAA =  1.00% (1)
                3: AGCTCGGAA =  1.00% (1)
                4: AGAGCGAAA =  1.00% (1)
                5: AGATCGGGA =  1.00% (1)


      --adapter2:  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
                   |||||||||||||||||||||||||||||||||
       Consensus:  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
         Quality:  525555555144141441430333303.2/22-2/-1..11111110--00000///..+....--*-),,,+++++++**(%'%%%$

        Top 5 most common 9-bp 5'-kmers:
                1: AGATCGGAA = 100.00% (100)

No files are generated from running the adapter identification step.

The consensus sequences inferred are compared to those specified using the ``--adapter1`` and ``--adapter2`` command-line options, or with the default values for these if no values have been given (as in this case). Pipes (|) indicate matches between the provided sequences and the consensus sequence, and "*" indicate the presence of unspecified bases (Ns).

Best practice is to compare the consensus with published `Illumina`_ or `BGI/MGI`_ adapter sequences and pick out the best matches. However, on occasion there will be differences between the published sequences and the observed adapter sequences. When using the consensus directly, it is not recommended to use the full consensus sequence, since the quality declines quickly towards the 3'.


Demultiplexing
-----------------------------------

AdapterRemoval supports simultaneous demultiplexing and adapter trimming; demultiplexing is carried out using a simple comparison between the specified barcode (a sequence of A, C, G, and T) and the first N bases of the mate 1 read, where N is the length of the barcode. Demultiplexing of double-indexed sequences is also supported, in which case two barcodes must be specified for each sample. The first barcode is then compared to first ``N_1`` bases of the mate 1 read, and the second barcode is compared to the first ``N_2`` bases of the mate 2 read. By default, this comparison requires a perfect match. Reads identified as containing a specific barcode(s) are then trimmed using adapter sequences including the barcode(s) as necessary. Reads for which no (pair of) barcodes matched are written to a separate file or pair of files (for paired end reads).

Demultiplexing is enabled by creating a table of barcodes, the first column of which species the sample name (using characters a-z, A-Z, 0-9, or _) and the second and (optional) third columns specifies the barcode sequences expected at the 5' termini of mate 1 and mate 2 reads, respectively.

For example, a table of barcodes from a double-indexed run might be as follows (see examples/barcodes.txt)::

    cat barcodes.txt
    sample_1 ATGCGGA TGAATCT
    sample_2 ATGGATT ATAGTGA
    sample_7 CAAAACT TCGCTGC

AdapterRemoval is invoked with the ``--barcode-list`` option, specifying the path to this table::

    adapterremoval3 --in-file1 demux_1.fastq --in-file2 demux_2.fastq --out-prefix output_demux --barcode-list barcodes.txt

This generates a set of output files for each sample specified in the barcode table, using ``output_demux`` as the prefix for output filenames, followed by a dot and the sample name, followed by a dot and the default name for a given file type. The reports generated by AdapterRemoval contains information about the number of reads identified for each sample and (in the JSON file) detailed per-sample quality metrics.

The maximum number of mismatches allowed when comparing barcodes is controlled using the options ``--barcode-mm``, ``--barcode-mm-r1``, and ``--barcode-mm-r2``, which specify the maximum number of mismatches total, and the maximum number of mismatches for the mate 1 and mate 2 barcodes respectively. Thus, if mm_1(i) and mm_2(i) represents the number of mismatches observed for barcode-pair i for a given pair of reads, these options require that

   1. mm_1(i) <= ``--barcode-mm-r1``
   2. mm_2(i) <= ``--barcode-mm-r2``
   3. mm_1(i) + mm_2(i) <= ``--barcode-mm``


If the ``--demultiplex-only`` option is used, then no trimming/processing is performed after the demultiplexing step::

    adapterremoval3 --in-file1 demux_1.fastq --in-file2 demux_2.fastq --out-prefix output_only_demux --barcode-list barcodes.txt --demultiplex-only

These reads will still contain adapters, and for paired reads/double indexed data these adapters will be prefixed by the barcode sequence(s). The adapter plus barcode sequences are reported for each sample in the JSON report file.


.. _this Wikipedia article: https://en.wikipedia.org/wiki/FASTQ_format#Encoding

.. _Illumina: https://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
.. _BGI/MGI: https://en.mgitech.cn/Download/download_file/id/71

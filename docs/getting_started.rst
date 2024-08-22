.. highlight:: Bash

Getting started
===============

To run AdapterRemoval on single-end FASTQ data, simply specify the location of FASTQ file(s) using the ``--file1`` command-line options::

	AdapterRemoval --file1 myreads_1.fastq.gz

To run AdapterRemoval on paired-end FASTQ data, specify the location of the mate 1 and mate 2 FASTQ files using the ``--file1`` and ``--file2`` command-line options::

    AdapterRemoval --file1 myreads_1.fastq.gz --file2 myreads_2.fastq.gz

The files may be uncompressed, gzip-compressed, or bzip2 compressed. When run in this manner, AdapterRemoval will save the trimmed reads in the current working directly, using filenames starting with 'your_output'. This behavior may be changed using the ``--basename`` option, or using specific options for each output file.

More examples of common usage may be found in the :doc:`examples` section of the documentation.


A note on specifying adapters
-----------------------------

AdapterRemoval relies on the user specifying the adapter sequences to be trimmed, using the ``--adapter1`` and ``--adapter2`` command-line options. By default, AdapterRemoval is setup to trim Illumina Truseq adapters, corresponding to the following command-line options::

    --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
    --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

It is therefore extremely important to specify the correct adapter sequences when running AdapterRemoval on a dataset that does not make use of these adapters. Failure to do so will result in the wrong sequences being trimmed, and actual adapter sequences being left in the resulting "trimmed" reads.

Adapter sequences are specified in the read orientation when using the ``--adapter1`` and ``--adapter2`` command-line options, directly corresponding to the sequence that is observed in the FASTQ files produced by the base calling software. If we were processing data generated using the above TrueSeq adapters, then we would therefore expect to find those sequences as-is in our FASTQ files (assuming that the read lengths are sufficiently long and that insert sizes are sufficiently short), typically followed by a low-quality A-tail::

    $ grep "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC......ATCTCGTATGCCGTCTTCTGCTTG" file1.fq
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAACAAGAAT
    CTGGAGTTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGAATCTCGTATGCCGTCTTCTGCTTGAAAAAAA
    GGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGAATCTCGTATGCCGTCTTCTGCTTGCAAATTGAAAACAC

    $ grep "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" file2.fq
    CAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTCAAAAAAAGAAAAACATCTTG
    GAACTCCAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTCAAAAAAAATAGA
    GAACTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTCAAAAACATAAGACCTA

The ambiguous bases representing the mate 1 barcode (the six Ns) have been replaced by single-character wildcards (dots) in the above grep commands, corresponding to how AdapterRemoval itself treats such characters.

For paired-end data, the ``--identify-adapters`` mode may be used to verify the choice of adapters, by attempting to reconstruct the adapter sequence directly from the FASTQ reads. See the :doc:`examples` section for a demonstration of this functionality.

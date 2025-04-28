###################
 Detailed overview
###################

This page describes the data-flow in AdapterRemoval, and list what command-line options affect which parts of the process. For a detailed overview of input and output files, see :doc:`input_and_output`.

************************
 1. Reading FASTQ files
************************

FASTQ data is read from one or more files specified with ``--in-file1`` and (optionally) ``--in-file2``. If either the ``--interleaved`` or the ``--interleaved-input`` options are specified, the files specified with ``--in-file1`` are expected to contain sequential pairs of FASTQ reads.

FASTQ quality scores are expected to be Phred+33_ encoded by default, but this can be changed using the ``--quality-format`` option. By default, AdapterRemoval will attempt to infer the character used to separate the read name from the mate number (for example ``/`` in ``@my_read/1`` and ``@my_read/2``), but this can be overridden using the ``--mate-separator`` option.

If enabled by the ``--mask-degenerate-bases`` option, IUPAC encoded degenerate bases (``B``, ``D``, ``H``, ``K``, ``M``, ``R``, ``S``, ``V``, ``W``, and ``Y``) are replaced with ``N``. If enabled by ``--convert-uracils`` option, positions containing uracils (``U``) are converted to thymine (``T``). All bases are furthermore converted to uppercase and any dots (``.``) are replaced with ``N``.

Basic statistics are collected based on this pre-processed data, including detailed per-nucleotide statistics are collected from a random selection of reads determined by the ``--report-sample-rate`` option. If enabled by the ``--report-duplication`` option, the amount of read duplication is estimated using the FastQC_ methodology per mate.

If specified, the ``--head`` option limits the total number of reads (in single ended mode) or pairs (in paired ended mode) that are read from the input.

*******************
 2. Demultiplexing
*******************

If a table of barcodes is specified using the ``--barcode-list``, then demultiplexing of the data is performed prior to trimming. By default, each barcode or pair of barcodes is expected to identify a unique sample, but this restriction can be loosened using the ``--multiple-barcodes`` option.

By default, these barcodes are considered to be in a ``unspecified`` orientation. An orientation may be specified with ``--barcode-orientation``, and if either ``forward`` or ``reverse`` is used, then the corresponding reverse or forward barcode pairs are automatically generated.

The ``--barcode-mm`` option defines a global maximum number of mismatches allowed across both barcodes (in paired end mode), while the ``--barcode-mm-r1`` and ``--barcode-mm-r2`` options sets limits for either barcode in addition to the global maximum.

If the ``--demultiplex-only`` option is set, the read processing and filtering steps are skipped and output is written as described in `5. Output`_.

********************
 3. Read processing
********************

If the ``--pre-trim3p`` option is specified, the reads are trimmed by the specified amounts of bases. Next poly-X tails are trimmed if the ``--pre-trim-polyx`` option is specified.

The reads and adapters (single end mode) or reads + adapters (paired end mode) are aligned, using the adapter sequences specified with ``--adapter1`` and ``--adapter2``, with a maximum error rate defined by ``--mm``. The alignment is performed as follows, for single end and paired end reads expectively:

.. code-block:: text

   SE: [Read 1 -------------]
                   [Adapter 1 -----]

   PE: [Adapter 2' ----][Read 1 ------------]
                 [Read 2' ------------][Adapter 1 ----]

Note that adapter 2 and read 2 are both reverse complemented in the above. If multiple adapters are specified using ``--adapter-list``, then alignment is carried out with each adapter or adapter pair to find the best alignment.

If the best alignment covers at least ``--merge-threshold`` bases, then the sequences are considered to be overlapping. For paired end reads the inferred adapter sequences are always trimmed, while for single end reads the overlap must be at least ``--min-adapter-overlap`` bases before the sequence is trimmed (if specified). If the ``--merge`` option is set and the best alignment is at least ``--merge-threshold`` long, then the two reads are merged into a single sequence (paired end mode only).

If the ``--post-trim5p`` or ``--post-trim3p`` options are specified, the reads are trimmed by the specified amounts of bases and poly-X tails are trimmed if the ``--post-trim-polyx`` option is specified. Only trimming with the ``--post-trim5p`` is performed on merged reads, since the 3p ends are typically located inside the reads or correspond to the 5p of the other mate.

Quality trimming is performed using the trimming algorithm specified using the ``--quality-trimming`` option and the related options. By default both ends of a read will be quality trimmed, but trimming can be limited to the 3p end via the ``--preserve5p`` option.

If the ``--prefix-read1``, ``--prefix-read2``, and ``--prefix-merged`` options are set, the mate 1, mate 2, and merged read names are prefixed using the corresponding values.

*******************
 4. Read filtering
*******************

Reads with more than ``--max-ns`` ambiguous bases (Ns), that are shorter than ``--min-length``, that are longer than ``--max-length``, that have a mean Phred encoded quality score less than ``--min-mean-quality``, or that have a complexity score less than ``--min-complexity`` are discarded. If only one read in a pair is filtered, then the read that was not filtered is classified as a "singleton" read.

***********
 5. Output
***********

Reads are written to files specified using ``--out-prefix`` and/or individual ``--out-*`` options. If ``--out-prefix`` is used, AdapterRemoval will automatically name the default output files, excluding discarded reads, which are not saved by default.

Output files may also be named using the corresponding ``--out-*`` options, for example ``--out-discarded`` for discarded reads. If ``--out-prefix`` is combined with ``--out-*`` options, then those options take priority.

If multiple ``--out-*`` options point to the same file, then all those reads are written to the same file in order. At least one ``--out-*`` option must be specified, but one can perform a dry run by for example using ``--out-prefix /dev/null`` (see below).

Basic statistics are collected for all output data, and detailed per-nucleotide statistics are collected from a random selection of reads determined by the ``--report-sample-rate`` option.

Supported output formats
========================

Output is written as gzip compressed FASTQ records by default (see below), but this can be controlled using either the ``--out-format`` option or by manually specifying a filename with the desired extension using one of the ``--out-*`` options. AdapterRemoval currently recognizes ``.fq.gz``, ``.fastq.gz``, ``.fq``, ``fastq``, ``.sam.gz``, ``.sam``, and ``.bam``.

If AdapterRemoval does not recognize the extension (for example when writing to STDOUT), then the format specified with ``--out-format`` is used.

FASTQ output
------------

When not performing demultiplexing, AdapterRemoval will generate the following files when using ``--out-prefix``, replacing the key ``${prefix}`` with the value of ``--out-prefix``:

-  ``${prefix}.json``
-  ``${prefix}.html``
-  ``${prefix}.r1.fastq.gz``
-  ``${prefix}.r2.fastq.gz`` (paired end mode)
-  ``${prefix}.singleton.gz`` (paired end mode, if any filtering is enabled)
-  ``${prefix}.merged.gz`` (paired end mode, if read merging is enabled)

If demultiplexing is enabled, the following files will be generated with ``--out-prefix``, with the key ``${sample}`` additionally replaced with each sample name:

-  ``${prefix}.json``
-  ``${prefix}.html``
-  ``${prefix}.unidentified.r1.fastq.gz``
-  ``${prefix}.unidentified.r2.fastq.gz`` (paired end mode)
-  ``${prefix}.${sample}.r1.fastq.gz``
-  ``${prefix}.${sample}.r2.fastq.gz`` (paired end mode)
-  ``${prefix}.${sample}.singleton.gz`` (paired end mode, if any filtering is enabled)
-  ``${prefix}.${sample}.merged.gz`` (paired end mode, if read merging is enabled)

The ``unidentified`` files contain any reads for which no sample could be identified using the supplied barcode sequences. Note that no trimming, filtering, or merging is performed if ``--demultiplex-only`` is used.

SAM/BAM output
--------------

When using ``--out-format`` with ``.sam``, ``.sam.gz``, or ``.bam``, AdapterRemoval defaults to writing a single file containing *all* output, including discarded reads. For example, for BAM output:

-  ``${prefix}.json``
-  ``${prefix}.html``
-  ``${prefix}.bam``

Or when demultiplexing is enabled, replacing the key ``${sample}`` with each sample name:

-  ``${prefix}.json``
-  ``${prefix}.html``
-  ``${prefix}.unidentified.bam``
-  ``${prefix}.${sample}.bam``

This can be overridden as described above, using the individual ``--out-*`` options.

STDIN, STDOUT, and /dev/null
============================

AdapterRemoval supports four special filenames: ``-``, ``/dev/stdin``, ``/dev/stdout``, and ``/dev/null``:

-  If ``-`` or ``/dev/stdin`` is used for either ``--in-file1`` or ``--in-file2``, then AdapterRemoval will read those reads from STDIN. Note that it is currently not possible to use ``--in-file1 - --in-file2 -`` for interleaved input, instead use ``--interleaved-input --in-file1 -``.
-  If ``-`` or ``/dev/stdout`` is used for any ``--out-*`` options, then that file is written to STDOUT. As noted above, multiple read types may be written to STDOUT in this manner, meaning that ``--interleaved-output --out-file1 -`` and ``--out-file1 - --out-file2 -`` give identical output.
-  If ``/dev/null`` is used with any ``--out-*`` option, then the corresponding files are not saved. Statistics are still collected for any processed reads.

As specifying ``/dev/null`` skips the compression/write steps, it is highly recommended to disable any undesired file types in this way when using the ``--out-prefix`` option. Alternatively, only specify the exact file-types that you want to save for the same effect.

Interleaved output
==================

Interleaved output, in which both mate 1 and mate 2 reads are written to the same file, can be obtained in several ways:

-  By using the ``--interleaved`` option to enabled interleaved input and output. The interleaved reads are written to ``--out-file1`` filename.
-  By using the ``--interleaved-output`` to enable interleaved output. The interleaved reads are written to ``--out-file1`` filename.
-  By specifying the same filename for ``--out-file1`` and ``--out-file2``. This includes special filenames ``-`` and ``/dev/stdout``, to write interleaved output to STDOUT.

Compressed output
=================

As described above, AdapterRemoval defaults to writing gzip compressed FASTQ files. The compression level used for these, as well as for gzip compressed SAM files and for BAM files, and be controlled via the ``--compression-level`` option.

Depending on the compression level, this option either enables block-based compression using ``libdeflate`` (for levels 2 and above), or stream-based compression using ``isa-l`` for level 1. Compression level 0 is uncompressed, but includes gzip headers and checksums for validation.

.. note::

   Most tools are compatible with block-based compression, but if a downstream tool only processes a few hundreds of your trimmed reads, then try running AdapterRemoval with ``--compression-level 1`` to enable stream-based compression.

.. tip::

   On fast storage, and if space is not a limiting factor, then ``--compression-level`` can be decreased to 1 for an up to 50-100% increase in throughput, at the cost of 10-20% larger files, depending on the data and hardware.

.. _fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

.. _phred+33: https://en.wikipedia.org/wiki/FASTQ_format#Quality

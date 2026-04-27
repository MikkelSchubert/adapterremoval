###################
 Detailed overview
###################

This page describes the data-flow in AdapterRemoval, and lists what command-line options affect which parts of the process.

************************
 1. Reading FASTQ files
************************

FASTQ data is read from one or more files specified with ``--in-file1`` and (optionally) ``--in-file2``. If either the ``--interleaved`` or the ``--interleaved-input`` options are specified, the file or files specified with ``--in-file1`` are expected to contain sequential pairs of FASTQ reads. Multiple files may be specified for ``--in-file1`` and ``--in-file2``, in which these sets of files are treated as if they were concatenated into one big file. See below for how to read data from STDIN.

FASTQ quality scores are expected to be Phred+33_ encoded by default, but this can be changed using the ``--quality-format`` option. By default, AdapterRemoval will attempt to infer the character used to separate the read name from the mate number (for example ``/`` in ``@my_read/1`` and ``@my_read/2``), but this can be overridden using the ``--mate-separator`` option.

If enabled by the ``--mask-degenerate-bases`` option, IUPAC encoded degenerate bases (``B``, ``D``, ``H``, ``K``, ``M``, ``R``, ``S``, ``V``, ``W``, and ``Y``) are replaced with ``N``. If enabled by the ``--convert-uracils`` option, positions containing uracils (``U``) are converted to thymine (``T``). All bases are furthermore converted to uppercase and any dots (``.``) are replaced with ``N``.

Basic statistics are collected based on this pre-processed data, including detailed per-nucleotide statistics, using a fraction of reads determined by the ``--report-sample-rate`` option. This sampling of data, for rates less than 1, represents the only random behavior in AdapterRemoval. If enabled by the ``--report-duplication`` option, the amount of read duplication is estimated using the FastQC_ methodology, individually for mate 1 and mate 2 reads.

If no adapters have been specified using ``--adapter1`` and ``--adapter2``, and ``--adapter-selection`` is set to ``auto``, AdapterRemoval will attempt to identify the adapters present in the first 100,000 reads or read pairs, using a built-in database of adapter sequences (see ``--adapter-database``). If no adapters could be identified in this manner, AdapterRemoval will fall back to the behavior specified via ``--adapter-fallback``, either aborting, assuming that the reads contain no adapter sequences, or attempting to infer adapters purely from overlap analyses (PE only).

If specified, the ``--head`` option limits the total number of reads (in single-end mode) or pairs (in paired-end mode) that are read from the input.

For a detailed description of supported input files, see :doc:`input_and_output`.

*******************
 2. Demultiplexing
*******************

If a table of barcodes is specified using the ``--barcode-table``, then demultiplexing of the data is performed prior to trimming of adapters and low-quality bases. By default, each barcode or pair of barcodes is expected to identify a unique sample, but this restriction can be loosened using the ``--multiple-barcodes`` option.

By default, these barcodes are considered to be in a ``unspecified`` orientation. An orientation may be specified with ``--barcode-orientation``, and if either ``forward`` or ``reverse`` is used, then the corresponding reverse or forward barcode pairs are automatically generated. When combined with ``--normalize-orientation``, this normalizes the orientation of merged reads (see below).

The ``--barcode-mm`` option defines a global maximum number of mismatches allowed across both barcodes (in paired-end mode), while the ``--barcode-mm-r1`` and ``--barcode-mm-r2`` options sets limits for either barcode in addition to the global maximum.

If the ``--demultiplex-only`` option is set, the read processing and filtering steps are skipped and output is written as described in `5. Output`_.

********************
 3. Read processing
********************

If the ``--pre-trim3p`` option is enabled, the reads are trimmed by the specified amounts of bases. Next, poly-X tails are trimmed if the ``--pre-trim-polyx`` option is enabled.

Each read and adapter sequence (single-end mode) or each read pair and adapter sequences (paired-end mode) are aligned, using the adapter sequences specified with ``--adapter1`` and ``--adapter2`` or selected via ``--adapter-selection`` (see above). The alignment is performed as follows, for single-end and paired-end reads respectively:

.. code-block:: text

    SE: [Read 1 -------------]
                    [Adapter 1 -----]

    PE: [Adapter 2' ----][Read 1 ------------]
                  [Read 2' ------------][Adapter 1 ----]

Note that adapter 2 and read 2 are both reverse complemented in the above figure. If multiple adapters are specified using ``--adapter-table``, then alignment is carried out with each adapter or adapter pair, in order to select the single best alignment, if any.

If a valid alignment is found, and if the error rate in the alignment does not exceed the maximum error rate defined by ``--mismatch-rate``, and if the alignment involves at least ``--min-overlap`` unambiguous bases (not 'N'), then the inferred adapter sequences are trimmed.

In paired-end mode only, if the ``--merge`` option is set and the best alignment is at least ``--merge-threshold`` long, then the two reads are merged into a single sequence. Ambiguous nucleotides (Ns) are not counted as part of the alignment, for the purpose of checking ``--merge-threshold``.

Depending on the ``--merge-strategy`` option, Phred scores of matching merged bases are either set to the maximum of the two observed Phred-scores, or to the sum of Phred scores (the probability that both reads are wrong). The output Phred score is capped by the ``--merge-quality-max``, to ensure that quality scores are compatible with downstream tools.

If the ``--post-trim5p`` or ``--post-trim3p`` options are specified, the reads are trimmed by the specified amounts of bases and poly-X tails are trimmed if the ``--post-trim-polyx`` option is specified. Only trimming with the ``--post-trim5p`` is performed on merged reads, since the 3' ends are typically located inside the reads or correspond to the 5' of the other mate (see alignment above).

Quality trimming is performed using the trimming algorithm specified using the ``--quality-trimming`` option and the related options. By default both ends of a read will be quality trimmed, but trimming can be limited to the 3' end via the ``--preserve5p`` option.

If the ``--prefix-read1``, ``--prefix-read2``, and ``--prefix-merged`` options are set, the mate 1, mate 2, and merged read names are prefixed using the corresponding values.

*******************
 4. Read filtering
*******************

Trimmed reads with more than ``--max-ns`` or ``--max-ns-fraction`` ambiguous bases (Ns), that are shorter than ``--min-length`` or longer than ``--max-length``, that have a mean Phred encoded quality score less than ``--min-mean-quality``, or that have a complexity score less than ``--min-complexity`` are discarded. If only one read in a pair is filtered, then the read that was not filtered is classified as a "singleton" read below.

***********
 5. Output
***********

Reads are written to files specified using any combination of ``--out-prefix`` and/or individual ``--out`` options, such as ``--out-discarded``. The ``--out-prefix`` automatically assigns filenames for all output using the specified prefix, excluding discarded reads, which are not saved by default. If ``--out-prefix`` is combined with individual ``--out`` options, then the individual options take priority.

If multiple ``--out`` options point to the same file, then all those reads are written to the same file in input-order. At least one ``--out`` option must be specified, but one can perform a dry run by setting any output option to ``/dev/null`` (see below).

Basic statistics are collected for all output data, and detailed per-nucleotide statistics are collected from a random selection of reads determined by the ``--report-sample-rate`` option.

See below for how to read data from STDIN. For a detailed overview output files, see :doc:`input_and_output`.

Supported output formats
========================

Output is written as gzip compressed FASTQ records by default (see below), but this can be controlled using either the ``--out-format`` option for filenames generated via the ``--out-prefix`` option, or by manually specifying a filename with the desired extension using one of the ``--out`` options. AdapterRemoval currently recognizes ``.fq.gz``, ``.fastq.gz``, ``.fq``, ``.fastq``, ``.sam.gz``, ``.sam``, and ``.bam``. If AdapterRemoval does not recognize the extension, then the format specified with ``--out-format`` is used.

The ``--stdout-format`` sets the format for reads written to STDOUT, which defaults to the same format as ``--out-format``, but uncompressed.

Interleaved output
==================

Interleaved output, in which both mate 1 and mate 2 reads are written to the same file, can be produced in several ways:

- By using the ``--interleaved`` or ``--interleaved-output`` option to enable interleaved input and output, or just interleaved output. The interleaved reads are written to ``--out-file1`` filename.
- By specifying the same filename for ``--out-file1`` and ``--out-file2``. This includes special filenames ``-`` and ``/dev/stdout``, to write interleaved output to STDOUT (see below).

Compressed output
=================

As described above, AdapterRemoval defaults to writing gzip compressed FASTQ files. The compression level used for these, as well as for gzip compressed SAM files and for BAM files, can be controlled via the ``--compression-level`` option.

Depending on the compression level, this option either enables block-based compression using ``libdeflate`` (for levels 2 and above), or stream-based compression using ``isa-l`` for level 1. Compression level 0 is uncompressed, but includes gzip headers and checksums for validation.

.. note::

    Most tools are compatible with block-based compression, but if a downstream tool only processes a few hundreds of your trimmed reads, then try running AdapterRemoval with ``--compression-level 1`` to enable stream-based compression. BAM output is always block-compressed.

.. tip::

    On fast storage, and if space is not a limiting factor, then ``--compression-level`` can be decreased to 1 for an up to 50-100% increase in throughput, at the cost of 10-20% larger files, depending on the data and hardware.

******************************
 STDIN, STDOUT, and /dev/null
******************************

AdapterRemoval supports four special filenames: ``-``, ``/dev/stdin``, ``/dev/stdout``, and ``/dev/null``:

- If ``-`` or ``/dev/stdin`` is used for either ``--in-file1`` or ``--in-file2``, then AdapterRemoval will read those reads from STDIN. Note that it is currently not possible to use ``--in-file1 - --in-file2 -`` for interleaved input, instead use ``--interleaved-input --in-file1 -``.
- If ``-`` or ``/dev/stdout`` is used for any ``--out-*`` options, then that file is written to STDOUT. As noted above, multiple read types may be written to STDOUT in this manner, meaning that ``--interleaved-output --out-file1 -`` and ``--out-file1 - --out-file2 -`` give identical output.
- If ``/dev/null`` is used with any ``--out-*`` option, then the corresponding files are not saved. Statistics are still collected for any processed reads.

As specifying ``/dev/null`` skips the compression/write steps, it is highly recommended to disable any undesired file types in this way when using the ``--out-prefix`` option. Alternatively, only specify the exact file-types that you want to save for the same effect.

.. _fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

.. _phred+33: https://en.wikipedia.org/wiki/FASTQ_format#Quality

Detailed overview
=================

The following page describes the data-flow in AdapterRemoval3, and list command-line options affect which parts of the process. For a detailed overview of input and output files, see :doc:`input_and_output`.

1. Reading FASTQ files
----------------------

FASTQ data is read from one or more files specified with `--in-file1` and (optionally) `--in-file2`. If either the `--interleaved` or the `--interleaved-input` options are specified, the files specified with `--in-file1` are expected to contain sequential pairs of FASTQ reads.

The quality is expected to be `Phred+33`_ by default, but this can be changed using the `--quality-format` option. By default AdapterRemoval3 will attempt to infer the mate-number separator (the character before the 1/2 mate number in read names), but this can be overridden using the `--mate-separator` option.

Basic statistics are collected for all input data, and detailed per-nucleotide statistics are collected from a random selection of reads determined by the `--report-sample-rate` option. If enabled by the `--report-duplication` option, the amount of read duplication is estimated using the `FastQC`_ methodology.

If specified, the `--head` option limits the total number of reads (in single ended mode) or pairs (in paired ended mode) that are read from the input.

2. Demultiplexing
-----------------

If a table of barcodes is specified using the `--barcode-list`, then demultiplexing of the data is performed prior to trimming. The `--barcode-mm` option defines a global maximum number of mismatches allowed across both barcodes (in paired end mode), while the `--barcode-mm-r1` and `--barcode-mm-r2` options sets limits for either barcode in addition to the global maximum.

If the `--demultiplex-only` option is set, the trimming steps are skipped and output is written as described in `5. Output`_.

3. Read processing
------------------

If the `--pre-trim3p` option is specified, the reads are trimmed by the specified amounts of bases. Next poly-X tails are trimmed if the `--pre-trim-polyx` option is specified.

The reads and adapters (single end mode) or reads + adapters (paired end mode) are aligned, using the adapter sequences specified with `--adapter1` and `--adapter2`, with a maximum error rate defined by `-mm`. The alignment is performed as follows, for single end and paired end reads expectively:

.. code-block:: text

    SE: [Read 1 -------------]
                    [Adapter 1 -----]

    PE: [Adapter 2' ----][Read 1 ------------]
                  [Read 2' ------------][Adapter 1 ----]

Note that adapter 2 and read 2 are both reverse complemented in the above. If multiple adapters are specified using `--adapter-list`, then alignment is carried out with each adapter or adapter pair to find the best alignment.

If the best alignment covers at least `--merge-threshold` bases, then the sequences are considered to be overlapping. For paired end reads the inferred adapter sequences are always trimmed, while for single end reads the overlap must be at least `--min-adapter-overlap` bases before the sequence is trimmed (if specified). If the `--merge` option is set and the best alignment is at least `--merge-threshold` long, then the two reads are merged into a single sequence (paired end mode only).

If the `--post-trim5p` or `--post-trim3p` options are specified, the reads are trimmed by the specified amounts of bases and poly-X tails are trimmed if the `--post-trim-polyx` option is specified. Only trimming with the `--post-trim5p` is performed on merged reads, since the 3p ends are typically located inside the reads or correspond to the 5p of the other mate.

Quality trimming is performed using the trimming algorithm specified using the `--trim-strategy` option and the related options. By default both ends of a read will be quality trimmed, but trimming can be limited to the 3p end via the `--preserve5p` option.

If the `--prefix-read1`, `--prefix-read2`, and `--prefix-merged` options are set, the mate 1, mate 2, and merged read names are prefixed using the corresponding values.

1. Read filtering
-----------------

Reads with more than `--max-ns` ambiguous bases (Ns), that are shorter than `--min-length`, that are longer than `--max-length`, or that have a complexity score less than `--min-complexity` are discarded. If one read in a pair is filtered, then the remaining read is classified as a "singleton" read.

5. Output
---------

Reads are written to the files specified using `--out-prefix` and/or the individual `--out-*` options, with the `--out-*` options taking priority. If multiple `--out-*` options point to the same file, then all of those reads are written to that file. Interleaved output can also be obtained by using the `--interleaved` or the `--interleaved-output` options, in which case mate 1 and mate 2 reads are both written to the file specified with `--in-file1` or `{prefix}.fastq` if `--out-prefix` is used.

If demultiplexing is enabled, the reads that did not match any barcodes will be written to `{prefix}.unidentified.fastq` in single-end mode and to `{prefix}.unidentified.r1.fastq` and `{prefix}.unidentified.r2.fastq` in paired end mode.

Gzip compression may be enabled using either the `--gzip` option for all files, or per file using the `--out-*` options by specifying a filename ending with a `.gz` extension. The `--gzip-level` option specifies the compression level used, with block-based compression using `libdeflate` used for compression levels above 3 and `isa-l` compression used for levels 3 and below.

The special filename `/dev/null` may be used with output options to signal that the corresponding reads should not be saved. Statistics are still collected for these reads, but compression/writing is not performed. Either `--out-prefix` or at least one `--out-*` option must be specified, but one can perform a dry run by for example using `--prefix /dev/null`.

Basic statistics are collected for all output data, and detailed per-nucleotide statistics are collected from a random selection of reads determined by the `--report-sample-rate` option.

.. _phred+33: https://en.wikipedia.org/wiki/FASTQ_format#Quality
.. _fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

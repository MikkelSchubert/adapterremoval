AdapterRemoval manpage
======================

Synopsis
--------

**adapterremoval3** [*options*...] --in-file1 <*filenames*> [--in-file2 <*filenames*>]


Description
-----------

:program:`adapterremoval3` removes residual adapter sequences from single-end (SE) or paired-end (PE) FASTQ reads, optionally trimming Ns and low qualities bases and/or merging overlapping paired-end mates into one read. Low quality reads are filtered based on the resulting length and the number of ambiguous nucleotides ('N') present following trimming. These operations may be combined with simultaneous demultiplexing using 5' barcode sequences. Reports containing statistics and plots in HTML and JSON format are generated after each run.

Alternatively, ``adapterremoval3`` may perform a read-only analysis of the input data, including the reconstruction of a consensus adapter sequences from paired-end data.

If you use this program, please cite the paper:

    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 12;9(1):88

    http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2


For detailed documentation, please see

    http://adapterremoval.readthedocs.io/


Options
-------

.. program:: adapterremoval3

.. option:: -h, --help

    Display summary of command-line options.

.. option:: -v, --version

    Print the version string.

.. option:: --licenses
    Print licenses for this software.

.. option:: --threads n

    Maximum number of threads. Defaults to 1.

.. option:: --simd name

    Selects the preferred SIMD instruction. Possible values are none, SSE2, AVX2, AVX512, and NEON. Options may be unavailable depending on the current system and depending on the compiler used to build AdapterRemoval. By default, AdapterRemoval will use the most advanced instruction available.


Input files
~~~~~~~~~~~

.. option:: --in-file1 filename [filenames...]
.. option:: --in-file2 filename [filenames...]

    Read FASTQ mate 1 and mate 2 reads from one or more files, either uncompressed or gzip compressed.  The ``--interleaved`` and ``--interleaved-input`` options may be used to enable reading of interleaved reads from the files specified with ``--in-file1``, in which case ``--in-file2`` is not used.

.. option:: --head n

    If set, AdapterRemoval will process only the first N reads/pairs of reads in the input. Accepts suffixes K (thousands), M (millions), and G (billions). By default, all data is processed.


Output file options
~~~~~~~~~~~~~~~~~~~

Output files for AdapterRemoval may be specified either via the ``--out-prefix`` option, which assigns default filenames, and/or via the individual ``--out-*`` to set or override the output filename for a given output read type. The same filename may be used for multiple ``--out-*`` options, in which case all of those reads are written to that file in input order.

Use ``-`` or ``/dev/stdin`` to read from standard input, use ``-`` or ``/dev/stdout`` to write to standard output, and use ``/dev/null`` to disable the generation of those types of reads. Statistics are still collected for read types written to ``/dev/null``, but the data itself is not serialized or compressed, in order to save time.

.. option:: --out-prefix path

    Prefix for the output files for which no filename was set using the corresponding options below, except for ``--out-discarded`` which is not saved by default. If this option is not used, then files for which no ``--out`` option was set will not be saved.

.. option:: --out-file1 filename
.. option:: --out-file2 filename

    Output files containing trimmed mate 1 reads and mate 2 reads. If interleaved output is enabled via the ``--interleaved`` or the ``--interleaved-output`` options, then only ``--out-file1`` is used, and this file will contain both mate 1 and mate 2 reads.

.. option:: --out-merged filename

    When used with ``--merged``, this file contains overlapping mate-pairs which have been merged into a single read. Setting this option in demultiplexing mode overrides ``--out-prefix`` for just this read type.

.. option:: --out-singleton filename

    Output file to which containing paired reads for which the mate has been discarded. Setting this option in demultiplexing mode overrides ``--out-prefix`` for just this read type.

.. option:: --out-unidentified1 filename
.. option:: --out-unidentified2 filename

    In demultiplexing mode, reads that could not be assigned to a single sample are written to these files. In interleaved mode, both mate 1 and mate 2 reads are written to ``--out-unidentified1``.

.. option:: --out-discarded filename

    Contains reads discarded due to the ``--min-length``, ``--max-length``, or ``--max-ns``, ``--min-mean-quality``, ``--min-complexity`` options. Setting this option in demultiplexing mode overrides ``--out-prefix`` for just this read type, but this option is not enabled by just setting ``--out-prefix``. This file will furthermore not be written if all filtering is disabled.

.. option:: --out-json filename
.. option:: --out-html filename

    Reports in JSON/HTML format containing information about the parameters used in the run, overall statistics on the reads before and after trimming, demultiplexing statistics, and the results of any analyses carried out during the runs. Analyses include insert inference of sizes, duplication levels, and any inferred consensus adapter sequences.

FASTQ options
~~~~~~~~~~~~~
.. option:: --quality-format name

    The Phred quality scores encoding used in input reads - either ``64`` for Phred+64 (Illumina 1.3+ and 1.5+) or ``33`` for Phred+33 (Illumina 1.8+, BGISeq, and more). In addition, the value ``solexa`` may be used to specify reads with Solexa encoded scores. The ``sam`` format may be used for Phred+33 encoded data with quality scores higher than that normally produced by squencing machines. Defaults to ``33``.

.. option:: --mate-separator separator

    Character separating the mate number (``1`` or ``2``) from the read name in FASTQ records, such as ``@my-read/1`` and ``@my-read/2``. This is typically either ``/`` or ``.``. By default, AdapterRemoval will attempt to infer this from the data.

.. option:: --interleaved-input

    Enable reading of interleaved FASTQ reads from the files specified with ``--in-file1``. Defaults to off.

.. option:: --interleaved-output

    Write paired-end reads to the file specified by ``--out-file1``, interleaving mate 1 and mate 2 reads. Defaults to off.

.. option:: --interleaved

    Enables ``--interleaved-input`` and ``--interleaved-output``. Defaults to off.

.. option:: --mask-degenerate-bases

    Mask degenerate/ambiguous IUPAC encoded bases (B/D/H/K/M/N/R/S/V/W/Y) in the input by replacing them with an 'N'; if this option is not used, AdapterRemoval will abort upon encountering degenerate bases in the input.

.. option:: --convert-uracils

    Convert uracils (U) to thymine (T) in input reads; if this option is not used, AdapterRemoval will abort upon encountering uracils in the input.


Output compression options
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: --out-format name

    Selects the default output format for files, one of either ``fastq`` for uncompressed FASTQ reads, ``fastq.gz`` for gzip compressed FASTQ reads, ``sam`` for uncompressed SAM records, ``sam.gz`` for gzip compressed SAM records, ``bam`` for BGZF compressed BAM records, and ``ubam`` for uncompressed BAM records. Setting a ``--out-*`` option overrides this option based on the extension used (except ``.ubam``). Defaults to ``fastq.gz``.

.. option:: --stdout-format name
    Selects the output format for data written to STDOUT, i.e. when writing to ``-`` or ``/dev/stdout``; choices are the same as for --out-format. By default, the uncompressed version of the current ``--out-format`` is used.

.. option:: --read-group [n, ...]

    Add read-group (RG) information to SAM/BAM output. Tags can either be specified as individual arguments, i.e. ``--read-group ID:foo SM:bar``, or as a string containing tags separated by tabs. Normally the latter can be written as ``--read-group $'ID:DS-1\tSM:TK-421\tPL:ILLUMINA'`` on the command-line. If the ID tag is not provided, the default ID ``1`` will be used.

.. option:: --compression-level n

    Sets the compression level for compressed output. Valid values are 0 to 13: Level 0 is uncompressed but includes gzip headers/checksums, level 1 is streamed for FASTQ output, which may be required in rare cases for compatibility. FASTQ output with compression levels 2 to 13, SAM, and BAM output is block compressed using the gzip compatible BGZF format. Defaults to compression level 5, roughly equivalent to ``gzip -4``.


FASTQ processing options
~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: --adapter1 adapter

    Adapter sequence expected to be found in mate 1 reads, specified in read direction. For a detailed description of how to provide the appropriate adapter sequences, see the "Adapters" section of the online documentation. Default is ``AGATCGGAAGAGCACACGTCTGAACTCCAGTCA``, intended for Illumina TruSeq and similar data.

.. option:: --adapter2 adapter

    Adapter sequence expected to be found in mate 2 reads, specified in read direction. For a detailed description of how to provide the appropriate adapter sequences, see the "Adapters" section of the online documentation. Default is ``AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT``, intended for Illumina TruSeq and similar data.

.. option:: --adapter-list filename

    Read one or more adapter sequences from a table. The first two columns (separated by whitespace) of each line in the file are expected to correspond to values passed to --adapter1 and --adapter2. In single-end mode, only column one is required. Lines starting with '#' are ignored. When multiple rows are found in the table, AdapterRemoval will try each adapter (pair), and select the best aligning adapters for each FASTQ read processed.

.. option:: --min-adapter-overlap length

    In single-end mode, reads are only trimmed if the overlap between read and the adapter is at least X bases long, not counting ambiguous nucleotides (N). Defaults to 1.

.. option:: --mismatch-rate rate

    The allowed fraction of mismatches allowed in the aligned region. If the value is less than 1, then the value is used directly. If ``--mismatch-rate`` is greater than 1, the rate is set to 1 / ``rate``. The default setting is 6 when trimming adapters, corresponding to a maximum mismatch rate of 1/6, and 10 when using ``--identify-adapters``.

.. option:: --shift n

    To allow for missing bases in the 5' end of the read, the program can let the alignment slip ``--shift`` bases in the 5' end. This corresponds to starting the alignment maximum ``--shift`` nucleotides into read 2 (for paired-end) or the adapter (for single-end). The default is 2.

.. option:: --merge

    In paired-end mode, merge overlapping mates into a single and recalculate the quality scores. The overlap needs to be at least ``--merge-threshold`` nucleotides, with a maximum number of mismatches determined by ``--mismatch-mate``. This option has no effect in single-end mode.

.. option:: --merge-threshold length

    The minimum overlap between mate 1 and mate 2 before the reads are merged into one, potentially longer sequence, when merging paired-end reads. Default is 11 bases.

.. option:: --merge-strategy name

    Determines how to assign quality scores to matching/mismatching bases during read merging. The ``maximum`` strategy uses ``Q=max(Q1, Q2)`` for matches while the ``additive`` strategy uses ``Q = Q1 + Q2``. Both strategies use ``Q = abs(Q1 - Q2)`` for mismatches, and picks the highest quality base, unless the qualities are the same in which case ``N`` is used. Setting this option implies --merge. Possible values are maximum, and additive. Defaults to ``maximum``.

.. option:: --merge-quality-max phred

    Sets the maximum Phred score for re-calculated quality scores when read merging is enabled with the 'additive' merging strategy. The value must be in the range 0 to 93, corresponding to Phred+33 encoded values of ``!`` to ``~``. Defaults to 41.

.. option:: --prefix-read1 X

    Adds the specified prefix to the names of mate 1 reads. Default to no prefix.

.. option:: --prefix-read2 X

    Adds the specified prefix to the names of mate 2 reads. Default to no prefix.

.. option:: --prefix-merged X

    Adds the specified prefix to merged read names. Default to no prefix.


Quality trimming options
~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: --pre-trim3p n [n]

    Trim the 3' of reads by a fixed amount after demultiplexing but before removing adapters. Specify one value to trim mate 1 and mate 2 reads the same amount, or two values separated by a space to trim each mate a different amount. Off by default.

.. option:: --post-trim5p n [n]

    Trim the 5' of reads by a fixed amount after removing adapters, but before carrying out quality based trimming. See ``--pre-trim3p``.

.. option:: --post-trim3p n [n]

    Trim the 3' of reads by a fixed amount after removing adapters, but before carrying out quality based trimming. See ``--pre-trim3p``.

.. option:: --quality-trimming method

    The method used for performing quality trimming; ``none`` to disable quality trimming, ``mott`` to enable trimming using the modified Mott's algorithm, ``window`` to perform window based quality trimming, and ``per-base`` to perform base-by-base trimming of low-quality bases and Ns (if enabled). Defaults to Mott's algorithm.

.. option:: --trim-mott-quality phred

    The threshold Phred score used when performing trimming quality based trimming using the modified Mott's algorithm. The value must be in the range 0 to 93, corresponding to Phred+33 encoded values of ``!`` to ``~``. Default to 13, which is an error rate of roughly 0.05.

.. option:: --trim-windows size

    Trim low quality bases using a sliding window based approach inspired by :program:`sickle` with the given window size. See the "Window based quality trimming" section of the manual page for a description of this algorithm. Applies when ``window`` based trimming is enabled using ``--quality-trimming``. Defaults to 0.1.

.. option:: --trim-min-quality minimum

    Inclusive minimum quality used when trimming low-quality bases with --trimming-strategy 'window' and 'per-base'. Applies when ``window`` based or ``per-base`` trimming is enabled using ``--quality-trimming``. Defaults to 2.

.. option:: --trim-qualities

    Trim consecutive stretches of low quality bases (threshold set by ``--trim-min-quality``) from the 5' and 3' termini. If trimming of Ns is also enabled (``--trim-ns``), then stretches of mixed low-quality bases and Ns are trimmed. Applies when ``per-base`` trimming is enabled using ``--quality-trimming``.

.. option:: --trim-ns

    Trim consecutive Ns from the 5' and 3' termini. If quality trimming is also enabled (``--trim-qualities``), then stretches of mixed low-quality bases and/or Ns are trimmed. Applies when ``window`` based or ``per-base`` trimming is enabled using ``--quality-trimming``.

.. option:: --pre-trim-polyx [nucleotides...]

    Enable trimming of poly-X tails prior to read alignment and adapter trimming. Zero or more nucleotides (A, C, G, T) may be specified, separated by spaces, with zero nucleotides corresponding to all of A, C, G, and T. Defaults to no trimming.

.. option:: --post-trim-polyx [nucleotides...]

    Enable trimming of poly-X tails after read alignment and adapter trimming/merging, but before trimming of low-quality bases. Merged reads are not trimmed by this option, since both ends are derived from the 5' of reads. See ``--pre-trim-polyx``. Off by default.

.. option:: --preserve5p

    If set, bases at the 5p will not be trimmed by ``mott``, ``window``, or ``per-base`` trimming, except if the entire read consists of low-quality bases. Merged reads will not be quality trimmed when this option is enabled due to the 3' ends being located inside the reads or overlapping the 5' of the source sequences.


Filtering options
~~~~~~~~~~~~~~~~~

.. option:: --max-ns n

    Discard reads containing more than ``--max`` ambiguous bases ('N') after trimming. Default is no maximum.


.. option:: --min-length length

    Reads shorter than this length are discarded following trimming. Defaults to 15.

.. option:: --max-length length

    Reads longer than this length are discarded following trimming. Defaults to no maximum.

.. option:: --min-mean-quality X

    Reads with a mean quality score less than this value (in the range 0 to 126) following trimming are discarded. Defaults to no minimum.

.. option:: min-complexity X

    Reads with a sequence quality less than this value after trimming are discarded. Complexity is measured as the fraction of positions that differ from the previous position. A suggested value is 0.3. Defaults to no minimum.


Demultiplexing options
~~~~~~~~~~~~~~~~~~~~~~

.. option:: --barcode-list filename

    Perform demultiplexing using table of one or two fixed-length barcodes for SE or PE reads. The table is expected to contain 2 or 3 white-space separated columns, the first of which represent the name of a given sample, and the second and third of which represent the mate 1 and (optionally) the mate 2 barcode sequence. For a detailed description, see the "Demultiplexing" section of the online documentation.

.. option:: --multiple-barcodes

    Allow for more than one barcode (pair) for each sample. If this option is not specified, AdapterRemoval will abort if multiple barcodes/barcode pairs identify the same sample.

.. option:: --barcode-orientation orientation

    Process barcodes sequences in both the barcode1-insert-barcode2 (``forward``) orientation and barcode2-insert-barcode1 (``reverse``) orientation. If ``forward`` or ``reverse`` is used, the barcode in the barcode table are assumed to be in that orientation, and the opposite sequence is generated. If ``explicit`` is used, the barcode table is expected to contain a 4th column specifying the orientation (``forward`` or ``reverse`` for each barcode), and only that orientation is used. Default is ``unspecified``.

.. option:: --normalize-orientation

    Reverse complement merged reads found to be in the reverse orientation, based on barcode orientation

.. option:: --barcode-mm n
    Maximum *total* number of mismatches allowed, when counting mismatches in both the mate 1 and the mate 2 barcode for paired reads.

.. option:: --barcode-mm-r1 n

    Maximum number of mismatches allowed for the mate 1 barcode; if not set, this value is equal to the ``--barcode-mm`` value; cannot be higher than the ``--barcode-mm`` value.

.. option:: --barcode-mm-r2 n

    Maximum number of mismatches allowed for the mate 2 barcode; if not set, this value is equal to the ``--barcode-mm`` value; cannot be higher than the ``--barcode-mm`` value.

.. option:: --demultiplex-only

    Only carry out demultiplexing using the list of barcodes supplied with --barcode-list. No other processing is done.


Reporting options
~~~~~~~~~~~~~~~~~~~~~~

.. option:: --report-only

    Write a report of the input data without performing any processing of the FASTQ reads. In addition, attempt to build a consensus adapter sequence from fully overlapping pairs of paired-end reads. The minimum overlap is controlled by ``--merge-threshold`` and the result will be compared with the values set using ``--adapter1`` and ``--adapter2``. Default is off.

.. option:: --report-title title

    Set the title used in the HTML report. Defaults to ``AdapterRemoval v3.0.0-alpha3``.

.. option:: --report-sample-rate X

    The fraction of reads to sample when generating base quality/composition curves for trimming reports. Using all data (``--report-sample-nth 1.0``) results in an 10-30% decrease in throughput and is typically not necessary, except for tiny datasets. Default is 0.10.

.. option:: --report-duplication <N>

    FastQC based duplicate detection, based on the frequency of the first N unique sequences observed. If the option is used without a value, an N of 100k is used, corresponding to FastQC defaults; a value of 0 disables the analysis. Accepts suffixes K, M, and G. Default is 0.


Logging options
---------------

--log-level name

    The minimum severity of messages to be written to stderr. Possible values are debug, info, warning, and error. Default is info.

--log-colors name

    Enable/disable the use of colors when writing log messages. If set to auto, colors will only be enabled if STDERR is a terminal and the NO_COLORS is environmental variable is not set. Possible values are auto, always, and never. Defaults to auto.

--log-progress name

    Specify the type of progress report used. If set to ``auto``, then a spinner will be used if STDERR is a terminal and the NO_COLORS environmental variable is not set, otherwise a log line will be written for every 1 million records processed. Possible values are ``auto``, ``spin``, ``log``, and ``never``. Default is ``auto``.

Window based quality trimming
-----------------------------

AdapterRemoval implements sliding window based approach to quality based base-trimming inspired by ``sickle``. If ``--trim-window-size`` is greater than or equal to 1, that number is used as the window size for all reads. If ``--trim-window-size`` is a number greater than or equal to 0 and less than 1, then that number is multiplied by the length of individual reads to determine the window size. If the window length is zero or is greater than the current read length, then the read length is used instead.

Reads are trimmed as follows for a given window size:

       1. The new 5' is determined by locating the first window where both the average quality and the quality of the first base in the window is greater than ``--trim-min-quality``.

       2. The new 3' is located by sliding the first window right, until the average quality becomes less than or equal to ``--trim-min-quality``. The new 3' is placed at the last base in that window where the quality is greater than or equal to ``--trim-min-quality``.

       3. If no 5' position could be determined, the read is discarded.


Exit status
-----------

AdapterRemoval exists with status 0 if the program ran successfully, and with a non-zero exit code if any errors were encountered. Output from AdapterRemoval should not be used if the program returned a non-zero exit code.


Reporting bugs
--------------

Please report any bugs using the AdapterRemoval issue-tracker:

https://github.com/MikkelSchubert/adapterremoval/issues


License
-------

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
at your option any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

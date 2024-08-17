AdapterRemoval manpage
======================

Synopsis
--------

**adapterremoval3** [*options*...] --in-file1 <*filenames*> [--in-file2 <*filenames*>]


Description
-----------

:program:`adapterremoval3` removes residual adapter sequences from single-end (SE) or paired-end (PE) FASTQ reads, optionally trimming Ns and low qualities bases and/or merging overlapping paired-end mates into one read. Low quality reads are filtered based on the resulting length and the number of ambiguous nucleotides ('N') present following trimming. These operations may be combined with simultaneous demultiplexing using 5' barcode sequences. Alternatively, ``adapterremoval3`` may attempt to reconstruct a consensus adapter sequences from paired-end data, in order to allow the identification of the adapter sequences originally used.

If you use this program, please cite the paper:

	Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 12;9(1):88

	http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2


For detailed documentation, please see

	http://adapterremoval.readthedocs.io/


Options
-------

.. program:: adapterremoval3

.. option:: --help

	Display summary of command-line options.

.. option:: --version

	Print the version string.

.. option:: --identify-adapters

	Attempt to build a consensus adapter sequence from fully overlapping pairs of paired-end reads. The minimum overlap is controlled by ``--merge-threshold``. The result will be compared with the values set using ``--adapter1`` and ``--adapter2``. No trimming is performed in this mode. Default is off.

.. option:: --threads n

	Maximum number of threads. Defaults to 1.


Input files
~~~~~~~~~~~

.. option:: --in-file1 filename [filenames...]

	Read FASTQ reads from one or more files, either uncompressed or gzip compressed. The ``--interleaved`` and ``--interleaved-input``  options may be used to enable reading of interleaved reads from these files.

.. option:: --in-file2 filename [filenames...]

	Read one or more FASTQ files containing mate 2 reads for a paired-end run. If specified, ``--in-file1`` must also be set.


Output file options
~~~~~~~~~~~~~~~~~~~

Output files for adapterremoval may be specified either via the ``--out-prefix`` option, which assigns default filenames, and/or via the individual ``--out-*`` to set or override the output filename for a given output read type. The same destination file may be used for multiple types of FASTQ reads, in which case the reads are written in input order.

.. option:: --out-prefix path

	Prefix for the output files for which no filename was set using the corresponding options below.

.. option:: --out-file1 filename
.. option:: --out-file2 filename

	Output files containing trimmed mate 1 reads and mate 2 reads. If interleaved output is enabled, then this file also contains mate 2 reads.

.. option:: --out-merged filename

	When used with --merged, this file contains overlapping mate-pairs which have been merged into a single read.

.. option:: --out-singleton filename

	Output file to which containing paired reads for which the mate has been discarded.

.. option:: --out-discarded filename

	Contains reads discarded due to the --min-length, --max-length or --max-ns options.

.. option:: --out-json filename
.. option:: --out-html filename

	Reports in JSON/HTML format containing information on the parameters used in the run as well as overall statistics on the reads before and after trimming.


FASTQ options
~~~~~~~~~~~~~
.. option:: --quality-format name

	The Phred quality scores encoding used in input reads - either '64' for Phred+64 (Illumina 1.3+ and 1.5+) or '33' for Phred+33 (Illumina 1.8+). In addition, the value 'solexa' may be used to specify reads with Solexa encoded scores. The 'sam' format may be used for Phred+33 data with very high quality scores. Default is 33.

.. option:: --mate-separator separator

	Character separating the mate number (1 or 2) from the read name in FASTQ records. This is typically either '/' or '.'. By default AdapterRemoval will attempt to infer this separator automatically.

.. option:: --interleaved-input

	Enable reading of interleaved FASTQ reads from the files specified with ``--in-file1``. Defaults to off.

.. option:: --interleaved-ouput

	Write paired-end reads to the file specified by ``--out-file1``, interleaving mate 1 and mate 2 reads. Defaults to off.

.. option:: --interleaved

	Enables ``--interleaved-input`` and ``--interleaved-output``. Defaults to off.


Output compression options
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: --gzip

	If set, all FASTQ files written by AdapterRemoval will be gzip compressed using the compression level specified using ``--gzip-level``. If ``--out-prefix`` is used then the ".gz" extension added to files for which no filename was specified. Gzip compression may also be enabled by manually specifying a ".gz" extension for a output files. Defaults to off.

.. option:: --gzip-level level

	Determines the compression level used when gzip'ing FASTQ files. Must be a value in the range 0 to 9, with 0 disabling compression and 9 being the best compression. For compression levels 4-9, block based compression is performed using libdeflate. This may cause compatibility issues in rare cases, which can be migitated by using a lower compression level. Defaults to 6.


FASTQ processing options
~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: --adapter1 adapter

	Adapter sequence expected to be found in mate 1 reads, specified in read direction. For a detailed description of how to provide the appropriate adapter sequences, see the "Adapters" section of the online documentation. Default is AGATCGGAAGAGCACACGTCTGAACTCCAGTCA.

.. option:: --adapter2 adapter

	Adapter sequence expected to be found in mate 2 reads, specified in read direction. For a detailed description of how to provide the appropriate adapter sequences, see the "Adapters" section of the online documentation. Default is AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT.

.. option:: --adapter-list filename

	Read one or more adapter sequences from a table. The first two columns (separated by whitespace) of each line in the file are expected to correspond to values passed to --adapter1 and --adapter2. In single-end mode, only column one is required. Lines starting with '#' are ignored. When multiple rows are found in the table, AdapterRemoval will try each adapter (pair), and select the best aligning adapters for each FASTQ read processed.

.. option:: --min-adapter-overlap length

	In single-end mode, reads are only trimmed if the overlap between read and the adapter is at least X bases long, not counting ambiguous nucleotides (N). Defaults to 0.

.. option:: --mismatch-rate rate

	The allowed fraction of mismatches allowed in the aligned region. If the value is less than 1, then the value is used directly. If ```--mismatchrate`` is greater than 1, the rate is set to 1 / ``--mismatchrate``. The default setting is 3 when trimming adapters, corresponding to a maximum mismatch rate of 1/3, and 10 when using ``--identify-adapters``.

.. option:: --shift n

	To allow for missing bases in the 5' end of the read, the program can let the alignment slip ``--shift`` bases in the 5' end. This corresponds to starting the alignment maximum ``--shift`` nucleotides into read2 (for paired-end) or the adapter (for single-end). The default is 2.

.. option:: --merge

	In paired-end mode, merge overlapping mates into a single and recalculate the quality scores. The overlap needs to be at least ``--merge-threshold`` nucleotides, with a maximum number of mismatches determined by ``--mismatch-mate``. This option has no effect in single-end mode.

.. option:: --merge-threshold length

	The minimum overlap between mate 1 and mate 2 before the reads are merged into one, when collapsing paired-end reads. Default is 11.

.. option:: --prefix-read1 X

	Adds the specified prefix to read 1 names. Default to no prefix.

.. option:: --prefix-read2 X

	Adds the specified prefix to read 2 names. Default to no prefix.

.. option:: --prefix-merged X

	Adds the specified prefix to merged read names. Default to no prefix.


Quality trimming options
~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: --pre-trim3p n [n]

	Trim the 3' of reads by a fixed amount after demultiplexing but before removing adapters. Specify one value to trim mate 1 and mate 2 reads the same amount, or two values separated by a space to trim each mate different amounts. Off by default.

.. option:: --post-trim5p n [n]

	Trim the 5' of reads by a fixed amount after removing adapters, but before carrying out quality based trimming. See ``--pre-trim3p``.

.. option:: --post-trim3p n [n]

	Trim the 3' of reads by a fixed amount after removing adapters, but before carrying out quality based trimming. See ``--pre-trim3p``.

.. option:: --trim-strategy name

	The strategy used for performing quality trimming; 'mott' to enable trimming using the modified Mott's algorithm, 'window' to perform window based quality trimming, and 'per-base' to perform base-by-base trimming of low-quality bases and Ns (if enabled). Defaults to Mott's algorithm.

.. option:: --trim-mott-rate rate

	The threshold value used when performing trimming quality based trimming using the modified Mott's algorithm. A value of zero or less disables trimming; a value greater than one is assumed to be a Phred encoded error rate (e.g. 13 ~= 0.05). Applies when Mott based trimming is enabled. Default to 0.05.

.. option:: --trim-windows size

	Trim low quality bases using a sliding window based approach inspired by :program:`sickle` with the given window size. See the "Window based quality trimming" section of the manual page for a description of this algorithm. Applies when window based trimming is enabled. Defaults to 0.1.

.. option:: --trim-min-quality minimum

	Inclusive minimum quality used when trimming low-quality bases with --trimming-strategy 'window' and 'per-base'. Applies when window based or per-base trimming is enabled. Default is 2.

.. option:: --trim-qualities

	Trim consecutive stretches of low quality bases (threshold set by ``--trim-min-quality``) from the 5' and 3' termini. If trimming of Ns is also enabled (``--trim-ns``), then stretches of mixed low-quality bases and Ns are trimmed. Applies when per-base trimming is enabled.

.. option:: --trim-ns

	Trim consecutive Ns from the 5' and 3' termini. If quality trimming is also enabled (``--trim-qualities``), then stretches of mixed low-quality bases and/or Ns are trimmed. Applies when window based or per-base trimming is enabled.

.. option:: --pre-trim-polyx nucleotides

	Enable trimming of poly-X tails prior to read alignment and adapter trimming. Zero or more nucleotides (A, C, G, T) may be specified. Zero or more nucleotides may be specified after the option seperated by spaces, with zero nucleotides corresponding to all of A, C, G, and T. Defaults to no trimming.

.. option:: --post-trim-polyx nucleotides

	Enable trimming of poly-X tails after read alignment and adapter trimming/merging, but before trimming of low-quality bases. Merged reads are not trimmed by this option (both ends are 5'). See `--pre-trim-polyx`. Off by default.

.. option:: --preserve5p

	If set, bases at the 5p will not be trimmed by ``--trim-mott-rate``. Merged reads will not be quality trimmed when this option is enabled due to the 3' ends being located inside the reads or overlapping the 5' of the source sequences.


Filtering options
~~~~~~~~~~~~~~~~~

.. option:: --max-ns n

	Discard reads containing more than ``--max`` ambiguous bases ('N') after trimming. Default is no maximum.


.. option:: --min-length length

	Reads shorter than this length are discarded following trimming. Defaults to 15.

.. option:: --max-length length

	Reads longer than this length are discarded following trimming. Defaults to no maximum.


Demultiplexing options
~~~~~~~~~~~~~~~~~~~~~~

.. option:: --barcode-list filename

	Perform demultiplxing using table of one or two fixed-length barcodes for SE or PE reads. The table is expected to contain 2 or 3 columns, the first of which represent the name of a given sample, and the second and third of which represent the mate 1 and (optionally) the mate 2 barcode sequence. For a detailed description, see the "Demultiplexing" section of the online documentation.

.. option:: --barcode-mm n
	Maximum number of mismatches allowed when counting mismatches in both the mate 1 and the mate 2 barcode for paired reads.

.. option:: --barcode-mm-r1 n

	Maximum number of mismatches allowed for the mate 1 barcode; if not set, this value is equal to the ``--barcode-mm`` value; cannot be higher than the ``--barcode-mm`` value.

.. option:: --barcode-mm-r2 n

	Maximum number of mismatches allowed for the mate 2 barcode; if not set, this value is equal to the ``--barcode-mm`` value; cannot be higher than the ``--barcode-mm`` value.

.. option:: --demultiplex-only

	Only carry out demultiplexing using the list of barcodes supplied with --barcode-list. No other processing is done.


Window based quality trimming
-----------------------------

AdapterRemoval implements sliding window based approach to quality based base-trimming inspired by ``sickle``. If ``window_size`` is greater than or equal to 1, that number is used as the window size for all reads. If ``window_size`` is a number greater than or equal to 0 and less than 1, then that number is multiplied by the length of individual reads to determine the window size. If the window length is zero or is greater than the current read length, then the read length is used instead.

Reads are trimmed as follows for a given window size:

       1. The new 5' is determined by locating the first window where both the average quality and the quality of the first base in the window is greater than ``--minquality``.

       2. The new 3' is located by sliding the first window right, until the average quality becomes less than or equal to ``--minquality``. The new 3' is placed at the last base in that window where the quality is greater than or equal to ``--minquality``.

       3. If no 5' position could be determined, the read is discarded.


Exit status
-----------

AdapterRemoval exists with status 0 if the program ran succesfully, and with a non-zero exit code if any errors were encountered. Do not use the output from AdapterRemoval if the program returned a non-zero exit code!


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
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

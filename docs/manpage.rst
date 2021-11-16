AdapterRemoval manpage
======================

Synopsis
--------

**AdapterRemoval** [*options*...] --file1 <*filenames*> [--file2 <*filenames*>]


Description
-----------

:program:`AdapterRemoval` removes residual adapter sequences from single-end (SE) or paired-end (PE) FASTQ reads, optionally trimming Ns and low qualities bases and/or collapsing overlapping paired-end mates into one read. Low quality reads are filtered based on the resulting length and the number of ambigious nucleotides ('N') present following trimming. These operations may be combined with simultaneous demultiplexing using 5' barcode sequences. Alternatively, ``AdapterRemoval`` may attempt to reconstruct a consensus adapter sequences from paired-end data, in order to allow the identification of the adapter sequences originally used.

If you use this program, please cite the paper:

	Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 12;9(1):88

	http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2


For detailed documentation, please see

	http://adapterremoval.readthedocs.io/


Options
-------

.. program:: AdapterRemoval

.. option:: --help

	Display summary of command-line options.

.. option:: --version

	Print the version string.

.. option:: --file1 filename [filenames...]

	Read FASTQ reads from one or more files, either uncompressed, bzip2 compressed, or gzip compressed. This contains either the single-end (SE) reads or, if paired-end, the mate 1 reads. If running in paired-end mode, both ``--file1`` and ``--file2`` must be set. See the primary documentation for a list of supported formats.

.. option:: --file2 filename [filenames...]

	Read one or more FASTQ files containing mate 2 reads for a paired-end run. If specified, ``--file1`` must also be set.

.. option:: --identify-adapters

	Attempt to build a consensus adapter sequence from fully overlapping pairs of paired-end reads. The minimum overlap is controlled by ``--minalignmentlength``. The result will be compared with the values set using ``--adapter1`` and ``--adapter2``. No trimming is performed in this mode. Default is off.

.. option:: --threads n

	Maximum number of threads. Defaults to 1.


FASTQ options
~~~~~~~~~~~~~
.. option:: --qualitybase base

	The Phred quality scores encoding used in input reads - either '64' for Phred+64 (Illumina 1.3+ and 1.5+) or '33' for Phred+33 (Illumina 1.8+). In addition, the value 'solexa' may be used to specify reads with Solexa encoded scores. Output is always written as Phred+33. Default is 33.

.. option:: --qualitymax base

	Specifies the maximum Phred score expected in input files, and used when writing output files. Possible values are 0 to 93 for Phred+33 encoded files, and 0 to 62 for Phred+64 encoded files. Defaults to 41.

.. option:: --mate-separator separator

	Character separating the mate number (1 or 2) from the read name in FASTQ records. Defaults to '/'.

.. option:: --interleaved

	Enables ``--interleaved-input`` and ``--interleaved-output``.

.. option:: --interleaved-input

	If set, input is expected to be a interleaved FASTQ files specified using ``--file1``, in which pairs of reads are written one after the other (e.g. read1/1, read1/2, read2/1, read2/2, etc.).

.. option:: --interleaved-ouput

	Write paired-end reads to a single file, interleaving mate 1 and mate 2 reads. By default, this file is named ``basename.paired.truncated``, but this may be changed using the ``--output1`` option.


Output file options
~~~~~~~~~~~~~~~~~~~
.. option:: --basename filename

	Prefix used for the naming output files, unless these names have been overridden using the corresponding command-line option (see below).

.. option:: --settings file

	Output file containing information on the parameters used in the run as well as overall statistics on the reads after trimming. Default filename is 'basename.settings'.

.. option:: --output1 file

	Output file containing trimmed mate1 reads. Default filename is 'basename.pair1.truncated' for paired-end reads, 'basename.truncated' for single-end reads, and 'basename.paired.truncated' for interleaved paired-end reads.

.. option:: --output2 file

	Output file containing trimmed mate 2 reads when ``--interleaved-output`` is not enabled. Default filename is 'basename.pair2.truncated' in paired-end mode.

.. option:: --singleton file

	Output file to which containing paired reads for which the mate has been discarded. Default filename is 'basename.singleton.truncated'.

.. option:: --outputmerged file

	If --merged is set, contains overlapping mate-pairs which have been merged into a single read. This does not include which have subsequently been trimmed due to low-quality or ambiguous nucleotides. Default filename is 'basename.merged'

.. option:: --discarded file

	Contains reads discarded due to the --minlength, --maxlength or --maxns options. Default filename is 'basename.discarded'.


Output compression options
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: --gzip

	If set, all FASTQ files written by AdapterRemoval will be gzip compressed using the compression level specified using ``--gzip-level``. The extension ".gz" is added to files for which no filename was given on the command-line. Defaults to off.

.. option:: --gzip-level level

	Determines the compression level used when gzip'ing FASTQ files. Must be a value in the range 0 to 9, with 0 disabling compression and 9 being the best compression. Defaults to 6.

.. option:: --bzip2

	If set, all FASTQ files written by AdapterRemoval will be bzip2 compressed using the compression level specified using ``--bzip2-level``. The extension ".bz2" is added to files for which no filename was given on the command-line. Defaults to off.

.. option:: --bzip2-level level

	Determines the compression level used when bzip2'ing FASTQ files. Must be a value in the range 1 to 9, with 9 being the best compression. Defaults to 9.


FASTQ trimming options
~~~~~~~~~~~~~~~~~~~~~~

.. option:: --adapter1 adapter

	Adapter sequence expected to be found in mate 1 reads, specified in read direction. For a detailed description of how to provide the appropriate adapter sequences, see the "Adapters" section of the online documentation. Default is AGATCGGAAGAGCACACGTCTGAACTCCAGTCA.

.. option:: --adapter2 adapter

	Adapter sequence expected to be found in mate 2 reads, specified in read direction. For a detailed description of how to provide the appropriate adapter sequences, see the "Adapters" section of the online documentation. Default is AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT.

.. option:: --adapter-list filename

	Read one or more adapter sequences from a table. The first two columns (separated by whitespace) of each line in the file are expected to correspond to values passed to --adapter1 and --adapter2. In single-end mode, only column one is required. Lines starting with '#' are ignored. When multiple rows are found in the table, AdapterRemoval will try each adapter (pair), and select the best aligning adapters for each FASTQ read processed.

.. option:: --minadapteroverlap length

	In single-end mode, reads are only trimmed if the overlap between read and the adapter is at least X bases long, not counting ambiguous nucleotides (N); this is independent of the ``--minalignmentlength`` when using ``--merge``, allowing a conservative selection of putative complete inserts in single-end mode, while ensuring that all possible adapter contamination is trimmed. The default is 0.

.. option:: --mm mismatchrate

	The allowed fraction of mismatches allowed in the aligned region. If the value is less than 1, then the value is used directly. If ```--mismatchrate`` is greater than 1, the rate is set to 1 / ``--mismatchrate``. The default setting is 3 when trimming adapters, corresponding to a maximum mismatch rate of 1/3, and 10 when using ``--identify-adapters``.

.. option:: --shift n

	To allow for missing bases in the 5' end of the read, the program can let the alignment slip ``--shift`` bases in the 5' end. This corresponds to starting the alignment maximum ``--shift`` nucleotides into read2 (for paired-end) or the adapter (for single-end). The default is 2.

.. option:: --trim5p n [n]

	Trim the 5' of reads by a fixed amount after removing adapters, but before carrying out quality based trimming. Specify one value to trim mate 1 and mate 2 reads the same amount, or two values separated by a space to trim each mate different amounts. Off by default.

.. option:: --trim3p n [n]

	Trim the 3' of reads by a fixed amount. See ``--trim5p``.

.. option:: --trimns

	Trim consecutive Ns from the 5' and 3' termini. If quality trimming is also enabled (``--trimqualities``), then stretches of mixed low-quality bases and/or Ns are trimmed.

.. option:: --maxns n

	Discard reads containing more than ``--max`` ambiguous bases ('N') after trimming. Default is 1000.

.. option:: --trimqualities

	Trim consecutive stretches of low quality bases (threshold set by ``--minquality``) from the 5' and 3' termini. If trimming of Ns is also enabled (``--trimns``), then stretches of mixed low-quality bases and Ns are trimmed.

.. option:: --trimwindows window_size

	Trim low quality bases using a sliding window based approach inspired by :program:`sickle` with the given window size. See the "Window based quality trimming" section of the manual page for a description of this algorithm.

.. option:: --minquality minimum

	Set the threshold for trimming low quality bases using ``--trimqualities`` and ``--trimwindows``. Default is 2.

.. option:: --preserve5p

	If set, bases at the 5p will not be trimmed by ``--trimns``, ``--trimqualities``, and ``--trimwindows``. Collapsed reads will not be quality/N trimmed when this option is enabled.

.. option:: --minlength length

	Reads shorter than this length are discarded following trimming. Defaults to 15.

.. option:: --maxlength length

	Reads longer than this length are discarded following trimming. Defaults to 4294967295.



FASTQ merging options
~~~~~~~~~~~~~~~~~~~~~

.. option:: --merge

	In paired-end mode, merge overlapping mates into a single and recalculate the quality scores. The overlap needs to be at least ``--minalignmentlength`` nucleotides, with a maximum number of mismatches determined by ``--mm``. This option has no effect in single-end mode.

.. option:: --minalignmentlength length

	The minimum overlap between mate 1 and mate 2 before the reads are merged into one, when collapsing paired-end reads, or when attempting to identify complete template sequences in single-end mode. Default is 11.

.. option:: --seed seed

	When collaping reads at positions where the two reads differ, and the quality of the bases are identical, AdapterRemoval will select a random base. This option specifies the seed used for the random number generator used by AdapterRemoval. This value is also written to the settings file. Note that setting the seed is not reliable in multithreaded mode, since the order of operations is non-deterministic.

.. option:: --merge-deterministic

	Enable deterministic mode; currently only affects --merge, different overlapping bases with equal quality are set to N quality 0, instead of being randomly sampled. Setting this option also sets --merge.

.. option:: --merge-conservatively

	Alternative merging algorithm inspired by FASTQ-join: For matching overlapping bases, the highest quality score is used. For mismatching overlapping bases, the highest quality base is used and the quality is set to the absolute difference in Phred-score between the two bases. For mismatching bases with identical quality scores, the base is set to 'N' and the quality score to 0 (Phred-encoded). Setting this option also sets --merge.


FASTQ demultiplexing options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

As of v2.2.2, AdapterRemoval implements sliding window based approach to quality based base-trimming inspired by ``sickle``. If ``window_size`` is greater than or equal to 1, that number is used as the window size for all reads. If ``window_size`` is a number greater than or equal to 0 and less than 1, then that number is multiplied by the length of individual reads to determine the window size. If the window length is zero or is greater than the current read length, then the read length is used instead.

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

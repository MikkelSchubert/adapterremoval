

### Version 2.2.2 - 2017-07-17

  * Made gzip and bzip2 support mandatory.
  * Added support for Intel compilers, courtesy of Kevin Murray (kdmurray91).


### Version 2.2.1a - 2017-05-17

  * Fixed compilation on OSX.


### Version 2.2.1 - 2017-05-15

  * Numerous spelling errors fixed courtesy of Andreas Tille.
  * Added support for specifying multiple filenames after --file1 and --file2,
    in which case the files are treated as if they were concatenated. This is
    supported for all operations. Special thanks to Stephen Clayton.
  * Added additional run-time checks to catch race-conditions.
  * Progress messages written to STDERR no longer cause subsequent error
    messages to be written to the same line.
  * Implemented quality trimming using a sliding window approach inspired by
    sickle (https://github.com/najoshi/sickle). Special thanks to Kevin Murray.
  * Fixed miscounting of the total number of retained nucleotides, where mate 1
    reads were being counted twice instead of counting both mate 1 and mate 2.


### Version 2.2.0 - 2016-10-27

  * AdapterRemoval now requires a C++11 compliant compiler; furthermore,
    multithreading is no longer an optional feature, as this is now
    implemented using the C++11 instead of directly calling pthreads.
  * Add explicit message not to use the results after failed runs.
  * Minor changes to .settings: Adapter numbers now 1-based; the 'Number of
    reads with adapters' is changed to 'Number of read pairs with adapters'
    when trimming PE reads; the 'Average read length of trimmed reads' is
    changed to 'Average read length of retained reads' for clarity.
  * Dropped the undocumented 'poor' classification for alignments; for
    statistical purposes, reads are either counted as aligned or not aligned.
    This ony changes how results are presented in the .settings files.
  * Rework selection of nucleotides at overlapping positions with the same
    quality, in order to prevent potential data-races during tie-breaking,
    when running in multi-threaded mode.
  * Added support for reading FASTQ files using Windows-style newlines (\r\n).
  * AdapterRemoval will not print a warning to STDERR if the same command-line
    option is specified multiple times.
  * Reworked handling of barcodes to avoid unnecessary memory allocations,
    which would cause problems when using long barcodes.
  * Added support for combining output files; this is enabled using the
    --combined-output option, and ensures that all reads are written to the
    same file, or pair of files (for non-interleaved PE reads). The sequence
    of reads that fail are replaced with a single 'N' with quality score 0.
  * Fixed bug in the counting of singleton reads used in '.settings' files.
  * Fixed mis-placement of underscore when pretty printing adapter sequences
    that included barcodes.
  * Fixed misprinting of mate 2 adapter sequences in the .settings file;
    these would be printed in the reverse complemented orientation, relative
    to how they were specified on the command-line.


### Version 2.1.7 - 2016-03-11

  * The mate number is now stripped from collapsed reads, where previously this
    would always be '\1' (if set). However, if meta-data is present in the
    reads, that found in the mate 1 read is retained.
  * The value used for --mate-separator is now written to the 'settings' file.
  * Improved 'make install'. This command now makes use of a PREFIX value to
    determine the installation destination (defaults to /usr/local), and
    includes the 'README.md' file and 'examples' folder in the installation.
  * Improved 'make test'. This command now attempts to download the required
    testing library automatically, using either wget or curl if available.


### Version 2.1.6 - 2016-03-03

  * Added support for reading / writing interleaved FASTQ files; this is
    enabled by the options --interleaved-input and --interleaved-output,
    respectively, or by setting --interleaved option which implies both of
    the former options. See the README for an example.
  * Fixed bug in a sanity check meant to detect if the mate 1 and mate 2 files
    were of unequal length. This is now correctly detected in all cases.
  * Expanded README with information about reading / writing FASTQ files with
    different PHRED encodings / maximum quality scores.

### Version 2.1.5 - 2016-02-19

  * Added the --mate-separator option, which specifies the character separating
    the mate number; by default this is '/', and AdapterRemoval will therefore
    identify mate numbers if read-names end with "/1" or "/2".
  * Fixed race condition which could result in premature termination when using
    --gzip or --bzip2 together with the --threads options.
  * Improved checks during compression and sanity checks following processing.


### Version 2.1.4 - 2016-02-09

  * Fixed bug which could occasionally result in failure when bzip2 compression
    was enabled, by attempting to compress empty buffer.
  * The following was contributed by Hannes Pétur Eggertsson:
    * Wrapped code in 'ar' namespace, and made it possible to compile
      AdapterRemoval as a static library (via the command 'make static'),
      allowing it to be used as part of other projects.
    * Updated instructions for installing GTest library using new repository.
    * Fixed typos.


### Version 2.1.3 - 2015-12-25

  * Added option --minadapteroverlap, which sets a minimum alignment length
    when carrying out trimming of single-end reads. The default (0) may result
    in an excess of false postiives around (1 - 2 bp long), which may be
    mitigated by running AdapterRemoval with '--minadapteroverlap 3'.
  * Greatly expanded README.md, adding several examples with test data included
    in the 'examples' folder, demonstrating common usage of the program.
  * Updated man-page with missing information and rewrote several parts.
  * Updated the help-text for several command-line options.
  * Avoid writing information to stdout, so that (SE) trimming can be piped.
    This can be accomplished by using the option --output1 /dev/stdout.
  * Fixed the --seed option, which was not properly applied during runtime.


### Version 2.1.2 - 2015-10-08

  * Changed the way "full-length" and "truncated collapsed" reads are counted
    in the .settings file; previously, all collapsed reads were counted, even
    if these were subsequently discarded (due to the length). Now only retained
    reads are counted, matching the behavior of AdapterRemoval v1.x.
  * Added setup instructions when running 'make test' for the first time.


### Version 2.1.1 - 2015-09-14

  * Fixed broken assert preventing the use of --adapter-list.
  * Fixed bug using --qualitybase-output for both input and output.


### Version 2.1.0 - 2015-09-08

Major changes:
  * Support for (transparently) reading and writing bzip2 files.
  * Parallelization of adapter trimming and identification using pthreads; the
    number of threads used is specified using --threads. Note that only one
    thread is allowed to perform IO (reads / writes) at a time, to prevent
    clobbering the disk, but compression (if enabled) is performed in parallel.
  * Support for combined demultiplexing and adapter-removal using the
    --barcode-list command-line option; when demultiplexing, the barcodes
    identified for a given read is added to the adapter sequence, in ordre to
    ensure correct trimming of the reads.
  * Features depending on external libraries (gzip, bzip2, and threading
    support) can be disabled in the Makefile on systems lacking these
    libraries.

Other changes / bug-fixes:
  * Display currently specified --adapter1 / adapter2 sequences for comparison
    when attempting to infer adapter sequences. Only the first pair is used, if
    multiple adapter pairs are specified.
  * Sites with no majority-base during adapter-identification are set to N.
  * Fixed failure to read of barcode sequences (--5prime / --5prime-list).
  * Progress report now shows total number of reads processes, for both single
    ended and pair ended analyses.
  * FASTQ reads with Solexa scores are now output as Solexa scores by default,
    rather than Phred+64. Note that the program represents quality scores using
    Phred scores internally, resulting in a lossy conversion. It is therefore
    recommended to convert to Phred scores rather than use Solexa scores.


### Version 2.0.0 - 2014-03-10

Version 2.0.0 of AdapterRemoval is a near complete rewrite, with the goal of
improved safety, increased speed, fixing a number of minor issues with
previous versions of AdapterRemoval, and adding a few new features:

Compatibility changes:
  * Command-line arguments --pcr1 and --pcr2 have been deprecated in favor of
    --adapter1 and --adapter2. While --pcr1 and --adapter1 are equivalent,
    --adapter2 expects the adapter sequence which may be observed in raw
    mate 2 reads, unlike --pcr2 which expected the sequence which could be
    observed in the reverse complement of mate 2 reads (cf. the README).
  * The use of --file1 and (optionally) --file2 is now required; reads will not
    be read from STDIN, nor written to STDOUT by default. To approximate the
    previous behavior, the following command may be used:
    $ AdapterRemoval --file1 /dev/stdin --output1 /dev/stdout
  * Per-read statistics of adapter / low-quality base trimming using --stats is
    no longer supported.

Major changes:
  * Strict validation of input FASTQ records, to ensure that records are well
    formed, that quality scores fall within the expected range given the
    specified format/offset, and more.
  * Limited support for Solexa quality scores; these are converted to and
    saved as Phred+33 or Phred+64 encoded scores.
  * Improved handling of asymmetric read-pairs, in which the length of the
    mate 1 read differs from the length of the mate 2 read.
  * Significant improvements in performance, resulting in a ~5x increase in the
    rate of adapter trimming in basic version, and a ~20x increase in the rate
    of adapter trimming in the SSE enabled version (the default).
  * Support for multiple adapter sequences as well as multiple barcode
    sequences; AdapterRemoval will favor the highest scoring alignment,
    favoring longer alignments over shorter alignments with the same score,
    and favoring alignments with the fewest ambiguous bases (N) involved if
    the score and length is identical.
  * If --collapse is set in single-ended mode, "collapsed" reads will be
    identified using the same criteria as for paired-ended mode, i.e. requiring
    that at least --minalignmentlen bases overlap, and written to .collapsed
    and .collapsed.truncated. This allows for the identification of reads
    that are complete inserts.
  * Added the ability to identify adapter sequences for paired-ended reads, by
    identifying reads which extends past the ends of the template sequence, and
    extracting the adapters from these.
  * Added support for reading / writing gzipped compressed FASTQ files; if
    enabled (using the --gzip flag), the ".gz" extension is added to filenames,
    unless the filenames are explicitly specified by the user.
  * Length distributions are now calculated per read-type post-trimming
    (mate 1, mate 2, collapsed, etc.) and written to the .settings file.


Other improvements / bug-fixes:
  * Barcodes may now contain Ns.
  * Fixed underestimation of error-probabilities during sequence collapse.
  * Fixed (futher) underestimation of error-probabilities of bases during
    collapsing, for conflicting base-calls with the same Phred score.
  * Fixed the maximum number of mismatches for alignments in the range of
    6 .. 9 bases always being 1, even if --mm was set to 0.
  * Fixed the maximum number of mismatches for alignments being calculated
    based on the length of the alignment including ambiguous bases (N),
    thereby inflating the number of mismatches allowed for poor alignments.
  * Replaced use of lower bits of rand() calls with random(), as the former
    generates low entropy bits in that range on some (non-Linux) platforms.
  * Fixed well-aligned reads being discarded due to the minimum-length
    requirement after trimming not being counted as well-aligned, resulting
    in the total number of alignments not matching the total number of reads.
  * Fixed bug in shifts for PE reads, which was causing some alignments to be
    missed for adapter-only (i.e. no insert sequence) sequences.
  * Improved input validation and sanity checks for command-line parameters.
  * It is now possible to explicitly specify the RNG seed, to allow individual
    runs to be reproduced; the seed is also written to the .settings file.
  * Seed is now initialized using a mix of seconds and microseconds, instead of
    the current time in seconds, to reduce the risk of multiple instances
    spawed within a short timespan from using the same seed.
  * An (optional) progress report is printed during usage, incidating the
    run-time and number of reads processed.


### Version 1.5.4 - 2014-04-23

  * Fixed bug in which collapsed reads would not be considered truncated if
    bases were trimmed from the 5' end.
  * Fixed bug in which the quality bases used for mate 2 during collapsing of
    overlapping read pairs made use of quality scores with a wrong orientation.
  * Reduced the amount of IO operations during trimming.


### Version 1.5.2 - 2013-10-22

Two changes to the program:
  * I have added a reference to the paper to both the man page and the help
     text.
  * I fixed a minor bug in the collapse code where two very low quality bases
    might give rise to a third low quality base being called. For example, a C
    with quality " and a T with quality ! would result in an A with quality #.
    This has been fixed so that the the result is now C with quality ".


### Version 1.5.0 - 2013-04-29

Small update: Due to user feedback, the program now outputs collapsed pairs in
two files: One contains full-length collapsed pairs constituting the full
insert, the other contains collapsed pairs that have been truncated due to low
qualities or Ns in the reads.


### Version 1.4.0 - 2013-03-24

I have made some fixes to the program:
  * The program can now handle the use of '.' instead of 'N' to encode
    undefined nucleotides.
  * There was a typo in the adapter sequence used for PCR2!
  * Some minor changes to output etc.


### Version 1.3.0 - 2013-02-10

I have updated AdapterRemoval and released version 1.3. These changes are based
on feedback from users of the program that had some very specific and well-
founded suggestions. Some of these changes are minor, others will have more
dramatic effects on the use of the program so please read these notes
carefully:

Minor changes:
  * I fixed an occasional segmentation fault.
  * Collapsed reads are now names "@M_...".
  * Collapsed reads are put in a separate file with extension ".collapsed".

Important changes:
  * The sequences PCR1 and PCR2 are now used as-is without reverse-
    complementation. You have to make sure that the sequences you search for
    are correct.
  * The default PCR1 and PCR2 sequences are now:
    PCR1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
    PCR2: AATGATACGGCGACCACCGAGATCACACTCTTTCCCTACACGACGCTCTTCCGATCT
  * I have changed the way PCR1 and PCR2 are used to make the program
    consistent. Now, you always search for the sequence PCR1 in READ1 (whether
    single end or paired end), and you search for PCR2 in READ2. In single end
    mode, this corresponds to having an empty READ2 and ignore PCR2 as
    illustrated below:

      * For paired end data, PCR2-READ1 aligned to READ2-PCR1.
      * For single end data, READ1 aligned to PCR1.

As always, please contact me with any questions or comments.

Stinus


### Version 1.1.0 - 2012-05-01

  * It is now possible to look for adapter in the 5' end of reads using the
    --5prime parameter.
  * Updated trimming of qualities.
  * Added option for discarding reads with too many gaps using --maxns max.
  * The programs handles lower vs upper case issues by translating all
    sequences to upper case.
  * The program now checks for inconsistent parameters.
  * Fixed some typographical issues with output.

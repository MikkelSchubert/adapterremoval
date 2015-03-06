=========
Changelog
=========

Version 2.0.0 - 2014-..-..
==========================
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
  * Improved handling of asymetric read-pairs, in which the length of the
    mate 1 read differs from the length of the mate 2 read.
  * Significant improvements in performance, resulting in a ~5x increase in the
    rate of adapter trimming in basic version, and a ~20x increase in the rate
    of adapter trimming in the SSE enabled version (the default).
  * Support for multiple adapter sequences as well as multiple barcode
    sequences; AdapterRemoval will favor the highest scoring alignment,
    favoring longer alignments over shorter alignments with the same score,
    and favoring alignments with the fewest ambigous bases (N) involved if
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
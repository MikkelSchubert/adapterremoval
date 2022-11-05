Migrating from older version
============================

Migrating AdapterRemoval to 3.x from 2.x
----------------------------------------

AdapterRemoval v3 makes a number of changes in default behavior, including the renaming of a number of options (old names may still be used, but are deprecated) and the renaming of all output filenames. See the changelog for an exhaustive list of changes in behavior.

**Adapter trimming**
 * Default adapters have been changed to the `recommended Illumina sequences`_. This makes the default settings more generally applicable, but may result in slightly different results. If nessesary, the old default adapter sequences can be specified using `--adapter1` and `--adapter2`.

**Read merging**
 * Read-merging is now assigns `N` to any mismatching position where both candidate bases have the same quality score, to make results deterministic, instead of picking one of the two at random, and picks the higest quality score for each merged position rather than calculate (roughly) the sum of quality scores, to prevent overestimation of base qualities.
 * Merged reads are no longer given a `M_` name prefix by default, but this may be re-enabled by using the new `--prefix-merged` option if this is needed.
 * It is no longer possible to "merge" unpaired reads.

**Qualty trimming and filtering**
 * The various quality trimming options `--trimwindows`, `--trimns`, `--trimqualities`, and `--minquality` have been removed in favor of the new `--trim-error-rate` option, which uses the modified Mott's algorithm. This is enabled by default with a threshold of `0.05`.
 * Filtering by sequence complexity using an method inspired by fastp has been added and enabled by default, requiring that 30% of positions differ from the previous position (`--min-complexity 0.3`).

**Reports**
 * The `.settings` file has been replaced by two reports, one report in HTML intended for human consumption and one report in JSON format intended for machine consumption.

**Misc**
 * FASTQ are now always written using Phred+33 quality scores. If you need FASTQ data in another format, then use a tool like `seqtk` to convert the output.
 * The `--qualitymax` option has been removed, and AdapterRemoval instead uses loser criteria to verify that input data is valid. Support for the full range valid of quality scores has been added via the `--quality-format sam` option.
 * The `--combined-output` option has been removed, but it is now possible to arbitrarily combine output files simply by specifying the same filename for multiple outputs. However, the ability to maintain pairing by writing 0-length trimmed reads as a single `N` was removed due to a negative impact on some downstream analyses.


Migrating to AdapterRemoval 2.x from 1.x
----------------------------------------

Command-line options mostly behave the same between AdapterRemoval v1 and AdapterRemoval v2, and scripts written with AdapterRemoval v1.x in mind should work with AdapterRemoval v2.x. A notable exception is the ``--pcr1`` and ``--pcr2`` options, which have been replaced by the ``--adapter1`` and ``--adapter2`` options described above. While the ``--pcr`` options are still supported for backwards compatibility, these should not be used going forward.

The difference between these two options is that ``--adapter2`` expects the mate 2 adapter sequence to be specified in the read orientation as described above, while the ``--pcr2`` expects the sequence to be in the same orientation as the mate 1 sequence, the reverse complement of the sequence observed in the mate 2 reads.

Using the common 13 bp Illumina adapter sequence (AGATCGGAAGAGC) as an example, this is how the options would be used in AdapterRemoval v2.x::

	adapterremoval3 --adapter1 AGATCGGAAGAGC --adapter2 AGATCGGAAGAGC ...

And in AdapterRemoval v1.x::

	adapterremoval3 --adapter1 AGATCGGAAGAGC --adapter2 GCTCTTCCGATCT ...


.. _recommended illumina sequences: https://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html

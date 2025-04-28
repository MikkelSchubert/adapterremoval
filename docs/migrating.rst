Migrating from older versions
=============================

Migrating AdapterRemoval to 3.x from 2.x
----------------------------------------

AdapterRemoval v3 makes a number of changes in default behavior, including the renaming of a number of options (old names may still be used, but are deprecated) and the renaming of all output filenames. See the changelog for an exhaustive list of changes in behavior.

**Output**
 * AdapterRemoval will no longer default to writing all files with the prefix ``your_output``. Instead, you should either set a prefix via the ``--out-prefix`` option (previously ``--basename``) or select just the reads you want to save via the various ``--out-*`` options.
 * Default filenames have been updated and all use the ``.fastq.gz`` extension. This can be overridden by using the various ``--out-*`` options on a per-file basis.
 * Discarded reads are no longer saved by default, even when using ``--out-prefix`` option. These can be saved by setting the ``--out-discarded`` option.
 * Output when using the ``--out-prefix`` option is now gzip compressed by default. This can be changed by using ``--out-format fastq``.
 * The files for discarded and singleton reads are not created if no filtering options are enabled.

**Adapter trimming**
 * Default adapters have been changed to the `recommended Illumina sequences`_. This makes the default settings more generally applicable, but may result in slightly different results. If necessary, the old default adapter sequences can be specified using ``--adapter1`` and ``--adapter2``.
 * The default ``--mm`` value has been increased to 6.

**Read merging**
 * Read-merging now assigns ``N`` to any mismatching position where both candidate bases have the same quality score, to make results deterministic, instead of picking one of the two at random, and picks the highest quality score for each merged position rather than calculate. A slightly more conservative version of the original score recalculation is available via ``--merge-strategy additive``, that sums the score of matches and subtracts the scores of mismatches.
 * Merged reads are no longer given a ``M_`` name prefix by default, but this may be re-enabled by using ``--prefix-merged M_``.
 * It is no longer possible to "merge" unpaired reads.

**Quality trimming and filtering**
 * The various quality trimming options ``--trimwindows``, ``--trimns``, ``--trimqualities``, and ``--minquality`` have been deprecated in favor of the new ``--trim-mott-rate`` option, which uses the modified Mott's algorithm. This is enabled by default with a threshold of ``0.05``, but the algorithm used may be changed back (or disabled) using the new ``--quality-trimming`` option.

**Reports**
 * The ``.settings`` file has been replaced by two reports, one report in HTML intended for human consumption and one report in JSON format intended for machine consumption. A JSON schema for the JSON report is included in the source distribution.

**Misc**
 * FASTQ are now always written using Phred+33 quality scores. If you need FASTQ data in another format, then use a tool like `seqtk`_ to convert the output.
 * The ``--qualitymax`` option has been removed, and AdapterRemoval instead uses looser criteria to verify that input data is valid. Support for the full range valid of quality scores has been added via the ``--quality-format sam`` option.
 * The ``--combined-output`` option has been removed, but it is now possible to arbitrarily combine output files simply by specifying the same filename for multiple outputs. However, the ability to maintain pairing by writing 0-length trimmed reads as a single `N` was removed due to a negative impact on some downstream analyses.
 * AdapterRemoval now defaults to using 2 threads. This can be changed using the ``--threads`` option.
 * The default (gzip) compression level has been decreased in order to increase throughput at the cost of output files being about 3% larger. This can be changed using the ``--compression-level`` option, that now gives access to a wider range of compression levels (0 to 13).

Migrating to AdapterRemoval 2.x from 1.x
----------------------------------------

Command-line options mostly behave the same between AdapterRemoval v1 and AdapterRemoval v2, and scripts written with AdapterRemoval v1.x in mind should work with AdapterRemoval v2.x. A notable exception is the ``--pcr1`` and ``--pcr2`` options, which have been replaced by the ``--adapter1`` and ``--adapter2`` options described above. While the ``--pcr`` options are still supported for backwards compatibility, these should not be used going forward.

The difference between these two options is that ``--adapter2`` expects the mate 2 adapter sequence to be specified in the read orientation as described above, while the ``--pcr2`` expects the sequence to be in the same orientation as the mate 1 sequence, the reverse complement of the sequence observed in the mate 2 reads.

Using the common 13 bp Illumina adapter sequence (AGATCGGAAGAGC) as an example, this is how the options would be used in AdapterRemoval v2.x::

    adapterremoval3 --adapter1 AGATCGGAAGAGC --adapter2 AGATCGGAAGAGC ...

And in AdapterRemoval v1.x::

    adapterremoval3 --adapter1 AGATCGGAAGAGC --adapter2 GCTCTTCCGATCT ...


.. _recommended illumina sequences: https://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
.. _seqtk: https://github.com/lh3/seqtk

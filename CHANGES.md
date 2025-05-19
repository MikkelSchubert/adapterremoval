# Changelog

## [3.0.0-alpha3] - 2025-05-19

This is the third alpha release of AdapterRemoval v3. As with the previous alpha
releases, changes that affect how AdapterRemoval is used (e.g. by removing
options) or that result in different output compared to previous versions are
marked with the label "[**BREAKING**]".

AdapterRemoval now uses `meson` for its build process, and `meson` is therefore
a build-time requirement. A `Makefile` is still provided to simplify setting up
and running the build. See the installation instructions in the documentation
for more information.

Major changes include support for hardware accelerated alignments using NEON on
modern Apple hardware, support for samples being identified by multiple
barcodes, support for handling barcodes in that may ligate in different
orientations, improved support for SAM/BAM output, and (optional) duplication
plot in HTML report.

### Added

- Multiple barcodes/barcode pairs may now be used to identify the same sample,
  via the `--multiple-barcodes` flag. The number of hits per barcode/barcode
  pair is reported in the HTML/JSON reports.
- Added support for handling barcodes that may ligate in different orientations
  (via `--barcode-orientation`) and for normalizing the orientation of merged
  reads (via `--normalize-orientation`).
- The `--use-colors` parameter may now be used to controls color output.
  Options are auto (default; enabled when run interactively), always, or never.
- The title of the HTML report can now be set via `--report-title`.
- Input files are now checked for duplicate filenames, in order to help
  prevents accidental data duplication.
- Alignments are now accelerated on Apple hardware using NEON instructions, for
  a roughly 3-fold increase in throughput.
- A duplication plot is now included in the HTML report if this is enabled,
  instead of only being reported in the JSON file.

### Changed

- [**BREAKING**] Changed  `CO` tags for read-groups in SAM/BAM files to `DS`
  (description) tags, in order to match the specification.
- [**BREAKING**] A number of changes have been made to the JSON report layout,
  including the moving, removal, and addition of sections. The layout is
  described in `schema.json`.
- [**BREAKING**] The minimum allowed/default value for `--min-adapter-overlap`
  was set to 1. In practice this has no effect, since length 0 alignments were
  never considered, but may break scripts running AdapterRemoval.
- [**BREAKING**] Drop support for raw error-rates to `--trim-mott-rate`, which
  was renamed to `--trim-mott-quality` to match other trimming options.
- [**BREAKING**] SAM/BAM output is now combined into a single file by default,
  including discarded reads. This can be overridden by setting the individual
  `--out-*` options.
- [**BREAKING**] Dropped `PG` tag from read-groups/records in SAM/BAM output.
- [**BREAKING**] Dropped (minimal) read-groups for SAM/BAM output. If desired,
  read-group information can be added with `--read-group`.
- [**BREAKING**] The `--report-duplication` option now supports k/m/g suffixes,
  and defaults to `100k` if used without an explicit value.
- [**BREAKING**] The `--read-group` option no longer attempts to unescape
  special characters. Instead, tags must be separated using embedded tabs
  (`--read-group $'ID:A\tSM:B'`) or provided as individual arguments
  (`--read-group 'ID:A' 'SM:B'`).
- Improved checks for conflicting command-line options.
- Barcodes are now recorded in FASTQ headers demultiplexing without trimming.
- The `$schema` URL is now included in the JSON report
- Makefile features are now enabled/disable with `true`/`false` instead of
  `yes`/`no`.
- Vega-lite is now loaded in the background, when opening the HTML reports,
  making the report readable before Vega-lite has loaded.
- Optimized alignments involving multiple possible adapter sequences, by only
  once performing the alignments that involve no adapter sequences.
- Optimized alignments involving multiple possible adapter sequences, by
  sorting the list of adapter sequences by hits. This increasing the odds that
  a good alignment is found early so that worse alignments can be skipped.
- The old Makefile was replaced with the Meson build system, but a wrapper
  Makefile is still provided/used as a convenience for setting the recommended
  build options.
- A number of small improvements were made to the `--help` text.
- Improved error messages when mismatching (paired) read names are detected.
- Singleton reads are now included in the overall summary statistics in
  JSON/HTML reports.
- Hardening flags are now enabled by default during compilation. This comes with
  a small performance cost, but most distros are also expected to enable similar
  flags by default.

### Fixed

- NA values were being written with '%' or 'bp' suffixes in HTML report.
- Some plots were omitted for merged reads in HTML report.
- Mate 2 adapters were reverse-complemented in JSON report when demultiplexing.
- SAM/BAM headers were not being written in demultiplexing mode.
- Mate 1/2 statistics were sampled independently, and thus potentially not
  derived from the same read pairs.
- The JSON/HTML reports would give different time-stamps for the run, since one
  gave the start time and the other the end time. Now start time is always used.
- Fixed failure when reading paired FASTQ files where read lengths differed
  between the two files
- Fixed report files and unidentified reads getting additional suffixes when
  filenames were specified manually during demultiplexing.
- Fixed `/dev/null` being listed as the path for some files when demultiplexing,
  and these outputs were disabled.
- Reverted the removal of support for '.' as equivalent to 'N' in FASTQ reads.
  This is found in some older data-sets (#112).
- Fixed misleading IO error messages, that would include descriptions of
  unrelated errors in some cases.

## [2.3.4] - 2024-08-24

This release adds a new couple of command-line options for handling non-ACGTN
bases in FASTQ data and back-ports a few minor fixes from the development
branch.

### Added

- Added support for converting Uracils (U) in input data to Thymine (T) via the
  `--convert-uracils` flag.
- Added support for replacing IUPAC-encoded degenerate bases with Ns via the
  `--mask-degenerate-bases` flag.
- Added DESTDIR support to `make install`.

### Fixed

- Improved progress timer accuracy, so updates occur closer to every 1M reads.

### Changed

- Minor improvements to `--help` text and documentation.

## [3.0.0-alpha2] - 2024-08-20

This is the second alpha release of AdapterRemoval v3. It is the intention that
a third alpha release, or the final 3.0 release, will follow within the next
couple of months.

As with alpha 1, changes that affect how AdapterRemoval is used (e.g. by
removing options) or that result in different output compared to AdapterRemoval
v2 are marked with the label "[**BREAKING**]".

In addition to changes listed below, this release includes increased throughput
thanks to improved parallelization of various steps in internal pipeline,
support for AVX512 and general improvements to the SIMD alignment algorithms,
loop unrolling of non-SIMD alignments to significantly increase throughput when
SIMD is not available, and a significant decrease in the number of allocations
to decrease overhead.

This release requires a compiler with support for c++17 and libdeflate is now a
mandatory dependency.

### Added

- Added support for converting Uracils (U) in input data to Thymine (T) via the
  `--convert-uracils` flag.
- Added support for replacing IUPAC-encoded degenerate bases with Ns via the
  `--mask-degenerate-bases` flag.
- Added support for writing output in SAM/BAM formats, with optional
  user-supplied read-group information.
- Added support for alignments using AVX512 instructions. AVX512 support only
  available when AdapterRemoval is compiled with GCC v11+ or Clang v8+.
- Added support selecting output file formats via the file extension and via
  the `--out-format` option. A corresponding option, `--stdout-format` was
  added to select the format for data written to STDOUT.
- Added support for reading from STDIN or writing to STDOUT when '-' is used as
  the filename, as an alternative to using `/dev/stdin` or `/dev/stdout`.
- Added dedicated threads solely for writing output data. This allows compute
  threads to work at full capacity, as long as the destination can consume
  written data fast enough. This may result in CPU utilization exceeding
  `--threads` by a couple of percent.
- Added support for setting DESTDIR when running `make install`.
- Added `--licenses` flag for displaying licenses of 3rd party code used by /
  incorporated into AdapterRemoval.
- Added `--simd` option allowing the user to select the specific SIMD
  instruction set they wish to use.
- Added `Containerfile` for building static binaries using alpine/musl.

### Changed

- [**BREAKING**] Changed the default `--mm`/`--mismatch-rate` from 1/3 to 1/6,
  in order to decrease the false positive rate, in particular for read merging.
- [**BREAKING**] Default to writing gzip-compressed FASTQ files; output written
  to STDOUT is uncompressed by default.
- [**BREAKING**] Discarded reads are no longer saved by default.
- [**BREAKING**] Output files for discarded reads and singleton (orphan)
  paired-end reads are only created if filtering is enabled.
- [**BREAKING**] The `--basename` / `--out-prefix` no longer defaults to
  `your_output`. Instead the user is required to set at least one `--out-*`
  option.
- [**BREAKING**] Merged `--identify-adapters` and `--report-only` commands. The
  adapter sequence is presently only reported in the HTML report, but will be
  added to the JSON report following some planned changes.
- [**BREAKING**] Reverted `--min-complexity` being enabled by default.
- Increased the default ``--threads`` value to 2.
- A number of command-line options were renamed for consistency; use of the old
  names is still supported, but will trigger a warning message.
- Re-organized compression: level 1 is streamed using isa-l, while levels 2-13
  correspond to libdeflate levels 1 to 12.
- Changed the default compression level to 5 on the new scale (libdeflate level
  4); this results in a ~40% increase in throughput at the cost of roughly ~3%
  larger output files.
- Setting an `--out-*` option in demultiplexing mode overrides the basename /
  prefix for that specific output type.
- Add smoothing to GC values calculated for the GC content curve, to account
  for the fact that possible GC% values are unevenly distributed depending on
  the read length.

### Removed

The following changes are all [**BREAKING**] as described above:

- Removed support for original merging algorithm has been removed. The
  `--merge-strategy additive` method produces very similar, but slightly more
  conservative scores.
- Removed the ability to randomly sample a base if no best base could be
  selected in case of mismatches. Such bases are now changed to `N`, while both
  methods assign a Phred score of 0 (`!`).

## [3.0.0-alpha1] - 2022-11-07

This is the first alpha release of AdapterRemoval v3. This is a major revision
of AdapterRemoval, with the goals of simplify usage by picking a sensible set of
default settings, adding new features to handle a wider range of data, providing
human/machine readable reports, and improving overall throughput.

This release features a number of breaking changes compared to AdapterRemoval v2
and it is therefore recommended that you carefully read the list of changes
below. Changes that affect how AdapterRemoval is used (e.g. by removing options)
or that result in different output compared to AdapterRemoval v2 are marked with
the label "[**BREAKING**]".

This is an alpha release; not all planned features are complete (more QC reports
are planned among other things), additional optimizations will be attempted, and
documentation is still needs to be expanded further before the final release.
Feedback is very welcome in the mean time.

### Added

- Reports are now available in JSON format for easy parsing and in HTML format
  for human consumptions. These replace the old `--settings` file.
- AVX2 enabled alignment algorithm for a significant performance boost (YMMV).
- Added support for detecting supported CPU extensions (SSE/AVX) at runtime.
- Support for combining output by simply by specifying the same filename for for
  multiple outputs types, e.g. `--output1 file.fq --output2 file.fq` will for
  example produce interleaved output.
- Added handling for `/dev/null` as a "magic" output filename. Read-types
  writing to this exact path will be discarded early in the pipeline, saving
  time previously spent processing, compressing, and writing FASTQ reads.
- Added read complexity filter inspired by [fastp].
- Added the ability to only processes the first `N` reads/read pairs via the
  newly added `--head N` command-line option.
- Added estimation of duplication rates based on the [FastQC] algorithm.
- Automatic detection of mate separators based on the first chunk of reads
  processed. The `--mate-separator` is therefore only required in cases where
  the results are ambiguous.
- Automatic gzip compression of output files with a `.gz` extension. This makes
  it possible to compress only a subset of files and removes the need for the
  `--gzip` option when manually specifying output files.
- Added options `--prefix-read1`, `--prefix-read2`, and `--prefix-merged` for
  adding custom prefixes to the names of FASTQ reads.
- Added support for trimming poly-X tails. Trimming can be done for a single
  nucleotide (e.g. poly-A) or for any combination of A, C, G, and T tails.

### Changed

- [**BREAKING**] Default adapters have been changed to the [recommended Illumina
  sequences], equivalent to the first 33 bp of the adapter sequences used by
  AdapterRemoval v2. This makes the default settings more generally applicable.
- [**BREAKING**] The trimming options `--trimwindows`, `--trimns`,
  `--trimqualities`, and `--minquality` have been deprecated in favor of a new
  the modified Mott's algorithm, which is enabled by default. The trimming
  algorithm used may be changed using the new `--trim-strategy` option.
- [**BREAKING**] Merging now defaults to using the conservative algorithm,
  meaning that matching quality scores are assigned `Q_match = max(Q_a, Q_b)`
  instead of `Q_match ~= Q_a + Q_b`, and that same-quality mismatches are
  assigned 'N' instead of one being picked at random. Motivated in part by
  `doi:10.1186/s12859-018-2579-2`. This can be changed using `--merge-strategy`.
- [**BREAKING**] The `--merge` option no longer has any effect when processing
  SE data; previously this option would treat reads with at
  `--minalignmentlength` adapter as pseudo-merged reads.
- [**BREAKING**] Merged reads are no longer given a `M_` name prefix and merged
  reads that have been trimmed after merging are no longer given an `MT_` name
  prefix. Instead, see the new option `--prefix-merged`.
- [**BREAKING**] Default filenames have all been revised and now include proper
  extensions to indicate the format.
- [**BREAKING**] The executable is now named `adapterremoval3`. This was done to
  allow v3 to coexist with AdapterRemoval v2 and to prevent accidental use of
  the wrong version.
- [**BREAKING**] Changed the default --maxns value from 1000 to "infinite"
- `--gzip` now defaults to compressing independent blocks of 64kb data using
  `libdeflate`. This significantly improves throughput in both single- and
  (especially) multi-threaded mode, but may be incompatible with a few programs.
  Compression levels of 3 and below use isa-l for compression and provides a
  more universally compatible output.
- The term "merging" is now used consistently instead of "collapsing", including
  for default output filenames. Options have been renamed, but old option names
  continue to work (except for `--outputcollapsedtruncated`).
- Improvements to alignment algorithm in order to terminate early if possible.
- Logging is now done more consistently and exposes options to increase or
  decrease the amount of messages printed (debug, info, warning, errors).

### Removed

The following changes are all [**BREAKING**] as described above:

- The `--outputcollapsedtruncated` has been removed and all merged reads
  (whether quality trimmed or not) are simply written to `--outputmerged`.
- The `--qualitybase-output` has been removed. Output is now always Phred+33.
- The `--combined-output` option has been removed in favor of allowing arbitrary
  merging of output files (see above).
- The `--settings` option has been replaced by `--out-json` and `--out-html` for
  machine and human readable reports, respectively.
- Removed support for guessing the intended command-line argument based on
  prefixes. I.e. `--th` will no longer be accepted for `--threads`. Due to the
  number of options added, removed, and renamed, this is no longer reliable.
- The deprecated `--pcr1` and `--pcr2` options have been removed.
- Dropped undocumented support for '.' as equivalent to 'N' in FASTQ reads.
- Support for reading and writing of bzip2 files has been removed.

## [2.3.3] - 2022-04-15

### Fixed

- Updated Catch2 to fix compilation with glibc 2.34, courtesy of loganrosen.

## [2.3.2] - 2021-03-17

### Added

- Added --collapse-conservatively merge algorithm based on FASTQ-join. See the
  man-page for more information

### Changed

- Improved error messages when AdapterRemoval failed to open or write FASTQ
  files (issue #42).
- Use 'install' in the Makefile; patch courtesy of Eric DEVEAUD.
- Added --collapse-deterministic to .settings file.

### Removed

- Removed benchmarking scripts which were included in the repo for the sake of
  making Schubert et al. 2016 reproducible. This is no longer relevant.

### Fixed

- Fixed build on some architectures. Patch courtesy of Andreas Tille/the Debian
  build team.
- Fixed display of max Phred scores in FASTQ validation error messages.
- Fixed --minadapteroverlap being misapplied in PE mode.

## [2.3.1] - 2019-06-23

### Added

- Added --preserve5p option. This option prevents AdapterRemoval from trimming
  the 5p of reads when the --trimqualities, --trimns, and --trimwindows options
  are used. Neither end of collapsed reads are trimmed when this option is used.

### Fixed

- Fixed Ns being miscounted as As when constructing consensus adapter sequences
  using --identify-adapters.

## [2.3.0] - 2019-03-12

### Added

- Added support for output files without a basename; to create such output
  files, use an empty basename (--basename "") or a basename ending with a slash
  (--basename path/).
- Added support for managing file handles to allow AdapterRemoval to run when
  the the number of output files exceeds the number of file handles, e.g. when
  demultiplexing large numbers of samples.

### Changed

- Reworked demultiplexing to improve performance for many paired barcodes.

### Fixed

- Fixed --collapse producing slightly different result on 32 bit and 64 bit
  architectures. Courtesy of Andreas Tille.

## [2.2.4] - 2019-02-10

### Fixed

- Fixed bug in --trim5p N which would AdapterRemoval to abort if N was greater
  than the pre-trimmed read length.

- Fixed --identify-adapters not respecting the --mate-separator option.

## [2.2.3] - 2019-01-22

### Added

- Added support for trimming reads by a fixed amount: --trim5p N --trim3p N.
  Different values may be given for each mate: --trim5p N1 N2. Trimming is
  carried out after adapters have been removed and reads have been collapsed, if
  enabled, but before quality trimming (Ns and low qualities).
- Added option for deterministic read merging (--collapse-deterministic). In
  this mode AdapterRemoval will set a merged base to 'N' with quality 0 if the
  corresponding bases on the two mates differ, and if both have the same quality
  score. The default behavior is to select one of the two bases at random.
- Added conda installation instructions, courtesy of Maxime Borry (maxibor).

### Fixed

- Fixed reporting of line numbers in error messages.
- Fixed reading mate 2 adapters specified via --adapter-list. Adapters would be
  used in the reverse orientation compared to --adapter2. Courtesy of Karolis
  (KarolisM).
- Fixed various typos and improved help/error messages.

## [2.2.2] - 2017-07-17

### Added

- Added support for Intel compilers, courtesy of Kevin Murray (kdmurray91).

### Changed

- Made gzip and bzip2 support mandatory.

## [2.2.1a] - 2017-05-17

### Fixed

- Fixed compilation on OSX.

## [2.2.1] - 2017-05-15

### Added

- Added support for specifying multiple filenames after --file1 and --file2, in
  which case the files are treated as if they were concatenated. This is
  supported for all operations. Special thanks to Stephen Clayton.
- Added additional run-time checks to catch race-conditions.
- Implemented quality trimming using a sliding window approach inspired by
  sickle (<https://github.com/najoshi/sickle>). Special thanks to Kevin Murray.

### Fixed

- Numerous spelling errors fixed courtesy of Andreas Tille.
- Progress messages written to STDERR no longer cause subsequent error messages
  to be written to the same line.
- Fixed miscounting of the total number of retained nucleotides, where mate 1
  reads were being counted twice instead of counting both mate 1 and mate 2.

## [2.2.0] - 2016-10-27

### Added

- Added support for reading FASTQ files using Windows-style newlines (\r\n).
- AdapterRemoval will now print a warning to STDERR if the same command-line
  option is specified multiple times.

### Changed

- AdapterRemoval now requires a C++11 compliant compiler; furthermore,
  multithreading is no longer an optional feature, as this is now implemented
  using the C++11 instead of directly calling pthreads.
- Add explicit message not to use the results after failed runs.
- Minor changes to .settings: Adapter numbers now 1-based; the 'Number of reads
  with adapters' is changed to 'Number of read pairs with adapters' when
  trimming PE reads; the 'Average read length of trimmed reads' is changed to
  'Average read length of retained reads' for clarity.
- Added support for combining output files; this is enabled using the
  --combined-output option, and ensures that all reads are written to the same
  file, or pair of files (for non-interleaved PE reads). The sequence o reads
  that fail are replaced with a single 'N' with quality score 0.

### Removed

- Dropped the undocumented 'poor' classification for alignments; for
  statistical purposes, reads are either counted as aligned or not aligned.
  This only changes how results are presented in the .settings files.

### Fixed

- Rework selection of nucleotides at overlapping positions with the same
  quality, in order to prevent potential data-races during tie-breaking, when
  running in multi-threaded mode.
- Reworked handling of barcodes to avoid unnecessary memory allocations, which
- would cause problems when using long barcodes.
- Fixed bug in the counting of singleton reads used in '.settings' files.
- Fixed mis-placement of underscore when pretty printing adapter sequences that
  included barcodes.
- Fixed misprinting of mate 2 adapter sequences in the .settings file; these
  would be printed in the reverse complemented orientation, relative to how
  they were specified on the command-line.

## [2.1.7] - 2016-03-11

### Added

- Improved 'make install'. This command now makes use of a PREFIX value to
  determine the installation destination (defaults to /usr/local), and includes
  the 'README.md' file and 'examples' folder in the installation.
- Improved 'make test'. This command now attempts to download the required
  testing library automatically, using either wget or curl if available.

### Changed

- The mate number is now stripped from collapsed reads, where previously this
  would always be '\1' (if set). However, if meta-data is present in the reads,
  that found in the mate 1 read is retained.
- The value used for --mate-separator is now written to the 'settings' file.

## [2.1.6] - 2016-03-03

### Added

- Added support for reading / writing interleaved FASTQ files; this is enabled
  by the options --interleaved-input and --interleaved-output, respectively, or
  by setting --interleaved option which implies both of the former options. See
  the README for an example.
- Expanded README with information about reading / writing FASTQ files with
  different PHRED encodings / maximum quality scores.

### Fixed

- Fixed bug in a sanity check meant to detect if the mate 1 and mate 2 files
  were of unequal length. This is now correctly detected in all cases.

## [2.1.5] - 2016-02-19

### Added

- Added the --mate-separator option, which specifies the character separating
  the mate number; by default this is '/', and AdapterRemoval will therefore
  identify mate numbers if read-names end with "/1" or "/2".

### Fixed

- Fixed race condition which could result in premature termination when using
  --gzip or --bzip2 together with the --threads options.
- Improved checks during compression and sanity checks following processing.

## [2.1.4] - 2016-02-09

### Changed

- Wrapped code in 'ar' namespace, and made it possible to compile
  AdapterRemoval as a static library (via the command 'make static'), allowing
  it to be used as part of other projects, courtesy of Hannes Pétur
  Eggertsson.
- Updated instructions for installing GTest library using new repository
  courtesy of Hannes Pétur Eggertsson.

### Fixed

- Fixed bug which could occasionally result in failure when bzip2 compression
  was enabled, by attempting to compress empty buffer.
- Fixed typos courtesy of Hannes Pétur Eggertsson

## [2.1.3] - 2015-12-25

### Added

- Added option --minadapteroverlap, which sets a minimum alignment length when
  carrying out trimming of single-end reads. The default (0) may result in an
  excess of false positives around (1 - 2 bp long), which may be mitigated by
  running AdapterRemoval with '--minadapteroverlap 3'.

### Changed

- Greatly expanded README.md, adding several examples with test data included in
  the 'examples' folder, demonstrating common usage of the program. ###
  Deprecated
- Updated man-page with missing information and rewrote several parts.
- Updated the help-text for several command-line options.

### Fixed

- Avoid writing information to stdout, so that (SE) trimming can be piped. This
  can be accomplished by using the option --output1 /dev/stdout.
- Fixed the --seed option, which was not properly applied during runtime.

## [2.1.2] - 2015-10-08

### Added

- Added setup instructions when running 'make test' for the first time.

### Changed

- Changed the way "full-length" and "truncated collapsed" reads are counted in
  the .settings file; previously, all collapsed reads were counted, even if
  these were subsequently discarded (due to the length). Now only retained reads
  are counted, matching the behavior of AdapterRemoval v1.x.

## [2.1.1] - 2015-09-14

### Fixed

- Fixed broken assert preventing the use of --adapter-list.
- Fixed bug using --qualitybase-output for both input and output.

## [2.1.0] - 2015-09-08

### Added

- Added support for (transparently) reading and writing bzip2 files.
- Added parallelization of adapter trimming and identification using pthreads;
  the number of threads used is specified using --threads. Note that only one
  thread is allowed to perform IO (reads / writes) at a time, to prevent
  clobbering the disk, but compression (if enabled) is performed in parallel.
- Added support for combined demultiplexing and adapter-removal using the
  --barcode-list command-line option; when demultiplexing, the barcodes
  identified for a given read is added to the adapter sequence, in order to
  ensure correct trimming of the reads.
- Features depending on external libraries (gzip, bzip2, and threading support)
  can now be disabled in the Makefile on systems lacking these libraries.
- Display currently specified --adapter1 / adapter2 sequences for comparison
  when attempting to infer adapter sequences. Only the first pair is used, if
  multiple adapter pairs are specified.

### Changed

- Sites with no majority-base during adapter-identification are set to N.
- FASTQ reads with Solexa scores are now output as Solexa scores by default,
  rather than Phred+64. Note that the program represents quality scores using
  Phred scores internally, resulting in a lossy conversion. It is therefore
  recommended to convert to Phred scores rather than use Solexa scores.
- Progress report now shows total number of reads processes, for both single
  ended and pair ended analyses.

### Fixed

- Fixed failure to read of barcode sequences (--5prime / --5prime-list).

## [2.0.0] - 2014-03-10

Version 2.0.0 of AdapterRemoval is a near complete rewrite, with the goal of
improved safety, increased speed, fixing a number of minor issues with previous
versions of AdapterRemoval, and adding a few new features.

### Breaking changes

- Command-line arguments --pcr1 and --pcr2 have been deprecated in favor of
  --adapter1 and --adapter2. While --pcr1 and --adapter1 are equivalent,
  --adapter2 expects the adapter sequence which may be observed in raw mate 2
  reads, unlike --pcr2 which expected the sequence which could be observed in
  the reverse complement of mate 2 reads (cf. the README).
- The use of --file1 and (optionally) --file2 is now required; reads will not be
  read from STDIN, nor written to STDOUT by default. To approximate the previous
  behavior, the following command may be used: $ AdapterRemoval --file1
  /dev/stdin --output1 /dev/stdout
- Per-read statistics of adapter / low-quality base trimming using --stats is no
  longer supported.
- Support for multiple adapter sequences as well as multiple barcode sequences;
  AdapterRemoval will favor the highest scoring alignment, favoring longer
  alignments over shorter alignments with the same score, and favoring
  alignments with the fewest ambiguous bases (N) involved if the score and
  length is identical.

### Added

- Limited support for Solexa quality scores; these are converted to and saved as
  Phred+33 or Phred+64 encoded scores.
- Added the ability to identify adapter sequences for paired-ended reads, by
  identifying reads which extends past the ends of the template sequence, and
  extracting the adapters from these.
- Added support for reading / writing gzipped compressed FASTQ files; if enabled
  (using the --gzip flag), the ".gz" extension is added to filenames, unless the
  filenames are explicitly specified by the user.
- It is now possible to explicitly specify the RNG seed, to allow individual
  runs to be reproduced; the seed is also written to the .settings file.
- An (optional) progress report is printed during usage, indicating the run-time
  and number of reads processed.

### Changed

- Strict validation of input FASTQ records, to ensure that records are well
  formed, that quality scores fall within the expected range given the specified
  format/offset, and more.
- Improved handling of asymmetric read-pairs, in which the length of the mate 1
  read differs from the length of the mate 2 read.
- Significant improvements in performance, resulting in a ~5x increase in the
  rate of adapter trimming in basic version, and a ~20x increase in the rate of
  adapter trimming in the SSE enabled version (the default).
- If --collapse is set in single-ended mode, "collapsed" reads will be
  identified using the same criteria as for paired-ended mode, i.e. requiring
  that at least --minalignmentlen bases overlap, and written to .collapsed and
  .collapsed.truncated. This allows for the identification of reads that are
  complete inserts.
- Length distributions are now calculated per read-type post-trimming (mate 1,
  mate 2, collapsed, etc.) and written to the .settings file.
- Barcodes may now contain Ns.
- Replaced use of lower bits of rand() calls with random(), as the former
  generates low entropy bits in that range on some (non-Linux) platforms.
- Seed is now initialized using a mix of seconds and microseconds, instead of
  the current time in seconds, to reduce the risk of multiple instances spawned
  within a short timespan from using the same seed.

### Fixed

- Fixed underestimation of error-probabilities during sequence collapse.
- Fixed (further) underestimation of error-probabilities of bases during
  collapsing, for conflicting base-calls with the same Phred score.
- Fixed the maximum number of mismatches for alignments in the range of 6 .. 9
  bases always being 1, even if --mm was set to 0.
- Fixed the maximum number of mismatches for alignments being calculated based
  on the length of the alignment including ambiguous bases (N), thereby
  inflating the number of mismatches allowed for poor alignments.
- Fixed well-aligned reads being discarded due to the minimum-length requirement
  after trimming not being counted as well-aligned, resulting in the total
  number of alignments not matching the total number of reads.
- Fixed bug in shifts for PE reads, which was causing some alignments to be
  missed for adapter-only (i.e. no insert sequence) sequences.
- Improved input validation and sanity checks for command-line parameters.

## [1.5.4] - 2014-04-23

### Changed

- Reduced the amount of IO operations during trimming.

### Fixed

- Fixed bug in which collapsed reads would not be considered truncated if bases
  were trimmed from the 5' end.
- Fixed bug in which the quality bases used for mate 2 during collapsing of
  overlapping read pairs made use of quality scores with a wrong orientation.

## [1.5.2] - 2013-10-22

### Changed

- Added a reference to the paper to both the man page and the help text.

### Fixed

- Fixed a minor bug in the collapse code where two very low quality bases might
  give rise to a third low quality base being called. For example, a C with
  quality " and a T with quality ! would result in an A with quality #. This has
  been fixed so that the the result is now C with quality ".

## [1.5.0] - 2013-04-29

### Changed

- Collapsed pairs are now written to two files: One contains full-length
  collapsed pairs constituting the full insert, the other contains collapsed
  pairs that have been truncated due to low qualities or Ns in the reads.

## [1.4.0] - 2013-03-24

### Added

- Support the use of '.' instead of 'N' to encode undefined nucleotides.
- Some minor changes to output etc.

### Fixed

- There was a typo in the adapter sequence used for PCR2!

## [1.3.0] - 2013-02-10

I have updated AdapterRemoval and released version 1.3. These changes are based
on feedback from users of the program that had some very specific and well-
founded suggestions. Some of these changes are minor, others will have more
dramatic effects on the use of the program so please read these notes carefully

### Changed

- The PCR1 and PCR2 sequences are used as-is (not reverse-complementation). You
  have to make sure that the sequences you search for are correct.
- The default PCR1 and PCR2 sequences are now:

- PCR1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
- PCR2: AATGATACGGCGACCACCGAGATCACACTCTTTCCCTACACGACGCTCTTCCGATCT

- Use of PCR1 and PCR2 changed to make the program consistent: Now, you always
  search for the sequence PCR1 in READ1 (whether single end or paired end), and
  you search for PCR2 in READ2. In single end mode, this corresponds to having
  an empty READ2 and ignore PCR2 as illustrated below:

- For paired end data, PCR2-READ1 aligned to READ2-PCR1.
- For single end data, READ1 aligned to PCR1.

- Collapsed reads are now names "@M\_...".
- Collapsed reads are put in a separate file with extension ".collapsed".

### Fixed

- Fixed an occasional segmentation fault.

## 1.1.0 - 2012-05-01

### Added

- It is now possible to look for adapter in the 5' end of reads using the
  --5prime parameter.
- Added option for discarding reads with too many gaps using --maxns max.

### Changed

- Updated trimming of qualities.
- The program handles lower vs upper case issues by translating all sequences
  to upper case.
- The program now checks for inconsistent parameters.

### Fixed

- Fixed some typographical issues with output.

## 1.0.0

- Initial release

[3.0.0-alpha3]: https://github.com/MikkelSchubert/adapterremoval/compare/3.0.0-alpha2...3.0.0-alpha3
[2.3.4]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.3.3...v2.3.4
[3.0.0-alpha2]: https://github.com/MikkelSchubert/adapterremoval/compare/3.0.0-alpha1...3.0.0-alpha2
[3.0.0-alpha1]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.3.3...3.0.0-alpha1
[2.3.3]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.3.2...v2.3.3
[2.3.2]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.3.1...v2.3.2
[2.3.1]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.3.0...v2.3.1
[2.3.0]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.2.4...v2.3.0
[2.2.4]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.2.3...v2.2.4
[2.2.3]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.2.2...v2.2.3
[2.2.2]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.2.1...v2.2.2a
[2.2.1a]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.2.1...v2.2.1a
[2.2.1]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.2.0...v2.2.1
[2.2.0]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.1.7...v2.2.0
[2.1.7]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.1.6...v2.1.7
[2.1.6]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.1.5...v2.1.6
[2.1.5]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.1.4...v2.1.5
[2.1.4]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.1.3...v2.1.4
[2.1.3]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.1.2...v2.1.3
[2.1.2]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.1.1...v2.1.2
[2.1.1]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.1.0...v2.1.1
[2.1.0]: https://github.com/MikkelSchubert/adapterremoval/compare/v2.0.0...v2.1.0
[2.0.0]: https://github.com/MikkelSchubert/adapterremoval/compare/v1.5.4...v2.0.0
[1.5.4]: https://github.com/MikkelSchubert/adapterremoval/compare/v1.5.2...v1.5.4
[1.5.2]: https://github.com/MikkelSchubert/adapterremoval/compare/v1.5.0...v1.5.2
[1.5.0]: https://github.com/MikkelSchubert/adapterremoval/compare/v1.4.0...v1.5.0
[1.4.0]: https://github.com/MikkelSchubert/adapterremoval/compare/v1.3.0...v1.4.0
[1.3.0]: https://github.com/MikkelSchubert/adapterremoval/compare/v1.1.0...v1.3.0
[fastp]: (https://github.com/OpenGene/fastp/)
[fastqc]: (https://github.com/s-andrews/FastQC)
[recommended illumina sequences]: https://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html

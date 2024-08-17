Input and output
================

Input files
-----------

FASTQ files
^^^^^^^^^^^

AdapterRemoval v3 expected standard FASTQ files as input specified using the `--in-file1` and `--in-file2` options. DNA sequences and quality scores are expected to be written on a single line, and DNA sequences are expected to contain only standard bases A, C, G, T, and N (ambiguous bases). AdapterRemoval v3 supports reading of FASTQ files with Phred+33, Phred+64, and Solexa encoded `quality scores`_. Input may optionally be gzip compressed.

When specifying paired reads, the FASTQ record names are expected to be identical excepting for the last two characters. The last two characters typically indicate the mate number, and consist of a matching separator character followed by either `1` or `2`. Examples include `/1` and `/2`, and `.1` and `.2`.


Table of adapters
-----------------

TODO

Table of barcodes
-----------------

TODO

Output files
------------

The default filenames for single end mode are as follows:

 * `{prefix}[.{sample}].fastq`
 * `{prefix}[.{sample}].discarded.fastq`
 * `{prefix}.unidentified.fastq`, if demultiplexing is enabled

The default filenames for paired end mode are as follows:

 * `{prefix}[.{sample}].r1.fastq`
 * `{prefix}[.{sample}].r2.fastq`
 * `{prefix}[.{sample}].merged.fastq`, if merging is enabled
 * `{prefix}[.{sample}].discarded.fastq`
 * `{prefix}[.{sample}].singleton.fastq`
 * `{prefix}.unidentified.r1.fastq`, if demultiplexing is enabled
 * `{prefix}.unidentified.r2.fastq`, if demultiplexing is enabled

The `{sample}` field is only added if demultiplexing is enabled.


FASTQ files
^^^^^^^^^^^

Output files with the `.fastq` extension correspond to that described for input FASTQ files, but are always Phred+33 encoded regardless of the encoding of the input data.


HTML report
^^^^^^^^^^^

TODO

JSON report
^^^^^^^^^^^

TODO


.. _quality scores: https://en.wikipedia.org/wiki/FASTQ_format#Quality

##################
 Input and output
##################

*************
 Input files
*************

FASTQ files
===========

AdapterRemoval v3 expected standard FASTQ files as input specified using the ``--in-file1`` and ``--in-file2`` options. DNA sequences and quality scores are expected to be written on a single line, and DNA sequences are expected to contain only standard bases A, C, G, T, and N (ambiguous bases). AdapterRemoval v3 supports reading of FASTQ files with Phred+33, Phred+64, and Solexa encoded `quality scores`_. Input may optionally be gzip compressed.

When specifying paired reads, the FASTQ record names are expected to be identical excepting for the last two characters. The last two characters typically indicate the mate number, and consist of a matching separator character followed by either ``1`` or ``2``. Examples include ``/1`` and ``/2``, and ``.1`` and ``.2``.

*******************
 Table of adapters
*******************

The table used with ``--adapter-table`` is expected to be a plain-text table containing one or two whitespace separated columns, with each line containing the same number of columns:

.. code::

   AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
   AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA   AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG

Column one corresponds to ``--adapter1`` and column two corresponds to ``--adapter2``. Only one column is required in single end mode, while both columns are required in paired end mode.

*******************
 Table of barcodes
*******************

The table used with ``--barcode-table`` is expected to be a plain-text table containing two or three whitespace separated columns, with each line containing the same number of columns:

.. code::

    sample_1 ATGCGGA TGAATCT
    sample_2 ATGGATT ATAGTGA
    sample_7 CAAAACT TCGCTGC

Sample names in the first column may only contain letters, numbers and underscores ``_``. The barcodes in the second column is expected to be found at the 5' of mate 1 reads and the barcode in the third column is expected to be found at the 5' of mate 2 reads.

If ``--barcode-orientation explicit`` is used, the table is expected to contain a fourth column specifying the orientation of each barcode pair:

.. code:: text

    sample_1 ATGCGGA TGAATCT forward
    sample_1 TGAATCT ATGCGGA reverse
    sample_2 ATGGATT ATAGTGA reverse
    sample_7 CAAAACT TCGCTGC forward

The following names for orientations are supported: ``forward``, ``reverse``, ``fwd``, ``rev``, ``+``, and ``-``. This table is used as is, and the "missing" orientations are not automatically generated.

**************
 Output files
**************

The default filenames for single end mode, when using ``--out-prefix`` are as follows:

   -  ``{prefix}[.{sample}].r1.fastq.gz``
   -  ``{prefix}[.{sample}].discarded.fastq.gz``
   -  ``{prefix}.unidentified.fastq.gz``, if demultiplexing is enabled

The default filenames for paired end mode, when using ``--out-prefix``, are as follows:

   -  ``{prefix}[.{sample}].r1.fastq.gz``
   -  ``{prefix}[.{sample}].r2.fastq.gz``
   -  ``{prefix}[.{sample}].merged.fastq.gz``, if merging is enabled
   -  ``{prefix}[.{sample}].singleton.fastq.gz``
   -  ``{prefix}.unidentified.r1.fastq.gz``, if demultiplexing is enabled
   -  ``{prefix}.unidentified.r2.fastq.gz``, if demultiplexing is enabled

Discarded reads are not written by default, and therefore does not have a default filename. The ``{sample}`` field is only added if demultiplexing is enabled.

For SAM and BAM output, read types are classified using SAM/BAM flags rather than by filename, and all reads are therefore written to a single file, depending on the chosen format:

   -  ``{prefix}[.{sample}].sam``
   -  ``{prefix}[.{sample}].sam.gz``
   -  ``{prefix}[.{sample}].bam``

Additionally, unlike FASTQ data, discarded reads are included by default, though this can be disabled by using ``--out-discarded /dev/null``.

FASTQ files
===========

Output files with the ``.fastq`` extension correspond to that described for input FASTQ files, but are always Phred+33 encoded regardless of the encoding of the input data.

SAM/BAM files
=============

TODO

HTML report
===========

TODO

Note that this file currently requires an internet connection, due to depending on external JavaScript files for interactive plotting.

JSON report
===========

TODO

A full `JSON schema`_ can be found in the ``schema.json`` file in the root of the source distribution of AdapterRemoval.

.. _json schema: https://json-schema.org/

.. _quality scores: https://en.wikipedia.org/wiki/FASTQ_format#Quality

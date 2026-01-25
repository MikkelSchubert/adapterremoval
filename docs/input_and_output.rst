##################
 Input and output
##################

*************
 Input files
*************

FASTQ files
===========

AdapterRemoval v3 expects standard FASTQ files as input specified using the ``--in-file1`` and ``--in-file2`` options. DNA sequences and quality scores are expected to be written on a single line, and DNA sequences are expected to contain only standard bases A, C, G, T, and N (ambiguous bases). AdapterRemoval v3 supports reading of FASTQ files with Phred+33, Phred+64, and Solexa encoded `quality scores`_. Input may optionally be gzip compressed.

When specifying paired reads, the FASTQ record names are expected to be identical excepting for the last two characters. The last two characters typically indicate the mate number, and consist of a matching separator character followed by either ``1`` or ``2``. Examples include ``/1`` and ``/2``, and ``.1`` and ``.2``.

*******************
 Table of adapters
*******************

The table used with ``--adapter-table`` is expected to be a plain-text table containing one or two whitespace-separated columns, with each line containing the same number of columns:

.. code-block:: text

    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA   AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG

Column one corresponds to ``--adapter1`` and column two corresponds to ``--adapter2``. Only one column is required in single-end mode, while both columns are required in paired-end mode.

*******************
 Table of barcodes
*******************

The table used with ``--barcode-table`` is expected to be a plain-text table containing two or three whitespace-separated columns, with each line containing the same number of columns:

.. code-block:: text

    sample_1 ATGCGGA TGAATCT
    sample_2 ATGGATT ATAGTGA
    sample_7 CAAAACT TCGCTGC

Sample names in the first column may only contain letters, numbers and underscores ``_``. The barcodes in the second column are expected to be found at the 5' of mate 1 reads and the barcode in the third column is expected to be found at the 5' of mate 2 reads.

If ``--barcode-orientation explicit`` is used, the table is expected to contain a fourth column specifying the orientation of each barcode pair:

.. code-block:: text

    sample_1 ATGCGGA TGAATCT forward
    sample_1 TGAATCT ATGCGGA reverse
    sample_2 ATGGATT ATAGTGA reverse
    sample_7 CAAAACT TCGCTGC forward

The following names for orientations are supported: ``forward``, ``reverse``, ``fwd``, ``rev``, ``+``, and ``-``. This table is used as is, and the "missing" orientations are not automatically generated.

**************
 Output files
**************

During a normal run, AdapterRemoval will generate a combination of trimmed reads and report files. For a description of how this is controlled, see :doc:`detailed_overview`. In the following examples, the value ``${prefix}`` is used to represent the path set using ``--out-prefix``.

FASTQ output
============

When not performing demultiplexing, AdapterRemoval will generate the following files:

- ``${prefix}.r1.fastq.gz``
- ``${prefix}.r2.fastq.gz`` (paired-end mode)
- ``${prefix}.singleton.fastq.gz`` (paired-end mode)
- ``${prefix}.merged.fastq.gz`` (paired-end mode, if read merging is enabled)

If demultiplexing is enabled, the following files will be generated with ``--out-prefix``, with the key ``${sample}`` additionally replaced with each sample name:

- ``${prefix}.unidentified.r1.fastq.gz``
- ``${prefix}.unidentified.r2.fastq.gz`` (paired-end mode)
- ``${prefix}.${sample}.r1.fastq.gz``
- ``${prefix}.${sample}.r2.fastq.gz`` (paired-end mode)
- ``${prefix}.${sample}.singleton.fastq.gz`` (paired-end mode)
- ``${prefix}.${sample}.merged.fastq.gz`` (paired-end mode, if read merging is enabled)

The ``unidentified`` files contain any reads for which no sample could be identified using the supplied barcode sequences. Note that no trimming, filtering, or merging is performed if ``--demultiplex-only`` is used.

SAM/BAM files
=============

When using ``--out-format`` with ``.sam``, ``.sam.gz``, or ``.bam``, AdapterRemoval defaults to writing a single file containing *all* output, including discarded reads. For example, for BAM output:

- ``${prefix}.bam``

Or when demultiplexing is enabled, replacing the key ``${sample}`` with each sample name:

- ``${prefix}.json``
- ``${prefix}.html``
- ``${prefix}.unidentified.bam``
- ``${prefix}.${sample}.bam``

This can be overridden as described above, using the individual ``--out`` options. To exclude discarded reads from the output, use ``--out-discarded /dev/null``.

SAM/BAM tags
------------

When demultiplexing, AdapterRemoval will assign read-groups based on the barcodes, using the provided sample names. For samples with multiple barcodes, a read-group will be created for each barcode or barcode pair. When ``--barcode-orientation`` is used, AdapterRemoval will additionally record the orientation of the barcodes using ``or`` tags in the corresponding ``@RG`` headers.

Additional read-group information can be set using ``--read-group``, which also adds read-group information when not demultiplexing.

HTML report
===========

When using the ``--out-html`` option or when using ``--out-prefix``, AdapterRemoval writes a single-file HTML report using the default filename ``${prefix}.html``.

An example of the HTML report generated by AdapterRemoval can be found here: `Example HTML report`_

.. warning::

    The HTML report currently requires an active internet connection, due to depending on external JavaScript files for drawing figures.

JSON report
===========

When using the ``--out-json`` option or when using ``--out-prefix``, AdapterRemoval writes a JSON report using the default filename ``${prefix}.json``.

An example of the JSON report generated by AdapterRemoval can be found here: `Example JSON report`_. A full `JSON schema`_ defining the content of the report can be found here: `Report JSON schema`_

.. _example html report: https://mikkelschubert.github.io/adapterremoval/examples/3.x.html

.. _example json report: https://mikkelschubert.github.io/adapterremoval/examples/3.x.json

.. _json schema: https://json-schema.org/

.. _quality scores: https://en.wikipedia.org/wiki/FASTQ_format#Quality

.. _report json schema: https://mikkelschubert.github.io/adapterremoval/schemas/3.0.0.json

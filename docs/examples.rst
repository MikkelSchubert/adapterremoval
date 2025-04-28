.. highlight:: Bash

###############
 Example usage
###############

The following examples all make use of the data included in the 'examples' folder.

*************
 Basic usage
*************

To trim single-end Illumina TruSeq or similar data, simplify specify one or more files using the --in-file1 option:

.. code::

   adapterremoval3 \
       --in-file1 reads_1.fastq \
       --out-prefix output_single

To trim paired-end FASTQ reads, instead specify two input files using the ``--in-file1`` and ``--in-file2`` options:

.. code::

   adapterremoval3 \
       --in-file1 reads_1.fastq \
       --in-file2 reads_2.fastq \
       --out-prefix output_paired \
       --merge

The ``--out-prefix`` option ensure that all output files start with the specified prefix. For the paired reads, the ``--merge`` option enables merging of reads overlapping at least 11bp (by default).

For non-Illumina/non-TruSeq data, it is important to :ref:`specify the correct adapters <specifying_adapters>` when running AdapterRemoval.

****************
 Output formats
****************

AdapterRemoval support writing processed reads as uncompressed FASTQ reads (``.fastq`` or ``.fq``), as gzip compressed FASTQ reads (``.fastq.gz`` or ``.fq.gz``), as uncompressed or compressed SAM records (``.sam`` and ``.sam.gz``), and as compressed BAM records (``.bam``) or uncompressed BAM records (``ubam``).

Formats are specified either using the ``--out-format`` option in conjunction with the ``-out-prefix`` option, or by specifying one of the above extensions when using one of the ``--out-*`` options (except for ``ubam``). If AdapterRemoval does not recognize the extension, then the format specified via ``--out-format`` is used, defaulting to gzip compressed FASTQ.

By default, the uncompressed version of the format specified via ``--out-format`` is used for records written to STDOUT, under the assumption that the output is to be processed by another program, but this may be overridden using the ``--stdout-format`` option.

************************************************************************
 Standard input (STDIN), standard output (STDOUT), and disabling output
************************************************************************

Reading from STDIN and writing to STDOUT can be accomplished *either* by using the special ``/dev/stdin`` and ``/dev/stdout`` files, *or* by using the filename ``-`` for either ``--in-file*`` or ``--out-*`` options:

.. code::

   some-command | adapterremoval3 --in-file1 - --out-file1 - | some-other-command

Input from STDIN and output to STDOUT can freely be interleaved as described in the *Interleaved input and output* below.

Finally, AdapterRemoval recognizes the special path ``/dev/null`` and will skip writing reads if this path is specified. That means that you can use the ``--prefix`` option to specify a default output path, and disable the read types that you do not care about. Meaning that instead of

.. code::

   adapterremoval3 \
       --in-file1 reads_1.fastq \
       --in-file2 reads_2.fastq \
       --out-html output.html \
       --out-json output.json \
       --out-file1 output.r1.fastq.gz \
       --out-file2 output.r2.fastq.gz \
       --out-merged output.merged.fastq.gz

you could write

.. code::

   adapterremoval3
      --in-file1 reads_1.fastq \
      --in-file2 reads_2.fastq \
      --out-prefix output \
      --out-singleton /dev/null

****************************
 Multiple input FASTQ files
****************************

More than one input file may be specified listed after the ``--in-file1`` and ``--in-file2`` options. Files are processed in the specified order, as if they had been concatenated using ``cat`` or ``zcat``:

.. code::

   adapterremoval3 \
       --in-file1 reads_1a.fastq reads_1b.fastq reads_1c.fastq
   adapterremoval3 \
       --in-file1 reads_1a.fastq reads_1b.fastq reads_1c.fastq \
       --in-file2 reads_2a.fastq reads_2b.fastq reads_2c.fastq

******************************
 Interleaved input and output
******************************

AdapterRemoval is able to read and write paired-end reads stored in a single, so-called interleaved FASTQ file (one pair at a time, first mate 1, then mate 2). This is accomplished by specifying the location of the file using ``--in-file1`` and *also* setting the ``--interleaved`` command-line option:

.. code::

   adapterremoval3 \
       --interleaved \
       --in-file1 interleaved.fastq \
       --out-prefix output_interleaved

Other than taking just a single input file, this mode operates almost exactly like paired end trimming (as described above); the mode differs only in that paired reads are not written to a 'r1' and a 'r2' file, but instead these are instead written to a single file. The location of this file is controlled using the ``--out-file1`` option.

Enabling either reading or writing of interleaved FASTQ files, both not both, can be accomplished by specifying either of the ``--interleaved-input`` and ``--interleaved-output`` options, both of which are enabled by the ``--interleaved`` option.

Alternatively, you can specify the same output file for multiple output types, in order to write all of those reads to a single file in interleaved mode:

.. code::

   adapterremoval3 \
       --in-file1 input_1.fastq.gz \
       --in-file2 input_2.fastq.gz \
       --out-file1 output_interleaved.fastq.gz \
       --out-file2 output_interleaved.fastq.gz

The ability to interleave output extends to all output types, except for the two reports (``--out-json`` and ``--out-html``), and one could for example write both discarded and singleton reads to the same file (``output_interleaved.discarded.fastq.gz``) using the following command:

.. code::

   adapterremoval3 \
        --in-file1 input_1.fastq.gz \
        --in-file2 input_2.fastq.gz \
        --out-prefix output_interleaved \
        --out-discarded output_interleaved.discarded.fastq.gz \
        --out-singleton output_interleaved.discarded.fastq.gz

***********************************
 Different quality score encodings
***********************************

By default, AdapterRemoval expects the quality scores in FASTQ reads to be Phred+33 encoded, meaning that the error probabilities are encoded as ``(char)('!' - 10 * log10(p))``. Most data will be encoded using Phred+33, but Phred+64 and 'Solexa' encoded quality scores are also supported. These are selected by specifying the ``--quality-format`` command-line option (specifying either '33', '64', or 'solexa'):

.. code::

   adapterremoval3 \
       --quality-format 64 \
       --in-file1 reads_q64.fastq \
       --out-prefix output_phred_64

Output is always saved as Phred+33. See `this Wikipedia article`_ for a detailed overview of Phred encoding schemes currently and previously in use.

*******************************************************
 Trimming paired-end reads with multiple adapter pairs
*******************************************************

It is possible to provide multiple, different sets of adapters for trimming, in which case AdapterRemoval will select the best match for each read (pair).  This is done by providing a one or two-column table, for SE and PE trimming, respectively, with each line containing one or two adapters, separated by whitespace.

For example, to specify both Illumina TruSeq and BIGSeq adapters, one might save the following in the file ``adapters.txt``:

.. code:: text

   AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
   AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA   AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG

This file is then specified using the ``--adapter-list`` option:

.. code::

   adapterremoval3 \
       --in-file1 reads_1.fastq \
       --in-file2 reads_2.fastq \
       --out-prefix output_multi \
       --adapter-list adapters.txt

Pairs of adapters are used exactly as written, and the resulting QC reports lists how frequently adapter or pairs of adapters were used.

Note that throughput decreases proportionally to the number of adapters, and it is therefore *not* recommended to use this functionality unless strictly necessary. When adapters differ only after the first N bases (for 20-30 bp), for example due to an embedded barcode, then is typically better to specify only the shared part of the adapter sequences on the command line.

It is also possibly to mask variable sites in an adapter sequence, such as barcodes, by setting sites to ``N``. These sites will then not be considered part of the alignments.

****************
 Demultiplexing
****************

AdapterRemoval supports simultaneous demultiplexing and adapter trimming; demultiplexing is carried out using a simple comparison between the specified barcode (a sequence of A, C, G, and T) and the first N bases of the mate 1 read, where N is the length of the barcode. Demultiplexing of double-indexed sequences is also supported, in which case two barcodes must be specified for each sample. The first barcode is then compared to first ``N_1`` bases of the mate 1 read, and the second barcode is compared to the first ``N_2`` bases of the mate 2 read. By default, this comparison requires a perfect match. Reads identified as containing a specific barcode(s) are then trimmed using adapter sequences including the barcode(s) as necessary. Reads for which no (pair of) barcodes matched are written to a separate file or a pair of files (for paired end reads).

Demultiplexing is enabled by creating a table of barcodes, the first column of which species the sample name (using characters a-z, A-Z, 0-9, or _) and the second and (optional) third columns specifies the barcode sequences expected at the 5' termini of mate 1 and mate 2 reads, respectively.

For example, a table of barcodes from a double-indexed run might be as follows (see examples/barcodes.txt):

.. code::

   cat barcodes.txt
   sample_1 ATGCGGA TGAATCT
   sample_2 ATGGATT ATAGTGA
   sample_7 CAAAACT TCGCTGC

AdapterRemoval is invoked with the ``--barcode-list`` option, specifying the path to this table:

.. code::

   adapterremoval3 --in-file1 demux_1.fastq --in-file2 demux_2.fastq --out-prefix output_demux --barcode-list barcodes.txt

This generates a set of output files for each sample specified in the barcode table, using ``output_demux`` as the prefix for output filenames, followed by a dot and the sample name, followed by a dot and the default name for a given file type. The reports generated by AdapterRemoval contains information about the number of reads identified for each sample and (in the JSON file) detailed per-sample quality metrics.

The maximum number of mismatches allowed when comparing barcodes is controlled using the options ``--barcode-mm``, ``--barcode-mm-r1``, and ``--barcode-mm-r2``, which specify the maximum number of mismatches total, and the maximum number of mismatches for the mate 1 and mate 2 barcodes respectively. Thus, if mm_1(i) and mm_2(i) represents the number of mismatches observed for barcode-pair i for a given pair of reads, these options require that

   #. mm_1(i) <= ``--barcode-mm-r1``
   #. mm_2(i) <= ``--barcode-mm-r2``
   #. mm_1(i) + mm_2(i) <= ``--barcode-mm``

If the ``--demultiplex-only`` option is used, then no trimming/processing is performed after the demultiplexing step:

.. code::

   adapterremoval3 --in-file1 demux_1.fastq --in-file2 demux_2.fastq --out-prefix output_only_demux --barcode-list barcodes.txt --demultiplex-only

These reads will still contain adapters, and for paired reads/double indexed data these adapters will be prefixed by the barcode sequence(s). The adapter plus barcode sequences are reported for each sample in the JSON report file.

***************************************************
 Quality reports and identifying adapter sequences
***************************************************

AdapterRemoval generates a detailed report of input and output data, as part of its operation. This report can additionally be run without performing read processing, meaning that statistics are only provided for the "raw" input data.

Additionally, when run without read processing enable,d AdapterRemoval attempts to infer a consensus adapter sequence. This is done based on fragments identified as belonging to the adapters through pairwise alignments of the reads.

By default, quality reports are written to ``prefix.html`` and ``prefix.json``, if the ``--out-prefix`` option is used. The files can either or both be written to specific files using ``--out-json`` and ``--out-html``:

This means that we can simply omit other output files to generate only the reports:

.. code::

   adapterremoval3 \
       --in-file1 reads_1.fastq \
       --in-file2 reads_2.fastq \
       --out-json my_report.json \
       --out-html my_report.html

To generate a report only for the input data, including inference of adapter sequences, use the ``--report-only`` mode. Since only the reports are generated in this mode, we can use the ``--out-prefix`` option to simplify the command:

.. code::

   adapterremoval3 \
       --report-only \
       --in-file1 reads_1.fastq \
       --in-file2 reads_2.fastq \
       --out-prefix my_report

The consensus sequences inferred are compared to those specified using the ``--adapter1`` and ``--adapter2`` command-line options, or with the default values for these if no values have been given (as in this case). Pipes (|) indicate matches between the provided sequences and the consensus sequence, and "*" indicate the presence of unspecified bases (Ns).

Best practice is to compare the consensus with published Illumina_ or `BGI/MGI`_ adapter sequences and pick out the best matches. However, on occasion there will be differences between the published sequences and the observed adapter sequences.

When using the consensus directly, it is not recommended to use the full consensus sequence, since the quality of the consensus sequence declines quickly towards the 3' end.

.. _bgi/mgi: https://en.mgitech.cn/Download/download_file/id/71

.. _illumina: https://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html

.. _this wikipedia article: https://en.wikipedia.org/wiki/FASTQ_format#Encoding

.. highlight:: Bash

#################
 Getting started
#################

To run AdapterRemoval on single-end FASTQ data, simply specify the location of FASTQ file(s) using the ``--in-file1`` command-line options:

.. code::

   adapterremoval3 --in-file1 myreads_1.fastq.gz

To run AdapterRemoval on paired-end FASTQ data, specify the location of the mate 1 and mate 2 FASTQ files using the ``--in-file1`` and ``--in-file2`` command-line options:

.. code::

   adapterremoval3 --in-file1 myreads_1.fastq.gz --in-file2 myreads_2.fastq.gz

The files may be uncompressed or gzip-compressed. When run in this manner, AdapterRemoval will save the trimmed reads in the current working directly, using filenames starting with 'your_output'. This behavior may be changed using the ``--out-prefix`` option, or using specific options for each output file. See the :doc:`input_and_output` page for more information about files generated by AdapterRemoval.

More examples of common usage may be found in the :doc:`examples` section of the documentation.

.. _specifying_adapters:

*******************************
 A note on specifying adapters
*******************************

AdapterRemoval uses the expected adapter sequences as part of the trimming/alignment process (see :doc:`detailed_overview`). It is therefore extremely important to specify the correct adapter sequences when running AdapterRemoval on a dataset that does not make use of these adapters. Failure to do so will result in the wrong sequences being trimmed, and actual adapter sequences being left in the resulting "trimmed" reads.

By default, AdapterRemoval is set up to trim the published `Illumina TruSeq sequences`_, which should be applicable to most Illumina data, corresponding to the following command-line options:

.. code::

   adapterremoval3 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

For BGISEQ/DNBSEQ/MGISEQ data, the published `BGI adapter sequences`_ should be used:

.. code::

   adapterremoval3 --adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG

Adapter sequences are specified in the read orientation when using the ``--adapter1`` and ``--adapter2`` command-line options, directly corresponding to the sequence that is observed in the FASTQ files produced by the base calling software. If we were processing data generated using the above TrueSeq adapters, then we would therefore expect to find those sequences as-is in our FASTQ files (assuming that the read lengths are sufficiently long and that insert sizes are sufficiently short):

.. code::

   $ grep "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" file1.fastq
   AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAACAAGAAT
   CTGGAGTTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGAATCTCGTATGCCGTCTTCTGCTTGAAAAAAA
   GGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGAATCTCGTATGCCGTCTTCTGCTTGCAAATTGAAAACAC

   $ grep "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" file2.fastq
   CAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTCAAAAAAAGAAAAACATCTTG
   GAACTCCAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTCAAAAAAAATAGA
   GAACTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTCAAAAACATAAGACCTA

For paired-end data, the ``--report-only`` mode may be used to verify the choice of adapters, by attempting to reconstruct the adapter sequence directly from the FASTQ reads. See the :doc:`examples` section for a demonstration of this functionality.

An 'N' in an adapter sequence is treated as a wildcard. An N will align against any other base, including Ns, but do not affect the score of the resulting alignment and are not counted as for the purpose of filters such as ``--min-adapter-overlap``.

.. _bgi adapter sequences: https://en.mgitech.cn/Download/download_file/id/71

.. _illumina truseq sequences: https://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html

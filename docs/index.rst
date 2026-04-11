################
 AdapterRemoval
################

AdapterRemoval is a multi-platform tool for processing High-Throughput Sequencing (HTS) data in FASTQ format. AdapterRemoval trims remnant adapter sequences, trims and filters low quality reads, merges overlapping paired-end reads, demultiplexes sequencing reads, and generates QC reports in human and machine-readable formats.

AdapterRemoval is designed for rapidly processing large datasets, and therefore accelerates operations using hardware specific instruction sets (SSE2, AVX, and NEON), via parallel processing of sequencing data, parallel compression of outputs, and other techniques.

- See the :doc:`getting_started` to get started processing your own data.
- See the :doc:`detailed_overview` page for a description of processing steps carried out by AdapterRemoval, and the options affecting each of these steps.
- See the :doc:`input_and_output` page for details on the input files read by AdapterRemoval and resulting output files.
- See the :doc:`migrating` page if you are upgrading from an older version of AdapterRemoval.
- For a complete list of command-line options, see the :doc:`manpage`.

For questions, suggestions, or bug reports, please `create an issue on GitHub`_.

If you use AdapterRemoval v3, then please cite the paper:

   Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 12;9(1):88 http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2

AdapterRemoval was originally published in Lindgreen 2012:

   Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation Sequencing Reads, BMC Research Notes, 5:337 http://www.biomedcentral.com/1756-0500/5/337/

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   getting_started
   detailed_overview
   input_and_output
   migrating
   manpage

####################
 Indices and tables
####################

-  :ref:`genindex`
-  :ref:`modindex`
-  :ref:`search`

.. _create an issue on github: https://github.com/MikkelSchubert/adapterremoval/issues/new

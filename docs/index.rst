AdapterRemoval
==============

AdapterRemoval searches for and removes remnant adapter sequences from High-Throughput Sequencing (HTS) data, trims low quality bases, merges merges overlapping paired-ended reads into (longer) consensus sequences, and generates quality QC reports (human and machine readable). Additionally, AdapterRemoval can construct a consensus adapter sequence from paired-ended reads.

See the :doc:`detailed_overview` page for a description of processing carried out by AdapterRemoval and the options affecting each step of the process.

If you use AdapterRemoval v3, then please cite the paper:

    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 12;9(1):88
    http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2

AdapterRemoval was originally published in Lindgreen 2012:

    Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation Sequencing Reads, BMC Research Notes, 5:337
    http://www.biomedcentral.com/1756-0500/5/337/


.. toctree::
    :maxdepth: 2
    :caption: Contents:

    installation
    getting_started
    examples
    detailed_overview
    input_and_output
    migrating
    manpage


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

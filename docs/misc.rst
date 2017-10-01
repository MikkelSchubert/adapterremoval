Miscellaneous
=============


Window-based quality trimming
-----------------------------

As of v2.2.2, AdapterRemoval implements sliding window based approach to quality based base-trimming inspired by `sickle`_. If ``window_size`` is greater than or equal to 1, that number is used as the window size for all reads. If ``window_size`` is a number greater than or equal to 0 and less than 1, then that number is multiplied by the length of individual reads to determine the window size. If the window length is zero or is greater than the current read length, then the read length is used instead.

Reads are trimmed as follows for a given window size:

  1. The new 5' is determined by locating the first window where both the average quality and the quality of the first base in the window is greater than ``--minquality``.

  2. The new 3' is located by sliding the first window right, until the average quality becomes less than or equal to ``--minquality``. The new 3' is placed at the last base in that window where the quality is greater than or equal to ``--minquality``.

  3. If no 5' position could be determined, the read is discarded.


Migrating from AdapterRemoval v1.x
----------------------------------

Command-line options mostly behave the same between AdapterRemoval v1 and AdapterRemoval v2, and scripts written with AdapterRemoval v1.x in mind should work with AdapterRemoval v2.x. A notable exception is the ``--pcr1`` and ``--pcr2`` options, which have been replaced by the ``--adapter1`` and ``--adapter2`` options described above. While the ``--pcr`` options are still supported for backwards compatibility, these should not be used going forward.

The difference between these two options is that ``--adapter2`` expects the mate 2 adapter sequence to be specified in the read orientation as described above, while the ``--pcr2`` expects the sequence to be in the same orientation as the mate 1 sequence, the reverse complement of the sequence observed in the mate 2 reads. 

Using the common 13 bp Illumina adapter sequence (AGATCGGAAGAGC) as an example, this is how the options would be used in AdapterRemoval v2.x::

	AdapterRemoval --adapter1 AGATCGGAAGAGC --adapter2 AGATCGGAAGAGC ...

And in AdapterRemoval v1.x::

	AdapterRemoval --adapter1 AGATCGGAAGAGC --adapter2 GCTCTTCCGATCT ...


.. _sickle: https://github.com/najoshi/sickle
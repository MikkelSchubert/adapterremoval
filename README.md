==============
AdapterRemoval
==============

This program searches for and removes remnant adapter sequences from your read data and (optionally) trims low quality bases from the 3' end of reads following adapter removal.  The program can analyze both single end and paired end data, and can be used to merge overlapping paired-ended reads into (longer) consensus sequences. Additionally, the the program may to recover a consensus adapter sequence for paired-ended data.

For detailed explanation of the parameters, please refer to the man page.  For comments, suggestions  and feedback please contact Mikkel Schubert (MikkelSch@gmail.com) and Stinus Lindgreen (stinus@binf.ku.dk).

If you use the program, please cite the paper:
    S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation Sequencing Reads, BMC Research Notes, 5:337
    http://www.biomedcentral.com/1756-0500/5/337/


Installation
============

Download and unpack the newest release from GitHub:

$ wget https://github.com/MikkelSchubert/adapterremoval/archive/v2.0.0.tar.gz
$ tar xvzf tar xvzf adapterremoval-2.0.0.tar.gz
$ cd adapterremoval-2.0.0

or

$ git clone https://github.com/MikkelSchubert/adapterremoval.git
$ cd adapterremoval

To compile, run

$ make

The resulting binary and man page is located in the "build" folder.

To install, run

$ sudo make install


Usage
=====

For program usage, please the the manual page. If AdapterRemoval has been installed, this may be accessed using the command "man AdapterRemoval". If AdapterRemoval has not been installed, the manual page may be read using the command "man build/AdapterRemoval.1" in the source folder once "make" has been run.